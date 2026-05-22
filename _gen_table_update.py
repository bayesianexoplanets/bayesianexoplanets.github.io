"""One-time generator: add Epoch/Duration + error bars to tois.csv & tois_new.csv,
and refresh NST mu/sigma/log10(p) for previously-empty rows.

Locked methods (verified in _calib.py):
  - known_planets NPZ cov = [Var(P), Var(phase), Var(tau), Cov(P,ph), Cov(ph,tau), Cov(P,tau)]
    -> err = sqrt(|cov[0..2]|).  (ttv.laplace, pipeline/ttv.py:96)
  - candidate CSV (tab-sep) has err_period/err_phase/err_tau + radiusp/radiusm.
  - Epoch = phase + t_start + 2457000,  t_start = StellarData_{TIC}.npy[0]  (BJD).
  - Duration = 2*tau ; err_Duration = 2*err_tau ; err_Epoch = err_phase.
  - radius err: new cands -> radiusp/radiusm ; known -> ExoFOP 'Planet Radius (R_Earth) err'.
  - NST: per-run max SNR (mask SNR==0), mu=mean, sigma=std(ddof=0),
    FullRun batch1..10 ; fall back to BigRun batch1..100 if no FullRun NST.
  - log10(p) = norm.logsf((SNR-mu)/sigma)/ln10.
"""
import os, glob
import numpy as np
import pandas as pd
from scipy.stats import norm

BASE = "/global/u2/j/julius/TESS corrected"
FULL = "/pscratch/sd/j/julius/exoprob/FullRun"
RECOV = "/pscratch/sd/j/julius/exoprob/RecoveryRun"
BIG = "/pscratch/sd/j/julius/exoprob/BigRun"
BULK = "/pscratch/sd/j/julius/Bulk Download/Data"
LN10 = np.log(10.0)
NEWCOLS = ["Epoch", "Duration", "err_Period", "err_Epoch", "err_Duration",
           "Radius_planet_errp", "Radius_planet_errm"]

# ---------- caches ----------
_tstart = {}
def t_start(tic):
    if tic not in _tstart:
        f = os.path.join(BULK, f"StellarData_{tic}.npy")
        _tstart[tic] = float(np.load(f)[0]) if os.path.exists(f) else np.nan
    return _tstart[tic]

_cand = {}
def cand_csv(run, tic):
    key = (run, tic)
    if key not in _cand:
        f = os.path.join(run, "candidates", "batch0", f"{tic}.0.csv")
        try:
            _cand[key] = pd.read_csv(f, sep="\t") if os.path.exists(f) else None
        except Exception:
            _cand[key] = None
    return _cand[key]

def npz_errors(tic, P, ph, tau, known_idx):
    """Return (err_P, err_ph, err_tau) from a known_planets NPZ whose params match, else None."""
    cands = []
    same = os.path.join(FULL, "known_planets", f"{tic}_{known_idx}.npz")
    if os.path.exists(same):
        cands.append(same)
    cands += [p for p in sorted(glob.glob(os.path.join(FULL, "known_planets", f"{tic}_*.npz"))) if p not in cands]
    for path in cands:
        try:
            z = np.load(path, allow_pickle=True)
            p = z["params"]; cov = np.asarray(z["cov"])
        except Exception:
            continue
        if not (np.isclose(p[0], P, rtol=1e-6, atol=1e-9) and
                np.isclose(p[1], ph, rtol=1e-6, atol=1e-8) and
                np.isclose(p[2], tau, rtol=1e-6, atol=1e-8)):
            continue
        if cov.shape == (6,):
            return (np.sqrt(abs(cov[0])), np.sqrt(abs(cov[1])), np.sqrt(abs(cov[2])))
        return None  # TTV cov -> can't use
    return None

def cand_errors(tic, P, ph, tau):
    """Return (err_P, err_ph, err_tau, radiusp, radiusm) from a matching candidate row, else None."""
    for run in (FULL, RECOV):
        df = cand_csv(run, tic)
        if df is None or "period" not in df.columns:
            continue
        m = df[np.isclose(df["period"], P, rtol=1e-4, atol=1e-6) &
               np.isclose(df["phase"], ph, rtol=1e-4, atol=1e-4) &
               np.isclose(df["tau"], tau, rtol=1e-4, atol=1e-5)]
        if len(m):
            r = m.iloc[0]
            return (float(r["err_period"]), float(r["err_phase"]), float(r["err_tau"]),
                    float(r["radiusp"]), float(r["radiusm"]))
    return None

def closest_errors(tic, P, ph, tau):
    """Fallback for rows with no exact match: use the error of the closest source
    (NPZ or candidate) by normalized P/phase/tau distance.
    Returns (err_P, err_ph, err_tau, radiusp, radiusm, kind) or None."""
    srcs = []
    for path in sorted(glob.glob(os.path.join(FULL, "known_planets", f"{tic}_*.npz"))):
        try:
            z = np.load(path, allow_pickle=True); p = z["params"]; cov = np.asarray(z["cov"])
        except Exception:
            continue
        if cov.shape == (6,):
            srcs.append((float(p[0]), float(p[1]), float(p[2]),
                         np.sqrt(abs(cov[0])), np.sqrt(abs(cov[1])), np.sqrt(abs(cov[2])),
                         None, None, "npz"))
    for run in (FULL, RECOV):
        df = cand_csv(run, tic)
        if df is not None and "period" in df.columns:
            for _, c in df.iterrows():
                srcs.append((float(c["period"]), float(c["phase"]), float(c["tau"]),
                             float(c["err_period"]), float(c["err_phase"]), float(c["err_tau"]),
                             float(c["radiusp"]), float(c["radiusm"]), "cand"))
    if not srcs:
        return None
    def dist(s):
        return (abs(s[0] - P) / max(P, 1e-6) + abs(s[1] - ph) / max(P, 1e-6)
                + abs(s[2] - tau) / max(tau, 1e-3))
    s = min(srcs, key=dist)
    return (s[3], s[4], s[5], s[6], s[7], s[8])

def nst_stats(tic):
    """per-run max SNR (mask SNR==0), mu=mean, sigma=std(ddof=0). FullRun batch1..10 else BigRun."""
    for run, nb in ((FULL, 11), (BIG, 101)):
        maxes = []
        for b in range(1, nb):
            f = os.path.join(run, "candidates", f"batch{b}", f"{tic}.0.csv")
            if os.path.exists(f):
                try:
                    df = pd.read_csv(f, sep="\t")
                except Exception:
                    continue
                if "SNR" in df.columns and len(df):
                    s = df["SNR"][df["SNR"] != 0]
                    if len(s):
                        maxes.append(float(s.max()))
        if maxes:
            a = np.array(maxes)
            return a.mean(), a.std(ddof=0), len(a)
    return None

# ---------- ExoFOP radius err ----------
exo = pd.read_csv("/global/u2/j/julius/TESS/tois.csv")
exo_rerr = {}
for _, e in exo[["TOI", "Planet Radius (R_Earth) err"]].dropna().iterrows():
    exo_rerr[round(float(e["TOI"]), 2)] = float(e["Planet Radius (R_Earth) err"])

# ======================= tois.csv (known planets) =======================
tois = pd.read_csv(os.path.join(BASE, "tois.csv"))
tois["_kidx"] = tois.groupby("TIC").cumcount()
mucol, sigcol, lpcol = "μ(SNR | null)", "σ(SNR | null)", "log10(p value)"

src_counts = {"npz": 0, "cand": 0, "closest": 0, "none": 0}
errP = np.full(len(tois), np.nan); errPh = np.full(len(tois), np.nan); errT = np.full(len(tois), np.nan)
rerrp = np.full(len(tois), np.nan); rerrm = np.full(len(tois), np.nan)
epoch = np.full(len(tois), np.nan); dur = np.full(len(tois), np.nan)

nst_filled = 0; nst_still_missing = 0; nst_missing_tics = []

for i, r in tois.iterrows():
    tic = int(r["TIC"]); P = r["Period"]; ph = r["Phase"]; tau = r["Tau"]
    # period/phase/tau errors via source cascade
    e = npz_errors(tic, P, ph, tau, int(r["_kidx"]))
    if e is not None:
        errP[i], errPh[i], errT[i] = e; src_counts["npz"] += 1
        rk = exo_rerr.get(round(float(r["TOI"]), 2)) if pd.notna(r["TOI"]) else None
        if rk is not None:
            rerrp[i] = rerrm[i] = rk
    else:
        c = cand_errors(tic, P, ph, tau)
        if c is not None:
            errP[i], errPh[i], errT[i], rerrp[i], rerrm[i] = c; src_counts["cand"] += 1
        else:
            cl = closest_errors(tic, P, ph, tau)
            rk = exo_rerr.get(round(float(r["TOI"]), 2)) if pd.notna(r["TOI"]) else None
            if cl is not None:
                errP[i], errPh[i], errT[i], rp, rm, kind = cl
                src_counts["closest"] += 1
                if kind == "cand" and rp is not None and np.isfinite(rp):
                    rerrp[i], rerrm[i] = rp, rm
                elif rk is not None:
                    rerrp[i] = rerrm[i] = rk
            else:
                src_counts["none"] += 1
                if rk is not None:
                    rerrp[i] = rerrm[i] = rk
    # Epoch / Duration
    ts = t_start(tic)
    if pd.notna(ts) and pd.notna(ph):
        epoch[i] = ph + ts + 2457000.0
    if pd.notna(tau):
        dur[i] = 2.0 * tau
    # NST refresh for previously-empty rows
    if pd.isna(r[mucol]):
        st = nst_stats(tic)
        if st is not None and st[1] > 0:
            mu, sig, n = st
            tois.at[i, mucol] = mu; tois.at[i, sigcol] = sig
            tois.at[i, lpcol] = norm.logsf((r["SNR"] - mu) / sig) / LN10
            nst_filled += 1
        else:
            nst_still_missing += 1; nst_missing_tics.append(tic)

tois["Epoch"] = epoch
tois["Duration"] = dur
tois["err_Period"] = errP
tois["err_Epoch"] = errPh
tois["err_Duration"] = 2.0 * errT
tois["Radius_planet_errp"] = rerrp
tois["Radius_planet_errm"] = rerrm
tois = tois.drop(columns=["_kidx"])
tois.to_csv(os.path.join(BASE, "tois.csv"), index=False)

print("=== tois.csv ===")
print("error source:", src_counts, "(of", len(tois), "rows)")
print(f"err_Period filled : {np.isfinite(errP).sum()}")
print(f"radius err filled : {np.isfinite(rerrp).sum()}")
print(f"Epoch filled      : {np.isfinite(epoch).sum()}  (missing StellarData: {np.isnan(epoch).sum()})")
print(f"NST newly filled  : {nst_filled}")
print(f"NST still missing : {nst_still_missing}  e.g. {nst_missing_tics[:15]}")

# ======================= tois_new.csv (new candidates) =======================
new = pd.read_csv(os.path.join(BASE, "tois_new.csv"))
nerrP = np.full(len(new), np.nan); nerrPh = np.full(len(new), np.nan); nerrT = np.full(len(new), np.nan)
nrerrp = np.full(len(new), np.nan); nrerrm = np.full(len(new), np.nan)
nepoch = np.full(len(new), np.nan); ndur = np.full(len(new), np.nan)
matched = 0
for i, r in new.iterrows():
    tic = int(r["TIC"]); cidx = int(r["_cand_idx"])
    df = cand_csv(FULL, tic)
    if df is not None and "event_id" in df.columns:
        m = df[df["event_id"] == cidx]
        if len(m):
            c = m.iloc[0]
            nerrP[i] = float(c["err_period"]); nerrPh[i] = float(c["err_phase"]); nerrT[i] = float(c["err_tau"])
            nrerrp[i] = float(c["radiusp"]); nrerrm[i] = float(c["radiusm"]); matched += 1
    ts = t_start(tic)
    if pd.notna(ts) and pd.notna(r["Phase"]):
        nepoch[i] = r["Phase"] + ts + 2457000.0
    if pd.notna(r["Tau"]):
        ndur[i] = 2.0 * r["Tau"]

new["Epoch"] = nepoch
new["Duration"] = ndur
new["err_Period"] = nerrP
new["err_Epoch"] = nerrPh
new["err_Duration"] = 2.0 * nerrT
new["Radius_planet_errp"] = nrerrp
new["Radius_planet_errm"] = nrerrm
new.to_csv(os.path.join(BASE, "tois_new.csv"), index=False)

print("\n=== tois_new.csv ===")
print(f"rows: {len(new)}  candidate-matched (errors+radius): {matched}")
print(f"Epoch filled: {np.isfinite(nepoch).sum()}  (missing StellarData: {np.isnan(nepoch).sum()})")

# spot checks
print("\n=== spot-check known rows ===")
for tic in [149603524, 69679391]:
    rr = tois[tois.TIC == tic].iloc[0]
    print(f"TIC {tic}: P={rr['Period']:.6f}+/-{rr['err_Period']:.2e}  "
          f"Epoch={rr['Epoch']:.4f}+/-{rr['err_Epoch']:.2e}  "
          f"Dur={rr['Duration']:.5f}+/-{rr['err_Duration']:.2e}  "
          f"Rp={rr['Radius_planet']}+/-{rr['Radius_planet_errp']}")
print("=== spot-check new candidates ===")
for tic in new["TIC"].head(2).astype(int):
    rr = new[new.TIC == tic].iloc[0]
    print(f"TIC {tic}: P={rr['Period']:.5f}+/-{rr['err_Period']:.2e}  "
          f"Epoch={rr['Epoch']:.4f}+/-{rr['err_Epoch']:.2e}  "
          f"Dur={rr['Duration']:.5f}+/-{rr['err_Duration']:.2e}  "
          f"Rp={rr['Radius_planet']} +{rr['Radius_planet_errp']}/-{rr['Radius_planet_errm']}")
