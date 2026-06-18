"""Phase 4-prep: for every validated catalog planet, find the RecoveryRun candidate that
matches the planet (period ratio ~1, NOT a harmonic, AND phase). That candidate's
(period, phase, tau) is the fiducial seed for the re-vet; planets with no clean match fall
back to an ExoFOP seed + grid search. Reports coverage and writes the seed table.

Reference (period, phase): ExoFOP if available, else the current catalog value.
  phase_ref = mod(Epoch(BJD) - (t_start+2457000), P_ref)   [ExoFOP], else catalog Phase.
Match: |P_cand/P_ref - 1| < PTOL  AND  circular |phase_cand - phase_ref| mod P_ref < phase_tol.
"""
import os
import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/TESS corrected"
TOIS = os.path.join(HERE, "tois.csv")
EXO = "/global/u2/j/julius/TESS/tois.csv"
RR = "/pscratch/sd/j/julius/exoprob/RecoveryRun/candidates/batch0"
BULK = "/pscratch/sd/j/julius/Bulk Download/Data"
PTOL = 0.01            # period within 1% (ratio ~ 1, rejects harmonics)

w = pd.read_csv(TOIS)
E = pd.read_csv(EXO)
_ts = {}


def t_start(tic):
    if tic not in _ts:
        try:
            _ts[tic] = float(np.load(f"{BULK}/StellarData_{tic}.npy")[0])
        except Exception:
            _ts[tic] = np.nan
    return _ts[tic]


def exo_ref(tic, toi):
    r = E[(E["TIC ID"] == tic) & (E["TOI"].astype(str) == str(toi))]
    if not len(r):
        return np.nan, np.nan
    P = float(r.iloc[0]["Period (days)"]); ep = float(r.iloc[0]["Epoch (BJD)"])
    if not np.isfinite(P) or P <= 0:
        return np.nan, np.nan
    ph = np.mod(ep - (t_start(tic) + 2457000.0), P) if np.isfinite(ep) and np.isfinite(t_start(tic)) else np.nan
    return P, ph


def circ_d(a, b, P):
    d = np.mod(a - b, P)
    return np.minimum(d, P - d)


rows = []
for tic, sub in w.groupby("TIC"):
    tic = int(tic)
    p = f"{RR}/{tic}.0.csv"
    rr = pd.read_csv(p, sep="\t") if os.path.exists(p) else None
    for _, r in sub.iterrows():
        toi = str(r["TOI"])
        P_ref, ph_ref = exo_ref(tic, toi)
        src = "exofop"
        if not np.isfinite(P_ref):                       # no ExoFOP -> use catalog as reference
            P_ref, ph_ref, src = float(r["Period"]), float(r["Phase"]), "catalog"
        match = None
        if rr is not None and np.isfinite(P_ref):
            cp = rr["period"].to_numpy(float)
            ok = np.abs(cp / P_ref - 1.0) < PTOL
            if np.isfinite(ph_ref):
                tau = rr["tau"].to_numpy(float)
                ptol = np.maximum(3 * tau, 0.03 * P_ref)
                ok &= circ_d(rr["phase"].to_numpy(float), ph_ref, P_ref) < ptol
            idx = np.where(ok)[0]
            if len(idx):
                match = rr.iloc[idx[int(np.argmax(rr.iloc[idx]["SNR"].to_numpy()))]]
        rows.append(dict(TIC=tic, TOI=toi, ref=src, P_ref=P_ref, ph_ref=ph_ref,
                         matched=match is not None,
                         rr_period=float(match["period"]) if match is not None else np.nan,
                         rr_phase=float(match["phase"]) if match is not None else np.nan,
                         rr_tau=float(match["tau"]) if match is not None else np.nan,
                         rr_SNR=float(match["SNR"]) if match is not None else np.nan,
                         rr_nVT=int(match["num_available_transits"]) if match is not None else -1))

m = pd.DataFrame(rows)
m.to_csv(os.path.join(HERE, "recovery_match.csv"), index=False)
print(f"planets: {len(m)} | matched RecoveryRun: {int(m.matched.sum())} "
      f"({100*m.matched.mean():.1f}%) | need grid-search: {int((~m.matched).sum())}")
print(f"  no RecoveryRun file: {int(sum(1 for t in m.TIC.unique() if not os.path.exists(f'{RR}/{t}.0.csv')))} TICs")
print("\nflagged TICs (should mostly match):")
for tic in [127530399, 333607525, 154618248, 156724719, 255704097, 38603673, 374829238, 298663873, 441420236, 29191596]:
    s = m[m.TIC == tic]
    for _, r in s.iterrows():
        print(f"  {tic}/{r['TOI']}: matched={r['matched']} rr(P={r['rr_period']:.4f} ph={r['rr_phase']:.3f} "
              f"tau={r['rr_tau']:.4f} SNR={r['rr_SNR']:.1f})" if r['matched'] else f"  {tic}/{r['TOI']}: NO MATCH -> grid search")
print("\nwrote recovery_match.csv")
