"""Fix-ups:
(i)   fill new-candidate stellar props from the known catalog (same TIC), Gaia fallback;
(ii)  replace non-finite (inf/NaN) new-candidate P/phase/tau (and radius) errors with the
      closest finite source on that star;
(iii) report Epoch in BTJD (BJD - 2457000) for both catalogs.
"""
import os, glob
import numpy as np
import pandas as pd

BASE = "/global/u2/j/julius/TESS corrected"
FULL = "/pscratch/sd/j/julius/exoprob/FullRun"
RECOV = "/pscratch/sd/j/julius/exoprob/RecoveryRun"

known = pd.read_csv(f"{BASE}/tois.csv")
new = pd.read_csv(f"{BASE}/tois_new.csv")

# ---------- (iii) Epoch -> BTJD ----------
known["Epoch"] = known["Epoch"] - 2457000.0
new["Epoch"] = new["Epoch"] - 2457000.0

# ---------- (i) stellar props for new candidates ----------
SCOLS = ["Mass", "Radius", "logg", "FEH", "Teff"]
host = known.groupby("TIC")[SCOLS].mean()
for c in SCOLS:
    new[c] = new[c].fillna(new["TIC"].map(host[c]))

# Gaia fallback for TICs not in the known catalog
gaia = pd.read_csv("/global/u2/j/julius/TESS/Gaia_TESS.csv").set_index("id_starname")
gmap = {"Mass": "iso_mass", "Radius": "iso_rad", "logg": "iso_logg",
        "FEH": "iso_feh", "Teff": "iso_teff"}
for i, r in new.iterrows():
    key = "tic" + str(int(r["TIC"]))
    if key in gaia.index:
        for c, gc in gmap.items():
            if pd.isna(new.at[i, c]):
                v = gaia.loc[key, gc]
                if np.isscalar(v) and pd.notna(v):
                    new.at[i, c] = float(v)
print("stellar fill remaining-empty:",
      {c: int(new[c].isna().sum()) for c in SCOLS})

# ---------- (ii) closest finite error fallback ----------
_cand = {}
def cand_csv(run, tic):
    k = (run, tic)
    if k not in _cand:
        f = os.path.join(run, "candidates", "batch0", f"{tic}.0.csv")
        try:
            _cand[k] = pd.read_csv(f, sep="\t") if os.path.exists(f) else None
        except Exception:
            _cand[k] = None
    return _cand[k]

def fin(*xs):
    return all(np.isfinite(x) for x in xs)

def closest_finite(tic, P, ph, tau):
    """nearest source whose P/phase/tau errors are all finite -> (eP,eph,etau,radp,radm,kind)"""
    srcs = []
    for path in sorted(glob.glob(os.path.join(FULL, "known_planets", f"{tic}_*.npz"))):
        try:
            z = np.load(path, allow_pickle=True); p = z["params"]; cov = np.asarray(z["cov"])
        except Exception:
            continue
        if cov.shape == (6,):
            eP, eph, et = np.sqrt(abs(cov[0])), np.sqrt(abs(cov[1])), np.sqrt(abs(cov[2]))
            if fin(eP, eph, et):
                srcs.append((float(p[0]), float(p[1]), float(p[2]), eP, eph, et, np.nan, np.nan, "npz"))
    for run in (FULL, RECOV):
        df = cand_csv(run, tic)
        if df is not None and "period" in df.columns:
            for _, c in df.iterrows():
                eP, eph, et = float(c["err_period"]), float(c["err_phase"]), float(c["err_tau"])
                if fin(eP, eph, et):
                    srcs.append((float(c["period"]), float(c["phase"]), float(c["tau"]),
                                 eP, eph, et, float(c["radiusp"]), float(c["radiusm"]), "cand"))
    if not srcs:
        return None
    s = min(srcs, key=lambda s: abs(s[0]-P)/max(P,1e-6) + abs(s[1]-ph)/max(P,1e-6)
            + abs(s[2]-tau)/max(tau,1e-3))
    return s[3], s[4], s[5], s[6], s[7], s[8]

def nonfinite(*xs):
    return any((pd.isna(x) or not np.isfinite(x)) for x in xs)

fixed = 0
for i, r in new.iterrows():
    bad_t = nonfinite(r["err_Period"], r["err_Epoch"], r["err_Duration"])
    bad_r = nonfinite(r["Radius_planet_errp"], r["Radius_planet_errm"])
    if not (bad_t or bad_r):
        continue
    cl = closest_finite(int(r["TIC"]), r["Period"], r["Phase"], r["Tau"])
    if cl is None:
        continue
    eP, eph, et, radp, radm, kind = cl
    if bad_t:
        new.at[i, "err_Period"] = eP
        new.at[i, "err_Epoch"] = eph
        new.at[i, "err_Duration"] = 2.0 * et
    if bad_r and kind == "cand" and fin(radp, radm):
        new.at[i, "Radius_planet_errp"] = radp
        new.at[i, "Radius_planet_errm"] = radm
    fixed += 1
print("new candidates with non-finite errors fixed via closest source:", fixed)

# also clean any non-finite errors left in the KNOWN catalog (write empty so JS shows none)
for col in ["err_Period", "err_Epoch", "err_Duration", "Radius_planet_errp", "Radius_planet_errm"]:
    n_inf = int((~np.isfinite(known[col].to_numpy(dtype=float))).sum() - known[col].isna().sum())
    known.loc[~np.isfinite(known[col].to_numpy(dtype=float)), col] = np.nan
    if n_inf:
        print(f"known {col}: {n_inf} non-finite -> blanked")

known.to_csv(f"{BASE}/tois.csv", index=False)
new.to_csv(f"{BASE}/tois_new.csv", index=False)

# spot checks
print("\n231077395:",
      new[new.TIC == 231077395][["Period", "Epoch", "err_Period", "err_Epoch",
                                  "err_Duration", "Radius_planet_errp"]].to_dict("records"))
print("38584799 stellar:",
      new[new.TIC == 38584799][["Mass", "Radius", "logg", "FEH", "Teff", "Epoch"]].to_dict("records"))
print("known Epoch range:", known["Epoch"].min(), "to", known["Epoch"].max())
