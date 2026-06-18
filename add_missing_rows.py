"""Add the 7 recoverable missing ExoFOP siblings to the website catalog (with ExoFOP
ephemerides) and reset the 2 period-mismatch rows to their ExoFOP periods, so a
--fix-period re-vet picks them up. Re-vet + merge fills SNR/radius/pass/plots afterward."""
import os, sys, shutil
os.chdir("/global/u2/j/julius")
import numpy as np, pandas as pd

HERE = "/global/u2/j/julius/TESS corrected"
TOIS = os.path.join(HERE, "tois.csv")
EXO = "/global/u2/j/julius/TESS/tois.csv"
BULK = "/pscratch/sd/j/julius/Bulk Download/Data"

NEW = [(32090583, "218.01"), (31852980, "487.02"), (64837857, "6650.02"),
       (342449055, "6704.01"), (234388232, "5088.02"), (119292328, "512.02"),
       (427320001, "2112.02")]
RESET = [(258514800, "1444.01"), (324609476, "5041.01")]

w = pd.read_csv(TOIS)
e = pd.read_csv(EXO)
RAD_COL = next((c for c in ["Planet Radius (R_Earth)", "Planet Radius (Rearth)"] if c in e.columns), None)


def exo(tic, toi):
    r = e[(e["TIC ID"] == tic) & (e["TOI"].astype(str) == toi)]
    return r.iloc[0] if len(r) else None


def ephem(tic, r):
    t_start = float(np.load(f"{BULK}/StellarData_{tic}.npy")[0])
    P = float(r["Period (days)"])
    phase = float(np.mod(float(r["Epoch (BJD)"]) - (t_start + 2457000.0), P))
    tau = float(r["Duration (hours)"]) * 0.5 / 24.0
    rp = float(r[RAD_COL]) if (RAD_COL and np.isfinite(r[RAD_COL])) else float("nan")
    return t_start, P, phase, tau, rp


shutil.copy(TOIS, TOIS + ".bak_before_missing")
print(f"backup -> {TOIS}.bak_before_missing  ({len(w)} rows)")

# --- reset the 2 mismatches to ExoFOP (so fix_period holds the correct period) ---
for tic, toi in RESET:
    r = exo(tic, toi)
    t_start, P, phase, tau, rp = ephem(tic, r)
    mask = (w["TIC"] == tic) & (w["TOI"].astype(str) == toi)
    if mask.sum() != 1:
        print(f"  RESET WARN {tic}/{toi}: {mask.sum()} rows"); continue
    i = w.index[mask][0]
    w.at[i, "Period"], w.at[i, "Phase"], w.at[i, "Tau"] = P, phase, tau
    w.at[i, "Epoch"], w.at[i, "Duration"] = t_start + phase, 2 * tau
    if np.isfinite(rp):
        w.at[i, "Radius_planet"] = rp
    print(f"  reset {tic}/{toi}: P->{P:.5f} phase->{phase:.4f} tau->{tau:.5f}")

# --- add the 7 new siblings (copy a same-TIC sibling row, override planet cols) ---
rows = []
for tic, toi in NEW:
    sib = w[w["TIC"] == tic]
    if not len(sib):
        print(f"  ADD WARN {tic}/{toi}: no sibling row"); continue
    r = exo(tic, toi)
    t_start, P, phase, tau, rp = ephem(tic, r)
    new = sib.iloc[0].copy()
    new["TOI"], new["Period"], new["Phase"], new["Tau"] = toi, P, phase, tau
    new["Epoch"], new["Duration"] = t_start + phase, 2 * tau
    new["Radius_planet"] = rp if np.isfinite(rp) else new["Radius_planet"]
    for c in ["SNR", "log10(p value)", "err_Period", "err_Epoch", "err_Duration",
              "Radius_planet_errp", "Radius_planet_errm", "Number of Valid Transits"]:
        new[c] = np.nan
    new["passed_all_tests"], new["failed_tests"] = False, "pending"
    rows.append(new)
    print(f"  add  {tic}/{toi}: P={P:.5f} phase={phase:.4f} tau={tau:.5f} rp={new['Radius_planet']}")

if rows:
    w = pd.concat([w, pd.DataFrame(rows)], ignore_index=True)
if "--apply" in sys.argv:
    w.to_csv(TOIS, index=False)
    print(f"WROTE {TOIS} ({len(w)} rows)")
else:
    print(f"DRY RUN — would write {len(w)} rows. pass --apply")
