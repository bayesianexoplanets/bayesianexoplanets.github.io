"""merge_catalog.py — merge the full CatalogRun re-vetting into the website catalog.

For every (TIC, TOI) the CatalogRun produced, replace the recomputed columns in
``TESS corrected/tois.csv``; keep the per-star NST (μ/σ/nst_samples) and stellar
(Mass/Radius/logg/FEH/Teff) + ``Has Visible TTVs`` columns; recompute
``log10(p value)`` from the new SNR and pass/fail with the standard rule; refresh
each planet's plot (``plots/{TIC}/{known_idx}.jpg`` ← ``CatalogRun/plots/{TIC}/{plot_index}.png``).
Rows with no CatalogRun match keep their current values (reported).

Conventions (verified): Duration=2·Tau, Epoch=t_start+Phase,
log10(p)=norm.logsf((SNR−μ)/σ)/ln(10), Radius_planet_err{p,m}=radius{p,m}.
TOI and plot_index come straight from the CatalogRun CSV (no read_known re-derivation).

Dry-run by default; pass --apply to write tois.csv + plots.
"""
import os
import sys
import glob

os.chdir("/global/u2/j/julius")
import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = "/global/u2/j/julius/TESS corrected"
TOIS = os.path.join(HERE, "tois.csv")
RUN = "/pscratch/sd/j/julius/exoprob/CatalogRun"
CAND = RUN + "/known_candidates"
PSRC = RUN + "/plots"
PDST = os.path.join(HERE, "plots")
STARS = "/pscratch/sd/j/julius/exoprob/FullRun/stars"
BULK = "/pscratch/sd/j/julius/Bulk Download/Data"
LN10 = np.log(10.0)
SPURIOUS_MAX, SNRD_MIN, SNR_OVERRIDE, NTRANSITS_MIN, HARMONICS_SNRD = 0.75, 1e-2, 50.0, 3, 0.95


def sharp(tic, _c={}):
    if tic not in _c:
        try:
            _c[tic] = bool(pd.read_csv(f"{STARS}/{tic}.csv", sep="\t")["has_sharp_peak"].iloc[0])
        except Exception:
            _c[tic] = False
    return _c[tic]


def t_start(tic, _c={}):
    if tic not in _c:
        _c[tic] = float(np.load(f"{BULK}/StellarData_{tic}.npy")[0])
    return _c[tic]


def failures(sp1, snrd, snr, ntr, sh):
    f = []
    if sp1 >= SPURIOUS_MAX:
        f.append("spurious")
    if snrd <= SNRD_MIN and snr <= SNR_OVERRIDE:
        f.append("snrd")
    if ntr < NTRANSITS_MIN:
        f.append("ntransits")
    if sh and snrd > HARMONICS_SNRD and snr <= SNR_OVERRIDE:
        f.append("sharp_freq")
    return f


def main(apply):
    w = pd.read_csv(TOIS)
    run = {}                      # (TIC, TOI) -> CatalogRun row
    for f in glob.glob(CAND + "/*.csv"):
        tic = int(os.path.basename(f).split(".")[0])
        try:
            d = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        for _, r in d.iterrows():
            run[(tic, str(r["TOI"]))] = r
    print(f"website rows={len(w)} | CatalogRun rows={len(run)} | TICs with CSV={len({k[0] for k in run})}")

    updated = plots = 0
    unmatched = []
    flips = {"pass->fail": 0, "fail->pass": 0}

    for tic, sub in w.groupby("TIC"):
        for known_idx, (wi, wrow) in enumerate(sub.iterrows()):
            toi = str(wrow["TOI"])
            r = run.get((int(tic), toi))
            if r is None:
                unmatched.append((int(tic), toi))
                continue
            snr = float(r["SNR"]); snrd = float(r["snrd_pvalue"])
            sp1 = float(r["spurious1"]); ntr = int(r["num_available_transits"])
            ph = float(r["phase"]); tauf = float(r["tau"])
            fl = failures(sp1, snrd, snr, ntr, sharp(int(tic)))
            passed = len(fl) == 0
            mu = float(wrow["μ(SNR | null)"]); sig = float(wrow["σ(SNR | null)"])
            log10p = float(norm.logsf((snr - mu) / sig) / LN10) if sig > 0 else wrow["log10(p value)"]

            was = bool(wrow["passed_all_tests"])
            if was and not passed:
                flips["pass->fail"] += 1
            elif passed and not was:
                flips["fail->pass"] += 1

            upd = {
                "Period": float(r["period"]), "Phase": ph, "Tau": tauf, "Duration": 2 * tauf,
                "Epoch": t_start(int(tic)) + ph, "SNR": snr,
                "err_Period": float(r["err_period"]), "err_Epoch": float(r["err_phase"]),
                "err_Duration": 2 * float(r["err_tau"]),
                "Radius_planet": float(r["radius"]), "Radius_planet_errp": float(r["radiusp"]),
                "Radius_planet_errm": float(r["radiusm"]), "Number of Valid Transits": ntr,
                "log10(p value)": log10p, "passed_all_tests": bool(passed), "failed_tests": "|".join(fl),
            }
            if apply:
                for k, v in upd.items():
                    w.at[wi, k] = v
            updated += 1

            src = os.path.join(PSRC, str(int(tic)), f"{int(r['plot_index'])}.png")
            dst = os.path.join(PDST, str(int(tic)), f"{known_idx}.jpg")
            if os.path.exists(src):
                if apply:
                    _to_jpg(src, dst)
                plots += 1
            else:
                if apply:
                    print(f"  WARN missing plot {src}")

    print(f"rows updated={updated} | plots refreshed={plots} | unmatched(kept)={len(unmatched)}")
    print(f"pass/fail flips vs current catalog: {flips}")
    if unmatched:
        print("sample unmatched (kept as-is):", unmatched[:15])
    if apply:
        w.to_csv(TOIS, index=False)
        print(f"WROTE {TOIS} ({len(w)} rows)")
    else:
        print("DRY RUN — pass --apply to write tois.csv + plots")


def _to_jpg(src, dst, max_w=900, quality=88):
    from PIL import Image
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    im = Image.open(src).convert("RGB")
    if im.width > max_w:
        im = im.resize((max_w, int(im.height * max_w / im.width)), Image.LANCZOS)
    im.save(dst, "JPEG", quality=quality)


if __name__ == "__main__":
    main(apply="--apply" in sys.argv)
