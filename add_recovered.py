"""add_recovered.py — add rerun-recovered planets missing from the website catalog.

Some rerun stars have planets in ExoFOP (and recovered by the rerun) that were
dropped when the catalog was first built from tois_corrected.csv. apply_rerun.py
only *updates* existing rows; this script *adds* the missing ones as new rows.

A new row reuses the per-star columns (μ/σ/nst_samples, Mass/Radius/logg/FEH/Teff,
Has Visible TTVs) from an existing sibling row of the same TIC, fills the
planet-specific columns from KnownRun/known_candidates/{tic}.csv, recomputes
log10(p) and passed_all_tests/failed_tests, inserts the row right after the TIC's
existing rows, and copies the planet's Stage-A plot to plots/{TIC}/{known_idx}.jpg.

Dry-run by default; pass --apply to write tois.csv + plots.
"""
import os
import sys
import glob

os.chdir("/global/u2/j/julius")
sys.path.insert(0, "/global/u2/j/julius")
import numpy as np
import pandas as pd
sys.path.insert(0, "/global/u2/j/julius/TESS corrected")
import sm_pvalue
from scipy.stats import norm
from TESS.load_tess import StarInfo_tess, read_known_planets_tess

HERE = "/global/u2/j/julius/TESS corrected"
TOIS = os.path.join(HERE, "tois.csv")
KNOWN = "/pscratch/sd/j/julius/exoprob/KnownRun/known_candidates"
PLOTS_SRC = "/pscratch/sd/j/julius/exoprob/KnownRun/plots"
PLOTS_DST = os.path.join(HERE, "plots")
STARS = "/pscratch/sd/j/julius/exoprob/FullRun/stars"
BULK = "/pscratch/sd/j/julius/Bulk Download/Data"
LN10 = np.log(10.0)
SPURIOUS_MAX, SNRD_MIN, SNR_OVERRIDE, NTRANSITS_MIN, HARMONICS_SNRD = 0.75, 1e-2, 50.0, 3, 0.95


def has_sharp_peak(tic):
    try:
        d = pd.read_csv(os.path.join(STARS, f"{tic}.csv"), sep="\t")
        return bool(d["has_sharp_peak"].iloc[0]) if not d.empty else False
    except Exception:
        return False


def planet_order(tic):
    """[(plot_index, period, TOI)] in read_known order, plus t_start."""
    sd = np.load(os.path.join(BULK, f"StellarData_{tic}.npy"))
    sdc = sd.copy()
    if np.isnan(sdc[3]):
        sdc[3] = 1
    if np.isnan(sdc[-1]):
        sdc[-1] = 1
        sdc[-2] = 1
    star = StarInfo_tess(*sdc[1:])
    pp = read_known_planets_tess(star, sd[0], False)
    e = pd.read_csv("TESS/tois.csv")
    e = e[e["TIC ID"] == tic]
    e = e[e["TFOPWG Disposition"].isin(["PC", "CP", "KP", "APC"])]
    eper = e["Period (days)"].to_numpy(float)
    etoi = e["TOI"].to_numpy()
    order = [(i, float(p.params[0]),
              str(etoi[int(np.nanargmin(np.abs(eper - float(p.params[0]))))]) if len(eper) else "")
             for i, p in enumerate(pp)]
    return order, float(sd[0])


def failed_list(sp1, snrd, snr, ntr, sharp):
    f = []
    if sp1 >= SPURIOUS_MAX:
        f.append("spurious")
    if snrd <= SNRD_MIN and snr <= SNR_OVERRIDE:
        f.append("snrd")
    if ntr < NTRANSITS_MIN:
        f.append("ntransits")
    if sharp and snrd > HARMONICS_SNRD and snr <= SNR_OVERRIDE:
        f.append("sharp_freq")
    return f


def main(apply):
    w = pd.read_csv(TOIS)
    cols = list(w.columns)
    add_by_tic = {}     # tic -> list of (new_row_dict, plot_index)
    plot_jobs = []      # (src_png, tic, known_idx)

    for f in sorted(glob.glob(os.path.join(KNOWN, "*.csv"))):
        tic = int(os.path.basename(f).split(".")[0])
        rr = pd.read_csv(f, sep="\t")
        wsub = w[w["TIC"] == tic]
        if len(wsub) == 0 or len(rr) == 0:
            continue
        wtoi = set(str(t) for t in wsub["TOI"])
        missing = rr[~rr["TOI"].astype(str).isin(wtoi)]
        if len(missing) == 0:
            continue
        sib = wsub.iloc[0]                       # per-star template (NST + stellar cols)
        order, ts = planet_order(tic)
        sharp = has_sharp_peak(tic)
        n_existing = len(wsub)
        new = []
        # order new planets by Stage-A plot index for determinism
        def plot_idx(cp):
            k = int(np.argmin([abs(per - cp) for (_, per, _) in order])) if order else -1
            return order[k][0] if k >= 0 else -1
        missing = missing.assign(_pidx=[plot_idx(float(r["catalog_period"])) for _, r in missing.iterrows()])
        missing = missing.sort_values("_pidx")
        for j, (_, r) in enumerate(missing.iterrows()):
            snr = float(r["SNR"]); snrd = float(r["snrd_pvalue"]); sp1 = float(r["spurious1"])
            ntr = int(r["num_available_transits"])
            fails = failed_list(sp1, snrd, snr, ntr, sharp)
            # Singh-Maddala posterior-mean-SF p-value: interpolate the sibling star's log10-SF grid.
            log10p = sm_pvalue.log10p_from_grid(sib["sm_sf_grid"], snr)
            if not np.isfinite(log10p):
                log10p = sib["log10(p value)"]
            row = {c: sib[c] for c in cols}     # start from sibling (per-star cols), then override
            row.update({
                "TIC": tic, "TOI": r["TOI"],
                "Period": float(r["period"]), "Phase": float(r["phase"]), "Tau": float(r["tau"]),
                "SNR": snr, "Radius_planet": float(r["radius"]),
                "Radius_planet_errp": float(r["radiusp"]), "Radius_planet_errm": float(r["radiusm"]),
                "Number of Valid Transits": ntr,
                "Epoch": ts + float(r["phase"]), "Duration": 2.0 * float(r["tau"]),
                "err_Period": float(r["err_period"]), "err_Epoch": float(r["err_phase"]),
                "err_Duration": 2.0 * float(r["err_tau"]),
                "log10(p value)": log10p,
                "passed_all_tests": bool(len(fails) == 0), "failed_tests": "|".join(fails),
            })
            known_idx = n_existing + j
            new.append(row)
            src = os.path.join(PLOTS_SRC, str(tic), f"{int(r['_pidx'])}.png")
            plot_jobs.append((src, tic, known_idx, os.path.exists(src)))
            print(f"  + TIC {tic} TOI {r['TOI']} -> new known_idx {known_idx} "
                  f"(plot {int(r['_pidx'])}{'' if os.path.exists(src) else ' MISSING'})  "
                  f"P={row['Period']:.4f} SNR={snr:.1f} {'PASS' if not fails else 'FAIL '+'|'.join(fails)}")
        add_by_tic[tic] = new

    n_new = sum(len(v) for v in add_by_tic.values())
    print(f"\nnew rows to add: {n_new}  across {len(add_by_tic)} TICs | plots to copy: {len(plot_jobs)}")
    missing_plots = [p for p in plot_jobs if not p[3]]
    if missing_plots:
        print("  WARNING missing plot sources:", [(p[1], p[0]) for p in missing_plots])

    if not apply:
        print("\nDRY RUN — pass --apply to write tois.csv and plots")
        return

    # splice new rows in right after each TIC's existing rows
    w = w.reset_index(drop=True)
    last_pos = {tic: w.index[w["TIC"] == tic].tolist()[-1] for tic in add_by_tic}
    parts = []
    for pos in range(len(w)):
        parts.append(w.iloc[[pos]])
        for tic, after in last_pos.items():
            if pos == after:
                parts.append(pd.DataFrame(add_by_tic[tic]).reindex(columns=cols))
    out = pd.concat(parts, ignore_index=True)
    out.to_csv(TOIS, index=False)
    print(f"\nWROTE {TOIS}  ({len(w)} -> {len(out)} rows)")

    from PIL import Image
    for src, tic, known_idx, ok in plot_jobs:
        if not ok:
            continue
        dst = os.path.join(PLOTS_DST, str(tic), f"{known_idx}.jpg")
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        im = Image.open(src).convert("RGB")
        if im.width > 900:
            im = im.resize((900, int(im.height * 900 / im.width)), Image.LANCZOS)
        im.save(dst, "JPEG", quality=88)
    print(f"copied {sum(1 for p in plot_jobs if p[3])} plots")


if __name__ == "__main__":
    main(apply="--apply" in sys.argv)
