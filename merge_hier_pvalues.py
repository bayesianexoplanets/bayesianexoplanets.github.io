"""merge_hier_pvalues.py -- put the rerun's hierarchical NST p-values into the website catalog.

For every row of tois.csv and tois_new.csv, the rerun's null unit of that TIC whose searched
period window contains the row's Period supplies the null columns: nst_samples (the unit's
per-batch max SNR), μ/σ(SNR | null) (their mean and population std), sm_sf_grid (the unit's
hierarchical posterior-mean log10 SF on sm_pvalue.XGRID) and log10(p value) evaluated at the
row's SNR through sm_pvalue.log10p_from_grid, exactly as apply_singh_maddala.py does. Rows
whose TIC has no unit containing their period (fewer than three transits fit the data, or a
star outside the rerun) have their null columns blanked and appear without a null run, the
convention merge_rerun23.py set for a deleted null. Dry run by default; --apply writes.
Usage: python merge_hier_pvalues.py [--apply] [--run Rerun_20260910]
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/exoplanets/TESS corrected"
sys.path.insert(0, HERE)
import sm_pvalue

NULL_COLS = ["log10(p value)", "μ(SNR | null)", "σ(SNR | null)", "sm_sf_grid", "nst_samples"]


def load_units(run):
    """Unit table (unit, TIC, period_min, period_max, samples, sf_grid) of a rerun."""
    outdir = f"/pscratch/sd/j/julius/exoprob/results/hierarchical_tess/{run.lower()}/"
    units = pd.read_csv(outdir + "nst_fullrun.csv")                    # unit, TIC, period_min, period_max, samples
    grid = pd.read_csv(outdir + "hier_grid.csv")[["unit", "sf_grid"]]
    return units.merge(grid, on="unit", how="left").set_index("unit")


def match_unit(units_of_tic, period):
    """Unit row whose window contains the period, or None."""
    if units_of_tic is None:
        return None
    hit = units_of_tic[(units_of_tic["period_min"] <= period) & (period < units_of_tic["period_max"])]
    return hit.iloc[0] if len(hit) else None


def merge(name, units, apply):
    """Rewrite the null columns of one catalog file from the rerun units; returns the summary row."""
    path = os.path.join(HERE, name)
    df = pd.read_csv(path)
    by_tic = {t: g for t, g in units.groupby("TIC")}
    old_lp = df["log10(p value)"].to_numpy(dtype=float).copy()    # a view would follow the rewrite

    n_matched, n_grid = 0, 0
    for i, row in df.iterrows():
        unit = match_unit(by_tic.get(int(row["TIC"])), float(row["Period"]))
        if unit is None or not isinstance(unit["sf_grid"], str):
            df.loc[i, NULL_COLS] = np.nan
            continue
        raw = unit["samples"] if isinstance(unit["samples"], str) else ""      # a window whose null searches found nothing
        samples = np.array([float(x) for x in raw.split("|") if x])
        df.loc[i, "nst_samples"] = "|".join(f"{x:.6f}" for x in samples) if len(samples) else np.nan
        df.loc[i, "μ(SNR | null)"] = round(float(samples.mean()), 4) if len(samples) else np.nan
        df.loc[i, "σ(SNR | null)"] = round(float(samples.std()), 6) if len(samples) else np.nan
        df.loc[i, "sm_sf_grid"] = unit["sf_grid"]
        df.loc[i, "log10(p value)"] = sm_pvalue.log10p_from_grid(unit["sf_grid"], float(row["SNR"]))
        n_matched += 1

    new_lp = df["log10(p value)"].to_numpy(dtype=float)
    both = np.isfinite(old_lp) & np.isfinite(new_lp)
    shift = new_lp[both] - old_lp[both]
    summary = dict(file=name, rows=len(df), matched=n_matched, unmatched=len(df) - n_matched,
                   old_significant=int((old_lp < -2).sum()), new_significant=int((new_lp < -2).sum()),
                   shift_median=float(np.median(shift)) if len(shift) else np.nan,
                   shift_p5=float(np.percentile(shift, 5)) if len(shift) else np.nan,
                   shift_p95=float(np.percentile(shift, 95)) if len(shift) else np.nan)
    if apply:
        df.to_csv(path, index=False)
    return summary


if __name__ == "__main__":
    apply = "--apply" in sys.argv
    run = sys.argv[sys.argv.index("--run") + 1] if "--run" in sys.argv else "Rerun_20260910"
    units = load_units(run)
    rows = [merge(name, units, apply) for name in ("tois.csv", "tois_new.csv")]
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    print("written" if apply else "dry run (pass --apply to write)")
