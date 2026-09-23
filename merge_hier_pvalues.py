"""merge_hier_pvalues.py -- put the rerun's hierarchical NST p-values into the website catalog.

For every row of tois.csv and tois_new.csv, the rerun's null unit of that TIC whose searched period
window contains the row's Period supplies the null columns (nst_samples, mu/sigma(SNR | null),
sm_sf_grid, log10(p value)), via the same sm_pvalue.log10p_from_grid path apply_singh_maddala.py
uses. Rows with no matching unit have their null columns blanked. Dry run by default; --apply writes.
Usage: python merge_hier_pvalues.py [RUN ...] [--apply]   (default: both runs serving the catalogs)
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/exoplanets/TESS corrected"
sys.path.insert(0, HERE)
import sm_pvalue

DEFAULT_RUNS = ["Rerun_20260910", "NewRun_20260919"]   # every run that serves these catalogs

NULL_COLS = ["log10(p value)", "μ(SNR | null)", "σ(SNR | null)", "sm_sf_grid", "nst_samples"]


def load_units(runs):
    """Unit table (unit, TIC, period_min, period_max, samples, sf_grid) pooled over several reruns.

    Several runs, because the catalog is no longer served by one. The 408 new hosts were searched in
    `NewRun_20260919` and everything else in `Rerun_20260910`, and a single-run call here does not
    simply miss the other run's rows: `merge` BLANKS every row it cannot match, so running this with
    one run would wipe the null columns of every row belonging to the other. Unit ids are only unique
    within a run, so they are namespaced as "run:unit" before pooling.
    """
    tables = []
    for run in runs:
        outdir = f"/pscratch/sd/j/julius/exoprob/results/hierarchical_tess/{run.lower()}/"
        units = pd.read_csv(outdir + "nst_fullrun.csv")                # unit, TIC, period_min, period_max, samples
        grid = pd.read_csv(outdir + "hier_grid.csv")[["unit", "sf_grid"]]
        merged = units.merge(grid, on="unit", how="left")
        merged["unit"] = [f"{run}:{u}" for u in merged["unit"]]
        merged["run"] = run
        tables.append(merged)
        print(f"  {run}: {len(merged)} units on {merged['TIC'].nunique()} stars")

    pooled = pd.concat(tables, ignore_index=True).set_index("unit")
    clashes = pooled.reset_index().duplicated(["TIC", "period_min", "period_max"]).sum()
    if clashes:
        print(f"  NOTE: {clashes} (TIC, window) pairs appear in more than one run; the first run listed wins")
    return pooled


def match_unit(units_of_tic, period):
    """Unit row whose window contains the period, or None.

    With several runs pooled a star can own two units covering the same period; the first row wins,
    which is the first run named on the command line.
    """
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
    runs = [a for a in sys.argv[1:] if not a.startswith("--")] or DEFAULT_RUNS
    print("pooling null units over", ", ".join(runs))
    units = load_units(runs)
    rows = [merge(name, units, apply) for name in ("tois.csv", "tois_new.csv")]
    table = pd.DataFrame(rows)
    print(table.round(3).to_string(index=False))

    # a jump in `unmatched` signals a missing run on the command line, not catalog growth
    if table["unmatched"].sum():
        print(f"\n{int(table['unmatched'].sum())} rows matched no unit and were BLANKED. "
              f"Check that every run serving these catalogs is listed.")
    print("written" if apply else "dry run (pass --apply to write)")
