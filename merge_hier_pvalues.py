"""merge_hier_pvalues.py -- put the rerun's hierarchical NST p-values into the website catalog.

For every row of tois.csv and tois_new.csv, the rerun's null unit of that TIC whose searched period
window contains the row's Period supplies the null columns (nst_samples, mu/sigma(SNR | null),
sm_sf_grid, log10(p value), log10(p local), n_period_bins). The local p is the exact quadrature value at
the row's SNR from row_pvalues.csv (sm_pvalue.log10p_from_grid interpolation only as the fallback); the
reported log10(p value) is GLOBAL, p = min(1, n_bins x p_local) with n_bins the star's number of
period-local null bins (results/star_nbins.csv; user decision 2026-09-24). sm_sf_grid stays the LOCAL
survival function, so the site's panel can draw it against the unit's own null samples. Rows with no
matching unit have their null columns blanked. Dry run by default; --apply writes.
Usage: python merge_hier_pvalues.py [RUN ...] [--apply]   (default: the pooled run directory)
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/exoplanets/TESS corrected"
sys.path.insert(0, HERE)
import sm_pvalue

DEFAULT_RUNS = ["Rerun_20260910"]   # the pooled fit over Rerun_20260910 + NewRun_20260919 writes all units here

NULL_COLS = ["log10(p value)", "log10(p local)", "n_period_bins", "μ(SNR | null)", "σ(SNR | null)", "sm_sf_grid", "nst_samples"]
NBINS = "/global/u2/j/julius/exoplanets/results/star_nbins.csv"


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
        units = pd.read_csv(outdir + "units.csv")                      # every null window, empty ones included
        samples = pd.read_csv(outdir + "nst_fullrun.csv")[["unit", "samples"]]
        grid = pd.read_csv(outdir + "hier_grid.csv")[["unit", "sf_grid"]]
        merged = units.merge(samples, on="unit", how="left").merge(grid, on="unit", how="left")
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


def load_exact(runs):
    """Exact log10 p per (TIC, Period, SNR): catalog rows and batch-0 candidates, from tess_rerun_pvalues."""
    exact = {}
    for run in runs:
        path = f"/pscratch/sd/j/julius/exoprob/results/hierarchical_tess/{run.lower()}/row_pvalues.csv"
        if os.path.exists(path):
            for r in pd.read_csv(path).itertuples(index=False):
                if np.isfinite(r.log10p):
                    exact.setdefault((int(r.TIC), round(float(r.Period), 6), round(float(r.SNR), 4)), r.log10p)
    return exact


def to_global(log10p_local, n_bins):
    """Bonferroni over the star's period bins: log10 min(1, n_bins x p_local)."""
    return min(0., log10p_local + np.log10(n_bins))


def merge(name, units, apply, exact=None, nbins=None):
    """Rewrite the null columns of one catalog file from the rerun units; returns the summary row.

    log10(p local) is the exact per-row quadrature value where one exists, else the grid interpolation;
    log10(p value) is its global (Bonferroni) version."""
    path = os.path.join(HERE, name)
    df = pd.read_csv(path)
    by_tic = {t: g for t, g in units.groupby("TIC")}
    old_lp = df["log10(p value)"].to_numpy(dtype=float).copy()    # a view would follow the rewrite
    exact = exact or {}
    for offset, column in enumerate(("log10(p local)", "n_period_bins"), start=1):
        if column not in df.columns:
            df.insert(df.columns.get_loc("log10(p value)") + offset, column, np.nan)

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
        n = int(nbins[int(row["TIC"])])
        df.loc[i, "sm_sf_grid"] = unit["sf_grid"]
        df.loc[i, "n_period_bins"] = n
        key = (int(row["TIC"]), round(float(row["Period"]), 6), round(float(row["SNR"]), 4))
        if key in exact:
            local = round(float(exact[key]), 4)
        else:
            local = sm_pvalue.log10p_from_grid(unit["sf_grid"], float(row["SNR"]))
            n_grid += 1
        df.loc[i, "log10(p local)"] = local
        df.loc[i, "log10(p value)"] = round(to_global(local, n), 4)
        n_matched += 1

    new_lp = df["log10(p value)"].to_numpy(dtype=float)
    both = np.isfinite(old_lp) & np.isfinite(new_lp)
    shift = new_lp[both] - old_lp[both]
    summary = dict(file=name, rows=len(df), matched=n_matched, unmatched=len(df) - n_matched, grid_fallback=n_grid,
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
    exact = load_exact(runs)
    nbins = pd.read_csv(NBINS).set_index("kepid")["n_bins"]
    print(f"  {len(exact)} exact per-row p-values; period-bin counts for {len(nbins)} stars")
    rows = [merge(name, units, apply, exact, nbins) for name in ("tois.csv", "tois_new.csv")]
    table = pd.DataFrame(rows)
    print(table.round(3).to_string(index=False))

    # a jump in `unmatched` signals a missing run on the command line, not catalog growth
    if table["unmatched"].sum():
        print(f"\n{int(table['unmatched'].sum())} rows matched no unit and were BLANKED. "
              f"Check that every run serving these catalogs is listed.")
    print("written" if apply else "dry run (pass --apply to write)")
