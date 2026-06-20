"""apply_singh_maddala.py — write the Bayesian Singh-Maddala p-values into the catalogs.

Reads the Phase-2 fit outputs (pvalue/sm_grid.csv, pvalue/sm_log10p_{tois,toisnew}.csv) and, for both
tois.csv and tois_new.csv:
  * sets 'log10(p value)' to the posterior-mean-SF value (per row; = grid interpolation at the SNR),
  * adds the per-star column 'sm_sf_grid' = log10(posterior-mean SF) tabulated on the shared SNR grid
    (sm_pvalue.XGRID), used by the website SF plot and by the merge-site recompute (both interpolate).
Legacy 'μ(SNR | null)', 'σ(SNR | null)', 'nst_samples' are kept (reference). Starts each catalog from
the pre-Singh-Maddala backup so re-runs are idempotent (no stale sm_* columns).
Dry-run by default; pass --apply to write.
"""
import os, sys, shutil
import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/TESS corrected"
FIT = "/pscratch/sd/j/julius/exoprob/pvalue"

grid = pd.read_csv(f"{FIT}/sm_grid.csv").set_index("TIC")["sf_grid"]
lp = {
    "tois": pd.read_csv(f"{FIT}/sm_log10p_tois.csv").set_index("row_index")["log10p"],
    "tois_new": pd.read_csv(f"{FIT}/sm_log10p_toisnew.csv").set_index("row_index")["log10p"],
}


def update(name, apply):
    path = f"{HERE}/{name}.csv"
    bak = f"{path}.bak_before_singh_maddala"
    src = bak if os.path.exists(bak) else path          # always rebuild from the pre-SM state
    df = pd.read_csv(src)
    new_lp = df.index.map(lambda i: lp[name].get(i, np.nan))
    sfg = df["TIC"].map(grid)
    print(f"{name}: rows={len(df)} | new log10p finite={int(np.isfinite(new_lp).sum())} | sm_sf_grid finite={int(sfg.notna().sum())}")
    if apply:
        df["log10(p value)"] = new_lp
        if "sm_sf_grid" in df.columns:
            df = df.drop(columns="sm_sf_grid")
        cols = list(df.columns)
        anchor = "σ(SNR | null)"
        pos = cols.index(anchor) + 1 if anchor in cols else len(cols)
        df.insert(pos, "sm_sf_grid", sfg.to_numpy())
        df.to_csv(path, index=False)
        print(f"   WROTE {path}")


if __name__ == "__main__":
    apply = "--apply" in sys.argv
    if apply:
        for n in ("tois", "tois_new"):
            if not os.path.exists(f"{HERE}/{n}.csv.bak_before_singh_maddala"):
                shutil.copy(f"{HERE}/{n}.csv", f"{HERE}/{n}.csv.bak_before_singh_maddala")
    for n in ("tois", "tois_new"):
        update(n, apply)
    if not apply:
        print("\nDRY RUN — pass --apply to write")
