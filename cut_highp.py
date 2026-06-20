"""cut_highp.py — apply the p <= 0.01 significance cut to the new-candidate catalog.

Run AFTER apply_singh_maddala.py --apply, which rebuilds tois_new.csv to the full pre-cut candidate
set with the new Singh-Maddala p-values. This keeps rows with log10(p value) <= -2 (p <= 0.01) in
tois_new.csv; the rest move to tois_new_dropped_highp.csv. Idempotent (re-running on the already-cut
file keeps all rows, since they already satisfy the cut). Dry-run by default; pass --apply to write.
"""
import sys
import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/TESS corrected"
LOGP_CUT = -2.0          # p <= 0.01


def main(apply):
    path = f"{HERE}/tois_new.csv"
    df = pd.read_csv(path)                               # the live (post-apply) catalog
    lp = pd.to_numeric(df["log10(p value)"], errors="coerce")
    keep = lp <= LOGP_CUT                                # significant; NaN -> dropped
    drop = ~keep
    print(f"input={len(df)} | keep (p<=0.01)={int(keep.sum())} | drop (p>0.01 or NaN)={int(drop.sum())}")
    if apply:
        df[keep].to_csv(path, index=False)
        df[drop].to_csv(f"{HERE}/tois_new_dropped_highp.csv", index=False)
        print(f"   WROTE {path} ({int(keep.sum())} rows) + tois_new_dropped_highp.csv ({int(drop.sum())} rows)")
    else:
        print("DRY RUN — pass --apply to write")


if __name__ == "__main__":
    main("--apply" in sys.argv)
