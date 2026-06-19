"""Phase 1: move TOIs that cannot be validated (< 3 valid transits) out of the website
catalog into a separate file, so their transits can still be masked when fitting the FGP
but they are no longer presented as validated candidates.

Removed = failed_tests contains 'ntransits'  OR  Number of Valid Transits < 3
          OR a manual mis-detection list (real planet has < 3 transits).
Rows are MOVED (full columns) to tois_removed_low_transit.csv, then dropped from tois.csv.
Dry-run by default; pass --apply to write (backs up tois.csv first)."""
import os
import sys
import shutil
import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/TESS corrected"
TOIS = os.path.join(HERE, "tois.csv")
OUT = os.path.join(HERE, "tois_removed_low_transit.csv")

# manual: real planet has <= 3 (signal-bearing) transits even though catalog
# num_available_transits over-counts; user-identified by visual inspection of the plots.
MANUAL = {
    (443666343, "3500.01"), (8988289, "316.01"), (399871011, "943.01"),
    (282794970, "6889.01"), (150907983, "6222.01"), (423670610, "850.01"),
    (437984134, "887.01"), (349488688, "2319.02"), (85293053, "1772.01"),
    (5592720, "4986.01"), (307232923, "6742.01"), (100757807, "811.01"),
    (46020827, "7471.01"), (356747847, "1958.01"),
    # long-ExoFOP: ExoFOP period gives <3 transits (unvalidatable); detected period is a
    # different/spurious signal -> remove (user-confirmed).
    (67395329, "6873.01"), (92443533, "6812.01"), (156514476, "6884.01"),
    (192415680, "6041.01"), (406684949, "5442.01"),
    # circumbinary TOI-1338 with huge TTVs: only a single planet transit is usable (the strong
    # short-period signal is the eclipsing binary, not the planet) -> unvalidatable (user-confirmed).
    (260128333, "1338.01"),
}

w = pd.read_csv(TOIS)
nvt = pd.to_numeric(w["Number of Valid Transits"], errors="coerce")
ft = w["failed_tests"].astype(str)

flag_ntr = ft.str.contains("ntransits", na=False)
flag_cnt = nvt < 3
flag_manual = w.apply(lambda r: (int(r["TIC"]), str(r["TOI"])) in MANUAL, axis=1)
remove = flag_ntr | flag_cnt | flag_manual

print(f"catalog rows: {len(w)}")
print(f"  flagged 'ntransits': {int(flag_ntr.sum())}")
print(f"  Number of Valid Transits < 3: {int(flag_cnt.sum())}")
print(f"  manual: {int(flag_manual.sum())}  {sorted(MANUAL)}")
print(f"  TOTAL to remove (union): {int(remove.sum())}")
print("\nremoved rows:")
for _, r in w[remove].iterrows():
    print(f"  TIC {int(r['TIC'])} TOI {r['TOI']}: P={r['Period']:.4f} nVT={r['Number of Valid Transits']} "
          f"fail={r['failed_tests']} pass={r['passed_all_tests']}")

removed = w[remove].copy()
kept = w[~remove].copy()

if "--apply" in sys.argv:
    shutil.copy(TOIS, TOIS + ".bak_before_lowtransit")
    # append to an existing removed file if present (idempotent-ish: dedupe on TIC+TOI)
    if os.path.exists(OUT):
        prev = pd.read_csv(OUT)
        removed = pd.concat([prev, removed], ignore_index=True).drop_duplicates(["TIC", "TOI"], keep="last")
    removed.to_csv(OUT, index=False)
    kept.to_csv(TOIS, index=False)
    print(f"\nWROTE {OUT} ({len(removed)} rows) and {TOIS} ({len(kept)} rows); backup .bak_before_lowtransit")
else:
    print(f"\nDRY RUN — would keep {len(kept)}, move {len(removed)} to {OUT}. pass --apply")
