"""add_failed_tests.py — add failed_tests column to tois.csv and tois_new.csv.

For tois_new.csv: looks up the matching row in batch0/{tic}.0.csv by event_id
(= _cand_idx) and applies the three measurable vetting thresholds.

For tois.csv: known TOIs don't have a direct batch0 link; set to '' if
passed_all_tests is True, 'see_pipeline' otherwise.
"""
import os
import pandas as pd

HERE = "/global/u2/j/julius/TESS corrected"
FULLRUN = "/pscratch/sd/j/julius/exoprob/FullRun/candidates"

TESTS = [
    ('spurious',  'spurious1',              lambda v: v >= 0.5),
    ('snrd',      'snrd_pvalue',            lambda v: v <= 0.01),
    ('ntransits', 'num_available_transits', lambda v: v < 3),
]


def batch0_row(tic, cand_idx):
    path = os.path.join(FULLRUN, "batch0", f"{int(tic)}.0.csv")
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, sep='\t')
        match = df[df['event_id'] == int(cand_idx)]
        return match.iloc[0] if len(match) > 0 else None
    except Exception:
        return None


def failed_tests_new(tic, cand_idx):
    row = batch0_row(tic, cand_idx)
    if row is None:
        return ''
    failed = []
    for name, col, test in TESTS:
        if col in row.index and pd.notna(row[col]):
            if test(row[col]):
                failed.append(name)
    return '|'.join(failed)


# ── tois_new.csv ──────────────────────────────────────────────────────────────
path_new = os.path.join(HERE, 'tois_new.csv')
df_new = pd.read_csv(path_new)
print(f"tois_new.csv: {len(df_new)} rows")

df_new['failed_tests'] = df_new.apply(
    lambda r: failed_tests_new(r['TIC'], r['_cand_idx']), axis=1
)
n_fail = (df_new['failed_tests'] != '').sum()
print(f"  {n_fail} rows have at least one failed test")
print("  Breakdown:")
for name, _, _ in TESTS:
    n = df_new['failed_tests'].str.contains(name, na=False).sum()
    print(f"    {name}: {n}")
df_new.to_csv(path_new, index=False)
print(f"  Wrote {path_new}")


# ── tois.csv ──────────────────────────────────────────────────────────────────
path_known = os.path.join(HERE, 'tois.csv')
df_known = pd.read_csv(path_known)
print(f"\ntois.csv: {len(df_known)} rows")
df_known['failed_tests'] = df_known['passed_all_tests'].apply(
    lambda v: '' if v is True or v == True else 'see_pipeline'
)
n_fail = (df_known['failed_tests'] != '').sum()
print(f"  {n_fail} rows with failed_tests='see_pipeline'")
df_known.to_csv(path_known, index=False)
print(f"  Wrote {path_known}")

print("\nDone.")
