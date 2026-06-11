"""assign_proposed_tois.py — populate proposed_toi column in tois_new.csv.

Numbers continue from the max integer TOI in tois.csv (currently 7475).
Rows are processed in existing CSV order (SNR descending). Each unique TIC
gets a new integer; multiple candidates from the same TIC share the integer
and get .01, .02, etc.
"""
import os
import pandas as pd

HERE = "/global/u2/j/julius/TESS corrected"

# Find max integer TOI in tois.csv
path_known = os.path.join(HERE, 'tois.csv')
df_known = pd.read_csv(path_known)
max_int = int(df_known['TOI'].dropna().astype(float).apply(lambda v: int(v)).max())
print(f"Max TOI integer in catalog: {max_int}")

# Assign proposed TOIs to new candidates
path_new = os.path.join(HERE, 'tois_new.csv')
df_new = pd.read_csv(path_new)
print(f"tois_new.csv: {len(df_new)} rows")

proposed = []
toi_counter = max_int
planet_per_tic = {}  # tic -> current planet number within this TOI group

for _, row in df_new.iterrows():
    tic = int(row['TIC'])
    if tic not in planet_per_tic:
        toi_counter += 1
        planet_per_tic[tic] = 1
    else:
        planet_per_tic[tic] += 1
    proposed.append(f"{toi_counter}.{planet_per_tic[tic]:02d}")

df_new['proposed_toi'] = proposed

# Sanity check
n_unique = df_new['proposed_toi'].nunique()
print(f"  Assigned {n_unique} unique proposed TOIs ({len(df_new)} candidates)")
print(f"  Range: {df_new['proposed_toi'].iloc[0]} .. {df_new['proposed_toi'].iloc[-1]}")
print(f"  Max TOI used: {toi_counter}")

df_new.to_csv(path_new, index=False)
print(f"  Wrote {path_new}")
print("\nDone.")
