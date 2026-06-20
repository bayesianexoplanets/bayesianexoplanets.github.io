"""assign_proposed_tois.py — populate the proposed_toi column in tois_new.csv.

Each new candidate is an ADDITIONAL planet at a star that already carries an ExoFOP TOI, so it
reuses that star's TOI integer and takes the next planet index:

    proposed_toi = <ExoFOP TOI integer>.0X ,  X = (# ExoFOP planets at that TIC) + 1   (always >= 2)

If a TIC ever has more than one new candidate they get consecutive indices (N+1, N+2, ...).
"""
import os
import pandas as pd

HERE = "/global/u2/j/julius/TESS corrected"
EXO = "/global/u2/j/julius/TESS/tois.csv"

# ExoFOP: TIC -> [TOI integer base, number of planets at that TIC]
E = pd.read_csv(EXO)
exo = {}
for _, r in E.iterrows():
    try:
        tic = int(r["TIC ID"]); base = int(float(r["TOI"]))
    except (ValueError, TypeError):
        continue
    if tic not in exo:
        exo[tic] = [base, 0]
    exo[tic][1] += 1

path_new = os.path.join(HERE, "tois_new.csv")
df = pd.read_csv(path_new)
print(f"tois_new.csv: {len(df)} rows")

proposed, missing = [], []
offset = {}                     # tic -> # new candidates already assigned for this TIC
for _, row in df.iterrows():
    tic = int(row["TIC"])
    if tic not in exo:
        missing.append(tic)
        proposed.append("")
        continue
    base, n_exo = exo[tic]
    k = offset.get(tic, 0)
    idx = n_exo + 1 + k         # X = N+1 (then N+2, ... for any further candidates on the same TIC)
    offset[tic] = k + 1
    proposed.append(f"{base}.{idx:02d}")

df["proposed_toi"] = proposed
df.to_csv(path_new, index=False)

if missing:
    print(f"  WARNING: {len(missing)} TIC(s) have no ExoFOP TOI (left blank): {missing[:10]}")
print(f"  assigned {sum(1 for p in proposed if p)} proposed TOIs (X = ExoFOP planet count + 1)")
print(f"  examples: {proposed[:6]}")
print(f"  wrote {path_new}")
