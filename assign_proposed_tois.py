"""assign_proposed_tois.py — populate the proposed_toi column in tois_new.csv.

A new candidate on a star that already carries an ExoFOP TOI takes that star's TOI integer and the planet
index after the highest one ExoFOP has for it, retired and false-positive entries included since their numbers
are taken (user rule 2026-09-24):

    proposed_toi = <TOI integer>.<highest ExoFOP index + 1>      e.g. TOI 10.01 known -> 10.02

Several candidates on one star take consecutive indices in the list's order (strongest first). A star with
no ExoFOP TOI has no integer to extend and is left blank.
"""
import os
import pandas as pd

HERE = "/global/u2/j/julius/exoplanets/TESS corrected"
EXO = "/global/u2/j/julius/exoplanets/TESS/tois.csv"

# ExoFOP: TIC -> [TOI integer, highest planet index at that TIC]
E = pd.read_csv(EXO)
exo = {}
for _, r in E.iterrows():
    try:
        tic, toi = int(r["TIC ID"]), float(r["TOI"])
    except (ValueError, TypeError):
        continue
    base, index = int(toi), int(round((toi - int(toi)) * 100))
    if tic not in exo or (base, index) > tuple(exo[tic]):
        exo[tic] = [base, index]

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
    base, highest = exo[tic]
    k = offset.get(tic, 0)
    idx = highest + 1 + k       # the index after ExoFOP's highest, then consecutive for further candidates
    offset[tic] = k + 1
    proposed.append(f"{base}.{idx:02d}")

df["proposed_toi"] = proposed
df.to_csv(path_new, index=False)

if missing:
    print(f"  WARNING: {len(missing)} TIC(s) have no ExoFOP TOI (left blank): {missing[:10]}")
print(f"  assigned {sum(1 for p in proposed if p)} proposed TOIs (highest ExoFOP index + 1)")
print(f"  examples: {proposed[:6]}")
print(f"  wrote {path_new}")
