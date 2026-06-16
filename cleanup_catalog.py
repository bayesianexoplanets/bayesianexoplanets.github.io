"""cleanup_catalog.py — drop specific (TIC, TOI) rows from the website catalog and
rebuild the affected TICs' per-planet plots so known_idx stays consistent.

The frontend assigns known_idx by per-TIC row order (index.html L573) and shows
plots/{TIC}/{known_idx}.jpg. After dropping a row, the surviving rows reindex, so
every surviving planet's plot must be re-emitted at its new known_idx (from the
KnownRun Stage-A png, matched via read_known order) and stale jpgs removed.

Edit DROP below. Dry-run by default; pass --apply.
"""
import os
import sys
import glob

os.chdir("/global/u2/j/julius")
sys.path.insert(0, "/global/u2/j/julius")
import numpy as np
import pandas as pd
from TESS.load_tess import StarInfo_tess, read_known_planets_tess

HERE = "/global/u2/j/julius/TESS corrected"
TOIS = os.path.join(HERE, "tois.csv")
PLOTS_SRC = "/pscratch/sd/j/julius/exoprob/KnownRun/plots"
PLOTS_DST = os.path.join(HERE, "plots")
BULK = "/pscratch/sd/j/julius/Bulk Download/Data"

# rows to drop:  TIC -> [TOI, ...]
DROP = {27491137: ["2076.01"]}


def planet_order(tic):
    sd = np.load(os.path.join(BULK, f"StellarData_{tic}.npy")).copy()
    if np.isnan(sd[3]):
        sd[3] = 1
    if np.isnan(sd[-1]):
        sd[-1] = 1
        sd[-2] = 1
    star = StarInfo_tess(*sd[1:])
    pp = read_known_planets_tess(star, sd[0], False)
    e = pd.read_csv("TESS/tois.csv")
    e = e[e["TIC ID"] == tic]
    e = e[e["TFOPWG Disposition"].isin(["PC", "CP", "KP", "APC"])]
    eper = e["Period (days)"].to_numpy(float)
    etoi = e["TOI"].to_numpy()
    # TOI -> Stage-A plot index
    out = {}
    for i, p in enumerate(pp):
        per = float(p.params[0])
        toi = str(etoi[int(np.nanargmin(np.abs(eper - per)))]) if len(eper) else ""
        out[toi] = i
    return out


def main(apply):
    from PIL import Image
    w = pd.read_csv(TOIS)
    drop_idx = []
    for tic, tois in DROP.items():
        sub = w[w["TIC"] == tic]
        toi2plot = planet_order(tic)
        survivors = [str(t) for t in sub["TOI"] if str(t) not in set(tois)]
        print(f"TIC {tic}: drop {tois} | survivors (new known_idx order): {survivors}")
        # rows to drop
        for ix, r in sub.iterrows():
            if str(r["TOI"]) in set(tois):
                drop_idx.append(ix)
        # plan plot rebuild: survivor j -> jpg j from KnownRun plot toi2plot[toi]
        for j, toi in enumerate(survivors):
            pidx = toi2plot.get(toi)
            src = os.path.join(PLOTS_SRC, str(tic), f"{pidx}.png") if pidx is not None else None
            dst = os.path.join(PLOTS_DST, str(tic), f"{j}.jpg")
            ok = src and os.path.exists(src)
            print(f"   known_idx {j} <- TOI {toi} (KnownRun plot {pidx}{'' if ok else ' MISSING'}) -> {dst}")
            if apply and ok:
                im = Image.open(src).convert("RGB")
                if im.width > 900:
                    im = im.resize((900, int(im.height * 900 / im.width)), Image.LANCZOS)
                im.save(dst, "JPEG", quality=88)
        # remove stale jpgs beyond survivor count
        for stale in range(len(survivors), len(sub)):
            p = os.path.join(PLOTS_DST, str(tic), f"{stale}.jpg")
            print(f"   remove stale jpg {p} (exists={os.path.exists(p)})")
            if apply and os.path.exists(p):
                os.remove(p)

    if apply:
        w2 = w.drop(index=drop_idx).reset_index(drop=True)
        w2.to_csv(TOIS, index=False)
        print(f"\nWROTE {TOIS}  ({len(w)} -> {len(w2)} rows; dropped {len(drop_idx)})")
    else:
        print(f"\nDRY RUN — would drop {len(drop_idx)} rows. Pass --apply.")


if __name__ == "__main__":
    main(apply="--apply" in sys.argv)
