"""Move re-detections out of the new-candidate list: to the catalog if they are known, off the site if they are duplicates.

A blind-search candidate that is not a new object does not belong on the new-planet page (user
decision 2026-09-18). Three cases, each already carried as a deterministic token by
`build_new_candidates.py`:

- `known_transit` -> PROMOTED to the catalog. These are not duplicates: the star has a single-transit
  TOI that ExoFOP lists with no period at all, and one of the candidate's transits falls on its epoch,
  so the blind search has MEASURED the period of a known single-transit TOI. That is a catalog result
  about a known object, not a discovery.
- `known_harmonic` -> DROPPED. An integer or half-integer ratio of a catalogued period of the same
  star, so it is the known planet's own comb, not a second object.
- `harmonic` -> DROPPED. A ratio of a STRONGER candidate of the same star, i.e. a duplicate of another
  row of this very list.
- `known_eb` -> DROPPED. The host is a known eclipsing binary.

A row that would be promoted but carries `ntransits` (fewer than three transits) is DROPPED instead: the
catalog removes such rows rather than flagging them (user decision 2026-09-21).

The frontend assigns a figure index by per-TIC row order (`index.html`), so a promoted row is appended
at the END of `tois.csv`: every existing row keeps its index and the new row takes the next one.

Dry run by default; pass --apply.
"""
import os
import shutil
import sys

import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/exoplanets/TESS corrected"
HOME = "/global/u2/j/julius/exoplanets/"
sys.path.insert(0, HERE)
sys.path.insert(0, HOME)
from merge_known_run import _to_jpg                                   # noqa: E402

TOIS = os.path.join(HERE, "tois.csv")
TOIS_NEW = os.path.join(HERE, "tois_new.csv")
DROPPED = os.path.join(HERE, "tois_new_dropped.csv")
PLOTS, PLOTS_NEW = os.path.join(HERE, "plots"), os.path.join(HERE, "plots_new")

PROMOTE = "known_transit"
DROP_TOKENS = ("known_harmonic", "harmonic", "known_eb")


def tokens(row):
    return [t for t in str(row.get("failed_tests") or "").split("|") if t]


def single_transit_toi(tic, period, epoch, duration, exofop):
    """The periodless TOI of the star whose epoch one of the candidate's transits covers."""
    rows = exofop[exofop["TIC ID"] == tic]
    period_col = [c for c in exofop.columns if c.startswith("Period")][0]
    epoch_col = [c for c in exofop.columns if c.startswith("Epoch")][0]
    for _, r in rows.iterrows():
        p, e = float(r[period_col]) if pd.notna(r[period_col]) else np.nan, float(r[epoch_col])
        if (np.isfinite(p) and p > 0) or not np.isfinite(e):
            continue
        if abs((e - 2457000. - epoch + period / 2.) % period - period / 2.) < duration:
            return r["TOI"]
    return None


def main(apply):
    website = pd.read_csv(TOIS)
    candidates = pd.read_csv(TOIS_NEW)
    exofop = pd.read_csv(HOME + "TESS/tois.csv")

    too_few = candidates.apply(lambda r: "ntransits" in tokens(r), axis=1)
    in_catalog = candidates.apply(lambda r: bool(np.any((website["TIC"] == r["TIC"])                 # never promote twice
                                                        & (np.abs(website["Period"] / r["Period"] - 1.) < 1e-3))), axis=1)
    promote = candidates[candidates.apply(lambda r: PROMOTE in tokens(r), axis=1) & ~too_few & ~in_catalog]
    drop = candidates[candidates.apply(lambda r: (any(t in tokens(r) for t in DROP_TOKENS) and PROMOTE not in tokens(r))
                                       or PROMOTE in tokens(r), axis=1).to_numpy() & ~candidates.index.isin(promote.index)]
    keep = candidates.drop(index=promote.index.union(drop.index))

    print("candidates %d -> promoted %d, dropped %d, kept %d" % (len(candidates), len(promote), len(drop), len(keep)))

    rows, copies = [], []
    for _, c in promote.iterrows():
        tic = int(c["TIC"])
        toi = single_transit_toi(tic, float(c["Period"]), float(c["Epoch"]), float(c["Duration"]), exofop)
        index_here = int((website["TIC"] == tic).sum())          # the frontend's per-TIC counter
        row = {col: c[col] for col in website.columns if col in c.index}
        row.update(TIC=tic, TOI=toi if toi is not None else c["TOI"],
                   passed_all_tests=False,
                   failed_tests="|".join(t for t in tokens(c) if t != PROMOTE))
        rows.append(row)
        copies.append((tic, int(c["_cand_idx"]), index_here))
        print("  promote TIC %-10d P=%9.5f SNR=%6.2f -> TOI %s, figure index %d (%s)"
              % (tic, c["Period"], c["SNR"], row["TOI"], index_here, row["failed_tests"] or "passes"))

    print("\ndropped:")
    for _, c in drop.iterrows():
        print("  TIC %-10d P=%9.5f SNR=%6.2f  %s" % (int(c["TIC"]), c["Period"], c["SNR"], "|".join(tokens(c))))

    if not apply:
        print("\nDRY RUN, pass --apply")
        return

    pd.concat([website, pd.DataFrame(rows)], ignore_index=True).to_csv(TOIS, index=False)
    keep.to_csv(TOIS_NEW, index=False)
    dropped = drop.assign(drop_reason=["|".join(t for t in tokens(c) if t in DROP_TOKENS + (PROMOTE, "ntransits"))
                                       for _, c in drop.iterrows()])
    previous = pd.read_csv(DROPPED) if os.path.exists(DROPPED) else pd.DataFrame()     # keeps earlier identifications
    pd.concat([previous, dropped], ignore_index=True).drop_duplicates(["TIC", "Period"], keep="last").to_csv(DROPPED, index=False)

    for tic, cand_idx, index_here in copies:
        source = os.path.join(PLOTS_NEW, str(tic), "%d.jpg" % cand_idx)
        target_dir = os.path.join(PLOTS, str(tic))
        os.makedirs(target_dir, exist_ok=True)
        if os.path.exists(source):
            shutil.copyfile(source, os.path.join(target_dir, "%d.jpg" % index_here))
    for tic in drop["TIC"].astype(int):
        stale = os.path.join(PLOTS_NEW, str(tic))
        if os.path.isdir(stale) and not (keep["TIC"].astype(int) == tic).any():
            shutil.rmtree(stale)

    print("\nWROTE %s (%d rows), %s (%d rows), %s (%d rows)"
          % (TOIS, len(website) + len(rows), TOIS_NEW, len(keep), DROPPED, len(drop)))


if __name__ == "__main__":
    main("--apply" in sys.argv)
