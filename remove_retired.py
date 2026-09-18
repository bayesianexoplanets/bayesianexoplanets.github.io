"""Remove signals the community has retired from the catalog and the website.

A row whose TFOPWG disposition is `FP` (false positive) or `FA` (false alarm) is no longer a planet
candidate, so it should not appear in a planet catalog at all. These rows are MOVED, with the
disposition and the TFOP comment that justifies the removal, to `TESS/tois_retired.csv` and
`TESS corrected/tois_retired.csv`, never deleted outright.

The frontend assigns a figure index by per-TIC row order (`index.html`), so removing a row reshuffles
every surviving row of that star. Each affected star therefore has all of its figures re-emitted from
the pretty renders of the known-planet run, matched by TOI, at the new indices.

Dry run by default; pass --apply.
"""
import glob
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
from flag_rules import RETIRED_DISPOSITIONS                           # noqa: E402
from constants import scratch                                         # noqa: E402

WEBSITE = os.path.join(HERE, "tois.csv")
CATALOG = HOME + "TESS/tois_corrected.csv"
EXOFOP = HOME + "TESS/tois.csv"
PLOTS = os.path.join(HERE, "plots")
PRETTY = scratch + "KnownRun_20260916/pretty/"


def retired_keys():
    """{(TIC, TOI rounded): (disposition, comment)} for every signal the community has retired."""
    exofop = pd.read_csv(EXOFOP)
    keys = {}
    for tic, toi, disposition, comment in zip(exofop["TIC ID"], exofop["TOI"],
                                              exofop["TFOPWG Disposition"], exofop["Comments"].fillna("")):
        if str(disposition).strip() in RETIRED_DISPOSITIONS:
            try:
                keys[(int(tic), round(float(toi), 2))] = (str(disposition).strip(), str(comment))
            except (ValueError, TypeError):
                continue
    return keys


def _pretty_for(tic, toi):
    """The rendered figure of one catalog row, or None."""
    path = PRETTY + "%d_%.2f.png" % (int(tic), float(toi))
    if os.path.exists(path):
        return path
    matches = glob.glob(PRETTY + "%d_*.png" % int(tic))
    return matches[0] if len(matches) == 1 else None


def reindex_star(tic, surviving):
    """Re-emit every figure of one star at its new per-TIC index."""
    target = os.path.join(PLOTS, str(int(tic)))
    for stale in glob.glob(os.path.join(target, "*.jpg")):
        os.remove(stale)
    written = 0
    for index, toi in enumerate(surviving):
        source = _pretty_for(tic, toi) if pd.notna(toi) else None
        if source:
            _to_jpg(source, os.path.join(target, "%d.jpg" % index))
            written += 1
    if written == 0 and os.path.isdir(target) and not os.listdir(target):
        os.rmdir(target)
    return written


def main(apply):
    keys = retired_keys()
    website = pd.read_csv(WEBSITE)
    catalog = pd.read_csv(CATALOG, sep="\t")
    for frame in (website, catalog):
        frame["_key"] = [(int(t), round(float(x), 2)) if pd.notna(x) else None
                         for t, x in zip(frame["TIC"], frame["TOI"])]

    web_drop = website["_key"].isin(keys)
    cat_drop = catalog["_key"].isin(keys)
    print("retired signals in the live ExoFOP table: %d" % len(keys))
    print("  website rows to remove : %d of %d (passing %d)"
          % (web_drop.sum(), len(website), int(website.loc[web_drop, "passed_all_tests"].sum())))
    print("  catalog rows to remove : %d of %d" % (cat_drop.sum(), len(catalog)))

    dropped = website[web_drop].copy()
    dropped["disposition"] = [keys[k][0] for k in dropped["_key"]]
    dropped["tfop_comment"] = [keys[k][1] for k in dropped["_key"]]
    print("\nhighest-SNR removals:")
    print(dropped.nlargest(8, "SNR")[["TIC", "TOI", "Period", "SNR", "disposition"]].to_string(index=False))

    affected = sorted({int(t) for t in dropped["TIC"]})
    print("\nstars affected: %d | of those keeping at least one row: %d"
          % (len(affected), sum(((website["TIC"] == t) & ~web_drop).any() for t in affected)))

    if not apply:
        print("\nDRY RUN, pass --apply")
        return

    dropped.drop(columns="_key").to_csv(os.path.join(HERE, "tois_retired.csv"), index=False)
    cat_dropped = catalog[cat_drop].copy()
    cat_dropped["disposition"] = [keys[k][0] for k in cat_dropped["_key"]]
    cat_dropped["tfop_comment"] = [keys[k][1] for k in cat_dropped["_key"]]
    cat_dropped.drop(columns="_key").to_csv(HOME + "TESS/tois_retired.csv", sep="\t", index=False)

    kept_web = website[~web_drop].drop(columns="_key")
    kept_web.to_csv(WEBSITE, index=False)
    catalog[~cat_drop].drop(columns="_key").to_csv(CATALOG, sep="\t", index=False)

    refreshed = 0
    for tic in affected:
        surviving = kept_web[kept_web["TIC"] == tic]["TOI"].tolist()
        refreshed += reindex_star(tic, surviving)

    print("\nWROTE %s (%d rows), %s (%d rows); %d figures re-emitted on %d stars"
          % (WEBSITE, len(kept_web), CATALOG, int((~cat_drop).sum()), refreshed, len(affected)))


if __name__ == "__main__":
    main("--apply" in sys.argv)
