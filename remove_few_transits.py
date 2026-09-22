"""Remove catalog rows with two or fewer valid transits, which our null method cannot validate.

A row's significance comes from a period-local null: the distribution of the best spurious signal in
the grid bin holding its period. That null needs three populated transit windows before it means
anything, which is why `pipeline/thresholds.py` carries `NTRANSITS_MIN = 3` as a STRUCTURAL cut. A
one- or two-transit signal therefore has no p-value we can stand behind, whatever its SNR.

Until now such rows stayed in the catalog carrying a `num_transits` failure flag. User decision of
2026-09-21: they are removed outright, as the retired false positives were, since we do not publish a
row we cannot independently verify. The same rule applies to the new candidate hosts at merge time,
where it drops 39 of the 360 fits (see `add_new_hosts.py`), and it subsumes every planet whose period
lies beyond its star's null grid entirely - all ten have two transits or fewer, because the grid stops
where the third populated window stops. It is the ONLY selection applied at merge time: nothing is
dropped on SNR, and significance is left to the p-value (user decision, 2026-09-21).

Rows are MOVED to `tois_few_transits.csv` in both trees, never deleted outright. The frontend assigns
a figure index by per-TIC row order, so every affected star has its figures re-emitted at the new
indices, exactly as `remove_retired.py` does.

Dry run by default; pass --apply.
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/exoplanets/TESS corrected"
HOME = "/global/u2/j/julius/exoplanets/"
sys.path.insert(0, HERE)
sys.path.insert(0, HOME)
from merge_known_run import _to_jpg                                   # noqa: E402
from constants import scratch                                         # noqa: E402
from pipeline.thresholds import NTRANSITS_MIN                         # noqa: E402

WEBSITE = os.path.join(HERE, "tois.csv")
CATALOG = HOME + "TESS/tois_corrected.csv"
PLOTS = os.path.join(HERE, "plots")
PRETTY = scratch + "KnownRun_20260916/pretty/"

WEB_COLUMN = "Number of Valid Transits"
CAT_COLUMN = "N_valid_transits"


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
    website = pd.read_csv(WEBSITE)
    catalog = pd.read_csv(CATALOG, sep="\t")

    web_counts = pd.to_numeric(website[WEB_COLUMN], errors="coerce")
    cat_counts = pd.to_numeric(catalog[CAT_COLUMN], errors="coerce")
    web_drop = web_counts < NTRANSITS_MIN
    cat_drop = cat_counts < NTRANSITS_MIN

    print("removing rows with fewer than %d valid transits" % NTRANSITS_MIN)
    print("  website rows : %d of %d  (0 transits %d, 1 %d, 2 %d)"
          % (web_drop.sum(), len(website), (web_counts == 0).sum(),
             (web_counts == 1).sum(), (web_counts == 2).sum()))
    print("  catalog rows : %d of %d  (0 transits %d, 1 %d, 2 %d)"
          % (cat_drop.sum(), len(catalog), (cat_counts == 0).sum(),
             (cat_counts == 1).sum(), (cat_counts == 2).sum()))
    if "passed_all_tests" in website.columns:
        passing = website["passed_all_tests"].astype(str).str.lower().isin(["true", "1"])
        print("  of the website rows removed, currently shown as passing: %d"
              % int((web_drop & passing).sum()))

    dropped = website[web_drop]
    if len(dropped):
        print("\nhighest-SNR removals:")
        print(dropped.nlargest(8, "SNR")[["TIC", "TOI", "Period", "SNR", WEB_COLUMN]].to_string(index=False))

    affected = sorted({int(t) for t in dropped["TIC"]})
    print("\nstars affected: %d | of those keeping at least one row: %d"
          % (len(affected), sum(((website["TIC"] == t) & ~web_drop).any() for t in affected)))

    if not apply:
        print("\nDRY RUN, pass --apply")
        return

    dropped.to_csv(os.path.join(HERE, "tois_few_transits.csv"), index=False)
    catalog[cat_drop].to_csv(HOME + "TESS/tois_few_transits.csv", sep="\t", index=False)

    kept_web = website[~web_drop]
    kept_web.to_csv(WEBSITE, index=False)
    catalog[~cat_drop].to_csv(CATALOG, sep="\t", index=False)

    refreshed = 0
    for tic in affected:
        surviving = kept_web[kept_web["TIC"] == tic]["TOI"].tolist()
        refreshed += reindex_star(tic, surviving)

    print("\nWROTE %s (%d rows), %s (%d rows); %d figures re-emitted on %d stars"
          % (WEBSITE, len(kept_web), CATALOG, int((~cat_drop).sum()), refreshed, len(affected)))


if __name__ == "__main__":
    main("--apply" in sys.argv)
