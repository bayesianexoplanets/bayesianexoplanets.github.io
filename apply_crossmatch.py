"""Apply the prior-art cross-match to the new-candidate list: known objects to the catalog, blends off the site.

The vetting tests ask whether a signal is real, never whether it is new. A cross-match against ExoFOP TOIs
and CTOIs, the NASA Exoplanet Archive, Gaia DR3 eclipsing binaries and variables, VSX, the TESS EB catalog
and SIMBAD (results/newcand_prior_verdicts.tsv for the 2026-09-18 list, tmp/hier_vi2/crossmatch_*.tsv for
the 2026-09-24 additions) sorts each candidate (user decision 2026-09-18, reaffirmed 2026-09-24):

- known_planet, known_toi, known_ctoi, known_toi_harmonic -> PROMOTED to tois.csv, appended at the end
  (the frontend's per-TIC figure index), with a `known_as` column naming the object: a confirmed planet
  absent from ExoFOP's TOI table, a TOI (most often one without a period that our search has measured) or
  a community TOI;
- known_eb, known_variable, blend_neighbour -> DROPPED to tois_new_dropped.csv with the identification:
  the signal is not a new planet of the target; so is a known object with fewer than three transits
  (`ntransits`), which the catalog removes rather than flags (user decision 2026-09-21);
- novel, unclear, or no cross-match row -> kept.
Names have commas stripped: the site's CSV parser splits on them.

Dry run by default; pass --apply.
"""
import os
import shutil
import sys

import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/exoplanets/TESS corrected"
HOME = "/global/u2/j/julius/exoplanets/"
TOIS = os.path.join(HERE, "tois.csv")
TOIS_NEW = os.path.join(HERE, "tois_new.csv")
DROPPED = os.path.join(HERE, "tois_new_dropped.csv")
PLOTS, PLOTS_NEW = os.path.join(HERE, "plots"), os.path.join(HERE, "plots_new")
SOURCES = [HOME + "results/newcand_prior_verdicts.tsv"] + \
    [f"/pscratch/sd/j/julius/exoprob/tmp/hier_vi2/crossmatch_{k}.tsv" for k in "abcd"]

PROMOTE = {"known_planet", "known_toi", "known_ctoi", "known_toi_harmonic"}
DROP = {"known_eb", "known_variable", "blend_neighbour"}
TOL = 0.001                        # relative period agreement of a cross-match row and a candidate


def load_crossmatch():
    """TIC -> list of (period, prior_status, prior_name) over every cross-match table."""
    out = {}
    for path in SOURCES:
        if not os.path.exists(path):
            print(f"  (missing {path})")
            continue
        table = pd.read_csv(path, sep="\t")
        for r in table.itertuples():
            out.setdefault(int(r.TIC), []).append((float(r.period), str(r.prior_status), str(r.prior_name)))
    return out


def status_of(tic, period, crossmatch):
    """(prior_status, prior_name) of the candidate, ('none', '') without a matching row."""
    for p, status, name in crossmatch.get(int(tic), []):
        if abs(p / period - 1.) < TOL:
            return status, name
    return "none", ""


def main(apply):
    website = pd.read_csv(TOIS)
    candidates = pd.read_csv(TOIS_NEW)
    crossmatch = load_crossmatch()
    status = [status_of(t, p, crossmatch) for t, p in zip(candidates["TIC"], candidates["Period"])]
    candidates["_status"] = [s for s, _ in status]
    candidates["_name"] = [n.replace(",", ";") for _, n in status]
    print("cross-match status of the %d candidates: %s" % (len(candidates), candidates["_status"].value_counts().to_dict()))

    too_few = candidates["failed_tests"].fillna("").astype(str).str.split("|").apply(lambda t: "ntransits" in t)
    promote = candidates[candidates["_status"].isin(PROMOTE) & ~too_few]
    drop = candidates[candidates["_status"].isin(DROP) | (candidates["_status"].isin(PROMOTE) & too_few)]
    keep = candidates.drop(index=promote.index.union(drop.index))

    rows, copies = [], []
    for _, c in promote.iterrows():
        tic = int(c["TIC"])
        index_here = int((website["TIC"] == tic).sum()) + sum(1 for t, _, _ in copies if t == tic)
        row = {col: c[col] for col in website.columns if col in c.index}
        row.update(TIC=tic, known_as=c["_name"])
        rows.append(row)
        copies.append((tic, int(c["_cand_idx"]), index_here))
        print("  promote TIC %-10d P=%9.5f SNR=%6.2f  %-16s %s" % (tic, c["Period"], c["SNR"], c["_status"], c["_name"]))
    for _, c in drop.iterrows():
        print("  drop    TIC %-10d P=%9.5f SNR=%6.2f  %-16s %s" % (int(c["TIC"]), c["Period"], c["SNR"], c["_status"], c["_name"]))
    print("candidates %d -> promoted %d, dropped %d, kept %d" % (len(candidates), len(promote), len(drop), len(keep)))

    if not apply:
        print("\nDRY RUN, pass --apply")
        return

    pd.concat([website, pd.DataFrame(rows)], ignore_index=True).to_csv(TOIS, index=False)
    keep.drop(columns=["_status", "_name"]).to_csv(TOIS_NEW, index=False)
    previous = pd.read_csv(DROPPED) if os.path.exists(DROPPED) else pd.DataFrame()
    reason = drop["_status"] + np.where(drop["_status"].isin(PROMOTE), "|ntransits", "")
    dropped = drop.assign(drop_reason=reason, prior_name=drop["_name"]).drop(columns=["_status", "_name"])
    pd.concat([previous, dropped], ignore_index=True).drop_duplicates(["TIC", "Period"], keep="last").to_csv(DROPPED, index=False)

    for tic, cand_idx, index_here in copies:
        source = os.path.join(PLOTS_NEW, str(tic), "%d.jpg" % cand_idx)
        os.makedirs(os.path.join(PLOTS, str(tic)), exist_ok=True)
        if os.path.exists(source):
            shutil.copyfile(source, os.path.join(PLOTS, str(tic), "%d.jpg" % index_here))
    for tic in set(promote["TIC"].astype(int)) | set(drop["TIC"].astype(int)):
        stale = os.path.join(PLOTS_NEW, str(tic))
        if os.path.isdir(stale) and not (keep["TIC"].astype(int) == tic).any():
            shutil.rmtree(stale)
    print("\nWROTE %s (%d rows), %s (%d rows), %s (+%d rows)" % (TOIS, len(website) + len(rows), TOIS_NEW, len(keep),
                                                                 DROPPED, len(drop)))


if __name__ == "__main__":
    main("--apply" in sys.argv)
