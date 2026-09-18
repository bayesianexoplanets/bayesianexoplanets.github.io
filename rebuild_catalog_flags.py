"""Rebuild the known-planet pass/fail flags of the website catalog from a fresh re-vetting run.

The flags `tois.csv` carried until 2026-09-16 were written in June from the CatalogRun summary, with
two of their three statistics still broken (61 % of rows had a per-transit consistency p-value of
exactly 1.0 and the spurious-transit statistic was inflated), so of 139 flagged rows 11 reproduced
as no failure, 17 borrowed another candidate's diagnostics, and 318 rows marked as passing would
fail the same rules today.

This script rewrites ONLY `failed_tests` and `passed_all_tests`, from the re-vetting run's own
diagnostics through `flag_rules.catalog_failures`. Every other column keeps its value: the
ephemeris, SNR, errors and radii come from the known-planet run, the null columns and the p-value
from the hierarchical merge, and the figures from the pretty-plot render.

Usage: python rebuild_catalog_flags.py [--apply] [--summary PATH]   (dry run by default)
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/exoplanets/TESS corrected"
HOME = "/global/u2/j/julius/exoplanets/"
sys.path.insert(0, HERE)
from flag_rules import catalog_failures, star_flags                # noqa: E402

TOIS = os.path.join(HERE, "tois.csv")
EXO = HOME + "TESS/tois.csv"
SUMMARY = HOME + "results/catalog_vetting_20260917.tsv"
FALLBACK_SUMMARY = "/pscratch/sd/j/julius/exoprob/CatalogRun_20260917/known_vetting_summary.csv"

MAX_MISSING = 40            # rows without a re-vetting row; above this the run is incomplete, so nothing is written.
                            # 18 of these are promoted candidate rows, which the catalog re-vet never covered by construction
MAX_CONFIRMED_FOLD = 0.02   # fold_shape may reject at most this fraction of the confirmed planets
DIAGNOSTICS = ["spurious1", "snrd_pvalue", "num_available_transits", "single_transit_ratio",
               "window_ok", "fold_absorbed", "fold_applicable"]


def load_summary(path):
    """The re-vetting summary, indexed by (TIC, TOI string) as the website rows are keyed."""
    for candidate in (path, FALLBACK_SUMMARY):
        if candidate and os.path.exists(candidate):
            summary = pd.read_csv(candidate, sep="\t")
            summary["_key"] = list(zip(summary["kepid"].astype(int), summary["TOI"].astype(str)))
            return summary.drop_duplicates("_key").set_index("_key"), candidate
    raise SystemExit(f"no re-vetting summary at {path} or {FALLBACK_SUMMARY}")


def rebuild(summary_path=SUMMARY, apply=False):
    """Recompute the flag columns; returns the website table and the per-row report."""
    website = pd.read_csv(TOIS)
    summary, used = load_summary(summary_path)
    exofop = pd.read_csv(EXO)
    flags = star_flags(exofop)
    dispositions = {(int(t), str(x)): str(d) for t, x, d           # the community's current call on each signal
                    in zip(exofop["TIC ID"], exofop["TOI"], exofop["TFOPWG Disposition"])}
    print(f"website rows={len(website)} | re-vetting rows={len(summary)} ({used})")

    old_failed = website["failed_tests"].fillna("").astype(str).to_numpy()
    new_failed, missing = [], []

    for index, row in website.iterrows():
        key = (int(row["TIC"]), str(row["TOI"]))
        merged = dict(row)
        if key in summary.index:
            revet = summary.loc[[key]].iloc[0]           # .loc[tuple] would read the key as (row, column)
            merged.update({column: revet[column] for column in DIAGNOSTICS if column in summary.columns})
        else:
            missing.append(key)
            new_failed.append(old_failed[index])          # keep what the row had
            continue
        merged["TFOPWG Disposition"] = dispositions.get((int(row["TIC"]), str(row["TOI"])), "")
        eb_host = flags.get(int(row["TIC"]), (False, []))[0]
        new_failed.append("|".join(catalog_failures(merged, eb_host)))

    new_failed = np.array(new_failed)
    report = pd.DataFrame({"TIC": website["TIC"], "TOI": website["TOI"], "SNR": website["SNR"],
                           "old": old_failed, "new": new_failed})
    report["flip"] = np.where((old_failed == "") & (new_failed != ""), "pass->fail",
                              np.where((old_failed != "") & (new_failed == ""), "fail->pass", ""))

    tokens = pd.Series([t for row in new_failed for t in row.split("|") if t])
    print(f"missing re-vetting rows: {len(missing)}" + (f" (first: {missing[:5]})" if missing else ""))
    print(f"passing: {int((old_failed == '').sum())} -> {int((new_failed == '').sum())} of {len(website)}")
    print("new failure tokens:", tokens.value_counts().to_dict())
    print("old failure tokens:", pd.Series([t for row in old_failed for t in row.split("|") if t]).value_counts().to_dict())
    print(f"flips: {int((report.flip == 'pass->fail').sum())} pass->fail, "
          f"{int((report.flip == 'fail->pass').sum())} fail->pass")

    confirmed = website["TOI"].notna() & (website.get("Disposition", pd.Series("", index=website.index)) == "confirmed")
    fold_rejected = pd.Series([("fold_shape" in row.split("|")) for row in new_failed], index=website.index)
    if confirmed.any():
        rate = float((fold_rejected & confirmed).sum()) / max(int(confirmed.sum()), 1)
        print(f"fold_shape rejects {rate:.3%} of the confirmed planets")
        if rate > MAX_CONFIRMED_FOLD:
            raise SystemExit(f"fold_shape rejects more than {MAX_CONFIRMED_FOLD:.0%} of confirmed planets; not writing")

    if len(missing) > MAX_MISSING:
        raise SystemExit(f"{len(missing)} rows have no re-vetting row (limit {MAX_MISSING}); re-run those TICs first")

    if apply:
        website["failed_tests"] = np.where(new_failed == "", np.nan, new_failed)
        website["passed_all_tests"] = new_failed == ""
        website.to_csv(TOIS, index=False)
        print(f"WROTE {TOIS}")
    else:
        print("dry run (pass --apply to write the flag columns)")

    return website, report


if __name__ == "__main__":
    summary_path = sys.argv[sys.argv.index("--summary") + 1] if "--summary" in sys.argv else SUMMARY
    _, report = rebuild(summary_path, apply="--apply" in sys.argv)
    flips = report[report.flip != ""].sort_values("SNR", ascending=False)
    if len(flips):
        print("\nflipped rows (strongest first):")
        print(flips.head(40).to_string(index=False))
