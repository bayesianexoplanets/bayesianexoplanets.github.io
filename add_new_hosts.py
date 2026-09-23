"""add_new_hosts.py -- bring the known planets of the new candidate hosts into both catalogs.

Carries forward the run's known-planet rows with at least NTRANSITS_MIN valid transits; no SNR cut
(significance is decided later by the null-based p-value). Writes both `TESS/tois_corrected.csv`
(pipeline catalog) and `TESS corrected/tois.csv` (website catalog); an existing (TIC, TOI) row is
updated in place rather than duplicated. Null columns are left empty for `merge_hier_pvalues.py`,
and flags are left failing (`failed_tests = 'not_vetted'`) until `rebuild_catalog_flags.py` runs.
Figures for touched stars are re-emitted from `NewRun_20260919/pretty/`.

Dry run by default; pass --apply.
Usage: python add_new_hosts.py [--apply] [--run NewRun_20260919]
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

from constants import scratch                                          # noqa: E402
from pipeline import catalog_from_known_planets as builder             # noqa: E402
from pipeline.thresholds import NTRANSITS_MIN                          # noqa: E402
from TESS.load_tess import load_stellar_data                           # noqa: E402
from merge_known_run import _to_jpg                                    # noqa: E402

WEBSITE = os.path.join(HERE, "tois.csv")
CATALOG = HOME + "TESS/tois_corrected.csv"
PLOTS = os.path.join(HERE, "plots")

def selected(folder, scratch_out):
    """The run's known-planet rows with at least NTRANSITS_MIN transits; no SNR cut."""
    table = builder.build(folder, out_path=scratch_out)
    table["TOI"] = table["TOI"].astype(float).round(2)

    enough = table["N_valid_transits"] >= NTRANSITS_MIN
    rows = table[enough].reset_index(drop=True)
    print("run %s: %d fits | >= %d transits %d on %d stars | dropped %d under %d transits"
          % (folder, len(table), NTRANSITS_MIN, len(rows), rows["TIC"].nunique(),
             int((~enough).sum()), NTRANSITS_MIN))
    print("  SNR of the carried rows: min %.2f, median %.2f, %d below 7.1 (kept; the p-value decides)"
          % (rows["SNR"].min(), rows["SNR"].median(), int((rows["SNR"] <= 7.1).sum())))
    return rows


def npz_by_planet(folder):
    """(TIC, TOI) -> the run's known-planet npz, for the covariance and the radius posterior."""
    table = {}
    for path in glob.glob(scratch + folder + "/known_planets/*.npz"):
        saved = np.load(path, allow_pickle=True)
        if "toi" in saved and np.isfinite(float(saved["toi"])):
            table[(int(os.path.basename(path).split("_")[0]), round(float(saved["toi"]), 2))] = saved
    return table


def merge_catalog(rows, apply):
    """Update or append the selected rows in the pipeline catalog."""
    catalog = pd.read_csv(CATALOG, sep="\t")
    catalog["TOI"] = catalog["TOI"].astype(float).round(2)
    known = set(zip(catalog["TIC"], catalog["TOI"]))

    keys = list(zip(rows["TIC"], rows["TOI"]))
    updates = [k for k in keys if k in known]
    fresh = rows[[k not in known for k in keys]]
    print("pipeline catalog: %d rows -> %d  (%d updated, %d appended)"
          % (len(catalog), len(catalog) + len(fresh), len(updates), len(fresh)))

    if apply:
        catalog = catalog.set_index(["TIC", "TOI"])
        for _, row in rows.iterrows():
            key = (row["TIC"], row["TOI"])
            if key in catalog.index:
                for column in catalog.columns:
                    if column in row.index:
                        catalog.at[key, column] = row[column]
        catalog = catalog.reset_index()
        catalog = pd.concat([catalog, fresh[catalog.columns]], ignore_index=True)
        catalog = catalog.sort_values(["TIC", "TOI"]).reset_index(drop=True)
        catalog.to_csv(CATALOG, sep="\t", index=False)
        print("  WROTE %s (%d rows)" % (CATALOG, len(catalog)))
    return len(fresh)


def website_values(row, saved, t_start):
    """One website row's values, in the conventions merge_known_run.py established."""
    cov = np.asarray(saved["cov"], dtype=float) if saved is not None and "cov" in saved else np.full(6, np.nan)
    radius = np.asarray(saved["radius"], dtype=float) if saved is not None and "radius" in saved else np.full(3, np.nan)
    return {
        "TIC": int(row["TIC"]), "TOI": float(row["TOI"]),
        "Period": float(row["Period"]), "Phase": float(row["Phase"]), "Tau": float(row["Tau"]),
        "SNR": float(row["SNR"]), "Duration": 2. * float(row["Tau"]), "Epoch": t_start + float(row["Phase"]),
        "err_Period": np.sqrt(cov[0]) if cov[0] > 0 else np.nan,
        "err_Epoch": np.sqrt(cov[1]) if cov[1] > 0 else np.nan,
        "err_Duration": 2. * np.sqrt(cov[2]) if cov[2] > 0 else np.nan,
        "Radius_planet": float(row["Radius_planet"]),
        "Radius_planet_errp": radius[1], "Radius_planet_errm": radius[2],
        "Mass": row["Mass"], "Radius": row["Radius"], "logg": row["logg"], "FEH": row["FEH"], "Teff": row["Teff"],
        "Number of Valid Transits": int(row["N_valid_transits"]),
        "Has Visible TTVs": False,
        "passed_all_tests": False, "failed_tests": "not_vetted",
        # blanked on updated rows too: merge_hier_pvalues.py rewrites all five from the run's own nulls
        "log10(p value)": np.nan, "nst_samples": np.nan, "sm_sf_grid": np.nan,
        "μ(SNR | null)": np.nan, "σ(SNR | null)": np.nan,
    }


def merge_website(rows, saved, apply):
    """Update or append the selected rows in the website catalog; returns the stars touched."""
    website = pd.read_csv(WEBSITE)
    website["TOI"] = website["TOI"].astype(float).round(2)
    position = {key: index for index, key in enumerate(zip(website["TIC"], website["TOI"]))}

    t_start, appended, updated = {}, [], 0
    for _, row in rows.iterrows():
        tic = int(row["TIC"])
        if tic not in t_start:
            t_start[tic] = float(load_stellar_data(tic)[0])
        values = website_values(row, saved.get((tic, row["TOI"])), t_start[tic])
        key = (tic, values["TOI"])
        if key in position:
            updated += 1
            if apply:
                for column, value in values.items():
                    website.at[website.index[position[key]], column] = value
        else:
            appended.append(values)

    print("website catalog: %d rows -> %d  (%d updated, %d appended)"
          % (len(website), len(website) + len(appended), updated, len(appended)))

    touched = sorted({int(row["TIC"]) for _, row in rows.iterrows()})
    merged = pd.concat([website, pd.DataFrame(appended)], ignore_index=True)
    merged = merged.sort_values(["TIC", "TOI"]).reset_index(drop=True)
    if apply:
        merged.to_csv(WEBSITE, index=False)
        print("  WROTE %s (%d rows)" % (WEBSITE, len(merged)))
    return merged, touched


def reindex_figures(folder, website, touched, apply):
    """Re-emit every figure of each touched star at its new per-TIC index.

    A star is reindexed only when EVERY one of its rows has a rendered source. Deleting first and
    re-emitting per row destroyed 37 existing figures when the run's pretty/ directory did not exist.
    """
    pretty = scratch + folder + "/pretty/"

    written = missing = skipped = 0
    for tic in touched:
        tois = website.loc[website["TIC"] == tic, "TOI"].tolist()
        sources = [pretty + "%d_%.2f.png" % (tic, toi) for toi in tois]
        if not all(os.path.exists(source) for source in sources):
            missing += sum(1 for source in sources if not os.path.exists(source))
            skipped += 1
            continue

        target = os.path.join(PLOTS, str(tic))
        if apply:
            for stale in glob.glob(os.path.join(target, "*.jpg")):
                os.remove(stale)
        for index, source in enumerate(sources):
            if apply:
                _to_jpg(source, os.path.join(target, "%d.jpg" % index))
            written += 1

    print("figures: %d re-emitted, %d sources missing, %d stars left untouched (render %s first)"
          % (written, missing, skipped, pretty))
    return written, missing


def main(folder, apply):
    scratch_out = scratch + "tmp/new_hosts_catalog.tsv"       # scratch, so the repo stays clean
    os.makedirs(os.path.dirname(scratch_out), exist_ok=True)
    rows = selected(folder, scratch_out)
    saved = npz_by_planet(folder)

    merge_catalog(rows, apply)
    website, touched = merge_website(rows, saved, apply)
    reindex_figures(folder, website, touched, apply)
    os.remove(scratch_out)

    if not apply:
        print("DRY RUN -- pass --apply to write both catalogs and the figures")


if __name__ == "__main__":
    run = "NewRun_20260919"
    if "--run" in sys.argv:
        run = sys.argv[sys.argv.index("--run") + 1]
    main(run, apply="--apply" in sys.argv)
