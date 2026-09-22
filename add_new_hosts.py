"""add_new_hosts.py -- bring the known planets of the new candidate hosts into both catalogs.

`NewRun_20260919` searched 408 hosts that carried no row of our own. Its known-planet stage produced
360 fits; the rows we publish are the **305** that are detected AND have at least three valid transits
(user decisions of 2026-09-21: a planet we did not recover ourselves is not ours to list, and a row with
two transits or fewer has no period-local null, so `pipeline.thresholds.NTRANSITS_MIN` removes it --
the same rule `remove_few_transits.py` applied to the existing catalog).

Detection is `SNR > SNR_DETECTED`. That constant is the campaign's threshold and it decides roughly
twenty-five rows here, so it is reported explicitly rather than applied silently.

Both trees are written, in the order the rest of the close-out expects:

  1. `TESS/tois_corrected.csv`  -- the pipeline catalog and the source of truth.
  2. `TESS corrected/tois.csv`  -- the website catalog, with the same value conventions
     `merge_known_run.py` uses (Duration = 2 Tau, Epoch = t_start + Phase, errors from the Laplace
     covariance of the run's npz, radius errors from its radius posterior).

A (TIC, TOI) already present is UPDATED in place rather than duplicated; 36 of the 305 are already in
the pipeline catalog and 32 in the website catalog, because the run covered a few stars that had a row
already. The newer fit wins, both catalogs coming from the same code version.

The null columns (`nst_samples`, mu/sigma, `sm_sf_grid`, `log10(p value)`) are left empty for
`merge_hier_pvalues.py`, and the flags are left failing (`failed_tests = 'not_vetted'`) so that a row
which never reaches `rebuild_catalog_flags.py` is conspicuous instead of being published as passing.

Figures are re-emitted per affected star, since the frontend indexes them by per-TIC row order; they
come from `NewRun_20260919/pretty/` and are simply reported missing until that render runs.

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

SNR_DETECTED = 7.1      # the campaign's detection threshold; see the module docstring


def selected(folder, scratch_out):
    """The run's publishable known-planet rows: detected, and with at least NTRANSITS_MIN transits."""
    table = builder.build(folder, out_path=scratch_out)
    table["TOI"] = table["TOI"].astype(float).round(2)

    detected = table["SNR"] > SNR_DETECTED
    enough = table["N_valid_transits"] >= NTRANSITS_MIN
    print("run %s: %d fits | detected %d | >= %d transits %d | publishable %d on %d stars"
          % (folder, len(table), int(detected.sum()), NTRANSITS_MIN, int(enough.sum()),
             int((detected & enough).sum()), table.loc[detected & enough, "TIC"].nunique()))
    print("  dropped: %d not detected, %d detected but under %d transits"
          % (int((~detected).sum()), int((detected & ~enough).sum()), NTRANSITS_MIN))
    return table[detected & enough].reset_index(drop=True)


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
        # blanked on updated rows too: the SNR has moved, so an old grid and its p-value no longer
        # belong together, and merge_hier_pvalues.py rewrites all five from the run's own nulls
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
    """Re-emit every figure of each touched star at its new per-TIC index."""
    pretty = scratch + folder + "/pretty/"

    written = missing = 0
    for tic in touched:
        tois = website.loc[website["TIC"] == tic, "TOI"].tolist()
        target = os.path.join(PLOTS, str(tic))
        if apply:
            for stale in glob.glob(os.path.join(target, "*.jpg")):
                os.remove(stale)
        for index, toi in enumerate(tois):
            source = pretty + "%d_%.2f.png" % (tic, toi)
            if os.path.exists(source):
                if apply:
                    _to_jpg(source, os.path.join(target, "%d.jpg" % index))
                written += 1
            else:
                missing += 1
    print("figures: %d available, %d still to render from %s" % (written, missing, pretty))
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
