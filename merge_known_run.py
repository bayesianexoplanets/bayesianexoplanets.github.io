"""merge_known_run.py — bring the website catalog up to date with the rebuilt curated catalog and the pretty plots.

For every website row (TIC, TOI) with a row in ``TESS/tois_corrected.csv`` (rebuilt from a known-planet run):
replace Period/Phase/Tau/SNR (Duration = 2 Tau, Epoch = t_start + Phase), the errors (Laplace covariance of the
run's npz), the planet radius with its errors (radius posterior of the run), the number of valid transits and the
Singh-Maddala log10(p value) at the new SNR; refresh ``plots/{TIC}/{known_idx}.jpg`` from the pretty figure
rendered from the same catalog values (``FOLDER/pretty/{TIC}_{TOI}.png``). The pass/fail flags are kept (the
known-planet stage carries no vetting statistics). Rows without a match keep their values (reported).

Usage: python merge_known_run.py FOLDER [--apply] [TIC ...]   (dry run by default)
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

os.chdir("/global/u2/j/julius/exoplanets")
sys.path.insert(0, "/global/u2/j/julius/exoplanets/TESS corrected"); sys.path.insert(0, "/global/u2/j/julius/exoplanets")
import sm_pvalue
from constants import scratch, home
from TESS.load_tess import load_stellar_data

HERE = "/global/u2/j/julius/exoplanets/TESS corrected"
TOIS = os.path.join(HERE, "tois.csv")
PDST = os.path.join(HERE, "plots")


def npz_by_planet(folder):
    """(TIC, TOI) -> the run's npz of that planet."""
    table = {}
    for path in glob.glob(scratch + folder + "/known_planets/*.npz"):
        saved = np.load(path, allow_pickle=True)
        if "toi" in saved and np.isfinite(float(saved["toi"])):
            tic = int(os.path.basename(path).split("_")[0])
            table[(tic, round(float(saved["toi"]), 2))] = saved
    return table


def _to_jpg(src, dst, max_w=900, quality=88):
    from PIL import Image
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    image = Image.open(src).convert("RGB")
    if image.width > max_w:
        image = image.resize((max_w, int(image.height * max_w / image.width)), Image.LANCZOS)
    image.save(dst, "JPEG", quality=quality)


def main(folder, apply, tics_filter=None):
    website = pd.read_csv(TOIS)
    catalog = pd.read_csv(home + "TESS/tois_corrected.csv", sep="\t")
    catalog["TOI"] = catalog["TOI"].astype(float).round(2)
    catalog = catalog.drop_duplicates(["TIC", "TOI"]).set_index(["TIC", "TOI"])
    saved = npz_by_planet(folder)
    print(f"website rows={len(website)} | catalog rows={len(catalog)} | npz planets={len(saved)}")

    updated = plots = 0
    unmatched, missing_plots, t_start = [], [], {}
    for tic, group in website.groupby("TIC", sort=False):
        if tics_filter and int(tic) not in tics_filter:
            continue
        for known_idx, (index, row) in enumerate(group.iterrows()):
            key = (int(tic), round(float(row["TOI"]), 2))
            if key not in catalog.index:
                unmatched.append(key)
                continue
            new = catalog.loc[key]
            planet = saved.get(key)
            cov = np.asarray(planet["cov"], dtype=float) if planet is not None and "cov" in planet else np.full(6, np.nan)
            radius = np.asarray(planet["radius"], dtype=float) if planet is not None and "radius" in planet else np.full(3, np.nan)
            if int(tic) not in t_start:
                t_start[int(tic)] = float(load_stellar_data(int(tic))[0])
            snr = float(new["SNR"])
            log10p = sm_pvalue.log10p_from_grid(row["sm_sf_grid"], snr)
            if not np.isfinite(log10p):
                log10p = row["log10(p value)"]
            values = {"Period": float(new["Period"]), "Phase": float(new["Phase"]), "Tau": float(new["Tau"]), "SNR": snr,
                      "Duration": 2. * float(new["Tau"]), "Epoch": t_start[int(tic)] + float(new["Phase"]),
                      "err_Period": np.sqrt(cov[0]) if cov[0] > 0 else np.nan, "err_Epoch": np.sqrt(cov[1]) if cov[1] > 0 else np.nan,
                      "err_Duration": 2. * np.sqrt(cov[2]) if cov[2] > 0 else np.nan,
                      "Radius_planet": float(new["Radius_planet"]), "Radius_planet_errp": radius[1], "Radius_planet_errm": radius[2],
                      "Number of Valid Transits": int(new["N_valid_transits"]), "log10(p value)": log10p}
            if apply:
                for column, value in values.items():
                    website.at[index, column] = value
            updated += 1
            src = scratch + folder + f"/pretty/{int(tic)}_{key[1]:.2f}.png"
            if os.path.exists(src):
                if apply:
                    _to_jpg(src, os.path.join(PDST, str(int(tic)), f"{known_idx}.jpg"))
                plots += 1
            else:
                missing_plots.append(key)

    new_only = len(catalog) - updated - len(unmatched)
    print(f"rows updated={updated} | unmatched website rows (kept)={len(unmatched)} | catalog rows not on the website={max(new_only, 0)} | "
          f"plots refreshed={plots} | plots missing={len(missing_plots)}")
    if unmatched:
        print("sample unmatched:", unmatched[:10])
    if apply:
        website.to_csv(TOIS, index=False)
        print(f"WROTE {TOIS} ({len(website)} rows)")
    else:
        print("DRY RUN — pass --apply to write tois.csv + plots")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    tics = {int(a) for a in args[1:] if a.isdigit()}
    main(args[0], apply="--apply" in sys.argv, tics_filter=tics or None)
