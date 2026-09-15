"""build_new_candidates.py -- rebuild tois_new.csv and plots_new/ from the rerun's p-value selection.

Rows: the candidates of `results/rerun_hier_top83_verdicts.csv` (local hierarchical p < 1e-4, vetting
flags passed, harmonics of known periods and within-star duplicates removed, visually reviewed).
Values come from the candidate's own batch0 row of the rerun (Period, Phase, Tau, SNR, errors, radius
posterior, number of transits, diagnostics), the star's stellar parameters from TESS/tois_corrected.csv,
Epoch = t_start + Phase and Duration = 2 Tau as in merge_known_run.py, and the pass/fail flags from the
rules of recompute_failed_tests.py applied to the candidate's own diagnostics. `Interesting` marks the
reviewers' `convincing` verdict; `verdict`, `confidence`, `note` carry the review. The null columns are
left empty for merge_hier_pvalues.py and proposed_toi for assign_proposed_tois.py. plots_new/{TIC}/{idx}.jpg
is the downscaled `{idx}_0.png` of the rerun (the FGP-subtracted diagnostic figure). Dry run by default.
Usage: python build_new_candidates.py [--apply]
"""
import os
import shutil
import sys

import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/exoplanets/TESS corrected"
HOME = "/global/u2/j/julius/exoplanets/"
RUN = "/pscratch/sd/j/julius/exoprob/Rerun_20260910/"
sys.path.insert(0, HERE)
from merge_known_run import _to_jpg
from recompute_failed_tests import SPURIOUS_MAX, SNRD_MIN, SNR_OVERRIDE, NTRANSITS_MIN, HARMONICS_SNRD

COLUMNS = ["TIC", "TOI", "Period", "Phase", "Tau", "SNR", "Radius_planet", "Mass", "Radius", "logg", "FEH", "Teff",
           "Number of Valid Transits", "Has Visible TTVs", "log10(p value)", "μ(SNR | null)", "σ(SNR | null)", "sm_sf_grid",
           "_cand_idx", "Multiplicity", "outlier_score", "Interesting", "ood_pvalue", "Epoch", "Duration", "err_Period",
           "err_Epoch", "err_Duration", "Radius_planet_errp", "Radius_planet_errm", "passed_all_tests", "failed_tests",
           "proposed_toi", "nst_samples", "verdict", "confidence", "note"]


def failed_tests(cand, sharp):
    """recompute_failed_tests.py rules on the candidate's own diagnostics."""
    failed = []
    snr, snrd = float(cand["SNR"]), float(cand["snrd_pvalue"])
    if float(cand["spurious1"]) >= SPURIOUS_MAX:
        failed.append("spurious")
    if snrd <= SNRD_MIN and snr <= SNR_OVERRIDE:
        failed.append("snrd")
    if int(cand["num_available_transits"]) < NTRANSITS_MIN:
        failed.append("ntransits")
    if sharp and snrd > HARMONICS_SNRD and snr <= SNR_OVERRIDE:
        failed.append("sharp_freq")
    return failed


def build(apply):
    selected = pd.read_csv(HOME + "results/rerun_hier_top83_verdicts.csv", sep="\t")
    catalog = pd.read_csv(HOME + "TESS/tois_corrected.csv", sep="\t")
    stellar = catalog.drop_duplicates("TIC").set_index("TIC")
    n_known = catalog.groupby("TIC").size()
    n_new = selected.groupby("tic").size()

    rows = []
    for _, s in selected.sort_values(["tic", "index"]).iterrows():
        tic, idx = int(s["tic"]), int(s["index"])
        cand = pd.read_csv(RUN + f"candidates/batch0/{tic}.csv", sep="\t")
        cand = cand[cand["event_id"] == idx].iloc[0]
        star = pd.read_csv(RUN + f"stars/{tic}.csv", sep="\t").iloc[0]
        fails = failed_tests(cand, bool(star["has_sharp_peak"]))
        st = stellar.loc[tic] if tic in stellar.index else None
        rows.append({
            "TIC": tic, "TOI": np.nan, "Period": float(cand["period"]), "Phase": float(cand["phase"]), "Tau": float(cand["tau"]),
            "SNR": float(cand["SNR"]), "Radius_planet": float(cand["radius"]),
            "Mass": st["Mass"] if st is not None else np.nan, "Radius": st["Radius"] if st is not None else np.nan,
            "logg": st["logg"] if st is not None else np.nan, "FEH": st["FEH"] if st is not None else np.nan,
            "Teff": st["Teff"] if st is not None else np.nan,
            "Number of Valid Transits": int(cand["num_available_transits"]), "Has Visible TTVs": np.nan,
            "log10(p value)": np.nan, "μ(SNR | null)": np.nan, "σ(SNR | null)": np.nan, "sm_sf_grid": np.nan,
            "_cand_idx": idx, "Multiplicity": int(n_known.get(tic, 0) + n_new.get(tic, 0)),
            "outlier_score": np.nan, "Interesting": bool(s["verdict"] == "convincing"), "ood_pvalue": np.nan,
            "Epoch": float(star["t_start"]) + float(cand["phase"]), "Duration": 2. * float(cand["tau"]),
            "err_Period": float(cand["err_period"]), "err_Epoch": float(cand["err_phase"]), "err_Duration": 2. * float(cand["err_tau"]),
            "Radius_planet_errp": float(cand["radiusp"]), "Radius_planet_errm": float(cand["radiusm"]),
            "passed_all_tests": len(fails) == 0, "failed_tests": "|".join(fails) if fails else np.nan,
            "proposed_toi": np.nan, "nst_samples": np.nan,
            "verdict": s["verdict"], "confidence": int(s["confidence"]), "note": s["note"]})
    table = pd.DataFrame(rows)[COLUMNS]
    print(f"{len(table)} candidates on {table.TIC.nunique()} stars | passed_all_tests {int(table.passed_all_tests.sum())} | "
          f"verdicts {table.verdict.value_counts().to_dict()} | stellar params missing {int(table.Teff.isna().sum())}")

    if apply:
        table.to_csv(os.path.join(HERE, "tois_new.csv"), index=False)
        dst_root = os.path.join(HERE, "plots_new")
        if os.path.isdir(dst_root):
            shutil.rmtree(dst_root)
        for _, r in table.iterrows():
            _to_jpg(RUN + f"plots/{int(r.TIC)}/{int(r._cand_idx)}_0.png", os.path.join(dst_root, str(int(r.TIC)), f"{int(r._cand_idx)}.jpg"))
        print(f"WROTE tois_new.csv and {len(table)} plots under plots_new/")
    else:
        print("dry run (pass --apply to write tois_new.csv and plots_new/)")
    return table


if __name__ == "__main__":
    build("--apply" in sys.argv)
