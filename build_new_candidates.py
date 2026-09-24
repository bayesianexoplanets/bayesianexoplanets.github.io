"""build_new_candidates.py -- rebuild tois_new.csv and plots_new/ from the p-value selection.

Rows: the `selected` candidates of `results/rerun_hier_new_candidates_clean.csv` (global hierarchical
p < 0.01, i.e. the local p times the star's number of period bins; vetting flags passed; harmonics of
known periods and within-star duplicates removed), from Rerun_20260910 or NewRun_20260919 (`run`),
each carrying the verdict of its visual review (2026-09-15 for the first 83, 2026-09-24 for the rest).
Values come from the candidate's own batch0 row of its run (Period, Phase, Tau, SNR, errors, radius
posterior, number of transits, diagnostics), the star's stellar parameters from TESS/tois_corrected.csv,
Epoch = t_start + Phase and Duration = 2 Tau as in merge_known_run.py, and the pass/fail flags from the
rules of recompute_failed_tests.py applied to the candidate's own diagnostics, plus four further
tests (user decision 2026-09-15): `known_eb` (the star has an ExoFOP row with TESS disposition EB,
or a false-positive row whose comment names an eclipsing binary, or is in KNOWN_EB_HOSTS),
`known_transit` (a single-transit TOI of the star without a period whose ExoFOP epoch falls in one
of the candidate's transit windows), `harmonic` (the period is an integer multiple 2..8 or fraction
of a stronger candidate of the same star, within 2 % in log). Leftovers of known planets (transit or
secondary-eclipse times, or the same period) never reach this list: the selection removes them
(analysis/known_overlap.py), which replaced the `known_harmonic` period-ratio token on 2026-09-24. The
visual review is NOT a test (user
decision 2026-09-15: tests must be deterministic); it stays in the verdict/confidence/note columns. `Interesting` marks the
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
SCRATCH = "/pscratch/sd/j/julius/exoprob/"
SELECTION = HOME + "results/rerun_hier_new_candidates_clean.csv"
VERDICTS = [HOME + "results/rerun_hier_top83_verdicts.csv"] + \
    [SCRATCH + f"tmp/hier_vi2/verdicts_{k}.csv" for k in range(9)]
sys.path.insert(0, HERE)
from merge_known_run import _to_jpg
from flag_rules import (SPURIOUS_MAX, SNRD_MIN, SNR_OVERRIDE, NTRANSITS_MIN, SINGLE_TRANSIT_MIN,
                        SIGNIFICANCE_MAX, KNOWN_EB_HOSTS, EB_COMMENT, star_flags)
sys.path.insert(0, HOME)
from pipeline.post import on_a_line
from false_alarms import fold_shape

FOLD_CACHES = [HOME + "results/rerun_fold_shape_83.tsv",   # the fold statistics of the selected candidates
               HOME + "results/newcand_fold_shape_20260924.tsv"]

SHARP_LINE_MIN_FREQ = 1.5   # c/d: notched lines below this are red-noise/window leakage, not coherent oscillations
                            # (on the 83 reviewed candidates every convincing/plausible match sits at 0.10-1.0 c/d,
                            # every pulsator/EB match at >= 1.85 c/d; 2026-09-15)
SHARP_LINE_TOL = 0.035      # c/d, the notch half width (post.on_a_line default) ...
SHARP_LINE_REL_TOL = 0.02   # ... or 2 % of the harmonic's frequency, whichever is larger (TIC 279769094: 3.06 vs 3.00 c/d)

HARMONIC_RATIOS = np.array([2, 3, 4, 5, 6, 7, 8, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8])
LOG_TOL = 0.02

COLUMNS = ["TIC", "TOI", "Period", "Phase", "Tau", "SNR", "Radius_planet", "Mass", "Radius", "logg", "FEH", "Teff",
           "Number of Valid Transits", "Has Visible TTVs", "log10(p value)", "μ(SNR | null)", "σ(SNR | null)", "sm_sf_grid",
           "_cand_idx", "Multiplicity", "outlier_score", "Interesting", "ood_pvalue", "Epoch", "Duration", "err_Period",
           "err_Epoch", "err_Duration", "Radius_planet_errp", "Radius_planet_errm", "passed_all_tests", "failed_tests",
           "proposed_toi", "nst_samples", "fold_absorbed", "fold_absorbed_pure", "fold_excess", "fold_duty",
           "verdict", "confidence", "note"]


def fold_stats_for(tic, period, phase, tau, cache):
    """Fold-shape statistics of one candidate, from the cache when it holds them."""
    key = (int(tic), round(float(period), 6))
    if key in cache:
        return cache[key]
    return fold_shape.fold_shape_stats(int(tic), float(period), float(phase), float(tau))


def load_fold_cache():
    """{(TIC, period): stats} from previous runs of the fold-shape test (submissions/newcand_assets.py)."""
    cache = {}
    for path in FOLD_CACHES:
        if os.path.exists(path):
            table = pd.read_csv(path, sep="\t")
            table = table[table["error"].isna()] if "error" in table.columns else table
            cache.update({(int(r["tic"]), round(float(r["period"]), 6)):
                          {c: r[c] for c in table.columns if c.startswith("fold_")} for _, r in table.iterrows()})
    return cache


def load_verdicts():
    """{(TIC, batch0 index): (verdict, confidence, note)} from every visual review, the latest winning."""
    out = {}
    for path in VERDICTS:
        if os.path.exists(path):
            table = pd.read_csv(path, sep="\t" if path.endswith(".csv") and "top83" in path else ",")
            for r in table.itertuples():
                confidence = int(r.confidence) if pd.notna(r.confidence) else 0
                out[(int(r.tic), int(r.index))] = (r.verdict, confidence, r.note)
    return out


def failed_tests(cand, line_freqs):
    """recompute_failed_tests.py rules on the candidate's own diagnostics; sharp_freq is the pipeline's
    frequency-matched cut (post.on_a_line: the candidate frequency or its first three harmonics within the
    notch half width of a notched line at or above SHARP_LINE_MIN_FREQ), not the legacy whole-star flag."""
    failed = []
    snr, snrd = float(cand["SNR"]), float(cand["snrd_pvalue"])
    if float(cand["spurious1"]) >= SPURIOUS_MAX:
        failed.append("spurious")
    if snrd <= SNRD_MIN and snr <= SNR_OVERRIDE:
        failed.append("snrd")
    if int(cand["num_available_transits"]) < NTRANSITS_MIN:
        failed.append("ntransits")
    single = float(cand.get("single_transit_ratio", np.nan))
    if np.isfinite(single) and single < SINGLE_TRANSIT_MIN:
        failed.append("single_transit")
    freq, lines = 1. / float(cand["period"]), line_freqs if isinstance(line_freqs, str) else ""
    if any(on_a_line(k * freq, lines, n_harm=1, tol=max(SHARP_LINE_TOL, SHARP_LINE_REL_TOL * k * freq), min_line_freq=SHARP_LINE_MIN_FREQ)
           for k in (1, 2, 3)):
        failed.append("sharp_freq")
    return failed


def star_flags(exo):
    """TIC -> (is an EB host, [epochs in BTJD of its single-transit TOIs without a period])."""
    period_col = [c for c in exo.columns if c.startswith("Period")][0]
    epoch_col = [c for c in exo.columns if c.startswith("Epoch")][0]
    comments = exo["Comments"].fillna("")
    eb = (exo["TESS Disposition"] == "EB") | ((exo["TFOPWG Disposition"] == "FP") & comments.str.contains(EB_COMMENT, case=False, regex=True))
    flags = {}
    for (tic, is_eb, period, epoch) in zip(exo["TIC ID"], eb, exo[period_col], exo[epoch_col]):
        try:
            tic = int(tic)
        except (ValueError, TypeError):
            continue
        host_eb, epochs = flags.get(tic, (False, []))
        if not (np.isfinite(period) and period > 0) and np.isfinite(epoch):
            epochs = epochs + [float(epoch) - 2457000.]
        flags[tic] = (host_eb or bool(is_eb), epochs)
    for tic in KNOWN_EB_HOSTS:
        host_eb, epochs = flags.get(tic, (False, []))
        flags[tic] = (True, epochs)
    return flags


def known_periods(exo, catalog):
    """TIC -> every catalogued (corrected catalog) and ExoFOP period of the star."""
    period_col = [c for c in exo.columns if c.startswith("Period")][0]
    periods = {}
    for tic, period in list(zip(catalog["TIC"], catalog["Period"])) + list(zip(exo["TIC ID"], exo[period_col])):
        try:
            tic, period = int(tic), float(period)
        except (ValueError, TypeError):
            continue
        if np.isfinite(period) and period > 0:
            periods.setdefault(tic, []).append(period)
    return {tic: np.array(v) for tic, v in periods.items()}


def is_harmonic(period, snr, star_candidates):
    """True when a stronger candidate of the star sits at an integer multiple or fraction of the period."""
    stronger = star_candidates[star_candidates["SNR"] > snr]
    if len(stronger) == 0:
        return False
    ratio = period / stronger["period"].to_numpy()[:, None] / HARMONIC_RATIOS[None, :]
    return bool((np.abs(np.log(ratio)) < LOG_TOL).any())


def build(apply):
    selected = pd.read_csv(SELECTION, sep="\t")
    selected = selected[selected["selected"]].rename(columns={"kepid": "tic", "event_id": "index"})
    verdicts = load_verdicts()
    catalog = pd.read_csv(HOME + "TESS/tois_corrected.csv", sep="\t")
    stellar = catalog.drop_duplicates("TIC").set_index("TIC")
    n_known = catalog.groupby("TIC").size()
    n_new = selected.groupby("tic").size()
    exo = pd.read_csv(HOME + "TESS/tois.csv")
    flags = star_flags(exo)
    fold_cache = load_fold_cache()

    rows = []
    for _, s in selected.sort_values(["tic", "index"]).iterrows():
        tic, idx = int(s["tic"]), int(s["index"])
        RUN = SCRATCH + s["run"] + "/"
        verdict, confidence, note = verdicts.get((tic, idx), ("not reviewed", 0, ""))
        star_candidates = pd.read_csv(RUN + f"candidates/batch0/{tic}.csv", sep="\t")
        cand = star_candidates[star_candidates["event_id"] == idx].iloc[0]
        star = pd.read_csv(RUN + f"stars/{tic}.csv", sep="\t").iloc[0]
        fails = failed_tests(cand, star["line_freqs"])
        host_eb, single_epochs = flags.get(tic, (False, []))
        epoch, period, duration = float(star["t_start"]) + float(cand["phase"]), float(cand["period"]), 2. * float(cand["tau"])
        if host_eb:
            fails.append("known_eb")
        if any(abs((e - epoch + period / 2.) % period - period / 2.) < duration for e in single_epochs):
            fails.append("known_transit")
        if is_harmonic(period, float(cand["SNR"]), star_candidates):
            fails.append("harmonic")
        fold = fold_stats_for(tic, period, float(cand["phase"]), float(cand["tau"]), fold_cache)
        if fold_shape.fold_shape_rejects(fold):
            fails.append("fold_shape")
        log10p = float(s.get("log10p_global", np.nan))
        if np.isfinite(log10p) and log10p >= SIGNIFICANCE_MAX:
            fails.append("significance")
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
            "outlier_score": np.nan, "Interesting": bool(verdict == "convincing"), "ood_pvalue": np.nan,
            "Epoch": float(star["t_start"]) + float(cand["phase"]), "Duration": 2. * float(cand["tau"]),
            "err_Period": float(cand["err_period"]), "err_Epoch": float(cand["err_phase"]), "err_Duration": 2. * float(cand["err_tau"]),
            "Radius_planet_errp": float(cand["radiusp"]), "Radius_planet_errm": float(cand["radiusm"]),
            "passed_all_tests": len(fails) == 0, "failed_tests": "|".join(fails) if fails else np.nan,
            "proposed_toi": np.nan, "nst_samples": np.nan,
            "fold_absorbed": fold.get("fold_absorbed", np.nan), "fold_absorbed_pure": fold.get("fold_absorbed_pure", np.nan),
            "fold_excess": fold.get("fold_excess", np.nan), "fold_duty": fold.get("fold_duty", np.nan),
            "verdict": verdict, "confidence": confidence, "note": note, "_run": s["run"]})
    runs = pd.DataFrame(rows)["_run"]
    table = pd.DataFrame(rows)[COLUMNS]
    print(f"{len(table)} candidates on {table.TIC.nunique()} stars | passed_all_tests {int(table.passed_all_tests.sum())} | "
          f"verdicts {table.verdict.value_counts().to_dict()} | stellar params missing {int(table.Teff.isna().sum())}")

    if apply:
        table.to_csv(os.path.join(HERE, "tois_new.csv"), index=False)
        dst_root = os.path.join(HERE, "plots_new")
        if os.path.isdir(dst_root):
            shutil.rmtree(dst_root)
        for run, (_, r) in zip(runs, table.iterrows()):
            _to_jpg(SCRATCH + f"{run}/plots/{int(r.TIC)}/{int(r._cand_idx)}_0_pretty.png",
                    os.path.join(dst_root, str(int(r.TIC)), f"{int(r._cand_idx)}.jpg"))
        print(f"WROTE tois_new.csv and {len(table)} plots under plots_new/")
    else:
        print("dry run (pass --apply to write tois_new.csv and plots_new/)")
    return table


if __name__ == "__main__":
    build("--apply" in sys.argv)
