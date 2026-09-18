"""The pass/fail tests of the website, in one place.

Both catalog pages use these rules, so a planet and a new candidate are judged by the same
statistics with the same thresholds. Every test is deterministic and computed from the run's own
diagnostics: the visual review is a column (`verdict`, `confidence`, `note`), never a test.

A missing statistic never fails its own test, the convention of `pipeline/post.py:standard_cuts`:
a row the run did not measure keeps its other tests rather than being hidden for a blank cell.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/exoplanets/TESS corrected"
HOME = "/global/u2/j/julius/exoplanets/"
sys.path.insert(0, HOME)
from false_alarms.fold_shape import fold_shape_rejects            # noqa: E402
from pipeline.thresholds import (SPURIOUS_MAX, SNRD_MIN, SNR_OVERRIDE,   # noqa: E402,F401  one source of truth, shared with pipeline/post.py
                                 SINGLE_TRANSIT_MIN, NTRANSITS_MIN)

SIGNIFICANCE_MAX = -2.      # website only: log10 p of the SNR against the star's own null in the planet's period bin

KNOWN_EB_HOSTS = {260128333: "TOI-1338: eclipsing binary host of a circumbinary planet (Kostov et al. 2020)"}
EB_COMMENT = r"\bEB\b|eclipsing binary|\bSB2\b"


def _number(row, key):
    """Row entry as a float, NaN when absent or unparsable."""
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError):
        return np.nan
    return value


def star_flags(exo):
    """TIC -> (is an EB host, [epochs in BTJD of its single-transit TOIs without a period])."""
    period_col = [c for c in exo.columns if c.startswith("Period")][0]
    epoch_col = [c for c in exo.columns if c.startswith("Epoch")][0]
    comments = exo["Comments"].fillna("")
    eb = (exo["TESS Disposition"] == "EB") | ((exo["TFOPWG Disposition"] == "FP")
                                              & comments.str.contains(EB_COMMENT, case=False, regex=True))
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


def catalog_failures(row, eb_host=False):
    """Names of the tests a known-planet row fails.

    Parameters
    ----------
    row : mapping
        The website row merged with its re-vetting diagnostics: `SNR`, `spurious1`, `snrd_pvalue`,
        `num_available_transits`, `single_transit_ratio`, `log10(p value)`, and the fold-shape
        columns (`fold_absorbed`, `fold_applicable`).
    eb_host : bool
        Whether the star is a known eclipsing binary.

    Returns
    -------
    list of str
    """
    failed = []
    snr = _number(row, "SNR")
    spurious, snrd = _number(row, "spurious1"), _number(row, "snrd_pvalue")
    ntransits, single = _number(row, "num_available_transits"), _number(row, "single_transit_ratio")
    log10p = _number(row, "log10(p value)")

    if np.isfinite(spurious) and spurious >= SPURIOUS_MAX:
        failed.append("spurious")
    if np.isfinite(snrd) and snrd <= SNRD_MIN and np.isfinite(snr) and snr <= SNR_OVERRIDE:
        failed.append("snrd")
    if np.isfinite(single) and single < SINGLE_TRANSIT_MIN:
        failed.append("single_transit")
    if np.isfinite(ntransits) and ntransits < NTRANSITS_MIN:
        failed.append("ntransits")
    if fold_shape_rejects(row):
        failed.append("fold_shape")
    if (log10p >= SIGNIFICANCE_MAX) if np.isfinite(log10p) else not (np.isfinite(snr) and snr > 0):
        failed.append("significance")
    if eb_host:
        failed.append("known_eb")

    return failed
