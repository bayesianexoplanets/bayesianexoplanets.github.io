"""recompute_failed_tests.py — replace 'see_pipeline' with real per-test failure reasons.

Reconstructs the per-test diagnostics behind tois.csv's passed_all_tests column
(commits 1693db8, 3c09181, 3b73018):

  spurious fail   : spurious1 >= 0.75
  snrd fail       : snrd_pvalue <= 1e-2  AND  SNR <= 50
  ntransits fail  : num_available_transits < 3
  sharp_freq fail : has_sharp_peak AND snrd_pvalue > 0.95 AND SNR <= 50

Diagnostics come from the RecoveryRun batch0 candidate matched by period
(1% relative, integer multiples up to 50; fallback = closest candidate),
plus per-star has_sharp_peak from FullRun/stars/{tic}.csv.

The script first verifies that the recomputed pass/fail agrees with the
existing passed_all_tests column, and only then writes failed_tests.
"""
import os
import numpy as np
import pandas as pd

HERE = "/global/u2/j/julius/TESS corrected"
RECOVERY = "/pscratch/sd/j/julius/exoprob/RecoveryRun/candidates/batch0"
STARS = "/pscratch/sd/j/julius/exoprob/FullRun/stars"

SPURIOUS_MAX = 0.75
SNRD_MIN = 1e-2
SNR_OVERRIDE = 50.0
NTRANSITS_MIN = 3
HARMONICS_SNRD = 0.95
REL_TOL = 0.01
MAX_MULT = 50

_batch_cache = {}
def batch0(tic):
    tic = int(tic)
    if tic not in _batch_cache:
        path = os.path.join(RECOVERY, f"{tic}.0.csv")
        try:
            _batch_cache[tic] = pd.read_csv(path, sep='\t')
        except Exception:
            _batch_cache[tic] = None
    return _batch_cache[tic]

_sharp_cache = {}
def sharp_peak(tic):
    tic = int(tic)
    if tic not in _sharp_cache:
        path = os.path.join(STARS, f"{tic}.csv")
        try:
            df = pd.read_csv(path, sep='\t')
            _sharp_cache[tic] = bool(df['has_sharp_peak'].iloc[0]) if not df.empty else False
        except Exception:
            _sharp_cache[tic] = False
    return _sharp_cache[tic]


def match_candidate(tic, period):
    """Return the batch0 row period-matched to `period`, else the closest one."""
    df = batch0(tic)
    if df is None or df.empty or not np.isfinite(period) or period <= 0:
        return None
    cp = df['period'].to_numpy(dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        x = np.where(cp > 0, np.maximum(cp / period, period / np.where(cp > 0, cp, np.nan)), np.inf)
    n = np.rint(x)
    rel = np.abs(x - n)
    ok = (rel < REL_TOL) & (n >= 1) & (n <= MAX_MULT)
    if ok.any():
        # prefer the smallest integer multiple, then the best relative error
        score = np.where(ok, 1.0 / n - rel * 1e-3, -np.inf)
        return df.iloc[int(np.argmax(score))]
    # fallback: closest in log-period
    dist = np.abs(np.log(np.where(cp > 0, cp, np.nan)) - np.log(period))
    if np.all(~np.isfinite(dist)):
        return None
    return df.iloc[int(np.nanargmin(dist))]


def failures(row):
    cand = match_candidate(row['TIC'], row['Period'])
    if cand is None:
        return None  # no diagnostics available at all
    # the SNR>50 override uses the matched candidate's own SNR (variant testing
    # against the existing passed_all_tests column: 96.6% agreement, the best of
    # the match-selection x SNR-source combinations)
    snr = float(cand['SNR'])
    failed = []
    if float(cand['spurious1']) >= SPURIOUS_MAX:
        failed.append('spurious')
    snrd = float(cand['snrd_pvalue'])
    if snrd <= SNRD_MIN and snr <= SNR_OVERRIDE:
        failed.append('snrd')
    if int(cand['num_available_transits']) < NTRANSITS_MIN:
        failed.append('ntransits')
    if sharp_peak(row['TIC']) and snrd > HARMONICS_SNRD and snr <= SNR_OVERRIDE:
        failed.append('sharp_freq')
    return failed


def main():
    path = os.path.join(HERE, 'tois.csv')
    df = pd.read_csv(path)
    print(f"tois.csv: {len(df)} rows, {df['passed_all_tests'].sum()} marked passed")

    recomputed, reasons = [], []
    for i, row in df.iterrows():
        f = failures(row)
        if f is None:
            recomputed.append(None)
            reasons.append(None)
        else:
            recomputed.append(len(f) == 0)
            reasons.append('|'.join(f))
        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{len(df)}", flush=True)

    rec = pd.Series(recomputed)
    have = rec.notna()
    agree = (rec[have] == df.loc[have, 'passed_all_tests']).sum()
    print(f"\nDiagnostics found for {have.sum()}/{len(df)} rows")
    print(f"Agreement with existing passed_all_tests: {agree}/{have.sum()} "
          f"({100*agree/have.sum():.2f}%)")

    dis = df.loc[have][rec[have] != df.loc[have, 'passed_all_tests']]
    if len(dis):
        print(f"Disagreements: {len(dis)} (existing=True but recomputed fail: "
              f"{(dis['passed_all_tests'] == True).sum()}, reverse: {(dis['passed_all_tests'] == False).sum()})")
        print(dis[['TIC', 'TOI', 'SNR', 'passed_all_tests']].head(15).to_string())

    # Failure breakdown among rows the catalog marks as failed
    failed_mask = df['passed_all_tests'] == False
    new_reasons = []
    for i in range(len(df)):
        if not failed_mask.iloc[i]:
            new_reasons.append('')          # passed -> no failure reason
        elif reasons[i]:                     # real reasons reconstructed
            new_reasons.append(reasons[i])
        else:
            new_reasons.append('see_pipeline')  # genuinely unavailable
    df['failed_tests'] = new_reasons

    n_real = sum(1 for i in range(len(df)) if failed_mask.iloc[i] and reasons[i])
    n_left = int(failed_mask.sum()) - n_real
    print(f"\nFailed rows: {failed_mask.sum()}  -> real reasons: {n_real}, see_pipeline left: {n_left}")
    from collections import Counter
    c = Counter()
    for i in range(len(df)):
        if failed_mask.iloc[i] and reasons[i]:
            for t in reasons[i].split('|'):
                c[t] += 1
    print("Breakdown:", dict(c))

    if agree / have.sum() > 0.95:
        df.to_csv(path, index=False)
        print(f"\nWrote {path}")
    else:
        print("\nAgreement below 95% — NOT writing. Rule reconstruction needs review.")


if __name__ == '__main__':
    main()
