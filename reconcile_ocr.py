#!/usr/bin/env python3
"""Reconcile tois.csv with the values displayed in the known-planet plots.

The plots and the NPZ fits diverged for some rows (NPZs were re-fit after the
plots were drawn). The website shows the plots, so the catalog must match the
plot. ocr_titles.tsv holds each plot's OCR'd snr/period/tau (title) and phase
(top-left x-axis center).

Logic per known-planet-plot row:
  - Period gate: trust OCR only if its period matches the table period within 5%
    (period is stable across re-fits; a mismatch means OCR misread or wrong plot).
  - Diverged: gate passes AND |OCR snr - table snr| > 2-sig-fig rounding.
  - Guard: if a diverged row's OCR snr claims a detection (>8) it must be backed
    by a candidate (best candidate within 2% period with SNR > 0.5*OCR snr); else
    the OCR snr is suspect -> manual review, no change.
  - Apply: overwrite Period/Phase/Tau/SNR from OCR for confirmed diverged rows.
"""
import os, sys
import numpy as np
import pandas as pd

ROOT = "/global/u2/j/julius/TESS corrected"
CSV  = f"{ROOT}/tois.csv"

_cc = {}
def best_cand_snr(TIC, P):
    if TIC not in _cc:
        d = []
        for base in ('/pscratch/sd/j/julius/exoprob/RecoveryRun',
                     '/pscratch/sd/j/julius/exoprob/FullRun'):
            cp = f'{base}/candidates/batch0/{TIC}.0.csv'
            if os.path.isfile(cp):
                try: d.append(pd.read_csv(cp, sep='\t')[['period', 'SNR']])
                except Exception: pass
        _cc[TIC] = pd.concat(d) if d else None
    c = _cc[TIC]
    if c is None: return np.nan
    cc = c[np.abs(c.period - P) / P < 0.02]
    return cc.SNR.max() if len(cc) else np.nan


def main(apply=False):
    df = pd.read_csv(CSV); df['known_idx'] = df.groupby('TIC').cumcount()
    ocr = pd.read_csv(f"{ROOT}/plots/ocr_titles.tsv", sep='\t')
    j = df.merge(ocr[['TIC', 'known_idx', 'snr', 'period', 'tau', 'phase']],
                 on=['TIC', 'known_idx'], how='left')

    ocr_rows = j['snr'].notna() | j['period'].notna()
    pgate = j['period'].notna() & (np.abs(j['period'] - j['Period']) / j['Period'] < 0.05)

    # gate-failed OCR rows -> manual review (wrong plot or OCR misread)
    review = j[ocr_rows & ~pgate]

    # diverged among gate-passing
    cand = j[pgate].copy()
    cand['snr_mismatch'] = (cand['snr'] - cand['SNR']).abs() > 0.15 * cand['snr'].abs() + 1.0
    div = cand[cand['snr_mismatch']].copy()

    # candidate-support guard: a claimed detection (OCR snr>8) needs a candidate
    div['best_cand'] = div.apply(lambda r: best_cand_snr(int(r.TIC), r.Period), axis=1)
    div['supported'] = (div['snr'] <= 8) | (div['best_cand'] > 0.5 * div['snr'])
    apply_set = div[div['supported']]
    suspect   = div[~div['supported']]

    print(f"OCR'd rows: {ocr_rows.sum()}")
    print(f"  gate pass: {pgate.sum()},  gate fail -> review: {(ocr_rows & ~pgate).sum()}")
    print(f"  diverged (gate-pass, snr mismatch): {len(div)}")
    print(f"    confirmed (apply): {len(apply_set)}")
    print(f"    suspect OCR (snr>8, no candidate) -> review: {len(suspect)}")

    pd.concat([review[['TIC','known_idx','Period','SNR','snr','period','tau','phase']],
               suspect[['TIC','known_idx','Period','SNR','snr','period','tau','phase']]]
              ).to_csv(f"{ROOT}/plots/_ocr_review.tsv", sep='\t', index=False)

    if len(apply_set):
        print("\nconfirmed diverged (table SNR -> plot SNR):")
        for _, r in apply_set.sort_values('SNR').head(40).iterrows():
            print(f"  TIC {int(r.TIC)} i{int(r.known_idx)}: "
                  f"SNR {r.SNR:7.2f}->{r['snr']:6.1f}  P {r.Period:.4f}->{r['period']}  "
                  f"tau {r.Tau:.4f}->{r['tau']}  phase {r.Phase:.3f}->{r['phase']}  "
                  f"(cand {r['best_cand']:.1f})")

    if apply and len(apply_set):
        base = pd.read_csv(CSV); base['known_idx'] = base.groupby('TIC').cumcount()
        upd = apply_set[['TIC', 'known_idx', 'snr', 'period', 'tau', 'phase']]
        merged = base.merge(upd, on=['TIC', 'known_idx'], how='left')
        m = merged['snr'].notna()
        merged.loc[m, 'SNR'] = merged.loc[m, 'snr']
        for col, src in [('Period', 'period'), ('Tau', 'tau'), ('Phase', 'phase')]:
            mm = m & merged[src].notna()
            merged.loc[mm, col] = merged.loc[mm, src]
        merged = merged.drop(columns=['known_idx', 'snr', 'period', 'tau', 'phase'])
        backup = CSV + '.preocr'
        import shutil
        if not os.path.isfile(backup):
            shutil.copy(CSV, backup); print(f"\nbacked up to {backup}")
        merged.to_csv(CSV, index=False)
        print(f"applied: updated {m.sum()} rows (Period/Phase/Tau/SNR) in tois.csv")


if __name__ == "__main__":
    main(apply=('--apply' in sys.argv))
