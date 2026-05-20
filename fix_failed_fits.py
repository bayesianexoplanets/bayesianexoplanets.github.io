#!/usr/bin/env python3
"""Fix rows where the known-planet fit FAILED but a candidate recovered the planet.

For some TOIs the per-planet ('known planet') fit converged to garbage (low or
negative SNR, broad tau) even though the blind candidate search (RecoveryRun /
FullRun batch0) clearly detected the transit. For these rows BOTH the catalog
value and the known-planet plot show the failed fit, so they are not caught by
reconcile_ocr.py (which only finds plot-vs-table divergence).

Here we detect them and, for each, switch the website's plot to the candidate
plot ({event_id}_0.png) and set the catalog Period/Phase/Tau/SNR to the
candidate's values, so table and plot both show the real detection.

Selection: table SNR < SNR_FAIL and the row is NOT a confirmed-diverged row
(plot also bad), AND a candidate exists within 2% of the catalog period with
SNR > SNR_CAND and >= MIN_TRANSITS transits.
"""
import os, sys
import numpy as np
import pandas as pd
from PIL import Image

ROOT  = "/global/u2/j/julius/TESS corrected"
CSV   = f"{ROOT}/tois.csv"
OUT   = f"{ROOT}/plots"
RUNS  = [('RecoveryRun', '/pscratch/sd/j/julius/exoprob/RecoveryRun'),
         ('FullRun',     '/pscratch/sd/j/julius/exoprob/FullRun')]

SNR_FAIL   = 5.0    # known-planet fit considered failed below this
SNR_CAND   = 10.0   # candidate considered a real recovery above this
MIN_TRANS  = 3      # NST needs >= 3 transits
PTOL       = 0.02   # candidate period within 2% of catalog period


def best_candidate(TIC, P):
    """Return (run, base, event_id, snr, period, phase, tau, ntr) or None."""
    best = None
    for run, base in RUNS:
        cp = f'{base}/candidates/batch0/{TIC}.0.csv'
        if not os.path.isfile(cp):
            continue
        try:
            c = pd.read_csv(cp, sep='\t')
        except Exception:
            continue
        cc = c[(np.abs(c.period - P) / P < PTOL) &
               (c.SNR > SNR_CAND) &
               (c.num_available_transits >= MIN_TRANS)]
        for _, r in cc.iterrows():
            png = f'{base}/plots/{TIC}/{int(r.event_id)}_0.png'
            if not os.path.isfile(png):
                continue
            if best is None or r.SNR > best[3]:
                best = (run, base, int(r.event_id), float(r.SNR),
                        float(r.period), float(r.phase), float(r.tau),
                        int(r.num_available_transits), png)
    return best


def main(apply=False):
    df = pd.read_csv(CSV); df['known_idx'] = df.groupby('TIC').cumcount()

    # confirmed-diverged rows are handled by reconcile_ocr.py; exclude them here
    diverged = set()
    rev = f"{ROOT}/plots/ocr_titles.tsv"
    if os.path.isfile(rev):
        ocr = pd.read_csv(rev, sep='\t')
        jj = df.merge(ocr[['TIC', 'known_idx', 'snr', 'period']],
                      on=['TIC', 'known_idx'], how='left')
        gate = jj['period'].notna() & (np.abs(jj['period'] - jj['Period']) / jj['Period'] < 0.05)
        divmask = gate & ((jj['snr'] - jj['SNR']).abs() > 0.15 * jj['snr'].abs() + 1.0) & (jj['snr'] > SNR_FAIL)
        diverged = set(zip(jj[divmask].TIC, jj[divmask].known_idx))

    rows = []
    for _, r in df[df.SNR < SNR_FAIL].iterrows():
        key = (r.TIC, r.known_idx)
        if key in diverged:
            continue  # plot is actually good; handled elsewhere
        bc = best_candidate(int(r.TIC), r.Period)
        if bc is None:
            continue
        run, base, ev, snr, P, phi, tau, ntr, png = bc
        rows.append(dict(TIC=int(r.TIC), known_idx=int(r.known_idx),
                         old_SNR=r.SNR, old_P=r.Period, old_phase=r.Phase, old_tau=r.Tau,
                         run=run, event_id=ev, SNR=snr, Period=P, Phase=phi, Tau=tau,
                         ntr=ntr, cand_png=png))
    res = pd.DataFrame(rows)
    res.to_csv(f"{OUT}/_failed_fits.tsv", sep='\t', index=False)
    print(f"failed-fit rows with a recovered candidate: {len(res)}")
    if len(res):
        print(res[['TIC','known_idx','old_SNR','SNR','old_P','Period','run','event_id','ntr']]
              .sort_values('SNR', ascending=False).to_string(index=False))

    if apply and len(res):
        # update plots and tois.csv + manifest
        man = pd.read_csv(f"{OUT}/manifest.tsv", sep='\t')
        base = pd.read_csv(CSV); base['known_idx'] = base.groupby('TIC').cumcount()
        backup = CSV + '.prefailed'
        import shutil
        if not os.path.isfile(backup):
            shutil.copy(CSV, backup); print(f"backed up to {backup}")
        for _, r in res.iterrows():
            # regenerate the website plot from the candidate plot
            img = Image.open(r.cand_png).convert('RGB')
            w, h = img.size
            outp = f"{OUT}/{int(r.TIC)}/{int(r.known_idx)}.jpg"
            img.resize((w // 2, h // 2), Image.LANCZOS).save(outp, 'jpeg', quality=85, optimize=True)
            # update catalog
            sel = (base.TIC == r.TIC) & (base.known_idx == r.known_idx)
            base.loc[sel, ['Period', 'Phase', 'Tau', 'SNR']] = [r.Period, r.Phase, r.Tau, r.SNR]
            # update manifest
            ms = (man.TIC == r.TIC) & (man.known_idx == r.known_idx)
            man.loc[ms, 'src'] = r.cand_png
            man.loc[ms, 'src_type'] = f"{r.run}_cand_failed_fit"
            man.loc[ms, 'snr_src'] = r.SNR
        base.drop(columns=['known_idx']).to_csv(CSV, index=False)
        man.to_csv(f"{OUT}/manifest.tsv", sep='\t', index=False)
        print(f"applied: switched {len(res)} rows to candidate plots + values")


if __name__ == "__main__":
    main(apply=('--apply' in sys.argv))
