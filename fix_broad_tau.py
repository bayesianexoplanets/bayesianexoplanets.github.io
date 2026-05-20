#!/usr/bin/env python3
"""Fix broad-tau fits (stellar-variability/EB contamination) using narrow candidates.

Some known-planet fits, especially on active stars (e.g. TOI-1136), lock the
matched filter onto stellar variability with an absurdly broad template
(tau >> a real transit half-duration of a few hours), giving an inflated SNR.
These evade the SNR-based failed-fit detector. Here we flag rows with
Tau > TAU_MAX days and, where a narrow candidate (tau < TAU_NARROW, >= MIN_TR
transits, within PTOL of the catalog period) recovered the real transit, switch
the website plot to the candidate plot and set Period/Phase/Tau/SNR to it.
"""
import os, sys
import numpy as np
import pandas as pd
from PIL import Image

ROOT = "/global/u2/j/julius/TESS corrected"
CSV  = f"{ROOT}/tois.csv"
OUT  = f"{ROOT}/plots"
RUNS = [('RecoveryRun', '/pscratch/sd/j/julius/exoprob/RecoveryRun'),
        ('FullRun',     '/pscratch/sd/j/julius/exoprob/FullRun')]
TAU_MAX, TAU_NARROW, MIN_TR, PTOL = 0.4, 0.3, 3, 0.02

def best_narrow(TIC, P):
    best = None
    for run, base in RUNS:
        cp = f'{base}/candidates/batch0/{TIC}.0.csv'
        if not os.path.isfile(cp):
            continue
        try: c = pd.read_csv(cp, sep='\t')
        except Exception: continue
        cc = c[(np.abs(c.period - P)/P < PTOL) & (c.tau < TAU_NARROW) &
               (c.num_available_transits >= MIN_TR)]
        for _, r in cc.iterrows():
            png = f'{base}/plots/{TIC}/{int(r.event_id)}_0.png'
            if not os.path.isfile(png): continue
            if best is None or r.SNR > best['SNR']:
                best = dict(run=run, event_id=int(r.event_id), SNR=float(r.SNR),
                            Period=float(r.period), Phase=float(r.phase), Tau=float(r.tau),
                            ntr=int(r.num_available_transits), png=png)
    return best

def main(apply=False):
    df = pd.read_csv(CSV); df['known_idx'] = df.groupby('TIC').cumcount()
    broad = df[df.Tau > TAU_MAX]
    rows, noc = [], []
    for _, r in broad.iterrows():
        bc = best_narrow(int(r.TIC), r.Period)
        if bc is None:
            noc.append(dict(TIC=int(r.TIC), known_idx=int(r.known_idx), TOI=r.TOI,
                            Period=r.Period, Tau=r.Tau, SNR=r.SNR)); continue
        rows.append(dict(TIC=int(r.TIC), known_idx=int(r.known_idx), TOI=r.TOI,
                         old_SNR=r.SNR, old_Tau=r.Tau, old_P=r.Period, **bc))
    res = pd.DataFrame(rows)
    res.to_csv(f"{OUT}/_broad_tau_fix.tsv", sep='\t', index=False)
    pd.DataFrame(noc).to_csv(f"{OUT}/_broad_tau_nocand.tsv", sep='\t', index=False)
    print(f"broad-tau rows (Tau>{TAU_MAX}d): {len(broad)}")
    print(f"  fixable with narrow candidate: {len(res)}")
    print(f"  no narrow candidate (review): {len(noc)}")

    if apply and len(res):
        man = pd.read_csv(f"{OUT}/manifest.tsv", sep='\t')
        base = pd.read_csv(CSV); base['known_idx'] = base.groupby('TIC').cumcount()
        import shutil
        bk = CSV + '.prebroadtau'
        if not os.path.isfile(bk): shutil.copy(CSV, bk); print(f"backup -> {bk}")
        for _, r in res.iterrows():
            sel = (base.TIC == r.TIC) & (base.known_idx == r.known_idx)
            base.loc[sel, ['Period','Phase','Tau','SNR']] = [r.Period, r.Phase, r.Tau, r.SNR]
            img = Image.open(r.png).convert('RGB'); w, h = img.size
            img.resize((w//2, h//2), Image.LANCZOS).save(f"{OUT}/{int(r.TIC)}/{int(r.known_idx)}.jpg", 'jpeg', quality=85, optimize=True)
            ms = (man.TIC == r.TIC) & (man.known_idx == r.known_idx)
            man.loc[ms, 'src'] = r.png; man.loc[ms, 'src_type'] = f"{r.run}_cand_broadtau"; man.loc[ms, 'snr_src'] = r.SNR
        base.drop(columns=['known_idx']).to_csv(CSV, index=False)
        man.to_csv(f"{OUT}/manifest.tsv", sep='\t', index=False)
        print(f"applied: switched {len(res)} broad-tau rows to narrow candidates")

if __name__ == "__main__":
    main(apply=('--apply' in sys.argv))
