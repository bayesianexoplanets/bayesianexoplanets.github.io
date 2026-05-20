#!/usr/bin/env python3
"""Unified reconciliation of tois.csv with the plots the website shows.

Per known-planet-plot row, using OCR of the plot (ocr_titles.tsv: snr/period/tau
= title, phase = top-left x-axis center) and the candidate CSVs:

  gate = OCR period within 5% of catalog period (else OCR unreliable -> review).

  Decision (first match wins):
   1) DIVERGED  - gate AND plot snr >= SNR_OK AND |plot snr - table snr| large:
        the plot already shows the real planet but the catalog (stale NPZ) does
        not. Trust the plot -> set Period/Phase/Tau/SNR from OCR (keep plot).
        Guard: a big jump with no supporting candidate -> review (OCR misread).
   2) FAILED_FIT - table snr < SNR_FAIL AND the plot does NOT show a good fit
        (plot snr < SNR_OK or gate fail) AND a candidate (snr>SNR_CAND, >=3
        transits, within 2% period) recovered the planet:
        switch the website plot to the candidate plot and set values to it.
   3) MINOR_SYNC - gate AND |plot snr - table snr| moderate: table<-OCR (cosmetic).
   4) leave unchanged.
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
SNR_OK, SNR_FAIL, SNR_CAND, MIN_TR, PTOL = 8.0, 8.0, 10.0, 3, 0.02

def diverged_amount(osnr, snr):
    return abs(osnr - snr) > 0.15 * abs(osnr) + 1.0

_cc = {}
def best_candidate(TIC, P):
    if TIC not in _cc:
        d = []
        for run, base in RUNS:
            cp = f'{base}/candidates/batch0/{TIC}.0.csv'
            if os.path.isfile(cp):
                try:
                    c = pd.read_csv(cp, sep='\t'); c['run'] = run; c['base'] = base
                    d.append(c)
                except Exception: pass
        _cc[TIC] = pd.concat(d) if d else None
    c = _cc[TIC]
    if c is None: return None
    cc = c[(np.abs(c.period - P)/P < PTOL) & (c.SNR > SNR_CAND) &
           (c.num_available_transits >= MIN_TR)]
    if not len(cc): return None
    r = cc.sort_values('SNR', ascending=False).iloc[0]
    png = f"{r.base}/plots/{TIC}/{int(r.event_id)}_0.png"
    if not os.path.isfile(png): return None
    return dict(run=r.run, event_id=int(r.event_id), SNR=float(r.SNR),
                Period=float(r.period), Phase=float(r.phase), Tau=float(r.tau),
                ntr=int(r.num_available_transits), png=png)

def best_cand_snr(TIC, P):
    bc = best_candidate(TIC, P)
    return bc['SNR'] if bc else np.nan

def classify(row, o):
    snr = row.SNR
    osnr, op, otau, oph = o['snr'], o['period'], o['tau'], o['phase']
    gate = (op is not None and not np.isnan(op) and abs(op - row.Period)/row.Period < 0.05)
    bc = best_candidate(int(row.TIC), row.Period)
    # 1) DIVERGED: plot shows the planet (snr >= SNR_OK) but the catalog (stale
    #    NPZ) is lower/negative. Restrict to osnr > snr: exponent-drop OCR errors
    #    make osnr SMALLER than the (correct) table value, so they cannot enter
    #    here -- they fall through to minor_sync (ignored).
    if gate and osnr is not None and osnr >= SNR_OK and osnr > snr and diverged_amount(osnr, snr):
        # exponent-drop sanity vs candidate: OCR much below candidate -> suspect
        if bc is not None and osnr < 0.4 * bc['SNR']:
            return ('review_suspect', None)
        if abs(osnr - snr) > 5 and bc is None:
            return ('review_suspect', None)
        return ('diverged', dict(Period=op, Phase=oph, Tau=otau, SNR=osnr))
    # 2) FAILED_FIT: plot does not show a good fit and a candidate recovered it
    plot_bad = (not gate) or (osnr is None) or (osnr < SNR_OK)
    if snr < SNR_FAIL and plot_bad and bc is not None:
        return ('failed_fit', bc)
    # 3) MINOR_SYNC: small low-snr mismatches; mostly OCR exponent-drop noise
    #    where the table is already correct -> IGNORE (do not apply).
    if gate and osnr is not None and diverged_amount(osnr, snr):
        return ('minor_sync', dict(Period=op, Phase=oph, Tau=otau, SNR=osnr))
    return ('ok', None)

def main(apply=False):
    df = pd.read_csv(CSV); df['known_idx'] = df.groupby('TIC').cumcount()
    ocr = pd.read_csv(f"{OUT}/ocr_titles.tsv", sep='\t')
    man = pd.read_csv(f"{OUT}/manifest.tsv", sep='\t')
    od = {(r.TIC, r.known_idx): dict(snr=r.snr, period=r.period, tau=r.tau, phase=r.phase)
          for r in ocr.itertuples()}

    # only known-planet-plot rows are OCR'd; others (candidate-sourced) already consistent
    kp = set(zip(man[man['src'].str.contains(r'/\d+\.png$', regex=True)].TIC,
                 man[man['src'].str.contains(r'/\d+\.png$', regex=True)].known_idx))

    buckets = {k: [] for k in ['diverged', 'failed_fit', 'minor_sync', 'review_suspect', 'ok']}
    for r in df.itertuples():
        key = (r.TIC, r.known_idx)
        if key not in kp or key not in od:
            continue
        cls, payload = classify(r, od[key])
        buckets[cls].append((r.TIC, r.known_idx, r.SNR, payload))

    for k in buckets:
        print(f"{k:16}: {len(buckets[k])}")

    if apply:
        base = pd.read_csv(CSV); base['known_idx'] = base.groupby('TIC').cumcount()
        import shutil
        bk = CSV + '.prereconcile'
        if not os.path.isfile(bk): shutil.copy(CSV, bk); print(f"backup -> {bk}")
        n_plot = 0
        for cls in ('diverged',):   # minor_sync intentionally ignored (OCR snr noise)
            for TIC, ki, _, p in buckets[cls]:
                sel = (base.TIC == TIC) & (base.known_idx == ki)
                for col in ('Period', 'Phase', 'Tau', 'SNR'):
                    if p[col] is not None and not (isinstance(p[col], float) and np.isnan(p[col])):
                        base.loc[sel, col] = p[col]
        for TIC, ki, _, p in buckets['failed_fit']:
            sel = (base.TIC == TIC) & (base.known_idx == ki)
            base.loc[sel, ['Period', 'Phase', 'Tau', 'SNR']] = [p['Period'], p['Phase'], p['Tau'], p['SNR']]
            img = Image.open(p['png']).convert('RGB'); w, h = img.size
            img.resize((w//2, h//2), Image.LANCZOS).save(f"{OUT}/{TIC}/{ki}.jpg", 'jpeg', quality=85, optimize=True)
            ms = (man.TIC == TIC) & (man.known_idx == ki)
            man.loc[ms, 'src'] = p['png']; man.loc[ms, 'src_type'] = f"{p['run']}_cand_failed_fit"; man.loc[ms, 'snr_src'] = p['SNR']
            n_plot += 1
        base.drop(columns=['known_idx']).to_csv(CSV, index=False)
        man.to_csv(f"{OUT}/manifest.tsv", sep='\t', index=False)
        print(f"applied: {len(buckets['diverged'])} diverged value syncs, "
              f"{n_plot} plot switches (failed_fit); minor_sync ignored")
    else:
        # dump review + samples
        rev = buckets['review_suspect']
        pd.DataFrame([(t, k, s) for t, k, s, _ in rev], columns=['TIC','known_idx','SNR']).to_csv(
            f"{OUT}/_review_suspect.tsv", sep='\t', index=False)

if __name__ == "__main__":
    main(apply=('--apply' in sys.argv))
