#!/usr/bin/env python3
"""OCR the snr/period/tau (title) and phase (top-left x-axis center) from the
known-planet diagnostic plots, so the catalog can be reconciled to exactly what
each plot displays.

In planet/transit.py:plot_lc the title is
    'snr = {snr:.2}, period = {period:.4}, tau = {tau:.2}'
and the top-left subplot uses xlim = (phase - span, phase + span), so
    phase = (min_xtick + max_xtick) / 2  of the top-left subplot.
"""
import os, re, sys
import numpy as np
import pandas as pd
from PIL import Image

ROOT = "/global/u2/j/julius/TESS corrected"
CSV  = f"{ROOT}/tois.csv"

_reader = None
def reader():
    global _reader
    if _reader is None:
        import easyocr
        _reader = easyocr.Reader(['en'], gpu=True, verbose=False)
    return _reader


def _num(s):
    """Parse a float, tolerating OCR artifacts (trailing _/,, sci notation)."""
    s = s.strip().rstrip('_,').replace(' ', '')
    s = s.replace('e+0', 'e+').replace('e-0', 'e-')
    try:
        return float(s)
    except ValueError:
        return None


def ocr_plot(path):
    """Return dict(snr, period, tau, phase) or None on failure."""
    im = Image.open(path).convert('RGB')
    W, H = im.size
    r = reader()

    # --- title strip (the snr/period/tau line, second title row) ---
    # Title format is always 'snr = {:.2}, period = {:.4}, tau = {:.2}'. Parse by
    # keyword: skip any non-digit junk (=, ', spaces, commas) between a keyword
    # and its value. No allowlist (forcing letters to digits corrupted the 'tau'
    # label into spurious numbers); natural OCR reads sci-notation snr correctly.
    title = im.crop((int(0.30*W), int(0.030*H), int(0.74*W), int(0.075*H)))
    boxes = r.readtext(np.array(title), detail=1)
    boxes.sort(key=lambda b: (b[0][0][0] + b[0][2][0]) / 2)  # left-to-right
    ttxt = ' '.join(t for _, t, _ in boxes).lower()
    NUM = r'[^\d\-]*(-?\d+\.?\d*(?:\s*e\s*[+-]?\s*\d+)?)'
    def grab(key):
        m = re.search(key + NUM, ttxt)
        return _num(m.group(1)) if m else None
    snr    = grab('snr')
    period = grab('period')
    tau    = grab('tau')

    # --- top-left x-ticks -> phase = (min+max)/2 ---
    xt = im.crop((int(0.11*W), int(0.475*H), int(0.52*W), int(0.508*H)))
    ticks = []
    for _, t, conf in r.readtext(np.array(xt), detail=1):
        v = _num(t)
        if v is not None and conf > 0.3:
            ticks.append(v)
    phase = (min(ticks) + max(ticks)) / 2 if len(ticks) >= 2 else None
    # sanity: phase must lie within [0, period] (allow small slack); else misread
    if phase is not None and period is not None:
        if phase < -0.05 * period or phase > 1.05 * period:
            phase = None

    return dict(snr=snr, period=period, tau=tau, phase=phase, n_ticks=len(ticks))


def main():
    df = pd.read_csv(CSV); df['known_idx'] = df.groupby('TIC').cumcount()
    m = pd.read_csv(f"{ROOT}/plots/manifest.tsv", sep='\t')
    # known-planet plots only: src ends with /{int}.png (not _0.png candidate plots)
    m = m[m['src'].str.contains(r'/\d+\.png$', regex=True)]
    work = m[['TIC', 'known_idx', 'src']].values.tolist()
    print(f"OCR over {len(work)} known-planet plots", flush=True)

    rows = []
    for n, (TIC, i, src) in enumerate(work):
        try:
            o = ocr_plot(src)
            o.update(TIC=int(TIC), known_idx=int(i), src=src)
            rows.append(o)
        except Exception as e:
            rows.append(dict(TIC=int(TIC), known_idx=int(i), src=src,
                             snr=None, period=None, tau=None, phase=None, err=str(e)))
        if (n + 1) % 250 == 0:
            print(f"  {n+1}/{len(work)}", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(f"{ROOT}/plots/ocr_titles.tsv", sep='\t', index=False)
    print(f"wrote plots/ocr_titles.tsv ({len(out)} rows, "
          f"{out['snr'].notna().sum()} with snr, {out['phase'].notna().sum()} with phase)",
          flush=True)


if __name__ == "__main__":
    main()
