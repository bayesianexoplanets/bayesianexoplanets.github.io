#!/usr/bin/env python3
"""Build plots/{TIC}/{i}.jpg tree for the GitHub Pages site.

For each row of tois.csv, find a source PNG whose underlying (Period, Phase)
matches the catalog within 1% (period) and 1% × P (phase, modulo P), then
downscale by 50% and re-encode as JPEG q=85.

Search order per row:
  1. FullRun/known_planets/{TIC}_{i}.npz at the canonical known_idx i.
  2. FullRun/known_planets/{TIC}_{j}.npz at any other j.
  3. RecoveryRun/candidates/batch0/{TIC}.0.csv  (event_id → {N}_0.png).
  4. FullRun/candidates/batch0/{TIC}.0.csv     (event_id → {N}_0.png).
"""
import os, sys, json, re, glob
import numpy as np
import pandas as pd
from PIL import Image
from multiprocessing import Pool

FULL   = "/pscratch/sd/j/julius/exoprob/FullRun"
RECOV  = "/pscratch/sd/j/julius/exoprob/RecoveryRun"
ROOT   = "/global/u2/j/julius/TESS corrected"
OUT    = f"{ROOT}/plots"
CSV    = f"{ROOT}/tois.csv"

PTOL_REL  = 0.01   # 1% period
PHTOL_REL = 0.01   # 1% × P phase (modulo)


def phase_dist(phi1, phi2, P):
    d = abs(phi1 - phi2) % P
    return min(d, P - d)


def check_match(P_csv, phi_csv, P_src, phi_src):
    if not (np.isfinite(P_src) and P_src > 0):
        return False, np.nan, np.nan
    pdrift  = abs(P_csv - P_src) / P_csv
    phdrift = phase_dist(phi_csv, phi_src, P_csv) / P_csv
    return (pdrift < PTOL_REL and phdrift < PHTOL_REL), pdrift, phdrift


def load_npz_params(path):
    try:
        z = np.load(path, allow_pickle=True)
        p = z['params']
        return float(p[0]), float(p[1]), float(p[2]), float(z['snr'])
    except Exception:
        return None, None, None, None


def _list_known_indices(TIC):
    """Return sorted list of j for which FullRun/known_planets/{TIC}_{j}.npz exists."""
    pat = re.compile(rf'^{TIC}_(\d+)\.npz$')
    out = []
    d = f"{FULL}/known_planets"
    try:
        for f in os.listdir(d):
            m = pat.match(f)
            if m: out.append(int(m.group(1)))
    except FileNotFoundError:
        pass
    return sorted(out)


def _load_candidate_csv(path):
    """Return list of (event_id, period, phase, tau, snr). Tab-separated, 27 cols."""
    try:
        df = pd.read_csv(path, sep='\t')
        return list(zip(df['event_id'].astype(int).tolist(),
                        df['period'].tolist(),
                        df['phase'].tolist(),
                        df['tau'].tolist(),
                        df['SNR'].tolist()))
    except Exception:
        return []


def find_source(TIC, i, P_csv, phi_csv):
    """Return dict with src, src_type, P_src, phi_src, p_drift, ph_drift or None."""
    # 1a. FullRun known-planet at j = i, full (P, ϕ) match
    npz_i = f"{FULL}/known_planets/{TIC}_{i}.npz"
    png_i = f"{FULL}/plots/{TIC}/{i}.png"
    P_src_i = phi_src_i = tau_src_i = snr_src_i = None
    if os.path.isfile(npz_i) and os.path.isfile(png_i):
        P_src_i, phi_src_i, tau_src_i, snr_src_i = load_npz_params(npz_i)
        if P_src_i is not None:
            ok, pd_, phd_ = check_match(P_csv, phi_csv, P_src_i, phi_src_i)
            if ok:
                return dict(src=png_i, src_type='FullRun_known_same_idx',
                            P_src=P_src_i, phi_src=phi_src_i, tau_src=tau_src_i, snr_src=snr_src_i,
                            p_drift=pd_, ph_drift=phd_)

    # 2. FullRun known-planet at any other j (full match)
    js = [j for j in _list_known_indices(TIC) if j != i]
    best = None
    for j in js:
        npz_j = f"{FULL}/known_planets/{TIC}_{j}.npz"
        png_j = f"{FULL}/plots/{TIC}/{j}.png"
        if not os.path.isfile(png_j): continue
        P_src, phi_src, tau_src, snr_src = load_npz_params(npz_j)
        if P_src is None: continue
        ok, pd_, phd_ = check_match(P_csv, phi_csv, P_src, phi_src)
        if ok and (best is None or pd_ < best['p_drift']):
            best = dict(src=png_j, src_type='FullRun_known_other_idx',
                        P_src=P_src, phi_src=phi_src, tau_src=tau_src, snr_src=snr_src,
                        p_drift=pd_, ph_drift=phd_, matched_j=j)
    if best is not None:
        return best

    # 3. RecoveryRun candidate
    cand_csv = f"{RECOV}/candidates/batch0/{TIC}.0.csv"
    if os.path.isfile(cand_csv):
        best = None
        for event_id, P_src, phi_src, tau_src, snr_src in _load_candidate_csv(cand_csv):
            ok, pd_, phd_ = check_match(P_csv, phi_csv, P_src, phi_src)
            if not ok: continue
            png = f"{RECOV}/plots/{TIC}/{event_id}_0.png"
            if not os.path.isfile(png): continue
            if best is None or pd_ < best['p_drift']:
                best = dict(src=png, src_type='RecoveryRun_cand',
                            P_src=P_src, phi_src=phi_src, tau_src=tau_src, snr_src=snr_src,
                            p_drift=pd_, ph_drift=phd_, matched_j=event_id)
        if best is not None:
            return best

    # 4. FullRun candidate
    cand_csv = f"{FULL}/candidates/batch0/{TIC}.0.csv"
    if os.path.isfile(cand_csv):
        best = None
        for event_id, P_src, phi_src, tau_src, snr_src in _load_candidate_csv(cand_csv):
            ok, pd_, phd_ = check_match(P_csv, phi_csv, P_src, phi_src)
            if not ok: continue
            png = f"{FULL}/plots/{TIC}/{event_id}_0.png"
            if not os.path.isfile(png): continue
            if best is None or pd_ < best['p_drift']:
                best = dict(src=png, src_type='FullRun_cand',
                            P_src=P_src, phi_src=phi_src, tau_src=tau_src, snr_src=snr_src,
                            p_drift=pd_, ph_drift=phd_, matched_j=event_id)
        if best is not None:
            return best

    # 5. FullRun known-planet at any j, period-only match within 2%.
    # Last resort — catches index swaps and rows where catalog phase/tau are
    # stale relative to the pipeline's actual fit. Phase/Tau/SNR in tois.csv
    # will be overwritten with the matched NPZ's values.
    best = None
    for j in _list_known_indices(TIC):
        npz_j = f"{FULL}/known_planets/{TIC}_{j}.npz"
        png_j = f"{FULL}/plots/{TIC}/{j}.png"
        if not os.path.isfile(png_j): continue
        P_src, phi_src, tau_src, snr_src = load_npz_params(npz_j)
        if P_src is None: continue
        pdrift = abs(P_csv - P_src) / P_csv
        if pdrift >= 0.02: continue
        d = abs(phi_csv - phi_src) % P_csv
        phdrift = min(d, P_csv - d) / P_csv
        label = ('FullRun_known_same_idx_period_only' if j == i
                 else 'FullRun_known_other_idx_period_only')
        if best is None or pdrift < best['p_drift']:
            best = dict(src=png_j, src_type=label,
                        P_src=P_src, phi_src=phi_src, tau_src=tau_src, snr_src=snr_src,
                        p_drift=pdrift, ph_drift=phdrift, matched_j=j)
    if best is not None:
        return best

    # 6. RecoveryRun candidate, period-only within 2%.
    cand_csv = f"{RECOV}/candidates/batch0/{TIC}.0.csv"
    if os.path.isfile(cand_csv):
        best = None
        for event_id, P_src, phi_src, tau_src, snr_src in _load_candidate_csv(cand_csv):
            if not np.isfinite(P_src) or P_src <= 0: continue
            pdrift = abs(P_csv - P_src) / P_csv
            if pdrift >= 0.02: continue
            png = f"{RECOV}/plots/{TIC}/{event_id}_0.png"
            if not os.path.isfile(png): continue
            d = abs(phi_csv - phi_src) % P_csv
            phdrift = min(d, P_csv - d) / P_csv
            if best is None or pdrift < best['p_drift']:
                best = dict(src=png, src_type='RecoveryRun_cand_period_only',
                            P_src=P_src, phi_src=phi_src, tau_src=tau_src, snr_src=snr_src,
                            p_drift=pdrift, ph_drift=phdrift, matched_j=event_id)
        if best is not None:
            return best

    # 7. FullRun candidate, period-only within 2%.
    cand_csv = f"{FULL}/candidates/batch0/{TIC}.0.csv"
    if os.path.isfile(cand_csv):
        best = None
        for event_id, P_src, phi_src, tau_src, snr_src in _load_candidate_csv(cand_csv):
            if not np.isfinite(P_src) or P_src <= 0: continue
            pdrift = abs(P_csv - P_src) / P_csv
            if pdrift >= 0.02: continue
            png = f"{FULL}/plots/{TIC}/{event_id}_0.png"
            if not os.path.isfile(png): continue
            d = abs(phi_csv - phi_src) % P_csv
            phdrift = min(d, P_csv - d) / P_csv
            if best is None or pdrift < best['p_drift']:
                best = dict(src=png, src_type='FullRun_cand_period_only',
                            P_src=P_src, phi_src=phi_src, tau_src=tau_src, snr_src=snr_src,
                            p_drift=pdrift, ph_drift=phdrift, matched_j=event_id)
        if best is not None:
            return best

    return None


def process_row(args):
    TIC, i, P_csv, phi_csv = args
    TIC, i = int(TIC), int(i)
    try:
        m = find_source(TIC, i, P_csv, phi_csv)
        if m is None:
            return dict(status='unmatched', TIC=TIC, known_idx=i,
                        P_csv=P_csv, phi_csv=phi_csv)
        img = Image.open(m['src']).convert('RGB')
        w, h = img.size
        img = img.resize((w // 2, h // 2), Image.LANCZOS)
        out_dir = f"{OUT}/{TIC}"
        os.makedirs(out_dir, exist_ok=True)
        out_path = f"{out_dir}/{i}.jpg"
        img.save(out_path, 'jpeg', quality=85, optimize=True)
        size_kb = os.path.getsize(out_path) / 1024
        return dict(status='ok', TIC=TIC, known_idx=i,
                    src=m['src'], src_type=m['src_type'],
                    P_csv=P_csv, P_src=m['P_src'],
                    phi_csv=phi_csv, phi_src=m['phi_src'],
                    tau_src=m.get('tau_src'),
                    snr_src=m.get('snr_src'),
                    p_drift=m['p_drift'], ph_drift=m['ph_drift'],
                    out=out_path, size_kb=size_kb)
    except Exception as e:
        return dict(status='error', TIC=TIC, known_idx=i, error=str(e),
                    P_csv=P_csv, phi_csv=phi_csv)


def main(limit=None, nproc=32):
    df = pd.read_csv(CSV)
    df['known_idx'] = df.groupby('TIC').cumcount()
    rows = df[['TIC', 'known_idx', 'Period', 'Phase']].astype(
        {'TIC': int, 'known_idx': int}).values.tolist()
    if limit is not None:
        rows = rows[:limit]
    print(f"processing {len(rows)} rows with {nproc} workers", flush=True)

    os.makedirs(OUT, exist_ok=True)

    with Pool(nproc) as pool:
        results = list(pool.imap_unordered(process_row, rows, chunksize=16))

    ok        = [r for r in results if r['status'] == 'ok']
    unmatched = [r for r in results if r['status'] == 'unmatched']
    errors    = [r for r in results if r['status'] == 'error']
    print(f"ok: {len(ok)}, unmatched: {len(unmatched)}, errors: {len(errors)}",
          flush=True)

    # manifest
    mdf = pd.DataFrame(ok)
    mdf = mdf.drop(columns=['status']).sort_values(['TIC', 'known_idx'])
    mdf.to_csv(f"{OUT}/manifest.tsv", sep='\t', index=False)

    udf = pd.DataFrame(unmatched + errors)
    udf.to_csv(f"{OUT}/unmatched.tsv", sep='\t', index=False)

    # summary
    if len(ok):
        print(f"source_type distribution:\n{mdf['src_type'].value_counts().to_string()}",
              flush=True)
        print(f"size: total {mdf['size_kb'].sum() / 1024:.1f} MB, "
              f"median {mdf['size_kb'].median():.1f} KB", flush=True)
        print(f"max p_drift: {mdf['p_drift'].max():.4%}, "
              f"max ph_drift: {mdf['ph_drift'].max():.4%}", flush=True)

    # Reconcile tois.csv with the plots the website shows: overwrite
    # Period/Phase/Tau/SNR of every matched row with the values from the plot it
    # links to, so the table value equals what is printed in the linked plot.
    if limit is None:
        backup = CSV + '.preplot_sync'
        if not os.path.isfile(backup):
            import shutil
            shutil.copy(CSV, backup)
            print(f"\nbacked up tois.csv to {backup}", flush=True)

        base = pd.read_csv(CSV)
        base['known_idx'] = base.groupby('TIC').cumcount()
        src = mdf[['TIC', 'known_idx', 'P_src', 'phi_src', 'tau_src', 'snr_src']]
        merged = base.merge(src, on=['TIC', 'known_idx'], how='left')
        mask = merged['P_src'].notna()
        merged.loc[mask, 'Period'] = merged.loc[mask, 'P_src']
        merged.loc[mask, 'Phase']  = merged.loc[mask, 'phi_src']
        merged.loc[mask, 'Tau']    = merged.loc[mask, 'tau_src']
        if 'SNR' in merged.columns:
            smask = mask & merged['snr_src'].notna()
            merged.loc[smask, 'SNR'] = merged.loc[smask, 'snr_src']
        merged = merged.drop(columns=['known_idx', 'P_src', 'phi_src', 'tau_src', 'snr_src'])
        merged.to_csv(CSV, index=False)
        print(f"reconciled tois.csv: Period/Phase/Tau/SNR synced to plot source "
              f"for {mask.sum()} rows", flush=True)


if __name__ == "__main__":
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 32
    main(limit=limit, nproc=nproc)
