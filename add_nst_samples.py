"""add_nst_samples.py — add nst_samples column to tois.csv and tois_new.csv.

For each TIC reads batch1..10 NST-run CSVs, takes the max SNR per run (the
null max-SNR for that simulation), and serialises as a pipe-separated string.
Rows for the same TIC share the same value. If fewer than 3 batches are found
the column is left empty and the frontend falls back to Gaussian-only CDF.
"""
import os
import pandas as pd

HERE = "/global/u2/j/julius/TESS corrected"
FULLRUN = "/pscratch/sd/j/julius/exoprob/FullRun/candidates"
NST_BATCHES = list(range(1, 11))


def nst_samples_for_tic(tic):
    samples = []
    for b in NST_BATCHES:
        path = os.path.join(FULLRUN, f"batch{b}", f"{tic}.0.csv")
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep='\t')
            if len(df) > 0 and 'SNR' in df.columns:
                samples.append(float(df['SNR'].max()))
        except Exception:
            pass
    return '|'.join(f"{s:.6f}" for s in samples) if len(samples) >= 3 else ''


for csv_name in ['tois.csv', 'tois_new.csv']:
    csv_path = os.path.join(HERE, csv_name)
    df = pd.read_csv(csv_path)

    # Build per-TIC cache to avoid redundant disk reads
    tics = df['TIC'].dropna().unique()
    print(f"\n{csv_name}: {len(tics)} unique TICs, {len(df)} rows")
    cache = {}
    for i, tic in enumerate(tics):
        cache[int(tic)] = nst_samples_for_tic(int(tic))
        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{len(tics)} TICs processed")

    df['nst_samples'] = df['TIC'].apply(lambda t: cache.get(int(t), ''))
    n_with = (df['nst_samples'] != '').sum()
    print(f"  {n_with}/{len(df)} rows have NST samples")
    df.to_csv(csv_path, index=False)
    print(f"  Wrote {csv_path}")

print("\nDone.")
