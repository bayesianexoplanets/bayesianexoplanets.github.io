"""Final verification: cov indices (0=P,1=phase,2=tau); per-run-max NST with SNR=0 mask."""
import os, glob
import numpy as np
import pandas as pd

FULL = "/pscratch/sd/j/julius/exoprob/FullRun"
tois = pd.read_csv("/global/u2/j/julius/TESS corrected/tois.csv")

print("--- cov indices (period=0, phase=1, tau=2) ---")
for tic in [149603524, 69679391, 229791084]:
    npzs = sorted(glob.glob(os.path.join(FULL, "known_planets", f"{tic}_*.npz")))
    z = np.load(npzs[0], allow_pickle=True)
    cov = np.asarray(z["cov"]); p = z["params"]
    if cov.shape == (6,):
        ep, eph, et = np.sqrt(np.abs(cov[0])), np.sqrt(np.abs(cov[1])), np.sqrt(np.abs(cov[2]))
        print(f"TIC {tic}: P={p[0]:.6f}+/-{ep:.3e}  phase={p[1]:.6f}+/-{eph:.3e}  tau={p[2]:.6f}+/-{et:.3e}")
    else:
        print(f"TIC {tic}: cov shape {cov.shape} (TTV)")

def per_run_maxes(tic, mask0):
    m = []
    for b in range(1, 11):
        f = os.path.join(FULL, "candidates", f"batch{b}", f"{tic}.0.csv")
        if os.path.exists(f):
            df = pd.read_csv(f, sep="\t")
            if "SNR" in df.columns and len(df):
                s = df["SNR"]
                if mask0:
                    s = s[s != 0]
                if len(s):
                    m.append(float(s.max()))
    return np.array(m)

print("\n--- per-run-max ddof=0, batch1..10: CSV vs no-mask vs mask(SNR!=0) ---")
print(f"{'TIC':>11} {'csv_mu':>8} {'csv_sig':>8} | {'nm_mu':>8} {'nm_sig':>8} | {'m_mu':>8} {'m_sig':>8}")
for tic in [59582240, 355800238, 149349867, 402319411, 443341002, 230292034, 141608198, 204496071, 129979528]:
    r = tois[tois.TIC == tic].iloc[0]
    a = per_run_maxes(tic, False); b = per_run_maxes(tic, True)
    print(f"{tic:>11} {r['μ(SNR | null)']:8.4f} {r['σ(SNR | null)']:8.4f} | "
          f"{a.mean():8.4f} {a.std(ddof=0):8.4f} | {b.mean():8.4f} {b.std(ddof=0):8.4f}")

# how often does SNR==0 appear in NST runs?
print("\n--- frequency of SNR==0 rows in NST runs (sample of 40 stars) ---")
nz = ntot = 0
for tic in tois.TIC.dropna().astype(int).head(40):
    for bb in range(1, 11):
        f = os.path.join(FULL, "candidates", f"batch{bb}", f"{tic}.0.csv")
        if os.path.exists(f):
            df = pd.read_csv(f, sep="\t")
            ntot += len(df); nz += int((df["SNR"] == 0).sum())
print(f"SNR==0 rows: {nz} of {ntot} total NST candidate rows")
