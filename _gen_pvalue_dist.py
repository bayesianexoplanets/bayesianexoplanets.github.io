"""Regenerate overview_pvalue_dist.png: KDE of -log10(p) for known vs new TOIs,
computed in log-space and normalized to P/P_max. Now includes the 129 known rows
whose NST stats (and hence log10(p)) were just filled in."""
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/global/u2/j/julius/TESS corrected"
known = pd.read_csv(f"{BASE}/tois.csv")
new = pd.read_csv(f"{BASE}/tois_new.csv")

def neglogp(df):
    x = -df["log10(p value)"].to_numpy(dtype=float)
    return x[np.isfinite(x) & (x > 0)]

xk, xn = neglogp(known), neglogp(new)

tg = np.linspace(np.log10(2), np.log10(1e7), 500)
xg = 10 ** tg
def kde_norm(x):
    y = gaussian_kde(np.log10(x))(tg)
    return y / y.max()
yk, yn = kde_norm(xk), kde_norm(xn)

PUR, ORG = "mediumpurple", "#ff7f0e"
fig, ax = plt.subplots(figsize=(9.5, 5.9))
ax.fill_between(xg, yk, color=PUR, alpha=0.28)
ax.fill_between(xg, yn, color=ORG, alpha=0.28)
ax.plot(xg, yk, color=PUR, lw=2.3, label="Known TOIs")
ax.plot(xg, yn, color=ORG, lw=2.3, label="New TOIs")
ax.set_xscale("log")
ax.set_xlim(2, 1e7); ax.set_ylim(0, 1.05)
ax.set_xlabel(r"$-\log_{10}(p)$")
ax.set_ylabel(r"$P/P_{\mathrm{max}}$")
ax.set_title("p-value distribution: known vs. new TOIs")
ax.grid(False)
ax.legend(loc="upper right", framealpha=0.95)
fig.tight_layout()
fig.savefig(f"{BASE}/overview_pvalue_dist.png", dpi=110)
print(f"wrote overview_pvalue_dist.png  | known n={len(xk)}  new n={len(xn)}")
