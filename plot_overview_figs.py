"""plot_overview_figs.py — regenerate the three Overview figures for the Singh-Maddala p-value.

  overview_nst_musigma.png : empirical NST PDF + SF (~5000 stars) vs the Singh-Maddala fit AND the
                             Gaussian fit -> shows the Gaussian misses the heavy tail.
  overview_nst_pvalue.png  : p-value estimator vs #NST samples -- Gaussian, no-prior Burr (MLE),
                             prior Burr (Bayes) -- vs the true p; shows the informed prior + the
                             Singh-Maddala model are needed.
  overview_pvalue_dist.png : -log10(p) distribution of known TOIs vs new candidates (new p-values).
Run after the catalogs carry the new 'log10(p value)' (Phase 3) for the third figure.
"""
import os, sys
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
import scipy.stats as stats
import fit_singh_maddala as sm

rc("mathtext", fontset="cm")
plt.rcParams["font.family"] = "serif"; plt.rcParams["font.serif"] = ["DejaVu Serif"]; plt.rcParams["font.size"] = 11
HERE = "/global/u2/j/julius/TESS corrected"
C, K, L = 8.309, 0.549, 7.386          # whole-sample Singh-Maddala params (NST Analysis.ipynb)
SM = stats.burr12(c=C, d=K, scale=L)


def load_all_nst():
    d = pd.read_csv(f"{HERE}/tois.csv")
    out = []
    for s in d["nst_samples"].dropna():
        v = sm.parse_nst(s)
        if len(v): out.append(v)
    return np.concatenate(out)


def fig_fit():
    nst = load_all_nst()
    mu, sig = nst.mean(), nst.std()
    x = np.linspace(0, 35, 300)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    # PDF
    ax[0].hist(nst, density=True, bins=100, range=(0, 35), histtype="step", color="black", lw=1.3,
               label=f"NST over {len(load_all_nst()):,}".replace(",", " ") + " samples")
    ax[0].plot(x, SM.pdf(x), lw=2, color="navy", label="Singh-Maddala")
    ax[0].plot(x, stats.norm(mu, sig).pdf(x), lw=2, color="crimson", ls="--", label="Gaussian")
    ax[0].set_yscale("log"); ax[0].set_xlim(0, 35); ax[0].set_ylim(1e-4, None)
    ax[0].set_xlabel("NST SNR"); ax[0].set_ylabel("pdf"); ax[0].legend(fontsize=9)
    # SF (the heavy-tail point)
    ax[1].hist(nst, density=True, bins=100, range=(0, 35), histtype="step", cumulative=-1,
               color="black", lw=1.3, label="NST (empirical SF)")
    ax[1].plot(x, SM.sf(x), lw=2, color="navy", label="Singh-Maddala SF")
    ax[1].plot(x, stats.norm(mu, sig).sf(x), lw=2, color="crimson", ls="--", label="Gaussian SF")
    ax[1].set_yscale("log"); ax[1].set_xlim(0, 35); ax[1].set_ylim(1e-4, 1)
    ax[1].set_xlabel("NST SNR"); ax[1].set_ylabel("sf  (= 1 - cdf)"); ax[1].legend(fontsize=9)
    fig.tight_layout(); fig.savefig(f"{HERE}/overview_nst_musigma.png", dpi=130, bbox_inches="tight")
    plt.close(fig); print("wrote overview_nst_musigma.png")


def fig_convergence(seeds=12, p_true=1e-2):
    """p-value estimator vs N samples: Gaussian, no-prior Burr (MLE), prior Burr (Bayes) vs true p."""
    ns = np.unique(np.geomspace(5, 1000, 22, dtype=int))
    rng = np.random.default_rng(0)
    # synthetic stars: per (seed, N) draw true params from the priors, generate N samples
    truth, datasets, testpts = [], [], []
    for s in range(seeds):
        c_t = 27.904 * np.exp(rng.normal(0, sm.CS)); k_t = 0.841 * np.exp(rng.normal(0, sm.KS))
        l_t = 4.864 + 1.816 * np.exp(rng.normal(0, sm.LS))
        d_full = stats.burr12(c=c_t, d=k_t, scale=l_t).rvs(size=ns.max(), random_state=rng)
        tp = stats.burr12(c=c_t, d=k_t, scale=l_t).isf(p_true)
        for n in ns:
            truth.append((c_t, k_t, l_t)); datasets.append(d_full[:n]); testpts.append(tp)
    # prior-Burr (Bayes) for ALL (seed,N) at once on the GPU
    c, k, lam = sm.fit_batch(datasets, num_warmup=800, num_samples=800, num_chains=4, seed=0, target_accept=0.95)
    P_prior = np.full((seeds, len(ns)), np.nan); P_mle = np.full_like(P_prior, np.nan); P_gauss = np.full_like(P_prior, np.nan)
    idx = 0
    for s in range(seeds):
        for j, n in enumerate(ns):
            tp = testpts[idx]; data = datasets[idx]
            P_prior[s, j] = 10 ** sm.log10p(c[:, idx], k[:, idx], lam[:, idx], tp)
            try:
                th = stats.burr12.fit(data, floc=0); P_mle[s, j] = stats.burr12(c=th[0], d=th[1], scale=th[3]).sf(tp)
            except Exception:
                pass
            P_gauss[s, j] = stats.norm(np.mean(data), np.std(data)).sf(tp)
            idx += 1

    def band(P):
        P = np.where(P > 0, P, np.nan)
        return np.nanpercentile(P, 16, 0), np.nanpercentile(P, 50, 0), np.nanpercentile(P, 84, 0)
    fig, axx = plt.subplots(figsize=(7.5, 5.2))
    for P, col, lab in [(P_gauss, "crimson", "Gaussian"), (P_mle, "forestgreen", "Singh-Maddala (no prior, MLE)"),
                        (P_prior, "navy", "Singh-Maddala (informed prior, Bayes)")]:
        lo, med, hi = band(P); axx.fill_between(ns, lo, hi, color=col, alpha=0.2); axx.plot(ns, med, color=col, lw=2, label=lab)
    axx.hlines(p_true, ns.min(), ns.max(), color="black", ls="dashed", lw=2, label=r"true $p = 10^{-2}$")
    axx.loglog(); axx.set_xlabel("number of NST samples"); axx.set_ylabel(r"$p$-value estimate"); axx.legend(fontsize=9, loc="lower right")
    axx.set_xlim(ns.min(), ns.max()); axx.set_ylim(1e-10, 3)   # clip the MLE's small-N excursions
    fig.tight_layout(); fig.savefig(f"{HERE}/overview_nst_pvalue.png", dpi=130, bbox_inches="tight")
    plt.close(fig); print("wrote overview_nst_pvalue.png")


def fig_pdist():
    from scipy.stats import gaussian_kde
    known = pd.read_csv(f"{HERE}/tois.csv")
    new = pd.read_csv(f"{HERE}/tois_new.csv")
    kk = -known.loc[known["passed_all_tests"] == True, "log10(p value)"].to_numpy(float)
    nn = -new["log10(p value)"].to_numpy(float)
    kk = kk[np.isfinite(kk)]; nn = nn[np.isfinite(nn)]
    xmax = max(kk.max(), nn.max())
    xg = np.linspace(0.01, xmax, 600)
    fig, axx = plt.subplots(figsize=(7.5, 5.2))
    for data, col, lab in [(kk, "navy", "known TOIs"), (nn, "darkorange", "new candidates")]:
        y = gaussian_kde(data)(xg); y = y / y.max()          # P / P_max
        axx.plot(xg, y, color=col, lw=2, label=lab)
        axx.fill_between(xg, y, color=col, alpha=0.12)
    axx.axvline(2, color="red", ls="dashed", lw=1.5, label=r"$p = 0.01$")
    axx.set_xlim(0.01, xmax); axx.set_ylim(0, None)
    axx.set_xlabel(r"$-\log p$"); axx.set_ylabel(r"$P/P_{\max}$"); axx.legend(fontsize=10)
    fig.tight_layout(); fig.savefig(f"{HERE}/overview_pvalue_dist.png", dpi=130, bbox_inches="tight")
    plt.close(fig); print("wrote overview_pvalue_dist.png")


if __name__ == "__main__":
    which = sys.argv[1:] or ["fit", "conv", "dist"]
    if "fit" in which: fig_fit()
    if "conv" in which: fig_convergence()
    if "dist" in which: fig_pdist()
