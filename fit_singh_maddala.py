"""fit_singh_maddala.py — per-star Bayesian Singh-Maddala (Burr XII) fit to the NST max-SNR samples.

The null (NST) max-SNR distribution is Singh-Maddala  f(x|c,k,lam) = (ck/lam)(x/lam)^(c-1)/(1+(x/lam)^c)^(k+1)
( = scipy.stats.burr12(c=c, d=k, scale=lam) ). We fit it per star to that star's <=10 NST samples by
NUTS, vectorized over stars, with whole-sample-informed log-normal priors (NST Analysis.ipynb):

    c   = 27.904 * LogNormal(0, 0.469)
    k   = 0.841  * LogNormal(0, 0.826)
    lam = 4.864  + 1.816 * LogNormal(0, 0.906)

The p-value at a candidate SNR x is the posterior-mean survival function
    p = E_post[ SF_burr12(x; c,k,lam) ] ,
computed stably (deep tail) with burr12.logsf + logsumexp:
    log10(p) = ( logsumexp_i burr12.logsf(x; c_i,k_i,lam_i) - log N ) / ln10 .
"""
import os
import numpy as np

os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS
from scipy.stats import burr12
from scipy.special import logsumexp

# prior hyperparameters (from NST Analysis.ipynb)
C0, CS = 27.904, 0.469
K0, KS = 0.841, 0.826
L_LOC, L0, LS = 4.864, 1.816, 0.906
WHOLE = (8.309, 0.549, 7.386)   # whole-sample (c,k,lam); fallback when a star has no usable samples
LN10 = np.log(10.0)


def _model(data, mask):
    """Vectorized over stars: data,mask shape (S, M); padded entries (mask=0) contribute nothing."""
    S = data.shape[0]
    with numpyro.plate("stars", S):
        c_raw = numpyro.sample("c_raw", dist.LogNormal(0.0, CS))
        k_raw = numpyro.sample("k_raw", dist.LogNormal(0.0, KS))
        l_raw = numpyro.sample("l_raw", dist.LogNormal(0.0, LS))
    c = numpyro.deterministic("c", C0 * c_raw)
    k = numpyro.deterministic("k", K0 * k_raw)
    lam = numpyro.deterministic("lam", L_LOC + L0 * l_raw)
    data_safe = jnp.where(mask, data, 1.0)              # avoid log(0) on padding
    logz = jnp.log(data_safe / lam[:, None])
    log1pzc = jnp.logaddexp(0.0, c[:, None] * logz)     # log(1 + (x/lam)^c), overflow-safe
    logp = ((jnp.log(c) + jnp.log(k) - jnp.log(lam))[:, None]
            + (c[:, None] - 1.0) * logz - (k[:, None] + 1.0) * log1pzc)
    numpyro.factor("lik", jnp.sum(jnp.where(mask, logp, 0.0)))


def fit_batch(samples_list, num_warmup=1000, num_samples=1000, num_chains=4,
              seed=0, target_accept=0.95):
    """samples_list: list of 1D arrays (each star's NST max-SNRs).
    Returns (c, k, lam) posterior draws, each shape (num_chains*num_samples, S)."""
    S = len(samples_list)
    M = max(1, max(len(s) for s in samples_list))
    data = np.ones((S, M)); mask = np.zeros((S, M), dtype=bool)
    for i, s in enumerate(samples_list):
        s = np.asarray(s, dtype=float)
        n = len(s)
        if n:
            data[i, :n] = s; mask[i, :n] = True
    mcmc = MCMC(NUTS(_model, target_accept_prob=target_accept),
                num_warmup=num_warmup, num_samples=num_samples, num_chains=num_chains,
                chain_method="vectorized", progress_bar=True)
    mcmc.run(jax.random.PRNGKey(seed), jnp.asarray(data), jnp.asarray(mask))
    post = mcmc.get_samples(group_by_chain=False)
    return np.asarray(post["c"]), np.asarray(post["k"]), np.asarray(post["lam"])


def log10p(c_draws, k_draws, lam_draws, x):
    """Posterior-mean SF at SNR x for one star, in log10. c/k/lam_draws: 1D arrays (N,).

    Uses the closed-form burr12.logsf = -k*log(1+(x/lam)^c), evaluated as
    -k*logaddexp(0, c*log(x/lam)) so it never overflows for deep-tail x (large SNR)."""
    logsf = -k_draws * np.logaddexp(0.0, c_draws * np.log(x / lam_draws))   # (N,)  == burr12.logsf
    return (logsumexp(logsf) - np.log(len(c_draws))) / LN10


def parse_nst(s):
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return np.array([])
    s = str(s).strip()
    if not s:
        return np.array([])
    return np.array([float(v) for v in s.split("|")], dtype=float)
