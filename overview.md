# Overview

We assign each candidate a significance by comparing its matched-filter SNR to a
per-star **null distribution**, built from null-simulation tests (NSTs). The
Gaussian p-value is

$$ p\text{-value} = 1 - \Phi\left(\frac{\mathrm{SNR} - \mu_\mathrm{null}}{\sigma_\mathrm{null}}\right), $$

where $\mu_\mathrm{null}$ and $\sigma_\mathrm{null}$ are the mean and standard
deviation of the maximum SNR recovered in each NST run on that star.

## How many NSTs are needed to estimate $\mu$ and $\sigma$?

Using ~100 NST runs for each of 100 randomly selected stars, we track the running
estimates of $\mu$ and $\sigma$ as a function of the number of NSTs used, $N$.

![Convergence of μ and σ with the number of NSTs](overview_nst_musigma.png)

## Is the Gaussian p-value reliable?

We pick, for each star, the SNR whose **true p-value is $0.01$** (the 99th
percentile of its null max-SNR samples) and estimate that p-value as a function
of $N$ from the Gaussian approximation.

![Gaussian p-value at the p ≈ 0.01 threshold](overview_nst_pvalue.png)

## How significant are the candidates?

Comparing the p-value distributions of the established (known) TOIs and our new
candidates — plotted as $-\log_{10}(p)$ on a log axis, so larger means more
significant:

![p-value distribution of known vs. new TOIs](overview_pvalue_dist.png)
