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

The left panel shows $\mu_\mathrm{null}$ and the right panel $\sigma_\mathrm{null}$
(each thin line is one star). The mean $\mu$ stabilises almost immediately, while the standard deviation
$\sigma$ is noisy for $N \lesssim 20$ and settles by $N \approx 30$ — so a few
tens of NSTs suffice to pin down the null distribution.

## Is the Gaussian p-value reliable?

We pick, for each star, the SNR whose **true p-value is $0.01$** (the 99th
percentile of its null max-SNR samples) and estimate that p-value as a function
of $N$ from the Gaussian approximation.

![Gaussian p-value at the p ≈ 0.01 threshold](overview_nst_pvalue.png)

The solid line is the median across the 20 stars and the shaded band is the
25th–75th percentile across them. The Gaussian estimate converges to the true
$p = 0.01$ (red dashed) by $N \approx 40$ and stays within a factor of ~2 of it,
so the per-star null max-SNR distribution is close to Gaussian and the Gaussian
p-value is a reliable significance estimate.
