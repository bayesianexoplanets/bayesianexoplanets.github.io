Transit timing variations (TTVs) for TESS Objects of Interest, measured with this catalog's matched-filter pipeline. Every transit of each candidate is fit individually with the single-transit matched filter (timing errors from the amplitude-marginalized Fisher matrix), and the O−C residuals against a weighted linear ephemeris are tested with a likelihood ratio between three nested noise models: measurement noise only, an added white timing jitter, and an added time-correlated (Matérn-5/2) TTV signal. The table lists the two resulting statistics: **TTV SNR** = √(2 Δlog L) for the coherent signal over white jitter, and **scatter SNR** = √(2 Δlog L) for excess white scatter over the measurement errors (which also absorbs error-bar miscalibration, so it should be read as a ranking statistic rather than a calibrated significance).

The entries are the best candidates from a visual inspection of the ~250 highest-ranked systems of the full ~5000-star catalog (ranked by the per-transit-normalized statistics): 25 systems with a clear, coherent O−C signal, plus the systems flagged as significant TTV detections by [Nabbie et al. (2026)](https://arxiv.org/abs/2606.17218) that our data can resolve (ticked in the last column; for those we show the planet with the strongest signal in our data). Known TTV systems recovered independently by this analysis include TOI-1136, TOI-2202 and TOI-282 (HD 28109). Amplitudes range from a few minutes (e.g. TOI-6101, TOI-1924) to days (TOI-4504).

Click a row to see the folded and per-epoch stacked transit (left) next to the O−C curve with its maximum-likelihood Gaussian-process overlay (right).

| Column | Description |
|---|---|
| TIC | TESS Input Catalog identifier |
| TOI | TESS Object of Interest number |
| Period [d] | Orbital period in days |
| τ [d] | Transit half-duration in days |
| SNR | Matched-filter signal-to-noise ratio of the periodic transit |
| TTV SNR | Coherent-TTV likelihood-ratio statistic |
| scatter SNR | Excess-white-scatter likelihood-ratio statistic |
| Nₜᵣ | Number of individually measured transits |
| Nabbie+2026 | System flagged as a significant TTV detection in Nabbie et al. (2026) |
