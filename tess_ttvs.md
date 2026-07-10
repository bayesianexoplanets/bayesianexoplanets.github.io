Transit timing variations (TTVs) for TESS Objects of Interest, measured with this catalog's matched-filter pipeline. Every transit of each candidate is fit individually with the single-transit matched filter (timing errors from the amplitude-marginalized Fisher matrix), giving the O−C residuals against the best-fit linear ephemeris. The O−C measurements are then interpreted with a physical model: a neural network (a normalizing flow) trained on 100,000 simulated two-planet systems maps each measured O−C sequence directly to the posterior distribution of the mass and orbit of the unseen planet perturbing the transiting one. The best two-planet model from this posterior, refined to the local maximum a posteriori solution, is the dashed line in the O−C panel, and the quoted **TTV SNR** is the matched-filter signal-to-noise of the O−C data against it, S = Σ m·y/σ² / √(Σ m²/σ²).

The entries are the best candidates from a visual inspection of the ~250 highest-ranked systems of the full ~5000-star catalog: 25 systems with a clear, coherent O−C signal, plus the systems flagged as significant TTV detections by [Nabbie et al. (2026)](https://arxiv.org/abs/2606.17218) that our data can resolve (ticked in the last column; for those we show the planet with the strongest signal in our data). Known TTV systems recovered independently by this analysis include TOI-1136, TOI-2202 and TOI-282 (HD 28109).

Click a row to see the folded and per-epoch stacked transit (left) next to the O−C measurements with the MAP two-planet model as the dashed line (top right), and the inferred mass and period of the perturbing planet (bottom right). There the filled contours are the neural posterior; where an independent N-body MCMC chain was run for the system, its posterior is overlaid as unfilled contours. Contours enclose 39% and 87% of the probability.

| Column | Description |
|---|---|
| TIC | TESS Input Catalog identifier |
| TOI | TESS Object of Interest number |
| Period [d] | Orbital period in days |
| τ [d] | Transit half-duration in days |
| SNR | Matched-filter signal-to-noise ratio of the periodic transit |
| TTV SNR | Matched-filter signal-to-noise of the O−C data against the MAP physical two-planet model |
| Nₜᵣ | Number of individually measured transits |
| Nabbie+2026 | System flagged as a significant TTV detection in Nabbie et al. (2026) |
