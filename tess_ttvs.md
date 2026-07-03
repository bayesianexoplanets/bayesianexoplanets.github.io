Transit timing variations (TTVs) for TESS Objects of Interest, measured with this catalog's matched-filter pipeline. Every transit of each candidate is fit individually with the single-transit matched filter (timing errors from the amplitude-marginalized Fisher matrix), giving the O−C residuals against the best-fit linear ephemeris. The TTV signal is then modeled as a sinusoid with a smoothly varying amplitude and phase, m(t) = f₁(t)·sin(ωt) + f₂(t)·cos(ωt), where the envelopes f₁, f₂ are penalized cubic splines with knots every ~3 orbital periods (penalty strength set by cross-validation) and the frequency ω is optimized starting from integer multiples of the orbital frequency. The quoted **TTV SNR** is the matched-filter signal-to-noise of the O−C data against this best-fit model, S = Σ m·y/σ² / √(Σ m²/σ²).

The entries are the best candidates from a visual inspection of the ~250 highest-ranked systems of the full ~5000-star catalog: 25 systems with a clear, coherent O−C signal, plus the systems flagged as significant TTV detections by [Nabbie et al. (2026)](https://arxiv.org/abs/2606.17218) that our data can resolve (ticked in the last column; for those we show the planet with the strongest signal in our data). Known TTV systems recovered independently by this analysis include TOI-1136, TOI-2202 and TOI-282 (HD 28109).

Click a row to see the folded and per-epoch stacked transit (left) next to the O−C measurements with the best-fit TTV model as the dashed line (right).

| Column | Description |
|---|---|
| TIC | TESS Input Catalog identifier |
| TOI | TESS Object of Interest number |
| Period [d] | Orbital period in days |
| τ [d] | Transit half-duration in days |
| SNR | Matched-filter signal-to-noise ratio of the periodic transit |
| TTV SNR | Matched-filter signal-to-noise of the O−C data against the best-fit TTV model |
| Nₜᵣ | Number of individually measured transits |
| Nabbie+2026 | System flagged as a significant TTV detection in Nabbie et al. (2026) |
