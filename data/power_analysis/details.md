# Power analysis: details

Generated 2026-10-29 by `power_analysis.py`. Summary: [report.md](report.md). Tables: `box_noise.csv`, `truncation.csv`, `points.csv`, `metrics.csv`, `bootstrap.csv`, `allocation.csv`, `equilibration.csv`, `thresholds.csv`, `literature.csv`.

## Method

1. **Run length.** Every Sage run is re-analysed with the benchmark's estimator (`mixing.observables` and `timeseries.estimate_mean`, as in `summarise-sage`), keeping 0.25, 0.5, 1, 1.5, 2 ns of production after 0.1 or 0.2 ns of equilibration. The analysis at the protocol length (1 ns) reproduces the Sage `runs.csv` exactly (tested).
2. **Single-run noise per box.** σ = k·sem, with sem the RMS within-run standard error over the 5 replicates and k the pooled ratio of replicate SD to within-run sem over all boxes at that length (k = 1.11 to 1.24 for h, 1.00 to 1.21 for density). This is smoother than the 5-replicate SD of each box and corrects the within-run sem for the slow modes it misses.
3. **Points.** ΔHmix and V^E are linear in the box means (ml-peg's `mixing_matrix`). Each pure box feeds all three compositions of its system, so noise is drawn per box and propagated, which keeps the correlation between points. The propagated ΔHmix noise at 1 ns (RMS over points 0.113 kJ/mol) matches the observed SD of Sage's replicate ΔHmix (RMS 0.130 kJ/mol, Sage `points.csv`).
4. **Metrics.** For each run length, 20000 draws of one run per box give the noise of every metric for a model whose true errors are 1, 0.5, 0.25, 0 × Sage's (replicate mean minus experiment). Two models run independently, so their metrics differ by noise with SD sqrt(2)·SD(metric). The **resolvable difference** is 2.80·sqrt(2)·SD(metric): a true difference that size is detected at the two-sided 5 % level 80 % of the time. A closed form sqrt(sᵀCs)/N (s the signs of the errors over σ_ref for MAZE, C the point covariance) is given as a check; it holds when errors are well above the noise.
5. **System choice.** A cluster bootstrap over systems (20000 resamples of the 10 mixtures, or of the 7 scored pure liquids) of Sage's noiseless metric measures how much the metric depends on which systems were picked.

## 1. Noise against run length

![noise scaling](figures/noise_scaling.png)

Single-run noise × sqrt(t / 1 ns), RMS over boxes and relative to 1 ns. A flat line would mean σ ∝ 1/sqrt(t). For h the scaled noise grows by 10% from 0.5 to 2 ns: short runs miss some slow fluctuations, so cutting the run length costs a little less precision than 1/sqrt(t) predicts.

| Production (ns) | σ ΔHmix (kJ/mol) | σ density (g/mL) | σ V^E (cm³/mol) | σ pure density (g/mL) | k (h) | k (density) |
|---|---|---|---|---|---|---|
| 0.25 | 0.187 | 0.00185 | 0.241 | 0.00158 | 1.14 | 1.14 |
| 0.5 | 0.148 | 0.00120 | 0.159 | 0.00114 | 1.24 | 1.00 |
| 1 | 0.099 | 0.00088 | 0.125 | 0.00092 | 1.17 | 1.06 |
| 1.5 | 0.078 | 0.00077 | 0.110 | 0.00078 | 1.11 | 1.11 |
| 2 | 0.068 | 0.00074 | 0.104 | 0.00075 | 1.13 | 1.21 |

Median single-run noise per point. k is the replicate SD / within-run sem ratio (1 if the within-run sem were right).

## 2. A model against experiment, per point

Minimum detectable error MDE = 2.80·sqrt(σ_sim² + σ_exp²): an error that size is significant at the 5 % level 80 % of the time. The fraction column shows how many of Sage's errors (a realistic model) exceed the MDE. Experimental errors are 6 to 15× smaller than the simulation noise (median σ_sim over median σ_exp, ΔHmix and density) and barely change these limits.

| Property | Production (ns) | Median MDE | Max MDE | Sage errors above MDE | Median σ_exp | Unit |
|---|---|---|---|---|---|---|
| ΔHmix | 0.5 | 0.417 | 0.665 | 37% | 0.0067 | kJ/mol |
| ΔHmix | 1 | 0.277 | 0.609 | 60% | 0.0067 | kJ/mol |
| ΔHmix | 2 | 0.192 | 0.523 | 63% | 0.0067 | kJ/mol |
| Mixture density | 0.5 | 0.00339 | 0.00506 | 67% | 0.00015 | g/mL |
| Mixture density | 1 | 0.00269 | 0.00385 | 70% | 0.00015 | g/mL |
| Mixture density | 2 | 0.00222 | 0.00353 | 78% | 0.00015 | g/mL |
| V^E | 0.5 | 0.447 | 0.803 | 22% | 0.0041 | cm³/mol |
| V^E | 1 | 0.35 | 0.57 | 22% | 0.0041 | cm³/mol |
| V^E | 2 | 0.292 | 0.478 | 30% | 0.0041 | cm³/mol |
| Pure density | 0.5 | 0.00384 | 0.00396 | 57% | 0.00015 | g/mL |
| Pure density | 1 | 0.00265 | 0.00393 | 57% | 0.00015 | g/mL |
| Pure density | 2 | 0.00214 | 0.00379 | 57% | 0.00015 | g/mL |

![per-point detection](figures/point_mde.png)

ΔHmix points above the diagonal have a Sage error that one run per box detects.

## 3. Ranking two models

Metric noise for a Sage-like model (SD of the score over repeat benchmark runs), the resolvable difference over all error profiles, and the SD of Sage's score over the choice of systems. Units as in the thresholds table of [report.md](report.md).

| Metric | SD (1 ns) | SD closed form | Resolvable (1 ns) | Resolvable (2 ns) | SD over systems |
|---|---|---|---|---|---|
| ΔHmix MAE | 0.024 | 0.026 | 0.053 to 0.094 | 0.041 to 0.073 | 0.17 |
| ΔHmix MAZE | 0.17 | 0.19 | 0.39 to 0.69 | 0.28 to 0.51 | 1.6 |
| Mixture density MAE | 0.00018 | 0.00018 | 0.00043 to 0.0007 | 0.00035 to 0.00057 | 0.0029 |
| Mixture density MAZE | 0.19 | 0.2 | 0.48 to 0.77 | 0.39 to 0.62 | 2.9 |
| V^E MAE | 0.025 | 0.03 | 0.066 to 0.1 | 0.054 to 0.084 | 0.078 |
| V^E MAZE | 0.19 | 0.22 | 0.48 to 0.74 | 0.39 to 0.61 | 0.64 |
| Pure density MAE | 0.00036 | 0.00037 | 0.00088 to 0.0014 | 0.00075 to 0.0012 | 0.0022 |
| Pure density MAZE | 0.35 | 0.36 | 0.85 to 1.4 | 0.7 to 1.1 | 2.2 |

Noise also biases each score upwards, most for good models: a model whose true ΔHmix MAE is 0.15 kJ/mol scores 0.18 on average at 1 ns, against 0.59 for a true 0.59. Differences between good models are slightly compressed, but their order is kept.

The two kinds of uncertainty answer different questions. The simulation noise says whether rerunning the benchmark would give the same ranking on these 10 systems; at 1 ns it would, for any gap above the resolvable difference. The system spread says how far a model's score depends on this set of systems, which no run length can reduce.

For a ranking, what matters is the system spread of the *difference* between two models, which depends on how alike their errors are. If model B's errors are model A's scaled by (1 + s), the gap is s·MAE_A and its system SD is s·0.17, i.e. 0.29× the gap at any size, so the choice of systems never threatens that ranking. If their errors are unrelated, the gap's system SD approaches sqrt(2)·0.17 = 0.24 kJ/mol, and gaps below about twice that are specific to this set of systems even though a rerun would reproduce them. Real model pairs (e.g. two sizes of one model family) are likely to sit between these cases.

## 4. Run length and cost

| Protocol (ns) | Cost | MACE-OFF23(M) GPU-days | Resolvable ΔHmix MAE (kJ/mol) | Resolvable density MAE (g/mL) | Perfect ΔHmix MAZE |
|---|---|---|---|---|---|
| 0.2 + 0.25 | 38% | 52 | 0.141 | 0.0014 | 1.07 |
| 0.2 + 0.5 | 58% | 81 | 0.122 | 0.0009 | 0.89 |
| 0.2 + 1 | 100% | 138 | 0.094 | 0.0007 | 0.62 |
| 0.2 + 1.5 | 142% | 196 | 0.078 | 0.0006 | 0.49 |
| 0.2 + 2 | 183% | 254 | 0.073 | 0.0006 | 0.45 |

Cost is Σ atoms × simulated time over the 50 boxes (76,065 atoms). GPU-days are an order of magnitude only: MACE-OFF23(M) runs 1.1×10⁶ steps/day for ~600-atom liquid boxes on one A100 ([arXiv:2312.15211](https://arxiv.org/abs/2312.15211), sec. III.3.2), scaled linearly to our 952 to 2,400-atom boxes. Faster and slower models scale the column, not the ratios.

**Equilibration.** The mean over the 100 to 200 ps window minus the mean over the 1 ns production, over all 250 Sage runs, is -0.044 kJ/mol for h (p = 0.004) and +0.0007 g/mL for density (p = 5e-04). Even for Sage, starting from Sage-built boxes, the window is shifted. By replicate, the density shift is -0.0009 (0), +0.0021 (1), +0.0019 (2), +0.0009 (3), -0.0007 (4) g/mL; replicate 0 is the distributed box that MLP runs start from. An MLP also has to relax from Sage's density to its own, which this test cannot see, so 0.2 ns stays. Table: [equilibration.csv](equilibration.csv).

**Per-box run lengths.** For a fixed variance of ΔHmix MAE, mixture density MAE and pure density MAE (equally weighted, each relative to uniform 1 ns), the cheapest allocation runs each box for t ∝ sqrt(A_b / n_b), where A_b is the box's variance contribution and n_b its atom count, with at least 0.5 ns. It ranges from 0.5 to 2.0 ns ([allocation.csv](allocation.csv)); the longest runs go to pure_n_methylaniline_298, pure_tetrahydropyran_298, pure_pyrrole_298, acetone_nmethylaniline_x230. It saves 14% of the cost. That assumes σ ∝ 1/sqrt(t), which section 1 shows is optimistic for short runs, so the real saving is smaller. Not worth the extra complexity.

## 5. Literature: how far apart are current MLPs?

| Model | Property | Result | Source |
|---|---|---|---|
| MACE-OFF23(M) | pure density, 109 organic liquids | MAE 0.09 g/mL | [Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2 |
| MACE-OFF23(S) | pure density, 109 organic liquids | MAE 0.23 g/mL (systematically too dense) | [Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2 |
| MACE-OFF23(L) | pure density, 109 organic liquids | no significant improvement over (M) | [Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2 |
| ANI-2x | pure density, 109 organic liquids | MAE and RMSE about twice MACE-OFF23(M) | [Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2 |
| MACE-OFF23(M) | ΔHvap, 109 organic liquids | systematic offset of about 2 kcal/mol (8 kJ/mol) | [Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2 |
| Organic_MPNICE | pure density, 62 solvents | about 4 % mean error (4.9 % in Table 12) | [Weber et al. 2025](https://arxiv.org/abs/2505.06462), sec. 3.1.4 |
| QRNN | pure density, 62 solvents | about 16 % mean error, systematically high | [Weber et al. 2025](https://arxiv.org/abs/2505.06462), sec. 3.1.4 |
| Hybrid / inorganic MPNICE | pure density, 62 solvents | 18 to 35 % mean absolute error | [Weber et al. 2025](https://arxiv.org/abs/2505.06462), Table 12 |
| ANI-2x, MACE-OFF23(S, M) vs CGenFF | density and ΔHvap of 6 liquids (water, methanol, acetone, benzene, hexane, NMA) | the MLPs miss in both directions (e.g. MACE-OFF23(M) too dense for water, too light for acetone and hexane); CGenFF best overall | [Picha et al. 2025](https://doi.org/10.1021/acs.jcim.5c01355) |
| BAMBOO | density of electrolyte solvent mixtures | about 0.01 g/mL average error, after a density alignment step fitted to experiment | [Gong et al. 2025](https://doi.org/10.1038/s42256-025-01009-7) ([arXiv:2404.07181](https://arxiv.org/abs/2404.07181)) |
| UMA-OMol vs SevenNet-OMat, Orb-OMat | density of ether and carbonate electrolytes | UMA-OMol much closer to experiment than models trained on inorganic data (qualitative in the abstract) | [arXiv:2603.20183](https://arxiv.org/abs/2603.20183) |

No study computing ΔHmix or V^E of molecular liquids with an MLP was found (searches of the web and of Consensus for MLP / neural network potential enthalpy of mixing, excess enthalpy and excess volume, October 2026). The closest work concerns metallic alloys ([arXiv:2406.17499](https://arxiv.org/abs/2406.17499)) or fits mixing data without MD.

- **Density.** Published MLPs differ in pure-liquid density MAE by about 0.05 to 0.25 g/mL. The benchmark resolves 0.0007 g/mL at 1 ns and 0.0014 g/mL at 0.25 ns. Density alone could be ranked from much shorter runs, but the runs serve ΔHmix too.
- **ΔHmix.** There is no MLP baseline. MLP ΔHvap errors are several kJ/mol (MACE-OFF23(M) offset about 8 kJ/mol; Picha et al. find large ΔHvap deviations for MACE-OFF23(M)). ΔHmix is a difference of such cohesive energies, so much of that error may cancel, but whatever remains is compared with |ΔHmix| near x = 0.5 of 0.6 to 2.6 kJ/mol. If MLPs differ in ΔHmix MAE by even a fifth of Sage's MAE (about 0.12 kJ/mol), 1 ns resolves it.

## Caveats

- **Sage stands in for MLP noise.** MLP fluctuations of h and volume should be similar (the same liquids at the same T and P), but correlation times may differ, especially for MLPs that change conformer barriers. One pilot MLP run per box, analysed with `mixing.box_result`, would check this: compare its sem and n_eff with the Sage `runs.csv`.
- **Slow modes.** k = 1.17 (h, 1 ns) is pooled over boxes. Single boxes have larger run-to-run spreads (for Sage, NMP + 1,2-dichloroethane; see the [Sage report](../sage_summary/report.md)), and an MLP may have slow modes of its own.
- **Systematic simulation errors** (finite size, starting structures) are not noise and are not included. Against the published 1000-molecule Sage runs, the median ΔHmix gap is 0.043 kJ/mol and the largest 0.22 kJ/mol (propanol_benzene_x246).
- **Error profiles.** The model pairs are built from Sage's errors. Real MLPs may spread their errors differently; the closed-form SD shows the result depends mainly on the noise, not on the profile.
