# Power analysis: details

Generated 2026-10-29 by `power_analysis.py`. Summary: [report.md](report.md). Tables: `box_noise.csv`, `truncation.csv`, `points.csv`, `metrics.csv`, `bootstrap.csv`, `allocation.csv`, `equilibration.csv`, `thresholds.csv`, `literature.csv`.

## Method

1. **Run length.** Every Sage run is re-analysed with the benchmark's estimator (`mixing.observables` and `timeseries.estimate_mean`, as in `summarise-sage`), keeping 0.25, 0.5, 1, 1.5, 2 ns of production after 0.1 or 0.2 ns of equilibration. The 2 ns analysis reproduces the Sage `runs.csv` exactly (tested).
2. **Single-run noise per box.** σ = k·sem, with sem the RMS within-run standard error over the 5 replicates and k the pooled ratio of replicate SD to within-run sem over all boxes at that length (k = 1.15 to 1.21 for h, 1.08 to 1.27 for density). This is smoother than the 5-replicate SD of each box and corrects the within-run sem for the slow modes it misses.
3. **Points.** ΔHmix and V^E are linear in the box means (ml-peg's `mixing_matrix`). Each pure box feeds all three compositions of its system, so noise is drawn per box and propagated, which keeps the correlation between points. The propagated ΔHmix noise at 2 ns (RMS over points 0.098 kJ/mol) matches the observed SD of Sage's replicate ΔHmix (RMS 0.101 kJ/mol, Sage `points.csv`).
4. **Metrics.** For each run length, 20000 draws of one run per box give the noise of every metric for a model whose true errors are 1, 0.5, 0.25, 0 × Sage's (replicate mean minus experiment). Two models run independently, so their metrics differ by noise with SD sqrt(2)·SD(metric). The **resolvable difference** is 2.80·sqrt(2)·SD(metric): a true difference that size is detected at the two-sided 5 % level 80 % of the time. A closed form sqrt(sᵀCs)/N (s the signs of the errors over σ_ref for MAZE, C the point covariance) is given as a check; it holds when errors are well above the noise.
5. **System choice.** A cluster bootstrap over systems (20000 resamples of the 13 mixtures, or of the 11 scored pure liquids) of Sage's noiseless metric measures how much the metric depends on which systems were picked.

## 1. Noise against run length

![noise scaling](figures/noise_scaling.png)

Single-run noise × sqrt(t / 2 ns), RMS over boxes and relative to 2 ns. A flat line would mean σ ∝ 1/sqrt(t). For h the scaled noise grows by 18% from 0.5 to 2 ns: short runs miss some slow fluctuations, so cutting the run length costs a little less precision than 1/sqrt(t) predicts.

| Production (ns) | σ ΔHmix (kJ/mol) | σ density (g/mL) | σ V^E (cm³/mol) | σ pure density (g/mL) | k (h) | k (density) |
|---|---|---|---|---|---|---|
| 0.25 | 0.200 | 0.00195 | 0.271 | 0.00170 | 1.18 | 1.23 |
| 0.5 | 0.150 | 0.00129 | 0.183 | 0.00119 | 1.21 | 1.08 |
| 1 | 0.109 | 0.00091 | 0.138 | 0.00094 | 1.19 | 1.11 |
| 1.5 | 0.088 | 0.00081 | 0.118 | 0.00076 | 1.15 | 1.17 |
| 2 | 0.080 | 0.00077 | 0.113 | 0.00071 | 1.19 | 1.27 |

Median single-run noise per point. k is the replicate SD / within-run sem ratio (1 if the within-run sem were right).

## 2. A model against experiment, per point

Minimum detectable error MDE = 2.80·sqrt(σ_sim² + σ_exp²): an error that size is significant at the 5 % level 80 % of the time. The fraction column shows how many of Sage's errors (a realistic model) exceed the MDE. Experimental errors are 5 to 12× smaller than the simulation noise (median σ_sim over median σ_exp, ΔHmix and density) and barely change these limits.

| Property | Production (ns) | Median MDE | Max MDE | Sage errors above MDE | Median σ_exp | Unit |
|---|---|---|---|---|---|---|
| ΔHmix | 0.5 | 0.42 | 0.795 | 38% | 0.0069 | kJ/mol |
| ΔHmix | 1 | 0.306 | 0.67 | 49% | 0.0069 | kJ/mol |
| ΔHmix | 2 | 0.225 | 0.551 | 67% | 0.0069 | kJ/mol |
| Mixture density | 0.5 | 0.00363 | 0.00578 | 69% | 0.00015 | g/mL |
| Mixture density | 1 | 0.00275 | 0.00432 | 75% | 0.00015 | g/mL |
| Mixture density | 2 | 0.00232 | 0.00358 | 83% | 0.00015 | g/mL |
| V^E | 0.5 | 0.515 | 0.866 | 17% | 0.0044 | cm³/mol |
| V^E | 1 | 0.388 | 0.639 | 17% | 0.0044 | cm³/mol |
| V^E | 2 | 0.318 | 0.537 | 22% | 0.0044 | cm³/mol |
| Pure density | 0.5 | 0.00345 | 0.00424 | 64% | 0.00015 | g/mL |
| Pure density | 1 | 0.00277 | 0.00409 | 64% | 0.00015 | g/mL |
| Pure density | 2 | 0.00233 | 0.00398 | 73% | 0.00015 | g/mL |

![per-point detection](figures/point_mde.png)

ΔHmix points above the diagonal have a Sage error that one run per box detects.

## 3. Ranking two models

Metric noise for a Sage-like model (SD of the score over repeat benchmark runs), the resolvable difference over all error profiles, and the SD of Sage's score over the choice of systems. Units as in the thresholds table of [report.md](report.md).

| Metric | SD (2 ns) | SD closed form | Resolvable (2 ns) | Resolvable (1 ns) | SD over systems |
|---|---|---|---|---|---|
| ΔHmix MAE | 0.019 | 0.02 | 0.041 to 0.075 | 0.054 to 0.096 | 0.13 |
| ΔHmix MAZE | 0.17 | 0.18 | 0.37 to 0.67 | 0.47 to 0.85 | 1.5 |
| Mixture density MAE | 0.00013 | 0.00013 | 0.00032 to 0.00053 | 0.0004 to 0.00066 | 0.0023 |
| Mixture density MAZE | 0.16 | 0.16 | 0.38 to 0.62 | 0.47 to 0.76 | 3.2 |
| V^E MAE | 0.02 | 0.024 | 0.053 to 0.08 | 0.064 to 0.095 | 0.06 |
| V^E MAZE | 0.21 | 0.25 | 0.55 to 0.81 | 0.67 to 0.98 | 0.72 |
| Pure density MAE | 0.00025 | 0.00026 | 0.00061 to 0.001 | 0.00071 to 0.0012 | 0.0019 |
| Pure density MAZE | 0.27 | 0.28 | 0.66 to 1.1 | 0.78 to 1.3 | 1.8 |

Noise also biases each score upwards, most for good models: a model whose true ΔHmix MAE is 0.14 kJ/mol scores 0.16 on average at 2 ns, against 0.55 for a true 0.55. Differences between good models are slightly compressed, but their order is kept.

The two kinds of uncertainty answer different questions. The simulation noise says whether rerunning the benchmark would give the same ranking on these 13 systems; at 2 ns it would, for any gap above the resolvable difference. The system spread says how far a model's score depends on this set of systems, which no run length can reduce.

For a ranking, what matters is the system spread of the *difference* between two models, which depends on how alike their errors are. If model B's errors are model A's scaled by (1 + s), the gap is s·MAE_A and its system SD is s·0.13, i.e. 0.24× the gap at any size, so the choice of systems never threatens that ranking. If their errors are unrelated, the gap's system SD approaches sqrt(2)·0.13 = 0.19 kJ/mol, and gaps below about twice that are specific to this set of systems even though a rerun would reproduce them. Real model pairs (e.g. two sizes of one model family) are likely to sit between these cases.

## 4. Shorter runs

| Protocol (ns) | Cost | MACE-OFF23(M) GPU-days | Resolvable ΔHmix MAE (kJ/mol) | Resolvable density MAE (g/mL) | Perfect ΔHmix MAZE |
|---|---|---|---|---|---|
| 0.2 + 0.25 | 20% | 71 | 0.139 | 0.0013 | 1.49 |
| 0.2 + 0.5 | 32% | 110 | 0.116 | 0.0009 | 1.17 |
| 0.2 + 1 | 55% | 189 | 0.096 | 0.0007 | 0.87 |
| 0.2 + 1.5 | 77% | 268 | 0.079 | 0.0006 | 0.70 |
| 0.2 + 2 | 100% | 347 | 0.075 | 0.0005 | 0.65 |

Cost is Σ atoms × simulated time over the 64 boxes (104,076 atoms). GPU-days are an order of magnitude only: MACE-OFF23(M) runs 1.1×10⁶ steps/day for ~600-atom liquid boxes on one A100 ([arXiv:2312.15211](https://arxiv.org/abs/2312.15211), sec. III.3.2), scaled linearly to our 952 to 2,800-atom boxes. Faster and slower models scale the column, not the ratios.

**Equilibration.** The mean over the 100 to 200 ps window minus the mean over the 2 ns production, over all 320 Sage runs, is -0.045 kJ/mol for h (p = 0.008) and +0.0008 g/mL for density (p = 2e-06). Even for Sage, starting from Sage-built boxes, the window is shifted. By replicate, the density shift is -0.0007 (0), +0.0017 (1), +0.0019 (2), +0.0013 (3), -0.0003 (4) g/mL; replicate 0 is the distributed box that MLP runs start from. An MLP also has to relax from Sage's density to its own, which this test cannot see, so 0.2 ns stays. Table: [equilibration.csv](equilibration.csv).

**Per-box run lengths.** For a fixed variance of ΔHmix MAE, mixture density MAE and pure density MAE (equally weighted, each relative to uniform 2 ns), the cheapest allocation runs each box for t ∝ sqrt(A_b / n_b), where A_b is the box's variance contribution and n_b its atom count, with at least 0.5 ns. It ranges from 0.5 to 4.7 ns ([allocation.csv](allocation.csv)); the longest runs go to pure_n_methylaniline_298, pure_2_dimethylamino_ethanol_298, pure_tetrahydropyran_298, pure_pyrrole_298. It saves 17% of the cost. That assumes σ ∝ 1/sqrt(t), which section 1 shows is optimistic for short runs, so the real saving is smaller. Not worth the extra complexity.

## 5. Literature: how far apart are current MLPs?

| Model | Property | Result | Source |
|---|---|---|---|
| MACE-OFF23(M) | pure density, 109 organic liquids | MAE 0.09 g/mL | [Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2 |
| MACE-OFF23(S) | pure density, 109 organic liquids | MAE 0.23 g/mL (systematically too dense) | [Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2 |
| MACE-OFF23(L) | pure density, 109 organic liquids | no significant improvement over (M) | [Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2 |
| ANI-2x | pure density, 109 organic liquids | MAE and RMSE about twice MACE-OFF23(M) | [Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2 |
| MACE-OFF23 | ΔHvap, 109 organic liquids | systematic offset of about 2 kcal/mol (8 kJ/mol) | [Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2 |
| Organic_MPNICE | pure density, 62 solvents | about 4 % mean error (4.9 % in Table 12) | [Weber et al. 2025](https://arxiv.org/abs/2505.06462), sec. 3.1.4 |
| QRNN | pure density, 62 solvents | about 16 % mean error, systematically high | [Weber et al. 2025](https://arxiv.org/abs/2505.06462), sec. 3.1.4 |
| Hybrid / inorganic MPNICE | pure density, 62 solvents | 18 to 35 % mean absolute error | [Weber et al. 2025](https://arxiv.org/abs/2505.06462), Table 12 |
| ANI-2x, MACE-OFF23(S, M) vs CGenFF | density and ΔHvap of 6 liquids (water, methanol, acetone, benzene, hexane, NMA) | the MLPs miss in both directions (e.g. MACE-OFF23(M) too dense for water, too light for acetone and hexane); CGenFF best overall | [Picha et al. 2025](https://doi.org/10.1021/acs.jcim.5c01355) |
| BAMBOO | density of electrolyte solvent mixtures | about 0.01 g/mL average error, after a density alignment step fitted to experiment | [Gong et al. 2025](https://doi.org/10.1038/s42256-025-01009-7) ([arXiv:2404.07181](https://arxiv.org/abs/2404.07181)) |
| UMA-OMol vs SevenNet-OMat, Orb-OMat | density of ether and carbonate electrolytes | UMA-OMol much closer to experiment than models trained on inorganic data (qualitative in the abstract) | [arXiv:2603.20183](https://arxiv.org/abs/2603.20183) |

No study computing ΔHmix or V^E of molecular liquids with an MLP was found (searches of the web and of Consensus for MLP / neural network potential enthalpy of mixing, excess enthalpy and excess volume, October 2026). The closest work concerns metallic alloys ([arXiv:2406.17499](https://arxiv.org/abs/2406.17499)) or fits mixing data without MD.

- **Density.** Published MLPs differ in pure-liquid density MAE by about 0.05 to 0.25 g/mL. The benchmark resolves 0.0005 g/mL at 2 ns and 0.0013 g/mL at 0.25 ns. Density alone could be ranked from much shorter runs, but the runs serve ΔHmix too.
- **ΔHmix.** There is no MLP baseline. MLP ΔHvap errors are several kJ/mol (MACE-OFF23 offset about 8 kJ/mol; Picha et al. find large ΔHvap deviations for MACE-OFF23(M)). ΔHmix is a difference of such cohesive energies, so much of that error may cancel, but whatever remains is compared with |ΔHmix| near x = 0.5 of 0.6 to 2.6 kJ/mol. If MLPs differ in ΔHmix MAE by even a fifth of Sage's MAE (about 0.11 kJ/mol), 2 ns resolves it.

## Caveats

- **Sage stands in for MLP noise.** MLP fluctuations of h and volume should be similar (the same liquids at the same T and P), but correlation times may differ, especially for MLPs that change conformer barriers. One pilot MLP run per box, analysed with `mixing.box_result`, would check this: compare its sem and n_eff with the Sage `runs.csv`.
- **Slow modes.** k = 1.19 (h, 2 ns) is pooled over boxes. Single boxes have larger run-to-run spreads (for Sage, NMP + 1,2-dichloroethane; see the [Sage report](../sage_summary/report.md)), and an MLP may have slow modes of its own.
- **Systematic simulation errors** (finite size, starting structures) are not noise and are not included. Against the published 1000-molecule Sage runs, the median ΔHmix gap is 0.047 kJ/mol and the largest 0.22 kJ/mol (dmae_methanol_x750).
- **Error profiles.** The model pairs are built from Sage's errors. Real MLPs may spread their errors differently; the closed-form SD shows the result depends mainly on the noise, not on the profile.
