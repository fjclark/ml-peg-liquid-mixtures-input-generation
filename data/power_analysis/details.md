# Power analysis: details

Generated 2026-11-06 by `power_analysis.py`. Summary: [report.md](report.md). Tables: `box_noise.csv`, `truncation.csv`, `points.csv`, `metrics.csv`, `bootstrap.csv`, `allocation.csv`, `equilibration.csv`, `thresholds.csv`, `literature.csv`.

## Method

1. **Run length.** Every Sage run is re-analysed with the benchmark's estimator (`mixing.observables` and `timeseries.estimate_mean`, as in `summarise-sage`), keeping 0.25, 0.5, 1, 1.5, 2 ns of production after 0.1 or 0.2 ns of equilibration. The analysis at the protocol length (1 ns) reproduces the Sage `runs.csv` exactly (tested).
2. **Single-run noise per box.** σ = k·sem, with sem the RMS within-run standard error over the 5 replicates and k the pooled ratio of replicate SD to within-run sem over all boxes at that length (k = 1.14 to 1.25 for h, 1.01 to 1.21 for density). This is smoother than the 5-replicate SD of each box and corrects the within-run sem for the slow modes it misses.
3. **Points.** ΔHmix and V^E are linear in the box means (ml-peg's `mixing_matrix`). Each pure box feeds all three compositions of its system, so noise is drawn per box and propagated, which keeps the correlation between points. The propagated ΔHmix noise at 1 ns (RMS over points 0.108 kJ/mol) matches the observed SD of Sage's replicate ΔHmix (RMS 0.122 kJ/mol, Sage `points.csv`).
4. **Metrics.** For each run length, 20000 draws of one run per box give the noise of every metric for a model whose true errors are 1, 0.5, 0.25, 0 × Sage's (replicate mean minus experiment). Two models run independently, so their metrics differ by noise with SD sqrt(2)·SD(metric). The **resolvable difference** is 2.80·sqrt(2)·SD(metric): a true difference that size is detected at the two-sided 5 % level 80 % of the time. A closed form sqrt(sᵀCs)/N (s the signs of the errors over σ_ref for MAZE, C the point covariance) is given as a check; it holds when errors are well above the noise.
5. **System choice.** A cluster bootstrap over systems (20000 resamples of the 12 mixtures, or of the 23 pure liquids) of Sage's noiseless metric measures how much the metric depends on which systems were picked.

## 1. Noise against run length

![noise scaling](figures/noise_scaling.png)

Single-run noise × sqrt(t / 1 ns), RMS over boxes and relative to 1 ns. A flat line would mean σ ∝ 1/sqrt(t). For h the scaled noise grows by 11% from 0.5 to 2 ns: short runs miss some slow fluctuations, so cutting the run length costs a little less precision than 1/sqrt(t) predicts.

| Production (ns) | σ ΔHmix (kJ/mol) | σ density (g/mL) | σ V^E (cm³/mol) | σ pure density (g/mL) | k (h) | k (density) |
|---|---|---|---|---|---|---|
| 0.25 | 0.177 | 0.00184 | 0.224 | 0.00163 | 1.16 | 1.13 |
| 0.5 | 0.144 | 0.00121 | 0.146 | 0.00112 | 1.25 | 1.01 |
| 1 | 0.095 | 0.00087 | 0.113 | 0.00086 | 1.18 | 1.05 |
| 1.5 | 0.079 | 0.00078 | 0.095 | 0.00075 | 1.14 | 1.12 |
| 2 | 0.068 | 0.00074 | 0.088 | 0.00070 | 1.18 | 1.21 |

Median single-run noise per point. k is the replicate SD / within-run sem ratio (1 if the within-run sem were right).

## 2. A model against experiment, per point

Minimum detectable error MDE = 2.80·sqrt(σ_sim² + σ_exp²): an error that size is significant at the 5 % level 80 % of the time. The fraction column shows how many of Sage's errors (a realistic model) exceed the MDE. Experimental errors are 4 to 11× smaller than the simulation noise (median σ_sim over median σ_exp, ΔHmix and density) and barely change these limits.

| Property | Production (ns) | Median MDE | Max MDE | Sage errors above MDE | Median σ_exp | Unit |
|---|---|---|---|---|---|---|
| ΔHmix | 0.5 | 0.404 | 0.671 | 44% | 0.0088 | kJ/mol |
| ΔHmix | 1 | 0.272 | 0.611 | 67% | 0.0088 | kJ/mol |
| ΔHmix | 2 | 0.196 | 0.542 | 67% | 0.0088 | kJ/mol |
| Mixture density | 0.5 | 0.00349 | 0.00508 | 79% | 0.0002 | g/mL |
| Mixture density | 1 | 0.00268 | 0.00383 | 82% | 0.0002 | g/mL |
| Mixture density | 2 | 0.00227 | 0.00353 | 88% | 0.0002 | g/mL |
| V^E | 0.5 | 0.41 | 0.807 | 27% | 0.0041 | cm³/mol |
| V^E | 1 | 0.316 | 0.564 | 27% | 0.0041 | cm³/mol |
| V^E | 2 | 0.246 | 0.478 | 36% | 0.0041 | cm³/mol |
| Pure density | 0.5 | 0.00331 | 0.0063 | 87% | 0.00015 | g/mL |
| Pure density | 1 | 0.0025 | 0.00501 | 91% | 0.00015 | g/mL |
| Pure density | 2 | 0.00214 | 0.00411 | 91% | 0.00015 | g/mL |

![per-point detection](figures/point_mde.png)

ΔHmix points above the diagonal have a Sage error that one run per box detects.

## 3. Ranking two models

Metric noise for a Sage-like model (SD of the score over repeat benchmark runs), the resolvable difference over all error profiles, and the SD of Sage's score over the choice of systems. Units as in the thresholds table of [report.md](report.md).

| Metric | SD (1 ns) | SD closed form | Resolvable (1 ns) | Resolvable (2 ns) | SD over systems |
|---|---|---|---|---|---|
| ΔHmix MAE | 0.021 | 0.023 | 0.047 to 0.082 | 0.037 to 0.065 | 0.15 |
| ΔHmix MAZE | 0.18 | 0.19 | 0.39 to 0.7 | 0.29 to 0.53 | 1.7 |
| Mixture density MAE | 0.00016 | 0.00016 | 0.00038 to 0.00063 | 0.00031 to 0.00051 | 0.0029 |
| Mixture density MAZE | 0.17 | 0.18 | 0.42 to 0.69 | 0.35 to 0.56 | 3 |
| V^E MAE | 0.021 | 0.026 | 0.055 to 0.085 | 0.046 to 0.072 | 0.078 |
| V^E MAZE | 0.17 | 0.2 | 0.43 to 0.67 | 0.35 to 0.57 | 0.72 |
| Pure density MAE | 0.0002 | 0.0002 | 0.00048 to 0.00079 | 0.00039 to 0.00066 | 0.0014 |
| Pure density MAZE | 0.2 | 0.21 | 0.5 to 0.81 | 0.41 to 0.66 | 1.7 |

Noise also biases each score upwards, most for good models: a model whose true ΔHmix MAE is 0.16 kJ/mol scores 0.19 on average at 1 ns, against 0.64 for a true 0.64. Differences between good models are slightly compressed, but their order is kept.

The two kinds of uncertainty answer different questions. The simulation noise says whether rerunning the benchmark would give the same ranking on these 12 systems; at 1 ns it would, for any gap above the resolvable difference. The system spread says how far a model's score depends on this set of systems, which no run length can reduce.

For a ranking, what matters is the system spread of the *difference* between two models, which depends on how alike their errors are. If model B's errors are model A's scaled by (1 + s), the gap is s·MAE_A and its system SD is s·0.15, i.e. 0.23× the gap at any size, so the choice of systems never threatens that ranking. If their errors are unrelated, the gap's system SD approaches sqrt(2)·0.15 = 0.21 kJ/mol, and gaps below about twice that are specific to this set of systems even though a rerun would reproduce them. Real model pairs (e.g. two sizes of one model family) are likely to sit between these cases.

## 4. Run length and cost

| Protocol (ns) | Cost | MACE-OFF23(M) GPU-days | Resolvable ΔHmix MAE (kJ/mol) | Resolvable density MAE (g/mL) | Perfect ΔHmix MAZE |
|---|---|---|---|---|---|
| 0.2 + 0.25 | 38% | 60 | 0.128 | 0.0013 | 1.22 |
| 0.2 + 0.5 | 58% | 94 | 0.112 | 0.0008 | 1.01 |
| 0.2 + 1 | 100% | 161 | 0.082 | 0.0006 | 0.69 |
| 0.2 + 1.5 | 142% | 228 | 0.070 | 0.0005 | 0.56 |
| 0.2 + 2 | 183% | 295 | 0.065 | 0.0005 | 0.50 |

Cost is Σ atoms × simulated time over the 59 boxes (88,640 atoms). GPU-days are an order of magnitude only: MACE-OFF23(M) runs 1.1×10⁶ steps/day for ~600-atom liquid boxes on one A100 ([arXiv:2312.15211](https://arxiv.org/abs/2312.15211), sec. III.3.2), scaled linearly to our 952 to 2,400-atom boxes. Faster and slower models scale the column, not the ratios.

**Equilibration.** The mean over the 100 to 200 ps window minus the mean over the 1 ns production, over all 295 Sage runs, is -0.042 kJ/mol for h (p = 0.002) and +0.0006 g/mL for density (p = 2e-04). Even for Sage, starting from Sage-built boxes, the window is shifted. By replicate, the density shift is -0.0006 (0), +0.0018 (1), +0.0016 (2), +0.0010 (3), -0.0007 (4) g/mL; replicate 0 is the distributed box that MLP runs start from. An MLP also has to relax from Sage's density to its own, which this test cannot see, so 0.2 ns stays. Table: [equilibration.csv](equilibration.csv).

**Per-box run lengths.** For a fixed variance of ΔHmix MAE, mixture density MAE and pure density MAE (equally weighted, each relative to uniform 1 ns), the cheapest allocation runs each box for t ∝ sqrt(A_b / n_b), where A_b is the box's variance contribution and n_b its atom count, with at least 0.5 ns. It ranges from 0.5 to 2.1 ns ([allocation.csv](allocation.csv)); the longest runs go to pure_1_2_dichloroethane_298, pure_n_methylaniline_298, acetone_nmethylaniline_x230, pure_1_propanol_298. It saves 7% of the cost. That assumes σ ∝ 1/sqrt(t), which section 1 shows is optimistic for short runs, so the real saving is smaller. Not worth the extra complexity.

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

- **Density.** Published MLPs differ in pure-liquid density MAE by about 0.05 to 0.25 g/mL. The benchmark resolves 0.0006 g/mL at 1 ns and 0.0013 g/mL at 0.25 ns. Density alone could be ranked from much shorter runs, but the runs serve ΔHmix too.
- **ΔHmix.** There is no MLP baseline. MLP ΔHvap errors are several kJ/mol (MACE-OFF23(M) offset about 8 kJ/mol; Picha et al. find large ΔHvap deviations for MACE-OFF23(M)). ΔHmix is a difference of such cohesive energies, so much of that error may cancel, but whatever remains is compared with |ΔHmix| near x = 0.5 of 0.6 to 2.6 kJ/mol. If MLPs differ in ΔHmix MAE by even a fifth of Sage's MAE (about 0.13 kJ/mol), 1 ns resolves it.

## Caveats

- **Sage stands in for MLP noise.** MLP fluctuations of h and volume should be similar (the same liquids at the same T and P), but correlation times may differ, especially for MLPs that change conformer barriers. One pilot MLP run per box, analysed with `mixing.box_result`, would check this: compare its sem and n_eff with the Sage `runs.csv`.
- **Slow modes.** k = 1.18 (h, 1 ns) is pooled over boxes. Single boxes have larger run-to-run spreads (for Sage, NMP + 1,2-dichloroethane; see the [Sage report](../sage_summary/report.md)), and an MLP may have slow modes of its own.
- **Systematic simulation errors** (finite size, starting structures) are not noise and are not included. Against the published 1000-molecule Sage runs, the median ΔHmix gap is 0.038 kJ/mol and the largest 0.22 kJ/mol (propanol_benzene_x246).
- **Error profiles.** The model pairs are built from Sage's errors. Real MLPs may spread their errors differently; the closed-form SD shows the result depends mainly on the noise, not on the profile.
