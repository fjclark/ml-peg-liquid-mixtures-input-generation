# Power analysis: can one run per box rank MLPs?

Generated 2026-10-29 by `power_analysis.py`. Method, full tables, literature and caveats: [details.md](details.md). Thresholds: [thresholds.csv](thresholds.csv), ready to paste into ml-peg as [metrics_snippet.yml](metrics_snippet.yml).

ml-peg runs every MLP once per box: 0.2 ns of NPT equilibration and 2 ns of production at a 1 fs time step. The 320 Sage runs (64 boxes × 5 replicates, [Sage report](../sage_summary/report.md)) stand in for MLP simulation noise. From them we estimate what one run per box can resolve, what it costs, and where ml-peg's 'good' and 'bad' thresholds should sit.

## Answer

- **Between models.** Two models are told apart (two-sided 5 % test, 80 % power) when their ΔHmix MAE differs by 0.04 to 0.075 kJ/mol, or their mixture density MAE by 0.0003 to 0.0005 g/mL. Published MLPs differ in pure-liquid density MAE by about 0.05 to 0.25 g/mL, about 90 to 500× more.
- **Against experiment.** One run per box detects an error at a single point above 0.22 kJ/mol for ΔHmix and 0.0023 g/mL for density (median over points). 67% of Sage's ΔHmix errors are larger.
- **The choice of systems matters more than noise.** Between repeat benchmark runs, a model's ΔHmix MAE moves by 0.019 kJ/mol (SD). Picking a different 13 systems moves Sage's by 0.13 kJ/mol, 7× more (17× for density). Longer runs cannot reduce this.
- **Cost.** 104,076 atoms over 64 boxes, about 347 MACE-OFF23(M) A100 GPU-days per model (order of magnitude). 0.2 + 1 ns costs 55% (saving about 158 GPU-days) and widens the resolvable ΔHmix MAE difference from 0.075 to 0.096 kJ/mol. Per-box run lengths would save only 17%.

We keep 0.2 + 2 ns, and will decide on 0.2 + 1 ns once MLP pilot runs show whether MLP correlation times match Sage's. Equilibration stays at 0.2 ns ([details](details.md#4-shorter-runs)).

![resolvable difference](figures/resolvable.png)

Blue: resolvable difference for a Sage-like error profile (band: 1 to 0 × Sage's errors). Orange: SD of Sage's metric over the choice of systems. Green: the score a model that matches experiment exactly gets from noise alone. Grey: the 'good' threshold.

## Thresholds

Every property is scored by MAZE (each error divided by a fixed σ_ref = sqrt(σ_Sage² + σ_exp²)), with the MAE shown at weight 0. Thresholds come from the error analysis: **good** is the score of a model that matches experiment exactly, with one run per box, plus 2 SD of its noise, so such a model scores 'good' about 98 % of the time. **Bad** is 5× Sage's score, as Sage was trained on these properties.

| Metric | Perfect model (2 ns) | Sage | Good | Bad | Unit |
|---|---|---|---|---|---|
| ΔHmix MAE | 0.071 ± 0.010 | 0.546 | 0.092 | 2.7 | kJ/mol |
| ΔHmix MAZE | 0.65 ± 0.09 | 5.37 | 0.84 | 27 |  |
| Mixture density MAE | 0.00063 ± 0.00008 | 0.00861 | 0.0008 | 0.043 | g/mL |
| Mixture density MAZE | 0.75 ± 0.10 | 10.7 | 0.94 | 54 |  |
| V^E MAE | 0.091 ± 0.013 | 0.224 | 0.12 | 1.1 | cm³/mol |
| V^E MAZE | 0.93 ± 0.14 | 2.38 | 1.2 | 12 |  |
| Pure density MAE | 0.00065 ± 0.00015 | 0.00847 | 0.00096 | 0.042 | g/mL |
| Pure density MAZE | 0.71 ± 0.17 | 8.43 | 1 | 42 |  |

The perfect-model floor rises for shorter runs (ΔHmix MAZE 0.65 at 2 ns, 0.87 at 1 ns), so the thresholds must be recomputed if the protocol changes.

## Main caveat

Sage stands in for MLP noise. MLP fluctuations of h and volume should be similar (the same liquids at the same T and P), but correlation times may differ, especially for MLPs that change conformer barriers. One pilot MLP run per box, analysed with `mixing.box_result`, would check this: compare its sem and n_eff with the Sage `runs.csv`. Other caveats are in [details.md](details.md#caveats).
