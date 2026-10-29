# Power analysis: can one run per box rank MLPs?

Generated 2026-10-29 by `power_analysis.py`. Method, full tables, literature and caveats: [details.md](details.md). Thresholds: [thresholds.csv](thresholds.csv), ready to paste into ml-peg as [metrics_snippet.yml](metrics_snippet.yml).

ml-peg runs every MLP once per box: 0.2 ns of NPT equilibration and 1 ns of production at a 1 fs time step. The 320 Sage runs (64 boxes × 5 replicates, [Sage report](../sage_summary/report.md)) stand in for MLP simulation noise. From them we estimate what one run per box can resolve, what it costs, and where ml-peg's 'good' and 'bad' thresholds should sit.

## Answer

- **Between models.** Two models are told apart (two-sided 5 % test, 80 % power) when their ΔHmix MAE differs by 0.05 to 0.096 kJ/mol, or their mixture density MAE by 0.0004 to 0.0007 g/mL. Published MLPs differ in pure-liquid density MAE by about 0.05 to 0.25 g/mL, about 80 to 400× more.
- **Against experiment.** One run per box detects an error at a single point above 0.31 kJ/mol for ΔHmix and 0.0028 g/mL for density (median over points). 62% of Sage's ΔHmix errors are larger.
- **The choice of systems matters more than noise.** Between repeat benchmark runs, a model's ΔHmix MAE moves by 0.024 kJ/mol (SD). Picking a different 13 systems moves Sage's by 0.13 kJ/mol, 6× more (14× for density). Longer runs cannot reduce this.
- **Cost.** 104,076 atoms over 64 boxes, about 189 MACE-OFF23(M) A100 GPU-days per model (order of magnitude). 0.2 + 2 ns would cost 1.8× as much (about 158 GPU-days more) and narrow the resolvable ΔHmix MAE difference from 0.096 to 0.075 kJ/mol. Per-box run lengths would save only 13%.

We run 0.2 + 1 ns for cost. The extra precision of 0.2 + 2 ns is small next to the spread over the choice of systems. Equilibration stays at 0.2 ns ([details](details.md#4-run-length-and-cost)).

![resolvable difference](figures/resolvable.png)

Blue: resolvable difference for a Sage-like error profile (band: 1 to 0 × Sage's errors). Orange: SD of Sage's metric over the choice of systems. Green: the score a model that matches experiment exactly gets from noise alone. Grey: the 'good' threshold.

## Thresholds

Every property is scored by MAZE (each error divided by a fixed σ_ref = sqrt(σ_Sage² + σ_exp²)), with the MAE shown at weight 0. Thresholds come from the error analysis: **good** is the score of a model that matches experiment exactly, with one run per box, plus 2 SD of its noise, so such a model scores 'good' about 98 % of the time. **Bad** is 5× Sage's score, as Sage was trained on these properties. Mixture and pure densities share one set of thresholds for simplicity, the larger of their two values, so both rules still hold for each.

| Metric | Perfect model (1 ns) | Sage | Good | Bad | Unit |
|---|---|---|---|---|---|
| ΔHmix MAE | 0.096 ± 0.014 | 0.55 | 0.12 | 2.8 | kJ/mol |
| ΔHmix MAZE | 0.64 ± 0.09 | 4.08 | 0.82 | 20 |  |
| Mixture density MAE | 0.00078 ± 0.00010 | 0.00862 | 0.0011 | 0.043 | g/mL |
| Mixture density MAZE | 0.78 ± 0.10 | 8.71 | 1.1 | 44 |  |
| V^E MAE | 0.11 ± 0.02 | 0.218 | 0.14 | 1.1 | cm³/mol |
| V^E MAZE | 0.76 ± 0.11 | 1.56 | 0.98 | 7.8 |  |
| Pure density MAE | 0.00077 ± 0.00018 | 0.00856 | 0.0011 | 0.043 | g/mL |
| Pure density MAZE | 0.72 ± 0.17 | 7.36 | 1.1 | 44 |  |

The perfect-model floor depends on the run length (ΔHmix MAZE 0.64 at 1 ns, 0.48 at 2 ns), so the thresholds must be recomputed if the protocol changes.

## Main caveat

Sage stands in for MLP noise. MLP fluctuations of h and volume should be similar (the same liquids at the same T and P), but correlation times may differ, especially for MLPs that change conformer barriers. One pilot MLP run per box, analysed with `mixing.box_result`, would check this: compare its sem and n_eff with the Sage `runs.csv`. Other caveats are in [details.md](details.md#caveats).
