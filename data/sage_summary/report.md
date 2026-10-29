# Sage (openff-2.3.0) reference runs

Generated 2026-10-29 by `summarise_sage.py`. Methods and statistics: [details.md](details.md). Every flagged box: [inspection.csv](inspection.csv).

Sage was run through the benchmark protocol (0.2 ns NPT equilibration, 1 ns production) on all 64 boxes, 5 independent replicates each (320 runs). The replicate spread gives σ_Sage, the MAZE denominator, and Sage's own score sets the 'bad' thresholds ([power analysis](../power_analysis/report.md)).

## Sage against experiment

| Property | Points | MAE | MAZE |
|---|---|---|---|
| ΔHmix | 39 | 0.55 kJ/mol | 4.08 |
| Density | 36 | 0.00862 g/mL | 8.71 |
| V^E | 36 | 0.218 cm³/mol | 1.56 |
| Pure density | 11 | 0.00856 g/mL | 7.36 |

The largest ΔHmix errors are for pyrrole_pyridine (1.90 kJ/mol) and water_pyridine (1.20 kJ/mol, wrong sign at 3 of 3 compositions). These are force-field misses, not run problems, and they are the systems most likely to separate models.

![Sage vs experiment](figures/sage_vs_exp.png)

![Sage V^E vs experiment](figures/sage_vs_exp_vex.png)

## Run-to-run noise

Pooled over all boxes, the SD of the replicate means is 1.19× the within-run standard error for h and 1.11× for density (p = 2e-05 and 1e-02). A single 1 ns run's standard error is therefore about 19% too small for h and 11% for density, because independently built boxes sample slow modes (conformers, local structure) that one run does not. This is why σ_Sage comes from the replicate spread rather than within-run errors. Per-box tests that survive Holm correction: none.

![within vs between](figures/within_between.png)

## Boxes to inspect

75 items over 47 boxes or systems, all in [inspection.csv](inspection.csv) with links to their plots. Most are places to look rather than problems:

| Category | Items | Boxes or systems |
|---|---|---|
| within/between mismatch | 33 | 22 |
| low sampling | 19 | 19 |
| start dependence | 1 | 1 |
| early density shift | 9 | 7 |
| microheterogeneity | 4 | 4 |
| published-Sage gap | 3 | 3 |
| force-field miss (not a run problem) | 6 | 6 |

Those that matter for the benchmark:

- **NMP + 1,2-dichloroethane** has slow Sage torsion sampling (see below). It is kept, with σ_Sage set to the median floor.
- **Start dependence**: pure_toluene_298 (density +0.0016 g/mL, t = +4.0). The distributed box (replicate 0) differs from the rebuilt replicates by more than their spread allows, so an MLP started from it may be slightly biased.
- **Published-Sage gap**: 3 of 51 matched points (next section).

## Box size

Our boxes (≥ 100 molecules) were compared with the published 1000-molecule openff-2.3.0 runs ([`fb-fit-v3-single-mean-k100`](https://github.com/openforcefield/ash-sage-rc2/blob/00c9524d38af53725a1a18feb3433f5c191cdb4a/04_benchmark/forcefields/fb-fit-v3-single-mean-k100.offxml?plain=1)) at 51 matched points. Points with |z| > 3: 3. Points that differ by more than 0.15 kJ/mol or 0.002 g/mL: 3, dmae_methanol_x247 (ΔHmix -0.15 kJ/mol, z = -4.0), dmae_methanol_x750 (ΔHmix -0.35 kJ/mol, z = -3.6) and propanol_benzene_x246 (ΔHmix +0.22 kJ/mol, z = +1.7). The median ΔHmix gap is 0.064 kJ/mol and the median density gap 0.0003 g/mL, well below Sage's errors against experiment, so box size does not bias the reference overall; the points above are in [inspection.csv](inspection.csv).

![finite size](figures/finite_size.png)

## NMP + 1,2-dichloroethane

In Sage, only 4% of 1,2-dichloroethane molecules switch between trans and gauche from one 100 ps frame to the next (49% if uncorrelated), a relaxation time of about 1.1 ns. The cause is Sage's barrier, 8 kJ/mol above MACE-OFF23(M), which makes isomerisation about 27× too slow. A model with a reference-like barrier should relax in about 41 ps, well inside the equilibration, so the system is a fair test for MLPs and was kept. Details and the torsion scans: [details.md](details.md#why-nmp--12-dichloroethane-was-kept).
