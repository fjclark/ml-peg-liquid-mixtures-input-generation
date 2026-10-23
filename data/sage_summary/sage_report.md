# Sage (openff-2.2.1) reference run

## Summary and systems to inspect

320 runs (64 boxes × 5 replicates) finished. Against experiment, Sage's replicate-mean ΔHmix MAE is 0.581 kJ/mol and its mixture density MAE 0.0120 g/mL.

**Within-run vs between-run uncertainty.** Within-run standard errors come from the autocorrelation (Chodera) estimator on every logged sample. Pooled over all boxes, the SD of the replicate means is 1.26× the within-run sem for h and 1.17× for density (1 expected; p = 2e-08 and 0.0002). 12 (h) and 9 (density) of 64 boxes have p_high < 0.05, against 3.2 expected by chance in each; 3 (h) and 2 (density) have p_low < 0.05. Tests still significant after Holm correction: 1. The replicates start from independently built boxes, so their spread includes slow modes a single run cannot see. MLIP runs are single runs, so their confidence intervals will be narrower than the true run-to-run spread by roughly this factor.

![within vs between](figures/within_between.png)

| Statistic | Observable | Tests | Median ratio | Pooled ratio | p | p_high < 0.05 | p_low < 0.05 |
|---|---|---|---|---|---|---|---|
| box | h | 64 | 1.11 | 1.26 | 2e-08 | 12 | 3 |
| box | density | 64 | 1.06 | 1.17 | 0.0002 | 9 | 2 |
| point | hmix | 39 | 1.08 | 1.20 | 0.0005 | 6 | 2 |
| point | vex | 39 | 0.90 | 0.96 | 0.6 | 3 | 5 |

**Multiple testing.** The pooled ratios are one test per row of the table (4 tests, uncorrected). With a Bonferroni correction over the 4, box h, box density, point hmix stay significant at p < 0.05. The per-box and per-point tests are two-sided and Holm-corrected together (206 tests, flagged at adjusted p < 0.05). The p_high < 0.05 counts above are unadjusted and should be read against the number expected by chance. Boxes listed under "within/between mismatch" below as "p < 0.05 unadjusted only" are places to look, not individually significant results; the evidence that within-run sems are too small comes from the pooled test.

**41 boxes or systems to inspect** (61 items). Categories are in rough order of concern for the benchmark.

### within/between mismatch (30)

Replicates (independently built boxes) scatter more than the within-run standard errors allow. Likely a slow mode (conformers, local structure) that one 2 ns run does not sample, so a single run's CI is too narrow here. Check the replicate time series for offsets that persist all run. Entries marked "unadjusted only" do not survive Holm correction.

| Box / system | Evidence | Plots |
|---|---|---|
| nmp_dichloroethane_x740 | h: replicate SD 2.6× the within-run sem (p_high = 3e-05, Holm-significant) | [replicates](figures/replicates/nmp_dichloroethane_x740.png), [convergence](figures/convergence/nmp_dichloroethane_x740.png) |
| methylpentanoate_octane_x730 | h: replicate SD 2.4× the within-run sem (p_high = 0.0002, p < 0.05 unadjusted only) | [replicates](figures/replicates/methylpentanoate_octane_x730.png), [convergence](figures/convergence/methylpentanoate_octane_x730.png) |
| methylpentanoate_octane_x490 | hmix: replicate SD 2.1× the within-run sem (p_high = 0.001, p < 0.05 unadjusted only) | [replicates](figures/replicates/methylpentanoate_octane_x490.png), [convergence](figures/convergence/methylpentanoate_octane_x490.png) |
| propionitrile_taee_x500 | h: replicate SD 2.1× the within-run sem (p_high = 0.001, p < 0.05 unadjusted only) | [replicates](figures/replicates/propionitrile_taee_x500.png), [convergence](figures/convergence/propionitrile_taee_x500.png) |
| nmp_dichloroethane_x740 | hmix: replicate SD 2.1× the within-run sem (p_high = 0.002, p < 0.05 unadjusted only) | [replicates](figures/replicates/nmp_dichloroethane_x740.png), [convergence](figures/convergence/nmp_dichloroethane_x740.png) |
| pure_n_methylaniline_298 | h: replicate SD 2.0× the within-run sem (p_high = 0.003, p < 0.05 unadjusted only) | [replicates](figures/replicates/pure_n_methylaniline_298.png), [convergence](figures/convergence/pure_n_methylaniline_298.png) |
| heptane_toluene_x250 | h: replicate SD 2.0× the within-run sem (p_high = 0.004, p < 0.05 unadjusted only) | [replicates](figures/replicates/heptane_toluene_x250.png), [convergence](figures/convergence/heptane_toluene_x250.png) |
| heptane_toluene_x250 | density: replicate SD 1.9× the within-run sem (p_high = 0.005, p < 0.05 unadjusted only) | [replicates](figures/replicates/heptane_toluene_x250.png), [convergence](figures/convergence/heptane_toluene_x250.png) |
| heptane_toluene_x750 | hmix: replicate SD 1.9× the within-run sem (p_high = 0.005, p < 0.05 unadjusted only) | [replicates](figures/replicates/heptane_toluene_x750.png), [convergence](figures/convergence/heptane_toluene_x750.png) |
| nmp_dichloroethane_x740 | density: replicate SD 1.9× the within-run sem (p_high = 0.006, p < 0.05 unadjusted only) | [replicates](figures/replicates/nmp_dichloroethane_x740.png), [convergence](figures/convergence/nmp_dichloroethane_x740.png) |
| methylpentanoate_octane_x230 | hmix: replicate SD 1.9× the within-run sem (p_high = 0.007, p < 0.05 unadjusted only) | [replicates](figures/replicates/methylpentanoate_octane_x230.png), [convergence](figures/convergence/methylpentanoate_octane_x230.png) |
| nmp_dichloroethane_x254 | density: replicate SD 1.9× the within-run sem (p_high = 0.008, p < 0.05 unadjusted only) | [replicates](figures/replicates/nmp_dichloroethane_x254.png), [convergence](figures/convergence/nmp_dichloroethane_x254.png) |
| butanone_dibutylamine_x220 | density: replicate SD 1.9× the within-run sem (p_high = 0.008, p < 0.05 unadjusted only) | [replicates](figures/replicates/butanone_dibutylamine_x220.png), [convergence](figures/convergence/butanone_dibutylamine_x220.png) |
| pure_octane_298 | h: replicate SD 1.8× the within-run sem (p_high = 0.009, p < 0.05 unadjusted only) | [replicates](figures/replicates/pure_octane_298.png), [convergence](figures/convergence/pure_octane_298.png) |
| methanol_water_x500 | h: replicate SD 1.8× the within-run sem (p_high = 0.01, p < 0.05 unadjusted only) | [replicates](figures/replicates/methanol_water_x500.png), [convergence](figures/convergence/methanol_water_x500.png) |
| pure_n_methylaniline_298 | density: replicate SD 1.8× the within-run sem (p_high = 0.01, p < 0.05 unadjusted only) | [replicates](figures/replicates/pure_n_methylaniline_298.png), [convergence](figures/convergence/pure_n_methylaniline_298.png) |
| propionitrile_taee_x200 | h: replicate SD 1.8× the within-run sem (p_high = 0.01, p < 0.05 unadjusted only) | [replicates](figures/replicates/propionitrile_taee_x200.png), [convergence](figures/convergence/propionitrile_taee_x200.png) |
| pure_tetrahydropyran_298 | density: replicate SD 1.8× the within-run sem (p_high = 0.01, p < 0.05 unadjusted only) | [replicates](figures/replicates/pure_tetrahydropyran_298.png), [convergence](figures/convergence/pure_tetrahydropyran_298.png) |
| methanol_water_x500 | density: replicate SD 1.8× the within-run sem (p_high = 0.01, p < 0.05 unadjusted only) | [replicates](figures/replicates/methanol_water_x500.png), [convergence](figures/convergence/methanol_water_x500.png) |
| propionitrile_taee_x500 | hmix: replicate SD 1.7× the within-run sem (p_high = 0.02, p < 0.05 unadjusted only) | [replicates](figures/replicates/propionitrile_taee_x500.png), [convergence](figures/convergence/propionitrile_taee_x500.png) |
| methanol_water_x199 | h: replicate SD 1.7× the within-run sem (p_high = 0.02, p < 0.05 unadjusted only) | [replicates](figures/replicates/methanol_water_x199.png), [convergence](figures/convergence/methanol_water_x199.png) |
| nmp_dichloroethane_x740 | vex: replicate SD 1.7× the within-run sem (p_high = 0.02, p < 0.05 unadjusted only) | [replicates](figures/replicates/nmp_dichloroethane_x740.png), [convergence](figures/convergence/nmp_dichloroethane_x740.png) |
| pure_heptane_298 | h: replicate SD 1.7× the within-run sem (p_high = 0.03, p < 0.05 unadjusted only) | [replicates](figures/replicates/pure_heptane_298.png), [convergence](figures/convergence/pure_heptane_298.png) |
| propionitrile_taee_x200 | vex: replicate SD 1.6× the within-run sem (p_high = 0.03, p < 0.05 unadjusted only) | [replicates](figures/replicates/propionitrile_taee_x200.png), [convergence](figures/convergence/propionitrile_taee_x200.png) |
| nmp_dichloroethane_x254 | vex: replicate SD 1.6× the within-run sem (p_high = 0.03, p < 0.05 unadjusted only) | [replicates](figures/replicates/nmp_dichloroethane_x254.png), [convergence](figures/convergence/nmp_dichloroethane_x254.png) |
| water_pyridine_x700 | h: replicate SD 1.6× the within-run sem (p_high = 0.04, p < 0.05 unadjusted only) | [replicates](figures/replicates/water_pyridine_x700.png), [convergence](figures/convergence/water_pyridine_x700.png) |
| pure_water_298 | density: replicate SD 1.6× the within-run sem (p_high = 0.05, p < 0.05 unadjusted only) | [replicates](figures/replicates/pure_water_298.png), [convergence](figures/convergence/pure_water_298.png) |
| pure_benzene_298 | h: replicate SD 1.6× the within-run sem (p_high = 0.05, p < 0.05 unadjusted only) | [replicates](figures/replicates/pure_benzene_298.png), [convergence](figures/convergence/pure_benzene_298.png) |
| water_pyridine_x700 | hmix: replicate SD 1.6× the within-run sem (p_high = 0.05, p < 0.05 unadjusted only) | [replicates](figures/replicates/water_pyridine_x700.png), [convergence](figures/convergence/water_pyridine_x700.png) |
| pure_toluene_298 | density: replicate SD 1.5× the within-run sem (p_high = 0.05, p < 0.05 unadjusted only) | [replicates](figures/replicates/pure_toluene_298.png), [convergence](figures/convergence/pure_toluene_298.png) |

### low sampling (5)

Long correlation times (viscous or slowly reorganising liquids). Check that block means still wander in the convergence plot.

| Box / system | Evidence | Plots |
|---|---|---|
| dmae_methanol_x500 | 1 replicate(s) with n_eff < 100 (minimum 64, τ up to 31 ps) | [replicates](figures/replicates/dmae_methanol_x500.png), [convergence](figures/convergence/dmae_methanol_x500.png) |
| dmae_methanol_x750 | 2 replicate(s) with n_eff < 100 (minimum 48, τ up to 42 ps) | [replicates](figures/replicates/dmae_methanol_x750.png), [convergence](figures/convergence/dmae_methanol_x750.png) |
| pure_2_dimethylamino_ethanol_298 | 3 replicate(s) with n_eff < 100 (minimum 89, τ up to 22 ps) | [replicates](figures/replicates/pure_2_dimethylamino_ethanol_298.png), [convergence](figures/convergence/pure_2_dimethylamino_ethanol_298.png) |
| pure_dibutylamine_298 | 1 replicate(s) with n_eff < 100 (minimum 70, τ up to 28 ps) | [replicates](figures/replicates/pure_dibutylamine_298.png), [convergence](figures/convergence/pure_dibutylamine_298.png) |
| pure_n_methylaniline_298 | 3 replicate(s) with n_eff < 100 (minimum 39, τ up to 51 ps) | [replicates](figures/replicates/pure_n_methylaniline_298.png), [convergence](figures/convergence/pure_n_methylaniline_298.png) |

### start dependence (1)

The distributed box (replicate 0) differs from the rebuilt replicates; the starting structure may bias MLIP runs.

| Box / system | Evidence | Plots |
|---|---|---|
| pure_methyl_pentanoate_298 | h_mean: replicate 0 offset -0.1655 (t = -4.9) | [replicates](figures/replicates/pure_methyl_pentanoate_298.png), [convergence](figures/convergence/pure_methyl_pentanoate_298.png) |

### early density shift (13)

Density still relaxing after the 0.2 ns discard, or a large fluctuation. Compare the replicates' time series.

| Box / system | Evidence | Plots |
|---|---|---|
| pyrrole_pyridine_x508 | rep3: 200–400 ps density differs from the rest by +0.0058 g/mL | [replicates](figures/replicates/pyrrole_pyridine_x508.png), [convergence](figures/convergence/pyrrole_pyridine_x508.png) |
| methanol_water_x500 | rep2: 200–400 ps density differs from the rest by -0.0051 g/mL | [replicates](figures/replicates/methanol_water_x500.png), [convergence](figures/convergence/methanol_water_x500.png) |
| nmp_dichloroethane_x254 | rep4: 200–400 ps density differs from the rest by +0.0059 g/mL | [replicates](figures/replicates/nmp_dichloroethane_x254.png), [convergence](figures/convergence/nmp_dichloroethane_x254.png) |
| dmae_methanol_x247 | rep3: 200–400 ps density differs from the rest by +0.0068 g/mL | [replicates](figures/replicates/dmae_methanol_x247.png), [convergence](figures/convergence/dmae_methanol_x247.png) |
| dmae_methanol_x500 | rep4: 200–400 ps density differs from the rest by +0.0062 g/mL | [replicates](figures/replicates/dmae_methanol_x500.png), [convergence](figures/convergence/dmae_methanol_x500.png) |
| pure_pyrrole_298 | rep3: 200–400 ps density differs from the rest by +0.0054 g/mL | [replicates](figures/replicates/pure_pyrrole_298.png), [convergence](figures/convergence/pure_pyrrole_298.png) |
| pure_methanol_298 | rep4: 200–400 ps density differs from the rest by +0.0052 g/mL | [replicates](figures/replicates/pure_methanol_298.png), [convergence](figures/convergence/pure_methanol_298.png) |
| pure_1_2_dichloroethane_298 | rep4: 200–400 ps density differs from the rest by +0.0117 g/mL | [replicates](figures/replicates/pure_1_2_dichloroethane_298.png), [convergence](figures/convergence/pure_1_2_dichloroethane_298.png) |
| pure_2_dimethylamino_ethanol_298 | rep0: 200–400 ps density differs from the rest by +0.0055 g/mL | [replicates](figures/replicates/pure_2_dimethylamino_ethanol_298.png), [convergence](figures/convergence/pure_2_dimethylamino_ethanol_298.png) |
| pure_n_methylaniline_298 | rep0: 200–400 ps density differs from the rest by -0.0060 g/mL | [replicates](figures/replicates/pure_n_methylaniline_298.png), [convergence](figures/convergence/pure_n_methylaniline_298.png) |
| pure_n_methylaniline_298 | rep3: 200–400 ps density differs from the rest by +0.0052 g/mL | [replicates](figures/replicates/pure_n_methylaniline_298.png), [convergence](figures/convergence/pure_n_methylaniline_298.png) |
| pure_1_propanol_298 | rep3: 200–400 ps density differs from the rest by +0.0056 g/mL | [replicates](figures/replicates/pure_1_propanol_298.png), [convergence](figures/convergence/pure_1_propanol_298.png) |
| pure_toluene_298 | rep3: 200–400 ps density differs from the rest by +0.0054 g/mL | [replicates](figures/replicates/pure_toluene_298.png), [convergence](figures/convergence/pure_toluene_298.png) |

Mean first-block shift over all boxes, by replicate. A shift common to most boxes of one replicate points at that replicate's starting structures rather than at a system:

| Replicate | Mean shift (g/mL) | Boxes with shift > 0 |
|---|---|---|
| 0 | -0.0001 | 50% |
| 1 | -0.0001 | 58% |
| 2 | -0.0003 | 38% |
| 3 | +0.0020 | 81% |
| 4 | +0.0018 | 81% |

### microheterogeneity (3)

Sage itself shows like-like clustering beyond random mixing; MLIP demixing flags for this box compare against this Sage value.

| Box / system | Evidence | Plots |
|---|---|---|
| water_pyridine_x296 | like-neighbour fraction 0.648 vs 0.580 for random mixing | [replicates](figures/replicates/water_pyridine_x296.png), [convergence](figures/convergence/water_pyridine_x296.png) |
| water_pyridine_x495 | like-neighbour fraction 0.604 vs 0.497 for random mixing | [replicates](figures/replicates/water_pyridine_x495.png), [convergence](figures/convergence/water_pyridine_x495.png) |
| water_pyridine_x700 | like-neighbour fraction 0.682 vs 0.579 for random mixing | [replicates](figures/replicates/water_pyridine_x700.png), [convergence](figures/convergence/water_pyridine_x700.png) |

### published-Sage gap (2)

Our replicate mean minus the published openff-2.2.1-ashgc mean (3 replicates, 1000 molecules) at our composition, beyond 0.15 kJ/mol or 0.002 g/mL. Both use openff-2.2.1 parameters with NAGL charges (published: model 0.1.0-rc.3, ours: 1.0.0), so the gap comes from box size, protocol (run length, equilibration, starting conformers) or the NAGL model version. z is the gap over the combined standard error of the two means: it tests whether the gap exceeds replicate noise, not what causes it. The plot shows every matched point (± 2 se), with the thresholds dotted.

| Box / system | Evidence | Plots |
|---|---|---|
| butanone_dibutylamine_x220 | hmix: ours − published -0.2103 (z = -2.5) | [system](figures/systems/butanone_dibutylamine.png) |
| butanone_dibutylamine_x500 | hmix: ours − published -0.1767 (z = -4.3) | [system](figures/systems/butanone_dibutylamine.png) |

![finite size](figures/finite_size.png)

### force-field miss (not a run problem) (7)

Sage is far from experiment; expected for some systems and useful for discriminating models.

| Box / system | Evidence | Plots |
|---|---|---|
| butanone_dibutylamine | mean |ΔHmix error| 0.47 kJ/mol, wrong sign at 1 of 3 x | [system](figures/systems/butanone_dibutylamine.png) |
| dipropylamine_dmf | mean |ΔHmix error| 0.52 kJ/mol | [system](figures/systems/dipropylamine_dmf.png) |
| dmae_methanol | mean |ΔHmix error| 0.61 kJ/mol | [system](figures/systems/dmae_methanol.png) |
| methanol_water | mean |ΔHmix error| 0.80 kJ/mol, wrong sign at 3 of 3 x | [system](figures/systems/methanol_water.png) |
| nmp_dichloroethane | mean |ΔHmix error| 0.53 kJ/mol | [system](figures/systems/nmp_dichloroethane.png) |
| pyrrole_pyridine | mean |ΔHmix error| 1.97 kJ/mol, wrong sign at 1 of 3 x | [system](figures/systems/pyrrole_pyridine.png) |
| water_pyridine | mean |ΔHmix error| 1.26 kJ/mol, wrong sign at 3 of 3 x | [system](figures/systems/water_pyridine.png) |

![Sage vs experiment](figures/sage_vs_exp.png)

Largest within/between mismatches:

![nmp_dichloroethane_x740](figures/replicates/nmp_dichloroethane_x740.png)
![methylpentanoate_octane_x730](figures/replicates/methylpentanoate_octane_x730.png)
![propionitrile_taee_x500](figures/replicates/propionitrile_taee_x500.png)

## Details

- Runs finished: 320 (64 boxes, up to 5 replicates).
- Mixture points with all replicates: 39 of 39.
- Production per run: 2000 ps after 200 ps equilibration.

## Replicate spread (MAZE denominator)

MAZE (mean absolute z-score error) is the ml-peg metric for ΔHmix and mixture density: |model − experiment| / σ_ref averaged over all 39 ΔHmix points (36 density points), giving one number per model. The denominator is σ_ref = sqrt(σ_Sage² + σ_exp²) (the σ_ref columns below, written to `systems.csv` as `sigma_ref_hmix` and `sigma_ref_density`), fixed per point and the same for every model. σ_Sage is one value per mixture: the SD of the replicate means, pooled over its compositions (RMS) and floored at the median over mixtures. σ_ref differs between a mixture's compositions only through σ_exp. The floor is there because model errors are mostly systematic and much larger than σ_ref, so σ_ref acts as a per-mixture weight: without the floor, the mixtures where Sage happens to sample most quietly would dominate MAZE, although an MLIP's noise there need not be as small. With the floor, mixtures are weighted roughly evenly and MAZE stays close to the MAE divided by a single σ.

![replicate spread](figures/replicate_spread.png)

| System | x | ΔHmix Sage | ΔHmix exp | SD | σ_Sage | σ_ref | ρ Sage | ρ exp | ρ SD | σ_ref ρ |
|---|---|---|---|---|---|---|---|---|---|---|
| acetone_nmethylaniline | 0.230 | -0.538 | -0.498 | 0.215 | 0.154 | 0.154 | 0.9554 | 0.9532 | 0.0010 | 0.0009 |
| acetone_nmethylaniline | 0.505 | -0.784 | -0.673 | 0.141 | 0.154 | 0.154 | 0.9045 | 0.9089 | 0.0005 | 0.0008 |
| acetone_nmethylaniline | 0.709 | -0.595 | -0.526 | 0.069 | 0.154 | 0.154 | 0.8559 | 0.8664 | 0.0009 | 0.0008 |
| butanone_dibutylamine | 0.220 | -0.101 | 0.423 | 0.172 | 0.119 | 0.120 | 0.7570 | 0.7604 | 0.0011 | 0.0008 |
| butanone_dibutylamine | 0.500 | 0.055 | 0.619 | 0.088 | 0.119 | 0.120 | 0.7612 | 0.7693 | 0.0007 | 0.0008 |
| butanone_dibutylamine | 0.790 | 0.094 | 0.428 | 0.074 | 0.119 | 0.120 | 0.7665 | 0.7837 | 0.0005 | 0.0008 |
| dipropylamine_dmf | 0.248 | 0.185 | 0.612 | 0.093 | 0.085 | 0.087 | 0.8528 | 0.8687 | 0.0006 | 0.0013 |
| dipropylamine_dmf | 0.500 | 0.123 | 0.748 | 0.058 | 0.085 | 0.088 | 0.8037 | 0.8115 | 0.0008 | 0.0012 |
| dipropylamine_dmf | 0.750 | 0.014 | 0.528 | 0.057 | 0.085 | 0.087 | 0.7639 | 0.7681 | 0.0005 | 0.0012 |
| dmae_methanol | 0.247 | -1.303 | -0.804 | 0.076 | 0.142 | 0.142 | 0.8275 | 0.8437 | 0.0004 | 0.0008 |
| dmae_methanol | 0.500 | -1.624 | -0.865 | 0.085 | 0.142 | 0.142 | 0.8610 | 0.8684 | 0.0010 | 0.0008 |
| dmae_methanol | 0.750 | -1.094 | -0.533 | 0.217 | 0.142 | 0.142 | 0.8804 | 0.8777 | 0.0009 | 0.0008 |
| heptane_toluene | 0.250 | 0.656 | 0.451 | 0.053 | 0.097 | 0.097 | 0.7873 | 0.8039 | 0.0011 | 0.0008 |
| heptane_toluene | 0.500 | 0.775 | 0.556 | 0.095 | 0.097 | 0.097 | 0.7399 | 0.7554 | 0.0007 | 0.0008 |
| heptane_toluene | 0.750 | 0.567 | 0.391 | 0.128 | 0.097 | 0.097 | 0.7003 | 0.7145 | 0.0007 | 0.0009 |
| methanol_water | 0.199 | 0.060 | -0.805 | 0.018 | 0.085 | 0.085 | 0.9135 | 0.9473 | 0.0005 | 0.0008 |
| methanol_water | 0.500 | 0.118 | -0.777 | 0.029 | 0.085 | 0.085 | 0.8388 | 0.8813 | 0.0012 | 0.0008 |
| methanol_water | 0.701 | 0.051 | -0.581 | 0.009 | 0.085 | 0.085 | 0.8031 | 0.8401 | 0.0005 | 0.0008 |
| methylpentanoate_octane | 0.230 | 1.066 | 0.769 | 0.156 | 0.147 | 0.147 | 0.7228 | 0.7323 | 0.0002 | 0.0009 |
| methylpentanoate_octane | 0.490 | 1.418 | 1.053 | 0.163 | 0.147 | 0.147 | 0.7685 | 0.7764 | 0.0003 | 0.0009 |
| methylpentanoate_octane | 0.730 | 0.936 | 0.789 | 0.117 | 0.147 | 0.147 | 0.8179 | 0.8234 | 0.0004 | 0.0009 |
| nmp_dichloroethane | 0.254 | -0.432 | -0.815 | 0.047 | 0.090 | 0.091 | 1.1597 | 1.1831 | 0.0016 | 0.0013 |
| nmp_dichloroethane | 0.500 | -0.504 | -1.264 | 0.039 | 0.090 | 0.090 | 1.1057 | 1.1278 | 0.0008 | 0.0013 |
| nmp_dichloroethane | 0.740 | -0.390 | -0.849 | 0.143 | 0.090 | 0.091 | 1.0533 | 1.0777 | 0.0013 | 0.0013 |
| propanol_benzene | 0.245 | 1.283 | 1.039 | 0.046 | 0.085 | 0.085 | 0.8422 | 0.8563 | 0.0009 | 0.0007 |
| propanol_benzene | 0.500 | 1.373 | 0.968 | 0.086 | 0.085 | 0.085 | 0.8263 | 0.8389 | 0.0004 | 0.0007 |
| propanol_benzene | 0.750 | 0.893 | 0.551 | 0.073 | 0.085 | 0.085 | 0.8098 | 0.8205 | 0.0006 | 0.0007 |
| propionitrile_taee | 0.200 | 0.177 | 0.501 | 0.091 | 0.085 | 0.085 | 0.7634 | 0.7638 | 0.0007 | 0.0008 |
| propionitrile_taee | 0.500 | 0.260 | 0.703 | 0.102 | 0.085 | 0.086 | 0.7676 | 0.7683 | 0.0007 | 0.0008 |
| propionitrile_taee | 0.700 | 0.189 | 0.588 | 0.055 | 0.085 | 0.086 | 0.7712 | 0.7715 | 0.0007 | 0.0008 |
| pyrrole_pyridine | 0.240 | -0.114 | -1.767 | 0.029 | 0.085 | 0.086 | 0.9748 | 0.9754 | 0.0007 | 0.0007 |
| pyrrole_pyridine | 0.508 | -0.071 | -2.650 | 0.028 | 0.085 | 0.087 | 0.9760 | 0.9712 | 0.0005 | 0.0007 |
| pyrrole_pyridine | 0.763 | 0.061 | -1.625 | 0.021 | 0.085 | 0.086 | 0.9769 | 0.9675 | 0.0003 | 0.0007 |
| tetrahydropyran_ethanol | 0.250 | 0.650 | 0.570 | 0.044 | 0.085 | 0.085 | 0.8125 | 0.8189 | 0.0007 | 0.0007 |
| tetrahydropyran_ethanol | 0.504 | 0.968 | 0.851 | 0.060 | 0.085 | 0.086 | 0.8351 | 0.8441 | 0.0008 | 0.0007 |
| tetrahydropyran_ethanol | 0.755 | 0.927 | 0.747 | 0.092 | 0.085 | 0.085 | 0.8519 | 0.8633 | 0.0005 | 0.0007 |
| water_pyridine | 0.296 | 0.133 | -1.127 | 0.034 | 0.085 | 0.092 | 0.9680 | nan | 0.0007 | 0.0007 |
| water_pyridine | 0.495 | 0.064 | -1.281 | 0.045 | 0.085 | 0.094 | 0.9672 | nan | 0.0004 | 0.0007 |
| water_pyridine | 0.700 | 0.074 | -1.097 | 0.045 | 0.085 | 0.091 | 0.9674 | nan | 0.0005 | 0.0007 |

## Finite-size check (our boxes vs published 1000-molecule runs)

Matched points within |Δx| ≤ 0.02: 52. The nearest published value is shifted to our composition along a curve fitted to the published points (two-term Redlich-Kister for ΔHmix, quadratic for density), as a composition offset alone can shift a density by more than the gaps tested. |z| > 3: 4. z uses the standard errors of the two replicate means, which are small, so small systematic offsets reach |z| > 3.

| Box | Property | x ours | x published (nearest) | Ours | Published (nearest) | Published (at our x) | Diff | z |
|---|---|---|---|---|---|---|---|---|
| acetone_nmethylaniline_x230 | hmix | 0.2300 | 0.2316 | -0.5381 | -0.6499 | -0.6473 | +0.1092 | +1.1 |
| acetone_nmethylaniline_x505 | hmix | 0.5050 | 0.5086 | -0.7842 | -0.7554 | -0.7571 | -0.0271 | -0.4 |
| acetone_nmethylaniline_x505 | density | 0.5050 | 0.4948 | 0.9045 | 0.9067 | 0.9045 | -0.0000 | -0.0 |
| acetone_nmethylaniline_x709 | hmix | 0.7091 | 0.7097 | -0.5955 | -0.5850 | -0.5859 | -0.0096 | -0.3 |
| butanone_dibutylamine_x220 | hmix | 0.2200 | 0.2154 | -0.1008 | 0.1071 | 0.1095 | -0.2103 | -2.5 |
| butanone_dibutylamine_x500 | hmix | 0.5000 | 0.5044 | 0.0547 | 0.2323 | 0.2314 | -0.1767 | -4.3 |
| butanone_dibutylamine_x500 | density | 0.5000 | 0.5018 | 0.7612 | 0.7607 | 0.7607 | +0.0006 | +1.8 |
| butanone_dibutylamine_x790 | hmix | 0.7900 | 0.7922 | 0.0943 | 0.1779 | 0.1791 | -0.0847 | -2.2 |
| dipropylamine_dmf_x248 | hmix | 0.2476 | 0.2482 | 0.1854 | 0.1702 | 0.1701 | +0.0153 | +0.3 |
| dipropylamine_dmf_x248 | density | 0.2476 | 0.2523 | 0.8528 | 0.8523 | 0.8533 | -0.0005 | -1.5 |
| dipropylamine_dmf_x500 | hmix | 0.5000 | 0.4995 | 0.1227 | 0.1762 | 0.1761 | -0.0534 | -1.6 |
| dipropylamine_dmf_x500 | density | 0.5000 | 0.5034 | 0.8037 | 0.8026 | 0.8032 | +0.0005 | +1.3 |
| dipropylamine_dmf_x750 | hmix | 0.7500 | 0.7496 | 0.0143 | 0.0921 | 0.0920 | -0.0777 | -2.1 |
| dipropylamine_dmf_x750 | density | 0.7500 | 0.7459 | 0.7639 | 0.7641 | 0.7635 | +0.0004 | +1.5 |
| dmae_methanol_x247 | hmix | 0.2470 | 0.2479 | -1.3032 | -1.3806 | -1.3778 | +0.0746 | +1.7 |
| dmae_methanol_x500 | hmix | 0.5000 | 0.4972 | -1.6240 | -1.6424 | -1.6401 | +0.0161 | +0.3 |
| dmae_methanol_x500 | density | 0.5000 | 0.5020 | 0.8610 | 0.8613 | 0.8611 | -0.0002 | -0.3 |
| dmae_methanol_x750 | hmix | 0.7500 | 0.7506 | -1.0942 | -1.0778 | -1.0799 | -0.0143 | -0.1 |
| methanol_water_x500 | hmix | 0.5000 | 0.5000 | 0.1183 | 0.1297 | 0.1297 | -0.0114 | -0.8 |
| methanol_water_x500 | density | 0.5000 | 0.4997 | 0.8388 | 0.8382 | 0.8382 | +0.0006 | +1.2 |
| methanol_water_x701 | hmix | 0.7014 | 0.7000 | 0.0509 | 0.0527 | 0.0523 | -0.0014 | -0.1 |
| methylpentanoate_octane_x230 | hmix | 0.2300 | 0.2255 | 1.0662 | 1.0335 | 1.0466 | +0.0196 | +0.2 |
| methylpentanoate_octane_x230 | density | 0.2300 | 0.2465 | 0.7228 | 0.7258 | 0.7231 | -0.0003 | -2.4 |
| methylpentanoate_octane_x490 | hmix | 0.4900 | 0.4887 | 1.4176 | 1.3921 | 1.3916 | +0.0260 | +0.4 |
| methylpentanoate_octane_x490 | density | 0.4900 | 0.4896 | 0.7685 | 0.7691 | 0.7692 | -0.0007 | -4.4 |
| methylpentanoate_octane_x730 | hmix | 0.7300 | 0.7323 | 0.9363 | 0.9599 | 0.9661 | -0.0298 | -0.5 |
| methylpentanoate_octane_x730 | density | 0.7300 | 0.7399 | 0.8179 | 0.8198 | 0.8176 | +0.0002 | +0.7 |
| nmp_dichloroethane_x254 | hmix | 0.2544 | 0.2531 | -0.4317 | -0.4647 | -0.4659 | +0.0342 | +1.4 |
| nmp_dichloroethane_x500 | hmix | 0.5000 | 0.5000 | -0.5040 | -0.5261 | -0.5261 | +0.0220 | +0.8 |
| nmp_dichloroethane_x500 | density | 0.5000 | 0.4975 | 1.1057 | 1.1060 | 1.1054 | +0.0003 | +0.6 |
| nmp_dichloroethane_x740 | hmix | 0.7404 | 0.7367 | -0.3904 | -0.3489 | -0.3448 | -0.0456 | -0.7 |
| propanol_benzene_x246 | hmix | 0.2455 | 0.2493 | 1.2830 | 1.2527 | 1.2431 | +0.0399 | +1.8 |
| propanol_benzene_x246 | density | 0.2455 | 0.2385 | 0.8422 | 0.8425 | 0.8421 | +0.0002 | +0.4 |
| propanol_benzene_x500 | hmix | 0.5000 | 0.4996 | 1.3732 | 1.3342 | 1.3338 | +0.0394 | +0.8 |
| propanol_benzene_x500 | density | 0.5000 | 0.5023 | 0.8263 | 0.8260 | 0.8261 | +0.0001 | +0.5 |
| propanol_benzene_x750 | hmix | 0.7500 | 0.7501 | 0.8935 | 0.8653 | 0.8656 | +0.0279 | +0.8 |
| propanol_benzene_x750 | density | 0.7500 | 0.7453 | 0.8098 | 0.8101 | 0.8098 | +0.0000 | +0.0 |
| propionitrile_taee_x200 | hmix | 0.2000 | 0.2009 | 0.1771 | 0.1685 | 0.1679 | +0.0092 | +0.2 |
| propionitrile_taee_x200 | density | 0.2000 | 0.2009 | 0.7634 | 0.7631 | 0.7630 | +0.0004 | +1.1 |
| propionitrile_taee_x500 | hmix | 0.5000 | 0.5025 | 0.2602 | 0.2641 | 0.2642 | -0.0040 | -0.1 |
| propionitrile_taee_x500 | density | 0.5000 | 0.5025 | 0.7676 | 0.7674 | 0.7674 | +0.0003 | +0.7 |
| propionitrile_taee_x700 | hmix | 0.7000 | 0.7012 | 0.1888 | 0.2065 | 0.2071 | -0.0183 | -0.7 |
| propionitrile_taee_x700 | density | 0.7000 | 0.7012 | 0.7712 | 0.7705 | 0.7705 | +0.0007 | +2.0 |
| pyrrole_pyridine_x240 | hmix | 0.2397 | 0.2395 | -0.1136 | -0.0647 | -0.0647 | -0.0489 | -1.7 |
| pyrrole_pyridine_x240 | density | 0.2397 | 0.2273 | 0.9748 | 0.9750 | 0.9751 | -0.0002 | -0.6 |
| pyrrole_pyridine_x508 | hmix | 0.5079 | 0.5095 | -0.0710 | -0.0110 | -0.0114 | -0.0596 | -4.7 |
| pyrrole_pyridine_x763 | hmix | 0.7634 | 0.7602 | 0.0612 | 0.0205 | 0.0207 | +0.0405 | +1.0 |
| pyrrole_pyridine_x763 | density | 0.7634 | 0.7554 | 0.9769 | 0.9778 | 0.9779 | -0.0010 | -3.6 |
| tetrahydropyran_ethanol_x250 | hmix | 0.2500 | 0.2500 | 0.6502 | 0.6490 | 0.6490 | +0.0012 | +0.0 |
| tetrahydropyran_ethanol_x504 | hmix | 0.5042 | 0.5000 | 0.9679 | 0.9609 | 0.9634 | +0.0045 | +0.1 |
| tetrahydropyran_ethanol_x504 | density | 0.5042 | 0.5040 | 0.8351 | 0.8339 | 0.8339 | +0.0012 | +2.6 |
| tetrahydropyran_ethanol_x755 | hmix | 0.7547 | 0.7500 | 0.9273 | 0.8775 | 0.8688 | +0.0585 | +1.2 |

## Convergence

Each run discards 200 ps of equilibration and keeps 2000 ps of production. h is the enthalpy per molecule, (E_pot + pV) / N in kJ/mol. The block tests split production into 10 blocks of 200 ps.

- **Number of independent samples.** The autocorrelation estimator gives n_eff, the number of effectively independent samples of h in a run; production time / n_eff is the correlation time, how long h takes to decorrelate. Median n_eff 548 (correlation time 3.6 ps), minimum 39. 9 runs fall below the ml-peg minimum of 100, where the standard error is less reliable.
- **Cross-check of the standard error.** A second, independent sem comes from the scatter of the 10 block means (SD / √10). Its ratio to the autocorrelation sem should be about 1; well above 1 would mean slow correlations the autocorrelation estimator missed, so the reported sem is too small. Median ratio 1.00 (h), 0.95 (density). With only 10 blocks, single-run ratios scatter by about ±24%, so only the median is meaningful.
- **Equilibration length.** If equilibration were too short, the density would still be relaxing in the first production block (200–400 ps). Mean density of that block minus the mean of the rest of the run: median |shift| 0.0017 g/mL, max 0.0117 g/mL. 13 runs exceed 0.005 g/mL (listed under systems to inspect).
- **Trends within a run.** Two tests on the block means of h and density: drift (linear fit of block mean against time, testing for a nonzero slope) and half-split (t-test of the first 5 blocks against the last 5). Over 1280 tests (runs × 2 observables × 2 tests), about 64 would reach p < 0.05 by chance, so p-values are Holm-corrected, which keeps the chance of any false flag at 0.05. Tests flagged after correction: 0.


### Start dependence (replicate 0 vs re-packed replicates)

Boxes with |t| > 4.0 (prediction interval from replicates 1+). This is a screen, not a flag: about 3.6 of the 128 tests are expected above the threshold by chance.

- pure_methyl_pentanoate_298 h_mean: offset -0.1655, t = -4.9

Plots: `figures/convergence/<box>.png`, `figures/systems/<system>.png` and `figures/replicates/<box>.png` (flagged boxes only).
