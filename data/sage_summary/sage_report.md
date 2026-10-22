# Sage (openff-2.2.1) reference run

## Summary and systems to inspect

320 runs (64 boxes × 5 replicates) finished. Against experiment, Sage's replicate-mean ΔHmix MAE is 0.581 kJ/mol and its mixture density MAE 0.0120 g/mL.

**Within-run vs between-run uncertainty.** Within-run standard errors come from the autocorrelation (Chodera) estimator on 1 ps bin averages. Pooled over all boxes, the SD of the replicate means is 1.23× the within-run sem for h and 1.16× for density (1 expected; p = 4e-07 and 0.0004). 5 (h) and 4 (density) of 64 boxes have p_high < 0.01, against 0.6 expected by chance; none is significantly too small. 0 tests remain significant after Holm correction. The replicates start from independently built boxes, so their spread includes slow modes a single run cannot see. MLIP runs are single runs, so their confidence intervals will be narrower than the true run-to-run spread by roughly this factor.

![within vs between](figures/within_between.png)

| Statistic | Observable | Tests | Median ratio | Pooled ratio | p | p_high < 0.01 | p_low < 0.01 |
|---|---|---|---|---|---|---|---|
| box | h | 64 | 1.08 | 1.23 | 4e-07 | 5 | 0 |
| box | density | 64 | 1.06 | 1.16 | 0.0004 | 4 | 0 |
| point | hmix | 39 | 1.07 | 1.17 | 0.002 | 3 | 0 |
| point | vex | 39 | 0.90 | 0.96 | 0.5 | 0 | 1 |

**39 boxes or systems to inspect** (54 items). Categories are in rough order of concern for the benchmark.

### within/between mismatch (12)

Replicates (independently built boxes) scatter more than the within-run standard errors allow. Likely a slow mode (conformers, local structure) that one 2 ns run does not sample, so a single run's CI is too narrow here. Check the replicate time series for offsets that persist all run.

| Box / system | Evidence | Plots |
|---|---|---|
| butanone_dibutylamine_x220 | density: replicate SD 1.8× the within-run sem (p_high = 0.009, p < 0.01 unadjusted only) | [replicates](figures/replicates/butanone_dibutylamine_x220.png), [convergence](figures/convergence/butanone_dibutylamine_x220.png) |
| heptane_toluene_x250 | h: replicate SD 1.9× the within-run sem (p_high = 0.006, p < 0.01 unadjusted only) | [replicates](figures/replicates/heptane_toluene_x250.png), [convergence](figures/convergence/heptane_toluene_x250.png) |
| heptane_toluene_x250 | density: replicate SD 1.9× the within-run sem (p_high = 0.005, p < 0.01 unadjusted only) | [replicates](figures/replicates/heptane_toluene_x250.png), [convergence](figures/convergence/heptane_toluene_x250.png) |
| methylpentanoate_octane_x730 | h: replicate SD 2.3× the within-run sem (p_high = 0.0002, p < 0.01 unadjusted only) | [replicates](figures/replicates/methylpentanoate_octane_x730.png), [convergence](figures/convergence/methylpentanoate_octane_x730.png) |
| nmp_dichloroethane_x254 | density: replicate SD 1.9× the within-run sem (p_high = 0.008, p < 0.01 unadjusted only) | [replicates](figures/replicates/nmp_dichloroethane_x254.png), [convergence](figures/convergence/nmp_dichloroethane_x254.png) |
| nmp_dichloroethane_x740 | h: replicate SD 2.5× the within-run sem (p_high = 8e-05, p < 0.01 unadjusted only) | [replicates](figures/replicates/nmp_dichloroethane_x740.png), [convergence](figures/convergence/nmp_dichloroethane_x740.png) |
| nmp_dichloroethane_x740 | density: replicate SD 1.9× the within-run sem (p_high = 0.006, p < 0.01 unadjusted only) | [replicates](figures/replicates/nmp_dichloroethane_x740.png), [convergence](figures/convergence/nmp_dichloroethane_x740.png) |
| propionitrile_taee_x500 | h: replicate SD 2.1× the within-run sem (p_high = 0.002, p < 0.01 unadjusted only) | [replicates](figures/replicates/propionitrile_taee_x500.png), [convergence](figures/convergence/propionitrile_taee_x500.png) |
| pure_n_methylaniline_298 | h: replicate SD 2.0× the within-run sem (p_high = 0.003, p < 0.01 unadjusted only) | [replicates](figures/replicates/pure_n_methylaniline_298.png), [convergence](figures/convergence/pure_n_methylaniline_298.png) |
| heptane_toluene_x750 | hmix: replicate SD 1.9× the within-run sem (p_high = 0.008, p < 0.01 unadjusted only) | [replicates](figures/replicates/heptane_toluene_x750.png), [convergence](figures/convergence/heptane_toluene_x750.png) |
| methylpentanoate_octane_x490 | hmix: replicate SD 2.1× the within-run sem (p_high = 0.002, p < 0.01 unadjusted only) | [replicates](figures/replicates/methylpentanoate_octane_x490.png), [convergence](figures/convergence/methylpentanoate_octane_x490.png) |
| nmp_dichloroethane_x740 | hmix: replicate SD 2.0× the within-run sem (p_high = 0.003, p < 0.01 unadjusted only) | [replicates](figures/replicates/nmp_dichloroethane_x740.png), [convergence](figures/convergence/nmp_dichloroethane_x740.png) |

### low sampling (10)

Long correlation times (viscous or slowly reorganising liquids). Check that block means still wander in the convergence plot.

| Box / system | Evidence | Plots |
|---|---|---|
| butanone_dibutylamine_x220 | 2 replicate(s) with n_eff < 100 (minimum 52, τ up to 38 ps) | [replicates](figures/replicates/butanone_dibutylamine_x220.png), [convergence](figures/convergence/butanone_dibutylamine_x220.png) |
| butanone_dibutylamine_x790 | 1 replicate(s) with n_eff < 100 (minimum 81, τ up to 25 ps) | [replicates](figures/replicates/butanone_dibutylamine_x790.png), [convergence](figures/convergence/butanone_dibutylamine_x790.png) |
| dmae_methanol_x247 | 1 replicate(s) with n_eff < 100 (minimum 93, τ up to 21 ps) | [replicates](figures/replicates/dmae_methanol_x247.png), [convergence](figures/convergence/dmae_methanol_x247.png) |
| dmae_methanol_x500 | 1 replicate(s) with n_eff < 100 (minimum 37, τ up to 54 ps) | [replicates](figures/replicates/dmae_methanol_x500.png), [convergence](figures/convergence/dmae_methanol_x500.png) |
| dmae_methanol_x750 | 4 replicate(s) with n_eff < 100 (minimum 31, τ up to 64 ps) | [replicates](figures/replicates/dmae_methanol_x750.png), [convergence](figures/convergence/dmae_methanol_x750.png) |
| propanol_benzene_x750 | 2 replicate(s) with n_eff < 100 (minimum 89, τ up to 22 ps) | [replicates](figures/replicates/propanol_benzene_x750.png), [convergence](figures/convergence/propanol_benzene_x750.png) |
| pure_2_dimethylamino_ethanol_298 | 4 replicate(s) with n_eff < 100 (minimum 56, τ up to 35 ps) | [replicates](figures/replicates/pure_2_dimethylamino_ethanol_298.png), [convergence](figures/convergence/pure_2_dimethylamino_ethanol_298.png) |
| pure_dibutylamine_298 | 2 replicate(s) with n_eff < 100 (minimum 44, τ up to 46 ps) | [replicates](figures/replicates/pure_dibutylamine_298.png), [convergence](figures/convergence/pure_dibutylamine_298.png) |
| pure_dipropylamine_298 | 1 replicate(s) with n_eff < 100 (minimum 96, τ up to 21 ps) | [replicates](figures/replicates/pure_dipropylamine_298.png), [convergence](figures/convergence/pure_dipropylamine_298.png) |
| pure_n_methylaniline_298 | 3 replicate(s) with n_eff < 100 (minimum 25, τ up to 81 ps) | [replicates](figures/replicates/pure_n_methylaniline_298.png), [convergence](figures/convergence/pure_n_methylaniline_298.png) |

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

### box-size gap (8)

Differs from the published 1000-molecule runs beyond the threshold. Could be box size or NAGL vs AM1-BCC charges.

| Box / system | Evidence | Plots |
|---|---|---|
| butanone_dibutylamine_x220 | hmix: ours − published -0.2629 (z = -3.1) | [system](figures/systems/butanone_dibutylamine.png) |
| butanone_dibutylamine_x500 | hmix: ours − published -0.2108 (z = -3.9) | [system](figures/systems/butanone_dibutylamine.png) |
| dipropylamine_dmf_x750 | hmix: ours − published -0.1560 (z = -6.1) | [system](figures/systems/dipropylamine_dmf.png) |
| dmae_methanol_x500 | hmix: ours − published -0.2461 (z = -4.5) | [system](figures/systems/dmae_methanol.png) |
| dmae_methanol_x750 | hmix: ours − published -0.2423 (z = -2.4) | [system](figures/systems/dmae_methanol.png) |
| methylpentanoate_octane_x230 | density: ours − published -0.0028 (z = -12.2) | [system](figures/systems/methylpentanoate_octane.png) |
| methylpentanoate_octane_x730 | density: ours − published -0.0021 (z = -9.3) | [system](figures/systems/methylpentanoate_octane.png) |
| nmp_dichloroethane_x500 | density: ours − published -0.0031 (z = -6.0) | [system](figures/systems/nmp_dichloroethane.png) |

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

σ_Sage is the SD of the replicate means, pooled per system (RMS over its compositions) and floored at the median over systems. σ_ref = sqrt(σ_Sage² + σ_exp²).

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

Matched points within |Δx| ≤ 0.02: 52. |z| > 3: 20. z uses the standard errors of the two replicate means, which are small, so small systematic offsets reach |z| > 3.

| Box | Property | Ours | Published | Diff | z |
|---|---|---|---|---|---|
| acetone_nmethylaniline_x230 | hmix | -0.5381 | -0.6150 | +0.0770 | +0.8 |
| acetone_nmethylaniline_x505 | hmix | -0.7842 | -0.8762 | +0.0920 | +1.4 |
| acetone_nmethylaniline_x505 | density | 0.9045 | 0.9061 | -0.0016 | -6.0 |
| acetone_nmethylaniline_x709 | hmix | -0.5955 | -0.7228 | +0.1274 | +3.7 |
| butanone_dibutylamine_x220 | hmix | -0.1008 | 0.1621 | -0.2629 | -3.1 |
| butanone_dibutylamine_x500 | hmix | 0.0547 | 0.2655 | -0.2108 | -3.9 |
| butanone_dibutylamine_x500 | density | 0.7612 | 0.7606 | +0.0007 | +2.0 |
| butanone_dibutylamine_x790 | hmix | 0.0943 | 0.1749 | -0.0805 | -1.9 |
| dipropylamine_dmf_x248 | hmix | 0.1854 | 0.1509 | +0.0345 | +0.7 |
| dipropylamine_dmf_x248 | density | 0.8528 | 0.8537 | -0.0009 | -2.8 |
| dipropylamine_dmf_x500 | hmix | 0.1227 | 0.1503 | -0.0275 | -1.1 |
| dipropylamine_dmf_x500 | density | 0.8037 | 0.8030 | +0.0007 | +1.8 |
| dipropylamine_dmf_x750 | hmix | 0.0143 | 0.1703 | -0.1560 | -6.1 |
| dipropylamine_dmf_x750 | density | 0.7639 | 0.7642 | -0.0003 | -1.0 |
| dmae_methanol_x247 | hmix | -1.3032 | -1.1534 | -0.1498 | -3.3 |
| dmae_methanol_x500 | hmix | -1.6240 | -1.3779 | -0.2461 | -4.5 |
| dmae_methanol_x500 | density | 0.8610 | 0.8596 | +0.0013 | +3.0 |
| dmae_methanol_x750 | hmix | -1.0942 | -0.8519 | -0.2423 | -2.4 |
| methanol_water_x500 | hmix | 0.1183 | 0.1242 | -0.0059 | -0.3 |
| methanol_water_x500 | density | 0.8388 | 0.8384 | +0.0004 | +0.8 |
| methanol_water_x701 | hmix | 0.0509 | 0.0492 | +0.0017 | +0.2 |
| methylpentanoate_octane_x230 | hmix | 1.0662 | 1.0873 | -0.0212 | -0.3 |
| methylpentanoate_octane_x230 | density | 0.7228 | 0.7256 | -0.0028 | -12.2 |
| methylpentanoate_octane_x490 | hmix | 1.4176 | 1.3225 | +0.0952 | +1.3 |
| methylpentanoate_octane_x490 | density | 0.7685 | 0.7693 | -0.0008 | -5.2 |
| methylpentanoate_octane_x730 | hmix | 0.9363 | 1.0154 | -0.0791 | -1.2 |
| methylpentanoate_octane_x730 | density | 0.8179 | 0.8200 | -0.0021 | -9.3 |
| nmp_dichloroethane_x254 | hmix | -0.4317 | -0.5099 | +0.0781 | +3.2 |
| nmp_dichloroethane_x500 | hmix | -0.5040 | -0.6079 | +0.1039 | +3.4 |
| nmp_dichloroethane_x500 | density | 1.1057 | 1.1088 | -0.0031 | -6.0 |
| nmp_dichloroethane_x740 | hmix | -0.3904 | -0.4535 | +0.0631 | +1.0 |
| propanol_benzene_x246 | hmix | 1.2830 | 1.2108 | +0.0723 | +1.6 |
| propanol_benzene_x246 | density | 0.8422 | 0.8411 | +0.0011 | +2.5 |
| propanol_benzene_x500 | hmix | 1.3732 | 1.2648 | +0.1084 | +2.0 |
| propanol_benzene_x500 | density | 0.8263 | 0.8255 | +0.0008 | +4.7 |
| propanol_benzene_x750 | hmix | 0.8935 | 0.7647 | +0.1288 | +2.4 |
| propanol_benzene_x750 | density | 0.8098 | 0.8106 | -0.0008 | -3.1 |
| propionitrile_taee_x200 | hmix | 0.1771 | 0.1845 | -0.0074 | -0.1 |
| propionitrile_taee_x200 | density | 0.7634 | 0.7632 | +0.0002 | +0.6 |
| propionitrile_taee_x500 | hmix | 0.2602 | 0.2827 | -0.0225 | -0.5 |
| propionitrile_taee_x500 | density | 0.7676 | 0.7684 | -0.0008 | -2.3 |
| propionitrile_taee_x700 | hmix | 0.1888 | 0.2811 | -0.0924 | -3.5 |
| propionitrile_taee_x700 | density | 0.7712 | 0.7717 | -0.0005 | -1.5 |
| pyrrole_pyridine_x240 | hmix | -0.1136 | -0.1439 | +0.0303 | +2.0 |
| pyrrole_pyridine_x240 | density | 0.9748 | 0.9760 | -0.0012 | -3.4 |
| pyrrole_pyridine_x508 | hmix | -0.0710 | -0.1193 | +0.0484 | +3.2 |
| pyrrole_pyridine_x763 | hmix | 0.0612 | -0.0080 | +0.0692 | +5.7 |
| pyrrole_pyridine_x763 | density | 0.9769 | 0.9772 | -0.0004 | -1.8 |
| tetrahydropyran_ethanol_x250 | hmix | 0.6502 | 0.6431 | +0.0071 | +0.1 |
| tetrahydropyran_ethanol_x504 | hmix | 0.9679 | 0.9341 | +0.0338 | +1.1 |
| tetrahydropyran_ethanol_x504 | density | 0.8351 | 0.8348 | +0.0003 | +0.8 |
| tetrahydropyran_ethanol_x755 | hmix | 0.9273 | 0.8482 | +0.0791 | +1.7 |

## Convergence

- Correlation time of h (production / n_eff): median 6.4 ps; minimum n_eff 25 (below 100: 21 runs).
- Block-average sem / autocorrelation sem: median 0.99 (h), 0.95 (density).
- Density shift of the first production block (200–400 ps) relative to the rest: median 0.00169, max 0.01166 g/mL.
- Drift and half-split tests flagged after Holm correction (p < 0.01): 0.


### Start dependence (replicate 0 vs re-packed replicates)

Boxes with |t| > 4.0 (prediction interval from replicates 1+). This is a screen, not a flag: about 3.6 of the 128 tests are expected above the threshold by chance.

- pure_methyl_pentanoate_298 h_mean: offset -0.1655, t = -4.9

Plots: `figures/convergence/<box>.png`, `figures/systems/<system>.png` and `figures/replicates/<box>.png` (flagged boxes only).
