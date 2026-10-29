# Decisions

The decisions behind the ml-peg `binary_mixtures` benchmark: 13 binary mixtures at 3 compositions (64 boxes with the pure liquids), built and run with Sage 2.3.0, scored by MAZE. The top section is the current state, one bullet per decision with its reason and evidence. The [history](#history) at the bottom lists every dated change, including those since superseded. Decisions were made by the benchmark author, many on Claude's recommendation, in planning sessions from 11 to 29 October 2026.

The rules are implemented in `src/liquid_mixtures_inputs/` and listed in [selection/report.md](selection/report.md). Data problems found along the way are in [docs/data_issues.md](docs/data_issues.md).

## Current decisions

### Scope and data

- **Properties.** ΔHmix, mixture density and V^E of binary liquid mixtures, plus pure densities not already scored by ml-peg. ml-peg had no ΔHmix benchmark and no excess properties beyond the ethanol-water density curve.
- **Three compositions per system** (near x = 0.25, 0.5 and 0.75), even at about 1.7× the liquid_densities cost, because composition dependence is the point of ΔHmix.
- **Maximally informative systems.** Every system needs ΔHmix and density at matching compositions, pure densities, and experimental uncertainties throughout, so each box scores several properties.
- **No overlap with ml-peg.** The ethanol-water density curve, water density and the 63 liquid_densities liquids are not re-scored; their pure boxes are still simulated, as ΔHmix needs them.
- **Sources.** The ash-sage-rc2 (Sage 2.3 refit) training and validation sets are the main pool, because they carry experimental uncertainties; the Sage 2.0 set has none. Every source is pinned to a commit.
- **Full isotherms from the source paper.** Isotherms and pure densities come from each mixture's own paper through the NIST ThermoML archive, and every pure density is cross-checked against ml-peg, Sage 2.0 and DIPPR/NIST (differences above 0.002 g/mL are flagged). The rc2 sets keep only 3 points per isotherm and include wrong pure densities (DMF, 1-propanol, piperidine) that flip the sign of V^E ([data_issues.md, 2.3](docs/data_issues.md#23-pure-densities-that-disagree-with-other-sources-verified-against-two-or-more-independent-papers)).
- **V^E uncertainties.** V^E is derived from the pure densities of the same paper. The density σ is treated as an error shared by every density in that paper (it cancels in V^E), and the paper's RMS residual about its Redlich-Kister fit gives the random part. Propagating the flat TRC σ as independent had made V^E error bars up to 40× the scatter (dipropylamine + DMF: 0.13 against 0.003 cm³/mol).
- **Methanol + water** densities from 10.1016/j.jct.2007.05.004 only, so its pure methanol is used for both V^E and the mixture densities; water is IAPWS-95.
- **Structure checks.** Every component must have zero formal charge, no enol or imidic acid motif, and match RDKit's canonical tautomer. The Sage data store 2-pyrrolidone as the wrong tautomer, and a wrong starting species cannot convert back in a 1 ns run (acetone's enol is ~10⁻⁸ of the equilibrium mixture, [Chiang et al. 1984](https://doi.org/10.1021/ja00314a055)).

### Filters

36 of the 248 systems with enough data pass every filter ([funnel](selection/report.md#screening-funnel)).

- **Elements** H, C, N, O and Cl only. Cl is covered by MACE-OFF, AIMNet2 and AceFF.
- **Size and dynamics.** ≤ 10 heavy atoms per molecule, and no viscous liquids (polyols, ethanolamines, polyamines), for cost and convergence in 1 ns of production.
- **No likely two-phase systems.** In a small box a model can stay mixed where the real system separates, or the reverse. Screens: H^E linear in x (the lever-rule signature, tested only when the largest |H^E| ≥ 0.4 kJ/mol, as ideal mixtures are also linear), H^E(x = 0.5) ≤ RT/2, and a literature miscibility check for every endothermic system.
- **Signal.** |ΔHmix(x ≈ 0.5)| ≥ 0.4 kJ/mol, well above the single-point detection limit.
- **Sage convergence.** Published Sage 3-replicate SD ≤ 0.10 kJ/mol (ΔHmix) and ≤ 0.001 g/mL (density), as a large spread suggests slow equilibration.
- **Uncertainties.** Relative ΔHmix uncertainty ≤ 15 % and pure-density uncertainties required, which removes placeholder errors (flat 25 %).

### Systems

- **The 13 systems** are listed with their data in [selection/report.md](selection/report.md#selected-systems). On review (2026-10-14) the set was reshaped to have fewer alcohols, a non-polar π control (heptane + toluene), a wider ΔHmix range (NMP + 1,2-dichloroethane, -1.26 kJ/mol), a ketone + amine pair that is endothermic in one case and exothermic in the other (2-butanone + dibutylamine against acetone + N-methylaniline), and O-H···π (1-propanol + benzene).
- **Aqueous systems.** Methanol + water (the smallest Sage spread in the pool, so it calibrates the noise) and water + pyridine at 303.15 K, ΔHmix only (the aromatic N-water motif of histidine and nucleobases; ThermoML mislabels the temperature, [data_issues.md, 2.2](docs/data_issues.md#22-pyridine--water-δhmix-temperature-mislabelled-verified)). Excluded: ethanol + water (already in ml-peg), 2-(dimethylamino)ethanol + water (large Sage SD), 2-pyrrolidone + water (viscous, wrong tautomer, placeholder errors), and acetophenone, acetylacetone and ethyl acetoacetate + water (likely two-phase).
- **Piperidine + water dropped.** Nitrogen inversion is frozen in Sage even at 500 K, so ΔHmix depended on the starting N-H orientation by ~0.3 kJ/mol.
- **NMP + 1,2-dichloroethane kept.** Its slow DCE torsion sampling is a Sage artefact: Sage's trans to gauche barrier is 25.4 kJ/mol against 17.3 kJ/mol for MACE-OFF23(M), which makes Sage's flips about 27× too slow ([Sage report](data/sage_summary/report.md#nmp--12-dichloroethane)). A model with a reference-like barrier relaxes in about 40 ps, inside the equilibration. ΔHmix comes from 10.1016/j.jct.2008.01.024 only, at its points nearest x = 0.2, 0.5 and 0.8, as the two papers disagree away from equimolar.
- **DMAE + methanol kept** for now: the anneal reaches the same conformer populations from different starts (see [open questions](#open-questions)).
- **Known gaps**, stated in the report: no N-H amide system and no amide-water system, and the set reflects ThermoML and Sage coverage. Reserves are kept in `selection.yml`.

### Boxes and protocol

- **Box size.** Edge ≥ 25 Å, ≥ 100 molecules and ≥ 20 of the minority component (952 to 2,800 atoms). Half of a 21.5 Å box is below the ~12 Å receptive field of a 2-layer, 6 Å MLP, precision per GPU-hour does not depend on box size, and 1000-molecule boxes would cost 10 to 15× more. Sage results in our boxes match the published 1000-molecule runs ([Sage report](data/sage_summary/report.md#box-size)).
- **Box build.** Packed at the experimental density with Sage 2.3.0, each molecule starting from the lowest-MMFF of 10 ETKDG conformers, annealed (450 K cooled to T over 1.6 ns) and equilibrated for 500 ps NVT. The first Sage runs showed slow conformer equilibria (tetrahydropyran, DCE, DMAE) and dependence on the starting conformer.
- **Protocol.** One run per model and box: 0.2 ns NPT equilibration and 1 ns production, for cost: 0.2 + 2 ns would cost 1.8× as much (about 190 against 350 MACE-OFF23(M) A100 GPU-days per model) and only narrow the resolvable ΔHmix MAE difference from 0.096 to 0.075 kJ/mol, small next to the spread over the choice of systems (0.13 kJ/mol). Equilibration stays at 0.2 ns, as the 100 to 200 ps window is still shifted for Sage ([power analysis](data/power_analysis/report.md)).
- **One estimator** for Sage and the MLPs: ml-peg's `timeseries.py` (autocorrelation) on every logged sample, without binning.

### Sage reference

- **Sage 2.3.0 throughout**, 5 replicates per box: replicate 0 is the distributed box, replicates 1 to 4 are independently re-packed. Each run is 2.2 ns, so the power analysis can compare the protocol with 2 ns of production, but every Sage summary (σ_Sage, the reference values and Sage's scores) uses only the protocol's window, 0.2 to 1.2 ns. The replicate spread gives σ_Sage. Mixing versions would mix parameters: 2.3.0 liquids are 0.01 to 0.03 g/mL denser than 2.2.1.
- **Finite-size check** against the published `fb-fit-v3-single-mean-k100` runs (2.3.0rc2, parameters identical to the released 2.3.0; only the NAGL model differs).
- **Significance.** Family-wise 0.05 (`FLAG_P`), Holm-adjusted for per-box tests and Bonferroni over the pooled tests.

### Metrics and thresholds

- **MAZE for every property** (ΔHmix, mixture density, V^E and pure density), with the MAE shown at weight 0. V^E uses MAE, not RMSE.
- **Denominator.** σ_ref = sqrt(σ_Sage² + σ_exp²). σ_Sage is the replicate SD pooled (RMS) over a system's compositions and floored at the median over systems, so that no system dominates MAZE because Sage happens to sample it quietly. NMP + DCE and pure DCE are set to the floor on all four properties, as their spread comes from Sage's slow torsion.
- **Thresholds from the error analysis.** 'Good' is the score of a model that matches experiment exactly, with one 0.2 + 1 ns run per box, plus 2 SD of its noise. 'Bad' is 5× Sage's score, as Sage was trained on these properties. Mixture and pure densities share one set of thresholds for simplicity, the larger of their two values. Values in [thresholds.csv](data/power_analysis/thresholds.csv), ready to paste as [metrics_snippet.yml](data/power_analysis/metrics_snippet.yml). They must be recomputed if the protocol changes.

### Open questions

- **DMAE + methanol.** Starting MLP runs from Sage-equilibrated conformer populations could hide a model that favours the wrong conformer, if that model also interconverts slowly. Options: keep it and add a conformer-population diagnostic, or drop it. (Tetrahydropyran + ethanol is unaffected; its equilibrium is essentially all chair.)
- **MLP correlation times.** The power analysis uses Sage as a proxy for MLP noise; MLP pilot runs should check this.
- **Within-run standard errors** of 1 ns runs are 19 % (h) and 11 % (density) too small on average; changing the estimator is not decided.
- **Fluorinated systems** (2,2,2-trifluoroethanol mixtures, which flip the sign of ΔHmix against their analogues) were suggested but not adopted; no decision is recorded.
- **Issue #936** (which asks for ~1000 molecules per box) is left to the benchmark author.

## History

Superseded decisions are marked ~~struck through~~ with what replaced them.

- 2026-10-11: Benchmark ΔHmix, mixture density and V^E from the Sage training and validation sets, starting from the recommended fast-converging list. About 12 systems, 3 compositions each. Elements H, C, N, O and Cl. Exclude likely two-phase systems.
- 2026-10-14: Choose maximally informative systems; avoid overlap with ml-peg; use the ash-sage-rc2 sets and pin all sources; take isotherms and pure densities from ThermoML and cross-check them; fix the 2-pyrrolidone tautomer.
- 2026-10-14: Hard filters: lever-rule screen and RT/2 bound, Sage 3-replicate SD, |ΔHmix(0.5)| ≥ 0.4 kJ/mol, no viscous liquids, ≤ 10 heavy atoms. Aqueous systems chosen (methanol + water, water + pyridine at 303.15 K) and others excluded. Set reshaped for diversity; reserves listed. Fluorinated systems suggested, not adopted.
- 2026-10-14: Box edge ≥ 25 Å, ~~≥ 150 molecules and ≥ 30 minority~~ (lowered on 2026-10-19). ~~Build boxes with Sage 2.2.1~~ (2.3.0 since 2026-10-28). Run Sage 5 times per box for σ_Sage, instead of scaling the published SDs.
- 2026-10-19: Lever-rule test only when the largest |H^E| ≥ 0.4 kJ/mol. Relative ΔHmix uncertainty ≤ 15 % and pure-density uncertainties required. Minority floor ≥ 20 and total ≥ 100 molecules. N-methylformamide + o-toluidine not added. Selection frozen at 14 systems.
- 2026-10-20: Drop piperidine + water (14 to 13 systems, 68 to 64 boxes). Build from the lowest-MMFF of 10 conformers with an anneal; rebuild every box and rerun Sage. Keep DMAE + methanol and NMP + DCE.
- 2026-10-22: One estimator for Sage and the MLPs (ml-peg's `timeseries.py`), ~~after 1 ps binning~~ (dropped on 2026-10-23). The within-run sem being low is left open.
- 2026-10-23: Check every curated component's charge and tautomer. V^E from each paper's own pure densities, with shared-error propagation. Methanol + water densities from one paper. Drop the 1 ps binning (median raw / binned sem 1.00). Keep NMP + DCE; its ΔHmix from one paper at x = 0.2, 0.5 and 0.8, with new boxes; σ_Sage of NMP + DCE at the floor. `FLAG_P` raised from 0.01 to 0.05. Keep 0.2 + 2 ns and report 0.2 + 1 ns as an option. ~~Flag the V^E RMSE and ΔHmix MAZE thresholds as below their noise floors~~ (thresholds set from the error analysis on 2026-10-29). NMF candidates moved to the report's reserves.
- 2026-10-28: Use Sage 2.3.0 throughout (boxes, reference runs and the finite-size comparison); remove every 2.2.1 run. Selection still uses published results from several Sage versions.
- 2026-10-29: ~~Raise the V^E RMSE 'good' threshold from 0.05 to 0.15 cm³/mol and keep ΔHmix MAZE 'good' at 1.0~~ (superseded the same day). Score every property by MAZE with MAE at weight 0; V^E moves from RMSE to MAE; σ_Sage and σ_ref added for V^E and pure densities, with NMP + DCE and pure DCE at the floor; 'good' = perfect-model score + 2 SD, 'bad' = 5× Sage. Mixture and pure densities share their thresholds (the larger of the two values).
- 2026-10-29: ~~Keep 0.2 + 2 ns~~: run 0.2 + 1 ns for cost. The Sage summaries, σ_Sage and the thresholds are recomputed from the first 0.2 + 1 ns of the existing 2.2 ns Sage runs.
