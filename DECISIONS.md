# Decisions: selecting the binary mixtures

A record of the decisions behind the 13 systems (64 boxes) in the ml-peg
`binary_mixtures` benchmark, with dates, who made each one, and why. It was
compiled on 2026-10-22 from the planning conversations (Claude Code sessions,
11 to 22 October 2026), `src/liquid_mixtures_inputs/selection.yml`,
`selection/report.md` and `data/sage_summary/sage_report.md`.

"User" means a decision stated or chosen by the benchmark author.
"Recommended" means a recommendation from Claude that the user accepted.
Quotes are the user's own words.

The rules themselves are implemented in `select_systems.py` and listed in
`selection/report.md`. Data problems found along the way are in
`docs/data_issues.md`.

## 1. Scope and data sources

| Date | Decision | By | Reason |
|---|---|---|---|
| 2026-10-11 | Benchmark condensed-phase binary mixture properties: ΔHmix and mixture density (and V^E), taken from the OpenFF Sage training/validation sets. The starting point was the recommended `dataset-fast-converging.csv` (< 10 ns systems). | User | The gap in ml-peg: there is no ΔHmix benchmark, and no excess properties beyond the ethanol-water density curve. |
| 2026-10-11 | Core set of about 12 systems. | User (recommended option) | Cost against coverage. |
| 2026-10-11 | Elements H, C, N, O and Cl only. | User (recommended option) | Cl is covered by MACE-OFF, AIMNet2 and AceFF; narrower models would otherwise miss systems. |
| 2026-10-11 | Keep 3 compositions per system, even at about 1.7× the liquid_densities cost. | User ("Keep 3 compositions") | Composition dependence is the point of ΔHmix. |
| 2026-10-14 | Choose "maximally informative" systems: both ΔHmix and mixture density at matching compositions, ideally with pure densities, and experimental errors throughout. | User | Each box should score several properties. |
| 2026-10-14 | Avoid overlap with ml-peg: the ethanol-water density curve, water density, and the 63 liquid_densities pure liquids. Pure densities already in ml-peg are simulated (ΔHmix needs them) but not re-scored. | User | "We don't want to duplicate effort". |
| 2026-10-14 | Use the ash-sage-rc2 (Sage 2.3 refit) training and validation sets as the main pool, because they carry experimental uncertainties (the Sage 2.0 set has none). Pin all sources to commit SHAs. | Recommended | Uncertainties are needed for MAZE and for the relative-error filter. |
| 2026-10-14 | Take pure densities and full isotherms from each mixture's own source paper through the NIST ThermoML archive, and cross-check every pure density (ml-peg, Sage 2.0, DIPPR/NIST), flagging differences > 0.002 g/mL. | Recommended | The rc2 pure densities include wrong values (DMF, 1-propanol, piperidine) that flip the sign of V^E. The rc2 sets also keep only 3 points per isotherm. |
| 2026-10-23 | Derive each V^E point from the pure densities of its own paper, and propagate the density σ as an error shared by every density in that paper (the larger of a common scale error and a common offset), plus the paper's RMS residual about the Redlich-Kister fit for the random part. Mixture density σ for scoring is unchanged. | User ("Can we combine the errors into Ve in a way which removes the systematic component"; recommended method) | V^E error bars were up to ~40× the scatter (dipropylamine + DMF: σ 0.13 against 0.003 cm³/mol RMS residual), because a flat ~0.12 % TRC-compiler σ, mostly calibration and purity, was propagated as independent although it cancels in V^E. |
| 2026-10-23 | Methanol + water densities from 10.1016/j.jct.2007.05.004 only (`density_dois`), so its pure methanol (0.78720 g/mL, ~0.0006 above other sources) is used for V^E and for ρ_exp. 10.1021/je700300y (12 points, no pure water) is dropped. Water stays IAPWS-95 (0.997047 g/mL; the paper gives 0.99705). ρ_exp moves from 0.94726 / 0.88130 / 0.84007 to 0.94723 / 0.88181 / 0.84068 g/mL at x = 0.2 / 0.5 / 0.7. | User ("Let's be consistent and use the same pure value throughout, from the paper"; chose the single paper over keeping both) | Pure and mixture densities from one paper share calibration, so one pure value applies throughout. |
| 2026-10-14 | Fix the 2-pyrrolidone tautomer (`OC1=NCCC1` → `O=C1CCCN1`) and guard every InChI round trip against imidic acid tautomers. | Recommended | The Sage data encode the wrong tautomer. |
| 2026-10-23 | Check every curated component before selection runs: zero formal charge, no enol or imidic acid motif, and a match to RDKit's canonical tautomer (unless listed in `TAUTOMER_REVIEWED`). | User (recommended option) | Guards `selection.yml` edits against the 2-pyrrolidone kind of error. A wrong starting species cannot be recovered in MD: acetone's enol is ~10⁻⁸ of the equilibrium mixture ([Chiang et al. 1984](https://doi.org/10.1021/ja00314a055)), so no correct model converts it back in 2 ns. |

## 2. Hard filters

| Date | Decision | By | Reason |
|---|---|---|---|
| 2026-10-11 | Exclude likely two-phase systems ("Flag, verify, likely drop"). | User | In a small MD box a model can stay mixed where the real system separates, or the reverse, so the comparison is not like for like. |
| 2026-10-14 | Two-phase screen: ΔHmix linear in x (H^E = a + b·x, the lever-rule signature), with a tolerance from σ_exp, instead of H^E/x constant. Equimolar ΔHmix must also be ≤ RT/2 (a conservative UCST screen), and each endothermic system gets a literature miscibility check. | User feedback | The line need not pass through the origin. False positives are unlikely, and full isotherms resolve them. Exclusions cite LLE or mutual-solubility data. |
| 2026-10-14 | Avoid systems with a large Sage 3-replicate SD: ≤ 0.10 kJ/mol for ΔHmix and ≤ 0.001 g/mL for density (max over 4 Sage versions). | User | A large repeat spread suggests long equilibration. |
| 2026-10-14 | Minimum signal: \|ΔHmix(x ≈ 0.5)\| ≥ 0.4 kJ/mol. | Recommended | The signal should sit well above the expected simulation noise (0.15 to 0.25 kJ/mol). |
| 2026-10-14 | No viscous liquids (polyols, ethanolamines, polyamines); ≤ 10 heavy atoms per molecule. | Recommended | Slow dynamics and cost. |
| 2026-10-19 | Only run the linearity (lever-rule) test when the largest \|H^E\| ≥ 0.4 kJ/mol. | User ("Yes, make that change") | Prompted by "Would the lever rule selection rule be satisfied by ideal liquids which mix?": ideal and near-ideal mixtures were being flagged as two-phase. This removed none of the proposed systems. |
| 2026-10-19 | Add a relative ΔHmix uncertainty ≤ 15 % filter, and make pure-density uncertainties a hard requirement. | User ("Add the 15 % relative uncertainty rule") | Removes placeholder or very loose experimental errors. |

Result: 248 systems with full data, of which 36 pass all hard filters.

## 3. Aqueous systems

| Date | Decision | By | Reason |
|---|---|---|---|
| 2026-10-14 | Include a few aqueous mixtures, relaxing rules 1 and 5 for them where needed. | User | "Relevant to biomolecular interactions". |
| 2026-10-14 | Drop 2-(dimethylamino)ethanol + water. | User | "Let's avoid 2-(dimethylamino)ethanol due to the large SD" (0.126 kJ/mol). |
| 2026-10-14 | Drop 2-pyrrolidone + water. | User | "Concerned about 2-pyrrolidone due to viscosity". It also has no valid Sage SD (wrong tautomer) and flat 25 % placeholder errors. |
| 2026-10-14 | Exclude ethanol + water, even as ΔHmix only. The user first chose to add it, then reversed. | User | "Let's leave out ethanol + water for now due to the already-present benchmark". |
| 2026-10-14 | Add water + pyridine at 303.15 K, ΔHmix only, with its own pure water and pyridine boxes at 303.15 K. | User (suggested it; accepted after review) | It is the only aromatic N-water system (the histidine/nucleobase hydration motif) and is non-ideal. ThermoML stores 298.15 K, but the paper (10.1021/je2002495) measured at 303.15 K. With no published Sage errors, MAZE uses our own Sage replicates. |
| 2026-10-14 | Keep methanol + water despite a composition mismatch of 0.06 in the rc2 subsets. | Recommended | Full ThermoML isotherms remove the mismatch, and it has the smallest Sage spread in the pool (0.03 kJ/mol), so it calibrates the noise. |
| 2026-10-14 | Exclude acetophenone, acetylacetone and ethyl acetoacetate + water (from the fast-converging list). | User (flag, verify, drop) | Linear H^E (the lever-rule signature), with LLE or solubility data cited in the report. |

## 4. Shaping the set (diversity)

After review on 2026-10-14 the user revised the set to fix: too many alcohols, no non-polar or π control, a narrow magnitude range, little biomolecular coverage, and no V^E diversity check.

| Change | By | Reason |
|---|---|---|
| isooctane + 1-propanol → 1-propanol + benzene | User (from review) | The same alcohol-association physics, plus O-H···π, from a single paper, and avoids the alcohol + alkane miscibility question. |
| methyl pentanoate + methanol → 2-butanone + dibutylamine | User (from review) | Fewer alcohol systems. It also pairs with acetone + N-methylaniline: a ketone with an N-H amine, endothermic in one and exothermic in the other. |
| NMP + 1,4-dichlorobutane → NMP + 1,2-dichloroethane | User (from review) | Widens the exothermic range (-1.26 vs -0.68 kJ/mol), at the cost of one new pure density. |
| Add heptane + toluene | User (from review) | A non-polar dispersion/π control that separates dispersion and D3 errors from hydrogen-bond errors. |
| Tabulate V^E sign and size per system in the report | User (from review) | V^E diversity was not a selection criterion, so the report states its coverage. |
| Reserves: pyrrole + benzene, cyclohexane + benzene, methyl pentanoate + methanol, NMP + 1,4-dichlorobutane, isooctane + 1-propanol, propyl acetate + 1,2-dichloroethane | Recommended | Replacements in case a literature check fails. |

Known gaps, accepted and stated in the report: there is no N-H (secondary) amide system and no amide-water system, and the set reflects ThermoML/Sage coverage (selection bias).

## 5. Candidates considered and not added

| Date | Candidate | Outcome | By |
|---|---|---|---|
| 2026-10-23 | N-methylformamide + o-toluidine, N-methylformamide + tetrahydropyran and N-methylformamide + 2-chlorotoluene at 308.15 K (the former "Candidate additions" report section) | Reviewed; not added yet. Moved to the report's Reserves section with their notes unchanged. | User |
| 2026-10-19 | N-methylformamide + o-toluidine at 308.15 K (would close the N-H amide gap; 5 extra boxes) | Not added: "Do not add NMA / o-toluidine for now". It stays in `selection.yml` under `reserves`. | User |
| 2026-10-14 | Fluorinated systems (2,2,2-trifluoroethanol with THF, water, THP, ethanol or methanol), suggested in a review session because fluorination flips the sign of ΔHmix against the non-fluorinated analogues | Not adopted. The element set stayed H/C/N/O/Cl. No explicit decision or reason is recorded in the logs. | Unrecorded |

## 6. Box sizes

| Date | Decision | By | Reason |
|---|---|---|---|
| 2026-10-14 | Box edge ≥ 25 Å (about 1.5 to 2k atoms) rather than 1000 molecules. Initially ≥ 150 molecules and ≥ 30 of the minority component. | User (recommended option) | Half of a 21.5 Å box is below the ~12 Å receptive field of a 2-layer, 6 Å MLIP. Precision per GPU-hour doesn't depend on box size, and 1000-molecule boxes would cost 10 to 15× more. |
| 2026-10-14 | Issue #936 (which says ~1000 molecules) is left to the user. | User ("No, I'll handle it") | |
| 2026-10-19 | Minority floor lowered to ≥ 20, after checking what Evaluator does (round(x·1000), with ≥ 1 of each). | User | Cost (the 150 total floor had pushed some boxes to 3.7k atoms). |
| 2026-10-19 | Total floor lowered to ≥ 100 molecules. | User ("drop the total floor to 100") | Cost; mixture boxes are now 1.1 to 2.5k atoms. |
| 2026-10-19 | "I'm happy with the selection and system sizes." Selection frozen at 14 systems. | User | |

## 7. Box build and the Sage reference

| Date | Decision | By | Reason |
|---|---|---|---|
| 2026-10-14 | Build pre-equilibrated boxes with Sage (openff-2.2.1) at the experimental density, with scripts in the repo. | User | Starts MLIP runs near equilibrium, which justifies the short 0.2 ns discard. |
| 2026-10-14 | Run Sage itself through the benchmark: 5 replicates on the benchmark boxes, giving σ_Sage for the MAZE denominator. This replaced scaling the published Sage SDs for box size. | User ("Actually, let's just run 5 x replicates with Sage with the proposed boxes") | The noise estimate then matches our box size and run length. |
| 2026-10-20 | Drop piperidine + water (14 → 13 systems, 68 → 64 boxes). | User ("Drop it", recommended option) | Nitrogen inversion is frozen in Sage, even in a 500 K anneal, so ΔHmix depended on the starting N-H orientation by ~0.3 kJ/mol. It also had the largest box-size gap against published Sage and is microheterogeneous for x between 0.08 and 0.5. |
| 2026-10-20 | Build from the lowest-MMFF of 10 ETKDG conformers, then anneal (450 K cooled to T over 1.6 ns) before 500 ps NVT. Rebuild all boxes and rerun all Sage replicates (v2). | User ("Rebuild all, rerun Sage", recommended option) | v1 showed slow conformer equilibria in Sage (THP twist-boat trap, DCE Cl-C-C-Cl, DMAE O-C-C-N) and start dependence. One clean protocol for every box. |
| 2026-10-20 | Keep DMAE + methanol and NMP + 1,2-dichloroethane. | User asked "Should we still keep DMAE and DCE?"; kept on the evidence | The anneal reaches the same conformer populations from different starts (DCE gauche fraction 0.48 vs 0.48; DMAE 0.69 vs 0.73). |
| 2026-10-22 | One estimator for Sage and the MLIPs (ml-peg's `timeseries.py`, applied after 1 ps binning), replacing pymbar. Sage is shown as a plot trace, not a table row. | User (recommended options) | Sage and MLIP uncertainties become directly comparable. |
| 2026-10-22 | Changing the estimator to close the within- vs between-run gap (within-run sem about 20 to 25 % low on average) is left to the user. | User | Not yet decided. |
| 2026-10-23 | Drop the 1 ps binning: Chodera is applied to every logged sample. | User ("Skip the binning for the sake of simplicity") | Over the 320 Sage runs, the raw sem / binned sem has a median of 1.00 (h and density). Only 5 % (h) and 2 % (density) of runs are more than 10 % lower, worst 0.67. The pooled replicate SD / within-run sem changes from 1.24 to 1.26 (h) and 1.21 to 1.22 (density). Not checked at the 10 fs MLIP logging interval. |
| 2026-10-23 | Keep NMP + 1,2-dichloroethane. | User ("keep it") | Its slow conformer sampling is a Sage artefact. In relaxed gas-phase scans the Sage 2.2.1 Cl-C-C-Cl trans→gauche barrier is ~6 kJ/mol above the MACE-OFF (M) reference (24.4 vs 18.2 kJ/mol), while the gauche/trans gap agrees. Sage flips are therefore ~12× too slow (relaxation ~0.9 ns, measured from 100 ps frames), and a reference-like model should relax in ~70 ps, inside the 0.2 ns equilibration. The scan is shown in `sage_report.md`. |
| 2026-10-23 | NMP + DCE ΔHmix from 10.1016/j.jct.2008.01.024 only (`hmix_dois`), at its points nearest x = 0.2, 0.5, 0.8 (`target_x`; boxes x = 0.200, 0.500, 0.804, replacing 0.254 and 0.740). New boxes built and Sage rerun; the old ones are in `data/superseded/nmp_dce_x254_x740/`. | User | The two papers disagree away from equimolar (10.1021/je020065c is 0.15 kJ/mol less exothermic at x ≈ 0.2); j.jct.2008.01.024 matches the fit better and looks more reliable. |
| 2026-10-23 | σ_Sage of NMP + DCE set to the median floor for ΔHmix and density (`SIGMA_SAGE_AT_FLOOR`). | User | Its replicate spread comes from Sage's slow torsion and would loosen the MAZE denominator for a model that does not share it. |
| 2026-10-23 | Significance threshold `FLAG_P` raised from 0.01 to 0.05 (Holm-adjusted for per-box tests, Bonferroni over the 4 pooled tests). The hard-coded 0.01 in the within/between counts, plot and inspection list now follow `FLAG_P`. | User ("make that change", recommended option) | 0.05 family-wise is the conventional level and still conservative after Holm. NMP + DCE x=0.74 (h) becomes Holm-significant; the pooled conclusions are unchanged. |
| 2026-10-23 | Keep the MLP protocol at 0.2 + 2 ns; report 0.2 + 1 ns as a cheaper option to decide on once MLP pilot runs exist. Equilibration stays at 0.2 ns. | User ("Keep 2 ns, note 1 ns option", recommended option) | Power analysis (`data/power_analysis/report.md`), with Sage as a proxy for MLP noise: at 2 ns two models are resolved when their ΔHmix MAE differs by 0.04 to 0.075 kJ/mol (density MAE 0.0005 g/mL); 1 ns costs 55 % and resolves 0.054 to 0.096 kJ/mol. Per-box run lengths would save only 17 %. The 100 to 200 ps window is still shifted for Sage (density p = 2e-6). Numbers from the Sage 2.3.0 runs (re-run 2026-10-29). |
| 2026-10-23 | Flag the V^E RMSE and ΔHmix MAZE 'good' thresholds as below their noise floors and propose new values (V^E 0.15 cm³/mol or a bias-corrected RMSE; MAZE 1.5), without editing ml-peg's `metrics.yml` yet. | User ("Flag + propose new values", recommended option) | A perfect model scores 0.119 cm³/mol on V^E RMSE (threshold 0.05) and 0.65 on ΔHmix MAZE (threshold 1.0) at 2 ns (Sage 2.3.0 runs). |
| 2026-10-28 | Use Sage openff-2.3.0 throughout: the boxes are rebuilt with it, the reference runs (σ_Sage, `systems.csv`) use it, and the finite-size check compares against the published `fb-fit-v3-single-mean-k100` (2.3.0rc2). All 2.2.1 runs and the 2.2.1 vs 2.3.0 comparison are removed. Supersedes the 2.2.1 parts of the 2026-10-14 box decision and the DCE barrier numbers of the 2026-10-23 NMP + DCE row (now in `sage_report.md`). System selection is unchanged and still uses published results from several Sage versions. Replicates 1 to 4 come from an earlier 2.3.0 run (re-packed and equilibrated with 2.3.0, independent of the distributed boxes); replicate 0 was rerun from the new boxes. | User ("consistently use Sage 2.3.0 throughout ... remove anything using Sage 2.2.1"; selection left as is, recommended option) | Comparing against 2.2.1 or 2.3.0rc1 mixed in different parameters (2.3.0 liquids are 0.01 to 0.03 g/mL denser than 2.2.1). The published fb-fit-v3 parameters are identical to the released openff-2.3.0 (checked handler by handler); only the NAGL model differs (0.1.0-rc.3 vs our 1.0.0). |

## 8. Open questions

- **Liquids with nanosecond conformer dynamics (raised 2026-10-22).** The user asked whether DCE and DMAE should stay at all. Starting MLIP runs from Sage-equilibrated conformer populations could hide a model that favours the wrong conformer, if that model also interconverts slowly.
  - Evidence: NMP + DCE is the most-flagged system in the v2 Sage within- vs between-run check (replicate SD 2.0 to 2.5× the within-run sem at x = 0.74), and pure DCE replicate 4 shows an early density shift.
  - Options under discussion:
    - keep both and add a conformer-population diagnostic (dihedral populations and transition counts from the saved frames), flagging unconverged boxes;
    - or drop NMP + DCE and DMAE + methanol.
  - THP + ethanol is not affected: its equilibrium is essentially all chair.
  - **Resolved for DCE (2026-10-23):** NMP + DCE is kept (section 7). DMAE + methanol is still open.
- **Fluorinated systems** (section 5): not decided either way.
