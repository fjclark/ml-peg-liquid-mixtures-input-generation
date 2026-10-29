# Data issues found while selecting the binary mixture systems

Issues found during the selection of the binary mixtures benchmark (October 2026).
Each entry says how it was found and whether it was verified. "Verified" means checked
against the primary source or several independent sources; "suspected" means a
consistent pattern that has not been confirmed against the paper.

Densities are in g/mL, enthalpies in kJ/mol, temperatures in K. "Sage sets" means the
OpenFF Sage 2.0 training set (openff-sage@a6a061b) and the Sage 2.3 refit training and
validation sets (ash-sage-rc2@00c9524).

## 1. Existing ml-peg data

### 1.1 liquid_densities: reference temperatures look about 5 K too high (suspected, strong pattern)

The liquid_densities inputs store `exp_temperature` per liquid, and the calc runs at that
temperature. Two separate patterns suggest that many references were measured about 5 K
below their label. Liquid densities fall by roughly 0.001 g/mL per K, so a 5 K offset
shifts the reference by about 0.005 g/mL.

**(a) Liquids labelled 298 K look like 20 °C values.** 21 liquids labelled 298 K also
have a 298.15 K pure density in the Sage sets. 17 of them differ by more than
0.002 g/mL, and every difference has the same sign (ml-peg higher):

| Liquid | ml-peg (labelled 298 K) | Sage sets, 298.15 K | Difference |
|---|---|---|---|
| tert-butanol | 0.7887 | 0.7807 | +0.0080 |
| methyl acetate | 0.9342 | 0.9271 | +0.0071 |
| 1,4-dioxane | 1.0337 | 1.0279 | +0.0058 |
| methanol | 0.7914 | 0.7863 | +0.0051 |
| chlorobenzene | 1.1058 | 1.1009 | +0.0049 |
| acetophenone | 1.0281 | 1.0234 | +0.0047 |
| anisole | 0.9940 | 0.9893 | +0.0047 |
| acetonitrile | 0.7857 | 0.7811 | +0.0046 |
| diglyme | 0.9434 | 0.9388 | +0.0046 |
| morpholine | 1.0005 | 0.9960 | +0.0045 |
| toluene | 0.8668 | 0.8623 | +0.0045 |
| triethylamine | 0.7275 | 0.7231 | +0.0044 |
| 1-butanol | 0.8095 | 0.8055 | +0.0040 |
| ethylene glycol | 1.1135 | 1.1096 | +0.0039 |
| ethanol | 0.7893 | 0.7864 | +0.0029 |
| glycerol | 1.2613 | 1.2585 | +0.0028 |
| benzene | 0.8765 | 0.8738 | +0.0027 |

- The four within 0.002 g/mL are 1,3-dioxolane, ethyl acetate, hexane and pyridine.
  Pyridine's two Sage values themselves disagree (see 2.3), and its ml-peg value
  (0.9819) is the 20 °C handbook value.
- **Water:** ml-peg uses 0.998 at 298 K. That is the 20 °C value (0.99820);
  [IAPWS-95](https://doi.org/10.1063/1.1461829) gives 0.99705 at 298.15 K.
- **Known 20 °C handbook values:** methanol 0.7914, ethanol 0.7893, benzene 0.8765 and
  toluene 0.8668 match common 20 °C handbook values exactly.
- **Caveats:**
  - Some Sage references are single points.
  - Acetonitrile's Sage value comes from a paper whose pure densities look high
    (see 2.3), so its true offset is larger.

**(b) Liquids labelled 303 K look like 25 °C values.** For six liquids labelled 303 K
we have several 298.15 K ThermoML measurements (fetched for this benchmark). Five
match the 298.15 K data almost exactly; at a true 303 K they should be about
0.005 g/mL lower:

| Liquid | ml-peg (labelled 303 K) | ThermoML, 298.15 K | Verdict |
|---|---|---|---|
| N,N-dimethylformamide | 0.9445 | 0.9438–0.9441 | 25 °C value |
| 1-propanol | 0.7997 | 0.7997 | 25 °C value |
| 2-butanone | 0.7999 | 0.7997–0.8001 | 25 °C value |
| 1,2-dichloroethane | 1.2454 | 1.2455–1.2474 | 25 °C value |
| acetone | 0.7845 | 0.7843–0.7853 | 25 °C value |
| N-methyl-2-pyrrolidone | 1.0230 | 1.0280–1.0286 | consistent with 303 K |

- **Impact:** an offset of 0.003 to 0.008 g/mL is comparable to the error of a good
  MLIP, so it biases the density MAE and can reorder models.
- **Suggested fix:** check the temperature of every reference against its source
  (the docs cite Perry's Handbook), then relabel `exp_temperature` (293 K or 298 K)
  or replace the values with data at the labelled temperature.
- **Scope:** the other liquids (other temperatures, or no Sage/ThermoML value at hand)
  were not checked.

### 1.2 liquid_densities docs count (verified)

`docs/source/user_guide/benchmarks/molecular_dynamics.rst` says "61 organic liquids",
but the S3 archive `liquid_densities.zip` holds 63 structures.

### 1.3 Cross-check tooling

`src/liquid_mixtures_inputs/mlpeg_liquids.yml` maps all 63 liquid_densities structure names to SMILES
(checked against the box formulas). It can be reused for a full check of section 1.1
against ThermoML.

## 2. Issues affecting the OpenFF (Sage) fits and benchmarks

### 2.1 2-pyrrolidone stored as the wrong tautomer (verified)

- **Error:** 2-pyrrolidone appears as `OC1=NCCC1` (the imidic acid) instead of
  `O=C1CCCN1` in:
  - the Sage 2.0 training set (31 rows, e.g. [line 145](https://github.com/openforcefield/openff-sage/blob/a6a061ba8ebc144e3fda46f5f7a9e7a352b9cd2f/data-set-curation/physical-property/optimizations/data-sets/sage-train-v1.csv?plain=1#L145));
  - the rc2 training set (30 rows, e.g. [line 769](https://github.com/openforcefield/ash-sage-rc2/blob/00c9524d38af53725a1a18feb3433f5c191cdb4a/01_download-data/physprop/final/output/training-set.csv?plain=1#L769));
  - the rc2 validation set (1 row, [line 45](https://github.com/openforcefield/ash-sage-rc2/blob/00c9524d38af53725a1a18feb3433f5c191cdb4a/01_download-data/physprop/final/output/validation-set.csv?plain=1#L45));
  - the rc2 benchmark outputs (`04_benchmark`, e.g. [line 2843](https://github.com/openforcefield/ash-sage-rc2/blob/00c9524d38af53725a1a18feb3433f5c191cdb4a/04_benchmark/phys-prop/output/training/summary-benchmarks-EnthalpyOfMixing.csv?plain=1#L2843)).
- **Consequence:** the force field was fitted to, and benchmarked on, a different
  molecule, so the published Sage errors and replicate SDs for 2-pyrrolidone systems are
  not meaningful.
- **Scope:** a scan of every SMILES in the three data sets for the `[OH]-C=N` pattern
  found no other cases.
- **Likely cause:** converting InChI to a molecule. Standard InChI does not fix the
  mobile hydrogen, and RDKit's `MolFromInchi` returns the imidic acid for N-H amides
  (verified for 2-pyrrolidone, N-methylformamide and formamide).
- **Not affected:** `dataset-fast-converging.csv` has the correct SMILES.

### 2.2 Pyridine + water ΔHmix temperature mislabelled (verified)

- **Error:** the [ThermoML record](https://trc.nist.gov/ThermoML/10.1021/je2002495.json)
  for 10.1021/je2002495 stores its ΔHmix at 298.15 K, but the
  [paper's](https://doi.org/10.1021/je2002495) title ("... Enthalpy of Mixing of
  Pyridine + Water at 303.15 K") and abstract say 303.15 K.
- **Where it propagates:** the Sage curation, including `dataset-fast-converging.csv`
  (three points at 298.15 K, [lines 14-16](../selection/inputs/dataset-fast-converging.csv)).
- **Consequence:** any fit or benchmark that used these points compared a 298.15 K
  simulation with 303.15 K data.
- **Suggested fix:** report it to NIST TRC.

### 2.3 Pure densities that disagree with other sources (verified against two or more independent papers)

Values used in the Sage sets that differ from the consensus by more than 0.002 g/mL:

| Liquid | Value in Sage sets | Source | Consensus 298.15 K | Note |
|---|---|---|---|---|
| DMF | 0.9501 | 10.1007/s10765-016-2096-3 | 0.9438–0.9441 (three papers) | +0.006 |
| N,N-dimethylacetamide | 0.9502 | same paper | ≈ 0.937 (literature) | Looks wrong; not cross-checked in our pool |
| acetonitrile | 0.7811 | same paper | ≈ 0.7766 (literature) | Same paper again; suspected |
| 1-propanol | 0.8032 | 10.1021/je900212d | 0.7997 | Equals the 20 °C value |
| piperidine | 0.8646 | 10.1016/j.fluid.2009.10.018 | 0.8565–0.8575 | +0.007 |
| NMP | 1.0258 | 10.1016/j.fluid.2008.02.023 | 1.0280–1.0286 (three papers) | −0.0025 |
| dibutylamine | 0.7578 | 10.1021/je800026p | 0.7555 (two papers) | +0.0023 |
| pyridine | 0.9829 | 10.1016/j.fluid.2015.04.031 | 0.9780–0.9783 (three papers) | +0.005 |

- **Pattern:** all three pure densities from 10.1007/s10765-016-2096-3 are high, which
  suggests a temperature or calibration problem in that paper rather than isolated
  typos.
- **Impact:** wrong pure densities shift the V^E derived from mixture densities. With
  the rc2 value for DMF, for example, the sign of V^E for dipropylamine + DMF flips.
- **Data checked:** the DMAc and acetonitrile rows were not checked against
  primary sources here; this benchmark uses neither.

### 2.4 Implausibly large pure-density uncertainties in the rc2 validation set (verified)

These do not affect the value, but they matter for any uncertainty-weighted objective.
Typical measured uncertainties are 0.0001 to 0.0005 g/mL, yet several papers carry
values 20 to 100 times larger:

| σ (g/mL) | Paper | Liquids |
|---|---|---|
| 0.02 | 10.1016/j.fluid.2013.11.008 | methanol, ethanol, 1-butanol, 2-propanol |
| 0.01 | 10.1016/j.jct.2012.07.006 | cyclohexane, heptane, octane, N-methylaniline, toluene (0.0122) |
| 0.01 | 10.1021/acs.jced.5b00542 | benzene, hexane |
| 0.008 | 10.1021/acs.jced.6b00019 | tert-butanol, methylcyclohexane, isooctane, p-xylene |
| 0.006 | 10.1007/s10765-016-2096-3 | DMF, DMAc, acetonitrile |

Mixture densities from 10.1016/j.jct.2005.10.011 (heptane + toluene) also carry
σ ≈ 0.012 g/mL.

### 2.5 Placeholder ΔHmix uncertainties (verified)

- **Placeholder values:** 10.1016/j.jct.2013.09.015 (2-pyrrolidone + water,
  2-pyrrolidone + pyridine, NMP + pyridine) has ΔHmix uncertainties of exactly 25 % of
  the value in the rc2 training set (e.g. σ = 0.1695 for ΔHmix = -0.678,
  [line 769](https://github.com/openforcefield/ash-sage-rc2/blob/00c9524d38af53725a1a18feb3433f5c191cdb4a/01_download-data/physprop/final/output/training-set.csv?plain=1#L769)). These are curation placeholders, not measured errors. (The
  Sage 2.0 training set has no ΔHmix uncertainties at all, e.g.
  [line 145](https://github.com/openforcefield/openff-sage/blob/a6a061ba8ebc144e3fda46f5f7a9e7a352b9cd2f/data-set-curation/physical-property/optimizations/data-sets/sage-train-v1.csv?plain=1#L145).)
- **Large but not placeholders:** several other systems have relative uncertainties
  above 15 % because their H^E is close to zero rather than because of placeholders:
  - ester + dibromoalkane systems (10.1016/j.jct.2005.07.009, 10.1016/j.jct.2006.05.004);
  - 10.1016/j.jct.2005.03.020;
  - 10.1016/j.jct.2005.11.003;
  - 10.1016/j.fluid.2007.03.017.

### 2.6 Likely two-phase systems in the training data

These systems are in the rc2 training set and in `dataset-fast-converging.csv`. All
three are aqueous mixtures with ΔHmix that is linear in x over x = 0.2 to 0.8, the
lever-rule signature of a two-phase region. A simulation of a mixed box compared with
data from a two-phase sample is not a like-for-like comparison, so these points are
poor fitting targets.

| System | Slopes between the three points | Evidence | Status |
|---|---|---|---|
| acetophenone + water | 0.677, 0.671 | Aqueous solubility 6.13 g/L at 25 °C (x ≈ 0.001; Southworth and Keller 1986, via [PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/7410#section=Solubility)) | Two-phase verified by linearity; organic-rich solubility limit not found |
| acetylacetone + water | 0.840, 0.867 | ThermoML mutual solubilities ([10.1016/j.fluid.2010.11.029](https://trc.nist.gov/ThermoML/10.1016/j.fluid.2010.11.029.json)): aqueous phase x ≈ 0.033, organic phase x ≈ 0.75–0.85 at 288–313 K | Verified |
| ethyl acetoacetate + water | 1.710, 1.469 | Aqueous solubility 2.86 g/100 mL at 20 °C (x ≈ 0.004; ICSC via [PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/8868#section=Solubility)); PubChem also lists 110 g/L at 17 °C (AQUASOL); no primary source checked | Suspected; x = 0.75 may be single-phase |

### 2.7 Disagreement between ΔHmix papers (verified)

- **NMP + 1,2-dichloroethane:** 10.1021/je020065c lies about 0.07 kJ/mol above
  10.1016/j.jct.2008.01.024. Measured against a joint Redlich-Kister fit:
  - je020065c: mean residual +0.066, RMS 0.087;
  - jct.2008.01.024: mean residual −0.013, RMS 0.019;
  - reported σ is 0.005 to 0.015, so the gap is several times larger than either
    paper's uncertainty.
- **Heptane + toluene** (10.1021/je400184t and 10.1016/j.jct.2008.12.004) agree within
  their uncertainties.

### 2.8 Slow conformer equilibration from single-conformer boxes (suspected for the fits)

Evaluator builds its boxes the same way the first version of this benchmark did:
`openff/evaluator/utils/packmol.py` (`_create_trajectory`) calls
`generate_conformers(n_conformers=1)`, and packmol copies that one conformer into
every molecule. Our first Sage runs (3.3, with openff-2.2.1) found
conformational equilibria that take nanoseconds: tetrahydropyran twist-boat to chair,
1,2-dichloroethane Cl-C-C-Cl, DMAE O-C-C-N and piperidine nitrogen inversion
(frozen even at 500 K). Properties of mixtures containing these molecules could
therefore depend on the starting conformer, especially with the Sage 2.0
protocol (200 ps equilibration, then 2 ns). Not checked:
- which conformer the OpenEye backend (used for the fits) returns; Omega
  usually gives a low-energy one, which avoids the tetrahydropyran twist-boat
  but still starts every molecule in the same torsional state;
- how far the rc2 pre-equilibration (until 300 uncorrelated samples) removes
  the bias.

## 3. Issues affecting this PR (binary_mixtures)

### 3.1 Handled in the selection code

| Issue | Handling |
|---|---|
| InChI → SMILES gives imidic-acid tautomers (2.1) | `thermoml.mol_from_inchi` converts them back to amides; covered by a unit test. `select_systems.TAUTOMER_FIXES` repairs the Sage CSVs, and Sage benchmark rows for the wrong tautomer are ignored. |
| Pyridine + water at 303.15 K, mislabelled 298.15 K (2.2) | `temperature_corrections` in `selection.yml`. Its mixture and pure boxes (water, pyridine) run at 303.15 K, with a 303.15 K pyridine reference density (10.1021/je900874z) and IAPWS-95 water. |
| Wrong pure densities in the Sage sets (2.3) | Pure densities come from the same paper as the mixture densities, and are cross-checked against every other source. The report lists sources that disagree by more than 0.002 g/mL. |
| ml-peg reference temperatures (1.1) | ml-peg pure densities are never used as references; water comes from IAPWS-95. |
| Noisy mixture densities (2.4) | Points with σ > 0.002 g/mL are left out of the V^E fits (heptane + toluene: 10 points from 10.1016/j.jct.2005.10.011). |
| Placeholder or near-zero ΔHmix errors (2.5) | Hard filter: relative ΔHmix uncertainty ≤ 15 %. |
| Two-phase systems (2.6) | Linear-H^E filter (applied when the largest \|H^E\| ≥ 0.4 kJ/mol) and the RT/2 screen; the three aqueous systems are excluded with the evidence above. |
| ThermoML expanded vs standard uncertainties | Expanded uncertainties are divided by the coverage factor 2. This matches the rc2 convention (checked). |

### 3.2 Open items and caveats for the selected systems

- **NMP + 1,2-dichloroethane (2.7):**
  - **Problem:** the reference values mix the two papers that disagree. The x = 0.25
    and 0.74 targets come from 10.1021/je020065c and the x = 0.5 target from
    10.1016/j.jct.2008.01.024, so the references carry an inter-paper offset of about
    0.07 kJ/mol.
  - **Option:** take all three targets from 10.1016/j.jct.2008.01.024, which has the
    smaller uncertainties.
  - **Size:** the offset is below the expected MAZE σ_ref (~0.1–0.2 kJ/mol), but it is
    not negligible.
- **Pyridine + water:** has no mixture density data, so it is scored on ΔHmix only, and
  it has no published Sage errors (σ_Sage comes from our own replicates).
- **Propionitrile + tert-amyl ethyl ether:** ΔHmix and densities are at 298.4 K
  (benchmark runs at 298.15 K). All data come from one paper.
- **Unverified literature values used only to explain exclusions:**
  - the ethyl acetoacetate aqueous solubility;
  - the organic-rich solubility limit for acetophenone + water.

  Neither affects the selected systems.
- **Composition matching:** NMP + 1,2-dichloroethane fails the 3-point composition
  match by 0.0003 (|Δx| = 0.0503). The full ThermoML isotherms remove the mismatch;
  recorded as a filter exception.

### 3.3 Slow conformer equilibration in the Sage reference run (fixed)

The first Sage reference run (v1, 2026-10-19/20, 68 boxes × 5 replicates) was
built from a single ETKDG conformer per molecule, copied into every molecule by
packmol, followed by 500 ps NVT. Checks on that run:

- **Tetrahydropyran:** the conformer was a twist-boat (27 kJ/mol above the chair
  in MMFF). About 20 % of rings were still non-chair in the distributed boxes and
  flipped during the first ~1 ns of NPT, biasing ΔHmix by up to +0.7 kJ/mol. Fixed
  by taking the lowest-MMFF-energy of 10 conformers.
- **Start dependence (replicate 0, the distributed box, against replicates 1-4,
  re-packed from the lowest-energy conformer):**
  - piperidine: every N-H stayed in its starting orientation (axial in replicate
    0, equatorial in replicates 1-4) for the whole 2.2 ns. Density offsets of up
    to 0.008 g/mL; pure h offset 0.69 kJ/mol; ΔHmix offsets ~0.3 kJ/mol. A 500 K
    anneal did not equilibrate it (0.87 vs 0.37 axial from opposite starts).
  - DMAE: density +0.0055 g/mL in replicate 0, relaxing over ~1.4 ns.
  - 1,2-dichloroethane: gauche fraction relaxes over > 2 ns from both sides
    (density 1.220 down and 1.200 up towards ~1.208 g/mL).
- **Fix:** an NVT anneal (450 K cooled to the target in 16 steps over 1.6 ns)
  before the 500 ps NVT. From opposite starting conformers it gives matching
  populations for DCE (gauche 0.48 vs 0.48; 0.65 vs 0.21 without). DMAE was only
  tested with a shorter anneal (500 K over 550 ps): 0.69 vs 0.73 (0.41 vs 0.73
  without). The v2 start-dependence screen checks every box.
- **Piperidine + water was dropped** (see `selection.yml`). It also had the largest
  box-size gap against the published 1000-molecule Sage run (ΔHmix −1.09 vs −0.74
  kJ/mol at x = 0.249, z = −4.8, [published value](https://github.com/openforcefield/ash-sage-rc2/blob/00c9524d38af53725a1a18feb3433f5c191cdb4a/04_benchmark/phys-prop/output/training/summary-benchmarks-EnthalpyOfMixing.csv?plain=1#L4806); the median gap
  over 35 points was 0.07 kJ/mol) and is microheterogeneous for x between 0.08 and
  0.5 ([arXiv:2203.01224](https://arxiv.org/abs/2203.01224)).
- All boxes were rebuilt with the anneal and the Sage reference run repeated (v2).
  `summarise_sage.py` now reports the replicate-0 start-dependence screen.
