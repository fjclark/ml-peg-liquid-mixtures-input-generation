# liquid_mixtures_inputs

Inputs for the ml-peg `binary_mixtures` benchmark, which tests how well machine-learned potentials (MLPs) predict the enthalpy of mixing (ΔHmix), mixture density and excess molar volume (V^E) of binary liquid mixtures. We screened 370 binary systems in the OpenFF Sage training and validation sets, kept the 36 that pass a set of hard filters, and chose 10 diverse mixtures at 3 compositions each (50 boxes including the pure liquids). Every box was built and equilibrated with Sage 2.3.0, and Sage was run 5 times per box and analysed with the benchmark protocol to measure simulation noise. That noise sets the MAZE denominators and ml-peg's 'good' thresholds. With one 0.2 + 1 ns run per box, two models can be told apart when their ΔHmix MAE differs by 0.05 to 0.1 kJ/mol, at a cost of about 140 MACE-OFF23(M) A100 GPU-days per model.

ml-peg does not need this repository to run the benchmark: it reads the pre-built boxes and reference data from S3. This repository records how they were made and regenerates them.

## Key decisions

Every decision, with its reason and history, is in [DECISIONS.md](DECISIONS.md).

- **Data.** ΔHmix and densities from the Sage training and validation sets, with full isotherms and pure densities taken from each mixture's own source paper through the NIST ThermoML archive. All sources are pinned. ([scope and data](DECISIONS.md#scope-and-data))
- **Filters.** H, C, N, O and Cl only; ≤ 10 heavy atoms; no viscous or likely two-phase systems; |ΔHmix(x ≈ 0.5)| ≥ 0.4 kJ/mol; a small published Sage replicate spread; experimental uncertainties throughout. 36 of 248 systems with enough data pass. ([filters](DECISIONS.md#filters))
- **Systems.** 10 mixtures chosen to cover association, π interactions, amines, two aqueous systems and a non-polar control, without overlapping ml-peg's existing density benchmarks. ([systems](DECISIONS.md#systems))
- **Boxes.** Edge ≥ 25 Å and ≥ 100 molecules (952 to 2,400 atoms), built from low-energy conformers and annealed with Sage 2.3.0, so runs start near equilibrium. ([boxes and protocol](DECISIONS.md#boxes-and-protocol))
- **Protocol.** One run per model and box: 0.2 ns NPT equilibration and 1 ns production, chosen for cost (2 ns would cost 1.8× as much for a modest gain in precision). ([boxes and protocol](DECISIONS.md#boxes-and-protocol))
- **Metrics.** Every property is scored by MAZE, dividing each error by σ_ref = sqrt(σ_Sage² + σ_exp²), with the MAE shown at weight 0. 'Good' is the score of a model that matches experiment exactly plus 2 SD of its noise; 'bad' is 5× Sage's score. ([metrics and thresholds](DECISIONS.md#metrics-and-thresholds))

## Outcomes

- **Resolution.** With one 0.2 + 1 ns run per box, two models are resolved (5 % test, 80 % power) when their ΔHmix MAE differs by 0.05 to 0.094 kJ/mol or their mixture density MAE by 0.0004 to 0.0007 g/mL. A single point detects an error above about 0.28 kJ/mol (ΔHmix) or 0.0027 g/mL (density). ([power analysis](data/power_analysis/report.md))
- **Cost.** 76,065 atoms over 50 boxes, about 138 MACE-OFF23(M) A100 GPU-days per model (order of magnitude). 0.2 + 2 ns would cost 1.8× as much (254 GPU-days) and narrow the resolvable ΔHmix MAE difference for a model with Sage-sized errors from 0.094 to 0.073 kJ/mol. ([details](data/power_analysis/details.md#4-run-length-and-cost))
- **Sage 2.3.0 against experiment.** ΔHmix MAE 0.59 kJ/mol (MAZE 5.0), mixture density MAE 0.0094 g/mL (MAZE 10.1), V^E MAE 0.24 cm³/mol (MAZE 1.9). Single 1 ns Sage runs understate their own standard errors by 17 % (h) and 6 % (density). ([Sage report](data/sage_summary/report.md))
- **Thresholds.** For example ΔHmix MAZE good 0.82, bad 25; all eight metrics are in [thresholds.csv](data/power_analysis/thresholds.csv), with a ready-to-paste [metrics_snippet.yml](data/power_analysis/metrics_snippet.yml).
- **Main caveat.** Sage stands in for MLP noise, which may be very different.

## Reports

| Report | What it answers |
|---|---|
| [selection/report.md](selection/report.md) | Which systems were chosen and why: filters, screening funnel, box plan, per-system data and sources. |
| [data/sage_summary/report.md](data/sage_summary/report.md) | How Sage does against experiment, how noisy single runs are, which boxes to inspect, and whether box size matters. Details in `details.md` and `inspection.csv`. |
| [data/power_analysis/report.md](data/power_analysis/report.md) | What one run per box can resolve, what it costs, and the thresholds. Details in `details.md`. |
| [DECISIONS.md](DECISIONS.md) | Every decision, its reason and evidence, and a dated history. |
| [docs/data_issues.md](docs/data_issues.md) | Data problems found along the way, in ml-peg's existing data and in the Sage training data. |

## Reproduce

Requires [pixi](https://pixi.sh) and, for `build-boxes` and `run-sage`, a CUDA GPU. There are two environments, so the conda-only OpenFF stack never has to solve against ml-peg's PyPI dependencies:

| Environment | Used for |
|---|---|
| `openff` | `build-boxes`, `run-sage`, `dce-scan` (OpenFF toolkit, NAGL, OpenMM, packmol) |
| `analysis` | `select`, `summarise-sage`, `power-analysis`, tests (ml-peg, RDKit) |

```bash
# 1. Select the systems (writes selection/)
pixi run -e analysis select --output-dir selection \
    --fast-set selection/inputs/dataset-fast-converging.csv

# 2. Build the boxes and run the Sage reference (GPU)
pixi run -e openff build-boxes --selection-dir selection \
    --output-dir data/binary_mixtures
pixi run -e openff run-sage --data-dir data/binary_mixtures \
    --output-dir data/sage_runs

# 3. Summarise the Sage runs (adds the Sage columns to systems.csv, pure.csv
#    and boxes.csv in place, once every replicate is complete)
pixi run -e analysis summarise-sage --data-dir data/binary_mixtures \
    --runs-dir data/sage_runs --output-dir data/sage_summary

# 4. Power analysis and thresholds (caches the re-analysis in truncation.csv)
pixi run -e analysis power-analysis --data-dir data/binary_mixtures \
    --runs-dir data/sage_runs --summary-dir data/sage_summary \
    --output-dir data/power_analysis
```

The ThermoML downloads are cached in `~/.cache/ml_peg/binary_mixtures_selection`. Both run scripts skip finished work, so they can be restarted. The DCE torsion scans behind the NMP + 1,2-dichloroethane decision are run with `dce-scan` (see `dce_torsion_scan.py`).

To publish the inputs, zip `data/binary_mixtures` and upload it to the S3 key `inputs/molecular_dynamics/binary_mixtures/binary_mixtures.zip`, with `binary_mixtures/` as the top-level directory.

```bash
cd data && zip -r binary_mixtures.zip binary_mixtures
```

Tests and linting:

```bash
pixi run -e analysis test
pixi run -e analysis pre-commit run --all-files
```

## Repository layout

`src/liquid_mixtures_inputs/`:

- `selection.yml`: the curated systems, with rationale, filter exceptions, pure-density references, dropped systems and reserves.
- `mlpeg_liquids.yml`: SMILES of ml-peg's liquid_densities liquids, so their pure densities are not scored twice.
- `thermoml.py`: fetches and parses NIST TRC ThermoML records by DOI.
- `select_systems.py`: screens the Sage data sets, checks the curated systems, fits the ThermoML isotherms and plans box sizes. `selection_report.py` writes the selection report.
- `build_boxes.py`: packs, parametrises (openff-2.3.0, NAGL charges, TIP3P water), anneals and equilibrates every box at the experimental density.
- `run_sage_benchmark.py`: runs every box with Sage for 2.2 ns (5 replicates). The summaries use the protocol's first 0.2 + 1 ns; the rest lets the power analysis compare 2 ns.
- `summarise_sage.py`: ΔHmix, densities and V^E per replicate, σ_Sage and σ_ref, convergence, within- vs between-run, finite-size and start-dependence checks. `sage_report.py` writes the Sage reports.
- `dce_torsion_scan.py`: Cl-C-C-Cl scans of 1,2-dichloroethane with Sage and MACE-OFF23(M), and trans/gauche flip rates in the Sage runs. The MACE step runs in the ml-peg virtual environment, as MACE is in neither pixi environment.
- `power_analysis.py`: re-analyses the Sage runs at shorter lengths to give detection limits, resolvable differences, cost and thresholds. `power_report.py` writes the power reports.
- `constants.py`: physical constants and the benchmark protocol.

Outputs: `selection/` (selection tables, report and figures), `data/sage_summary/` and `data/power_analysis/` (tables, reports and figures). `data/binary_mixtures/` and `data/sage_runs/` are not committed.

## Data provenance

- OpenFF Sage 2.0 training set: `openforcefield/openff-sage` (MIT), pinned commit in `select_systems.py`.
- Sage 2.3 refit training and validation sets and published Sage benchmark results: `openforcefield/ash-sage-rc2`, pinned commit.
- Experimental data: NIST TRC ThermoML archive, cited by DOI per point.

## AI use

Claude Opus 5.5 was used to help plan, draft and iteratively improve this code, with human verification and approval of all code and decisions.
