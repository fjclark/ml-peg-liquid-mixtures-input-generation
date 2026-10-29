# liquid_mixtures_inputs

System selection, Sage box building and the Sage reference runs for the ml-peg
`binary_mixtures` benchmark (enthalpy of mixing, mixture density and excess
volume of 13 binary liquid mixtures).

ml-peg does not need this repository to run the benchmark: it reads the
pre-built boxes and the reference data from S3. This repository records how
the systems were chosen and lets those inputs be regenerated.

See [`DECISIONS.md`](DECISIONS.md) for a dated record of decisions made
while selecting the systems (filters, exceptions, dropped and replaced
systems, box sizes, the box build), with the reasons and the open questions.

## Requirements

- [pixi](https://pixi.sh).
- A CUDA GPU for `build-boxes` and `run-sage`

There are two environments, so the conda-only OpenFF stack never has to solve
against ml-peg's PyPI dependencies:

| Environment | Used for |
|---|---|
| `openff` | `build-boxes`, `run-sage` (OpenFF toolkit, NAGL, OpenMM, packmol) |
| `analysis` | `select`, `summarise-sage`, tests (ml-peg, RDKit) |

## Contents

- `src/liquid_mixtures_inputs/`:
  - `selection.yml`: the curated systems, with rationale, filter exceptions,
    pure-density references and dropped systems.
  - `mlpeg_liquids.yml`: SMILES for ml-peg's liquid_densities liquids, so their
    pure densities are not scored twice.
  - `thermoml.py`: fetches and parses NIST TRC ThermoML records by DOI.
  - `select_systems.py`: screens the OpenFF Sage physical-property sets, checks the curated systems, pulls full isotherms from ThermoML, cross-checks pure densities and plans box sizes.
  - `report.py`: writes the selection report and figures.
  - `build_boxes.py`: packs, parametrises (openff-2.3.0, NAGL AM1-BCC charges, TIP3P water) and NVT-equilibrates every benchmark box at the experimental density. Each molecule starts from the lowest-MMFF-energy of 10 ETKDG conformers, and the NVT stage includes an anneal (450 K cooled to the target over 1.6 ns) before 500 ps at the target temperature. Packmol copies one conformer into every molecule, and an initial Sage reference run (v1) showed that several conformational equilibria take nanoseconds in Sage: tetrahydropyran twist-boat to chair, the Cl-C-C-Cl torsion of 1,2-dichloroethane, the O-C-C-N torsion of DMAE and nitrogen inversion in piperidine. Without the anneal, results depended on the starting conformer.  Piperidine inversion does not equilibrate even at 500 K, so piperidine + water was dropped (see `selection.yml`).
  - `run_sage_benchmark.py`: runs every box with Sage as a reference model
    (5 independent replicates, the benchmark's 0.2 + 2 ns NPT protocol).
  - `summarise_sage.py`: ΔHmix, densities and V^E per replicate, the pooled replicate SD used as the MAZE denominator, convergence checks, a within- vs between-run comparison, a finite-size comparison with the published 1000-molecule Sage results, and a start-dependence screen (replicate 0, the distributed box, against the re-packed replicates).
  - `dce_torsion_scan.py`: relaxed Cl-C-C-Cl scans of 1,2-dichloroethane with Sage and MACE-OFF23(M), and trans/gauche flip rates in the Sage runs, behind the report's case for keeping NMP + DCE. The MACE step runs in the ml-peg virtual environment, as MACE is in neither pixi environment.
- `DECISIONS.md`: the selection decisions and their reasons.
- `docs/data_issues.md`: potential data problems (to be verified) found during selection (ml-peg references, OpenFF training data, and this benchmark).
- `selection/`: the committed selection outputs (`candidates.csv`, `systems.csv`, `pure.csv`, `isotherms.csv`, `report.md`, `figures/`), and `selection/inputs/dataset-fast-converging.csv`, the `--fast-set` input.
- `data/sage_summary/`: the committed Sage reference summary.

## Usage

### 1. Select the systems

```bash
pixi run -e analysis select --output-dir selection
```

Add `--fast-set selection/inputs/dataset-fast-converging.csv` to also assess
the recommended fast-converging list in `report.md`.

This writes, in `--output-dir`:

- `candidates.csv`: every binary with ≥ 3 enthalpy and ≥ 3 density points
  near 298 K, with all screening metrics and failures.
- `systems.csv`: one row per mixture box (molecule counts, box composition,
  experimental ΔHmix, density and V^E at that composition, uncertainties and
  DOIs).
- `pure.csv`: chosen pure densities with all cross-check values.
- `isotherms.csv`: the ThermoML points behind each fit.
- `report.md` and `figures/`.

### 2. Build the boxes and run the Sage reference

```bash
pixi run -e openff build-boxes --selection-dir selection \
    --output-dir data/binary_mixtures
pixi run -e openff run-sage --data-dir data/binary_mixtures \
    --output-dir data/sage_runs
pixi run -e analysis summarise-sage --data-dir data/binary_mixtures \
    --runs-dir data/sage_runs --output-dir data/sage_summary
```

`build-boxes` writes `boxes/<box_id>.xyz` (extxyz with the SMILES, mapped
SMILES, molecule counts and `exp_temperature` in `atoms.info`), `boxes.csv`,
and copies `systems.csv` and `pure.csv`. `summarise-sage` adds the Sage
reference columns to `systems.csv` and `boxes.csv` **in place**. It does so
only when every replicate of every system is complete; otherwise it leaves
the existing columns as they are. Both run scripts skip finished work, so they
can be restarted.

### 3. Upload to S3

```bash
cd data && zip -r binary_mixtures_v2.zip binary_mixtures
```

ml-peg expects it at the S3 key
`inputs/molecular_dynamics/binary_mixtures/binary_mixtures_v2.zip`, with
`binary_mixtures/` as the top-level directory. Bump the zip name whenever
columns change, to avoid stale ml-peg cached copies.

## Tests and linting

```bash
pixi run -e analysis test
pixi run -e analysis pre-commit install
```

## Data provenance

- OpenFF Sage 2.0 training set: `openforcefield/openff-sage` (MIT), pinned commit in `select_systems.py`.
- Sage 2.3 refit training/validation sets and published Sage benchmark
  results: `openforcefield/ash-sage-rc2`, pinned commit.
- Experimental data: NIST TRC ThermoML archive, cited by DOI per point.

## AI Use

Claude Opus 5.5 was used to help plan, draft, and iteratively improve this code, with human verification and approval of all code and decisions.
