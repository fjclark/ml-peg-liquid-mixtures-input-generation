"""
Build Sage-equilibrated starting boxes for the binary mixtures benchmark.

Run in the ``openff`` pixi environment (the OpenFF stack is conda-only)::

    pixi run -e openff build-boxes --selection-dir selection

Each box is packed with packmol at the experimental density, parametrised with
openff-2.3.0 (NAGL AM1-BCC charges, TIP3P water from the force field library
charges), minimised, annealed from 450 K to the target temperature over 1.6 ns
and run for 500 ps NVT at the experimental density. Every molecule starts from
the lowest-MMFF-energy of 10 ETKDG conformers. See DECISIONS.md (box build) for
why the anneal is needed.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import json
import logging
import multiprocessing
from pathlib import Path
import shutil
import subprocess
import tempfile

from ase import Atoms
from ase.io import read, write
import numpy as np
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import Quantity
import openmm
from openmm import unit
import pandas as pd
from rdkit.Chem import AllChem
import yaml

from liquid_mixtures_inputs.constants import AVOGADRO

HERE = Path(__file__).parent
FORCE_FIELD = "openff-2.3.0.offxml"
NAGL_MODEL = "openff-gnn-am1bcc-1.0.0.pt"
WATER = "O"
PACK_TOLERANCE = 2.0  # Angstrom
TIMESTEP_FS = 2.0
FRICTION_PER_PS = 1.0
WARMUP_PS = 20.0  # 1 fs warm-up after minimisation
NVT_PS = 500.0
# NVT anneal from ANNEAL_T_HOT down to the target temperature in equal steps,
# so slow conformer populations (e.g. Cl-C-C-Cl in 1,2-dichloroethane, O-C-C-N in
# DMAE) start near equilibrium instead of all in the packed conformer.
ANNEAL_T_HOT = 450.0  # K
ANNEAL_STEPS = 16
ANNEAL_STEP_PS = 100.0

LOGGER = logging.getLogger("build_boxes")


@dataclass
class BoxSpec:
    """
    Composition and target state of one simulation box.

    Attributes
    ----------
    box_id
        Unique identifier, used as the file stem.
    kind
        ``"mixture"`` or ``"pure"``.
    smiles
        SMILES of each component (one for pure boxes).
    counts
        Number of molecules of each component.
    temperature
        Temperature in K.
    density
        Density used for packing and NVT equilibration, in g/mL.
    density_source
        Where ``density`` comes from.
    system_id
        Mixture system identifier, empty for pure boxes.
    x
        Mole fraction of the first component.
    """

    box_id: str
    kind: str
    smiles: list[str]
    counts: list[int]
    temperature: float
    density: float
    density_source: str
    system_id: str = ""
    x: float = 1.0


def file_slug(name: str) -> str:
    """
    Turn a molecule name into a file-name-safe slug.

    Parameters
    ----------
    name
        Molecule name.

    Returns
    -------
    str
        Lower-case slug.
    """
    characters = [c if c.isalnum() else "_" for c in name.lower()]
    return "_".join(part for part in "".join(characters).split("_") if part)


def molar_mass(smiles: str) -> float:
    """
    Return the molar mass of a molecule in g/mol.

    Parameters
    ----------
    smiles
        SMILES string.

    Returns
    -------
    float
        Molar mass.
    """
    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    return sum(atom.mass.m_as("dalton") for atom in molecule.atoms)


def load_specs(selection_dir: Path, selection_yml: Path) -> list[BoxSpec]:
    """
    Read the mixture and pure box plans written by ``select_systems.py``.

    Parameters
    ----------
    selection_dir
        Directory holding ``systems.csv`` and ``pure.csv``.
    selection_yml
        Curated ``selection.yml`` (for component names).

    Returns
    -------
    list[BoxSpec]
        Mixture boxes followed by pure boxes.
    """
    selection = yaml.safe_load(selection_yml.read_text())
    names = {}
    for system in selection["systems"]:
        names[system["a"]] = system["names"][0]
        names[system["b"]] = system["names"][1]
    systems = pd.read_csv(selection_dir / "systems.csv")
    pure = pd.read_csv(selection_dir / "pure.csv")
    pure_density = {
        (row.smiles, round(row.temperature, 2)): row.density
        for row in pure.itertuples()
    }

    specs = []
    for row in systems.itertuples():
        density, source = row.density_exp, "experiment (Redlich-Kister V^E fit)"
        if pd.isna(density):
            # Ideal mixing from the pure densities (V^E is small next to the
            # equilibration tolerance; NPT finds the model's own density anyway).
            temperature = round(row.temperature, 2)
            mass_a, mass_b = molar_mass(row.a), molar_mass(row.b)
            ideal_volume = (
                row.x_box * mass_a / pure_density[(row.a, temperature)]
                + (1 - row.x_box) * mass_b / pure_density[(row.b, temperature)]
            )
            density = (row.x_box * mass_a + (1 - row.x_box) * mass_b) / ideal_volume
            source = "ideal mixing of pure densities (no mixture density data)"
        specs.append(
            BoxSpec(
                box_id=f"{row.system_id}_x{round(row.x_box * 1000):03d}",
                kind="mixture",
                smiles=[row.a, row.b],
                counts=[int(row.n_a), int(row.n_b)],
                temperature=float(row.temperature),
                density=float(density),
                density_source=source,
                system_id=row.system_id,
                x=float(row.x_box),
            )
        )
    for row in pure.itertuples():
        specs.append(
            BoxSpec(
                box_id=f"pure_{file_slug(names[row.smiles])}_{round(row.temperature)}",
                kind="pure",
                smiles=[row.smiles],
                counts=[int(row.n_molecules)],
                temperature=float(row.temperature),
                density=float(row.density),
                density_source=str(row.source),
            )
        )
    return specs


def set_low_energy_conformer(molecule: Molecule, n_conformers: int = 10) -> None:
    """
    Keep the lowest-energy of several conformers as the molecule's only conformer.

    A single ETKDG conformer can sit in a high-energy basin that is slow to
    leave (for tetrahydropyran it is often a twist-boat, 27 kJ/mol above the
    chair in MMFF, which flips to chair in Sage over ~1 ns). Packmol copies the
    conformer into every molecule, so the whole box would share the basin.

    Parameters
    ----------
    molecule
        Molecule, modified in place.
    n_conformers
        Number of ETKDG conformers to compare after MMFF minimisation. They
        are embedded with RDKit directly, because ``generate_conformers``
        prunes by RMSD and returns a single conformer for small molecules.
    """
    rdkit_molecule = molecule.to_rdkit()
    AllChem.EmbedMultipleConfs(rdkit_molecule, n_conformers, randomSeed=1)
    results = AllChem.MMFFOptimizeMoleculeConfs(rdkit_molecule, maxIters=5000)
    energies = [energy for _, energy in results]
    best = rdkit_molecule.GetConformer(int(np.argmin(energies))).GetPositions()
    molecule._conformers = None
    molecule.add_conformer(Quantity(best, "angstrom"))


def prepare_molecule(smiles: str) -> Molecule:
    """
    Build a molecule with a conformer and (except water) NAGL charges.

    Parameters
    ----------
    smiles
        SMILES string.

    Returns
    -------
    Molecule
        OpenFF molecule.
    """
    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    set_low_energy_conformer(molecule)
    # Water takes the force field's TIP3P library charges instead.
    if smiles != WATER:
        molecule.assign_partial_charges(NAGL_MODEL)
    return molecule


def box_edge(spec: BoxSpec) -> float:
    """
    Return the cubic box edge (Angstrom) at the target density.

    Parameters
    ----------
    spec
        Box specification.

    Returns
    -------
    float
        Edge length.
    """
    mass = sum(n * molar_mass(s) for s, n in zip(spec.smiles, spec.counts, strict=True))
    volume_cm3 = mass / AVOGADRO / spec.density
    return volume_cm3 ** (1 / 3) * 1e8


def pack(molecules: list[Molecule], counts: list[int], edge: float, seed: int):
    """
    Pack molecules into a cubic box with packmol.

    Parameters
    ----------
    molecules
        One molecule (with a conformer) per component.
    counts
        Number of copies of each.
    edge
        Box edge in Angstrom.
    seed
        Packmol random seed.

    Returns
    -------
    np.ndarray
        Packed coordinates in Angstrom, in topology order.
    """
    packmol = shutil.which("packmol")
    if packmol is None:
        raise RuntimeError("packmol not found; run inside the pixi environment")
    margin = PACK_TOLERANCE / 2
    with tempfile.TemporaryDirectory() as temporary:
        work_dir = Path(temporary)
        lines = [
            f"tolerance {PACK_TOLERANCE}",
            "filetype xyz",
            "output packed.xyz",
            f"seed {seed}",
            "nloop 400",
        ]
        low, high = margin, edge - margin
        for i, (molecule, count) in enumerate(zip(molecules, counts, strict=True)):
            xyz = molecule.conformers[0].m_as("angstrom")
            symbols = [atom.symbol for atom in molecule.atoms]
            body = "\n".join(
                f"{symbol} {x:.4f} {y:.4f} {z:.4f}"
                for symbol, (x, y, z) in zip(symbols, xyz, strict=True)
            )
            (work_dir / f"mol{i}.xyz").write_text(f"{len(symbols)}\nmol{i}\n{body}\n")
            lines += [
                f"structure mol{i}.xyz",
                f"  number {count}",
                f"  inside box {low:.3f} {low:.3f} {low:.3f} "
                f"{high:.3f} {high:.3f} {high:.3f}",
                "end structure",
            ]
        (work_dir / "pack.inp").write_text("\n".join(lines) + "\n")
        with open(work_dir / "pack.inp") as stdin:
            result = subprocess.run(
                [packmol], stdin=stdin, cwd=work_dir, capture_output=True, text=True
            )
        packed = work_dir / "packed.xyz"
        if not packed.exists():
            raise RuntimeError(f"packmol failed:\n{result.stdout[-2000:]}")
        if "Success!" not in result.stdout:
            LOGGER.warning("packmol did not reach full convergence; minimising anyway")
        return read(packed).get_positions()


def make_simulation(
    molecules: list[Molecule],
    counts: list[int],
    positions: np.ndarray,
    edge: float | np.ndarray,
    temperature: float,
    barostat: bool,
    seed: int,
    platform: str = "CUDA",
):
    """
    Parametrise a box with Sage and build an OpenMM simulation.

    Parameters
    ----------
    molecules
        One charged molecule per component.
    counts
        Number of copies of each.
    positions
        Coordinates in Angstrom, in topology order.
    edge
        Cubic box edge in Angstrom, or a 3x3 cell matrix.
    temperature
        Temperature in K.
    barostat
        Add a 1 atm Monte Carlo barostat (NPT) if ``True``.
    seed
        Seed for the integrator and barostat.
    platform
        OpenMM platform name.

    Returns
    -------
    tuple[openmm.app.Simulation, int]
        The simulation and the number of degrees of freedom.
    """
    from openmm import app

    copies = [
        molecule
        for molecule, count in zip(molecules, counts, strict=True)
        for _ in range(count)
    ]
    topology = Topology.from_molecules(copies)
    cell = np.diag([edge] * 3) if np.ndim(edge) == 0 else np.asarray(edge)
    topology.box_vectors = Quantity(cell, "angstrom")
    force_field = ForceField(FORCE_FIELD)
    # Water has no NAGL charges and takes the force field's TIP3P library charges.
    charged = [m for m in molecules if m.partial_charges is not None]
    interchange = force_field.create_interchange(
        topology, charge_from_molecules=charged
    )
    interchange.positions = Quantity(positions, "angstrom")
    system = interchange.to_openmm_system()
    if barostat:
        force = openmm.MonteCarloBarostat(1.0 * unit.atmosphere, temperature, 25)
        force.setRandomNumberSeed(seed)
        system.addForce(force)
    integrator = openmm.LangevinMiddleIntegrator(
        temperature * unit.kelvin,
        FRICTION_PER_PS / unit.picosecond,
        TIMESTEP_FS * unit.femtosecond,
    )
    integrator.setRandomNumberSeed(seed)
    properties = {"Precision": "mixed"} if platform == "CUDA" else {}
    simulation = app.Simulation(
        interchange.to_openmm_topology(),
        system,
        integrator,
        openmm.Platform.getPlatformByName(platform),
        properties,
    )
    simulation.context.setPositions((positions * 0.1) * unit.nanometer)
    simulation.context.setPeriodicBoxVectors(*(cell * 0.1 * unit.nanometer))
    dof = 3 * system.getNumParticles() - system.getNumConstraints()
    if any(isinstance(f, openmm.CMMotionRemover) for f in system.getForces()):
        dof -= 3
    return simulation, dof


def equilibrate_nvt(simulation, temperature: float, seed: int) -> None:
    """
    Minimise, warm up at 1 fs, anneal, then run NVT at the packing density.

    Packmol copies one conformer into every molecule. The anneal (450 K cooled
    to the target over 1.6 ns) lets slow torsions reach their equilibrium
    populations: for 1,2-dichloroethane, opposite starting conformers gave
    matching gauche fractions afterwards (0.48 vs 0.48, against 0.65 vs 0.21
    without the anneal).

    Parameters
    ----------
    simulation
        OpenMM simulation (without barostat).
    temperature
        Temperature in K.
    seed
        Velocity seed.
    """
    simulation.minimizeEnergy(tolerance=10 * unit.kilojoule_per_mole / unit.nanometer)
    simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin, seed)
    integrator = simulation.integrator
    integrator.setStepSize(1.0 * unit.femtosecond)
    simulation.step(round(WARMUP_PS * 1000))
    integrator.setStepSize(TIMESTEP_FS * unit.femtosecond)
    for anneal_temperature in np.linspace(ANNEAL_T_HOT, temperature, ANNEAL_STEPS):
        integrator.setTemperature(anneal_temperature * unit.kelvin)
        simulation.step(round(ANNEAL_STEP_PS * 1000 / TIMESTEP_FS))
    integrator.setTemperature(temperature * unit.kelvin)
    simulation.step(round(NVT_PS * 1000 / TIMESTEP_FS))


def whole_molecule_positions(
    positions: np.ndarray, counts: list[int], sizes: list[int], edge: float
) -> np.ndarray:
    """
    Make molecules whole and put each molecule's centre inside the box.

    Parameters
    ----------
    positions
        Coordinates in Angstrom (molecules may be split across the boundary).
    counts
        Number of molecules of each component.
    sizes
        Atoms per molecule of each component.
    edge
        Cubic box edge in Angstrom.

    Returns
    -------
    np.ndarray
        Coordinates with whole molecules.
    """
    whole = positions.copy()
    start = 0
    for count, size in zip(counts, sizes, strict=True):
        for _ in range(count):
            molecule = whole[start : start + size]
            # Unwrap every atom next to the first, then wrap the centre into the box.
            molecule -= edge * np.round((molecule - molecule[0]) / edge)
            molecule -= edge * np.floor(molecule.mean(axis=0) / edge)
            whole[start : start + size] = molecule
            start += size
    return whole


def build_box(
    spec: BoxSpec, seed: int, platform: str = "CUDA", molecules=None
) -> tuple[Atoms, list[Molecule]]:
    """
    Pack, parametrise and NVT-equilibrate one box.

    Parameters
    ----------
    spec
        Box specification.
    seed
        Seed for packmol, velocities and the thermostat.
    platform
        OpenMM platform name.
    molecules
        Pre-built molecules (one per component), to reuse charges.

    Returns
    -------
    tuple[Atoms, list[Molecule]]
        The equilibrated structure and the molecules used.
    """
    molecules = molecules or [prepare_molecule(s) for s in spec.smiles]
    for molecule in molecules:
        if not molecule.n_conformers:
            set_low_energy_conformer(molecule)
    edge = box_edge(spec)
    positions = pack(molecules, spec.counts, edge, seed)
    simulation, _ = make_simulation(
        molecules, spec.counts, positions, edge, spec.temperature, False, seed, platform
    )
    equilibrate_nvt(simulation, spec.temperature, seed)
    state = simulation.context.getState(getPositions=True)
    final = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    sizes = [molecule.n_atoms for molecule in molecules]
    final = whole_molecule_positions(np.asarray(final), spec.counts, sizes, edge)
    symbols = [
        atom.symbol
        for molecule, count in zip(molecules, spec.counts, strict=True)
        for _ in range(count)
        for atom in molecule.atoms
    ]
    atoms = Atoms(symbols, positions=final, cell=[edge] * 3, pbc=True)
    atoms.info.update(box_info(spec, molecules, seed))
    return atoms, molecules


def box_info(spec: BoxSpec, molecules: list[Molecule], seed: int) -> dict:
    """
    Metadata stored in ``atoms.info`` of each box.

    Parameters
    ----------
    spec
        Box specification.
    molecules
        Molecules used, one per component.
    seed
        Seed used to build the box.

    Returns
    -------
    dict
        Scalar metadata (strings, ints and floats only, for extxyz).
    """
    info = {
        "box_id": spec.box_id,
        "kind": spec.kind,
        "system_id": spec.system_id or "none",
        "n_components": len(spec.smiles),
        "exp_temperature": spec.temperature,
        "pack_density": round(spec.density, 6),
        "density_source": spec.density_source,
        "x_a": round(spec.x, 6),
        "n_molecules": int(sum(spec.counts)),
        "charge": 0,
        "spin": 1,
        "seed": seed,
        "builder": f"{FORCE_FIELD} {NAGL_MODEL} packmol NVT {NVT_PS:g} ps",
    }
    for key, smiles, count, molecule in zip(
        "ab", spec.smiles, spec.counts, molecules, strict=False
    ):
        info[f"smiles_{key}"] = smiles
        info[f"n_{key}"] = int(count)
        info[f"n_atoms_{key}"] = molecule.n_atoms
        info[f"mapped_smiles_{key}"] = molecule.to_smiles(mapped=True)
    return info


def molecules_from_info(info: dict) -> tuple[list[Molecule], list[int]]:
    """
    Rebuild the charged molecules and counts of a box from ``atoms.info``.

    Parameters
    ----------
    info
        ``atoms.info`` written by :func:`box_info`.

    Returns
    -------
    tuple[list[Molecule], list[int]]
        Molecules in topology atom order, and their counts.
    """
    molecules, counts = [], []
    for key in "ab"[: int(info["n_components"])]:
        molecule = Molecule.from_mapped_smiles(
            info[f"mapped_smiles_{key}"], allow_undefined_stereo=True
        )
        if info[f"smiles_{key}"] != WATER:
            molecule.assign_partial_charges(NAGL_MODEL)
        molecules.append(molecule)
        counts.append(int(info[f"n_{key}"]))
    return molecules, counts


def _build_and_write(spec: BoxSpec, seed: int, platform: str, box_dir: Path) -> int:
    """
    Build one box and write it as extxyz (worker process entry point).

    Parameters
    ----------
    spec
        Box specification.
    seed
        Build seed.
    platform
        OpenMM platform name.
    box_dir
        Output directory.

    Returns
    -------
    int
        Number of atoms in the box.
    """
    atoms, _ = build_box(spec, seed, platform)
    # Write to a hidden file first, so an interrupted build is never mistaken
    # for a finished box.
    partial = box_dir / f".{spec.box_id}.xyz"
    write(partial, atoms, format="extxyz")
    partial.rename(box_dir / f"{spec.box_id}.xyz")
    return len(atoms)


def main() -> None:
    """Build every box listed by ``select_systems.py``."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-dir", type=Path, required=True)
    parser.add_argument("--selection-yml", type=Path, default=HERE / "selection.yml")
    parser.add_argument("--output-dir", type=Path, default=Path("data/binary_mixtures"))
    parser.add_argument("--box-id", nargs="*", help="Only build these boxes")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--platform", default="CUDA")
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    specs = load_specs(args.selection_dir, args.selection_yml)
    box_dir = args.output_dir / "boxes"
    box_dir.mkdir(parents=True, exist_ok=True)
    for name in ("systems.csv", "pure.csv"):
        shutil.copy(args.selection_dir / name, args.output_dir / name)
    pd.DataFrame(
        [
            vars(spec)
            | {
                "smiles": ";".join(spec.smiles),
                "counts": ";".join(map(str, spec.counts)),
            }
            for spec in specs
        ]
    ).to_csv(args.output_dir / "boxes.csv", index=False)

    todo = [
        spec
        for spec in specs
        if (not args.box_id or spec.box_id in args.box_id)
        and not (box_dir / f"{spec.box_id}.xyz").exists()
    ]
    LOGGER.info("%d boxes to build with %d workers", len(todo), args.workers)
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(args.workers, mp_context=context) as pool:
        futures = {
            pool.submit(
                _build_and_write, spec, args.seed, args.platform, box_dir
            ): spec.box_id
            for spec in todo
        }
        for future in as_completed(futures):
            LOGGER.info("built %s (%d atoms)", futures[future], future.result())
    (args.output_dir / "README.json").write_text(
        json.dumps(
            {
                "force_field": FORCE_FIELD,
                "charges": NAGL_MODEL,
                "nvt_ps": NVT_PS,
                "timestep_fs": TIMESTEP_FS,
                "seed": args.seed,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
