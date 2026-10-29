"""
Cl-C-C-Cl torsion of 1,2-dichloroethane (DCE): gas-phase scans and Sage sampling.

Supports the report's case for keeping NMP + DCE. Four steps::

    pixi run -e openff dce-scan sage
    ../ml-peg/.venv/bin/python src/liquid_mixtures_inputs/dce_torsion_scan.py mace
    pixi run -e openff dce-scan flips
    pixi run -e openff dce-scan plot

``sage`` runs a relaxed scan with the benchmark force field and writes the
optimised structures. ``mace`` re-optimises them with MACE-OFF23(M) at the same
fixed dihedrals. MACE is in neither pixi environment, so this step uses the
ml-peg virtual environment next to this repo. ``flips`` measures how often DCE
molecules switch between trans and gauche in the Sage reference runs, and
``plot`` draws both scans and prints the numbers quoted in the report.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

# Constants are defined here rather than imported from the package, because the
# MACE step runs this file directly in the ml-peg virtual environment.
SMILES = "ClCCCl"
DIHEDRAL = [0, 1, 2, 3]  # Cl-C-C-Cl in the OpenFF atom order of SMILES
ANGLES = np.arange(0, 181, 15)
GAUCHE_MAX = 120.0  # |dihedral| below this counts as gauche
KJ_PER_MOL_PER_EV = 96.48533212
RT_KJ_PER_MOL = 8.314462618e-3 * 298.15
OUT = Path("data/sage_summary")


def moving_atoms(n_atoms: int) -> list[int]:
    """
    Find the atoms on the second carbon's side of the C-C bond (Cl and two H).

    Parameters
    ----------
    n_atoms
        Number of atoms in DCE.

    Returns
    -------
    list[int]
        Atom indices to rotate when setting the dihedral.
    """
    from openff.toolkit import Molecule

    molecule = Molecule.from_smiles(SMILES)
    assert molecule.n_atoms == n_atoms
    hydrogens = {
        neighbour.molecule_atom_index
        for neighbour in molecule.atom(2).bonded_atoms
        if neighbour.atomic_number == 1
    }
    return sorted({3} | hydrogens)


def scan_sage(out: Path) -> None:
    """
    Relaxed Sage scan, restraining the dihedral with a stiff harmonic term.

    Parameters
    ----------
    out
        Output directory for ``dce_scan_sage.xyz`` and ``dce_torsion_scan.csv``.
    """
    from ase import Atoms
    from ase.io import write
    from openff.toolkit import ForceField, Molecule
    import openmm
    from openmm import unit
    import pandas as pd

    from liquid_mixtures_inputs.build_boxes import FORCE_FIELD, NAGL_MODEL

    molecule = Molecule.from_smiles(SMILES)
    molecule.generate_conformers(n_conformers=1)
    molecule.assign_partial_charges(NAGL_MODEL)
    interchange = ForceField(FORCE_FIELD).create_interchange(
        molecule.to_topology(), charge_from_molecules=[molecule]
    )
    system = interchange.to_openmm_system()
    restraint = openmm.CustomTorsionForce(
        "0.5*k*d^2; d = min(abs(theta-theta0), 2*pi-abs(theta-theta0)); pi = 3.14159265"
    )
    restraint.addGlobalParameter("k", 1e5)
    restraint.addPerTorsionParameter("theta0")
    restraint.addTorsion(*DIHEDRAL, [0.0])
    restraint.setForceGroup(1)
    system.addForce(restraint)
    context = openmm.Context(
        system,
        openmm.VerletIntegrator(0.001),
        openmm.Platform.getPlatformByName("Reference"),
    )
    symbols = [atom.symbol for atom in molecule.atoms]
    atoms = Atoms(symbols, positions=molecule.conformers[0].m_as("angstrom"))
    mask = np.zeros(len(atoms), bool)
    mask[moving_atoms(len(atoms))] = True
    frames, rows = [], []
    for angle in ANGLES[::-1]:  # from trans down to cis
        atoms.set_dihedral(*DIHEDRAL, angle, mask=mask)
        restraint.setTorsionParameters(0, *DIHEDRAL, [np.radians(angle)])
        restraint.updateParametersInContext(context)
        context.setPositions(atoms.get_positions() * 0.1 * unit.nanometer)
        openmm.LocalEnergyMinimizer.minimize(context, 1e-4, 0)
        state = context.getState(getEnergy=True, getPositions=True, groups={0})
        atoms.set_positions(
            state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        )
        frame = atoms.copy()
        frame.info["dihedral"] = float(angle)
        frames.append(frame)
        rows.append(
            {
                "dihedral": float(angle),
                "sage": state.getPotentialEnergy().value_in_unit(
                    unit.kilojoule_per_mole
                ),
                "sage_achieved": atoms.get_dihedral(*DIHEDRAL),
            }
        )
    write(out / "dce_scan_sage.xyz", frames[::-1], format="extxyz")
    table = pd.DataFrame(rows[::-1])
    table["sage"] -= table.sage.iloc[-1]
    table.to_csv(out / "dce_torsion_scan.csv", index=False)


def scan_mace(out: Path) -> None:
    """
    Re-optimise the Sage scan structures with MACE-OFF23(M) at fixed dihedrals.

    Parameters
    ----------
    out
        Directory holding ``dce_scan_sage.xyz`` and ``dce_torsion_scan.csv``.
    """
    from ase.constraints import FixInternals
    from ase.io import read
    from ase.optimize import BFGS
    from mace.calculators import mace_off
    import pandas as pd

    calc = mace_off(model="medium", device="cpu", default_dtype="float64")
    table = pd.read_csv(out / "dce_torsion_scan.csv")
    energies = []
    for atoms in read(out / "dce_scan_sage.xyz", index=":"):
        angle = atoms.info["dihedral"]
        atoms.set_constraint(FixInternals(dihedrals_deg=[[angle, DIHEDRAL]]))
        atoms.calc = calc
        BFGS(atoms, logfile=None).run(fmax=1e-3, steps=2000)
        energies.append(atoms.get_potential_energy() * KJ_PER_MOL_PER_EV)
    table["mace_off_m"] = np.array(energies) - energies[-1]
    table.to_csv(out / "dce_torsion_scan.csv", index=False)


def flip_statistics(data_dir: Path, runs_dir: Path, out: Path) -> None:
    """
    Fraction of DCE molecules changing trans/gauche state between frames.

    For a two-state process the probability of a change after a lag Δt is
    2·p_t·p_g·(1 - exp(-Δt/τ)), which gives the relaxation time τ.

    Parameters
    ----------
    data_dir
        Benchmark data directory with the distributed boxes.
    runs_dir
        Sage runs (``<box>/rep*_frames.xyz``, one frame per 100 ps).
    out
        Output directory for ``dce_flips.csv``.
    """
    from ase.io import read
    import pandas as pd

    from liquid_mixtures_inputs.constants import EQUILIBRATION_PS
    from liquid_mixtures_inputs.run_sage_benchmark import FRAME_EVERY_PS

    rows = []
    for box_path in sorted((data_dir / "boxes").glob("*.xyz")):
        info = read(box_path).info
        keys = "ab"[: int(info["n_components"])]
        if SMILES not in [info[f"smiles_{k}"] for k in keys]:
            continue
        offset = 0
        for k in keys:
            if info[f"smiles_{k}"] == SMILES:
                start, n_mol = offset, int(info[f"n_{k}"])
                size = int(info[f"n_atoms_{k}"])
            offset += int(info[f"n_{k}"]) * int(info[f"n_atoms_{k}"])
        for path in sorted((runs_dir / info["box_id"]).glob("rep*_frames.xyz")):
            frames = [
                f
                for f in read(path, index=":")
                if f.info["time_ps"] >= EQUILIBRATION_PS
            ]
            gauche = []
            for frame in frames:
                states = []
                for molecule_index in range(n_mol):
                    first = start + molecule_index * size
                    assert frame[first].symbol == "Cl" == frame[first + 3].symbol
                    angle = frame.get_dihedral(*(first + i for i in DIHEDRAL), mic=True)
                    states.append(abs(((angle + 180) % 360) - 180) < GAUCHE_MAX)
                gauche.append(states)
            gauche = np.array(gauche)
            gauche_fraction = gauche.mean()
            changed = (gauche[1:] != gauche[:-1]).mean()
            limit = 2 * gauche_fraction * (1 - gauche_fraction)
            ratio = changed / limit
            rows.append(
                {
                    "box_id": info["box_id"],
                    "replicate": path.name.split("_")[0],
                    "n_dce": n_mol,
                    "gauche_fraction": gauche_fraction,
                    "changed_per_frame": changed,
                    "uncorrelated_limit": limit,
                    "tau_ps": -FRAME_EVERY_PS / np.log(1 - ratio)
                    if ratio < 1
                    else np.nan,
                }
            )
    pd.DataFrame(rows).to_csv(out / "dce_flips.csv", index=False)


def torsion_summary(table, column: str) -> tuple[float, float]:
    """
    Return the gauche minimum and trans-to-gauche barrier of one scan.

    Parameters
    ----------
    table
        ``dce_torsion_scan.csv`` as a data frame (energies relative to trans).
    column
        Scan column: ``sage`` or ``mace_off_m``.

    Returns
    -------
    tuple[float, float]
        Gauche energy and barrier height above trans, in kJ/mol.
    """
    gauche = table[(table.dihedral > 30) & (table.dihedral < 100)][column].min()
    barrier = table[(table.dihedral > 90) & (table.dihedral < 160)][column].max()
    return float(gauche), float(barrier)


def plot(out: Path) -> None:
    """
    Plot both scans and print the numbers used in the report.

    Parameters
    ----------
    out
        Directory with ``dce_torsion_scan.csv`` and ``dce_flips.csv``.
    """
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    table = pd.read_csv(out / "dce_torsion_scan.csv")
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for column, label in (("mace_off_m", "MACE-OFF23(M)"), ("sage", "Sage 2.3.0")):
        # The scan is symmetric about 0 and 180 degrees; show -180 to 180.
        x = np.r_[-table.dihedral[::-1], table.dihedral[1:]]
        y = np.r_[table[column][::-1], table[column][1:]]
        ax.plot(x, y, "o-", ms=3, label=label)
    ax.set_xlabel("Cl-C-C-Cl dihedral (°)")
    ax.set_ylabel("energy relative to trans (kJ/mol)")
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "figures" / "dce_torsion_scan.png", dpi=110)
    plt.close(fig)

    g_sage, b_sage = torsion_summary(table, "sage")
    g_ref, b_ref = torsion_summary(table, "mace_off_m")
    slowdown = np.exp((b_sage - b_ref) / RT_KJ_PER_MOL)
    print(f"gauche - trans: Sage {g_sage:.1f}, MACE-OFF23(M) {g_ref:.1f} kJ/mol")
    print(f"barrier above trans: Sage {b_sage:.1f}, MACE-OFF23(M) {b_ref:.1f} kJ/mol")
    print(f"Arrhenius slowdown at 298 K: {slowdown:.1f}x")
    max_dev = (table.sage_achieved - table.dihedral).abs()
    print(f"max |achieved - target| Sage dihedral: {max_dev.max():.3f} deg")
    flips = pd.read_csv(out / "dce_flips.csv")
    by_box = flips.groupby("box_id")[
        [
            "n_dce",
            "gauche_fraction",
            "changed_per_frame",
            "uncorrelated_limit",
            "tau_ps",
        ]
    ].mean()
    print(by_box.round(3).to_string())
    print(
        f"tau median {flips.tau_ps.median():.0f} ps; reference-like barrier: "
        f"{flips.tau_ps.median() / slowdown:.0f} ps"
    )


def main() -> None:
    """Run one step."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=["sage", "mace", "flips", "plot"])
    parser.add_argument("--output-dir", type=Path, default=OUT)
    parser.add_argument("--data-dir", type=Path, default=Path("data/binary_mixtures"))
    parser.add_argument("--runs-dir", type=Path, default=Path("data/sage_runs"))
    args = parser.parse_args()
    if args.step == "sage":
        scan_sage(args.output_dir)
    elif args.step == "mace":
        scan_mace(args.output_dir)
    elif args.step == "flips":
        flip_statistics(args.data_dir, args.runs_dir, args.output_dir)
    else:
        plot(args.output_dir)


if __name__ == "__main__":
    main()
