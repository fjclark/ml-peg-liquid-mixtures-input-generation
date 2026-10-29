"""
Run the binary mixtures benchmark with Sage (openff-2.3.0) as a reference model.

Run in the ``openff`` pixi environment::

    pixi run -e openff run-sage --data-dir data/binary_mixtures \
        --output-dir data/sage_runs

Every box is run as several independent replicates with the benchmark protocol
(0.2 ns NPT equilibration followed by 2 ns NPT production, 298.15 K or the box's
``exp_temperature``, 1 atm). Replicate 0 starts from the distributed box.
Later replicates are re-packed with a new packmol seed and re-equilibrated
(500 ps NVT), so their starting configurations are independent.

Each run writes ``<box_id>/rep<r>.csv`` in the same scalar format as the MLIP
calc (time in ps, energies in eV, volume in Angstrom^3, density in g/mL), plus
``rep<r>_frames.xyz`` with one frame every 100 ps for structural checks.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import multiprocessing
from pathlib import Path
import time

from ase import Atoms
from ase.io import read, write
import numpy as np
from openmm import unit

from liquid_mixtures_inputs.build_boxes import (
    TIMESTEP_FS,
    BoxSpec,
    build_box,
    make_simulation,
    molecules_from_info,
    prepare_molecule,
)

LOGGER = logging.getLogger("run_sage_benchmark")

EQUILIBRATION_PS = 200.0
PRODUCTION_PS = 2000.0
LOG_EVERY_PS = 0.1
FRAME_EVERY_PS = 100.0
KJ_PER_MOL_TO_EV = 1 / 96.48533212
AMU_PER_A3_TO_G_ML = 1.66053906660
CSV_HEADER = "step,time_ps,epot_eV,ekin_eV,temperature_K,volume_A3,density_g_cm3"


class ScalarReporter:
    """
    OpenMM reporter writing the benchmark's scalar CSV format.

    Parameters
    ----------
    path
        CSV path.
    interval
        Reporting interval in steps.
    dof
        Degrees of freedom, for the instantaneous temperature.
    total_mass
        Total mass in amu.
    """

    def __init__(self, path: Path, interval: int, dof: int, total_mass: float):
        """
        Open the CSV and write its header.

        Parameters
        ----------
        path
            CSV path.
        interval
            Reporting interval in steps.
        dof
            Degrees of freedom, for the instantaneous temperature.
        total_mass
            Total mass in amu.
        """
        self._file = open(path, "w")  # noqa: SIM115
        self._file.write(CSV_HEADER + "\n")
        self._interval = interval
        self._dof = dof
        self._mass = total_mass

    def describeNextReport(self, simulation):  # noqa: N802
        """
        Tell OpenMM when the next report is due and what it needs.

        Parameters
        ----------
        simulation
            The simulation.

        Returns
        -------
        dict
            Report request.
        """
        steps = self._interval - simulation.currentStep % self._interval
        return {"steps": steps, "periodic": None, "include": ["energy"]}

    def report(self, simulation, state) -> None:
        """
        Append one line of scalars.

        Parameters
        ----------
        simulation
            The simulation.
        state
            OpenMM state with energies.
        """
        epot = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        ekin = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        volume = state.getPeriodicBoxVolume().value_in_unit(unit.angstrom**3)
        temperature = 2 * ekin / (self._dof * 0.00831446261815324)
        step = simulation.currentStep
        self._file.write(
            f"{step},{step * TIMESTEP_FS / 1000:.4f},"
            f"{epot * KJ_PER_MOL_TO_EV:.6f},{ekin * KJ_PER_MOL_TO_EV:.6f},"
            f"{temperature:.3f},{volume:.4f},"
            f"{self._mass / volume * AMU_PER_A3_TO_G_ML:.6f}\n"
        )

    def close(self) -> None:
        """Close the CSV."""
        self._file.close()


def spec_from_info(info: dict) -> BoxSpec:
    """
    Rebuild a box specification from the ``atoms.info`` of a distributed box.

    Parameters
    ----------
    info
        ``atoms.info`` written by ``build_boxes.box_info``.

    Returns
    -------
    BoxSpec
        Specification used to re-pack independent replicates.
    """
    keys = "ab"[: int(info["n_components"])]
    return BoxSpec(
        box_id=info["box_id"],
        kind=info["kind"],
        smiles=[info[f"smiles_{k}"] for k in keys],
        counts=[int(info[f"n_{k}"]) for k in keys],
        temperature=float(info["exp_temperature"]),
        density=float(info["pack_density"]),
        density_source=str(info["density_source"]),
        system_id=info["system_id"],
        x=float(info["x_a"]),
    )


def run_replicate(
    box_path: Path,
    replicate: int,
    out_dir: Path,
    platform: str,
    run_ps: float = EQUILIBRATION_PS + PRODUCTION_PS,
):
    """
    Run one replicate of one box with the benchmark protocol.

    Parameters
    ----------
    box_path
        Distributed extxyz box.
    replicate
        Replicate index; 0 starts from the distributed box.
    out_dir
        Output directory for this box.
    platform
        OpenMM platform name.
    run_ps
        Total NPT length (equilibration plus production) in ps.

    Returns
    -------
    float
        Wall time in seconds.
    """
    start = time.time()
    atoms = read(box_path)
    info = atoms.info
    seed = 1000 + replicate
    molecules, counts = molecules_from_info(info)
    if replicate == 0:
        positions, cell = atoms.get_positions(), np.asarray(atoms.get_cell())
    else:
        # Independent start: new packing and NVT pre-equilibration.
        spec = spec_from_info(info)
        built, _ = build_box(spec, seed, platform, molecules)
        positions, cell = built.get_positions(), np.asarray(built.get_cell())
    temperature = float(info["exp_temperature"])
    simulation, dof = make_simulation(
        molecules, counts, positions, cell, temperature, True, seed, platform
    )
    simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin, seed)

    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_csv = out_dir / f".rep{replicate}.csv"
    tmp_frames = out_dir / f".rep{replicate}_frames.xyz"
    tmp_frames.unlink(missing_ok=True)
    log_steps = round(LOG_EVERY_PS * 1000 / TIMESTEP_FS)
    frame_steps = round(FRAME_EVERY_PS * 1000 / TIMESTEP_FS)
    reporter = ScalarReporter(tmp_csv, log_steps, dof, float(atoms.get_masses().sum()))
    simulation.reporters.append(reporter)
    total_steps = round(run_ps * 1000 / TIMESTEP_FS)
    symbols = atoms.get_chemical_symbols()
    while simulation.currentStep < total_steps:
        simulation.step(min(frame_steps, total_steps - simulation.currentStep))
        state = simulation.context.getState(getPositions=True)
        frame = Atoms(
            symbols,
            positions=state.getPositions(asNumpy=True).value_in_unit(unit.angstrom),
            cell=state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom),
            pbc=True,
        )
        frame.info["time_ps"] = simulation.currentStep * TIMESTEP_FS / 1000
        write(tmp_frames, frame, format="extxyz", append=True)
    reporter.close()
    tmp_frames.rename(out_dir / f"rep{replicate}_frames.xyz")
    tmp_csv.rename(out_dir / f"rep{replicate}.csv")
    return time.time() - start


def _job(box_path: Path, replicate: int, out_dir: Path, platform: str, run_ps: float):
    """
    Worker entry point that also reports which job finished.

    Parameters
    ----------
    box_path
        Distributed extxyz box.
    replicate
        Replicate index.
    out_dir
        Output directory for this box.
    platform
        OpenMM platform name.
    run_ps
        Total NPT length in ps.

    Returns
    -------
    tuple[str, int, float]
        Box id, replicate and wall time.
    """
    wall = run_replicate(box_path, replicate, out_dir, platform, run_ps)
    return box_path.stem, replicate, wall


def main() -> None:
    """Run every box and replicate that has no finished output yet."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=5)
    parser.add_argument("--box-id", nargs="*", help="Only run these boxes")
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--platform", default="CUDA")
    parser.add_argument(
        "--run-ps",
        type=float,
        default=EQUILIBRATION_PS + PRODUCTION_PS,
        help="Total NPT length (shorten only for testing)",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    boxes = sorted((args.data_dir / "boxes").glob("*.xyz"))
    if args.box_id:
        boxes = [b for b in boxes if b.stem in args.box_id]
    # Replicate-major order, so complete replicate sets finish first.
    jobs = [
        (box, r)
        for r in range(args.replicates)
        for box in boxes
        if not (args.output_dir / box.stem / f"rep{r}.csv").exists()
    ]
    LOGGER.info("%d runs to do with %d workers", len(jobs), args.workers)
    # Warm the NAGL model cache once before forking workers.
    prepare_molecule("CCO")
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(args.workers, mp_context=context) as pool:
        futures = [
            pool.submit(
                _job, box, r, args.output_dir / box.stem, args.platform, args.run_ps
            )
            for box, r in jobs
        ]
        for done, future in enumerate(as_completed(futures), 1):
            try:
                box_id, replicate, wall = future.result()
                LOGGER.info(
                    "[%d/%d] %s rep%d done in %.0f s",
                    done,
                    len(jobs),
                    box_id,
                    replicate,
                    wall,
                )
            except Exception:  # keep the other runs going
                LOGGER.exception("run failed")


if __name__ == "__main__":
    main()
