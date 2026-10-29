"""
Select binary liquid mixtures for the binary_mixtures benchmark.

1. Download the OpenFF Sage physical-property sets (pinned commits) and screen
   every binary system against the hard filters, writing ``candidates.csv``.
2. Check the curated systems in ``selection.yml``.
3. Fetch the source papers' full isotherms and pure densities from the NIST
   TRC ThermoML archive, cross-check pure densities, fit Redlich-Kister curves
   and plan box sizes, writing ``systems.csv``, ``pure.csv`` and
   ``isotherms.csv``.
4. Write the selection report (``selection_report.py``).

Run with::

    pixi run -e analysis select --output-dir selection \
        --fast-set selection/inputs/dataset-fast-converging.csv

The reasons for every threshold below are in ``DECISIONS.md`` and in the
hard-filter table of ``selection/report.md``.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import math
from pathlib import Path
import statistics

from ase.io import read
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
import yaml

from liquid_mixtures_inputs.constants import AVOGADRO, GAS_CONSTANT_KJ_PER_MOL_K
from liquid_mixtures_inputs.thermoml import (
    Measurement,
    canonical_smiles,
    download,
    fetch_record,
    parse_record,
)
from ml_peg.calcs.utils.utils import BENCHMARK_DATA_DIR, download_s3_data

HERE = Path(__file__).parent
CACHE_DIR = BENCHMARK_DATA_DIR / "binary_mixtures_selection"

SAGE_COMMIT = "a6a061ba8ebc144e3fda46f5f7a9e7a352b9cd2f"
ASH_COMMIT = "00c9524d38af53725a1a18feb3433f5c191cdb4a"
RAW = "https://raw.githubusercontent.com/openforcefield"
SOURCES = {
    "sage-2.0-train": f"{RAW}/openff-sage/{SAGE_COMMIT}/data-set-curation/"
    "physical-property/optimizations/data-sets/sage-train-v1.csv",
    "rc2-train": f"{RAW}/ash-sage-rc2/{ASH_COMMIT}/01_download-data/physprop/final/"
    "output/training-set.csv",
    "rc2-validation": f"{RAW}/ash-sage-rc2/{ASH_COMMIT}/01_download-data/physprop/"
    "final/output/validation-set.csv",
}
SAGE_BENCHMARKS = {
    f"{split}-{property_name}": f"{RAW}/ash-sage-rc2/{ASH_COMMIT}/04_benchmark/"
    f"phys-prop/output/{split}/summary-benchmarks-{file_name}.csv"
    for split in ("training", "validation")
    for property_name, file_name in (
        ("hmix", "EnthalpyOfMixing"),
        ("density", "Density"),
    )
}
SAGE_FORCE_FIELDS = ("openff-2.0.0", "openff-2.1.0", "openff-2.2.1", "openff-2.3.0rc1")
# Sage 2.3.0rc2 is benchmarked under its fit name: 04_benchmark/forcefields/
# sage-2.3.0rc2-v3-k100.offxml is byte-identical to fb-fit-v3-single-mean-k100.offxml.
# Tabulated in the report only, not used by the repeat-SD filter.
SAGE_RC2 = "fb-fit-v3-single-mean-k100"

# The Sage sets encode 2-pyrrolidone as its imidic acid tautomer.
TAUTOMER_FIXES = {"OC1=NCCC1": "O=C1CCCN1"}

# Structure checks on every curated component. An OH on a non-aromatic carbon
# double-bonded to C or N is an enol or imidic acid, never the dominant form
# of our neutral liquids.
ENOL_OR_IMIDIC_ACID = Chem.MolFromSmarts("[OX2H1]-[#6;!a]=[#6,#7]")
TAUTOMER_ENUMERATOR = rdMolStandardize.TautomerEnumerator()
# RDKit's canonical tautomer follows scoring rules, not solution stability.
# Canonical SMILES reviewed and kept despite differing from it.
TAUTOMER_REVIEWED: set[str] = set()

# Hard filters
SCREEN_T_RANGE = (297.0, 299.5)  # K
MAX_COMPOSITION_MISMATCH = 0.05
ALLOWED_ELEMENTS = {"C", "H", "N", "O", "Cl"}
MAX_HEAVY_ATOMS = 10
VISCOUS_SMARTS = {
    "polyol": "[OX2H][CX4][CX4][OX2H]",
    "1,3-diol": "[OX2H][CX4][CX4][CX4][OX2H]",
    "glycerol-like": "[OX2H][CX4][CX4]([OX2H])",
    "polyamine": "[NX3;!$(N-C=O)][CX4][CX4][NX3;!$(N-C=O)]",
    "polyamine-3": "[NX3;!$(N-C=O)][CX4][CX4][CX4][NX3;!$(N-C=O)]",
    "ethanolamine (2+ OH)": "[OX2H][CX4][CX4][NX3][CX4][CX4][OX2H]",
}
UCST_LIMIT = 0.5 * GAS_CONSTANT_KJ_PER_MOL_K * 298.15  # RT/2 in kJ/mol
MAX_REPEAT_SD = {"hmix": 0.10, "density": 0.001}  # kJ/mol, g/mL
MIN_SIGNAL = 0.4  # kJ/mol, |ΔHmix| at x ~ 0.5
# Relative experimental ΔHmix uncertainties among systems passing the other
# filters are 0.5-3.5 %. Larger values are either flat 25 % uncertainties
# assigned during dataset curation or come from H^E close to zero.
MAX_RELATIVE_SIGMA_HMIX = 0.15
# Uncertainty floors used when weighting fits: 2 % of the value for the
# two-phase test, 0.001 in absolute units for the Redlich-Kister fits.
RELATIVE_SIGMA_FLOOR = 0.02
FIT_SIGMA_FLOOR = 1e-3

# Box planning
MIN_EDGE = 25.0  # Angstrom
MIN_MOLECULES = 100
MIN_MINORITY = 20
MAX_X_SHIFT = 0.01
DEFAULT_TARGET_X = (0.25, 0.5, 0.75)
MAX_TARGET_DISTANCE = 0.1  # a measured point must lie this close to each target

PURE_DENSITY_TOLERANCE = 0.002  # g/mL, cross-check flag
MAX_DENSITY_SIGMA = 0.002  # g/mL; noisier mixture densities are not used for V^E
STATE_T_TOLERANCE = 0.7  # K, window around a system's target temperature
MLPEG_T_TOLERANCE = 1.0  # K, for using an ml-peg pure density
PRESSURE_RANGE_KPA = (90.0, 110.0)


@dataclass
class Point:
    """
    One experimental point for a binary system, in the system's orientation.

    Attributes
    ----------
    property_name
        ``hmix``, ``density`` or ``vex``.
    x
        Mole fraction of component a.
    temperature
        Temperature in K.
    value
        Value in kJ/mol, g/mL or cm3/mol.
    sigma
        Standard uncertainty, or ``None``.
    doi
        Source DOI.
    source
        Where the point came from (data-set name, ``thermoml``, ``derived``,
        ``reference`` or ``ml-peg``).
    """

    property_name: str
    x: float
    temperature: float
    value: float
    sigma: float | None
    doi: str
    source: str


@dataclass
class RedlichKisterFit:
    """
    Redlich-Kister expansion y = x(1-x) sum_k a_k (2x-1)^k.

    Attributes
    ----------
    coefficients
        The fitted a_k.
    """

    coefficients: np.ndarray

    def __call__(self, x: float) -> float:
        """
        Evaluate the fitted curve.

        Parameters
        ----------
        x
            Mole fraction of component a.

        Returns
        -------
        float
            Fitted excess property.
        """
        return float(
            sum(
                coefficient * x * (1 - x) * (2 * x - 1) ** k
                for k, coefficient in enumerate(self.coefficients)
            )
        )


def fix_smiles(smiles: str) -> str:
    """
    Canonicalise SMILES and repair the known Sage tautomer error.

    Parameters
    ----------
    smiles
        SMILES from a Sage data set.

    Returns
    -------
    str
        Canonical SMILES of the intended molecule.
    """
    return canonical_smiles(TAUTOMER_FIXES.get(smiles, smiles))


def structure_problems(smiles: str) -> list[str]:
    """
    Check a component for charges and unlikely tautomers.

    Parameters
    ----------
    smiles
        Molecule SMILES.

    Returns
    -------
    list[str]
        Problems found (empty if the structure passes).
    """
    molecule = Chem.MolFromSmiles(smiles)
    problems = []
    if any(atom.GetFormalCharge() for atom in molecule.GetAtoms()):
        problems.append("non-zero formal charge")
    if molecule.HasSubstructMatch(ENOL_OR_IMIDIC_ACID):
        problems.append("enol or imidic acid tautomer")
    stored = Chem.MolToSmiles(molecule)
    canonical = Chem.MolToSmiles(TAUTOMER_ENUMERATOR.Canonicalize(molecule))
    if canonical != stored and stored not in TAUTOMER_REVIEWED:
        problems.append(
            f"RDKit canonical tautomer is {canonical} (review, then add "
            "to TAUTOMER_REVIEWED if correct)"
        )
    return problems


def check_structures(systems: list[dict]) -> None:
    """
    Raise if any curated component fails the structure checks.

    Parameters
    ----------
    systems
        ``systems`` entries from ``selection.yml``.

    Raises
    ------
    ValueError
        Listing every failing component.
    """
    components = sorted({system[key] for system in systems for key in ("a", "b")})
    failures = []
    for smiles in components:
        problems = structure_problems(smiles)
        if problems:
            failures.append(f"{smiles}: {'; '.join(problems)}")
    if failures:
        raise ValueError("Structure checks failed:\n" + "\n".join(failures))


def molecule_info(smiles: str) -> dict:
    """
    Compute composition and size information for a molecule.

    Parameters
    ----------
    smiles
        Molecule SMILES.

    Returns
    -------
    dict
        Heavy-atom count, element set, molar mass (g/mol), total atoms and
        the names of any matched viscous-liquid patterns.
    """
    molecule = Chem.MolFromSmiles(smiles)
    with_hydrogens = Chem.AddHs(molecule)
    viscous = [
        name
        for name, smarts in VISCOUS_SMARTS.items()
        if molecule.HasSubstructMatch(Chem.MolFromSmarts(smarts))
    ]
    return {
        "heavy_atoms": molecule.GetNumHeavyAtoms(),
        "elements": {atom.GetSymbol() for atom in with_hydrogens.GetAtoms()},
        "molar_mass": Descriptors.MolWt(molecule),
        "n_atoms": with_hydrogens.GetNumAtoms(),
        "viscous": viscous,
    }


def parse_sage_row(row: dict, source: str) -> tuple[tuple[str, str], Point]:
    """
    Parse one row of a Sage physical-property data set.

    Parameters
    ----------
    row
        CSV row in the OpenFF data-set format.
    source
        Name of the data set, stored on the point.

    Returns
    -------
    tuple[tuple[str, str], Point]
        The sorted SMILES pair (``(smiles, "")`` for a pure liquid) and the
        point, with x of the first SMILES of the pair.
    """
    first = fix_smiles(row["Component 1"])
    second = fix_smiles(row["Component 2"]) if row["Component 2"] else ""
    x_first = float(row["Mole Fraction 1"]) if row["Mole Fraction 1"] else 1.0
    if second and second < first:
        first, second, x_first = second, first, 1.0 - x_first
    hmix = row.get("EnthalpyOfMixing Value (kJ / mol)")
    if hmix:
        property_name, value = "hmix", float(hmix)
        sigma = row.get("EnthalpyOfMixing Uncertainty (kJ / mol)")
    else:
        property_name = "density"
        value = float(row["Density Value (g / ml)"])
        sigma = row.get("Density Uncertainty (g / ml)")
    point = Point(
        property_name=property_name,
        x=x_first,
        temperature=float(row["Temperature (K)"]),
        value=value,
        sigma=float(sigma) if sigma else None,
        doi=row["Source"],
        source=source,
    )
    return (first, second), point


def load_points(cache: Path) -> dict[tuple[str, str], list[Point]]:
    """
    Load binary and pure points from the Sage sets, keyed by sorted SMILES pair.

    Pure-liquid densities are stored under the key ``(smiles, "")``. Mole
    fractions refer to the first SMILES of the key.

    Parameters
    ----------
    cache
        Cache directory for downloads.

    Returns
    -------
    dict[tuple[str, str], list[Point]]
        Points per system.
    """
    points = defaultdict(list)
    for name, url in SOURCES.items():
        path = download(url, cache / f"{name}.csv")
        with path.open() as handle:
            for row in csv.DictReader(handle):
                key, point = parse_sage_row(row, name)
                points[key].append(point)
    return points


def load_sage_benchmarks(cache: Path) -> dict[tuple[str, str], list[dict]]:
    """
    Load published Sage benchmark results (3-replicate means and SDs).

    Parameters
    ----------
    cache
        Cache directory for downloads.

    Returns
    -------
    dict[tuple[str, str], list[dict]]
        Per sorted SMILES pair: force field, property, value, reference value
        and replicate SD for each benchmarked point.
    """
    results = defaultdict(list)
    for name, url in SAGE_BENCHMARKS.items():
        property_name = name.split("-")[-1]
        path = download(url, cache / f"{name}.csv")
        with path.open() as handle:
            for row in csv.DictReader(handle):
                if (
                    row["forcefield"] not in (*SAGE_FORCE_FIELDS, SAGE_RC2)
                    or not row["smiles_2"]
                ):
                    continue
                key = tuple(
                    sorted((fix_smiles(row["smiles_1"]), fix_smiles(row["smiles_2"])))
                )
                results[key].append(
                    {
                        "forcefield": row["forcefield"],
                        "property_name": property_name,
                        "value": float(row["value"]),
                        "reference": float(row["reference_value"]),
                        "sd": float(row["uncertainty"]) if row["uncertainty"] else None,
                        "substance": row["substance"],
                        "tautomer_error": "OC1=NCCC1"
                        in (row["smiles_1"], row["smiles_2"]),
                    }
                )
    return results


def load_mlpeg_liquids() -> dict[str, dict]:
    """
    Read ml-peg's liquid_densities pure liquids and their reference densities.

    Returns
    -------
    dict[str, dict]
        Canonical SMILES to name, reference density (g/mL) and temperature (K).
    """
    names = yaml.safe_load((HERE / "mlpeg_liquids.yml").read_text())
    data_dir = (
        download_s3_data(
            filename="liquid_densities.zip",
            key="inputs/molecular_dynamics/liquid_densities/liquid_densities.zip",
        )
        / "liquid_densities"
        / "equilibrated_structures_xyz"
    )
    liquids = {}
    for stem, smiles in names.items():
        atoms = read(data_dir / f"{stem}.xyz")
        liquids[canonical_smiles(smiles)] = {
            "name": stem,
            "density": float(atoms.info["exp_density"]),
            "temperature": float(atoms.info["exp_temperature"]),
        }
    return liquids


def interior_hmix_arrays(points: list[Point]) -> tuple[np.ndarray, ...]:
    """
    Return sorted x, ΔHmix and floored σ of the points away from the pure ends.

    Parameters
    ----------
    points
        Enthalpy of mixing points of one system at one temperature.

    Returns
    -------
    tuple[np.ndarray, ...]
        The x, value and σ arrays, where σ has a 2 % relative floor so that tiny
        reported uncertainties do not dominate.
    """
    interior = sorted((p for p in points if 0.05 < p.x < 0.95), key=lambda p: p.x)
    x = np.array([p.x for p in interior])
    value = np.array([p.value for p in interior])
    sigma = np.array(
        [
            max(p.sigma or 0.0, RELATIVE_SIGMA_FLOOR * abs(p.value), 1e-3)
            for p in interior
        ]
    )
    return x, value, sigma


def linear_fit_chi2(points: list[Point]) -> float | None:
    """
    Reduced chi^2 of a straight-line fit to interior excess-enthalpy points.

    A single-phase H^E curve gives a large value; a two-phase region, where
    H^E = a + b x, gives a value near 1.

    Parameters
    ----------
    points
        Enthalpy of mixing points of one system at one temperature.

    Returns
    -------
    float | None
        Reduced chi^2, or ``None`` with fewer than four interior points.
    """
    x, value, sigma = interior_hmix_arrays(points)
    if len(x) < 4:
        return None
    coefficients = np.polyfit(x, value, 1, w=1 / sigma)
    residuals = (np.polyval(coefficients, x) - value) / sigma
    return float(np.sum(residuals**2)) / (len(x) - 2)


def two_phase_signature(points: list[Point]) -> bool:
    """
    Test whether excess enthalpies are linear in x (lever-rule signature).

    In a two-phase region H^E = a + b x. With full isotherms this fits a
    straight line through the interior points and flags reduced chi^2 <= 2.
    With exactly three points it compares the two consecutive slopes against
    their propagated uncertainty.

    An ideal or near-ideal mixture (H^E close to zero everywhere) is also a
    straight line, so the test is only applied when the largest measured
    |H^E| reaches ``MIN_SIGNAL``. The largest value rather than the x = 0.5
    value is used because a two-phase line can be small at x = 0.5 and large
    at one end (acetophenone + water: 0.34 kJ/mol at x = 0.5, 0.54 at 0.8).

    Parameters
    ----------
    points
        Enthalpy of mixing points of one system at one temperature.

    Returns
    -------
    bool
        ``True`` if the data are consistent with a straight line.
    """
    x, value, sigma = interior_hmix_arrays(points)
    if len(x) < 3 or np.max(np.abs(value)) < MIN_SIGNAL:
        return False
    if len(x) == 3:
        slopes = np.diff(value) / np.diff(x)
        slope_variances = [
            (sigma[i] ** 2 + sigma[i + 1] ** 2) / (x[i + 1] - x[i]) ** 2
            for i in range(2)
        ]
        return bool(abs(slopes[0] - slopes[1]) < 2 * math.sqrt(sum(slope_variances)))
    return linear_fit_chi2(points) <= 2.0


def max_composition_mismatch(hmix: list[Point], density: list[Point]) -> float:
    """
    Largest distance from an enthalpy composition to the nearest density one.

    Parameters
    ----------
    hmix
        Enthalpy of mixing points.
    density
        Mixture density points.

    Returns
    -------
    float
        Maximum |delta x|.
    """
    return max(min(abs(h.x - d.x) for d in density) for h in hmix)


def largest_sage_sd(sage_rows: list[dict], property_name: str) -> float | None:
    """
    Return the largest published Sage 3-replicate SD of one property.

    Parameters
    ----------
    sage_rows
        Published Sage results of one system.
    property_name
        ``hmix`` or ``density``.

    Returns
    -------
    float | None
        The largest SD, or ``None`` if none is published.
    """
    return max(
        (
            row["sd"]
            for row in sage_rows
            if row["property_name"] == property_name
            and row["sd"]
            and row["forcefield"] in SAGE_FORCE_FIELDS
        ),
        default=None,
    )


def filter_failures(row: dict, info_a: dict, info_b: dict) -> list[str]:
    """
    List the hard filters a screened system fails.

    Parameters
    ----------
    row
        Screening metrics from :func:`screen_candidates`.
    info_a, info_b
        Output of :func:`molecule_info` for the two components.

    Returns
    -------
    list[str]
        Failed filters, in the order of the report's filter table.
    """
    failures = []
    if row["max_dx"] > MAX_COMPOSITION_MISMATCH:
        failures.append("composition mismatch")
    if not (row["exp_sigma_hmix"] and row["exp_sigma_density"]):
        failures.append("missing experimental uncertainty")
    if (row["hmix_max_rel_sigma"] or 0.0) > MAX_RELATIVE_SIGMA_HMIX:
        failures.append(
            f"relative ΔHmix uncertainty > {100 * MAX_RELATIVE_SIGMA_HMIX:.0f} %"
        )
    if not row["pure_density_sigma"]:
        failures.append("missing pure density with uncertainty")
    if not (info_a["elements"] | info_b["elements"]) <= ALLOWED_ELEMENTS:
        failures.append("elements")
    if max(info_a["heavy_atoms"], info_b["heavy_atoms"]) > MAX_HEAVY_ATOMS:
        failures.append("too many heavy atoms")
    if row["viscous"]:
        failures.append("viscous")
    if row["linear_two_phase"]:
        failures.append("two-phase (linear H^E)")
    if row["ucst_screen"]:
        failures.append("two-phase (H^E > RT/2)")
    for property_name in ("hmix", "density"):
        sd = row[f"sage_rep_sd_{property_name}"]
        if sd is None:
            failures.append(f"no Sage repeat SD ({property_name})")
        elif sd > MAX_REPEAT_SD[property_name]:
            failures.append(f"Sage repeat SD ({property_name})")
    if abs(row["hmix_eq"]) < MIN_SIGNAL:
        failures.append("small signal")
    return failures


def in_screen_window(point: Point) -> bool:
    """
    Test whether a point lies in the screening temperature window.

    Parameters
    ----------
    point
        Experimental point.

    Returns
    -------
    bool
        ``True`` inside ``SCREEN_T_RANGE``.
    """
    return SCREEN_T_RANGE[0] <= point.temperature <= SCREEN_T_RANGE[1]


def screen_candidates(
    points: dict[tuple[str, str], list[Point]],
    sage: dict[tuple[str, str], list[dict]],
    mlpeg: dict[str, dict],
) -> list[dict]:
    """
    Apply the hard filters to every binary system with enough data.

    Parameters
    ----------
    points
        Experimental points from :func:`load_points`.
    sage
        Published Sage benchmark results.
    mlpeg
        The ml-peg pure liquids.

    Returns
    -------
    list[dict]
        One row per candidate system with every metric and filter outcome.
    """
    pure_in_window = {
        key[0]: point
        for key, system_points in points.items()
        if not key[1]
        for point in system_points
        if in_screen_window(point)
    }
    rows = []
    for (a, b), system_points in points.items():
        if not b:
            continue
        window = [p for p in system_points if in_screen_window(p)]
        # One point per distinct composition (the last one listed wins).
        hmix = list(
            {round(p.x, 3): p for p in window if p.property_name == "hmix"}.values()
        )
        density = list(
            {round(p.x, 3): p for p in window if p.property_name == "density"}.values()
        )
        if len(hmix) < 3 or len(density) < 3:
            continue
        info_a, info_b = molecule_info(a), molecule_info(b)
        hmix_equimolar = min(hmix, key=lambda p: abs(p.x - 0.5)).value
        hmix_by_x = sorted(hmix, key=lambda p: p.x)
        sage_rows = [r for r in sage.get((a, b), []) if not r["tautomer_error"]]
        sage_221_errors = [
            abs(r["value"] - r["reference"])
            for r in sage_rows
            if r["forcefield"] == "openff-2.2.1" and r["property_name"] == "hmix"
        ]
        row = {
            "a": a,
            "b": b,
            "sources": "+".join(sorted({p.source for p in window})),
            "heavy_atoms": f"{info_a['heavy_atoms']}/{info_b['heavy_atoms']}",
            "elements": "".join(sorted(info_a["elements"] | info_b["elements"])),
            "n_hmix": len(hmix),
            "n_density": len(density),
            "hmix_x": " ".join(f"{p.x:.3f}" for p in hmix_by_x),
            "hmix_values": " ".join(f"{p.value:.4f}" for p in hmix_by_x),
            "hmix_eq": round(hmix_equimolar, 4),
            "max_dx": round(max_composition_mismatch(hmix, density), 3),
            "exp_sigma_hmix": all(p.sigma for p in hmix),
            "hmix_max_rel_sigma": max(
                (p.sigma / abs(p.value) for p in hmix if p.sigma and p.value),
                default=None,
            ),
            "exp_sigma_density": all(p.sigma for p in density),
            "pure_density_sigma": all(
                c in pure_in_window and pure_in_window[c].sigma for c in (a, b)
            ),
            "in_mlpeg": "+".join(c for c in (a, b) if c in mlpeg),
            "n_in_mlpeg": sum(c in mlpeg for c in (a, b)),
            "sage_rep_sd_hmix": largest_sage_sd(sage_rows, "hmix"),
            "sage_rep_sd_density": largest_sage_sd(sage_rows, "density"),
            "sage_2.2.1_hmix_mae": round(statistics.mean(sage_221_errors), 4)
            if sage_221_errors
            else None,
            "viscous": "+".join(info_a["viscous"] + info_b["viscous"]),
            "linear_two_phase": two_phase_signature(hmix),
            "ucst_screen": hmix_equimolar > UCST_LIMIT,
        }
        failures = filter_failures(row, info_a, info_b)
        row["failures"] = "; ".join(failures)
        row["passes"] = not failures
        rows.append(row)
    return sorted(
        rows, key=lambda r: (not r["passes"], r["n_in_mlpeg"], -abs(r["hmix_eq"]))
    )


def thermoml_point(measurement: Measurement, x: float) -> Point:
    """
    Convert a ThermoML measurement into a point at mole fraction ``x`` of a.

    Parameters
    ----------
    measurement
        Parsed ThermoML measurement.
    x
        Mole fraction of component a (1.0 for a pure liquid).

    Returns
    -------
    Point
        The point, with ``source`` "thermoml".
    """
    return Point(
        measurement.property_name,
        x,
        measurement.temperature,
        measurement.value,
        measurement.uncertainty,
        measurement.doi,
        "thermoml",
    )


def thermoml_points(
    a: str,
    b: str,
    dois: set[str],
    temperature: float,
    window: float,
    cache: Path,
    corrections: dict[str, float] | None = None,
) -> tuple[list[Point], dict[str, list[Point]], dict[str, str]]:
    """
    Collect full isotherms and pure densities from ThermoML for one system.

    Parameters
    ----------
    a, b
        Canonical SMILES of the two components; x refers to ``a``.
    dois
        Source DOIs to fetch.
    temperature
        Target temperature in K.
    window
        Accepted |T - target| in K.
    cache
        Cache directory for ThermoML records.
    corrections
        DOIs whose recorded mixture temperatures are wrong, mapped to the
        correct temperature. Pure densities from those records are skipped.

    Returns
    -------
    tuple[list[Point], dict[str, list[Point]], dict[str, str]]
        Mixture points, pure-density points per component, and any DOI that
        could not be fetched or parsed mapped to the error message.
    """
    mixture, pure, errors = [], defaultdict(list), {}
    for doi in sorted(dois):
        try:
            measurements = parse_record(doi, fetch_record(doi, cache / "thermoml"))
        except Exception as error:  # noqa: BLE001 - report and continue
            errors[doi] = str(error)
            continue
        corrected_temperature = (corrections or {}).get(doi)
        for measurement in measurements:
            if corrected_temperature is not None:
                is_pure = (
                    len(measurement.smiles) == 1 or 1.0 in measurement.mole_fractions
                )
                if is_pure:
                    continue
                measurement.temperature = corrected_temperature
            if abs(measurement.temperature - temperature) > window:
                continue
            pressure = measurement.pressure
            if pressure is not None and not (
                PRESSURE_RANGE_KPA[0] <= pressure <= PRESSURE_RANGE_KPA[1]
            ):
                continue
            is_density = measurement.property_name == "density"
            if set(measurement.smiles) == {a, b}:
                fractions = dict(
                    zip(measurement.smiles, measurement.mole_fractions, strict=True)
                )
                x = fractions[a]
                if x in (0.0, 1.0) and is_density:
                    pure[a if x == 1.0 else b].append(thermoml_point(measurement, 1.0))
                elif 0.0 < x < 1.0:
                    mixture.append(thermoml_point(measurement, x))
            elif len(measurement.smiles) == 1 and measurement.smiles[0] in (a, b):
                if is_density:
                    pure[measurement.smiles[0]].append(thermoml_point(measurement, 1.0))
    return mixture, dict(pure), errors


def redlich_kister(x: np.ndarray, y: np.ndarray, sigma: np.ndarray) -> RedlichKisterFit:
    """
    Fit a weighted Redlich-Kister expansion y = x(1-x) sum_k a_k (2x-1)^k.

    The number of terms is min(4, n - 1) so that at least one degree of
    freedom remains when n >= 2.

    Parameters
    ----------
    x
        Mole fractions.
    y
        Excess property values.
    sigma
        Standard uncertainties used as weights.

    Returns
    -------
    RedlichKisterFit
        The fitted curve.
    """
    n_terms = max(1, min(4, len(x) - 1))
    basis = np.stack([x * (1 - x) * (2 * x - 1) ** k for k in range(n_terms)], axis=1)
    weights = 1.0 / sigma
    coefficients, *_ = np.linalg.lstsq(
        basis * weights[:, None], y * weights, rcond=None
    )
    return RedlichKisterFit(coefficients)


def fit_points(points: list[Point]) -> RedlichKisterFit | None:
    """
    Fit a Redlich-Kister curve to excess-property points.

    Parameters
    ----------
    points
        Points to fit.

    Returns
    -------
    RedlichKisterFit | None
        Fitted curve, or ``None`` if there are no points.
    """
    if not points:
        return None
    x = np.array([p.x for p in points])
    y = np.array([p.value for p in points])
    sigma = np.array([max(p.sigma or 0.0, FIT_SIGMA_FLOOR) for p in points])
    return redlich_kister(x, y, sigma)


def plan_box(x: float, density: float, info_a: dict, info_b: dict) -> dict:
    """
    Choose molecule counts for a box (minimum edge, molecule and minority counts).

    Parameters
    ----------
    x
        Target mole fraction of component a (1.0 for a pure box of a).
    density
        Expected density in g/mL.
    info_a, info_b
        Output of :func:`molecule_info` for the two components.

    Returns
    -------
    dict
        Molecule counts, achieved x, atom count and expected edge length.
    """
    mean_mass = x * info_a["molar_mass"] + (1 - x) * info_b["molar_mass"]
    volume_cm3 = (MIN_EDGE * 1e-8) ** 3
    n_total = max(MIN_MOLECULES, math.ceil(density * volume_cm3 * AVOGADRO / mean_mass))
    is_pure = x in (0.0, 1.0)
    # Add molecules until x is close enough and the minority is large enough.
    while True:
        n_a = round(x * n_total)
        n_b = n_total - n_a
        if is_pure or (
            abs(n_a / n_total - x) <= MAX_X_SHIFT and min(n_a, n_b) >= MIN_MINORITY
        ):
            break
        n_total += 1
    mass = (n_a * info_a["molar_mass"] + n_b * info_b["molar_mass"]) / AVOGADRO
    edge = (mass / density) ** (1 / 3) * 1e8
    return {
        "n_a": n_a,
        "n_b": n_b,
        "x_box": n_a / n_total,
        "n_atoms": n_a * info_a["n_atoms"] + n_b * info_b["n_atoms"],
        "edge": edge,
    }


def choose_pure_density(
    smiles: str,
    temperature: float,
    thermoml_pure: list[Point],
    sage_pure: list[Point],
    mlpeg: dict[str, dict],
    preferred_dois: list[str] | None = None,
    reference: dict | None = None,
) -> dict:
    """
    Pick the reference pure density and cross-check it against other sources.

    Preference order: an explicit reference from ``selection.yml``; the value
    from the paper that supplies most of the mixture densities (so V^E is
    self-consistent); the median ThermoML value; rc2/Sage data-set values;
    ml-peg. The ml-peg value is used only when its temperature is within 1 K.

    Parameters
    ----------
    smiles
        Component SMILES.
    temperature
        Target temperature in K.
    thermoml_pure
        Pure densities from the fetched ThermoML records.
    sage_pure
        Pure densities from the Sage/rc2 data sets.
    mlpeg
        The ml-peg pure liquids.
    preferred_dois
        Mixture-density DOIs, most points first.
    reference
        Explicit reference value (keys ``density``, ``sigma``, ``source``).

    Returns
    -------
    dict
        Chosen value, uncertainty and source, all cross-check values and the
        sources that disagree by more than ``PURE_DENSITY_TOLERANCE``.
    """
    near_thermoml = [
        p
        for p in thermoml_pure
        if abs(p.temperature - temperature) <= STATE_T_TOLERANCE
    ]
    near_sage = [
        p for p in sage_pure if abs(p.temperature - temperature) <= STATE_T_TOLERANCE
    ]
    checks = {f"thermoml:{p.doi}": p.value for p in near_thermoml}
    checks.update({f"{p.source}:{p.doi}": p.value for p in near_sage})
    mlpeg_liquid = mlpeg.get(smiles)
    mlpeg_usable = (
        mlpeg_liquid is not None
        and abs(mlpeg_liquid["temperature"] - temperature) <= MLPEG_T_TOLERANCE
    )
    if mlpeg_usable:
        checks["ml-peg"] = mlpeg_liquid["density"]
    same_paper = [
        p for doi in (preferred_dois or []) for p in near_thermoml if p.doi == doi
    ]
    if reference:
        chosen = Point(
            "density",
            1.0,
            temperature,
            reference["density"],
            reference["sigma"],
            reference["source"],
            "reference",
        )
        checks["reference"] = reference["density"]
    elif same_paper:
        chosen = same_paper[0]
    elif near_thermoml:
        median = statistics.median(p.value for p in near_thermoml)
        chosen = min(near_thermoml, key=lambda p: abs(p.value - median))
    elif near_sage:
        chosen = min(near_sage, key=lambda p: p.sigma or 1.0)
    elif mlpeg_usable:
        chosen = Point(
            "density",
            1.0,
            mlpeg_liquid["temperature"],
            mlpeg_liquid["density"],
            None,
            "ml-peg",
            "ml-peg",
        )
    else:
        chosen = None

    values = list(checks.values())
    spread = max(values) - min(values) if len(values) > 1 else 0.0
    if chosen is None:
        outliers, source = [], ""
    else:
        outliers = [
            name
            for name, value in checks.items()
            if abs(value - chosen.value) > PURE_DENSITY_TOLERANCE
        ]
        is_reference = chosen.source == "reference"
        source = chosen.doi if is_reference else f"{chosen.source}:{chosen.doi}"
    return {
        "density": chosen.value if chosen else None,
        "sigma": chosen.sigma if chosen else None,
        "source": source,
        "checks": checks,
        "spread": spread,
        "outliers": outliers,
        "mlpeg_density": mlpeg_liquid["density"] if mlpeg_liquid else None,
        "mlpeg_temperature": mlpeg_liquid["temperature"] if mlpeg_liquid else None,
    }


def excess_volumes(
    mixture_densities: list[Point],
    thermoml_pure: dict[str, list[Point]],
    chosen_pure: dict[str, dict],
    a: str,
    b: str,
    molar_mass_a: float,
    molar_mass_b: float,
) -> list[Point]:
    """
    Derive V^E from mixture densities, with correlated uncertainties.

    V^E = M/rho - x M_a/rho_a - (1 - x) M_b/rho_b. Each point uses the pure
    densities of its own paper (closest temperature in the state window)
    where it reports them, so calibration and purity errors common to the
    paper largely cancel; otherwise it falls back to the chosen reference
    pure density.

    ThermoML gives one combined uncertainty with no split into systematic and
    random parts. All of it is treated as systematic: one error shared by the
    densities from the same paper. Its V^E uncertainty is the larger of two
    shared errors, a common scale error of the largest relative sigma and a
    common offset of the largest sigma among the densities used. The random
    part is added from the scatter about the fit by :func:`add_fit_scatter`.
    A pure density from another source is uncorrelated with the mixture and
    adds in quadrature.

    Parameters
    ----------
    mixture_densities
        Mixture densities in g/mL, x of component a.
    thermoml_pure
        Pure densities per component from the fetched ThermoML records.
    chosen_pure
        Chosen reference pure densities per component.
    a, b
        Component SMILES.
    molar_mass_a, molar_mass_b
        Molar masses in g/mol.

    Returns
    -------
    list[Point]
        V^E points in cm3/mol with their correlated standard uncertainty.
    """
    excess = []
    for point in mixture_densities:
        x = point.x
        # (density, sigma, from the same paper) for the mixture, a and b.
        densities = [(point.value, point.sigma, True)]
        for component in (a, b):
            same_paper = [
                pure
                for pure in thermoml_pure.get(component, [])
                if pure.doi == point.doi
                and abs(pure.temperature - point.temperature) <= STATE_T_TOLERANCE
            ]
            if same_paper:
                pure = min(
                    same_paper, key=lambda q: abs(q.temperature - point.temperature)
                )
                densities.append((pure.value, pure.sigma, True))
            else:
                chosen = chosen_pure[component]
                densities.append((chosen["density"], chosen["sigma"], False))
        (rho_mixture, _, _), (rho_a, _, _), (rho_b, _, _) = densities
        mean_mass = x * molar_mass_a + (1 - x) * molar_mass_b
        value = (
            mean_mass / rho_mixture
            - x * molar_mass_a / rho_a
            - (1 - x) * molar_mass_b / rho_b
        )
        # Sensitivities dV^E/drho for the mixture, a and b.
        sensitivities = [
            -mean_mass / rho_mixture**2,
            x * molar_mass_a / rho_a**2,
            (1 - x) * molar_mass_b / rho_b**2,
        ]
        # (sensitivity, density, sigma) of the densities sharing the paper's error.
        shared = [
            (sensitivity, rho, sigma_rho or 0.0)
            for sensitivity, (rho, sigma_rho, same) in zip(
                sensitivities, densities, strict=True
            )
            if same
        ]
        if not any(sigma_rho for _, _, sigma_rho in shared):
            sigma_vex = None
        else:
            # Common scale error and common offset, each as large as the
            # paper's largest relative or absolute sigma at this point.
            largest_relative_sigma = max(
                sigma_rho / rho for _, rho, sigma_rho in shared
            )
            largest_sigma = max(sigma_rho for _, _, sigma_rho in shared)
            scale_error = (
                abs(sum(sensitivity * rho for sensitivity, rho, _ in shared))
                * largest_relative_sigma
            )
            offset_error = (
                abs(sum(sensitivity for sensitivity, _, _ in shared)) * largest_sigma
            )
            shared_error = max(scale_error, offset_error)
            independent_errors = [
                sensitivity * (sigma_rho or 0.0)
                for sensitivity, (_, sigma_rho, same) in zip(
                    sensitivities, densities, strict=True
                )
                if not same
            ]
            sigma_vex = math.sqrt(
                shared_error**2 + sum(error**2 for error in independent_errors)
            )
        excess.append(
            Point("vex", x, point.temperature, value, sigma_vex, point.doi, "derived")
        )
    return excess


def add_fit_scatter(vex: list[Point], curve: RedlichKisterFit) -> list[Point]:
    """
    Add the random part of the derived V^E uncertainty from the fit residuals.

    The random scatter of each paper's derived points is estimated as its
    RMS residual about the Redlich-Kister fit, corrected for the degrees of
    freedom the fit uses, and added in quadrature to the correlated part.
    Directly measured V^E points keep their reported uncertainty.

    Parameters
    ----------
    vex
        V^E points, derived and direct.
    curve
        Redlich-Kister fit to ``vex``.

    Returns
    -------
    list[Point]
        The points with updated uncertainties for derived points.
    """
    n_points = len(vex)
    degrees_of_freedom = n_points - max(1, min(4, n_points - 1))
    if degrees_of_freedom <= 0:
        return vex
    squared_residuals = defaultdict(list)
    for p in vex:
        if p.source == "derived":
            squared_residuals[p.doi].append((p.value - curve(p.x)) ** 2)
    scatter = {
        doi: math.sqrt(sum(squares) / len(squares) * n_points / degrees_of_freedom)
        for doi, squares in squared_residuals.items()
    }
    updated = []
    for p in vex:
        if p.source == "derived":
            sigma = math.sqrt((p.sigma or 0.0) ** 2 + scatter[p.doi] ** 2)
            p = Point(
                p.property_name, p.x, p.temperature, p.value, sigma, p.doi, p.source
            )
        updated.append(p)
    return updated


def orient_to_component_a(point: Point, flip: bool) -> Point:
    """
    Re-express a Sage-set point (keyed by the sorted pair) with x of component a.

    Parameters
    ----------
    point
        Point whose x refers to the first SMILES of the sorted pair.
    flip
        ``True`` if component a is the second SMILES of the sorted pair.

    Returns
    -------
    Point
        Point with x of component a.
    """
    return Point(
        point.property_name,
        1 - point.x if flip else point.x,
        point.temperature,
        point.value,
        point.sigma,
        point.doi,
        point.source,
    )


def find_reference(
    references: list[dict], smiles: str, temperature: float
) -> dict | None:
    """
    Find an explicit pure-density reference for one component in ``selection.yml``.

    Parameters
    ----------
    references
        ``pure_density_references`` from ``selection.yml``.
    smiles
        Canonical SMILES.
    temperature
        Temperature in K.

    Returns
    -------
    dict | None
        The matching reference.
    """
    for reference in references:
        if (
            canonical_smiles(reference["smiles"]) == smiles
            and abs(reference["temperature"] - temperature) < 0.1
        ):
            return reference
    return None


def target_points(hmix: list[Point], target_x) -> list[Point]:
    """
    Pick the measured ΔHmix points nearest each target composition.

    Parameters
    ----------
    hmix
        ΔHmix points sorted by x.
    target_x
        Target compositions of component a.

    Returns
    -------
    list[Point]
        One point per target that has a measurement within
        ``MAX_TARGET_DISTANCE``, without duplicates.
    """
    targets = []
    for goal in target_x:
        if hmix:
            nearest = min(hmix, key=lambda p: abs(p.x - goal))
            if abs(nearest.x - goal) <= MAX_TARGET_DISTANCE and nearest not in targets:
                targets.append(nearest)
    return targets


def box_row(
    spec: dict,
    target: Point,
    context: dict,
) -> dict:
    """
    Build one ``systems.csv`` row: box plan and experimental values at box x.

    Parameters
    ----------
    spec
        Entry from ``selection.yml``.
    target
        The measured ΔHmix point the box is built for.
    context
        The system's components, molecule info, pure densities, fits and
        mixture densities, as assembled by :func:`build_system`.

    Returns
    -------
    dict
        The row.
    """
    a, b = context["a"], context["b"]
    info, pures = context["info"], context["pures"]
    hmix_fit, vex_fit = context["hmix_fit"], context["vex_fit"]
    densities, density_dois = context["densities"], context["density_dois"]
    molar_mass_a, molar_mass_b = info[a]["molar_mass"], info[b]["molar_mass"]
    have_pure_densities = bool(pures[a]["density"] and pures[b]["density"])

    density_guess = (
        pures[a]["density"] * target.x + pures[b]["density"] * (1 - target.x)
        if have_pure_densities
        else 1.0
    )
    box = plan_box(target.x, density_guess, info[a], info[b])
    x_box = box["x_box"]
    row = {
        "system_id": spec["id"],
        "a": a,
        "b": b,
        "temperature": context["temperature"],
        "x_target": round(target.x, 4),
        "x_box": round(x_box, 4),
        "n_a": box["n_a"],
        "n_b": box["n_b"],
        "n_atoms": box["n_atoms"],
        "hmix_measured": target.value,
        "hmix_measured_sigma": target.sigma,
        "hmix_measured_T": target.temperature,
        # The measured value, moved from the measured x to the box x along the fit.
        "hmix_exp": round(hmix_fit(x_box) + (target.value - hmix_fit(target.x)), 5),
        "hmix_sigma": target.sigma,
        "hmix_doi": target.doi,
    }
    if not (vex_fit and have_pure_densities):
        row.update(
            vex_exp=None,
            density_exp=None,
            density_sigma=None,
            density_doi="",
            density_nearest_dx=None,
            edge=round(box["edge"], 2),
        )
        return row

    mean_mass = x_box * molar_mass_a + (1 - x_box) * molar_mass_b
    ideal_volume = (
        x_box * molar_mass_a / pures[a]["density"]
        + (1 - x_box) * molar_mass_b / pures[b]["density"]
    )
    vex_box = vex_fit(x_box)
    row["vex_exp"] = round(vex_box, 5)
    row["density_exp"] = round(mean_mass / (ideal_volume + vex_box), 5)
    nearest = min(densities, key=lambda p: abs(p.x - x_box)) if densities else None
    # Uncertainty from the main density paper near x_box (other papers
    # sometimes report implausible uncertainties).
    local_sigmas = [
        p.sigma
        for p in densities
        if abs(p.x - x_box) <= 0.05 and p.sigma and p.doi == density_dois[0]
    ] or [p.sigma for p in densities if abs(p.x - x_box) <= 0.05 and p.sigma]
    if local_sigmas:
        row["density_sigma"] = statistics.median(local_sigmas)
    else:
        row["density_sigma"] = nearest.sigma if nearest else None
    row["density_doi"] = nearest.doi if nearest else ""
    row["density_nearest_dx"] = round(abs(nearest.x - x_box), 3) if nearest else None
    row["edge"] = round(
        plan_box(target.x, row["density_exp"], info[a], info[b])["edge"], 2
    )
    return row


def build_system(
    spec: dict,
    points: dict[tuple[str, str], list[Point]],
    sage: dict[tuple[str, str], list[dict]],
    mlpeg: dict[str, dict],
    cache: Path,
    references: list[dict] | None = None,
) -> dict:
    """
    Assemble reference data, fits and box plans for one curated system.

    Parameters
    ----------
    spec
        Entry from ``selection.yml``.
    points
        Sage data-set points.
    sage
        Published Sage benchmark results.
    mlpeg
        The ml-peg pure liquids.
    cache
        Cache directory.
    references
        Explicit pure-density references from ``selection.yml``.

    Returns
    -------
    dict
        Everything the outputs and report need for this system.
    """
    references = references or []
    a, b = canonical_smiles(spec["a"]), canonical_smiles(spec["b"])
    temperature = float(spec["temperature"])
    key = tuple(sorted((a, b)))
    flip = key[0] != a

    # Fetch every paper that the Sage sets cite for this system or its pure liquids.
    set_points = [orient_to_component_a(p, flip) for p in points.get(key, [])]
    dois = {
        p.doi
        for p in set_points
        if abs(p.temperature - temperature) <= STATE_T_TOLERANCE
    }
    sage_pure = {c: points.get((c, ""), []) for c in (a, b)}
    dois |= {
        p.doi
        for c in (a, b)
        for p in sage_pure[c]
        if abs(p.temperature - temperature) <= STATE_T_TOLERANCE
    }
    dois |= set(spec.get("extra_dois", []))
    mixture, thermoml_pure, errors = thermoml_points(
        a,
        b,
        dois,
        temperature,
        STATE_T_TOLERANCE,
        cache,
        spec.get("temperature_corrections"),
    )

    hmix = sorted((p for p in mixture if p.property_name == "hmix"), key=lambda p: p.x)
    # Restrict ΔHmix to chosen papers where the sources disagree.
    if spec.get("hmix_dois"):
        hmix = [p for p in hmix if p.doi in spec["hmix_dois"]]
    all_densities = sorted(
        (p for p in mixture if p.property_name == "density"), key=lambda p: p.x
    )
    # Restrict densities to one paper, so its pure densities serve throughout.
    if spec.get("density_dois"):
        all_densities = [p for p in all_densities if p.doi in spec["density_dois"]]
    densities = [
        p for p in all_densities if p.sigma is None or p.sigma <= MAX_DENSITY_SIGMA
    ]
    noisy_density = Counter(p.doi for p in all_densities if p not in densities)
    density_dois = [doi for doi, _ in Counter(p.doi for p in densities).most_common()]

    info = {a: molecule_info(a), b: molecule_info(b)}
    pures = {
        c: choose_pure_density(
            c,
            temperature,
            thermoml_pure.get(c, []),
            sage_pure[c],
            mlpeg,
            density_dois,
            find_reference(references, c, temperature),
        )
        for c in (a, b)
    }

    vex = sorted((p for p in mixture if p.property_name == "vex"), key=lambda p: p.x)
    if pures[a]["density"] and pures[b]["density"]:
        vex += excess_volumes(
            densities,
            thermoml_pure,
            pures,
            a,
            b,
            info[a]["molar_mass"],
            info[b]["molar_mass"],
        )

    hmix_fit = fit_points(hmix)
    vex_fit = fit_points(vex)
    if vex_fit:
        # Refit once the random scatter is part of each point's σ.
        vex = add_fit_scatter(vex, vex_fit)
        vex_fit = fit_points(vex)

    context = {
        "a": a,
        "b": b,
        "temperature": temperature,
        "info": info,
        "pures": pures,
        "hmix_fit": hmix_fit,
        "vex_fit": vex_fit,
        "densities": densities,
        "density_dois": density_dois,
    }
    targets = target_points(hmix, spec.get("target_x", DEFAULT_TARGET_X))
    rows = [box_row(spec, target, context) for target in targets]

    sage_rows = [r for r in sage.get(key, []) if not r["tautomer_error"]]
    return {
        "spec": spec,
        "a": a,
        "b": b,
        "temperature": temperature,
        "info": info,
        "pures": pures,
        "hmix": hmix,
        "density": densities,
        "vex": vex,
        "hmix_fit": hmix_fit,
        "vex_fit": vex_fit,
        "rows": rows,
        "set_points": set_points,
        "sage": sage_rows,
        "errors": errors,
        "noisy_density": dict(noisy_density),
        "linear_chi2": linear_fit_chi2(hmix),
        "x_range": (hmix[0].x, hmix[-1].x) if hmix else None,
    }


def pure_rows(systems: list[dict], mlpeg: dict[str, dict]) -> list[dict]:
    """
    Build one ``pure.csv`` row per pure liquid and temperature.

    Parameters
    ----------
    systems
        Output of :func:`build_system` for every system.
    mlpeg
        The ml-peg pure liquids.

    Returns
    -------
    list[dict]
        Rows with the chosen density, its cross-checks and the box plan.
    """
    rows, seen = [], set()
    for system in systems:
        temperature = system["temperature"]
        for smiles in (system["a"], system["b"]):
            if (smiles, temperature) in seen:
                continue
            seen.add((smiles, temperature))
            pure = system["pures"][smiles]
            info = system["info"][smiles]
            box = plan_box(1.0, pure["density"] or 1.0, info, info)
            rows.append(
                {
                    "smiles": smiles,
                    "temperature": temperature,
                    "density": pure["density"],
                    "sigma": pure["sigma"],
                    "source": pure["source"],
                    "in_mlpeg": smiles in mlpeg,
                    # Liquids already in ml-peg's liquid_densities are run
                    # (ΔHmix needs them) but not scored twice.
                    "scored": smiles not in mlpeg,
                    "mlpeg_density": pure["mlpeg_density"],
                    "mlpeg_temperature": pure["mlpeg_temperature"],
                    "spread": round(pure["spread"], 5),
                    "outliers": "; ".join(pure["outliers"]),
                    "checks": "; ".join(
                        f"{source}={density:.5f}"
                        for source, density in pure["checks"].items()
                    ),
                    "n_molecules": box["n_a"],
                    "n_atoms": box["n_atoms"],
                    "edge": round(box["edge"], 2),
                }
            )
    return rows


def isotherm_rows(systems: list[dict]) -> list[dict]:
    """
    List every experimental point behind the fits, one row each.

    Parameters
    ----------
    systems
        Output of :func:`build_system` for every system.

    Returns
    -------
    list[dict]
        Rows of ``isotherms.csv``.
    """
    return [
        {
            "system_id": system["spec"]["id"],
            "prop": point.property_name,
            "x": point.x,
            "temperature": point.temperature,
            "value": point.value,
            "sigma": point.sigma,
            "doi": point.doi,
            "source": point.source,
        }
        for system in systems
        for point in system["hmix"] + system["density"] + system["vex"]
    ]


def write_csv(path: Path, rows: list[dict]) -> None:
    """
    Write a list of dictionaries to CSV, using the union of keys as header.

    Parameters
    ----------
    path
        Output path.
    rows
        Rows to write.
    """
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    """Run the selection and write all outputs."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output-dir", type=Path, default=Path("selection"))
    parser.add_argument("--cache-dir", type=Path, default=CACHE_DIR)
    parser.add_argument(
        "--fast-set",
        type=Path,
        default=None,
        help="Optional CSV of recommended fast-converging systems to assess.",
    )
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    selection = yaml.safe_load((HERE / "selection.yml").read_text())
    check_structures(selection["systems"])
    points = load_points(args.cache_dir)
    sage = load_sage_benchmarks(args.cache_dir)
    mlpeg = load_mlpeg_liquids()

    candidates = screen_candidates(points, sage, mlpeg)
    write_csv(out / "candidates.csv", candidates)

    references = selection.get("pure_density_references", [])
    systems = [
        build_system(spec, points, sage, mlpeg, args.cache_dir, references)
        for spec in selection["systems"]
    ]
    write_csv(
        out / "systems.csv", [row for system in systems for row in system["rows"]]
    )
    write_csv(out / "pure.csv", pure_rows(systems, mlpeg))
    write_csv(out / "isotherms.csv", isotherm_rows(systems))

    # Imported here because the report imports this module.
    from liquid_mixtures_inputs.selection_report import write_report

    write_report(
        out, selection, candidates, systems, points, sage, mlpeg, args.fast_set
    )


if __name__ == "__main__":
    main()
