"""
Select binary liquid mixtures for the binary_mixtures benchmark.

Stage 1 of the benchmark build:

1. Download the OpenFF Sage physical-property sets (pinned commits) and screen
   every binary system against the hard filters, writing ``candidates.csv``.
2. Check the curated systems in ``selection.yml`` against those filters.
3. Fetch the source papers' full isotherms and pure densities from the NIST
   TRC ThermoML archive, cross-check pure densities, fit Redlich-Kister curves
   and plan box sizes, writing ``systems.csv``, ``pure.csv`` and
   ``isotherms.csv``.
4. Write the selection report (``report.py``).

Run with::

    pixi run -e analysis select --output-dir selection \
        --fast-set selection/inputs/dataset-fast-converging.csv
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import math
from pathlib import Path
import statistics
import urllib.request

from ase.io import read
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import yaml

from liquid_mixtures_inputs.thermoml import (
    Measurement,
    canonical_smiles,
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
    f"{split}-{prop}": f"{RAW}/ash-sage-rc2/{ASH_COMMIT}/04_benchmark/phys-prop/"
    f"output/{split}/summary-benchmarks-{name}.csv"
    for split in ("training", "validation")
    for prop, name in (("hmix", "EnthalpyOfMixing"), ("density", "Density"))
}
SAGE_FORCE_FIELDS = ("openff-2.0.0", "openff-2.1.0", "openff-2.2.1", "openff-2.3.0rc1")
# Sage 2.3.0rc2 is benchmarked under its fit name: 04_benchmark/forcefields/
# sage-2.3.0rc2-v3-k100.offxml is byte-identical to fb-fit-v3-single-mean-k100.offxml.
# Tabulated in the report only, not used by the repeat-SD filter.
SAGE_RC2 = "fb-fit-v3-single-mean-k100"

# The Sage sets encode 2-pyrrolidone as its imidic acid tautomer.
TAUTOMER_FIXES = {"OC1=NCCC1": "O=C1CCCN1"}

# Hard-filter settings (see the plan / report for the rationale).
SCREEN_T_RANGE = (297.0, 299.5)
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
UCST_LIMIT = 0.5 * 8.314462618e-3 * 298.15  # RT/2 in kJ/mol
MAX_REPEAT_SD = {"hmix": 0.10, "density": 0.001}
MIN_SIGNAL = 0.4  # kJ/mol at x ~ 0.5
# Relative experimental ΔHmix uncertainties among systems passing the other
# filters are 0.5-3.5 %. Larger values are either flat 25 % uncertainties
# assigned during dataset curation or come from H^E close to zero.
MAX_REL_SIGMA_HMIX = 0.15

# Box planning (stage 2 uses these counts).
MIN_EDGE = 25.0  # Angstrom
MIN_MOLECULES = 100
MIN_MINORITY = 20
MAX_X_SHIFT = 0.01
AVOGADRO = 6.02214076e23

PURE_DENSITY_TOLERANCE = 0.002  # g/mL, cross-check flag
MAX_DENSITY_SIGMA = 0.002  # g/mL; noisier mixture densities are not used for V^E
STATE_T_TOLERANCE = 0.7  # K, window around a system's target temperature
# Systems whose enthalpy data sit just outside the window (see selection.yml).
STATE_T_EXTRA: dict[str, float] = {}


@dataclass
class Point:
    """
    One experimental point for a binary system, in the system's orientation.

    Attributes
    ----------
    prop
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
        Where the point came from (data-set name or ``thermoml``).
    """

    prop: str
    x: float
    temperature: float
    value: float
    sigma: float | None
    doi: str
    source: str


def fetch(url: str, path: Path) -> Path:
    """
    Download a file once and cache it.

    Parameters
    ----------
    url
        URL to download.
    path
        Local cache path.

    Returns
    -------
    Path
        Path to the cached file.
    """
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(url, headers={"User-Agent": "ml-peg"})
        with urllib.request.urlopen(request, timeout=120) as resp:
            path.write_bytes(resp.read())
    return path


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
    mol = Chem.MolFromSmiles(smiles)
    mol_h = Chem.AddHs(mol)
    viscous = [
        name
        for name, smarts in VISCOUS_SMARTS.items()
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts))
    ]
    return {
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "elements": {a.GetSymbol() for a in mol_h.GetAtoms()},
        "molar_mass": Descriptors.MolWt(mol),
        "n_atoms": mol_h.GetNumAtoms(),
        "viscous": viscous,
    }


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
        path = fetch(url, cache / f"{name}.csv")
        with path.open() as handle:
            for row in csv.DictReader(handle):
                c1 = fix_smiles(row["Component 1"])
                c2 = fix_smiles(row["Component 2"]) if row["Component 2"] else ""
                x1 = float(row["Mole Fraction 1"]) if row["Mole Fraction 1"] else 1.0
                if c2 and c2 < c1:
                    c1, c2, x1 = c2, c1, 1.0 - x1
                hmix = row.get("EnthalpyOfMixing Value (kJ / mol)")
                if hmix:
                    prop, value = "hmix", float(hmix)
                    sigma = row.get("EnthalpyOfMixing Uncertainty (kJ / mol)")
                else:
                    prop, value = "density", float(row["Density Value (g / ml)"])
                    sigma = row.get("Density Uncertainty (g / ml)")
                points[(c1, c2)].append(
                    Point(
                        prop=prop,
                        x=x1,
                        temperature=float(row["Temperature (K)"]),
                        value=value,
                        sigma=float(sigma) if sigma else None,
                        doi=row["Source"],
                        source=name,
                    )
                )
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
        prop = name.split("-")[-1]
        path = fetch(url, cache / f"{name}.csv")
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
                        "prop": prop,
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


def linear_fit_chi2(points: list[Point]) -> float | None:
    """
    Reduced chi^2 of a straight-line fit to interior excess-enthalpy points.

    A single-phase H^E curve gives a large value; a two-phase region, where
    H^E = a + b x, gives a value near 1. A relative uncertainty floor of 2 %
    stops tiny reported uncertainties dominating.

    Parameters
    ----------
    points
        Enthalpy of mixing points of one system at one temperature.

    Returns
    -------
    float | None
        Reduced chi^2, or ``None`` with fewer than four interior points.
    """
    pts = sorted((p for p in points if 0.05 < p.x < 0.95), key=lambda p: p.x)
    if len(pts) < 4:
        return None
    x = np.array([p.x for p in pts])
    h = np.array([p.value for p in pts])
    s = np.array([max(p.sigma or 0.0, 0.02 * abs(p.value), 1e-3) for p in pts])
    coeffs = np.polyfit(x, h, 1, w=1 / s)
    return float(np.sum(((np.polyval(coeffs, x) - h) / s) ** 2)) / (len(pts) - 2)


def two_phase_signature(points: list[Point]) -> bool:
    """
    Test whether excess enthalpies are linear in x (lever-rule signature).

    In a two-phase region H^E = a + b x. With full isotherms this fits a
    straight line through the interior points and flags reduced chi^2 <= 2.
    With exactly three points it compares the two consecutive slopes against
    their propagated uncertainty. A relative floor of 2 % stops tiny reported
    uncertainties from hiding curvature.

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
    pts = sorted((p for p in points if 0.05 < p.x < 0.95), key=lambda p: p.x)
    if len(pts) < 3 or max(abs(p.value) for p in pts) < MIN_SIGNAL:
        return False
    x = np.array([p.x for p in pts])
    h = np.array([p.value for p in pts])
    s = np.array([max(p.sigma or 0.0, 0.02 * abs(p.value), 1e-3) for p in pts])
    if len(pts) == 3:
        slopes = np.diff(h) / np.diff(x)
        var = [(s[i] ** 2 + s[i + 1] ** 2) / (x[i + 1] - x[i]) ** 2 for i in range(2)]
        return bool(abs(slopes[0] - slopes[1]) < 2 * math.sqrt(var[0] + var[1]))
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
    pure_with_sigma = {
        key[0]: p
        for key, pts in points.items()
        if not key[1]
        for p in pts
        if SCREEN_T_RANGE[0] <= p.temperature <= SCREEN_T_RANGE[1]
    }
    rows = []
    for (a, b), pts in points.items():
        if not b:
            continue
        window = [
            p for p in pts if SCREEN_T_RANGE[0] <= p.temperature <= SCREEN_T_RANGE[1]
        ]
        hmix = list({round(p.x, 3): p for p in window if p.prop == "hmix"}.values())
        dens = list({round(p.x, 3): p for p in window if p.prop == "density"}.values())
        if len(hmix) < 3 or len(dens) < 3:
            continue
        info_a, info_b = molecule_info(a), molecule_info(b)
        h_eq = min(hmix, key=lambda p: abs(p.x - 0.5)).value
        sage_rows = [r for r in sage.get((a, b), []) if not r["tautomer_error"]]
        sd = {
            prop: max(
                (
                    r["sd"]
                    for r in sage_rows
                    if r["prop"] == prop
                    and r["sd"]
                    and r["forcefield"] in SAGE_FORCE_FIELDS
                ),
                default=None,
            )
            for prop in ("hmix", "density")
        }
        mae_221 = [
            abs(r["value"] - r["reference"])
            for r in sage_rows
            if r["forcefield"] == "openff-2.2.1" and r["prop"] == "hmix"
        ]
        row = {
            "a": a,
            "b": b,
            "sources": "+".join(sorted({p.source for p in window})),
            "heavy_atoms": f"{info_a['heavy_atoms']}/{info_b['heavy_atoms']}",
            "elements": "".join(sorted(info_a["elements"] | info_b["elements"])),
            "n_hmix": len(hmix),
            "n_density": len(dens),
            "hmix_x": " ".join(f"{p.x:.3f}" for p in sorted(hmix, key=lambda p: p.x)),
            "hmix_values": " ".join(
                f"{p.value:.4f}" for p in sorted(hmix, key=lambda p: p.x)
            ),
            "hmix_eq": round(h_eq, 4),
            "max_dx": round(max_composition_mismatch(hmix, dens), 3),
            "exp_sigma_hmix": all(p.sigma for p in hmix),
            "hmix_max_rel_sigma": max(
                (p.sigma / abs(p.value) for p in hmix if p.sigma and p.value),
                default=None,
            ),
            "exp_sigma_density": all(p.sigma for p in dens),
            "pure_density_sigma": all(
                c in pure_with_sigma and pure_with_sigma[c].sigma for c in (a, b)
            ),
            "in_mlpeg": "+".join(c for c in (a, b) if c in mlpeg),
            "n_in_mlpeg": sum(c in mlpeg for c in (a, b)),
            "sage_rep_sd_hmix": sd["hmix"],
            "sage_rep_sd_density": sd["density"],
            "sage_2.2.1_hmix_mae": round(statistics.mean(mae_221), 4)
            if mae_221
            else None,
            "viscous": "+".join(info_a["viscous"] + info_b["viscous"]),
            "linear_two_phase": two_phase_signature(hmix),
            "ucst_screen": h_eq > UCST_LIMIT,
        }
        failures = []
        if row["max_dx"] > MAX_COMPOSITION_MISMATCH:
            failures.append("composition mismatch")
        if not (row["exp_sigma_hmix"] and row["exp_sigma_density"]):
            failures.append("missing experimental uncertainty")
        if (row["hmix_max_rel_sigma"] or 0.0) > MAX_REL_SIGMA_HMIX:
            failures.append(
                f"relative ΔHmix uncertainty > {100 * MAX_REL_SIGMA_HMIX:.0f} %"
            )
        if not row["pure_density_sigma"]:
            failures.append("missing pure density with uncertainty")
        if not set(row["elements"].replace("Cl", "X")) <= set("CHNOX"):
            failures.append("elements")
        if max(info_a["heavy_atoms"], info_b["heavy_atoms"]) > MAX_HEAVY_ATOMS:
            failures.append("too many heavy atoms")
        if row["viscous"]:
            failures.append("viscous")
        if row["linear_two_phase"]:
            failures.append("two-phase (linear H^E)")
        if row["ucst_screen"]:
            failures.append("two-phase (H^E > RT/2)")
        for prop in ("hmix", "density"):
            if sd[prop] is None:
                failures.append(f"no Sage repeat SD ({prop})")
            elif sd[prop] > MAX_REPEAT_SD[prop]:
                failures.append(f"Sage repeat SD ({prop})")
        if abs(h_eq) < MIN_SIGNAL:
            failures.append("small signal")
        row["failures"] = "; ".join(failures)
        row["passes"] = not failures
        rows.append(row)
    return sorted(
        rows, key=lambda r: (not r["passes"], r["n_in_mlpeg"], -abs(r["hmix_eq"]))
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
            measurements: list[Measurement] = parse_record(
                doi, fetch_record(doi, cache / "thermoml")
            )
        except Exception as err:  # noqa: BLE001 - report and continue
            errors[doi] = str(err)
            continue
        fixed_t = (corrections or {}).get(doi)
        for m in measurements:
            if fixed_t is not None:
                if len(m.smiles) == 1 or 1.0 in m.mole_fractions:
                    continue
                m.temperature = fixed_t
            if abs(m.temperature - temperature) > window:
                continue
            if m.pressure is not None and not 90.0 <= m.pressure <= 110.0:
                continue
            fractions = dict(zip(m.smiles, m.mole_fractions, strict=True))
            if set(m.smiles) == {a, b}:
                x = fractions[a]
                if x in (0.0, 1.0) and m.prop == "density":
                    comp = a if x == 1.0 else b
                    pure[comp].append(
                        Point(
                            "density",
                            1.0,
                            m.temperature,
                            m.value,
                            m.uncertainty,
                            doi,
                            "thermoml",
                        )
                    )
                elif 0.0 < x < 1.0:
                    mixture.append(
                        Point(
                            m.prop,
                            x,
                            m.temperature,
                            m.value,
                            m.uncertainty,
                            doi,
                            "thermoml",
                        )
                    )
            elif len(m.smiles) == 1 and m.smiles[0] in (a, b) and m.prop == "density":
                pure[m.smiles[0]].append(
                    Point(
                        "density",
                        1.0,
                        m.temperature,
                        m.value,
                        m.uncertainty,
                        doi,
                        "thermoml",
                    )
                )
    return mixture, dict(pure), errors


def redlich_kister(x: np.ndarray, y: np.ndarray, sigma: np.ndarray):
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
    Callable[[float], float]
        The fitted curve.
    """
    n_terms = max(1, min(4, len(x) - 1))
    basis = np.stack([x * (1 - x) * (2 * x - 1) ** k for k in range(n_terms)], axis=1)
    weights = 1.0 / sigma
    coeffs, *_ = np.linalg.lstsq(basis * weights[:, None], y * weights, rcond=None)

    def curve(value: float) -> float:
        """
        Evaluate the fitted curve.

        Parameters
        ----------
        value
            Mole fraction.

        Returns
        -------
        float
            Fitted excess property.
        """
        return float(
            sum(
                c * value * (1 - value) * (2 * value - 1) ** k
                for k, c in enumerate(coeffs)
            )
        )

    return curve


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
    pure = x in (0.0, 1.0)
    while True:
        n_a = round(x * n_total)
        n_b = n_total - n_a
        ok = pure or (
            abs(n_a / n_total - x) <= MAX_X_SHIFT and min(n_a, n_b) >= MIN_MINORITY
        )
        if ok:
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
        Chosen value, uncertainty and source, all cross-check values and a
        flag if they disagree by more than 0.002 g/mL.
    """
    near = [
        p
        for p in thermoml_pure
        if abs(p.temperature - temperature) <= STATE_T_TOLERANCE
    ]
    near_sage = [
        p for p in sage_pure if abs(p.temperature - temperature) <= STATE_T_TOLERANCE
    ]
    checks = {f"thermoml:{p.doi}": p.value for p in near}
    checks.update({f"{p.source}:{p.doi}": p.value for p in near_sage})
    ml = mlpeg.get(smiles)
    if ml and abs(ml["temperature"] - temperature) <= 1.0:
        checks["ml-peg"] = ml["density"]
    same_paper = [p for doi in (preferred_dois or []) for p in near if p.doi == doi]
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
    elif near:
        best = statistics.median(p.value for p in near)
        chosen = min(near, key=lambda p: abs(p.value - best))
    elif near_sage:
        chosen = min(near_sage, key=lambda p: p.sigma or 1.0)
    elif ml and abs(ml["temperature"] - temperature) <= 1.0:
        chosen = Point(
            "density", 1.0, ml["temperature"], ml["density"], None, "ml-peg", "ml-peg"
        )
    else:
        chosen = None
    values = list(checks.values())
    spread = max(values) - min(values) if len(values) > 1 else 0.0
    outliers = (
        [k for k, v in checks.items() if abs(v - chosen.value) > PURE_DENSITY_TOLERANCE]
        if chosen
        else []
    )
    return {
        "density": chosen.value if chosen else None,
        "sigma": chosen.sigma if chosen else None,
        "source": (
            chosen.doi
            if chosen and chosen.source == "reference"
            else f"{chosen.source}:{chosen.doi}"
            if chosen
            else ""
        ),
        "checks": checks,
        "spread": spread,
        "outliers": outliers,
        "mlpeg_density": ml["density"] if ml else None,
        "mlpeg_temperature": ml["temperature"] if ml else None,
    }


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
    window = STATE_T_EXTRA.get(spec["id"], STATE_T_TOLERANCE)
    key = tuple(sorted((a, b)))
    flip = key[0] != a

    def oriented(p: Point) -> Point:
        """
        Re-express a Sage-set point with x referring to component a.

        Parameters
        ----------
        p
            Point keyed by the sorted pair.

        Returns
        -------
        Point
            Point with x of component a.
        """
        return Point(
            p.prop,
            1 - p.x if flip else p.x,
            p.temperature,
            p.value,
            p.sigma,
            p.doi,
            p.source,
        )

    set_points = [oriented(p) for p in points.get(key, [])]
    dois = {p.doi for p in set_points if abs(p.temperature - temperature) <= window}
    sage_pure = {c: points.get((c, ""), []) for c in (a, b)}
    dois |= {
        p.doi
        for c in (a, b)
        for p in sage_pure[c]
        if abs(p.temperature - temperature) <= STATE_T_TOLERANCE
    }
    dois |= set(spec.get("extra_dois", []))
    mixture, thermo_pure, errors = thermoml_points(
        a, b, dois, temperature, window, cache, spec.get("temperature_corrections")
    )

    hmix = sorted((p for p in mixture if p.prop == "hmix"), key=lambda p: p.x)
    all_dens = sorted((p for p in mixture if p.prop == "density"), key=lambda p: p.x)
    dens = [p for p in all_dens if p.sigma is None or p.sigma <= MAX_DENSITY_SIGMA]
    noisy_density = Counter(p.doi for p in all_dens if p not in dens)
    density_dois = [doi for doi, _ in Counter(p.doi for p in dens).most_common()]

    info = {a: molecule_info(a), b: molecule_info(b)}
    pures = {
        c: choose_pure_density(
            c,
            temperature,
            thermo_pure.get(c, []),
            sage_pure[c],
            mlpeg,
            density_dois,
            next(
                (
                    r
                    for r in references
                    if canonical_smiles(r["smiles"]) == c
                    and abs(r["temperature"] - temperature) < 0.1
                ),
                None,
            ),
        )
        for c in (a, b)
    }
    vex_direct = sorted((p for p in mixture if p.prop == "vex"), key=lambda p: p.x)

    # Excess volumes from mixture densities and the chosen pure densities.
    ma, mb = info[a]["molar_mass"], info[b]["molar_mass"]
    vex = list(vex_direct)
    if pures[a]["density"] and pures[b]["density"]:
        for p in dens:
            mean_mass = p.x * ma + (1 - p.x) * mb
            ideal = (
                p.x * ma / pures[a]["density"] + (1 - p.x) * mb / pures[b]["density"]
            )
            sigma = mean_mass / p.value**2 * p.sigma if p.sigma else None
            vex.append(
                Point(
                    "vex",
                    p.x,
                    p.temperature,
                    mean_mass / p.value - ideal,
                    sigma,
                    p.doi,
                    "derived",
                )
            )

    def fit(pts: list[Point], floor: float):
        """
        Fit a Redlich-Kister curve to excess-property points.

        Parameters
        ----------
        pts
            Points to fit.
        floor
            Minimum uncertainty used for weighting.

        Returns
        -------
        Callable[[float], float] | None
            Fitted curve, or ``None`` if there are no points.
        """
        if not pts:
            return None
        x = np.array([p.x for p in pts])
        y = np.array([p.value for p in pts])
        s = np.array([max(p.sigma or 0.0, floor) for p in pts])
        return redlich_kister(x, y, s)

    hmix_fit = fit(hmix, 1e-3)
    vex_fit = fit(vex, 1e-3)

    # Target compositions: measured enthalpy points nearest 0.25, 0.5 and 0.75.
    targets = []
    for goal in (0.25, 0.5, 0.75):
        if hmix:
            nearest = min(hmix, key=lambda p: abs(p.x - goal))
            if abs(nearest.x - goal) <= 0.1 and nearest not in targets:
                targets.append(nearest)

    rows = []
    for target in targets:
        rho_guess = (
            pures[a]["density"] * target.x + pures[b]["density"] * (1 - target.x)
            if (pures[a]["density"] and pures[b]["density"])
            else 1.0
        )
        box = plan_box(target.x, rho_guess, info[a], info[b])
        x_box = box["x_box"]
        row = {
            "system_id": spec["id"],
            "a": a,
            "b": b,
            "temperature": temperature,
            "x_target": round(target.x, 4),
            "x_box": round(x_box, 4),
            "n_a": box["n_a"],
            "n_b": box["n_b"],
            "n_atoms": box["n_atoms"],
            "hmix_measured": target.value,
            "hmix_measured_sigma": target.sigma,
            "hmix_measured_T": target.temperature,
            "hmix_exp": round(hmix_fit(x_box) + (target.value - hmix_fit(target.x)), 5),
            "hmix_sigma": target.sigma,
            "hmix_doi": target.doi,
        }
        if vex_fit and pures[a]["density"] and pures[b]["density"]:
            mean_mass = x_box * ma + (1 - x_box) * mb
            ideal = (
                x_box * ma / pures[a]["density"]
                + (1 - x_box) * mb / pures[b]["density"]
            )
            vex_box = vex_fit(x_box)
            row["vex_exp"] = round(vex_box, 5)
            row["density_exp"] = round(mean_mass / (ideal + vex_box), 5)
            nearest_d = min(dens, key=lambda p: abs(p.x - x_box)) if dens else None
            # Uncertainty from the main density paper near x_box (other papers
            # sometimes report implausible uncertainties).
            local = [
                p.sigma
                for p in dens
                if abs(p.x - x_box) <= 0.05 and p.sigma and p.doi == density_dois[0]
            ] or [p.sigma for p in dens if abs(p.x - x_box) <= 0.05 and p.sigma]
            row["density_sigma"] = (
                statistics.median(local)
                if local
                else (nearest_d.sigma if nearest_d else None)
            )
            row["density_doi"] = nearest_d.doi if nearest_d else ""
            row["density_nearest_dx"] = (
                round(abs(nearest_d.x - x_box), 3) if nearest_d else None
            )
            row["edge"] = round(
                plan_box(target.x, row["density_exp"], info[a], info[b])["edge"], 2
            )
        else:
            row.update(
                vex_exp=None,
                density_exp=None,
                density_sigma=None,
                density_doi="",
                density_nearest_dx=None,
                edge=round(box["edge"], 2),
            )
        rows.append(row)

    sage_rows = [r for r in sage.get(key, []) if not r["tautomer_error"]]
    return {
        "spec": spec,
        "a": a,
        "b": b,
        "temperature": temperature,
        "info": info,
        "pures": pures,
        "hmix": hmix,
        "density": dens,
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
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    """Run the selection and write all stage 1 outputs."""
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

    pure_rows, seen = [], set()
    for system in systems:
        temperature = system["temperature"]
        for smiles in (system["a"], system["b"]):
            if (smiles, temperature) in seen:
                continue
            seen.add((smiles, temperature))
            pure = system["pures"][smiles]
            info = system["info"][smiles]
            box = plan_box(1.0, pure["density"] or 1.0, info, info)
            pure_rows.append(
                {
                    "smiles": smiles,
                    "temperature": temperature,
                    "density": pure["density"],
                    "sigma": pure["sigma"],
                    "source": pure["source"],
                    "in_mlpeg": smiles in mlpeg,
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
    write_csv(out / "pure.csv", pure_rows)

    iso_rows = [
        {
            "system_id": system["spec"]["id"],
            "prop": point.prop,
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
    write_csv(out / "isotherms.csv", iso_rows)

    from liquid_mixtures_inputs.report import write_report

    write_report(
        out, selection, candidates, systems, points, sage, mlpeg, args.fast_set
    )


if __name__ == "__main__":
    main()
