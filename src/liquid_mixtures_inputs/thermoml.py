"""Fetch and parse NIST TRC ThermoML records for binary liquid mixtures."""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
import urllib.request

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

# Standard InChI does not fix mobile hydrogens, so RDKit can return the imidic
# acid tautomer of an amide (e.g. CN=CO for N-methylformamide, the same error
# that gives OC1=NCCC1 for 2-pyrrolidone in the Sage sets). Move the H back.
IMIDIC_ACID = Chem.MolFromSmarts("[OX2H1]-[#6]=[#7]")
TO_AMIDE = AllChem.ReactionFromSmarts("[O;H1:1]-[#6:2]=[#7:3]>>[O;H0:1]=[#6:2]-[#7:3]")

THERMOML_URL = "https://trc.nist.gov/ThermoML/{doi}.json"

# ThermoML property names we use, mapped to (short name, factor to our units).
# Densities are converted to g/mL, excess volumes to cm3/mol.
PROPERTIES = {
    "Mass density, kg/m3": ("density", 1e-3),
    "Excess molar enthalpy (molar enthalpy of mixing), kJ/mol": ("hmix", 1.0),
    "Excess molar volume, m3/mol": ("vex", 1e6),
}


@dataclass
class Measurement:
    """
    One liquid-phase measurement from a ThermoML record.

    Attributes
    ----------
    doi
        DOI of the source paper.
    property_name
        Short property name: ``density``, ``hmix`` or ``vex``.
    smiles
        Canonical SMILES of each component, in record order.
    mole_fractions
        Mole fraction of each component, aligned with ``smiles``.
    temperature
        Temperature in K.
    pressure
        Pressure in kPa, or ``None`` if not reported.
    value
        Property value in g/mL, kJ/mol or cm3/mol.
    uncertainty
        Standard uncertainty in the same units, or ``None`` if not reported.
    """

    doi: str
    property_name: str
    smiles: tuple[str, ...]
    mole_fractions: tuple[float, ...]
    temperature: float
    pressure: float | None
    value: float
    uncertainty: float | None = field(default=None)


def canonical_smiles(smiles: str) -> str:
    """
    Return the RDKit canonical SMILES of a molecule.

    Parameters
    ----------
    smiles
        Input SMILES.

    Returns
    -------
    str
        Canonical SMILES.
    """
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def mol_from_inchi(inchi: str) -> Chem.Mol | None:
    """
    Build a molecule from an InChI, restoring amide tautomers.

    Parameters
    ----------
    inchi
        Standard InChI string.

    Returns
    -------
    Chem.Mol | None
        The molecule, or ``None`` if RDKit cannot parse the InChI.
    """
    mol = Chem.MolFromInchi(inchi)
    while mol is not None and mol.HasSubstructMatch(IMIDIC_ACID):
        product = TO_AMIDE.RunReactants((mol,))[0][0]
        for atom in product.GetAtoms():
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)
        Chem.SanitizeMol(product)
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(product))
    return mol


def download(url: str, path: Path) -> Path:
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
        with urllib.request.urlopen(request, timeout=120) as response:
            path.write_bytes(response.read())
    return path


def fetch_record(doi: str, cache_dir: Path) -> dict:
    """
    Download (or read from cache) the ThermoML JSON record for a DOI.

    Parameters
    ----------
    doi
        DOI of the paper, e.g. ``10.1021/je301238f``.
    cache_dir
        Directory used to cache downloaded records.

    Returns
    -------
    dict
        Parsed JSON record.
    """
    path = cache_dir / f"{doi.replace('/', '_')}.json"
    return json.loads(download(THERMOML_URL.format(doi=doi), path).read_text())


def _standard_uncertainty(container: dict) -> float | None:
    """
    Extract a standard uncertainty from a ThermoML property value.

    Expanded uncertainties are divided by the coverage factor 2 (ThermoML
    expanded uncertainties are reported at 95 % confidence).

    Parameters
    ----------
    container
        ``PropertyValue`` entry of a ThermoML record.

    Returns
    -------
    float | None
        Standard uncertainty in ThermoML units, or ``None`` if absent.
    """
    for key in ("CombinedUncertainty", "PropUncertainty"):
        uncertainty = container.get(key)
        if isinstance(uncertainty, list):
            uncertainty = uncertainty[0] if uncertainty else None
        if not uncertainty:
            continue
        for standard_key in ("nCombStdUncertValue", "nStdUncertValue"):
            if standard_key in uncertainty:
                return float(uncertainty[standard_key])
        for expanded_key in ("nCombExpandUncertValue", "nExpandUncertValue"):
            if expanded_key in uncertainty:
                return float(uncertainty[expanded_key]) / 2.0
    return None


def _property_name(property_entry: dict) -> str | None:
    """
    Return the ThermoML property name of a ``Property`` entry.

    Parameters
    ----------
    property_entry
        ``Property`` entry of a ``PureOrMixtureData`` block.

    Returns
    -------
    str | None
        The ``ePropName`` string, or ``None`` if it cannot be found.
    """
    group = property_entry["Property-MethodID"]["PropertyGroup"]
    for value in group.values():
        if isinstance(value, dict) and "ePropName" in value:
            return value["ePropName"]
    return None


def _parse_state(entries: list[dict], values: dict[int, float], org_to_smiles):
    """
    Read temperature, pressure and mole fractions from variables or constraints.

    Parameters
    ----------
    entries
        ``Variable`` or ``Constraint`` entries of a ThermoML block.
    values
        Mapping from variable/constraint number to its value.
    org_to_smiles
        Mapping from ThermoML compound number to canonical SMILES.

    Returns
    -------
    tuple[float | None, float | None, dict[str, float]]
        Temperature (K), pressure (kPa) and mole fractions by SMILES.
    """
    temperature = pressure = None
    fractions: dict[str, float] = {}
    for entry in entries:
        identifier = entry.get("VariableID") or entry.get("ConstraintID")
        number = entry.get("nVarNumber", entry.get("nConstraintNumber"))
        if number not in values:
            continue
        kind = identifier.get("VariableType") or identifier.get("ConstraintType")
        if "eTemperature" in kind:
            temperature = values[number]
        elif "ePressure" in kind:
            pressure = values[number]
        elif kind.get("eComponentComposition") == "Mole fraction":
            compound = identifier["RegNum"]["nOrgNum"]
            fractions[org_to_smiles[compound]] = values[number]
    return temperature, pressure, fractions


def parse_record(doi: str, record: dict) -> list[Measurement]:
    """
    Extract liquid density, excess enthalpy and excess volume measurements.

    Parameters
    ----------
    doi
        DOI of the source paper.
    record
        Parsed ThermoML JSON record.

    Returns
    -------
    list[Measurement]
        All pure-liquid and binary liquid measurements of supported properties.
    """
    org_to_smiles = {}
    for compound in record["Compound"]:
        mol = mol_from_inchi(compound["sStandardInChI"])
        if mol is not None:
            org_to_smiles[compound["RegNum"]["nOrgNum"]] = Chem.MolToSmiles(mol)

    measurements = []
    for block in record.get("PureOrMixtureData", []):
        compounds = [c["RegNum"]["nOrgNum"] for c in block["Component"]]
        phases = {p["ePhase"] for p in block.get("PhaseID", [])}
        if (
            phases != {"Liquid"}
            or len(compounds) > 2
            or any(c not in org_to_smiles for c in compounds)
        ):
            continue
        smiles = tuple(org_to_smiles[c] for c in compounds)
        # ThermoML property number -> (our property name, unit factor)
        wanted = {}
        for property_entry in block["Property"]:
            name = _property_name(property_entry)
            if name in PROPERTIES:
                wanted[property_entry["nPropNumber"]] = PROPERTIES[name]
        if not wanted:
            continue

        # The state is set by constraints (fixed for the block) and variables
        # (per row); a variable overrides a constraint.
        constraint_values = {
            c["nConstraintNumber"]: float(c["nConstraintValue"])
            for c in block.get("Constraint", [])
        }
        fixed_temperature, fixed_pressure, fixed_fractions = _parse_state(
            block.get("Constraint", []), constraint_values, org_to_smiles
        )
        for row in block["NumValues"]:
            variable_values = {
                v["nVarNumber"]: float(v["nVarValue"])
                for v in row.get("VariableValue", [])
            }
            row_temperature, row_pressure, row_fractions = _parse_state(
                block.get("Variable", []), variable_values, org_to_smiles
            )
            temperature = (
                row_temperature if row_temperature is not None else fixed_temperature
            )
            pressure = row_pressure if row_pressure is not None else fixed_pressure
            fractions = {**fixed_fractions, **row_fractions}
            if temperature is None:
                continue
            if len(smiles) == 1:
                fractions = {smiles[0]: 1.0}
            elif len(fractions) == 1:
                ((known, x),) = fractions.items()
                other = next(s for s in smiles if s != known)
                fractions[other] = 1.0 - x
            if set(fractions) != set(smiles):
                continue
            for property_value in row["PropertyValue"]:
                if property_value["nPropNumber"] not in wanted:
                    continue
                property_name, factor = wanted[property_value["nPropNumber"]]
                uncertainty = _standard_uncertainty(property_value)
                measurements.append(
                    Measurement(
                        doi=doi,
                        property_name=property_name,
                        smiles=smiles,
                        mole_fractions=tuple(fractions[s] for s in smiles),
                        temperature=temperature,
                        pressure=pressure,
                        value=float(property_value["nPropValue"]) * factor,
                        uncertainty=None
                        if uncertainty is None
                        else uncertainty * factor,
                    )
                )
    return measurements
