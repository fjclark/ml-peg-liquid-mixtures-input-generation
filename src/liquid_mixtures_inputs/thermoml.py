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
    prop
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
    prop: str
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
    if not path.exists():
        cache_dir.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(
            THERMOML_URL.format(doi=doi), headers={"User-Agent": "ml-peg"}
        )
        with urllib.request.urlopen(request, timeout=60) as resp:
            path.write_bytes(resp.read())
    return json.loads(path.read_text())


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
        unc = container.get(key)
        if isinstance(unc, list):
            unc = unc[0] if unc else None
        if not unc:
            continue
        for std_key in ("nCombStdUncertValue", "nStdUncertValue"):
            if std_key in unc:
                return float(unc[std_key])
        for exp_key in ("nCombExpandUncertValue", "nExpandUncertValue"):
            if exp_key in unc:
                return float(unc[exp_key]) / 2.0
    return None


def _property_name(prop: dict) -> str | None:
    """
    Return the ThermoML property name of a ``Property`` entry.

    Parameters
    ----------
    prop
        ``Property`` entry of a ``PureOrMixtureData`` block.

    Returns
    -------
    str | None
        The ``ePropName`` string, or ``None`` if it cannot be found.
    """
    group = prop["Property-MethodID"]["PropertyGroup"]
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
        ident = entry.get("VariableID") or entry.get("ConstraintID")
        number = entry.get("nVarNumber", entry.get("nConstraintNumber"))
        if number not in values:
            continue
        kind = ident.get("VariableType") or ident.get("ConstraintType")
        if "eTemperature" in kind:
            temperature = values[number]
        elif "ePressure" in kind:
            pressure = values[number]
        elif kind.get("eComponentComposition") == "Mole fraction":
            fractions[org_to_smiles[ident["RegNum"]["nOrgNum"]]] = values[number]
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
        orgs = [c["RegNum"]["nOrgNum"] for c in block["Component"]]
        phases = {p["ePhase"] for p in block.get("PhaseID", [])}
        if (
            phases != {"Liquid"}
            or len(orgs) > 2
            or any(o not in org_to_smiles for o in orgs)
        ):
            continue
        smiles = tuple(org_to_smiles[o] for o in orgs)
        props = {}
        for prop in block["Property"]:
            name = _property_name(prop)
            if name in PROPERTIES:
                props[prop["nPropNumber"]] = PROPERTIES[name]
        if not props:
            continue

        constraint_values = {
            c["nConstraintNumber"]: float(c["nConstraintValue"])
            for c in block.get("Constraint", [])
        }
        t_c, p_c, x_c = _parse_state(
            block.get("Constraint", []), constraint_values, org_to_smiles
        )
        for row in block["NumValues"]:
            var_values = {
                v["nVarNumber"]: float(v["nVarValue"])
                for v in row.get("VariableValue", [])
            }
            t_v, p_v, x_v = _parse_state(
                block.get("Variable", []), var_values, org_to_smiles
            )
            temperature = t_v if t_v is not None else t_c
            pressure = p_v if p_v is not None else p_c
            fractions = {**x_c, **x_v}
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
            for pv in row["PropertyValue"]:
                if pv["nPropNumber"] not in props:
                    continue
                prop, factor = props[pv["nPropNumber"]]
                unc = _standard_uncertainty(pv)
                measurements.append(
                    Measurement(
                        doi=doi,
                        prop=prop,
                        smiles=smiles,
                        mole_fractions=tuple(fractions[s] for s in smiles),
                        temperature=temperature,
                        pressure=pressure,
                        value=float(pv["nPropValue"]) * factor,
                        uncertainty=None if unc is None else unc * factor,
                    )
                )
    return measurements
