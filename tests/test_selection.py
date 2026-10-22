"""Test the binary_mixtures system-selection helpers (offline)."""

from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem

from liquid_mixtures_inputs import select_systems as ss
from liquid_mixtures_inputs.thermoml import (
    mol_from_inchi,
    parse_record,
)


@pytest.mark.parametrize(
    "smiles", ["CNC=O", "NC=O", "O=C1CCCN1", "CC(=O)NC", "CN(C)C=O", "c1ccncc1"]
)
def test_inchi_round_trip_keeps_amides(smiles):
    """Test standard-InChI round trips do not produce imidic acid tautomers."""
    inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles))
    assert Chem.MolToSmiles(mol_from_inchi(inchi)) == Chem.MolToSmiles(
        Chem.MolFromSmiles(smiles)
    )


def _points(xs, values, sigma=0.005):
    """
    Build enthalpy points for a test isotherm.

    Parameters
    ----------
    xs
        Mole fractions.
    values
        Excess enthalpies.
    sigma
        Uncertainty applied to every point.

    Returns
    -------
    list[ss.Point]
        Points.
    """
    return [
        ss.Point("hmix", x, 298.15, v, sigma, "doi", "test")
        for x, v in zip(xs, values, strict=True)
    ]


def test_two_phase_signature():
    """Test the linear-H^E test flags lever-rule data but not a curved isotherm."""
    # Acetophenone + water from the fast-converging list (two-phase).
    assert ss.two_phase_signature(_points([0.2, 0.5, 0.8], [0.1379, 0.3413, 0.5414]))
    # Pyrrole + pyridine (single phase, strongly curved).
    assert not ss.two_phase_signature(
        _points([0.239, 0.509, 0.76], [-1.7662, -2.649, -1.6452])
    )
    # Full isotherms: a straight line versus a parabola.
    xs = np.linspace(0.1, 0.9, 9)
    assert ss.two_phase_signature(_points(xs, 0.05 + 0.6 * xs))
    assert not ss.two_phase_signature(_points(xs, 2.0 * xs * (1 - xs)))
    # Ideal and near-ideal mixtures are straight lines but not two-phase.
    assert not ss.two_phase_signature(_points([0.25, 0.5, 0.75], [0.0, 0.0, 0.0]))
    assert not ss.two_phase_signature(_points(xs, 0.2 * xs * (1 - xs)))


def test_redlich_kister_recovers_curve():
    """Test the Redlich-Kister fit reproduces a known two-term curve."""
    xs = np.linspace(0.05, 0.95, 12)
    ys = xs * (1 - xs) * (-3.0 + 0.8 * (2 * xs - 1))
    curve = ss.redlich_kister(xs, ys, np.full_like(xs, 1e-3))
    assert curve(0.37) == pytest.approx(0.37 * 0.63 * (-3.0 + 0.8 * (2 * 0.37 - 1)))


def test_plan_box_constraints():
    """Test box plans respect the edge, molecule and composition rules."""
    info_a = ss.molecule_info("CCCCC(=O)OC")
    info_b = ss.molecule_info("CCCCCCCC")
    box = ss.plan_box(0.2255, 0.73, info_a, info_b)
    assert box["edge"] >= ss.MIN_EDGE - 0.5
    assert box["n_a"] + box["n_b"] >= ss.MIN_MOLECULES
    assert min(box["n_a"], box["n_b"]) >= ss.MIN_MINORITY
    assert abs(box["x_box"] - 0.2255) <= ss.MAX_X_SHIFT


def test_parse_record_binary_enthalpy():
    """Test parsing of a minimal ThermoML enthalpy-of-mixing block."""
    water = Chem.MolToInchi(Chem.MolFromSmiles("O"))
    methanol = Chem.MolToInchi(Chem.MolFromSmiles("CO"))
    record = {
        "Compound": [
            {"RegNum": {"nOrgNum": 1}, "sStandardInChI": methanol},
            {"RegNum": {"nOrgNum": 2}, "sStandardInChI": water},
        ],
        "PureOrMixtureData": [
            {
                "Component": [{"RegNum": {"nOrgNum": 1}}, {"RegNum": {"nOrgNum": 2}}],
                "PhaseID": [{"ePhase": "Liquid"}],
                "Property": [
                    {
                        "nPropNumber": 1,
                        "Property-MethodID": {
                            "PropertyGroup": {
                                "ExcessPartialApparentEnergyProp": {
                                    "ePropName": "Excess molar enthalpy (molar "
                                    "enthalpy of mixing), kJ/mol"
                                }
                            }
                        },
                    }
                ],
                "Constraint": [
                    {
                        "nConstraintNumber": 1,
                        "nConstraintValue": 298.15,
                        "ConstraintID": {
                            "ConstraintType": {"eTemperature": "Temperature, K"}
                        },
                    }
                ],
                "Variable": [
                    {
                        "nVarNumber": 1,
                        "VariableID": {
                            "RegNum": {"nOrgNum": 1},
                            "VariableType": {"eComponentComposition": "Mole fraction"},
                        },
                    }
                ],
                "NumValues": [
                    {
                        "VariableValue": [{"nVarNumber": 1, "nVarValue": 0.3}],
                        "PropertyValue": [
                            {
                                "nPropNumber": 1,
                                "nPropValue": -0.85,
                                "CombinedUncertainty": {"nCombExpandUncertValue": 0.01},
                            }
                        ],
                    }
                ],
            }
        ],
    }
    (m,) = parse_record("10.0/test", record)
    assert m.prop == "hmix"
    assert dict(zip(m.smiles, m.mole_fractions, strict=True)) == pytest.approx(
        {"CO": 0.3, "O": 0.7}
    )
    assert m.temperature == 298.15
    assert m.value == -0.85
    assert m.uncertainty == pytest.approx(0.005)
