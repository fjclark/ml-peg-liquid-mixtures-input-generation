"""Test the binary_mixtures system-selection helpers (offline)."""

from __future__ import annotations

import numpy as np
import pytest
from rdkit import Chem

from liquid_mixtures_inputs import select_systems as ss
from liquid_mixtures_inputs.thermoml import (
    Measurement,
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


# Analytic binary mixture for the build_system tests: ethanol (A) + heptane (B),
# with x the mole fraction of A. Pure densities in g/mL, V^E in cm3/mol and
# H^E in kJ/mol. Both excess curves are two-term Redlich-Kister expansions, so
# the fits in build_system can reproduce them exactly.
ETHANOL, HEPTANE = "CCO", "CCCCCCC"
RHO_A, RHO_B = 0.785, 0.680
# Offset from round values, so every box x differs from its measured target
GRID = np.round(np.arange(0.05, 0.951, 0.05) + 0.0037, 4)
DOI = "10.0/analytic"


def _vex(x):
    """
    Return the analytic excess molar volume.

    Parameters
    ----------
    x
        Mole fraction of A.

    Returns
    -------
    float
        V^E in cm3/mol.
    """
    return x * (1 - x) * (-1.2 + 0.3 * (2 * x - 1))


def _hmix(x):
    """
    Return the analytic excess molar enthalpy.

    Parameters
    ----------
    x
        Mole fraction of A.

    Returns
    -------
    float
        H^E in kJ/mol.
    """
    return x * (1 - x) * (2.0 - 0.5 * (2 * x - 1))


def _density(x):
    """
    Return the analytic mixture density.

    Parameters
    ----------
    x
        Mole fraction of A.

    Returns
    -------
    float
        Density in g/mL.
    """
    m_a = ss.molecule_info(ETHANOL)["molar_mass"]
    m_b = ss.molecule_info(HEPTANE)["molar_mass"]
    ideal = x * m_a / RHO_A + (1 - x) * m_b / RHO_B
    return (x * m_a + (1 - x) * m_b) / (ideal + _vex(x))


def _analytic_record(record_order):
    """
    Build ThermoML measurements of the analytic mixture.

    Parameters
    ----------
    record_order
        Component order of the record, ``(A, B)`` or ``(B, A)``.

    Returns
    -------
    list[Measurement]
        Pure densities, mixture densities and excess enthalpies at 298.15 K.
    """
    measurements = [
        Measurement(DOI, "density", (c,), (1.0,), 298.15, 101.325, rho, 1e-4)
        for c, rho in ((ETHANOL, RHO_A), (HEPTANE, RHO_B))
    ]
    for x in GRID:
        fractions = {ETHANOL: x, HEPTANE: 1 - x}
        mole_fractions = tuple(fractions[c] for c in record_order)
        measurements += [
            Measurement(
                DOI,
                "density",
                record_order,
                mole_fractions,
                298.15,
                101.325,
                _density(x),
                1e-4,
            ),
            Measurement(
                DOI,
                "hmix",
                record_order,
                mole_fractions,
                298.15,
                101.325,
                _hmix(x),
                5e-3,
            ),
        ]
    return measurements


@pytest.mark.parametrize("record_order", [(ETHANOL, HEPTANE), (HEPTANE, ETHANOL)])
@pytest.mark.parametrize("spec_order", [(ETHANOL, HEPTANE), (HEPTANE, ETHANOL)])
def test_build_system_reproduces_analytic_mixture(
    monkeypatch, tmp_path, record_order, spec_order
):
    """
    Test the reference values equal the analytic mixture at each box x.

    The box compositions differ from the measured ones by up to
    ``MAX_X_SHIFT``, so this also tests the shift of ΔHmix to x_box. Swapping
    the components in the record or in the spec must give x -> 1 - x and the
    same values.
    """
    monkeypatch.setattr(ss, "fetch_record", lambda doi, cache: None)
    monkeypatch.setattr(
        ss, "parse_record", lambda doi, record: _analytic_record(record_order)
    )
    spec = {
        "id": "analytic",
        "a": spec_order[0],
        "b": spec_order[1],
        "temperature": 298.15,
        "extra_dois": [DOI],
    }
    system = ss.build_system(spec, {}, {}, {}, tmp_path)

    assert system["pures"][ETHANOL]["density"] == RHO_A
    assert system["pures"][HEPTANE]["density"] == RHO_B
    assert len(system["rows"]) == 3
    assert all(row["x_box"] != row["x_target"] for row in system["rows"])
    for row in system["rows"]:
        # Exact box composition (systems.csv rounds x_box to 4 decimals)
        x_box = row["n_a"] / (row["n_a"] + row["n_b"])
        x_a = x_box if spec_order[0] == ETHANOL else 1 - x_box
        # Values are rounded to 5 decimals in build_system
        assert row["hmix_exp"] == pytest.approx(_hmix(x_a), abs=1e-5)
        assert row["vex_exp"] == pytest.approx(_vex(x_a), abs=1e-5)
        assert row["density_exp"] == pytest.approx(_density(x_a), abs=1e-5)


@pytest.mark.parametrize(
    ("name", "raw", "expected"),
    [
        ("Mass density, kg/m3", 1000.0, 1.0),
        ("Excess molar volume, m3/mol", 1e-6, 1.0),
        ("Excess molar enthalpy (molar enthalpy of mixing), kJ/mol", 1.0, 1.0),
    ],
)
def test_parse_record_unit_factors(name, raw, expected):
    """Test ThermoML values convert to g/mL, cm3/mol and kJ/mol."""
    record = {
        "Compound": [
            {
                "RegNum": {"nOrgNum": 1},
                "sStandardInChI": Chem.MolToInchi(Chem.MolFromSmiles("CO")),
            },
            {
                "RegNum": {"nOrgNum": 2},
                "sStandardInChI": Chem.MolToInchi(Chem.MolFromSmiles("O")),
            },
        ],
        "PureOrMixtureData": [
            {
                "Component": [{"RegNum": {"nOrgNum": 1}}, {"RegNum": {"nOrgNum": 2}}],
                "PhaseID": [{"ePhase": "Liquid"}],
                "Property": [
                    {
                        "nPropNumber": 1,
                        "Property-MethodID": {
                            "PropertyGroup": {"Group": {"ePropName": name}}
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
                        "VariableValue": [{"nVarNumber": 1, "nVarValue": 0.5}],
                        "PropertyValue": [
                            {
                                "nPropNumber": 1,
                                "nPropValue": raw,
                                "CombinedUncertainty": {
                                    "nCombStdUncertValue": raw / 100
                                },
                            }
                        ],
                    }
                ],
            }
        ],
    }
    (m,) = parse_record("10.0/test", record)
    assert m.value == pytest.approx(expected)
    assert m.uncertainty == pytest.approx(expected / 100)
