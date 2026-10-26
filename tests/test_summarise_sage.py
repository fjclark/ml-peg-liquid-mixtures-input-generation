"""Test the Sage summary statistics."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from liquid_mixtures_inputs.summarise_sage import (
    SIZE_GAP,
    compare_versions,
    finite_size_check,
    mixture_points,
    pooled_sigma,
    version_mae,
    within_between,
)
from ml_peg.analysis.utils.timeseries import estimate_mean


def ar1(n_samples: int, phi: float, rng: np.random.Generator) -> np.ndarray:
    """
    Generate an AR(1) series with unit noise.

    Parameters
    ----------
    n_samples
        Number of samples.
    phi
        Autoregressive coefficient.
    rng
        Random generator.

    Returns
    -------
    np.ndarray
        Time series.
    """
    noise = rng.standard_normal(n_samples)
    series = np.zeros(n_samples)
    for i in range(1, n_samples):
        series[i] = phi * series[i - 1] + noise[i]
    return series


def test_within_between():
    """Correct within-run errors pass; per-replicate offsets are flagged."""
    rng = np.random.default_rng(1)
    offsets = [-0.4, -0.2, 0.0, 0.2, 0.4]
    rows = []
    for box in range(20):
        for rep in range(5):
            offset = offsets[rep] if box == 0 else 0.0
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                estimate = estimate_mean(ar1(5000, 0.8, rng) + offset)
            rows.append(
                {
                    "box_id": f"box{box}",
                    "replicate": rep,
                    "h_mean": estimate.mean,
                    "h_sigma": estimate.sem,
                    "density_mean": estimate.mean,
                    "density_sigma": estimate.sem,
                    "hmix": estimate.mean,
                    "hmix_sigma_run": estimate.sem,
                    "vex": estimate.mean,
                    "vex_sigma_run": estimate.sem,
                }
            )
    runs = pd.DataFrame(rows)
    table = within_between(runs, runs)
    flagged = set(table[table.flag == "high"].box_id)
    assert flagged == {"box0"}
    assert (table.flag == "low").sum() == 0


def _runs(x: float, excess_h: float, excess_v: float) -> tuple[pd.DataFrame, ...]:
    """
    Build per-run statistics for one mixture box and its two pure boxes.

    Parameters
    ----------
    x
        Mole fraction of component a in the mixture box.
    excess_h
        Mixture h minus the mole-fraction average of the pure h (kJ/mol).
    excess_v
        Mixture molar volume minus the mole-fraction average (cm3/mol).

    Returns
    -------
    tuple[pd.DataFrame, ...]
        ``runs`` and ``systems`` tables for :func:`mixture_points`.
    """
    h = {"a": -45.3, "b": -61.7}
    v = {"a": 58.4, "b": 147.5}
    common = {"temperature": 298.15, "h_sigma": 0.01, "v_molar_sigma": 0.02}
    rows = []
    for rep in range(3):
        rows += [
            {
                "box_id": "pure_a",
                "kind": "pure",
                "smiles": "CCO",
                "x": 1.0,
                "replicate": rep,
                "h_mean": h["a"],
                "v_molar_mean": v["a"],
                "density_mean": 0.785,
                "density_sigma": 1e-4,
                **common,
            },
            {
                "box_id": "pure_b",
                "kind": "pure",
                "smiles": "CCCCCCC",
                "x": 1.0,
                "replicate": rep,
                "h_mean": h["b"],
                "v_molar_mean": v["b"],
                "density_mean": 0.680,
                "density_sigma": 1e-4,
                **common,
            },
            {
                "box_id": "test_x300",
                "kind": "mixture",
                "smiles": "CCO;CCCCCCC",
                "x": x,
                "replicate": rep,
                "h_mean": x * h["a"] + (1 - x) * h["b"] + excess_h,
                "v_molar_mean": x * v["a"] + (1 - x) * v["b"] + excess_v,
                "density_mean": 0.72,
                "density_sigma": 1e-4,
                **common,
            },
        ]
    systems = pd.DataFrame(
        [
            {
                "system_id": "test",
                "x_box": x,
                "a": "CCO",
                "b": "CCCCCCC",
                "hmix_exp": 0.0,
                "hmix_sigma": 0.01,
                "density_exp": 0.72,
                "density_sigma": 1e-4,
                "vex_exp": 0.0,
            }
        ]
    )
    return pd.DataFrame(rows), systems


@pytest.mark.parametrize(("excess_h", "excess_v"), [(0.0, 0.0), (0.7, -0.35)])
def test_mixture_points_recovers_excess(excess_h, excess_v):
    """Ideal mixing gives zero ΔHmix and V^E; a known excess is recovered."""
    runs, systems = _runs(0.3, excess_h, excess_v)
    points = mixture_points(runs, systems)
    assert len(points) == 3
    assert points.hmix.to_numpy() == pytest.approx(excess_h, abs=1e-12)
    assert points.vex.to_numpy() == pytest.approx(excess_v, abs=1e-12)


def _published_density(x):
    """
    Return a smooth, curved density for the finite-size test.

    Parameters
    ----------
    x
        Mole fraction of component a.

    Returns
    -------
    float
        Density in g/mL (ideal mixing of 0.885 and 0.699 g/mL liquids with
        molar masses 116 and 114 g/mol, plus a V^E term).
    """
    mass = 116.16 * x + 114.23 * (1 - x)
    volume = 116.16 * x / 0.885 + 114.23 * (1 - x) / 0.699 + x * (1 - x) * 0.8
    return mass / volume


def _published_hmix(x):
    """
    Return a skewed excess enthalpy for the finite-size test.

    Parameters
    ----------
    x
        Mole fraction of component a.

    Returns
    -------
    float
        ΔHmix in kJ/mol.
    """
    return x * (1 - x) * (5.6 - 1.2 * (2 * x - 1))


def test_finite_size_check_corrects_composition_offset():
    """
    Published values at nearby compositions are moved to our composition.

    Our values equal the published curve exactly at our compositions, which
    are up to 0.0165 from the published ones (as for methylpentanoate +
    octane), so every difference must vanish. Without the shift the density
    gap would exceed ``SIZE_GAP``.
    """
    x_pub = np.array([0.2465, 0.4896, 0.7399])
    x_ours = np.array([0.23, 0.49, 0.73])
    published = pd.DataFrame(
        [
            {
                "smiles_1": "CCCCC(=O)OC",
                "smiles_2": "CCCCCCCC",
                "x_1": x,
                "prop": prop,
                "value": f(x),
                "uncertainty": 0.001,
            }
            for x in x_pub
            for prop, f in (("density", _published_density), ("hmix", _published_hmix))
        ]
    )
    # Published in the other orientation, to test the swap as well
    published.loc[published.prop == "hmix", ["smiles_1", "smiles_2"]] = [
        "CCCCCCCC",
        "CCCCC(=O)OC",
    ]
    published.loc[published.prop == "hmix", "x_1"] = 1 - x_pub
    systems = pd.DataFrame(
        {"system_id": "mo", "a": "CCCCC(=O)OC", "b": "CCCCCCCC", "x_box": x_ours}
    )
    summary = pd.DataFrame(
        {
            "box_id": [f"mo_x{round(x * 1000):03d}" for x in x_ours],
            "x": x_ours,
            "n_rep": 5,
            "hmix_sage": [_published_hmix(x) for x in x_ours],
            "hmix_sd": 0.05,
            "density_sage": [_published_density(x) for x in x_ours],
            "density_sd": 0.0005,
        }
    )
    finite = finite_size_check(summary, systems, published)
    assert len(finite) == 6
    density = finite[finite.prop == "density"]
    assert density["diff"].abs().max() < 1e-4
    assert finite[finite.prop == "hmix"]["diff"].abs().max() < 5e-3
    raw_gap = (density.ours - density.published_nearest).abs().max()
    assert raw_gap > SIZE_GAP["density"]


def test_pooled_sigma_floor_exception():
    """Listed systems get the floor even when their own spread is larger."""
    rng = np.random.default_rng(0)
    rows = []
    for system_id, scale in (
        ("quiet", 0.05),
        ("noisy", 0.1),
        ("nmp_dichloroethane", 0.5),
    ):
        for replicate in range(5):
            value = scale * rng.standard_normal()
            rows.append(
                {
                    "system_id": system_id,
                    "box_id": f"{system_id}_x500",
                    "x": 0.5,
                    "replicate": replicate,
                    "hmix": value,
                    "hmix_sigma_run": 0.01,
                    "density": 1.0 + value / 100,
                    "vex": 0.0,
                    "hmix_exp": 0.0,
                    "hmix_exp_sigma": 0.0,
                    "density_exp": 1.0,
                    "density_exp_sigma": 0.0,
                    "vex_exp": 0.0,
                }
            )
    summary = pooled_sigma(pd.DataFrame(rows), 5).set_index("system_id")
    for prop in ("hmix", "density"):
        sd = summary[f"{prop}_sd"]
        floor = float(np.median(sd))
        assert summary.loc["nmp_dichloroethane", f"sigma_sage_{prop}"] == (
            pytest.approx(floor, abs=1e-6)
        )
        assert sd["nmp_dichloroethane"] > floor
        assert summary.loc["quiet", f"sigma_sage_{prop}"] == pytest.approx(
            floor, abs=1e-6
        )


def _version_summary(shift: float, sd: float) -> pd.DataFrame:
    """
    Two-point summary in the format of ``pooled_sigma``.

    Parameters
    ----------
    shift
        Added to every Sage value.
    sd
        Replicate SD of every property.

    Returns
    -------
    pd.DataFrame
        Summary rows for systems ``a`` and ``b``.
    """
    rows = []
    for system_id, exp in (("a", 1.0), ("b", -1.0)):
        row = {"system_id": system_id, "box_id": f"{system_id}_x500", "x": 0.5}
        row["n_rep"] = 5
        for prop in ("hmix", "density", "vex"):
            row[f"{prop}_sage"] = exp + 0.5 + shift
            row[f"{prop}_sd"] = sd
            row[f"{prop}_exp"] = exp
        rows.append(row)
    return pd.DataFrame(rows)


def test_compare_versions():
    """Differences, Welch t and MAE follow from the replicate means and SDs."""
    ref = _version_summary(0.0, 0.1)
    other = _version_summary(-0.3, 0.2)
    table = compare_versions(ref, other)
    assert len(table) == 6
    np.testing.assert_allclose(table["diff"], -0.3)
    se = np.sqrt((0.1**2 + 0.2**2) / 5)
    np.testing.assert_allclose(table.se, se)
    np.testing.assert_allclose(table.t, -0.3 / se)
    assert (table.p_adj >= table.p).all()

    mae = version_mae(table, {"a": "training", "b": "validation"})
    total = mae[mae.set == "all"].set_index("prop")
    np.testing.assert_allclose(total.mae_ref, 0.5)
    np.testing.assert_allclose(total.mae_other, 0.2)
    assert (total.n_better == 2).all()
    assert set(mae.set) == {"training", "validation", "all"}
