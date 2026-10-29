"""Test the power analysis of the binary mixtures benchmark."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from liquid_mixtures_inputs.power_analysis import (
    REFERENCE,
    TOTAL_PS,
    allocate,
    metric_noise,
    point_sigma,
    truncated_estimates,
)
from ml_peg.analysis.molecular_dynamics.binary_mixtures.mixing import (
    box_result,
    read_production,
)


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


@pytest.fixture
def scalar_log(tmp_path, request):
    """
    Write a synthetic scalar log in the calc's CSV format.

    Parameters
    ----------
    tmp_path
        Temporary directory.
    request
        Request; ``request.param`` is the AR(1) coefficient (default 0.9).

    Returns
    -------
    Path
        CSV path.
    """
    phi = getattr(request, "param", 0.9)
    rng = np.random.default_rng(0)
    time = np.arange(1, int(TOTAL_PS * 10) + 1) * 0.1
    log = pd.DataFrame(
        {
            "step": np.arange(1, len(time) + 1) * 100,
            "time_ps": time,
            "epot_eV": -40.0 + 0.5 * ar1(len(time), phi, rng),
            "ekin_eV": 50.0,
            "temperature_K": 298.0,
            "volume_A3": 16000.0 + 20 * ar1(len(time), phi, rng),
            "density_g_cm3": 0.9,
        }
    )
    path = tmp_path / "rep0.csv"
    log.to_csv(path, index=False)
    return path


def test_reference_window_matches_benchmark(scalar_log):
    """The 0.2 + 2 ns window gives the benchmark analysis's means and sems."""
    rows = pd.DataFrame(truncated_estimates(("box", 0, 120, scalar_log)))
    ours = rows[
        (rows.equilibration_ps == REFERENCE[0]) & (rows.production_ps == REFERENCE[1])
    ].set_index("observable")
    prod = read_production(scalar_log, REFERENCE[0], TOTAL_PS)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        bench = box_result(prod, 120)
    for name, estimate in bench.estimates.items():
        assert ours.loc[name, "mean"] == pytest.approx(estimate.mean, rel=1e-12)
        assert ours.loc[name, "sem"] == pytest.approx(estimate.sem, rel=1e-12)


@pytest.mark.parametrize("scalar_log", [0.5], indirect=True)
def test_sem_scales_as_inverse_sqrt_time(scalar_log):
    """For a short-memory series, sem(0.25 ns) / sem(2 ns) ≈ sqrt(8)."""
    rows = pd.DataFrame(truncated_estimates(("box", 0, 120, scalar_log)))
    h = rows[(rows.observable == "h") & (rows.equilibration_ps == REFERENCE[0])]
    sem = h.set_index("production_ps")["sem"]
    assert sem[250.0] / sem[2000.0] == pytest.approx(np.sqrt(8), rel=0.25)


def test_incomplete_run_is_skipped(tmp_path):
    """A log that stops early gives no rows."""
    path = tmp_path / "rep0.csv"
    pd.DataFrame({"time_ps": [0.1, 1000.0]}).to_csv(path, index=False)
    assert truncated_estimates(("box", 0, 10, path)) == []


def test_point_sigma_matches_sampled_covariance():
    """sqrt(diag(M Σ Mᵀ)) matches the SD of sampled point noise."""
    rng = np.random.default_rng(2)
    matrix = np.array([[1.0, -0.2, -0.8, 0.0], [0.0, -0.5, -0.5, 1.0]])
    sigma_box = np.array([0.1, 0.3, 0.2, 0.15])
    draws = (rng.standard_normal((200_000, 4)) * sigma_box) @ matrix.T
    np.testing.assert_allclose(
        point_sigma(matrix, sigma_box), draws.std(axis=0), rtol=0.01
    )


def test_metric_noise_matches_closed_form():
    """With errors far above the noise, MAE noise matches sqrt(sᵀCs)/N."""
    rng = np.random.default_rng(3)
    n_points, n_boxes = 30, 40
    matrix = np.zeros((n_points, n_boxes))
    matrix[np.arange(n_points), np.arange(n_points)] = 1.0
    # Points share pure boxes in groups of three, as in the benchmark
    for i in range(n_points):
        matrix[i, n_points + i // 3] = -0.5
    prop = {
        "matrix": matrix,
        "error": rng.choice([-1.0, 1.0], n_points) * rng.uniform(5, 10, n_points),
        "metric": "mae",
    }
    sigma_box = rng.uniform(0.05, 0.2, n_boxes)
    row = metric_noise(prop, sigma_box, 1.0, rng)
    assert row["sd"] == pytest.approx(row["sd_closed_form"], rel=0.03)
    assert abs(row["bias"]) < 3 * row["sd"] / np.sqrt(1000)


def test_metric_noise_floor_for_perfect_model():
    """A perfect model's MAE is the mean of σ·sqrt(2/π) over points."""
    rng = np.random.default_rng(4)
    sigma_box = np.array([0.1, 0.2, 0.3])
    prop = {"matrix": np.eye(3), "error": np.ones(3), "metric": "mae"}
    row = metric_noise(prop, sigma_box, 0.0, rng)
    assert row["observed_mean"] == pytest.approx(
        np.mean(sigma_box) * np.sqrt(2 / np.pi), rel=0.02
    )


def test_allocation_reaches_target_at_lowest_cost():
    """The allocation hits the target variance, keeps the floor and beats uniform."""
    coefficient = np.array([1.0, 4.0, 9.0, 0.01])
    n_atoms = np.array([1000.0, 2000.0, 1500.0, 1200.0])
    uniform_t = 2000.0
    target = np.sum(coefficient / uniform_t)
    t = allocate(coefficient, n_atoms, 500.0, target)
    assert np.sum(coefficient / t) == pytest.approx(target, rel=1e-9)
    assert t.min() >= 500.0 - 1e-9
    assert np.sum(n_atoms * t) < np.sum(n_atoms * uniform_t)
    free = t > 500.0 + 1e-9
    ratio = t[free] / np.sqrt(coefficient[free] / n_atoms[free])
    np.testing.assert_allclose(ratio, ratio[0], rtol=1e-9)
