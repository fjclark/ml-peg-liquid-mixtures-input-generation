"""Test the Sage summary statistics."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from liquid_mixtures_inputs.summarise_sage import within_between
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
