"""
Test whether one MLP run per box can rank models against experiment.

Run in the ``analysis`` pixi environment after ``summarise-sage``::

    pixi run -e analysis power-analysis \
        --data-dir D --runs-dir R --output-dir O

The Sage replicates stand in for MLP simulation noise. Every Sage run is
re-analysed at shorter production lengths with the benchmark's estimator, giving
the single-run noise of each box against run length. That noise is propagated
through the linear map from box means to ΔHmix and V^E (``mixing_matrix``; the
pure boxes are shared between a system's three compositions, so their noise is
correlated across points) to give, per run length:

- the smallest error against experiment a single point can detect;
- the noise in each benchmark metric (MAE, MAZE) and the smallest
  difference between two models it can resolve;
- the spread of Sage's metric over the choice of systems (cluster bootstrap),
  which shows how far a score depends on this set of systems;
- a GPU-time allocation across boxes that keeps the metric noise of uniform
  runs of the protocol length at a lower cost;
- ml-peg's 'good' and 'bad' thresholds (``thresholds.csv`` and
  ``metrics_snippet.yml``).

The reports are written by ``power_report.py``.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import warnings

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
import yaml  # noqa: E402

from liquid_mixtures_inputs.constants import (  # noqa: E402
    EQUILIBRATION_PS,
    PRODUCTION_PS,
    SAGE_RUN_PS,
)
from ml_peg.analysis.molecular_dynamics.binary_mixtures.mixing import (  # noqa: E402
    mixing_matrix,
    observables,
)
from ml_peg.analysis.utils.timeseries import estimate_mean  # noqa: E402

REFERENCE = (EQUILIBRATION_PS, PRODUCTION_PS)  # the benchmark protocol
TOTAL_PS = EQUILIBRATION_PS + PRODUCTION_PS
ALTERNATIVE_PS = 2000.0  # the longer production compared with the protocol
# Every Sage run is re-analysed for each pair of these lengths.
EQUILIBRATION_LENGTHS_PS = (100.0, 200.0)
PRODUCTION_LENGTHS_PS = (250.0, 500.0, 1000.0, 1500.0, 2000.0)
EARLY_WINDOW = (100.0, 100.0)  # 100 to 200 ps, for the equilibration test
OBSERVABLES = ("h", "density", "v_molar")
# Two-sided 5 % test with 80 % power
Z_DETECT = stats.norm.ppf(0.975) + stats.norm.ppf(0.8)
# Error profiles: multiples of Sage's replicate-mean errors
PROFILES = (1.0, 0.5, 0.25, 0.0)
N_DRAWS = 20000
N_BOOTSTRAP = 20000
MIN_ALLOCATION_PS = 500.0
# Categorical slots 1-3 of the dataviz reference palette (pass all-pairs CVD checks)
COLORS = ("#2a78d6", "#eb6834", "#1baf7a")
# ml-peg thresholds (DECISIONS.md, metrics): 'good' is what a model that matches
# experiment exactly scores with one run per box, plus GOOD_NOISE_SDS standard
# deviations of that score; 'bad' is BAD_SAGE_MULTIPLE times Sage's score, as
# Sage was trained on these properties.
GOOD_NOISE_SDS = 2.0
BAD_SAGE_MULTIPLE = 5.0
# Properties that share one set of thresholds (DECISIONS.md, metrics): mixture and
# pure densities, for simplicity. The shared value is the larger of the two, so
# the rules above still hold for each.
THRESHOLD_GROUPS = {"pure_density": "density"}
LABELS = {
    "hmix": "ΔHmix",
    "density": "Mixture density",
    "vex": "V^E",
    "pure_density": "Pure density",
}


def truncated_estimates(task: tuple) -> list[dict]:
    """
    Re-analyse one Sage run at every equilibration and production length.

    Parameters
    ----------
    task
        ``(box_id, replicate, n_molecules, csv_path)``.

    Returns
    -------
    list[dict]
        One row per (equilibration, production, observable), or none if the run
        is incomplete.
    """
    box_id, replicate, n_molecules, csv_path = task
    data = pd.read_csv(csv_path)
    if data.time_ps.iloc[-1] < SAGE_RUN_PS - 0.5:
        return []
    rows = []
    windows = [
        (equilibration, production)
        for equilibration in EQUILIBRATION_LENGTHS_PS
        for production in PRODUCTION_LENGTHS_PS
    ]
    for equilibration, production in [*windows, EARLY_WINDOW]:
        window = data[
            (data.time_ps > equilibration)
            & (data.time_ps <= equilibration + production + 1e-6)
        ]
        for name, series in observables(window, n_molecules).items():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                estimate = estimate_mean(series)
            rows.append(
                {
                    "box_id": box_id,
                    "replicate": replicate,
                    "equilibration_ps": equilibration,
                    "production_ps": production,
                    "observable": name,
                    "mean": estimate.mean,
                    "sem": estimate.sem,
                    "n_eff": estimate.n_eff,
                }
            )
    return rows


def load_truncated(runs: pd.DataFrame, runs_dir: Path, workers: int) -> pd.DataFrame:
    """
    Re-analyse every Sage run at every run length.

    Parameters
    ----------
    runs
        ``runs.csv`` from ``summarise-sage`` (box, replicate, molecule count).
    runs_dir
        Output directory of ``run_sage_benchmark.py``.
    workers
        Number of processes.

    Returns
    -------
    pd.DataFrame
        One row per run, length and observable.
    """
    tasks = [
        (
            r.box_id,
            r.replicate,
            r.n_molecules,
            runs_dir / r.box_id / f"rep{r.replicate}.csv",
        )
        for r in runs.itertuples()
    ]
    with ProcessPoolExecutor(workers) as pool:
        return pd.DataFrame(
            [row for rows in pool.map(truncated_estimates, tasks) for row in rows]
        )


def box_noise(truncated: pd.DataFrame) -> pd.DataFrame:
    """
    Single-run noise of every box at every run length.

    ``sem`` is the RMS within-run standard error over replicates, ``sd`` the SD
    of the replicate means. The pooled ratio sqrt(Σ(n-1)·SD²/sem² / Σ(n-1))
    over boxes, per length and observable, is the inflation ``k`` that turns a
    within-run sem into the run-to-run spread; ``sigma`` = k·sem is the
    single-run noise used downstream (smoother than the 5-replicate SD).

    Parameters
    ----------
    truncated
        Output of :func:`load_truncated`.

    Returns
    -------
    pd.DataFrame
        One row per box, length and observable.
    """
    keys = ["equilibration_ps", "production_ps", "observable"]
    noise = (
        truncated.groupby([*keys, "box_id"])
        .agg(
            n=("mean", "count"),
            mean=("mean", "mean"),
            sd=("mean", lambda s: s.std(ddof=1)),
            sem=("sem", lambda s: float(np.sqrt(np.mean(s**2)))),
        )
        .reset_index()
    )
    chi2 = (noise["n"] - 1) * (noise["sd"] / noise["sem"]) ** 2
    pooled = (
        noise.assign(chi=chi2, dof=noise["n"] - 1)
        .groupby(keys)
        .apply(
            lambda group: float(np.sqrt(group.chi.sum() / group.dof.sum())),
            include_groups=False,
        )
        .rename("k")
    )
    noise = noise.join(pooled, on=keys)
    noise["sigma"] = noise["k"] * noise["sem"]
    return noise


def early_shifts(truncated: pd.DataFrame) -> pd.Series:
    """
    Compute each run's 100 to 200 ps window mean minus its production mean.

    Parameters
    ----------
    truncated
        Output of :func:`load_truncated`.

    Returns
    -------
    pd.Series
        Shift indexed by box, replicate and observable.
    """
    index = ["box_id", "replicate", "observable"]
    early = truncated[
        (truncated.equilibration_ps == EARLY_WINDOW[0])
        & (truncated.production_ps == EARLY_WINDOW[1])
    ].set_index(index)
    full = truncated[
        (truncated.equilibration_ps == REFERENCE[0])
        & (truncated.production_ps == REFERENCE[1])
    ].set_index(index)
    return (early["mean"] - full["mean"]).dropna()


def equilibration_shift(truncated: pd.DataFrame) -> pd.DataFrame:
    """
    Test whether the 100 to 200 ps window differs from the production after it.

    If it does not, 100 ps of equilibration would do for Sage. The two windows
    do not overlap. The test is a one-sample t-test of the per-run shifts
    (window mean minus production mean) over all runs, so it does not rely on
    the within-run sem of a 100 ps window, which is too small.

    Parameters
    ----------
    truncated
        Output of :func:`load_truncated`.

    Returns
    -------
    pd.DataFrame
        One row per observable: number of runs, mean and RMS shift, and the
        two-sided p-value of the mean shift.
    """
    shift = early_shifts(truncated)
    rows = []
    for observable, values in shift.groupby(level="observable"):
        rows.append(
            {
                "observable": observable,
                "n_runs": len(values),
                "mean_shift": float(values.mean()),
                "rms_shift": float(np.sqrt((values**2).mean())),
                "p": float(stats.ttest_1samp(values, 0.0).pvalue),
            }
        )
    return pd.DataFrame(rows)


def metric_sets(data_dir: Path, runs: pd.DataFrame) -> dict[str, dict]:
    """
    Points, box map and Sage errors of every scored property.

    Parameters
    ----------
    data_dir
        Benchmark data directory with the Sage columns (``systems.csv``,
        ``boxes.csv``, ``pure.csv``).
    runs
        ``runs.csv`` from ``summarise-sage``.

    Returns
    -------
    dict[str, dict]
        Per property: ``matrix`` (points × boxes), ``observable`` (box
        observable it maps), ``error`` (Sage replicate mean minus experiment),
        ``sigma_exp``, ``sigma_ref`` (MAZE denominator), ``system`` (cluster
        label per point) and ``unit``.
    """
    systems = pd.read_csv(data_dir / "systems.csv")
    boxes = pd.read_csv(data_dir / "boxes.csv").sort_values("box_id")
    pure = pd.read_csv(data_dir / "pure.csv")
    box_ids = list(boxes.box_id)
    column = {box_id: i for i, box_id in enumerate(box_ids)}
    mixing = mixing_matrix(systems, boxes, box_ids)

    with_density = systems[systems.density_exp.notna()]
    density = np.zeros((len(with_density), len(box_ids)))
    density[np.arange(len(with_density)), [column[b] for b in with_density.box_id]] = 1

    sage_density = (
        runs[runs.kind == "pure"]
        .groupby(["smiles", "temperature"])["density_mean"]
        .mean()
    )
    pure_ids = boxes[boxes.kind == "pure"].set_index(["smiles", "temperature"]).box_id
    scored = pure[pure.scored].reset_index(drop=True)
    pure_matrix = np.zeros((len(scored), len(box_ids)))
    for i, row in enumerate(scored.itertuples()):
        pure_matrix[i, column[pure_ids.loc[(row.smiles, row.temperature)]]] = 1
    pure_error = np.array(
        [
            sage_density.loc[(r.smiles, r.temperature)] - r.density
            for r in scored.itertuples()
        ]
    )
    return {
        "box_ids": box_ids,
        "hmix": {
            "matrix": mixing,
            "observable": "h",
            "error": (systems.hmix_sage - systems.hmix_exp).to_numpy(),
            "sigma_exp": systems.hmix_sigma.to_numpy(),
            "sigma_ref": systems.sigma_ref_hmix.to_numpy(),
            "system": systems.system_id.to_numpy(),
            "unit": "kJ/mol",
        },
        "density": {
            "matrix": density,
            "observable": "density",
            "error": (with_density.density_sage - with_density.density_exp).to_numpy(),
            "sigma_exp": with_density.density_sigma.to_numpy(),
            "sigma_ref": with_density.sigma_ref_density.to_numpy(),
            "system": with_density.system_id.to_numpy(),
            "unit": "g/mL",
        },
        "vex": {
            "matrix": mixing[systems.density_exp.notna().to_numpy()],
            "observable": "v_molar",
            "error": (with_density.vex_sage - with_density.vex_exp).to_numpy(),
            "sigma_exp": with_density.vex_sigma.fillna(0.0).to_numpy(),
            "sigma_ref": with_density.sigma_ref_vex.to_numpy(),
            "system": with_density.system_id.to_numpy(),
            "unit": "cm³/mol",
        },
        "pure_density": {
            "matrix": pure_matrix,
            "observable": "density",
            "error": pure_error,
            "sigma_exp": scored.sigma.to_numpy(),
            "sigma_ref": scored.sigma_ref_density.to_numpy(),
            "system": scored.smiles.to_numpy(),
            "unit": "g/mL",
        },
    }


def metric(errors: np.ndarray, scale: np.ndarray | None = None) -> np.ndarray:
    """
    Benchmark metric (MAE, or MAZE with ``scale``) over the last axis.

    Parameters
    ----------
    errors
        Model minus experiment, shape (..., n_points).
    scale
        Per-point denominator σ_ref (MAZE), or None (MAE).

    Returns
    -------
    np.ndarray
        Metric with the last axis reduced.
    """
    if scale is not None:
        errors = errors / scale
    return np.mean(np.abs(errors), axis=-1)


def point_sigma(matrix: np.ndarray, sigma_box: np.ndarray) -> np.ndarray:
    """
    Single-run noise of each point, from independent box noise.

    Parameters
    ----------
    matrix
        Points × boxes linear map.
    sigma_box
        Single-run noise of every box.

    Returns
    -------
    np.ndarray
        Point SDs sqrt(diag(M Σ Mᵀ)) with Σ = diag(sigma_box²).
    """
    return np.sqrt((matrix**2) @ sigma_box**2)


def metric_noise(
    property_set: dict,
    sigma_box: np.ndarray,
    profile: float,
    rng: np.random.Generator,
    scale: np.ndarray | None = None,
) -> dict:
    """
    Monte Carlo noise of a metric for one model with one run per box.

    The model's true errors are ``profile`` × Sage's errors; each draw adds
    correlated point noise M·(σ_box·z).

    Parameters
    ----------
    property_set
        One entry of :func:`metric_sets`.
    sigma_box
        Single-run noise of every box for ``prop["observable"]``.
    profile
        Multiple of Sage's errors used as the model's true errors.
    rng
        Random generator.
    scale
        Per-point denominator (MAZE), or None.

    Returns
    -------
    dict
        True metric, mean and SD of the observed metric, bias, the closed-form
        SD (sqrt(sᵀCs)/N with s = sign(error) / scale, valid when errors
        are well above the noise) and the resolvable difference between two
        independent models, Z_DETECT·√2·SD.
    """
    matrix, error = property_set["matrix"], profile * property_set["error"]
    noise = (rng.standard_normal((N_DRAWS, len(sigma_box))) * sigma_box) @ matrix.T
    observed = metric(error + noise, scale)
    true = float(metric(error, scale))
    sd = float(observed.std(ddof=1))
    closed = np.nan
    if profile > 0:
        weights = np.sign(error) / (scale if scale is not None else 1.0)
        covariance = matrix @ np.diag(sigma_box**2) @ matrix.T
        closed = float(np.sqrt(weights @ covariance @ weights) / len(error))
    return {
        "true": true,
        "observed_mean": float(observed.mean()),
        "bias": float(observed.mean()) - true,
        "sd": sd,
        "sd_closed_form": closed,
        "resolvable": float(Z_DETECT * np.sqrt(2) * sd),
    }


def system_bootstrap(
    property_set: dict, rng: np.random.Generator, scale: np.ndarray | None = None
) -> float:
    """
    SD of Sage's noiseless metric over resampled systems (cluster bootstrap).

    Parameters
    ----------
    property_set
        One entry of :func:`metric_sets`.
    rng
        Random generator.
    scale
        Per-point denominator (MAZE), or None.

    Returns
    -------
    float
        Bootstrap SD of the metric.
    """
    labels = property_set["system"]
    clusters = [np.flatnonzero(labels == system) for system in np.unique(labels)]
    error = property_set["error"]
    if scale is not None:
        error = error / scale
    values = np.empty(N_BOOTSTRAP)
    for i in range(N_BOOTSTRAP):
        picked = rng.integers(len(clusters), size=len(clusters))
        values[i] = metric(np.concatenate([error[clusters[j]] for j in picked]))
    return float(values.std(ddof=1))


def allocate(
    coefficient: np.ndarray,
    n_atoms: np.ndarray,
    min_ps: float,
    target: float,
) -> np.ndarray:
    """
    Production lengths minimising GPU time for a target metric variance.

    The variance is Σ_b A_b / t_b (single-run variance ∝ 1/t) and the cost
    Σ_b n_b·t_b. Without bounds the optimum is t_b ∝ sqrt(A_b / n_b); boxes that
    would fall below ``min_ps`` are fixed there and the rest re-solved.

    Parameters
    ----------
    coefficient
        A_b, the variance contribution of box b times its run length.
    n_atoms
        Atoms per box (cost per ps).
    min_ps
        Shortest allowed production.
    target
        Variance to reach.

    Returns
    -------
    np.ndarray
        Production length per box in ps.
    """
    fixed = np.zeros(len(coefficient), dtype=bool)
    while True:
        remaining = target - np.sum(coefficient[fixed] / min_ps)
        free = ~fixed
        shape = np.sqrt(coefficient[free] / n_atoms[free])
        # Σ A/(λ·shape) = remaining  ⇒  λ = Σ (A/shape) / remaining
        multiplier = np.sum(coefficient[free] / shape) / remaining
        lengths = np.full(len(coefficient), min_ps)
        lengths[free] = multiplier * shape
        below_minimum = free & (lengths < min_ps)
        if not below_minimum.any():
            return lengths
        fixed |= below_minimum


def sigma_at(
    noise: pd.DataFrame, box_ids: list[str], observable: str, production: float
) -> np.ndarray:
    """
    Single-run noise of every box for one observable and production length.

    Parameters
    ----------
    noise
        Output of :func:`box_noise`.
    box_ids
        Box order.
    observable
        Observable.
    production
        Production length in ps (with the reference equilibration).

    Returns
    -------
    np.ndarray
        Single-run noise σ per box, in ``box_ids`` order.
    """
    rows = noise[
        (noise.equilibration_ps == REFERENCE[0])
        & (noise.production_ps == production)
        & (noise.observable == observable)
    ].set_index("box_id")
    return rows.sigma.loc[box_ids].to_numpy()


def analyse(
    sets: dict, noise: pd.DataFrame, n_atoms: np.ndarray, seed: int = 0
) -> dict:
    """
    Per-point detection limits, metric noise and the time allocation.

    Parameters
    ----------
    sets
        Output of :func:`metric_sets`.
    noise
        Output of :func:`box_noise`.
    n_atoms
        Atoms per box, in ``sets["box_ids"]`` order.
    seed
        Random seed.

    Returns
    -------
    dict
        Tables ``points``, ``metrics``, ``bootstrap`` and ``allocation``.
    """
    rng = np.random.default_rng(seed)
    box_ids = sets["box_ids"]
    property_sets = {name: value for name, value in sets.items() if name != "box_ids"}
    point_rows, metric_rows, bootstrap_rows = [], [], []
    for name, property_set in property_sets.items():
        # Every property is scored as MAE (physical units) and MAZE (σ_ref).
        variants = [("MAE", None), ("MAZE", property_set["sigma_ref"])]
        for label, scale in variants:
            bootstrap_rows.append(
                {
                    "property": name,
                    "metric": label,
                    "sage": float(metric(property_set["error"], scale)),
                    "bootstrap_sd": system_bootstrap(property_set, rng, scale),
                }
            )
        for production in PRODUCTION_LENGTHS_PS:
            sigma_box = sigma_at(noise, box_ids, property_set["observable"], production)
            sigma_sim = point_sigma(property_set["matrix"], sigma_box)
            sigma_total = np.sqrt(sigma_sim**2 + property_set["sigma_exp"] ** 2)
            for i in range(len(sigma_sim)):
                point_rows.append(
                    {
                        "property": name,
                        "point": i,
                        "system_id": property_set["system"][i],
                        "production_ps": production,
                        "sigma_sim": sigma_sim[i],
                        "sigma_exp": property_set["sigma_exp"][i],
                        "mde": Z_DETECT * sigma_total[i],
                        "sage_error": property_set["error"][i],
                    }
                )
            for label, scale in variants:
                for profile in PROFILES:
                    row = metric_noise(property_set, sigma_box, profile, rng, scale)
                    metric_rows.append(
                        {
                            "property": name,
                            "metric": label,
                            "production_ps": production,
                            "profile": profile,
                            **row,
                        }
                    )
    return {
        "points": pd.DataFrame(point_rows),
        "metrics": pd.DataFrame(metric_rows),
        "bootstrap": pd.DataFrame(bootstrap_rows),
        "allocation": allocation_table(property_sets, noise, box_ids, n_atoms),
    }


def allocation_table(
    property_sets: dict, noise: pd.DataFrame, box_ids: list[str], n_atoms: np.ndarray
) -> pd.DataFrame:
    """
    Choose per-box production lengths that keep the noise of uniform runs.

    ΔHmix MAE, mixture density MAE and pure density MAE are weighted together,
    each variance normalised by its value for uniform runs of the protocol
    length. With random error signs, E[(Σ_i s_i M_ib)²] = Σ_i M_ib², so
    A_b = Σ_i M_ib² c_b / N².

    Parameters
    ----------
    property_sets
        Entries of :func:`metric_sets` by property.
    noise
        Output of :func:`box_noise`.
    box_ids
        Box order.
    n_atoms
        Atoms per box.

    Returns
    -------
    pd.DataFrame
        Weight and production length per box.
    """
    total = np.zeros(len(box_ids))
    for name in ("hmix", "density", "pure_density"):
        property_set = property_sets[name]
        # Single-run variance × run length, which is constant if σ ∝ 1/sqrt(t).
        variance_times_length = (
            sigma_at(noise, box_ids, property_set["observable"], REFERENCE[1]) ** 2
            * REFERENCE[1]
        )
        contribution = (
            (property_set["matrix"] ** 2).sum(axis=0)
            * variance_times_length
            / len(property_set["error"]) ** 2
        )
        total += contribution / np.sum(contribution / REFERENCE[1])
    used = total > 0
    lengths = np.full(len(box_ids), np.nan)
    lengths[used] = allocate(
        total[used], n_atoms[used], MIN_ALLOCATION_PS, np.sum(total) / REFERENCE[1]
    )
    return pd.DataFrame(
        {
            "box_id": box_ids,
            "n_atoms": n_atoms,
            "weight": total,
            "production_ps": lengths,
        }
    )


def round_significant(value: float, digits: int = 2) -> float:
    """
    Round a number to a few significant figures.

    Parameters
    ----------
    value
        Number to round.
    digits
        Significant figures.

    Returns
    -------
    float
        The rounded number.
    """
    return float(f"{value:.{digits}g}")


def thresholds(metrics: pd.DataFrame, bootstrap: pd.DataFrame) -> pd.DataFrame:
    """
    Set ml-peg's 'good' and 'bad' thresholds of every metric from the noise.

    good = perfect-model mean + ``GOOD_NOISE_SDS`` × its SD with one run per
    box of the benchmark protocol, so a model that matches experiment scores
    'good' about 98 % of the time. bad = ``BAD_SAGE_MULTIPLE`` × Sage's score.
    Properties in ``THRESHOLD_GROUPS`` share the larger of the group's values.

    Parameters
    ----------
    metrics
        ``metrics`` table from :func:`analyse`.
    bootstrap
        ``bootstrap`` table from :func:`analyse` (holds Sage's score).

    Returns
    -------
    pd.DataFrame
        One row per metric with the perfect-model mean and SD, Sage's score
        and the two thresholds (2 significant figures).
    """
    perfect = metrics[
        (metrics.production_ps == REFERENCE[1]) & (metrics.profile == 0.0)
    ].set_index(["property", "metric"])
    rows = []
    for row in bootstrap.itertuples():
        noise = perfect.loc[(row.property, row.metric)]
        rows.append(
            {
                "property": row.property,
                "metric": row.metric,
                "perfect_mean": noise.observed_mean,
                "perfect_sd": noise.sd,
                "sage": row.sage,
                "good": round_significant(
                    noise.observed_mean + GOOD_NOISE_SDS * noise.sd
                ),
                "bad": round_significant(BAD_SAGE_MULTIPLE * row.sage),
            }
        )
    table = pd.DataFrame(rows)
    group = table.property.replace(THRESHOLD_GROUPS)
    for column in ("good", "bad"):
        table[column] = table.groupby([group, "metric"])[column].transform("max")
    return table


# ml-peg metric names, units and what each point set holds.
MLPEG_NAMES = {
    "hmix": ("ΔHmix", "kJ mol^-1", "enthalpy of mixing"),
    "density": ("Density", "g cm^-3", "mixture densities"),
    "vex": ("V^E", "cm^3 mol^-1", "excess molar volume"),
    "pure_density": (
        "Pure density",
        "g cm^-3",
        "pure-liquid densities not already scored by the liquid densities benchmark",
    ),
}


def metrics_snippet(table: pd.DataFrame, sets: dict) -> str:
    """
    Write the thresholds in the format of ml-peg's ``metrics.yml``.

    MAZE is scored; MAE is shown with weight 0, so each score keeps its
    physical units without counting twice.

    Parameters
    ----------
    table
        Output of :func:`thresholds`.
    sets
        Output of :func:`metric_sets` (for the number of points).

    Returns
    -------
    str
        YAML text for the ``metrics`` block.
    """
    entries = {}
    for row in table.itertuples():
        name, unit, description = MLPEG_NAMES[row.property]
        n_points = len(sets[row.property]["error"])
        if row.metric == "MAZE":
            tooltip = (
                f"Mean absolute z-score of the {description} ({n_points} points), "
                "dividing each error by a fixed sqrt(σ_exp² + σ_Sage²), where σ_Sage "
                "is the spread of 5 independent Sage runs"
            )
        else:
            tooltip = f"Mean absolute error of the {description} ({n_points} points)"
        entry = {"good": row.good, "bad": row.bad}
        entry["unit"] = None if row.metric == "MAZE" else unit
        entry["tooltip"] = tooltip
        entry["level_of_theory"] = "Experimental"
        if row.metric == "MAE":
            entry["weight"] = 0
        entries[f"{name} {row.metric}"] = entry
    header = (
        "# Generated by `pixi run -e analysis power-analysis` (thresholds.csv).\n"
        f"# good: perfect-model score + {GOOD_NOISE_SDS:g} SD with one "
        f"{REFERENCE[0] / 1000:g} + {REFERENCE[1] / 1000:g} ns run per box; "
        f"bad: {BAD_SAGE_MULTIPLE:g} x Sage (openff-2.3.0).\n"
    )
    return header + yaml.safe_dump(
        {"metrics": entries}, allow_unicode=True, sort_keys=False
    )


def box_atoms(data_dir: Path, box_ids: list[str]) -> np.ndarray:
    """
    Atom count of every box.

    Parameters
    ----------
    data_dir
        Benchmark data directory (``systems.csv`` and ``pure.csv``).
    box_ids
        Box order.

    Returns
    -------
    np.ndarray
        Atoms per box.
    """
    systems = pd.read_csv(data_dir / "systems.csv").set_index("box_id").n_atoms
    boxes = pd.read_csv(data_dir / "boxes.csv").set_index("box_id")
    pure = (
        pd.read_csv(data_dir / "pure.csv").set_index(["smiles", "temperature"]).n_atoms
    )
    return np.array(
        [
            systems[b]
            if b in systems.index
            else pure.loc[(boxes.smiles[b], boxes.temperature[b])]
            for b in box_ids
        ],
        dtype=float,
    )


# Published MLP condensed-phase results (pure liquids; no MLP ΔHmix or V^E
# study of molecular liquids was found). Each row cites its source.
LITERATURE = [
    {
        "model": "MACE-OFF23(M)",
        "property": "pure density, 109 organic liquids",
        "result": "MAE 0.09 g/mL",
        "source": (
            "[Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2"
        ),
    },
    {
        "model": "MACE-OFF23(S)",
        "property": "pure density, 109 organic liquids",
        "result": "MAE 0.23 g/mL (systematically too dense)",
        "source": (
            "[Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2"
        ),
    },
    {
        "model": "MACE-OFF23(L)",
        "property": "pure density, 109 organic liquids",
        "result": "no significant improvement over (M)",
        "source": (
            "[Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2"
        ),
    },
    {
        "model": "ANI-2x",
        "property": "pure density, 109 organic liquids",
        "result": "MAE and RMSE about twice MACE-OFF23(M)",
        "source": (
            "[Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2"
        ),
    },
    {
        "model": "MACE-OFF23(M)",
        "property": "ΔHvap, 109 organic liquids",
        "result": "systematic offset of about 2 kcal/mol (8 kJ/mol)",
        "source": (
            "[Kovács et al. 2025](https://doi.org/10.1021/jacs.4c07099), sec. III.3.2"
        ),
    },
    {
        "model": "Organic_MPNICE",
        "property": "pure density, 62 solvents",
        "result": "about 4 % mean error (4.9 % in Table 12)",
        "source": "[Weber et al. 2025](https://arxiv.org/abs/2505.06462), sec. 3.1.4",
    },
    {
        "model": "QRNN",
        "property": "pure density, 62 solvents",
        "result": "about 16 % mean error, systematically high",
        "source": "[Weber et al. 2025](https://arxiv.org/abs/2505.06462), sec. 3.1.4",
    },
    {
        "model": "Hybrid / inorganic MPNICE",
        "property": "pure density, 62 solvents",
        "result": "18 to 35 % mean absolute error",
        "source": "[Weber et al. 2025](https://arxiv.org/abs/2505.06462), Table 12",
    },
    {
        "model": "ANI-2x, MACE-OFF23(S, M) vs CGenFF",
        "property": (
            "density and ΔHvap of 6 liquids (water, methanol, acetone, benzene, "
            "hexane, NMA)"
        ),
        "result": (
            "the MLPs miss in both directions (e.g. MACE-OFF23(M) too dense for water, "
            "too light for acetone and hexane); CGenFF best overall"
        ),
        "source": "[Picha et al. 2025](https://doi.org/10.1021/acs.jcim.5c01355)",
    },
    {
        "model": "BAMBOO",
        "property": "density of electrolyte solvent mixtures",
        "result": (
            "about 0.01 g/mL average error, after a density alignment step fitted to "
            "experiment"
        ),
        "source": (
            "[Gong et al. 2025](https://doi.org/10.1038/s42256-025-01009-7) "
            "([arXiv:2404.07181](https://arxiv.org/abs/2404.07181))"
        ),
    },
    {
        "model": "UMA-OMol vs SevenNet-OMat, Orb-OMat",
        "property": "density of ether and carbonate electrolytes",
        "result": (
            "UMA-OMol much closer to experiment than models trained on inorganic data "
            "(qualitative in the abstract)"
        ),
        "source": "[arXiv:2603.20183](https://arxiv.org/abs/2603.20183)",
    },
]


def plot_noise_scaling(noise: pd.DataFrame, path: Path) -> None:
    """
    Plot single-run noise × sqrt(t) against production length.

    Flat lines mean σ ∝ 1/sqrt(t); a rising line means slow modes that only
    longer runs pick up.

    Parameters
    ----------
    noise
        Output of :func:`box_noise`.
    path
        Output PNG.
    """
    fig, ax = plt.subplots(figsize=(5, 3.4))
    ref = noise[noise.equilibration_ps == REFERENCE[0]]
    for color, obs, label in zip(
        COLORS, ("h", "density"), ("h", "density"), strict=False
    ):
        rms = (
            ref[ref.observable == obs]
            .groupby("production_ps")
            .sigma.apply(lambda s: float(np.sqrt(np.mean(s**2))))
        )
        scaled = rms * np.sqrt(rms.index / REFERENCE[1]) / rms[REFERENCE[1]]
        ax.plot(rms.index / 1000, scaled, "o-", color=color, lw=2, ms=6, label=label)
    ax.axhline(1, color="0.5", lw=0.8, ls="--")
    ax.set_ylim(0, 1.2)
    ax.set_xlabel("Production length (ns)")
    reference_ns = REFERENCE[1] / 1000
    ax.set_ylabel(f"σ·sqrt(t / {reference_ns:g} ns) / σ({reference_ns:g} ns)")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_resolvable(
    metrics: pd.DataFrame,
    threshold_table: pd.DataFrame,
    path: Path,
) -> None:
    """
    Plot the resolvable metric difference against production length.

    One panel per metric. The line is the Sage-like error profile, the band
    spans all profiles; dotted is the observed metric of a perfect model.

    Parameters
    ----------
    metrics
        ``metrics`` table from :func:`analyse`.
    threshold_table
        Output of :func:`thresholds`.
    path
        Output PNG.
    """
    good_values = threshold_table.set_index(["property", "metric"]).good
    panels = [
        ("hmix", "MAE", "kJ/mol"),
        ("density", "MAE", "g/mL"),
        ("hmix", "MAZE", ""),
        ("vex", "MAE", "cm³/mol"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.4))
    for ax, (prop, label, unit) in zip(axes.flat, panels, strict=True):
        rows = metrics[(metrics.property == prop) & (metrics.metric == label)]
        t = np.array(PRODUCTION_LENGTHS_PS) / 1000
        band = rows.groupby("production_ps").resolvable.agg(["min", "max"])
        ax.fill_between(t, band["min"], band["max"], color=COLORS[0], alpha=0.2, lw=0)
        sage = rows[rows.profile == 1.0].set_index("production_ps").resolvable
        ax.plot(
            t,
            sage.loc[list(PRODUCTION_LENGTHS_PS)],
            "o-",
            color=COLORS[0],
            lw=2,
            ms=6,
            label="resolvable difference",
        )
        floor = rows[rows.profile == 0.0].set_index("production_ps").observed_mean
        ax.plot(
            t,
            floor.loc[list(PRODUCTION_LENGTHS_PS)],
            "s:",
            color=COLORS[2],
            lw=2,
            ms=6,
            label="perfect-model score",
        )
        good = good_values[(prop, label)]
        ax.axhline(good, color="0.4", lw=0.8)
        ax.annotate(
            "ml-peg 'good'",
            (t[-1], good),
            textcoords="offset points",
            xytext=(0, 3),
            ha="right",
            fontsize=8,
            color="0.3",
        )
        ax.set_title(f"{LABELS[prop]} {label}", fontsize=10)
        ax.set_xlabel("Production length (ns)")
        ax.set_ylabel(unit or "dimensionless")
        ax.set_ylim(bottom=0)
        ax.spines[["top", "right"]].set_visible(False)
    axes.flat[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_point_mde(points: pd.DataFrame, path: Path) -> None:
    """
    Plot each ΔHmix point's detection limit against Sage's error.

    Points above the diagonal are errors one run per box would detect.

    Parameters
    ----------
    points
        ``points`` table from :func:`analyse`.
    path
        Output PNG.
    """
    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    hmix = points[points.property == "hmix"]
    for color, prod in zip(COLORS, (REFERENCE[1], 500.0), strict=False):
        rows = hmix[hmix.production_ps == prod]
        ax.scatter(
            rows.mde,
            rows.sage_error.abs(),
            s=36,
            color=color,
            edgecolors="white",
            linewidths=1,
            label=f"{prod / 1000:g} ns",
        )
    top = float(max(hmix.mde.max(), hmix.sage_error.abs().max())) * 1.05
    ax.plot([0, top], [0, top], color="0.5", lw=0.8, ls="--")
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_xlabel("Minimum detectable error (kJ/mol)")
    ax.set_ylabel("|Sage ΔHmix error| (kJ/mol)")
    ax.legend(frameon=False, title="Production")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    """Run the power analysis and write the tables, figures and report."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data-dir", type=Path, default=Path("data/binary_mixtures"))
    parser.add_argument("--runs-dir", type=Path, default=Path("data/sage_runs"))
    parser.add_argument("--summary-dir", type=Path, default=Path("data/sage_summary"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/power_analysis"))
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    runs = pd.read_csv(args.summary_dir / "runs.csv")
    cache = out / "truncation.csv"
    if cache.exists():
        truncated = pd.read_csv(cache)
    else:
        truncated = load_truncated(runs, args.runs_dir, args.workers)
        truncated.to_csv(cache, index=False)
    noise = box_noise(truncated)
    noise.to_csv(out / "box_noise.csv", index=False)
    equilibration = equilibration_shift(truncated)
    equilibration.to_csv(out / "equilibration.csv", index=False)

    sets = metric_sets(args.data_dir, runs)
    n_atoms = box_atoms(args.data_dir, sets["box_ids"])
    results = analyse(sets, noise, n_atoms)
    for name, table in results.items():
        table.to_csv(out / f"{name}.csv", index=False)
    pd.DataFrame(LITERATURE).to_csv(out / "literature.csv", index=False)
    threshold_table = thresholds(results["metrics"], results["bootstrap"])
    threshold_table.to_csv(out / "thresholds.csv", index=False)
    (out / "metrics_snippet.yml").write_text(metrics_snippet(threshold_table, sets))

    figures = out / "figures"
    figures.mkdir(exist_ok=True)
    plot_noise_scaling(noise, figures / "noise_scaling.png")
    plot_resolvable(
        results["metrics"],
        threshold_table,
        figures / "resolvable.png",
    )
    plot_point_mde(results["points"], figures / "point_mde.png")

    from liquid_mixtures_inputs.power_report import write_reports

    write_reports(
        out,
        {
            "noise": noise,
            "equilibration": equilibration,
            "shifts": early_shifts(truncated),
            **results,
            "thresholds": threshold_table,
            "n_atoms": n_atoms,
            "sage_points": pd.read_csv(args.summary_dir / "points.csv"),
            "finite": pd.read_csv(args.summary_dir / "finite_size.csv"),
        },
    )


if __name__ == "__main__":
    main()
