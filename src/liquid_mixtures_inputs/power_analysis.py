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
- the noise in each benchmark metric (MAE, MAZE, RMSE) and the smallest
  difference between two models it can resolve;
- the spread of the metric over the choice of systems (cluster bootstrap), for
  comparison;
- a GPU-time allocation across boxes that keeps the metric noise of uniform
  2 ns runs at a lower cost.
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

from ml_peg.analysis.molecular_dynamics.binary_mixtures.mixing import (  # noqa: E402
    mixing_matrix,
    observables,
)
from ml_peg.analysis.utils.timeseries import estimate_mean  # noqa: E402

TOTAL_PS = 2200.0
EQUILIBRATION_PS = (100.0, 200.0)
PRODUCTION_PS = (250.0, 500.0, 1000.0, 1500.0, 2000.0)
REFERENCE = (200.0, 2000.0)  # the benchmark protocol
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
# MACE-OFF23(M) on one A100 via OpenMM: 1.1e6 steps/day for ~600-atom liquid boxes
# (arXiv:2312.15211, section III.3.2). Used only for an order-of-magnitude cost.
MACE_OFF_M_ATOM_STEPS_PER_DAY = 1.1e6 * 600
# Benchmark metric thresholds in ml-peg (binary_mixtures/metrics.yml)
GOOD = {
    ("hmix", "MAE"): 0.2,
    ("hmix", "MAZE"): 1.0,
    ("density", "MAE"): 0.005,
    ("density", "MAZE"): 3.0,
    ("vex", "RMSE"): 0.05,
    ("pure_density", "MAE"): 0.005,
}
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
    if data.time_ps.iloc[-1] < TOTAL_PS - 0.5:
        return []
    rows = []
    windows = [(eq, prod) for eq in EQUILIBRATION_PS for prod in PRODUCTION_PS]
    for eq, prod in [*windows, EARLY_WINDOW]:
        window = data[(data.time_ps > eq) & (data.time_ps <= eq + prod + 1e-6)]
        for name, y in observables(window, n_molecules).items():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                estimate = estimate_mean(y)
            rows.append(
                {
                    "box_id": box_id,
                    "replicate": replicate,
                    "equilibration_ps": eq,
                    "production_ps": prod,
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
    chi = (noise["n"] - 1) * (noise["sd"] / noise["sem"]) ** 2
    pooled = (
        noise.assign(chi=chi, dof=noise["n"] - 1)
        .groupby(keys)
        .apply(
            lambda g: float(np.sqrt(g.chi.sum() / g.dof.sum())), include_groups=False
        )
        .rename("k")
    )
    noise = noise.join(pooled, on=keys)
    noise["sigma"] = noise["k"] * noise["sem"]
    return noise


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
    index = ["box_id", "replicate", "observable"]
    early = truncated[
        (truncated.equilibration_ps == EARLY_WINDOW[0])
        & (truncated.production_ps == EARLY_WINDOW[1])
    ].set_index(index)
    full = truncated[
        (truncated.equilibration_ps == REFERENCE[0])
        & (truncated.production_ps == REFERENCE[1])
    ].set_index(index)
    shift = (early["mean"] - full["mean"]).dropna()
    rows = []
    for obs, values in shift.groupby(level="observable"):
        rows.append(
            {
                "observable": obs,
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
        ``sigma_exp``, ``sigma_ref`` (MAZE denominator, or None), ``system``
        (cluster label per point), ``metric`` ("mae" or "rmse") and ``unit``.
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
            "metric": "mae",
            "unit": "kJ/mol",
        },
        "density": {
            "matrix": density,
            "observable": "density",
            "error": (with_density.density_sage - with_density.density_exp).to_numpy(),
            "sigma_exp": with_density.density_sigma.to_numpy(),
            "sigma_ref": with_density.sigma_ref_density.to_numpy(),
            "system": with_density.system_id.to_numpy(),
            "metric": "mae",
            "unit": "g/mL",
        },
        "vex": {
            "matrix": mixing[systems.density_exp.notna().to_numpy()],
            "observable": "v_molar",
            "error": (with_density.vex_sage - with_density.vex_exp).to_numpy(),
            "sigma_exp": np.zeros(len(with_density)),
            "sigma_ref": None,
            "system": with_density.system_id.to_numpy(),
            "metric": "rmse",
            "unit": "cm³/mol",
        },
        "pure_density": {
            "matrix": pure_matrix,
            "observable": "density",
            "error": pure_error,
            "sigma_exp": scored.sigma.to_numpy(),
            "sigma_ref": None,
            "system": scored.smiles.to_numpy(),
            "metric": "mae",
            "unit": "g/mL",
        },
    }


def metric(
    errors: np.ndarray, kind: str, scale: np.ndarray | None = None
) -> np.ndarray:
    """
    Benchmark metric over the last axis of an error array.

    Parameters
    ----------
    errors
        Model minus experiment, shape (..., n_points).
    kind
        "mae" or "rmse".
    scale
        Per-point denominator (MAZE), or None.

    Returns
    -------
    np.ndarray
        Metric with the last axis reduced.
    """
    if scale is not None:
        errors = errors / scale
    if kind == "rmse":
        return np.sqrt(np.mean(errors**2, axis=-1))
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
    prop: dict,
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
    prop
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
        SD (MAE only: sqrt(sᵀCs)/N with s = sign(error), valid when errors
        are well above the noise) and the resolvable difference between two
        independent models, Z_DETECT·√2·SD.
    """
    matrix, error = prop["matrix"], profile * prop["error"]
    noise = (rng.standard_normal((N_DRAWS, len(sigma_box))) * sigma_box) @ matrix.T
    observed = metric(error + noise, prop["metric"], scale)
    true = float(metric(error, prop["metric"], scale))
    sd = float(observed.std(ddof=1))
    closed = np.nan
    if prop["metric"] == "mae" and profile > 0:
        weights = np.sign(error) / (scale if scale is not None else 1.0)
        cov = matrix @ np.diag(sigma_box**2) @ matrix.T
        closed = float(np.sqrt(weights @ cov @ weights) / len(error))
    return {
        "true": true,
        "observed_mean": float(observed.mean()),
        "bias": float(observed.mean()) - true,
        "sd": sd,
        "sd_closed_form": closed,
        "resolvable": float(Z_DETECT * np.sqrt(2) * sd),
    }


def system_bootstrap(
    prop: dict, rng: np.random.Generator, scale: np.ndarray | None = None
) -> float:
    """
    SD of Sage's noiseless metric over resampled systems (cluster bootstrap).

    Parameters
    ----------
    prop
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
    labels = prop["system"]
    clusters = [np.flatnonzero(labels == s) for s in np.unique(labels)]
    error = prop["error"] if scale is None else prop["error"] / scale
    values = np.empty(N_BOOTSTRAP)
    for i in range(N_BOOTSTRAP):
        pick = rng.integers(len(clusters), size=len(clusters))
        values[i] = metric(
            np.concatenate([error[clusters[j]] for j in pick]), prop["metric"]
        )
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
        lam = np.sum(coefficient[free] / shape) / remaining
        t = np.full(len(coefficient), min_ps)
        t[free] = lam * shape
        newly = free & (t < min_ps)
        if not newly.any():
            return t
        fixed |= newly


def sigma_at(
    noise: pd.DataFrame, box_ids: list[str], obs: str, prod: float
) -> np.ndarray:
    """
    Single-run noise of every box for one observable and production length.

    Parameters
    ----------
    noise
        Output of :func:`box_noise`.
    box_ids
        Box order.
    obs
        Observable.
    prod
        Production length in ps (with the reference equilibration).

    Returns
    -------
    np.ndarray
        Single-run noise σ per box, in ``box_ids`` order.
    """
    rows = noise[
        (noise.equilibration_ps == REFERENCE[0])
        & (noise.production_ps == prod)
        & (noise.observable == obs)
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
    props = {k: v for k, v in sets.items() if k != "box_ids"}
    point_rows, metric_rows, boot_rows = [], [], []
    for name, prop in props.items():
        variants = [(prop["metric"].upper(), None)]
        if prop["sigma_ref"] is not None:
            variants.append(("MAZE", prop["sigma_ref"]))
        for label, scale in variants:
            boot_rows.append(
                {
                    "property": name,
                    "metric": label,
                    "sage": float(metric(prop["error"], prop["metric"], scale)),
                    "bootstrap_sd": system_bootstrap(prop, rng, scale),
                }
            )
        for prod in PRODUCTION_PS:
            sigma_box = sigma_at(noise, box_ids, prop["observable"], prod)
            sigma_sim = point_sigma(prop["matrix"], sigma_box)
            sigma_tot = np.sqrt(sigma_sim**2 + prop["sigma_exp"] ** 2)
            for i in range(len(sigma_sim)):
                point_rows.append(
                    {
                        "property": name,
                        "point": i,
                        "system_id": prop["system"][i],
                        "production_ps": prod,
                        "sigma_sim": sigma_sim[i],
                        "sigma_exp": prop["sigma_exp"][i],
                        "mde": Z_DETECT * sigma_tot[i],
                        "sage_error": prop["error"][i],
                    }
                )
            for label, scale in variants:
                for profile in PROFILES:
                    row = metric_noise(prop, sigma_box, profile, rng, scale)
                    metric_rows.append(
                        {
                            "property": name,
                            "metric": label,
                            "production_ps": prod,
                            "profile": profile,
                            **row,
                        }
                    )

    # Allocation for ΔHmix MAE, mixture density MAE and pure density MAE together,
    # each variance normalised by its value for uniform 2 ns runs. With random
    # error signs, E[(Σ_i s_i M_ib)²] = Σ_i M_ib², so A_b = Σ_i M_ib² c_b / N².
    total = np.zeros(len(box_ids))
    for name in ("hmix", "density", "pure_density"):
        prop = props[name]
        c = (
            sigma_at(noise, box_ids, prop["observable"], REFERENCE[1]) ** 2
            * REFERENCE[1]
        )
        a = (prop["matrix"] ** 2).sum(axis=0) * c / len(prop["error"]) ** 2
        total += a / np.sum(a / REFERENCE[1])
    used = total > 0
    uniform_cost = np.sum(n_atoms * (REFERENCE[0] + REFERENCE[1]))
    t = np.full(len(box_ids), np.nan)
    t[used] = allocate(
        total[used], n_atoms[used], MIN_ALLOCATION_PS, np.sum(total) / REFERENCE[1]
    )
    allocation = pd.DataFrame(
        {
            "box_id": box_ids,
            "n_atoms": n_atoms,
            "weight": total,
            "production_ps": t,
        }
    )
    allocation.attrs["uniform_cost"] = uniform_cost
    allocation.attrs["allocated_cost"] = float(
        np.nansum(n_atoms[used] * (REFERENCE[0] + t[used]))
    )
    return {
        "points": pd.DataFrame(point_rows),
        "metrics": pd.DataFrame(metric_rows),
        "bootstrap": pd.DataFrame(boot_rows),
        "allocation": allocation,
    }


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
        "model": "MACE-OFF23",
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
    ax.set_ylabel("σ·sqrt(t / 2 ns) / σ(2 ns)")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_resolvable(metrics: pd.DataFrame, bootstrap: pd.DataFrame, path: Path) -> None:
    """
    Plot the resolvable metric difference against production length.

    One panel per metric. The line is the Sage-like error profile, the band
    spans all profiles; dashed is the system-bootstrap SD of Sage's metric and
    dotted the observed metric of a perfect model.

    Parameters
    ----------
    metrics
        ``metrics`` table from :func:`analyse`.
    bootstrap
        ``bootstrap`` table from :func:`analyse`.
    path
        Output PNG.
    """
    panels = [
        ("hmix", "MAE", "kJ/mol"),
        ("density", "MAE", "g/mL"),
        ("hmix", "MAZE", ""),
        ("vex", "RMSE", "cm³/mol"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.4))
    for ax, (prop, label, unit) in zip(axes.flat, panels, strict=True):
        rows = metrics[(metrics.property == prop) & (metrics.metric == label)]
        t = np.array(PRODUCTION_PS) / 1000
        band = rows.groupby("production_ps").resolvable.agg(["min", "max"])
        ax.fill_between(t, band["min"], band["max"], color=COLORS[0], alpha=0.2, lw=0)
        sage = rows[rows.profile == 1.0].set_index("production_ps").resolvable
        ax.plot(
            t,
            sage.loc[list(PRODUCTION_PS)],
            "o-",
            color=COLORS[0],
            lw=2,
            ms=6,
            label="resolvable difference",
        )
        boot = bootstrap[(bootstrap.property == prop) & (bootstrap.metric == label)]
        ax.axhline(
            float(boot.bootstrap_sd.iloc[0]),
            color=COLORS[1],
            lw=2,
            ls="--",
            label="SD over choice of systems",
        )
        floor = rows[rows.profile == 0.0].set_index("production_ps").observed_mean
        ax.plot(
            t,
            floor.loc[list(PRODUCTION_PS)],
            "s:",
            color=COLORS[2],
            lw=2,
            ms=6,
            label="perfect-model score",
        )
        good = GOOD.get((prop, label))
        if good is not None:
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
    for color, prod in zip(COLORS, (2000.0, 500.0), strict=False):
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


def md_table(frame: pd.DataFrame) -> str:
    """
    Format a data frame as a Markdown table.

    Parameters
    ----------
    frame
        Table with already formatted cells.

    Returns
    -------
    str
        Markdown table.
    """
    lines = [
        "| " + " | ".join(frame.columns) + " |",
        "|" + "---|" * len(frame.columns),
    ]
    lines += ["| " + " | ".join(str(v) for v in row) + " |" for row in frame.to_numpy()]
    return "\n".join(lines)


def gpu_days(atom_ps: float) -> float:
    """
    Convert atom·ps at a 1 fs time step into MACE-OFF23(M) A100 GPU-days.

    Parameters
    ----------
    atom_ps
        Sum over boxes of atoms × simulated ps.

    Returns
    -------
    float
        GPU-days, assuming cost linear in the number of atoms.
    """
    return atom_ps * 1000 / MACE_OFF_M_ATOM_STEPS_PER_DAY


def write_report(out: Path, noise, equilibration, results, n_atoms) -> None:
    """
    Write ``report.md`` from the computed tables.

    Parameters
    ----------
    out
        Output directory (figures are in ``out / "figures"``).
    noise
        Output of :func:`box_noise`.
    equilibration
        Output of :func:`equilibration_shift`.
    results
        Output of :func:`analyse`.
    n_atoms
        Atoms per box.
    """
    metrics, points = results["metrics"], results["points"]
    boot = results["bootstrap"].set_index(["property", "metric"])
    allocation = results["allocation"]

    def m(prop, label, prod, profile, col):
        """
        Look up one value of the metrics table.

        Parameters
        ----------
        prop
            Property.
        label
            Metric name.
        prod
            Production length in ps.
        profile
            Error profile.
        col
            Column.

        Returns
        -------
        float
            The value.
        """
        row = metrics[
            (metrics.property == prop)
            & (metrics.metric == label)
            & (metrics.production_ps == prod)
            & (metrics.profile == profile)
        ]
        return float(row[col].iloc[0])

    def res_range(prop, label, prod):
        """
        Range of the resolvable difference over error profiles.

        Parameters
        ----------
        prop
            Property.
        label
            Metric name.
        prod
            Production length in ps.

        Returns
        -------
        tuple[float, float]
            Smallest and largest resolvable difference.
        """
        rows = metrics[
            (metrics.property == prop)
            & (metrics.metric == label)
            & (metrics.production_ps == prod)
        ].resolvable
        return float(rows.min()), float(rows.max())

    def pt(prop, prod):
        """
        Select the per-point rows of one property and run length.

        Parameters
        ----------
        prop
            Property.
        prod
            Production length in ps.

        Returns
        -------
        pd.DataFrame
            Matching rows of the points table.
        """
        return points[(points.property == prop) & (points.production_ps == prod)]

    ref = REFERENCE[1]
    h_sd = m("hmix", "MAE", ref, 1.0, "sd")
    h_boot = float(boot.loc[("hmix", "MAE"), "bootstrap_sd"])
    d_sd = m("density", "MAE", ref, 1.0, "sd")
    d_boot = float(boot.loc[("density", "MAE"), "bootstrap_sd"])
    total_atoms = float(n_atoms.sum())
    cost = {p: total_atoms * (REFERENCE[0] + p) for p in PRODUCTION_PS}
    alloc_cost = float(
        np.nansum(allocation.n_atoms * (REFERENCE[0] + allocation.production_ps))
    )
    eq = equilibration.set_index("observable")
    k = (
        noise[(noise.equilibration_ps == REFERENCE[0])]
        .groupby(["production_ps", "observable"])
        .k.first()
    )

    lines = [
        "# Can one run per box rank MLPs? (power analysis)",
        "",
        "ml-peg runs every MLP once per box: 0.2 ns of NPT equilibration and 2 ns of "
        "production at a 1 fs time step. This report asks whether that is enough to "
        "tell MLPs apart from each other and from experiment, given the simulation "
        "and experimental error bars, and whether shorter runs would do.",
        "",
        "The 320 Sage runs (64 boxes × 5 replicates, `data/sage_summary/`) stand in "
        "for MLP simulation noise. This assumes an MLP's fluctuations and correlation "
        "times are similar to Sage's, which no MLP run has checked yet.",
        "",
        "## Summary",
        "",
        "- **Simulation noise is small next to the choice of systems.** At 2 ns, the "
        f"ΔHmix MAE of one model moves by {h_sd:.3f} kJ/mol (SD) between repeat "
        "benchmark runs. Which 13 systems were chosen moves Sage's MAE by "
        f"{h_boot:.2f} kJ/mol (cluster bootstrap), {h_boot / h_sd:.0f}× more. For "
        f"mixture density the factor is {d_boot / d_sd:.0f}× ({d_sd:.5f} vs "
        f"{d_boot:.4f} g/mL). How much the choice of systems matters for a *ranking* "
        "depends on how alike two models' errors are (section 3).",
        "- **Two models can be told apart** (95 % test, 80 % power) when their ΔHmix "
        "MAE "
        f"differs by {res_range('hmix', 'MAE', ref)[0]:.2f} to "
        f"{res_range('hmix', 'MAE', ref)[1]:.3f} kJ/mol, or their mixture density MAE "
        "by "
        f"{res_range('density', 'MAE', ref)[1]:.4f} g/mL. Published MLPs differ in "
        "pure-liquid density MAE by about 0.05 to 0.25 g/mL (section 6), 100 to 500× "
        "more.",
        "- **Against experiment, per point,** one run per box detects a ΔHmix error "
        "above "
        f"{pt('hmix', ref).mde.median():.2f} kJ/mol (median over points; "
        f"{(pt('hmix', ref).sage_error.abs() > pt('hmix', ref).mde).mean():.0%} of "
        "Sage's "
        "errors are larger) and a density error above "
        f"{pt('density', ref).mde.median():.4f} "
        "g/mL. Experimental errors are 5 to 12× smaller than simulation noise (median "
        "σ_exp vs σ_sim) and barely change these limits.",
        "- **Two metric thresholds sit below the noise floor.** A perfect model "
        "scores "
        f"{m('vex', 'RMSE', ref, 0.0, 'observed_mean'):.3f} cm³/mol on V^E RMSE "
        "(ml-peg "
        f"'good' is {GOOD[('vex', 'RMSE')]}) and "
        f"{m('hmix', 'MAZE', ref, 0.0, 'observed_mean'):.2f} "
        f"on ΔHmix MAZE ('good' is {GOOD[('hmix', 'MAZE')]}). See section 5.",
        f"- **Run length.** Production of 1 ns costs {cost[1000.0] / cost[ref]:.0%} of "
        "the "
        f"current protocol and widens the resolvable ΔHmix MAE difference from "
        f"{m('hmix', 'MAE', ref, 1.0, 'resolvable'):.3f} to "
        f"{m('hmix', 'MAE', 1000.0, 1.0, 'resolvable'):.3f} kJ/mol. Noise falls a "
        "little "
        "more slowly than 1/sqrt(t), because longer runs sample slow modes that short "
        "runs miss. Per-box run lengths optimised for cost save only "
        f"{1 - alloc_cost / cost[ref]:.0%}. Equilibration should stay at 0.2 ns.",
        "",
        "## Recommendation",
        "",
        "- Keep 0.2 + 2 ns as the default. At that length a rerun reproduces any "
        "ranking whose gap exceeds the resolvable difference, and simulation noise is "
        "a small part of each metric's uncertainty except for V^E RMSE.",
        f"- 0.2 + 1 ns is a reasonable cheaper option (about "
        f"{gpu_days(cost[ref]) - gpu_days(cost[1000.0]):.0f} MACE-OFF23(M) A100 "
        "GPU-days "
        f"saved per model, of {gpu_days(cost[ref]):.0f}; see section 4). Decide once "
        "MLP "
        "pilot runs show whether MLP correlation times match Sage's.",
        "- Do not shorten equilibration: the 100 to 200 ps window is still shifted "
        "for Sage "
        "(section 4), and an MLP must also relax from Sage's structure to its own "
        "density.",
        "- Raise the ΔHmix MAZE and V^E RMSE 'good' thresholds above their noise "
        "floors, "
        "or bias-correct V^E RMSE (section 5).",
        "",
        "## Method",
        "",
        "1. **Run length.** Every Sage run is re-analysed with the benchmark's "
        "estimator "
        "(`mixing.observables` and `timeseries.estimate_mean`, as in `summarise-sage`) "
        f"keeping {', '.join(f'{p / 1000:g}' for p in PRODUCTION_PS)} ns of production "
        "after "
        "0.2 ns of equilibration. The 2 ns analysis reproduces `runs.csv` exactly "
        "(tested).",
        "2. **Single-run noise per box.** σ = k·sem, with sem the RMS within-run "
        "standard "
        "error over the 5 replicates and k the pooled ratio of replicate SD to "
        "within-run "
        "sem over all boxes at that length (k ≈ 1.2 for h, 1.0 to 1.2 for density). "
        "This "
        "is smoother than the 5-replicate SD of each box and corrects the within-run "
        "sem "
        "for the slow modes it misses.",
        "3. **Points.** ΔHmix and V^E are linear in the box means (ml-peg's "
        "`mixing_matrix`). Each pure box feeds all three compositions of its system, "
        "so "
        "noise is drawn per box and propagated, which keeps the correlation between "
        "points. The propagated ΔHmix noise at 2 ns (RMS over points "
        f"{float(np.sqrt((pt('hmix', ref).sigma_sim ** 2).mean())):.3f} kJ/mol) "
        "matches "
        "the observed SD of Sage's replicate ΔHmix "
        "(RMS 0.095 kJ/mol, `data/sage_summary/points.csv`).",
        f"4. **Metrics.** For each run length, {N_DRAWS} draws of one run per box give "
        "the "
        "noise of every benchmark metric for a model whose true errors are 1, 0.5, "
        "0.25 "
        "or 0 × Sage's (replicate-mean minus experiment). Two models run "
        "independently, "
        "so their metrics differ by noise with SD sqrt(2)·SD(metric). The **resolvable "
        f"difference** is {Z_DETECT:.2f}·sqrt(2)·SD(metric): a true difference that "
        "size "
        "is detected at the two-sided 5 % level 80 % of the time. For MAE, a closed "
        "form "
        "sqrt(sᵀCs)/N (s the signs of the errors, C the point covariance) is given as "
        "a "
        "check; it applies when errors are well above the noise.",
        f"5. **System choice.** A cluster bootstrap over systems ({N_BOOTSTRAP} "
        "resamples of "
        "the 13 mixtures, or of the 11 scored pure liquids) of Sage's noiseless metric "
        "measures how much the metric depends on which systems were picked.",
        "",
        "## 1. Noise against run length",
        "",
        "![noise scaling](figures/noise_scaling.png)",
        "",
        "Single-run noise × sqrt(t / 2 ns), RMS over boxes and relative to 2 ns. A "
        "flat "
        "line would mean σ ∝ 1/sqrt(t). For h the scaled noise grows by about 15 % "
        "from 0.5 to 2 ns: short runs miss some slow fluctuations, so cutting the run "
        "length costs a little less precision than 1/sqrt(t) predicts. The density "
        "line wanders by a similar amount without a trend, which is within the "
        "uncertainty of k from 5 replicates.",
        "",
    ]
    rows = []
    for p in PRODUCTION_PS:
        rows.append(
            [
                f"{p / 1000:g}",
                f"{pt('hmix', p).sigma_sim.median():.3f}",
                f"{pt('density', p).sigma_sim.median():.5f}",
                f"{pt('vex', p).sigma_sim.median():.3f}",
                f"{pt('pure_density', p).sigma_sim.median():.5f}",
                f"{k.loc[(p, 'h')]:.2f}",
                f"{k.loc[(p, 'density')]:.2f}",
            ]
        )
    lines += [
        md_table(
            pd.DataFrame(
                rows,
                columns=[
                    "Production (ns)",
                    "σ ΔHmix (kJ/mol)",
                    "σ density (g/mL)",
                    "σ V^E (cm³/mol)",
                    "σ pure density (g/mL)",
                    "k (h)",
                    "k (density)",
                ],
            )
        ),
        "",
        "Median single-run noise per point. k is the replicate SD / within-run sem "
        "ratio "
        "(1 if the within-run sem were right).",
        "",
        "## 2. A model against experiment, per point",
        "",
        "Minimum detectable error MDE = "
        f"{Z_DETECT:.2f}·sqrt(σ_sim² + σ_exp²): an error that size is significant at "
        "the 5 % level 80 % of the time. The fraction column shows how many of Sage's "
        "errors (a realistic model) exceed the MDE.",
        "",
    ]
    rows = []
    for prop in ("hmix", "density", "vex", "pure_density"):
        for p in (500.0, 1000.0, 2000.0):
            q = pt(prop, p)
            rows.append(
                [
                    LABELS[prop],
                    f"{p / 1000:g}",
                    f"{q.mde.median():.4g}",
                    f"{q.mde.max():.4g}",
                    f"{(q.sage_error.abs() > q.mde).mean():.0%}",
                    f"{q.sigma_exp.median():.2g}",
                ]
            )
    lines += [
        md_table(
            pd.DataFrame(
                rows,
                columns=[
                    "Property",
                    "Production (ns)",
                    "Median MDE",
                    "Max MDE",
                    "Sage errors above MDE",
                    "Median σ_exp",
                ],
            )
        ),
        "",
        "Units: kJ/mol (ΔHmix), g/mL (densities), cm³/mol (V^E).",
        "",
        "![per-point detection](figures/point_mde.png)",
        "",
        "ΔHmix points above the diagonal have a Sage error that one run per box "
        "detects.",
        "",
        "## 3. Ranking two models",
        "",
        "![resolvable difference](figures/resolvable.png)",
        "",
        "Blue: resolvable difference for a Sage-like error profile (band: all "
        "profiles, "
        "1 to 0 × Sage's errors). Orange: SD of Sage's metric over the choice of "
        "systems. "
        "Green: the score a perfect model gets from noise alone. Grey: ml-peg's 'good' "
        "threshold.",
        "",
    ]
    rows = []
    for prop, label in (
        ("hmix", "MAE"),
        ("hmix", "MAZE"),
        ("density", "MAE"),
        ("density", "MAZE"),
        ("vex", "RMSE"),
        ("pure_density", "MAE"),
    ):
        lo2, hi2 = res_range(prop, label, ref)
        lo1, hi1 = res_range(prop, label, 1000.0)
        rows.append(
            [
                f"{LABELS[prop]} {label}",
                f"{m(prop, label, ref, 1.0, 'true'):.4g}",
                f"{m(prop, label, ref, 1.0, 'sd'):.2g}",
                f"{m(prop, label, ref, 1.0, 'sd_closed_form'):.2g}",
                f"{lo2:.2g} to {hi2:.2g}",
                f"{lo1:.2g} to {hi1:.2g}",
                f"{float(boot.loc[(prop, label), 'bootstrap_sd']):.2g}",
                f"{m(prop, label, ref, 0.0, 'observed_mean'):.3g}",
                f"{GOOD[(prop, label)]}",
            ]
        )
    lines += [
        md_table(
            pd.DataFrame(
                rows,
                columns=[
                    "Metric",
                    "Sage",
                    "SD (sim, 2 ns)",
                    "SD closed form",
                    "Resolvable (2 ns)",
                    "Resolvable (1 ns)",
                    "SD over systems",
                    "Perfect-model score (2 ns)",
                    "ml-peg 'good'",
                ],
            )
        ),
        "",
        "Noise also biases each score upwards, most for good models: a model whose "
        "true "
        f"ΔHmix MAE is {m('hmix', 'MAE', ref, 0.25, 'true'):.2f} kJ/mol scores "
        f"{m('hmix', 'MAE', ref, 0.25, 'observed_mean'):.2f} on average at 2 ns, "
        "against "
        f"{m('hmix', 'MAE', ref, 1.0, 'observed_mean'):.2f} for a true "
        f"{m('hmix', 'MAE', ref, 1.0, 'true'):.2f}. Differences between good models "
        "are "
        "therefore slightly compressed, but their order is kept.",
        "",
        "The two kinds of uncertainty answer different questions. The simulation noise "
        "says whether rerunning the benchmark would give the same ranking on these 13 "
        "systems; at 2 ns it would, for any gap above the resolvable difference. The "
        "system spread says how far a model's score depends on this particular set of "
        "systems, which no run length can reduce.",
        "",
        "For a ranking, what matters is the system spread of the *difference* "
        "between two models, which depends on how alike their errors are. If model B's "
        "errors are model A's scaled by (1 + s), the gap is s·MAE_A and its system SD "
        f"is s·{h_boot:.2f}, i.e. {h_boot / m('hmix', 'MAE', ref, 1.0, 'true'):.2f}× "
        "the gap at any size, so the choice of systems never threatens that ranking "
        "and simulation noise is the only limit. If their errors are unrelated, the "
        f"gap's system SD approaches sqrt(2)·{h_boot:.2f} = "
        f"{np.sqrt(2) * h_boot:.2f} kJ/mol, and gaps below about twice that are "
        "specific to this set of systems even though a rerun would reproduce them. "
        "Real model pairs (e.g. two sizes of one model family) are likely to sit "
        "between these cases.",
        "",
        "## 4. Shorter runs",
        "",
    ]
    rows = []
    for p in PRODUCTION_PS:
        rows.append(
            [
                f"0.2 + {p / 1000:g}",
                f"{cost[p] / cost[ref]:.0%}",
                f"{gpu_days(cost[p]):.0f}",
                f"{m('hmix', 'MAE', p, 1.0, 'resolvable'):.3f}",
                f"{m('density', 'MAE', p, 1.0, 'resolvable'):.4f}",
                f"{m('hmix', 'MAZE', p, 0.0, 'observed_mean'):.2f}",
                f"{m('vex', 'RMSE', p, 0.0, 'observed_mean'):.3f}",
            ]
        )
    lines += [
        md_table(
            pd.DataFrame(
                rows,
                columns=[
                    "Protocol (ns)",
                    "Cost",
                    "MACE-OFF23(M) GPU-days",
                    "Resolvable ΔHmix MAE",
                    "Resolvable density MAE",
                    "Perfect ΔHmix MAZE",
                    "Perfect V^E RMSE",
                ],
            )
        ),
        "",
        f"Cost is Σ atoms × simulated time over the 64 boxes ({total_atoms:.0f} "
        "atoms). "
        "GPU-days are an order of magnitude only: MACE-OFF23(M) runs 1.1×10⁶ steps/day "
        "for ~600-atom liquid boxes on one A100 "
        "([arXiv:2312.15211](https://arxiv.org/abs/2312.15211), sec. III.3.2), scaled "
        "linearly to our 950 to 2,800-atom boxes. Faster and slower models scale the "
        "column, not the ratios.",
        "",
        "**Equilibration.** The mean over the 100 to 200 ps window minus the mean over "
        "the 2 ns production, over all 320 Sage runs:",
        "",
        md_table(
            pd.DataFrame(
                [
                    [
                        o,
                        f"{eq.loc[o, 'mean_shift']:+.4g}",
                        f"{eq.loc[o, 'rms_shift']:.4g}",
                        f"{eq.loc[o, 'p']:.1g}",
                    ]
                    for o in ("h", "density")
                ],
                columns=["Observable", "Mean shift", "RMS shift", "p (mean = 0)"],
            )
        ),
        "",
        "Units: kJ/mol per molecule (h), g/mL (density). Even for Sage, starting from "
        "Sage-built boxes, the window is shifted on average, mostly in replicates 1 "
        "to 3 "
        "(re-packed boxes; density +0.0013 to +0.0016 g/mL) rather than replicate 0, "
        "the "
        "distributed box that MLP runs start from (-0.0004 g/mL). An MLP also has to "
        "relax from Sage's density "
        "to its own, which this test cannot see, so 0.2 ns stays.",
        "",
        "**Per-box run lengths.** For a fixed variance of ΔHmix MAE, mixture density "
        "MAE "
        "and pure density MAE (equally weighted, each relative to uniform 2 ns), the "
        "cheapest allocation runs each box for t ∝ sqrt(A_b / n_b), where A_b is the "
        "box's variance contribution and n_b its atom count, with at least "
        f"{MIN_ALLOCATION_PS / 1000:g} ns. It ranges from "
        f"{allocation.production_ps.min() / 1000:.1f} to "
        f"{allocation.production_ps.max() / 1000:.1f} ns (`allocation.csv`): pure "
        "water "
        "and the small-weight pure boxes get less, slow pure liquids (N-methylaniline, "
        "DMAE, pyrrole, THP) more. It saves "
        f"{1 - alloc_cost / cost[ref]:.0%} of the cost. That assumes σ ∝ 1/sqrt(t), "
        "which section 1 shows is optimistic for short runs, so the real saving is "
        "smaller. Not worth the extra complexity.",
        "",
        "## 5. Thresholds below the noise floor",
        "",
        "- **V^E RMSE.** A perfect model scores "
        f"{m('vex', 'RMSE', ref, 0.0, 'observed_mean'):.3f} "
        "cm³/mol at 2 ns (median per-point σ "
        f"{pt('vex', ref).sigma_sim.median():.3f}): "
        "V^E is a small difference of molar volumes, and RMSE adds the noise in "
        f"quadrature. The 'good' threshold of {GOOD[('vex', 'RMSE')]} cannot be "
        "reached. "
        "Options: raise 'good' to 0.15; or report a bias-corrected "
        "sqrt(max(0, MSE − mean σ_i²)), with σ_i each point's own propagated standard "
        "error inflated by k ≈ 1.2.",
        "- **ΔHmix MAZE.** A perfect model scores "
        f"{m('hmix', 'MAZE', ref, 0.0, 'observed_mean'):.2f} "
        f"at 2 ns and {m('hmix', 'MAZE', 1000.0, 0.0, 'observed_mean'):.2f} at 1 ns, "
        "because σ_ref "
        "is close to the single-run noise. The 'good' threshold of "
        f"{GOOD[('hmix', 'MAZE')]} is barely reachable. Proposal: 'good' = 1.5.",
        f"- ΔHmix MAE ('good' {GOOD[('hmix', 'MAE')]} kJ/mol, floor "
        f"{m('hmix', 'MAE', ref, 0.0, 'observed_mean'):.3f}), density MAE and MAZE, "
        "and pure "
        "density MAE are well above their floors.",
        "",
        "## 6. Literature: how far apart are current MLPs?",
        "",
        md_table(
            pd.DataFrame(LITERATURE)[["model", "property", "result", "source"]].rename(
                columns=str.capitalize
            )
        ),
        "",
        "No study computing ΔHmix or V^E of molecular liquids with an MLP was found "
        "(searches of the web and of Consensus for MLP / neural network potential "
        "enthalpy of mixing, excess enthalpy and excess volume, October 2026). The "
        "closest work concerns metallic alloys "
        "([arXiv:2406.17499](https://arxiv.org/abs/2406.17499)) or fits mixing data "
        "without MD. For comparison, Sage's ΔHmix MAE on this set is "
        f"{m('hmix', 'MAE', ref, 1.0, 'true'):.2f} kJ/mol "
        "(`data/sage_summary/sage_report.md`).",
        "",
        "**What this means for the benchmark.**",
        "",
        "- **Density.** Published MLPs differ in pure-liquid density MAE by about "
        "0.05 to "
        "0.25 g/mL (MACE-OFF23 S vs M, ANI-2x vs MACE-OFF23(M), MPNICE variants). The "
        f"benchmark resolves {res_range('density', 'MAE', ref)[1]:.4f} g/mL at 2 ns "
        "and "
        f"{res_range('density', 'MAE', 250.0)[1]:.4f} g/mL at 0.25 ns. Density alone "
        "could be ranked from runs a tenth as long, but the runs serve ΔHmix too.",
        "- **ΔHmix.** There is no MLP baseline. MLP ΔHvap errors are several kJ/mol "
        "(MACE-OFF23 offset about 8 kJ/mol; Picha et al. find large ΔHvap deviations "
        "for MACE-OFF23(M)). ΔHmix is a difference of such cohesive energies, so much "
        "of that error may cancel, but whatever remains is compared with |ΔHmix| at "
        "x ≈ 0.5 of 0.4 to 2.6 kJ/mol. If MLPs differ in ΔHmix MAE by even a fifth of "
        "Sage's MAE "
        f"(about {0.2 * m('hmix', 'MAE', ref, 1.0, 'true'):.2f} kJ/mol), 2 ns resolves "
        "it. "
        "This benchmark would be the first systematic MLP ΔHmix data, which is itself "
        "an argument for it.",
        "",
        "## Caveats",
        "",
        "- **Sage stands in for MLP noise.** MLP fluctuations of h and volume should "
        "be "
        "similar (the same physical liquids at the same T and P), but correlation "
        "times "
        "may differ, especially for MLPs that change conformer barriers. One pilot "
        "MLP run "
        "per box, analysed with `mixing.box_result`, would check this: compare its sem "
        "and n_eff with `runs.csv`.",
        "- **Slow modes.** k ≈ 1.2 is pooled over boxes. Single boxes (e.g. NMP + DCE "
        "at x = 0.74 for Sage) have larger run-to-run spreads, and an MLP may have "
        "slow "
        "modes of its own.",
        "- **Systematic simulation errors** (finite size, starting structures) are not "
        "noise and are not included. In the finite-size check in `sage_report.md` most "
        "ΔHmix gaps to the published 1000-molecule runs are below 0.06 kJ/mol; the "
        "largest are acetone + N-methylaniline (0.11) and butanone + dibutylamine (up "
        "to "
        "0.21 kJ/mol).",
        "- **Error profiles.** The model pairs are built from Sage's errors. Real "
        "MLPs may "
        "spread their errors differently; the closed-form SD shows the result depends "
        "mainly on the noise, not on the profile.",
        "",
        "Generated by `pixi run -e analysis power-analysis`. Tables: `box_noise.csv`, "
        "`truncation.csv`, `points.csv`, `metrics.csv`, `bootstrap.csv`, "
        "`allocation.csv`, `equilibration.csv`, `literature.csv`.",
        "",
    ]
    (out / "report.md").write_text("\n".join(lines))


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

    figures = out / "figures"
    figures.mkdir(exist_ok=True)
    plot_noise_scaling(noise, figures / "noise_scaling.png")
    plot_resolvable(
        results["metrics"], results["bootstrap"], figures / "resolvable.png"
    )
    plot_point_mde(results["points"], figures / "point_mde.png")
    write_report(out, noise, equilibration, results, n_atoms)


if __name__ == "__main__":
    main()
