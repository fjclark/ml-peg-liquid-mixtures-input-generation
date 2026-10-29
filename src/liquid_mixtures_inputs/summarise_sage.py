"""
Summarise the Sage reference runs of the binary mixtures benchmark.

Run in the ``analysis`` pixi environment after ``run_sage_benchmark.py``::

    pixi run -e analysis summarise-sage \
        --data-dir D --runs-dir R --output-dir O

It computes ΔHmix, mixture density and V^E for every replicate, the replicate
spread used as the MAZE denominator (σ_Sage, pooled per system), convergence
diagnostics, a comparison of within-run and between-run uncertainties, and a
finite-size check against the published 1000-molecule openff-2.3.0 results.
Per-run statistics use the same estimators as the benchmark analysis
(``ml_peg/analysis/molecular_dynamics/binary_mixtures/mixing.py``).
``systems.csv``, ``pure.csv`` and ``boxes.csv`` in ``--data-dir`` are updated
in place with the Sage reference columns. The reports are written by
``sage_report.py``.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import re
import warnings

import matplotlib as mpl

mpl.use("Agg")
from ase.io import read  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from liquid_mixtures_inputs.constants import (  # noqa: E402
    AVOGADRO,
    EQUILIBRATION_PS,
    KJ_PER_MOL_PER_EV,
    PRODUCTION_PS,
)
from ml_peg.analysis.molecular_dynamics.binary_mixtures.mixing import (  # noqa: E402
    MIN_N_EFF,
    box_result,
    frame_checks,
    holm,
    mixture_box_id,
    read_production,
)

# Family-wise significance level for every test in this summary (Holm-adjusted
# per box and per point, Bonferroni over the pooled tests).
FLAG_P = 0.05
START_T = 4.0  # |t| above which replicate 0 is reported as start-dependent
# Sage 2.3.0rc2, benchmarked under its fit name. Its bonded, vdW and library-charge
# parameters are identical to the released openff-2.3.0.offxml; its NAGL charges
# use openff-gnn-am1bcc-0.1.0-rc.3 where we use 1.0.0.
PUBLISHED_FF = "fb-fit-v3-single-mean-k100"
PUBLISHED_FF_URL = (
    "https://github.com/openforcefield/ash-sage-rc2/blob/"
    "00c9524d38af53725a1a18feb3433f5c191cdb4a/04_benchmark/forcefields/"
    "fb-fit-v3-single-mean-k100.offxml?plain=1"
)
PUBLISHED_X_TOLERANCE = 0.02
# Systems whose σ_Sage is set to the floor: their replicate spread comes from a
# Sage artefact that a model need not share (the slow Cl-C-C-Cl torsion of
# 1,2-dichloroethane; see DECISIONS.md).
SIGMA_SAGE_AT_FLOOR = {"nmp_dichloroethane"}
PURE_SIGMA_SAGE_AT_FLOOR = {"ClCCCl"}  # pure 1,2-dichloroethane, same reason
# Thresholds for the "systems to inspect" summary
SIZE_GAP = {"hmix": 0.15, "density": 0.002}  # kJ/mol, g/mL
EARLY_SHIFT = 0.005  # g/mL
EARLY_WINDOW_PS = 200.0  # start of production tested for an early density shift
DEMIX_EXCESS = 0.05  # like-neighbour fraction above random mixing
FF_HMIX_ERROR = 0.4  # kJ/mol, system-mean |Sage - experiment|


def load_run(task: tuple) -> dict | None:
    """
    Calculate statistics and frame checks of one Sage replicate.

    Parameters
    ----------
    task
        ``(box, box_path, csv_path)``: the ``boxes.csv`` row as a dict, the
        input box and the replicate's scalar CSV.

    Returns
    -------
    dict | None
        One row of per-run statistics, or ``None`` if the run is incomplete.
    """
    box, box_path, csv_path = task
    production = read_production(
        csv_path, EQUILIBRATION_PS, EQUILIBRATION_PS + PRODUCTION_PS
    )
    if production is None:
        return None
    # The Sage runs are longer than the protocol; keep only its window.
    end = EQUILIBRATION_PS + PRODUCTION_PS + 1e-6
    production = production[production.time_ps <= end].reset_index(drop=True)
    reference = read(box_path)
    n_molecules = int(reference.info["n_molecules"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        result = box_result(production, n_molecules)
    row = {
        "box_id": box["box_id"],
        "replicate": int(csv_path.stem[3:]),
        "kind": box["kind"],
        "system_id": box["system_id"],
        "smiles": box["smiles"],
        "temperature": box["temperature"],
        "x": box["x"],
        "n_molecules": n_molecules,
        "run_ps": result.run_ps,
        "mean_T": float(production.temperature_K.mean()),
    }
    for name, estimate in result.estimates.items():
        row[f"{name}_mean"] = estimate.mean
        row[f"{name}_sigma"] = estimate.sem
        row[f"{name}_n_eff"] = estimate.n_eff
        row[f"{name}_tau_ps"] = PRODUCTION_PS / estimate.n_eff
        row[f"{name}_p_drift"] = result.p_drift[name]
        row[f"{name}_p_half"] = result.p_half[name]
        row[f"{name}_block_ratio"] = result.block_ratio[name]
    # First 200 ps of production against the rest: a sign that the density was
    # still relaxing after equilibration.
    time = production.time_ps
    first = production.density_g_cm3[time < EQUILIBRATION_PS + EARLY_WINDOW_PS]
    rest = production.density_g_cm3[time >= EQUILIBRATION_PS + EARLY_WINDOW_PS]
    row["density_first_block_shift"] = float(first.mean() - rest.mean())

    frames_path = csv_path.with_name(f"{csv_path.stem}_frames.xyz")
    row["reacted"], row["like_fraction"] = np.nan, np.nan
    if frames_path.exists():
        frames = [
            frame
            for frame in read(frames_path, ":")
            if EQUILIBRATION_PS - 1e-6
            <= frame.info["time_ps"]
            <= EQUILIBRATION_PS + PRODUCTION_PS + 1e-6
        ]
        reacted, fraction = frame_checks(frames, reference, box["kind"] == "mixture")
        row["reacted"], row["like_fraction"] = reacted, fraction
    return row


def load_runs(data_dir: Path, runs_dir: Path, workers: int) -> pd.DataFrame:
    """
    Read every finished replicate and compute per-run statistics.

    Parameters
    ----------
    data_dir
        Benchmark data directory (``boxes.csv`` and ``boxes/``).
    runs_dir
        Output directory of ``run_sage_benchmark.py``.
    workers
        Number of processes.

    Returns
    -------
    pd.DataFrame
        One row per box and replicate.
    """
    boxes = pd.read_csv(data_dir / "boxes.csv")
    tasks = [
        (box, data_dir / "boxes" / f"{box['box_id']}.xyz", path)
        for box in boxes.to_dict("records")
        for path in sorted((runs_dir / box["box_id"]).glob("rep*.csv"))
    ]
    with ProcessPoolExecutor(workers) as pool:
        rows = [row for row in pool.map(load_run, tasks) if row is not None]
    return pd.DataFrame(rows)


def mixture_points(runs: pd.DataFrame, systems: pd.DataFrame) -> pd.DataFrame:
    """
    Compute ΔHmix, density and V^E for every mixture box and replicate.

    Each replicate's ΔHmix uses the same replicate of the pure boxes.

    Parameters
    ----------
    runs
        Output of :func:`load_runs`.
    systems
        ``systems.csv`` with experimental values.

    Returns
    -------
    pd.DataFrame
        One row per mixture box and replicate.
    """
    pure = runs[runs.kind == "pure"].set_index(["smiles", "temperature", "replicate"])
    rows = []
    for experiment in systems.itertuples():
        box_id = mixture_box_id(experiment.system_id, experiment.x_box)
        for run in runs[runs.box_id == box_id].itertuples():
            try:
                a = pure.loc[(experiment.a, run.temperature, run.replicate)]
                b = pure.loc[(experiment.b, run.temperature, run.replicate)]
            except KeyError:
                continue
            x = run.x
            point = {
                "system_id": experiment.system_id,
                "box_id": box_id,
                "replicate": run.replicate,
                "x": x,
                "density": run.density_mean,
                "density_sigma_run": run.density_sigma,
                "hmix_exp": experiment.hmix_exp,
                "hmix_exp_sigma": experiment.hmix_sigma,
                "density_exp": experiment.density_exp,
                "density_exp_sigma": experiment.density_sigma,
                "vex_exp": experiment.vex_exp,
                "vex_exp_sigma": getattr(experiment, "vex_sigma", np.nan),
            }
            # Excess property = mixture - x·pure a - (1 - x)·pure b, per molecule.
            for name, observable in (("hmix", "h"), ("vex", "v_molar")):
                point[name] = (
                    getattr(run, f"{observable}_mean")
                    - x * a[f"{observable}_mean"]
                    - (1 - x) * b[f"{observable}_mean"]
                )
                point[f"{name}_sigma_run"] = np.sqrt(
                    getattr(run, f"{observable}_sigma") ** 2
                    + (x * a[f"{observable}_sigma"]) ** 2
                    + ((1 - x) * b[f"{observable}_sigma"]) ** 2
                )
            rows.append(point)
    return pd.DataFrame(rows)


def pooled_sigma(points: pd.DataFrame, n_replicates: int) -> pd.DataFrame:
    """
    Compute σ_Sage and σ_ref (the MAZE denominator) of every mixture point.

    σ_Sage is the replicate SD, pooled (RMS) over a system's compositions and
    floored at the median over systems, so that no system dominates MAZE
    because Sage happens to sample it quietly. Systems in
    ``SIGMA_SAGE_AT_FLOOR`` get the floor itself. σ_ref = sqrt(σ_Sage² + σ_exp²).

    Parameters
    ----------
    points
        Output of :func:`mixture_points`.
    n_replicates
        Required number of replicates per point.

    Returns
    -------
    pd.DataFrame
        One row per mixture box with means, SDs, pooled σ_Sage and σ_ref.
    """
    grouped = points.groupby(["system_id", "box_id"])
    summary = grouped.agg(
        x=("x", "first"),
        n_rep=("replicate", "count"),
        hmix_sage=("hmix", "mean"),
        hmix_sd=("hmix", "std"),
        hmix_sigma_run=("hmix_sigma_run", "mean"),
        density_sage=("density", "mean"),
        density_sd=("density", "std"),
        vex_sage=("vex", "mean"),
        vex_sd=("vex", "std"),
        hmix_exp=("hmix_exp", "first"),
        hmix_exp_sigma=("hmix_exp_sigma", "first"),
        density_exp=("density_exp", "first"),
        density_exp_sigma=("density_exp_sigma", "first"),
        vex_exp=("vex_exp", "first"),
        vex_exp_sigma=("vex_exp_sigma", "first"),
    ).reset_index()
    summary["complete"] = summary.n_rep >= n_replicates
    for property_name in ("hmix", "density", "vex"):
        # RMS of the replicate SDs over a system's compositions.
        pooled = summary.groupby("system_id")[f"{property_name}_sd"].apply(
            lambda sds: float(np.sqrt(np.nanmean(sds**2)))
        )
        floor = float(np.nanmedian(pooled))
        sigma = summary.system_id.map(pooled).clip(lower=floor)
        sigma[summary.system_id.isin(SIGMA_SAGE_AT_FLOOR)] = floor
        summary[f"sigma_sage_{property_name}"] = sigma.round(6)
        summary[f"sigma_ref_{property_name}"] = np.sqrt(
            summary[f"sigma_sage_{property_name}"] ** 2
            + summary[f"{property_name}_exp_sigma"].fillna(0.0) ** 2
        ).round(6)
    return summary


def pure_sigma(runs: pd.DataFrame, pure: pd.DataFrame) -> pd.DataFrame:
    """
    Compute σ_Sage and σ_ref (the MAZE denominator) of every pure density.

    The same rule as :func:`pooled_sigma`, with one box per liquid and
    temperature: the replicate SD, floored at the median over pure boxes.
    Liquids in ``PURE_SIGMA_SAGE_AT_FLOOR`` get the floor itself.

    Parameters
    ----------
    runs
        Per-run statistics.
    pure
        ``pure.csv`` (experimental density and ``sigma``).

    Returns
    -------
    pd.DataFrame
        One row per pure box with the Sage mean, SD, σ_Sage and σ_ref.
    """
    table = (
        runs[runs.kind == "pure"]
        .groupby(["smiles", "temperature"])
        .density_mean.agg(n_rep="count", density_sage="mean", density_sd="std")
        .reset_index()
        .merge(
            pure[["smiles", "temperature", "density", "sigma", "scored"]],
            on=["smiles", "temperature"],
            how="left",
        )
        .rename(columns={"density": "density_exp", "sigma": "density_exp_sigma"})
    )
    floor = float(np.nanmedian(table.density_sd))
    sigma = table.density_sd.clip(lower=floor)
    sigma[table.smiles.isin(PURE_SIGMA_SAGE_AT_FLOOR)] = floor
    table["sigma_sage_density"] = sigma.round(6)
    table["sigma_ref_density"] = np.sqrt(
        table.sigma_sage_density**2 + table.density_exp_sigma.fillna(0.0) ** 2
    ).round(6)
    return table


def published_x_first(row: pd.Series) -> float:
    """
    Parse the mole fraction of ``smiles_1`` from a published substance string.

    Parameters
    ----------
    row
        Row of the published results, e.g. substance ``CO{solv}{x=0.25}|O{...}``.

    Returns
    -------
    float
        Mole fraction, or NaN if it cannot be parsed.
    """
    pattern = re.escape(str(row.smiles_1)) + r"\{solv\}\{x=([0-9.]+)\}"
    match = re.search(pattern, row.substance)
    return float(match.group(1)) if match else np.nan


def published_sage(cache_dir: Path) -> pd.DataFrame:
    """
    Load the published openff-2.3.0 results (3 replicates, 1000 molecules).

    Parameters
    ----------
    cache_dir
        Cache of ``select_systems.py`` holding the ash-sage-rc2 summary CSVs.

    Returns
    -------
    pd.DataFrame
        Published values with the parsed mole fraction of ``smiles_1``.
    """
    frames = []
    for name in (
        "training-hmix",
        "validation-hmix",
        "training-density",
        "validation-density",
    ):
        path = cache_dir / f"{name}.csv"
        if path.exists():
            frames.append(pd.read_csv(path))
    data = pd.concat(frames)
    data = data[data.forcefield == PUBLISHED_FF].copy()
    data["x_1"] = data.apply(published_x_first, axis=1)
    data["prop"] = data.type.map({"EnthalpyOfMixing": "hmix", "Density": "density"})
    return data


def redlich_kister_basis(x) -> np.ndarray:
    """
    Evaluate the first two Redlich-Kister basis functions, x(1-x) and x(1-x)(2x-1).

    Parameters
    ----------
    x
        Mole fraction(s).

    Returns
    -------
    np.ndarray
        Basis values with shape (len(x), 2).
    """
    x = np.atleast_1d(x)
    return np.stack([x * (1 - x), x * (1 - x) * (2 * x - 1)], axis=1)


def shift_to_composition(
    property_name: str, xs: np.ndarray, values: np.ndarray, x: float
) -> float:
    """
    Move the published value nearest ``x`` to ``x`` along a fitted curve.

    The curve only supplies the change between the nearest published
    composition and ``x`` (at most ``PUBLISHED_X_TOLERANCE``), so its shape
    matters little. For ΔHmix it is a two-term Redlich-Kister fit, which is
    zero at x = 0 and 1; a straight line to those limits misses the curvature
    near the ends. For density it is a polynomial of degree up to 2.

    Parameters
    ----------
    property_name
        ``hmix`` or ``density``.
    xs
        Published mole fractions of component a.
    values
        Published values.
    x
        Composition to shift to.

    Returns
    -------
    float
        Published value shifted to ``x``.
    """
    nearest = int(np.argmin(np.abs(xs - x)))
    if property_name == "hmix":
        n_terms = min(2, len(xs))
        coefficients, *_ = np.linalg.lstsq(
            redlich_kister_basis(xs)[:, :n_terms], values, rcond=None
        )
        fitted_at_x = float(redlich_kister_basis(x)[0, :n_terms] @ coefficients)
        fitted_at_nearest = float(
            redlich_kister_basis(xs[nearest])[0, :n_terms] @ coefficients
        )
    else:
        polynomial = np.polyfit(xs, values, min(2, len(xs) - 1))
        fitted_at_x = float(np.polyval(polynomial, x))
        fitted_at_nearest = float(np.polyval(polynomial, xs[nearest]))
    return float(values[nearest] + fitted_at_x - fitted_at_nearest)


def finite_size_check(summary: pd.DataFrame, systems, published) -> pd.DataFrame:
    """
    Compare our replicate means with published 1000-molecule Sage values.

    The published compositions differ from ours by up to
    ``PUBLISHED_X_TOLERANCE``, which can shift a density by more than the
    gaps being tested (about 0.2 g/mL per unit x for some systems), so the
    nearest published value is shifted to our composition with
    :func:`shift_to_composition`. Published pure densities are not used, as
    the published table has no temperatures.

    Parameters
    ----------
    summary
        Output of :func:`pooled_sigma`.
    systems
        ``systems.csv``.
    published
        Output of :func:`published_sage`.

    Returns
    -------
    pd.DataFrame
        One row per matched point.
    """
    rows = []
    components = systems.assign(
        box_id=[
            mixture_box_id(system_id, x_box)
            for system_id, x_box in zip(systems.system_id, systems.x_box, strict=True)
        ]
    ).set_index("box_id")
    for point in summary.itertuples():
        a, b = components.loc[point.box_id, ["a", "b"]]
        for property_name in ("hmix", "density"):
            values = published[published.prop == property_name]
            same_order = values[(values.smiles_1 == a) & (values.smiles_2 == b)]
            swapped = values[(values.smiles_1 == b) & (values.smiles_2 == a)]
            candidates = pd.concat(
                [
                    same_order.assign(xa=same_order.x_1),
                    swapped.assign(xa=1 - swapped.x_1),
                ]
            )
            if candidates.empty:
                continue
            nearest = candidates.iloc[(candidates.xa - point.x).abs().argmin()]
            if abs(nearest.xa - point.x) > PUBLISHED_X_TOLERANCE:
                continue
            published_at_x = shift_to_composition(
                property_name,
                candidates.xa.to_numpy(),
                candidates.value.to_numpy(),
                point.x,
            )
            ours = getattr(point, f"{property_name}_sage")
            # Standard errors of the two replicate means (published: 3 replicates).
            se_ours = getattr(point, f"{property_name}_sd") / np.sqrt(point.n_rep)
            se_published = nearest.uncertainty / np.sqrt(3)
            rows.append(
                {
                    "box_id": point.box_id,
                    "prop": property_name,
                    "x_ours": point.x,
                    "x_published": nearest.xa,
                    "ours": ours,
                    "published_nearest": nearest.value,
                    "published": published_at_x,
                    "diff": ours - published_at_x,
                    "z": (ours - published_at_x) / np.hypot(se_ours, se_published),
                }
            )
    return pd.DataFrame(rows)


def plot_convergence(runs_dir: Path, box_id: str, path: Path) -> None:
    """
    Plot block-averaged density and enthalpy against time for all replicates.

    Parameters
    ----------
    runs_dir
        Output directory of ``run_sage_benchmark.py``.
    box_id
        Box to plot.
    path
        Output PNG path.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
    for csv in sorted((runs_dir / box_id).glob("rep*.csv")):
        data = pd.read_csv(csv)
        block = ((data.time_ps - 1e-6) // 100).astype(int)  # (0, 100], (100, 200], ...
        means = data.groupby(block).agg(
            t=("time_ps", "mean"), rho=("density_g_cm3", "mean"), e=("epot_eV", "mean")
        )
        axes[0].plot(means.t, means.rho, marker=".", label=csv.stem)
        axes[1].plot(means.t, means.e - means.e.iloc[-10:].mean(), marker=".")
    for ax, label in zip(
        axes, ("density (g/mL), 100 ps blocks", "Epot - late mean (eV)"), strict=True
    ):
        ax.axvline(EQUILIBRATION_PS, color="k", ls="--", lw=0.8)
        ax.set_xlabel("time (ps)")
        ax.set_ylabel(label)
    axes[0].legend(fontsize=7)
    fig.suptitle(box_id)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_system(points: pd.DataFrame, summary: pd.DataFrame, path: Path) -> None:
    """
    Plot Sage ΔHmix and density against experiment for one system.

    Parameters
    ----------
    points
        Per-replicate rows of one system.
    summary
        Per-box rows of one system.
    path
        Output PNG path.
    """
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    for ax, prop, unit_label in (
        (axes[0], "hmix", "ΔHmix (kJ/mol)"),
        (axes[1], "density", "density (g/mL)"),
    ):
        ax.scatter(points.x, points[prop], s=10, alpha=0.5, label="Sage replicates")
        ax.errorbar(
            summary.x,
            summary[f"{prop}_exp"],
            yerr=summary[f"{prop}_exp_sigma"],
            fmt="ks",
            ms=4,
            label="experiment",
        )
        ax.set_xlabel("x (component a)")
        ax.set_ylabel(unit_label)
    axes[0].legend(fontsize=7)
    fig.suptitle(summary.system_id.iloc[0])
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def start_dependence(runs: pd.DataFrame) -> pd.DataFrame:
    """
    Compare replicate 0 (the distributed box) with the re-packed replicates.

    Replicates 1 and up are re-packed and re-equilibrated, so a significant
    offset of replicate 0 means the result depends on the starting
    configuration (for example a conformer population that does not relax
    within the run). The statistic is a prediction-interval t value:
    (r0 - mean) / (SD * sqrt(1 + 1/n)).

    Parameters
    ----------
    runs
        Per-run statistics.

    Returns
    -------
    pd.DataFrame
        One row per box and observable with the offset and its t value.
    """
    rows = []
    for box_id, group in runs.groupby("box_id"):
        others = group[group.replicate > 0]
        if group.replicate.min() != 0 or len(others) < 2:
            continue
        for column in ("h_mean", "density_mean"):
            first = group.loc[group.replicate == 0, column].iloc[0]
            spread = others[column].std(ddof=1) * np.sqrt(1 + 1 / len(others))
            offset = first - others[column].mean()
            rows.append(
                {
                    "box_id": box_id,
                    "observable": column,
                    "offset": offset,
                    "t": offset / spread,
                }
            )
    return pd.DataFrame(rows)


def within_between(runs: pd.DataFrame, points: pd.DataFrame) -> pd.DataFrame:
    """
    Compare the replicate spread with the within-run standard errors.

    If the within-run (autocorrelation) standard errors are right, the SD of
    the n replicate means follows (n - 1)·SD²/sem² ~ χ²(n - 1). The ratio
    SD_between / sem_within is tested two-sided, with Holm across all tests.
    sem_within is the RMS of the replicates' standard errors.

    Parameters
    ----------
    runs
        Per-run statistics.
    points
        Per-replicate mixture points.

    Returns
    -------
    pd.DataFrame
        One row per box and observable (``level`` "box") and per mixture point
        and property (``level`` "point").
    """
    rows = []
    # (level, table, observables, mean column, standard error column)
    tests = (
        ("box", runs, ("h", "density"), "{}_mean", "{}_sigma"),
        ("point", points, ("hmix", "vex"), "{}", "{}_sigma_run"),
    )
    for level, table, observables, mean_column, sem_column in tests:
        for box_id, group in table.groupby("box_id"):
            for observable in observables:
                means = group[mean_column.format(observable)]
                sems = group[sem_column.format(observable)]
                n = len(means)
                if n < 3:
                    continue
                sd = float(means.std(ddof=1))
                sem = float(np.sqrt(np.mean(sems**2)))
                chi2 = (n - 1) * (sd / sem) ** 2
                rows.append(
                    {
                        "level": level,
                        "box_id": box_id,
                        "observable": observable,
                        "n": n,
                        "sd_between": sd,
                        "sem_within": sem,
                        "ratio": sd / sem,
                        "chi2": chi2,
                        "p_high": float(stats.chi2.sf(chi2, n - 1)),
                        "p_low": float(stats.chi2.cdf(chi2, n - 1)),
                    }
                )
    table = pd.DataFrame(rows)
    p_two = np.minimum(1.0, 2 * np.minimum(table.p_high, table.p_low))
    table["p_adj"] = holm(p_two)
    table["flag"] = np.where(
        table.p_adj < FLAG_P, np.where(table.ratio > 1, "high", "low"), ""
    )
    return table


def pooled_ratio(table: pd.DataFrame) -> pd.DataFrame:
    """
    Pool the within/between test over all boxes or points of one observable.

    Parameters
    ----------
    table
        Output of :func:`within_between`.

    Returns
    -------
    pd.DataFrame
        Pooled ratio sqrt(Σχ² / Σdof), median ratio and the two-sided p-value.
    """
    rows = []
    for (level, observable), group in table.groupby(
        ["level", "observable"], sort=False
    ):
        degrees_of_freedom = float((group.n - 1).sum())
        chi2 = float(group.chi2.sum())
        p = 2 * min(
            stats.chi2.sf(chi2, degrees_of_freedom),
            stats.chi2.cdf(chi2, degrees_of_freedom),
        )
        rows.append(
            {
                "level": level,
                "observable": observable,
                "n_tests": len(group),
                "median_ratio": float(group.ratio.median()),
                "pooled_ratio": float(np.sqrt(chi2 / degrees_of_freedom)),
                "p": float(min(1.0, p)),
                "n_sig_high": int((group.p_high < FLAG_P).sum()),
                "n_sig_low": int((group.p_low < FLAG_P).sum()),
            }
        )
    return pd.DataFrame(rows)


def random_like_fraction(box: pd.Series) -> float:
    """
    Like-neighbour fraction expected for an ideally random mixture.

    Parameters
    ----------
    box
        Row of ``boxes.csv``.

    Returns
    -------
    float
        Molecule-averaged probability that a neighbour is the same component.
    """
    n_a, n_b = (int(c) for c in str(box.counts).split(";"))
    n = n_a + n_b
    return (n_a * (n_a - 1) + n_b * (n_b - 1)) / (n * (n - 1))


def plot_within_between(table: pd.DataFrame, path: Path) -> None:
    """
    Plot between-run SD against within-run standard error for every box.

    Parameters
    ----------
    table
        Output of :func:`within_between`.
    path
        Output PNG path.
    """
    boxes = table[table.level == "box"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, obs, unit in (
        (axes[0], "h", "h (kJ/mol)"),
        (axes[1], "density", "density (g/mL)"),
    ):
        data = boxes[boxes.observable == obs]
        n = int(data.n.median())
        lo, hi = data.sem_within.min() / 1.5, data.sem_within.max() * 1.5
        grid = np.geomspace(lo, hi, 50)
        ax.plot(grid, grid, "k-", lw=0.8, label="SD = sem")
        for q in (FLAG_P / 2, 1 - FLAG_P / 2):
            factor = np.sqrt(stats.chi2.ppf(q, n - 1) / (n - 1))
            ax.plot(grid, grid * factor, "k--", lw=0.6)
        flagged = data.p_high < FLAG_P
        ax.scatter(data.sem_within[~flagged], data.sd_between[~flagged], s=12)
        ax.scatter(
            data.sem_within[flagged],
            data.sd_between[flagged],
            s=16,
            c="tab:red",
            label=f"p_high < {FLAG_P}",
        )
        for r in data[flagged].itertuples():
            ax.annotate(r.box_id, (r.sem_within, r.sd_between), fontsize=6)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(f"within-run sem, RMS over replicates: {unit}")
        ax.set_ylabel(f"SD of replicate means: {unit}")
        ax.set_title(f"{obs}: dashed lines are the {1 - FLAG_P:.0%} χ²({n - 1}) range")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_replicates(runs_dir: Path, runs: pd.DataFrame, box_id: str, path: Path):
    """
    Plot replicate means and block-mean time series of one box.

    Parameters
    ----------
    runs_dir
        Output directory of ``run_sage_benchmark.py``.
    runs
        Per-run statistics.
    box_id
        Box to plot.
    path
        Output PNG path.
    """
    group = runs[runs.box_id == box_id].sort_values("replicate")
    n_mol = int(group.n_molecules.iloc[0])
    fig, axes = plt.subplots(2, 2, figsize=(10, 6))
    for col, obs, unit in ((0, "h", "h (kJ/mol)"), (1, "density", "density (g/mL)")):
        ax = axes[0, col]
        ax.errorbar(
            group.replicate,
            group[f"{obs}_mean"],
            yerr=2 * group[f"{obs}_sigma"],
            fmt="o",
            capsize=3,
        )
        ax.axhline(group[f"{obs}_mean"].mean(), color="k", lw=0.8)
        ax.set_xlabel("replicate")
        ax.set_ylabel(f"{unit}, mean ± 2 sem")
    for csv in sorted((runs_dir / box_id).glob("rep*.csv")):
        data = pd.read_csv(csv)
        block = ((data.time_ps - 1e-6) // 100).astype(int)
        # Enthalpy per molecule, (E_pot + pV) / N in kJ/mol, at 1 atm.
        h = (
            data.epot_eV * KJ_PER_MOL_PER_EV
            + 101325.0 * data.volume_A3 * 1e-30 * AVOGADRO / 1000
        ) / n_mol
        means = (
            data.assign(h=h)
            .groupby(block)
            .agg(t=("time_ps", "mean"), h=("h", "mean"), rho=("density_g_cm3", "mean"))
        )
        axes[1, 0].plot(means.t, means.h, marker=".", label=csv.stem)
        axes[1, 1].plot(means.t, means.rho, marker=".")
    for ax, unit in ((axes[1, 0], "h (kJ/mol)"), (axes[1, 1], "density (g/mL)")):
        ax.axvline(EQUILIBRATION_PS, color="k", ls="--", lw=0.8)
        ax.set_xlabel("time (ps)")
        ax.set_ylabel(f"{unit}, 100 ps blocks")
    axes[1, 0].legend(fontsize=7)
    fig.suptitle(box_id)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_finite_size(finite: pd.DataFrame, path: Path) -> None:
    """
    Plot our Sage values minus the published 1000-molecule values.

    Parameters
    ----------
    finite
        Output of ``finite_size_check``.
    path
        Output PNG path.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, prop, unit in (
        (axes[0], "hmix", "ΔHmix (kJ/mol)"),
        (axes[1], "density", "density (g/mL)"),
    ):
        data = finite[finite.prop == prop].reset_index(drop=True)
        se = (data["diff"] / data.z).abs()
        ax.errorbar(data["diff"], np.arange(len(data)), xerr=2 * se, fmt="o", ms=3)
        ax.axvline(0, color="k", lw=0.8)
        for sign in (-1, 1):
            ax.axvline(sign * SIZE_GAP[prop], color="tab:red", ls=":", lw=0.8)
        ax.set_yticks(np.arange(len(data)), data.box_id, fontsize=6)
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
        ax.set_xlabel(f"ours - published: {unit} (± 2 se)")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_replicate_spread(summary: pd.DataFrame, path: Path) -> None:
    """
    Plot per-point replicate SDs against the pooled σ_Sage and σ_ref per system.

    Parameters
    ----------
    summary
        Output of :func:`pooled_sigma`.
    path
        Output PNG path.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    systems = summary.groupby("system_id").sigma_sage_hmix.first().sort_values()
    rows = summary.system_id.map({s: i for i, s in enumerate(systems.index)})
    for ax, prop, unit in (
        (axes[0], "hmix", "ΔHmix (kJ/mol)"),
        (axes[1], "density", "density (g/mL)"),
    ):
        ax.scatter(summary[f"{prop}_sd"], rows, s=12, label="replicate SD at one x")
        ax.scatter(
            summary[f"sigma_sage_{prop}"],
            rows,
            marker="|",
            s=150,
            c="k",
            label="σ_Sage",
        )
        ax.scatter(
            summary[f"sigma_ref_{prop}"],
            rows,
            s=30,
            facecolors="none",
            edgecolors="tab:red",
            label="σ_ref",
        )
        ax.axvline(
            summary[f"sigma_sage_{prop}"].min(),
            color="k",
            ls=":",
            lw=0.8,
            label="floor (median system)",
        )
        ax.set_yticks(np.arange(len(systems)), systems.index, fontsize=7)
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
        ax.set_xlabel(f"SD of replicate means: {unit}")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_sage_vs_exp(summary: pd.DataFrame, path: Path) -> None:
    """
    Parity plot of Sage replicate-mean ΔHmix against experiment.

    Parameters
    ----------
    summary
        Per-point summary.
    path
        Output PNG path.
    """
    fig, ax = plt.subplots(figsize=(5.5, 5))
    for system_id, group in summary.groupby("system_id"):
        ax.errorbar(
            group.hmix_exp,
            group.hmix_sage,
            xerr=group.hmix_exp_sigma,
            yerr=group.hmix_sd,
            fmt="o",
            ms=4,
            label=system_id,
        )
    lims = [
        min(summary.hmix_exp.min(), summary.hmix_sage.min()) - 0.2,
        max(summary.hmix_exp.max(), summary.hmix_sage.max()) + 0.2,
    ]
    ax.plot(lims, lims, "k-", lw=0.8)
    ax.set_xlabel("experimental ΔHmix (kJ/mol)")
    ax.set_ylabel("Sage ΔHmix (kJ/mol), mean ± replicate SD")
    ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def inspection_items(
    runs, summary, finite, run_flags, start, within, boxes
) -> pd.DataFrame:
    """
    Collect boxes and systems worth inspecting, by rule.

    Parameters
    ----------
    runs
        Per-run statistics.
    summary
        Per-point summary.
    finite
        Finite-size comparison.
    run_flags
        Holm-flagged drift and half-split tests as (box_id, text).
    start
        Start-dependence table.
    within
        Within/between table.
    boxes
        ``boxes.csv``.

    Returns
    -------
    pd.DataFrame
        Items with ``category``, ``target``, ``evidence`` and ``box`` (the box
        to plot, or ``None`` for a system-level item).
    """
    items = []
    mismatch = within[(within.flag != "") | (within.p_high < FLAG_P)]
    # Holm-significant first, then by unadjusted p (many adjusted p tie at 1)
    mismatch = mismatch.assign(unflagged=mismatch.flag == "").sort_values(
        ["unflagged", "p_high"]
    )
    for r in mismatch.itertuples():
        strength = "Holm-significant" if r.flag else f"p < {FLAG_P} unadjusted only"
        items.append(
            {
                "category": "within/between mismatch",
                "target": r.box_id,
                "evidence": f"{r.observable}: replicate SD {r.ratio:.1f}× the "
                f"within-run sem (p_high = {r.p_high:.1g}, {strength})",
                "box": r.box_id,
            }
        )
    low = runs[(runs.h_n_eff < MIN_N_EFF) | (runs.density_n_eff < MIN_N_EFF)]
    for box_id, group in low.groupby("box_id"):
        n_eff = group[["h_n_eff", "density_n_eff"]].min(axis=1)
        items.append(
            {
                "category": "low sampling",
                "target": box_id,
                "evidence": f"{len(group)} replicate(s) with n_eff < {MIN_N_EFF} "
                f"(minimum {n_eff.min():.0f}, τ up to "
                f"{PRODUCTION_PS / n_eff.min():.0f} ps)",
                "box": box_id,
            }
        )
    for box_id, text in run_flags:
        items.append(
            {
                "category": "drift or half-split",
                "target": box_id,
                "evidence": text,
                "box": box_id,
            }
        )
    for r in start[start.t.abs() > START_T].itertuples():
        items.append(
            {
                "category": "start dependence",
                "target": r.box_id,
                "evidence": f"{r.observable}: replicate 0 offset {r.offset:+.4f} "
                f"(t = {r.t:+.1f})",
                "box": r.box_id,
            }
        )
    for r in runs[runs.density_first_block_shift.abs() > EARLY_SHIFT].itertuples():
        items.append(
            {
                "category": "early density shift",
                "target": r.box_id,
                "evidence": f"rep{r.replicate}: 200–400 ps density differs from "
                f"the rest by {r.density_first_block_shift:+.4f} g/mL",
                "box": r.box_id,
            }
        )
    if "reacted" in runs:
        for r in runs[runs.reacted.fillna(False).astype(bool)].itertuples():
            items.append(
                {
                    "category": "bonding changed",
                    "target": r.box_id,
                    "evidence": f"rep{r.replicate}: molecules differ from the "
                    "input box",
                    "box": r.box_id,
                }
            )
        random = (
            boxes[boxes.kind == "mixture"]
            .set_index("box_id")
            .apply(random_like_fraction, axis=1)
        )
        mixed = runs[runs.kind == "mixture"].groupby("box_id").like_fraction.mean()
        for box_id, value in mixed.items():
            if value - random[box_id] > DEMIX_EXCESS:
                items.append(
                    {
                        "category": "microheterogeneity",
                        "target": box_id,
                        "evidence": f"like-neighbour fraction {value:.3f} vs "
                        f"{random[box_id]:.3f} for random mixing",
                        "box": box_id,
                    }
                )
    for r in finite.itertuples():
        if abs(r.diff) > SIZE_GAP[r.prop]:
            items.append(
                {
                    "category": "published-Sage gap",
                    "target": r.box_id,
                    "evidence": f"{r.prop}: ours − published {r.diff:+.4f} "
                    f"(z = {r.z:+.1f})",
                    "box": None,
                }
            )
    for system_id, group in summary.groupby("system_id"):
        error = (group.hmix_sage - group.hmix_exp).abs().mean()
        wrong_sign = (np.sign(group.hmix_sage) != np.sign(group.hmix_exp)).sum()
        if error > FF_HMIX_ERROR or wrong_sign:
            sign_text = f", wrong sign at {wrong_sign} of {len(group)} x" * bool(
                wrong_sign
            )
            items.append(
                {
                    "category": "force-field miss (not a run problem)",
                    "target": system_id,
                    "evidence": f"mean |ΔHmix error| {error:.2f} kJ/mol{sign_text}",
                    "box": None,
                }
            )
    return pd.DataFrame(items)


def main() -> None:
    """Summarise the Sage runs and update ``systems.csv`` and ``boxes.csv``."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--runs-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path.home() / ".cache/ml_peg/binary_mixtures_selection",
    )
    parser.add_argument("--replicates", type=int, default=5)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    out = args.output_dir
    (out / "figures" / "convergence").mkdir(parents=True, exist_ok=True)
    (out / "figures" / "systems").mkdir(parents=True, exist_ok=True)

    systems = pd.read_csv(args.data_dir / "systems.csv")
    boxes = pd.read_csv(args.data_dir / "boxes.csv")
    runs = load_runs(args.data_dir, args.runs_dir, args.workers)
    runs.to_csv(out / "runs.csv", index=False)
    points = mixture_points(runs, systems)
    points.to_csv(out / "replicate_points.csv", index=False)
    summary = pooled_sigma(points, args.replicates)
    summary.to_csv(out / "points.csv", index=False)
    pure = pd.read_csv(args.data_dir / "pure.csv")
    pure_summary = pure_sigma(runs, pure)
    pure_summary.to_csv(out / "pure_points.csv", index=False)
    finite = finite_size_check(summary, systems, published_sage(args.cache_dir))
    finite.to_csv(out / "finite_size.csv", index=False)
    within = within_between(runs, points)
    within.to_csv(out / "within_between.csv", index=False)
    pooled = pooled_ratio(within)

    tests = [
        (
            r.box_id,
            f"{r.box_id} rep{r.replicate} {obs} {kind}",
            getattr(r, f"{obs}_p_{kind}"),
        )
        for r in runs.itertuples()
        for obs in ("h", "density")
        for kind in ("drift", "half")
    ]
    adjusted = holm(np.array([p for *_, p in tests]))
    flags = [
        (box_id, f"{name} (p_adj = {p:.2g})")
        for (box_id, name, _), p in zip(tests, adjusted, strict=True)
        if p < FLAG_P
    ]
    for box_id in runs.box_id.unique():
        plot_convergence(
            args.runs_dir, box_id, out / "figures/convergence" / f"{box_id}.png"
        )
    for system_id, group in points.groupby("system_id"):
        plot_system(
            group,
            summary[summary.system_id == system_id],
            out / "figures/systems" / f"{system_id}.png",
        )
    plot_within_between(within, out / "figures" / "within_between.png")
    plot_finite_size(finite, out / "figures" / "finite_size.png")
    plot_sage_vs_exp(summary, out / "figures" / "sage_vs_exp.png")
    plot_replicate_spread(summary, out / "figures" / "replicate_spread.png")
    start = start_dependence(runs)
    start.to_csv(out / "start_dependence.csv", index=False)
    items = inspection_items(runs, summary, finite, flags, start, within, boxes)
    (out / "figures" / "replicates").mkdir(exist_ok=True)
    for box_id in sorted(items.box.dropna().unique()):
        plot_replicates(
            args.runs_dir, runs, box_id, out / "figures/replicates" / f"{box_id}.png"
        )
    from liquid_mixtures_inputs.sage_report import write_reports

    write_reports(
        out,
        {
            "runs": runs,
            "summary": summary,
            "pure_summary": pure_summary,
            "finite": finite,
            "within": within,
            "pooled": pooled,
            "start": start,
            "flags": flags,
            "items": items,
        },
    )

    if not (summary.complete.all() and len(summary) == len(systems)):
        return
    cols = [
        "box_id",
        "hmix_sage",
        "density_sage",
        "vex_sage",
        "sigma_sage_hmix",
        "sigma_sage_density",
        "sigma_sage_vex",
        "sigma_ref_hmix",
        "sigma_ref_density",
        "sigma_ref_vex",
    ]
    updated = systems.assign(
        box_id=[
            mixture_box_id(s, x)
            for s, x in zip(systems.system_id, systems.x_box, strict=True)
        ]
    )
    updated = updated.drop(columns=[c for c in cols[1:] if c in updated]).merge(
        summary[cols], on="box_id", how="left"
    )
    updated.to_csv(args.data_dir / "systems.csv", index=False)

    pure_columns = ["density_sage", "sigma_sage_density", "sigma_ref_density"]
    pure = pure.drop(columns=[c for c in pure_columns if c in pure]).merge(
        pure_summary[["smiles", "temperature", *pure_columns]],
        on=["smiles", "temperature"],
        how="left",
    )
    pure.to_csv(args.data_dir / "pure.csv", index=False)

    like = (
        runs[runs.kind == "mixture"]
        .groupby("box_id")
        .like_fraction.agg(like_fraction_sage="mean", like_fraction_sage_sd="std")
        .reset_index()
    )
    boxes = boxes.drop(
        columns=[
            c for c in ("like_fraction_sage", "like_fraction_sage_sd") if c in boxes
        ]
    ).merge(like, on="box_id", how="left")
    boxes.to_csv(args.data_dir / "boxes.csv", index=False)


if __name__ == "__main__":
    main()
