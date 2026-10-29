"""
Write the power-analysis reports from the tables of ``power_analysis.py``.

``report.md`` gives the answer, the cost, the thresholds and the main caveat;
``details.md`` holds the method, the full tables, the literature and every
caveat.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from liquid_mixtures_inputs.power_analysis import (
    ALTERNATIVE_PS,
    BAD_SAGE_MULTIPLE,
    GOOD_NOISE_SDS,
    LABELS,
    LITERATURE,
    MIN_ALLOCATION_PS,
    N_BOOTSTRAP,
    N_DRAWS,
    PRODUCTION_LENGTHS_PS,
    PROFILES,
    REFERENCE,
    Z_DETECT,
    round_significant,
)

# MACE-OFF23(M) on one A100 via OpenMM: 1.1e6 steps/day for ~600-atom liquid boxes
# (arXiv:2312.15211, section III.3.2). Used only for an order-of-magnitude cost.
MACE_OFF_M_ATOM_STEPS_PER_DAY = 1.1e6 * 600
# Spread of published MLP pure-density MAEs (LITERATURE: MACE-OFF23 S vs M,
# ANI-2x vs MACE-OFF23(M), MPNICE variants), in g/mL.
LITERATURE_DENSITY_SPREAD = (0.05, 0.25)
METRIC_ORDER = [
    ("hmix", "MAE"),
    ("hmix", "MAZE"),
    ("density", "MAE"),
    ("density", "MAZE"),
    ("vex", "MAE"),
    ("vex", "MAZE"),
    ("pure_density", "MAE"),
    ("pure_density", "MAZE"),
]
UNITS = {
    "hmix": "kJ/mol",
    "density": "g/mL",
    "vex": "cm³/mol",
    "pure_density": "g/mL",
}


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


def metric_value(
    metrics: pd.DataFrame,
    property_name: str,
    metric: str,
    production: float,
    profile: float,
    column: str,
) -> float:
    """
    Look up one value of the metrics table.

    Parameters
    ----------
    metrics
        ``metrics`` table from ``power_analysis.analyse``.
    property_name
        Property, e.g. ``hmix``.
    metric
        ``MAE`` or ``MAZE``.
    production
        Production length in ps.
    profile
        Error profile (multiple of Sage's errors).
    column
        Column to return, e.g. ``resolvable``.

    Returns
    -------
    float
        The value.
    """
    row = metrics[
        (metrics.property == property_name)
        & (metrics.metric == metric)
        & (metrics.production_ps == production)
        & (metrics.profile == profile)
    ]
    return float(row[column].iloc[0])


def resolvable_range(
    metrics: pd.DataFrame, property_name: str, metric: str, production: float
) -> tuple[float, float]:
    """
    Find the range of the resolvable difference over the error profiles.

    Parameters
    ----------
    metrics
        ``metrics`` table from ``power_analysis.analyse``.
    property_name
        Property, e.g. ``hmix``.
    metric
        ``MAE`` or ``MAZE``.
    production
        Production length in ps.

    Returns
    -------
    tuple[float, float]
        Smallest and largest resolvable difference.
    """
    rows = metrics[
        (metrics.property == property_name)
        & (metrics.metric == metric)
        & (metrics.production_ps == production)
    ].resolvable
    return float(rows.min()), float(rows.max())


def point_rows(
    points: pd.DataFrame, property_name: str, production: float
) -> pd.DataFrame:
    """
    Select the per-point rows of one property and run length.

    Parameters
    ----------
    points
        ``points`` table from ``power_analysis.analyse``.
    property_name
        Property, e.g. ``hmix``.
    production
        Production length in ps.

    Returns
    -------
    pd.DataFrame
        Matching rows.
    """
    return points[
        (points.property == property_name) & (points.production_ps == production)
    ]


def costs(n_atoms: np.ndarray) -> dict[float, float]:
    """
    Compute the cost (atoms × simulated ps) of each protocol.

    Parameters
    ----------
    n_atoms
        Atoms per box.

    Returns
    -------
    dict[float, float]
        Cost per production length, with the standard equilibration.
    """
    return {
        production: float(n_atoms.sum()) * (REFERENCE[0] + production)
        for production in PRODUCTION_LENGTHS_PS
    }


def allocation_cost(allocation: pd.DataFrame) -> float:
    """
    Compute the cost of the per-box allocation.

    Parameters
    ----------
    allocation
        ``allocation`` table from ``power_analysis.analyse``.

    Returns
    -------
    float
        Atoms × simulated ps.
    """
    return float(
        np.nansum(allocation.n_atoms * (REFERENCE[0] + allocation.production_ps))
    )


def mean_and_sd(mean: float, sd: float) -> str:
    """
    Format a mean and SD to the same decimal places (2 significant in the mean).

    Parameters
    ----------
    mean
        Mean, e.g. 0.0713.
    sd
        Standard deviation, e.g. 0.0105.

    Returns
    -------
    str
        For example ``"0.071 ± 0.010"``.
    """
    decimals = max(0, 1 - int(np.floor(np.log10(abs(mean)))))
    return f"{mean:.{decimals}f} ± {sd:.{decimals}f}"


def threshold_rows(threshold_table: pd.DataFrame) -> list[str]:
    """
    Format the thresholds as a Markdown table.

    Parameters
    ----------
    threshold_table
        Output of ``power_analysis.thresholds``.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    table = threshold_table.set_index(["property", "metric"])
    rows = []
    for property_name, metric in METRIC_ORDER:
        row = table.loc[(property_name, metric)]
        unit = UNITS[property_name] if metric == "MAE" else ""
        rows.append(
            [
                f"{LABELS[property_name]} {metric}",
                mean_and_sd(row.perfect_mean, row.perfect_sd),
                f"{row.sage:.3g}",
                f"{row.good:g}",
                f"{row.bad:g}",
                unit,
            ]
        )
    frame = pd.DataFrame(
        rows,
        columns=[
            "Metric",
            f"Perfect model ({REFERENCE[1] / 1000:g} ns)",
            "Sage",
            "Good",
            "Bad",
            "Unit",
        ],
    )
    return [md_table(frame)]


def write_summary(out: Path, tables: dict) -> None:
    """
    Write the short ``report.md``.

    Parameters
    ----------
    out
        Output directory (figures are in ``out / "figures"``).
    tables
        The tables named in :func:`write_reports`.
    """
    metrics, points = tables["metrics"], tables["points"]
    reference = REFERENCE[1]
    cost = costs(tables["n_atoms"])
    hmix_low, hmix_high = resolvable_range(metrics, "hmix", "MAE", reference)
    density_low, density_high = resolvable_range(metrics, "density", "MAE", reference)
    hmix_points = point_rows(points, "hmix", reference)
    density_points = point_rows(points, "density", reference)
    saved = 1 - allocation_cost(tables["allocation"]) / cost[reference]
    floor = metric_value(metrics, "hmix", "MAZE", reference, 0.0, "observed_mean")
    floor_longer = metric_value(
        metrics, "hmix", "MAZE", ALTERNATIVE_PS, 0.0, "observed_mean"
    )
    reference_ns, longer_ns = reference / 1000, ALTERNATIVE_PS / 1000
    n_boxes = len(tables["n_atoms"])
    n_runs = len(tables["shifts"].xs("density", level="observable"))
    protocol = f"{REFERENCE[0] / 1000:g} + {reference_ns:g} ns"
    longer = f"{REFERENCE[0] / 1000:g} + {longer_ns:g} ns"

    lines = [
        "# Power analysis: can one run per box rank MLPs?",
        "",
        f"Generated {date.today()} by `power_analysis.py`. Method, full tables, "
        "literature and caveats: [details.md](details.md). Thresholds: "
        "[thresholds.csv](thresholds.csv), ready to paste into ml-peg as "
        "[metrics_snippet.yml](metrics_snippet.yml).",
        "",
        f"ml-peg runs every MLP once per box: {REFERENCE[0] / 1000:g} ns of NPT "
        f"equilibration and {reference / 1000:g} ns of production at a 1 fs time "
        f"step. The {n_runs} Sage runs ({n_boxes} boxes × {n_runs // n_boxes} "
        "replicates, "
        "[Sage report](../sage_summary/report.md)) stand in for MLP simulation "
        "noise. From them we estimate what one run per box can resolve, what it "
        "costs, and where ml-peg's 'good' and 'bad' thresholds should sit.",
        "",
        "## Answer",
        "",
        "- **Between models.** Two models are told apart (two-sided 5 % test, "
        f"80 % power) when their ΔHmix MAE differs by {hmix_low:.2f} to "
        f"{hmix_high:.3f} kJ/mol, or their mixture density MAE by "
        f"{density_low:.4f} to {density_high:.4f} g/mL. Published MLPs differ in "
        f"pure-liquid density MAE by about {LITERATURE_DENSITY_SPREAD[0]} to "
        f"{LITERATURE_DENSITY_SPREAD[1]} g/mL, about "
        f"{round_significant(LITERATURE_DENSITY_SPREAD[0] / density_high, 1):g} to "
        f"{round_significant(LITERATURE_DENSITY_SPREAD[1] / density_high, 1):g}× "
        "more.",
        "- **Against experiment.** One run per box detects an error at a single "
        f"point above {hmix_points.mde.median():.2f} kJ/mol for ΔHmix and "
        f"{density_points.mde.median():.4f} g/mL for density (median over "
        "points). "
        f"{(hmix_points.sage_error.abs() > hmix_points.mde).mean():.0%} of Sage's "
        "ΔHmix errors are larger.",
        f"- **Cost.** {tables['n_atoms'].sum():,.0f} atoms over {n_boxes} boxes, about "
        f"{gpu_days(cost[reference]):.0f} MACE-OFF23(M) A100 GPU-days per model "
        f"(order of magnitude). {longer} would cost "
        f"{cost[ALTERNATIVE_PS] / cost[reference]:.1f}× as much (about "
        f"{gpu_days(cost[ALTERNATIVE_PS] - cost[reference]):.0f} GPU-days more) "
        "and narrow the resolvable ΔHmix MAE difference from "
        f"{hmix_high:.3f} to "
        f"{resolvable_range(metrics, 'hmix', 'MAE', ALTERNATIVE_PS)[1]:.3f} "
        f"kJ/mol. Per-box run lengths would save only {saved:.0%}.",
        "",
        f"We run {protocol} for cost: it already resolves density gaps far "
        f"smaller than those between published MLPs, and {longer} would narrow "
        "the resolvable ΔHmix MAE difference by only "
        f"{1 - resolvable_range(metrics, 'hmix', 'MAE', ALTERNATIVE_PS)[1] / hmix_high:.0%}"
        f" for {cost[ALTERNATIVE_PS] / cost[reference]:.1f}× the cost. "
        "Equilibration stays at "
        f"{REFERENCE[0] / 1000:g} ns ([details](details.md#4-run-length-and-cost)).",
        "",
        "![resolvable difference](figures/resolvable.png)",
        "",
        "Blue: resolvable difference for a Sage-like error profile (band: 1 to 0 "
        "× Sage's errors). Green: the score a model that matches experiment exactly gets from noise "
        "alone. Grey: the 'good' threshold.",
        "",
        "## Thresholds",
        "",
        "Every property is scored by MAZE (each error divided by a fixed σ_ref = "
        "sqrt(σ_Sage² + σ_exp²)), with the MAE shown at weight 0. Thresholds come "
        "from the error analysis: **good** is the score of a model that matches "
        f"experiment exactly, with one run per box, plus {GOOD_NOISE_SDS:g} SD of "
        "its noise, so such a model scores 'good' about 98 % of the time. **Bad** "
        f"is {BAD_SAGE_MULTIPLE:g}× Sage's score, as Sage was trained on these "
        "properties. Mixture and pure densities share one set of thresholds for "
        "simplicity, the larger of their two values, so both rules still hold "
        "for each.",
        "",
        *threshold_rows(tables["thresholds"]),
        "",
        "The perfect-model floor depends on the run length (ΔHmix MAZE "
        f"{floor:.2f} at {reference_ns:g} ns, {floor_longer:.2f} at {longer_ns:g} "
        "ns), so the thresholds must be recomputed if the protocol changes.",
        "",
        "## Main caveat",
        "",
        "Sage stands in for MLP noise. MLP fluctuations of h and volume should be "
        "similar (the same liquids at the same T and P), but correlation times may "
        "differ, especially for MLPs that change conformer barriers. One pilot MLP "
        "run per box, analysed with `mixing.box_result`, would check this: compare "
        "its sem and n_eff with the Sage `runs.csv`. Other caveats are in "
        "[details.md](details.md#caveats).",
        "",
    ]
    (out / "report.md").write_text("\n".join(lines))


def method_section(tables: dict) -> list[str]:
    """
    Write the method section of ``details.md``.

    Parameters
    ----------
    tables
        The tables named in :func:`write_reports`.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    reference = REFERENCE[1]
    n_mixtures = point_rows(tables["points"], "hmix", reference).system_id.nunique()
    n_pure = len(point_rows(tables["points"], "pure_density", reference))
    noise = tables["noise"]
    standard = noise[noise.equilibration_ps == REFERENCE[0]]
    k = standard.groupby("observable").k.agg(["min", "max"])
    propagated = point_rows(tables["points"], "hmix", reference).sigma_sim
    observed = tables["sage_points"].hmix_sd
    lengths = ", ".join(f"{p / 1000:g}" for p in PRODUCTION_LENGTHS_PS)
    profiles = ", ".join(f"{p:g}" for p in PROFILES)
    return [
        "## Method",
        "",
        "1. **Run length.** Every Sage run is re-analysed with the benchmark's "
        "estimator (`mixing.observables` and `timeseries.estimate_mean`, as in "
        f"`summarise-sage`), keeping {lengths} ns of production after 0.1 or "
        f"0.2 ns of equilibration. The analysis at the protocol length "
        f"({reference / 1000:g} ns) reproduces the Sage `runs.csv` exactly "
        "(tested).",
        "2. **Single-run noise per box.** σ = k·sem, with sem the RMS within-run "
        "standard error over the 5 replicates and k the pooled ratio of replicate "
        "SD to within-run sem over all boxes at that length (k = "
        f"{k.loc['h', 'min']:.2f} to {k.loc['h', 'max']:.2f} for h, "
        f"{k.loc['density', 'min']:.2f} to {k.loc['density', 'max']:.2f} for "
        "density). This is smoother than the 5-replicate SD of each box and "
        "corrects the within-run sem for the slow modes it misses.",
        "3. **Points.** ΔHmix and V^E are linear in the box means (ml-peg's "
        "`mixing_matrix`). Each pure box feeds all three compositions of its "
        "system, so noise is drawn per box and propagated, which keeps the "
        "correlation between points. The propagated ΔHmix noise at "
        f"{reference / 1000:g} ns (RMS over "
        f"points {np.sqrt((propagated**2).mean()):.3f} kJ/mol) matches the "
        "observed SD of Sage's replicate ΔHmix (RMS "
        f"{np.sqrt((observed**2).mean()):.3f} kJ/mol, Sage `points.csv`).",
        f"4. **Metrics.** For each run length, {N_DRAWS} draws of one run per box "
        "give the noise of every metric for a model whose true errors are "
        f"{profiles} × Sage's (replicate mean minus experiment). Two models run "
        "independently, so their metrics differ by noise with SD "
        "sqrt(2)·SD(metric). The **resolvable difference** is "
        f"{Z_DETECT:.2f}·sqrt(2)·SD(metric): a true difference that size is "
        "detected at the two-sided 5 % level 80 % of the time. A closed form "
        "sqrt(sᵀCs)/N (s the signs of the errors over σ_ref for MAZE, C the "
        "point covariance) is given as a check; it holds when errors are well "
        "above the noise.",
        f"5. **System choice.** A cluster bootstrap over systems ({N_BOOTSTRAP} "
        f"resamples of the {n_mixtures} mixtures, or of the {n_pure} scored pure "
        "liquids) of "
        "Sage's noiseless metric measures how much the metric depends on which "
        "systems were picked.",
        "",
    ]


def noise_section(tables: dict) -> list[str]:
    """
    Write the noise-against-run-length section of ``details.md``.

    Parameters
    ----------
    tables
        The tables named in :func:`write_reports`.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    noise, points = tables["noise"], tables["points"]
    standard = noise[noise.equilibration_ps == REFERENCE[0]]
    k = standard.groupby(["production_ps", "observable"]).k.first()
    h_noise = standard[standard.observable == "h"].groupby("production_ps").sigma
    h_rms = h_noise.apply(lambda sigma: float(np.sqrt(np.mean(sigma**2))))
    scaled = h_rms * np.sqrt(h_rms.index / REFERENCE[1])
    growth = scaled[ALTERNATIVE_PS] / scaled[500.0] - 1
    reference_ns, longer_ns = REFERENCE[1] / 1000, ALTERNATIVE_PS / 1000
    rows = []
    for production in PRODUCTION_LENGTHS_PS:
        median_sigma = {
            name: point_rows(points, name, production).sigma_sim.median()
            for name in LABELS
        }
        rows.append(
            [
                f"{production / 1000:g}",
                f"{median_sigma['hmix']:.3f}",
                f"{median_sigma['density']:.5f}",
                f"{median_sigma['vex']:.3f}",
                f"{median_sigma['pure_density']:.5f}",
                f"{k.loc[(production, 'h')]:.2f}",
                f"{k.loc[(production, 'density')]:.2f}",
            ]
        )
    columns = [
        "Production (ns)",
        "σ ΔHmix (kJ/mol)",
        "σ density (g/mL)",
        "σ V^E (cm³/mol)",
        "σ pure density (g/mL)",
        "k (h)",
        "k (density)",
    ]
    return [
        "## 1. Noise against run length",
        "",
        "![noise scaling](figures/noise_scaling.png)",
        "",
        f"Single-run noise × sqrt(t / {reference_ns:g} ns), RMS over boxes and "
        f"relative to {reference_ns:g} ns. A flat line would mean σ ∝ 1/sqrt(t). "
        f"For h the scaled noise grows by {growth:.0%} from 0.5 to {longer_ns:g} "
        "ns: short runs miss some slow fluctuations, "
        "so cutting the run length costs a little less precision than 1/sqrt(t) "
        "predicts.",
        "",
        md_table(pd.DataFrame(rows, columns=columns)),
        "",
        "Median single-run noise per point. k is the replicate SD / within-run sem "
        "ratio (1 if the within-run sem were right).",
        "",
    ]


def point_section(tables: dict) -> list[str]:
    """
    Write the per-point detection section of ``details.md``.

    Parameters
    ----------
    tables
        The tables named in :func:`write_reports`.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    points = tables["points"]
    reference = REFERENCE[1]
    ratios = [
        point_rows(points, name, reference).sigma_sim.median()
        / point_rows(points, name, reference).sigma_exp.median()
        for name in ("hmix", "density")
    ]
    rows = []
    for property_name in ("hmix", "density", "vex", "pure_density"):
        for production in (500.0, REFERENCE[1], ALTERNATIVE_PS):
            selected = point_rows(points, property_name, production)
            rows.append(
                [
                    LABELS[property_name],
                    f"{production / 1000:g}",
                    f"{selected.mde.median():.3g}",
                    f"{selected.mde.max():.3g}",
                    f"{(selected.sage_error.abs() > selected.mde).mean():.0%}",
                    f"{selected.sigma_exp.median():.2g}",
                    UNITS[property_name],
                ]
            )
    columns = [
        "Property",
        "Production (ns)",
        "Median MDE",
        "Max MDE",
        "Sage errors above MDE",
        "Median σ_exp",
        "Unit",
    ]
    return [
        "## 2. A model against experiment, per point",
        "",
        f"Minimum detectable error MDE = {Z_DETECT:.2f}·sqrt(σ_sim² + σ_exp²): an "
        "error that size is significant at the 5 % level 80 % of the time. The "
        "fraction column shows how many of Sage's errors (a realistic model) "
        "exceed the MDE. Experimental errors are "
        f"{min(ratios):.0f} to {max(ratios):.0f}× smaller than the simulation "
        "noise (median σ_sim over median σ_exp, ΔHmix and density) and barely "
        "change these limits.",
        "",
        md_table(pd.DataFrame(rows, columns=columns)),
        "",
        "![per-point detection](figures/point_mde.png)",
        "",
        "ΔHmix points above the diagonal have a Sage error that one run per box "
        "detects.",
        "",
    ]


def ranking_section(tables: dict) -> list[str]:
    """
    Write the section on ranking two models for ``details.md``.

    Parameters
    ----------
    tables
        The tables named in :func:`write_reports`.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    metrics = tables["metrics"]
    bootstrap = tables["bootstrap"].set_index(["property", "metric"])
    reference = REFERENCE[1]
    n_systems = point_rows(tables["points"], "hmix", reference).system_id.nunique()
    rows = []
    for property_name, metric in METRIC_ORDER:
        low, high = resolvable_range(metrics, property_name, metric, reference)
        low_longer, high_longer = resolvable_range(
            metrics, property_name, metric, ALTERNATIVE_PS
        )
        sage_like = metrics[
            (metrics.property == property_name)
            & (metrics.metric == metric)
            & (metrics.production_ps == reference)
            & (metrics.profile == 1.0)
        ].iloc[0]
        systems_sd = bootstrap.loc[(property_name, metric), "bootstrap_sd"]
        rows.append(
            [
                f"{LABELS[property_name]} {metric}",
                f"{sage_like.sd:.2g}",
                f"{sage_like.sd_closed_form:.2g}",
                f"{low:.2g} to {high:.2g}",
                f"{low_longer:.2g} to {high_longer:.2g}",
                f"{systems_sd:.2g}",
            ]
        )
    columns = [
        "Metric",
        f"SD ({reference / 1000:g} ns)",
        "SD closed form",
        f"Resolvable ({reference / 1000:g} ns)",
        f"Resolvable ({ALTERNATIVE_PS / 1000:g} ns)",
        "SD over systems",
    ]
    hmix_systems = float(bootstrap.loc[("hmix", "MAE"), "bootstrap_sd"])
    sage_mae = metric_value(metrics, "hmix", "MAE", reference, 1.0, "true")
    quarter_true = metric_value(metrics, "hmix", "MAE", reference, 0.25, "true")
    quarter_observed = metric_value(
        metrics, "hmix", "MAE", reference, 0.25, "observed_mean"
    )
    full_observed = metric_value(
        metrics, "hmix", "MAE", reference, 1.0, "observed_mean"
    )
    return [
        "## 3. Ranking two models",
        "",
        "Metric noise for a Sage-like model (SD of the score over repeat "
        "benchmark runs), the resolvable difference over all error profiles, and "
        "the SD of Sage's score over the choice of systems. Units as in the "
        "thresholds table of [report.md](report.md).",
        "",
        md_table(pd.DataFrame(rows, columns=columns)),
        "",
        "Noise also biases each score upwards, most for good models: a model whose "
        f"true ΔHmix MAE is {quarter_true:.2f} kJ/mol scores "
        f"{quarter_observed:.2f} on average at {reference / 1000:g} ns, against "
        f"{full_observed:.2f} "
        f"for a true {sage_mae:.2f}. Differences between good models are slightly "
        "compressed, but their order is kept.",
        "",
        "The two kinds of uncertainty answer different questions. The simulation "
        "noise says whether rerunning the benchmark would give the same ranking on "
        f"these {n_systems} systems; at {reference / 1000:g} ns it would, for any gap "
        "above the resolvable "
        "difference. The system spread says how far a model's score depends on "
        "this set of systems, which no run length can reduce.",
        "",
        "For a ranking, what matters is the system spread of the *difference* "
        "between two models, which depends on how alike their errors are. If "
        "model B's errors are model A's scaled by (1 + s), the gap is s·MAE_A and "
        f"its system SD is s·{hmix_systems:.2f}, i.e. "
        f"{hmix_systems / sage_mae:.2f}× the gap at any size, so the choice of "
        "systems never threatens that ranking. If their errors are unrelated, the "
        f"gap's system SD approaches sqrt(2)·{hmix_systems:.2f} = "
        f"{np.sqrt(2) * hmix_systems:.2f} kJ/mol, and gaps below about twice that "
        "are specific to this set of systems even though a rerun would reproduce "
        "them. Real model pairs (e.g. two sizes of one model family) are likely "
        "to sit between these cases.",
        "",
    ]


def run_length_section(tables: dict) -> list[str]:
    """
    Write the section on run length, cost and per-box runs for ``details.md``.

    Parameters
    ----------
    tables
        The tables named in :func:`write_reports`.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    metrics, allocation, n_atoms = (
        tables["metrics"],
        tables["allocation"],
        tables["n_atoms"],
    )
    reference = REFERENCE[1]
    cost = costs(n_atoms)
    rows = []
    for production in PRODUCTION_LENGTHS_PS:
        hmix = metric_value(metrics, "hmix", "MAE", production, 1.0, "resolvable")
        density = metric_value(metrics, "density", "MAE", production, 1.0, "resolvable")
        floor = metric_value(metrics, "hmix", "MAZE", production, 0.0, "observed_mean")
        rows.append(
            [
                f"0.2 + {production / 1000:g}",
                f"{cost[production] / cost[reference]:.0%}",
                f"{gpu_days(cost[production]):.0f}",
                f"{hmix:.3f}",
                f"{density:.4f}",
                f"{floor:.2f}",
            ]
        )
    columns = [
        "Protocol (ns)",
        "Cost",
        "MACE-OFF23(M) GPU-days",
        "Resolvable ΔHmix MAE (kJ/mol)",
        "Resolvable density MAE (g/mL)",
        "Perfect ΔHmix MAZE",
    ]
    equilibration = tables["equilibration"].set_index("observable")
    density_shifts = tables["shifts"].xs("density", level="observable")
    by_replicate = density_shifts.groupby(level="replicate").mean()
    longest = allocation.nlargest(4, "production_ps").box_id
    saved = 1 - allocation_cost(allocation) / cost[reference]
    return [
        "## 4. Run length and cost",
        "",
        md_table(pd.DataFrame(rows, columns=columns)),
        "",
        f"Cost is Σ atoms × simulated time over the {len(n_atoms)} boxes "
        f"({n_atoms.sum():,.0f} atoms). GPU-days are an order of magnitude only: "
        "MACE-OFF23(M) runs 1.1×10⁶ steps/day for ~600-atom liquid boxes on one A100 "
        "([arXiv:2312.15211](https://arxiv.org/abs/2312.15211), sec. III.3.2), "
        f"scaled linearly to our {n_atoms.min():,.0f} to {n_atoms.max():,.0f}-atom "
        "boxes. Faster and slower models scale the column, not the ratios.",
        "",
        "**Equilibration.** The mean over the 100 to 200 ps window minus the mean "
        f"over the {reference / 1000:g} ns production, over all {len(density_shifts)} "
        "Sage runs, is "
        f"{equilibration.loc['h', 'mean_shift']:+.3f} kJ/mol for h (p = "
        f"{equilibration.loc['h', 'p']:.1g}) and "
        f"{equilibration.loc['density', 'mean_shift']:+.4f} g/mL for density "
        f"(p = {equilibration.loc['density', 'p']:.0e}). Even for Sage, starting "
        "from Sage-built boxes, the window is shifted. By replicate, the density "
        "shift is "
        + ", ".join(
            f"{shift:+.4f} ({replicate})" for replicate, shift in by_replicate.items()
        )
        + " g/mL; replicate 0 is the distributed box that MLP runs start from. An "
        "MLP also has to relax from Sage's density to its own, which this test "
        "cannot see, so 0.2 ns stays. Table: [equilibration.csv](equilibration.csv).",
        "",
        "**Per-box run lengths.** For a fixed variance of ΔHmix MAE, mixture "
        "density MAE and pure density MAE (equally weighted, each relative to "
        f"uniform {reference / 1000:g} ns), the cheapest allocation runs each box "
        "for "
        "t ∝ sqrt(A_b / n_b), where A_b is the box's variance contribution and n_b "
        f"its atom count, with at least {MIN_ALLOCATION_PS / 1000:g} ns. It ranges "
        f"from {allocation.production_ps.min() / 1000:.1f} to "
        f"{allocation.production_ps.max() / 1000:.1f} ns "
        f"([allocation.csv](allocation.csv)); the longest runs go to "
        f"{', '.join(longest)}. It saves {saved:.0%} of the cost. That assumes "
        "σ ∝ 1/sqrt(t), which section 1 shows is optimistic for short runs, so the "
        "real saving is smaller. Not worth the extra complexity.",
        "",
    ]


def literature_section(tables: dict) -> list[str]:
    """
    Write the literature section of ``details.md``.

    Parameters
    ----------
    tables
        The tables named in :func:`write_reports`.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    metrics, sage_points = tables["metrics"], tables["sage_points"]
    reference = REFERENCE[1]
    sage_mae = metric_value(metrics, "hmix", "MAE", reference, 1.0, "true")
    # |ΔHmix| at the composition nearest x = 0.5 in each system.
    middle = sage_points.loc[
        (sage_points.x - 0.5).abs().groupby(sage_points.system_id).idxmin()
    ].hmix_exp.abs()
    density_protocol = resolvable_range(metrics, "density", "MAE", reference)[1]
    density_short = resolvable_range(
        metrics, "density", "MAE", PRODUCTION_LENGTHS_PS[0]
    )[1]
    literature = pd.DataFrame(LITERATURE)[["model", "property", "result", "source"]]
    return [
        "## 5. Literature: how far apart are current MLPs?",
        "",
        md_table(literature.rename(columns=str.capitalize)),
        "",
        "No study computing ΔHmix or V^E of molecular liquids with an MLP was "
        "found (searches of the web and of Consensus for MLP / neural network "
        "potential enthalpy of mixing, excess enthalpy and excess volume, October "
        "2026). The closest work concerns metallic alloys "
        "([arXiv:2406.17499](https://arxiv.org/abs/2406.17499)) or fits mixing "
        "data without MD.",
        "",
        "- **Density.** Published MLPs differ in pure-liquid density MAE by about "
        f"{LITERATURE_DENSITY_SPREAD[0]} to {LITERATURE_DENSITY_SPREAD[1]} g/mL. "
        f"The benchmark resolves {density_protocol:.4f} g/mL at "
        f"{reference / 1000:g} ns and "
        f"{density_short:.4f} g/mL at {PRODUCTION_LENGTHS_PS[0] / 1000:g} ns. "
        "Density alone could be ranked from much shorter runs, but the runs serve "
        "ΔHmix too.",
        "- **ΔHmix.** There is no MLP baseline. MLP ΔHvap errors are several "
        "kJ/mol (MACE-OFF23(M) offset about 8 kJ/mol; Picha et al. find large ΔHvap "
        "deviations for MACE-OFF23(M)). ΔHmix is a difference of such cohesive "
        "energies, so much of that error may cancel, but whatever remains is "
        f"compared with |ΔHmix| near x = 0.5 of {middle.min():.1f} to "
        f"{middle.max():.1f} kJ/mol. If MLPs differ in ΔHmix MAE by even a fifth "
        f"of Sage's MAE (about {0.2 * sage_mae:.2f} kJ/mol), "
        f"{reference / 1000:g} ns resolves it.",
        "",
    ]


def caveats_section(tables: dict) -> list[str]:
    """
    Write the caveats section of ``details.md``.

    Parameters
    ----------
    tables
        The tables named in :func:`write_reports`.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    finite = tables["finite"]
    hmix_gaps = finite[finite.prop == "hmix"].set_index("box_id")["diff"].abs()
    k_h = tables["noise"]
    k_h = k_h[
        (k_h.observable == "h")
        & (k_h.equilibration_ps == REFERENCE[0])
        & (k_h.production_ps == REFERENCE[1])
    ].k.iloc[0]
    return [
        "## Caveats",
        "",
        "- **Sage stands in for MLP noise.** MLP fluctuations of h and volume "
        "should be similar (the same liquids at the same T and P), but correlation "
        "times may differ, especially for MLPs that change conformer barriers. One "
        "pilot MLP run per box, analysed with `mixing.box_result`, would check "
        "this: compare its sem and n_eff with the Sage `runs.csv`.",
        f"- **Slow modes.** k = {k_h:.2f} (h, {REFERENCE[1] / 1000:g} ns) is "
        "pooled over boxes. Single "
        "boxes have larger run-to-run spreads (for Sage, NMP + "
        "1,2-dichloroethane; see the [Sage report](../sage_summary/report.md)), "
        "and an MLP may have slow modes of its own.",
        "- **Systematic simulation errors** (finite size, starting structures) are "
        "not noise and are not included. Against the published 1000-molecule Sage "
        f"runs, the median ΔHmix gap is {hmix_gaps.median():.3f} kJ/mol and the "
        f"largest {hmix_gaps.max():.2f} kJ/mol ({hmix_gaps.idxmax()}).",
        "- **Error profiles.** The model pairs are built from Sage's errors. Real "
        "MLPs may spread their errors differently; the closed-form SD shows the "
        "result depends mainly on the noise, not on the profile.",
        "",
    ]


def write_details(out: Path, tables: dict) -> None:
    """
    Write ``details.md``.

    Parameters
    ----------
    out
        Output directory.
    tables
        The tables named in :func:`write_reports`.
    """
    lines = [
        "# Power analysis: details",
        "",
        f"Generated {date.today()} by `power_analysis.py`. Summary: "
        "[report.md](report.md). Tables: `box_noise.csv`, `truncation.csv`, "
        "`points.csv`, `metrics.csv`, `bootstrap.csv`, `allocation.csv`, "
        "`equilibration.csv`, `thresholds.csv`, `literature.csv`.",
        "",
        *method_section(tables),
        *noise_section(tables),
        *point_section(tables),
        *ranking_section(tables),
        *run_length_section(tables),
        *literature_section(tables),
        *caveats_section(tables),
    ]
    (out / "details.md").write_text("\n".join(lines))


def write_reports(out: Path, tables: dict) -> None:
    """
    Write ``report.md`` and ``details.md``.

    Parameters
    ----------
    out
        Output directory.
    tables
        ``noise``, ``equilibration``, ``shifts`` (per-run early-window shifts),
        ``points``, ``metrics``, ``bootstrap``, ``allocation``, ``thresholds``,
        ``n_atoms``, ``sage_points`` (the Sage ``points.csv``) and ``finite``
        (the Sage ``finite_size.csv``).
    """
    write_summary(out, tables)
    write_details(out, tables)
