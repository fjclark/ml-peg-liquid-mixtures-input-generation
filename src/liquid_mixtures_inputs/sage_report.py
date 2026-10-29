"""
Write the Sage reference-run reports from the tables of ``summarise_sage.py``.

``report.md`` is the short summary; ``details.md`` holds the methods, the
statistics behind each claim and the NMP + 1,2-dichloroethane case; and
``inspection.csv`` lists every box or system worth a closer look.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from liquid_mixtures_inputs.constants import EQUILIBRATION_PS, PRODUCTION_PS
from liquid_mixtures_inputs.dce_torsion_scan import dce_summary
from liquid_mixtures_inputs.summarise_sage import (
    EARLY_SHIFT,
    EARLY_WINDOW_PS,
    FLAG_P,
    PUBLISHED_FF,
    PUBLISHED_FF_URL,
    PUBLISHED_X_TOLERANCE,
    SIGMA_SAGE_AT_FLOOR,
    SIZE_GAP,
    START_T,
)
from ml_peg.analysis.molecular_dynamics.binary_mixtures.mixing import (
    MIN_N_EFF,
    N_BLOCKS,
)

# What each inspection category means, in rough order of concern.
CATEGORY_NOTES = {
    "within/between mismatch": (
        "Replicates (independently built boxes) scatter more than the within-run "
        "standard errors allow: likely a slow mode (conformers, local structure) "
        f"that one {PRODUCTION_PS / 1000:g} ns run does not sample, so a single "
        "run's confidence interval "
        "is too narrow there. Most entries are p < 0.05 before Holm correction "
        "only, so they are places to look, not individual results."
    ),
    "low sampling": (
        f"A run with fewer than {MIN_N_EFF} independent samples (long correlation "
        "times, viscous or slowly reorganising liquids)."
    ),
    "drift or half-split": "A trend within a run, significant after Holm correction.",
    "start dependence": (
        "The distributed box (replicate 0) differs from the rebuilt replicates, "
        "so the starting structure may bias MLP runs."
    ),
    "early density shift": (
        f"The density of the first {EARLY_WINDOW_PS:.0f} ps of production "
        f"differs from the rest of the run by more than {EARLY_SHIFT} g/mL: still "
        f"relaxing after the {EQUILIBRATION_PS / 1000:g} ns discard, or a large "
        "fluctuation."
    ),
    "bonding changed": "Molecules changed during a classical run: a bug.",
    "microheterogeneity": (
        "Sage itself shows like-like clustering beyond random mixing; MLP "
        "demixing flags for this box should be compared against this Sage value."
    ),
    "published-Sage gap": (
        f"Our replicate mean differs from the published 1000-molecule run by "
        f"more than {SIZE_GAP['hmix']} kJ/mol or {SIZE_GAP['density']} g/mL."
    ),
    "force-field miss (not a run problem)": (
        "Sage is far from experiment. Expected for some systems, and useful for "
        "separating models."
    ),
}


START_UNITS = {"h_mean": ("h", "kJ/mol"), "density_mean": ("density", "g/mL")}
GAP_UNITS = {"hmix": ("ΔHmix", "kJ/mol", 2), "density": ("density", "g/mL", 4)}


def sage_scores(summary: pd.DataFrame, pure_summary: pd.DataFrame) -> pd.DataFrame:
    """
    Score Sage's replicate means against experiment with the benchmark metrics.

    Parameters
    ----------
    summary
        Per mixture point, with Sage means, experiment and σ_ref.
    pure_summary
        Per pure box, with Sage mean, experiment, σ_ref and ``scored``.

    Returns
    -------
    pd.DataFrame
        One row per property: number of points, MAE (with unit) and MAZE.
    """
    with_density = summary[summary.density_exp.notna()]
    scored_pure = pure_summary[pure_summary.scored.fillna(False).astype(bool)]
    sets = [
        ("ΔHmix", "kJ/mol", summary, "hmix_sage", "hmix_exp", "sigma_ref_hmix"),
        (
            "Density",
            "g/mL",
            with_density,
            "density_sage",
            "density_exp",
            "sigma_ref_density",
        ),
        ("V^E", "cm³/mol", with_density, "vex_sage", "vex_exp", "sigma_ref_vex"),
        (
            "Pure density",
            "g/mL",
            scored_pure,
            "density_sage",
            "density_exp",
            "sigma_ref_density",
        ),
    ]
    rows = []
    for name, unit, table, sage, experiment, sigma_ref in sets:
        error = (table[sage] - table[experiment]).abs()
        rows.append(
            {
                "property": name,
                "points": len(table),
                "mae": float(error.mean()),
                "unit": unit,
                "maze": float((error / table[sigma_ref]).mean()),
            }
        )
    return pd.DataFrame(rows)


def system_errors(summary: pd.DataFrame) -> pd.DataFrame:
    """
    Compute Sage's mean |ΔHmix error| per system and its wrong-sign count.

    Parameters
    ----------
    summary
        Per mixture point, with Sage means and experiment.

    Returns
    -------
    pd.DataFrame
        One row per system, largest error first.
    """
    error = (summary.hmix_sage - summary.hmix_exp).abs()
    wrong_sign = np.sign(summary.hmix_sage) != np.sign(summary.hmix_exp)
    table = pd.DataFrame(
        {"system_id": summary.system_id, "error": error, "wrong_sign": wrong_sign}
    )
    return (
        table.groupby("system_id")
        .agg(
            error=("error", "mean"),
            wrong_sign=("wrong_sign", "sum"),
            points=("error", "size"),
        )
        .sort_values("error", ascending=False)
        .reset_index()
    )


def join_words(items: list[str]) -> str:
    """
    Join items as ``"a"``, ``"a and b"`` or ``"a, b and c"``.

    Parameters
    ----------
    items
        Texts to join.

    Returns
    -------
    str
        The joined text, or ``"none"`` if there are no items.
    """
    if not items:
        return "none"
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


def miss_text(row) -> str:
    """
    Describe one system's ΔHmix error, e.g. ``water_pyridine (1.20 kJ/mol)``.

    Parameters
    ----------
    row
        Row of :func:`system_errors`.

    Returns
    -------
    str
        System, mean |error| and, if any, the wrong-sign count.
    """
    text = f"{row.system_id} ({row.error:.2f} kJ/mol"
    if row.wrong_sign:
        text += f", wrong sign at {row.wrong_sign:.0f} of {row.points} compositions"
    return text + ")"


def start_text(row) -> str:
    """
    Describe one start-dependent box, e.g. ``x (density -0.0019 g/mL, t = -4.8)``.

    Parameters
    ----------
    row
        Row of the start-dependence table.

    Returns
    -------
    str
        Box, observable, offset of replicate 0 and its t value.
    """
    name, unit = START_UNITS[row.observable]
    return f"{row.box_id} ({name} {row.offset:+.4f} {unit}, t = {row.t:+.1f})"


def gap_text(row) -> str:
    """
    Describe one published-Sage gap, e.g. ``x (ΔHmix -0.22 kJ/mol, z = -1.8)``.

    Parameters
    ----------
    row
        Row of the finite-size table.

    Returns
    -------
    str
        Box, property, our mean minus the published mean, and z.
    """
    name, unit, digits = GAP_UNITS[row.prop]
    return f"{row.box_id} ({name} {row.diff:+.{digits}f} {unit}, z = {row.z:+.1f})"


def inspection_table(items: pd.DataFrame) -> pd.DataFrame:
    """
    Add the plot paths to the inspection items, for ``inspection.csv``.

    Parameters
    ----------
    items
        Output of ``summarise_sage.inspection_items``.

    Returns
    -------
    pd.DataFrame
        Columns ``category``, ``target``, ``evidence`` and ``plots``.
    """
    plots = []
    for item in items.itertuples():
        if isinstance(item.box, str):
            plots.append(
                f"figures/replicates/{item.box}.png; figures/convergence/{item.box}.png"
            )
        else:
            system = item.target.rsplit("_x", 1)[0]
            plots.append(f"figures/systems/{system}.png")
    return items[["category", "target", "evidence"]].assign(plots=plots)


def holm_significant(within: pd.DataFrame) -> str:
    """
    List the within/between tests that survive Holm correction.

    Parameters
    ----------
    within
        Within/between table.

    Returns
    -------
    str
        For example ``"pure_1_2_dichloroethane_298 (density)"``, or ``"none"``.
    """
    flagged = within[within.flag != ""]
    names = [f"{row.box_id} ({row.observable})" for row in flagged.itertuples()]
    return ", ".join(names) or "none"


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
    runs, items, finite = tables["runs"], tables["items"], tables["finite"]
    pooled = tables["pooled"].set_index(["level", "observable"])
    scores = sage_scores(tables["summary"], tables["pure_summary"])
    misses = system_errors(tables["summary"])
    dce = dce_summary(out)
    h_ratio = pooled.loc[("box", "h"), "pooled_ratio"]
    density_ratio = pooled.loc[("box", "density"), "pooled_ratio"]
    finite_gap = finite["diff"].abs() > finite.prop.map(SIZE_GAP)
    gaps = finite[finite_gap]
    start = tables["start"][tables["start"].t.abs() > START_T]
    counts = items.groupby("category", sort=False).target.agg(["size", "nunique"])
    worst, second = misses.iloc[0], misses.iloc[1]

    lines = [
        "# Sage (openff-2.3.0) reference runs",
        "",
        f"Generated {date.today()} by `summarise_sage.py`. Methods and statistics: "
        "[details.md](details.md). Every flagged box: "
        "[inspection.csv](inspection.csv).",
        "",
        f"Sage was run through the benchmark protocol ({EQUILIBRATION_PS / 1000:g} "
        f"ns NPT equilibration, {PRODUCTION_PS / 1000:g} ns production) on all "
        f"{runs.box_id.nunique()} boxes, {runs.replicate.nunique()} independent "
        f"replicates each ({len(runs)} runs). The replicate spread gives σ_Sage, "
        "the MAZE denominator, and Sage's own score sets the 'bad' thresholds "
        "([power analysis](../power_analysis/report.md)).",
        "",
        "## Sage against experiment",
        "",
        "| Property | Points | MAE | MAZE |",
        "|---|---|---|---|",
        *[
            f"| {row.property} | {row.points} | {row.mae:.3g} {row.unit} | "
            f"{row.maze:.3g} |"
            for row in scores.itertuples()
        ],
        "",
        f"The largest ΔHmix errors are for {miss_text(worst)} and "
        f"{miss_text(second)}. These are force-field misses, not run problems, "
        "and they are the systems most likely to separate models.",
        "",
        "![Sage vs experiment](figures/sage_vs_exp.png)",
        "",
        "## Run-to-run noise",
        "",
        "Pooled over all boxes, the SD of the replicate means is "
        f"{h_ratio:.2f}× the within-run standard error for h and "
        f"{density_ratio:.2f}× for density (p = "
        f"{pooled.loc[('box', 'h'), 'p']:.0e} and "
        f"{pooled.loc[('box', 'density'), 'p']:.0e}). A single "
        f"{PRODUCTION_PS / 1000:g} ns run's "
        f"standard error is therefore about {(h_ratio - 1):.0%} too small for h "
        f"and {(density_ratio - 1):.0%} for density, because independently "
        "built boxes sample slow modes (conformers, local structure) that one run "
        "does not. This is why σ_Sage comes from the replicate spread rather than "
        "within-run errors. Per-box tests that survive Holm correction: "
        f"{holm_significant(tables['within'])}.",
        "",
        "![within vs between](figures/within_between.png)",
        "",
        "## Boxes to inspect",
        "",
        f"{len(items)} items over {items.target.nunique()} boxes or systems, all "
        "in [inspection.csv](inspection.csv) with links to their plots. Most are "
        "places to look rather than problems:",
        "",
        "| Category | Items | Boxes or systems |",
        "|---|---|---|",
        *[
            f"| {category} | {counts.loc[category, 'size']} | "
            f"{counts.loc[category, 'nunique']} |"
            for category in CATEGORY_NOTES
            if category in counts.index
        ],
        "",
        "Those that matter for the benchmark:",
        "",
        "- **NMP + 1,2-dichloroethane** has slow Sage torsion sampling (see "
        "below). It is kept, with σ_Sage set to the median floor.",
        "- **Start dependence**: "
        + join_words([start_text(row) for row in start.itertuples()])
        + ". The distributed box (replicate 0) differs from the rebuilt "
        "replicates by more than their spread allows, so an MLP started from it "
        "may be slightly biased.",
        f"- **Published-Sage gap**: {len(gaps)} of {len(finite)} matched points "
        "(next section).",
        "",
        "## Box size",
        "",
        f"Our boxes (≥ 100 molecules) were compared with the published "
        f"1000-molecule openff-2.3.0 runs ([`{PUBLISHED_FF}`]({PUBLISHED_FF_URL})) "
        f"at {len(finite)} matched points. Points with |z| > 3: "
        f"{int((finite.z.abs() > 3).sum())}. Points that differ by more than "
        f"{SIZE_GAP['hmix']} kJ/mol or "
        f"{SIZE_GAP['density']} g/mL: {len(gaps)}, "
        + join_words([gap_text(row) for row in gaps.itertuples()])
        + ". The median ΔHmix gap is "
        f"{finite[finite.prop == 'hmix']['diff'].abs().median():.3f} kJ/mol and "
        "the median density gap "
        f"{finite[finite.prop == 'density']['diff'].abs().median():.4f} g/mL, well "
        "below Sage's errors against experiment, so box size does not bias the "
        "reference overall; the points above are in "
        "[inspection.csv](inspection.csv).",
        "",
        "![finite size](figures/finite_size.png)",
        "",
        "## NMP + 1,2-dichloroethane",
        "",
        "In Sage, only "
        f"{dce['changed_per_frame']:.0%} of 1,2-dichloroethane molecules switch "
        "between trans and gauche from one 100 ps frame to the next "
        f"({dce['uncorrelated_limit']:.0%} if uncorrelated), a relaxation time of "
        f"about {dce['tau_sage_ps'] / 1000:.1f} ns. The cause is Sage's barrier, "
        f"{dce['barrier_sage'] - dce['barrier_reference']:.0f} kJ/mol above "
        "MACE-OFF23(M), which makes isomerisation about "
        f"{dce['slowdown']:.0f}× too slow. A model with a reference-like barrier "
        f"should relax in about {dce['tau_reference_ps']:.0f} ps, well inside the "
        "equilibration, so the system is a fair test for MLPs and was kept. "
        "Details and the torsion scans: "
        "[details.md](details.md#why-nmp--12-dichloroethane-was-kept).",
        "",
    ]
    (out / "report.md").write_text("\n".join(lines))


def within_between_section(tables: dict) -> list[str]:
    """
    Write the within- vs between-run section of ``details.md``.

    Parameters
    ----------
    tables
        The tables named in :func:`write_reports`.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    pooled, within, items = tables["pooled"], tables["within"], tables["items"]
    passed = pooled[pooled.p * len(pooled) < FLAG_P]
    bonferroni = (
        ", ".join(f"{row.level} {row.observable}" for row in passed.itertuples())
        or "none"
    )
    # The three boxes with the largest replicate SD / within-run sem.
    mismatched = items[items.category == "within/between mismatch"].box.unique()
    box_level = within[(within.level == "box") & within.box_id.isin(mismatched)]
    largest = box_level.groupby("box_id").ratio.max().nlargest(3).index
    return [
        "## Within-run vs between-run uncertainty",
        "",
        "Within-run standard errors come from the autocorrelation (Chodera) "
        "estimator on every logged sample. If they were right, the SD of the "
        "replicate means would equal them on average (ratio 1). Each box and "
        "mixture point gets a chi-squared test of the replicate means against "
        "their within-run errors; the pooled ratio combines all of them.",
        "",
        "| Level | Observable | Tests | Median ratio | Pooled ratio | p | "
        f"p_high < {FLAG_P} | p_low < {FLAG_P} |",
        "|---|---|---|---|---|---|---|---|",
        *[
            f"| {row.level} | {row.observable} | {row.n_tests} | "
            f"{row.median_ratio:.2f} | {row.pooled_ratio:.2f} | {row.p:.1g} | "
            f"{row.n_sig_high} | {row.n_sig_low} |"
            for row in pooled.itertuples()
        ],
        "",
        f"**Multiple testing.** The pooled ratios are {len(pooled)} tests; with a "
        f"Bonferroni correction over them, {bonferroni} stay significant at "
        f"p < {FLAG_P}. The per-box and per-point tests are two-sided and "
        f"Holm-corrected together ({len(within)} tests). The p_high counts above "
        f"are unadjusted: about {FLAG_P * pooled.n_tests.iloc[0]:.1f} per row are "
        "expected by chance. Surviving correction: "
        f"{holm_significant(within)}. The evidence that within-run errors are too "
        "small therefore comes from the pooled test, not from any one box. MLP "
        "runs are single runs, so their confidence intervals will be too narrow "
        "by roughly the pooled ratio.",
        "",
        "Largest mismatches:",
        "",
        *[f"![{box_id}](figures/replicates/{box_id}.png)" for box_id in largest],
        "",
    ]


def maze_section(tables: dict) -> list[str]:
    """
    Write the MAZE denominator section of ``details.md``.

    Parameters
    ----------
    tables
        The tables named in :func:`write_reports`.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    summary = tables["summary"]
    return [
        "## Replicate spread (MAZE denominator)",
        "",
        "MAZE divides each error by a fixed σ_ref = sqrt(σ_Sage² + σ_exp²), the "
        "same for every model, and averages |error| / σ_ref over the points. "
        "σ_Sage is one value per mixture: the SD of the replicate means, pooled "
        "over its compositions (RMS) and floored at the median over mixtures. "
        "Pure densities use the same rule, with one box per liquid. The floor is "
        "there because model errors are mostly systematic and much larger than "
        "σ_ref, so σ_ref acts as a weight: without it, the mixtures where Sage "
        "happens to sample most quietly would dominate MAZE, although an MLP's "
        "noise there need not be as small. With the floor, MAZE stays close to "
        "the MAE divided by a single σ. "
        f"{', '.join(sorted(SIGMA_SAGE_AT_FLOOR))} and pure 1,2-dichloroethane "
        "are set to the floor, as their spread comes from a Sage artefact (last "
        "section).",
        "",
        "Median σ_Sage over mixtures: "
        f"{summary.sigma_sage_hmix.median():.3f} kJ/mol (ΔHmix), "
        f"{summary.sigma_sage_density.median():.4f} g/mL (density), "
        f"{summary.sigma_sage_vex.median():.3f} cm³/mol (V^E). The values are "
        "written to `systems.csv` and `pure.csv` (`sigma_sage_*`, `sigma_ref_*`); "
        "per-point tables are in [points.csv](points.csv) and "
        "[pure_points.csv](pure_points.csv).",
        "",
        "![replicate spread](figures/replicate_spread.png)",
        "",
    ]


def finite_size_section(tables: dict) -> list[str]:
    """
    Write the finite-size section of ``details.md``.

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
    return [
        "## Finite-size check",
        "",
        f"Our replicate means are compared with the published runs "
        f"[`{PUBLISHED_FF}`]({PUBLISHED_FF_URL}) (openff-2.3.0rc2, 3 replicates of "
        "1000 molecules), whose bonded, vdW and library-charge parameters are "
        "identical to the released openff-2.3.0. Both use NAGL charges "
        f"(published: model 0.1.0-rc.3, [line 487]({PUBLISHED_FF_URL}#L487); ours: "
        "1.0.0), so any gap comes from box size, protocol or the NAGL version. "
        f"{len(finite)} points match within |Δx| ≤ {PUBLISHED_X_TOLERANCE}. The "
        "nearest published value is shifted to our composition along a curve "
        "fitted to the published points (two-term Redlich-Kister for ΔHmix, "
        "quadratic for density), as a composition offset alone can shift a "
        "density by more than the gaps tested. z is the gap over the combined "
        "standard error of the two means: it tests whether the gap exceeds "
        "replicate noise, not what causes it. "
        f"|z| > 3: {int((finite.z.abs() > 3).sum())}. "
        "All points: [finite_size.csv](finite_size.csv).",
        "",
    ]


def convergence_section(tables: dict) -> list[str]:
    """
    Write the convergence section of ``details.md``.

    Parameters
    ----------
    tables
        The tables named in :func:`write_reports`.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    runs, flags, start = tables["runs"], tables["flags"], tables["start"]
    shift = runs.density_first_block_shift
    block_ps = PRODUCTION_PS / N_BLOCKS
    n_tests = 4 * len(runs)
    expected_start = len(start) * 2 * stats.t.sf(START_T, runs.replicate.max() - 1)
    by_replicate = runs.groupby("replicate").density_first_block_shift
    return [
        "## Convergence",
        "",
        f"h is the enthalpy per molecule, (E_pot + pV) / N in kJ/mol. The block "
        f"tests split production into {N_BLOCKS} blocks of {block_ps:.0f} ps.",
        "",
        f"- **Independent samples.** Median n_eff {runs.h_n_eff.median():.0f} "
        f"(correlation time {runs.h_tau_ps.median():.1f} ps), minimum "
        f"{runs.h_n_eff.min():.0f}. {int((runs.h_n_eff < MIN_N_EFF).sum())} runs "
        f"fall below ml-peg's minimum of {MIN_N_EFF}.",
        "- **Cross-check of the standard error.** The scatter of the block means "
        f"(SD / √{N_BLOCKS}) over the autocorrelation sem has a median of "
        f"{runs.h_block_ratio.median():.2f} (h) and "
        f"{runs.density_block_ratio.median():.2f} (density). Single-run ratios "
        f"scatter by about ±{1 / np.sqrt(2 * (N_BLOCKS - 1)):.0%}, so only the "
        "median is meaningful.",
        f"- **Equilibration.** The first {EARLY_WINDOW_PS:.0f} ps of production "
        f"({EQUILIBRATION_PS:.0f} to {EQUILIBRATION_PS + EARLY_WINDOW_PS:.0f} ps) "
        "minus the "
        f"rest of the run: median |shift| {shift.abs().median():.4f} g/mL, max "
        f"{shift.abs().max():.4f} g/mL. {int((shift.abs() > EARLY_SHIFT).sum())} "
        f"runs exceed {EARLY_SHIFT} g/mL. Mean shift by replicate (g/mL, with "
        "the share of boxes above 0): "
        + ", ".join(
            f"{replicate}: {shift.mean():+.4f} ({(shift > 0).mean():.0%})"
            for replicate, shift in by_replicate
        )
        + ". A shift shared by most boxes of one replicate points at that "
        "replicate's starting structures rather than at a system.",
        "- **Trends within a run.** Drift (slope of the block means) and "
        f"half-split (first {N_BLOCKS // 2} blocks against the last "
        f"{N_BLOCKS // 2}) tests on h and density: {n_tests} tests, "
        f"Holm-corrected. Flagged: {len(flags)}.",
        *[f"  - {text}" for _, text in flags],
        f"- **Start dependence.** Replicate 0 (the distributed box) against a "
        f"prediction interval from the rebuilt replicates; |t| > {START_T} is "
        f"reported. About {expected_start:.1f} of the {len(start)} tests are "
        "expected above this by chance. All tests: "
        "[start_dependence.csv](start_dependence.csv).",
        "",
        "Every box has a convergence plot in `figures/convergence/` and every "
        "system a plot in `figures/systems/`.",
        "",
    ]


def inspection_section(tables: dict) -> list[str]:
    """
    Write the inspection-category section of ``details.md``.

    Parameters
    ----------
    tables
        The tables named in :func:`write_reports`.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    items = tables["items"]
    lines = [
        "## Inspection categories",
        "",
        "Every item is in [inspection.csv](inspection.csv), with its evidence "
        "and plots. The categories, in rough order of concern:",
        "",
    ]
    for category, note in CATEGORY_NOTES.items():
        count = int((items.category == category).sum())
        if count:
            lines.append(f"- **{category}** ({count}). {note}")
    return [*lines, ""]


def dce_section(out: Path) -> list[str]:
    """
    Write the NMP + 1,2-dichloroethane section of ``details.md``.

    Parameters
    ----------
    out
        Output directory with ``dce_flips.csv`` and ``dce_torsion_scan.csv``.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    dce = dce_summary(out)
    return [
        "## Why NMP + 1,2-dichloroethane was kept",
        "",
        "In Sage the Cl-C-C-Cl torsion of 1,2-dichloroethane (DCE) relaxes "
        "slowly. In every DCE box (the three mixtures and pure DCE), only about "
        f"{dce['changed_per_frame']:.0%} of DCE molecules switch between trans "
        "and gauche from one 100 ps frame to the next "
        f"({dce['uncorrelated_limit']:.0%} if uncorrelated), a relaxation time of "
        f"about {dce['tau_sage_ps'] / 1000:.1f} ns ([dce_flips.csv](dce_flips.csv)). "
        f"With {dce['fewest_molecules']:.0f} DCE molecules in the DCE-poor box, a "
        f"{PRODUCTION_PS / 1000:g} ns run samples this poorly, so the replicates "
        "scatter more than their "
        "within-run errors. The cause is Sage's barrier. In relaxed gas-phase "
        "scans (below; energies relative to trans, "
        "[dce_torsion_scan.csv](dce_torsion_scan.csv)), the gauche/trans gap "
        f"agrees with MACE-OFF23(M) ({dce['gauche_sage']:.1f} vs "
        f"{dce['gauche_reference']:.1f} kJ/mol), but the trans to gauche barrier "
        f"is {dce['barrier_sage']:.1f} vs {dce['barrier_reference']:.1f} kJ/mol. "
        f"At 298 K that makes Sage's isomerisation about {dce['slowdown']:.0f}× "
        "too slow (Arrhenius, equal prefactors). A model with a reference-like "
        f"barrier should relax in about {dce['tau_reference_ps']:.0f} ps, inside "
        f"the {EQUILIBRATION_PS:.0f} ps equilibration, and sample many flips per "
        "run, so the system is a fair test for MLPs. σ_Sage of the mixture and of "
        "pure DCE is set to the median floor so that Sage's slow sampling does "
        "not loosen the MAZE denominator.",
        "",
        "![DCE torsion scan](figures/dce_torsion_scan.png)",
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
    runs, summary = tables["runs"], tables["summary"]
    lines = [
        "# Sage reference runs: details",
        "",
        f"Generated {date.today()} by `summarise_sage.py`. Summary: "
        "[report.md](report.md).",
        "",
        f"{len(runs)} runs ({runs.box_id.nunique()} boxes, "
        f"{runs.replicate.nunique()} replicates), each {EQUILIBRATION_PS:.0f} ps "
        f"of equilibration and {PRODUCTION_PS:.0f} ps of production. "
        f"{int(summary.complete.sum())} of {len(summary)} mixture points have all "
        "replicates. Per-run statistics: [runs.csv](runs.csv).",
        "",
        *within_between_section(tables),
        *maze_section(tables),
        *finite_size_section(tables),
        *convergence_section(tables),
        *inspection_section(tables),
        *dce_section(out),
    ]
    (out / "details.md").write_text("\n".join(lines))


def write_reports(out: Path, tables: dict) -> None:
    """
    Write ``report.md``, ``details.md`` and ``inspection.csv``.

    Parameters
    ----------
    out
        Output directory.
    tables
        ``runs``, ``summary`` (mixture points), ``pure_summary``, ``finite``,
        ``within`` (within/between tests), ``pooled``, ``start``, ``flags``
        (Holm-flagged run tests as (box_id, text)) and ``items`` (inspection
        items).
    """
    inspection_table(tables["items"]).to_csv(out / "inspection.csv", index=False)
    write_summary(out, tables)
    write_details(out, tables)
