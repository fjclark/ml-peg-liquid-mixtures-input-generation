"""Write the selection report (``report.md`` plus PNG figures)."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from datetime import date
import math
from pathlib import Path
import statistics

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import Draw  # noqa: E402

from liquid_mixtures_inputs import select_systems  # noqa: E402
from liquid_mixtures_inputs.thermoml import canonical_smiles  # noqa: E402

AMIDE_SMILES = {canonical_smiles(s) for s in ("CNC=O", "CC(=O)NC", "CC(=O)N", "NC=O")}
FORCE_FIELD_DIR = (
    f"{select_systems.RAW}/ash-sage-rc2/{select_systems.ASH_COMMIT}/"
    "04_benchmark/forcefields"
)
RC2_OFFXML = f"{FORCE_FIELD_DIR}/sage-2.3.0rc2-v3-k100.offxml"
FB_FIT_OFFXML = f"{FORCE_FIELD_DIR}/{select_systems.SAGE_RC2}.offxml"
# The rc2 training-set row with a placeholder 25 % ΔHmix uncertainty.
PLACEHOLDER_SIGMA_LINE = 769


def format_number(value, digits: int = 3) -> str:
    """
    Format a number for a Markdown table.

    Parameters
    ----------
    value
        Number or ``None``.
    digits
        Decimal places.

    Returns
    -------
    str
        Formatted value, or an en dash for missing values.
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "–"
    return f"{value:.{digits}f}"


def github_url(raw_url: str) -> str:
    """
    Convert a raw GitHub file URL into a line-linkable GitHub page URL.

    Parameters
    ----------
    raw_url
        ``raw.githubusercontent.com`` URL pinned to a commit.

    Returns
    -------
    str
        The ``github.com/.../blob/...?plain=1`` URL for the same file.
    """
    repository, path = raw_url.removeprefix(select_systems.RAW + "/").split("/", 1)
    return f"https://github.com/openforcefield/{repository}/blob/{path}?plain=1"


def display_name(system: dict) -> str:
    """
    Return the display name of a system.

    Parameters
    ----------
    system
        System dictionary from :func:`select_systems.build_system`.

    Returns
    -------
    str
        ``"A + B"``.
    """
    return " + ".join(system["spec"]["names"])


def sorted_pair(a: str, b: str) -> tuple[str, str]:
    """
    Return the canonical, sorted SMILES pair used as a system key.

    Parameters
    ----------
    a, b
        SMILES of the two components.

    Returns
    -------
    tuple[str, str]
        The sorted canonical SMILES.
    """
    return tuple(sorted((canonical_smiles(a), canonical_smiles(b))))


def plot_repeat_sd(candidates: list[dict], path: Path) -> None:
    """
    Plot histograms of published Sage replicate SDs over all candidates.

    Parameters
    ----------
    candidates
        Candidate rows from :func:`select_systems.screen_candidates`.
    path
        Output PNG path.
    """
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
    for ax, (property_name, unit) in zip(
        axes, (("hmix", "kJ/mol"), ("density", "g/mL")), strict=True
    ):
        column = f"sage_rep_sd_{property_name}"
        values = [c[column] for c in candidates if c[column]]
        ax.hist(values, bins=40, color="#4a7ab5")
        ax.axvline(
            select_systems.MAX_REPEAT_SD[property_name],
            color="#c0392b",
            ls="--",
            label="limit",
        )
        ax.set_xlabel(f"max Sage 3-replicate SD ({unit})")
        ax.set_ylabel("systems")
        ax.set_title("ΔHmix" if property_name == "hmix" else "density")
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_system(system: dict, path: Path) -> None:
    """
    Plot experimental ΔHmix(x) and V^E(x) with fits and box compositions.

    Parameters
    ----------
    system
        System dictionary from :func:`select_systems.build_system`.
    path
        Output PNG path.
    """
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    grid = np.linspace(0.0, 1.0, 201)
    for ax, key, fit_key, label in (
        (axes[0], "hmix", "hmix_fit", "ΔHmix (kJ/mol)"),
        (axes[1], "vex", "vex_fit", "V^E (cm³/mol)"),
    ):
        by_doi = defaultdict(list)
        for point in system[key]:
            by_doi[point.doi].append(point)
        # The six papers with the most points.
        for doi, points in sorted(by_doi.items(), key=lambda item: -len(item[1]))[:6]:
            ax.errorbar(
                [p.x for p in points],
                [p.value for p in points],
                yerr=[p.sigma or 0.0 for p in points],
                fmt="o",
                ms=3,
                label=doi,
            )
        if system[fit_key]:
            ax.plot(grid, [system[fit_key](x) for x in grid], color="k", lw=1)
        box_x = [row["x_box"] for row in system["rows"]]
        box_values = [row[f"{key}_exp"] for row in system["rows"]]
        if all(value is not None for value in box_values):
            ax.plot(
                box_x, box_values, "s", mfc="none", mec="#c0392b", ms=9, label="box x"
            )
        ax.axhline(0, color="0.7", lw=0.5)
        ax.set_xlim(0, 1)
        ax.set_xlabel(f"x ({system['spec']['names'][0]})")
        ax.set_ylabel(label)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=6, frameon=False)
    fig.suptitle(
        f"{display_name(system)} at {system['temperature']:.2f} K", fontsize=10
    )
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def hard_filters_section() -> list[str]:
    """
    Describe the hard filters, using the thresholds from ``select_systems``.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    t_low, t_high = select_systems.SCREEN_T_RANGE
    max_relative = 100 * select_systems.MAX_RELATIVE_SIGMA_HMIX
    rows = [
        (
            "(entry) ≥ 3 ΔHmix and ≥ 3 density points",
            f"At {t_low}–{t_high} K, counted per distinct composition.",
            "Three compositions per system need a measured ΔHmix and density.",
        ),
        (
            "composition mismatch",
            "Every ΔHmix composition has a density point within "
            f"\\|Δx\\| ≤ {select_systems.MAX_COMPOSITION_MISMATCH}.",
            "Both properties are scored on the same box.",
        ),
        (
            "missing experimental uncertainty",
            "Every ΔHmix and density point reports an uncertainty.",
            "Needed for the MAZE denominator and to judge agreement.",
        ),
        (
            f"relative ΔHmix uncertainty > {max_relative:.0f} %",
            f"Largest relative ΔHmix uncertainty ≤ {max_relative:.0f} %.",
            "Relative experimental ΔHmix uncertainties among systems passing the "
            "other filters are 0.5–3.5 %. Larger values are flat 25 % uncertainties "
            "assigned during dataset curation (e.g. the rc2 entries for "
            "10.1016/j.jct.2013.09.015) or come from H^E near zero; either way they "
            "would inflate the MAZE denominator.",
        ),
        (
            "missing pure density with uncertainty",
            "Both pure liquids have a density with an uncertainty in the window.",
            "Needed for V^E and for scoring new pure densities.",
        ),
        (
            "elements",
            "Only C, H, N, O and Cl.",
            "Chosen scope: common organic and biomolecular elements.",
        ),
        (
            "too many heavy atoms",
            f"≤ {select_systems.MAX_HEAVY_ATOMS} heavy atoms per molecule.",
            "Keeps boxes small at MLIP cost while holding enough molecules.",
        ),
        (
            "viscous",
            "No diol, polyol, polyamine or bis-hydroxyethyl amine motifs "
            "(SMARTS in `VISCOUS_SMARTS`).",
            "Slow dynamics would not converge in 1 ns of production.",
        ),
        (
            "two-phase (linear H^E)",
            "ΔHmix is not consistent with a straight line in x (three points: "
            "consecutive slopes differ by > 2σ; more points: reduced χ² of a linear "
            "fit > 2). σ has a 2 % relative floor. Applied only when the largest "
            f"measured \\|H^E\\| ≥ {select_systems.MIN_SIGNAL} kJ/mol.",
            "A two-phase region gives H^E = a + b·x (lever rule). An ideal or "
            "near-ideal mixture is also a straight line (H^E ≈ 0), so the test is "
            "skipped for small H^E; the small-signal filter handles those.",
        ),
        (
            "two-phase (H^E > RT/2)",
            f"ΔHmix(x≈0.5) ≤ RT/2 = {select_systems.UCST_LIMIT:.2f} kJ/mol.",
            "Regular-solution critical value: H^E = βx(1 − x) splits into two "
            "phases below T_c = β/(2R) ([discussion]"
            "(https://scipython.com/blog/the-miscibility-of-a-regular-solution/)), "
            "and H^E(0.5) = β/4 = RT/2 at T_c. A conservative screen for an "
            "upper critical solution temperature above 298 K.",
        ),
        (
            "Sage repeat SD (hmix / density)",
            "Largest published Sage 3-replicate SD over openff-2.0.0 to 2.3.0rc1 "
            f"≤ {select_systems.MAX_REPEAT_SD['hmix']} kJ/mol (ΔHmix) and "
            f"≤ {select_systems.MAX_REPEAT_SD['density']} g/mL (density).",
            "A large spread between replicates suggests slow equilibration.",
        ),
        (
            "no Sage repeat SD (hmix / density)",
            "A published Sage SD exists for the property.",
            "Without it the convergence check above cannot be made.",
        ),
        (
            "small signal",
            f"\\|ΔHmix(x≈0.5)\\| ≥ {select_systems.MIN_SIGNAL} kJ/mol.",
            "Our estimate of the simulation σ is 0.15–0.25 kJ/mol (see "
            "`DECISIONS.md`), so smaller values might discriminate between "
            "models less well.",
        ),
    ]
    return [
        "## Hard filters",
        "",
        "A system must pass every filter below to be a candidate. The filter names "
        "match the labels in the funnel table and in the `failures` column of "
        "`candidates.csv`. A few systems are admitted as listed exceptions; the "
        "Selected systems table gives the filters each one fails and why it was kept.",
        "",
        "| Filter | Rule | Why |",
        "|---|---|---|",
        *[f"| {name} | {rule} | {why} |" for name, rule, why in rows],
        "",
        "Systems are also checked for literature evidence of miscibility at the run "
        "temperature (see Per-system details); that check is manual.",
        "",
    ]


def data_sources_section() -> list[str]:
    """
    List the data sources with pinned links.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    sources = select_systems.SOURCES
    return [
        "## Data sources",
        "",
        "- OpenFF Sage 2.0 training set "
        f"([openff-sage@{select_systems.SAGE_COMMIT[:7]}]"
        f"({github_url(sources['sage-2.0-train'])})): "
        "no experimental uncertainties.",
        "- Sage 2.3 refit "
        f"[training]({github_url(sources['rc2-train'])}) and "
        f"[validation]({github_url(sources['rc2-validation'])}) sets "
        f"(ash-sage-rc2@{select_systems.ASH_COMMIT[:7]}): "
        "experimental uncertainties for every point.",
        "- Published Sage benchmark results (ash-sage-rc2 "
        f"[`04_benchmark`]({github_url(select_systems.SAGE_BENCHMARKS['training-hmix'])})"
        "): mean "
        "and SD over 3 independent replicates for openff-2.0.0, 2.1.0, 2.2.1, "
        "2.3.0rc1 and 2.3.0rc2 (labelled `fb-fit-v3-single-mean-k100` in the "
        "benchmark files; "
        f"[`sage-2.3.0rc2-v3-k100.offxml`]({github_url(RC2_OFFXML)}) is identical to "
        f"[`fb-fit-v3-single-mean-k100.offxml`]({github_url(FB_FIT_OFFXML)})). "
        "2.3.0rc2 is shown in the Sage baseline table only and is not used by "
        "the repeat-SD filter.",
        "- Full isotherms and pure densities from the NIST TRC ThermoML archive, "
        "fetched by DOI. ThermoML stores either a standard uncertainty (one "
        "standard deviation) or an expanded uncertainty (a 95 % interval, "
        "coverage factor k = 2). Standard uncertainties are used as given; "
        "expanded ones are divided by 2. Every σ in this report is therefore a "
        "standard uncertainty. The ash-sage-rc2 sets use the same convention: "
        "their σ equals half the ThermoML expanded value.",
        "- V^E is derived from mixture densities and the pure densities of the "
        "same paper. ThermoML gives one combined density σ, often a flat "
        "~0.1 % assessed by the TRC data-file compiler, which is mostly "
        "calibration and purity error shared by every density in the paper. "
        "It is treated as fully correlated within a paper, so it largely "
        "cancels in V^E (the larger of a common scale error and a common "
        "offset is kept). The random part is each paper's RMS residual about "
        "the Redlich-Kister fit, added in quadrature. A pure density from "
        "another source is uncorrelated and adds in quadrature. Mixture "
        "density σ is used unchanged, since there the shared error counts.",
        "- ml-peg's liquid_densities inputs, to avoid re-scoring pure densities.",
        "",
    ]


def largest_percent(candidate: dict) -> float | None:
    """
    Return a candidate's largest relative ΔHmix uncertainty in percent.

    Parameters
    ----------
    candidate
        Candidate row.

    Returns
    -------
    float | None
        The percentage, or ``None`` if no uncertainty was reported.
    """
    relative = candidate["hmix_max_rel_sigma"]
    return relative and 100 * relative


def funnel_section(candidates: list[dict], points: dict) -> list[str]:
    """
    Count the systems at each screening step and per failed filter.

    Parameters
    ----------
    candidates
        Candidate rows.
    points
        Sage data-set points.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    n_binary = sum(1 for key in points if key[1])
    failures = Counter(
        failure
        for candidate in candidates
        for failure in candidate["failures"].split("; ")
        if failure
    )
    linear_flags = ", ".join(
        f"`{c['a']}` + `{c['b']}` ({format_number(largest_percent(c), 0)} %)"
        for c in candidates
        if c["linear_two_phase"]
    )
    placeholder_url = (
        f"{github_url(select_systems.SOURCES['rc2-train'])}#L{PLACEHOLDER_SIGMA_LINE}"
    )
    return [
        "## Screening funnel",
        "",
        "- Distinct binary systems across the Sage 2.0 training set and the "
        f"Sage 2.3 rc2 training and validation sets: {n_binary}.",
        "- With ≥ 3 ΔHmix and ≥ 3 density points at "
        f"{select_systems.SCREEN_T_RANGE[0]}–"
        f"{select_systems.SCREEN_T_RANGE[1]} K: {len(candidates)}.",
        f"- Passing every hard filter: {sum(c['passes'] for c in candidates)}.",
        "",
        "Failures by filter (a system can fail several):",
        "",
        "| Filter | Systems failing |",
        "|---|---|",
        *[f"| {name} | {count} |" for name, count in failures.most_common()],
        "",
        "The H^E-linearity test flags systems whose enthalpies are consistent with a "
        "straight line in x (the lever-rule signature of a two-phase region, H^E = a + "
        f"b·x) and whose largest |H^E| is at least {select_systems.MIN_SIGNAL} kJ/mol. "
        "Systems it "
        "flags in the pool, with the largest relative ΔHmix uncertainty: "
        + linear_flags
        + ". Large relative uncertainties can hide real curvature, so a flag on data "
        "with errors of ~25 % (placeholder values in the rc2 sets, e.g. "
        f"[σ = 0.1695 for ΔHmix = -0.678]({placeholder_url})) "
        "is inconclusive "
        "rather than evidence of two phases. Such systems are also rejected by the "
        "relative-uncertainty filter.",
        "",
        "![Sage repeat SD](figures/sage_repeat_sd.png)",
        "",
    ]


def read_fast_set(path: Path) -> dict[tuple[str, str], list]:
    """
    Read the recommended fast-converging list (Sage data-set format).

    Parameters
    ----------
    path
        CSV in Sage data-set format.

    Returns
    -------
    dict[tuple[str, str], list]
        Points per sorted SMILES pair, in file order.
    """
    systems = defaultdict(list)
    with path.open() as handle:
        for row in csv.DictReader(handle):
            key, point = select_systems.parse_sage_row(row, "fast-set")
            systems[key].append(point)
    return systems


def fast_set_outcome(key, points: list, candidates_by_pair: dict, selection: dict):
    """
    Describe what happened to one system of the fast-converging list, and why.

    Parameters
    ----------
    key
        Sorted SMILES pair.
    points
        The list's points for this system.
    candidates_by_pair
        Candidate rows keyed by sorted SMILES pair.
    selection
        Parsed ``selection.yml``.

    Returns
    -------
    str
        "selected" (with a note on wrong temperatures) or "excluded: <reason>".
    """
    chosen = {sorted_pair(s["a"], s["b"]): s for s in selection["systems"]}
    dropped = {
        sorted_pair(d["a"], d["b"]): d["reason"]
        for d in selection.get("dropped", [])
        if "a" in d
    }
    two_phase = {
        sorted_pair(t["a"], t["b"]): t for t in selection.get("two_phase_fast_set", [])
    }
    if key in chosen:
        outcome = "selected"
        corrections = chosen[key].get("temperature_corrections", {})
        for doi in sorted({p.doi for p in points if p.doi in corrections}):
            listed = sorted({p.temperature for p in points if p.doi == doi})
            if listed != [corrections[doi]]:
                outcome += (
                    f". Note: the list gives "
                    f"{', '.join(f'{t:.2f}' for t in listed)} K for "
                    f"[{doi}](https://doi.org/{doi}), which is incorrect; the "
                    f"paper reports measurements at {corrections[doi]:.2f} K"
                )
        return outcome
    if key in dropped:
        return f"excluded: {dropped[key]}"
    if key in two_phase:
        return f"excluded: two-phase. {two_phase[key].get('lle', '')}"
    candidate = candidates_by_pair.get(key)
    if candidate and candidate["failures"]:
        return f"excluded: {candidate['failures']}"
    return "excluded (see dropped systems)"


def fast_set_section(path: Path, candidates: list[dict], selection: dict) -> list[str]:
    """
    Assess a recommended fast-converging system list.

    Parameters
    ----------
    path
        CSV in Sage data-set format.
    candidates
        Candidate rows.
    selection
        Parsed ``selection.yml``.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    candidates_by_pair = {(c["a"], c["b"]): c for c in candidates}
    lines = [
        "## Recommended fast-converging systems",
        "",
        "| System | Points (ΔHmix/ρ) | Linear H^E (two-phase signature) | "
        "Slopes between points | Outcome |",
        "|---|---|---|---|---|",
    ]
    for key, points in read_fast_set(path).items():
        hmix = sorted(
            (p for p in points if p.property_name == "hmix"), key=lambda p: p.x
        )
        slopes = (
            np.diff([p.value for p in hmix]) / np.diff([p.x for p in hmix])
            if len(hmix) > 1
            else []
        )
        linear = select_systems.two_phase_signature(hmix)
        outcome = fast_set_outcome(key, points, candidates_by_pair, selection)
        lines.append(
            f"| `{key[0]}` + `{key[1]}` | {len(hmix)}/{len(points) - len(hmix)} | "
            f"{'yes' if linear else 'no'} | "
            f"{', '.join(f'{slope:.3f}' for slope in slopes)} | {outcome} |"
        )
    lines.append("")
    for system in selection.get("two_phase_fast_set", []):
        lines.append(
            f"- **{' + '.join(system['names'])}**: {system.get('lle', 'TODO')}"
        )
    lines.append("")
    return lines


def selected_systems_section(systems: list[dict], candidates: list[dict]) -> list[str]:
    """
    Tabulate the selected systems with their filter status.

    Parameters
    ----------
    systems
        System dictionaries.
    candidates
        Candidate rows.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    candidates_by_pair = {(c["a"], c["b"]): c for c in candidates}
    lines = [
        "## Selected systems",
        "",
        "| # | System | T (K) | ΔHmix at box x (kJ/mol) | V^E at box x (cm³/mol) | "
        "Hard filters | Interaction |",
        "|---|---|---|---|---|---|---|",
    ]
    for number, system in enumerate(systems, 1):
        candidate = candidates_by_pair.get(tuple(sorted((system["a"], system["b"]))))
        if candidate and candidate["passes"]:
            status = "pass"
        else:
            failed = candidate["failures"] if candidate else "not in the screened pool"
            reason = system["spec"].get("exception", "TODO").replace("|", "\\|")
            status = f"exception. Fails: {failed}. {reason}"
        hmix = " / ".join(format_number(r["hmix_exp"], 2) for r in system["rows"])
        vex = " / ".join(format_number(r["vex_exp"], 2) for r in system["rows"])
        lines.append(
            f"| {number} | {display_name(system)} | {system['temperature']:.2f} | "
            f"{hmix} | {vex} | {status} | {system['spec']['interaction']} |"
        )
    n_mixture = sum(len(s["rows"]) for s in systems)
    pure_keys = {(c, s["temperature"]) for s in systems for c in (s["a"], s["b"])}
    lines += [
        "",
        f"{len(systems)} systems, {n_mixture} mixture boxes and "
        f"{len(pure_keys)} pure boxes "
        f"({n_mixture + len(pure_keys)} in total).",
        "",
    ]
    return lines


def structures_section(systems: list[dict], figure_dir: Path) -> list[str]:
    """
    Draw both components of every system and tabulate the images.

    Parameters
    ----------
    systems
        System dictionaries.
    figure_dir
        Figure directory; structures are written to ``structures/`` inside it.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    structure_dir = figure_dir / "structures"
    structure_dir.mkdir(exist_ok=True)
    lines = [
        "### Structures",
        "",
        "Component a is the one whose mole fraction x is quoted.",
        "",
        "| # | Component a | Component b |",
        "|---|---|---|",
    ]
    for number, system in enumerate(systems, 1):
        cells = []
        for key, name in zip(("a", "b"), system["spec"]["names"], strict=True):
            path = structure_dir / f"{system['spec']['id']}_{key}.png"
            Draw.MolToFile(Chem.MolFromSmiles(system[key]), str(path), size=(240, 180))
            cells.append(
                f"![{name}](figures/structures/{path.name})<br>{name}<br>`{system[key]}`"
            )
        lines.append(f"| {number} | {cells[0]} | {cells[1]} |")
    return [*lines, ""]


def vex_coverage_section(systems: list[dict]) -> list[str]:
    """
    Count the systems by the sign of V^E at their box compositions.

    Parameters
    ----------
    systems
        System dictionaries.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    signs = Counter()
    for system in systems:
        values = [r["vex_exp"] for r in system["rows"] if r["vex_exp"] is not None]
        if not values:
            signs["no density data"] += 1
        elif all(v < 0 for v in values):
            signs["negative"] += 1
        elif all(v > 0 for v in values):
            signs["positive"] += 1
        else:
            signs["sign change"] += 1
    return [
        "### V^E coverage",
        "",
        "V^E diversity was not a selection criterion. Sign at the three box "
        "compositions:",
        "",
        *[f"- {sign}: {count}" for sign, count in signs.most_common()],
        "",
    ]


def box_plan_section(systems: list[dict]) -> list[str]:
    """
    Tabulate box sizes and the rule that sets each one.

    Parameters
    ----------
    systems
        System dictionaries.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    binding = Counter()
    box_rows = []
    for system in systems:
        for row in system["rows"]:
            n_total = row["n_a"] + row["n_b"]
            if row["edge"] < select_systems.MIN_EDGE + 0.15:
                reason = "edge"
            elif n_total == select_systems.MIN_MOLECULES:
                reason = "total"
            else:
                reason = "minority / x rounding"
            binding[reason] += 1
            box_rows.append(
                f"| {display_name(system)} | {row['x_box']:.3f} | "
                f"{row['n_a']} / {row['n_b']} | {n_total} | {row['n_atoms']} | "
                f"{row['edge']:.1f} | {reason} |"
            )
    atoms = [row["n_atoms"] for system in systems for row in system["rows"]]
    return [
        "### Box plan",
        "",
        f"Boxes have edge ≥ {select_systems.MIN_EDGE} Å, ≥ "
        f"{select_systems.MIN_MOLECULES} molecules and ≥ "
        f"{select_systems.MIN_MINORITY} minority molecules, with x within "
        f"{select_systems.MAX_X_SHIFT} of the "
        "measured composition. Evaluator, by comparison, uses a fixed total "
        "(1000 molecules for Sage), rounds x·N per component and only requires "
        "at least one molecule of each, so at the compositions used here its "
        "minority counts are about 200 or more.",
        "",
        f"Mixture boxes: {min(atoms)}–{max(atoms)} atoms (median "
        f"{int(statistics.median(atoms))}); total {sum(atoms)} atoms across "
        f"{len(atoms)} boxes. Constraint that sets the size: "
        + ", ".join(f"{reason} {count}" for reason, count in binding.most_common())
        + ".",
        "",
        "| System | x (box) | Molecules a / b | Total | Atoms | Edge (Å) | Set by |",
        "|---|---|---|---|---|---|---|",
        *box_rows,
    ]


def published_mean_error(
    rows: list[dict], force_field: str, property_name: str, digits: int
) -> str:
    """
    Compute the MAE of one published Sage version for one property.

    Parameters
    ----------
    rows
        Published Sage results of one system.
    force_field
        Force-field label in the benchmark files.
    property_name
        ``hmix`` or ``density``.
    digits
        Decimal places.

    Returns
    -------
    str
        Formatted MAE, or an en dash if there are no results.
    """
    errors = [
        abs(r["value"] - r["reference"])
        for r in rows
        if r["forcefield"] == force_field and r["property_name"] == property_name
    ]
    return format_number(statistics.mean(errors), digits) if errors else "–"


def sage_baseline_section(systems: list[dict]) -> list[str]:
    """
    Tabulate published Sage errors per system.

    Parameters
    ----------
    systems
        System dictionaries.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    columns = {ff: ff for ff in select_systems.SAGE_FORCE_FIELDS} | {
        select_systems.SAGE_RC2: "openff-2.3.0rc2"
    }

    lines = [
        "",
        "## Sage baseline (published, 1000-molecule boxes)",
        "",
        "| System | Max rep. SD ΔHmix | Max rep. SD ρ | "
        + " | ".join(f"{label} ΔHmix MAE" for label in columns.values())
        + " | openff-2.2.1 ρ MAE |",
        "|---" * (4 + len(columns)) + "|",
    ]
    for system in systems:
        rows = system["sage"]
        if not rows:
            lines.append(
                f"| {display_name(system)} | – | – | "
                + " | ".join("–" for _ in columns)
                + " | – |"
            )
            continue
        # SDs over the same force fields as the repeat-SD filter.
        sd_hmix = select_systems.largest_sage_sd(rows, "hmix")
        sd_density = select_systems.largest_sage_sd(rows, "density")
        hmix_maes = [published_mean_error(rows, ff, "hmix", 3) for ff in columns]
        lines.append(
            f"| {display_name(system)} | {format_number(sd_hmix)} | "
            f"{format_number(sd_density, 4)} | "
            + " | ".join(hmix_maes)
            + f" | {published_mean_error(rows, 'openff-2.2.1', 'density', 4)} |"
        )
    return lines


def pure_densities_section(systems: list[dict], mlpeg: dict) -> list[str]:
    """
    Tabulate the chosen pure densities and any disagreeing sources.

    Parameters
    ----------
    systems
        System dictionaries.
    mlpeg
        The ml-peg pure liquids.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    lines = [
        "",
        "## Pure densities",
        "",
        "Each pure density is taken from the paper that supplies most of the mixture "
        "densities, so V^E is self-consistent, unless an explicit reference is given "
        "(water from [IAPWS-95](https://doi.org/10.1063/1.1461829)). Other sources "
        "disagreeing by more than "
        f"{select_systems.PURE_DENSITY_TOLERANCE} g/mL are listed.",
        "",
        "| Component | T (K) | ρ (g/mL) | σ | Source | In ml-peg | ml-peg value | "
        "Disagreeing sources |",
        "|---|---|---|---|---|---|---|---|",
    ]
    seen = set()
    for system in systems:
        for smiles in (system["a"], system["b"]):
            if (smiles, system["temperature"]) in seen:
                continue
            seen.add((smiles, system["temperature"]))
            pure = system["pures"][smiles]
            liquid = mlpeg.get(smiles)
            mlpeg_text = (
                f"{liquid['density']:.4f} @ {liquid['temperature']:.0f} K"
                if liquid
                else "–"
            )
            lines.append(
                f"| `{smiles}` | {system['temperature']:.2f} | "
                f"{format_number(pure['density'], 5)} | "
                f"{format_number(pure['sigma'], 5)} | {pure['source']} | "
                f"{'yes' if liquid else 'no'} | "
                f"{mlpeg_text} | {', '.join(pure['outliers']) or '–'} |"
            )
    lines += [
        "",
        "Many ml-peg liquid_densities references appear to be measured about 5 K "
        "below their labelled temperature (20 °C values labelled 298 K, 25 °C values "
        "labelled 303 K), so they are not used here; see "
        "[`docs/data_issues.md`](../docs/data_issues.md) section 1.1. "
        "Pure densities of components already in ml-peg are simulated (needed for "
        "ΔHmix) but not re-scored.",
        "",
    ]
    return lines


def system_details(number: int, system: dict, candidate: dict | None) -> list[str]:
    """
    Describe the data sources, fit quality and boxes of one system.

    Parameters
    ----------
    number
        Position in the selection.
    system
        System dictionary.
    candidate
        Its screening row, or ``None`` if it was not screened.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    spec = system["spec"]
    hmix_dois = Counter(p.doi for p in system["hmix"])
    density_dois = Counter(p.doi for p in system["density"])
    lines = [
        f"### {number}. {display_name(system)} (`{system['a']}` + `{system['b']}`)",
        "",
        f"- **Interaction:** {spec['interaction']}",
        f"- **Rationale:** {spec['rationale']}",
        f"- **Miscibility:** {spec.get('miscibility', 'TODO')}",
    ]
    if spec.get("exception"):
        lines.append(f"- **Filter exception:** {spec['exception']}")
    if candidate and candidate["failures"]:
        lines.append(f"- **Failed filters:** {candidate['failures']}")
    lines += [
        f"- **ΔHmix data:** {sum(hmix_dois.values())} points from "
        + ", ".join(f"{doi} ({n})" for doi, n in hmix_dois.most_common()),
        f"- **Density data:** {sum(density_dois.values())} points from "
        + (
            ", ".join(f"{doi} ({n})" for doi, n in density_dois.most_common()) or "none"
        ),
    ]
    if system["x_range"]:
        chi2 = system["linear_chi2"]
        lines.append(
            f"- **Single-phase evidence:** ΔHmix measured at {len(system['hmix'])} "
            "compositions spanning "
            f"x = {system['x_range'][0]:.2f}–{system['x_range'][1]:.2f}; "
            "straight-line fit reduced χ² = "
            + (f"{chi2:.0f}" if chi2 is not None else "n/a (too few points)")
            + " (a two-phase region would give ≈ 1)."
        )
    derived = [p for p in system["vex"] if p.source == "derived"]
    if derived and system["vex_fit"]:
        rms = math.sqrt(
            sum((p.value - system["vex_fit"](p.x)) ** 2 for p in derived) / len(derived)
        )
        lines.append(
            f"- **V^E fit:** {len(derived)} points derived from densities; RMS "
            f"residual {rms:.4f} cm³/mol, median σ "
            f"{statistics.median(p.sigma or 0.0 for p in derived):.4f} cm³/mol."
        )
    if system["noisy_density"]:
        lines.append(
            "- **Density points not used for V^E** (σ > "
            f"{select_systems.MAX_DENSITY_SIGMA} g/mL): "
            + ", ".join(f"{doi} ({n})" for doi, n in system["noisy_density"].items())
        )
    if system["errors"]:
        lines.append(f"- **ThermoML errors:** {system['errors']}")
    lines += [
        "",
        "| x box | n_a / n_b | atoms | edge (Å) | ΔHmix ± σ (kJ/mol) | "
        "ρ ± σ (g/mL) | V^E (cm³/mol) |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in system["rows"]:
        lines.append(
            f"| {row['x_box']:.3f} | {row['n_a']} / {row['n_b']} | {row['n_atoms']} | "
            f"{format_number(row['edge'], 1)} | {format_number(row['hmix_exp'])} ± "
            f"{format_number(row['hmix_sigma'], 4)} | "
            f"{format_number(row['density_exp'], 4)} ± "
            f"{format_number(row['density_sigma'], 5)} | "
            f"{format_number(row['vex_exp'])} |"
        )
    lines += ["", f"![{spec['id']}](figures/{spec['id']}.png)", ""]
    return lines


def known_gaps_section(points: dict) -> list[str]:
    """
    List gaps in chemical coverage, with the amide systems the Sage sets have.

    Parameters
    ----------
    points
        Sage data-set points.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    amide_systems = [
        (key, system_points)
        for key, system_points in points.items()
        if key[1] and set(key) & AMIDE_SMILES
    ]
    lines = [
        "## Known gaps",
        "",
        "- **No N-H (secondary) amide.** Systems with N-methylformamide, "
        "N-methylacetamide, acetamide or formamide in the Sage sets:",
        "",
    ]
    if amide_systems:
        for key, system_points in amide_systems:
            with_sigma = sum(1 for p in system_points if p.sigma)
            lines.append(
                f"  - `{key[0]}` + `{key[1]}`: {len(system_points)} points "
                f"({with_sigma} with uncertainties) from "
                + ", ".join(sorted({p.source for p in system_points}))
            )
    else:
        lines.append("  - none")
    lines += [
        "",
        "  The Sage 2.0 set drops their uncertainties and has enthalpies only at "
        "308.15 K, but the ThermoML records carry full uncertainties (and the "
        "standard-InChI conversion would otherwise turn N-methylformamide into "
        "`CN=CO`; `thermoml.py` repairs amide tautomers). See reserves.",
        "- **No amide-water system.** 2-pyrrolidone + water was dropped (see below).",
        "- **No ester.** Methyl pentanoate + octane was dropped for cost (see below).",
        "- **Selection bias.** Every system comes from the Sage training/validation "
        "pools, so the set reflects ThermoML coverage and Sage's curation choices "
        "rather than a designed spread.",
        "- **Endothermic range is capped** at RT/2 by the two-phase screen.",
        "",
    ]
    return lines


def dropped_and_reserves_section(selection: dict) -> list[str]:
    """
    List dropped systems, reserves and the command to reproduce the report.

    Parameters
    ----------
    selection
        Parsed ``selection.yml``.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    return [
        "## Dropped systems",
        "",
        *[f"- **{d['system']}**: {d['reason']}" for d in selection.get("dropped", [])],
        "",
        "## Reserves",
        "",
        *[
            f"- **{r['system']}**: {r['note']}" if "note" in r else f"- {r['system']}"
            for r in selection.get("reserves", [])
        ],
        "",
        "## Reproduce",
        "",
        "```bash",
        "pixi run -e analysis select --output-dir selection \\",
        "    [--fast-set selection/inputs/dataset-fast-converging.csv]",
        "```",
        "",
    ]


def write_report(
    out: Path,
    selection: dict,
    candidates: list[dict],
    systems: list[dict],
    points: dict,
    sage: dict,
    mlpeg: dict,
    fast_set: Path | None,
) -> Path:
    """
    Write ``report.md`` and its figures.

    Parameters
    ----------
    out
        Output directory.
    selection
        Parsed ``selection.yml``.
    candidates
        Candidate rows.
    systems
        System dictionaries.
    points
        Sage data-set points.
    sage
        Published Sage benchmark results (kept for the call signature; the
        per-system results are in ``systems``).
    mlpeg
        The ml-peg pure liquids.
    fast_set
        Optional fast-converging system CSV.

    Returns
    -------
    Path
        Path to the report.
    """
    figure_dir = out / "figures"
    figure_dir.mkdir(exist_ok=True)
    plot_repeat_sd(candidates, figure_dir / "sage_repeat_sd.png")
    for system in systems:
        plot_system(system, figure_dir / f"{system['spec']['id']}.png")

    candidates_by_pair = {(c["a"], c["b"]): c for c in candidates}
    lines = [
        "# Binary mixtures benchmark: system selection report",
        "",
        f"Generated {date.today().isoformat()} by `select_systems.py`.",
        "",
        *hard_filters_section(),
        *data_sources_section(),
        *funnel_section(candidates, points),
    ]
    if fast_set:
        lines += fast_set_section(fast_set, candidates, selection)
    lines += [
        *selected_systems_section(systems, candidates),
        *structures_section(systems, figure_dir),
        *vex_coverage_section(systems),
        *box_plan_section(systems),
        *sage_baseline_section(systems),
        *pure_densities_section(systems, mlpeg),
        "## Per-system details",
        "",
    ]
    for number, system in enumerate(systems, 1):
        candidate = candidates_by_pair.get(tuple(sorted((system["a"], system["b"]))))
        lines += system_details(number, system, candidate)
    lines += [
        *known_gaps_section(points),
        *dropped_and_reserves_section(selection),
    ]
    path = out / "report.md"
    path.write_text("\n".join(lines))
    return path
