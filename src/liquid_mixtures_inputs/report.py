"""Write the stage 1 selection report (Markdown plus PNG figures)."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from datetime import date
from pathlib import Path
import statistics

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import Draw  # noqa: E402

from liquid_mixtures_inputs import select_systems as ss  # noqa: E402
from liquid_mixtures_inputs.thermoml import canonical_smiles  # noqa: E402

AMIDE_SMILES = {canonical_smiles(s) for s in ("CNC=O", "CC(=O)NC", "CC(=O)N", "NC=O")}
_FF_DIR = f"{ss.RAW}/ash-sage-rc2/{ss.ASH_COMMIT}/04_benchmark/forcefields"
RC2_OFFXML = f"{_FF_DIR}/sage-2.3.0rc2-v3-k100.offxml"
FB_FIT_OFFXML = f"{_FF_DIR}/{ss.SAGE_RC2}.offxml"


def _fmt(value, digits: int = 3) -> str:
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


def _github(raw_url: str) -> str:
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
    repo, path = raw_url.removeprefix(ss.RAW + "/").split("/", 1)
    return f"https://github.com/openforcefield/{repo}/blob/{path}?plain=1"


def _names(system: dict) -> str:
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
    for ax, (prop, unit) in zip(
        axes, (("hmix", "kJ/mol"), ("density", "g/mL")), strict=True
    ):
        values = [
            c[f"sage_rep_sd_{prop}"] for c in candidates if c[f"sage_rep_sd_{prop}"]
        ]
        ax.hist(values, bins=40, color="#4a7ab5")
        ax.axvline(ss.MAX_REPEAT_SD[prop], color="#c0392b", ls="--", label="limit")
        ax.set_xlabel(f"max Sage 3-replicate SD ({unit})")
        ax.set_ylabel("systems")
        ax.set_title("ΔHmix" if prop == "hmix" else "density")
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
        for p in system[key]:
            by_doi[p.doi].append(p)
        for doi, pts in sorted(by_doi.items(), key=lambda kv: -len(kv[1]))[:6]:
            ax.errorbar(
                [p.x for p in pts],
                [p.value for p in pts],
                yerr=[p.sigma or 0.0 for p in pts],
                fmt="o",
                ms=3,
                label=doi,
            )
        if system[fit_key]:
            ax.plot(grid, [system[fit_key](x) for x in grid], color="k", lw=1)
        targets = [r["x_box"] for r in system["rows"]]
        values = [r[f"{key}_exp"] for r in system["rows"]]
        if all(v is not None for v in values):
            ax.plot(
                targets, values, "s", mfc="none", mec="#c0392b", ms=9, label="box x"
            )
        ax.axhline(0, color="0.7", lw=0.5)
        ax.set_xlim(0, 1)
        ax.set_xlabel(f"x ({system['spec']['names'][0]})")
        ax.set_ylabel(label)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=6, frameon=False)
    fig.suptitle(f"{_names(system)} at {system['temperature']:.2f} K", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def _sage_table(systems: list[dict]) -> list[str]:
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
    columns = {ff: ff for ff in ss.SAGE_FORCE_FIELDS} | {ss.SAGE_RC2: "openff-2.3.0rc2"}
    lines = [
        "| System | Max rep. SD ΔHmix | Max rep. SD ρ | "
        + " | ".join(f"{label} ΔHmix MAE" for label in columns.values())
        + " | openff-2.2.1 ρ MAE |",
        "|---" * (4 + len(columns)) + "|",
    ]
    for s in systems:
        rows = s["sage"]
        if not rows:
            lines.append(
                f"| {_names(s)} | – | – | "
                + " | ".join("–" for _ in columns)
                + " | – |"
            )
            continue
        # SDs over the same force fields as the repeat-SD filter.
        filtered = [r for r in rows if r["forcefield"] in ss.SAGE_FORCE_FIELDS]
        sd_h = max(
            (r["sd"] for r in filtered if r["prop"] == "hmix" and r["sd"]),
            default=None,
        )
        sd_d = max(
            (r["sd"] for r in filtered if r["prop"] == "density" and r["sd"]),
            default=None,
        )
        maes = []
        for ff in columns:
            errs = [
                abs(r["value"] - r["reference"])
                for r in rows
                if r["forcefield"] == ff and r["prop"] == "hmix"
            ]
            maes.append(_fmt(statistics.mean(errs)) if errs else "–")
        d_errs = [
            abs(r["value"] - r["reference"])
            for r in rows
            if r["forcefield"] == "openff-2.2.1" and r["prop"] == "density"
        ]
        lines.append(
            f"| {_names(s)} | {_fmt(sd_h)} | {_fmt(sd_d, 4)} | "
            + " | ".join(maes)
            + f" | {_fmt(statistics.mean(d_errs), 4) if d_errs else '–'} |"
        )
    return lines


def _fast_set_section(path: Path, candidates: list[dict], selection: dict) -> list[str]:
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
    systems = defaultdict(list)
    with path.open() as handle:
        for row in csv.DictReader(handle):
            a, b = ss.fix_smiles(row["Component 1"]), ss.fix_smiles(row["Component 2"])
            x = float(row["Mole Fraction 1"])
            if b < a:
                a, b, x = b, a, 1 - x
            prop = "hmix" if row["EnthalpyOfMixing Value (kJ / mol)"] else "density"
            value = (
                row["EnthalpyOfMixing Value (kJ / mol)"]
                or row["Density Value (g / ml)"]
            )
            sigma = (
                row["EnthalpyOfMixing Uncertainty (kJ / mol)"]
                or row["Density Uncertainty (g / ml)"]
            )
            systems[(a, b)].append(
                ss.Point(
                    prop,
                    x,
                    float(row["Temperature (K)"]),
                    float(value),
                    float(sigma) if sigma else None,
                    row["Source"],
                    "fast-set",
                )
            )
    by_pair = {(c["a"], c["b"]): c for c in candidates}
    chosen = {
        tuple(sorted((canonical_smiles(s["a"]), canonical_smiles(s["b"])))): s
        for s in selection["systems"]
    }
    dropped = {
        tuple(sorted((canonical_smiles(d["a"]), canonical_smiles(d["b"])))): d["reason"]
        for d in selection.get("dropped", [])
        if "a" in d
    }
    lle = {
        tuple(sorted((canonical_smiles(t["a"]), canonical_smiles(t["b"])))): t
        for t in selection.get("two_phase_fast_set", [])
    }
    lines = [
        "| System | Points (ΔHmix/ρ) | Linear H^E (two-phase signature) | "
        "Slopes between points | Outcome |",
        "|---|---|---|---|---|",
    ]
    for key, pts in systems.items():
        h = sorted((p for p in pts if p.prop == "hmix"), key=lambda p: p.x)
        slopes = (
            np.diff([p.value for p in h]) / np.diff([p.x for p in h])
            if len(h) > 1
            else []
        )
        linear = ss.two_phase_signature(h)
        if key in chosen:
            outcome = "selected"
            corrections = chosen[key].get("temperature_corrections", {})
            for doi in sorted({p.doi for p in pts if p.doi in corrections}):
                listed = sorted({p.temperature for p in pts if p.doi == doi})
                if listed != [corrections[doi]]:
                    outcome += (
                        f". Note: the list gives "
                        f"{', '.join(f'{t:.2f}' for t in listed)} K for "
                        f"[{doi}](https://doi.org/{doi}), which is incorrect; the "
                        f"paper reports measurements at {corrections[doi]:.2f} K"
                    )
        elif key in dropped:
            outcome = f"excluded: {dropped[key]}"
        elif key in lle:
            outcome = f"excluded: two-phase. {lle[key].get('lle', '')}"
        else:
            cand = by_pair.get(key)
            outcome = (
                f"excluded: {cand['failures']}"
                if cand and cand["failures"]
                else "excluded (see dropped systems)"
            )
        lines.append(
            f"| `{key[0]}` + `{key[1]}` | {len(h)}/{len(pts) - len(h)} | "
            f"{'yes' if linear else 'no'} | "
            f"{', '.join(f'{s:.3f}' for s in slopes)} | {outcome} |"
        )
    return lines


def _structures_section(systems: list[dict], fig_dir: Path) -> list[str]:
    """
    Draw both components of every system and tabulate the images.

    Parameters
    ----------
    systems
        System dictionaries.
    fig_dir
        Figure directory; structures are written to ``structures/`` inside it.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    struct_dir = fig_dir / "structures"
    struct_dir.mkdir(exist_ok=True)
    lines = [
        "### Structures",
        "",
        "Component a is the one whose mole fraction x is quoted.",
        "",
        "| # | Component a | Component b |",
        "|---|---|---|",
    ]
    for i, s in enumerate(systems, 1):
        cells = []
        for key, name in zip(("a", "b"), s["spec"]["names"], strict=True):
            path = struct_dir / f"{s['spec']['id']}_{key}.png"
            Draw.MolToFile(Chem.MolFromSmiles(s[key]), str(path), size=(240, 180))
            cells.append(
                f"![{name}](figures/structures/{path.name})<br>{name}<br>`{s[key]}`"
            )
        lines.append(f"| {i} | {cells[0]} | {cells[1]} |")
    return [*lines, ""]


def _hard_filters_section() -> list[str]:
    """
    Describe the hard filters, using the thresholds from ``select_systems``.

    Returns
    -------
    list[str]
        Markdown lines.
    """
    t_low, t_high = ss.SCREEN_T_RANGE
    rows = [
        (
            "(entry) ≥ 3 ΔHmix and ≥ 3 density points",
            f"At {t_low}–{t_high} K, counted per distinct composition.",
            "Three compositions per system need a measured ΔHmix and density.",
        ),
        (
            "composition mismatch",
            "Every ΔHmix composition has a density point within "
            f"\\|Δx\\| ≤ {ss.MAX_COMPOSITION_MISMATCH}.",
            "Both properties are scored on the same box.",
        ),
        (
            "missing experimental uncertainty",
            "Every ΔHmix and density point reports an uncertainty.",
            "Needed for the MAZE denominator and to judge agreement.",
        ),
        (
            f"relative ΔHmix uncertainty > {100 * ss.MAX_REL_SIGMA_HMIX:.0f} %",
            "Largest relative ΔHmix uncertainty "
            f"≤ {100 * ss.MAX_REL_SIGMA_HMIX:.0f} %.",
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
            f"≤ {ss.MAX_HEAVY_ATOMS} heavy atoms per molecule.",
            "Keeps boxes small at MLIP cost while holding enough molecules.",
        ),
        (
            "viscous",
            "No diol, polyol, polyamine or bis-hydroxyethyl amine motifs "
            "(SMARTS in `VISCOUS_SMARTS`).",
            "Slow dynamics would not converge in 2 ns.",
        ),
        (
            "two-phase (linear H^E)",
            "ΔHmix is not consistent with a straight line in x (three points: "
            "consecutive slopes differ by > 2σ; more points: reduced χ² of a linear "
            "fit > 2). σ has a 2 % relative floor. Applied only when the largest "
            f"measured \\|H^E\\| ≥ {ss.MIN_SIGNAL} kJ/mol.",
            "A two-phase region gives H^E = a + b·x (lever rule). An ideal or "
            "near-ideal mixture is also a straight line (H^E ≈ 0), so the test is "
            "skipped for small H^E; the small-signal filter handles those.",
        ),
        (
            "two-phase (H^E > RT/2)",
            f"ΔHmix(x≈0.5) ≤ RT/2 = {ss.UCST_LIMIT:.2f} kJ/mol.",
            "Regular-solution critical value: H^E = βx(1 − x) splits into two "
            "phases below T_c = β/(2R) ([derivation]"
            "(https://scipython.com/blog/the-miscibility-of-a-regular-solution/)), "
            "and H^E(0.5) = β/4 = RT/2 at T_c. A conservative screen for an "
            "upper critical solution temperature above 298 K.",
        ),
        (
            "Sage repeat SD (hmix / density)",
            "Largest published Sage 3-replicate SD over openff-2.0.0 to 2.3.0rc1 "
            f"≤ {ss.MAX_REPEAT_SD['hmix']} kJ/mol (ΔHmix) and "
            f"≤ {ss.MAX_REPEAT_SD['density']} g/mL (density).",
            "A large spread between replicates suggests slow equilibration.",
        ),
        (
            "no Sage repeat SD (hmix / density)",
            "A published Sage SD exists for the property.",
            "Without it the convergence check above cannot be made.",
        ),
        (
            "small signal",
            f"\\|ΔHmix(x≈0.5)\\| ≥ {ss.MIN_SIGNAL} kJ/mol.",
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
        Published Sage benchmark results.
    mlpeg
        The ml-peg pure liquids.
    fast_set
        Optional fast-converging system CSV.

    Returns
    -------
    Path
        Path to the report.
    """
    fig_dir = out / "figures"
    fig_dir.mkdir(exist_ok=True)
    plot_repeat_sd(candidates, fig_dir / "sage_repeat_sd.png")
    for s in systems:
        plot_system(s, fig_dir / f"{s['spec']['id']}.png")

    n_binary = sum(1 for k in points if k[1])
    failures = Counter(f for c in candidates for f in c["failures"].split("; ") if f)
    by_pair = {(c["a"], c["b"]): c for c in candidates}

    lines = [
        "# Binary mixtures benchmark: system selection report",
        "",
        f"Generated {date.today().isoformat()} by `select_systems.py`.",
        "",
        *_hard_filters_section(),
        "## Data sources",
        "",
        "- OpenFF Sage 2.0 training set "
        f"([openff-sage@{ss.SAGE_COMMIT[:7]}]({_github(ss.SOURCES['sage-2.0-train'])})): "
        "no experimental uncertainties.",
        "- Sage 2.3 refit "
        f"[training]({_github(ss.SOURCES['rc2-train'])}) and "
        f"[validation]({_github(ss.SOURCES['rc2-validation'])}) sets "
        f"(ash-sage-rc2@{ss.ASH_COMMIT[:7]}): "
        "experimental uncertainties for every point.",
        "- Published Sage benchmark results (ash-sage-rc2 "
        f"[`04_benchmark`]({_github(ss.SAGE_BENCHMARKS['training-hmix'])})): mean "
        "and SD over 3 independent replicates for openff-2.0.0, 2.1.0, 2.2.1, "
        "2.3.0rc1 and 2.3.0rc2 (labelled `fb-fit-v3-single-mean-k100` in the "
        "benchmark files; "
        f"[`sage-2.3.0rc2-v3-k100.offxml`]({_github(RC2_OFFXML)}) is identical to "
        f"[`fb-fit-v3-single-mean-k100.offxml`]({_github(FB_FIT_OFFXML)})). "
        "2.3.0rc2 is shown in the Sage baseline table only and is not used by "
        "the repeat-SD filter.",
        "- Full isotherms and pure densities from the NIST TRC ThermoML archive, "
        "fetched by DOI. ThermoML stores either a standard uncertainty (one "
        "standard deviation) or an expanded uncertainty (a 95 % interval, "
        "coverage factor k = 2). Standard uncertainties are used as given; "
        "expanded ones are divided by 2. Every σ in this report is therefore a "
        "standard uncertainty. The ash-sage-rc2 sets use the same convention: "
        "their σ equals half the ThermoML expanded value.",
        "- ml-peg's liquid_densities inputs, to avoid re-scoring pure densities.",
        "",
        "## Screening funnel",
        "",
        "- Distinct binary systems across the Sage 2.0 training set and the "
        f"Sage 2.3 rc2 training and validation sets: {n_binary}.",
        f"- With ≥ 3 ΔHmix and ≥ 3 density points at {ss.SCREEN_T_RANGE[0]}–"
        f"{ss.SCREEN_T_RANGE[1]} K: {len(candidates)}.",
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
        f"b·x) and whose largest |H^E| is at least {ss.MIN_SIGNAL} kJ/mol. Systems it "
        "flags in the pool, with the largest relative ΔHmix uncertainty: "
        + ", ".join(
            f"`{c['a']}` + `{c['b']}` "
            f"({_fmt(c['hmix_max_rel_sigma'] and 100 * c['hmix_max_rel_sigma'], 0)} %)"
            for c in candidates
            if c["linear_two_phase"]
        )
        + ". Large relative uncertainties can hide real curvature, so a flag on data "
        "with errors of ~25 % (placeholder values in the rc2 sets, e.g. "
        f"[σ = 0.1695 for ΔHmix = -0.678]({_github(ss.SOURCES['rc2-train'])}#L769)) "
        "is inconclusive "
        "rather than evidence of two phases. Such systems are also rejected by the "
        "relative-uncertainty filter.",
        "",
        "![Sage repeat SD](figures/sage_repeat_sd.png)",
        "",
    ]

    if fast_set:
        lines += [
            "## Recommended fast-converging systems",
            "",
            *_fast_set_section(fast_set, candidates, selection),
            "",
        ]
        for t in selection.get("two_phase_fast_set", []):
            lines.append(f"- **{' + '.join(t['names'])}**: {t.get('lle', 'TODO')}")
        lines.append("")

    lines += [
        "## Selected systems",
        "",
        "| # | System | T (K) | ΔHmix at box x (kJ/mol) | V^E at box x (cm³/mol) | "
        "Hard filters | Interaction |",
        "|---|---|---|---|---|---|---|",
    ]
    for i, s in enumerate(systems, 1):
        key = tuple(sorted((s["a"], s["b"])))
        cand = by_pair.get(key)
        if cand and cand["passes"]:
            status = "pass"
        else:
            failed = cand["failures"] if cand else "not in the screened pool"
            reason = s["spec"].get("exception", "TODO").replace("|", "\\|")
            status = f"exception. Fails: {failed}. {reason}"
        h = " / ".join(_fmt(r["hmix_exp"], 2) for r in s["rows"])
        v = " / ".join(_fmt(r["vex_exp"], 2) for r in s["rows"])
        lines.append(
            f"| {i} | {_names(s)} | {s['temperature']:.2f} | {h} | {v} | {status} | "
            f"{s['spec']['interaction']} |"
        )
    n_mix = sum(len(s["rows"]) for s in systems)
    pure_keys = {(c, s["temperature"]) for s in systems for c in (s["a"], s["b"])}
    lines += [
        "",
        f"{len(systems)} systems, {n_mix} mixture boxes and "
        f"{len(pure_keys)} pure boxes "
        f"({n_mix + len(pure_keys)} in total).",
        "",
        *_structures_section(systems, fig_dir),
        "### V^E coverage",
        "",
        "V^E diversity was not a selection criterion. Sign at the three box "
        "compositions:",
        "",
    ]
    signs = Counter()
    for s in systems:
        vals = [r["vex_exp"] for r in s["rows"] if r["vex_exp"] is not None]
        if not vals:
            signs["no density data"] += 1
        elif all(v < 0 for v in vals):
            signs["negative"] += 1
        elif all(v > 0 for v in vals):
            signs["positive"] += 1
        else:
            signs["sign change"] += 1
    lines += [f"- {k}: {v}" for k, v in signs.most_common()]
    lines += [
        "",
        "### Box plan",
        "",
        f"Boxes have edge ≥ {ss.MIN_EDGE} Å, ≥ {ss.MIN_MOLECULES} molecules and ≥ "
        f"{ss.MIN_MINORITY} minority molecules, with x within {ss.MAX_X_SHIFT} of the "
        "measured composition. Evaluator, by comparison, uses a fixed total "
        "(1000 molecules for Sage), rounds x·N per component and only requires "
        "at least one molecule of each, so at the compositions used here its "
        "minority counts are about 200 or more.",
        "",
    ]
    binding = Counter()
    box_rows = []
    for s in systems:
        for r in s["rows"]:
            n_total = r["n_a"] + r["n_b"]
            if r["edge"] < ss.MIN_EDGE + 0.15:
                reason = "edge"
            elif n_total == ss.MIN_MOLECULES:
                reason = "total"
            else:
                reason = "minority / x rounding"
            binding[reason] += 1
            box_rows.append(
                f"| {_names(s)} | {r['x_box']:.3f} | {r['n_a']} / {r['n_b']} | "
                f"{n_total} | {r['n_atoms']} | {r['edge']:.1f} | {reason} |"
            )
    atoms = [r["n_atoms"] for s in systems for r in s["rows"]]
    lines += [
        f"Mixture boxes: {min(atoms)}–{max(atoms)} atoms (median "
        f"{int(statistics.median(atoms))}); total {sum(atoms)} atoms across "
        f"{len(atoms)} boxes. Constraint that sets the size: "
        + ", ".join(f"{k} {v}" for k, v in binding.most_common())
        + ".",
        "",
        "| System | x (box) | Molecules a / b | Total | Atoms | Edge (Å) | Set by |",
        "|---|---|---|---|---|---|---|",
        *box_rows,
    ]
    lines += ["", "## Sage baseline (published, 1000-molecule boxes)", ""]
    lines += _sage_table(systems)
    lines += ["", "## Pure densities", ""]
    lines += [
        "Each pure density is taken from the paper that supplies most of the mixture "
        "densities, so V^E is self-consistent, unless an explicit reference is given "
        "(water from [IAPWS-95](https://doi.org/10.1063/1.1461829)). Other sources "
        "disagreeing by more than "
        f"{ss.PURE_DENSITY_TOLERANCE} g/mL are listed.",
        "",
        "| Component | T (K) | ρ (g/mL) | σ | Source | In ml-peg | ml-peg value | "
        "Disagreeing sources |",
        "|---|---|---|---|---|---|---|---|",
    ]
    seen = set()
    for s in systems:
        for c in (s["a"], s["b"]):
            if (c, s["temperature"]) in seen:
                continue
            seen.add((c, s["temperature"]))
            p = s["pures"][c]
            ml = mlpeg.get(c)
            ml_txt = f"{ml['density']:.4f} @ {ml['temperature']:.0f} K" if ml else "–"
            lines.append(
                f"| `{c}` | {s['temperature']:.2f} | {_fmt(p['density'], 5)} | "
                f"{_fmt(p['sigma'], 5)} | {p['source']} | {'yes' if ml else 'no'} | "
                f"{ml_txt} | {', '.join(p['outliers']) or '–'} |"
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
        "## Per-system details",
        "",
    ]
    for i, s in enumerate(systems, 1):
        spec = s["spec"]
        key = tuple(sorted((s["a"], s["b"])))
        cand = by_pair.get(key)
        dois_h = Counter(p.doi for p in s["hmix"])
        dois_d = Counter(p.doi for p in s["density"])
        lines += [
            f"### {i}. {_names(s)} (`{s['a']}` + `{s['b']}`)",
            "",
            f"- **Interaction:** {spec['interaction']}",
            f"- **Rationale:** {spec['rationale']}",
            f"- **Miscibility:** {spec.get('miscibility', 'TODO')}",
        ]
        if spec.get("exception"):
            lines.append(f"- **Filter exception:** {spec['exception']}")
        if cand and cand["failures"]:
            lines.append(f"- **Failed filters:** {cand['failures']}")
        lines += [
            f"- **ΔHmix data:** {sum(dois_h.values())} points from "
            + ", ".join(f"{d} ({n})" for d, n in dois_h.most_common()),
            f"- **Density data:** {sum(dois_d.values())} points from "
            + (", ".join(f"{d} ({n})" for d, n in dois_d.most_common()) or "none"),
        ]
        if s["x_range"]:
            chi2 = s["linear_chi2"]
            lines.append(
                f"- **Single-phase evidence:** ΔHmix measured at {len(s['hmix'])} "
                "compositions spanning "
                f"x = {s['x_range'][0]:.2f}–{s['x_range'][1]:.2f}; "
                "straight-line fit reduced χ² = "
                + (f"{chi2:.0f}" if chi2 is not None else "n/a (too few points)")
                + " (a two-phase region would give ≈ 1)."
            )
        if s["noisy_density"]:
            lines.append(
                "- **Density points not used for V^E** (σ > "
                f"{ss.MAX_DENSITY_SIGMA} g/mL): "
                + ", ".join(f"{d} ({n})" for d, n in s["noisy_density"].items())
            )
        if s["errors"]:
            lines.append(f"- **ThermoML errors:** {s['errors']}")
        lines += [
            "",
            "| x box | n_a / n_b | atoms | edge (Å) | ΔHmix ± σ (kJ/mol) | "
            "ρ ± σ (g/mL) | V^E (cm³/mol) |",
            "|---|---|---|---|---|---|---|",
        ]
        for r in s["rows"]:
            lines.append(
                f"| {r['x_box']:.3f} | {r['n_a']} / {r['n_b']} | {r['n_atoms']} | "
                f"{_fmt(r['edge'], 1)} | {_fmt(r['hmix_exp'])} ± "
                f"{_fmt(r['hmix_sigma'], 4)} | "
                f"{_fmt(r['density_exp'], 4)} ± {_fmt(r['density_sigma'], 5)} | "
                f"{_fmt(r['vex_exp'])} |"
            )
        lines += ["", f"![{spec['id']}](figures/{spec['id']}.png)", ""]

    amide_systems = [
        (k, pts) for k, pts in points.items() if k[1] and set(k) & AMIDE_SMILES
    ]
    lines += [
        "## Known gaps",
        "",
        "- **No N-H (secondary) amide.** Systems with N-methylformamide, "
        "N-methylacetamide, acetamide or formamide in the Sage sets:",
        "",
    ]
    if amide_systems:
        for k, pts in amide_systems:
            with_sigma = sum(1 for p in pts if p.sigma)
            lines.append(
                f"  - `{k[0]}` + `{k[1]}`: {len(pts)} points "
                f"({with_sigma} with uncertainties) from "
                + ", ".join(sorted({p.source for p in pts}))
            )
    else:
        lines.append("  - none")
    lines += [
        "",
        "  The Sage 2.0 set drops their uncertainties and has enthalpies only at "
        "308.15 K, but the ThermoML records carry full uncertainties (and the "
        "standard-InChI conversion would otherwise turn N-methylformamide into "
        "`CN=CO`; `thermoml.py` repairs amide tautomers). See candidate additions.",
        "- **No amide-water system.** 2-pyrrolidone + water was dropped (see below).",
        "- **Selection bias.** Every system comes from the Sage training/validation "
        "pools, so the set reflects ThermoML coverage and Sage's curation choices "
        "rather than a designed spread.",
        "- **Endothermic range is capped** at RT/2 by the two-phase screen.",
        "",
        "## Candidate additions (pending decision)",
        "",
        *[
            f"- **{c['system']}**: {c['note']}"
            for c in selection.get("candidate_additions", [])
        ],
        "",
        "## Dropped systems",
        "",
        *[f"- **{d['system']}**: {d['reason']}" for d in selection.get("dropped", [])],
        "",
        "## Reserves",
        "",
        *[f"- {r}" for r in selection.get("reserves", [])],
        "",
        "## Reproduce",
        "",
        "```bash",
        "pixi run -e analysis select --output-dir selection \\",
        "    [--fast-set selection/inputs/dataset-fast-converging.csv]",
        "```",
        "",
    ]
    path = out / "report.md"
    path.write_text("\n".join(lines))
    return path
