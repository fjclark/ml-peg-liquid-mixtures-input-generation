"""
Screen the NIST ThermoML archive for systems with under-represented groups.

A one-off search beyond the Sage sets: every binary system in the archive
snapshot (papers from about 2003 to 2019) with a component in one of the
``functional_groups`` of ``selection.yml`` is screened at each temperature
from 297 to 315 K with the hard filters that need no Sage results. The ranked
shortlist is written to ``archive_candidates.csv``. A chosen system enters
``selection.yml`` through ``extra_dois`` with an ``exception`` (it has no
published Sage repeat SDs), so the main selection never reads the archive.

Pure densities are taken from the same records as the mixture data, as
binary-mixture papers report their pure liquids.

Run with::

    pixi run -e analysis screen-archive --selection-dir selection
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from functools import cache
import json
from pathlib import Path
import re
import tarfile

from rdkit import Chem
import yaml

from liquid_mixtures_inputs import select_systems as sel
from liquid_mixtures_inputs.thermoml import download, mol_from_inchi, parse_record

ARCHIVE_URL = "https://data.nist.gov/od/ds/mds2-2422/ThermoML.v2020-09-30.tgz"
T_RANGE = (297.0, 315.0)  # K; N-methylacetamide melts near 301 K
FORMULA = re.compile(r"([A-Z][a-z]?)(\d*)")


def formula_passes(inchi: str) -> bool:
    """
    Apply the element and heavy-atom filters to an InChI's formula layer.

    Parameters
    ----------
    inchi
        Standard InChI.

    Returns
    -------
    bool
        ``True`` for a single neutral species within the filters.
    """
    layers = inchi.split("/")
    if len(layers) < 2 or "." in layers[1] or "/q" in inchi or "/p" in inchi:
        return False
    counts = {e: int(n or 1) for e, n in FORMULA.findall(layers[1])}
    heavy = sum(n for e, n in counts.items() if e != "H")
    return set(counts) <= sel.ALLOWED_ELEMENTS and heavy <= sel.MAX_HEAVY_ATOMS


@cache
def component_passes(smiles: str) -> bool:
    """
    Apply the element, heavy-atom and neutral-molecule filters to a component.

    Parameters
    ----------
    smiles
        Canonical SMILES.

    Returns
    -------
    bool
        ``True`` for a single neutral molecule within the filters.
    """
    if "." in smiles or "+" in smiles or "-]" in smiles:
        return False
    info = sel.molecule_info(smiles)
    return (
        info["elements"] <= sel.ALLOWED_ELEMENTS
        and info["heavy_atoms"] <= sel.MAX_HEAVY_ATOMS
    )


@cache
def smiles_from_inchi(inchi: str) -> str | None:
    """
    Convert an InChI to canonical SMILES, as ``parse_record`` does.

    Parameters
    ----------
    inchi
        Standard InChI.

    Returns
    -------
    str | None
        Canonical SMILES, or ``None`` if RDKit cannot parse it.
    """
    mol = mol_from_inchi(inchi)
    return None if mol is None else Chem.MolToSmiles(mol)


def archive_measurements(archive: Path, groups: list[dict]) -> list:
    """
    Parse every archive record with a component in one of ``groups``.

    Parameters
    ----------
    archive
        The archive tarball.
    groups
        ``functional_groups`` entries from ``selection.yml``.

    Returns
    -------
    list[Measurement]
        Liquid measurements from those records, within ``T_RANGE``.
    """
    measurements = []
    with tarfile.open(archive) as tar:
        for member in tar:
            if not member.name.endswith(".json"):
                continue
            record = json.load(tar.extractfile(member))
            inchis = [c.get("sStandardInChI", "") for c in record.get("Compound", [])]
            if not any(
                formula_passes(i)
                and (smiles := smiles_from_inchi(i))
                and sel.functional_groups(smiles, groups)
                for i in inchis
            ):
                continue
            doi = member.name.removesuffix(".json")
            measurements += [
                m
                for m in parse_record(doi, record)
                if T_RANGE[0] <= m.temperature <= T_RANGE[1]
            ]
    return measurements


def screen(
    measurements: list,
    groups: list[dict],
    selected_pure: set[tuple[str, float]],
    two_phase: set[tuple[str, str]] = frozenset(),
) -> list[dict]:
    """
    Apply the hard filters that need no Sage results to each system and T.

    Parameters
    ----------
    measurements
        Output of :func:`archive_measurements`.
    groups
        ``functional_groups`` entries from ``selection.yml``.
    selected_pure
        (SMILES, T) of the pure boxes already in the set.
    two_phase
        Sorted SMILES pairs that the literature shows to be two-phase
        (``two_phase_fast_set`` in ``selection.yml``).

    Returns
    -------
    list[dict]
        One row per system and temperature with ≥ 3 ΔHmix and ≥ 3 density
        points, passing rows first.
    """
    mixtures = defaultdict(list)
    pure = defaultdict(list)
    for m in measurements:
        temperature = round(m.temperature, 2)
        if len(m.smiles) == 1:
            if m.property_name == "density" and m.uncertainty:
                pure[(m.smiles[0], temperature)].append(m)
            continue
        if m.pressure is not None and not (
            sel.PRESSURE_RANGE_KPA[0] <= m.pressure <= sel.PRESSURE_RANGE_KPA[1]
        ):
            continue
        a, b = sorted(m.smiles)
        x = m.mole_fractions[m.smiles.index(a)]
        if 0.0 < x < 1.0:
            mixtures[(a, b, temperature)].append(sel.thermoml_point(m, x))

    rows = []
    for (a, b, temperature), points in mixtures.items():
        hmix = list(
            {round(p.x, 3): p for p in points if p.property_name == "hmix"}.values()
        )
        density = list(
            {round(p.x, 3): p for p in points if p.property_name == "density"}.values()
        )
        if len(hmix) < 3 or len(density) < 3:
            continue
        if not (component_passes(a) and component_passes(b)):
            continue
        info_a, info_b = sel.molecule_info(a), sel.molecule_info(b)
        hmix_eq = min(hmix, key=lambda p: abs(p.x - 0.5))
        row = {
            "a": a,
            "b": b,
            "temperature": temperature,
            "groups": "; ".join(
                sorted({g for c in (a, b) for g in sel.functional_groups(c, groups)})
            ),
            "hmix_eq": round(hmix_eq.value, 4),
            "x_eq": round(hmix_eq.x, 3),
            "n_hmix": len(hmix),
            "n_density": len(density),
            "new_pure_boxes": sum(
                (c, temperature) not in selected_pure for c in (a, b)
            ),
            "heavy_atoms": f"{info_a['heavy_atoms']}/{info_b['heavy_atoms']}",
            "dois": " ".join(sorted({p.doi for p in hmix + density})),
        }
        failures = []
        if not all(p.sigma for p in hmix + density):
            failures.append("missing experimental uncertainty")
        elif max(p.sigma / abs(p.value) for p in hmix if p.value) > (
            sel.MAX_RELATIVE_SIGMA_HMIX
        ):
            failures.append("relative ΔHmix uncertainty > 15 %")
        if not all(pure.get((c, temperature)) for c in (a, b)):
            failures.append("missing pure density with uncertainty")
        if (
            len(sel.target_points(sorted(hmix, key=lambda p: p.x), (0.25, 0.5, 0.75)))
            < 3
        ):
            failures.append("no ΔHmix near x = 0.25, 0.5 and 0.75")
        if info_a["viscous"] or info_b["viscous"]:
            failures.append("viscous")
        if (info_a["acid"] and info_b["base"]) or (info_a["base"] and info_b["acid"]):
            failures.append("acid-base pair")
        if (a, b) in two_phase:
            failures.append("two-phase (literature)")
        if sel.two_phase_signature(hmix):
            failures.append("two-phase (linear H^E)")
        if hmix_eq.value > 0.5 * sel.GAS_CONSTANT_KJ_PER_MOL_K * temperature:
            failures.append("two-phase (H^E > RT/2)")
        if abs(hmix_eq.value) < sel.MIN_SIGNAL:
            failures.append("small signal")
        row["failures"] = "; ".join(failures)
        row["passes"] = not failures
        rows.append(row)
    return sorted(
        rows,
        key=lambda r: (
            not r["passes"],
            r["groups"],
            r["new_pure_boxes"],
            -abs(r["hmix_eq"]),
        ),
    )


def main() -> None:
    """Screen the archive and write ``archive_candidates.csv``."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--selection-dir", type=Path, default=Path("selection"))
    parser.add_argument("--cache-dir", type=Path, default=sel.CACHE_DIR)
    args = parser.parse_args()

    selection = yaml.safe_load((sel.HERE / "selection.yml").read_text())
    groups = selection["functional_groups"]
    two_phase = {
        tuple(sorted((sel.canonical_smiles(s["a"]), sel.canonical_smiles(s["b"]))))
        for s in selection["two_phase_fast_set"]
    }
    with (args.selection_dir / "pure.csv").open() as handle:
        selected_pure = {
            (row["smiles"], round(float(row["temperature"]), 2))
            for row in csv.DictReader(handle)
        }
    archive = download(ARCHIVE_URL, args.cache_dir / Path(ARCHIVE_URL).name)
    rows = screen(
        archive_measurements(archive, groups), groups, selected_pure, two_phase
    )
    sel.write_csv(args.selection_dir / "archive_candidates.csv", rows)
    print(f"{sum(r['passes'] for r in rows)} of {len(rows)} archive systems pass")


if __name__ == "__main__":
    main()
