"""
Write citations for every experimental data source.

Each DOI behind the reference data (the isotherm points of ``isotherms.csv``
and the pure densities of ``pure.csv``) is resolved at doi.org as CSL-JSON
(cached), and written to ``dois.txt`` (one DOI per line, to copy) and
``sources.csv`` (one row per paper, with what it supplies).
"""

from __future__ import annotations

from collections import defaultdict
import csv
import html
import json
from pathlib import Path
import re

from liquid_mixtures_inputs.thermoml import download

DOI_URL = "https://doi.org/{doi}"
CSL_JSON = "application/vnd.citationstyles.csl+json"
DOI_PATTERN = re.compile(r"10\.\d{4,9}/\S+")
PROPERTY_LABELS = {"hmix": "ΔHmix", "density": "density", "vex": "V^E"}


def fetch_metadata(doi: str, cache_dir: Path) -> dict:
    """
    Download (or read from cache) the CSL-JSON metadata of a DOI.

    Parameters
    ----------
    doi
        DOI, e.g. ``10.1021/je301238f``.
    cache_dir
        Directory used to cache downloaded records.

    Returns
    -------
    dict
        CSL-JSON metadata.
    """
    path = cache_dir / "doi" / f"{doi.replace('/', '_')}.json"
    url = DOI_URL.format(doi=doi)
    return json.loads(download(url, path, {"Accept": CSL_JSON}).read_text())


def clean(text: str) -> str:
    """
    Strip HTML tags and entities and collapse whitespace.

    Parameters
    ----------
    text
        Text from a CSL-JSON field, e.g. a title with ``<i>`` tags.

    Returns
    -------
    str
        Plain text.
    """
    return " ".join(html.unescape(re.sub(r"<[^>]+>", "", text)).split())


def year(metadata: dict) -> str:
    """
    Return the publication year of a CSL-JSON record.

    The print year comes first, as journals cite the volume's year; ``issued``
    is the earliest date, often the online one of the year before.

    Parameters
    ----------
    metadata
        CSL-JSON metadata.

    Returns
    -------
    str
        Year, or an empty string if none is given.
    """
    for field in ("published-print", "issued"):
        parts = metadata.get(field, {}).get("date-parts", [[None]])[0]
        if parts and parts[0]:
            return str(parts[0])
    return ""


def short_citation(metadata: dict) -> str:
    """
    Format a one-line citation.

    Parameters
    ----------
    metadata
        CSL-JSON metadata.

    Returns
    -------
    str
        E.g. ``Alonso et al. (2010), J. Chem. Eng. Data 55, 2505-2511``.
    """
    authors = metadata.get("author", [])
    names = [a.get("family") or a.get("literal", "") for a in authors]
    if len(names) > 2:
        who = f"{names[0]} et al."
    else:
        who = " and ".join(names)
    journal = clean(
        metadata.get("container-title-short") or metadata.get("container-title", "")
    )
    where = " ".join(part for part in (journal, metadata.get("volume", "")) if part)
    page = metadata.get("page", "")
    return f"{who} ({year(metadata)}), {where}" + (f", {page}" if page else "")


def source_uses(
    isotherms: list[dict], pure: list[dict], specs: list[dict]
) -> dict[str, list[str]]:
    """
    List what each DOI supplies to the benchmark.

    Parameters
    ----------
    isotherms
        Rows of ``isotherms.csv``.
    pure
        Rows of ``pure.csv``.
    specs
        Curated systems from ``selection.yml``.

    Returns
    -------
    dict[str, list[str]]
        Uses per DOI, e.g. ``["pyrrole + pyridine (ΔHmix, density)"]``.

    Raises
    ------
    ValueError
        If a pure density has a source without a DOI.
    """
    systems = {spec["id"]: " + ".join(spec["names"]) for spec in specs}
    molecules = {
        smiles: name
        for spec in specs
        for smiles, name in zip((spec["a"], spec["b"]), spec["names"], strict=True)
    }
    properties = defaultdict(lambda: defaultdict(list))
    for row in isotherms:
        labels = properties[row["doi"]][systems[row["system_id"]]]
        if PROPERTY_LABELS[row["prop"]] not in labels:
            labels.append(PROPERTY_LABELS[row["prop"]])
    uses = {
        doi: [f"{system} ({', '.join(labels)})" for system, labels in by_system.items()]
        for doi, by_system in properties.items()
    }
    for row in pure:
        if not row["source"]:
            continue
        match = DOI_PATTERN.search(row["source"])
        if match is None:
            raise ValueError(f"pure density source without a DOI: {row['source']}")
        name = molecules.get(row["smiles"], row["smiles"])
        uses.setdefault(match.group(), []).append(
            f"pure {name} density ({float(row['temperature']):g} K)"
        )
    return uses


def write_sources(
    out: Path,
    isotherms: list[dict],
    pure: list[dict],
    specs: list[dict],
    cache_dir: Path,
) -> None:
    """
    Write ``dois.txt`` and ``sources.csv`` for every data source.

    Parameters
    ----------
    out
        Output directory.
    isotherms
        Rows of ``isotherms.csv``.
    pure
        Rows of ``pure.csv``.
    specs
        Curated systems from ``selection.yml``.
    cache_dir
        Directory used to cache downloaded records.
    """
    uses = source_uses(isotherms, pure, specs)
    metadata = {doi: fetch_metadata(doi, cache_dir) for doi in uses}
    citations = {doi: short_citation(metadata[doi]) for doi in uses}
    ordered = sorted(uses, key=lambda doi: (citations[doi], doi))
    with (out / "sources.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["doi", "url", "citation", "title", "uses"])
        for doi in ordered:
            writer.writerow(
                [
                    doi,
                    DOI_URL.format(doi=doi),
                    citations[doi],
                    clean(metadata[doi].get("title", "")),
                    "; ".join(uses[doi]),
                ]
            )
    (out / "dois.txt").write_text("".join(f"{doi}\n" for doi in ordered))
