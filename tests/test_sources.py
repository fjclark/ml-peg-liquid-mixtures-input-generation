"""Test the data-source citation helpers (offline)."""

from __future__ import annotations

import pytest

from liquid_mixtures_inputs import sources

METADATA = {
    "DOI": "10.1021/je900874z",
    "author": [
        {"given": "Iván", "family": "Alonso"},
        {"given": "Víctor", "family": "Alonso"},
        {"given": "Isaías García", "family": "de la Fuente"},
    ],
    "title": "Ketone + Amine Mixtures for <i>N-</i>Methylaniline",
    "container-title": "Journal of Chemical &amp; Engineering Data",
    "container-title-short": "J. Chem. Eng. Data",
    "volume": "55",
    "issue": "7",
    "page": "2505-2511",
    "issued": {"date-parts": [[2009, 12, 1]]},
    "published-print": {"date-parts": [[2010, 7, 8]]},
}

SPECS = [
    {
        "id": "acetone_nma",
        "a": "CC(C)=O",
        "b": "CNc1ccccc1",
        "names": ["acetone", "N-methylaniline"],
    },
]


def test_short_citation():
    """Test the one-line citation uses et al. and the print year."""
    assert (
        sources.short_citation(METADATA)
        == "Alonso et al. (2010), J. Chem. Eng. Data 55, 2505-2511"
    )


def test_source_uses_lists_properties_and_pure_densities():
    """Test every DOI is listed with the properties and pure densities it gives."""
    isotherms = [
        {"system_id": "acetone_nma", "prop": "hmix", "doi": "10.1000/h"},
        {"system_id": "acetone_nma", "prop": "density", "doi": "10.1000/d"},
        {"system_id": "acetone_nma", "prop": "density", "doi": "10.1000/d"},
        {"system_id": "acetone_nma", "prop": "vex", "doi": "10.1000/d"},
    ]
    pure = [
        {"smiles": "CC(C)=O", "temperature": 298.15, "source": "thermoml:10.1000/d"},
        {"smiles": "O", "temperature": 298.15, "source": "iapws-95:10.1000/w"},
    ]
    assert sources.source_uses(isotherms, pure, SPECS) == {
        "10.1000/h": ["acetone + N-methylaniline (ΔHmix)"],
        "10.1000/d": [
            "acetone + N-methylaniline (density, V^E)",
            "pure acetone density (298.15 K)",
        ],
        "10.1000/w": ["pure O density (298.15 K)"],
    }


def test_source_uses_rejects_a_pure_source_without_doi():
    """Test a pure density that cannot be cited stops the selection."""
    pure = [{"smiles": "O", "temperature": 298.15, "source": "ml-peg"}]
    with pytest.raises(ValueError, match="without a DOI"):
        sources.source_uses([], pure, SPECS)
