"""Stage 4 — open-data metadata registry.

Queries dataset metadata (NOT the data itself) across:
  - GEO       NCBI Entrez esearch/esummary on 'gds'
  - SRA       NCBI Entrez esearch/esummary on 'sra'
  - cBioPortal /api/studies
  - UK Biobank Showcase (field metadata only — no participant data)
  - GTEx, CellxGene, Open Targets — placeholders

Every fetcher accepts a ``live`` flag. Default offline returns stub entries.
"""
from __future__ import annotations

import asyncio
from dataclasses import dataclass, field
from typing import Sequence


@dataclass
class DatasetEntry:
    accession: str          # e.g. "GSE286332", "PRJEB11591", "thca_tcga_pub"
    source: str             # "GEO" | "SRA" | "cBioPortal" | "UKB" | ...
    title: str = ""
    organism: str = "Homo sapiens"
    n_samples: int = 0
    assay: str = ""         # "RNA-seq", "scRNA-seq", "WGS", ...
    keywords: list[str] = field(default_factory=list)
    url: str = ""
    access: str = "open"    # "open" | "controlled" | "DUA"

    def all_terms(self) -> set[str]:
        return {t.lower() for t in (self.keywords + [self.assay]) if t}


async def _query_geo(terms: Sequence[str]) -> list[DatasetEntry]:
    """Live GEO esearch + esummary. Stub here; wire NCBI Entrez when going live."""
    return []


async def _query_sra(terms: Sequence[str]) -> list[DatasetEntry]:
    return []


async def _query_cbioportal(terms: Sequence[str]) -> list[DatasetEntry]:
    """https://www.cbioportal.org/api/studies — keyword filter client-side."""
    return []


async def _query_ukb(terms: Sequence[str]) -> list[DatasetEntry]:
    """UKB Showcase: field metadata is open; participant data is DUA-gated."""
    return []


async def query_datasets(
    terms: Sequence[str],
    sources: Sequence[str] = ("GEO", "SRA", "cBioPortal", "UKB"),
    *,
    live: bool = False,
    stub_entries: Sequence[DatasetEntry] | None = None,
) -> list[DatasetEntry]:
    """Tier-2 fan-out across sources. Returns flat DatasetEntry list."""
    if not live:
        return list(stub_entries or [])

    coros = []
    if "GEO" in sources:
        coros.append(_query_geo(terms))
    if "SRA" in sources:
        coros.append(_query_sra(terms))
    if "cBioPortal" in sources:
        coros.append(_query_cbioportal(terms))
    if "UKB" in sources:
        coros.append(_query_ukb(terms))

    results = await asyncio.gather(*coros, return_exceptions=True)
    out: list[DatasetEntry] = []
    for r in results:
        if isinstance(r, Exception):
            continue
        out.extend(r)
    return out
