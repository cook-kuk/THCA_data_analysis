"""Stage 3 — keyword-indexed method paper library.

A flat in-memory index over method papers. Each entry has a canonical name
('ComBat-seq', 'CellRanger', 'Scanpy', 'cookHLA', 'arcasHLA', ...) plus
keywords + the application domains it touches.

For v0 we use bag-of-keywords + Jaccard scoring. The interface is stable so
v1 can swap in sentence-transformer embeddings without touching callers.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Sequence


@dataclass
class MethodEntry:
    name: str
    paper_title: str = ""
    year: int = 0
    domains: list[str] = field(default_factory=list)        # e.g. "RNA-seq", "scRNA-seq", "HLA imputation"
    keywords: list[str] = field(default_factory=list)       # technique terms
    code_url: str = ""

    def all_terms(self) -> set[str]:
        return {t.lower() for t in (self.domains + self.keywords + [self.name]) if t}


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


class MethodIndex:
    """Lookup methods relevant to a query (e.g. researcher's keywords)."""

    def __init__(self, entries: Sequence[MethodEntry] | None = None):
        self._entries: list[MethodEntry] = list(entries or [])

    def add(self, entry: MethodEntry) -> None:
        self._entries.append(entry)

    def extend(self, entries: Iterable[MethodEntry]) -> None:
        self._entries.extend(entries)

    def __len__(self) -> int:
        return len(self._entries)

    def search(self, terms: Iterable[str], top_k: int = 10) -> list[tuple[MethodEntry, float]]:
        """Return up to top_k (entry, score) tuples sorted by Jaccard overlap."""
        q = {t.lower() for t in terms if t}
        if not q:
            return []
        scored = [(e, _jaccard(q, e.all_terms())) for e in self._entries]
        scored = [(e, s) for e, s in scored if s > 0]
        scored.sort(key=lambda x: x[1], reverse=True)
        return scored[:top_k]


def default_seed_entries() -> list[MethodEntry]:
    """A tiny seed registry — enough to demo. Real deployments load from JSON/CSV."""
    return [
        MethodEntry(name="ComBat-seq", paper_title="Zhang 2020", year=2020,
                    domains=["RNA-seq", "batch correction"],
                    keywords=["batch effect", "negative binomial", "LODO"]),
        MethodEntry(name="arcasHLA", paper_title="Orenbuch 2020", year=2020,
                    domains=["HLA imputation", "RNA-seq"],
                    keywords=["HLA typing", "class I", "class II", "MHC"]),
        MethodEntry(name="CIBERSORTx", paper_title="Newman 2019", year=2019,
                    domains=["deconvolution", "bulk RNA-seq"],
                    keywords=["cell type fractions", "immune", "signature matrix"]),
        MethodEntry(name="CLAM", paper_title="Lu 2021", year=2021,
                    domains=["pathology", "weakly supervised"],
                    keywords=["WSI", "attention MIL", "H&E"]),
        MethodEntry(name="DESeq2", paper_title="Love 2014", year=2014,
                    domains=["RNA-seq", "differential expression"],
                    keywords=["count model", "shrinkage", "GLM"]),
    ]
