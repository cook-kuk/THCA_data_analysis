"""patent_priorart_mcp — prior-art collision risk screening.

IMPORTANT: this is *research planning* prior-art screening, NOT a legal
freedom-to-operate (FTO) opinion. Final FTO requires patent-attorney review.
The MCP only surfaces risk signals so research-topic claims can be
differentiated before submission.

Offline stub matches against a small in-memory corpus of representative
prior-art entries; v0.3 wires Google Patents / KIPRIS / WIPO PATENTSCOPE.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence


@dataclass
class PriorArtHit:
    accession: str          # e.g. "US 11,xxx,xxx" or "KR-10-xxxx-xxxx"
    title: str
    assignee: str
    similarity: float = 0.0  # 0..1, surface keyword overlap
    notes: str = ""


_REP_PRIOR_ART = [
    PriorArtHit("US-2023-DT-ONC", "Oncology digital twin platform",
                "Anonymous Pharma", 0.0, "treatment response simulation"),
    PriorArtHit("US-2024-AI-COSCI", "AI co-scientist for hypothesis generation",
                "Big Tech Co", 0.0, "multi-agent research planning"),
    PriorArtHit("KR-10-XXXX-1", "수술영상 phase recognition AI",
                "국내 의료AI 기업", 0.0, "surgical phase recognition"),
]


def screen_collision_risk(
    topic_title: str,
    topic_keywords: Sequence[str],
    *,
    live: bool = False,
    stub_corpus: Sequence[PriorArtHit] | None = None,
) -> list[PriorArtHit]:
    """Return prior-art hits ranked by surface keyword overlap.

    NOT a legal FTO. For research planning only.
    """
    corpus = list(stub_corpus or _REP_PRIOR_ART)
    if live:
        raise NotImplementedError("Live patent search deferred to v0.3")

    q = {k.lower() for k in topic_keywords if k}
    q.update(t.lower() for t in topic_title.split())
    hits: list[PriorArtHit] = []
    for entry in corpus:
        ent_terms = {t.lower() for t in (entry.title + " " + entry.notes).split() if t}
        overlap = len(q & ent_terms) / max(1, len(q | ent_terms))
        if overlap > 0:
            hits.append(PriorArtHit(
                accession=entry.accession, title=entry.title,
                assignee=entry.assignee, similarity=round(overlap, 3),
                notes=entry.notes,
            ))
    hits.sort(key=lambda h: h.similarity, reverse=True)
    return hits
