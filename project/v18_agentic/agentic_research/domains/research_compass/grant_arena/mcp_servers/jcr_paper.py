"""jcr_paper_mcp — JCR Q1 / top-tier biomed filter (offline stub).

User preference (2026-05-08): publication target floor = Nature Communications
or higher. Sub-NC venues are filtered OUT, not reported as Q1. Replace with
a live JCR/InCites fetcher in v0.3.

The authoritative source is Clarivate's JCR. Offline we hand-curate a small
allowlist of top-tier biomedical journals so demos and tests work without API.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class JCRRecord:
    journal: str
    quartile: str           # "Q1" .. "Q4"
    category_top_pct: float # 0..1, lower = better (e.g. 0.05 = top 5%)
    tier: str = "Q1"        # "TOP" (NC and above) | "Q1" | "below_NC"


_TOP_TIER_BIOMED = {
    "nature", "nature medicine", "nature cancer", "nature communications",
    "nature methods", "nature biotechnology", "nature genetics", "nature reviews",
    "cell", "cell reports medicine", "cancer cell", "cancer discovery",
    "lancet", "lancet oncology", "nejm", "new england journal of medicine",
    "science", "science advances", "science translational medicine",
    "jama", "jama oncology",
}


_BELOW_NC_BUT_Q1 = {
    "jci insight", "npj digital medicine", "npj quantum information",
    "bioinformatics", "genome biology", "genome research",
    "scientific reports", "plos one", "frontiers",
}


def is_jcr_q1(journal: str, *, min_tier: str = "TOP") -> JCRRecord | None:
    """Return a JCRRecord if the journal meets the floor.

    min_tier="TOP" (default per user preference): only Nature Communications
    and above qualify. min_tier="Q1": include below-NC Q1 venues too.
    """
    j = journal.strip().lower()
    if j in _TOP_TIER_BIOMED:
        return JCRRecord(journal=journal, quartile="Q1",
                         category_top_pct=0.05, tier="TOP")
    if min_tier == "Q1" and j in _BELOW_NC_BUT_Q1:
        return JCRRecord(journal=journal, quartile="Q1",
                         category_top_pct=0.10, tier="Q1")
    return None


def journal_tier(journal: str) -> str:
    """Return 'TOP' | 'Q1' | 'below_Q1'."""
    j = journal.strip().lower()
    if j in _TOP_TIER_BIOMED:
        return "TOP"
    if j in _BELOW_NC_BUT_Q1:
        return "Q1"
    return "below_Q1"
