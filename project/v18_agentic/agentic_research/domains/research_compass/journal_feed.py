"""Stage 2 — recent-paper feed.

Pulls the last N days of papers from:
  - BRIC 한빛사 (한국 연구자 highlight feed)
  - bioRxiv API (search-based)
  - Top-tier journal RSS (Nature, Cell, NEJM, Science) — optional

All fetchers are async and respect a ``live`` flag. Default offline mode
returns whatever is in ``stub_items`` so tests/demos work without network.
"""
from __future__ import annotations

import asyncio
from dataclasses import dataclass, field
from datetime import datetime
from typing import Sequence


@dataclass
class JournalItem:
    title: str
    journal: str           # e.g. "Nature", "BRIC-한빛사", "bioRxiv"
    date: str              # ISO yyyy-mm-dd
    abstract: str = ""
    url: str = ""
    keywords: list[str] = field(default_factory=list)
    source: str = ""       # "bric" | "biorxiv" | "rss"


_DEFAULT_TOP_JOURNALS = ("Nature", "Cell", "NEJM", "Science", "Lancet", "Nat Med", "Cancer Cell", "Cell Rep Med", "JCI Insight")


async def _fetch_bric_hanbitsa(days: int) -> list[JournalItem]:
    """Live BRIC 한빛사 fetcher. Stubbed in offline mode — wire HTTP here when going live."""
    # Real implementation: GET https://www.ibric.org/myboard/list.php?Board=hanbitsa
    # then parse the listing. Left as a TODO behind the live flag.
    return []


async def _fetch_biorxiv(query: str, days: int) -> list[JournalItem]:
    """Live bioRxiv API fetcher.

    bioRxiv exposes /details/biorxiv/{from}/{to} for date-window listings;
    keyword filtering happens client-side. Wire when live.
    """
    return []


async def _fetch_top_journal_rss(journals: Sequence[str], days: int) -> list[JournalItem]:
    """Live RSS fan-out. Wire feedparser or httpx when live."""
    return []


async def fetch_journal_feed(
    days: int = 90,
    query: str = "",
    journals: Sequence[str] = _DEFAULT_TOP_JOURNALS,
    *,
    live: bool = False,
    stub_items: Sequence[JournalItem] | None = None,
) -> list[JournalItem]:
    """Tier-2 fan-out across BRIC + bioRxiv + RSS, returning a flat list.

    Parameters
    ----------
    days : look-back window
    query : free-text query for the search-based sources (bioRxiv)
    journals : RSS feed sources to include
    live : if False, returns stub_items only (test/demo path)
    stub_items : injected items for offline runs
    """
    if not live:
        return list(stub_items or [])

    results = await asyncio.gather(
        _fetch_bric_hanbitsa(days),
        _fetch_biorxiv(query, days),
        _fetch_top_journal_rss(journals, days),
        return_exceptions=True,
    )
    items: list[JournalItem] = []
    for r in results:
        if isinstance(r, Exception):
            continue
        items.extend(r)
    items.sort(key=lambda x: x.date, reverse=True)
    return items


def filter_recent(items: Sequence[JournalItem], cutoff_iso: str) -> list[JournalItem]:
    """Keep items dated >= cutoff_iso. Items without a parseable date are dropped."""
    keep: list[JournalItem] = []
    for it in items:
        try:
            if datetime.fromisoformat(it.date) >= datetime.fromisoformat(cutoff_iso):
                keep.append(it)
        except ValueError:
            continue
    return keep
