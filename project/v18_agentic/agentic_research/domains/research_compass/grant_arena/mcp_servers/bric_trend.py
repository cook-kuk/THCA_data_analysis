"""bric_trend_mcp — Korean biomedical trend radar (BRIC 한빛사 + Bio뉴스).

BRIC 한빛사 is a curated channel of Korean researchers' high-impact biomed
papers. Offline mode returns injected stubs; v0.3 will scrape the listing
page (static HTML, low-friction).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence


@dataclass
class BRICItem:
    title: str
    korean_pi: str
    affiliation: str
    journal: str
    year: int
    keywords: list[str] = field(default_factory=list)
    url: str = ""


async def search_bric_hanbitsa(
    query: str,
    *,
    live: bool = False,
    stub_items: Sequence[BRICItem] | None = None,
) -> list[BRICItem]:
    if not live:
        return list(stub_items or [])
    raise NotImplementedError("Live BRIC scraper deferred to v0.3")
