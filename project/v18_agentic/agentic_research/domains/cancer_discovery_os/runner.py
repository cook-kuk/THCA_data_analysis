"""CancerDiscoveryOS — top-level runner that dispatches a query to subsystems.

Decision tree (rule-based v0; LLM-classified v0.3):

    intent: GRANT_DISCOVERY    → research_compass.grant_arena.GrantArena.run(...)
    intent: TOPIC_RECOMMEND    → research_compass.ResearchCompass.run(cv)
    intent: WET_LAB_PROTOCOL   → BiomniAdapter.go(task, capability=...)
    intent: ASSET_LOOKUP       → LabManifest.search(query)
    intent: FRONTIER_SCAN      → grant_arena.frontier_scout.scout_watchlist(...)
    intent: UNKNOWN            → return list of available subsystems + manifest summary

Pure dispatcher. No LLM call. The OS object owns pre-built compass +
arena instances so callers don't re-initialize state per query.
"""
from __future__ import annotations

import asyncio
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from ..research_compass import ResearchCompass
from ..research_compass.grant_arena import GrantArena, parse_rfp
from .biomni_adapter import BiomniAdapter
from .manifest import LabManifest, default_yu_cook_manifest


class IntentTag(str, Enum):
    GRANT_DISCOVERY = "grant_discovery"
    TOPIC_RECOMMEND = "topic_recommend"
    WET_LAB_PROTOCOL = "wet_lab_protocol"
    ASSET_LOOKUP = "asset_lookup"
    FRONTIER_SCAN = "frontier_scan"
    UNKNOWN = "unknown"


@dataclass
class OSQuery:
    text: str
    intent: IntentTag | None = None          # None → auto-classify
    payload: dict[str, Any] = field(default_factory=dict)


@dataclass
class OSResponse:
    intent: IntentTag
    handler: str
    result: Any = None
    notes: list[str] = field(default_factory=list)


_INTENT_KEYWORDS = {
    IntentTag.GRANT_DISCOVERY: ("grant", "rfp", "과제", "제안서", "samsung", "삼성"),
    IntentTag.TOPIC_RECOMMEND: ("topic", "research direction", "what should i work on", "주제 추천"),
    IntentTag.WET_LAB_PROTOCOL: ("crispr", "scrna", "single cell", "admet", "gwas",
                                  "wet lab", "experiment", "protocol", "perturbation"),
    IntentTag.ASSET_LOOKUP: ("asset", "manifest", "what do we have", "lab inventory",
                              "내 자산", "내 데이터"),
    IntentTag.FRONTIER_SCAN: ("watchlist", "what's hot", "frontier", "trending",
                                "최신 트렌드", "hot methods"),
}


def _classify(text: str) -> IntentTag:
    t = text.lower()
    for tag, kws in _INTENT_KEYWORDS.items():
        if any(kw in t for kw in kws):
            return tag
    return IntentTag.UNKNOWN


@dataclass
class CancerDiscoveryOS:
    manifest: LabManifest = field(default_factory=default_yu_cook_manifest)
    biomni: BiomniAdapter = field(default_factory=BiomniAdapter)
    compass: ResearchCompass | None = None
    arena: GrantArena | None = None

    def __post_init__(self):
        if self.compass is None:
            self.compass = ResearchCompass(live=False)
        if self.arena is None:
            self.arena = GrantArena(topics_per_generator=2, evolve_top_k=2)

    def handle(self, query: OSQuery) -> OSResponse:
        intent = query.intent or _classify(query.text)

        if intent is IntentTag.ASSET_LOOKUP:
            term = query.payload.get("search", query.text)
            hits = self.manifest.search(term) or self.manifest.entries
            return OSResponse(intent=intent, handler="LabManifest.search", result=hits)

        if intent is IntentTag.WET_LAB_PROTOCOL:
            return OSResponse(
                intent=intent, handler="BiomniAdapter.go",
                result=self.biomni.go(query.text, capability=query.payload.get("capability")),
                notes=[f"biomni.available={self.biomni.available}"],
            )

        if intent is IntentTag.TOPIC_RECOMMEND:
            cv = query.payload.get("cv", "")
            result = asyncio.get_event_loop().run_until_complete(self.compass.run(cv)) if cv else None
            return OSResponse(intent=intent, handler="ResearchCompass.run", result=result,
                              notes=[] if cv else ["no cv supplied; pass payload['cv']"])

        if intent is IntentTag.GRANT_DISCOVERY:
            needed = ("profile", "datasets", "method_index")
            missing = [k for k in needed if k not in query.payload]
            if missing:
                return OSResponse(
                    intent=intent, handler="GrantArena.run",
                    notes=[f"missing payload keys: {missing}; "
                           "provide profile/datasets/method_index/bric_items/rfp_text"],
                )
            rfp = query.payload.get("rfp") or parse_rfp(query.payload.get("rfp_text", ""))
            arena_res = self.arena.run(
                profile=query.payload["profile"],
                rfp=rfp,
                datasets=query.payload["datasets"],
                method_index=query.payload["method_index"],
                bric_items=query.payload.get("bric_items", ()),
            )
            return OSResponse(intent=intent, handler="GrantArena.run", result=arena_res)

        if intent is IntentTag.FRONTIER_SCAN:
            return OSResponse(
                intent=intent, handler="frontier_scout.scout_watchlist",
                notes=["live scout deferred to v0.3; pass payload['stub_by_name'] to run offline"],
            )

        return OSResponse(
            intent=IntentTag.UNKNOWN, handler="manifest_summary",
            result={
                "manifest_size": len(self.manifest),
                "subsystems": ["LabManifest", "ResearchCompass", "GrantArena",
                               "frontier_scout", "BiomniAdapter"],
                "biomni_available": self.biomni.available,
            },
            notes=["intent classification fell through; here is the OS overview"],
        )
