"""Offline demo for cancer_discovery_os.

Shows the OS handling 5 different intents — asset lookup, wet-lab protocol
(Biomni stub), grant discovery (arena), topic recommendation (compass),
and unknown (overview). No network, no LLM, no Biomni install.

Run:
    cd project/v18_agentic && python3 -m examples.cancer_discovery_os_demo
"""
from __future__ import annotations

import asyncio

from agentic_research.domains.research_compass import (
    DatasetEntry, MethodIndex, ResearchProfile,
)
from agentic_research.domains.research_compass.method_index import default_seed_entries
from agentic_research.domains.research_compass.grant_arena.mcp_servers.bric_trend import BRICItem
from agentic_research.domains.cancer_discovery_os import (
    CancerDiscoveryOS, OSQuery, IntentTag,
)


def main():
    os_ = CancerDiscoveryOS()
    print(f"Manifest: {len(os_.manifest)} entries (PI={os_.manifest.pi}, Co-PI={os_.manifest.co_pi})")
    print(f"Biomni capabilities: {len(os_.biomni.capabilities)} (available={os_.biomni.available})")

    print("\n--- 1) ASSET_LOOKUP: 'PTC' ---")
    r = os_.handle(OSQuery("what do we have on PTC?"))
    print(f"  intent={r.intent.value} handler={r.handler} hits={len(r.result)}")
    for e in r.result[:3]:
        print(f"    [{e.kind.value}] {e.id} — {e.title}")

    print("\n--- 2) WET_LAB_PROTOCOL: CRISPR screen ---")
    r = os_.handle(OSQuery("plan a CRISPR screen for thyroid dedifferentiation"))
    print(f"  intent={r.intent.value} handler={r.handler}")
    print(f"  result={r.result}")
    print(f"  notes={r.notes}")

    print("\n--- 3) GRANT_DISCOVERY ---")
    profile = ResearchProfile(
        name="Cook", skills=["agentic AI", "spatial transcriptomics"],
        prior_topics=["thyroid cancer", "active surveillance"],
        methods_used=["target trial emulation", "ComBat-seq"],
    )
    datasets = [
        DatasetEntry(accession="TCGA-THCA", source="GDC", assay="RNA-seq",
                     keywords=["thyroid", "PTC"], access="open"),
        DatasetEntry(accession="GSE286332", source="GEO", assay="RNA-seq",
                     keywords=["thyroid", "Korean"], access="open"),
    ]
    method_index = MethodIndex(default_seed_entries())
    bric = [BRICItem(title="MAeSTro-EXP outcomes", korean_pi="Yu HW",
                     affiliation="SNUBH", journal="Thyroid", year=2025,
                     keywords=["active surveillance", "Korean"])]
    r = os_.handle(OSQuery(
        "삼성 30억 grant topic 발굴",
        payload={
            "profile": profile, "datasets": datasets,
            "method_index": method_index, "bric_items": bric,
            "rfp_text": "Advanced AI, 디지털 헬스, 정밀의료, 30억 3년",
        },
    ))
    arena_res = r.result
    print(f"  intent={r.intent.value} handler={r.handler}")
    print(f"  topics={len(arena_res.topics)} battles={len(arena_res.battle_log)}")
    print("  top 3:")
    for t in arena_res.ranking[:3]:
        sc = arena_res.final_scores[t.title]
        print(f"    [{sc.tier} {sc.total}] {t.title}  by={t.proposed_by}")

    print("\n--- 4) TOPIC_RECOMMEND (no CV → notes only) ---")
    r = os_.handle(OSQuery("recommend a research direction"))
    print(f"  intent={r.intent.value} handler={r.handler} notes={r.notes}")

    print("\n--- 5) UNKNOWN → OS overview ---")
    r = os_.handle(OSQuery("hello world"))
    print(f"  intent={r.intent.value}")
    print(f"  result={r.result}")


if __name__ == "__main__":
    main()
