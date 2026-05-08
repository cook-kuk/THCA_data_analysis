"""Smoke tests for cancer_discovery_os."""
from __future__ import annotations

from agentic_research.domains.research_compass import (
    DatasetEntry, MethodIndex, ResearchProfile,
)
from agentic_research.domains.research_compass.method_index import default_seed_entries
from agentic_research.domains.research_compass.grant_arena.mcp_servers.bric_trend import BRICItem
from agentic_research.domains.cancer_discovery_os import (
    BiomniAdapter, CancerDiscoveryOS, EntryKind, IntentTag, LabManifest,
    OSQuery, default_biomni_capabilities, default_yu_cook_manifest,
)


def test_default_manifest_has_core_entries():
    m = default_yu_cook_manifest()
    assert m.pi and m.co_pi
    ids = {e.id for e in m.entries}
    for required in ("paper1_dm1", "paper2_image_dm1", "paper11_pancancer",
                     "cohort_korean_k2", "cohort_baba", "cohort_ctc_emt",
                     "platform_v18", "platform_v19_grant_arena", "grant_causal_kthyro"):
        assert required in ids, f"missing manifest entry: {required}"


def test_manifest_search_and_by_kind():
    m = default_yu_cook_manifest()
    grants = m.by_kind(EntryKind.GRANT)
    assert grants and grants[0].id == "grant_causal_kthyro"
    thyroid_hits = m.search("thyroid")
    assert any("DM1" in e.title or "PTC" in e.title or "thyroid" in e.title.lower() for e in thyroid_hits)


def test_biomni_adapter_stubs_responses():
    a = BiomniAdapter()
    assert not a.available
    caps = a.list_capabilities()
    assert any(c.name == "crispr_screen_plan" for c in caps)
    out = a.go("plan a CRISPR screen for thyroid dedifferentiation")
    assert out["status"] == "stubbed"
    assert out["matched_capability"] == "crispr_screen_plan"


def test_biomni_adapter_falls_back_to_lab_qa():
    out = BiomniAdapter().go("what is the Hashimoto-overlap effect on PTC?")
    assert out["matched_capability"] == "lab_qa"


def test_os_intent_classifies_asset_lookup():
    os_ = CancerDiscoveryOS()
    r = os_.handle(OSQuery("what assets do we have on PTC?"))
    assert r.intent is IntentTag.ASSET_LOOKUP
    assert r.handler == "LabManifest.search"
    assert isinstance(r.result, list) and len(r.result) > 0


def test_os_intent_classifies_wet_lab():
    r = CancerDiscoveryOS().handle(OSQuery("plan a CRISPR screen"))
    assert r.intent is IntentTag.WET_LAB_PROTOCOL
    assert r.result["status"] == "stubbed"


def test_os_grant_discovery_runs_arena():
    profile = ResearchProfile(
        name="t", skills=["agentic AI"],
        prior_topics=["thyroid cancer"],
        methods_used=["target trial emulation"],
    )
    datasets = [DatasetEntry(accession="TCGA-THCA", source="GDC", assay="RNA-seq",
                             keywords=["thyroid"], access="open")]
    bric = [BRICItem(title="K", korean_pi="Y", affiliation="S", journal="J", year=2025,
                     keywords=["thyroid"])]
    r = CancerDiscoveryOS().handle(OSQuery(
        "삼성 grant 발굴",
        payload={
            "profile": profile, "datasets": datasets,
            "method_index": MethodIndex(default_seed_entries()),
            "bric_items": bric,
            "rfp_text": "Advanced AI, 디지털 헬스, 30억",
        },
    ))
    assert r.intent is IntentTag.GRANT_DISCOVERY
    assert r.result is not None
    assert r.result.topics


def test_os_unknown_returns_overview():
    r = CancerDiscoveryOS().handle(OSQuery("hello world"))
    assert r.intent is IntentTag.UNKNOWN
    assert "subsystems" in r.result
    assert r.result["biomni_available"] is False


def test_default_biomni_capabilities_complete():
    caps = default_biomni_capabilities()
    names = {c.name for c in caps}
    for required in ("crispr_screen_plan", "scrna_annotation",
                     "admet_prediction", "gwas_causal_gene", "rare_disease_diagnosis"):
        assert required in names
