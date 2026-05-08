"""cancer_discovery_os — top-level integration of Yu-Cook lab work + Biomni.

Combines:
  - v18 framework (orchestration patterns, six-component spine)
  - research_compass (CV → topic recommendation, 5-stage pipeline)
  - grant_arena (competitive multi-agent grant-topic discovery + Elo)
  - frontier_scout (top-researcher paper monitoring)
  - the lab's own assets (Paper 1-3 + Paper 11, neoantigen hub, K2/BABA/CTC,
    v17 thyroid sprint discoveries, K-Thyro Foundation grant design)
  - Biomni (Stanford SNAP general-purpose biomedical AI agent) — adapter stub

Design contract:
  - Offline-first. No LLM calls, no Biomni install required, no network.
  - Manifest-driven: every work asset is a Manifest entry; the runner
    dispatches a natural-language query to the right sub-system based on
    intent classification (rule-based v0; LLM-classifier v0.3).
  - Marathon-mode safe: scaffolding/infra only. Live LLM + Biomni wiring
    deferred to post-2026-06-13.
"""
from __future__ import annotations

from .manifest import (
    LabManifest, ManifestEntry, EntryKind, default_yu_cook_manifest,
)
from .biomni_adapter import BiomniAdapter, BiomniCapability, default_biomni_capabilities
from .runner import CancerDiscoveryOS, OSQuery, OSResponse, IntentTag

__all__ = [
    "LabManifest", "ManifestEntry", "EntryKind", "default_yu_cook_manifest",
    "BiomniAdapter", "BiomniCapability", "default_biomni_capabilities",
    "CancerDiscoveryOS", "OSQuery", "OSResponse", "IntentTag",
]
