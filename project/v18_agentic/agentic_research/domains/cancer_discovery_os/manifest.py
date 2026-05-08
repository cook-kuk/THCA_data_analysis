"""LabManifest — registry of all Yu-Cook lab assets in one machine-readable place.

Why a manifest matters: the agent can only orchestrate what it knows about.
Encoding papers, datasets, cohorts, prior findings, and platform modules as
a typed registry lets every downstream agent reason over them without
re-deriving structure each time.

This is the source-of-truth registry for the OS. v0.3 will sync from
project/CLAUDE.md, MEMORY.md, and the papers_hub_2026_05_04 directory.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class EntryKind(str, Enum):
    PAPER = "paper"                  # in-progress or shipped manuscript
    COHORT = "cohort"                # patient/sample collection
    DATASET = "dataset"              # public or internal data asset
    METHOD = "method"                # analytical method / pipeline / signature
    PLATFORM = "platform"            # software module (v18, research_compass, ...)
    GRANT = "grant"                  # grant proposal or proposal seed
    DISCOVERY = "discovery"          # validated finding (DM1, NBNR, etc.)


@dataclass
class ManifestEntry:
    id: str
    kind: EntryKind
    title: str
    summary: str = ""
    tags: list[str] = field(default_factory=list)
    status: str = "active"           # "active" | "shipped" | "frozen" | "backlog"
    path: str = ""                   # repo-relative location, if any
    n_samples: int | None = None
    journal_target: str = ""
    notes: str = ""


@dataclass
class LabManifest:
    pi: str = ""
    co_pi: str = ""
    entries: list[ManifestEntry] = field(default_factory=list)

    def by_kind(self, kind: EntryKind) -> list[ManifestEntry]:
        return [e for e in self.entries if e.kind == kind]

    def search(self, query: str) -> list[ManifestEntry]:
        q = query.lower()
        return [
            e for e in self.entries
            if q in e.title.lower() or q in e.summary.lower()
            or any(q in t.lower() for t in e.tags)
        ]

    def __len__(self) -> int:
        return len(self.entries)


def default_yu_cook_manifest() -> LabManifest:
    """Snapshot of Yu-Cook lab assets as of 2026-05-08.

    Sourced from MEMORY.md / CLAUDE.md / papers_hub_2026_05_04. Update by
    re-running the manifest sync (deferred). Until then, treat this as a
    point-in-time registry — verify against current code before acting.
    """
    return LabManifest(
        pi="Hyeong Won Yu",
        co_pi="Seungho Cook",
        entries=[
            ManifestEntry(
                id="paper1_dm1",
                kind=EntryKind.PAPER,
                title="DM1 molecular dark matter (8-gene + driver-exclusion)",
                summary="Pillar 1 STRONG. TIERA67 ARI=0.90, GSE286332 d=-1.6, "
                        "Hashimoto-like generalization OR 5x.",
                tags=["thyroid", "DM1", "main", "v17"],
                status="active",
                journal_target="Cell Rep Med / JCI Insight",
            ),
            ManifestEntry(
                id="paper2_image_dm1",
                kind=EntryKind.PAPER,
                title="H&E → DM1 / pathology projection",
                summary="N=59 mean AUC 0.83±0.14, ViT-L+CLAM PASS_LAUNCH 2026-05-08",
                tags=["pathology", "image", "DM1"],
                status="active",
            ),
            ManifestEntry(
                id="paper3_ici_track_a",
                kind=EntryKind.PAPER,
                title="ICI vulnerability dark thyroid cancer (Track A FROZEN)",
                summary="Design bundle FROZEN; Track B-lite ATC d=+2.53 myeloid / -2.51 thyroid_diff",
                tags=["ICI", "ATC", "DIAL"],
                status="frozen",
            ),
            ManifestEntry(
                id="paper4_korean_gd_hla",
                kind=EntryKind.PAPER,
                title="Korean GD HLA Pan-Asian (backlog)",
                summary="Renumbered Paper 3 → Paper 4 backlog (4/4 gating)",
                tags=["HLA", "Korean", "Graves"],
                status="backlog",
            ),
            ManifestEntry(
                id="paper11_pancancer",
                kind=EntryKind.PAPER,
                title="Pan-cancer DM1 prognostic / DepMap druggable",
                summary="TCGA pancan 11,069 × 33 lineages; 10 cancers Cox FDR<0.1; "
                        "MYC d=-0.50, NAMPT d=-0.45 in DM1-high",
                tags=["pan-cancer", "DepMap", "prognostic"],
                status="active",
                journal_target="Nat Commun",
            ),
            ManifestEntry(
                id="cohort_korean_k2",
                kind=EntryKind.COHORT,
                title="Korean K2 (PRJEB11591, Yoo 2016 SNU-GMI)",
                summary="n=260 RNA-seq; arcasHLA imputation done; DPB1*05:01 56% replicates Kim 2014",
                tags=["Korean", "RNA-seq", "HLA"],
                status="active",
                n_samples=260,
            ),
            ManifestEntry(
                id="cohort_baba",
                kind=EntryKind.COHORT,
                title="BABA robotic thyroidectomy archive (~5,000 cases)",
                summary="Yu group SNUH archive — used for surgical physical AI / video AI",
                tags=["surgical-video", "BABA", "robotic"],
                status="active",
                n_samples=5000,
            ),
            ManifestEntry(
                id="cohort_ctc_emt",
                kind=EntryKind.COHORT,
                title="PTC perioperative CTC/EMT prospective cohort",
                summary="62 PTC patients, pre/2-week/3-month CTC; 87% pre-op detection",
                tags=["CTC", "EMT", "liquid biopsy"],
                status="active",
                n_samples=62,
            ),
            ManifestEntry(
                id="cohort_maestro",
                kind=EntryKind.COHORT,
                title="MAeSTro / MASTER / KoMPASS active surveillance",
                summary="Multicenter prospective AS vs surgery for low-risk PTMC",
                tags=["active surveillance", "PTMC", "Korean"],
                status="active",
            ),
            ManifestEntry(
                id="dataset_gse286332",
                kind=EntryKind.DATASET,
                title="GSE286332 (Korean PTC vs PTC+HT, n=18)",
                summary="10,380 DEGs; 8-gene d=-1.6, HLA-II d=+3.65, IFN-γ FDR=2e-4",
                tags=["GEO", "Korean", "Hashimoto"],
                status="active", n_samples=18,
            ),
            ManifestEntry(
                id="discovery_dm1_subB_nbnr",
                kind=EntryKind.DISCOVERY,
                title="DM1 sub-B = NBNR cluster (TCGA equivalent of K2)",
                summary="n=56, 96% mutation-negative, 4× Hashimoto-like rate; bridges Paper 1 to K2",
                tags=["DM1", "NBNR", "bridge"],
                status="active",
            ),
            ManifestEntry(
                id="discovery_proteogenomic_v1",
                kind=EntryKind.DISCOVERY,
                title="Proteogenomic v1 corroboration (Mun 2025)",
                summary="thyroid_diff d=-1.91, myeloid d=+1.52; 5/7 8-gene members r<-0.5 along PTC→PDTC→ATC",
                tags=["proteomics", "Mun2025", "DM1"],
                status="active",
            ),
            ManifestEntry(
                id="discovery_rai_dediff",
                kind=EntryKind.DISCOVERY,
                title="RAI-dediff axis (GSE151179 post vs pre-RAI)",
                summary="thyroid_diff d=-1.01 p=0.0001 + HLA-II +0.62 + myeloid +0.55",
                tags=["RAI", "dedifferentiation"],
                status="active",
            ),
            ManifestEntry(
                id="platform_v18",
                kind=EntryKind.PLATFORM,
                title="v18 agentic_research framework",
                summary="5 patterns + 6-component spine; 7 tests pass; 22-slide Korean talk",
                tags=["framework", "patterns"],
                status="active",
                path="project/v18_agentic/",
            ),
            ManifestEntry(
                id="platform_v19_compass",
                kind=EntryKind.PLATFORM,
                title="v19 research_compass (CV → topic recommendation)",
                summary="5-stage pipeline (cv_profile, journal_feed, method_index, data_registry, topic_ranker)",
                tags=["compass", "v19"],
                status="active",
                path="project/v18_agentic/agentic_research/domains/research_compass/",
            ),
            ManifestEntry(
                id="platform_v19_grant_arena",
                kind=EntryKind.PLATFORM,
                title="v19 grant_arena (competitive grant-topic discovery)",
                summary="8 generators × 5 critics × Elo battle × 10-factor scoring; "
                        "frontier_scout for top-researcher paper mining",
                tags=["grant", "agent-arena", "v19"],
                status="active",
                path="project/v18_agentic/agentic_research/domains/research_compass/grant_arena/",
            ),
            ManifestEntry(
                id="platform_neoantigen_hub",
                kind=EntryKind.PLATFORM,
                title="Neoantigen vaccine hub 2026-05-07",
                summary="Post-TESLA platform · 119,237 master / 13 sources / 31 analyses / "
                        "7 live pages at port 80; top KRAS G12D GADGVGKSAL XGBoost proba=0.871",
                tags=["neoantigen", "vaccine", "TESLA"],
                status="active",
            ),
            ManifestEntry(
                id="grant_causal_kthyro",
                kind=EntryKind.GRANT,
                title="Causal K-Thyro Foundation (Samsung 30억 / 3yr)",
                summary="Public-biodata pretrained PTC state model + AS-vs-surgery target trial "
                        "emulation + Korean clinical validation; 8 WPs, 5 weakness audits done",
                tags=["grant", "samsung", "30억", "causal", "thyroid"],
                status="active",
                path="project/grant_samsung_2026/causal_kthyro_weakness_audit.md",
            ),
            ManifestEntry(
                id="method_dial_lite",
                kind=EntryKind.METHOD,
                title="DIAL-lite (4/7 module FLIP audit)",
                summary="Hugo+Riaz audit; n=415 pooled 8 cohorts — DIAL audit necessity proven",
                tags=["DIAL", "ICI", "audit"],
                status="active",
            ),
        ],
    )
