"""BiomniAdapter — interface stub for Stanford SNAP's Biomni agent.

Biomni (https://github.com/snap-stanford/biomni) is an LLM-backed
biomedical research agent with code execution + retrieval-augmented
planning + 11GB datalake. Public API:

    from biomni.agent import A1
    agent = A1(path='./data', llm='claude-sonnet-4-5')
    agent.go("plan a CRISPR screen for X")

This adapter does NOT install or invoke Biomni. It defines the capability
surface so our runner can route queries to Biomni *as if* it were
available, and stub responses for offline tests/demos. v0.3 wires the
real `biomni.agent.A1.go(...)` and Biomni-R0 reasoning model.

Marathon-mode rationale: Biomni runs full-system code execution under an
LLM. Wiring it live before paper ship is a security + cost + distraction
hazard. Interface now, install later.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence


@dataclass
class BiomniCapability:
    name: str                        # e.g. "crispr_screen_plan"
    description: str
    inputs: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)


def default_biomni_capabilities() -> list[BiomniCapability]:
    """Capability list distilled from Biomni's README (2026-05 snapshot)."""
    return [
        BiomniCapability(
            name="crispr_screen_plan",
            description="Plan and execute CRISPR screens to identify genes regulating a process",
            inputs=["phenotype", "cell_system", "screen_format"],
            outputs=["gene_panel", "guide_RNA_design", "readout_protocol"],
        ),
        BiomniCapability(
            name="scrna_annotation",
            description="Single-cell RNA-seq cell-type annotation + hypothesis generation",
            inputs=["scRNA_matrix", "tissue_context"],
            outputs=["cell_type_labels", "marker_set", "hypothesis_list"],
        ),
        BiomniCapability(
            name="admet_prediction",
            description="Predict ADMET properties for chemical compounds",
            inputs=["smiles_or_inchi"],
            outputs=["admet_table", "liability_flags"],
        ),
        BiomniCapability(
            name="gwas_causal_gene",
            description="GWAS causal gene identification + variant prioritization",
            inputs=["gwas_summary_stats", "phenotype"],
            outputs=["ranked_gene_list", "variant_effect_table"],
        ),
        BiomniCapability(
            name="rare_disease_diagnosis",
            description="Diagnose rare diseases from patient phenotype + variant data",
            inputs=["hpo_terms", "vcf_or_variant_list"],
            outputs=["differential_diagnosis", "evidence_per_candidate"],
        ),
        BiomniCapability(
            name="lab_qa",
            description="Lab-bench question answering over PubMed + protocol library",
            inputs=["natural_language_question"],
            outputs=["answer", "citations"],
        ),
    ]


@dataclass
class BiomniAdapter:
    """Thin facade over a (future) Biomni A1 agent.

    Use ``available=True`` once Biomni is installed and the LLM key is
    configured. While ``available=False`` (default), every call returns a
    deterministic stub describing what *would* happen.
    """

    available: bool = False
    capabilities: list[BiomniCapability] = field(default_factory=default_biomni_capabilities)
    llm: str = "claude-sonnet-4-5"
    datalake_path: str = "./data"

    def list_capabilities(self) -> list[BiomniCapability]:
        return list(self.capabilities)

    def supports(self, capability_name: str) -> bool:
        return any(c.name == capability_name for c in self.capabilities)

    def go(self, task: str, *, capability: str | None = None) -> dict:
        """Dispatch a free-text task; in stub mode, return a structured no-op.

        Live mode (v0.3) calls ``A1(path=..., llm=...).go(task)`` and parses
        the trajectory. We intentionally do not import biomni here so the
        repo runs without the 11GB datalake or biomni install.
        """
        if not self.available:
            matched = capability or self._best_capability(task)
            return {
                "status": "stubbed",
                "matched_capability": matched,
                "task": task,
                "note": ("Biomni not installed/enabled. Set adapter.available=True "
                         "and `pip install biomni` plus configure an LLM API key to run live."),
            }
        raise NotImplementedError("Live Biomni wiring deferred to v0.3")

    def _best_capability(self, task: str) -> str | None:
        t = task.lower()
        keyword_to_cap = {
            "crispr": "crispr_screen_plan",
            "scrna": "scrna_annotation",
            "single cell": "scrna_annotation",
            "admet": "admet_prediction",
            "compound": "admet_prediction",
            "gwas": "gwas_causal_gene",
            "rare disease": "rare_disease_diagnosis",
            "phenotype": "rare_disease_diagnosis",
        }
        for k, cap in keyword_to_cap.items():
            if k in t:
                return cap
        return "lab_qa"
