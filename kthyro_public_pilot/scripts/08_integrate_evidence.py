#!/usr/bin/env python3
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from pilot_utils import ensure_standard_dirs, pilot_root, setup_logging


ROWS = [
    "RAI-restorable / redifferentiation axis",
    "HLA/APM immune visibility axis",
    "HLA-low invisible tumor axis",
    "CD8 exclusion / myeloid-CAF barrier axis",
    "Drug-delivery failure proxy axis",
    "Aggressive dedifferentiation / dark-lineage axis",
    "AI companion diagnostic feasibility",
    "Drug/perturbation candidate feasibility",
]


def exists_table(name: str) -> bool:
    return (pilot_root() / "results" / "tables" / name).exists()


def load(name: str) -> pd.DataFrame:
    p = pilot_root() / "results" / "tables" / name
    return pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame()


def main() -> None:
    parser = argparse.ArgumentParser(description="Integrate public-pilot evidence and make GO/NO-GO call.")
    parser.parse_args()
    ensure_standard_dirs()
    logger = setup_logging("08_integrate_evidence")
    tcga = load("tcga_thca_patient_vulnerability_scores.tsv")
    spatial = load("spatial_slide_niche_summary.tsv")
    coherence = load("spatial_coherence_statistics.tsv")
    scrna = load("scrna_celltype_module_scores.tsv")
    drug = load("drug_pilot_candidate_rankings.tsv")
    driver = load("tcga_thca_driver_summary.tsv")

    axis_cols = ["rai_differentiation_score", "immune_visibility_score", "cd8_exclusion_proxy", "drug_delivery_failure_proxy", "aggressive_dedifferentiation_score"]
    supported_axes = 0
    if not tcga.empty:
        for c in axis_cols:
            if c in tcga and tcga[c].notna().sum() >= 100 and tcga[c].std() > 0.1:
                supported_axes += 1
    spatial_coherent = False
    if not coherence.empty:
        niche = coherence[coherence["label_set"].eq("niche_label")]
        spatial_coherent = bool((niche["z"].fillna(0) > 2).sum() >= max(1, len(niche) // 3))
    drug_supported = not drug.empty and "evidence_category" in drug and drug["evidence_category"].astype(str).str.startswith(("A", "B", "C")).any()
    scrna_supported = not scrna.empty and "status" not in scrna.columns

    rows = []
    for row in ROWS:
        if "RAI" in row:
            tcga_msg = f"Observed RAI score variability across {len(tcga)} TCGA tumors." if not tcga.empty else "not available"
            spatial_msg = "RAI-high/low spot maps and niche fractions generated." if not spatial.empty else "not available"
            validation = "Iodide uptake after MEK/BRAF/RET/NTRK perturbation; NIS/TG/TPO/PAX8 mIHC/GeoMx."
            strength = "moderate" if not tcga.empty and not spatial.empty else "weak"
            slide = "RAI state is an axis, not a mutation average."
        elif "HLA/APM immune visibility" in row:
            tcga_msg = "HLA-I/APM, IFN, cytotoxic modules scored separately." if not tcga.empty else "not available"
            spatial_msg = "HLA/APM maps and HLA-visible niche labels generated." if not spatial.empty else "not available"
            validation = "Protein HLA-I/B2M/TAP1 plus CD8/GZMB mIHC; perturbation rescue with IFN/epigenetic candidates."
            strength = "moderate" if not tcga.empty and not spatial.empty else "weak"
            slide = "Immune visibility can be quantified as a separable therapeutic state."
        elif "HLA-low" in row:
            tcga_msg = "HLA-low invisible score combines low APM with epithelial context." if not tcga.empty else "not available"
            spatial_msg = "HLA-low tumor-like niche labels generated." if not spatial.empty else "not available"
            validation = "Confirm HLA-I/B2M loss in PAX8+ tumor regions; no peptide-presentation claims from RNA alone."
            strength = "moderate" if not tcga.empty and not spatial.empty else "weak"
            slide = "HLA-low territories nominate immune-invisibility rescue experiments."
        elif "CD8 exclusion" in row:
            tcga_msg = "CAF/myeloid minus cytotoxic exclusion proxy computed." if not tcga.empty else "not available"
            spatial_msg = "CAF/myeloid/CD8-excluded niche maps and KNN coherence generated." if not spatial.empty else "not available"
            validation = "mIHC for CD8, GZMB, CD68, CD163, ACTA2, FAP, COL1A1 with distance-to-tumor analysis."
            strength = "moderate" if not tcga.empty and spatial_coherent else "weak"
            slide = "Spatial data are needed to distinguish inflamed from excluded immune states."
        elif "Drug-delivery" in row:
            tcga_msg = "ECM/hypoxia minus vascular proxy computed; proxy only." if not tcga.empty else "not available"
            spatial_msg = "Drug-delivery failure proxy maps generated." if not spatial.empty else "not available"
            validation = "Fresh tissue fluorescent drug/liposome/nanoparticle distribution imaging."
            strength = "weak" if not spatial.empty else "not supported"
            slide = "Delivery failure is testable by imaging, but public RNA gives only a proxy."
        elif "Aggressive" in row:
            tcga_msg = "Dedifferentiation/proliferation/MAPK axis computed; DM1 optional from local artifacts." if not tcga.empty else "not available"
            spatial_msg = "Aggressive dedifferentiation score mapped per spot." if not spatial.empty else "not available"
            validation = "Pathology grade, Ki-67/TOP2A, phospho-ERK, CNV/fusion/WES where available."
            strength = "moderate" if not tcga.empty and not spatial.empty else "weak"
            slide = "Aggressive state is measurable as expression program, not only driver mutation."
        elif "AI companion" in row:
            tcga_msg = "Patient-level score table and subtypes generated." if not tcga.empty else "not available"
            spatial_msg = "Spot-level score table and slide-level niche summary generated." if not spatial.empty else "not available"
            validation = "Train pathology-spatial model on paired FFPE, ROI GeoMx, mIHC, perturbation assays."
            strength = "strong" if not tcga.empty and not spatial.empty else "weak"
            slide = "Outputs are already shaped as companion-diagnostic training labels."
        else:
            tcga_msg = "Expression states linked to candidate axes." if not tcga.empty else "not available"
            spatial_msg = "Spatial niches nominate where to test perturbations." if not spatial.empty else "not available"
            validation = "Organoid/slice perturbation panel with iodide uptake, HLA/APM rescue, CD8 co-culture, delivery imaging."
            strength = "moderate" if drug_supported else "weak"
            slide = "Public drug resources nominate candidates; functional restoration must be validated."
        rows.append(
            {
                "axis": row,
                "TCGA expression evidence": tcga_msg,
                "TCGA clinical/driver association": "Driver/clinical summary generated." if not driver.empty and "status" not in driver.columns else "limited or unavailable",
                "spatial evidence": spatial_msg,
                "spatial coherence": "Non-random same-niche structure in multiple slides." if spatial_coherent else "weak/not established or pending",
                "scRNA/cell-type support": "scRNA cell-type module sanity check generated." if scrna_supported else "not available; deconvolution plan required",
                "drug/resource support": "Candidate table generated from DepMap/PRISM-derived public resources." if drug_supported else "limited",
                "experimental validation needed": validation,
                "claim strength": strength,
                "recommended Samsung slide message": slide,
                "risk": "RNA-derived score may conflate cell composition and activation; public data not paired with perturbation.",
                "mitigation": "Validate with FFPE mIHC/GeoMx and fresh tissue functional perturbation/imaging assays.",
            }
        )
    matrix = pd.DataFrame(rows)
    out = pilot_root() / "results" / "tables" / "integrated_public_pilot_evidence_matrix.tsv"
    matrix.to_csv(out, sep="\t", index=False)

    testable_hypotheses = [
        "MEK/BRAF/RET/NTRK-axis inhibition can restore iodide uptake in RAI-low/MAPK-high tumor niches.",
        "HLA/APM-low PAX8+ territories can be rescued by IFN/epigenetic perturbation and validated by HLA-I/B2M/TAP1 protein.",
        "CAF/myeloid-rich regions spatially exclude CD8/GZMB and can be nominated for barrier-modulating combinations.",
        "ECM/hypoxia-high vascular-low territories predict poor fluorescent drug/nanoparticle penetration ex vivo.",
    ]
    if supported_axes >= 3 and spatial_coherent and len(testable_hypotheses) >= 3:
        decision = "GO"
        reason = f"{supported_axes}/5 main axes have TCGA support and spatial same-niche coherence is non-random in cancer slides."
    elif supported_axes >= 3:
        decision = "CONDITIONAL GO"
        reason = f"{supported_axes}/5 main axes have TCGA support, but spatial coherence is weak/missing and requires validation."
    else:
        decision = "NO-GO"
        reason = "Insufficient thyroid-relevant public-data support was loaded."
    report = pilot_root() / "results" / "reports" / "GO_NO_GO_DECISION.md"
    report.write_text(
        f"# GO/NO-GO Decision\n\nDecision: **{decision}**\n\nReason: {reason}\n\n"
        "Testable hypotheses:\n"
        + "\n".join(f"- {h}" for h in testable_hypotheses)
        + "\n\nClaim boundary: this is hypothesis-generating public-data evidence, not clinical deployment evidence.\n"
    )
    logger.info("Decision: %s. Evidence matrix written to %s", decision, out)


if __name__ == "__main__":
    main()
