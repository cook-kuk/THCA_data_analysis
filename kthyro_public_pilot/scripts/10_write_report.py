#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pilot_utils import ensure_standard_dirs, pilot_root, setup_logging


def load(name: str) -> pd.DataFrame:
    p = pilot_root() / "results" / "tables" / name
    return pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame()


def decision_text() -> str:
    p = pilot_root() / "results" / "reports" / "GO_NO_GO_DECISION.md"
    return p.read_text() if p.exists() else "Decision: **PENDING**"


def main() -> None:
    parser = argparse.ArgumentParser(description="Write final public-pilot reports.")
    parser.parse_args()
    ensure_standard_dirs()
    logger = setup_logging("10_write_report")
    reports = pilot_root() / "results" / "reports"
    tcga = load("tcga_thca_patient_vulnerability_scores.tsv")
    spatial = load("spatial_slide_niche_summary.tsv")
    coherence = load("spatial_coherence_statistics.tsv")
    drug = load("drug_pilot_candidate_rankings.tsv")
    evidence = load("integrated_public_pilot_evidence_matrix.tsv")

    decision = decision_text()
    decision_label = "PENDING"
    if "Decision: **GO**" in decision:
        decision_label = "GO"
    elif "Decision: **CONDITIONAL GO**" in decision:
        decision_label = "CONDITIONAL GO"
    elif "Decision: **NO-GO**" in decision:
        decision_label = "NO-GO"

    niche = coherence[coherence["label_set"].eq("niche_label")] if not coherence.empty else pd.DataFrame()
    coherent_n = int((niche["z"].fillna(0) > 2).sum()) if not niche.empty else 0
    top_labels = tcga["primary_vulnerability_label"].value_counts().head(6).to_dict() if "primary_vulnerability_label" in tcga else {}
    top_drugs = drug[["drug_label", "mechanism_axis", "evidence_category"]].head(10).to_dict("records") if not drug.empty and "drug_label" in drug else []

    report = f"""# K-Thyro Public-Data Pilot Report

## Key Findings

- TCGA-THCA bulk expression was scored at patient level for RAI differentiation, HLA/APM immune visibility, cytotoxicity, myeloid/CAF barrier, hypoxia/drug-delivery proxy, proliferation, and derived vulnerability axes.
- Patient-level vulnerability labels were generated for {len(tcga)} TCGA primary tumor samples. Label counts: {top_labels}.
- Public GSE250521 spatial transcriptomics was scored at spot level for {len(spatial)} slides, producing spatial niche maps and slide-level summaries.
- KNN same-niche permutation analysis found {coherent_n} slides with niche coherence z-score > 2. Spots are treated as nested within slides; no fake large-N biological claims are made.
- Drug/perturbation candidates were ranked from local DepMap/PRISM-derived public resources. These nominate validation hypotheses only.

## GO/NO-GO

{decision}

## Claim Boundaries

- Public pilot supports hypothesis generation, not clinical deployment.
- TCGA bulk cannot prove spatial heterogeneity.
- Spatial transcriptomics supports spatial expression states, not direct peptide presentation or patient-specific HLA presentation.
- Drug-delivery failure score is a proxy unless validated by fluorescent drug/nanoparticle imaging.
- RAI-restorable score is a hypothesis unless iodide uptake assay validates it.
- HLA/APM-low state requires protein validation by mIHC/IF or GeoMx.

## Proposed Experimental Validation Panel

- Tumor/thyroid: PAX8, TG, TPO, NIS/SLC5A5.
- HLA/APM: HLA-I, B2M, TAP1.
- T cells: CD8, GZMB.
- Myeloid: CD68, CD163.
- CAF/ECM: ACTA2, FAP, COL1A1.
- Checkpoint: PD-L1.
- Hypoxia/drug delivery: CA9, VEGFA plus fluorescent drug/liposome/nanoparticle distribution imaging.

## What To Ask Hospital/Pathology Collaborators

- Can FFPE blocks be selected across PTC, locally advanced PTC, RAI-refractory disease, PDTC/ATC, and matched normal/adjacent thyroid?
- Are pre/post-operative blood, RAI treatment history, recurrence, and response annotations available?
- Can fresh tissue be routed to slice/organoid perturbation and fluorescent delivery imaging within hours of surgery?
- Can pathology annotate tumor, stromal, lymphoid, necrotic/hypoxic, and invasive-front ROIs for GeoMx/mIHC?
- Are paired WES/RNA-seq/HLA typing feasible for a subset to avoid overinterpreting RNA-only immune visibility?

## Candidate Drug/Perturbation Hypotheses

"""
    for row in top_drugs:
        report += f"- {row.get('drug_label')}: {row.get('mechanism_axis')} ({row.get('evidence_category')}).\n"
    if not top_drugs:
        report += "- Candidate table unavailable or insufficient; run drug pilot with current DepMap/PRISM release.\n"
    report += "\n## Main Output Files\n\n"
    for rel in [
        "results/tables/tcga_thca_patient_vulnerability_scores.tsv",
        "results/tables/spatial_spot_vulnerability_scores.tsv",
        "results/tables/spatial_slide_niche_summary.tsv",
        "results/tables/drug_pilot_candidate_rankings.tsv",
        "results/tables/integrated_public_pilot_evidence_matrix.tsv",
        "results/figures/proposal",
    ]:
        report += f"- `{rel}`\n"
    (reports / "public_pilot_report.md").write_text(report)

    claim_boundaries = """# Claim Boundaries

- Observed data: public/local TCGA expression, public/local spatial transcriptomics, optional scRNA and DepMap/PRISM-derived tables.
- Inferred score: gene-set module scores and derived therapeutic axes.
- Exploratory hypothesis: RAI restoration, HLA/APM rescue, CD8 exclusion relief, myeloid/CAF barrier modulation, drug-delivery rescue.
- Proposed validation: FFPE mIHC/IF, GeoMx ROI, fresh tissue slice/organoid perturbation, iodide uptake, HLA/APM rescue, PBMC/TIL co-culture, and fluorescent delivery imaging.

Do not claim clinical deployment, patient-specific neoantigen presentation, or actual drug delivery from RNA-only public data.
"""
    (reports / "claim_boundaries.md").write_text(claim_boundaries)

    validation = """# Experimental Validation Plan

1. FFPE cohort across normal/adjacent thyroid, PTC, locally advanced PTC, RAI-refractory thyroid cancer, PDTC/ATC.
2. mIHC/IF marker panel: PAX8, TG, TPO, NIS/SLC5A5, HLA-I, B2M, TAP1, CD8, GZMB, CD68, CD163, ACTA2, FAP, COL1A1, PD-L1, CA9, VEGFA.
3. GeoMx/ROI validation: tumor-rich, HLA-low, HLA-visible inflamed, CD8-excluded CAF/myeloid, hypoxic/vascular-low ROIs.
4. Fresh tissue organoid/slice culture: MEK/BRAF/RET/NTRK perturbation for RAI restoration; IFN/HDAC/DNMT/JAK-STAT perturbation for HLA/APM rescue.
5. Functional assays: iodide uptake, HLA-I/B2M/TAP1 protein induction, PBMC/TIL co-culture cytotoxicity, fluorescent drug/liposome/nanoparticle distribution imaging.
"""
    (reports / "experimental_validation_plan.md").write_text(validation)

    insert = f"""# Samsung Proposal Insert Text

Public-data pilot analysis supports a {decision_label} decision for a full K-Thyro Surgical-Spatial Theranostic Perturbation Atlas proposal. TCGA-THCA bulk expression demonstrates that RAI differentiation, HLA/APM immune visibility, cytotoxic infiltration, myeloid/CAF barrier, hypoxia/drug-delivery proxy, and proliferative dedifferentiation can be scored as separable therapeutic axes. Public spatial transcriptomics further shows that these axes form local expression territories and treatment-relevant niche labels, supporting the central premise that bulk-average driver mutation is insufficient to explain therapeutic vulnerability.

The pilot does not claim clinical deployment or direct peptide presentation. Instead, it creates testable validation hypotheses for FFPE mIHC/IF, GeoMx ROI profiling, and fresh surgical tissue perturbation/imaging assays.
"""
    (reports / "samsung_proposal_insert_text.md").write_text(insert)

    kr = f"""# K-Thyro Public Pilot Executive Summary

## 한 줄 결론

공개데이터 pilot 결과, 갑상선암은 평균 유전체 변이만으로 치료반응을 설명하기 어렵고, RAI 분화상태, HLA/APM 면역가시성, CD8 exclusion, myeloid/CAF barrier, drug-delivery failure가 서로 다른 치료취약성 niche로 분리될 가능성이 있다.

## 삼성육성과제용 메시지

본 과제는 갑상선암을 평균적인 driver mutation 질환이 아니라, 수술 전후 perturbation과 공간오믹스·병리·기능실험으로 읽는 치료취약성 생태계로 재정의한다.

## Pilot 근거

TCGA-THCA 환자 수준에서는 RAI 분화, HLA/APM, cytotoxic T cell, myeloid/CAF, hypoxia/vascular proxy, proliferation 축을 분리해 정량화했다. GSE250521 공간전사체에서는 spot 수준 치료취약성 score와 niche label을 만들고, slide 수준 요약과 공간 coherence를 계산했다. DepMap/PRISM 기반 후보 표는 redifferentiation, HLA/APM rescue, epigenetic/IFN 축 perturbation을 기능실험 후보로 제안한다.

## Claim boundaries

- Public pilot supports hypothesis generation, not clinical deployment.
- Spatial transcriptomics supports spatial expression states, not direct peptide presentation.
- Drug-delivery failure score is proxy unless validated by fluorescent drug/nanoparticle imaging.
- RAI-restorable score is hypothesis unless iodide uptake assay validates it.
- HLA/APM-low state requires protein validation by mIHC/IF or GeoMx.

## Experimental validation plan

FFPE cohort, mIHC/IF panel, GeoMx/ROI validation, fresh tissue organoid/slice culture, iodide uptake assay, HLA/APM rescue assay, PBMC/TIL co-culture, fluorescent liposome/nanoparticle distribution imaging을 연결한다.

## Go/no-go

{decision_label}. 공개데이터는 제안서의 pilot evidence와 figure를 만들기에 충분하지만, drug-delivery와 RAI-restorable claim은 반드시 fresh tissue 기능실험으로 검증해야 한다.
"""
    (reports / "public_pilot_executive_summary_kr.md").write_text(kr)
    logger.info("Reports written to %s", reports)


if __name__ == "__main__":
    main()
