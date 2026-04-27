#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from v17p3_common import RPT, TAB


def j(name: str) -> dict:
    p = TAB / name
    return json.loads(p.read_text()) if p.exists() else {}


def t(name: str) -> pd.DataFrame:
    p = TAB / name
    return pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame()


def write(path: Path, text: str) -> None:
    path.write_text(text.strip() + "\n", encoding="utf-8")


def main() -> None:
    f1, f2, f3 = j("F1_summary.json"), j("F2_summary.json"), j("F3_summary.json")
    a1, a2, a3 = j("A1_summary.json"), j("A2_summary.json"), j("A3_summary.json")
    a4, a5, a6 = j("A4_summary.json"), j("A5_summary.json"), j("A6_summary.json")
    ext = t("F1_external_5cohort_recovery.tsv")
    gsea = t("F2_gsea_hallmark_proper.tsv")
    drugs = t("A3_top_drugs_dm1_selective.tsv")
    outline = f"""
# v17p3 Paper Outline FINAL

## Title 후보
1. A driver-orthogonal transcriptional axis predicts differentiation and therapeutic vulnerability in thyroid cancer
2. DM1/DM2 define an actionable transcriptional taxonomy beyond BRAF/RAS in papillary thyroid carcinoma
3. From dark matter to clinical actionability: multi-cohort validation of thyroid cancer transcriptional states
4. Recovered external validation and proper GSEA define DM1/DM2 as an actionable thyroid cancer axis
5. RAI-relevant thyroid TF circuitry stratifies DM1/DM2 beyond canonical drivers

## Abstract angle
- Biology: DM1/DM2는 driver-negative residual bucket이 아니라 thyroid differentiation/MAPK axis를 반영하는 subtype이다.
- Method: external gene-ID recovery와 proper GSEA로 Phase 2의 약점을 보정했다.
- Clinical: BRAF/RAS와 직교하는 DM score가 RAI-related score와 drug sensitivity를 stratify한다.

## Main figures
1. DM1/DM2 overview and Phase 1-3 timeline
2. Robustness + external recovery
3. Proper GSEA and TF/thyroid differentiation axis
4. Clinical reframe and age-independent signature
5. BRAF/RAS-positive internal stratification
6. Drug actionability
7. RAI prediction / PDTC validation
"""
    write(RPT / "v17p3_paper_outline_FINAL.md", outline)
    resp = f"""
# v17p3 Response To Reviewers V2

1. TCGA-only? → `F1_5cohort_forest.html`, recovered robust cohorts={f1.get('robust_cohort_count_dia_auc_gt_0_85')}
2. Biology proxy only? → proper GSEA hallmark FDR<0.01 count={f2.get('hallmark_fdr_lt_0_01')}
3. Clinical relevance weak? → alternative endpoints p<0.05 count={f3.get('n_endpoints_p_lt_0_05')}
4. Driver groups은? → `A2_dm_score_by_driver.html`
5. Actionability? → `A3_drug_volcano_dm1_vs_dm2.html`, `A6_rai_score_distribution.html`
6. Pan-cancer implication? → `A4_pancancer_heatmap.html`
7. Genomic depth? → `A5_immune_landscape.html` with honest CNV missing caveat
"""
    write(RPT / "v17p3_response_to_reviewers_v2.md", resp)
    dump = f"""
# v17p3 External Review Dump V3

## 1. F1 external recovery
- recovered robust cohorts: {f1.get('robust_cohort_count_dia_auc_gt_0_85')}
- GSE213647 TierA overlap: {f1.get('gse213647_gene_overlap_tiera')}

## 2. F2 proper GSEA
- Hallmark FDR<0.01: {f2.get('hallmark_fdr_lt_0_01')}
- Reactome FDR<0.05: {f2.get('reactome_fdr_lt_0_05')}
- proxy/proper sign agreement: {f2.get('proxy_proper_sign_agreement')}

## 3. F3 alternative endpoints
- significant endpoint count: {f3.get('n_endpoints_p_lt_0_05')}
- best endpoint: {f3.get('best_endpoint')}

## 4. A1 scRNA
- patients: {a1.get('patients')}
- mixed patients: {a1.get('mixed_patients')}
- dominance range: {a1.get('dominance_ratio_range')}

## 5. A2 BRAF/RAS internal stratification
- corr(BRAF, probDM2): {a2.get('driver_dm_correlation_braf_vs_probdm2')}
- BRAF/DM2-like n: {a2.get('braf_dm2_like_n')}
- RAS/DM1-like n: {a2.get('ras_dm1_like_n')}

## 6. A3 drug response
- FDR<0.1 selective drugs: {a3.get('n_drugs_fdr_lt_0_1')}
- best DM1-selective drug: {a3.get('best_dm1_drug')}

## 7. A4 pan-cancer
- applicable cancers: {a4.get('applicable_cancers')}

## 8. A5 genomic/immune
- TMB DM2-DM1: {a5.get('tmb_dm2_minus_dm1')}

## 9. A6 RAI prediction
- PDTC low-RAI AUC: {a6.get('pdtc_auc_using_low_rai')}

## 10. Venue update
- Bioinformatics: 80%
- npj Precision Oncology: 65%
- Genome Medicine: 35%
- Nature Communications: 10%
"""
    write(RPT / "v17p3_external_review_dump_v3.md", dump)


if __name__ == "__main__":
    main()

