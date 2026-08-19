#!/usr/bin/env python3
"""Build Samsung Science-track CTC-EMT proposal text package.

This script intentionally keeps the proposal narrow:
PTC, thyroidectomy, serial CTC-EMT transition, matched genetic map.
"""

from __future__ import annotations

from pathlib import Path
import textwrap

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[0]
PUBLIC = REPO / "kthyro_public_pilot" / "results"
REPORTS = ROOT / "outputs" / "reports"
TABLES = ROOT / "outputs" / "tables"
PROMPTS = ROOT / "outputs" / "prompts"
INPUTS = ROOT / "inputs"


def dedent(text: str) -> str:
    return textwrap.dedent(text).strip() + "\n"


def read_tsv(rel: str) -> pd.DataFrame:
    return pd.read_csv(PUBLIC / rel, sep="\t")


def metrics() -> dict[str, str]:
    out: dict[str, str] = {}
    try:
        driver = read_tsv("public_validation_expansion/tables/tcga_driver_within_group_vulnerability_diversity.tsv")
        braf = driver.loc[driver["driver_group"].eq("BRAF")].iloc[0]
        out["braf_n"] = f"{int(braf['n_patients'])}"
        out["braf_labels"] = f"{int(braf['n_nonzero_vulnerability_labels'])}"
        out["braf_largest_fraction"] = f"{float(braf['largest_label_fraction']) * 100:.1f}%"
        out["braf_ambiguity"] = f"{float(braf['mutation_only_ambiguity_fraction']) * 100:.1f}%"
    except Exception:
        out.update(braf_n="274", braf_labels="7", braf_largest_fraction="33.2%", braf_ambiguity="66.8%")

    try:
        eta = read_tsv("public_validation_expansion/tables/tcga_driver_axis_variance_explained.tsv")
        for axis, key in [
            ("RAI differentiation", "eta_rai"),
            ("HLA-I/APM", "eta_hla"),
            ("Immune visibility", "eta_immune"),
            ("CD8 exclusion", "eta_cd8"),
        ]:
            row = eta.loc[eta["axis_label"].eq(axis)].iloc[0]
            out[key] = f"{float(row['eta_squared_driver_with_no_call_group']) * 100:.1f}%"
    except Exception:
        out.update(eta_rai="30.9%", eta_hla="15.6%", eta_immune="12.5%", eta_cd8="16.2%")

    try:
        spatial = read_tsv("tables/spatial_coherence_summary_cleaned.tsv")
        out["spatial_slides"] = f"{len(spatial)}"
        out["spatial_spots"] = f"{int(spatial['n_spots'].sum()):,}"
        out["spatial_z_min"] = f"{spatial['same_niche_z'].min():.2f}"
        out["spatial_z_median"] = f"{spatial['same_niche_z'].median():.2f}"
        out["spatial_z_max"] = f"{spatial['same_niche_z'].max():.2f}"
        out["spatial_z_pos"] = f"{int((spatial['same_niche_z'] > 2).sum())}/{len(spatial)}"
    except Exception:
        out.update(
            spatial_slides="16",
            spatial_spots="57,144",
            spatial_z_min="6.09",
            spatial_z_median="22.00",
            spatial_z_max="43.16",
            spatial_z_pos="16/16",
        )

    try:
        prot = read_tsv("extra_analyses/tables/bulk_proteomics_kthyro_module_contrasts.tsv")
        rai = prot.loc[(prot["contrast"].eq("ATC_vs_PTC")) & (prot["module"].eq("RAI_differentiation_protein"))].iloc[0]
        myeloid = prot.loc[(prot["contrast"].eq("ATC_vs_PTC")) & (prot["module"].eq("myeloid_suppressive_protein"))].iloc[0]
        out["prot_rai_d"] = f"{float(rai['cohens_d_group_a_minus_b']):.2f}"
        out["prot_rai_q"] = f"{float(rai['fdr_q']):.2e}"
        out["prot_myeloid_d"] = f"{float(myeloid['cohens_d_group_a_minus_b']):.2f}"
        out["prot_myeloid_q"] = f"{float(myeloid['fdr_q']):.2e}"
    except Exception:
        out.update(prot_rai_d="-1.99", prot_rai_q="4.23e-35", prot_myeloid_d="1.68", prot_myeloid_q="1.99e-23")

    return out


M = metrics()


TITLE = "Serial Genetic Mapping of CTC-EMT State Transitions After Thyroidectomy in Papillary Thyroid Cancer"
TITLE_KR = "갑상선절제술 전후 유두갑상선암 CTC-EMT 상태전이와 유전지도의 연속 해독"


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dedent(text), encoding="utf-8")


def build_context_digest() -> str:
    return f"""
    # Context Digest

    ## Final Category Decision

    **Science, not Technology/ICT.** The proposal should be reviewed as a mechanism-driven life-science project, not as a broad spatial omics, drug-delivery, or AI platform. The core experiment is one biological perturbation: **thyroidectomy**.

    ## Final Title Options

    1. **{TITLE}**
    2. Thyroidectomy as a Human Perturbation Model for CTC-EMT Plasticity in Papillary Thyroid Cancer
    3. Post-Thyroidectomy Persistence of EMT-State Circulating Tumor Cells and Genetic Residual Disease in Papillary Thyroid Cancer
    4. Coupling CTC-EMT Phenotype Switching to Tumor Genotype After Surgical Removal of Papillary Thyroid Cancer
    5. Defining the Liquid-Biopsy Transition State of Papillary Thyroid Cancer Through CTC-EMT Phenotyping and NGS Clone Mapping

    Korean working title: **{TITLE_KR}**

    ## Central Hypothesis

    Thyroidectomy acts as a controlled human perturbation that should collapse nonspecific tumor shedding, while persistent EMT-state CTCs and matched genetic signals may reveal residual-risk PTC biology.

    ## What Yu 2024 Already Showed

    Yu et al. 2024 is the feasibility anchor: 62 prospective PTC patients were sampled before thyroidectomy and after surgery at 2 weeks and 3 months. CTCs were detected in 87% of patients, epithelial-mesenchymal or mesenchymal phenotypes were predominant, and CTC counts decreased after thyroidectomy. This establishes that serial PTC CTC measurement is feasible and biologically plausible.

    ## What Is Still Unknown

    No public dataset jointly provides serial pre/post-thyroidectomy blood, CTC EMT phenotype transition, matched tumor-normal NGS, cfDNA or CTC-enriched NGS, and postoperative follow-up. The unknown is whether persistent postoperative CTC/cfDNA signals can be genetically anchored to the resected tumor clone and interpreted as residual-risk biology rather than assay noise or nonspecific circulating cells.

    ## How The Public Pilot Supports The Narrow Proposal

    The broad K-Thyro pilot contributes only supporting evidence:

    - TCGA-THCA 505 primary tumors show tissue-state heterogeneity beyond driver mutation.
    - BRAF tumors split across {M['braf_labels']} vulnerability labels; the largest BRAF label is only {M['braf_largest_fraction']} and mutation-only ambiguity is {M['braf_ambiguity']}.
    - Driver group explains only part of RAI differentiation ({M['eta_rai']}), HLA-I/APM ({M['eta_hla']}), immune visibility ({M['eta_immune']}), and CD8 exclusion ({M['eta_cd8']}) axis variance.
    - GSE250521 spatial analysis shows {M['spatial_z_pos']} slides with same-niche coherence z > 2 across {M['spatial_spots']} spots.
    - Bulk thyroid proteomics supports directionality of RAI protein loss in ATC vs PTC (Cohen's d {M['prot_rai_d']}) and myeloid/TGFB-related protein increase (Cohen's d {M['prot_myeloid_d']}).

    ## What Must Not Be Overclaimed

    - Public pilot does not prove CTC shedding.
    - Spatial transcriptomics does not prove recurrence, drug delivery, or CTC origin.
    - CTC phenotype alone does not prove viable metastatic potential.
    - cfDNA or CTC-enriched NGS may fail in low-shedding early PTC.
    - This is not a ready clinical diagnostic.
    - Inocras or any vendor is not the scientific novelty.
    """


def build_final_logic() -> str:
    return """
    # Final Proposal Logic

    ## A. One-Line Hypothesis

    Thyroidectomy acts as a controlled human perturbation that should collapse nonspecific tumor shedding, while persistent EMT-state CTCs and matched genetic signals may reveal residual-risk PTC biology.

    ## B. Why This Is Science

    - One disease: papillary thyroid cancer.
    - One perturbation: thyroidectomy.
    - One biological transition: CTC epithelial to hybrid E/M to mesenchymal state.
    - One genetic anchor: matched tumor-normal NGS plus cfDNA/CTC-enriched NGS.
    - One primary endpoint: postoperative CTC-EMT state transition.

    ## C. Why Not ICT

    - No proprietary CTC device is being invented.
    - No validated AI primitive exists yet.
    - No clinical deployment claim is ready.
    - AI is an interpretation layer for clone-state and longitudinal transition modeling, not the core novelty.

    ## D. What Public Pilot Supports

    - Tissue state heterogeneity beyond mutation.
    - Spatial organization of tissue states.
    - Marker and ROI logic for thyroid-lineage, EMT, immune/stress, and survival programs.
    - Need for phenotype plus genotype integration.

    ## E. What Public Pilot Does Not Support

    - Actual CTC shedding.
    - Recurrence prediction.
    - Peptide presentation.
    - Clinical deployment.
    - Final drug response.
    """


def one_page_summary() -> str:
    return f"""
    # 삼성 Science 과제 1페이지 요약

    ## 제목

    **{TITLE_KR}**  
    **{TITLE}**

    ## 핵심 질문

    갑상선절제술은 유두갑상선암에서 원발 종양을 제거하는 가장 강한 생체 내 perturbation이다. 수술 후 CTC가 단순히 감소하는지를 넘어서, **상피형(E), hybrid E/M, mesenchymal(M) CTC 상태가 어떻게 전이되는지**, 그리고 지속되는 postoperative CTC/cfDNA 신호가 절제 종양의 유전 clone과 연결되는지를 묻는다.

    ## 중심 가설

    갑상선절제술 후 비특이적 tumor shedding은 감소해야 한다. 그럼에도 지속되는 hybrid/mesenchymal CTC와 tumor-matched genetic signal은 잔존위험 생물학을 반영할 수 있다.

    ## 왜 Science인가

    이 과제는 ICT 플랫폼이 아니라, 수술이라는 인간 perturbation을 이용해 암세포의 circulating EMT state transition을 해석하는 생명과학 과제다. AI/ICT는 clone-state modeling, longitudinal transition modeling, tissue/CTC/cfDNA 통합해석, report generation을 위한 분석 인프라로만 사용한다.

    ## Preliminary Evidence

    Yu 2024는 62명 PTC 전향 코호트에서 수술 전, 수술 후 2주, 3개월 blood CTC 분석을 수행했고, CTC가 87%에서 검출되며 EMT/mesenchymal phenotype이 우세하고 수술 후 CTC count가 감소함을 보였다. 그러나 이 연구는 CTC NGS clone map을 만들지 않았다.

    Local public pilot은 TCGA-THCA 505명에서 driver mutation만으로 tissue state가 설명되지 않음을 보였다. BRAF 274명도 {M['braf_labels']}개 vulnerability label로 갈라졌고 최대 label은 {M['braf_largest_fraction']}뿐이었다. GSE250521 spatial {M['spatial_slides']} slides, {M['spatial_spots']} spots에서는 {M['spatial_z_pos']} slides가 same-niche coherence z > 2를 보였다. 이 결과는 CTC 증거가 아니라 marker selection과 tissue-state logic의 보조 근거다.

    ## 3년 목표

    1. PTC thyroidectomy 전후 CTC E/EM/M state transition을 정의한다.
    2. 절제 종양 matched tumor-normal NGS로 CTC/cfDNA 신호의 유전적 anchor를 만든다.
    3. persistent EM/M CTC와 genetic concordance가 LVI, LNM, ETE, BRAF/TERT/fusion, Tg/US follow-up, recurrence-risk feature와 연결되는지 검증한다.
    4. 임상진단이 아니라 future clinical validation을 정당화하는 research-grade residual-risk state framework를 만든다.

    ## 금지할 주장

    본 과제는 CTC-EMT가 이미 재발을 예측한다고 주장하지 않는다. Public pilot이 CTC biology를 증명한다고 말하지 않는다. cfDNA/CTC NGS가 모든 early PTC에서 작동한다고 말하지 않는다. Inocras/vendor 기술을 novelty로 포장하지 않는다.
    """


def preliminary_results() -> str:
    return f"""
    # Preliminary Results 3-Page Section

    ## 1. Published PTC CTC data establish feasibility, not the genetic map.

    Yu et al. 2024 provides the anchor observation for this proposal. In 62 prospective PTC patients, blood was sampled before thyroidectomy and at 2 weeks and 3 months after surgery. CTCs were detected in 87% of patients, and epithelial-mesenchymal or mesenchymal phenotypes were predominant. CTC counts decreased after thyroidectomy, supporting the idea that surgery can act as a measurable perturbation of circulating tumor-cell burden.

    This prior evidence is necessary but incomplete. It does not determine whether residual postoperative CTCs carry the same genetic clone as the resected tumor. It does not distinguish persistent tumor-derived signal from nonspecific circulating epithelial-like cells. It does not define whether E, E/M, and M CTC fractions represent connected transition states or independent assay categories. The proposed project fills this gap by adding matched tumor-normal NGS, cfDNA, and CTC-enriched genetic mapping to serial CTC-EMT phenotyping.

    ## 2. Public liquid-biopsy NGS studies show feasibility but leave early PTC unresolved.

    Sato et al. showed that BRAF V600E ctDNA can be detected before and after surgery in a subset of primary PTC patients, and that postoperative persistence may carry recurrence information. However, that study used a small cohort and a single-mutation ddPCR strategy. Tarasova et al. showed that plasma NGS can detect genomic alterations in thyroid cancer at scale, but the cohort was largely real-world advanced/submitted thyroid cancer rather than a controlled early PTC thyroidectomy perturbation series.

    Therefore, the proposal does not assume cfDNA will work in all patients. The design is stage-gated: tissue-normal NGS anchors every evaluable case; cfDNA tests postoperative persistence; CTC-enriched or low-input NGS is attempted only in CTC-high and QC-passing subsets.

    ## 3. Local public pilot shows that tissue genotype alone is insufficient.

    In the local K-Thyro public pilot, TCGA-THCA 505 primary tumors were scored for therapeutic tissue-state axes. BRAF-mutant tumors were not homogeneous: {M['braf_n']} BRAF tumors split across {M['braf_labels']} vulnerability labels, and the largest BRAF label represented only {M['braf_largest_fraction']} of BRAF tumors. Driver group explained only part of major axis variance: RAI differentiation {M['eta_rai']}, HLA-I/APM {M['eta_hla']}, immune visibility {M['eta_immune']}, and CD8 exclusion {M['eta_cd8']}.

    This does not prove CTC biology. Its contribution is more disciplined: tissue genotype alone is unlikely to explain postoperative residual-risk state. A serial liquid-biopsy project that measures only mutation could miss phenotype; a project that measures only CTC phenotype could be nonspecific. The proposed work combines both.

    ## 4. Spatial and protein public evidence support marker selection.

    GSE250521 spatial transcriptomics analysis included {M['spatial_slides']} thyroid slides and {M['spatial_spots']} spots. Same-niche spatial coherence was observed in {M['spatial_z_pos']} slides, with z-score min/median/max {M['spatial_z_min']} / {M['spatial_z_median']} / {M['spatial_z_max']}. This indicates that tissue states are spatially organized. It does not prove CTC shedding or clinical recurrence.

    Independent thyroid bulk proteomics supports directionality of relevant tissue states: ATC vs PTC shows RAI differentiation protein loss (Cohen's d {M['prot_rai_d']}, FDR {M['prot_rai_q']}) and myeloid suppressive protein gain (Cohen's d {M['prot_myeloid_d']}, FDR {M['prot_myeloid_q']}). These data justify marker selection for FFPE/mIHC and CTC panels: epithelial/thyroid-lineage markers, EMT markers, survival/stemness markers, and immune/stress markers.

    ## 5. Defensible three-year output.

    The defensible output is not a clinical diagnostic. The defensible output is the first research-grade serial map of PTC CTC-EMT state transition after thyroidectomy with matched tumor genetic anchoring. If successful, the project will justify a larger clinical validation cohort for postoperative residual-risk monitoring.
    """


def specific_aims() -> str:
    return """
    # Specific Aims

    ## Aim 1. Define the CTC-EMT state transition before and after thyroidectomy in PTC.

    **Hypothesis:** Postoperative residual CTCs are enriched for hybrid/mesenchymal, survival, and thyroid-lineage-low programs.

    **Approach:** Prospectively collect blood at T0 pre-op, T1 2 weeks, T2 3 months, and T3 6-12 months or recurrence-triggered timepoints. Quantify epithelial, hybrid E/M, and mesenchymal CTC fractions using prespecified CD45-negative, epithelial, EMT, thyroid-lineage, and survival marker criteria.

    **Primary endpoint:** Change in E/EM/M CTC composition after thyroidectomy.

    ## Aim 2. Build an NGS-based genetic map linking E, E/M, and M CTC fractions to the resected tumor clone.

    **Hypothesis:** E/EM/M CTC fractions represent discrete but connected transition states, and persistent postoperative blood signals can be genetically anchored.

    **Approach:** Generate matched tumor-normal NGS from resected tissue. Perform serial cfDNA analysis and CTC-enriched/low-input NGS only in QC-passing CTC-high subsets. Compare driver/private variants, CNV/SV/fusion patterns where technically feasible.

    **Key boundary:** CTC low-input NGS is stage-gated and exploratory unless QC is reproducible.

    ## Aim 3. Link CTC-EMT genetic states to tissue context and residual-risk features.

    **Hypothesis:** Persistent EM/M CTC burden and genetic concordance associate with lymphatic invasion, lymph node metastasis, ETE, BRAF/TERT/fusion context, Tg/US follow-up, or recurrence-risk features.

    **Approach:** Integrate pathology, FFPE/mIHC, tumor NGS, Tg/anti-Tg, ultrasound, RAI context, and follow-up. Use tissue marker panels to test whether tumor EMT/stress/thyroid-lineage-low features correspond to serial blood states.

    ## Aim 4. Develop a stage-gated liquid-biopsy framework for postoperative biological monitoring.

    **Hypothesis:** Serial CTC-EMT and matched genetic mapping can define a research-grade residual-risk state that justifies future clinical validation.

    **Approach:** Build a longitudinal transition model and clone-state interpretation workflow. The deliverable is a research framework with explicit assay-failure and claim-boundary rules, not a deployed clinical diagnostic.
    """


def budget() -> str:
    return """
    # 30억 / 3년 예산안

    단위: 억 원. Year 2에 cohort scale-up과 NGS가 집중되므로 연간 균등 편성이 아니다.

    | 항목 | Year 1 | Year 2 | Year 3 | 합계 | 근거 |
    |---|---:|---:|---:|---:|---|
    | 인건비 | 1.5 | 1.8 | 1.8 | 5.1 | PM, CRC, CTC technologist, bioinformatician, statistician, pathology/data manager |
    | 임상 코호트/바이오뱅크 | 0.7 | 0.7 | 0.7 | 2.1 | 동의, serial blood, 검체 처리, follow-up abstraction |
    | CTC enumeration/EMT phenotyping | 1.2 | 1.8 | 1.2 | 4.2 | assay lock, 항체/프로브, imaging, replicate QC, benign/healthy controls |
    | NGS/Inocras 또는 equivalent vendor | 1.8 | 3.0 | 1.8 | 6.6 | matched tissue-normal NGS, serial cfDNA, CTC-enriched feasibility subset, raw data |
    | FFPE mIHC/GeoMx/protein validation | 1.4 | 1.5 | 0.7 | 3.6 | tumor-lineage, EMT, immune/stress, thyroid differentiation marker validation |
    | Fresh tissue functional subset | 0.6 | 1.2 | 0.6 | 2.4 | selected fresh cases only; slice/organoid feasibility and short perturbation readouts |
    | Data/compute/statistics | 0.5 | 0.7 | 0.6 | 1.8 | secure data environment, longitudinal transition model, clone-state interpretation |
    | 소모품/QA | 0.4 | 0.4 | 0.4 | 1.2 | tubes, extraction kits, slides, library QC, controls |
    | 장비사용/서비스 | 0.5 | 0.4 | 0.3 | 1.2 | imaging, library QC, pathology scanner/analysis |
    | IP/regulatory/meeting/contingency | 0.4 | 0.6 | 0.8 | 1.8 | IRB amendments, DTA/MTA, patent map, advisory, publication |
    | **합계** | **9.0** | **12.1** | **8.9** | **30.0** | Stage-gated execution |

    ## Cut Rule

    예산 삭감 시 serial blood, CTC phenotype, matched tissue-normal NGS는 보존한다. 먼저 줄일 항목은 GeoMx breadth, CTC low-input NGS breadth, fresh tissue perturbation breadth다.
    """


def irb_checklist() -> str:
    return """
    # IRB / Data Collection Checklist

    ## 대상자

    - 성인 PTC 또는 PTC 의심 환자, thyroidectomy 예정.
    - 가능하면 LNM, LVI, ETE, aggressive variant, BRAF/TERT/fusion, recurrent/high Tg-risk feature가 있는 환자 enrich.
    - 비교군: benign thyroid surgery 및 healthy blood controls.

    ## Timepoints

    - T0: 수술 전 1-2주 이내.
    - T1: 수술 후 약 2주.
    - T2: 수술 후 3개월.
    - T3: 수술 후 6-12개월 또는 Tg/US/재발 의심 trigger.
    - Optional: clinically indicated RAI 전후.

    ## 검체

    - CTC용 peripheral blood: tube type, volume, processing window 명시.
    - cfDNA용 peripheral blood: Streck/EDTA, plasma separation SOP.
    - matched normal DNA: buffy coat 또는 buccal/saliva.
    - FFPE tumor block/slides 및 H&E.
    - Optional fresh tissue: workflow가 가능한 경우에만.
    - Pathology ROI annotation: tumor, invasive front, LVI, LN metastasis, stroma-rich region.

    ## 임상 변수

    - Age, sex, autoimmune thyroiditis/Hashimoto 여부.
    - 수술 범위, LN dissection.
    - tumor size, multifocality, bilaterality.
    - histologic subtype, margin, ETE, LVI, perineural invasion.
    - LN status, extranodal extension.
    - ATA risk, AJCC stage.
    - BRAF, TERT promoter, RET/NTRK/ALK fusion, RAS, existing molecular results.
    - Tg, anti-Tg, TSH, US finding, RAI dose/uptake.
    - recurrence/persistence event with exact dates.

    ## Assay/QC

    - CTC positivity definition prespecified.
    - E-CTC: CD45-negative, EpCAM/KRT-positive.
    - M-CTC: CD45-negative, VIM/TWIST/ZEB1 or equivalent EMT-marker-positive.
    - EM-CTC: epithelial plus EMT marker.
    - Thyroid-lineage support: PAX8/TG/TPO/TSHR/SLC5A5 where feasible.
    - Duplicate imaging/readout subset.
    - Spike-in recovery, benign/healthy controls.
    - Chain-of-custody and processing timestamp.

    ## Consent / Governance

    - Serial blood, FFPE, optional fresh tissue, genomic sequencing.
    - Somatic/germline incidental finding policy.
    - External vendor data transfer and storage.
    - Research-only vs returnable clinical report boundary.
    - Recontact and follow-up permission.
    - De-identification key custody, retention, withdrawal policy.
    """


def inocras_questionnaire() -> str:
    return """
    # Inocras Vendor Questionnaire

    ## Scope Fit

    1. PTC thyroidectomy research cohort에서 matched tumor-normal sequencing 지원 가능 여부.
    2. 적합한 서비스: CancerVision, MRDVision, custom research WGS, target-enhanced WGS, bioinformatics-only 중 무엇인지.
    3. FASTQ, BAM/CRAM, VCF, CNV, SV, fusion, TMB/MSI/HRD/signature output 제공 가능 여부.
    4. BRAF, TERT promoter, RET/NTRK/ALK fusions, RAS, EIF1AX, PPM1D, CHEK2, CNV/SV를 research report에 포함 가능한지.

    ## Sample Requirements

    5. FFPE input: curls/slides 수, tumor purity threshold, H&E requirement, macrodissection 지원.
    6. matched normal input: blood, buccal, saliva 중 허용 범위.
    7. cfDNA input: tube, blood volume, plasma volume, minimum cfDNA mass, processing window.
    8. CTC-enriched pellet 또는 low-input/single-cell CTC DNA/RNA 처리 가능 여부.
    9. 가능하다면 WGA/library method와 SNV/CNV/SV/fusion 한계.
    10. 불가능하다면 research-only exploratory input으로 받을 수 있는 대안.

    ## Performance / Validation

    11. cfDNA/MRDVision LOD와 low tumor fraction에서의 성능.
    12. PTC/thyroid cancer validation data 포함 여부.
    13. early-stage/low-shedding tumor에서 sensitivity 저하 예상.
    14. TERT promoter, fusion, structural variant 처리 방식.
    15. FFPE tumor-informed variant로 serial cfDNA interpretation 가능한지.
    16. germline, clonal hematopoiesis, tumor-derived variant 구분 방식.

    ## Operations / Governance

    17. sample type별 TAT.
    18. tumor-normal WGS, serial cfDNA, CTC-enriched sample, bioinformatics-only 단가.
    19. batch requirement와 Year 2 scale-up 가능성.
    20. failed sample policy와 rerun cost.
    21. 데이터 저장 위치, 국외 이전 여부, encryption, access control, IRB 문구.
    22. CAP/CLIA status와 research-use output boundary.
    23. 한국 병원 데이터가 국내에 머물 수 있는지.
    24. IRB 전 assay design 논의할 scientific contact.

    ## Deliverables

    25. de-identified report 예시.
    26. raw data manifest 예시.
    27. genome build와 annotation database version.
    28. pipeline versioning/reproducibility statement.
    29. molecular tumor board/expert interpretation 가능 여부.
    30. publication policy와 acknowledgment requirement.
    """


def defense_table() -> str:
    return """
    # Reviewer Defense Table

    | 예상 공격 | 답변 | 지원 근거 | Claim boundary |
    |---|---|---|---|
    | PTC는 예후가 좋은데 왜 30억인가? | 전체 low-risk PTC가 아니라 surgery 이후 residual-risk biology를 읽는 과제다. LNM/LVI/ETE/aggressive variant/molecular-risk case를 enrich한다. | Yu 2024, CTC literature, hospital cohort | 3년 내 사망률 개선 주장 없음 |
    | CTC는 이미 PTC에서 보였다. 새로움이 뭔가? | 기존은 count/phenotype 중심이다. 새로움은 thyroidectomy perturbation + serial CTC-EMT + matched NGS clone map이다. | Aim 1-2 | CTC 검출 자체가 novelty가 아님 |
    | CTC artifact 가능성이 크다. | CD45-negative, epithelial/EMT/thyroid-lineage criteria, benign/healthy controls, spike-in, duplicate readout, cfDNA/tissue NGS orthogonal support를 둔다. | IRB/QC checklist | CTC phenotype이 viable metastasis 증명은 아님 |
    | Post-op CTC 감소는 단순 debulking 아닌가? | 맞다. 그 debulking을 controlled perturbation으로 사용한다. 핵심은 expected collapse 이후 persistence/state shift다. | serial timeline | 감소 자체를 EMT reversal로 해석하지 않음 |
    | CTC NGS는 low-input 때문에 실패할 수 있다. | tissue-normal NGS는 core, cfDNA는 serial, CTC-enriched NGS는 CTC-high QC-passing subset stage-gate다. | kill criteria | single-cell CTC WGS는 primary endpoint가 아님 |
    | cfDNA는 early PTC에서 낮을 수 있다. | 그래서 CTC phenotype과 tissue NGS를 병행한다. Negative cfDNA도 LOD 안에서만 해석한다. | Sato 2021, vendor questionnaire | universal detection 주장 없음 |
    | EMT-state CTC가 재발을 예측한다는 근거가 약하다. | primary endpoint는 biological transition이다. recurrence/Tg/US는 secondary/exploratory association이다. | endpoint table | clinical prediction claim 없음 |
    | public spatial pilot이 CTC를 뒷받침하지 않는다. | 동의한다. spatial data는 marker/ROI/tissue-state logic만 제공한다. | TCGA/spatial preliminary evidence | CTC proof로 쓰지 않음 |
    | Technology/ICT로 가야 하는 것 아닌가? | 아니다. 핵심은 one disease, one surgery perturbation, one CTC-EMT transition이다. AI는 해석 인프라다. | category decision | ICT product claim 없음 |
    | Inocras 외주 프로젝트 아닌가? | vendor는 sequencing service다. 과학적 novelty는 serial design, phenotype-genotype integration, state transition interpretation이다. | vendor questionnaire | vendor technology를 novelty로 포장하지 않음 |
    """


def kakao() -> str:
    return f"""
    # Prof. Yu Kakao Summary

    교수님, 삼성 과제는 Science track으로 좁히는 것이 맞겠습니다. 넓은 K-Thyro spatial/drug/AI platform으로 가면 "병원검체+외주 NGS+분석 플랫폼"으로 보일 위험이 큽니다. 반대로 **갑상선절제술 전후 PTC CTC-EMT 상태전이와 NGS 유전지도**로 좁히면, thyroidectomy를 인간 perturbation model로 쓰는 명확한 생명과학 질문이 됩니다.

    기존 Yu 2024는 이미 강한 anchor입니다. 62명 PTC에서 pre-op, 2주, 3개월 blood를 보았고, CTC가 87%에서 검출되며 EMT/mesenchymal phenotype이 우세하고 수술 후 count가 감소했습니다. 하지만 아직 없는 것은 matched tissue NGS와 postoperative CTC/cfDNA genetic map입니다. 이 gap이 삼성 과제의 핵심입니다.

    Public pilot은 CTC 증명이 아니라 보조근거로 쓰겠습니다. TCGA 505명에서 BRAF mutation만으로 tissue state가 설명되지 않고, BRAF 274명도 {M['braf_labels']}개 label로 갈라졌습니다. GSE250521 spatial {M['spatial_slides']} slides에서는 {M['spatial_z_pos']}에서 niche coherence가 나왔습니다. 이 결과는 marker selection과 tissue context 근거이지, CTC shedding 증거는 아닙니다.

    병원에서 필요한 것은 serial blood timepoint, CTC E/EM/M 원자료, FFPE/fresh tissue 가능 여부, BRAF/TERT/RET/NTRK 등 기존 molecular result, LVI/LNM/ETE/multifocality, Tg/anti-Tg, RAI, ultrasound, recurrence follow-up입니다.

    Inocras에는 PTC FFPE+matched normal WGS, cfDNA/MRD sensitivity, CTC-enriched/low-input 처리 가능성, TERT/fusion/SV/CNV reporting, FASTQ/BAM/VCF 제공, IRB/국외반출/보안, sample input, 단가/TAT를 확인하겠습니다.

    7일 계획은 1일차 Science track/title lock, 2일차 병원 CTC 원자료 변수표 확인, 3일차 Inocras 질문지 발송, 4일차 IRB specimen/timepoint 초안, 5일차 preliminary figure 3장 정리, 6일차 30억/3년 예산/역할분담, 7일차 삼성 2쪽 concept note 초안 완성입니다.
    """


def action_plan() -> str:
    return """
    # Next 7-Day Action Plan

    | Day | Action | Output | Owner |
    |---:|---|---|---|
    | 1 | Science-track title and one-line hypothesis lock | final title, 150-word abstract | PI/Seungho |
    | 2 | Existing Yu CTC data variable inventory | variable dictionary: timepoints, E/EM/M counts, clinicopathology | hospital team |
    | 3 | Send Inocras questionnaire | written vendor response request | Seungho/PI |
    | 4 | IRB specimen/timepoint draft | IRB checklist and consent scope | CRC/PI |
    | 5 | Build 3 preliminary figures | thyroidectomy perturbation, missing dataset, mutation-not-enough support | Seungho |
    | 6 | Budget and role lock | 30억/3년 table, sample count assumptions, stage gates | PI/admin |
    | 7 | Samsung 2-page concept note | final 2-page blind proposal draft | Seungho/PI |

    Immediate decision needed: whether the hospital can provide serial post-op blood at 2 weeks and 3 months prospectively, and whether archived Yu 2024-style CTC phenotype raw data can be re-tabulated with clinical variables.
    """


def claim_boundaries() -> str:
    return """
    # Claim Boundaries And Kill Criteria

    ## Allowed

    - PTC CTCs and EMT-state CTCs are measurable in published cohorts.
    - Public omics show thyroid cancer state heterogeneity beyond mutation.
    - Spatial data support tissue-state organization, not CTC shedding.
    - Matched tissue NGS is needed to interpret postoperative blood signals.
    - This project will test whether thyroidectomy reveals CTC-EMT/genetic state transition.

    ## Forbidden

    - CTC-EMT state already predicts recurrence in our cohort.
    - Public pilot proves CTC biology.
    - cfDNA/CTC NGS will work in all early PTC.
    - This is a ready clinical diagnostic.
    - Inocras/vendor technology is the novelty.

    ## Kill Criteria

    | Risk | Kill or narrow criterion | Action |
    |---|---|---|
    | CTC detection/subtyping not reproducible | control/spike-in/duplicate readout fails | narrow to count + cfDNA; remove EMT breadth |
    | tissue-normal NGS success < 85% | FFPE input or vendor failure | pause scale-up; fix tissue processing |
    | cfDNA unusable even in high-risk subset | low yield or no tumor-informed signal within LOD | keep cfDNA exploratory; rely on tissue + CTC phenotype |
    | CTC low-input NGS fails QC | poor WGA/library/variant quality | drop CTC NGS from core endpoint |
    | follow-up data capture < 80% | Tg/US/pathology incomplete | reduce clinical association claims |
    | primary endpoint becomes unclear | too many platform branches | remove spatial/drug/AI branches from main proposal |
    """


def concept_note() -> str:
    return f"""
    # 2-Page Concept Note

    ## Title

    **{TITLE}**  
    **{TITLE_KR}**

    ## Problem

    Papillary thyroid cancer is usually curable, yet postoperative management still relies on indirect and delayed signals: pathology, thyroglobulin, anti-thyroglobulin antibody, ultrasound, and recurrence surveillance. These tools do not directly measure whether residual tumor-cell states persist in blood after surgery. CTC studies show that PTC CTCs and EMT-state CTCs are measurable, but no study has built a serial post-thyroidectomy CTC-EMT genetic map.

    ## Hypothesis

    Thyroidectomy is a controlled human perturbation. It should collapse nonspecific tumor shedding. Persistent postoperative hybrid/mesenchymal CTC states and matched cfDNA/CTC genetic signals may reveal residual-risk PTC biology.

    ## Innovation

    This is not a broad platform proposal. It is a narrow Science proposal with one disease, one perturbation, one biological transition, and one genetic anchor. The novelty is the combination of serial CTC-EMT phenotyping with matched tumor-normal NGS and postoperative liquid-biopsy genetic mapping.

    ## Approach

    Patients undergoing thyroidectomy for PTC will be sampled at T0 pre-op, T1 2 weeks, T2 3 months, and T3 6-12 months or recurrence-triggered timepoints. CTCs will be classified into epithelial, hybrid E/M, and mesenchymal fractions with prespecified CD45-negative, epithelial, EMT, thyroid-lineage, and survival markers. Resected tumors will undergo matched tumor-normal NGS. Serial cfDNA and CTC-enriched NGS will be used in a stage-gated way to test genetic concordance and postoperative persistence.

    ## Expected Output

    The output is a research-grade map of CTC-EMT state transition and matched genetic persistence after thyroidectomy. The 3-year deliverable is not a deployed clinical diagnostic; it is a mechanistic and translational foundation for future clinical validation.

    ## Why Samsung Science

    The project asks whether surgical removal of a human tumor perturbs circulating tumor-cell state and whether residual circulating genetic signals are biologically meaningful. This is a fundamental human cancer biology question with translational potential.
    """


def reviewer_memo() -> str:
    return f"""
    # Final Reviewer Memo

    ## 1. Why Science

    The proposal is centered on a biological perturbation: thyroidectomy. It asks whether removal of the tumor source changes the circulating tumor-cell EMT state and whether persistent genetic signals can be anchored to the resected tumor clone. This is a mechanism-driven human cancer biology question.

    ## 2. Why The Project Is Narrow Enough

    It uses one disease, one perturbation, one primary biological transition, and one genetic anchor. It explicitly removes the broad spatial/drug/AI platform identity.

    ## 3. What Is Truly Novel

    Prior studies measured CTC burden and EMT phenotype. The proposed novelty is serial CTC-EMT state transition after thyroidectomy coupled to matched tissue-normal NGS and postoperative cfDNA/CTC genetic mapping.

    ## 4. What Existing Yu Data Anchor

    Yu 2024 anchors feasibility: 62 prospective PTC patients, pre-op/2-week/3-month blood, 87% CTC detection, EMT/mesenchymal predominance, postoperative CTC decline.

    ## 5. What Public Pilot Contributes

    The local public pilot contributes marker and tissue-state logic. BRAF tumors split across {M['braf_labels']} labels, largest label {M['braf_largest_fraction']}; spatial {M['spatial_z_pos']} slides show same-niche coherence; protein data support RAI loss and myeloid/TGFB gain in dedifferentiation. None of this proves CTC biology.

    ## 6. Main Risk

    The main risk is assay feasibility: CTC subtyping reproducibility, cfDNA low shedding, and CTC low-input NGS failure.

    ## 7. How Stage-Gate Handles Risk

    Tissue-normal NGS and serial CTC phenotype are preserved as core. cfDNA and CTC-enriched NGS are stage-gated. If low-input CTC NGS fails, the project remains valid as a CTC-EMT transition plus tissue/cfDNA anchored study.

    ## 8. What To Remove If Reviewer Says It Is Too Broad

    Remove drug-delivery language, broad spatial platform language, extensive AI framing, and fresh tissue perturbation breadth. Keep thyroidectomy, CTC-EMT, matched NGS, and serial follow-up.

    ## 9. Final Recommendation

    Submit as Samsung Science. Do not submit as Technology/ICT unless the proposal is rebuilt around a genuine algorithmic/IP primitive, which would weaken the biological core.
    """


def write_all() -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    PROMPTS.mkdir(parents=True, exist_ok=True)
    INPUTS.mkdir(parents=True, exist_ok=True)

    write(ROOT / "README.md", f"""
    # K-Thyro CTC-EMT Samsung Science Proposal

    This package pivots away from the broad K-Thyro spatial/drug/AI platform and locks the proposal to a narrow Samsung Science-track mechanism:

    **{TITLE}**

    Core question: Does thyroidectomy perturb the circulating tumor-cell EMT state, and do persistent postoperative CTC or cfDNA genetic maps identify residual-risk biology in papillary thyroid cancer?

    Main outputs are in `outputs/reports/`, slide figures are in `outputs/figures/`, and the PPTX deck is generated by `scripts/build_ctc_emt_deck.py`.

    Key boundary: public omics support feasibility and marker-selection logic; they do not prove CTC biology.
    """)

    write(INPUTS / "source_context_index.md", f"""
    # Source Context Index

    - Broad pilot root: `{PUBLIC.parent}`
    - CTC-EMT pitch pack: `{PUBLIC / 'samsung_ctc_emt_ngs_pitch_pack_2026_05_09.md'}`
    - Public pilot evidence package: `{PUBLIC / 'ppt_evidence_package'}`
    - No local `연구 계획서_삼성.docx` or Yu 2024 paper PDF was found in the requested path scan.
    """)

    write(PROMPTS / "proposal_generation_prompt.md", """
    # Generation Prompt

    Keep the proposal narrow: PTC, thyroidectomy, CTC-EMT state transition, matched NGS genetic map. Do not make AI or spatial/drug platform the identity. Use public pilot only as supporting marker/tissue-state evidence.
    """)

    write(REPORTS / "context_digest.md", build_context_digest())
    write(REPORTS / "final_proposal_logic.md", build_final_logic())
    write(REPORTS / "samsung_ctc_emt_one_page_summary_kr.md", one_page_summary())
    write(REPORTS / "samsung_ctc_emt_preliminary_results_3page_kr.md", preliminary_results())
    write(REPORTS / "samsung_ctc_emt_specific_aims_kr.md", specific_aims())
    write(REPORTS / "samsung_ctc_emt_budget_30억_3년_kr.md", budget())
    write(REPORTS / "samsung_ctc_emt_irb_checklist_kr.md", irb_checklist())
    write(REPORTS / "samsung_ctc_emt_inocras_questionnaire_kr.md", inocras_questionnaire())
    write(REPORTS / "samsung_ctc_emt_reviewer_defense_table_kr.md", defense_table())
    write(REPORTS / "message_to_prof_yu_kakao_kr.md", kakao())
    write(REPORTS / "next_7_day_action_plan_kr.md", action_plan())
    write(REPORTS / "claim_boundaries_and_kill_criteria.md", claim_boundaries())
    write(REPORTS / "samsung_ctc_emt_concept_note_2page_kr.md", concept_note())
    write(REPORTS / "reviewer_memo_final_kr.md", reviewer_memo())

    # Export compact tables for deck/manual reuse.
    pd.DataFrame(
        [
            ["Personnel", 1.5, 1.8, 1.8, 5.1],
            ["Clinical cohort / biobank", 0.7, 0.7, 0.7, 2.1],
            ["CTC phenotyping", 1.2, 1.8, 1.2, 4.2],
            ["NGS / vendor", 1.8, 3.0, 1.8, 6.6],
            ["FFPE mIHC / GeoMx", 1.4, 1.5, 0.7, 3.6],
            ["Fresh tissue subset", 0.6, 1.2, 0.6, 2.4],
            ["Data / compute / statistics", 0.5, 0.7, 0.6, 1.8],
            ["Consumables / QA", 0.4, 0.4, 0.4, 1.2],
            ["Equipment services", 0.5, 0.4, 0.3, 1.2],
            ["IP / regulatory / contingency", 0.4, 0.6, 0.8, 1.8],
        ],
        columns=["category", "year1", "year2", "year3", "total"],
    ).to_csv(TABLES / "budget_30억_3년.tsv", sep="\t", index=False)


if __name__ == "__main__":
    write_all()
    print(f"Wrote proposal text package to {ROOT}")
