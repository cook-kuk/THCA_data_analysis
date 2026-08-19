# Samsung CTC-EMT-NGS Pitch Pack

Date: 2026-05-09  
Working directory: `/home/seungho/personal/THCA_data_analysis/kthyro_public_pilot`  
Core constraint: **PTC CTC-EMT state transition and NGS genetic map after thyroidectomy**  

## Executive Verdict

**Submit to Science, not Technology/ICT.**

Brutal reviewer read: the broad K-Thyro spatial/drug platform is interesting but too diffuse for a 30억/3년 Samsung pitch. If framed as Technology/ICT, it looks like a hospital cohort plus vendor sequencing plus analysis dashboard. That is not enough as an ICT primitive. The defensible center is a biological question:

> Does thyroidectomy perturb the circulating tumor-cell EMT state, and do persistent/postoperative CTC or cfDNA genetic maps identify the residual-risk clone in papillary thyroid cancer?

Use ICT/AI as analysis infrastructure, not as the category identity. Use the public pilot as feasibility and marker-selection evidence, not as direct proof of CTC biology.

## 0. Public Data And Literature Scan

### Direct CTC/EMT Thyroid Evidence

| Source | What it supports | What it does not support | Use in proposal |
|---|---|---|---|
| Yu et al., Int J Surg 2024, prospective PTC thyroidectomy CTC study, 62 PTC patients; blood before surgery, 2 weeks, 3 months; CTC detected in 87%; EMT/mesenchymal phenotypes predominant; counts decreased after thyroidectomy | PTC has measurable CTCs; CTC-EMT state changes after thyroidectomy are biologically and operationally plausible | No CTC NGS genetic map; no clone tracking; not a mechanism proof | Opening rationale and hospital-data justification |
| Li et al., Mol Clin Oncol 2022, 394 thyroid cancer patients; 270 PTC; CanPatrol CTC subtyping; CTC positivity about 95%; >6 CTCs and mesenchymal CTCs associated with poorer outcomes in differentiated thyroid cancer | EMT-state CTCs have prognostic signal in thyroid cancer | Mixed histologies, observational, not thyroidectomy-serial NGS | Reviewer defense: CTC/EMT signal is not invented |
| Lin et al., Thyroid 2018, CEC/CTC markers EpCAM/PDPN/TSHR correlated with remission/survival | Thyroid-origin circulating epithelial cells have clinical association | Not NGS, not EMT transition, not postoperative clone map | Background only |
| TUMORFISHER PTC 2026 prospective report, 210 PTC patients, preoperative CTC risk stratification | Preoperative CTC can stratify PTC risk in a recent prospective cohort | Not serial thyroidectomy transition; platform-specific | Secondary support, not core |

### Liquid Biopsy / NGS Evidence

| Source | What it supports | What it does not support | Use in proposal |
|---|---|---|---|
| Sato et al., Thyroid 2021, 22 primary PTC patients, BRAF V600E ctDNA before/after surgery by ddPCR; presurgery ctDNA detected in 5/16 BRAF-mutant tumors; one postsurgery-positive patient recurred | Postoperative blood mutation persistence can carry residual-risk signal | Single mutation/ddPCR only; small N; not CTC; not genome-wide | Justifies serial blood NGS and postoperative timepoints |
| Tarasova et al., Thyroid 2024, Guardant360 plasma NGS in 1094 thyroid cancer samples | Plasma NGS can detect actionable thyroid cancer alterations across subtypes | Mostly advanced/real-world submitted samples; not early PTC thyroidectomy kinetics | Supports feasibility of plasma NGS, not proof of early PTC monitoring |
| TCGA-THCA Cell 2014, 496 PTC genomic landscape | PTC has driver classes, fusions, differentiation states; tissue NGS map is necessary | No CTC or postoperative liquid biopsy | Tissue-genetic reference map |

### Public Omics Used In Local Pilot

| Local/public resource | Local result | Use after narrowing |
|---|---|---|
| TCGA-THCA bulk, 505 primary tumors | 7 vulnerability labels; BRAF tumors split across all 7 labels; driver group explains only part of axis variance | Marker/axis selection; argument that mutation alone is insufficient |
| GSE250521 spatial transcriptomics, 16 slides, 57,144 spots | 16/16 slides same-niche coherence z > 2; therapeutic states are spatially organized | Tissue context for where CTC-shedding states may originate; ROI selection |
| External GPL570 thyroid cohorts, 206 samples | 168 contrast tests; 65 strong, 12 moderate; ATC/advanced contrasts reproduce RAI-low/aggressive/barrier-high directions | Cross-cohort reproducibility of tissue-state axes |
| Bulk thyroid proteomics proxy, 461 samples | RAI protein axis down in ATC; myeloid/TGFB barrier axis up in ATC | Protein-level directionality, not spatial protein proof |

### Gap That Becomes The Proposal

There is **no strong public dataset** that jointly provides:

- PTC serial blood before and after thyroidectomy,
- CTC EMT phenotype transition,
- matched tissue NGS,
- CTC-enriched or cfDNA NGS,
- postoperative residual/recurrent outcome.

That absence is the proposal’s opening. Do not pretend the pilot already proves it.

## 1. Science vs Technology/ICT Decision

| Criterion | Science | Technology/ICT |
|---|---:|---:|
| Originality of biological question | High: surgery as a controlled perturbation of CTC-EMT state | Medium: analysis platform is derivative unless new assay/algorithm is invented |
| Fit to Samsung language | Strong: life-science mechanism and new principle | Weak-to-medium: healthcare ICT fusion possible, but reviewer will ask what is the core ICT invention |
| Defensibility with current evidence | Strong: public CTC and local omics support feasibility | Weak: no proprietary CTC device, no validated AI product |
| Risk of looking like service integration | Lower | High |
| Best reviewer hook | "Thyroidectomy reveals the liquid-biopsy state transition of residual PTC clones" | "Multimodal liquid-biopsy platform" sounds overbroad |

**Final category: Science.**  
Use Technology/ICT only as a secondary capability: data integration, clone-state modeling, and clinical-grade reporting.

## 2. Title Rewrites

### Science-Style Narrow Mechanism Titles

1. **Serial Genetic Mapping of CTC-EMT State Transitions After Thyroidectomy in Papillary Thyroid Cancer**
2. Thyroidectomy as a Human Perturbation Model for CTC-EMT Plasticity in Papillary Thyroid Cancer
3. Post-Thyroidectomy Persistence of EMT-State Circulating Tumor Cells and Genetic Residual Disease in Papillary Thyroid Cancer
4. Coupling CTC-EMT Phenotype Switching to Tumor Genotype After Surgical Removal of Papillary Thyroid Cancer
5. Defining the Liquid-Biopsy Transition State of Papillary Thyroid Cancer Through CTC-EMT Phenotyping and NGS Clone Mapping

### Technology/ICT-Style Platform Titles

1. Surgical Liquid-Biopsy Atlas Platform for CTC-EMT and NGS Clone Tracking in Papillary Thyroid Cancer
2. A Multimodal CTC-EMT and Genome-Mapping Platform for Post-Thyroidectomy PTC Monitoring
3. AI-Integrated Serial Liquid Biopsy for CTC State Transition and Genetic Residual Mapping in PTC
4. CTC-to-Genome Digital Twin Platform for Thyroidectomy-Guided Papillary Thyroid Cancer Stratification
5. Precision Thyroidectomy Monitoring Platform Using CTC-EMT Phenotyping, cfDNA, and Target-Enhanced WGS

### Best Final Title

**Serial Genetic Mapping of CTC-EMT State Transitions After Thyroidectomy in Papillary Thyroid Cancer**

Korean working title:

**갑상선절제술 전후 유두갑상선암 CTC-EMT 상태전이와 유전지도의 연속 해독**

Why this title wins: it is narrow, testable, and honest. It does not claim prediction, cure, AI product, or platform completion before data exist.

## 3. 12-Slide Samsung Pitch Deck Outline

| Slide | Title | One-line message | Main support / figure |
|---:|---|---|---|
| 1 | The Narrow Problem | PTC is usually curable, but current postoperative risk tools cannot read the real-time residual tumor-cell state. | Clinical schematic; no overclaim |
| 2 | Core Hypothesis | Thyroidectomy is a controlled perturbation that should collapse benign tumor shedding, while persistent EMT-state CTCs/genetic clones mark residual-risk biology. | Concept diagram |
| 3 | Why Now | Prior PTC CTC studies show CTCs and EMT phenotypes are measurable, but no one has paired serial CTC-EMT transition with NGS clone mapping. | Public literature table |
| 4 | Public Pilot GO | Local pilot shows PTC/thyroid cancer states are not captured by driver mutation alone. | `results/ppt_evidence_package/figures/tcga_mutation_not_enough_driver_label_fraction.png` |
| 5 | Tissue States Are Spatial | GSE250521: 16/16 slides show non-random therapeutic niche coherence. | `results/ppt_evidence_package/figures/spatial_coherence_barplot_dark.png` |
| 6 | Mechanistic Bridge To CTC | EMT/high-risk tissue states nominate markers for CTC phenotyping: epithelial, hybrid, mesenchymal, thyroid-lineage, immune-barrier context. | Marker table: EpCAM/KRT/PAX8/TG/VIM/TWIST/ZEB1/CD45-negative |
| 7 | Study Design | Prospective serial sampling around thyroidectomy: T0 pre-op, T1 2 weeks, T2 3 months, T3 6-12 months/high-risk recurrence trigger. | Cohort timeline |
| 8 | NGS Genetic Map | Matched tumor-normal NGS anchors the clone; cfDNA and CTC-enriched/low-input NGS test postoperative persistence. | Tissue-blood-clone map |
| 9 | Primary Endpoints | Primary: CTC-EMT state transition after surgery. Secondary: concordance with tissue genotype/cfDNA, Tg/US/recurrence-risk features. | Endpoint matrix |
| 10 | What We Will Not Claim | This is not a clinical deployment test, not proof of metastasis, not a vendor product, and not a replacement for pathology/Tg. | Reviewer boundary slide |
| 11 | 30억/3년 Execution | Year 1 assay lock and pilot cohort; Year 2 scale serial cohort and NGS; Year 3 validation model and mechanistic closure. | Budget/milestone waterfall |
| 12 | Final Ask | Fund a Science project that turns thyroidectomy into a human perturbation experiment for liquid-biopsy cancer-state transitions. | One-sentence close |

Design direction: dense scientific deck, not marketing. Use a restrained clinical palette, one clear schematic per slide, and keep local public-pilot figures as "supporting feasibility" rather than the title story.

## 4. Three-Page Preliminary Results Section

### Preliminary Result 1. Public evidence supports measurable CTC-EMT states in thyroid cancer, but not yet genetic clone mapping.

Published thyroid cancer CTC studies provide a realistic entry point for the proposed work. In a prospective PTC thyroidectomy cohort, peripheral blood CTCs were measured before surgery and after thyroidectomy, and CTCs were detected in most patients. The key biological signal was not just CTC count but the predominance of epithelial-mesenchymal or mesenchymal phenotypes. CTC burden decreased after thyroidectomy, particularly in clinicopathologic high-risk contexts such as lymphatic invasion, lymph node metastasis, and BRAF V600E mutation. A separate thyroid cancer cohort using CTC subtype analysis also reported high CTC positivity and associated mesenchymal CTC burden with worse outcomes in differentiated thyroid cancer.

These studies justify the feasibility of serial CTC phenotyping in PTC. They do not answer the proposed question. Existing data do not show whether postoperative CTC persistence represents the same genetic clone as the resected tumor, whether EMT-state switching is concordant with tumor genotype, or whether a residual liquid-biopsy clone can be distinguished from nonspecific postoperative cell release. This is the exact gap addressed here.

The proposal therefore starts from a conservative premise: CTC-EMT state is a measurable biological phenotype, not yet a validated clinical decision test. The proposed advance is to convert a count-and-marker observation into a serial perturbation map: before thyroidectomy, after surgical removal, and during early postoperative surveillance, with matched tissue NGS anchoring the genetic identity of the signal.

### Preliminary Result 2. Local public-data pilot shows that PTC therapeutic states cannot be reduced to driver mutation alone.

The local K-Thyro public-data pilot analyzed TCGA-THCA primary tumors and public spatial transcriptomics to test whether thyroid cancer states are separable beyond average mutation or bulk expression. TCGA-THCA analysis included 505 primary tumors and produced 7 patient-level vulnerability labels: RAI-readable differentiated (127), Mixed/Other (125), Drug-delivery barrier-high (105), HLA-visible inflamed (44), CD8-excluded myeloid/CAF-high (42), HLA-low immune-invisible (32), and RAI-low dedifferentiated (30).

The driver-mutation audit is the most important defense for the narrowed CTC proposal. BRAF-mutant tumors were not homogeneous: 274 BRAF tumors split across all 7 vulnerability labels, and the largest BRAF label represented only 33.2% of BRAF cases. Driver group explained only part of the variance across axes: 30.9% for RAI differentiation, 15.6% for HLA-I/APM, 12.5% for immune visibility, and 16.2% for CD8 exclusion. Even after descriptive adjustment for driver, stage, epithelial-lineage score, age, and sex, large unexplained fractions remained, and residual axis values still separated vulnerability labels.

This result does not prove CTC biology. Its value is narrower and stronger: it shows that tissue genotype alone is unlikely to capture the state information we want from serial blood. A postoperative NGS map without CTC-EMT state could miss biology; CTC-EMT phenotyping without tissue genotype could be nonspecific. The proposed project combines both because the public pilot says neither layer is sufficient alone.

### Preliminary Result 3. Spatial and protein-level public data support tissue-state context for selecting CTC markers and validation ROIs.

GSE250521 spatial transcriptomics analysis included 16 thyroid slides and 57,144 spots across normal, PTC, locally advanced PTC, and ATC conditions. Applying the same K-Thyro state-scoring framework, all 16 slides showed same-niche spatial coherence with z-score > 2. The minimum, median, and maximum z-scores were 6.09, 22.00, and 43.16. This means the tissue states are not random spot noise. However, spots are nested within slides, and the result must be interpreted as spatial organization of expression states, not as clinical outcome prediction.

Additional analyses strengthened the tissue-state interpretation. Spatial hotspot analysis supported biologically interpretable patterns such as delivery-failure proxy co-localizing with hypoxia-high territories and HLA-high co-localizing with IFN/CD8-high territories. Independent bulk thyroid proteomics across 461 samples showed that the RAI differentiation protein axis decreases in ATC, while myeloid/TGFB barrier protein axes increase with dedifferentiation. This is not public spatial proteomics proof. It is protein-level directionality that helps select the markers to validate in FFPE and CTC assays.

For the CTC-EMT project, this preliminary evidence provides the marker logic: epithelial/thyroid-lineage markers for tumor identity, EMT markers for state transition, and tissue NGS/cfDNA for genetic anchoring. It also defines the claim boundary. The project will not claim that spatial transcriptomics proves CTC shedding. Instead, tissue-state maps will guide where tumor EMT/barrier/dedifferentiation states may originate and which FFPE/mIHC features should be correlated with serial CTC changes.

### Preliminary Result 4. The strongest proposal is a stage-gated serial cohort, not a broad platform claim.

The public pilot supports a GO decision only if the proposal is narrowed. A broad theranostic spatial-drug platform would invite reviewer attacks: too many axes, no direct CTC data, no clinical deployment evidence, and vendor-dependent NGS. The narrowed project can be defended because it has one perturbation, one disease, one serial liquid-biopsy phenomenon, and one genetic anchoring strategy.

The proposed preliminary workflow is:

1. Enroll PTC patients undergoing thyroidectomy with enriched high-risk features where possible.
2. Collect blood at preoperative baseline, 2 weeks, 3 months, and later surveillance or recurrence-triggered timepoints.
3. Quantify epithelial, hybrid epithelial-mesenchymal, and mesenchymal CTC states using a prespecified marker panel and strict CD45-negative/thyroid-lineage criteria.
4. Generate matched tumor-normal NGS maps from resected tissue.
5. Test whether postoperative CTC/cfDNA signals match the tissue clone and whether persistent EMT-state CTCs are associated with residual-risk features.
6. Validate tissue origin and state with FFPE/mIHC and selected spatial/ROI assays.

The defensible claim after 3 years is not "clinical test ready." The defensible claim is:

> Thyroidectomy reveals a measurable liquid-biopsy state transition in PTC, and matched NGS can distinguish postoperative genetic persistence from nonspecific CTC/CTC-like signals.

That is a Science claim.

## 5. Kakao Message To Prof. Yu

교수님, 삼성 과제 방향을 다시 좁혀보니 넓은 "공간오믹스+drug platform"보다 **"갑상선절제술 전후 PTC CTC-EMT 상태전이와 NGS 유전지도"**가 훨씬 방어 가능해 보입니다.

이유는 명확합니다. 넓게 가면 심사위원이 "기존 병원검체+외주 NGS+분석 플랫폼 아닌가?"라고 공격할 가능성이 큽니다. 반대로 CTC-EMT를 수술 전후 serial로 보고, 절제 조직 NGS와 postoperative blood signal을 매칭하면 **thyroidectomy를 인간 perturbation model로 쓰는 생명과학 질문**이 됩니다.

Public pilot에서 확인한 것은 직접 CTC 증명이 아니라 보조 근거입니다. TCGA 505명에서 BRAF mutation만으로는 RAI/immune/barrier/aggressive state가 설명되지 않았고, BRAF 274명도 7개 vulnerability label로 갈라졌습니다. GSE250521 spatial 16개 slide에서는 16/16에서 niche coherence가 나와서 tissue state가 무작위가 아니라는 점은 지지됩니다. 다만 이건 CTC를 증명하지 않으므로 claim boundary를 낮추겠습니다.

병원에서 필요한 데이터는 다음입니다: 수술 전/후 혈액 timepoint(가능하면 pre-op, 2주, 3개월, 6-12개월), CTC count와 E/EM/M phenotype 원자료, FFPE/fresh tumor block, BRAF/TERT/RET/NTRK 등 기존 병리/분자결과, lymphatic invasion/LN metastasis/ETE/multifocality, Tg/anti-Tg, RAI 여부, 초음파/재발 follow-up입니다.

Inocras에는 다음을 확인하려고 합니다: PTC FFPE+matched normal WGS/target-enhanced WGS 가능 범위, cfDNA/MRDVision의 PTC low-shedding sensitivity, CTC-enriched fraction 또는 single/low-input CTC DNA/RNA 처리 가능 여부, fusion/SV/CNV/TERT/BRAF/RET/NTRK 보고 범위, FASTQ/BAM/VCF 제공 여부, IRB/국외반출/데이터보안 조건, sample input과 단가/TAT입니다.

7일 계획은 1일차 제목/Science 트랙 확정, 2일차 교수님 기존 CTC 원자료 변수표 확인, 3일차 Inocras 질문지 발송, 4일차 IRB specimen/timepoint 초안, 5일차 preliminary figure 3장 정리, 6일차 30억/3년 예산표와 역할분담, 7일차 삼성 2쪽 제안서 초안으로 가겠습니다.

핵심 문장은 이렇게 잡겠습니다:  
**"갑상선절제술 전후 유두갑상선암 CTC-EMT 상태전이와 유전지도의 연속 해독"**

## 6. Reviewer Defense Table

| Expected attack | Answer | Figure/table support | Claim boundary |
|---|---|---|---|
| PTC is too indolent. Why spend 30억? | The target is not all low-risk PTC; it is residual-risk biology after surgery and avoidance of blind overtreatment/undertreatment. Enrich for LN+, LVI+, ETE+, aggressive histology, BRAF/TERT/fusion or Tg-risk cases. | Public CTC studies; hospital cohort table to be built | We do not claim mortality reduction in 3 years |
| CTC in PTC has already been shown. What is new? | Prior work measured count and EMT phenotype. The new layer is serial thyroidectomy perturbation plus matched NGS clone map. | Yu 2024 + proposed NGS workflow | Novelty is not CTC detection alone |
| CTC assay artifacts are common. | Use strict CD45-negative, thyroid-lineage positive, epithelial/EMT marker rules; include healthy/benign thyroid controls, spike-in controls, duplicate readers, and orthogonal cfDNA/tissue NGS. | IRB checklist; assay QC slide | CTC phenotype is a measured assay state, not automatically viable metastatic cell |
| Postoperative CTC decrease is just tumor debulking. | Correct. The project uses debulking as a controlled perturbation; the key signal is persistence or state shift after the expected decline. | Serial timepoint design | We do not claim all decline is mechanistic EMT reversal |
| NGS from CTCs may fail due to low input. | Stage-gate: tissue-normal NGS for all selected cases; cfDNA for serial blood; CTC-enriched/low-input NGS only in CTC-high subset after feasibility lock. | Inocras questionnaire; budget stage gates | Single-cell CTC WGS is exploratory, not guaranteed primary endpoint |
| cfDNA may be low in early PTC. | That is why the design combines CTC phenotype, cfDNA, and tissue NGS. Negative cfDNA is informative only with assay LOD and tumor-informed context. | Sato 2021 small PTC ctDNA result | We do not promise universal ctDNA detection |
| EMT-state CTC does not prove recurrence. | The endpoint is state transition and genetic concordance first; recurrence/Tg/US are secondary exploratory associations. | Endpoint matrix | No clinical deployment claim |
| Public spatial/drug pilot does not support CTC. | Agreed. It supports tissue-state marker logic and driver insufficiency, not CTC shedding. | TCGA driver insufficiency; GSE250521 coherence | Spatial evidence is auxiliary |
| This is Technology/ICT with no technology. | We submit as Science. The central discovery is biological state transition after surgery. | Track decision slide | AI is analysis support, not the novelty claim |
| Inocras makes it a vendor project. | Inocras is a sequencing service option. The original science is serial sample design, phenotype-genotype integration, and interpretation. | Vendor role table | Vendor output is raw/genomic data, not the scientific conclusion |
| Budget is too high. | Serial clinical sampling, CTC phenotyping, matched tissue-normal NGS, cfDNA/CTC feasibility, FFPE/mIHC, and data integration are the cost drivers. Stage gates prevent open-ended spend. | Budget table | If Samsung trims, reduce CTC-low-input NGS and GeoMx first, not serial blood |

## 7. Detailed 30억 / 3년 Budget

Unit: 억 KRW. Direct research budget basis; institutional overhead, if required by local rules, must be handled separately or absorbed by reducing lower-priority assays.

| Category | Year 1 | Year 2 | Year 3 | Total | Rationale |
|---|---:|---:|---:|---:|---|
| Personnel | 1.5 | 1.8 | 1.8 | 5.1 | PM, clinical CRC, CTC technologist, bioinformatician, statistician, part-time pathology/data manager |
| Clinical cohort / biobank ops | 0.7 | 0.7 | 0.7 | 2.1 | Consent, serial blood draw, sample processing, freezer/cold chain, follow-up abstraction |
| CTC enumeration and EMT phenotyping | 1.2 | 1.8 | 1.2 | 4.2 | Assay lock, antibodies/probes, imaging, replicate QC, benign/healthy controls |
| NGS / Inocras or equivalent vendor | 1.8 | 3.0 | 1.8 | 6.6 | Matched tissue-normal NGS, cfDNA/MRD pilot, CTC-enriched/low-input feasibility subset, raw data delivery |
| FFPE mIHC / GeoMx / protein validation | 1.4 | 1.5 | 0.7 | 3.6 | Tumor-lineage, EMT, immune/barrier, NIS/RAI markers; ROI validation in selected cases |
| Fresh tissue functional validation | 0.6 | 1.2 | 0.6 | 2.4 | Slice/organoid feasibility, short perturbation assays, viability/iodide or marker response only in selected fresh cases |
| Data/compute/statistics | 0.5 | 0.7 | 0.6 | 1.8 | Secure data environment, workflow, clone-state model, reproducible analysis |
| Consumables and QA | 0.4 | 0.4 | 0.4 | 1.2 | Tubes, Streck/EDTA, slides, extraction kits, library QC, controls |
| Equipment use / service contracts | 0.5 | 0.4 | 0.3 | 1.2 | Imaging, sequencer/library QC instruments, pathology scanner/analysis service |
| IP/regulatory/meetings/contingency | 0.4 | 0.6 | 0.8 | 1.8 | IRB amendments, data transfer review, patent map, external advisory, publication/meeting, contingency |
| **Total** | **9.0** | **12.1** | **8.9** | **30.0** | Non-flat budget because Year 2 is the sequencing/scale-up peak |

Budget cut hierarchy if needed:

1. Preserve serial blood, CTC phenotype, matched tissue-normal NGS.
2. Reduce GeoMx breadth before reducing cohort timepoints.
3. Make CTC low-input NGS a CTC-high subset if feasibility is limited.
4. Defer fresh tissue perturbation breadth to a small mechanistic subset.

## 8. IRB / Data Collection Checklist

### Study Population

- Adult patients with suspected or confirmed PTC undergoing thyroidectomy.
- Enrich high-risk features where possible: clinically suspicious LN, large tumor, ETE, tall-cell/hobnail/diffuse sclerosing features, BRAF/TERT/fusion if available, recurrent disease.
- Comparator groups: benign thyroid surgery and/or healthy blood controls for assay specificity.

### Timepoints

- T0: preoperative baseline, ideally within 1-2 weeks before surgery.
- T1: early postoperative, around 2 weeks.
- T2: 3 months.
- T3: 6-12 months or clinically triggered by Tg/anti-Tg/ultrasound concern.
- Optional: before/after RAI when clinically indicated.

### Specimens

- Peripheral blood for CTC: specify tube type, volume, processing time, storage condition.
- Peripheral blood for cfDNA: Streck or EDTA protocol; plasma separation SOP.
- Matched normal DNA: buffy coat or saliva/buccal, depending vendor requirements.
- FFPE tumor block/slides and H&E reference.
- Fresh tissue optional: only if surgical workflow supports rapid handoff.
- Pathology ROI annotation for tumor, invasive front, lymphovascular invasion, LN metastasis, stroma-rich areas.

### Clinical Variables

- Age, sex, prior thyroid disease, Hashimoto/autoimmune status.
- Surgery type, extent, LN dissection.
- Tumor size, multifocality, bilaterality.
- Histologic subtype, margin, ETE, LVI, perineural invasion.
- LN status: central/lateral, number examined, number positive, extranodal extension.
- ATA risk group, AJCC stage.
- BRAF V600E, TERT promoter, RET/NTRK/ALK fusions, RAS if already tested.
- Tg, anti-Tg, TSH, ultrasound findings, RAI dose/uptake if applicable.
- Recurrence/persistence events with absolute dates.

### Assay/QC Requirements

- Prespecified CTC positivity definition.
- Epithelial CTC: EpCAM/KRT positive, CD45 negative.
- Mesenchymal CTC: VIM/TWIST/ZEB1 or equivalent EMT marker positive, CD45 negative, thyroid-lineage support where possible.
- Hybrid CTC: epithelial plus EMT marker.
- Thyroid-lineage support: PAX8/TG/TTF1/NIS where feasible, acknowledging sensitivity limits.
- Duplicate imaging/readout for a subset.
- Spike-in recovery and benign/healthy controls.
- Chain of custody and sample processing timestamp.

### Consent / Governance

- Serial blood and tissue genomic testing.
- Somatic plus potential germline finding handling.
- Data transfer to external sequencing vendor.
- Return-of-results policy: research-only vs clinical report.
- Recontact for additional blood/follow-up.
- De-identification key management.
- Data retention and withdrawal policy.

## 9. Inocras Vendor Questionnaire

### Scope Fit

1. Can you support a research cohort of PTC thyroidectomy patients with matched tumor-normal sequencing?
2. Which product/service is appropriate: CancerVision, MRDVision, custom research WGS, targeted-enhanced WGS, or bioinformatics-only?
3. Can the project receive FASTQ, BAM/CRAM, VCF, CNV, SV, fusion, TMB/MSI/HRD/signature outputs?
4. Can reports be research-only and customized for BRAF, TERT promoter, RET/NTRK/ALK fusions, RAS, EIF1AX, PPM1D, CHEK2, CNV/SV?

### Sample Requirements

5. Minimum FFPE tumor input: curls/slides, tumor purity threshold, H&E requirement, macrodissection support.
6. Matched normal requirement: blood, buccal, or both.
7. cfDNA input requirement: tube type, blood volume, plasma volume, minimum cfDNA mass, processing window.
8. Can you process CTC-enriched cell pellets or low-input/single-cell CTC DNA/RNA?
9. If yes, what WGA or library method is used, and what are expected SNV/CNV/SV/fusion limitations?
10. If no, can you accept CTC-enriched DNA as exploratory research input with no clinical report?

### Performance / Validation

11. What is the validated limit of detection for cfDNA/MRDVision under low tumor fraction?
12. Is PTC or thyroid cancer represented in your validation data?
13. How does assay sensitivity change for low-shedding, low-stage tumors?
14. How are TERT promoter, fusions, and structural variants handled?
15. Can tumor-informed variants from FFPE be used to design or interpret serial cfDNA?
16. Can you distinguish germline, clonal hematopoiesis, and tumor-derived variants?

### Operations

17. Turnaround time for each sample type and batch size.
18. Cost per sample for tumor-normal WGS, cfDNA serial sample, CTC-enriched exploratory sample, and bioinformatics-only analysis.
19. Batch requirements and whether Year 2 scale-up can be scheduled.
20. Failed sample policy and re-run cost.
21. Data storage location, cross-border transfer, encryption, access control, and IRB language needed.
22. CAP/CLIA status for clinical services and boundary for research-use outputs.
23. Whether Korean hospital data can remain within Korea or must be transferred abroad.
24. Named scientific contact for assay design before IRB submission.

### Deliverables

25. Example de-identified report and raw data manifest.
26. Variant annotation database versions.
27. Genome build.
28. Pipeline versioning and reproducibility statement.
29. Molecular Tumor Board or expert interpretation availability.
30. Publication policy and vendor acknowledgment requirements.

## 10. Claim Boundaries And Kill Criteria

### Allowed Claims

- PTC CTCs and EMT-state CTCs are measurable in published cohorts.
- Local public omics show thyroid cancer state heterogeneity beyond driver mutation.
- Spatial transcriptomics supports tissue-state organization, not CTC shedding.
- Matched tissue NGS is needed to interpret postoperative blood signals.
- The project will test whether thyroidectomy reveals a CTC-EMT/genetic state transition.

### Forbidden Claims

- CTC-EMT state already predicts recurrence in our cohort.
- Public pilot proves CTC biology.
- cfDNA/CTC NGS will work in all early-stage PTC.
- This is a ready clinical diagnostic.
- Inocras or any vendor technology is the scientific novelty.

### Stage-Gated Kill Criteria

| Gate | Pass criterion | If fail |
|---|---|---|
| Assay feasibility | Reproducible CTC detection/subtyping in controls and early patient set | Narrow to CTC count + cfDNA, drop EMT subtyping breadth |
| NGS feasibility | Tissue-normal NGS success >85%; cfDNA usable in high-risk subset | Keep tissue map; treat cfDNA as exploratory |
| CTC low-input feasibility | Genotype signal in CTC-enriched subset above prespecified QC | Do not spend Year 2/3 money on single-cell CTC NGS |
| Clinical follow-up | Tg/US/pathology data captured for >80% of enrolled patients | Reduce outcome claims to biological transition only |
| Reviewer defense | One primary endpoint remains clean and testable | Cut spatial/drug appendices from main proposal |

## Sources Used

- Samsung Science program page: https://www.samsungstf.org/ssrfPr/program/basic.do
- Samsung ICT program page: https://samsungstf.org/ssrfPr/program/ict.do
- Samsung program overview: https://www.samsungstf.org/ssrfPr/intro/business.do
- Yu et al. PTC CTC thyroidectomy prospective cohort, PubMed/PMC record: https://pubmed.ncbi.nlm.nih.gov/38445526/
- Li et al. CTC EMT/CD133 thyroid cancer: https://www.spandidos-publications.com/10.3892/mco.2022.2574/abstract
- Sato et al. BRAF V600E ctDNA before/after PTC surgery: https://journals.sagepub.com/doi/abs/10.1089/thy.2021.0267
- Tarasova et al. plasma NGS thyroid cancer: https://journals.sagepub.com/doi/abs/10.1089/thy.2023.0204
- TCGA-THCA integrated genomic characterization: https://gdc.cancer.gov/about-data/publications/thca_2014
- GSE184362 PTC scRNA: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184362
- GSE250521 thyroid spatial transcriptomics: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250521
- Inocras CancerVision provider page: https://inocras.com/cancer-provider/
- Inocras researcher services: https://inocras.com/ko/researchers/
