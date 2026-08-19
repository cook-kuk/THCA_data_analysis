# Moderna V940 Public Reverse Engineering KR

## One-line take
Moderna/Merck의 V940(mRNA-4157)에서 공개적으로 배울 수 있는 핵심은 **“환자별 tumor mutanome + HLA + RNA 정보를 이용해 제한된 수의 neoantigen을 고르고, 단일 mRNA/LNP 제품으로 제조한 뒤 PD-1 blockade와 병용해 adjuvant recurrence endpoint를 본다”**입니다.

우리는 proprietary mRNA 제조/selection algorithm을 복제하지 않습니다. 공개 원리를 CLEAN-Neo++/NeoImmune-Stack의 **candidate selection operating layer**로 번역합니다.

## What Moderna Publicly Does

| Layer | Publicly visible Moderna/V940 pattern | What we copy legally | What we do not copy |
|---|---|---|---|
| Disease setting | resected high-risk melanoma adjuvant setting; expanding to NSCLC, cSCC, RCC, urothelial | surgical/adjuvant residual-risk framing | proprietary trial operations |
| Input | tumor DNA sequence, HLA type, and public reports mention WES/RNA-seq/HLA profiling | tumor-normal WES/targeted sequencing, RNA expression, HLA typing as required fields | private sample logistics/manufacturing details |
| Selection | algorithmically derived patient-specific neoantigens | multi-branch ranking with leakage audit | internal Moderna algorithm |
| Payload | single synthetic mRNA encoding up to 34 neoantigens in LNP | “top-N capped list” logic, e.g. top 20/34 ranked candidates | mRNA construct design/manufacturing |
| Mechanism | endogenous translation, natural antigen processing/presentation, T-cell response | rank candidates by presentation + immunogenicity + patient context | claim that ranking equals in vivo efficacy |
| Combination | pembrolizumab/anti-PD-1 | position output as ICI-compatible prioritization | drug-combination clinical claims |
| Endpoint | RFS/DMFS in trials; immune response in phase 1 | retrospective patient-level hit rate@20/@34; wet-lab validation bridge | clinical efficacy claim |

## Key public evidence
- KEYNOTE-942 phase 2b tested V940 plus pembrolizumab versus pembrolizumab alone in resected high-risk melanoma; 157 patients were randomized 2:1, and the combination showed longer recurrence-free survival, with 18-month RFS 79% versus 62% in the Lancet report.
- The Lancet article states mRNA-4157 encodes up to 34 neoantigens in an LNP formulation and is tailored to the patient's tumor mutanome and HLA type.
- Merck/Moderna describe V940 as a single synthetic mRNA coding for up to 34 neoantigens selected from the patient's tumor DNA mutational signature.
- Phase 1 mechanistic work reports de novo and strengthened pre-existing T-cell responses to targeted neoantigens, but this remains treatment-specific immune-response evidence, not a generic proof for our candidates.
- Phase 3 INTerpath-001 is active/not recruiting in high-risk melanoma with estimated enrollment 1,089 and completion around 2030; Moderna/Merck have expanded INTerpath programs into NSCLC, cSCC, RCC, and urothelial cancer.

## CLEAN-Neo++ translation

### Moderna-like top-N cap
Use `top20` for wet-lab feasibility and `top34` as the Moderna-analog cap. A strong report should show both:
- patient_hit_rate@20 / Recall@20
- patient_hit_rate@34 / Recall@34
- local-rescued public-missed positives inside top34
- binding-only false-positive traps inside top34

### Moderna-like input contract
Required:
- patient_id
- tumor-normal mutation calls
- mutant peptide and WT peptide
- HLA class I typing
- RNA expression
- VAF/clonality
- HLA LOH / B2M / APM context
- assay labels when available

Optional but high-value:
- MS immunopeptidomics
- T-cell assay
- ctDNA/MRD
- treatment context, especially PD-1/ICI exposure

### Moderna-like scoring architecture
1. Generate candidate universe from mutation calls.
2. Add frozen public presentation predictors.
3. Add local immunogenicity branches:
   - Wave8_TCR_SelfSim_full
   - Structure_LR
   - ESM2_Bayesian
   - quantum_kernel_no_anchor_gamma1.0 / W7A_QK_only
4. Add patient context gates:
   - expression
   - clonality/VAF
   - HLA LOH/APM/B2M
   - tumor type and treatment setting
5. Run leakage audit.
6. Emit top20/top34 candidates with evidence and rejection reasons.

## What makes this “대박” if executed correctly
The Moderna insight is not “predict binding.” It is the full operating chain:

patient tumor data → ranked patient-specific neoantigens → capped manufacturable/testable set → ICI-compatible immune activation hypothesis → recurrence endpoint.

NeoImmune-Stack can own the upstream selection layer:

patient tumor data → leakage-aware top-N candidate queue → wet-lab-ready rationale → MS/T-cell validation bridge.

## Reviewer-safe claim
“Inspired by the public design logic of individualized neoantigen therapy programs such as V940, CLEAN-Neo++ builds a leakage-aware, patient-level candidate selection layer that ranks patient-specific neoantigens for downstream validation. It does not claim vaccine manufacturing, clinical efficacy, or Moderna-equivalent proprietary selection.”

## Sources
- Lancet KEYNOTE-942 summary: https://www.sciencedirect.com/science/article/pii/S0140673623022687
- Merck V940 description / PRIME designation: https://www.merck.com/news/mrna-4157-v940-an-investigational-personalized-mrna-cancer-vaccine-in-combination-with-keytruda-pembrolizumab-receives-prime-scheme-designation-from-the-european-medicines-agency-for-adjuva/
- Merck phase 3 melanoma initiation: https://www.merck.com/news/merck-and-moderna-initiate-phase-3-study-evaluating-v940-mrna-4157-in-combination-with-keytruda-pembrolizumab-for-adjuvant-treatment-of-patients-with-resected-high-riskstage-iib-iv-melanom/
- ClinicalTrials INTerpath-001: https://clinicaltrials.gov/study/NCT05933577
- Merck/Moderna INTerpath expansion in NSCLC and other tumors: https://www.merck.com/news/merck-and-moderna-initiate-phase-3-trial-evaluating-adjuvant-v940-mrna-4157-in-combination-with-keytruda-pembrolizumab-after-neoadjuvant-keytruda-and-chemotherapy-in-patients-with-certain-ty/
- Phase 1 immune-response paper PubMed: https://pubmed.ncbi.nlm.nih.gov/39115419/
