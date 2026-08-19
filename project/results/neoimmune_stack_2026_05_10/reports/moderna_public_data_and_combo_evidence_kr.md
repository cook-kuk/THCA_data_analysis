# Moderna/V940 Public Data and Combination-Therapy Evidence KR

## Bottom line
모더나/Merck V940(mRNA-4157) 관련 공개 데이터는 있습니다. 다만 **환자별 neoantigen sequence/candidate raw table, 모델 feature table, proprietary selection score**는 공개되어 있지 않은 것으로 봐야 합니다.

공개된 것은 크게 네 종류입니다.

1. **임상 aggregate efficacy/safety data**: KEYNOTE-942 phase 2b, V940 + pembrolizumab vs pembrolizumab.
2. **면역반응 실험 데이터/논문**: KEYNOTE-603 phase 1, T-cell response characterization.
3. **병용치료 임상시험 registry/protocol data**: INTerpath melanoma/NSCLC/cSCC/RCC/urothelial programs.
4. **회사 공시/press release**: up to 34 neoantigens, NGS + machine-learning/proprietary algorithm, single mRNA, LNP, automated workflow.

## 1. Clinical efficacy/safety public data

| Trial | Public source | What is public | Useful for us | What is not public |
|---|---|---|---|---|
| KEYNOTE-942 / mRNA-4157-P201 / NCT03897881 | Lancet 2024 | randomized phase 2b design, n=157, V940+pembro vs pembro, RFS/DMFS, safety aggregate | endpoint framing, adjuvant/post-resection logic, PD-1 combination rationale | per-patient mutation/neoantigen list, individual vaccine sequence, selection score |
| INTerpath-001 / NCT05933577 | ClinicalTrials.gov/NCI/Merck | phase 3 melanoma design, enrollment target, schedule, comparator | phase 3 endpoint design and top-N adjuvant framing | results not posted yet; no raw neoantigen table |

### Key KEYNOTE-942 numbers from public Lancet summary
- 157 randomized patients.
- 107 received combination assignment, 50 pembrolizumab monotherapy assignment.
- 18-month recurrence-free survival: 79% vs 62%.
- HR for recurrence/death: 0.561 with two-sided p=0.053.
- Grade >=3 treatment-related adverse events: 25% combination vs 18% monotherapy.
- No V940-related grade 4-5 events reported in the summary.

Use in our paper: **adjuvant recurrence-risk framing**, not as our performance benchmark.

## 2. Experimental/immune-response public data

| Study | Public source | What is public | Useful for us | What is not public |
|---|---|---|---|---|
| KEYNOTE-603 / NCT03313778 | Cancer Discovery 2024 / PubMed abstract | phase 1 T-cell response characterization, NSCLC n=4 monotherapy, melanoma n=12 with pembrolizumab, safety and immunogenicity summary | mechanistic proof that individualized neoantigen mRNA can induce/strengthen T-cell responses | full patient-level neoantigen sequences and raw assay matrices may be restricted or in journal supplements only |

### Key public immune-response points
- V940 targets up to 34 patient-specific tumor neoantigens.
- Phase 1 examined resected NSCLC monotherapy and resected melanoma plus pembrolizumab.
- No grade 4/5 adverse events or dose-limiting toxicities in the PubMed abstract.
- mRNA-4157 alone induced de novo and strengthened pre-existing T-cell responses.
- Combination therapy showed sustained neoantigen-specific T-cell responses and expansion of cytotoxic CD8/CD4 T cells in public abstract summaries.

Use in our paper: **why immunogenicity/T-cell response matters beyond binding**.

## 3. Combination-therapy public data

| Program | Disease | Phase/status public | Combination | Public use |
|---|---|---|---|---|
| KEYNOTE-942 | resected high-risk melanoma | phase 2b published | V940 + pembrolizumab vs pembrolizumab | proof-of-concept clinical frame |
| INTerpath-001 / V940-001 | high-risk stage II-IV melanoma | phase 3, active/not recruiting, estimated n=1089 | V940 + pembrolizumab vs placebo + pembrolizumab | registration-grade endpoint design |
| INTerpath-002 / NCT06077760 | NSCLC | phase 3 mentioned in company program | V940 + pembrolizumab | extension to lung cancer |
| INTerpath-009 | resectable NSCLC post neoadjuvant therapy | phase 3 initiated Oct 2024 | V940 + pembrolizumab after neoadjuvant pembro+chemo | minimal-residual-risk setting |
| INTerpath-007 / NCT06295809 | cutaneous squamous cell carcinoma | phase 2/3 | V940 + pembrolizumab, neoadjuvant/adjuvant | event-free survival framing |
| INTerpath-004 / NCT06307431 | RCC | phase 2 | V940 + pembrolizumab | adjuvant high-risk setting |
| INTerpath-005 / NCT06305767 | muscle-invasive urothelial carcinoma | phase 2 | V940 + pembrolizumab | adjuvant post-resection setting |

Use in our paper/business: the winning clinical pattern is **minimal residual/adjuvant disease + PD-1 + personalized top-N neoantigen set**.

## 4. Company public technical disclosure

From Moderna SEC filings and Merck releases, the public technical pattern is:
- NGS identifies patient-specific HLA type and tumor neoantigens.
- Algorithm/proprietary or machine-learning based selection designs an mRNA.
- Up to 34 neoantigens are encoded.
- Class I/CD8 and class II/CD4 responses are targeted/predicted.
- Neoantigens are encoded in a single mRNA sequence.
- Formulated in proprietary LNPs.
- Intramuscular administration.
- Automated manufacturing workflow for rapid turnaround.

What we can legally mirror:
- top34 cap,
- patient-specific HLA/mutation/RNA input contract,
- class I + class II awareness,
- local immunogenicity branch,
- top-N wet-lab queue,
- adjuvant/minimal residual risk framing.

What we cannot claim/copy:
- Moderna's proprietary algorithm,
- Moderna's LNP,
- Moderna's mRNA construct,
- Moderna's manufacturing workflow,
- Moderna-equivalent efficacy.

## How to use this immediately in NeoImmune-Stack

1. Add `top34` next to top20.
2. Add class II candidate branch later, because V940 is described as targeting both CD8 and CD4 responses.
3. Add `treatment_context` and `ICI_context` fields.
4. Add MRD/ctDNA optional field because adjuvant recurrence-risk trials increasingly use MRD biology.
5. Add assay-confirmed immune response endpoints:
   - ELISpot/ICS/tetramer/TCR expansion if available.
   - MS presentation if presentation is claimed.
6. Benchmark public patient-level data first; use hospital data for prospective validation later.

## Direct answer
**있다.** 공개된 모더나 데이터는 “raw patient neoantigen table”이 아니라, 임상 결과/안전성/면역반응/병용치료 설계 데이터입니다. 우리에게 가장 유용한 것은 patient-level raw peptide가 아니라 **top34 cap + adjuvant/ICI setting + T-cell response endpoint + manufacturing-testable candidate selection logic**입니다.

## Sources
- KEYNOTE-942 Lancet: https://www.sciencedirect.com/science/article/pii/S0140673623022687
- KEYNOTE-603 phase 1 immune-response PubMed: https://pubmed.ncbi.nlm.nih.gov/39115419/
- ClinicalTrials INTerpath-001 / NCT05933577: https://clinicaltrials.gov/study/NCT05933577
- Merck V940 PRIME description: https://www.merck.com/news/mrna-4157-v940-an-investigational-personalized-mrna-cancer-vaccine-in-combination-with-keytruda-pembrolizumab-receives-prime-scheme-designation-from-the-european-medicines-agency-for-adjuva/
- Merck phase 3 melanoma program: https://www.merck.com/news/merck-and-moderna-initiate-phase-3-study-evaluating-v940-mrna-4157-in-combination-with-keytruda-pembrolizumab-for-adjuvant-treatment-of-patients-with-resected-high-riskstage-iib-iv-melanom/
- Merck/Moderna INTerpath expansion: https://www.merck.com/news/merck-and-moderna-initiate-phase-3-trial-evaluating-adjuvant-v940-mrna-4157-in-combination-with-keytruda-pembrolizumab-after-neoadjuvant-keytruda-and-chemotherapy-in-patients-with-certain-ty/
- Moderna 2024 SEC 10-K public technical description: https://www.sec.gov/Archives/edgar/data/1682852/000168285225000022/mrna-20241231.htm
