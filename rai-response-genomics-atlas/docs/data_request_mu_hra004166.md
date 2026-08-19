# Data request — Mu et al. JCEM 2024 / NGDC HRA004166

## Target

Mu et al. **Characterizing Genetic Alterations Related to Radioiodine Avidity in Metastatic Thyroid Cancer.** J Clin Endocrinol Metab 2024 Apr 19;109(5):1231-1240. PMID 38060243 · PMC11031230. 220 patients enrolled, **214 analyzed** (6 missed RAI scan), with four RAI uptake patterns:

1. **I-RAIR** (initially RAI refractory) — **n = 80**.
2. **C-RAIA** (continually RAI avid) — **n = 48**.
3. **P-RAIR** (partly RAI refractory) — **n = 10**.
4. **G-RAIR** (gradually RAI refractory) — **n = 19**.

Data accession: **HRA004166** at NGDC (https://ngdc.cncb.ac.cn/), access restricted for patient confidentiality.

## Verified mutation patterns (PMC fetch 2026-05-21)

- **BRAF V600E**: enriched in I-RAIR (61.1% of mutated cases).
- **TERT promoter**: 50.7% of mutated I-RAIR cases.
- **RAS family**: more frequent in I-RAIA (4.5% I-RAIR vs higher in I-RAIA).
- **TP53**: associated with I-RAIR pattern.
- **Late-hit mutations (TERT/TP53/PIK3CA combined)**: 50.0% of I-RAIR vs 26.9% of I-RAIA.

This is **already enough for Figure 4's "molecular gray zone" framing** without per-patient access — the published frequencies show driver mutations alone do not partition the four uptake classes cleanly. P-RAIR (partly refractory) and G-RAIR (gradually refractory) are the gray zone.

## Why we need this data

This is the **single best argument** for the manuscript's "molecular gray zone" framing — Mu defines four ordered uptake patterns rather than a binary. Even at the mutation-frequency level (already in the paper's supplementary), it lets us show that driver mutations alone do not partition the four patterns cleanly; an expression-based panel is required.

If we can secure per-patient mutation data, we can:

- Re-fit a multinomial logistic regression with the panel proxy.
- Quantify the panel's added value beyond BRAF/RAS/TERT/TP53/PIK3CA + fusions.
- Replicate Mu 2024's mutation-frequency table with our 8-gene Tier-4 proxy overlaid.

## What to extract from public sources first (no DAC required)

1. **Paper supplementary tables** — JCEM typically publishes per-class mutation-frequency tables. These alone enable the driver-overlap argument.
2. **NGDC HRA004166 access page** — check public summary statistics, sample counts, sequencing coverage.
3. **GitHub / Zenodo** — Chinese groups sometimes post analysis code separately; check the paper's data availability statement.

## Data access conditions (NGDC HRA)

HRA accessions usually require:

- Institutional Data Access Committee (DAC) application via NGDC's GSA-Human portal.
- PI signature on a Data Use Agreement.
- IRB approval at the requester's institution.
- A research proposal (1–2 pages).
- Sometimes a fee or local collaborator co-PI.

Timeline expectation: 4–12 weeks from application to data access, with variable success rate for non-Chinese requesters.

## Email draft — corresponding author (parallel to DAC application)

Subject: Data request and possible collaboration: RAI uptake pattern cohort (Mu et al. JCEM 2024, HRA004166)

> Dear Dr. Mu (or current corresponding author),
>
> I am Seungho Cook (orcid id), a computational researcher in the Department of [department], Seoul National University Bundang Hospital. With Dr. Yu Hyeong-won we are validating a parsimonious 8-gene thyroid differentiation / iodide-handling panel (Cook, Yu et al. 2026, preprint at bioRxiv link) across multiple RAI-response cohorts.
>
> Your 2024 JCEM paper introducing the four RAI uptake patterns in 220 distant-metastatic DTC patients (HRA004166) is the most informative dataset for the "molecular gray zone" framing — uptake patterns clearly stratify beyond the binary refractory/avid distinction. We would like to:
>
> (a) Validate our 8-gene panel against your four uptake-pattern classes (or their published per-patient mutation summary, if direct expression / NGS data are restricted).
> (b) Where possible, fit a multinomial model with the panel score as an added covariate alongside BRAF / RAS / TERT promoter / TP53 / PIK3CA / fusion status.
>
> We have submitted (or plan to submit) a Data Access Committee request through NGDC for HRA004166. Beyond the formal channel, we would welcome any direct collaboration arrangement: co-authorship on the resulting paper, joint statistical analysis, or sharing of per-patient mutation tables only (if expression data are not available).
>
> Even if direct data sharing is not feasible, summary tables you would be willing to share (per-class panel-gene expression means, or per-class mutation frequencies broken down by BRAF/RAS/TERT/TP53/PIK3CA/fusion) would significantly strengthen our manuscript.
>
> Happy to send the preprint and analysis plan ahead of any commitment.
>
> Best regards,
>
> Seungho Cook (국승호) · Co-PI Dr. Yu Hyeong-won · Seoul National University Bundang Hospital
> [email] · ORCID [link]

## Tracking

- [ ] PubMed citation verified
- [ ] Supplementary tables downloaded → `data/external/mu_2024_supp/`
- [ ] NGDC HRA004166 page reviewed; access conditions noted
- [ ] DAC application submitted (date: ____)
- [ ] Email to corresponding author (date: ____)
- [ ] Response received (date: ____)
- [ ] Data access granted / declined / pending
