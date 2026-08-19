---
title: Submission readiness assessment
audit_date: 2026-07-30
auditor: manuscript-audit-agent
---

# 09 — Submission readiness

## One-sentence verdict

The manuscript is NOT submission-ready: three fatal scientific issues (HR direction unverified, figure architecture unresolved, abstract overclaims), six unfilled author-voice sections, two missing figures (3C and 3D), two TODO placeholder citations, and terminology violations throughout require resolution before the package can be submitted to Nature Communications.

---

## Dimension scores (0 = not ready, 1 = repairable, 2 = ready)

| Dimension | Score | Rationale |
|---|---|---|
| central_claim_validity | 1 | Core biology (DM1 axis, methylation, external validation) is scientifically solid and well-supported. Fatal issues are in how the clinical claims are framed (title, abstract, IHC "immediate deployment"). Repairable with targeted language edits. |
| statistical_validity | 1 | Survival analysis direction (reference group) must be verified from source code before claims about which state is adverse can be made. All other statistical methods are appropriate. Single-cohort interaction with 33 events is marginal but disclosed. Repairable once direction is confirmed. |
| clinical_restraint | 0 | Abstract uses "immediate translational deployment," "candidate treatment-selection framework," and "clinical-grade IHC 3-plex derivative" without adequate qualification. The IHC result is a methylation proxy. The title uses "predicts." Limitations are entirely absent. Multiple overclaims present. Requires substantive revision. |
| external_validation | 1 | Multi-cohort validation is genuinely strong (14-19 entries, 80/80 direction, Lee p=2.7e-7, Mun 7/7, single-cell intrinsic). Issues: denominator definitions ("80 cells," "19 entries" vs "14 entries") need precise definition; GPL570 fourth cohort accession missing; Lee n discrepancy (370 vs 632) needs resolution. Repairable. |
| figures | 0 | Fig 3C is missing entirely. Fig 3D needs rebuild. Fig 5D (ATA mosaic) needs build. Fig 7 and Fig 8 are in a web directory rather than the manuscript figure directory, and have no matching captions file. Three competing figure architectures exist without resolution. "Extended Data" terminology must be corrected to "Supplementary Figure" throughout. Not ready. |
| methods | 1 | Cohort accessions are mostly present. Key gaps: KMeans seed unstated, HM450 probe aggregation unstated, fourth GPL570 cohort missing, GSE126698 missing from v2 cohort list, ATA risk-tier operationalization absent. Two TODO placeholders for Zenodo/repo URL remain. Repairable with targeted additions. |
| reproducibility | 0 | Random seed unspecified for clustering. No figure-to-source-data manifest. Code repository not established. Zenodo DOI not obtained. Two TODO placeholders in submitted Methods draft. These are direct NC policy requirements that cannot be waived. Not ready. |
| citations | 1 | Landa 2016 correctly cited in v2 (Krishnamoorthy misattribution corrected). Liu 2018 PFI cited correctly. Chakravarty 2011 and Rothenberg 2015 cited but full references need author confirmation. Mun 2025 has no accession or DOI. Lee n discrepancy needs resolution. Repairable with author input. |
| journal_formatting | 0 | "Extended Data" used throughout (must be "Supplementary Figure" for NC). Life Sciences Reporting Summary absent. Statistics and Reproducibility subsection needs formatting. Data Availability needs accession numbers inline. Not ready. |
| administrative_completeness | 0 | Six voice-protected author sections unfilled (Introduction hook, Introduction aims, Discussion §3.1, Limitations, Cover letter, Reviewer Q9). Author Contributions incomplete. Acknowledgements/Funding absent. Institutional corresponding email unverified. Competing interests needs verification. Not ready. |

---

## Total: 5/20 (25% ready)

---

## Priority action queue (estimated effort)

### Tier 0 — Blocking everything (must be done before any other step)

| Action | Owner | Estimated effort |
|---|---|---|
| Verify OS and PFI Cox model reference direction from source code | Author / analyst | 30 min |
| Make editorial decision on figure architecture (5 vs 6 vs 8 figures) | Author | 15 min |
| Correct title ("predicts" → "marks" or equivalent) | Author decision | 5 min |

### Tier 1 — Author voice (uniquely requires author; cannot be delegated)

| Action | Owner | Estimated effort |
|---|---|---|
| Draft Introduction hook (V1) | Author | 1-2 hr |
| Draft Introduction aims paragraph (V2) | Author | 30 min |
| Draft Discussion §3.1 mechanism interpretation (V3) | Author | 2-3 hr |
| Draft Limitations paragraph (V4) | Author | 1-2 hr |
| Draft Cover letter paragraph 1 (V5) | Author | 30 min |
| Draft Reviewer Q9 (V6) | Author | 30 min |

### Tier 2 — Computational (can be done in parallel with Tier 1)

| Action | Owner | Estimated effort |
|---|---|---|
| Build Fig 3C (beta vs expression scatter, 4-gene grid) | Analyst | 1-2 hr |
| Rebuild/confirm Fig 3D (per-driver-class beta bar) | Analyst | 30 min |
| Build Fig 5D (ATA mosaic) | Analyst | 30 min |
| Resolve Lee n=370 vs n=632 and label both analyses | Analyst | 1 hr |
| Identify fourth GPL570 cohort and add accession | Analyst | 30 min |
| Add GSE126698 to cohort list | Analyst | 15 min |
| Obtain and add Mun 2025 accession | Author/Analyst | 30 min |
| Specify KMeans seed and add to Methods | Analyst | 10 min |
| Specify HM450 probe aggregation method | Analyst | 15 min |
| Define ATA risk-tier operationalization rules | Analyst | 30 min |

### Tier 3 — Editorial cleanup (can be done at any time)

| Action | Owner | Estimated effort |
|---|---|---|
| Replace all "Extended Data" / "ED" with "Supplementary Figure" | Analyst | 30 min |
| Rewrite abstract (design disclosure, IHC proxy, treatment-selection qualifier) | Analyst with author review | 1 hr |
| Correct MAR language throughout | Analyst | 15 min |
| Correct "silences" to "epigenetic correlate" in Fig 3 title | Analyst | 5 min |
| Replace GSE151179 approximate values with exact values | Analyst | 15 min |
| Copy/regenerate Fig 7 and 8 to NC submission directory at 600 dpi | Analyst | 1 hr |
| Reconcile v2 draft with chosen figure architecture | Analyst | 2-3 hr |

### Tier 4 — Administration (can be done in parallel)

| Action | Owner | Estimated effort |
|---|---|---|
| Establish public GitHub/code repository | Analyst | 2-4 hr |
| Obtain Zenodo DOI | Analyst | 1 hr |
| Confirm institutional corresponding email | Author | 5 min |
| Complete Author Contributions | Author | 30 min |
| Complete Acknowledgements and Funding | Author | 30 min |
| Verify Competing Interests | Author | 15 min |
| Download and complete Life Sciences Reporting Summary | Analyst | 2-4 hr |
| Create figure-to-source-data manifest | Analyst | 2-3 hr |

---

## Earliest possible submission estimate

Assuming all Tier 0 decisions are made today (2026-07-30), Tier 1 voice sections are drafted within 1 week, and Tier 2-4 tasks are completed in parallel:

**Optimistic: 2-3 weeks from today (by ~2026-08-20)**

**Conservative: 4-5 weeks from today (by ~2026-09-01)**

The rate-limiting steps are:
1. Author voice sections (6 sections; cannot be parallelized by analyst)
2. Fig 3C build (depends on confirmed HR direction from Tier 0)
3. Life Sciences Reporting Summary (time-consuming form)
4. Code repository and Zenodo DOI (requires clean code organization)
