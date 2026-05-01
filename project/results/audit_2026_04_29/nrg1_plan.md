---
title: "NRG1 Korean germline × somatic — separate paper trajectory plan"
date: 2026-04-29
purpose: PROMPT I deliverable. Plan only — no germline data is currently in the project.
status: BLOCKED on data access. Recommend defer to dedicated trajectory after 분당 SNUH or KoGES germline data arrives.
---

# 1. Why NRG1?

**Memory `v17_dark_matter_pivot_2026_04_29.md` (PM portfolio review)** lists NRG1 germline×somatic as the **third separate trajectory** (after main 8-gene Dark Matter paper and Graves' parallel paper). Idea 3 in the rThyroid portfolio.

**Biological premise.** Korean GWAS (KARE / KoGES) reports NRG1 SNPs as thyroid cancer risk loci. NRG1 is a Schwann-cell / RAS-pathway-adjacent ligand. Hypothesis: germline NRG1 risk-allele carriers may have systematically different somatic mutation landscape (BRAF / RAS / RET-fusion frequency).

# 2. Current data gap

| Resource | In project? | Source needed |
|----------|-------------|---------------|
| Korean germline genotypes (NRG1 SNPs, ≥10K subjects) | NO | KoGES (gov.kr request), 분당 biobank (post-collaboration), CUKB (Catholic University Korea Biobank) |
| Korean somatic mutations matched to germline | NO | Tied to same biobank |
| 1000G EAS reference for NRG1 allele frequency | NO (downloadable) | Ensembl REST or 1000G FTP |
| TCGA-THCA Asian subgroup matched germline | partial | TCGA controlled-access (dbGaP) |

# 3. Recommended scope for separate paper

| Section | Content | Source |
|---------|---------|--------|
| Discovery | KoGES NRG1 lead SNP × somatic BRAF / RAS / TERT | KoGES + matched WGS/NGS |
| Replication | 분당 cohort or Yoo 262 if germline genotyped | post-collaboration |
| Functional | TCGA-THCA Asian subgroup interaction term | dbGaP request |
| Mechanism | NRG1 / ERBB axis literature review | PubMed |

Target venue (small paper, separate from main): **J Med Genet** (IF ~3) or **EJHG / Cancer Genet** (IF 3-5). Not in main 8-gene Dark Matter paper.

# 4. Pre-flight checklist before starting

- [ ] Confirm Korean GWAS lead SNP for NRG1 thyroid cancer risk (rs### — find in PMID + LDlink)
- [ ] Verify dbGaP access status for TCGA-THCA controlled germline data
- [ ] Outreach to KoGES (Korean CDC) for genotype access — 6-12 month review
- [ ] Outreach to 분당 SNUH for biobank access (already drafted, 2026-04-27)

# 5. What this paper does NOT need from the main 8-gene paper

It does not depend on the 8-gene panel results. It uses driver mutation status (BRAF / RAS / RET / TERT) only. So it can run in parallel without blocking the main Dark Matter paper.

# 6. Decision

**Defer. Plan only.** No germline data → no analysis possible today. Pre-flight checklist above. When germline data arrives via 분당 outreach or KoGES request, this trajectory becomes a 3-month side-paper executable in parallel with main paper revisions.

# 7. Files

- This plan: `project/results/audit_2026_04_29/nrg1_plan.md`
- Pre-existing related: `project/notebooks_or_scripts/v17_external_K601_TERT.py` (TERT external work, not NRG1)
- Outreach email (Bundang): `project/outreach/email_bundang_KR_v2.md` (2026-04-27)
