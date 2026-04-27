# Landa collaboration email — v1 polished draft

**To:** landai@mskcc.org (Iñigo Landa, MSKCC — Fagin lab)
**Subject:** GSE76039 + your PDTC/ATC cohort — strongest external validation of a TCGA-THCA transcriptomic axis

---

Dear Dr. Landa,

I am Seungho Cook, a computational researcher in Seoul. I am writing because your published PDTC/ATC cohort (Landa et al., *JCI* 2016) is the strongest external validator of a transcriptomic axis we identified in TCGA-THCA, and I wanted to share the result with you directly before submission.

**The result.** In 513 TCGA-THCA primary tumours, unsupervised clustering yields a binary axis (DM1 / DM2) orthogonal to BRAF/RAS status. The axis maps onto a clear dedifferentiation gradient PTC → PDTC → ATC. When tested on GSE76039 (n=37, the public expression subset of your cohort), it returns a direction-correct AUC of **0.974** — by far the cleanest external signal we have. Adding TERT promoter status (36 TCGA-THCA carriers) yields 4-group OS stratification (BRAF / RAS / TERT+ / Triple-neg) with logrank p = 4.92e-06.

This means your data is not background validation — it is the reason the dedifferentiation framing is defensible at all. We would like to acknowledge that explicitly in the paper, and to ask two things.

**Two asks.**

1. *Clinical metadata beyond cBioPortal.* The cBioPortal mirror of MSK-PDTC-ATC (n=81) gives genomic calls but not RAI status, time from PTC primary to PDTC/ATC progression, or treatment trajectory. If any of these are available alongside the GSE76039 expression data, even a partial join would strengthen the trajectory analysis substantially.

2. *Manuscript feedback.* I would be honoured to send you the current draft (npj Precision Oncology target) for review. Given that your PDTC/ATC framing is foundational to our dedifferentiation argument, I would welcome any form of involvement that fits your interest — from a comment on the framing through to co-authorship if the contribution warrants it.

A short technical summary follows. I am happy to share the full manuscript and supplementary materials whenever convenient.

Thank you for the foundational work — without GSE76039 we would not have been able to test the dedifferentiation hypothesis externally at all.

Best regards,

Seungho Cook
kukshomr@gmail.com
[Affiliation line — pending]

---

## Attachment — one-page summary

**Discovery cohort.** TCGA-THCA n=513, RNA-seq (STAR Counts), unsupervised consensus clustering on a 67-gene THCA-relevant feature set → DM1 (n=109) / DM2 (n=69) within the dark-matter-positive subset.

**External validation cohorts.**

| Cohort | n | Source | Result |
|---|---:|---|---|
| GSE76039 | 37 | Landa 2016 PDTC/ATC subset | **AUC 0.974** (direction-correct dedifferentiation projection) |
| GSE126698, GSE213647 | combined ~120 | technical replication | AUC 0.81–0.89 |
| GSE27155 | 79 | classical PTC array | direction-correct |

**4-group integration with TERT.** Recovered 36 TCGA-THCA TERT promoter mutations from cBioPortal `thca_tcga_pub` (Sanger-validated 2014 *Cell* paper). 4-group OS logrank p = 4.92e-06; TERT+ events 6/36 (16.7%) vs WT 10/468 (2.1%).

**8-gene RAI-response panel.** 5-fold CV AUC 0.954, ΔAUC vs BRAF V600E status alone = +0.132.

**Limitations declared.** (a) GSE76039 n=37 is small; we present bootstrap CI bands. (b) PDTC/ATC trajectory is currently inferred from histology, not longitudinal — your clinical metadata would directly address this. (c) TERT calls are Sanger-mirror, not re-called.

---

## Internal notes (not part of email)

- Landa is a working scientist, not a senior emeritus — tone should be peer-to-peer, not deferential.
- "Direction-correct AUC 0.974" is the hook — concrete and unambiguous. He'll recognize whether the number is plausible immediately.
- Defer co-authorship ask: "any form of involvement" is the soft form. Don't propose authorship in first contact.
- If he replies, the natural follow-up is a data-sharing discussion + 30-min Zoom. Do not propose Zoom in first contact.
- Reference: Landa et al. *J Clin Invest* 2016, PMID 27613709 ("Genomic and transcriptomic hallmarks of poorly differentiated and anaplastic thyroid cancers"). Verify GSE76039 = the expression subset before sending.
