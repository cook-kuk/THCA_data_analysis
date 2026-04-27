# Xing collaboration email — v1 polished draft

**To:** mxing1@jhmi.edu (Mingzhao Xing, MD, Johns Hopkins)
**Subject:** Transcriptomic axis (DM1/DM2) extending the BRAF/TERT genotype framework — TCGA-THCA reanalysis, npj submission

---

Dear Professor Xing,

I am Seungho Cook, an independent computational researcher in Seoul. I am writing to share a TCGA-THCA reanalysis that builds directly on your BRAF/TERT genotype framework (Xing et al., *JCO* 2014; Liu R. et al., *JAMA Oncol* 2017), and to ask whether your group has expression data we could use to externally validate it.

**Headline finding.** In 513 TCGA-THCA primary tumours, an unsupervised transcriptomic clustering yields a binary axis (DM1 / DM2) that is statistically orthogonal to BRAF/RAS mutation status — only 17 of 335 mutation-carrying tumours violate the canonical mapping. Anchoring this axis to an 8-gene RAI-response panel produces 5-fold CV AUC 0.954 vs 0.822 for BRAF V600E status alone (ΔAUC = +0.132).

**TERT integration.** We recovered 36 TERT promoter mutations for TCGA-THCA from the cBioPortal `thca_tcga_pub` study (Sanger-validated by the original 2014 *Cell* publication MAF) at canonical positions chr5:1295228 (C228T) and 1295250 (C250T). Adding these as a fourth genotype yields a 4-group survival stratification (BRAF / RAS / TERT+ / Triple-neg) with logrank p = 4.92e-06 — a direct transcriptomic-era extension of your 4-genotype framework in *JAMA Oncol* 2017.

**Two specific asks.**

1. *External validation.* Our DM1/DM2 axis would benefit substantially from cross-validation in your 1,051-patient JHU cohort (Liu R. 2017) or the 507-patient cohort from *JCO* 2014. Even a small subset with available expression data (microarray or RNA-seq) would let us report a true external validation of the 8-gene panel.

2. *Manuscript review.* I would be honoured to send you the current manuscript draft (npj Precision Oncology target) for any feedback you are willing to provide. Your input on the framing — particularly how the transcriptomic axis is positioned relative to your published genotype framework — would be invaluable, and I would welcome any form of involvement that fits your interest, from acknowledgement to co-authorship.

I have attached a one-page summary of the methodology and key results below. I am of course happy to share the full manuscript and supplementary materials at your convenience.

Thank you for considering this — your foundational work is the reason we believe a transcriptomic 4-group framework is worth proposing in the first place.

Best regards,

Seungho Cook
kukshomr@gmail.com
[Affiliation line: independent / formal affiliation pending]

---

## Attachment — one-page summary (paste below the signature)

**Cohorts.** TCGA-THCA n=513 (discovery); GSE76039 n=37 PDTC/ATC (Landa et al. *JCI* 2016, dedifferentiation gradient validation, AUC 0.974); GSE126698, GSE213647, GSE27155 (technical replication).

**8-gene RAI panel.** Selected by Lasso-stable signature on the DM1/DM2 axis. Genes: [list 8 — fill in from sample_master]. 5-fold CV AUC 0.954.

**Orthogonality.** BRAF×DM cross-tab shows DM1 = 79% BRAF, DM2 = 31% BRAF; conditional independence test rejects nested-axis hypothesis (p < 1e-4).

**TERT integration.** 36 / 513 (7.0%) TCGA-THCA samples carry promoter mutations. 73% are stage III/IV. TDS16 differentiation score: TERT+ mean 6.24 vs WT 6.88 (Δ = -0.64). 4-group OS logrank p = 4.92e-06.

**Limitations** (declared up-front). (a) TERT calls are Sanger-validated mirror data, not re-called from BAMs. (b) 6 events in TERT+ → wide HR CI; we present Firth-corrected estimates. (c) Single-cohort transcriptomic discovery — external validation in your cohort is the gap this email aims to close.

---

## Internal notes (not part of email)

- Xing's tone preference: concise, mutation-fluent, no hedging on numbers. He will ignore vague openers.
- The "co-author collaboration" ask is softened to "any form of involvement that fits your interest" — gives him room to decline cleanly.
- Lead with the orthogonality finding because that's what extends his framework; the TERT integration confirms it without competing.
- If he replies positively, the natural follow-up is a 30-min Zoom; do not propose this in first contact.
- Liu R. 2017 paper is *JAMA Oncology*, not JCO — the user's draft mixed this up. Verify before sending: PMID 28199474.
- Xing M. 2014 *JCO* 4-genotype paper is correct: PMID 24463735.
