[Date: 2026-04-27]

Editor-in-Chief
_Bioinformatics (Oxford University Press)_

Dear Editor,

We submit our manuscript "**Recovered external validation and proper GSEA define DM1/DM2 as an actionable thyroid cancer transcriptional axis**" for consideration.

We present (i) the **DIAL framework** (Direction Identifiability After Leave-out), a leak-detection safeguard for batch-correction protocols, and (ii) a multi-layer 5-cohort validation pipeline that recovers a robust transcriptomic axis (DM1/DM2) for thyroid cancer sub-classification. The DIAL framework caught its own data-leak artefact in our prior v5.1 protocol — a methodological self-audit we transparently report as v5.2. We also present a gene-ID recovery pipeline that improved cross-platform transfer in our 5-cohort analysis, and a proper preranked GSEA workflow with explicit pathway-name normalisation that achieves 100 % sign-agreement with our prior proxy GSEA on the overlapping pathway subset.

The manuscript matches _Bioinformatics_'s scope on **methodology and reproducibility**. The pipeline scripts (`notebooks_or_scripts/v17p35_*.py`) implement deterministic seeds (random_state = 42), version-pinned dependencies (sklearn 1.8.0, gseapy 1.13, scanpy 1.12.1), and explicit per-fold ComBat-seq batch-correction protocol that prevents leakage. The DIAL framework is generalisable beyond thyroid cancer to any cohort-pooling scenario in genomics.

The clinical relevance — an 8-gene RAI panel that outperforms BRAF V600E status by ΔAUC = +0.132 — is presented as a downstream application of the methodological framework rather than as the headline; this matches _Bioinformatics_' methods-paper voice better than a clinical-translation venue.

**Suggested reviewers** (no conflicts of interest declared):
- **Aedin Culhane, PhD** — University of Limerick — genomic data integration and reproducibility.
- **Lior Pachter, PhD** — Caltech — single-cell and bulk transcriptomic methods.
- **Casey Greene, PhD** — University of Colorado — biomedical data science and ML reproducibility.

All raw outputs, pipeline scripts, and reproducibility archives are available in the submission package. The authors declare no competing interests.

Sincerely,

Seungho Cook  
kukshomr@gmail.com

[PI signature TBD]
