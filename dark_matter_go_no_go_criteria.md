# Go / no-go criteria for the BRAF/RAS-negative thyroid cancer taxonomy paper

Working title: **Molecular taxonomy of BRAF/RAS-negative thyroid cancer reveals noncanonical driver classes and radioiodine-related differentiation states**

## Go decision

Proceed if the paper is framed as:

- molecular taxonomy,
- noncanonical driver-class enrichment,
- differentiation and RAI-gene expression states,
- external validation of expression signatures,
- exploratory clinical outcome only.

Do not proceed if the intended frame requires:

- survival prediction,
- proven recurrence stratification,
- proven RAI resistance,
- therapeutic target discovery,
- fusion-driven subtype discovery from current TCGA counts alone.

## Minimum go criteria

| Criterion | Required threshold | Current status |
|---|---|---|
| BRAF/RAS-negative prevalence | At least 5% of TCGA mutation-callable cohort | Pass: 137/482 = 28.42% |
| Subclass sample size | At least one interpretable subclass with n >= 5 | Pass: DICER1/EIF1AX/PPM1D n=7; TERT-only n=5; true driver-negative n=125 |
| Noncanonical driver enrichment | Directionally coherent and statistically supported | Pass: DICER1/EIF1AX/PPM1D 1/81 DM1 vs 6/55 DM2, Fisher p=0.0175 |
| Differentiation/RAI-gene expression | Clear effect between subclasses | Pass: RAI score p=4.57e-16; TDS16 p=3.96e-16 |
| Single-cell cellular-state support | Replication in at least one public scRNA dataset | Pass: GSE184362 pooled r=0.889; local GSE241184 r=0.905 |
| Survival limitation acknowledged | Survival not used as headline | Pass required; local PFI events only 9 |
| RAI-response validation path | Direct public dataset available | Pass path: GSE151179/GSE151180 |

## No-go triggers

| Trigger | Consequence |
|---|---|
| GSE151179 shows no direction-consistent RAI-gene/8-gene signal | Remove RAI-response language; keep "differentiation state" only |
| GSE184362/GSE232237 fail epithelial-cell gradient replication after QC | Downgrade single-cell validation to local-only; target lower-impact journal |
| DICER1/EIF1AX/PPM1D enrichment disappears under final mutation curation | Remove "driver-enriched subclass" from title; keep taxonomy/expression-state paper |
| Fusion counts remain limited to 1-3 cases without validated calls | Keep fusion as overlay only |
| Any analysis requires Paper 1 edits or reopens Paper 1 figures | Stop; keep this as separate advance paper |
| Team wants survival/recurrence as primary endpoint | No-go with current data |
| Team wants therapeutic target discovery claim | No-go unless drug-response or functional dependency data are added |

## Required analyses before manuscript

1. Freeze discovery denominator:
   - Primary: 482 TCGA mutation-callable tumors.
   - Report BRAF/RAS-negative count as 137/482.
   - Mention v17 178-row broader triple-negative count only as an extended integration convention.

2. Rebuild final subclass table:
   - BRAF/RAS-negative total.
   - DICER1/EIF1AX/PPM1D.
   - TERT-only.
   - fusion/alternative MAPK.
   - true driver-negative.

3. Recompute six-gene RAI score:
   - SLC5A5, TPO, TSHR, TG, PAX8, NKX2-1.
   - Compare with TDS16 and existing 8-gene score.

4. Validate clinical RAI signature:
   - GSE151179 mRNA.
   - GSE151180 miRNA optional.
   - Use exact RAI-avid/refractory metadata labels.

5. Validate cellular gradient:
   - GSE184362 malignant/thyrocyte cells.
   - GSE232237 malignant/thyrocyte cells.

6. Add limitations table:
   - TCGA outcome event scarcity.
   - no TCGA RAI-response endpoint.
   - incomplete fusion ascertainment.
   - controlled Korean advanced cohort pending EGA access.

## Claim rules

Allowed wording:

- "RAI-related differentiation state"
- "RAI-handling gene expression"
- "candidate noncanonical driver class"
- "fusion overlay"
- "true driver-negative residual class"
- "survival analysis was underpowered and exploratory"

Avoid wording:

- "RAI resistance subtype"
- "predicts recurrence"
- "improves survival stratification"
- "actionable therapeutic targets"
- "fusion-driven thyroid cancer class"
- "drug target discovery"

## Final go/no-go verdict

**Go for feasibility and first-pass manuscript planning.**

Current data are strong enough for a rigorous taxonomy paper if the paper stays honest: driver-negative prevalence, DICER1/EIF1AX/PPM1D enrichment, TERT-only small class, sparse fusion overlay, and differentiation/RAI-gene state. Outcome and therapy claims are not ready.

