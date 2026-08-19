# AUDIT 05 — external cohort audit

Generated 2026-08-06. Sources are the open supplementary tables of the three published
cohorts; nothing here required a data-access application.

---

## 1. Boucai 2023 — resolution of the n = 23 versus n = 24 discrepancy

**Verdict: n = 23 is correct for the driver analysis, and the excluded patient is ERAI_3.**

### Patient flow

| Stage | n | Source |
|---|---|---|
| Published clinical cohort, exceptional responders | 8 | Table S1, rows ER# 1–8 |
| Published clinical cohort, non-responders | 16 | 1:2 matched design stated in the paper |
| **Published total** | **24** | |
| Exceptional responders with somatic mutation data | **7** | Table S2, sheet "Somatic mutations in ER and NR" |
| Non-responders with somatic mutation data | **16** | same sheet |
| **Included in the driver meta-analysis** | **23** | |
| Excluded | 1 | **ERAI_3 — present in Table S1 clinical listing, absent from the mutation sheet** |

### Sample identifiers, verbatim

```
ER (7):  ERAI_1_M  ERAI_2  ERAI_4_M  ERAI_5  ERAI_6_M  ERAI_7  ERAI_8
NR (16): NRAI_1A  NRAI_1B_M  NRAI_2A  NRAI_2B  NRAI_3A  NRAI_3B  NRAI_4A  NRAI_4B
         NRAI_5A  NRAI_5B  NRAI_6A  NRAI_6B  NRAI_7A  NRAI_7B  NRAI_8A  NRAI_8B
```

The ER series runs 1, 2, 4, 5, 6, 7, 8 — **ERAI_3 is the missing identifier**, which matches
the 8 → 7 drop exactly.

### Are NRAI_1A and NRAI_1B one patient or two?

This determines whether the non-responder arm contributes 8 patients or 16, and therefore
whether the analysis treated within-patient replicates as independent. Resolved by comparing
truncal driver mutations within each A/B pair:

| Pair | A driver | B driver | Verdict |
|---|---|---|---|
| NRAI_1A / NRAI_1B_M | NRAS Q61R | **BRAF V600E** | different patients |
| NRAI_2A / NRAI_2B | NRAS/HRAS Q61K | Q61R | different patients |
| NRAI_3A / NRAI_3B | N486S | **KRAS G12C** | different patients |
| NRAI_5A / NRAI_5B | BRAF V600E | BRAF V600E | same driver, different passenger set |
| NRAI_8A / NRAI_8B | BRAF V600E | BRAF V600E | same driver, different passenger set |

Three of five pairs carry **discordant truncal drivers**. Two tumours from one thyroid
carcinoma patient essentially never differ at the truncal driver, so A and B are distinct
patients. This is also what the paper's stated 1:2 matched design predicts: each exceptional
responder is matched to two non-responders, labelled A and B.

**Consequence:** the meta-analysis unit of analysis was correct — 23 patients, each
contributing one sample. No within-patient replication was treated as independent. The
`_M` suffix denotes a metastasis specimen and no patient appears twice.

### Required wording

Any figure caption or Methods sentence using n = 23 must state:
"one of eight exceptional responders (ERAI_3) is listed in the clinical table but has no
somatic mutation data in the supplementary mutation table and was therefore excluded from
the driver classification."

---

## 2. Driver meta-analysis — strength of claim

Three cohorts, 294 analysable patients, random-effects pooling of the BRAF/RAS-negative
compartment versus all other drivers:

| Cohort | n | refractory in BRAF/RAS-negative | χ² P (3 driver classes) |
|---|---|---|---|
| Siraj 2022 | 158 | 26/56 (46%) | 0.675 |
| Zhang 2026 | 113 | 27/56 (48%) | 0.109 |
| Boucai 2023 | 23 | 2/6 (33%) | 0.034 |
| **Pooled (random effects)** | **294** | **OR 0.75, 95% CI 0.30–1.87** | **P = 0.53, I² = 63%** |

**The pooled estimate is null, but this does not license a strong claim.** The confidence
interval spans 0.30 to 1.87, so moderate associations in either direction remain compatible
with the data. Heterogeneity is substantial (I² ≈ 63%) on only three studies, where τ² is
itself poorly estimated. Endpoint definitions differ across cohorts: Siraj uses a
seven-criterion refractoriness definition, Zhang uses scan-based avidity, Boucai uses RECIST
exceptional response. And the smallest cohort, Boucai at n = 23, is the only one reaching
nominal significance and points in the **opposite** direction.

### Wording that must be used

> Across three clinically heterogeneous cohorts, driver class did not consistently
> discriminate radioiodine-refractory disease, with substantial between-study heterogeneity.

### Wording that must not be used

- "the driver axis is dead" / "드라이버 축은 확실히 죽었습니다"
- "all three cohorts were consistent" / "세 코호트 모두 일관"
- "driver mutations do not matter"
- "the meta-analysis disproves driver biology"
- any phrasing implying the meta-analysis compared driver class **against** the eight-gene
  score in the same patients — it did not, and no such head-to-head exists

---

## 3. Zhang 2026 — strength of claim

Within 113 advanced-DTC patients and one clinical endpoint, the proteomic consensus subtype
association (χ² P = 1.4 × 10⁻⁷; CC1 21% → CC2 57% → CC3 87% refractory) is far stronger than
the driver-class association (P = 0.109) in the same patients.

Using one cohort and one endpoint removes cohort-level confounding. It does **not** remove
patient-level confounding: disease burden, metastatic site, histology, stage, prior
treatment, cumulative activity, tumour cellularity and proteomic batch could all differ by
subtype, and the supplementary table does not carry enough of these to adjust properly.

### Required wording

> Within the same advanced-DTC cohort and clinical endpoint, proteomic subtype showed a
> substantially stronger unadjusted association with refractory disease than driver class.

### Prohibited wording

- "cannot be explained by confounding"
- "proves expression state is causal"
- "directly validates the eight-gene panel"
- "demonstrates superiority of our panel"

The CC subtypes are whole-proteome consensus classes. They are **orthogonal molecular-state
evidence**, not a validation of the eight-gene transcriptomic score.

---

## 4. GSE173248 — inventory classification

Discovered by phenotype-based GEO re-mining rather than by searching for "radioiodine".

| Field | Value |
|---|---|
| n | 22 samples |
| Modality | circulating microRNA |
| Phenotype field found | `rhtsh stimulation: Complete Remission` versus `Persistent/Recurrent Disease` |
| Words "radioiodine", "RAI", "iodine" in the record | absent |
| Endpoint | post-treatment clinical disease status under rhTSH stimulation |
| Radioiodine directness | indirect |
| Eight-gene panel compatibility | **none** — miRNA, not mRNA |

**Classification:** post-treatment outcome-labelled circulating-miRNA cohort.

It may appear in the data-landscape inventory as an example of a phenotype hidden from
keyword search. It must **not** appear in the figure of eight-gene panel tests. Because the
GEO record does not state that all patients received radioiodine ablation, it must not be
described as a "RAI ablation cohort" until that is verified from the source publication.
