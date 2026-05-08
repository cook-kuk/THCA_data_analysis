# Track 28 — DM1-high × IFN-γ-low "fail-to-induce" THCA subgroup

**Boundary.** HLA-I/II are **gene-expression modules** (transcript-level) only.
Per `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`, no allele
genotyping is performed on cancer cohorts in this track. All figures carry
the boilerplate "HLA-I/II = gene-expression module — not allele genotype.
Cancer-cohort allele genotyping out of scope."

**Setup.** Track 6′ (pan-cancer) found that DM1-high × HLA-I-low predicts BRCA
OS (HR=2.20, FDR=1.6e-4). Tracks 5/12/21/25 ruled out HLA-I LoH/LoF in
TCGA-THCA — DM1×HLA-I correlates **positively**, consistent with inflamed
induction biology. We therefore re-ask whether a residual subgroup of
**DM1-high samples that fail to mount the inflamed signature** (low IFN-γ)
exists in TCGA-THCA, and whether they are biologically or clinically distinct.

---

## 1. Quadrant definition

Median splits on DM1 (Track 5) and IFN-γ Hallmark (Track 26) yield four
quadrants on **n = 527 TCGA-THCA primary tumors**:

| Quadrant | Definition | n | % |
|---|---|---:|---:|
| Q_inflamed_dark | DM1 > med, IFN-γ > med | 166 | 31.5 |
| Q_quiet_bright | DM1 < med, IFN-γ < med | 167 | 31.7 |
| **Q_fail_to_induce** | **DM1 > med, IFN-γ < med** | **97** | **18.4** |
| Q_inflamed_bright | DM1 < med, IFN-γ > med | 97 | 18.4 |

The fail-to-induce subgroup is **97/527 (18.4%)**. By construction DM1 is
matched across the two DM1-high quadrants (mean DM1 ≈ −0.22 vs −0.21,
Cohen's d=−0.07, n.s.), so any difference vs Q_inflamed_dark is attributable
to the IFN-γ axis, not DM1 magnitude.

→ Figure F01 (quadrant scatter), Tables T01–T02.

---

## 2. Driver enrichment per quadrant

Fisher exact (each driver vs all-others, FDR-BH across all driver × quadrant
tests):

- **Q_fail_to_induce is BRAF-enriched** — BRAF 77/97 vs 217/430 outside,
  OR=3.78, p=1.3e-7, FDR=3.2e-7.
- Q_fail_to_induce is **RAS-depleted** — RAS 1/97 vs 53/430 outside,
  OR=0.07, p=2.8e-4.
- Q_fail_to_induce is also TripleNeg-depleted (OR=0.43, p=1.6e-3) and
  Other-driver-neutral.

Q_inflamed_dark: even stronger BRAF skew (141/166, OR=7.67, p=2.9e-21).
Q_quiet_bright: RAS-enriched (OR=38, p=3.6e-23).

→ Interpretation: the fail-to-induce subgroup is **BRAF-driven but
immune-silent** — not a different driver class, but BRAF tumours that fail to
recruit/activate IFN-γ. RAS tumours essentially never appear in the
DM1-high band at all (consistent with Round 5 finding RAS=98% DM1).

→ Figure F02, Table T03.

---

## 3. Clinical features per quadrant

Mann-Whitney fail vs inflamed-dark (FDR-BH):

| Feature | mean fail | mean inflamed-dark | Cohen's d | p | FDR |
|---|---:|---:|---:|---:|---:|
| age_at_dx | 46.8 | 48.9 | −0.13 | 0.31 | 0.36 |
| RAI score (v17 proxy) | 6.76 | 6.59 | +0.33 | 2.2e-3 | 3.1e-3 |
| TDS-16 score | 6.18 | 5.98 | +0.42 | 1.9e-4 | 3.3e-4 |
| HLA-I module | −0.29 | +0.46 | −1.70 | 1.9e-27 | 6.5e-27 |
| HLA-II module | −0.30 | +0.70 | −1.71 | 7.2e-26 | 1.7e-25 |

Categorical breakdowns (T04_clin_cat_*): stage / sex / histology subtype /
molecular subtype / TDS group / dediff flag / ATA risk proxy distributions
are available; fail-to-induce shows no extreme stage skew vs inflamed-dark.

Notable: fail-to-induce has **higher residual differentiation** (TDS-16
+0.42, RAI score +0.33) than inflamed-dark — i.e., these BRAF tumours have
not lost lineage as much as their inflamed counterparts.

→ Tables T04 (categorical), T05 (continuous summary), T06 (tests).

---

## 4. TERT promoter mutation per quadrant

| Quadrant | TERT-mut / total | % | Fisher OR vs other quadrants | p |
|---|:---:|---:|---:|---:|
| Q_inflamed_dark | 19/166 | 11.4 | 2.62 | **0.008** |
| Q_fail_to_induce | 7/97 | 7.2 | 1.08 | 0.83 |
| Q_quiet_bright | 9/167 | 5.4 | 0.70 | 0.46 |
| Q_inflamed_bright | 1/97 | 1.0 | 0.12 | **0.007** |

TERT promoter mutations cluster in **inflamed-dark** (the canonical
aggressive DM1 + inflamed compartment). Fail-to-induce does **not** carry an
excess TERT burden — consistent with this being a less-progressed BRAF state
that has not acquired secondary TERT yet.

→ Figure F03, Tables T07–T08.

---

## 5. HLA / cytolytic / APM modules + cell composition per quadrant

Per construction HLA-I and IFN-γ correlate positively in TCGA-THCA, so the
fail-to-induce quadrant is also low on most immune-induced modules. Top
features (fail vs inflamed-dark, MW Cohen's d):

| Feature | d | FDR |
|---|---:|---:|
| INTERFERON_GAMMA_RESPONSE | −2.02 | 2.4e-34 |
| INFLAMMATORY_RESPONSE | −2.07 | 6.3e-34 |
| CD4_T_helper marker | −1.91 | 3.5e-31 |
| ALLOGRAFT_REJECTION | −1.69 | 1.2e-27 |
| APM_score (B2M/TAP/PSMB/NLRC5/IRF1) | −1.72 | 1.6e-26 |
| HLA-II recomp | −1.67 | 4.2e-24 |
| myeloid_suppressive | −1.52 | 7.4e-24 |
| Neutrophil marker | −1.49 | 6.0e-23 |
| Macrophage_M2 | −1.52 | 1.1e-22 |
| T_reg | −1.54 | 2.8e-22 |
| Cytolytic activity (GZMA+PRF1) | −1.25 | 4.1e-18 |

→ Fail-to-induce is uniformly **immune-deserted**, not selectively
HLA-low — every adaptive- and most innate-immune marker drops together.

→ Figure F04 (HLA/cyt/APM boxplots), Figure F08 (cell-type heatmap),
Tables T10–T11.

---

## 6. Methylation per quadrant (RAI-lineage 8-gene)

Per-gene HM450 β (mean):

| Quadrant | mean_8g_β | TG | TPO | SLC5A5 |
|---|---:|---:|---:|---:|
| Q_quiet_bright | 0.287 | 0.529 | 0.560 | 0.583 |
| Q_fail_to_induce | 0.353 | 0.664 | 0.891 | 0.578 |
| Q_inflamed_bright | 0.370 | 0.616 | 0.800 | 0.573 |
| Q_inflamed_dark | 0.398 | 0.687 | 0.889 | 0.580 |

Fail-to-induce sits **between bright and inflamed-dark** on global
RAI-lineage promoter β. TG and TPO promoter methylation in fail-to-induce
already match inflamed-dark (TG 0.66 vs 0.69; TPO 0.89 vs 0.89), but mean β
is lower because other loci have not closed. Consistent with the earlier
clinical finding that fail-to-induce is a **less-advanced DM1 state**.

→ Figure F05, Table T09.

---

## 7. Survival per quadrant

Multivariate logrank across the four quadrants:

| Endpoint | n | events | χ² | p |
|---|---:|---:|---:|---:|
| OS | 526 | 16 | 0.28 | 0.96 |
| DSS | 520 | 7 | 1.90 | 0.59 |
| **PFI** | 526 | 56 | 8.69 | **0.034** |
| **DFI** | 368 | 27 | 14.88 | **0.0019** |

Cox HR (Q_fail_to_induce vs Q_inflamed_dark reference):

| Endpoint | Univariate HR (95% CI), p | Multivariate (stage_high + age + BRAF) HR, p |
|---|---|---|
| OS | 0.87 (0.22–3.50), p=0.85 | 0.71 (0.20–2.55), p=0.60 |
| DSS | 1.65 (0.23–11.7), p=0.62 | 0.99 (0.20–5.01), p=0.99 |
| PFI | 0.79 (0.41–1.54), p=0.49 | 0.79 (0.41–1.54), p=0.50 |
| **DFI** | **0.28 (0.08–0.95), p=0.041** | 0.36 (0.12–1.05), p=0.061 |

→ Fail-to-induce does **not** carry worse survival in TCGA-THCA. The only
significant univariate term, DFI, goes the *opposite* direction (HR=0.28,
fewer disease-free events, i.e., **less recurrence**), and attenuates after
stage/age/BRAF adjustment. The pan-cancer Track 6′ BRCA mortality signal
**does not transfer to TCGA-THCA**.

The TCGA-THCA event base is small (OS=16, DSS=7) and the cohort is heavily
indolent PTC, so power for OS/DSS is limited; the result is consistent with
the published TCGA-THCA outcome envelope.

→ Figures F06 (KM 4-panel), F07 (Cox forest), Tables T12–T13.

---

## 8. Gene-level DEG + Hallmark GSEA

Fail-to-induce vs inflamed-dark (Welch's t on log-norm pancan expression,
20,053 genes, FDR-BH):

**Top down (depleted in fail-to-induce):** ICOS, TIGIT, GBP5, IL2RG,
TNFSF13B (BAFF), NOD2, DAPP1, CTLA4, LTA, CSF2RB. All FDR < 1e-35. Pure
adaptive-immune / IFN-γ-axis signature loss.

**Top up (enriched in fail-to-induce):** mostly metabolic / mitochondrial
genes (see GSEA below).

Hallmark-proxy GSEA (signed MW on per-gene t):

| Hallmark | direction | n_in_set | mean_t_in_set | FDR |
|---|---|---:|---:|---:|
| INTERFERON_GAMMA_RESPONSE | down_in_fail | 31 | −10.84 | 6.5e-19 |
| INTERFERON_ALPHA_RESPONSE | down_in_fail | 17 | −8.48 | 1.3e-9 |
| ALLOGRAFT_REJECTION | down_in_fail | 14 | −10.64 | 4.6e-9 |
| INFLAMMATORY_RESPONSE | down_in_fail | 12 | −9.41 | 1.3e-7 |
| TNFA_SIGNALING_VIA_NFKB | down_in_fail | 18 | −6.20 | 1.8e-6 |
| EMT | down_in_fail | 14 | −6.29 | 5.1e-6 |
| G2M_CHECKPOINT | down_in_fail | 14 | −4.90 | 8.3e-6 |
| P53_PATHWAY | down_in_fail | 10 | −4.56 | 4.3e-4 |
| E2F_TARGETS | down_in_fail | 13 | −3.19 | 2.8e-3 |
| **FATTY_ACID_METABOLISM** | up_in_fail | 14 | +1.63 | **3.1e-3** |
| **OXIDATIVE_PHOSPHORYLATION** | up_in_fail | 16 | +1.15 | **0.015** |

→ Fail-to-induce loses inflamed-, EMT-, and proliferation-axis programs and
**gains lipid-oxidative metabolic programs**, consistent with a quieter,
less-dedifferentiated thyrocyte phenotype. This matches Round 5 GSEA's
"HT-overlap immune-up + fatty-acid-down" pattern reading in the
opposite direction.

→ Figure F09, Tables T14 (full DEG) / T15 (top-50 up) / T16 (top-50 down) /
T17 (GSEA).

---

## 9. RAI-refractoriness (post-RAI dediff) signature per quadrant

Composite from `v19_paper3_rai_dediff_axis_2026_05_06`:
RAI_dediff = −thyroid_diff_z + 0.5·myeloid_z + 0.5·HLA-II_z. Higher means
more like the post-RAI dediff state.

| Quadrant | mean RAI_dediff | n |
|---|---:|---:|
| Q_quiet_bright | −1.20 | 167 |
| Q_inflamed_bright | +0.06 | 97 |
| Q_fail_to_induce | +0.09 | 97 |
| Q_inflamed_dark | +1.16 | 166 |

Fail-to-induce vs inflamed-dark: d = −1.39, p = 7.8e-24 (T19). Inflamed-dark
is the THCA quadrant most resembling the post-RAI refractory state;
fail-to-induce sits roughly mid-axis. The pattern closes the loop: a BRAF
DM1 tumour is closest to RAI-refractory phenotype only **once it has also
recruited the inflamed/myeloid program**. BRAF DM1 alone (fail-to-induce) is
not yet RAI-dediff-like.

→ Figure F10, Tables T18–T19.

---

## 10. Cell-type composition per quadrant

(Marker-gene module proxies; no precomputed CIBERSORTx for full TCGA-THCA in
this repo.)

Inflamed-dark > fail-to-induce on every adaptive lineage marker
(CD8_T, CD4_T_helper, T_reg, B_cell, NK_cell, Dendritic) and every myeloid
marker (Macrophage_M1, Macrophage_M2, Neutrophil). Stromal_fibroblast and
Endothelial markers also drop in fail-to-induce. There is no marker module
where fail-to-induce is enriched.

→ Figure F08 (heatmap), Table T10. The CIBERSORT-based file
`project/results/v17p2/tables/cibersort_immune_fractions.tsv` covers only
n=178 of 527 with 3 cell columns (Treg / M2 / Exh-T) so we did not
substitute it; it can be cross-checked downstream.

---

## 11. Limitations

1. **Median-split definition.** "Fail-to-induce" is a median × median
   construct, not a clinically validated subtype. Tertile splits (T02)
   show similar gradients.
2. **Small triple-negative n.** Only 18 triple-negative tumours fall in
   fail-to-induce, limiting power for that driver class.
3. **Endpoint power.** TCGA-THCA OS events = 16 across 526 samples; DSS = 7.
   Survival null findings here cannot rule out a small effect; they only
   establish that there is no large mortality penalty in this cohort,
   unlike pan-cancer Track 6′ BRCA.
4. **Observational, not interventional.** The cell composition / GSEA /
   methylation comparisons are descriptive cross-section.
5. **Module-level only.** Per boundary, no allele genotyping; HLA-low here
   is RNA expression, not allele LoH or LoF.
6. **Marker-gene cell composition.** Used in absence of full CIBERSORTx
   coverage; the direction is reproducible across markers but absolute
   fractions are not estimated.

---

## 12. Implication for Paper 1

A **DM1-high × IFN-γ-low "fail-to-induce" subgroup exists** in TCGA-THCA at
roughly 18% (n=97/527), is overwhelmingly **BRAF-driven** (OR=3.8 vs
elsewhere), is **NOT TERT-mutation enriched** (7.2% vs inflamed-dark 11.4%),
and is **NOT associated with worse OS / DSS / PFI / DFI** in TCGA-THCA. The
DEG/GSEA picture is a clean immune-deserted, less-dedifferentiated, more
metabolically quiet BRAF state that has not engaged the inflamed-induction
loop.

For Paper 1 framing this means:

- **Do not promote fail-to-induce to a separate primary axis.** It is a
  state along the DM1 axis, not an independent prognostic stratum. The
  pan-cancer Track 6′ DM1 × HLA-low BRCA mortality signal does **not**
  carry over to THCA.
- The finding does **strengthen** the existing Paper 1 narrative that DM1
  in THCA is induction-biology, not LoF: if HLA-I-low DM1 were dangerous
  via immune-escape, fail-to-induce would carry the worst outcomes — it
  doesn't.
- Optional supplementary use: report fail-to-induce as evidence that the
  IFN-γ axis (not DM1 per se) is what links DM1 to the inflamed/dediff/
  TERT/RAI-refractory continuum (Tracks 26, 33, 34).
- Reviewer Q-bank value: pre-empts the predictable "is your DM1 axis just
  an inflammation surrogate?" question — Track 28 shows that ~18% of
  DM1-high tumours are uninflamed, that they are BRAF, and that they have
  no clinical penalty in THCA.

---

## Output paths

- Tables: `project/results/hla_deepdive_2026_05_08/track28_immune_cold_thca/tables/T01–T20`
- Figures: `project/results/hla_deepdive_2026_05_08/track28_immune_cold_thca/figs/F01–F10`
- Summary JSON: `project/results/hla_deepdive_2026_05_08/track28_immune_cold_thca/track28_summary.json`
- Per-sample export: `tables/T20_per_sample_full.tsv` (n=527, 51 cols incl. quadrant, modules, drivers, RAI proxy, cell-type modules, hallmarks).
- Script: `scripts/hla_deepdive_2026_05_08/track28/run_track28.py`
