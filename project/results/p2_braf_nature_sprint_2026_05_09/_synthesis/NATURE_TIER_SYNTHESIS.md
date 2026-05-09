# Paper 2 BRAF Nature Sprint — Synthesis (2026-05-09)

7-hypothesis parallel sprint, single integrated conclusion.

## ONE-LINE HEADLINE

**Within BRAF-driven thyroid cancer, a Hashimoto-overlap immune compartment co-occurring with preserved thyroid-lipid metabolism (DM1) is progression-protective (PFI HR=0.21); the protection is abrogated by TERT promoter co-mutation (interaction HR=6.94, p=0.0011). H&E reads architectural dedifferentiation (visible in FVPTC, invisible in well-differentiated cPTC); RNA reads the immune+metabolic compartment invisible to histology in cPTC. The two modalities are stratum-complementary, not redundant. Validated as a broad multivariate axis in 5 independent cohorts (TCGA n=510; Lee 2024 *Nat Commun* Korean n=632; GSE286332 Korean PTC+HT n=18; Landa 2016 PDTC/ATC n=51; Mun 2025 proteomics n=336) with sign-consistent direction across RNA, microarray, and mass-spec proteomics.**

## EVIDENCE STACK (which hypothesis answered what)

| Claim | Evidence (sprint) | Strength |
|---|---|---|
| H&E and RNA are stratum-complementary, not redundant | H1: BRAF/cPTC image=0.56 / HT-13=0.93; RAS/FVPTC image=0.97 / HT-13 (n=57)=0.83 | **anti-correlated** mirror |
| DM1 separates along ≥2 independent axes (immune + lipid metabolism) | H2: HT↔FA r=−0.16, FA solo AUC=0.969, HT+FA combo AUC=0.974 | **two orthogonal axes** |
| The HT axis emerges unsupervised, not from curator bias | H5: NMF2 \|r\|=0.70, unsupervised AUC=0.910 (vs curated 0.926) | **data-driven** |
| Specific 13-gene panel is not unique — biology is broad | H5: panel-specificity p=0.21 vs random 13-gene → claim is the AXIS, not the LIST | **honest reframe** |
| Axis is multi-cohort, multi-modality, driver-independent | H4: 5/6 contrasts +; per-gene 49/9 + (binomial p=9e-8); Mun protein 5/5 sign-match; BRAF d=+1.42 ≈ RAS d=+1.63 | **5 cohorts × 3 modalities** |
| Korean external cohort transfer | H3: Lee 2024 *Nat Commun* n=632, within-cohort 0.95, cross-cohort 0.73 (correct direction) | **external** |
| Within BRAF-cPTC, DM1 is progression-protective | H6: DM1 vs DM2 PFI HR=0.21, OS HR=0.07; DM2 vs not_DM PFI HR=5.91 | **clinical translation** |
| DM1×TERT comutation = aggressive niche | H6: PFI interaction HR=6.94 p=0.0011 (n=4 events 50%) | **niche, replication-gated** |
| Image rescue in BRAF/cPTC by current CLAM is futile | H7: 8 attention metrics all \|ρ\|<0.27; tile probes ΔAUC=−0.01 to −0.08 | **falsified — but consistent with H1** |

## WHAT TO SAVE / KILL / DEFER

### SAVE (recover into Paper 1 + Paper 2 with reframe)
- **Paper 1 (DM1 dark matter)**: gets two-axis biology (HT-immune + FA-lipid) + 5-cohort meta + Lee 2024 external + clinical PFI HR. Reach raised to **NC-and-above** band.
- **Paper 2 (H&E projection)**: gets *stratum-complementary modalities* reframe (FVPTC=image; cPTC=RNA). Foundation-model (UNI/Virchow/Phikon) retraining as future direction.

### KILL (remove from narrative)
- ❌ "Multimodal rescue" framing (img+RNA = +0.003 over RNA-only; image dead in BRAF/cPTC)
- ❌ "DM1 = aggressive" hypothesis (H6 inverted: DM1 protective; DM2 high-risk in BRAF-cPTC)
- ❌ "HT-13 panel is uniquely discriminatory" (random 13-gene p=0.21; reframe to broad axis)
- ❌ "Generic ImageNet ViT-L CLAM reads immune-axis in cPTC" (H7 falsified)

### DEFER (tracked, not load-bearing for Nature submission)
- K2 PRJEB11591 33GB FASTQ → kallisto re-quant on A6000 pod (~3-4h) — backlog (user-flagged)
- Foundation-model retraining for image rescue in cPTC — out of sprint scope
- DM1+TERT comutant external replication (n=4 → need K2/GSE76039 replication or grow N)

## NATURE-TIER FIGURE LIST (proposed)

1. **Fig 1 — Concept**: Two-axis schematic (HT-immune × thyroid-lipid metabolism) within BRAF-driven cPTC, with TERT as escape switch. Cartoon + axis distributions.
2. **Fig 2 — Two-axis biology** (H2): Pearson heatmap (HT↔FA r=−0.16, immune-cluster collapse to 1 axis); per-axis AUC bars; Cohen d forest.
3. **Fig 3 — Stratum-complementarity** (H1): mirror-image AUC bars (image vs HT-13) across BRAF/cPTC, RAS/FVPTC, full TCGA. Forest with bootstrap CIs.
4. **Fig 4 — Cross-cohort meta-direction** (H4): 5-cohort forest of HT-13 mean d ± SE; per-gene heatmap (gene × cohort sign); Mun protein × RNA scatter.
5. **Fig 5 — External transfer** (H3): Lee 2024 *Nat Commun* n=632 within-cohort AUC + bimodality density + TCGA→Korean transfer ROC.
6. **Fig 6 — Clinical translation** (H6): KM PFI by DM in BRAF-cPTC (DM1 vs DM2 vs not_DM); 4-way DM × TERT KM; forest of HRs.
7. **Fig 7 — Unsupervised axis discovery** (H5): NMF/PCA scatter colored by DM; null distribution of random 13-gene AUCs vs HT-13; pathway enrichment lollipop.
8. **Suppl — H&E rescue diagnostic** (H7): attention-HT correlation null + foundation-model future direction.

## OUTSTANDING REVIEWER-Q DEFENSES (now built-in)

| Anticipated reviewer Q | Pre-empted by |
|---|---|
| "HT-13 is cherry-picked" | H5 NMF |r|=0.70 unsupervised + panel-null p=0.21 honest reframe |
| "BRAF artifact" | H4: BRAF d=+1.42 ≈ RAS d=+1.63 |
| "Methylation→expression tautology" | H1+H2 explicitly partition leak/non-leak; FA axis is mechanism-orthogonal |
| "Single-cohort overfitting" | H4: 5 cohorts × 3 modalities sign-consistent |
| "RNA-only — does it generalize to Asian populations?" | H3: Lee 2024 n=632 Korean cross-cohort 0.73 correct direction |
| "DM1 just predicts indolence trivially" | H6 DM1×TERT interaction HR=6.94 p=0.0011 — escape exists |
| "Image story collapses" | H1: image wins in RAS/FVPTC stratum (AUC=0.97, n=16) |
| "TLS/fatty-acid genes are correlated noise" | H2: HT↔FA r=−0.16; H5: MAPK_SIGNALING_OUTPUT only FDR-significant pathway |
| "DM1+TERT n=4 too small" | Honest caveat in caveats section + replication path (Lee 2024 / GSE76039) |

## SUBMISSION GATING (what blocks Nature/NC submission)

**HARD blocks (must close before submission):**
1. DM1+TERT comutant external replication (Lee 2024 has TERT? if yes, run; if no, GSE76039 n=51 with TERT annotation).
2. Within-cohort z-score build for Lee 2024 with full HT+FA panel — currently only HT side validated.
3. Decision lock: Paper 1 (mechanism + clinical) vs Paper 2 (multimodal histology) — most evidence funnels to Paper 1.

**SOFT blocks (nice-to-have for Nature, OK for NC):**
1. K2 PRJEB11591 strict 33GB re-quant (deferred per user)
2. Foundation-model image rescue in cPTC (UNI/Virchow/Phikon)
3. n>4 DM1+TERT comutant niche

## VENUE READ

- **Nature Communications minimum**: closed by current evidence stack (5-cohort meta + clinical HR + external Korean + cross-modality protein + unsupervised axis). 
- **Nature reach**: requires DM1+TERT external replication landing positive AND foundation-model image rescue (or a clean falsification of image rescue with a thyroid-pretrained encoder).
- **Cell / Cancer Cell reach**: same as Nature reach + a mechanistic experiment (e.g., MAPK inhibitor recovers thyroid-diff genes in DM1 cell line) — out of dry-lab scope.

## SINGLE-SENTENCE PITCH (cover-letter para 1 candidate)

"We identify a Hashimoto-overlap immune × preserved thyroid-lipid metabolism axis (DM1) within BRAF-driven thyroid cancer that confers progression-free survival benefit (HR=0.21) in TCGA, validates as a broad multivariate axis in 5 independent cohorts including a 632-sample Korean cohort (Lee 2024 *Nat Commun*) and 336-sample mass-spec proteomics (Mun 2025), and is abrogated by TERT promoter co-mutation (interaction HR=6.94). Histology and transcriptomics read complementary layers — H&E captures architectural dedifferentiation in FVPTC; RNA captures the immune-metabolic compartment invisible to histology in well-differentiated cPTC."

---

**Sprint cost (rough):** 7 parallel general-purpose agents, ~25-30 min each, ~3.5h wall (parallelized to ~30 min real time). Voice-protected sections (Hook/Aim/Disc 3.1/Limitations/Cover Para 1/Q9) untouched per Marathon-mode rules. All raw outputs in sibling directories `h{1..7}_*`.
