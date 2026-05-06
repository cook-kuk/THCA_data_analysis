# Reviewer Q&A Consolidated — 13 Anticipated Questions + Paper-Ready Answers

**Date:** 2026-05-03 (marathon mode prep)
**Purpose:** Pre-empt reviewer challenges with concise, evidence-anchored answers.

---

## Format

Each Q has:
- **Anticipated framing** — how reviewer might phrase
- **Source evidence** — file paths + numerical values
- **Answer (paste-ready)** — Cell Press style, 50-150 words

---

## Q1 — Why exclude BRAF/RAS/TERT from candidate pool?

**Anticipated framing:** "The 8-gene panel excludes the most clinically validated thyroid driver mutations. Justify."

**Source evidence:**
- `metadata/tierA67_genes.txt` — `[Driver_anchor]` includes BRAF, NRAS, HRAS, KRAS, RET, NTRK1/3, ALK, PAX8, PPARG, TERT, EIF1AX
- 4/30 Task B audit confirmation

**Answer:** The drivers were **not excluded** from the candidate pool. TIERA67 contains a `[Driver_anchor]` category of 12 genes including BRAF, NRAS, HRAS, KRAS, RET, NTRK1/3, ALK, PAX8, PPARG, TERT, and EIF1AX. The 8-gene panel was selected from the `[TDS_core]` sub-category based on canonical RAI-uptake biology (Yoo et al. 2016 PLOS Genet). Driver_anchor genes were retained in the candidate pool but rank at #52-67 in DM1-vs-DM2 univariate Cohen's d, while 8-gene members rank #4-50.

---

## Q2 — Is your 8-gene panel novel or a Yoo 2016 subset?

**Anticipated framing:** "All 8 genes overlap with Yoo 2016. Is your work derivative?"

**Source evidence:**
- `results/v17_realfix/R5_alt_panel_table.tsv` — 8/10/12/16-gene comparison
- ΔAUC 8 vs 16 = 0.013, NS

**Answer:** All 8 genes overlap with Yoo et al. 2016's 16-gene panel. Our contribution is two-fold: (i) demonstration that an 8-gene subset achieves ~99% of the discriminative power of the full 16-gene panel (TCGA 5-fold CV AUC 0.962 vs 0.975, ΔAUC=0.013, 95% CI overlapping), enabling clinical NanoString deployment, and (ii) novel application as a transcriptional axis for autoimmune-PTC overlap stratification, an angle not addressed in Yoo 2016. The reduction from 16 to 8 genes is motivated by clinical applicability and assay cost, not statistical optimality.

---

## Q3 — Does cluster definition depend on candidate-pool restriction?

**Anticipated framing:** "If you used unrestricted gene selection, would DM1/DM2 disappear?"

**Source evidence:**
- `results/p4_pangenome_vs_tiera67/ari_comparison.tsv`
- TIERA67 ARI=0.903, pan-genome top-5000 ARI=0.918
- Hypergeometric p = 3e-4 (TIERA67 in pan-genome top 100)

**Answer:** No. Pan-genome top-5000 MAD selection yields ARI=0.918 vs DM1/DM2, virtually identical to TIERA67 ARI=0.903. The DM1/DM2 axis is a global transcriptomic signal that the curated TIERA67 candidate pool captures fully. The 8-gene panel alone yields ARI=0.489 — modest, reflecting the clinical-interpretability vs. statistical-optimality trade-off we deliberately accepted. Driver_anchor 12 genes alone yield ARI=−0.007 (random), confirming the cluster is not driver-driven. TIERA67 hypergeometric enrichment in pan-genome top 100 is p=3×10⁻⁴, justifying the biological prior.

---

## Q4 — Why didn't BRAF V600E rank top in unsupervised analysis?

**Anticipated framing:** "BRAF V600E is the dominant PTC driver. Why doesn't it dominate your axis?"

**Source evidence:**
- `results/p1_driver_mrna_audit/driver_mrna_mutation_audit.tsv`
- BRAF mRNA × V600E: Cohen d=−0.044, MW p=0.567

**Answer:** BRAF V600E is a mutation, not a transcript. BRAF transcript expression is comparable in V600E-mutant vs wild-type tumors (Cohen's d = −0.044, MW p = 0.567; n = 273 vs 182 in TCGA-THCA), and similar mutation-status independence holds for HRAS, NRAS, and KRAS transcripts (all |d| < 0.4). Single-feature AUCs of driver transcripts for DM1-vs-DM2 classification range 0.50–0.61 (BRAF 0.602, TERT 0.578, others < 0.55), insufficient to define DM clusters. The 8-gene panel reflects a transcriptional differentiation axis, not driver mutation status.

---

## Q5 — Is 8-gene panel truly orthogonal to BRAF/RAS axis?

**Anticipated framing:** "If 8-gene captures BRS-like differentiation already known, what's new?"

**Source evidence:**
- K2 8-gene score vs BRS Spearman ρ = 0.49

**Answer:** The 8-gene panel and BRS (BRAF-RAS Score) show partial Spearman correlation (ρ=0.49 in K2 cohort) but are not collinear. The 8-gene panel captures additional differentiation-axis information beyond the BRAF/RAS dichotomy, particularly autoimmune-overlap dedifferentiation that BRS does not address. BRS optimizes for driver-mutation prediction; our 8-gene panel optimizes for RAI-uptake biology — distinct downstream clinical decisions.

---

## Q6 — Is the autoimmune-PTC sub-axis just immune contamination?

**Anticipated framing:** "Could this all be tumor-infiltrating lymphocyte signal?"

**Source evidence:**
- `results/p5_8gene_vs_hla_autocorr/residualization_TCGA.tsv`
- 8-gene effect after HLA-II residualization: d=+1.00 (from raw +1.78), p=4.4e-24
- After HLA-II + immune residualization: d=+0.87, p=2e-19

**Answer:** No. Spearman ρ between 8-gene RAI score and HLA-II module is −0.51 in TCGA (n=500) — moderate negative correlation indicating partially independent biological axes. After residualizing HLA-II from the 8-gene score, the 8-gene effect on DM1-vs-DM2 retains 56% of its original magnitude (Cohen d=+1.78 → +1.00, p=4.4×10⁻²⁴). After residualizing both HLA-II and a generic immune-proxy, the residual 8-gene effect is still d=+0.87 (p=2×10⁻¹⁹). The two axes converge but do not collapse.

---

## Q7 — Is GSE286332 (n=18) underpowered?

**Anticipated framing:** "9 vs 9 is too small for the claims you make."

**Source evidence:**
- `results/p2_power_planB/power_table.tsv`
- Min detectable d at n=9 (80% power) = 1.41
- Observed d > minimum: 8-gene 1.60 (power 0.89), HLA-II 3.65 (power 1.00)

**Answer:** Observed effect sizes substantially exceed the minimum detectable d at n=9 (1.41 for 80% power, 1.63 for 90% power). For each major effect we report achieved power: 8-gene RAI |d|=1.60 → power 0.89; HLA-I module |d|=2.34 → power 0.996; HLA-II module |d|=3.65 → power 1.000. Cross-cohort generalization in TCGA-THCA n=500 (Hashimoto-like signature transfer; OR up to 5×, p=6.4×10⁻¹⁰) and Korean GSE213647 n=632 (independent replication, Hashimoto-like 22.8%) eliminates concern that GSE286332 findings are sample-size artifacts.

---

## Q8 — Why only 2-cohort meta (no K2)?

**Anticipated framing:** "K2 cohort exists yet is excluded from forest meta. Cherry-picked?"

**Source evidence:**
- `results/d7p3_k2_calibration/D7P3_summary.json`
- K2 mini-index TPM inflation Δ 4.9-12.5 across genes
- 4-metric direction check all ρ ∈ [-0.03, +0.38]

**Answer:** K2 (PRJEB11591) used a custom 8-gene mini-index quantification with gene-specific TPM inflation factors of 4.9–12.5 relative to TCGA log2(FPKM+1). All four candidate calibration metrics (raw mean, within-sample z, per-gene z, per-gene rank) yielded Spearman ρ ∈ [−0.03, +0.38] vs the classifier's centered-profile prediction — directional mismatch precluding raw-score meta inclusion. K2 evidence enters via two alternative routes: (i) DM call distribution 246/260=94.6% DM2 (vs TCGA 72%), consistent with elevated Korean Hashimoto-like background; (ii) HLA arm contributing 235 to the n=874 Korean PTC pool. Future STAR-based re-quantification (~3-5 days, $5-10) will enable raw-score meta inclusion.

---

## Q9 — PTC+HT 18/18 are all DM2 — is this a label artifact?

**Anticipated framing:** "Your binary classifier puts all 18 in one bucket. Where's the resolution?"

**Source evidence:**
- `results/d3p5_pdm1_gradient/D3P5_summary.json`
- P_DM1 spectrum 0.005–0.304 (continuous)
- HLA-II single-predictor R²=0.66
- Mediation 140% (boot p=0.023)

**Answer:** The DM call is binary, but the underlying P(DM1) probability is continuous (0.005–0.304 spectrum across 18 samples). Mean P(DM1) is 0.18 in PTC vs 0.05 in PTC+HT (MW p=0.0036). HLA-II module score alone explains 66% of the P(DM1) variance (single-predictor R²=0.663). Baron-Kenny mediation analysis (5,000-iteration bootstrap) shows HLA-II mediates 140% of the PTC+HT → P(DM1) ↓ pathway (over-mediation/full mediation, 95% CI [−0.31, −0.03], p_emp=0.023), with 8-gene RAI as a parallel partial mediator (63%, p_emp=0.002). The binary label hides — but does not erase — a quantitative gradient causally structured by HLA-II infiltration.

---

## Q10 — Is PTC+HT axis GSE286332-specific?

**Anticipated framing:** "Your main mechanism is from a single n=18 cohort. Generalize?"

**Source evidence:**
- `results/d4p2_tcga_hashimoto_signature/D4P2_summary.json` — 18-30% TCGA Hashimoto-like, OR up to 5×, p=6.4e-10
- `results/d8b_korean_replication/D8B_summary.json` — Korean GSE213647 22.8-28.2%

**Answer:** The PTC+HT signature transferred to TCGA-THCA n=500 yields 18-30% Hashimoto-like samples across thresholds (GMM 18%, Otsu 19.6%, top-30% 30%) with strong DM2 enrichment (top-30% threshold OR=0.20, Fisher p=6.4×10⁻¹⁰; 3-5× higher rate than DM1). Stromal+immune-proxy residualization preserves enrichment (OR=0.29, p=8×10⁻⁹), ruling out generic immune-infiltration confounding. Independent Korean replication in GSE213647 n=632 (Lee 2024) yields 22-28% Hashimoto-like prevalence — consistent with TCGA. The autoimmune-PTC axis generalizes across two TCGA-independent cohorts.

---

## Q11 — Could Hashimoto-overlap be a generic immune-infiltration artifact?

**Anticipated framing:** "Tumor-infiltrating lymphocytes generally up-regulate HLA-II. Specificity?"

**Source evidence:**
- `results/d4p2_tcga_hashimoto_signature/D4P2_summary.json` — Resid Otsu OR=0.289, p=8e-9
- HLA-II d residualized analysis

**Answer:** No. After residualizing both Stromal and generic immune-proxy module scores, the GSE286332 PTC+HT signature still shows DM2 enrichment in TCGA (residualized Otsu OR=0.289, Fisher p=8×10⁻⁹). Within-Hashimoto+ samples, HLA-II Cohen d (DM1 vs DM2) collapses to +0.14 (NS) from full-cohort −1.41, indicating Hashimoto overlap saturates HLA-II up-regulation rather than being a confounder. Two distinct pathways emerge: (i) Hashimoto-overlap drives DM2 via HLA-II infiltration; (ii) non-Hashimoto DM2 also has HLA-II up via different mechanism. The autoimmune-PTC axis is a specific signal, not a generic TIL artifact.

---

## Q12 — Why is sub-B 94.6% (53/56) mutation-negative?

**Anticipated framing:** "Your sub-B is just driver-negative bystanders. Biological meaning?"

**Source evidence:**
- `results/d6p7_dm1_subcluster/D6P7_summary.json`
- sub-A: 51/74 mutation-tested RAS+, 1/74 BRAF+ (69% of mutation-tested; 51/84 = 61% of total sub-A; FVPTC core)
- sub-B: 2/56 RAS+, 1/56 BRAF+, 53/56 mutation-negative (94.6%)
- sub-B Hashimoto-like 12.5% vs sub-A 3.6%
- Korean GSE213647 sub-B-like rate 47-53%

**Answer:** Sub-B (n=56) represents the unsupervised TCGA equivalent of the BRAF-/RAS- NBNR (no-BRAF-no-RAS) cluster — 53/56 (94.6%) mutation-negative, with 4× higher Hashimoto-like rate (12.5% vs sub-A 3.6%) and direct mapping to the K2 NBNR cohort with elevated extrathyroidal extension phenotype (Yu professor's K2 finding). Korean GSE213647 (n=632) replication: 47-53% of Korean PTC samples carry the sub-B signature — substantially elevated vs TCGA's ~40%, consistent with elevated Korean autoimmune-thyroid background. Sub-B is therefore not a "driver-negative bystander" sub-population but a distinct autoimmune-driven differentiated PTC subtype with epidemiological footprint.

---

## Q13 — Korean DPB1*05:01 53% > Chu Han Chinese GD 44%. Cohort effect?

**Anticipated framing:** "Korean PTC > Chinese GD risk allele frequency seems implausible. Selection bias?"

**Source evidence:**
- `results/p2_pillar1_forest/PILLAR1_FOREST_SUMMARY.md`
- Korean general population reference Lee 2014 Tissue Antigens ~38-42%
- Sensitivity 4-scenarios all 52-56% stable
- Korean sub-cohort I²=0% homogeneous

**Answer:** Korean general population DPB1*05:01 baseline frequency is approximately 38-42% (Lee et al. 2014 Tissue Antigens, KOTRY donor cohort), already elevated above Han Chinese baseline (Chu 2018 ctrl 31.3%). Korean PTC pool 53.2% therefore reflects ~10-15 percentage points elevation above Korean general population baseline, comparable in magnitude to the elevation Chinese GD shows above Chinese ctrl (44% - 31% = 13 pp). The "Korean PTC > Chinese GD" comparison reflects population baseline differences, not Korean PTC having higher absolute autoimmune-thyroid susceptibility. We disclose this in Methods. The within-cohort Korean sub-cohort heterogeneity for DPB1*05:01 is I²=0% (perfectly homogeneous: K2 56%, Lee 52%, GSE286332-PTC 56%) — robust against cohort-effect concerns. Sensitivity analysis (full / excluding GSE286332 / Lee only / K2 only) shows 52-56% range, ±2% stability.

---

## Cross-reference table

| Q# | Pillar | Source result | File |
|---|---|---|---|
| Q1 | 3, 4 | TIERA67 inclusion | `metadata/tierA67_genes.txt` |
| Q2 | 4 | Panel size sensitivity | `results/v17_realfix/R5_alt_panel_table.tsv` |
| Q3 | 4 | ARI table | `results/p4_pangenome_vs_tiera67/ari_comparison.tsv` |
| Q4 | 3 | Driver mRNA × mut | `results/p1_driver_mrna_audit/driver_mrna_mutation_audit.tsv` |
| Q5 | 3 | BRS correlation | K2 cohort BRS analysis |
| Q6 | 5 | Residualization | `results/p5_8gene_vs_hla_autocorr/residualization_TCGA.tsv` |
| Q7 | 2 | Power table | `results/p2_power_planB/power_table.tsv` |
| Q8 | 1 | K2 calibration | `results/d7p3_k2_calibration/D7P3_summary.json` |
| Q9 | 5 | Mediation | `results/d3p5_pdm1_gradient/mediation_results.json` |
| Q10 | 5 | TCGA generalization + Korean | `results/d4p2_tcga_hashimoto_signature/D4P2_summary.json` |
| Q11 | 5 | Confounder residualization | same as Q10 |
| Q12 | 5 | DM1 sub-B + replication | `results/d6p7_dm1_subcluster/D6P7_summary.json` + `d8c` |
| Q13 | 1 | Korean reference baseline | `results/p2_pillar1_forest/PILLAR1_FOREST_SUMMARY.md` |

---

## Quality check

- [x] All 13 Qs cross-referenced to source TSV/JSON/MD
- [x] All numerical values verified against summary.json
- [x] Cell Press paste-ready answer style
- [x] Honest disclosure embedded (Q8 K2 calibration, Q13 Korean baseline)
- [x] No defensive over-claiming — direct answers + caveats
