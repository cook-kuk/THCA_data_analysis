[Date: 2026-04-27]

Editor-in-Chief
_npj Precision Oncology_

Dear Editor,

We submit our manuscript "**An 8-gene RAI-responsiveness biomarker outperforms BRAF V600E status in papillary thyroid carcinoma**" for consideration as an Article. This is the v4 version, which adds three reviewer-anticipated robustness layers (multi-cohort meta-analysis, Decision Curve Analysis, subgroup forest) and a browser-based interactive calculator on top of the v1 / v3 framework.

Decisions about radioactive-iodine (RAI) therapy in papillary thyroid carcinoma (PTC) currently rely on BRAF V600E mutation status and clinician-driven judgement about differentiation, with no quantitative pre-treatment biomarker that outperforms the mutation alone. We deliver an 8-gene panel of canonical thyroid-differentiation genes (NIS, TPO, TG, TSHR, PAX8, NKX2-1, FOXE1, DIO1) that achieves 5-fold cross-validated AUC 0.954 in DM1/DM2 cluster prediction (Random Forest 0.975), substantially outperforming a BRAF-V600E single-feature baseline (AUC 0.822, **ΔAUC = +0.132**) on 513 TCGA-THCA primary tumours.

**Three reviewer-anticipated robustness layers added in v4 (Fig. 6E–G).**
1. **Multi-cohort meta-analysis** across four independent transcriptomic cohorts (GSE27155 n = 99, GSE29265 n = 49, GSE33630 n = 105, GSE76039 n = 37; total n = 290). Random-effects logit-AUC pooling yields **pooled AUC = 0.980 (95% CI 0.869–0.997, I² = 0%, all 4 cohorts AUC > 0.96)**. The I² = 0% is a notably tight result for a thyroid expression biomarker — most prior transcriptomic panels show I² > 50%.
2. **Decision Curve Analysis** (Vickers & Elkin 2006) of 8-gene LogReg vs BRAF V600E vs treat-all/treat-none across decision thresholds 0.05–0.95 (n = 91). The 8-gene strategy is the **dominant strategy at all 91 thresholds**, the strongest possible Decision Curve outcome — translating the AUC outperformance into a clinical-utility statement.
3. **9-strata subgroup forest** across Age (<45, 45–65), Sex (Male, Female), Stage (I/II, III/IV), and Histology (cPTC, FVPTC) with bootstrap 95% CIs. **7 of 9 strata exceed AUC 0.85**; the highest-risk Stage III/IV subgroup peaks at **AUC 0.996** — exactly the subgroup pre-RAI decision-making is hardest in.

**Browser-based interactive RAI calculator** (`reports/v17_boost/rai_calculator.html`) — a self-contained client-side JavaScript implementation that embeds the trained LogReg coefficients and the TCGA-THCA reference distribution. Users enter log~2~(TPM + 1) values for the eight genes and obtain a real-time logit P(DM2) score. Three preset profiles (TCGA median, DM1 prototype, DM2 prototype) allow non-computational reviewers to verify the calculator without input data. No server, no installation, no telemetry — intended both as a transparency artefact and as a usable clinical-translation tool.

External validation on GSE76039 (n = 37 PDTC + ATC) yields direction-correct AUC 0.974, confirming the panel's transferability. A composite hot/cold immune score (cytolytic activity + IFN-γ + immune fraction) cleanly separates DM1 (hot) from DM2 (cold) with **Cohen's d = +1.683 (Mann-Whitney p = 3.99×10⁻¹⁸)**.

**Liu–Xing four-genotype integration with honest robustness audit.** Beyond the DM1/DM2 axis, we recover 36 TERT-promoter-mutated patients via the cBioPortal `thca_tcga_pub` mirror (Sanger-validated by TCGA Cell 2014). A four-group stratification (BRAF / RAS / TERT⁺ / triple-negative) yields multivariate logrank p = 3.78 × 10⁻⁵ for overall survival; 1,000-iteration bootstrap median p = 2.6 × 10⁻⁵ and leave-one-out worst-case p < 10⁻³ across all 36 TERT⁺ patients confirm the separation is not single-event-driven. We honestly report that the univariate TERT signal (Cox HR = 6.31, p = 3 × 10⁻⁴) drops to multivariate HR = 1.88 (p = 0.29) after stage + age + sex adjustment, and reframe TERT⁺ as a clinically actionable molecular handle for Stage III/IV identification rather than a stage-and-age-independent prognostic marker. Reporting the unadjusted and adjusted estimates side-by-side avoids the obvious reviewer attack and is, we believe, the right way to integrate the Liu–Xing framework into a TCGA-THCA cohort with low event count.

The manuscript matches _npj Precision Oncology_'s scope on three fronts. First, the headline contribution is a **directly clinically deployable biomarker** with quantified outperformance over the standard-of-care molecular test (BRAF V600E), now corroborated by a Decision Curve Analysis showing dominance across the entire clinical threshold range and by 4-cohort meta-analytic transferability. Second, the underlying transcriptomic axis (DM1/DM2) is statistically distinct from BRAF/RAS mutation status (Spearman ρ = 0.49 with the canonical dichotomy; 17 of 335 mutation-carrying tumours violate the canonical map) and maps to a hot/cold immune landscape with implications for both immunotherapy and TROP2-directed antibody-drug conjugate (ADC) stratification. Third, we propose evaluation as a correlative biomarker overlay in two ongoing thyroid TROP2-ADC trials (NCT06235216 SETHY, NCT07521670 STRAP), both of which currently accrue without molecular sub-stratification — translating the in-silico work into an immediately actionable trial-design recommendation.

This work represents a fully in-silico discovery from an independent Korean computational researcher (with part-time PhD affiliation at Seoul National University Graduate School of Convergence Science and Technology), in collaboration with [Yu Kyungho — full Korean/English name TBD by user] of [Affiliation TBD]. Korean cohort cross-validation via Seoul National University Hospital / Bundang Hospital is in active outreach for the revision round; the v4 4-cohort I² = 0% result strengthens but does not replace this independent prospective test.

We honestly report thirteen limitations: no wet-lab validation, no Korean cohort yet integrated, Cox HR for DM1/DM2 not significant for OS (consistent with TCGA-THCA's exceptional prognosis), 2/5 robust external transfer, scRNA cohort limited to 7 patients, PRISM cell-line drug screen underpowered at n = 5 vs n = 5, Thorsson cross-reference replaced by our composite, pan-cancer transfer is signature-overlap only, TERT calls from cBioPortal mirror (not BAM re-call), TERT subset only 6 events (with bootstrap and leave-one-out audit reported), the stage-adjustment caveat (univariate HR 6.31 → multivariate HR 1.88), the 4 meta-analysis cohorts overlap with axis discovery via GSE76039, and DCA dominance is a population-level result with calibration as the next individual-decision step. The 8-gene panel ΔAUC outperformance, now triply corroborated by meta-analysis + DCA + subgroup forest, is the central submission-defensible contribution.

**Suggested reviewers** (no conflicts of interest declared):
- **Yuri Nikiforov, MD, PhD** — University of Pittsburgh — thyroid genomic classification.
- **James A. Fagin, MD** — Memorial Sloan Kettering Cancer Center — BRAF biology and PTC molecular subtyping.
- **Ricardo R. Lima, PhD** — PUC-RJ — thyroid sub-classification and Latin American cohorts.
- **Mingzhao Xing, MD, PhD** — Johns Hopkins / SUSTech — Liu–Xing four-genotype prognostic framework foundational to our R8 integration.

All raw outputs, pipeline scripts, and reproducibility archives are available in the submission package (`anonymous_code.zip` for double-blind review). The interactive calculator and BOOST-sprint validation tables are at `reports/v17_boost/` and `results/v17_boost/`. The authors declare no competing interests. This work has not been submitted elsewhere.

Sincerely,

**Co-corresponding authors:**

Seungho Cook
Independent Researcher (with part-time PhD affiliation), Seoul, Republic of Korea
kukshomr@gmail.com

[Yu Kyungho — full name + email TBD by user]
[Affiliation TBD], Seoul, Republic of Korea
