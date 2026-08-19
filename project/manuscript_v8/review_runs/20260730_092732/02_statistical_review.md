---
title: Statistical review
audit_date: 2026-07-30
auditor: manuscript-audit-agent (statistical reviewer mode)
---

# 02 — Statistical review

## 1. Survival analysis

### 1.1 Overall survival — TCGA and MSK

**Design:** Cox proportional-hazards with DM1 status (binary) as predictor; outcome = overall survival.

**Issues identified:**

- **Event count (TCGA):** TCGA-THCA primary-tumour OS has only 16 events in n=504 (3.2% event rate). With one binary predictor, the rule of 10 events per variable is barely satisfied. The wide confidence interval (HR 2.30, 95% CI 0.77–6.88) reflects this: the CI crosses 1.0, meaning the TCGA-alone estimate is not statistically significant. The manuscript correctly notes this.
- **MSK cohort:** n=117, 38 events (~32%). This cohort is advanced-disease-enriched (predominantly PDTC and ATC). Pooling TCGA (primary PTC) with MSK (advanced thyroid cancer) via meta-analysis introduces population heterogeneity that the DerSimonian-Laird method is designed to handle, but with only 2 studies, I^2=0% is trivially achieved even when true heterogeneity exists (low power of Cochran Q with k=2). The I^2=0% claim should not be interpreted as confirming homogeneity; it simply means the two point estimates are numerically close.
- **Reference group direction (CRITICAL):** The manuscript does not explicitly state in either the v2 draft or the STAR methods which group is the reference in the Cox model (DM2 as reference with DM1 as test, or vice versa). The AUDIT_LOCKED_RESULTS.md flags this as unresolved (rows S-01 through S-03). If HR=2.53 means DM1 patients have 2.53× the hazard of DM2 patients, the direction is correct and consistent with biological interpretation (DM1 = poor differentiation = adverse). This must be confirmed from the lifelines model summary output before submission.
- **Multivariable adjustment:** STAR methods (line 204) mentions "multivariable models included age, stage and selected driver covariates where sample size allowed." The v2 figure legend (Fig 5b) references a multivariate Cox forest. However, the main narrative (v2 Results §2) reports only the unadjusted estimate in the pooled meta-analysis. It is unclear whether the pooled HR 2.53 is adjusted or unadjusted. This must be clarified.

**Verdict:** The OS meta-analysis is structurally sound for a two-cohort preliminary survival analysis but the direction claim requires source-code verification, and I^2=0% should not be presented as proof of homogeneity.

---

### 1.2 Progression-free interval and BRAF interaction (S-04, S-05)

**Design:** Cox PH with continuous standardized DM1 score as predictor; endpoint = PFI per Liu 2018; interaction tested as multiplicative term (score × BRAF).

**Issues identified:**

- **Reference direction (CRITICAL — flagged in AUDIT_LOCKED_RESULTS.md S-04):** The HR per +1 SD DM1 score in BRAF+ patients is 0.66 (95% CI 0.48–0.92). Since DM1 = iodine-handling-LOW state, a higher DM1 score could mean MORE lineage loss (i.e., closer to DM1 pole) OR could represent the z-standardized value within the BRAF+ cohort where scoring direction is arbitrary. If HR=0.66 means a higher DM1 score (more lineage loss) is associated with LOWER PFI hazard (better prognosis), this is biologically counterintuitive and could indicate a sign error in score direction or reference coding. Alternatively, if the DM1 score is coded as the iodine-handling component (higher = more intact), then HR=0.66 per +1 SD means more intact lineage = lower hazard, which IS biologically plausible. The manuscript must resolve this from the source code.

- **Interaction term (S-05):** The interaction HR=0.47 (p=0.022) indicates that the DM1 score × BRAF interaction modifies the DM1 effect on PFI. With 49 total PFI events in n=472, and only 33 events in the BRAF+ subset (n=287), the interaction is estimated from a small number of events. Standard guidance (e.g., Peduzzi 1995) suggests at minimum 10 events per variable in the full interaction model; this model includes score, BRAF, and score×BRAF, requiring ~30 events for adequate estimation — marginally met. The result should be clearly labeled as hypothesis-generating.

- **Multiplicity:** The manuscript tests DM1 interaction with BRAF and with RAS separately (p=0.022 and p=0.41, respectively). Two interaction tests on the same dataset without multiplicity correction should be noted. Even at unadjusted α=0.05, the BRAF interaction p=0.022 is borderline. Bonferroni correction would require p<0.025 to remain significant — barely satisfied.

- **Subgroup vs interaction (CRITICAL):** v2 Results §3 (line 78) describes stratified analysis within BRAF+ as a separate Cox model, not as a subgroup of the interaction model. The stratified HR=0.66 is from a separate model, not from the interaction model output. These are two different quantities. Both are valid, but they must be clearly distinguished and not conflated.

- **IHC 3-plex vs Full-8 comparison (IHC-01 vs C-S04):** The v2 abstract states the Full-8 panel HR = 0.66 (log-rank p=7.6e-3) and the IHC 3-plex HR = 0.65 (log-rank p=1.7e-4). The IHC 3-plex has a smaller Cohen's d but a more significant log-rank p. This is mathematically possible if the 3-gene mean methylation proxy is a better single-variable separator than the 8-gene score for PFI in BRAF+. However, comparison of multiple panel combinations on the same dataset inflates type I error. The manuscript does acknowledge this is exploratory but the figure structure (14 combinations tested) requires explicit multiplicity warning.

---

### 1.3 ARI computation

**Question:** Are the two partitions compared correctly?

The manuscript reports ARI values for several comparisons:
- Pan-genome top-5000 partition vs 8-gene-derived DM1/DM2 labels: ARI = 0.92
- TIERA67 partition vs 8-gene labels: ARI = 0.90
- Driver-anchor partition vs 8-gene labels: ARI = -0.007
- 8-gene panel alone vs itself (circular): should be 1.0; manuscript reports 0.49 (this likely means 8-gene KMeans partition vs some other reference, not self-comparison)

**Issue:** The 8-gene panel ARI = 0.49 (04_results.md line 26) is compared against "the panel-driven DM1/DM2 partition." If the 8-gene panel IS the reference partition, comparing it to itself yields ARI=1.0 by definition. ARI=0.49 likely means the 8-gene-alone clustering was compared to the pan-genome-derived labels, not the reverse. This needs clarification in the methods.

---

### 1.4 Fisher OR for fusion enrichment

**Claim:** DM1 76.8% fusion-positive vs DM2 30.9%; OR = 7.41 (95% CI 4.38-12.55, p=1.9e-13).

**Structural question:** The manuscript never states the 2×2 contingency table explicitly. Given DM1 fusion+ = 63/82, DM1 fusion- = 19/82. If DM2 fusion+ rate is 30.9%, we need the DM2 denominator (SV-tested DM2 tumors). From n=504 total, 28.4% DM1 = ~143 DM1 tumors; ~361 DM2 tumors. If 97.3% of 557 primary tumors have SV data, and 82 SV-tested DM1 tumors are available, then DM2 SV-tested = approximately 542-82=460. DM2 fusion+ at 30.9% of 460 = ~142. This would give a 2×2 of: DM1 fusion+ 63, DM1 fusion- 19; DM2 fusion+ 142, DM2 fusion- 318. The OR from this table: (63×318)/(19×142) = 20034/2698 = 7.42. This matches the reported OR of 7.41, suggesting the numbers are internally consistent, but the actual denominators should be stated explicitly in the manuscript.

**Missingness:** Chi-square p=0.56 for MAR is correctly presented, but the claim "consistent with missing-at-random" should be softened per CLAIM_LANGUAGE_MATRIX.md. Chi-square p=0.56 rules out strong MNAR but does not confirm MAR.

---

### 1.5 Methylation Cohen's d direction

**Claim:** DM1 has higher mean beta than DM2 (0.385 vs 0.253); per-gene Cohen's d values are all positive for 7 of 8 genes.

**Direction verification:** DM1 = iodine-handling-low = expected to have higher promoter methylation (silencing). Mean beta DM1 > DM2 is biologically consistent. Cohen's d = (DM1_mean - DM2_mean)/pooled_SD would be positive (DM1 higher). This is internally consistent and biologically correct. 

**Arithmetic (M-11):** (0.385-0.253)/0.253 = 52.2% increase. The manuscript claims "52% higher" — accurate.

**SLC5A5 exception:** d=0.22 (p=0.42, NS) is correctly presented as non-significant and attributed to post-translational/enhancer regulation of NIS. This is a responsible negative finding.

---

### 1.6 Meta-analysis: fixed vs random effects

**Issue:** The v2 draft text (line 72) says "Random-effects meta-analysis yielded a pooled hazard ratio of 2.53." The STAR methods (line 204) correctly specify DerSimonian-Laird random-effects. Figure 5 caption (line 254) says "random-effects overall-survival meta-analysis." However, with only two studies (k=2), DerSimonian-Laird becomes numerically equivalent to the weighted fixed-effects estimate when I^2=0%, because the between-study variance tau^2 is estimated as zero. The manuscript should note that with k=2, the distinction between fixed and random effects is minimal and the pooled estimate should be interpreted with caution.

**I^2=0%:** Not meaningful with k=2 (Cochran Q statistic with df=1 has only ~30% power to detect true I^2=50%). The claim of "no detectable between-cohort heterogeneity" is accurate but should not imply homogeneity confirmation.

---

### 1.7 External validation statistical summary

**"80/80 cells in expected direction" (V-05):** Definition of "cell" is undefined in the v2 text. The auditor cannot verify whether this is 8 genes × 10 cohorts = 80, or 10 genes × 8 cohorts, or some other decomposition. The specific cohorts contributing to the 80×80 matrix must be listed.

**Lee 2024 n=370 Cohen's d=0.24:** This is a surprisingly small effect size (d=0.24) compared to the overall validation claim of mean d=2.81. The p=2.7e-7 significance despite d=0.24 is explained by the large sample size (n=370). The manuscript correctly reports both; readers should be aware these are different analyses using different comparison variables (8-gene score vs subB dedifferentiation score).

---

### 1.8 Power simulation

**Monte Carlo design:** 500 simulations per n; exponential hazard; TCGA BRAF+ priors (HR=1.49, event rate 11.7%, 70% completeness). Empirical power = fraction reaching p<0.05 in Cox regression.

**Issues:**
- HR=1.49 is the point estimate from TCGA. The 95% CI for the BRAF+ HR (0.66 with CI 0.48-0.92 per +1 SD; converting: HR per +1 SD of 0.66 in one direction implies HR in the other direction) is needed to set sensitivity bounds. The manuscript does report sensitivity analyses at HR=1.40 and HR=1.22, which is appropriate.
- The simulation uses TCGA BRAF+ priors for a prospective cohort. If the prospective cohort has a different BRAF+ prevalence, event rate, or follow-up distribution, the power estimate may not apply.
- The power of >90% at n=200 is correctly labeled as "support for feasibility" with the appropriate caveat about TCGA priors.

---

## Summary table

| Analysis | Critical issue | Severity |
|---|---|---|
| OS pooled meta-analysis | Reference direction unverified from source | FATAL |
| PFI BRAF+ HR=0.66 | Direction counterintuitive — requires source code | FATAL |
| DM1 × BRAF interaction | Single-cohort; 33 events; multiplicity not corrected | MAJOR |
| I^2=0% with k=2 | Statistically misleading to interpret as homogeneity | MAJOR |
| ARI=0.49 description | Comparison groups not clearly defined | MINOR |
| Fisher OR (fusion) | 2×2 table not shown in manuscript | MINOR |
| 80/80 cell matrix | "Cell" definition absent | MAJOR |
| Lee 2024 n=370 vs n=632 | Two separate analyses not distinguished | MAJOR |
| MAR language | Chi-square p=0.56 ≠ MAR confirmed | MINOR |
| Power simulation | TCGA priors may not generalize; properly caveated | MINOR |
