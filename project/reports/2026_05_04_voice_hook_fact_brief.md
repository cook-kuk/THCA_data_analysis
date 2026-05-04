# Paper 1 Hook ¶1 — Voice Fact Brief (본인 키보드용)

**Date:** 2026-05-04 · **Owner:** Seungho Cook
**Mode:** Marathon scaffolding — fact-only briefing for user's voice-protected Hook ¶1.
**This file does NOT contain a Hook draft.** No prose, no sentences, no rhetorical phrasing. Bullets only.
**Sources verified:** OUTLINE v8 §Intro/§Abstract; PAPER1_MARATHON_AUDIT_BUNDLE 2026-05-04; molecular-only LOCK §2/§4/§5; HM closure report; R1–R5 prose `2026-05-03`.

---

## 1. Clinical problem (3 bullets)

- DTC (differentiated thyroid carcinoma) shows 15–35% recurrence; ATA 2015 risk stratification spans 3–5% (low risk) to 50–75% (high risk) with anatomic/pathologic features only — no transcriptional axis.
- Within ATA tiers, no molecular axis predicts dedifferentiation trajectory or RAI-uptake collapse.
- Korean / Asian-specific autoimmune-thyroid context (Hashimoto's overlap, DPB1\*05:01) is not represented in TCGA-derived frameworks.

## 2. Why BRAF/RAS is exhausted (3 bullets)

- BRAF transcript V600E vs WT (TCGA-THCA n = 273 vs 182): Cohen's d = −0.044, MW p = 0.567 — mutation status does not propagate to BRAF mRNA level.
- Single-feature DM-discrimination AUC: BRAF 0.602, TERT 0.578, KRAS 0.525, NRAS 0.521, HRAS 0.500 — all below 0.7 biomarker threshold.
- Driver_anchor 12-gene KMeans k=2 clustering vs DM1/DM2: ARI = −0.007 (random). Drivers alone cannot define the DM transcriptional axis. (TCGA Cancer Network 2014 BRS dichotomy + BRAF V600E + TERT promoter mutation = "marginal benefit" gap.)

## 3. What dark-matter gap remains (3 bullets)

- BRAF-neg ∩ RAS-neg subset n = 156 in TCGA-THCA has no driver-based risk stratification (this paper's molecular axis target).
- DM2 sub-population mechanism (less-differentiated cluster) was not mechanistically explained before this work.
- Existing landscapes — TCGA 2014 BRS dichotomy, Yoo 2016 16-gene panel, Krishnamoorthy / Landa 2016 PDTC/ATC genomic hallmarks, Pan/Ge 2025 ATC proteogenomic subtypes — none address autoimmune-PTC overlap as a distinct mechanism.

## 4. What this paper resolves (5 bullets)

- 8-gene RAI panel (SLC5A5/TPO/TG/TSHR/PAX8/NKX2-1/FOXE1/DIO1) defines DM1/DM2 transcriptional axis on TCGA-THCA n = 500 (DM1 = 140, DM2 = 360); centered-profile classifier 5-fold CV AUC = 0.962. Cluster definition robust to candidate-pool: TIERA67 ARI = 0.903 ≈ pan-genome top-5000 MAD ARI = 0.918; Driver_anchor alone ARI = −0.007.
- Pan-Asian HLA cohort n = 874 (K2 PRJEB11591 235 + Lee 2024 GSE213647 630 + GSE286332 PTC arm 9): DPB1\*05:01 carrier 53.2% Korean PTC vs 31.3% Han Chinese ctrl (Chu 2018) vs 44.0% Han Chinese GD; Korean-vs-ctrl OR = 2.50 (p = 4×10⁻²⁶); Cochran's Q = 1.18, I² = 0% across 3 Korean sub-cohorts.
- GSE286332 PTC vs PTC+HT (n = 18: 9 + 9): 10,380 DEGs (BH-FDR < 0.05); 8-gene RAI Cohen d = −1.60 (PAX8 d = −2.32, NKX2-1 d = −1.92, FOXE1 d = −1.75, SLC5A5 d = +0.25 preserved); HLA-II module d = +3.65 (p = 4×10⁻⁴); GSEA Hallmark IFN-γ NES = +1.80 FDR = 2×10⁻⁴, KEGG Type I diabetes NES = +1.92 FDR = 0.
- TCGA Hashimoto-like signature (top 150 up + 50 down GSE286332 DEGs, padj < 0.01, |LFC| > 1) classifies 18–30% of TCGA n = 500 (GMM 18%, Otsu 19.6%, top-30% = 30%); DM2-enriched OR = 0.20, p = 6.4×10⁻¹⁰; confounder-residualized OR = 0.29, p = 8×10⁻⁹; Korean GSE213647 replication 22.8% (Otsu 28.2%).
- DM1 sub-B (n = 56) = 53/56 (94.6%) mutation-negative NBNR cluster; sub-A (n = 84) = 51/74 mut-tested RAS+ (69%; 51/84 = 61% of total) classical FVPTC core. Korean GSE213647 sub-B-like rate 47–53% (GMM 47.2%, Otsu 52.5%) — bridges Yu professor K2 NBNR observation.
- Survival (TCGA): PFI dichotomized HR = 2.04 [1.15–3.61] p = 0.015; DFI multivariate HR = 1.41 p = 0.025 (BRAF/RAS/age/stage adjusted).

## 5. What NOT to overclaim (5 bullets)

- Verbs: use **supports** / **consistent with** / **stratifies orthogonally to**. Forbidden: `proves`, `establishes`, `demonstrates the framework`, `definitively`.
- "All risks resolved" / "모든 risk 해소" is false. Honest open negatives remain: dark-matter subset N1 (BRAF-neg ∩ RAS-neg n = 156) PFI HR = 1.20 p = 0.80 (NS); spot-level TROP2 × DM1_resid ρ = −0.016 (Q3 NEG); H&E → DM1 closure battery NO-GO (r = 0.022, AUROC = 0.511). These belong in Limitations / honest-negative supp, NOT in Hook.
- No image / pathology / morphology language for DM1 inference. Closure battery 2026-05-04 = full 폐기 (verdict D). Forbidden in Hook: `H&E-inferable`, `pathology-AI triage`, `WSI-validated`, `tile-level DM1 inference`, `morphology-derived DM1`, `H&E predicts DM1`.
- TROP2 framing must stay **tumor-level**: "tumor-population vulnerability" — NOT "DM1-high spots are TROP2-high" (Q3 spot-level NEG). Functional validation (TROP2 IHC, sacituzumab IC50, organoid) "remains required" — must appear in Limitations.
- No venue-probability percentages (`Cancer Cell 도전권 20–30%`, `Nat Cancer 50–65% reach`); no Bundang FFPE numbers (outreach stage, no data); no `RunPod G1/G2/G3` / `image-DM1 Phase B/C` references in any active-claim form.

## 6. Candidate opening angles (label-only, no prose)

- **A. clinical unmet need angle** — ATA 2015 anatomic-only stratification + 15–35% recurrence range + no transcriptional axis frame.
- **B. molecular dark matter angle** — BRAF/RAS-negative n = 156 subset with no driver-based risk anchor + DM2 unexplained mechanism frame.
- **C. RAI-lineage silencing angle** — TF backbone collapse (FOXE1, NKX2-1, PAX8) + iodide-uptake transcriptional axis frame, anchored in canonical RAI biology.
- **D. fusion/epigenetic actionability angle** — TF collapse → DNMT silencing → STAT3/AP-1 activation → TROP2 re-expression 4-step framework, "consistent with" actionability hypothesis (tumor-level only).

(Selection is user's voice decision. No combination prescribed.)

## 7. Exact numbers SAFE to use in Hook

| domain | number | source |
|---|---|---|
| DTC recurrence range | 15–35% | OUTLINE Abstract Background |
| ATA 2015 recurrence span | 3–5% (low) → 50–75% (high) | OUTLINE Intro ¶1, ATA cheatsheet |
| TCGA-THCA cohort | n = 500 (DM1 = 140, DM2 = 360) | Methods M1 |
| BRAF V600E vs WT transcript | Cohen d = −0.044, p = 0.567, n = 273 vs 182 | R3, audit table #3 |
| Driver single-feature AUC ceiling | BRAF 0.602, TERT 0.578 | R3 |
| Driver_anchor cluster ARI | −0.007 | R3, R4 |
| Pan-Asian HLA pooled cohort | n = 874 | R1 |
| Chu 2018 GD vs ctrl | n = 1,468 vs 1,490 (total 2,958) | M1 |
| DPB1\*05:01 carrier | 53.2% Korean / 31.3% ctrl / 44.0% GD | R1 |
| DPB1\*05:01 Korean-vs-ctrl OR | 2.50 (95% CI 2.07–3.02), p = 4×10⁻²⁶ | R1 |
| GSE286332 cohort | n = 18 (9 PTC + 9 PTC+HT) | R2 |
| GSE286332 DEGs | 10,380 (BH-FDR < 0.05) | R2 |
| 8-gene RAI Cohen d (PTC+HT) | −1.60 (MW p = 0.008) | R2 |
| HLA-II Cohen d (PTC+HT) | +3.65 (p = 4×10⁻⁴) | R2 |
| GSEA IFN-γ | NES = +1.80, FDR = 2×10⁻⁴ | R2 |
| TIERA67 vs pan-genome ARI | 0.903 vs 0.918 | R4 |
| TCGA Hashimoto-like rate | 18–30% (GMM 18%, Otsu 19.6%, top-30% = 30%) | R5b |
| TCGA Hashimoto-like × DM2 | OR = 0.20, p = 6.4×10⁻¹⁰ | R5b |
| DM1 sub-B mutation-negative | 53/56 (94.6%) | R5d, post-HM-closure value |
| DM1 sub-A RAS+ | 51/74 mut-tested (69%); 51/84 = 61% of total | R5d, post-HM-closure phrasing |
| Korean GSE213647 sub-B-like | 47.2% GMM / 52.5% Otsu | R5d |
| TCGA PFI HR (dichotomized) | 2.04 [1.15–3.61], p = 0.015 | LOCK §2 #3 |
| TCGA DFI multivariate HR | 1.41, p = 0.025 (BRAF/RAS/age/stage adjusted) | LOCK §2 #3 |
| TF collapse (sample-level) | FOXE1 −1.21 / NKX2-1 −0.65 / STAT3 +1.55 / FOSL1 +0.98 / DNMT1 +0.74 | LOCK §2 #4 |

## 8. Exact numbers NOT safe in Hook

| number | reason | proper home |
|---|---|---|
| H&E LOSO Spearman r = 0.022, AUROC = 0.511 | image-DM1 closure NO-GO; Hook is positive frame | (omit; Methods supp negative-feasibility paragraph only) |
| Spot-level TROP2 × DM1_resid ρ = −0.016 | Q3 NEG; cannot be a Hook claim | Limitations / honest-negative supp |
| Dark-matter subset N1 PFI HR = 1.20, p = 0.80 (NS) | honest negative; underpowered | Limitations / supp |
| "96% mutation-negative" (sub-B) | superseded 2026-05-04; correct = 53/56 (94.6%) | use 94.6% only |
| "DM1-high spots are TROP2-high" | Q3 NEG | not citeable |
| Cancer Cell 20–30% / Nat Cancer 50–65% (venue probability) | speculative internal planning | not submission-facing |
| Bundang FFPE n > 50 | outreach stage, no data delivered | omit |
| RunPod G1/G2/G3 / image-DM1 Phase B/C | superseded by closure NO-GO | omit |
| TROP2 Δmean +4.33, FDR = 1.2×10⁻³² (tumor-level) | safe in Results / Discussion as **tumor-level**; **NOT** in Hook spot-level form | Results panel + tumor-level framing |
| Foundation-model expected gain (+0.05–0.15 Spearman) | derived from closure battery; not Hook material | Methods supp |
| 8-gene single-cluster ARI 0.49 | Pillar 4 honest framing; "modest, clinically interpretable" — Hook overload risk | R4 / Discussion |
| Bootstrap CI [−0.997, −0.916] (12-slide ST) | Supp-only spatial detail | Supp |
| Moran's I 0.36 (28 ST slides) | Supp spatial structure | Supp |

---

## 9. Recent multimodal AI citation candidate (Intro / Discussion context, NOT Hook material)

**Reference:** Valanarasu JMJ, Xu H, Usuyama N, ..., Poon H. *Multimodal AI generates virtual population for tumor microenvironment modeling.* **Cell** 189(2), Jan 2026. **DOI: 10.1016/j.cell.2025.11.016**. (BibTeX `@Valanarasu2026` in `project/manuscript_v8/03_intro_references.bib`.)

What it does (one bullet each):
- **GigaTIME** = cross-modal H&E → virtual multiplex-IF translator (image-to-image, NOT image-to-RNA)
- Trained on **40 × 10⁶ paired cells**, 21 proteins
- Applied to 14,256 patients × 51 Providence hospitals × 24 cancers × 306 subtypes → 299,376 virtual mIF datasets
- Independent **10,200 TCGA validation**
- 1,234 protein-biomarker-staging-survival associations recovered

Why it matters for THIS Paper 1 (citation logic, NOT Hook prose):
- **Strengthens our negative result**: GigaTIME = SOTA proof that H&E carries information when paired with molecular ground-truth at scale; our negative result on H&E-alone → DM1_resid prediction (closure battery 2026-05-04) is therefore not a "wrong-model" failure — it is a different task class (continuous molecular regression) at radically different scale (3.2K tiles vs 40M cells).
- **Bounds the image-DM1 angle**: future H&E-DM1 work would require RNA-paired training at GigaTIME-class scale, not ImageNet-frozen tile embedding. Already encoded in `2026_05_04_paper1_supp_methods_he_negative_feasibility.md` (STAR Methods supp draft, commit `cf14656`).
- **Frames the molecular axis**: our 8-gene RNA axis IS the molecular ground-truth that GigaTIME-class systems would need for thyroid; the paper establishes the axis on RNA, the multimodal-AI bridge can be future work.

Where author may cite (suggestions only — voice-protected):
- **Introduction ¶ on multimodal AI frontier** (1 sentence): "Recent multimodal AI for the tumor microenvironment (Valanarasu et al., 2026) translates H&E into virtual multiplex-IF at population scale, but H&E-only molecular subtype inference for thyroid cancer dedifferentiation has not been established."
- **Discussion ¶ on multimodal AI bridge as future work** (1–2 sentences): "GigaTIME-class multimodal AI (Valanarasu et al., 2026) demonstrates the population-scale feasibility of H&E ↔ molecular-layer translation given sufficient paired training; our 8-gene RNA molecular ground-truth supplies one such target axis for future H&E-paired training in thyroid."
- **Limitations ¶** (1 sentence, optional): "We did not pursue H&E-based triage of the molecular subtype here (closure battery 2026-05-04 NO-GO); RNA-paired multimodal training at the scale of Valanarasu et al. (2026) would be required to revisit this angle."
- **STAR Methods H&E negative-feasibility supp** (already integrated in `cf14656` follow-up — commit imminent)

What NOT to do with this citation in Hook:
- ❌ Do NOT use as "image-DM1 is feasible" claim — GigaTIME does protein translation, not transcriptional axis regression
- ❌ Do NOT use as "we did the same thing GigaTIME does" — task class and scale are different
- ❌ Do NOT cite in any Hook prose claiming H&E-DM1 success (closure battery NO-GO standing)

---

## Voice-protected reminder (do NOT auto-write)

- Hook ¶1 = 본인 키보드 strict. Default Claude 권한 = scaffolding/infra만 (`v17_sprint_vs_marathon_violation` 메모리).
- "고고" / "다 해줘" / "faster" 명령은 voice-protected sprint generate 권한 NOT 부여.
- Suggested user-decision toggles (per OUTLINE ¶1 line 37): "no mechanistic compass" vs "without molecular guidance" — user's pick.
- Landa 2016 cite save 적용은 Discussion §3.1 책임. Hook ¶1에서 Krishnamoorthy 2025 인용 시 사실 4개 모두 misattribution 위험 → Landa 2016 JCI 126(3):1052–1066 사용 (`v17_landa2016_cite_save` 메모리).

---

Hook fact brief ready. User writes the paragraph.
