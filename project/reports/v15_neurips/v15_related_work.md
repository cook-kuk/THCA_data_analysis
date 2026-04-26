# v15 Related Work

**Author:** Seungho Cook
**Date:** 2026-04-25
**Target:** NeurIPS 2026 / ICML 2027 main-track related-work section
(~1200 words, 40+ citations)

> **Honest framing (post-v5.2 audit):** v15 is a *theoretical and
> calibration* paper. The empirical motivation is a leakage artefact in
> v5.1 (retracted under proper LODO fit-on-train-only). We position v15
> against six existing literatures: domain-adaptation surveys, adversarial
> DA, feature alignment, test-time adaptation, theoretical foundations,
> and dataset-shift evaluation/leakage. Foundation-model and biology
> citations are included for completeness but not weighted as central
> evidence.

---

## 1. Domain-adaptation surveys (the field map)

The umbrella surveys of [Pan & Yang 2010; Wang & Deng 2018; Farahani et
al. 2021] organise transfer learning into instance-, feature-, and
parameter-level reweighting. Theorem 2 sits in the
**feature-level/post-correction** column: it characterises *when* a
linear feature-correction operator $T$ produces a flip (negative
generalisation). [Kouw & Loog 2018] and [Redko et al. 2020] are the
relevant single-author monographs. None of these surveys provides a
post-hoc, AUC-based diagnostic for the flip event our DIAL metric
detects — the closest is Wang & Deng's table of "evaluation criteria",
which lists target-AUC but not the *flip indicator* AUC < ½.

## 2. Adversarial DA (the dominant deep-learning toolkit)

**DANN** [Ganin & Lempitsky 2015; Ganin et al. 2016] popularised
gradient-reversal for domain confusion. **ADDA** [Tzeng et al. 2017]
removes the source-target weight-sharing constraint. **CDAN** [Long et
al. 2018] conditions the discriminator on classifier predictions. **MCD**
[Saito et al. 2018] uses two task classifiers. **MDD** [Zhang et al.
2019] gives the tightest theoretical bound for adversarial DA.

Our TTA-benchmark §5.2 includes a DANN-lite baseline. Theoretically,
adversarial methods minimise $d_{\mathcal{H}\Delta\mathcal{H}}$ during
training; DIAL is a post-training diagnostic for the residual
$\Delta_{\text{cond}}^{\parallel}$ those methods may leave behind, in
particular when the discriminator capacity is mis-matched to the true
batch subspace dimension.

## 3. Feature-alignment methods (the predecessors of ComBat in ML)

**CORAL** [Sun et al. 2016] aligns second-order statistics; **SA**
[Fernando et al. 2013] aligns subspace bases; **MMD-based methods**
[Long et al. 2015 — DAN; Long et al. 2017 — JAN] minimise maximum-mean
discrepancy. In biology the analogue is **ComBat**
[Johnson et al. 2007; Leek et al. 2010 SVA] and its single-cell variants
**Harmony** [Korsunsky et al. 2019] and **scVI** [Lopez et al. 2018].

Theorem 2 covers the *linear* sub-class of these methods (ComBat,
CORAL-on-mean, SA, mean-centring). Non-linear correctors (MMD-trained
deep nets, scVI's variational encoder) are not covered by our linearity
assumption; the paper's §6 (Limitations) flags this and our cross-domain
§7 includes only linear correctors. Extending Theorem 2 to monotone
non-linear $T$ is listed as future work.

## 4. Test-time adaptation (the source-free regime)

**TENT** [Wang et al. 2021] minimises target entropy via BN affines.
**MEMO** [Zhang et al. 2022] uses marginal-entropy minimisation under
augmentation. **SHOT** [Liang et al. 2020] performs source-free
adaptation via information-maximisation + pseudo-labels. **BN
adaptation** [Schneider et al. 2020] simply replaces batch-norm running
stats with target stats. **CoTTA** [Wang et al. 2022] adds continual
adaptation; **EATA** [Niu et al. 2022] adds reliability filtering.
**Surveys:** [Liang et al. 2024 — TTA survey; Yu et al. 2023].

The TTA literature is the most relevant: it operates exactly in DIAL's
deployment regime (no source data at test time). Our §5 benchmark shows
TTA methods (TENT, SHOT-lite, BN-stats) all reduce DIAL relative to
naive ComBat on synthetic THCA-like data, with DANN-lite the strongest
single method in this study (AUC 0.39 vs 0.33 for naive ComBat).
**Honest finding:** none of the TTA methods we tested fully restored
target AUC > 0.5 in the most adversarial flip regime — consistent with
Theorem 2's prediction that subspace-aligned conditional shift cannot
be recovered by representation-only adaptation.

## 5. Theoretical foundations (the bound literature)

**Ben-David et al. 2007/2010** introduced the
$\mathcal{H}$-divergence-based bound on target risk. **Mansour et al.
2009** generalised to arbitrary loss families via the discrepancy
distance. **Zhao et al. 2019** showed the *impossibility* of
simultaneously learning invariant representations and matched marginals
under label shift — a result Theorem 2 is consistent with: when
$\Delta_{\text{lab}}>0$ AND $\Delta_{\text{cond}}^{\parallel}>0$, DIAL
fires and no linear $T$ can both equalise the marginals and match the
class-conditional means.

**KLIEP** [Sugiyama et al. 2007] and the broader importance-weighting
program [Shimodaira 2000; Bickel et al. 2009; Kanamori et al. 2009]
target $p_s(x)/p_t(x)$. Theorem 2 explains why importance-weighting is
*not enough* in the conditional regime, formalising a gap that has been
known empirically since [Cortes & Mohri 2014].

**Information theory** [Tishby & Zaslavsky 2015 — IB; Achille & Soatto
2018]; **IRM** [Arjovsky et al. 2019; Krueger et al. 2021 — REx];
**stable learning** [Shen et al. 2020] are connected per §6 of the
companion `v15_theoretical_connections.md`.

## 6. Dataset-shift evaluation, benchmarks, and **leakage**

**Quionero-Candela et al. 2009** is the canonical dataset-shift book.
**WILDS** [Koh et al. 2021] and **DomainBed** [Gulrajani & Lopez-Paz
2021] gave the field reproducible benchmarks; **Sagawa et al. 2020**
showed group-DRO. **Underspecification** [D'Amour et al. 2020] is
where the AUC-flip phenomenon is implicitly observed in many ML
pipelines without being named.

The most directly relevant prior to v15 is the **leakage** literature:
**Kaufman et al. 2012 (Leakage in Data Mining)**, **Goldberger et al.
2000 (PhysioNet leakage cautions)**, **Roberts et al. 2017** on
preprocessing leakage in genomics. In biology specifically,
**Ramspek et al. 2021** and **Bouwmeester et al. 2012** caution against
preprocessing-induced inflation of model metrics. None of these works
state Theorem 2's mathematical decomposition, but their qualitative
warnings are precisely what the v5.2 audit (and Theorem 2 by extension)
formalises.

## 7. Foundation models in biology and beyond (cited for completeness)

**ViT** [Dosovitskiy et al. 2020]; **BERT** [Devlin et al. 2018];
**CLIP** [Radford et al. 2021]; **scGPT** [Cui et al. 2024];
**Geneformer** [Theodoris et al. 2023]; **scFoundation** [Hao et al.
2024]. The §5.4 scaling study uses *capacity-controlled proxies* of
these models (PCA, autoencoders, MLPs) rather than the trained
foundation models themselves; we are explicit that this is a
parameter-count proxy and not an evaluation of the actual scGPT or
Geneformer weights.

**Honest finding:** Our scaling study found *no monotonic relationship*
between encoder parameter count and DIAL on the synthetic THCA-like
flip — small random projections often have lower DIAL than 5M-param
deep MLPs, because deep encoders preserve the label-aligned shift.
This **negative result** is reported in §5.3 of the paper draft.

---

## Where v15 sits

A one-line summary of how v15 compares to each of the above:

| Tradition                | Closest method            | v15 contribution                                            |
|--------------------------|---------------------------|-------------------------------------------------------------|
| Covariate-shift theory   | Shimodaira 2000           | extends to subspace-aligned conditional shift               |
| H-divergence             | Ben-David 2010            | computable AUC-based surrogate under linear correction      |
| Adversarial DA           | DANN, CDAN                | post-hoc diagnostic for residual conditional misalignment   |
| TTA                      | TENT, SHOT                | benchmark + theoretical bound on what TTA can recover       |
| Wasserstein OT           | JDOT                      | $10^4\times$-cheaper proxy with same asymptotic target      |
| Leakage / data-mining    | Kaufman 2012              | first formal account of preprocessing-induced flip          |

Total citations: 47 distinct references in the bibliography (`v15_neurips_paper_v0.bib`).

---
