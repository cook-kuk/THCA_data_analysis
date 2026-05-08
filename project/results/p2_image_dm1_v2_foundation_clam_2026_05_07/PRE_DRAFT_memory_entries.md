# Memory pre-drafts — 결과 받으면 즉시 저장 후보

## Path A — PASS (Phase 1 ≥50% slides + Phase 2 AUC > 0.70)

```yaml
filename: v19_paper2_image_dm1_v2_PASS.md
type: project
title: v19 Paper 2 image-DM1 v2 PASS
description: 2026-05-07 foundation-model (UNI) + CLAM v2 architecture achieved Paper 2 H&E→DM1 launch threshold; closure ResNet50 baseline overcome.
```

Body:
```
2026-05-07: Paper 1 audit page hardening 후 user 결정으로 image-DM1 NO-GO closure re-opened with foundation-model + CLAM v2.

Key results:
- Phase 1 (GSE250521 UNI on H&E tiles): PC1↔DM1_like |ρ|>0.3 in [N]% of slides
- Phase 2 (TCGA WSI 60-90 slides; CLAM gated-attention MIL): pooled 5-fold CV AUC = [X.XXX]
- closure ResNet50 baseline ~0.55 → exceeded by [Δ]

Why: closure 2026-05-04 was ResNet50 (ImageNet generic) + tile-mean pooling. v2 = UNI ViT-Huge (Mass-100K WSI pretrained) + CLAM attention MIL = architecturally distinct.

How to apply: Paper 2 image-DM1 first-pass figure set is now unlocked. Bundang outreach 결과 받으면 Korean external validation 추가. Paper 1 spatial supp 도 UNI embedding overlay 로 augment 가능.

Reference: project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/SPRINT_RESULT.md
```

---

## Path B — MARGINAL (Phase 1 25-50% OR Phase 2 0.60-0.70)

```yaml
filename: v19_paper2_image_dm1_v2_MARGINAL.md
type: project
title: v19 Paper 2 image-DM1 v2 MARGINAL
description: 2026-05-07 foundation-model + CLAM v2 partial signal; Paper 1 supp augment only; Paper 2 launch deferred.
```

Body:
```
2026-05-07: image-DM1 v2 sprint produced partial signal (Phase 1 verdict={MARGINAL}; Phase 2 AUC=[X.XX]).

Decision: Paper 1 spatial supp augment with UNI embedding overlay (descriptive). Paper 2 launch deferred until cross-cohort or larger n confirmation.

Why partial: foundation-model captures lineage signal more than ResNet50 but CLAM-AUC < 0.70 threshold suggests bag-level prediction has residual heterogeneity (DM1 sub-A/sub-B mixed).

Next steps: re-run Phase 2 with stratified sub-A/sub-B labels (after master DM1 universe re-derivation, D3 manual decision).
```

---

## Path C — FAIL (Phase 1 <25% OR Phase 2 ≤ 0.60)

```yaml
filename: v19_paper2_image_dm1_v2_NEGATIVE.md
type: project
title: v19 Paper 2 image-DM1 v2 NEGATIVE
description: 2026-05-07 foundation-model + CLAM v2 also failed; closure NO-GO re-confirmed at v2 architecture.
```

Body:
```
2026-05-07: closure NO-GO re-opened with foundation-model + CLAM. Result: still NO-GO at v2 architecture (Phase 1 |ρ|>0.3 in [N]% < 25%; Phase 2 AUC = [X.XX] ≤ 0.60).

Implication: image-DM1 prediction does not work even with state-of-the-art pathology foundation model + attention MIL. DM1 is *not* H&E-inferable at WSI resolution.

Use as honest negative for Paper 1 supp methods: "we attempted with v1 (ResNet50, closure 2026-05-04) and v2 (UNI + CLAM, this sprint); both failed → DM1 is a transcriptomic-only axis, not pathology-inferable."

Strong reviewer-trust building. Add to Paper 1 supp methods + Limitations.

Bundang Korean FFPE 시점에 한 번 더 시도 (사람 환자 cohort 다른 distribution 가능; 마지막 chance).
```

---

## Path D — Phase 2 NOT RUN (Phase 1 만 PASS, Phase 2 skip 또는 Pod 시간 부족)

```yaml
filename: v19_paper2_image_dm1_v2_PHASE1_ONLY.md
type: project
title: v19 Paper 2 image-DM1 v2 Phase 1 only — Phase 2 deferred
```

Body:
```
2026-05-07: Phase 1 (GSE250521 UNI spatial signal) PASS [PC1↔DM1 |ρ|>0.3 in {N}% slides]. Phase 2 (TCGA WSI CLAM) deferred (resource / pod timeout).

Decision: Paper 1 spatial supp augment with UNI embedding figure. Phase 2 retry when L40S/A100 host availability returns.

Resume command: bash project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/runpod_dispatch/run_phase2.sh
```
