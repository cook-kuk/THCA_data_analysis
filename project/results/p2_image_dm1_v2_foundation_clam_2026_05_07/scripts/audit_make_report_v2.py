#!/usr/bin/env python3
"""Final audit report — phase 1 + 2 + 3 combined."""
from __future__ import annotations
import json
from pathlib import Path

import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
AUD = ROOT / "analysis_supp" / "audit_ras_auc100"


def main():
    p1 = json.loads((AUD / "AUDIT_SUMMARY.json").read_text())
    p2 = json.loads((AUD / "PHASE2_AUDIT_SUMMARY.json").read_text())
    rs = (json.loads((AUD / "RETRAIN_SUMMARY.json").read_text())
          if (AUD / "RETRAIN_SUMMARY.json").exists() else None)
    p3 = (json.loads((AUD / "PHASE3_AUDIT_SUMMARY.json").read_text())
          if (AUD / "PHASE3_AUDIT_SUMMARY.json").exists() else None)
    p4 = (json.loads((AUD / "phase4" / "PHASE4_SOLUTIONS_SUMMARY.json").read_text())
          if (AUD / "phase4" / "PHASE4_SOLUTIONS_SUMMARY.json").exists() else None)

    md = []
    md.append("# Paper 2 image-DM1 — RAS-like AUC=1.000 audit (v2 final)")
    md.append("")
    md.append("## Verdict")
    md.append("")
    md.append("**RAS-like AUC=1.000 is fragile AND the broader image-DM1 "
              "claim has multiple paper-blocking confounds.**")
    md.append("")
    md.append("| Issue | Severity | Evidence |")
    md.append("|---|:---:|---|")
    md.append("| RAS_like ∩ FVPTC = 16/16 (same slides, not independent) | "
              "**MAJOR** | C1 |")
    md.append("| RAS-like 1.000 separation gap = 0.017 prob | "
              "**MAJOR** | C2 |")
    md.append("| Multi-seed RAS-like AUC range [0.44, 0.92] | "
              "**MAJOR** | A multi-seed |")
    md.append("| CPU rerun seed=42 gives 0.79 not 1.00 (irreproducible) | "
              "**MAJOR** | A seed=42 |")
    md.append("| Male AUC = 0.35 (worse than random) | "
              "**PAPER-KILLING** | E1 |")
    md.append("| All 3 RAS-like DM1 cases from 2 TSS centers (DJ, FK); "
              "all 10 TSS=EM cases are DM2 | **PAPER-KILLING** | E2 |")
    md.append("| 0 off-diagonal (cPTC, RAS_like) cases — cannot "
              "disentangle histology from molecular subtype | "
              "**PAPER-KILLING** | E4 |")
    md.append("| Clinical-only LR (histology + sex + TSS) AUC=0.768 "
              "**> CLAM 0.746** | **PAPER-KILLING** | E5 |")
    md.append("")

    # --- A. data structure facts ---
    md.append("## A. Data structure (no retraining)")
    md.append("")
    md.append("- 59 slides = 59 unique cases (no patient-level leakage).")
    md.append("- Histology × molecular subtype crosstab (in our 59):")
    md.append("")
    md.append("| | RAS_like | BRAF_like | unknown |")
    md.append("|---|---:|---:|---:|")
    md.append("| FVPTC | 16 | 0 | 0 |")
    md.append("| cPTC | 0 | 41 | 0 |")
    md.append("| unknown | 0 | 0 | 2 |")
    md.append("")
    md.append("→ **Histology and molecular subtype are 100% co-linear in "
              "our subset.** RAS_like ∩ FVPTC = 16/16. The two subgroup "
              "AUC=1.000 rows are one finding labelled twice.")
    md.append("")

    # --- B. RAS-like 16-slide raw distribution
    raw = pd.read_csv(AUD / "c2_RAS_like_16_predictions.tsv", sep="\t")
    raw_sorted = raw.sort_values("prob_DM1")
    md.append("## B. Raw OOF predictions for the 16 RAS-like slides")
    md.append("")
    md.append("| rank | submitter_id | TSS | label | prob_DM1 | fold |")
    md.append("|---:|---|:---:|:---:|---:|---:|")
    for i, row in enumerate(raw_sorted.itertuples(), 1):
        tss = row.submitter_id.split("-")[1]
        md.append(f"| {i} | {row.submitter_id} | {tss} | "
                  f"{'**DM1**' if row.label==1 else 'DM2'} | "
                  f"{row.prob_DM1:.3f} | {int(row.fold)} |")
    md.append("")
    md.append("**Separation gap = 0.017** (lowest DM1 0.549 − highest DM2 "
              "0.532). 1 misrank → AUC = 0.974. Bootstrap CI [1.000, "
              "1.000] is deterministic resampling, not generalization.")
    md.append("")

    # --- C. Statistical (perm + bootstrap)
    md.append("## C. Permutation null + bootstrap CI")
    md.append("")
    md.append("| group | n (pos/neg) | observed AUC | perm p (2-sided) | "
              "boot 95% CI |")
    md.append("|---|---:|---:|---:|---|")
    for k, v in p1["C3_permutation"].items():
        md.append(f"| {k} | {v['n']} ({v['n_pos']}/{v['n_neg']}) | "
                  f"{v['observed_auc']:.3f} | "
                  f"{v['permutation_p_two_sided']:.4f} | "
                  f"[{v['bootstrap_ci95_lo']:.3f}, "
                  f"{v['bootstrap_ci95_hi']:.3f}] |")
    md.append("")
    md.append("Combinatorial baseline P(perfect | random)=1/C(16,3)="
              "0.0018. Perm p=0.003 means model ranking is barely above "
              "the random-shuffle null (1 in 333).")
    md.append("")

    # --- D. Split-stress retrain
    md.append("## D. Split-stress retrain")
    md.append("")
    if rs is not None:
        md.append("| split strategy | overall AUC | RAS-like AUC | "
                  "BRAF-like AUC |")
        md.append("|---|---:|---:|---:|")
        md.append("| original (GPU, KFold seed=42) | 0.746 | **1.000** | "
                  "0.592 |")
        for r in rs.get("strategy_A", []):
            md.append(f"| KFold seed={r['seed']} (CPU) | "
                      f"{r['overall_auc']:.3f} | {r['ras_like_auc']:.3f} | "
                      f"{r['braf_like_auc']:.3f} |")
        b = rs.get("strategy_B", {})
        if b:
            md.append(f"| **StratifiedKFold(3)** | "
                      f"{b['overall_auc']:.3f} | "
                      f"**{b['ras_like_auc']:.3f}** | "
                      f"{b['braf_like_auc']:.3f} |")
        c = rs.get("strategy_C", {})
        if c:
            md.append(f"| **LOO 16 RAS-like** | n/a | "
                      f"**{c['loo_auc']:.3f}** | n/a |")
        md.append("")
        md.append("- Original 1.000 does not reproduce on CPU (seed=42 "
                  "gives 0.795). Multi-seed range [0.44, 0.92].")
        md.append("- Stratified-CV / LOO converge to **0.92** — the "
                  "honest estimate.")
        md.append("")

    # --- E. Phase 2 confounds (paper-killers)
    md.append("## E. Phase-2 confounders — paper-killing findings")
    md.append("")

    # E1
    md.append("### E1. Sex-stratified AUC")
    md.append("")
    md.append("| sex | n (pos/neg) | AUC | 95% boot CI |")
    md.append("|---|---:|---:|---|")
    for r in p2["E1_sex"]:
        md.append(f"| {r['group']} | {r['n']} ({r['n_pos']}/{r['n_neg']}) | "
                  f"{r['auc']:.3f} | [{r['ci_lo']:.3f}, "
                  f"{r['ci_hi']:.3f}] |")
    md.append("")
    md.append("→ Model is **worse than random in n=13 Males**. CI lower "
              "bound 0.000. Sex-stratified review will flag this.")
    md.append("")

    # E2
    md.append("### E2. TSS distribution among 16 RAS-like slides")
    md.append("")
    md.append("| TSS | n | n DM1 | mean prob_DM1 |")
    md.append("|---|---:|---:|---:|")
    for r in p2["E2_tss"]:
        md.append(f"| {r['tss']} | {int(r['n'])} | {int(r['n_pos'])} | "
                  f"{r['mean_prob']:.3f} |")
    md.append("")
    md.append("→ **All 3 RAS-like DM1 cases concentrate in TSS=DJ (×2) + "
              "FK (×1).** All 10 RAS-like cases from TSS=EM are DM2. "
              "The model could be classifying **scanning center / "
              "staining batch**, not biology.")
    md.append("")

    # E5
    e5 = p2["E5_clinical_baseline"]
    md.append("### E5. Clinical-covariate baseline beats CLAM")
    md.append("")
    md.append("| Predictor | AUC |")
    md.append("|---|---:|")
    md.append("| Histology only (FVPTC indicator) | 0.680 |")
    md.append("| Histology + sex + TSS dummies (5-fold OOF LR) | "
              f"**{e5['clinical_only_oof_auc']:.3f}** |")
    md.append(f"| CLAM-UNI overall (5-fold OOF) | 0.746 |")
    md.append("")
    md.append(f"→ **CLAM image gain over clinical = "
              f"{e5['clam_image_gain']:+.3f}** — image model **does "
              "not add over a 3-feature logistic regression** on "
              "histology + sex + TSS.")
    md.append("")

    # --- F. Phase 3 decisive retraining
    md.append("## F. Phase-3 decisive retraining")
    md.append("")
    if p3 is not None:
        # E6
        md.append("### E6. Init-only variance (StratifiedKFold split fixed)")
        md.append("")
        md.append("| init seed | overall AUC | RAS-like AUC | "
                  "BRAF-like AUC |")
        md.append("|---:|---:|---:|---:|")
        for r in p3["E6_init_variance"]:
            md.append(f"| {r['init_seed']} | {r['overall_auc']:.3f} | "
                      f"{r['ras_like_auc']:.3f} | "
                      f"{r['braf_like_auc']:.3f} |")
        ras_vals = [r["ras_like_auc"] for r in p3["E6_init_variance"]]
        md.append("")
        md.append(f"RAS-like across 6 init seeds at fixed split: median "
                  f"{sorted(ras_vals)[len(ras_vals)//2]:.3f}, range "
                  f"[{min(ras_vals):.3f}, {max(ras_vals):.3f}]")
        md.append("")

        # E7
        md.append("### E7. Label-shuffle null (3 reps)")
        md.append("")
        md.append("| rep | overall AUC | RAS-like AUC |")
        md.append("|---:|---:|---:|")
        for r in p3["E7_label_shuffle_null"]:
            md.append(f"| {r['rep']} | {r['overall_shuf_auc']:.3f} | "
                      f"{r['ras_like_shuf_auc']:.3f} |")
        md.append("")
        md.append("→ Shuffled-label AUC ≈ 0.5 = no memorization shortcut. "
                  "Healthy null (model isn't fitting noise structure).")
        md.append("")

        # E8
        e8 = p3["E8_leave_one_tss_out"]
        md.append("### E8. Leave-one-TSS-out — TSS batch-effect probe")
        md.append("")
        md.append("| TSS held out | n test | n DM1 | LOTO AUC |")
        md.append("|---|---:|---:|---:|")
        for r in e8["per_tss"]:
            auc_str = (f"{r['tss_holdout_auc']:.3f}"
                       if r["tss_holdout_auc"] is not None
                       and not pd.isna(r["tss_holdout_auc"])
                       else "n/a")
            md.append(f"| {r['tss_held_out']} | {r['n_test']} | "
                      f"{r['n_pos']} | {auc_str} |")
        md.append("")
        md.append(f"- Pooled OOF overall AUC under LOTO = "
                  f"**{e8['pooled_overall_auc']:.3f}**")
        md.append(f"- Pooled OOF RAS-like AUC under LOTO = "
                  f"**{e8['pooled_ras_like_auc']:.3f}**")
        md.append("")
        md.append("→ If pooled LOTO RAS-like AUC drops sharply vs the "
                  "0.92 from random / stratified split, the AUC was "
                  "TSS-batch-driven. If it stays ≥ 0.85, the signal "
                  "survives center-holdout.")
        md.append("")

    # --- Phase 4 SOLUTIONS attempted
    if p4 is not None:
        md.append("## G. Phase-4 — solutions attempted")
        md.append("")
        md.append("Tried 4 fixes for the TSS / sex / image-gain failures:")
        md.append("")
        # S1
        md.append("### S1. TSS-ComBat on UNI features "
                  "(per-tile mean centering by TSS)")
        md.append("")
        md.append("| init seed | overall | RAS-like | BRAF-like | Male | "
                  "Female |")
        md.append("|---:|---:|---:|---:|---:|---:|")
        for r in p4["S1_tss_combat"]:
            md.append(f"| {r['init_seed']} | {r['overall']:.3f} | "
                      f"{r['ras']:.3f} | {r['braf']:.3f} | "
                      f"{r['male']:.3f} | {r['female']:.3f} |")
        ras1 = [r["ras"] for r in p4["S1_tss_combat"]]
        male1 = [r["male"] for r in p4["S1_tss_combat"]]
        md.append("")
        md.append(f"- RAS-like median {sorted(ras1)[len(ras1)//2]:.3f} "
                  f"(was 0.92 stratified, 1.000 original) — TSS effect "
                  "partially removed → AUC drops to honest range.")
        md.append(f"- Male median {sorted(male1)[len(male1)//2]:.3f} "
                  f"(was 0.35 raw) — sex asymmetry partially fixed.")
        md.append(f"- Overall median ~0.66 (was 0.74) — image gain "
                  "smaller after batch correction.")
        md.append("")

        # S2
        md.append("### S2. TSS-balanced StratifiedKFold on (label × TSS-bucket)")
        md.append("")
        md.append("| init seed | overall | RAS-like | BRAF-like |")
        md.append("|---:|---:|---:|---:|")
        for r in p4["S2_tss_balanced_split"]:
            md.append(f"| {r['init_seed']} | {r['overall']:.3f} | "
                      f"{r['ras']:.3f} | {r['braf']:.3f} |")
        ras2 = [r["ras"] for r in p4["S2_tss_balanced_split"]]
        md.append("")
        md.append(f"- RAS-like median "
                  f"{sorted(ras2)[len(ras2)//2]:.3f} — "
                  "explicit TSS stratification preserves the apparent "
                  "RAS-like advantage. Combined with the LOTO 0.308 "
                  "result, this confirms: RAS-like AUC stays high "
                  "**only when training and test contain the same TSS "
                  "centers**.")
        md.append("")

        # S3
        md.append("### S3. Multimodal LR (CLAM_prob + clinical features)")
        md.append("")
        md.append("| Predictor | OOF AUC |")
        md.append("|---|---:|")
        for r in p4["S3_multimodal_LR"]:
            md.append(f"| {r['predictor']} | {r['auc']:.3f} |")
        md.append("")
        clin = next(r for r in p4["S3_multimodal_LR"]
                    if r["predictor"] == "clinical_only")
        multi = next(r for r in p4["S3_multimodal_LR"]
                     if r["predictor"] == "multimodal")
        md.append(f"→ **Image channel net gain over clinical = "
                  f"{multi['auc'] - clin['auc']:+.3f}**. Adding CLAM "
                  "to clinical features does not improve AUC; "
                  "the image model is **redundant with clinical metadata**.")
        md.append("")

        # S5
        md.append("### S5. EM-out subsample (drop 10 RAS-like-EM-DM2 anchor)")
        md.append("")
        md.append("| init seed | overall | RAS-like | BRAF-like |")
        md.append("|---:|---:|---:|---:|")
        for r in p4["S5_em_out"]:
            md.append(f"| {r['init_seed']} | {r['overall']:.3f} | "
                      f"{r['ras']:.3f} | {r['braf']:.3f} |")
        ras5 = [r["ras"] for r in p4["S5_em_out"]]
        md.append("")
        md.append(f"- After dropping the 10 (TSS=EM, RAS_like, all DM2) "
                  f"slides: RAS-like median = "
                  f"{sorted(ras5)[len(ras5)//2]:.3f} "
                  "(was 0.92). Signal degrades meaningfully without "
                  "the EM-DM2 anchor — confirms that the apparent "
                  "RAS-like advantage was driven by the EM-DM2 "
                  "majority.")
        md.append("")

    # --- H. Bottom line
    md.append("## H. Bottom line for Paper 2")
    md.append("")
    md.append("**Cannot ship 'image-DM1 AUC=0.83 (UNI)' as the headline "
              "as currently framed.** Three independent confounds:")
    md.append("")
    md.append("1. **Sex asymmetry** — Female 0.83 vs Male 0.35. "
              "Reviewers will ask for sex-balanced training.")
    md.append("2. **TSS batch confound** — molecular subtype is "
              "co-linear with scanning center. Need leave-one-TSS-out "
              "as primary reporting metric, OR ComBat-equivalent "
              "feature normalization.")
    md.append("3. **No image gain over clinical** — histology + sex + "
              "TSS dummies already give AUC=0.768. CLAM at 0.746 is "
              "*below* this baseline; the image channel is not "
              "the load-bearing signal in TCGA.")
    md.append("")
    md.append("**Path forward options:**")
    md.append("")
    md.append("- **Defer Paper 2 launch** until Korean K2 H&E (or "
              "FFPE multi-site) cohort is in. K2 will have different "
              "TSS structure → if RAS-like AUC=0.92 holds on K2 with "
              "TSS-balanced split, the signal is real.")
    md.append("- **Reframe to clinical-augmented model**: report "
              "(histology + sex + TSS + CLAM-prob) as the predictor; "
              "compare against (clinical only) baseline; main claim "
              "becomes \"CLAM adds X% on top of clinical\". Currently "
              "X = 0%.")
    md.append("- **Restrict to TSS-balanced subsample** for the main "
              "result; report TSS-stratified analysis as primary.")
    md.append("")
    md.append("**Do NOT publish RAS-like AUC=1.000 as currently framed.**")
    md.append("")
    md.append("---")
    md.append("")
    md.append("## I. Final synthesis — what each fix bought us")
    md.append("")
    md.append("| Fix | What it solves | What it doesn't |")
    md.append("|---|---|---|")
    md.append("| S1 ComBat | TSS center batch effect; sex asymmetry "
              "partially | n_pos=3 RAS-like instability remains |")
    md.append("| S2 TSS-balanced split | confirms TSS-balanced AUC ≈ 0.91 | "
              "RAS-like AUC still seed-dependent [0.80, 0.92]; "
              "**LOTO 0.308 still the truth** |")
    md.append("| S3 multimodal LR | quantifies image gain | image gain "
              "= NEGATIVE; clinical alone (0.768) > clinical + "
              "image (0.756) |")
    md.append("| S5 EM-out | confirms EM-DM2 anchor is load-bearing | "
              "RAS-like AUC drops to 0.72; signal partially survives |")
    md.append("")
    md.append("**No fix recovers the 1.000.** The honest TCGA-only "
              "ceiling is overall AUC ≈ 0.65–0.80, RAS-like AUC ≈ "
              "0.70–0.90 (highly seed-dependent), and the image "
              "channel adds ZERO over (histology + sex + TSS) "
              "logistic regression.")
    md.append("")
    md.append("**The fundamental fix is more data, not more analysis.**")
    md.append("")
    md.append("- TCGA-THCA full WSI = ~500 slides (10× current 59)")
    md.append("- K2 cohort (planned) = ~320 slides with different TSS "
              "structure → independent test of any signal that survives "
              "TCGA training")
    md.append("- FFPE multi-site (planned) = 600 slides → adequate "
              "power for stratified analysis")
    md.append("")
    md.append("Until then: **defer Paper 2 main-text image-DM1 claim**, "
              "or reframe as \"clinical-only + image as auxiliary\" "
              "with image gain reported truthfully (0% on TCGA).")
    md.append("")

    out = AUD / "AUDIT_REPORT_v2.md"
    out.write_text("\n".join(md))
    print(f"WROTE {out}")


if __name__ == "__main__":
    main()
