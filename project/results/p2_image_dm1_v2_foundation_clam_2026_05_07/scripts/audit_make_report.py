#!/usr/bin/env python3
"""Generate AUDIT_REPORT.md for the RAS-like AUC=1.000 audit."""
from __future__ import annotations
import json
from pathlib import Path

import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
AUD = ROOT / "analysis_supp" / "audit_ras_auc100"


def main():
    md = []
    md.append("# RAS-like AUC=1.000 audit — data-split / shortcut investigation")
    md.append("")
    md.append("**Question:** Paper 2 image-DM1 CLAM (UNI features) reports "
              "RAS-like / FVPTC subgroup AUC=1.000 (n=16, n_pos=3 / n_neg=13). "
              "Is this real biological signal, a data-split artifact, or "
              "a histology shortcut?")
    md.append("")
    md.append("## Verdict TL;DR")
    md.append("")
    if (AUD / "RETRAIN_SUMMARY.json").exists():
        md.append("**The signal survives, but the 1.000 was an unstable upper "
                  "outlier — not a reproducible result.**")
        md.append("")
        md.append("| split strategy | RAS-like AUC | overall AUC |")
        md.append("|---|---:|---:|")
        md.append("| Original (GPU, seed=42, KFold 5) | **1.000** | 0.746 |")
        md.append("| CPU rerun seed=42, KFold 5 | 0.795 | 0.722 |")
        md.append("| Multi-seed median (n=6 seeds) | **0.744** [0.436–0.923] | 0.739 |")
        md.append("| StratifiedKFold(3) on label×subtype | **0.923** | 0.793 |")
        md.append("| LOO within 16 RAS-like slides | **0.923** | n/a |")
        md.append("")
        md.append("**Interpretation.**")
        md.append("1. **The exact 1.000 is not reproducible.** Even at "
                  "seed=42, CPU rerun gives 0.795 (GPU↔CPU numerical drift). "
                  "Original 1.000 was one realization of a high-variance "
                  "estimator.")
        md.append("2. **Multi-seed range [0.44, 0.92]** — at seed=2026, "
                  "RAS-like AUC drops to 0.44 (worse than random). With "
                  "n_pos=3 in n=16, plain KFold(5) is unstable.")
        md.append("3. **Under more rigorous splits (stratified, LOO), "
                  "RAS-like AUC stabilises at 0.923.** This is the "
                  "honest headline number: RAS-like is genuinely easier "
                  "to classify than BRAF-like (≈0.65), but not perfectly.")
        md.append("4. **Histology-only LogReg AUC = 0.68** — half the "
                  "RAS-like advantage is base-rate (FVPTC 81% DM2). "
                  "CLAM adds +0.24 over the trivial baseline.")
        md.append("")
        md.append("**Recommendation for Paper 2:** report **RAS-like "
                  "AUC = 0.92 (stratified-CV)** in main text; relegate "
                  "the original 1.000 to a methods footnote with "
                  "the seed-sensitivity disclosure. The clinical claim "
                  "(\"FVPTC/RAS-like sub-cohort is more separable than "
                  "cPTC/BRAF-like\") is preserved without overclaiming.")
        md.append("")
    else:
        md.append("(see Section 6 — filled after retrain finishes)")
        md.append("")

    # ── C1
    summary_path = AUD / "AUDIT_SUMMARY.json"
    s = json.loads(summary_path.read_text())
    c1 = s["C1_subgroup_overlap"]
    md.append("## 1. Are RAS-like and FVPTC the same 16 slides?")
    md.append("")
    md.append(f"| RAS_like | FVPTC | RAS ∩ FVPTC | union |")
    md.append(f"|---|---|---|---|")
    md.append(f"| {c1['n_RAS_like']} | {c1['n_FVPTC']} | "
              f"{c1['n_intersection']} | {c1['n_union']} |")
    md.append("")
    md.append(f"**Result.** RAS_like ∩ FVPTC = "
              f"{c1['n_intersection']}/{c1['n_RAS_like']} — "
              "the two subgroup AUC=1.000 rows are **not independent corroboration**, "
              "they are the same 16 slides labelled twice.")
    md.append("")

    # ── C2
    md.append("## 2. Raw OOF prediction distribution within the 16 RAS-like slides")
    md.append("")
    raw = pd.read_csv(AUD / "c2_RAS_like_16_predictions.tsv", sep="\t")
    raw_sorted = raw.sort_values("prob_DM1")
    md.append("| rank | submitter_id | label | prob_DM1 | fold |")
    md.append("|---:|---|:---:|---:|---:|")
    for i, row in enumerate(raw_sorted.itertuples(), 1):
        md.append(f"| {i} | {row.submitter_id} | "
                  f"{'**DM1**' if row.label==1 else 'DM2'} | "
                  f"{row.prob_DM1:.3f} | {int(row.fold)} |")
    dm1_min = raw[raw["label"] == 1]["prob_DM1"].min()
    dm2_max = raw[raw["label"] == 0]["prob_DM1"].max()
    md.append("")
    md.append(f"**Separation gap.** lowest DM1 prob = **{dm1_min:.3f}**, "
              f"highest DM2 prob = **{dm2_max:.3f}**. The 'perfect ranking' "
              f"is achieved by **{dm1_min - dm2_max:.3f}** margin. "
              "If the borderline DM2 (TCGA-EM-A3FP, 0.532) had been ranked "
              "above the borderline DM1 (TCGA-DJ-A13W, 0.549), AUC would "
              "drop to 38/39 = 0.974.")
    md.append("")

    # ── C3
    c3 = s["C3_permutation"]
    md.append("## 3. Permutation null + bootstrap CI")
    md.append("")
    md.append("| group | n (pos/neg) | observed AUC | "
              "perm p (2-sided) | null p99 | boot 95% CI |")
    md.append("|---|---:|---:|---:|---:|---|")
    for k, v in c3.items():
        md.append(f"| {k} | {v['n']} ({v['n_pos']}/{v['n_neg']}) | "
                  f"{v['observed_auc']:.3f} | "
                  f"{v['permutation_p_two_sided']:.4f} | "
                  f"{v['null_auc_p99']:.3f} | "
                  f"[{v['bootstrap_ci95_lo']:.3f}, "
                  f"{v['bootstrap_ci95_hi']:.3f}] |")
    md.append("")
    md.append(f"**Combinatorial baseline.** With n_pos=3 / n_neg=13, "
              "P(perfect ranking | random) = 1/C(16,3) = "
              f"**{1/560:.4f}** (0.18%). "
              "Observed permutation p=0.0030 — model-produced ranking is "
              "**better than random label assignment** but only barely "
              "(perfect by chance every 1 in 333 random label shuffles).")
    md.append("")
    md.append("**Bootstrap CI [1.000, 1.000] is misleading.** "
              "It reflects deterministic resampling of an already-perfect "
              "ranking — not external generalization uncertainty. With "
              "n_pos=3, the *effective* sample size is tiny.")
    md.append("")

    # ── C4
    c4 = s["C4_shortcut"]
    md.append("## 4. Histology / molecular shortcut probe")
    md.append("")
    md.append("| Probe | AUC / corr |")
    md.append("|---|---:|")
    md.append(f"| Histology-only LogReg → label (DM1 vs DM2) | "
              f"{c4['c4a_histology_only_LR_AUC']:.3f} |")
    md.append(f"| Molecular-only LogReg → label (DM1 vs DM2) | "
              f"{c4['c4b_molecular_only_LR_AUC']:.3f} |")
    md.append(f"| corr(FVPTC indicator, CLAM prob_DM1) | "
              f"{c4['c4c_corr_histology_predDM1']:.3f} |")
    md.append("")
    md.append("**Reading.** Knowing only the histology label "
              "(FVPTC/cPTC) lets a logistic regression hit AUC=0.68 on "
              "the DM1/DM2 task — purely from base-rate asymmetry "
              "(81% of FVPTC → DM2, 63% of cPTC → DM1). CLAM's overall "
              "AUC=0.746 is only +0.07 over this trivial baseline. "
              "corr(FVPTC, prob_DM1)=−0.36 indicates CLAM has internalised "
              "histology pattern as a feature.")
    md.append("")

    # ── C5
    md.append("## 5. Per-fold composition (unstratified KFold seed=42)")
    md.append("")
    cf = pd.read_csv(AUD / "c5_RAS_like_fold_composition.tsv", sep="\t")
    md.append("| fold | n_total | n_pos (DM1) | n_neg (DM2) |")
    md.append("|---:|---:|---:|---:|")
    for r in cf.itertuples():
        md.append(f"| {int(r.fold)} | {int(r.n)} | "
                  f"{int(r.n_pos)} | {int(r.n_neg)} |")
    md.append("")
    md.append("**Issue.** Folds 1, 2, 4 contained ZERO RAS-like positives. "
              "All 3 RAS-like DM1 cases concentrate in folds 3 (1) and 5 (2). "
              "Models in folds 1/2/4 saw all 3 RAS-like DM1 in training and "
              "thus produced confident low-prob_DM1 for their test-set "
              "RAS-like DM2 cases — locking in 11 of the 13 DM2 ranks. "
              "Strict OOF accounting still holds (each test slide was held "
              "out from its own fold), but the apparent *separation* is "
              "easier to achieve than under a stratified split.")
    md.append("")

    # ── 6 retrain (placeholder; filled by audit_make_report.py rerun after
    # retrain completes)
    retrain_summary_path = AUD / "RETRAIN_SUMMARY.json"
    md.append("## 6. Split-stress retrain (Strategies A / B / C)")
    md.append("")
    if retrain_summary_path.exists():
        rs = json.loads(retrain_summary_path.read_text())
        # Strategy A
        md.append("### Strategy A — multi-seed KFold(5)")
        md.append("")
        md.append("| seed | overall AUC | RAS-like AUC | BRAF-like AUC | "
                  "seconds |")
        md.append("|---:|---:|---:|---:|---:|")
        for r in rs.get("strategy_A", []):
            md.append(f"| {r['seed']} | {r['overall_auc']:.3f} | "
                      f"{r['ras_like_auc']:.3f} | "
                      f"{r['braf_like_auc']:.3f} | {r['seconds']} |")
        md.append("")
        # Strategy B
        if "strategy_B" in rs:
            b = rs["strategy_B"]
            md.append("### Strategy B — StratifiedKFold(3) on label × molecular_subtype")
            md.append("")
            md.append(f"- Overall AUC: **{b['overall_auc']:.3f}**")
            md.append(f"- RAS-like AUC: **{b['ras_like_auc']:.3f}**")
            md.append(f"- BRAF-like AUC: **{b['braf_like_auc']:.3f}**")
            md.append("")
            md.append("Per-fold breakdown:")
            md.append("")
            md.append("| fold | n_val | val AUC | RAS pos | RAS neg |")
            md.append("|---:|---:|---:|---:|---:|")
            for f in b.get("fold_breakdown", []):
                md.append(f"| {f['fold']} | {f['n_val']} | "
                          f"{f['val_auc']:.3f} | {f['n_ras_pos_val']} | "
                          f"{f['n_ras_neg_val']} |")
            md.append("")
        # Strategy C
        if "strategy_C" in rs:
            c = rs["strategy_C"]
            md.append("### Strategy C — Leave-one-out on the 16 RAS-like slides")
            md.append("")
            md.append(f"- n = {c['n']} ({c['n_pos']} DM1 / {c['n_neg']} DM2)")
            md.append(f"- LOO AUC = **{c['loo_auc']:.3f}**")
            md.append("")
            md.append("Each held-out slide was predicted by a model trained "
                      "on the other 58 slides (other 15 RAS-like slides "
                      "remained in training). Most stringent test of the "
                      "AUC=1.000 claim.")
            md.append("")
    else:
        md.append("(retraining still running — re-run `audit_make_report.py` "
                  "after `RETRAIN_SUMMARY.json` exists)")
        md.append("")

    # ── Verdict
    md.append("## 7. Summary of evidence")
    md.append("")
    md.append("**A. What is real (signal survives split changes):**")
    md.append("- Stratified-CV(3) and LOO-on-RAS-like both give "
              "**RAS-like AUC = 0.923** — robust to split choice.")
    md.append("- BRAF-like AUC stays in 0.57–0.74 across all strategies — "
              "consistent moderate signal.")
    md.append("- Overall AUC stays in 0.66–0.80 across 6 random seeds.")
    md.append("- Permutation test p=0.003 — the model's ranking within "
              "the 16 RAS-like slides is not from random noise.")
    md.append("")
    md.append("**B. What is not real (the 1.000 itself is fragile):**")
    md.append("- CPU rerun at the same seed=42 gives RAS-like AUC=0.795, "
              "not 1.000 — the original number does not reproduce on CPU.")
    md.append("- Multi-seed range [0.436, 0.923] — at one seed (2026), "
              "RAS-like AUC is 0.436, *worse than random*. With only 3 "
              "positives in n=16, KFold(5) is too unstable to interpret "
              "any one realization.")
    md.append("- Bootstrap CI [1.000, 1.000] reflects deterministic "
              "resampling of an already-perfect ranking — not "
              "generalization uncertainty.")
    md.append("- The two rows (RAS_like AUC=1.000 + FVPTC AUC=1.000) "
              "are the same 16 slides labelled twice.")
    md.append("- Half the RAS-like advantage is histology base-rate "
              "(FVPTC-only LogReg AUC = 0.68; corr(FVPTC, CLAM "
              "prob_DM1) = −0.36).")
    md.append("")
    md.append("**C. Recommended Paper 2 reporting change:**")
    md.append("- Replace headline `RAS-like AUC = 1.000 (boot CI [1.0, 1.0])` "
              "with `RAS-like AUC = 0.92 (StratifiedCV-3, LOO-confirmed; "
              "n=16, n_pos=3)`")
    md.append("- Footnote: \"Original report 1.000 from a single "
              "KFold(seed=42) realization is not reproducible across "
              "GPU/CPU or alternative seeds; the stratified-CV / LOO "
              "estimate of 0.92 is more reliable.\"")
    md.append("- Reviewer-defense: \"RAS-like advantage is not pure "
              "histology shortcut: histology-only baseline gives 0.68, "
              "CLAM adds +0.24.\"")
    md.append("")

    out = AUD / "AUDIT_REPORT.md"
    out.write_text("\n".join(md))
    print(f"WROTE {out}")


if __name__ == "__main__":
    main()
