#!/usr/bin/env python3
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from neoimmune_common import (
    ensure_run_dir,
    metric_binary,
    norm_series,
    patient_topn_metrics,
    peptide_features,
    precision_recall_at_k,
    safe_read_table,
    simple_train_predict,
    write_md,
    write_tsv,
)


CLEAN_MODEL_TOKENS = [
    "Structure_LR",
    "Wave8",
    "TCR",
    "SelfSim",
    "ESM2_Bayesian",
    "GP_quantum",
    "VQC",
    "W7A",
    "W7B",
    "quantum_kernel",
    "Stack_mean",
    "Stack_median",
    "Stack_LR",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    label_col = "label_immunogenicity"
    canon[label_col] = pd.to_numeric(canon[label_col], errors="coerce")
    local = safe_read_table(outdir / "local_scores" / "all_local_scores.tsv.gz")
    if local.empty:
        local_wide = pd.DataFrame({"candidate_id": canon["candidate_id"]})
    else:
        clean = local[local["local_model_name"].astype(str).str.contains("|".join(CLEAN_MODEL_TOKENS), case=False, na=False)].copy()
        clean["normalized_score"] = pd.to_numeric(clean["normalized_score"], errors="coerce")
        local_wide = clean.pivot_table(index="candidate_id", columns="local_model_name", values="normalized_score", aggfunc="mean").reset_index()
    d = canon.merge(local_wide, on="candidate_id", how="left")
    feat = peptide_features(d)
    score_cols = [c for c in local_wide.columns if c != "candidate_id"]
    for c in score_cols:
        feat[c] = pd.to_numeric(d[c], errors="coerce")
    use = d[label_col].notna()
    d.loc[use, "clean_science_score"] = simple_train_predict(pd.concat([d.loc[use, ["candidate_id", label_col, "source_study", "patient_id"]], feat.loc[use]], axis=1), list(feat.columns), label_col, "source_study")
    if (~use).any():
        d.loc[~use, "clean_science_score"] = norm_series(feat.loc[~use].mean(axis=1), True).fillna(0.5)
    d["clean_science_score"] = pd.to_numeric(d["clean_science_score"], errors="coerce").fillna(0.5).clip(0, 1)
    flags_path = outdir / "metrics" / "leakage_candidate_flags.tsv"
    if flags_path.exists():
        flags = safe_read_table(flags_path)[["candidate_id", "any_existing_overlap_flag"]]
        d = d.merge(flags, on="candidate_id", how="left")
        d["any_existing_overlap_flag"] = d["any_existing_overlap_flag"].fillna(False).astype(bool)
    else:
        d["any_existing_overlap_flag"] = False

    pred_cols = ["candidate_id", "patient_id", "dataset_source", "peptide_mut", "hla_allele", label_col, "clean_science_score"]
    write_tsv(d[pred_cols].sort_values("clean_science_score", ascending=False), outdir / "clean_track" / "clean_science_predictions.tsv")

    rows = []
    eval_df = d[d[label_col].notna()].copy()
    if not eval_df.empty:
        base = metric_binary(eval_df[label_col], eval_df["clean_science_score"])
        for k in [10, 20, 50]:
            base.update(precision_recall_at_k(eval_df, "clean_science_score", label_col, k))
        base.update(patient_topn_metrics(eval_df, "clean_science_score", label_col, 20))
        base.update({"track": "clean_science_track", "model": "CLEAN-Neo++ clean local ranker", "split": "source_grouped_cv_or_existing_labels", "features": "peptide/WT/HLA/tumor context/local-only scores"})
        rows.append(base)
        for c in score_cols:
            tmp = metric_binary(eval_df[label_col], eval_df[c])
            for k in [10, 20]:
                tmp.update(precision_recall_at_k(eval_df, c, label_col, k))
            tmp.update(patient_topn_metrics(eval_df, c, label_col, 20))
            tmp.update({"track": "clean_science_track", "model": c, "split": "existing_local_artifact", "features": "single existing local score"})
            rows.append(tmp)
    lb = pd.DataFrame(rows)
    ordered = ["track", "model", "split", "features", "n", "positives", "AUROC", "AUPRC", "Precision@10", "Precision@20", "Recall@10", "Recall@20", "Recall@50", "patient_hit_rate@20", "patient_recall@20", "patients_evaluated", "calibration_brier"]
    for c in ordered:
        if c not in lb.columns:
            lb[c] = np.nan
    lb = lb[ordered + [c for c in lb.columns if c not in ordered]].sort_values(["AUPRC", "patient_hit_rate@20"], ascending=False)
    write_tsv(lb, outdir / "metrics" / "leaderboard_clean_track.tsv")

    strict_rows = []
    strict_df = eval_df[~eval_df["any_existing_overlap_flag"]].copy() if "any_existing_overlap_flag" in eval_df.columns else eval_df.copy()
    if not strict_df.empty:
        base = metric_binary(strict_df[label_col], strict_df["clean_science_score"])
        for k in [10, 20, 50]:
            base.update(precision_recall_at_k(strict_df, "clean_science_score", label_col, k))
        base.update(patient_topn_metrics(strict_df, "clean_science_score", label_col, 20))
        base.update({"track": "clean_science_track", "model": "CLEAN-Neo++ clean local ranker", "split": "strict_no_existing_overlap", "features": "peptide/WT/HLA/tumor context/local-only scores"})
        strict_rows.append(base)
        for c in score_cols:
            tmp_df = strict_df[strict_df[c].notna()].copy()
            if tmp_df.empty:
                continue
            tmp = metric_binary(tmp_df[label_col], tmp_df[c])
            for k in [10, 20]:
                tmp.update(precision_recall_at_k(tmp_df, c, label_col, k))
            tmp.update(patient_topn_metrics(tmp_df, c, label_col, 20))
            tmp.update({"track": "clean_science_track", "model": c, "split": "strict_no_existing_overlap", "features": "single existing local score"})
            strict_rows.append(tmp)
    strict_lb = pd.DataFrame(strict_rows)
    for c in ordered:
        if c not in strict_lb.columns:
            strict_lb[c] = np.nan
    strict_lb = strict_lb[ordered + [c for c in strict_lb.columns if c not in ordered]].sort_values(["AUPRC", "patient_hit_rate@20"], ascending=False)
    write_tsv(strict_lb, outdir / "metrics" / "leaderboard_clean_track_strict_no_overlap.tsv")
    md = [
        "# Clean science track",
        "",
        "External public predictor scores are excluded from training features.",
        "",
        f"- Evaluated labeled rows: {len(eval_df):,}",
        f"- Strict no-existing-overlap labeled rows: {len(strict_df):,}",
        f"- Local score columns available: {len(score_cols):,}",
        "",
        "## Boundary",
        "This is a leakage-aware integration scaffold. Any apparent gain still needs source-heldout/no-overlap confirmation before a paper claim.",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "clean_science_track.md")


if __name__ == "__main__":
    main()
