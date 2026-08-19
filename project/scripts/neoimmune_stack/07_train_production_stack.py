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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    canon["label_immunogenicity"] = pd.to_numeric(canon["label_immunogenicity"], errors="coerce")

    ext = safe_read_table(outdir / "external_scores" / "all_external_scores.tsv.gz")
    ext = ext[~ext["candidate_id"].astype(str).str.startswith("__")].copy()
    ext["normalized_score"] = pd.to_numeric(ext["normalized_score"], errors="coerce")
    ext_wide = ext.pivot_table(index="candidate_id", columns="model_name", values="normalized_score", aggfunc="mean").reset_index() if not ext.empty else pd.DataFrame({"candidate_id": canon["candidate_id"]})

    local = safe_read_table(outdir / "local_scores" / "all_local_scores.tsv.gz")
    local["normalized_score"] = pd.to_numeric(local["normalized_score"], errors="coerce")
    local_wide = local.pivot_table(index="candidate_id", columns="local_model_name", values="normalized_score", aggfunc="mean").reset_index() if not local.empty else pd.DataFrame({"candidate_id": canon["candidate_id"]})

    d = canon.merge(ext_wide, on="candidate_id", how="left").merge(local_wide, on="candidate_id", how="left", suffixes=("", "_local"))
    feat = peptide_features(d)
    score_cols = [c for c in list(ext_wide.columns) + list(local_wide.columns) if c != "candidate_id" and c in d.columns]
    for c in score_cols:
        feat[c] = pd.to_numeric(d[c], errors="coerce")
    use = d["label_immunogenicity"].notna()
    d.loc[use, "production_stack_score"] = simple_train_predict(pd.concat([d.loc[use, ["candidate_id", "label_immunogenicity", "source_study", "patient_id"]], feat.loc[use]], axis=1), list(feat.columns), "label_immunogenicity", "source_study")
    if (~use).any():
        d.loc[~use, "production_stack_score"] = norm_series(feat.loc[~use].mean(axis=1), True).fillna(0.5)
    d["production_stack_score"] = pd.to_numeric(d["production_stack_score"], errors="coerce").fillna(0.5).clip(0, 1)
    d["model_disagreement_score"] = feat[score_cols].std(axis=1, skipna=True) if score_cols else np.nan
    d["abstention_coverage"] = feat[score_cols].notna().mean(axis=1) if score_cols else 0.0
    flags_path = outdir / "metrics" / "leakage_candidate_flags.tsv"
    if flags_path.exists():
        flags = safe_read_table(flags_path)[["candidate_id", "any_existing_overlap_flag"]]
        d = d.merge(flags, on="candidate_id", how="left")
        d["any_existing_overlap_flag"] = d["any_existing_overlap_flag"].fillna(False).astype(bool)
    else:
        d["any_existing_overlap_flag"] = False

    pred_cols = ["candidate_id", "patient_id", "dataset_source", "peptide_mut", "hla_allele", "label_immunogenicity", "production_stack_score", "model_disagreement_score", "abstention_coverage"]
    write_tsv(d[pred_cols].sort_values("production_stack_score", ascending=False), outdir / "production_track" / "production_stack_predictions.tsv")

    rows = []
    eval_df = d[d["label_immunogenicity"].notna()].copy()
    if not eval_df.empty:
        base = metric_binary(eval_df["label_immunogenicity"], eval_df["production_stack_score"])
        for k in [10, 20, 50]:
            base.update(precision_recall_at_k(eval_df, "production_stack_score", "label_immunogenicity", k))
        base.update(patient_topn_metrics(eval_df, "production_stack_score", "label_immunogenicity", 20))
        base.update({"track": "production_stack_track", "model": "NeoImmune-Stack production stack", "split": "source_grouped_cv_or_existing_labels", "features": "local + frozen external + patient context"})
        rows.append(base)
        for c in score_cols:
            tmp = metric_binary(eval_df["label_immunogenicity"], eval_df[c])
            for k in [10, 20]:
                tmp.update(precision_recall_at_k(eval_df, c, "label_immunogenicity", k))
            tmp.update(patient_topn_metrics(eval_df, c, "label_immunogenicity", 20))
            tmp.update({"track": "production_stack_track", "model": c, "split": "existing_artifact", "features": "single frozen/local score"})
            rows.append(tmp)
    lb = pd.DataFrame(rows)
    ordered = ["track", "model", "split", "features", "n", "positives", "AUROC", "AUPRC", "Precision@10", "Precision@20", "Recall@10", "Recall@20", "Recall@50", "patient_hit_rate@20", "patient_recall@20", "patients_evaluated", "calibration_brier"]
    for c in ordered:
        if c not in lb.columns:
            lb[c] = np.nan
    lb = lb[ordered + [c for c in lb.columns if c not in ordered]].sort_values(["patient_hit_rate@20", "AUPRC"], ascending=False)
    write_tsv(lb, outdir / "metrics" / "leaderboard_production_track.tsv")

    strict_rows = []
    strict_df = eval_df[~eval_df["any_existing_overlap_flag"]].copy() if "any_existing_overlap_flag" in eval_df.columns else eval_df.copy()
    if not strict_df.empty:
        base = metric_binary(strict_df["label_immunogenicity"], strict_df["production_stack_score"])
        for k in [10, 20, 50]:
            base.update(precision_recall_at_k(strict_df, "production_stack_score", "label_immunogenicity", k))
        base.update(patient_topn_metrics(strict_df, "production_stack_score", "label_immunogenicity", 20))
        base.update({"track": "production_stack_track", "model": "NeoImmune-Stack production stack", "split": "strict_no_existing_overlap", "features": "local + frozen external + patient context"})
        strict_rows.append(base)
        for c in score_cols:
            tmp_df = strict_df[strict_df[c].notna()].copy()
            if tmp_df.empty:
                continue
            tmp = metric_binary(tmp_df["label_immunogenicity"], tmp_df[c])
            for k in [10, 20]:
                tmp.update(precision_recall_at_k(tmp_df, c, "label_immunogenicity", k))
            tmp.update(patient_topn_metrics(tmp_df, c, "label_immunogenicity", 20))
            tmp.update({"track": "production_stack_track", "model": c, "split": "strict_no_existing_overlap", "features": "single frozen/local score"})
            strict_rows.append(tmp)
    strict_lb = pd.DataFrame(strict_rows)
    for c in ordered:
        if c not in strict_lb.columns:
            strict_lb[c] = np.nan
    strict_lb = strict_lb[ordered + [c for c in strict_lb.columns if c not in ordered]].sort_values(["patient_hit_rate@20", "AUPRC"], ascending=False)
    write_tsv(strict_lb, outdir / "metrics" / "leaderboard_production_track_strict_no_overlap.tsv")
    md = [
        "# Production stack track",
        "",
        "Frozen external scores are allowed here and are explicitly barred from clean-science claims.",
        "",
        f"- Evaluated labeled rows: {len(eval_df):,}",
        f"- Strict no-existing-overlap labeled rows: {len(strict_df):,}",
        f"- Feature score columns available: {len(score_cols):,}",
        "",
        "## Boundary",
        "Production-stack performance is practical prioritization evidence, not proof of a new clean algorithm and not proof of clinical vaccine efficacy.",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "production_stack_track.md")


if __name__ == "__main__":
    main()
