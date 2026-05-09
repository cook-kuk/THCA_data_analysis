#!/usr/bin/env python3
"""Fold-safe ESM2 + QK/source-gated fusion tests.

The important guardrail: for every outer fold, model-pair and weight selection
uses only out-of-fold predictions for the outer-train rows. Outer-test labels
are used only once for final evaluation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import OUT, SEED, ensure_dirs, metrics


PRIMARY = ["exact_peptide_hla_holdout", "near_peptide_cluster_holdout", "hla_stratified_group_5fold", "hla_supertype_heldout"]
BASE_CANDIDATES = [
    "prespecified_rf_qk_no_anchor_w0.5",
    "rule_gate_rf_qk_fallback_train_selected",
    "v2_cf_plm_lr",
    "v2_multimodal_lr",
    "v2_rule_moe_cf_plm_structure",
    "anchor_rf",
]
ESM_CANDIDATES = [
    "v2_cf_esm2_150m_lr_fast",
    "v2_multimodal_esm2_150m_lr_fast",
    "v2_cf_esm2_650m_lr_fast",
    "v2_multimodal_esm2_650m_lr_fast",
    "v2_cf_esm2_35m_lr_fast",
    "v2_multimodal_esm2_35m_lr_fast",
]
WEIGHTS = [0.0, 0.25, 0.5, 0.75, 1.0]  # base weight; 1.0 means no ESM2 contribution.


def load_predictions() -> pd.DataFrame:
    parts = [pd.read_csv(OUT / "predictions/all_predictions.tsv", sep="\t")]
    for p in sorted((OUT / "predictions").glob("esm2_*_fast_predictions.tsv")):
        parts.append(pd.read_csv(p, sep="\t"))
    df = pd.concat(parts, ignore_index=True, sort=False)
    df["row_id"] = df["row_id"].astype(str)
    df["label"] = pd.to_numeric(df["label"], errors="coerce").fillna(0).astype(int)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score"])
    return df


def percentile_against_train(train_scores: np.ndarray, scores: np.ndarray) -> np.ndarray:
    train_scores = np.asarray(train_scores, dtype=float)
    scores = np.asarray(scores, dtype=float)
    if len(train_scores) < 10 or np.nanstd(train_scores) <= 1e-12:
        lo, hi = np.nanmin(scores), np.nanmax(scores)
        return (scores - lo) / (hi - lo + 1e-9)
    return np.asarray([(train_scores <= s).mean() for s in scores], dtype=float)


def model_scores(pred: pd.DataFrame, split: str, model: str, row_ids: list[str]) -> pd.DataFrame:
    sub = pred[(pred["split_name"].eq(split)) & (pred["model_name"].eq(model)) & (pred["row_id"].isin(row_ids))].copy()
    if sub.empty:
        return sub
    sub = sub.sort_values(["row_id", "fold_id"]).drop_duplicates("row_id", keep="first")
    return sub[["row_id", "label", "score"]]


def candidate_train_score(pred: pd.DataFrame, split: str, train_ids: list[str], base: str, esm: str, w: float) -> tuple[float, float, int]:
    b = model_scores(pred, split, base, train_ids).rename(columns={"score": "base_score"})
    e = model_scores(pred, split, esm, train_ids).rename(columns={"score": "esm_score", "label": "label_esm"})
    m = b.merge(e, on="row_id", how="inner")
    if len(m) < 30 or m["label"].nunique() < 2:
        return (-np.inf, -np.inf, len(m))
    b_pct = percentile_against_train(m["base_score"].values, m["base_score"].values)
    e_pct = percentile_against_train(m["esm_score"].values, m["esm_score"].values)
    score = w * b_pct + (1 - w) * e_pct
    met = metrics(m["label"].values, score, ks=(10,))
    return (float(met.get("AUPRC", -np.inf)), float(met.get("top10_precision", -np.inf)), len(m))


def fuse_fold(pred: pd.DataFrame, split: str, fold_id: str, fold: pd.DataFrame) -> list[dict]:
    train_ids = fold.loc[fold["role"].eq("train"), "row_id"].astype(str).tolist()
    test_ids = fold.loc[fold["role"].eq("test"), "row_id"].astype(str).tolist()
    if not test_ids:
        return []
    grid = []
    for base in BASE_CANDIDATES:
        for esm in ESM_CANDIDATES:
            for w in WEIGHTS:
                auprc, top10, n_train = candidate_train_score(pred, split, train_ids, base, esm, w)
                grid.append({"base": base, "esm": esm, "base_weight": w, "train_AUPRC": auprc, "train_top10": top10, "n_train_pred": n_train})
    grid_df = pd.DataFrame(grid)
    usable = grid_df[np.isfinite(grid_df["train_AUPRC"])].copy()
    if usable.empty:
        return []
    # Deterministic nested selection: primary train AUPRC, then train top10,
    # then stronger base contribution for stability.
    sel = usable.sort_values(["train_AUPRC", "train_top10", "base_weight"], ascending=[False, False, False]).iloc[0]
    rows = []
    variants = [
        ("nested_train_selected_esm2_qk_gate", sel["base"], sel["esm"], float(sel["base_weight"]), "reviewer_safe_internal_locked", "nested_oof_train_selected"),
        ("prespecified_qk75_esm25_exact150m", "prespecified_rf_qk_no_anchor_w0.5", "v2_cf_esm2_150m_lr_fast", 0.75, "exploratory_descriptive_only", "posthoc_prespecified_for_audit"),
        ("prespecified_qk50_esm650m_near", "rule_gate_rf_qk_fallback_train_selected", "v2_cf_esm2_650m_lr_fast", 0.50, "exploratory_descriptive_only", "posthoc_prespecified_for_audit"),
        ("prespecified_v2base75_esm650m", "v2_multimodal_lr", "v2_multimodal_esm2_650m_lr_fast", 0.75, "exploratory_descriptive_only", "posthoc_prespecified_for_audit"),
    ]
    for name, base, esm, w, status, selection in variants:
        b_train = model_scores(pred, split, base, train_ids).rename(columns={"score": "base_train_score"})
        e_train = model_scores(pred, split, esm, train_ids).rename(columns={"score": "esm_train_score"})
        b_test = model_scores(pred, split, base, test_ids).rename(columns={"score": "base_score"})
        e_test = model_scores(pred, split, esm, test_ids).rename(columns={"score": "esm_score", "label": "label_esm"})
        m = b_test.merge(e_test, on="row_id", how="inner")
        if m.empty:
            continue
        base_train_scores = b_train["base_train_score"].values if not b_train.empty else b_test["base_score"].values
        esm_train_scores = e_train["esm_train_score"].values if not e_train.empty else e_test["esm_score"].values
        b_pct = percentile_against_train(base_train_scores, m["base_score"].values)
        e_pct = percentile_against_train(esm_train_scores, m["esm_score"].values)
        score = w * b_pct + (1 - w) * e_pct
        for rid, lab, sc in zip(m["row_id"].values, m["label"].astype(int).values, score):
            rows.append({
                "split_name": split,
                "fold_id": fold_id,
                "row_id": rid,
                "label": lab,
                "score": float(sc),
                "model_name": name,
                "model_family": "esm2_qk_foldsafe_gate",
                "claim_status": status,
                "source_script": "cross_neo_v2_esm2_qk_gate_fusion",
                "selected_base": base,
                "selected_esm": esm,
                "base_weight": w,
                "selection_status": selection,
                "train_selected_AUPRC": float(sel["train_AUPRC"]),
                "train_selected_top10": float(sel["train_top10"]),
            })
    return rows


def main() -> None:
    ensure_dirs()
    pred = load_predictions()
    rows = []
    selected_rows = []
    for split in PRIMARY:
        split_file = OUT / "splits" / f"{split}.tsv"
        if not split_file.exists():
            continue
        sdf = pd.read_csv(split_file, sep="\t")
        for fold_id, fold in sdf.groupby("fold_id"):
            out = fuse_fold(pred, split, str(fold_id), fold)
            rows.extend(out)
            for r in out:
                if r["model_name"] == "nested_train_selected_esm2_qk_gate":
                    selected_rows.append({k: r[k] for k in ["split_name", "fold_id", "selected_base", "selected_esm", "base_weight", "train_selected_AUPRC", "train_selected_top10"]})
                    break
    fused = pd.DataFrame(rows)
    if fused.empty:
        raise SystemExit("no fused predictions produced")
    fused["rank"] = fused.groupby(["split_name", "fold_id", "model_name"])["score"].rank(method="first", ascending=False)
    fused["rank_pct"] = fused["rank"] / fused.groupby(["split_name", "fold_id", "model_name"])["row_id"].transform("size")
    fused.to_csv(OUT / "predictions/esm2_qk_gate_fusion_predictions.tsv", sep="\t", index=False)
    pd.DataFrame(selected_rows).to_csv(OUT / "metrics/esm2_qk_gate_selected_weights.tsv", sep="\t", index=False)

    mets = []
    for (split, model), g in fused.groupby(["split_name", "model_name"]):
        mets.append({"split_name": split, "model_name": model, "model_family": g["model_family"].iloc[0], "claim_status": g["claim_status"].iloc[0], **metrics(g["label"].values, g["score"].values)})
    met = pd.DataFrame(mets)
    met.to_csv(OUT / "metrics/esm2_qk_gate_fusion_metrics.tsv", sep="\t", index=False)

    base = pd.read_csv(OUT / "metrics/all_model_all_split_metrics_plus_runpod_esm2.tsv", sep="\t") if (OUT / "metrics/all_model_all_split_metrics_plus_runpod_esm2.tsv").exists() else pd.read_csv(OUT / "metrics/all_model_all_split_metrics.tsv", sep="\t")
    combined = pd.concat([base, met], ignore_index=True, sort=False)
    combined.to_csv(OUT / "metrics/all_model_all_split_metrics_plus_runpod_esm2_gate.tsv", sep="\t", index=False)
    comp_rows = []
    for split in PRIMARY:
        b = base[base["split_name"].eq(split)].sort_values(["AUPRC", "top10_precision"], ascending=False).head(1)
        f = met[(met["split_name"].eq(split)) & (met["claim_status"].eq("reviewer_safe_internal_locked"))].sort_values(["AUPRC", "top10_precision"], ascending=False).head(1)
        c = combined[combined["split_name"].eq(split)].sort_values(["AUPRC", "top10_precision"], ascending=False).head(1)
        comp_rows.append({
            "split_name": split,
            "pre_gate_best_model": b["model_name"].iloc[0] if not b.empty else "",
            "pre_gate_best_AUPRC": b["AUPRC"].iloc[0] if not b.empty else np.nan,
            "pre_gate_best_top10": b["top10_precision"].iloc[0] if not b.empty else np.nan,
            "gate_best_model": f["model_name"].iloc[0] if not f.empty else "",
            "gate_best_AUPRC": f["AUPRC"].iloc[0] if not f.empty else np.nan,
            "gate_best_top10": f["top10_precision"].iloc[0] if not f.empty else np.nan,
            "combined_best_model": c["model_name"].iloc[0] if not c.empty else "",
            "combined_best_AUPRC": c["AUPRC"].iloc[0] if not c.empty else np.nan,
            "combined_best_top10": c["top10_precision"].iloc[0] if not c.empty else np.nan,
        })
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(OUT / "metrics/esm2_qk_gate_primary_comparison.tsv", sep="\t", index=False)

    xlsx = OUT / "CROSS_Neo_v2_all_results_summary_plus_runpod_esm2_gate.xlsx"
    with pd.ExcelWriter(xlsx) as writer:
        combined.head(30000).to_excel(writer, sheet_name="all_plus_gate", index=False)
        comp.to_excel(writer, sheet_name="primary_gate_compare", index=False)
        met.to_excel(writer, sheet_name="gate_metrics", index=False)
        pd.DataFrame(selected_rows).to_excel(writer, sheet_name="selected_weights", index=False)

    lines = [
        "# ESM2 + QK Gate Fusion Report",
        "",
        "Nested gate uses only same-split out-of-fold predictions from outer-train rows to select base model, ESM2 expert, and base/ESM2 weight.",
        "",
        "## Primary Locked Split Comparison",
        "",
        comp.to_markdown(index=False),
        "",
        "## Gate Metrics",
        "",
        met.sort_values(["split_name", "AUPRC", "top10_precision"], ascending=[True, False, False]).to_markdown(index=False),
        "",
        "## Selected Weights",
        "",
        pd.DataFrame(selected_rows).to_markdown(index=False),
        "",
        f"Workbook: `{xlsx}`",
    ]
    (OUT / "ESM2_QK_gate_fusion_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v2-esm2-qk-gate] fused_rows={len(fused)} metrics={len(met)} xlsx={xlsx}")
    print(comp.to_string(index=False))


if __name__ == "__main__":
    main()
