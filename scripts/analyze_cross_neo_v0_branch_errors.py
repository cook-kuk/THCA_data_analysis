#!/usr/bin/env python3
"""CROSS-Neo v0 branch error decomposition for v1 design."""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

from cross_neo_v1_common import V0, V1, add_group_ranks, ensure_v1_dirs, metrics, read_v0_expert_predictions


CORE_EXPERTS = [
    "C_counterfactual_rf",
    "D_structure_geometry_rf",
    "E_quantum_fixed_rf",
    "F_CROSS_all_rf",
    "qk_no_anchor_gamma1",
    "qk_quantum_only_gamma1",
    "late_prespecified_equal_weight_Cw0.5_qk_no_anchor_gamma1",
    "late_prespecified_equal_weight_Cw0.5_qk_quantum_only_gamma1",
]


def top_flags(df: pd.DataFrame) -> pd.DataFrame:
    out = add_group_ranks(df, ["split_name", "fold_id", "expert"], "score")
    out["top5_flag"] = (out["rank_desc"] <= 5).astype(int)
    out["top10_flag"] = (out["rank_desc"] <= 10).astype(int)
    out["positive_low_rank"] = ((out["label"] == 1) & (out["rank_pct_desc"] > 0.5)).astype(int)
    out["negative_top10"] = ((out["label"] == 0) & (out["top10_flag"] == 1)).astype(int)
    return out


def rank_correlation(ranked: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for split, sub in ranked[ranked["expert"].isin(CORE_EXPERTS)].groupby("split_name"):
        wide = sub.pivot_table(index=["fold_id", "sample_id"], columns="expert", values="rank_pct_desc")
        for a, b in itertools.combinations([c for c in CORE_EXPERTS if c in wide.columns], 2):
            x = wide[[a, b]].dropna()
            rows.append(
                {
                    "split_name": split,
                    "expert_a": a,
                    "expert_b": b,
                    "n": len(x),
                    "spearman_rank_corr": float(x[a].corr(x[b], method="spearman")) if len(x) >= 3 else np.nan,
                }
            )
    return pd.DataFrame(rows)


def case_tables(ranked: pd.DataFrame, meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    wide = ranked[ranked["expert"].isin(CORE_EXPERTS)].pivot_table(
        index=["split_name", "fold_id", "sample_id", "label"],
        columns="expert",
        values=["score", "rank_desc", "rank_pct_desc"],
        aggfunc="mean",
    )
    wide.columns = [f"{a}__{b}" for a, b in wide.columns]
    wide = wide.reset_index().merge(meta, on="sample_id", how="left")
    c_rank = "rank_pct_desc__C_counterfactual_rf"
    q_cols = [c for c in wide.columns if c.startswith("rank_pct_desc__qk_")]
    q_score_cols = [c for c in wide.columns if c.startswith("score__qk_")]
    wide["best_qk_rank_pct"] = wide[q_cols].min(axis=1)
    wide["best_qk_score"] = wide[q_score_cols].max(axis=1)
    wide["qk_delta_vs_c"] = wide["best_qk_score"] - wide["score__C_counterfactual_rf"]
    qk_rescue = wide[(wide["label"] == 1) & (wide[c_rank] > 0.5) & (wide["best_qk_rank_pct"] <= 0.25)].copy()
    qk_harm = wide[(wide["label"] == 0) & (wide[c_rank] > 0.25) & (wide["best_qk_rank_pct"] <= 0.15)].copy()
    s_rank = "rank_pct_desc__D_structure_geometry_rf"
    structure_harm = wide[(wide["label"] == 0) & (wide.get(s_rank, np.inf) <= 0.15) & (wide[c_rank] > 0.25)].copy()
    return (
        qk_rescue.sort_values(["best_qk_rank_pct", "qk_delta_vs_c"], ascending=[True, False]),
        qk_harm.sort_values(["best_qk_rank_pct", "qk_delta_vs_c"], ascending=[True, False]),
        structure_harm.sort_values(s_rank),
    )


def oracle_upper_bound(ranked: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for split, sub in ranked[ranked["expert"].isin(CORE_EXPERTS)].groupby("split_name"):
        wide = sub.pivot_table(index=["fold_id", "sample_id", "label"], columns="expert", values="score", aggfunc="mean").reset_index()
        expert_cols = [c for c in CORE_EXPERTS if c in wide.columns]
        for e in expert_cols:
            mm = metrics(wide["label"].to_numpy(), wide[e].to_numpy())
            rows.append({"split_name": split, "strategy": e, "label_informed_oracle": 0, **mm})
        oracle = wide[expert_cols].max(axis=1).where(wide["label"] == 1, wide[expert_cols].min(axis=1))
        mm = metrics(wide["label"].to_numpy(), oracle.to_numpy())
        rows.append({"split_name": split, "strategy": "oracle_label_informed_best_branch", "label_informed_oracle": 1, **mm})
    return pd.DataFrame(rows)


def complementarity(ranked: pd.DataFrame) -> pd.DataFrame:
    rows = []
    core = ranked[ranked["expert"].isin(CORE_EXPERTS)]
    for (split, fold), sub in core.groupby(["split_name", "fold_id"]):
        by = {}
        for expert, g in sub.groupby("expert"):
            top = g.sort_values("score", ascending=False).head(10)
            by[expert] = set(top["sample_id"])
            rows.append(
                {
                    "split_name": split,
                    "fold_id": fold,
                    "expert": expert,
                    "top10_positive_n": int(top["label"].sum()),
                    "top10_negative_n": int((1 - top["label"]).sum()),
                    "top10_precision": float(top["label"].mean()) if len(top) else np.nan,
                }
            )
        for a, b in itertools.combinations(sorted(by), 2):
            rows.append(
                {
                    "split_name": split,
                    "fold_id": fold,
                    "expert": f"overlap::{a}::{b}",
                    "top10_positive_n": np.nan,
                    "top10_negative_n": np.nan,
                    "top10_precision": np.nan,
                    "top10_overlap_n": len(by[a] & by[b]),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    ensure_v1_dirs()
    expert = read_v0_expert_predictions()
    ranked = top_flags(expert)
    master = pd.read_csv(V0 / "master_table.tsv", sep="\t")
    retr = pd.read_csv(V0 / "retrieval_features_by_fold.tsv", sep="\t")
    geom = pd.read_csv(V0 / "structure_geometry_features.tsv", sep="\t")
    ood = pd.read_csv(V0 / "ood_scores.tsv", sep="\t") if (V0 / "ood_scores.tsv").exists() else pd.DataFrame()
    meta_cols = ["sample_id", "peptide_mut", "hla", "hla_supertype", "study"]
    meta = master[meta_cols].copy()
    meta["peptide_length"] = meta["peptide_mut"].astype(str).str.len()
    meta = meta.merge(geom[["sample_id", "structure_missing", "structure_low_confidence", "mean_pLDDT_peptide"]], on="sample_id", how="left")
    if len(ood):
        meta = meta.merge(ood[["sample_id", "ood_score"]].drop_duplicates("sample_id"), on="sample_id", how="left")
    ranked = ranked.merge(meta, on="sample_id", how="left")
    ranked = ranked.merge(
        retr[["split_name", "fold_id", "sample_id", "retrieval_leakage_risk", "near_peptide_similarity_train"]],
        on=["split_name", "fold_id", "sample_id"],
        how="left",
    )
    ranked.to_csv(V1 / "branch_error_decomposition.tsv", sep="\t", index=False)
    corr = rank_correlation(ranked)
    corr.to_csv(V1 / "branch_rank_correlation.tsv", sep="\t", index=False)
    qk_rescue, qk_harm, structure_harm = case_tables(ranked, meta)
    qk_rescue.to_csv(V1 / "qk_rescue_cases.tsv", sep="\t", index=False)
    qk_harm.to_csv(V1 / "qk_harm_cases.tsv", sep="\t", index=False)
    structure_harm.to_csv(V1 / "structure_harm_cases.tsv", sep="\t", index=False)
    oracle = oracle_upper_bound(ranked)
    oracle.to_csv(V1 / "oracle_branch_upper_bound.tsv", sep="\t", index=False)
    comp = complementarity(ranked)
    comp.to_csv(V1 / "branch_complementarity.tsv", sep="\t", index=False)
    dist_rows = []
    for col in ["label", "hla_supertype", "study", "retrieval_leakage_risk", "structure_missing"]:
        for keys, sub in ranked[ranked["expert"].isin(CORE_EXPERTS)].groupby(["split_name", "expert", col], dropna=False):
            dist_rows.append(
                {
                    "split_name": keys[0],
                    "expert": keys[1],
                    "stratum": col,
                    "value": keys[2],
                    "n": len(sub),
                    "mean_score": float(sub["score"].mean()),
                    "median_rank_pct": float(sub["rank_pct_desc"].median()),
                    "label_rate": float(sub["label"].mean()),
                }
            )
    pd.DataFrame(dist_rows).to_csv(V1 / "branch_score_distribution_by_strata.tsv", sep="\t", index=False)
    best_oracle = oracle.sort_values("AUPRC", ascending=False).head(8)
    lines = [
        "# CROSS-Neo v0 Branch Error Decomposition",
        "",
        "Branch ranks are computed inside each v0 split/fold. Late-fusion rows are reconstructed from fold-safe v0 predictions.",
        "",
        f"QK rescue cases: {len(qk_rescue)}",
        f"QK harm cases: {len(qk_harm)}",
        f"Structure harm cases: {len(structure_harm)}",
        "",
        "## Rank Correlation Snapshot",
        corr.sort_values("spearman_rank_corr").head(12).to_markdown(index=False),
        "",
        "## Oracle Upper Bound Snapshot",
        best_oracle.to_markdown(index=False),
        "",
        "Interpretation: oracle rows are label-informed upper bounds for complementarity diagnosis only, not deployable models.",
    ]
    (V1 / "branch_error_decomposition.md").write_text("\n".join(lines) + "\n")
    print(f"[v1-branch] ranked={len(ranked)} qk_rescue={len(qk_rescue)} qk_harm={len(qk_harm)}")


if __name__ == "__main__":
    main()
