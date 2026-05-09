#!/usr/bin/env python3
"""Collate all local neoantigen algorithm/test-set metrics into one XLSX.

The workbook is intentionally strict: it does not just rank by best score, but
also records leakage-control tier, public-predictor caveats, source-heldout
status, exploratory flags, and whether a row is reviewer-safe for a main claim.
"""

from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
OUTDIR = REPO / "project/results/cross_neo_v1_lockdown"
XLSX = OUTDIR / "NEOANTIGEN_ALL_ALGORITHMS_ALL_TESTSETS_FAIR_COMPARISON.xlsx"
SUMMARY = OUTDIR / "NEOANTIGEN_ALL_ALGORITHMS_ALL_TESTSETS_FAIR_COMPARISON_SUMMARY.md"

PUBLIC_PREDICTORS = [
    "MHCflurry",
    "NetMHCpan",
    "NetMHCstabpan",
    "BigMHC",
    "PRIME",
    "MixMHCpred",
    "DeepImmuno",
    "TransPHLA",
]

PRIMARY_SPLITS = [
    "exact_peptide_hla_holdout",
    "near_peptide_cluster_holdout",
    "hla_stratified_group_5fold",
    "hla_supertype_heldout",
]

CANONICAL_FILES = [
    ("v1_lockdown", "v1_lockdown_comparison", "project/results/cross_neo_v1_lockdown/v1_lockdown_comparison.tsv"),
    ("v1_lockdown", "locked_anchor_metrics", "project/results/cross_neo_v1_lockdown/locked_anchor_metrics.tsv"),
    ("v1_lockdown", "foldsafe_fusion_metrics", "project/results/cross_neo_v1_lockdown/foldsafe_fusion_metrics.tsv"),
    ("v1_lockdown", "source_topk_rescue_metrics", "project/results/cross_neo_v1_lockdown/source_topk_rescue_metrics.tsv"),
    ("v1_lockdown", "primary_delta_vs_anchor", "project/results/cross_neo_v1_lockdown/primary_delta_vs_anchor.tsv"),
    ("v1_lockdown", "qk_net_benefit", "project/results/cross_neo_v1_lockdown/qk_net_benefit_by_split.tsv"),
    ("v1_lockdown", "source_root_cause", "project/results/cross_neo_v1_lockdown/source_root_cause_matrix.tsv"),
    ("v1_lockdown", "public_overlap_manifest", "project/results/cross_neo_v1_lockdown/public_overlap_manifest.tsv"),
    ("v1_diagnostic", "v1_model_comparison", "project/results/cross_neo_v1/v1_model_comparison.tsv"),
    ("v1_diagnostic", "gated_moe_metrics", "project/results/cross_neo_v1/gated_moe_metrics.tsv"),
    ("v1_diagnostic", "source_balanced_metrics", "project/results/cross_neo_v1/source_balanced_metrics.tsv"),
    ("v1_diagnostic", "pu_corrected_metrics", "project/results/cross_neo_v1/pu_corrected_metrics.tsv"),
    ("v1_diagnostic", "pu_ranking_metrics", "project/results/cross_neo_v1/pu_ranking_metrics.tsv"),
    ("v1_diagnostic", "decoy_focal_metrics", "project/results/cross_neo_v1/decoy_focal_metrics.tsv"),
    ("v0", "metrics_by_split", "project/results/cross_neo_v0/metrics_by_split.tsv"),
    ("v0", "late_fusion_oof_metrics", "project/results/cross_neo_v0/late_fusion_oof_metrics.tsv"),
    ("v0", "qk_fallback_metrics", "project/results/cross_neo_v0/qk_fallback_metrics.tsv"),
    ("v0", "source_heldout_metrics", "project/results/cross_neo_v0/source_heldout_metrics.tsv"),
    ("v0", "evaluation_contract_metrics", "project/results/cross_neo_v0/evaluation_contract_metrics.tsv"),
    ("clean_neobench", "clean_neobench_leaderboard", "project/results/clean_neobench_barneo_2026_05_09/clean_neobench_leaderboard.tsv"),
    ("clean_neobench", "clean_neobench_split_metrics", "project/results/clean_neobench_barneo_2026_05_09/clean_neobench_split_metrics.tsv"),
    ("clean_neobench", "barneo_bma_method_weights", "project/results/clean_neobench_barneo_2026_05_09/barneo_bma_method_weights.tsv"),
    ("legacy_wave", "wave3_algorithm_sweep_results", "project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/algorithm_sweep_results.tsv"),
    ("legacy_wave", "wave3_algorithm_sweep_pivot", "project/results/p_neo_bayesian_2026_05_09/wave3_algorithm_sweep/algorithm_sweep_pivot.tsv"),
    ("legacy_wave", "wave2_combined_results", "project/results/p_neo_bayesian_2026_05_09/wave2/wave2_combined_results.tsv"),
    ("legacy_wave", "wave2_structure_ablation", "project/results/p_neo_bayesian_2026_05_09/wave2/structure_ablation_results.tsv"),
    ("legacy_wave", "wave2_qk_results", "project/results/p_neo_bayesian_2026_05_09/wave2/qk_results.tsv"),
    ("legacy_wave", "wave2_vqc_results", "project/results/p_neo_bayesian_2026_05_09/wave2/vqc_results.tsv"),
    ("legacy_wave", "wave1_results_summary", "project/results/p_neo_bayesian_2026_05_09/wave1/results_summary.tsv"),
    ("legacy_wave", "root_results_summary", "project/results/p_neo_bayesian_2026_05_09/results_summary.tsv"),
    ("legacy_wave", "fallback_sweep_summary", "project/results/p_neo_bayesian_2026_05_09/clean_neo_v0_2026_05_09/fallback_algorithm_sweep_summary.tsv"),
]


def clean_sheet_name(name: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_]+", "_", name)[:31]
    return name or "sheet"


def read_table(path: Path) -> pd.DataFrame | None:
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        if path.suffix.lower() == ".csv":
            return pd.read_csv(path)
        return pd.read_csv(path, sep="\t")
    except Exception:
        return None


def first_present(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c in df.columns:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def combine_algorithm(df: pd.DataFrame, table_name: str) -> pd.Series:
    cols = df.columns
    if {"feature_group", "model"}.issubset(cols):
        return df["feature_group"].astype(str) + " / " + df["model"].astype(str)
    if {"fusion", "status"}.issubset(cols):
        return df["fusion"].astype(str)
    for candidates in [
        ["method"],
        ["algorithm"],
        ["model"],
        ["branch"],
        ["fusion"],
        ["predictor"],
        ["feature_group"],
        ["expert"],
    ]:
        col = first_present(df, candidates)
        if col:
            return df[col].astype(str)
    return pd.Series([table_name] * len(df), index=df.index)


def combine_testset(df: pd.DataFrame, table_name: str) -> pd.Series:
    if "heldout_study" in df.columns:
        return "source_heldout_" + df["heldout_study"].astype(str)
    for candidates in [
        ["split_name"],
        ["split"],
        ["testset"],
        ["test_set"],
        ["dataset"],
        ["cohort"],
        ["subset"],
        ["evaluation_set"],
        ["condition"],
        ["heldout"],
    ]:
        col = first_present(df, candidates)
        if col:
            return df[col].astype(str)
    return pd.Series([table_name] * len(df), index=df.index)


def numeric_col(df: pd.DataFrame, candidates: list[str]) -> pd.Series:
    col = first_present(df, candidates)
    if col:
        return pd.to_numeric(df[col], errors="coerce")
    return pd.Series([np.nan] * len(df), index=df.index)


def infer_leakage_tier(phase: str, table_name: str, source_path: str) -> str:
    s = f"{phase} {table_name} {source_path}".lower()
    if "cross_neo_v1_lockdown" in s:
        return "locked_foldsafe_audited"
    if "clean_neobench" in s:
        return "clean_neobench_audited"
    if "cross_neo_v1" in s:
        return "v1_diagnostic_not_final_lock"
    if "cross_neo_v0" in s:
        return "v0_locked_internal"
    if "wave3_algorithm_sweep" in s:
        return "public_comparator_or_legacy_sweep_caveated"
    return "legacy_exploratory"


def infer_testset_severity(testset: str) -> str:
    t = str(testset).lower()
    if "source_heldout" in t or "study" in t:
        return "hard_source_or_study_shift"
    if "near_peptide" in t or "cluster" in t:
        return "hard_near_neighbor_holdout"
    if "hla_supertype" in t or "hla_stratified" in t or "hla" in t:
        return "hard_hla_holdout"
    if "exact_peptide" in t or "pmhc" in t:
        return "hard_exact_pmhc_holdout"
    if "repeated" in t or "internal" in t or "cv" in t:
        return "internal_cv"
    return "other_or_legacy"


def public_predictor_flag(algorithm: str) -> bool:
    a = str(algorithm).lower()
    return any(p.lower() in a for p in PUBLIC_PREDICTORS)


def infer_claim_status(row: pd.Series) -> str:
    alg = str(row.get("algorithm", ""))
    tier = str(row.get("leakage_control_tier", ""))
    testset = str(row.get("testset", ""))
    source_table = str(row.get("table_name", ""))
    status = str(row.get("status", ""))
    family = str(row.get("method_family", ""))
    if public_predictor_flag(alg):
        return "public_pretrained_comparator_caveated_not_feature"
    if "exploratory" in status.lower() or "exploratory" in alg.lower():
        return "exploratory_descriptive_only"
    if "qk_" in alg.lower() and not any(x in alg.lower() for x in ["fusion", "prespecified", "nested", "rule_gate"]):
        return "qk_diagnostic_or_fallback_not_main_claim"
    if "source_heldout" in testset.lower() and pd.to_numeric(pd.Series([row.get("top10_precision")]), errors="coerce").iloc[0] == 0:
        return "source_shift_failure"
    if tier == "locked_foldsafe_audited" and any(x in family.lower() for x in ["prespecified", "nested", "rule_gate", "clean_anchor"]):
        return "reviewer_safe_internal_locked"
    if tier == "clean_neobench_audited":
        return "clean_benchmark_caveated"
    if "legacy" in tier or "wave" in source_table:
        return "legacy_or_exploratory"
    return "diagnostic"


def normalize_metrics(df: pd.DataFrame, phase: str, table_name: str, rel_path: str) -> pd.DataFrame:
    out = pd.DataFrame()
    out["phase"] = phase
    out["table_name"] = table_name
    out["source_file"] = rel_path
    out["algorithm"] = combine_algorithm(df, table_name)
    out["testset"] = combine_testset(df, table_name)
    out["n"] = numeric_col(df, ["n", "N", "n_test", "num_samples", "local_n"])
    out["n_pos"] = numeric_col(df, ["n_pos", "positives", "positive", "pos", "n_positive"])
    out["prevalence"] = numeric_col(df, ["prevalence", "positive_rate", "pos_rate"])
    out["AUPRC"] = numeric_col(df, ["AUPRC", "auprc", "average_precision", "AP", "AveragePrecision"])
    out["AUROC"] = numeric_col(df, ["AUROC", "auroc", "ROC_AUC", "roc_auc", "AUC", "auc"])
    out["Brier"] = numeric_col(df, ["Brier", "brier", "brier_score"])
    out["ECE"] = numeric_col(df, ["ECE", "ece"])
    for k in (5, 10, 20):
        out[f"top{k}_precision"] = numeric_col(df, [f"top{k}_precision", f"precision_at_{k}", f"precision@{k}", f"top{k}", f"p@{k}"])
        out[f"recall_at_{k}"] = numeric_col(df, [f"recall_at_{k}", f"recall@{k}"])
        out[f"enrichment_at_{k}"] = numeric_col(df, [f"enrichment_at_{k}", f"enrichment@{k}"])
    for col in ["status", "method_family", "fusion_family", "split_name", "heldout_study", "feature_group", "model", "branch", "fusion", "method"]:
        if col in df.columns:
            out[col] = df[col].astype(str)
    out["leakage_control_tier"] = infer_leakage_tier(phase, table_name, rel_path)
    out["testset_severity"] = out["testset"].map(infer_testset_severity)
    out["public_predictor_forbidden_as_feature"] = out["algorithm"].map(public_predictor_flag)
    out["primary_rank_metric"] = np.where(out["AUPRC"].notna(), out["AUPRC"], out["AUROC"])
    out["is_primary_locked_split"] = out["testset"].isin(PRIMARY_SPLITS)
    out["claim_status"] = out.apply(infer_claim_status, axis=1)
    metric_cols = ["AUPRC", "AUROC", "top5_precision", "top10_precision", "top20_precision", "Brier", "ECE"]
    has_metric = out[metric_cols].notna().any(axis=1)
    return out[has_metric].copy()


def discover_extra_metric_files() -> list[tuple[str, str, str]]:
    roots = [
        REPO / "project/results/p_neo_bayesian_2026_05_09",
        REPO / "project/results/clean_neobench_barneo_2026_05_09",
        REPO / "project/results/cross_neo_v1",
    ]
    known = {str(REPO / p) for _, _, p in CANONICAL_FILES}
    out = []
    patterns = ["*metrics*.tsv", "*results*.tsv", "*summary*.tsv", "*leaderboard*.tsv", "*auroc*.tsv"]
    deny = ["prediction", "predictions", "candidate_scores", "cases", "features", "input", "output", "master", "manifest", "weights"]
    for root in roots:
        if not root.exists():
            continue
        for pat in patterns:
            for path in root.rglob(pat):
                if str(path) in known:
                    continue
                name = path.name.lower()
                if any(d in name for d in deny):
                    continue
                if path.stat().st_size > 8_000_000:
                    continue
                rel = str(path.relative_to(REPO))
                phase = "auto_discovered"
                out.append((phase, path.stem[:40], rel))
    seen = set()
    unique = []
    for item in out:
        if item[2] not in seen:
            seen.add(item[2])
            unique.append(item)
    return unique


def build_all_metrics() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    manifest_rows = []
    raw_sheets: dict[str, pd.DataFrame] = {}
    normalized = []
    file_specs = CANONICAL_FILES + discover_extra_metric_files()
    for phase, table_name, rel in file_specs:
        path = REPO / rel
        df = read_table(path)
        status = "loaded" if df is not None else "missing_or_unreadable"
        manifest_rows.append(
            {
                "phase": phase,
                "table_name": table_name,
                "source_file": rel,
                "status": status,
                "rows": 0 if df is None else len(df),
                "columns": "" if df is None else ",".join(map(str, df.columns)),
            }
        )
        if df is None or df.empty:
            continue
        sheet = clean_sheet_name(table_name)
        raw_sheets[sheet] = df.head(5000).copy()
        norm = normalize_metrics(df, phase, table_name, rel)
        if not norm.empty:
            normalized.append(norm)
    all_metrics = pd.concat(normalized, ignore_index=True) if normalized else pd.DataFrame()
    if not all_metrics.empty:
        sort_cols = ["testset_severity", "testset", "AUPRC", "top10_precision", "AUROC"]
        all_metrics = all_metrics.sort_values(sort_cols, ascending=[True, True, False, False, False])
        all_metrics["rank_within_testset_by_AUPRC"] = all_metrics.groupby("testset")["AUPRC"].rank(ascending=False, method="min")
        all_metrics["rank_within_testset_by_top10"] = all_metrics.groupby("testset")["top10_precision"].rank(ascending=False, method="min")
    return all_metrics, pd.DataFrame(manifest_rows), raw_sheets


def best_by_testset(all_metrics: pd.DataFrame) -> pd.DataFrame:
    if all_metrics.empty:
        return all_metrics
    d = all_metrics.copy()
    d["safe_sort"] = np.where(d["claim_status"].eq("reviewer_safe_internal_locked"), 0, 1)
    return (
        d.sort_values(["testset", "safe_sort", "AUPRC", "top10_precision", "AUROC"], ascending=[True, True, False, False, False])
        .groupby("testset")
        .head(5)
        .reset_index(drop=True)
    )


def dedupe_metrics(all_metrics: pd.DataFrame) -> pd.DataFrame:
    if all_metrics.empty:
        return all_metrics
    priority = {
        "locked_foldsafe_audited": 0,
        "clean_neobench_audited": 1,
        "v1_diagnostic_not_final_lock": 2,
        "v0_locked_internal": 3,
        "public_comparator_or_legacy_sweep_caveated": 4,
        "legacy_exploratory": 5,
    }
    d = all_metrics.copy()
    d["_tier_priority"] = d["leakage_control_tier"].map(priority).fillna(9)
    metric_key = [
        "algorithm",
        "testset",
        "n",
        "n_pos",
        "AUPRC",
        "AUROC",
        "top5_precision",
        "top10_precision",
        "top20_precision",
    ]
    present = [c for c in metric_key if c in d.columns]
    d = d.sort_values(["_tier_priority", "source_file"]).drop_duplicates(present, keep="first")
    return d.drop(columns=["_tier_priority"])


def fair_leaderboard(all_metrics: pd.DataFrame) -> pd.DataFrame:
    if all_metrics.empty:
        return all_metrics
    allowed = all_metrics[
        all_metrics["claim_status"].isin(
            [
                "reviewer_safe_internal_locked",
                "clean_benchmark_caveated",
                "qk_diagnostic_or_fallback_not_main_claim",
                "source_shift_failure",
            ]
        )
    ].copy()
    return allowed.sort_values(["is_primary_locked_split", "testset", "AUPRC", "top10_precision"], ascending=[False, True, False, False])


def primary_matrix(all_metrics: pd.DataFrame, value: str) -> pd.DataFrame:
    if all_metrics.empty:
        return pd.DataFrame()
    d = all_metrics[all_metrics["testset"].isin(PRIMARY_SPLITS)].copy()
    if d.empty:
        return pd.DataFrame()
    d = d.sort_values([value, "AUPRC", "top10_precision"], ascending=False).drop_duplicates(["algorithm", "testset"])
    return d.pivot(index="algorithm", columns="testset", values=value).reset_index()


def write_summary(all_metrics: pd.DataFrame, manifest: pd.DataFrame) -> None:
    dedup = dedupe_metrics(all_metrics)
    primary = dedup[dedup["testset"].isin(PRIMARY_SPLITS)]
    best_primary = (
        primary.sort_values(["testset", "AUPRC", "top10_precision"], ascending=[True, False, False])
        .groupby("testset")
        .head(3)
        if not primary.empty
        else pd.DataFrame()
    )
    lines = [
        "# Neoantigen All-Algorithm / All-Testset Fair Comparison",
        "",
        f"Workbook: `{XLSX.relative_to(REPO)}`",
        "",
        "This workbook is strict: it includes legacy/exploratory/public-comparator rows, but flags claim status and leakage-control tier so they are not mixed with reviewer-safe locked rows.",
        "",
        "## Counts",
        "",
        f"- metric rows raw: {len(all_metrics)}",
        f"- metric rows deduplicated for ranking: {len(dedup)}",
        f"- source tables scanned: {len(manifest)}",
        f"- loaded source tables: {int(manifest['status'].eq('loaded').sum()) if len(manifest) else 0}",
        "",
        "## Best Primary Locked Rows",
        "",
        best_primary[
            [
                "testset",
                "algorithm",
                "n",
                "n_pos",
                "prevalence",
                "AUPRC",
                "AUROC",
                "top5_precision",
                "top10_precision",
                "claim_status",
            ]
        ].to_markdown(index=False)
        if not best_primary.empty
        else "No primary rows found.",
        "",
        "## Use Rules",
        "",
        "- Public pretrained predictors are comparator-only unless row-level public overlap is resolved.",
        "- QK standalone rows are diagnostic/fallback only; no quantum advantage claim.",
        "- Source-heldout failures remain visible and should not be filtered out.",
        "- Primary ranking should use AUPRC and top-k before AUROC.",
    ]
    SUMMARY.write_text("\n".join(lines) + "\n")


def write_workbook() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    all_metrics, manifest, raw_sheets = build_all_metrics()
    readme = pd.DataFrame(
        [
            {
                "field": "purpose",
                "value": "Fair all-algorithm/all-testset neoantigen comparison with leakage and claim-status flags.",
            },
            {"field": "primary_metrics", "value": "AUPRC, top5/top10/top20 precision, enrichment over prevalence."},
            {"field": "secondary_metrics", "value": "AUROC, Brier, ECE."},
            {"field": "forbidden_feature_rule", "value": "MHCflurry/NetMHCpan/NetMHCstabpan/BigMHC/PRIME/MixMHCpred/DeepImmuno scores are comparator-only unless explicitly audited."},
            {"field": "current_decision", "value": "CROSS-Neo v1 remains HOLD: internal locked fusion improves, source-heldout/public overlap block KEEP."},
        ]
    )
    dedup = dedupe_metrics(all_metrics)
    with pd.ExcelWriter(XLSX, engine="xlsxwriter") as writer:
        readme.to_excel(writer, index=False, sheet_name="README")
        all_metrics.to_excel(writer, index=False, sheet_name="all_metrics_long")
        dedup.to_excel(writer, index=False, sheet_name="all_metrics_dedup")
        fair_leaderboard(dedup).to_excel(writer, index=False, sheet_name="fair_leaderboard")
        best_by_testset(dedup).to_excel(writer, index=False, sheet_name="best5_by_testset")
        primary_matrix(dedup, "AUPRC").to_excel(writer, index=False, sheet_name="primary_AUPRC_matrix")
        primary_matrix(dedup, "top10_precision").to_excel(writer, index=False, sheet_name="primary_top10_matrix")
        manifest.to_excel(writer, index=False, sheet_name="source_manifest")
        for sheet_name, rel in [
            ("v1_primary_delta", "project/results/cross_neo_v1_lockdown/primary_delta_vs_anchor.tsv"),
            ("v1_qk_net", "project/results/cross_neo_v1_lockdown/qk_net_benefit_by_split.tsv"),
            ("v1_source_root", "project/results/cross_neo_v1_lockdown/source_root_cause_matrix.tsv"),
            ("v1_reviewer_matrix", "project/results/cross_neo_v1_lockdown/reviewer_response_matrix.tsv"),
            ("v1_public_overlap", "project/results/cross_neo_v1_lockdown/public_overlap_manifest.tsv"),
            ("clean_neobench_board", "project/results/clean_neobench_barneo_2026_05_09/clean_neobench_leaderboard.tsv"),
        ]:
            df = read_table(REPO / rel)
            if df is not None:
                df.to_excel(writer, index=False, sheet_name=sheet_name)
        for name, df in list(raw_sheets.items())[:40]:
            if name in writer.sheets:
                name = clean_sheet_name(f"raw_{name}")
            df.to_excel(writer, index=False, sheet_name=name)

        workbook = writer.book
        header_fmt = workbook.add_format({"bold": True, "bg_color": "#D9EAF7", "border": 1})
        warn_fmt = workbook.add_format({"bg_color": "#FCE4D6"})
        good_fmt = workbook.add_format({"bg_color": "#E2F0D9"})
        for sheet_name, ws in writer.sheets.items():
            ws.freeze_panes(1, 0)
            ws.autofilter(0, 0, 0, 30)
            ws.set_row(0, None, header_fmt)
            ws.set_column(0, 0, 22)
            ws.set_column(1, 3, 32)
            ws.set_column(4, 20, 14)
            if sheet_name in {"all_metrics_long", "fair_leaderboard", "best5_by_testset"}:
                ws.conditional_format("A2:AZ20000", {"type": "text", "criteria": "containing", "value": "source_shift_failure", "format": warn_fmt})
                ws.conditional_format("A2:AZ20000", {"type": "text", "criteria": "containing", "value": "reviewer_safe_internal_locked", "format": good_fmt})
    write_summary(all_metrics, manifest)
    print(f"[neo-xlsx] wrote {XLSX}")
    print(f"[neo-xlsx] summary {SUMMARY}")
    print(f"[neo-xlsx] metric_rows_raw={len(all_metrics)} metric_rows_dedup={len(dedup)} source_tables={len(manifest)}")


if __name__ == "__main__":
    write_workbook()
