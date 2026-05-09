#!/usr/bin/env python3
"""Write reviewer-safe CLEAN-NeoBench and BAR-Neo reports."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, value_counts_md


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="Output directory")
    return parser.parse_args()


def load(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def write(path: Path, text: str) -> None:
    path.write_text(text.rstrip() + "\n")


def flag_count_table(flags: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if flags.empty:
        return pd.DataFrame()
    for col in flags.drop(columns=["candidate_id"], errors="ignore").columns:
        vc = flags[col].fillna("NA").astype(str).value_counts(dropna=False).head(8)
        for value, count in vc.items():
            rows.append({"flag": col, "value": value, "count": int(count)})
    return pd.DataFrame(rows)


def missing_or_empty_count(df: pd.DataFrame, col: str) -> int:
    s = df[col]
    return int((s.isna() | s.fillna("").astype(str).eq("")).sum())


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dir(output_root)
    master = load(output_root / "clean_neobench_master.tsv")
    scores = load(output_root / "clean_neobench_method_scores.tsv")
    flags = load(output_root / "clean_neobench_overlap_flags.tsv")
    leaderboard = load(output_root / "clean_neobench_leaderboard.tsv")
    metrics = load(output_root / "clean_neobench_split_metrics.tsv")
    barneo = load(output_root / "barneo_candidate_scores.tsv")
    reasons = load(output_root / "barneo_abstention_reasons.tsv")

    clean_board = leaderboard[leaderboard.get("method_role", "").isin(["anchor", "internal_candidate", "bounded_fallback"])] if len(leaderboard) else pd.DataFrame()
    public_board = leaderboard[leaderboard.get("method_role", "").eq("caveated_public_comparator")] if len(leaderboard) else pd.DataFrame()
    source_board = metrics[metrics.get("split_contract", "").eq("source_heldout")].sort_values(["split_group", "AUPRC"], ascending=[True, False]) if len(metrics) else pd.DataFrame()
    hla_board = metrics[metrics.get("split_contract", "").isin(["hla_heldout", "korean_hla_focus"])].sort_values(["split_contract", "split_group", "AUPRC"], ascending=[True, True, False]) if len(metrics) else pd.DataFrame()
    low_prev = metrics[metrics.get("split_contract", "").eq("low_prevalence_heldout")].sort_values("AUPRC", ascending=False) if len(metrics) else pd.DataFrame()
    calibration = metrics.sort_values("calibration_ece", ascending=True) if len(metrics) and "calibration_ece" in metrics else pd.DataFrame()

    method_card = f"""
# CLEAN-NeoBench Method Card

## What This Framework Does

CLEAN-NeoBench is a leakage-aware AI neoantigen predictor benchmarking framework. It normalizes candidate identity, HLA metadata, labels, public/internal predictor outputs, overlap flags, split contracts, and reviewer-facing metrics into one rerunnable benchmark surface.

## What It Does Not Claim

- No clinical vaccine selection.
- No new SOTA predictor claim.
- No external validation proven.
- No standalone QK/quantum superiority claim.
- No statement that public pretrained tools are clean baselines without row-level training-corpus overlap audit.
- No unified Class I/Class II predictor.

## Benchmark Contracts

- Exact peptide-HLA holdout.
- Near peptide cluster holdout using existing clusters or k-mer/identity fallback.
- Source protein/window holdout when metadata exist.
- Source/study heldout, including CEDAR, NEPdb, and TESLA-like stress groups when present.
- HLA allele and HLA supertype heldout.
- Patient/study heldout when metadata exist.
- Low-prevalence heldout for TESLA-like settings.
- Korean-HLA focus board for A*24:02, A*11:01, A*02:01, B*15:01, B*40:01, C*01:02, C*03:03, C*07:02.

## Method Zoo

{value_counts_md(scores["method_role"] if len(scores) else pd.Series(dtype=str))}

## Leakage Audit Design

Overlap flags include exact peptide, exact peptide-HLA, near peptide cluster, source protein/window, study, patient, and public-tool training overlap status. Public pretrained methods are caveated unless official training corpora are row-audited.

## Metrics

CLEAN-NeoBench emphasizes AUPRC, top-k precision/recall/hit, positive rank, calibration Brier/ECE, confidence coverage, risk at high confidence, and abstention rate. AUROC is reported but is not the sole decision metric.

## Limitations

Missing WT peptide, source window, patient, expression, clonality, HLA LOH, and B2M metadata reduce the strength of source-window and patient-gated claims. Missing metrics are reported as unavailable rather than fabricated.

## Claim Boundary

Allowed claims: leakage-aware benchmark, benchmark-adaptive reliability ranking, calibrated candidate prioritization with abstention, reviewer-safe comparison of internal and public AI predictors, and research triage framework.
"""
    write(output_root / "CLEAN_NEOBENCH_METHOD_CARD.md", method_card)

    leaderboard_md = f"""
# CLEAN-NeoBench Leaderboard

## Overall Leaderboard

{dataframe_to_markdown(leaderboard.head(25) if len(leaderboard) else leaderboard)}

## Clean Internal Leaderboard

{dataframe_to_markdown(clean_board.head(25))}

## Caveated Public Comparator Leaderboard

{dataframe_to_markdown(public_board.head(25))}

## Source-Heldout Leaderboard

{dataframe_to_markdown(source_board.head(40))}

## HLA-Heldout and Korean-HLA Focus Board

{dataframe_to_markdown(hla_board.head(50))}

## Low-Prevalence Board

{dataframe_to_markdown(low_prev.head(25))}

## Calibration Board

{dataframe_to_markdown(calibration[["split_contract", "split_group", "method_name", "method_role", "n_total", "n_pos", "AUPRC", "calibration_ece", "calibration_brier"]].head(30) if len(calibration) else calibration)}

## Abstention Board

{dataframe_to_markdown(barneo[["candidate_id", "barneo_score", "patient_gated_score", "confidence_score", "confidence_bin", "abstain", "abstention_reason_primary"]].head(30) if len(barneo) else barneo)}

## Interpretation Rule

Prefer methods that perform across strict/no-overlap, source-heldout, HLA-heldout, low-prevalence, and calibration boards. Do not promote a method that only wins in one split or hides source-heldout collapse.
"""
    write(output_root / "CLEAN_NEOBENCH_LEADERBOARD.md", leaderboard_md)

    split_audit = f"""
# CLEAN-NeoBench Split Audit

## Candidate Counts

- Candidates: {len(master)}
- Labeled candidates: {int(master["label"].notna().sum()) if len(master) else 0}
- Method score rows: {len(scores)}
- Metric rows: {len(metrics)}

## Source Distribution

{value_counts_md(master["source_name"] if len(master) else pd.Series(dtype=str))}

## Positive Counts by Source

{dataframe_to_markdown(master.groupby("source_name", dropna=False)["label"].agg(["count", "sum", "mean"]).reset_index().sort_values("count", ascending=False).head(30) if len(master) else pd.DataFrame())}

## Overlap Flag Summary

{dataframe_to_markdown(flag_count_table(flags).head(80))}

## Leakage Risk Distribution

{value_counts_md(master["leakage_risk_level"] if len(master) else pd.Series(dtype=str))}

## Missing Metadata Summary

{dataframe_to_markdown(pd.DataFrame({"column": master.columns, "missing_or_empty": [missing_or_empty_count(master, c) for c in master.columns]}).sort_values("missing_or_empty", ascending=False).head(40) if len(master) else pd.DataFrame())}

## Split Contract Rows

{value_counts_md(metrics["split_contract"] if len(metrics) else pd.Series(dtype=str))}
"""
    write(output_root / "CLEAN_NEOBENCH_SPLIT_AUDIT.md", split_audit)

    public = scores[scores["uses_public_pretraining"].astype(bool)] if len(scores) else pd.DataFrame()
    public_caveat = f"""
# CLEAN-NeoBench Public Tool Caveat

## Included Public or Public-Pretrained Tools

{dataframe_to_markdown(public[["method_name", "method_family", "method_role", "uses_public_pretraining", "training_overlap_audited", "clean_comparator_allowed", "caveat"]].drop_duplicates().sort_values("method_name").head(80) if len(public) else public)}

## Caveat Rule

Public pretrained tools are caveated comparators until their official training corpora are row-audited against benchmark candidates. A no-overlap flag against our local master table is not the same as a no-overlap guarantee against MHCflurry, NetMHCpan, BigMHC, PRIME, MixMHCpred, NetMHCstabpan, or other public training data.

## Allowed Comparisons

- Public tools may be shown as caveated upper-bound or field-context comparators.
- Internal models may be compared against them only with the training-overlap caveat visible.
- Clean internal anchor status is reserved for locally audited methods such as Structure_LR.

## Forbidden Claims

- Public pretrained tools are clean external baselines without row-level training-corpus audit.
- A public tool win proves external validation.
- A public score can be used as a clean training feature before overlap audit.
"""
    write(output_root / "CLEAN_NEOBENCH_PUBLIC_TOOL_CAVEAT.md", public_caveat)

    barneo_card = f"""
# BAR-Neo Method Card

## Algorithm

BAR-Neo is a Benchmark-Adaptive Reliability Neoantigen Ranking layer. It is not a giant deep model. It combines candidate intrinsic features, method rank percentiles, predictor disagreement, overlap/leakage flags, HLA/source support, historical source-heldout behavior, and transparent patient gates.

## Feature Groups

- Candidate intrinsic features: peptide length, hydrophobic/aromatic/charge proxy, cysteine count, glycine/proline count, anchor residue proxies.
- HLA features: HLA gene and supertype one-hot, allele support count, Korean-HLA focus flags.
- Method-derived features: score percentiles, number of methods available, clean internal score, caveated public score, public-vs-internal disagreement.
- Reliability features: exact/near/source/study/patient overlap flags, leakage risk level, source prevalence, source-heldout historical top-k behavior.
- Uncertainty features: predictor disagreement, low method availability, underrepresented HLA, uncertainty-only model support where available.

## Calibration Logic

Method scores are direction-normalized and calibrated with fold-safe isotonic/logistic calibration when feasible. Small methods fall back to rank-percentile calibration and are marked through abstention logic when unstable.

## Abstention Logic

BAR-Neo abstains or lowers confidence for high leakage risk, underrepresented HLA, source-shift risk, high method disagreement, caveated-public-only evidence, MHC-II rows in a Class-I benchmark, failed presentation gates, low-prevalence source collapse, and sparse method coverage.

## Patient-Gated Scoring Logic

`patient_gated_score = barneo_score * disease_gate * presentation_gate * antigen_gate * immune_context_gate * model_confidence_gate`

PAAD high-priority context is resected/MRD/low-burden disease. THCA high-priority research context is ATC, PDTC, progressive radioiodine-refractory DTC, or high-risk recurrence. Unknown patient metadata are marked as uncertainty rather than inferred.

## Limitations

BAR-Neo currently reflects available public/local benchmark labels and sparse patient metadata. It is a research triage and reliability layer, not a clinical treatment selector.
"""
    write(output_root / "BAR_NEO_METHOD_CARD.md", barneo_card)

    reason_counts = reasons["abstention_reason"].replace("", "no_abstention_reason").value_counts().reset_index() if len(reasons) else pd.DataFrame()
    reason_counts.columns = ["abstention_reason", "count"] if len(reason_counts) else []
    public_weak = barneo[barneo["abstention_reason_all"].fillna("").str.contains("public", case=False, na=False)].head(20) if len(barneo) else pd.DataFrame()
    source_weak = barneo[barneo["abstention_reason_all"].fillna("").str.contains("source", case=False, na=False)].head(20) if len(barneo) else pd.DataFrame()
    abstention_report = f"""
# BAR-Neo Abstention Report

## Abstention Counts

- Candidates scored: {len(barneo)}
- Abstain true: {int(barneo["abstain"].sum()) if len(barneo) else 0}
- Abstain false: {int((~barneo["abstain"].astype(bool)).sum()) if len(barneo) else 0}

## Top Abstention Reasons

{dataframe_to_markdown(reason_counts.head(30))}

## High Score / Low Confidence Cases

{dataframe_to_markdown(barneo.sort_values(["barneo_score", "confidence_score"], ascending=[False, True])[["candidate_id", "barneo_score", "patient_gated_score", "confidence_score", "confidence_bin", "abstain", "abstention_reason_primary"]].head(20) if len(barneo) else barneo)}

## Public-Strong / Internal-Weak Caveat Examples

{dataframe_to_markdown(public_weak[["candidate_id", "barneo_score", "confidence_score", "abstention_reason_primary", "abstention_reason_all"]] if len(public_weak) else public_weak)}

## Source-Shift Abstention Examples

{dataframe_to_markdown(source_weak[["candidate_id", "barneo_score", "confidence_score", "abstention_reason_primary", "abstention_reason_all"]] if len(source_weak) else source_weak)}
"""
    write(output_root / "BAR_NEO_ABSTENTION_REPORT.md", abstention_report)

    update_manifest(
        output_root,
        "write_clean_neobench_reports",
        {
            "reports": [
                "CLEAN_NEOBENCH_METHOD_CARD.md",
                "CLEAN_NEOBENCH_LEADERBOARD.md",
                "CLEAN_NEOBENCH_SPLIT_AUDIT.md",
                "CLEAN_NEOBENCH_PUBLIC_TOOL_CAVEAT.md",
                "BAR_NEO_METHOD_CARD.md",
                "BAR_NEO_ABSTENTION_REPORT.md",
            ],
            "warnings": [],
        },
    )
    manifest_path = output_root / "run_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    else:
        manifest = {}
    manifest["output_files"] = sorted(p.name for p in output_root.iterdir() if p.is_file())
    manifest["summary"] = {
        "n_candidates": int(len(master)),
        "n_methods": int(scores["method_name"].nunique()) if len(scores) else 0,
        "n_metric_rows": int(len(metrics)),
        "n_barneo_scores": int(len(barneo)),
    }
    if len(master):
        manifest["missing_metadata"] = {
            col: missing_or_empty_count(master, col)
            for col in [
                "wt_peptide",
                "patient_id",
                "source_protein_window",
                "expression_tpm",
                "mutant_expression",
                "vaf",
                "clonality",
                "hla_loh",
                "b2m_status",
                "antigen_processing_status",
                "immune_context_score",
                "tls_score",
                "ifng_score",
                "cytolytic_score",
            ]
            if col in master.columns
        }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"[clean-neobench-reports] wrote reports to {output_root}")


if __name__ == "__main__":
    main()
