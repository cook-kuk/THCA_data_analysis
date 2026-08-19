#!/usr/bin/env python3
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from neoimmune_common import ensure_run_dir, safe_read_table, write_md, write_tsv


def annotate(row: pd.Series) -> list[str]:
    reasons = []
    expr = pd.to_numeric(row.get("expression_tpm"), errors="coerce")
    vaf = pd.to_numeric(row.get("vaf"), errors="coerce")
    if pd.notna(expr) and expr < 1:
        reasons.append("low expression")
    if pd.notna(vaf) and vaf < 0.05:
        reasons.append("low VAF/subclonal risk")
    if str(row.get("hla_loh_status", "")).lower() in {"true", "loh", "loss", "1"}:
        reasons.append("HLA LOH/APM risk")
    if str(row.get("b2m_status", "")).lower() in {"loss", "mutated", "deficient"}:
        reasons.append("B2M/APM loss risk")
    dis = pd.to_numeric(row.get("model_disagreement_score"), errors="coerce")
    if pd.notna(dis) and dis > 0.25:
        reasons.append("model disagreement")
    cov = pd.to_numeric(row.get("abstention_coverage"), errors="coerce")
    if pd.notna(cov) and cov < 0.25:
        reasons.append("sparse method coverage")
    if str(row.get("leakage_risk_level", "")).lower().find("high") >= 0:
        reasons.append("source-study/leakage risk")
    if pd.isna(row.get("label_immunogenicity")):
        reasons.append("unlabeled candidate")
    return reasons or ["no obvious metadata-level failure flag"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    prod = safe_read_table(outdir / "production_track" / "production_stack_predictions.tsv")
    clean = safe_read_table(outdir / "clean_track" / "clean_science_predictions.tsv")
    d = canon.merge(prod[["candidate_id", "production_stack_score", "model_disagreement_score", "abstention_coverage"]], on="candidate_id", how="left")
    d = d.merge(clean[["candidate_id", "clean_science_score"]], on="candidate_id", how="left")
    d["label_immunogenicity"] = pd.to_numeric(d["label_immunogenicity"], errors="coerce")
    d["production_stack_score"] = pd.to_numeric(d["production_stack_score"], errors="coerce").fillna(0.5)
    d["clean_science_score"] = pd.to_numeric(d["clean_science_score"], errors="coerce").fillna(0.5)
    d["production_rank"] = d["production_stack_score"].rank(ascending=False, method="first")
    d["clean_rank"] = d["clean_science_score"].rank(ascending=False, method="first")

    false_pos = d[(d["label_immunogenicity"] == 0) & (d["production_rank"] <= 200)].copy()
    rescued = d[(d["label_immunogenicity"] == 1) & (d["production_rank"] <= 200) & (d["clean_rank"] > d["production_rank"] + 50)].copy()
    if false_pos.empty:
        false_pos = d[(d["label_immunogenicity"] == 0)].sort_values("production_stack_score", ascending=False).head(50).copy()
    if rescued.empty:
        rescued = d[(d["label_immunogenicity"] == 1)].sort_values("production_stack_score", ascending=False).head(50).copy()

    false_pos["case_type"] = "false_positive_high_rank"
    rescued["case_type"] = "rescued_positive_or_top_positive"
    cases = pd.concat([false_pos.head(100), rescued.head(100)], ignore_index=True, sort=False)
    if not cases.empty:
        cases["audit_annotation"] = cases.apply(lambda r: "; ".join(annotate(r)), axis=1)
        cases["rescue_annotation"] = np.where(
            cases["case_type"].eq("rescued_positive_or_top_positive"),
            "candidate should be reviewed for rescue by local/TCR/structure/quantum/patient-context branch; inspect component attributions before wet-lab discussion",
            "",
        )
    cols = [
        "case_type",
        "candidate_id",
        "patient_id",
        "dataset_source",
        "peptide_mut",
        "peptide_wt",
        "hla_allele",
        "label_immunogenicity",
        "production_stack_score",
        "clean_science_score",
        "model_disagreement_score",
        "abstention_coverage",
        "expression_tpm",
        "vaf",
        "hla_loh_status",
        "b2m_status",
        "leakage_risk_level",
        "audit_annotation",
        "rescue_annotation",
    ]
    for c in cols:
        if c not in cases.columns:
            cases[c] = np.nan
    write_tsv(cases[cols], outdir / "predictions" / "failure_cases.tsv")

    md = [
        "# Failure case audit",
        "",
        "## False-positive annotations",
        "- low expression",
        "- low VAF/subclonal risk",
        "- weak or missing presentation evidence",
        "- high WT/self similarity when available",
        "- poor TCR/self-similarity branch evidence when available",
        "- HLA LOH/APM/B2M risk",
        "- model disagreement",
        "- source-study/leakage risk",
        "- external predictor conflict",
        "",
        "## Rescued positive annotations",
        "- missed by binding-only/public comparators when available",
        "- rescued by TCR/self-similarity branch",
        "- rescued by structure branch",
        "- rescued by quantum kernel branch",
        "- rescued by production stack",
        "- rescued by patient context",
        "",
        f"- Failure/rescue rows written: {len(cases):,}",
        "",
        "## Boundary",
        "This audit explains ranking behavior; it does not assert biological validation.",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "failure_case_audit.md")


if __name__ == "__main__":
    main()

