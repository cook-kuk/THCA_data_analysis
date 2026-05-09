#!/usr/bin/env python3
"""Prepare CROSS-Neo-TCR structure job manifests.

The manifest is intentionally conservative: rows without full TCR/MHC sequences
are retained and flagged as blocked instead of being discarded.
"""

from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path

import pandas as pd

from common import OUT


TCR_OUT = OUT / "tcr_extension"
JOB_OUT = TCR_OUT / "structure_jobs"


def clean(value: object) -> str:
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value or "").strip()
    return "" if text.lower() in {"", "na", "nan", "none", "<na>"} else text


def available_tool_status() -> dict[str, str]:
    checks = {
        "TCRdock": shutil.which("TCRdock") or shutil.which("tcrdock"),
        "TCRmodel2": shutil.which("TCRmodel2") or shutil.which("tcrmodel2"),
        "AlphaFold-Multimer_colabfold_batch": shutil.which("colabfold_batch"),
        "AlphaFold_python_module": "available" if importlib.util.find_spec("alphafold") else "",
        "AlphaFold3": shutil.which("alphafold3") or ("available" if importlib.util.find_spec("alphafold3") else ""),
        "Boltz": shutil.which("boltz") or ("available" if importlib.util.find_spec("boltz") else ""),
        "Chai-1": shutil.which("chai-lab") or ("available" if importlib.util.find_spec("chai_lab") else ""),
        "tFold-TCR": shutil.which("tfold") or ("available" if importlib.util.find_spec("tfold") else ""),
        "ImmuneBuilder": "available" if importlib.util.find_spec("ImmuneBuilder") else "",
    }
    return {k: str(v) if v else "not_found" for k, v in checks.items()}


def read_optional(path: Path) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path, sep="\t", low_memory=False)
    return pd.DataFrame()


def row_set(path: Path) -> set[str]:
    df = read_optional(path)
    return set(df["row_id"].astype(str)) if not df.empty and "row_id" in df.columns else set()


def choose_tool(has_full_tcr: bool, has_mhc_sequence: bool, status: dict[str, str]) -> tuple[str, str, str]:
    if has_full_tcr and has_mhc_sequence and status["AlphaFold-Multimer_colabfold_batch"] != "not_found":
        return "AlphaFold-Multimer_colabfold_batch", "gpu_hours", "runnable_local_colabfold_multimer"
    if has_full_tcr and has_mhc_sequence and status["Boltz"] != "not_found":
        return "Boltz_or_Boltz2", "gpu_hours", "runnable_if_boltz_configured"
    if has_full_tcr and has_mhc_sequence and status["Chai-1"] != "not_found":
        return "Chai-1", "gpu_hours", "runnable_if_chai_configured"
    if has_full_tcr and status["TCRdock"] != "not_found":
        return "TCRdock", "minutes_to_gpu_hours", "runnable_if_tcrdock_configured"
    if has_full_tcr and status["TCRmodel2"] != "not_found":
        return "TCRmodel2", "minutes_to_gpu_hours", "runnable_if_tcrmodel2_configured"
    if not has_full_tcr:
        return "blocked_missing_full_tcr_sequence", "blocked", "diagnostic_manifest_only_requires_full_tcr"
    if not has_mhc_sequence:
        return "blocked_missing_mhc_sequence", "blocked", "diagnostic_manifest_only_requires_mhc_sequence"
    return "blocked_no_supported_structure_tool", "blocked", "diagnostic_manifest_only_no_local_tool"


def priority_for(row: pd.Series, rescue_rows: set[str], fp_rows: set[str], candidate_rows: set[str]) -> str:
    rid = str(row.get("neo_row_id", row.get("row_id", "")))
    label = pd.to_numeric(row.get("label_binary", 0), errors="coerce")
    paired = bool(row.get("tcr_paired_available", False))
    exact = clean(row.get("link_category", row.get("tcr_link_category", ""))) == "exact_tcr_pmhc_match"
    if paired and exact and (label == 1 or rid in rescue_rows or rid in candidate_rows or rid in fp_rows):
        return "P0"
    if paired and exact:
        return "P1"
    if paired:
        return "P1"
    return "P2"


def main() -> None:
    JOB_OUT.mkdir(parents=True, exist_ok=True)
    status = available_tool_status()
    linked = pd.read_csv(TCR_OUT / "tcr_neo_linked_registry.tsv", sep="\t", low_memory=False)
    evidence = pd.read_csv(TCR_OUT / "tcr_neo_linkage_evidence_examples.tsv", sep="\t", low_memory=False)

    rescue_rows = row_set(TCR_OUT / "tcr_rescue_cases.tsv")
    fp_rows = row_set(TCR_OUT / "tcr_false_positive_audit.tsv")
    candidate_rows = row_set(TCR_OUT / "tcr_wetlab_candidate_prioritization.tsv")
    top_candidates = read_optional(TCR_OUT / "tcr_wetlab_candidate_prioritization.tsv")
    top_candidate_rank: dict[str, int] = {}
    if not top_candidates.empty:
        for i, rid in enumerate(top_candidates["row_id"].astype(str).head(500), start=1):
            top_candidate_rank[rid] = i

    neo_cols = [
        "row_id",
        "source_dataset",
        "peptide",
        "mutant_peptide",
        "wildtype_peptide",
        "hla_4digit",
        "label_binary",
        "tcr_link_category",
        "tcr_evidence_count",
        "paired_tcr_evidence_count",
        "cancer_context_evidence_count",
        "pathogen_context_evidence_count",
    ]
    neo_cols = [c for c in neo_cols if c in linked.columns]
    jobs_src = evidence.merge(linked[neo_cols], left_on="neo_row_id", right_on="row_id", how="left", suffixes=("", "_neo"))
    jobs_src = jobs_src[jobs_src["tcr_paired_available"].fillna(False).astype(bool)].copy()
    jobs_src["cdr3_alpha"] = jobs_src["tcr_cdr3_alpha"].map(clean)
    jobs_src["cdr3_beta"] = jobs_src["tcr_cdr3_beta"].map(clean)
    jobs_src = jobs_src[jobs_src["cdr3_alpha"].ne("") & jobs_src["cdr3_beta"].ne("")].copy()
    jobs_src = jobs_src.drop_duplicates(["neo_row_id", "tcr_row_id", "neo_peptide", "neo_hla"])

    rows: list[dict[str, object]] = []
    seen: set[tuple[str, str, str]] = set()
    for rec in jobs_src.itertuples(index=False):
        row = pd.Series(rec._asdict())
        rid = clean(row.get("neo_row_id"))
        priority = priority_for(row, rescue_rows, fp_rows, candidate_rows)
        mode_list = ["mutant"]
        if clean(row.get("wildtype_peptide")):
            mode_list.append("wildtype")
        has_full_tcr = bool(clean(row.get("tcr_alpha_full_sequence")) and clean(row.get("tcr_beta_full_sequence")))
        has_mhc_sequence = bool(clean(row.get("mhc_sequence")))
        tool, runtime, claim_status = choose_tool(has_full_tcr, has_mhc_sequence, status)
        for mode in mode_list:
            key = (rid, clean(row.get("tcr_row_id")), mode)
            if key in seen:
                continue
            seen.add(key)
            peptide = clean(row.get("wildtype_peptide")) if mode == "wildtype" else clean(row.get("mutant_peptide")) or clean(row.get("neo_peptide"))
            rows.append(
                {
                    "job_id": f"TCRSTRUCT_{len(rows):06d}",
                    "neo_row_id": rid,
                    "tcr_row_id": clean(row.get("tcr_row_id")),
                    "source_dataset": clean(row.get("source_dataset")) or clean(row.get("tcr_source_dataset")),
                    "peptide": peptide,
                    "hla_4digit": clean(row.get("hla_4digit")) or clean(row.get("neo_hla")),
                    "mhc_sequence": clean(row.get("mhc_sequence")),
                    "cdr3_alpha": clean(row.get("cdr3_alpha")),
                    "cdr3_beta": clean(row.get("cdr3_beta")),
                    "tcr_alpha_full_sequence": clean(row.get("tcr_alpha_full_sequence")),
                    "tcr_beta_full_sequence": clean(row.get("tcr_beta_full_sequence")),
                    "mutant_peptide": clean(row.get("mutant_peptide")) or clean(row.get("neo_peptide")),
                    "wildtype_peptide": clean(row.get("wildtype_peptide")),
                    "mode": mode,
                    "tool": tool,
                    "priority": priority,
                    "expected_runtime_class": runtime,
                    "claim_status": claim_status,
                    "label_binary": row.get("label_binary", ""),
                    "link_category": clean(row.get("link_category")),
                    "candidate_rank": top_candidate_rank.get(rid, ""),
                    "is_rescue_case": rid in rescue_rows,
                    "is_false_positive_audit_case": rid in fp_rows,
                    "tcr_source_dataset": clean(row.get("tcr_source_dataset")),
                    "tcr_antigen_source": clean(row.get("tcr_antigen_source")),
                    "tcr_binding_label_binary": row.get("tcr_binding_label_binary", ""),
                }
            )

    jobs = pd.DataFrame(rows)
    if not jobs.empty:
        priority_order = {"P0": 0, "P1": 1, "P2": 2}
        jobs["_p"] = jobs["priority"].map(priority_order).fillna(9)
        jobs = jobs.sort_values(["_p", "candidate_rank", "neo_row_id", "tcr_row_id", "mode"]).drop(columns=["_p"])
    jobs.to_csv(JOB_OUT / "tcr_structure_job_manifest.tsv", sep="\t", index=False, na_rep="NA")
    jobs.head(500).to_csv(JOB_OUT / "tcr_structure_job_manifest_preview.tsv", sep="\t", index=False, na_rep="NA")

    tool_df = pd.DataFrame([{"tool": k, "local_status": v} for k, v in status.items()])
    tool_df.to_csv(JOB_OUT / "tcr_structure_tool_availability.tsv", sep="\t", index=False)
    summary = (
        jobs.groupby(["priority", "tool", "claim_status"], dropna=False)
        .size()
        .reset_index(name="n_jobs")
        .sort_values(["priority", "n_jobs"], ascending=[True, False])
        if not jobs.empty
        else pd.DataFrame(columns=["priority", "tool", "claim_status", "n_jobs"])
    )
    summary.to_csv(JOB_OUT / "tcr_structure_job_summary.tsv", sep="\t", index=False, na_rep="NA")

    lines = [
        "# CROSS-Neo-TCR Structure Job Manifest Report",
        "",
        f"Candidate paired TCR-pMHC jobs emitted: {len(jobs):,}",
        "",
        "## Tool Availability",
        "",
        "| Tool | Local status |",
        "|---|---|",
    ]
    lines.extend([f"| {k} | {v} |" for k, v in status.items()])
    lines.extend(
        [
            "",
            "## Priority Summary",
            "",
            "| Priority | Tool | Claim status | n jobs |",
            "|---|---|---|---:|",
        ]
    )
    for r in summary.itertuples(index=False):
        lines.append(f"| {r.priority} | {r.tool} | {r.claim_status} | {int(r.n_jobs)} |")
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "- Jobs are retained even when full TCR or MHC sequences are missing.",
            "- CDR3-only rows are useful for prioritization but are not sufficient for TCR-pMHC structure prediction with most tools.",
            "- Local runnable path currently favors ColabFold/AlphaFold-Multimer only when full chain sequences are available.",
        ]
    )
    (TCR_OUT / "tcr_structure_pipeline_report.md").write_text("\n".join(lines) + "\n")
    print(f"[tcr-structure-jobs] jobs={len(jobs)} p0={(jobs['priority'].eq('P0').sum() if not jobs.empty else 0)}")


if __name__ == "__main__":
    main()
