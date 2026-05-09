#!/usr/bin/env python3
"""Create structure-preparation jobs for counterfactual MD controls."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
DESIGN = OUT / "counterfactual_design"


def exists_path(x: object) -> bool:
    if pd.isna(x):
        return False
    s = str(x)
    return bool(s and s.lower() not in {"nan", "none", "na"} and Path(s).exists())


def choose_route(row: pd.Series) -> tuple[str, str, bool, str]:
    readiness = str(row.get("readiness", ""))
    blocked = "" if pd.isna(row.get("blocked_by", "")) else str(row.get("blocked_by", ""))
    control = str(row.get("control_type", ""))
    complex_kind = str(row.get("complex_kind", ""))
    seq_status = str(row.get("sequence_status", ""))
    has_template = exists_path(row.get("template_pdb_path", ""))

    if readiness.startswith("ready_TCR_pMHC_positive_control") and has_template:
        return "openmm_prepare_existing_extracted_positive_control", "", True, "run_openmm_pilot_after_pdbfixer"
    if readiness.startswith("ready_existing_tcr_pmhc_template") and has_template:
        return "openmm_prepare_existing_candidate_template", "", True, "run_openmm_pilot_after_pdbfixer"
    if readiness.startswith("ready_template_or_structure_generation_needed") and has_template:
        return "openmm_prepare_existing_candidate_template", "", True, "run_openmm_pilot_after_pdbfixer"
    if control == "wildtype":
        if "manual_confirmation_required" in readiness or "requires_manual_confirmation" in seq_status:
            return "confirm_wt_then_thread_or_predict_structure", "confirm_wt_sequence_before_simulation", False, "manual_confirm_then_build_structure"
        return "blocked_missing_wt_sequence", blocked or "missing_wt_sequence", False, "curate_wt_sequence"
    if control == "anchor_preserved_scrambled_decoy":
        if has_template:
            return "thread_decoy_on_candidate_template_then_pdbfixer", "prepare_decoy_structure_from_template", False, "thread_peptide_then_minimize"
        return "predict_decoy_pmhc_structure", blocked or "prepare_decoy_structure", False, "run_structure_prediction_first"
    if control == "same_or_similar_hla_positive_control":
        if readiness.startswith("ready_existing_positive_control_pdb") and has_template:
            return "fetch_extract_positive_control_then_openmm_prepare", "", True, "run_openmm_pilot_after_pdbfixer"
        if "structure_generation_required" in readiness:
            return "predict_or_fetch_positive_control_structure", blocked or "prepare_positive_control_structure", False, "fetch_or_predict_control_structure"
    if has_template:
        return "openmm_prepare_existing_template", blocked, not bool(blocked), "run_openmm_pilot_after_pdbfixer"
    return "blocked_missing_template_or_structure", blocked or "missing_template_or_structure", False, "prepare_structure_first"


def main() -> None:
    manifest_path = DESIGN / "counterfactual_md_batch_manifest.tsv"
    if not manifest_path.exists():
        raise SystemExit("Run 14_prepare_md_counterfactual_batch.py first")
    manifest = pd.read_csv(manifest_path, sep="\t")
    rows = []
    for i, r in manifest.iterrows():
        route, blocker, openmm_ready, next_action = choose_route(r)
        input_pdb = "" if pd.isna(r.get("template_pdb_path", "")) else str(r.get("template_pdb_path", ""))
        if not exists_path(input_pdb):
            input_pdb = ""
        rows.append(
            {
                "prep_id": f"MDPREP_{i:04d}",
                "batch_id": r.get("batch_id", ""),
                "priority_rank": r.get("priority_rank", ""),
                "row_id": r.get("row_id", ""),
                "peptide": r.get("peptide", ""),
                "hla_4digit": r.get("hla_4digit", ""),
                "control_type": r.get("control_type", ""),
                "complex_kind": r.get("complex_kind", ""),
                "sequence": r.get("sequence", ""),
                "sequence_status": r.get("sequence_status", ""),
                "input_template_pdb": input_pdb,
                "structure_route": route,
                "openmm_input_ready": bool(openmm_ready),
                "blocked_by": blocker,
                "next_action": next_action,
                "replicates": r.get("replicates", 3),
                "ns_per_replicate": r.get("ns_per_replicate", 10),
                "claim_status": "structure_prep_manifest_only_not_immunogenicity_evidence",
            }
        )
    jobs = pd.DataFrame(rows)
    jobs.to_csv(DESIGN / "structure_prep_job_manifest.tsv", sep="\t", index=False)
    ready = jobs[jobs["openmm_input_ready"]].copy()
    blocked = jobs[~jobs["openmm_input_ready"]].copy()
    ready.to_csv(DESIGN / "openmm_ready_structure_jobs.tsv", sep="\t", index=False)
    blocked.to_csv(DESIGN / "blocked_structure_jobs.tsv", sep="\t", index=False)

    command_rows = []
    for _, r in ready.iterrows():
        outdir = f"{r['batch_id']}_rep{{REP}}"
        command_rows.append(
            {
                "prep_id": r["prep_id"],
                "batch_id": r["batch_id"],
                "input_pdb": r["input_template_pdb"],
                "command_template": (
                    "python run_openmm_pilot.py "
                    f"--input {r['input_template_pdb']} "
                    f"--outdir {outdir} "
                    "--mode explicit --platform CUDA --ns 10 --report-steps 50000 "
                    "--minimize-iterations 1000 --timestep-fs 2.0 --temperature-k 300"
                ),
            }
        )
    commands = pd.DataFrame(command_rows)
    commands.to_csv(DESIGN / "openmm_ready_command_templates.tsv", sep="\t", index=False)

    report = [
        "# Counterfactual Structure Preparation Jobs",
        "",
        "## Boundary",
        "",
        "This file prepares future MD inputs. It does not launch new simulations and does not establish immunogenicity.",
        "",
        "## Summary",
        "",
        f"- total prep jobs: {len(jobs)}",
        f"- OpenMM-ready existing/extracted template jobs: {len(ready)}",
        f"- blocked or structure-generation-required jobs: {len(blocked)}",
        "",
        "## Ready Jobs",
        "",
        ready[[c for c in ["prep_id", "batch_id", "peptide", "hla_4digit", "control_type", "complex_kind", "sequence", "structure_route", "input_template_pdb"] if c in ready.columns]].head(30).to_markdown(index=False)
        if not ready.empty
        else "No ready jobs.",
        "",
        "## Blocked Job Counts",
        "",
        blocked.groupby(["structure_route", "blocked_by"], dropna=False).size().reset_index(name="n").to_markdown(index=False)
        if not blocked.empty
        else "No blocked jobs.",
    ]
    (DESIGN / "structure_prep_job_report.md").write_text("\n".join(report) + "\n")
    print(f"[md-structure-prep] jobs={len(jobs)} openmm_ready={len(ready)} blocked={len(blocked)}")
    print(DESIGN / "structure_prep_job_manifest.tsv")


if __name__ == "__main__":
    main()
