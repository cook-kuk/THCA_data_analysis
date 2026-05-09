#!/usr/bin/env python3
"""Parse local TCR-pMHC structure outputs into QC and feature tables."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import OUT, safe_parquet


TCR_OUT = OUT / "tcr_extension"
JOB_OUT = TCR_OUT / "structure_jobs"
STRUCT_OUT = TCR_OUT / "structure_outputs"
LOCAL_AF_PILOT = Path("/data/neoantigen_vaccine_hub/experiments/alphafold_pilot")


def clean(value: object) -> str:
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value or "").strip()
    return "" if text.lower() in {"", "na", "nan", "none", "<na>"} else text


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def mean_plddt(scores: dict) -> float:
    vals = scores.get("plddt", [])
    return float(np.mean(vals)) if vals else np.nan


def parse_local_af2_pilots() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for score_path in sorted(LOCAL_AF_PILOT.glob("results*/**/*scores*.json")):
        scores = read_json(score_path)
        stem = score_path.name.split("_scores_", 1)[0]
        pdbs = list(score_path.parent.glob(f"{stem}*unrelaxed*.pdb"))
        pae = list(score_path.parent.glob(f"{stem}*predicted_aligned_error*.json"))
        rows.append(
            {
                "job_id": f"LOCAL_AF2_PILOT_{len(rows):03d}",
                "neo_row_id": "",
                "mode": "local_existing_pilot",
                "tool": "AlphaFold-Multimer_colabfold_batch",
                "output_found": True,
                "pdb_path": str(pdbs[0]) if pdbs else "",
                "score_json_path": str(score_path),
                "pae_json_path": str(pae[0]) if pae else "",
                "plddt_mean": mean_plddt(scores),
                "iptm": float(scores.get("iptm", np.nan)) if "iptm" in scores else np.nan,
                "ptm": float(scores.get("ptm", np.nan)) if "ptm" in scores else np.nan,
                "max_pae": float(scores.get("max_pae", np.nan)) if "max_pae" in scores else np.nan,
                "interface_contact_count": np.nan,
                "peptide_tcr_contact_count": np.nan,
                "peptide_mhc_contact_count": np.nan,
                "cdr3a_peptide_contact_count": np.nan,
                "cdr3b_peptide_contact_count": np.nan,
                "crossing_angle": np.nan,
                "peptide_bulge_proxy": np.nan,
                "mutation_residue_contact_flag": np.nan,
                "mutant_wt_interface_delta": np.nan,
                "wt_cross_reactivity_risk": np.nan,
                "claim_status": "local_pmhc_or_incomplete_pilot_not_tcr_claim",
                "qc_warning": "Existing local AF2 pilot was parsed for parser validation; do not claim TCR-pMHC interface unless full TCR output is present.",
            }
        )
    return rows


def main() -> None:
    STRUCT_OUT.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(JOB_OUT / "tcr_structure_job_manifest.tsv", sep="\t", low_memory=False)
    rows: list[dict[str, object]] = []
    for rec in manifest.itertuples(index=False):
        row = rec._asdict()
        rows.append(
            {
                "job_id": clean(row.get("job_id")),
                "neo_row_id": clean(row.get("neo_row_id")),
                "mode": clean(row.get("mode")),
                "tool": clean(row.get("tool")),
                "priority": clean(row.get("priority")),
                "output_found": False,
                "pdb_path": "",
                "score_json_path": "",
                "pae_json_path": "",
                "plddt_mean": np.nan,
                "iptm": np.nan,
                "ptm": np.nan,
                "max_pae": np.nan,
                "interface_contact_count": np.nan,
                "peptide_tcr_contact_count": np.nan,
                "peptide_mhc_contact_count": np.nan,
                "cdr3a_peptide_contact_count": np.nan,
                "cdr3b_peptide_contact_count": np.nan,
                "crossing_angle": np.nan,
                "peptide_bulge_proxy": np.nan,
                "mutation_residue_contact_flag": np.nan,
                "mutant_wt_interface_delta": np.nan,
                "wt_cross_reactivity_risk": np.nan,
                "claim_status": clean(row.get("claim_status")) or "no_structure_output",
                "qc_warning": "No parsed structure output found for this manifest job.",
            }
        )
    rows.extend(parse_local_af2_pilots())

    features = pd.DataFrame(rows)
    safe_parquet(features, STRUCT_OUT / "tcr_structure_features.parquet")
    qc_cols = [
        "job_id",
        "neo_row_id",
        "mode",
        "tool",
        "priority",
        "output_found",
        "pdb_path",
        "score_json_path",
        "plddt_mean",
        "iptm",
        "ptm",
        "max_pae",
        "claim_status",
        "qc_warning",
    ]
    features[[c for c in qc_cols if c in features.columns]].to_csv(STRUCT_OUT / "tcr_structure_qc.tsv", sep="\t", index=False, na_rep="NA")

    n_found = int(features["output_found"].fillna(False).astype(bool).sum())
    n_jobs = int(len(manifest))
    lines = [
        "# CROSS-Neo-TCR Structure Output Parse Report",
        "",
        f"Manifest jobs checked: {n_jobs:,}",
        f"Parsed structure-like outputs: {n_found:,}",
        "",
        "## Interpretation",
        "",
        "- No manifest TCR-pMHC job currently has a parsed local output.",
        "- Existing AlphaFold/ColabFold pilot files under `/data/neoantigen_vaccine_hub/experiments/alphafold_pilot` were parsed only to validate the parser and confidence-field extraction.",
        "- Interface/contact/delta features remain missing until actual TCR-pMHC mutant/WT jobs are generated.",
    ]
    (STRUCT_OUT / "tcr_structure_parse_report.md").write_text("\n".join(lines) + "\n")
    with (TCR_OUT / "tcr_structure_pipeline_report.md").open("a") as fh:
        fh.write("\n## Parse Status\n\n")
        fh.write(f"- Manifest jobs checked: {n_jobs:,}\n")
        fh.write(f"- Parsed structure-like outputs: {n_found:,}\n")
        fh.write("- Interface features are not claimable until actual TCR-pMHC outputs exist.\n")
    print(f"[tcr-structure-parse] manifest_jobs={n_jobs} parsed_outputs={n_found}")


if __name__ == "__main__":
    main()
