#!/usr/bin/env python3
"""Select same/similar-HLA positive controls for the next MD batch.

Positive controls are for runtime and analysis sanity checks. They are not
cancer-specificity evidence and must not be used as neoantigen ground truth.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
CROSS = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
DESIGN = OUT / "counterfactual_design"

TARGETS = {
    ("HMTEVVRHC", "HLA-A*02:01"): "A02",
    ("GADGVGKSAL", "HLA-C*08:02"): "C08",
}


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", **kwargs) if path.exists() else pd.DataFrame()


def has_value(s: pd.Series) -> pd.Series:
    return s.notna() & ~s.astype(str).str.lower().isin(["", "nan", "na", "none"])


def select_controls() -> pd.DataFrame:
    cols = [
        "row_id",
        "source_dataset",
        "species",
        "antigen_source",
        "peptide",
        "peptide_length",
        "hla_4digit",
        "hla_supertype",
        "mhc_class",
        "cdr3_alpha",
        "cdr3_beta",
        "paired_tcr_available",
        "binding_label_binary",
        "structure_pdb_id",
        "is_cancer_context",
        "is_pathogen_context",
    ]
    reg = read_tsv(CROSS / "tcr_extension/tcr_registry.tsv", usecols=lambda c: c in cols, low_memory=False)
    if reg.empty:
        return pd.DataFrame()
    reg["peptide"] = reg["peptide"].astype(str).str.upper()
    reg["has_structure"] = has_value(reg["structure_pdb_id"])
    reg["paired_tcr_available"] = reg["paired_tcr_available"].fillna(False).astype(bool)
    reg["binding_label_binary"] = pd.to_numeric(reg["binding_label_binary"], errors="coerce").fillna(0).astype(int)
    reg["peptide_length"] = pd.to_numeric(reg["peptide_length"], errors="coerce")
    reg["is_cancer_context"] = reg["is_cancer_context"].fillna(False).astype(bool)
    reg["is_pathogen_context"] = reg["is_pathogen_context"].fillna(False).astype(bool)

    rows: list[dict] = []
    for (target_peptide, target_hla), target_supertype in TARGETS.items():
        candidates = reg[
            (reg["binding_label_binary"] == 1)
            & (reg["peptide"] != target_peptide)
            & (reg["peptide_length"].between(8, 11, inclusive="both"))
            & (reg["mhc_class"].astype(str).str.upper().isin(["I", "1", "CLASS I"]))
        ].copy()
        candidates["match_class"] = np.select(
            [
                candidates["hla_4digit"].astype(str).eq(target_hla),
                candidates["hla_supertype"].astype(str).eq(target_supertype),
            ],
            ["exact_hla", "same_supertype"],
            default="other_hla",
        )
        candidates = candidates[candidates["match_class"] != "other_hla"].copy()
        if candidates.empty:
            continue
        candidates["selection_score"] = (
            candidates["match_class"].map({"exact_hla": 100, "same_supertype": 55}).fillna(0)
            + candidates["has_structure"].astype(int) * 25
            + candidates["paired_tcr_available"].astype(int) * 20
            + candidates["is_cancer_context"].astype(int) * 8
            - candidates["is_pathogen_context"].astype(int) * 3
        )
        candidates = candidates.sort_values(
            ["selection_score", "match_class", "has_structure", "paired_tcr_available", "is_cancer_context", "peptide"],
            ascending=[False, True, False, False, False, True],
        )
        seen_peptides: set[str] = set()
        rank = 0
        for _, r in candidates.iterrows():
            pep = str(r["peptide"])
            if pep in seen_peptides:
                continue
            seen_peptides.add(pep)
            rank += 1
            rows.append(
                {
                    "target_peptide": target_peptide,
                    "target_hla_4digit": target_hla,
                    "control_rank": rank,
                    "control_peptide": pep,
                    "control_hla_4digit": r.get("hla_4digit", ""),
                    "control_hla_supertype": r.get("hla_supertype", ""),
                    "match_class": r["match_class"],
                    "selection_score": r["selection_score"],
                    "control_tcr_row_id": r.get("row_id", ""),
                    "source_dataset": r.get("source_dataset", ""),
                    "antigen_source": r.get("antigen_source", ""),
                    "paired_tcr_available": bool(r.get("paired_tcr_available", False)),
                    "has_structure": bool(r.get("has_structure", False)),
                    "structure_pdb_id": r.get("structure_pdb_id", ""),
                    "is_cancer_context": bool(r.get("is_cancer_context", False)),
                    "is_pathogen_context": bool(r.get("is_pathogen_context", False)),
                    "claim_status": "positive_control_only_not_cancer_specificity_evidence",
                }
            )
            if rank >= 10:
                break
    out = pd.DataFrame(rows)
    return out


def update_batch_manifest(controls: pd.DataFrame) -> pd.DataFrame:
    manifest_path = DESIGN / "counterfactual_md_batch_manifest.tsv"
    manifest = read_tsv(manifest_path)
    if manifest.empty or controls.empty:
        return manifest
    updated = manifest.copy()
    pos_mask = updated["control_type"].astype(str).eq("same_or_similar_hla_positive_control")
    for idx, row in updated[pos_mask].iterrows():
        sub = controls[
            (controls["target_peptide"].astype(str) == str(row["peptide"]))
            & (controls["target_hla_4digit"].astype(str) == str(row["hla_4digit"]))
        ].copy()
        if sub.empty:
            continue
        # Prefer structure-backed control for pMHC/TCR-pMHC if present, otherwise exact-HLA sequence control.
        if str(row["complex_kind"]) == "TCR-pMHC":
            preferred = sub[(sub["has_structure"]) & (sub["paired_tcr_available"])]
            if preferred.empty:
                preferred = sub[sub["paired_tcr_available"]]
        else:
            preferred = sub[sub["has_structure"]]
        if preferred.empty:
            preferred = sub
        chosen = preferred.sort_values("selection_score", ascending=False).iloc[0]
        updated.loc[idx, "sequence"] = chosen["control_peptide"]
        updated.loc[idx, "sequence_status"] = (
            f"selected_positive_control_{chosen['match_class']}; "
            f"tcr_row={chosen['control_tcr_row_id']}; pdb={chosen['structure_pdb_id']}"
        )
        updated.loc[idx, "readiness"] = "ready_existing_positive_control_pdb" if bool(chosen["has_structure"]) else "structure_generation_required"
        updated.loc[idx, "blocked_by"] = "" if bool(chosen["has_structure"]) else "prepare_positive_control_structure"
        updated.loc[idx, "rationale"] = "positive control for MD runtime/contact-analysis sanity; not cancer specificity evidence"
    updated.to_csv(manifest_path, sep="\t", index=False)
    return updated


def write_report(controls: pd.DataFrame, manifest: pd.DataFrame) -> None:
    DESIGN.mkdir(parents=True, exist_ok=True)
    controls.to_csv(DESIGN / "positive_control_candidates.tsv", sep="\t", index=False)
    selected = manifest[manifest.get("control_type", pd.Series(dtype=str)).astype(str).eq("same_or_similar_hla_positive_control")].copy()
    selected.to_csv(DESIGN / "positive_control_selected_manifest_rows.tsv", sep="\t", index=False)
    lines = [
        "# MD Positive Control Selection",
        "",
        "## Boundary",
        "",
        "Positive controls are sanity controls for MD runtime, contact analysis, and visualization. They are not cancer specificity or neoantigen validation evidence.",
        "",
        "## Summary",
        "",
        f"- candidate control rows: {len(controls)}",
        f"- selected manifest rows updated: {len(selected)}",
        "",
        "## Selected Rows",
        "",
        selected[
            [
                c
                for c in [
                    "batch_id",
                    "peptide",
                    "hla_4digit",
                    "complex_kind",
                    "sequence",
                    "sequence_status",
                    "readiness",
                    "blocked_by",
                ]
                if c in selected.columns
            ]
        ].head(20).to_markdown(index=False)
        if not selected.empty
        else "No selected rows.",
    ]
    (DESIGN / "positive_control_selection_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    DESIGN.mkdir(parents=True, exist_ok=True)
    controls = select_controls()
    manifest = update_batch_manifest(controls)
    write_report(controls, manifest)
    exact = int((controls.get("match_class", pd.Series(dtype=str)) == "exact_hla").sum()) if not controls.empty else 0
    with_structure = int(controls.get("has_structure", pd.Series(dtype=bool)).sum()) if not controls.empty else 0
    print(f"[md-positive-controls] candidates={len(controls)} exact_hla={exact} structure_backed={with_structure}")
    print(DESIGN / "positive_control_candidates.tsv")


if __name__ == "__main__":
    main()
