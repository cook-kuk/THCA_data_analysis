#!/usr/bin/env python3
"""Step 3 — build cohort-level harmonized sample + feature index.

Reads `processed/mtbls3339_*` and `processed/cellrepmed2026_*` outputs
from steps 01 and 02. Constructs:
  - harmonized_thyroid_multiomics_sample_index.tsv
  - harmonized_thyroid_multiomics_feature_index.tsv

**Important**: this is **cohort-level**, not sample-level cross-merge.
Each row stays inside its source cohort; we just normalize column
names so downstream code can iterate uniformly. Every row carries a
`mapping_confidence` column:
  - high   = direct named field (e.g. an explicit "Subtype" column)
  - medium = inferred from filename / sheet category
  - low    = heuristic only (column-pattern guess)
HLA validation data is intentionally excluded from this index per
project boundary (`v17_arcasHLA_korean_k2`).
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path
import pandas as pd

SAMPLE_FIELDS = ["cohort", "sample_id", "patient_id", "disease_stage",
                 "histology", "recurrence_risk", "subtype", "modality",
                 "source_file", "mapping_confidence"]


def lower_lookup(df: pd.DataFrame, *candidates: str) -> str | None:
    cols_l = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_l:
            return cols_l[cand.lower()]
    return None


def harmonize_mtbls(processed_dir: Path) -> pd.DataFrame:
    smp_path = processed_dir / "mtbls3339_sample_metadata.tsv"
    if not smp_path.exists():
        return pd.DataFrame(columns=SAMPLE_FIELDS)
    df = pd.read_csv(smp_path, sep="\t", dtype=str, keep_default_na=False)
    out = []
    sid_col = lower_lookup(df, "Sample Name", "Sample_Name", "sample_id", "Source Name")
    pid_col = lower_lookup(df, "Subject", "Subject ID", "patient", "patient_id")
    hist_col = lower_lookup(df, "Characteristics[disease]", "Characteristics[Disease]",
                            "Characteristics[organism part]", "histology")
    stage_col = lower_lookup(df, "Factor Value[stage]", "Factor Value[Stage]", "stage")
    recur_col = lower_lookup(df, "Factor Value[recurrence risk]",
                             "Factor Value[Recurrence Risk]", "Recurrence Risk",
                             "recurrence_risk")
    for _, r in df.iterrows():
        sid = r[sid_col] if sid_col else ""
        pid = r[pid_col] if pid_col else ""
        hist = r[hist_col] if hist_col else ""
        stage = r[stage_col] if stage_col else ""
        recur = r[recur_col] if recur_col else ""
        confidence = "high" if sid_col else "low"
        out.append({
            "cohort": "MTBLS3339",
            "sample_id": sid, "patient_id": pid,
            "disease_stage": stage, "histology": hist,
            "recurrence_risk": recur, "subtype": "",
            "modality": "metabolomics",
            "source_file": smp_path.name,
            "mapping_confidence": confidence,
        })
    return pd.DataFrame(out, columns=SAMPLE_FIELDS)


def harmonize_crm_clinical(processed_dir: Path) -> pd.DataFrame:
    """Prefer the authoritative PDF-extracted Table S1 (FUSCC-1…FUSCC-113 with CC subtype)
    when present; fall back to xlsx-derived clinical heuristic otherwise."""
    out = []
    s1_path = processed_dir / "cellrepmed2026_table_S1_clinical.tsv"
    if s1_path.exists():
        df = pd.read_csv(s1_path, sep="\t", dtype=str, keep_default_na=False)
        # Strict patient row filter — drop trailing abbreviation/legend rows.
        df = df[df["Patient ID"].astype(str).str.match(r"^FUSCC-\d+$")].copy()
        for _, r in df.iterrows():
            out.append({
                "cohort": "CellRepMed2026",
                "sample_id": r["Patient ID"],
                "patient_id": r["Patient ID"],
                "disease_stage": r.get("Status of advanced disease", ""),
                "histology": r.get("Histological Type", "advanced_DTC"),
                "recurrence_risk": r.get("RAI sensitivity", ""),
                "subtype": r.get("CC subtype", ""),
                "modality": "proteomic+NGS+pathology",
                "source_file": s1_path.name,
                "mapping_confidence": "high",
            })
        return pd.DataFrame(out, columns=SAMPLE_FIELDS)

    # Fallback: heuristic xlsx clinical sheet concat from 02_ingest.
    clin_path = processed_dir / "cellrepmed2026_clinical.tsv"
    if not clin_path.exists():
        return pd.DataFrame(columns=SAMPLE_FIELDS)
    df = pd.read_csv(clin_path, sep="\t", dtype=str, keep_default_na=False)
    sid_col = lower_lookup(df, "sample_id", "sample", "case_id", "caseid", "id",
                           "patient_id", "patient")
    sub_col = lower_lookup(df, "subtype", "cluster", "consensus", "consensus_cluster")
    for _, r in df.iterrows():
        sid = r[sid_col] if sid_col else ""
        if not sid:
            cands = [c for c in df.columns if not c.startswith("__source")]
            if cands:
                sid = str(r[cands[0]])
        confidence = "high" if (sid_col and sub_col) else ("medium" if sid_col else "low")
        out.append({
            "cohort": "CellRepMed2026",
            "sample_id": str(sid), "patient_id": str(sid),
            "disease_stage": "",
            "histology": "advanced_DTC",
            "recurrence_risk": "",
            "subtype": str(r[sub_col]) if sub_col else "",
            "modality": "proteomic",
            "source_file": f'{r.get("__source_file__","")}::{r.get("__source_sheet__","")}',
            "mapping_confidence": confidence,
        })
    return pd.DataFrame(out, columns=SAMPLE_FIELDS)


def harmonize_crm_subtypes(processed_dir: Path) -> pd.DataFrame:
    sub_path = processed_dir / "cellrepmed2026_subtypes.tsv"
    if not sub_path.exists():
        return pd.DataFrame(columns=SAMPLE_FIELDS)
    df = pd.read_csv(sub_path, sep="\t", dtype=str, keep_default_na=False)
    out = []
    sid_col = lower_lookup(df, "sample_id", "sample", "case_id", "id", "patient")
    sub_col = lower_lookup(df, "subtype", "cluster", "consensus", "cc1", "cc2", "cc3")
    for _, r in df.iterrows():
        sid = r[sid_col] if sid_col else ""
        sub = r[sub_col] if sub_col else ""
        if not sid:
            continue
        out.append({
            "cohort": "CellRepMed2026",
            "sample_id": str(sid), "patient_id": str(sid),
            "disease_stage": "", "histology": "advanced_DTC",
            "recurrence_risk": "",
            "subtype": str(sub),
            "modality": "proteomic_subtype_call",
            "source_file": f'{r.get("__source_file__","")}::{r.get("__source_sheet__","")}',
            "mapping_confidence": "high" if (sid_col and sub_col) else "medium",
        })
    return pd.DataFrame(out, columns=SAMPLE_FIELDS)


def build_feature_index(processed_dir: Path) -> pd.DataFrame:
    rows = []
    feat_anno = processed_dir / "mtbls3339_feature_annotation.tsv"
    if feat_anno.exists():
        df = pd.read_csv(feat_anno, sep="\t", dtype=str, keep_default_na=False)
        rows.append(pd.DataFrame({
            "cohort": "MTBLS3339",
            "feature_type": "metabolite",
            "feature_id": df["__feature_id__"] if "__feature_id__" in df.columns
                          else df.iloc[:, 0],
            "annotation_summary": df.astype(str).agg(" | ".join, axis=1).str[:200],
            "source_file": df.get("source_maf", ""),
        }))
    crm_feat = processed_dir / "cellrepmed2026_feature_annotation.tsv"
    if crm_feat.exists():
        df = pd.read_csv(crm_feat, sep="\t", dtype=str, keep_default_na=False)
        rows.append(pd.DataFrame({
            "cohort": "CellRepMed2026",
            "feature_type": "protein_or_gene",
            "feature_id": df.iloc[:, 0],
            "annotation_summary": df.astype(str).agg(" | ".join, axis=1).str[:200],
            "source_file": df.get("__source_file__", ""),
        }))
    crm_prot = processed_dir / "cellrepmed2026_protein_matrix_wide.tsv"
    if crm_prot.exists():
        df = pd.read_csv(crm_prot, sep="\t", dtype=str, keep_default_na=False, nrows=200000)
        rows.append(pd.DataFrame({
            "cohort": "CellRepMed2026",
            "feature_type": "protein_matrix_row",
            "feature_id": df.iloc[:, 0],
            "annotation_summary": "",
            "source_file": df.get("__source_file__", ""),
        }))
    if not rows:
        return pd.DataFrame(columns=["cohort", "feature_type", "feature_id",
                                     "annotation_summary", "source_file"])
    return pd.concat(rows, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--processed-dir", required=True, type=Path)
    args = ap.parse_args()
    if not args.processed_dir.exists():
        sys.exit(f"ERROR: --processed-dir not found: {args.processed_dir}")

    pieces = []
    pieces.append(harmonize_mtbls(args.processed_dir))
    pieces.append(harmonize_crm_clinical(args.processed_dir))
    pieces.append(harmonize_crm_subtypes(args.processed_dir))
    sample_idx = pd.concat(pieces, ignore_index=True)
    sample_idx = sample_idx[sample_idx["sample_id"].astype(str).str.len() > 0]

    sample_path = args.processed_dir / "harmonized_thyroid_multiomics_sample_index.tsv"
    sample_idx.to_csv(sample_path, sep="\t", index=False)

    feat_idx = build_feature_index(args.processed_dir)
    feat_path = args.processed_dir / "harmonized_thyroid_multiomics_feature_index.tsv"
    feat_idx.to_csv(feat_path, sep="\t", index=False)

    print(f"[03] sample index: {len(sample_idx)} rows → {sample_path}")
    print(f"[03] feature index: {len(feat_idx)} rows → {feat_path}")
    if len(sample_idx) == 0:
        print("[03] WARNING: empty sample index — upstream ingest produced no usable rows.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
