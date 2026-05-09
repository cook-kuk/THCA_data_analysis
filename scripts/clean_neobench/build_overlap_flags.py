#!/usr/bin/env python3
"""Add strict leakage and overlap flags to the CLEAN-NeoBench master table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import ensure_dir, near_similarity, normalize_empty, update_manifest, write_tsv


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="Output directory")
    parser.add_argument("--near-threshold", type=float, default=0.80, help="Near-peptide similarity threshold")
    return parser.parse_args()


def as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    s = normalize_empty(value).lower()
    return s in {"1", "true", "yes", "y"}


def train_mask(master: pd.DataFrame) -> pd.Series:
    if "_original_split" in master.columns:
        return master["_original_split"].fillna("").astype(str).str.lower().eq("train")
    if "source_name" in master.columns:
        return master["source_name"].fillna("").astype(str).str.lower().isin({"cedar", "nepdb", "tesla_mmc4", "tesla_mmc7"})
    return pd.Series(False, index=master.index)


def leakage_level(row: pd.Series) -> str:
    if as_bool(row.get("exact_peptide_hla_train_overlap")) or as_bool(row.get("public_tool_training_overlap_any")):
        return "high"
    if as_bool(row.get("exact_peptide_train_overlap")) or as_bool(row.get("near_peptide_train_overlap")):
        return "medium"
    if as_bool(row.get("study_train_overlap")) or as_bool(row.get("patient_train_overlap")):
        return "medium"
    return "low"


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    ensure_dir(output_root)
    master_path = output_root / "clean_neobench_master.tsv"
    master = pd.read_csv(master_path, sep="\t")
    warnings: list[str] = []
    mask_train = train_mask(master)
    if not mask_train.any():
        warnings.append("no explicit train rows found; overlap flags use empty train sets")

    train = master[mask_train].copy()
    train_peptides = set(train["peptide"].fillna("").astype(str))
    train_phla = set(zip(train["peptide"].fillna("").astype(str), train["hla_allele_4digit"].fillna("").astype(str)))
    train_clusters = set(train["split_near_peptide_cluster"].fillna("").astype(str))
    train_windows = set(train["source_protein_window"].fillna("").astype(str)) - {""}
    train_studies = set(train["study_id"].fillna("").astype(str)) - {""}
    train_patients = set(train["patient_id"].fillna("").astype(str)) - {""}

    flags = master[["candidate_id"]].copy()
    flags["exact_peptide_train_overlap"] = master["peptide"].fillna("").astype(str).isin(train_peptides)
    flags["exact_peptide_hla_train_overlap"] = [
        (pep, hla) in train_phla for pep, hla in zip(master["peptide"].fillna("").astype(str), master["hla_allele_4digit"].fillna("").astype(str))
    ]
    flags["near_peptide_train_overlap"] = master["split_near_peptide_cluster"].fillna("").astype(str).isin(train_clusters)
    if not flags["near_peptide_train_overlap"].any() and len(train_peptides):
        # Expensive fallback only when cluster metadata are not useful.
        train_list = list(train_peptides)
        near = []
        for pep in master["peptide"].fillna("").astype(str):
            near.append(any(near_similarity(pep, tpep) >= args.near_threshold for tpep in train_list))
        flags["near_peptide_train_overlap"] = near
        warnings.append("near-peptide overlap used k-mer/identity fallback")
    flags["source_protein_window_train_overlap"] = master["source_protein_window"].fillna("").astype(str).isin(train_windows)
    flags["study_train_overlap"] = master["study_id"].fillna("").astype(str).isin(train_studies)
    flags["patient_train_overlap"] = master["patient_id"].fillna("").astype(str).isin(train_patients)
    flags["public_tool_training_overlap_any"] = master["public_tool_training_overlap_any"].map(as_bool)
    detail = master["public_tool_training_overlap_detail"].fillna("").astype(str)
    flags["public_tool_training_overlap_detail"] = detail.where(detail.str.len() > 0, "not_audited_or_unavailable")
    flags["leakage_risk_level"] = flags.apply(leakage_level, axis=1)

    for col in flags.columns:
        if col == "candidate_id":
            continue
        master[col] = flags[col]

    write_tsv(flags, output_root / "clean_neobench_overlap_flags.tsv")
    write_tsv(master, master_path)
    update_manifest(
        output_root,
        "build_overlap_flags",
        {
            "n_candidates": int(len(master)),
            "n_train_reference": int(mask_train.sum()),
            "overlap_counts": {col: int(flags[col].sum()) for col in flags.columns if col.endswith("_overlap") or col.endswith("_any")},
            "leakage_risk_counts": flags["leakage_risk_level"].value_counts().to_dict(),
            "warnings": warnings,
        },
    )
    print(f"[clean-neobench-overlap] rows={len(flags)} output={output_root / 'clean_neobench_overlap_flags.tsv'}")


if __name__ == "__main__":
    main()
