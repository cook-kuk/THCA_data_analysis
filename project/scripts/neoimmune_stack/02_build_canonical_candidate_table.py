#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from neoimmune_common import (
    CANONICAL_COLUMNS,
    HUB_ROOT,
    REPO_ROOT,
    ensure_run_dir,
    first_existing,
    safe_read_table,
    write_md,
    write_tsv_gz,
)


def coalesce(df: pd.DataFrame, names: list[str], default=np.nan) -> pd.Series:
    out = pd.Series(default, index=df.index)
    for n in names:
        if n in df.columns:
            out = out.where(out.notna() & (out.astype(str) != ""), df[n])
    return out


def col_as_str(df: pd.DataFrame, name: str, default: str = "") -> pd.Series:
    if name in df.columns:
        return df[name].fillna(default).astype(str)
    return pd.Series(default, index=df.index, dtype=str)


def canonicalize_clean_neobench(df: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    out["candidate_id"] = coalesce(df, ["candidate_id", "sample_id"])
    out["patient_id"] = coalesce(df, ["patient_id"]).fillna("unknown_patient")
    out["sample_id"] = coalesce(df, ["sample_id", "candidate_id"])
    out["dataset_source"] = coalesce(df, ["source_dataset", "source_name", "study_id"])
    out["tumor_type"] = coalesce(df, ["cancer_type", "disease_context"])
    out["gene"] = coalesce(df, ["gene"])
    out["mutation_id"] = coalesce(df, ["mutation_id"])
    out["source_type"] = "clean_neobench"
    out["peptide_mut"] = coalesce(df, ["mut_peptide", "peptide"])
    out["peptide_wt"] = coalesce(df, ["wt_peptide"])
    out["peptide_length"] = coalesce(df, ["peptide_length"])
    out["flank_left"] = np.nan
    out["flank_right"] = np.nan
    out["hla_allele"] = coalesce(df, ["hla", "hla_allele_4digit"])
    out["hla_class"] = coalesce(df, ["mhc_class"])
    out["expression_tpm"] = coalesce(df, ["expression_tpm", "mutant_expression"])
    out["vaf"] = coalesce(df, ["vaf"])
    out["clonality"] = coalesce(df, ["clonality"])
    out["hla_loh_status"] = coalesce(df, ["hla_loh"])
    out["apm_score"] = coalesce(df, ["antigen_processing_status", "immune_context_score"])
    out["b2m_status"] = coalesce(df, ["b2m_status"])
    out["tap1_expr"] = np.nan
    out["tap2_expr"] = np.nan
    label_type = col_as_str(df, "label_type")
    label = coalesce(df, ["label"])
    is_presentation = label_type.str.contains("presentation|ligand|ms", case=False, na=False)
    out["label_presentation"] = np.where(is_presentation, label, np.nan)
    out["label_immunogenicity"] = np.where(~is_presentation, label, np.nan)
    out["assay_type"] = coalesce(df, ["label_type", "assay_type"])
    out["assay_result"] = coalesce(df, ["label"])
    out["validation_level"] = coalesce(df, ["label_type", "leakage_risk_level"])
    out["train_test_group"] = coalesce(df, ["split_source_heldout", "split_study_heldout", "split_patient_heldout"])
    out["source_study"] = coalesce(df, ["study_id", "source_name"])
    for c in df.columns:
        if c.startswith("split_") or "overlap" in c or c in ["leakage_risk_level", "hla_supertype", "source_protein_window"]:
            out[c] = df[c]
    return out


def canonicalize_cross_neo(df: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    out["candidate_id"] = coalesce(df, ["sample_id"])
    miss = out["candidate_id"].isna() | out["candidate_id"].astype(str).eq("")
    out.loc[miss, "candidate_id"] = [f"CNV0_{i:05d}" for i in range(int(miss.sum()))]
    out["patient_id"] = coalesce(df, ["patient_id"]).fillna("unknown_patient")
    out["sample_id"] = coalesce(df, ["sample_id"])
    out["dataset_source"] = coalesce(df, ["study"])
    out["tumor_type"] = "cancer"
    out["gene"] = coalesce(df, ["source_protein"])
    out["mutation_id"] = np.nan
    out["source_type"] = "cross_neo_v0"
    out["peptide_mut"] = coalesce(df, ["peptide_mut"])
    out["peptide_wt"] = coalesce(df, ["peptide_wt"])
    out["peptide_length"] = out["peptide_mut"].fillna("").astype(str).str.len()
    out["flank_left"] = np.nan
    out["flank_right"] = np.nan
    out["hla_allele"] = coalesce(df, ["hla"])
    out["hla_class"] = "I"
    out["expression_tpm"] = np.nan
    out["vaf"] = np.nan
    out["clonality"] = np.nan
    out["hla_loh_status"] = np.nan
    out["apm_score"] = np.nan
    out["b2m_status"] = np.nan
    out["tap1_expr"] = np.nan
    out["tap2_expr"] = np.nan
    out["label_presentation"] = np.nan
    out["label_immunogenicity"] = coalesce(df, ["label"])
    out["assay_type"] = coalesce(df, ["assay_type"])
    out["assay_result"] = coalesce(df, ["label"])
    out["validation_level"] = coalesce(df, ["strict_set_flag", "public_overlap_flags"])
    out["train_test_group"] = coalesce(df, ["split"])
    out["source_study"] = coalesce(df, ["study"])
    for c in ["strict_set_flag", "public_overlap_flags", "split", "near_peptide_cluster", "tcr_motif_score", "self_exact_match", "self_hamming1_count_log", "self_hamming2_count_log", "self_blosum_max", "source_window_15aa", "source_window_30aa"]:
        if c in df.columns:
            out[c] = df[c]
    return out


def canonicalize_hub(df: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    out["candidate_id"] = coalesce(df, ["canonical_id", "source_record_id"])
    out["patient_id"] = coalesce(df, ["patient_id"]).fillna("unknown_patient")
    out["sample_id"] = coalesce(df, ["canonical_id", "source_record_id"])
    out["dataset_source"] = coalesce(df, ["source_dataset"])
    out["tumor_type"] = coalesce(df, ["cancer_type", "tumor_type_raw"])
    out["gene"] = coalesce(df, ["gene"])
    out["mutation_id"] = coalesce(df, ["protein_change", "dna_change"])
    out["source_type"] = "neoantigen_hub"
    out["peptide_mut"] = coalesce(df, ["mt_peptide", "peptide_core"])
    out["peptide_wt"] = coalesce(df, ["wt_peptide"])
    out["peptide_length"] = coalesce(df, ["peptide_length"])
    out["flank_left"] = np.nan
    out["flank_right"] = np.nan
    out["hla_allele"] = coalesce(df, ["hla_normalized", "hla_allele"])
    out["hla_class"] = coalesce(df, ["mhc_class"])
    out["expression_tpm"] = coalesce(df, ["expression_tpm"])
    out["vaf"] = coalesce(df, ["variant_allele_fraction"])
    out["clonality"] = coalesce(df, ["clonality_if_available"])
    out["hla_loh_status"] = np.nan
    out["apm_score"] = np.nan
    out["b2m_status"] = np.nan
    out["tap1_expr"] = np.nan
    out["tap2_expr"] = np.nan
    out["label_presentation"] = coalesce(df, ["ms_presented"])
    out["label_immunogenicity"] = coalesce(df, ["immunogenic_positive", "y"])
    out["assay_type"] = coalesce(df, ["assay_type", "assay_platform"])
    out["assay_result"] = coalesce(df, ["validation_label"])
    out["validation_level"] = coalesce(df, ["validation_level", "test_set_safety"])
    out["train_test_group"] = coalesce(df, ["test_set_safety"])
    out["source_study"] = coalesce(df, ["source_dataset", "pmid"])
    for c in ["binding_affinity_nm", "binding_affinity_nM", "binding_rank", "el_score", "presentation_score", "prime_score", "bigmhc_score", "foreignness_score", "self_similarity_score", "tcr_recognition_score", "test_set_safety", "leakage_reason", "conflict_flag"]:
        if c in df.columns:
            out[c] = df[c]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)

    sources = []
    p = first_existing([REPO_ROOT / "project/results/clean_neobench_barneo_2026_05_09/clean_neobench_master.tsv"])
    if p:
        sources.append((p, canonicalize_clean_neobench(safe_read_table(p))))
    p = first_existing([REPO_ROOT / "project/results/cross_neo_v0/master_table.tsv"])
    if p:
        sources.append((p, canonicalize_cross_neo(safe_read_table(p))))
    p = first_existing([HUB_ROOT / "data_processed/benchmark_ready_immunogenicity.csv", HUB_ROOT / "data_processed/neoantigen_master.csv"])
    if p:
        sources.append((p, canonicalize_hub(safe_read_table(p))))

    if not sources:
        raise SystemExit("No candidate source tables found.")

    df = pd.concat([d for _, d in sources], ignore_index=True, sort=False)
    for c in CANONICAL_COLUMNS:
        if c not in df.columns:
            df[c] = np.nan
    df["peptide_mut"] = df["peptide_mut"].fillna("").astype(str).str.strip()
    df["hla_allele"] = df["hla_allele"].fillna("").astype(str).str.strip()
    df = df[df["peptide_mut"].ne("")]
    df["candidate_id"] = df["candidate_id"].fillna("").astype(str)
    missing = df["candidate_id"].eq("") | df["candidate_id"].eq("nan")
    df.loc[missing, "candidate_id"] = [f"NIS_{i:07d}" for i in range(int(missing.sum()))]
    df = df.drop_duplicates(subset=["candidate_id"], keep="first")
    df = df[CANONICAL_COLUMNS + [c for c in df.columns if c not in CANONICAL_COLUMNS]]

    write_tsv_gz(df, outdir / "data" / "canonical_candidates.tsv.gz")
    try:
        parquet_df = df.copy()
        for c in parquet_df.columns:
            if parquet_df[c].dtype == "object":
                parquet_df[c] = parquet_df[c].astype("string")
        parquet_df.to_parquet(outdir / "data" / "canonical_candidates.parquet", index=False)
        parquet_msg = "parquet written"
    except Exception as e:
        parquet_msg = f"parquet unavailable: {e!r}"
        (outdir / "data" / "canonical_candidates.parquet.ERROR.txt").write_text(parquet_msg, encoding="utf-8")

    md = [
        "# Canonical candidate table",
        "",
        f"- Rows: {len(df):,}",
        f"- Columns: {len(df.columns):,}",
        f"- Parquet status: {parquet_msg}",
        "",
        "## Sources used",
    ]
    for p, d in sources:
        md.append(f"- `{p}`: {len(d):,} rows before merge")
    md += [
        "",
        "## Claim boundary",
        "This table standardizes candidate records and labels; it does not validate vaccine efficacy or antigen presentation.",
    ]
    write_md("\n".join(md) + "\n", outdir / "data" / "canonical_candidates_README.md")


if __name__ == "__main__":
    main()
