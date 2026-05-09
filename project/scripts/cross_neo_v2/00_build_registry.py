#!/usr/bin/env python3
"""Build CROSS-Neo 2.0 canonical registry."""

from __future__ import annotations

import pandas as pd

from common import OUT, V0, ensure_dirs, hla_gene, hla_supertype, normalize_hla, peptide_cluster, safe_parquet


def main() -> None:
    ensure_dirs()
    src = pd.read_csv(V0 / "master_table.tsv", sep="\t")
    reg = pd.DataFrame()
    reg["row_id"] = src["sample_id"].astype(str)
    reg["source_dataset"] = src["study"].fillna("unknown").astype(str)
    reg["study_id"] = src["study"].fillna("unknown").astype(str)
    for col in ["patient_id", "assay_type"]:
        reg[col] = src[col] if col in src else pd.NA
    reg["cancer_type"] = pd.NA
    reg["pathogen_type"] = pd.NA
    reg["peptide"] = src["peptide_mut"].astype(str)
    reg["mutant_peptide"] = src["peptide_mut"].astype(str)
    reg["wildtype_peptide"] = src["peptide_wt"] if "peptide_wt" in src else pd.NA
    reg["hla_raw"] = src["hla"].astype(str)
    reg["hla_4digit"] = reg["hla_raw"].map(normalize_hla)
    reg["hla_gene"] = reg["hla_4digit"].map(hla_gene)
    reg["hla_supertype"] = reg["hla_4digit"].map(hla_supertype)
    reg["peptide_length"] = reg["peptide"].astype(str).str.len()
    reg["mutation_position"] = pd.NA
    for col in ["source_protein", "source_window_15aa", "source_window_30aa"]:
        reg[col] = src[col] if col in src else pd.NA
    reg["left_flank"] = pd.NA
    reg["right_flank"] = pd.NA
    reg["source_window"] = reg["source_window_30aa"].combine_first(reg["source_window_15aa"])
    reg["label_raw"] = src["label"]
    reg["label_binary"] = pd.to_numeric(src["label"], errors="coerce").astype("Int64")
    reg["label_strength"] = pd.NA
    for col in ["expression", "TPM", "RNA_evidence", "clonality", "VAF", "train_test_original", "public_overlap_flags"]:
        reg[col] = src[col] if col in src else pd.NA
    reg["train_test_original"] = src["split"] if "split" in src else reg["train_test_original"]
    reg["strict_set_flag"] = src["strict_set_flag"].astype(bool) if "strict_set_flag" in src else False
    reg["near_peptide_cluster"] = src["near_peptide_cluster"] if "near_peptide_cluster" in src else reg["peptide"].map(peptide_cluster)
    reg["exact_peptide_hla_key"] = reg["peptide"].astype(str) + "|" + reg["hla_4digit"].astype(str)
    reg["peptide_only_key"] = reg["peptide"].astype(str)
    reg["mutant_wt_hla_key"] = reg["mutant_peptide"].astype(str) + "|" + reg["wildtype_peptide"].astype(str) + "|" + reg["hla_4digit"].astype(str)
    reg["near_peptide_cluster_key"] = reg["near_peptide_cluster"].astype(str)

    for col in [
        "patient_id", "assay_type", "wildtype_peptide", "source_protein", "source_window", "expression",
        "TPM", "RNA_evidence", "clonality", "VAF", "public_overlap_flags",
    ]:
        reg[f"missing_{col}"] = reg[col].isna() | reg[col].astype(str).isin(["", "nan", "None", "NA"])

    dup_cols = ["exact_peptide_hla_key", "peptide_only_key", "mutant_wt_hla_key", "near_peptide_cluster_key"]
    for col in dup_cols:
        reg[f"duplicate_count_{col}"] = reg.groupby(col)["row_id"].transform("size")

    safe_parquet(reg, OUT / "canonical_registry.parquet")
    reg.to_csv(OUT / "canonical_registry.tsv", sep="\t", index=False, na_rep="NA")

    miss = []
    for c in reg.columns:
        miss.append({"column": c, "missing_n": int(reg[c].isna().sum()), "missing_frac": float(reg[c].isna().mean())})
    pd.DataFrame(miss).to_csv(OUT / "registry_missingness_report.tsv", sep="\t", index=False, na_rep="NA")
    reg.groupby("source_dataset").agg(n=("row_id", "size"), positives=("label_binary", "sum"), prevalence=("label_binary", "mean")).reset_index().to_csv(OUT / "registry_source_counts.tsv", sep="\t", index=False, na_rep="NA")
    reg.groupby(["source_dataset", "label_binary"]).size().reset_index(name="n").to_csv(OUT / "registry_label_balance.tsv", sep="\t", index=False, na_rep="NA")

    lines = [
        "# CROSS-Neo 2.0 Registry Report",
        "",
        f"Rows: {len(reg)}",
        f"Positive prevalence: {reg['label_binary'].mean():.4f}",
        "",
        "## Source Counts",
        "",
        reg.groupby("source_dataset").agg(n=("row_id", "size"), positives=("label_binary", "sum"), prevalence=("label_binary", "mean")).reset_index().to_markdown(index=False),
        "",
        "## Red Flags",
        "",
        f"- WT peptide missing fraction: {reg['missing_wildtype_peptide'].mean():.3f}",
        f"- Source window missing fraction: {reg['missing_source_window'].mean():.3f}",
        f"- Exact peptide-HLA duplicated rows: {int((reg['duplicate_count_exact_peptide_hla_key'] > 1).sum())}",
        "- Public predictor scores are not included as features.",
    ]
    (OUT / "00_registry_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v2-registry] rows={len(reg)} out={OUT}")


if __name__ == "__main__":
    main()
