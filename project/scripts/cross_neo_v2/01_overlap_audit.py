#!/usr/bin/env python3
"""CROSS-Neo 2.0 overlap/leakage audit."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from common import OUT, V1, ensure_dirs, kmer_set, levenshtein, read_table, seq_similarity


def jaccard(a: set[str], b: set[str]) -> float:
    return 0.0 if not a or not b else len(a & b) / len(a | b)


def split_rows(reg: pd.DataFrame) -> list[dict[str, object]]:
    rows = []
    strict = reg[reg["strict_set_flag"].astype(bool)]
    folds_path = OUT.parent / "cross_neo_v0/folds.tsv"
    if folds_path.exists():
        folds = pd.read_csv(folds_path, sep="\t")
        for (split, fold), g in folds.groupby(["split_name", "fold_id"]):
            test_ids = set(g["sample_id"].astype(str))
            train_ids = set(strict["row_id"]) - test_ids
            rows.append({"split_name": split, "fold_id": fold, "train_ids": train_ids, "test_ids": test_ids})
    for src in ["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"]:
        test_ids = set(reg.loc[reg["source_dataset"].eq(src), "row_id"])
        train_ids = set(reg.loc[reg["source_dataset"].isin(["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"]) & ~reg["source_dataset"].eq(src), "row_id"])
        if test_ids and train_ids:
            rows.append({"split_name": f"source_heldout_{src}", "fold_id": src, "train_ids": train_ids, "test_ids": test_ids})
    return rows


def main() -> None:
    ensure_dirs()
    reg = read_table(OUT / "canonical_registry.tsv")
    public_manifest = pd.read_csv(V1 / "public_overlap_manifest.tsv", sep="\t") if (V1 / "public_overlap_manifest.tsv").exists() else pd.DataFrame()
    by_split, by_row = [], []
    for item in split_rows(reg):
        train = reg[reg["row_id"].isin(item["train_ids"])]
        test = reg[reg["row_id"].isin(item["test_ids"])]
        train_peps = train["peptide"].astype(str).tolist()
        train_pep_set = set(train_peps)
        train_pair_set = set(train["exact_peptide_hla_key"].astype(str))
        train_mwh = set(train["mutant_wt_hla_key"].astype(str))
        train_hla = set(train["hla_4digit"].astype(str))
        train_super = set(train["hla_supertype"].astype(str))
        train_source = set(train["source_dataset"].astype(str))
        train_prot = set(train["source_protein"].dropna().astype(str))
        row_stats = []
        for _, r in test.iterrows():
            pep = str(r["peptide"])
            lev1 = any(levenshtein(pep, p, 1) <= 1 for p in train_peps[:5000])
            lev2 = any(levenshtein(pep, p, 2) <= 2 for p in train_peps[:5000])
            km = kmer_set(pep, 4)
            max_j = max((jaccard(km, kmer_set(p, 4)) for p in train_peps[:5000]), default=0.0)
            rec = {
                "split_name": item["split_name"], "fold_id": item["fold_id"], "row_id": r["row_id"],
                "exact_peptide_overlap": int(pep in train_pep_set),
                "exact_peptide_hla_overlap": int(str(r["exact_peptide_hla_key"]) in train_pair_set),
                "mutant_wt_hla_overlap": int(str(r["mutant_wt_hla_key"]) in train_mwh),
                "levenshtein_le1_overlap": int(lev1),
                "levenshtein_le2_overlap": int(lev2),
                "max_kmer_jaccard": max_j,
                "near_kmer_jaccard_ge_0_5": int(max_j >= 0.5),
                "source_protein_overlap": int(str(r["source_protein"]) in train_prot) if pd.notna(r["source_protein"]) else 0,
                "hla_allele_overlap": int(str(r["hla_4digit"]) in train_hla),
                "hla_supertype_overlap": int(str(r["hla_supertype"]) in train_super),
                "study_source_overlap": int(str(r["source_dataset"]) in train_source),
            }
            risk_score = rec["exact_peptide_hla_overlap"] * 3 + rec["exact_peptide_overlap"] * 2 + rec["levenshtein_le1_overlap"] + rec["near_kmer_jaccard_ge_0_5"]
            rec["leakage_risk"] = "high" if risk_score >= 3 else ("medium" if risk_score else "clean_or_shift")
            row_stats.append(rec)
        rr = pd.DataFrame(row_stats)
        by_row.append(rr)
        agg = {"split_name": item["split_name"], "fold_id": item["fold_id"], "train_n": len(train), "test_n": len(test)}
        for c in rr.columns:
            if c not in {"split_name", "fold_id", "row_id", "leakage_risk"}:
                agg[c + "_rate"] = float(rr[c].mean()) if len(rr) else 0.0
        agg["high_risk_n"] = int((rr["leakage_risk"] == "high").sum()) if len(rr) else 0
        by_split.append(agg)
    row_df = pd.concat(by_row, ignore_index=True)
    split_df = pd.DataFrame(by_split)
    row_df.to_csv(OUT / "overlap_audit_by_row.tsv", sep="\t", index=False, na_rep="NA")
    split_df.to_csv(OUT / "overlap_audit_by_split.tsv", sep="\t", index=False, na_rep="NA")
    split_df.groupby("split_name").mean(numeric_only=True).reset_index().to_csv(OUT / "leakage_risk_summary.tsv", sep="\t", index=False, na_rep="NA")
    missing_public = public_manifest[~public_manifest.get("status", "").astype(str).eq("local_audit_available")] if len(public_manifest) else pd.DataFrame()
    missing_public.to_csv(OUT / "public_training_corpus_missing_manifest.tsv", sep="\t", index=False, na_rep="NA")

    heat = split_df.groupby("split_name").mean(numeric_only=True)[["exact_peptide_hla_overlap_rate", "levenshtein_le1_overlap_rate", "near_kmer_jaccard_ge_0_5_rate", "study_source_overlap_rate"]]
    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(heat.values, aspect="auto", cmap="magma")
    ax.set_xticks(range(len(heat.columns))); ax.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(heat.index))); ax.set_yticklabels(heat.index, fontsize=7)
    fig.colorbar(im, ax=ax, label="overlap rate")
    fig.tight_layout(); fig.savefig(OUT / "overlap_heatmap.png", dpi=180); plt.close(fig)

    lines = ["# CROSS-Neo 2.0 Overlap Audit", "", split_df.groupby("split_name").mean(numeric_only=True).reset_index().to_markdown(index=False), "", "## Public Corpus Gaps", "", missing_public.to_markdown(index=False) if len(missing_public) else "No missing public corpus manifest rows."]
    (OUT / "01_overlap_audit_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v2-overlap] rows={len(row_df)} splits={split_df['split_name'].nunique()}")


if __name__ == "__main__":
    main()
