#!/usr/bin/env python3
"""Fold-safe retrieval features for CROSS-Neo v0."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cross_neo_v0_common import OUT, ensure_dirs, load_master, load_or_build_folds, seq_similarity


def max_similarity(query: str, refs: list[str]) -> float:
    if not refs:
        return 0.0
    return float(max(seq_similarity(query, r) for r in refs))


def best_source_similarity(row: pd.Series, train: pd.DataFrame) -> float:
    for col in ("source_window_30aa", "source_window_15aa", "source_protein"):
        query = str(row.get(col, "") or "")
        refs = [str(x) for x in train.get(col, pd.Series(dtype=str)).fillna("").tolist() if str(x)]
        if query and refs:
            return max_similarity(query, refs)
    return np.nan


def best_tcr_similarity(row: pd.Series, train: pd.DataFrame) -> float:
    vals = []
    for col in ("tcr_beta", "tcr_alpha"):
        query = str(row.get(col, "") or "")
        refs = [str(x) for x in train.get(col, pd.Series(dtype=str)).fillna("").tolist() if str(x)]
        if query and refs:
            vals.append(max_similarity(query, refs))
    return float(max(vals)) if vals else 0.0


def main() -> None:
    ensure_dirs()
    master = load_master()
    folds = load_or_build_folds(master)
    strict_ids = set(master.loc[master["strict_set_flag"].astype(bool), "sample_id"])
    rows = []
    for (split_name, fold_id), fold in folds.groupby(["split_name", "fold_id"], sort=False):
        test_ids = set(fold["sample_id"])
        train = master[master["sample_id"].isin(strict_ids - test_ids)].copy()
        test = master[master["sample_id"].isin(test_ids)].copy()
        train_pep = set(train["peptide_mut"])
        train_pair = set(zip(train["peptide_mut"], train["hla"]))
        train_pos = train[train["label"] == 1]
        train_pos_peps = train_pos["peptide_mut"].tolist()
        for _, r in test.iterrows():
            pep = str(r["peptide_mut"])
            hla = str(r["hla"])
            same_hla = train[train["hla"] == hla]
            same_hla_pos = same_hla[same_hla["label"] == 1]
            exact_pair = int((pep, hla) in train_pair)
            exact_pep = int(pep in train_pep)
            near_all = max_similarity(pep, train["peptide_mut"].tolist())
            near_pos = max_similarity(pep, train_pos_peps)
            same_hla_near = max_similarity(pep, same_hla["peptide_mut"].tolist())
            wt = str(r.get("peptide_wt", "") or "")
            wt_self = seq_similarity(pep, wt) if wt else float(r.get("self_hamming1_count_log", 0.0)) + 0.5 * float(r.get("self_hamming2_count_log", 0.0))
            risk = "clean_no_reference"
            if exact_pair or exact_pep:
                risk = "exact_hit"
            elif max(near_all, same_hla_near) >= 0.75:
                risk = "near_hit"
            if len(train) == 0:
                risk = "unknown"
            rows.append(
                {
                    "split_name": split_name,
                    "fold_id": fold_id,
                    "sample_id": r["sample_id"],
                    "exact_peptide_hla_hit_train": exact_pair,
                    "exact_peptide_hit_train": exact_pep,
                    "near_peptide_similarity_train": near_all,
                    "near_positive_peptide_similarity_train": near_pos,
                    "same_hla_near_peptide_similarity_train": same_hla_near,
                    "source_protein_window_similarity_train": best_source_similarity(r, train),
                    "wt_self_similarity": wt_self,
                    "tcr_motif_similarity_if_available": best_tcr_similarity(r, train),
                    "same_hla_train_density": int(len(same_hla)),
                    "same_hla_positive_rate_train": float(same_hla["label"].mean()) if len(same_hla) else np.nan,
                    "train_fold_prevalence": float(train["label"].mean()) if len(train) else np.nan,
                    "retrieval_leakage_risk": risk,
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "retrieval_features_by_fold.tsv", sep="\t", index=False)
    audit = out.groupby(["split_name", "retrieval_leakage_risk"]).size().reset_index(name="n")
    audit.to_csv(OUT / "retrieval_leakage_audit.tsv", sep="\t", index=False)
    lines = [
        "# Retrieval Leakage Audit",
        "",
        "Retrieval features were computed using only each fold's strict-set train rows.",
        "",
        audit.to_markdown(index=False),
    ]
    (OUT / "retrieval_leakage_audit.md").write_text("\n".join(lines) + "\n")
    print(f"[retrieval] wrote rows={len(out)} folds={out[['split_name','fold_id']].drop_duplicates().shape[0]}")


if __name__ == "__main__":
    main()
