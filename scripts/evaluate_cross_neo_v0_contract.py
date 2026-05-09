#!/usr/bin/env python3
"""Evaluation contract tables for CROSS-Neo v0."""

from __future__ import annotations

import pandas as pd

from cross_neo_v0_common import OUT, ece_score, ensure_dirs, load_master, metrics


def main() -> None:
    ensure_dirs()
    master = load_master()
    preds = pd.read_csv(OUT / "oof_predictions.tsv", sep="\t")
    retr = pd.read_csv(OUT / "retrieval_features_by_fold.tsv", sep="\t")
    merged = preds.merge(master[["sample_id", "hla", "hla_supertype", "study", "public_overlap_flags"]], on="sample_id", how="left")
    rows = []
    for (split, group, model), sub in merged.groupby(["split_name", "feature_group", "model"]):
        mm = metrics(sub["label"].to_numpy(), sub["score"].to_numpy())
        mm["ECE"] = ece_score(sub["label"].to_numpy(), sub["score"].to_numpy())
        rows.append({"contract_split": split, "feature_group": group, "model": model, **mm})
    out = pd.DataFrame(rows).sort_values(["contract_split", "AUPRC"], ascending=[True, False])
    out.to_csv(OUT / "evaluation_contract_metrics.tsv", sep="\t", index=False)

    best = out.sort_values("AUPRC", ascending=False).iloc[0]
    bpred = merged[
        (merged["split_name"] == best["contract_split"])
        & (merged["feature_group"] == best["feature_group"])
        & (merged["model"] == best["model"])
    ]
    per_hla = []
    for hla, sub in bpred.groupby("hla"):
        if len(sub) >= 3:
            per_hla.append({"hla": hla, **metrics(sub["label"].to_numpy(), sub["score"].to_numpy())})
    pd.DataFrame(per_hla).to_csv(OUT / "per_hla_metrics.tsv", sep="\t", index=False)

    per_study = []
    for study, sub in bpred.groupby("study"):
        per_study.append({"study": study, **metrics(sub["label"].to_numpy(), sub["score"].to_numpy())})
    pd.DataFrame(per_study).to_csv(OUT / "per_study_metrics.tsv", sep="\t", index=False)

    fail_rows = []
    b = bpred.merge(retr[["split_name", "fold_id", "sample_id", "retrieval_leakage_risk"]], on=["split_name", "fold_id", "sample_id"], how="left")
    for risk, sub in b.groupby("retrieval_leakage_risk", dropna=False):
        fail_rows.append({"failure_mode": f"retrieval_{risk}", **metrics(sub["label"].to_numpy(), sub["score"].to_numpy())})
    for flag, desc in [
        ("source_protein_window_holdout", "not available: source windows missing"),
        ("study_heldout", "not available: strict rows all ITSNdb"),
        ("time_heldout_if_date_available", "not available: only ITSNdb publication year currently populated"),
        ("public_corpus_overlap_audit", "partial only: in-house master overlap flag exists; public tool train rows not downloaded/audited"),
    ]:
        fail_rows.append({"failure_mode": flag, "note": desc})
    pd.DataFrame(fail_rows).to_csv(OUT / "failure_modes.tsv", sep="\t", index=False)
    print(f"[eval] best={best['feature_group']} {best['model']} AUPRC={best['AUPRC']:.3f}")


if __name__ == "__main__":
    main()
