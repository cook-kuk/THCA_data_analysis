#!/usr/bin/env python3
"""Score public industrial locked testset with transparent CROSS-lite adapters.

This is not the full KG-GA model. The public industrial v0 testset currently
contains peptide, HLA and response labels only. Full KG-GA needs additional
features such as TCR expert, MD/control, impact portfolio, foreignness and
CROSS_BMA. This script therefore creates explicitly labeled minimal adapters
using only features available at scoring time.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "project/results/public_industrial_neoantigen_benchmark_scout_2026_05_10"


def minmax(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    lo = x.min(skipna=True)
    hi = x.max(skipna=True)
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        return pd.Series(np.where(x.notna(), 0.5, np.nan), index=x.index)
    return (x - lo) / (hi - lo)


def peptide_quality_proxy(peptide: str) -> float:
    pep = str(peptide or "").upper()
    if not pep:
        return 0.0
    length_bonus = 1.0 if len(pep) in {9, 10, 11} else 0.65 if 8 <= len(pep) <= 15 else 0.0
    aromatic = sum(aa in "FWY" for aa in pep) / len(pep)
    charged = sum(aa in "KRDE" for aa in pep) / len(pep)
    hydrophobic = sum(aa in "AILMVPFWY" for aa in pep) / len(pep)
    # Mild non-label heuristic: reward plausible class-I length and moderate biochemical diversity.
    diversity = min(1.0, len(set(pep)) / 8)
    return float(0.45 * length_bonus + 0.20 * diversity + 0.15 * hydrophobic + 0.10 * aromatic + 0.10 * (1 - abs(charged - 0.22)))


def metric_row(df: pd.DataFrame, score_col: str, name: str) -> dict[str, object]:
    y = pd.to_numeric(df["label"], errors="coerce")
    s = pd.to_numeric(df[score_col], errors="coerce")
    mask = y.notna() & s.notna()
    d = df.loc[mask].copy()
    y = y[mask].astype(int)
    s = s[mask]
    row: dict[str, object] = {
        "algorithm": name,
        "score_column": score_col,
        "n": int(len(d)),
        "positives": int(y.sum()),
        "positive_rate": float(y.mean()),
        "AUPRC": float(average_precision_score(y, s)) if y.nunique() == 2 else np.nan,
        "AUROC": float(roc_auc_score(y, s)) if y.nunique() == 2 else np.nan,
    }
    ordered = d.assign(_score=s).sort_values("_score", ascending=False)
    for k in [3, 5, 10, 15, 25]:
        top = ordered.head(k)
        row[f"top{k}_hits"] = int(top["label"].sum())
        row[f"top{k}_precision"] = float(top["label"].mean())
    return row


def main() -> None:
    locked = pd.read_csv(BASE / "public_industrial_locked_testset_v0.tsv", sep="\t")
    im = pd.read_csv(BASE / "locked_testset_bigmhc_im_predictions.csv")
    el = pd.read_csv(BASE / "locked_testset_bigmhc_el_predictions.csv")
    if "candidate_id" not in im.columns:
        im.insert(0, "candidate_id", locked["industrial_candidate_id"].to_numpy())
    if "candidate_id" not in el.columns:
        el.insert(0, "candidate_id", locked["industrial_candidate_id"].to_numpy())
    df = locked.merge(
        im[["candidate_id", "BigMHC_IM"]],
        left_on="industrial_candidate_id",
        right_on="candidate_id",
        how="left",
    ).drop(columns=["candidate_id"])
    df = df.merge(
        el[["candidate_id", "BigMHC_EL"]],
        left_on="industrial_candidate_id",
        right_on="candidate_id",
        how="left",
    ).drop(columns=["candidate_id"])

    df["bigmhc_im_norm"] = minmax(df["BigMHC_IM"]).fillna(0.0)
    df["bigmhc_el_norm"] = minmax(df["BigMHC_EL"]).fillna(0.0)
    df["peptide_quality_proxy"] = df["peptide"].map(peptide_quality_proxy)
    df["assay_endpoint_is_tetramer"] = df["assay_endpoint_inferred"].astype(str).eq("tetramer_or_response").astype(float)

    # Public-only minimal adapters. They intentionally do not use labels, source id, or response text.
    df["cross_industrial_minimal_score"] = (
        0.55 * df["bigmhc_im_norm"]
        + 0.25 * df["bigmhc_el_norm"]
        + 0.20 * df["peptide_quality_proxy"]
    )
    df["cross_claimsafe_minimal_score"] = (
        0.45 * df["bigmhc_im_norm"]
        + 0.20 * df["bigmhc_el_norm"]
        + 0.25 * df["peptide_quality_proxy"]
        + 0.10 * (1.0 - df["assay_endpoint_is_tetramer"] * 0.0)
    )
    # KG-GA-lite keeps the full KG-GA idea but uses only present modalities.
    # Missing TCR/MD/impact/foreignness are not imputed as signal.
    df["kg_ga_public_minimal_score"] = (
        0.70 * df["cross_industrial_minimal_score"]
        + 0.20 * df["bigmhc_im_norm"]
        + 0.10 * df["peptide_quality_proxy"]
    )

    metrics = pd.DataFrame(
        [
            metric_row(df, "BigMHC_IM", "BigMHC_IM"),
            metric_row(df, "BigMHC_EL", "BigMHC_EL"),
            metric_row(df, "cross_industrial_minimal_score", "CROSS_industrial_minimal"),
            metric_row(df, "cross_claimsafe_minimal_score", "CROSS_claimsafe_minimal"),
            metric_row(df, "kg_ga_public_minimal_score", "KG_GA_public_minimal_adapter"),
        ]
    ).sort_values(["AUPRC", "AUROC"], ascending=False)

    df.to_csv(BASE / "public_industrial_locked_testset_cross_adapter_scores_v0.tsv", sep="\t", index=False)
    metrics.to_csv(BASE / "public_industrial_locked_testset_cross_adapter_metrics_v0.tsv", sep="\t", index=False)
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "input": str(BASE / "public_industrial_locked_testset_v0.tsv"),
        "n": int(len(df)),
        "positives": int(df["label"].sum()),
        "best_algorithm": metrics.iloc[0].to_dict(),
        "claim_boundary": "CROSS/KG-GA minimal adapter only; full KG-GA requires TCR/MD/impact/foreignness/CROSS_BMA fields and manual QA",
        "output_files": [
            str(BASE / "public_industrial_locked_testset_cross_adapter_scores_v0.tsv"),
            str(BASE / "public_industrial_locked_testset_cross_adapter_metrics_v0.tsv"),
        ],
    }
    (BASE / "public_industrial_locked_testset_cross_adapter_summary_v0.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
