#!/usr/bin/env python3
"""Prediction-only audit for the UNI-final TCGA-THCA CLAM OOF table.

This does not retrain UNI under LOTO. It records the available OOF-level
checks so the 0.874 UNI result is not left completely unaudited.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold


ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/"
            "p2_image_dm1_v2_foundation_clam_2026_05_07")
META = Path("/home/seungho/personal/THCA_data_analysis/project/metadata/"
            "sample_master_v3.tsv")
OUT = ROOT / "analysis_supp" / "audit_uni_predictions"
OUT.mkdir(parents=True, exist_ok=True)


def load_merged() -> pd.DataFrame:
    pred = pd.read_csv(ROOT / "phase2_tcga_clam_UNI"
                       / "clam_per_slide_predictions.tsv", sep="\t")
    man = pd.read_csv(ROOT / "phase2_tcga_clam"
                      / "slide_manifest.tsv", sep="\t")
    df = pred.merge(
        man[["file_id", "case_id", "submitter_id"]],
        left_on="slide", right_on="file_id", how="left"
    )
    meta = pd.read_csv(META, sep="\t", low_memory=False)
    meta = meta[meta["dataset"] == "TCGA-THCA"].copy()
    meta["case_short"] = meta["sample_id"].str.extract(r"^(TCGA-[^-]+-[^-]+)")
    meta_case = (meta.sort_values("sample_id")
                     .drop_duplicates("case_short", keep="first")
                     [["case_short", "histology_subtype", "molecular_subtype",
                       "sex", "age", "ajcc_stage_group"]])
    df = df.merge(meta_case, left_on="submitter_id",
                  right_on="case_short", how="left")
    df["tss"] = df["submitter_id"].str.extract(r"TCGA-([^-]+)-")
    return df


def subgroup_auc(df: pd.DataFrame, col: str, val: str) -> dict:
    sub = df[df[col] == val]
    auc = np.nan
    if len(sub) and sub["label"].nunique() > 1:
        auc = float(roc_auc_score(sub["label"], sub["prob_DM1"]))
    return {
        "group": f"{col}={val}",
        "n": int(len(sub)),
        "n_pos": int(sub["label"].sum()) if len(sub) else 0,
        "auc": auc,
    }


def clinical_oof_auc(df: pd.DataFrame) -> float:
    sub = df.dropna(subset=["histology_subtype", "sex", "tss"]).copy()
    sub = sub[sub["histology_subtype"].isin(["cPTC", "FVPTC"])]
    x = pd.get_dummies(sub[["histology_subtype", "sex", "tss"]],
                       drop_first=True).astype(float).values
    y = sub["label"].values
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    oof = np.full(len(y), np.nan)
    for tr, va in kf.split(x):
        lr = LogisticRegression(max_iter=1000, C=1.0)
        lr.fit(x[tr], y[tr])
        oof[va] = lr.predict_proba(x[va])[:, 1]
    return float(roc_auc_score(y, oof))


def tss_only_ras_auc(df: pd.DataFrame) -> float:
    ras = df[df["molecular_subtype"] == "RAS_like"].copy()
    x = pd.get_dummies(ras[["tss"]], drop_first=False).astype(float).values
    y = ras["label"].values
    lr = LogisticRegression(C=1e6, max_iter=1000)
    lr.fit(x, y)
    return float(roc_auc_score(y, lr.predict_proba(x)[:, 1]))


def main() -> None:
    df = load_merged()
    df.to_csv(OUT / "uni_merged_preds_subgroups.tsv", sep="\t", index=False)
    ras = df[df["molecular_subtype"] == "RAS_like"].copy()
    ras.sort_values("prob_DM1").to_csv(
        OUT / "uni_RAS_like_16_predictions.tsv", sep="\t", index=False
    )

    summary = {
        "n_slides": int(len(df)),
        "n_pos": int(df["label"].sum()),
        "n_neg": int(len(df) - df["label"].sum()),
        "overall_auc": float(roc_auc_score(df["label"], df["prob_DM1"])),
        "subgroup_auc": [
            subgroup_auc(df, "molecular_subtype", "RAS_like"),
            subgroup_auc(df, "histology_subtype", "FVPTC"),
            subgroup_auc(df, "molecular_subtype", "BRAF_like"),
            subgroup_auc(df, "histology_subtype", "cPTC"),
        ],
        "sex_auc": [
            subgroup_auc(df, "sex", "Female"),
            subgroup_auc(df, "sex", "Male"),
        ],
        "clinical_only_oof_auc": clinical_oof_auc(df),
        "tss_only_within_raslike_auc": tss_only_ras_auc(df),
        "raslike_tss_counts": (
            ras.groupby("tss")["label"]
               .agg(n="count", n_pos="sum")
               .reset_index()
               .to_dict("records")
        ),
        "note": (
            "Prediction-only audit: no UNI LOTO retraining was run. "
            "TSS-label determinism in the RAS-like subset remains identical "
            "to the ViT-L audit."
        ),
    }
    (OUT / "UNI_PREDICTION_ONLY_AUDIT_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, default=str)
    )

    lines = [
        "# UNI-final prediction-only audit",
        "",
        "This is an OOF-table audit only. It does not replace UNI LOTO retraining.",
        "",
        "| Metric | Value |",
        "|---|---:|",
        f"| n slides | {summary['n_slides']} |",
        f"| overall AUC | {summary['overall_auc']:.3f} |",
        f"| clinical-only OOF AUC | {summary['clinical_only_oof_auc']:.3f} |",
        f"| TSS-only within-RAS-like AUC | {summary['tss_only_within_raslike_auc']:.3f} |",
    ]
    for rec in summary["subgroup_auc"]:
        lines.append(f"| {rec['group']} AUC | {rec['auc']:.3f} |")
    for rec in summary["sex_auc"]:
        lines.append(f"| {rec['group']} AUC | {rec['auc']:.3f} |")
    lines.extend([
        "",
        "## Interpretation",
        "",
        "UNI-final OOF predictions improve overall performance relative to the "
        "ViT-L audit table and do not reproduce the RAS-like AUC=1.000. "
        "However, the RAS-like TSS-label structure is unchanged: DJ/FK contain "
        "the 3 RAS-like DM1 cases and EM contains 10 RAS-like DM2 cases. "
        "Therefore K2 H&E or UNI-specific LOTO retraining is still required "
        "before any center-generalized image-DM1 subgroup claim.",
        "",
    ])
    (OUT / "UNI_PREDICTION_ONLY_AUDIT.md").write_text("\n".join(lines))
    print(f"WROTE {OUT / 'UNI_PREDICTION_ONLY_AUDIT_SUMMARY.json'}")
    print(f"WROTE {OUT / 'UNI_PREDICTION_ONLY_AUDIT.md'}")


if __name__ == "__main__":
    main()
