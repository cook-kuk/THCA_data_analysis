#!/usr/bin/env python3
"""BioDarwin RF-anchor expansion v1.

Train the in-house RF biophys baseline on the master pool, score the BioDarwin
academic/industrial candidate tables, and compare:

1. RF only
2. simple heuristic RF + BigMHC + public-anchor mix
3. logistic stacking model
4. current BioDarwin mode bank baseline

This is an experiment driver, not a claim generator.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
BIO = ROOT / "project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10"
WAVE = ROOT / "project/results/p_neo_bayesian_2026_05_09"
OUT = ROOT / "project/results/biodarwin_rf_anchor_expansion_2026_05_11"
OUT.mkdir(parents=True, exist_ok=True)

import sys
sys.path.insert(0, str(ROOT / "project/results/cancer_vaccine_robustness_2026_05_09"))
from _common import biophys, build_hla_onehot_factory, fit_rf, normalize_hla  # noqa: E402


def safe_metrics(y_true, y_score):
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    if len(np.unique(y_true)) < 2:
        return {
            "n": int(len(y_true)),
            "positives": int(y_true.sum()),
            "positive_rate": float(y_true.mean()) if len(y_true) else np.nan,
            "AUPRC": np.nan,
            "AUROC": np.nan,
            "top5_precision": np.nan,
            "top5_hits": np.nan,
            "top10_precision": np.nan,
            "top10_hits": np.nan,
            "top24_precision": np.nan,
            "top24_hits": np.nan,
        }
    order = np.argsort(-y_score)
    out = {
        "n": int(len(y_true)),
        "positives": int(y_true.sum()),
        "positive_rate": float(y_true.mean()),
        "AUPRC": float(average_precision_score(y_true, y_score)),
        "AUROC": float(roc_auc_score(y_true, y_score)),
    }
    for k in [5, 10, 24]:
        kk = min(k, len(y_true))
        hits = int(y_true[order[:kk]].sum())
        out[f"top{k}_precision"] = float(hits / kk)
        out[f"top{k}_hits"] = hits
    return out


def score_rf(train_df, eval_df):
    train_df = train_df.copy()
    eval_df = eval_df.copy()
    train_df["hla_norm"] = train_df["hla_allele_4digit"].apply(normalize_hla)
    eval_df["hla_norm"] = eval_df["hla_allele_4digit"].apply(normalize_hla)
    train_df = train_df.dropna(subset=["hla_norm"]).reset_index(drop=True)
    eval_df = eval_df.dropna(subset=["hla_norm"]).reset_index(drop=True)
    encode, _, _ = build_hla_onehot_factory(train_df["hla_norm"].tolist())

    def feats(df):
        bp = biophys(df["peptide"].tolist())
        hla = encode(df["hla_norm"].tolist())
        return np.hstack([bp, hla]).astype(np.float32)

    Xtr = feats(train_df)
    ytr = train_df["label"].astype(int).to_numpy()
    rf = fit_rf(Xtr, ytr)
    return rf, feats(eval_df), eval_df


def normalize01(x):
    x = np.asarray(x, dtype=float)
    lo = np.nanmin(x)
    hi = np.nanmax(x)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.zeros_like(x, dtype=float)
    return (x - lo) / (hi - lo)


def enrich(df, rf_score):
    out = df.copy().reset_index(drop=True)
    out["rf_biophys_score"] = rf_score
    out["rf_biophys_norm"] = normalize01(rf_score)
    out["rf_bigmhc_gap"] = np.abs(out["rf_biophys_norm"] - out["bigmhc_im_norm"])
    out["rf_public_anchor_gap"] = np.abs(out["rf_biophys_norm"] - out["biodarwin_public_anchor_v3_score"])
    out["rf_anchor_heuristic_v1"] = (
        0.45 * out["rf_biophys_norm"]
        + 0.25 * out["bigmhc_im_norm"]
        + 0.20 * out["biodarwin_public_anchor_v3_score"]
        + 0.10 * (1.0 - out["rf_bigmhc_gap"])
    )
    return out


def main():
    t0 = time.time()
    summary = {"generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"), "out_dir": str(OUT)}

    train_pool = pd.read_csv(WAVE / "bundle.tsv", sep="\t")
    train_pool = train_pool[train_pool["split"] == "train"].dropna(subset=["HLA_norm"]).copy()
    train_pool["hla_allele_4digit"] = train_pool["HLA_norm"]
    train_pool["peptide"] = train_pool["peptide"].astype(str).str.upper()

    mb = BIO / "mode_bank_ensemble_v4"
    academic = pd.read_csv(mb / "biodarwin_mode_bank_academic_scores.tsv", sep="\t")
    industrial = pd.read_csv(mb / "biodarwin_mode_bank_industrial_scores.tsv", sep="\t")
    academic["hla_allele_4digit"] = academic["hla_allele_4digit"].astype(str)
    industrial["hla_allele_4digit"] = industrial["hla_allele"].astype(str)

    rf, X_acad, acad_eval = score_rf(train_pool, academic)
    _, X_ind, ind_eval = score_rf(train_pool, industrial)

    acad = enrich(acad_eval, rf.predict_proba(X_acad)[:, 1])
    ind = enrich(ind_eval, rf.predict_proba(X_ind)[:, 1])

    # train a simple stacker only on academic_train_sources to avoid touching the validation rows
    train_mask = acad["split_biodarwin"].eq("academic_train_sources")
    stack_features = ["rf_biophys_norm", "bigmhc_im_norm", "biodarwin_public_anchor_v3_score", "rf_bigmhc_gap", "seq_cytotoxic_prior", "seq_helper_prior"]
    X_stack = acad.loc[train_mask, stack_features].to_numpy(dtype=float)
    y_stack = acad.loc[train_mask, "label"].astype(int).to_numpy()
    stacker = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(class_weight="balanced", max_iter=500, random_state=42)),
    ])
    stacker.fit(X_stack, y_stack)
    acad["rf_anchor_stack_v1"] = stacker.predict_proba(acad[stack_features].to_numpy(dtype=float))[:, 1]
    ind["rf_anchor_stack_v1"] = stacker.predict_proba(ind[stack_features].to_numpy(dtype=float))[:, 1]

    model_cols = {
        "RF_biophys": "rf_biophys_score",
        "RF_anchor_heuristic_v1": "rf_anchor_heuristic_v1",
        "RF_anchor_stack_v1": "rf_anchor_stack_v1",
        "Current_mode_bank_v4": "biodarwin_mode_bank_v4_score",
        "Public_anchor_v3": "biodarwin_public_anchor_v3_score",
        "BigMHC_IM": "bigmhc_im_norm",
    }

    rows = []
    for split_name, df in [("academic_train_sources", acad[acad["split_biodarwin"].eq("academic_train_sources")]),
                           ("academic_validation_like", acad[acad["split_biodarwin"].eq("academic_validation_like")]),
                           ("industrial_locked_v0", ind[ind["split_biodarwin"].eq("industrial_locked_v0")])]:
        y = df["label"].astype(int).to_numpy()
        for model_name, col in model_cols.items():
            met = safe_metrics(y, df[col].to_numpy(dtype=float))
            rows.append({"split": split_name, "algorithm": model_name, **met})

    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "rf_anchor_expansion_metrics.tsv", sep="\t", index=False)
    acad.to_csv(OUT / "rf_anchor_expansion_academic_scores.tsv", sep="\t", index=False)
    ind.to_csv(OUT / "rf_anchor_expansion_industrial_scores.tsv", sep="\t", index=False)

    # focus summary for the user-facing report
    focus = metrics[metrics["split"].isin(["academic_validation_like", "industrial_locked_v0"])].copy()
    pivot = focus.pivot(index="algorithm", columns="split", values=["AUPRC", "AUROC", "top10_hits"])
    pivot.to_csv(OUT / "rf_anchor_expansion_pivot.tsv", sep="\t")

    summary.update({
        "elapsed_s": round(time.time() - t0, 2),
        "n_academic": int(len(acad)),
        "n_industrial": int(len(ind)),
        "train_rows": int(len(train_pool)),
        "train_pos": int(train_pool["label"].sum()),
        "best_validation_like": focus[focus["split"].eq("academic_validation_like")].sort_values("AUPRC", ascending=False).head(1).to_dict(orient="records"),
        "best_industrial": focus[focus["split"].eq("industrial_locked_v0")].sort_values("AUPRC", ascending=False).head(1).to_dict(orient="records"),
        "outputs": {
            "metrics": str(OUT / "rf_anchor_expansion_metrics.tsv"),
            "academic_scores": str(OUT / "rf_anchor_expansion_academic_scores.tsv"),
            "industrial_scores": str(OUT / "rf_anchor_expansion_industrial_scores.tsv"),
            "pivot": str(OUT / "rf_anchor_expansion_pivot.tsv"),
        }
    })

    (OUT / "rf_anchor_expansion_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    report = []
    report.append("# BioDarwin RF anchor expansion v1")
    report.append("")
    report.append(f"Generated: {summary['generated_at']}")
    report.append("")
    report.append("## Best validation-like")
    report.append(pd.DataFrame(summary["best_validation_like"]).to_markdown(index=False))
    report.append("")
    report.append("## Best industrial locked")
    report.append(pd.DataFrame(summary["best_industrial"]).to_markdown(index=False))
    report.append("")
    report.append("## Focus metrics")
    report.append(focus.to_markdown(index=False))
    report.append("")
    report.append("## Claim boundary")
    report.append("- RF is used as an AUPRC anchor and routing feature.")
    report.append("- This is an experiment driver, not a locked SOTA claim.")
    report.append("- Overlap-aware interpretation still applies.")
    (OUT / "RF_ANCHOR_EXPANSION_REPORT.md").write_text("\n".join(report), encoding="utf-8")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
