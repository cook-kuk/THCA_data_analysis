#!/usr/bin/env python3
"""BioDarwin regime router v1.

Hybrid router:
- internal-feature rows -> current mode bank
- external OOD rows -> public anchor + BigMHC with RF disagreement penalty

This is the next step after RF anchor expansion: let RF influence routing,
not dominate the final score.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
MB = ROOT / "project/results/biodarwin_pan_vaccine_ga_rl_2026_05_10/mode_bank_ensemble_v4"
RF = ROOT / "project/results/biodarwin_rf_anchor_expansion_2026_05_11"
OUT = ROOT / "project/results/biodarwin_regime_router_2026_05_11"
OUT.mkdir(parents=True, exist_ok=True)


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


def main():
    t0 = time.time()
    acad = pd.read_csv(MB / "biodarwin_mode_bank_academic_scores.tsv", sep="\t")
    ind = pd.read_csv(MB / "biodarwin_mode_bank_industrial_scores.tsv", sep="\t")
    acad_rf = pd.read_csv(RF / "rf_anchor_expansion_academic_scores.tsv", sep="\t")
    ind_rf = pd.read_csv(RF / "rf_anchor_expansion_industrial_scores.tsv", sep="\t")

    acad = acad.merge(acad_rf[["candidate_id", "rf_biophys_norm", "rf_bigmhc_gap", "rf_anchor_heuristic_v1"]], on="candidate_id", how="left")
    ind = ind.merge(ind_rf[["industrial_candidate_id", "rf_biophys_norm", "rf_bigmhc_gap", "rf_anchor_heuristic_v1"]], on="industrial_candidate_id", how="left")

    for df in (acad, ind):
        for c in ["biodarwin_mode_bank_v4_score", "biodarwin_public_anchor_v3_score", "bigmhc_im_norm", "rf_bigmhc_gap", "rf_anchor_heuristic_v1", "rf_biophys_norm"]:
            if c in df.columns:
                df[c] = df[c].fillna(0.0)

    # External fallback: public anchor + BigMHC with an RF disagreement penalty.
    acad["regime_router_v1_score"] = np.where(
        acad["internal_feature_available"].astype(bool),
        acad["biodarwin_mode_bank_v4_score"],
        0.70 * acad["biodarwin_public_anchor_v3_score"] + 0.30 * acad["bigmhc_im_norm"] - 0.10 * acad["rf_bigmhc_gap"],
    )
    ind["regime_router_v1_score"] = np.where(
        ind["internal_feature_available"].astype(bool),
        ind["biodarwin_mode_bank_v4_score"],
        0.70 * ind["biodarwin_public_anchor_v3_score"] + 0.30 * ind["bigmhc_im_norm"] - 0.10 * ind["rf_bigmhc_gap"],
    )

    rows = []
    for split_name, df in [
        ("academic_train_sources", acad[acad["split_biodarwin"].eq("academic_train_sources")]),
        ("academic_validation_like", acad[acad["split_biodarwin"].eq("academic_validation_like")]),
        ("industrial_locked_v0", ind[ind["split_biodarwin"].eq("industrial_locked_v0")]),
    ]:
        y = df["label"].astype(int).to_numpy()
        for model_name, col in [
            ("Current_mode_bank_v4", "biodarwin_mode_bank_v4_score"),
            ("Public_anchor_v3", "biodarwin_public_anchor_v3_score"),
            ("RF_biophys", "rf_biophys_norm"),
            ("RF_anchor_heuristic_v1", "rf_anchor_heuristic_v1"),
            ("Regime_router_v1", "regime_router_v1_score"),
        ]:
            rows.append({"split": split_name, "algorithm": model_name, **safe_metrics(y, df[col].to_numpy(dtype=float))})

    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "regime_router_metrics.tsv", sep="\t", index=False)
    acad.to_csv(OUT / "regime_router_academic_scores.tsv", sep="\t", index=False)
    ind.to_csv(OUT / "regime_router_industrial_scores.tsv", sep="\t", index=False)
    pivot = metrics[metrics["split"].isin(["academic_validation_like", "industrial_locked_v0"])].pivot(index="algorithm", columns="split", values=["AUPRC", "AUROC", "top10_hits"])
    pivot.to_csv(OUT / "regime_router_pivot.tsv", sep="\t")

    best_val = metrics[metrics["split"].eq("academic_validation_like")].sort_values("AUPRC", ascending=False).head(1).to_dict(orient="records")
    best_ind = metrics[metrics["split"].eq("industrial_locked_v0")].sort_values("AUPRC", ascending=False).head(1).to_dict(orient="records")

    summary = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "elapsed_s": round(time.time() - t0, 2),
        "out_dir": str(OUT),
        "best_validation_like": best_val,
        "best_industrial": best_ind,
        "outputs": {
            "metrics": str(OUT / "regime_router_metrics.tsv"),
            "academic_scores": str(OUT / "regime_router_academic_scores.tsv"),
            "industrial_scores": str(OUT / "regime_router_industrial_scores.tsv"),
            "pivot": str(OUT / "regime_router_pivot.tsv"),
        },
    }
    (OUT / "regime_router_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    report = []
    report.append("# BioDarwin regime router v1")
    report.append("")
    report.append(f"Generated: {summary['generated_at']}")
    report.append("")
    report.append("## Best validation-like")
    report.append(pd.DataFrame(best_val).to_markdown(index=False))
    report.append("")
    report.append("## Best industrial locked")
    report.append(pd.DataFrame(best_ind).to_markdown(index=False))
    report.append("")
    report.append("## Claim boundary")
    report.append("- Internal-feature rows use current mode bank.")
    report.append("- External rows use public anchor + BigMHC with RF disagreement penalty.")
    report.append("- This is a router, not a locked SOTA claim.")
    (OUT / "REGIME_ROUTER_REPORT.md").write_text("\n".join(report), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
