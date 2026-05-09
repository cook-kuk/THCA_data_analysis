#!/usr/bin/env python3
"""Evaluate strict no-TCR/self-exact ESMFold proxy features."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parent
CURATION = ROOT.parent / "curation_2026_05_09"

FEATURES = [
    "mean_pLDDT_peptide",
    "min_pLDDT_peptide",
    "mean_pLDDT_HLA",
    "anchor_pLDDT",
    "interface_contacts_8A",
    "interface_contacts_10A",
    "n_buried_residues_8A",
    "mean_min_pep_to_hla_CA_dist",
    "max_min_pep_to_hla_CA_dist",
    "radius_of_gyration_peptide",
    "end_to_end_CA_dist",
    "peptide_helicity_proxy",
]

BASELINES = [
    "Wave8_TCR_SelfSim_full",
    "Wave8_TCR_SelfSim_no_exact",
    "Wave8_TCR_motif_only",
    "Wave8_TCR_only",
    "MHCflurry",
    "Structure_LR",
    "BigMHC_IM",
    "PRIME",
    "NetMHCpan_4.1",
    "ESM2_Bayesian",
]


def metric_ci(y: np.ndarray, score: np.ndarray, metric: str, n_boot: int = 400) -> tuple[float, float, float]:
    rng = np.random.default_rng(7)
    mask = np.isfinite(score)
    y = y[mask]
    score = score[mask]
    if y.min() == y.max():
        return np.nan, np.nan, np.nan
    fn = roc_auc_score if metric == "auroc" else average_precision_score
    point = float(fn(y, score))
    vals: list[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if y[idx].min() == y[idx].max():
            continue
        vals.append(float(fn(y[idx], score[idx])))
    lo, hi = np.quantile(vals, [0.025, 0.975])
    return point, float(lo), float(hi)


def fit_predict_cv(x: np.ndarray, y: np.ndarray, repeated: bool = False) -> np.ndarray:
    clf = make_pipeline(
        StandardScaler(),
        LogisticRegression(C=0.5, class_weight="balanced", max_iter=5000, solver="liblinear"),
    )
    cv = (
        RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=17)
        if repeated
        else StratifiedKFold(n_splits=5, shuffle=True, random_state=17)
    )
    if repeated:
        pred = np.zeros(len(y), dtype=float)
        counts = np.zeros(len(y), dtype=float)
        for train_idx, test_idx in cv.split(x, y):
            clf.fit(x[train_idx], y[train_idx])
            pred[test_idx] += clf.predict_proba(x[test_idx])[:, 1]
            counts[test_idx] += 1
        return pred / counts
    return cross_val_predict(clf, x, y, cv=cv, method="predict_proba")[:, 1]


def add_metric(rows: list[dict], name: str, y: np.ndarray, score: np.ndarray, category: str, note: str = "") -> None:
    auroc, auroc_lo, auroc_hi = metric_ci(y, score, "auroc")
    auprc, auprc_lo, auprc_hi = metric_ci(y, score, "auprc")
    mask = np.isfinite(score)
    rows.append(
        {
            "method": name,
            "category": category,
            "n": int(mask.sum()),
            "n_pos": int(y[mask].sum()),
            "auroc": auroc,
            "auroc_ci_lo": auroc_lo,
            "auroc_ci_hi": auroc_hi,
            "auprc": auprc,
            "auprc_ci_lo": auprc_lo,
            "auprc_ci_hi": auprc_hi,
            "note": note,
        }
    )


def main() -> None:
    strict = pd.read_csv(ROOT / "strict_esmfold_features.tsv", sep="\t")
    wide = pd.read_csv(CURATION / "curated_predictions_itsndb_wide.tsv", sep="\t")
    merged = strict.merge(
        wide[["peptide", "hla", "label", *[c for c in BASELINES if c in wide.columns]]],
        on=["peptide", "hla", "label"],
        how="left",
        validate="one_to_one",
    )
    if len(merged) != len(strict):
        raise RuntimeError(f"merge row mismatch: {len(strict)} -> {len(merged)}")

    y = merged["label"].to_numpy(int)
    rows: list[dict] = []

    for method in BASELINES:
        if method not in merged:
            continue
        add_metric(rows, method, y, merged[method].to_numpy(float), "existing_method")

    for feat in FEATURES:
        score = merged[feat].to_numpy(float)
        raw_auc = roc_auc_score(y, score)
        inv_auc = roc_auc_score(y, -score)
        direction = "raw" if raw_auc >= inv_auc else "inverse"
        add_metric(
            rows,
            f"ESMFold_univariate__{feat}",
            y,
            score if direction == "raw" else -score,
            "structure_univariate",
            f"direction={direction}; raw_auroc={raw_auc:.3f}",
        )

    x3d = merged[FEATURES].to_numpy(float)
    pred_3d = fit_predict_cv(x3d, y, repeated=True)
    add_metric(rows, "ESMFold_3D_LR_repeated5x5_CV", y, pred_3d, "structure_cv", "strict-set internal CV")

    combo_methods = ["Structure_LR", "MHCflurry", "Wave8_TCR_SelfSim_no_exact", "Wave8_TCR_SelfSim_full"]
    for method in combo_methods:
        if method not in merged:
            continue
        mask = merged[method].notna().to_numpy()
        x_combo = np.column_stack([merged.loc[mask, FEATURES].to_numpy(float), merged.loc[mask, method].to_numpy(float)])
        pred = fit_predict_cv(x_combo, y[mask], repeated=True)
        add_metric(
            rows,
            f"{method}+ESMFold_3D_LR_repeated5x5_CV",
            y[mask],
            pred,
            "combo_cv",
            "strict-set internal CV; includes existing score as one feature",
        )

    out = pd.DataFrame(rows).sort_values(["auroc", "auprc"], ascending=False)
    out.to_csv(ROOT / "strict_esmfold_method_comparison.tsv", sep="\t", index=False)
    merged.assign(ESMFold_3D_CV_score=pred_3d).to_csv(
        ROOT / "strict_esmfold_merged_predictions.tsv", sep="\t", index=False
    )

    feature_rows = [
        {
            "feature": feat,
            "mean": float(merged[feat].mean()),
            "sd": float(merged[feat].std()),
            "positive_mean": float(merged.loc[merged["label"] == 1, feat].mean()),
            "negative_mean": float(merged.loc[merged["label"] == 0, feat].mean()),
        }
        for feat in FEATURES
    ]
    pd.DataFrame(feature_rows).to_csv(ROOT / "strict_esmfold_feature_summary.tsv", sep="\t", index=False)

    top = out.head(12)
    lines = [
        "# Wave8B Strict ESMFold Structure Summary",
        "",
        f"Strict set: n={len(merged)}, positives={int(y.sum())}. Filters: no in-house overlap, no TCR exact hit, no self exact hit, HLA pseudo-sequence available.",
        "",
        "## Top Same-Row Results",
        "",
        "| method | category | n / pos | AUROC | AUPRC | note |",
        "|---|---|---:|---:|---:|---|",
    ]
    for _, r in top.iterrows():
        lines.append(
            f"| {r['method']} | {r['category']} | {int(r['n'])} / {int(r['n_pos'])} | "
            f"{r['auroc']:.3f} [{r['auroc_ci_lo']:.3f}, {r['auroc_ci_hi']:.3f}] | "
            f"{r['auprc']:.3f} [{r['auprc_ci_lo']:.3f}, {r['auprc_ci_hi']:.3f}] | {r['note']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- The ESMFold branch is a pseudo-complex proxy: peptide + `GGGGS` + 34-aa HLA pseudo-sequence folded as one chain.",
            "- It is useful as a fast structure-derived feature generator, not as a physical pMHC structure claim.",
            "- Internal CV on n=89 is hypothesis-generating. The paper-grade claim still requires a leakage-controlled external/time-split benchmark.",
        ]
    )
    (ROOT / "WAVE8B_STRICT_STRUCTURE_SUMMARY.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
