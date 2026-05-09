#!/usr/bin/env python3
"""Source-bias, prevalence-shift and PU correction diagnostics for CROSS-Neo v1."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

from cross_neo_v1_common import V0, V1, build_base_feature_table, ensure_v1_dirs, fill_by_train_median, load_master, metrics


SOURCES = ["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation", "ITSNdb strict"]


def feature_cols(base: pd.DataFrame) -> list[str]:
    cols = [c for c in base.columns if c.startswith("cf_")]
    cols += [c for c in ["GP_quantum", "VQC", "W7A_QK_only", "W7A_full", "mean_pLDDT_peptide", "structure_missing", "structure_low_confidence", "peptide_len"] if c in base.columns]
    return cols


def train_rf(train: pd.DataFrame, test: pd.DataFrame, cols: list[str], sample_weight=None) -> np.ndarray:
    if train["label"].nunique() < 2:
        return np.full(len(test), train["label"].mean() if len(train) else 0.5)
    xtr, xte = fill_by_train_median(train, test, cols)
    clf = RandomForestClassifier(
        n_estimators=220,
        max_depth=3,
        min_samples_leaf=5,
        class_weight="balanced",
        random_state=20260509,
        n_jobs=-1,
    )
    clf.fit(xtr, train["label"].to_numpy(int), sample_weight=sample_weight)
    return clf.predict_proba(xte)[:, 1]


def source_weights(train: pd.DataFrame) -> np.ndarray:
    src_counts = train["study"].value_counts()
    w = train["study"].map(lambda s: 1.0 / max(1, src_counts.get(s, 1))).to_numpy(float)
    return w * len(w) / w.sum()


def pu_weights(train: pd.DataFrame) -> np.ndarray:
    w = np.ones(len(train), dtype=float)
    low_capture = train["study"].isin(["NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"])
    w[(train["label"].to_numpy(int) == 0) & low_capture.to_numpy()] = 0.25
    w[(train["label"].to_numpy(int) == 1)] = 1.5
    return w


def main() -> None:
    ensure_v1_dirs()
    master = load_master()
    base = build_base_feature_table(master)
    cols = feature_cols(base)
    strict_ids = set(master.loc[master["strict_set_flag"].astype(bool), "sample_id"])
    diag_rows = []
    for source, sub in master.groupby("study"):
        if source.startswith("ITSNdb"):
            continue
        b = base[base["sample_id"].isin(sub["sample_id"])].copy()
        diag_rows.append(
            {
                "source": source,
                "n": len(sub),
                "n_pos": int(sub["label"].sum()),
                "prevalence": float(sub["label"].mean()),
                "median_peptide_len": float(sub["peptide_mut"].astype(str).str.len().median()),
                "hla_nunique": int(sub["hla"].nunique()),
                "hla_top_fraction": float(sub["hla"].value_counts(normalize=True).iloc[0]),
                "wt_available_rate": float(sub["peptide_wt"].fillna("").astype(str).str.len().gt(0).mean()),
                "source_protein_available_rate": float(sub["source_protein"].fillna("").astype(str).str.len().gt(0).mean()),
                "mean_counterfactual_feature": float(b[[c for c in cols if c.startswith("cf_")]].mean(numeric_only=True).mean()),
                "mean_qk_feature": float(b[[c for c in ["GP_quantum", "VQC", "W7A_QK_only", "W7A_full"] if c in b]].mean(numeric_only=True).mean()),
            }
        )
    strict = master[master["sample_id"].isin(strict_ids)]
    diag_rows.append(
        {
            "source": "ITSNdb strict",
            "n": len(strict),
            "n_pos": int(strict["label"].sum()),
            "prevalence": float(strict["label"].mean()),
            "median_peptide_len": float(strict["peptide_mut"].astype(str).str.len().median()),
            "hla_nunique": int(strict["hla"].nunique()),
            "hla_top_fraction": float(strict["hla"].value_counts(normalize=True).iloc[0]),
            "wt_available_rate": float(strict["peptide_wt"].fillna("").astype(str).str.len().gt(0).mean()),
            "source_protein_available_rate": float(strict["source_protein"].fillna("").astype(str).str.len().gt(0).mean()),
        }
    )
    pd.DataFrame(diag_rows).to_csv(V1 / "source_shift_diagnostics.tsv", sep="\t", index=False)

    # Source predictability sanity check.
    pool = base[base["study"].isin(["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"])].copy()
    pred_rows = []
    if pool["study"].nunique() >= 2:
        tmp = pool[["sample_id", "study", *cols]].copy()
        for c in cols:
            tmp[c] = pd.to_numeric(tmp[c], errors="coerce").fillna(tmp[c].median() if tmp[c].notna().any() else 0.0)
        y = tmp["study"].astype(str)
        n_splits = min(5, y.value_counts().min())
        if n_splits >= 2:
            clf = RandomForestClassifier(n_estimators=160, max_depth=4, min_samples_leaf=8, random_state=20260509, n_jobs=-1)
            cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=20260509)
            yp = cross_val_predict(clf, tmp[cols].to_numpy(float), y, cv=cv)
            pred_rows.append({"task": "predict_source_from_allowed_features", "n": len(y), "accuracy": accuracy_score(y, yp), "macro_f1": f1_score(y, yp, average="macro")})
    pd.DataFrame(pred_rows).to_csv(V1 / "source_predictability.tsv", sep="\t", index=False)

    # Source-heldout corrected training on local train pool.
    train_pool = base[base["study"].isin(["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"])].copy()
    metric_rows, pred_out = [], []
    for heldout, test in train_pool.groupby("study"):
        if test["label"].nunique() < 2 or len(test) < 20:
            continue
        train = train_pool[train_pool["study"] != heldout].copy()
        strategies = {
            "source_balanced_rf": source_weights(train),
            "pu_weighted_rf": pu_weights(train),
            "source_balanced_plus_pu_rf": source_weights(train) * pu_weights(train),
        }
        for name, w in strategies.items():
            score = train_rf(train, test, cols, sample_weight=w)
            # Prior correction uses train prevalence only; it affects cross-source pooled calibration, not within-source rank.
            train_prev = np.clip(train["label"].mean(), 1e-4, 1 - 1e-4)
            global_prev = np.clip(train_pool["label"].mean(), 1e-4, 1 - 1e-4)
            logit = np.log(np.clip(score, 1e-5, 1 - 1e-5) / np.clip(1 - score, 1e-5, 1))
            prior_score = 1 / (1 + np.exp(-(logit + np.log(train_prev / (1 - train_prev)) - np.log(global_prev / (1 - global_prev)))))
            for variant, pred in [(name, score), (name + "_train_prior_calibrated", prior_score)]:
                mm = metrics(test["label"].to_numpy(int), pred)
                metric_rows.append({"heldout_study": heldout, "model": variant, **mm})
                for sid, y, p in zip(test["sample_id"], test["label"], pred):
                    pred_out.append({"heldout_study": heldout, "model": variant, "sample_id": sid, "label": int(y), "score": float(p)})
    met = pd.DataFrame(metric_rows).sort_values(["heldout_study", "AUPRC"], ascending=[True, False])
    met.to_csv(V1 / "source_balanced_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(pred_out).to_csv(V1 / "source_balanced_predictions.tsv", sep="\t", index=False)
    pu = met[met["model"].str.contains("pu", na=False)].copy()
    pu.to_csv(V1 / "pu_corrected_metrics.tsv", sep="\t", index=False)

    v0_source = pd.read_csv(V0 / "source_heldout_metrics.tsv", sep="\t") if (V0 / "source_heldout_metrics.tsv").exists() else pd.DataFrame()
    lines = [
        "# CROSS-Neo v1 Source Bias Report",
        "",
        "This is a source-heldout stress analysis, not external validation.",
        "",
        "## Source Diagnostics",
        pd.DataFrame(diag_rows).to_markdown(index=False),
        "",
        "## Source Predictability",
        pd.DataFrame(pred_rows).to_markdown(index=False) if pred_rows else "Not enough sources for CV source prediction.",
        "",
        "## Corrected Source-Heldout Metrics",
        met.sort_values("AUPRC", ascending=False).head(18).to_markdown(index=False) if len(met) else "No corrected source-heldout metrics.",
        "",
        "## v0 Source-Heldout Reference",
        v0_source.head(12).to_markdown(index=False) if len(v0_source) else "Missing.",
        "",
        "Interpretation: source labels are highly imbalanced and source identity is partly predictable from allowed features. Any source-heldout top-k collapse should be treated as assay/source shift unless corrected models show stable top-k rescue.",
    ]
    (V1 / "source_bias_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v1-source] corrected_metrics={len(met)}")


if __name__ == "__main__":
    main()
