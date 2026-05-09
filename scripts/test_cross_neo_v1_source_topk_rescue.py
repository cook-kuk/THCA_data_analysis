#!/usr/bin/env python3
"""Leakage-safe source-heldout top-k rescue tests for CROSS-Neo v1."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler

from cross_neo_v1_lockdown_common import OUT, V0, ensure_dirs, load_master, metric_row, peptide_cluster, sigmoid, summarize_predictions


SOURCES = ["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"]


def feature_table(master: pd.DataFrame) -> pd.DataFrame:
    emb = np.load(V0 / "counterfactual_embeddings.npy")
    idx = pd.read_csv(V0 / "counterfactual_feature_index.tsv", sep="\t")
    n = min(80, emb.shape[1])
    df = pd.DataFrame(emb[:, :n], columns=[f"cf_{i:03d}" for i in range(n)])
    df.insert(0, "sample_id", idx["sample_id"].values)
    meta = master[["sample_id", "peptide_mut", "hla", "hla_supertype", "study", "label", "near_peptide_cluster"]].copy()
    meta["peptide_length"] = meta["peptide_mut"].astype(str).str.len()
    return meta.merge(df, on="sample_id", how="left")


def design_matrix(train: pd.DataFrame, test: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, list[str]]:
    cf_cols = [c for c in train.columns if c.startswith("cf_")]
    cat_cols = ["hla_supertype"]
    tr_cat = pd.get_dummies(train[cat_cols].fillna("UNKNOWN").astype(str), prefix=cat_cols)
    te_cat = pd.get_dummies(test[cat_cols].fillna("UNKNOWN").astype(str), prefix=cat_cols)
    tr_cat, te_cat = tr_cat.align(te_cat, join="outer", axis=1, fill_value=0)
    num_cols = cf_cols + ["peptide_length"]
    xtr = np.hstack([train[num_cols].fillna(0).to_numpy(float), tr_cat.to_numpy(float)])
    xte = np.hstack([test[num_cols].fillna(0).to_numpy(float), te_cat.to_numpy(float)])
    scaler = StandardScaler()
    xtr = scaler.fit_transform(xtr)
    xte = scaler.transform(xte)
    return xtr, xte, num_cols + tr_cat.columns.tolist()


def fit_rf(xtr: np.ndarray, ytr: np.ndarray, sample_weight: np.ndarray | None, seed: int) -> RandomForestClassifier:
    clf = RandomForestClassifier(
        n_estimators=260,
        max_depth=4,
        min_samples_leaf=6,
        class_weight="balanced",
        random_state=seed,
        n_jobs=-1,
    )
    clf.fit(xtr, ytr, sample_weight=sample_weight)
    return clf


def train_score_standardize(train_scores: np.ndarray, test_scores: np.ndarray) -> np.ndarray:
    mu = float(np.mean(train_scores))
    sd = float(np.std(train_scores))
    if sd < 1e-6:
        return test_scores
    return sigmoid((test_scores - mu) / sd)


def hla_rank_normalize(train: pd.DataFrame, train_scores: np.ndarray, test: pd.DataFrame, test_scores: np.ndarray) -> np.ndarray:
    train = train.copy()
    test = test.copy()
    train["_score"] = train_scores
    out = []
    global_scores = np.sort(train_scores)
    for (_, r), s in zip(test.iterrows(), test_scores):
        pool = train.loc[train["hla_supertype"].astype(str).eq(str(r.get("hla_supertype"))), "_score"].to_numpy(float)
        if len(pool) < 15:
            pool = global_scores
        out.append(float((np.asarray(pool) <= s).mean()))
    return np.asarray(out, dtype=float)


def diversify_scores(test: pd.DataFrame, scores: np.ndarray, hla_limit: int = 3, cluster_limit: int = 2) -> np.ndarray:
    d = test[["sample_id", "hla", "near_peptide_cluster"]].copy()
    d["score"] = scores
    d = d.sort_values("score", ascending=False).reset_index(drop=True)
    selected = []
    hla_counts: dict[str, int] = {}
    cluster_counts: dict[str, int] = {}
    for _, r in d.iterrows():
        hla = str(r["hla"])
        cluster = str(r["near_peptide_cluster"])
        penalty = 0.0
        if hla_counts.get(hla, 0) >= hla_limit:
            penalty += 0.10 * (hla_counts[hla] - hla_limit + 1)
        if cluster_counts.get(cluster, 0) >= cluster_limit:
            penalty += 0.12 * (cluster_counts[cluster] - cluster_limit + 1)
        selected.append((r["sample_id"], float(r["score"] - penalty)))
        hla_counts[hla] = hla_counts.get(hla, 0) + 1
        cluster_counts[cluster] = cluster_counts.get(cluster, 0) + 1
    return test["sample_id"].map(dict(selected)).to_numpy(float)


def ood_abstain_scores(train: pd.DataFrame, test: pd.DataFrame, scores: np.ndarray) -> np.ndarray:
    train_hla = set(train["hla"].astype(str))
    train_cluster = set(train["near_peptide_cluster"].astype(str))
    out = scores.copy()
    for i, (_, r) in enumerate(test.iterrows()):
        penalty = 0.0
        if str(r["hla"]) not in train_hla:
            penalty += 0.20
        if str(r["near_peptide_cluster"]) not in train_cluster:
            penalty += 0.15
        if int(r["peptide_length"]) not in set(train["peptide_length"].astype(int)):
            penalty += 0.08
        out[i] = max(0.0, out[i] - penalty)
    return out


def main() -> None:
    ensure_dirs()
    master = load_master()
    feat = feature_table(master)
    train_pool = feat[feat["study"].isin(SOURCES)].copy()
    existing = pd.read_csv(V0 / "source_heldout_predictions.tsv", sep="\t")
    qk_existing = (
        existing[existing["method"].eq("sourceheld_qk_compact_gamma1")][["heldout_study", "sample_id", "score"]]
        .groupby(["heldout_study", "sample_id"], as_index=False)["score"]
        .mean()
        .rename(columns={"score": "qk_score"})
    )

    pred_rows = []
    for heldout in SOURCES:
        test = train_pool[train_pool["study"].eq(heldout)].copy().reset_index(drop=True)
        train = train_pool[~train_pool["study"].eq(heldout)].copy().reset_index(drop=True)
        if test.empty or train["label"].nunique() < 2:
            continue
        xtr, xte, _ = design_matrix(train, test)
        ytr = train["label"].astype(int).to_numpy()
        yte = test["label"].astype(int).to_numpy()

        source_counts = train["study"].value_counts().to_dict()
        source_w = train["study"].map(lambda s: 1.0 / source_counts.get(s, 1)).to_numpy(float)
        source_w = source_w / np.mean(source_w)
        pos_rate = max(1e-6, float(ytr.mean()))
        pos_w = np.where(ytr == 1, min(10.0, 1.0 / pos_rate), 1.0)
        pos_w = pos_w / np.mean(pos_w)
        pu_w = np.where(ytr == 1, 1.0, 0.35)

        variants = {
            "A_score_standardization_train_only": None,
            "B_rank_normalization_hla_supertype_train_only": None,
            "C_source_balanced_training": source_w,
            "D_positive_class_reweighting": pos_w,
            "E_pu_style_conservative_ranker": pu_w,
        }

        base_clf = fit_rf(xtr, ytr, None, 20260509)
        train_base = base_clf.predict_proba(xtr)[:, 1]
        test_base = base_clf.predict_proba(xte)[:, 1]

        scores: dict[str, np.ndarray] = {
            "anchor_rf_retrained_cf": test_base,
            "A_score_standardization_train_only": train_score_standardize(train_base, test_base),
            "B_rank_normalization_hla_supertype_train_only": hla_rank_normalize(train, train_base, test, test_base),
        }
        for method, sw in variants.items():
            if method in scores:
                continue
            clf = fit_rf(xtr, ytr, sw, 20260519 + len(scores))
            scores[method] = clf.predict_proba(xte)[:, 1]

        qk_series = test[["sample_id"]].merge(qk_existing[qk_existing["heldout_study"].eq(heldout)], on="sample_id", how="left")["qk_score"]
        qk = qk_series.to_numpy(float)
        qk = np.where(np.isfinite(qk), qk, test_base)
        fusion = 0.5 * test_base + 0.5 * qk
        scores["H_anchor_rf_plus_qk_compact_w0.5"] = fusion
        scores["F_topk_diversification_hla_cluster"] = diversify_scores(test, fusion)
        scores["G_abstention_extreme_ood"] = ood_abstain_scores(train, test, fusion)

        for method, score in scores.items():
            for sid, y, s in zip(test["sample_id"], yte, score):
                pred_rows.append(
                    {
                        "split_name": f"source_heldout_{heldout}",
                        "fold_id": heldout,
                        "heldout_study": heldout,
                        "method": method,
                        "sample_id": sid,
                        "label": int(y),
                        "score": float(s),
                    }
                )

    pred = pd.DataFrame(pred_rows)
    pred.to_csv(OUT / "source_topk_rescue_predictions.tsv", sep="\t", index=False)
    metrics = summarize_predictions(pred, "method")
    metrics.to_csv(OUT / "source_topk_rescue_metrics.tsv", sep="\t", index=False)

    best = metrics.sort_values(["split_name", "top10_precision", "AUPRC"], ascending=[True, False, False]).groupby("split_name").head(3)
    lines = [
        "# CROSS-Neo v1 Source Top-k Rescue Report",
        "",
        "All variants hold out the target source during fitting. Heldout source labels are used only for final evaluation.",
        "",
        "Methods A/B are train-score normalization/rank normalization, C/D/E are anchor-RF weighting variants, F is label-free top-k diversification, G is OOD score suppression, and H is fixed anchor+QK compact fusion.",
        "",
        "## Best Rows Per Source-Heldout Split",
        "",
        best.to_markdown(index=False) if len(best) else "No rescue rows.",
        "",
        "## Full Metrics",
        "",
        metrics.sort_values(["split_name", "method"]).to_markdown(index=False) if len(metrics) else "No metrics.",
    ]
    (OUT / "source_topk_rescue_report.md").write_text("\n".join(lines) + "\n")
    print(f"[source-topk-rescue] predictions={len(pred)} metrics={len(metrics)}")


if __name__ == "__main__":
    main()
