#!/usr/bin/env python3
"""Source-heldout stress test for CROSS-Neo v0 branches.

This uses the local train pool sources as held-out domains. It is not external
validation, but it is a harsher shortcut check than random internal CV.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from cross_neo_v0_common import INPUT, OUT, ensure_dirs, load_master, metrics, seq_similarity


STRUCT_COLS = [
    "pwm_p2",
    "pwm_pOmega",
    "pwm_mean",
    "pwm_sum",
    "blosum_max",
    "blosum_mean",
    "kd_mean",
    "kd_std",
    "chg_mean",
    "chg_sum",
    "hel_mean",
    "she_mean",
    "pep_len",
    "aromatic_frac",
]


def qkernel(a: np.ndarray, b: np.ndarray, gamma: float = 1.0, chunk: int = 256) -> np.ndarray:
    out = np.zeros((len(a), len(b)), dtype=float)
    for start in range(0, len(a), chunk):
        aa = a[start : start + chunk]
        diff = aa[:, None, :] - b[None, :, :]
        log_k = np.log(np.clip(np.cos(gamma * diff) ** 2, 1e-8, 1.0)).mean(axis=2)
        out[start : start + chunk] = np.exp(log_k)
    return out


def fit_qk(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray) -> np.ndarray:
    rng = np.random.default_rng(20260509)
    if len(train_x) > 800:
        pos = np.where(train_y == 1)[0]
        neg = np.where(train_y == 0)[0]
        n_pos = min(len(pos), 300)
        n_neg = min(len(neg), 800 - n_pos)
        keep = np.concatenate([rng.choice(pos, n_pos, replace=False), rng.choice(neg, n_neg, replace=False)])
        train_x = train_x[keep]
        train_y = train_y[keep]
    scaler = StandardScaler()
    xtr = scaler.fit_transform(train_x)
    xte = scaler.transform(test_x)
    n_comp = min(12, xtr.shape[1], len(xtr) - 1)
    pca = PCA(n_components=n_comp, random_state=20260509)
    xtr = pca.fit_transform(xtr)
    xte = pca.transform(xte)
    clf = SVC(kernel="precomputed", C=1.0, class_weight="balanced")
    clf.fit(qkernel(xtr, xtr), train_y)
    raw = clf.decision_function(qkernel(xte, xtr))
    return 1.0 / (1.0 + np.exp(-raw))


def train_only_retrieval(test: pd.DataFrame, train: pd.DataFrame) -> pd.DataFrame:
    pep_list = train["peptide_mut"].astype(str).tolist()
    pair_set = set(zip(train["peptide_mut"], train["hla"]))
    rows = []
    for _, r in test.iterrows():
        pep = str(r["peptide_mut"])
        hla = str(r["hla"])
        same = train[train["hla"] == hla]
        near = max((seq_similarity(pep, q) for q in pep_list), default=0.0)
        same_near = max((seq_similarity(pep, q) for q in same["peptide_mut"].astype(str)), default=0.0)
        rows.append(
            {
                "sample_id": r["sample_id"],
                "exact_pair": int((pep, hla) in pair_set),
                "near_similarity": near,
                "same_hla_near_similarity": same_near,
                "same_hla_density": len(same),
                "same_hla_pos_rate": float(same["label"].mean()) if len(same) else 0.0,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    ensure_dirs()
    master = load_master()
    train_master = master[master["split"] == "train"].copy().reset_index(drop=True)

    struct = pd.read_csv(INPUT / "wave2/structure_features.tsv", sep="\t")
    cf = np.load(OUT / "counterfactual_embeddings.npy")
    cf_names = pd.read_csv(OUT / "counterfactual_feature_names.tsv", sep="\t")["feature"].tolist()
    cf_idx = pd.read_csv(OUT / "counterfactual_feature_index.tsv", sep="\t")
    cf_df = pd.DataFrame(cf[:, :80], columns=[f"cf_{n}" for n in cf_names[:80]])
    cf_df.insert(0, "sample_id", cf_idx["sample_id"])

    df = train_master[["sample_id", "peptide_mut", "hla", "label", "study"]].merge(
        struct[["peptide", "HLA_norm", *STRUCT_COLS]],
        left_on=["peptide_mut", "hla"],
        right_on=["peptide", "HLA_norm"],
        how="left",
    )
    df = df.merge(cf_df, on="sample_id", how="left")
    feat_cols = [c for c in df.columns if c.startswith("cf_")] + STRUCT_COLS
    for c in feat_cols:
        df[c] = df[c].astype(float).fillna(df[c].median() if df[c].notna().any() else 0.0)

    rows = []
    pred_rows = []
    studies = [s for s, g in df.groupby("study") if g["label"].nunique() == 2 and len(g) >= 30]
    for study in studies:
        test = df[df["study"] == study].copy()
        train = df[df["study"] != study].copy()
        if train["label"].nunique() < 2 or test["label"].nunique() < 2:
            continue
        retr = train_only_retrieval(test, train)
        ytr = train["label"].to_numpy(int)
        yte = test["label"].to_numpy(int)

        # Counterfactual RF branch.
        cf_cols = [c for c in feat_cols if c.startswith("cf_")]
        rf = RandomForestClassifier(
            n_estimators=200,
            max_depth=3,
            min_samples_leaf=5,
            class_weight="balanced",
            random_state=20260509,
            n_jobs=-1,
        )
        rf.fit(train[cf_cols].to_numpy(float), ytr)
        p_cf = rf.predict_proba(test[cf_cols].to_numpy(float))[:, 1]

        # Quantum-kernel fallback on compact local features.
        q_cols = STRUCT_COLS + cf_cols[:24]
        p_qk = fit_qk(train[q_cols].to_numpy(float), ytr, test[q_cols].to_numpy(float))
        p_fuse = 0.5 * p_cf + 0.5 * p_qk

        for method, pred in [
            ("sourceheld_counterfactual_rf", p_cf),
            ("sourceheld_qk_compact_gamma1", p_qk),
            ("sourceheld_prespecified_late_fusion_w0.5", p_fuse),
        ]:
            mm = metrics(yte, pred)
            rows.append({"heldout_study": study, "method": method, **mm})
            for sid, y, p in zip(test["sample_id"], yte, pred):
                pred_rows.append({"heldout_study": study, "method": method, "sample_id": sid, "label": int(y), "score": float(p)})

    out = pd.DataFrame(rows).sort_values(["method", "AUPRC"], ascending=[True, False])
    out.to_csv(OUT / "source_heldout_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(pred_rows).to_csv(OUT / "source_heldout_predictions.tsv", sep="\t", index=False)
    lines = [
        "# Source-Heldout Stress Test",
        "",
        "This is not external validation. It holds out local train-pool sources to test source shortcut sensitivity.",
        "",
        out.to_markdown(index=False),
    ]
    (OUT / "source_heldout_report.md").write_text("\n".join(lines) + "\n")
    print(f"[source-heldout] rows={len(out)} studies={len(studies)}")


if __name__ == "__main__":
    main()
