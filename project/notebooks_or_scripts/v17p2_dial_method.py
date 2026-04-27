#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
from inmoose.pycombat import pycombat_norm
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score, davies_bouldin_score, roc_auc_score, silhouette_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from v17p2_common import FIG, LOG, PHASE1_TAB, TAB, V5_CC, dia_auc, identifiability_auc, load_dark_cohort, log_line, zscore_pca, write_json
from v17_common import align_dataset_expr_meta


def combat_blend(train: pd.DataFrame, test: pd.DataFrame, lam: float) -> pd.DataFrame:
    joined = pd.concat([train, test], axis=0)
    batch = np.array(["TCGA"] * len(train) + ["EXT"] * len(test))
    corrected = pycombat_norm(joined.T, batch=batch, covar_mod=None)
    corrected = corrected.T if isinstance(corrected, pd.DataFrame) else pd.DataFrame(np.asarray(corrected).T, index=joined.index, columns=joined.columns)
    return (1.0 - lam) * joined + lam * corrected


def pancancer_table(cancer: str) -> pd.DataFrame:
    base = V5_CC / cancer
    X = np.load(base / "X_combined.npz")["X"]
    Y = pd.read_csv(base / "Y.tsv", sep="\t", header=None)
    B = pd.read_csv(base / "B.tsv", sep="\t", header=None)
    n = min(X.shape[0], len(Y), len(B))
    X = X[:n]
    y = Y.iloc[:n, 0].to_numpy()
    b = B.iloc[:n, 0].astype(str).to_numpy()
    clf = make_pipeline(StandardScaler(), PCA(n_components=min(50, X.shape[0] - 1, X.shape[1]), random_state=42), LogisticRegression(max_iter=3000, random_state=42))
    prob = cross_val_predict(clf, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=42), method="predict_proba")[:, 1]
    Xp = PCA(n_components=min(50, X.shape[0] - 1, X.shape[1]), random_state=42).fit_transform(StandardScaler().fit_transform(X))
    return pd.DataFrame([{
        "cancer": cancer,
        "n": int(n),
        "label_auc": float(roc_auc_score(y, prob)),
        "dia_auc": float(dia_auc(y, prob)),
        "identifiability": float(identifiability_auc(pd.DataFrame(Xp), b)),
    }])


def main() -> None:
    log_line(LOG, "v17p2 layer5 start")
    markers = pd.read_csv(PHASE1_TAB / "dark_matter_cluster_markers.tsv", sep="\t")
    sig = sorted(markers["gene"].unique())[:30]
    tcga_expr, meta = load_dark_cohort()
    y = (meta["v17_dark_cluster"] == "DM2").astype(int).to_numpy()
    ext_expr, ext_meta = align_dataset_expr_meta("GSE27155", genes=sig, tumor_only=True)
    ext_meta = ext_meta[~ext_meta["driver_anchor"].isin(["BRAF", "RAS"])].copy()
    ext_expr = ext_expr.loc[ext_meta["sample_id"]]
    common = [g for g in sig if g in tcga_expr.columns and g in ext_expr.columns]
    centroid_dm1 = tcga_expr.loc[y == 0, common].mean(axis=0)
    centroid_dm2 = tcga_expr.loc[y == 1, common].mean(axis=0)
    d1 = ((ext_expr[common] - centroid_dm1) ** 2).sum(axis=1)
    d2 = ((ext_expr[common] - centroid_dm2) ** 2).sum(axis=1)
    pseudo = (d2 < d1).astype(int).to_numpy()
    if len(np.unique(pseudo)) < 2:
        pseudo[0] = 1 - pseudo[0]

    sens_rows = []
    for lam in np.linspace(0.0, 1.0, 11):
        mixed = combat_blend(tcga_expr[common], ext_expr[common], float(lam))
        tc = mixed.iloc[: len(tcga_expr)]
        ex = mixed.iloc[len(tcga_expr):]
        clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000, random_state=42))
        clf.fit(tc, y)
        prob = clf.predict_proba(ex)[:, 1]
        joined = pd.concat([tc, ex], axis=0)
        ident = identifiability_auc(joined, np.array(["TCGA"] * len(tc) + ["EXT"] * len(ex)))
        sens_rows.append({"lambda": lam, "DIA_AUC": dia_auc(pseudo, prob), "identifiability": ident})
    sens = pd.DataFrame(sens_rows)
    sens.to_csv(TAB / "dial_combat_sensitivity_curve.tsv", sep="\t", index=False)

    pc_parts = []
    for cancer in ["THCA", "LUAD", "COAD"]:
        try:
            pc_parts.append(pancancer_table(cancer))
        except Exception as e:
            log_line(LOG, f"v17p2 layer5 pancancer skip {cancer}: {e}")
    pc_df = pd.concat(pc_parts, ignore_index=True) if pc_parts else pd.DataFrame(columns=["cancer", "n", "label_auc", "dia_auc", "identifiability"])
    pc_df.to_csv(TAB / "dial_pancancer_application.tsv", sep="\t", index=False)

    alt_rows = []
    X, _, _ = zscore_pca(tcga_expr[common], n_components=10)
    km = KMeans(n_clusters=2, random_state=42, n_init=25).fit(X)
    sil = silhouette_score(X, km.labels_)
    db = davies_bouldin_score(X, km.labels_)
    ari = adjusted_rand_score(y, km.labels_)
    alt_rows.append({"dataset": "TCGA_DarkMatter", "metric": "silhouette", "value": float(sil)})
    alt_rows.append({"dataset": "TCGA_DarkMatter", "metric": "davies_bouldin", "value": float(db)})
    alt_rows.append({"dataset": "TCGA_DarkMatter", "metric": "ARI_vs_v17", "value": float(ari)})
    alt_rows.append({"dataset": "TCGA_DarkMatter", "metric": "DIAL_identifiability", "value": float(sens.loc[sens['lambda'] == 1.0, 'identifiability'].iloc[0])})
    alt_rows.append({"dataset": "TCGA_DarkMatter", "metric": "DIA_AUC_ext", "value": float(sens.loc[sens['lambda'] == 1.0, 'DIA_AUC'].iloc[0])})
    alt = pd.DataFrame(alt_rows)
    alt.to_csv(TAB / "dial_vs_alternatives.tsv", sep="\t", index=False)

    sim_rows = []
    rng = np.random.default_rng(42)
    for n in [40, 80, 160]:
        for batch_shift in [0.0, 0.5, 1.0, 1.5, 2.0]:
            for imbalance in [0.5, 0.7, 0.85]:
                n1 = int(n * imbalance)
                n0 = n - n1
                x0 = rng.normal(0, 1, size=(n0, 10))
                x1 = rng.normal(0.8, 1, size=(n1, 10))
                Xs = np.vstack([x0, x1])
                ysim = np.array([0] * n0 + [1] * n1)
                batch = np.array(["A"] * (n // 2) + ["B"] * (n - n // 2))
                Xs[batch == "B"] += batch_shift
                clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, random_state=42))
                prob = cross_val_predict(clf, Xs, ysim, cv=StratifiedKFold(5, shuffle=True, random_state=42), method="predict_proba")[:, 1]
                sim_rows.append({
                    "n": n,
                    "batch_shift": batch_shift,
                    "imbalance": imbalance,
                    "label_auc": float(roc_auc_score(ysim, prob)),
                    "identifiability": float(identifiability_auc(pd.DataFrame(Xs), batch)),
                })
    sim = pd.DataFrame(sim_rows)
    sim.to_csv(TAB / "dial_failure_modes.tsv", sep="\t", index=False)

    px.line(sens, x="lambda", y=["DIA_AUC", "identifiability"], markers=True, title="ComBat sensitivity threshold").write_html(FIG / "dial_combat_threshold.html", include_plotlyjs="cdn")
    if not pc_df.empty:
        px.scatter(pc_df, x="dia_auc", y="cancer", size="n", color="identifiability", title="DIAL pan-cancer feasibility").write_html(FIG / "dial_pancancer_forest.html", include_plotlyjs="cdn")
    else:
        px.scatter(title="Pan-cancer feasibility unavailable").write_html(FIG / "dial_pancancer_forest.html", include_plotlyjs="cdn")
    px.line_polar(alt, r="value", theta="metric", line_close=True, title="DIAL vs alternative metrics").write_html(FIG / "dial_vs_alternatives_radar.html", include_plotlyjs="cdn")
    px.scatter(sim, x="batch_shift", y="identifiability", color="imbalance", size="n", title="DIAL failure-mode simulation").write_html(FIG / "dial_failure_mode_simulation.html", include_plotlyjs="cdn")

    threshold = sens.loc[sens["identifiability"] < 0.5, "lambda"]
    summary = {
        "combat_threshold_lambda_ident_lt_0_5": float(threshold.min()) if not threshold.empty else None,
        "pancancer_tested": int(len(pc_df)),
        "max_pancancer_dia_auc": float(pc_df["dia_auc"].max()) if not pc_df.empty else None,
    }
    write_json(TAB / "dial_method_summary.json", summary)
    log_line(LOG, f"v17p2 layer5 done {summary}")


if __name__ == "__main__":
    main()
