#!/usr/bin/env python3
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from sklearn.cluster import AgglomerativeClustering, KMeans, SpectralClustering
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from sklearn.mixture import GaussianMixture

from v17_idea1_dark_matter import leiden_target_k
from v17p2_common import FIG, LOG, SEED, TAB, dm_binary_labels, load_dark_cohort, log_line, write_json, zscore_pca

try:
    import hdbscan  # type: ignore
except Exception:
    hdbscan = None


N_BOOT = int(os.environ.get("V17P2_BOOTSTRAPS", "300"))
N_PERM = int(os.environ.get("V17P2_PERMUTATIONS", "300"))
ABLATION_REP = int(os.environ.get("V17P2_ABLATION_REP", "40"))
SUBSAMPLE_REP = int(os.environ.get("V17P2_SUBSAMPLE_REP", "60"))


def method_labels(X: np.ndarray) -> dict[str, np.ndarray]:
    out = {
        "Leiden": leiden_target_k(X, 2, SEED),
        "KMeans": KMeans(n_clusters=2, random_state=SEED, n_init=25).fit_predict(X),
        "GMM": GaussianMixture(n_components=2, random_state=SEED).fit_predict(X),
        "Spectral": SpectralClustering(n_clusters=2, random_state=SEED, affinity="nearest_neighbors", assign_labels="kmeans").fit_predict(X),
    }
    if hdbscan is not None:
        lab = hdbscan.HDBSCAN(min_cluster_size=18).fit_predict(X)
        if len(np.unique(lab[lab >= 0])) >= 2:
            if (lab < 0).any():
                fill = pd.Series(lab).replace(-1, pd.Series(lab[lab >= 0]).mode().iloc[0]).to_numpy()
                out["HDBSCAN"] = fill.astype(int)
            else:
                out["HDBSCAN"] = lab.astype(int)
        else:
            out["HDBSCAN"] = AgglomerativeClustering(n_clusters=2, linkage="ward").fit_predict(X)
    else:
        out["HDBSCAN_fallback"] = AgglomerativeClustering(n_clusters=2, linkage="ward").fit_predict(X)
    return out


def bootstrap_consensus(X: np.ndarray, labels_full: np.ndarray, n_boot: int) -> tuple[np.ndarray, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    n = X.shape[0]
    same = np.zeros((n, n), dtype=float)
    seen = np.zeros((n, n), dtype=float)
    for b in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        uniq = np.unique(idx)
        Xi = X[uniq]
        labs = KMeans(n_clusters=2, random_state=SEED + b, n_init=20).fit_predict(Xi)
        pos = {u: i for i, u in enumerate(uniq)}
        for i in uniq:
            for j in uniq:
                seen[i, j] += 1.0
                if labs[pos[i]] == labs[pos[j]]:
                    same[i, j] += 1.0
    consensus = np.divide(same, seen, out=np.zeros_like(same), where=seen > 0)
    rows = []
    for i in range(n):
        peer = np.where(labels_full == labels_full[i])[0]
        peer = peer[peer != i]
        persistence = float(consensus[i, peer].mean()) if len(peer) else np.nan
        rows.append({"sample_id": i, "cluster": int(labels_full[i]), "bootstrap_persistence": persistence})
    return consensus, pd.DataFrame(rows)


def permuted_silhouette(X: np.ndarray, n_perm: int) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    obs_lab = KMeans(n_clusters=2, random_state=SEED, n_init=25).fit_predict(X)
    obs = float(silhouette_score(X, obs_lab))
    for i in range(n_perm):
        Xp = X.copy()
        for j in range(X.shape[1]):
            Xp[:, j] = rng.permutation(Xp[:, j])
        lab = KMeans(n_clusters=2, random_state=SEED + i, n_init=10).fit_predict(Xp)
        rows.append({"iter": i, "metric": "null_silhouette", "value": float(silhouette_score(Xp, lab))})
    null = pd.DataFrame(rows)
    p_emp = (1.0 + float((null["value"] >= obs).sum())) / (len(null) + 1.0)
    out = pd.concat([pd.DataFrame([{"iter": -1, "metric": "observed_silhouette", "value": obs, "empirical_p": p_emp}]), null], ignore_index=True)
    return out


def feature_ablation(expr: pd.DataFrame, full_labels: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    sizes = [30, 40, 50]
    rows = []
    genes = expr.columns.to_list()
    for size in sizes:
        for rep in range(ABLATION_REP):
            keep = sorted(rng.choice(genes, size=size, replace=False).tolist())
            X, _, _ = zscore_pca(expr[keep], n_components=10)
            lab = KMeans(n_clusters=2, random_state=SEED + rep, n_init=15).fit_predict(X)
            rows.append({"panel_size": size, "rep": rep, "ARI_vs_full": adjusted_rand_score(full_labels, lab)})
    imp = []
    for i, gene in enumerate(genes):
        keep = [g for g in genes if g != gene]
        X, _, _ = zscore_pca(expr[keep], n_components=10)
        lab = KMeans(n_clusters=2, random_state=SEED + i, n_init=15).fit_predict(X)
        ari = adjusted_rand_score(full_labels, lab)
        imp.append({"gene": gene, "leave_one_out_ARI": ari, "importance_drop": 1.0 - ari})
    return pd.DataFrame(rows), pd.DataFrame(imp).sort_values("importance_drop", ascending=False)


def subsample_stability(expr: pd.DataFrame, full_labels: np.ndarray) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    n = len(expr)
    for frac in [0.5, 0.7, 0.9]:
        k = max(20, int(round(n * frac)))
        for rep in range(SUBSAMPLE_REP):
            idx = np.sort(rng.choice(n, size=k, replace=False))
            sub = expr.iloc[idx]
            X, _, _ = zscore_pca(sub, n_components=10)
            lab = KMeans(n_clusters=2, random_state=SEED + rep, n_init=15).fit_predict(X)
            rows.append({"subsample_fraction": frac, "rep": rep, "n_sub": k, "ARI_vs_full": adjusted_rand_score(full_labels[idx], lab)})
    return pd.DataFrame(rows)


def main() -> None:
    log_line(LOG, "v17p2 layer1 start")
    expr, meta = load_dark_cohort()
    X, _, _ = zscore_pca(expr, n_components=10)
    full_labels = dm_binary_labels(meta)

    methods = method_labels(X)
    metric_rows = []
    for a, la in methods.items():
        for b, lb in methods.items():
            metric_rows.append({
                "method_a": a,
                "method_b": b,
                "ARI": adjusted_rand_score(la, lb),
                "NMI": normalized_mutual_info_score(la, lb),
            })
    metric_df = pd.DataFrame(metric_rows)
    metric_df.to_csv(TAB / "robustness_5method_ARI.tsv", sep="\t", index=False)

    consensus, persistence = bootstrap_consensus(X, full_labels, N_BOOT)
    persistence["sample_id"] = meta["sample_id"].values
    persistence["cluster_name"] = meta["v17_dark_cluster"].values
    persistence.to_csv(TAB / "bootstrap_persistence.tsv", sep="\t", index=False)

    perm = permuted_silhouette(X, N_PERM)
    perm.to_csv(TAB / "permutation_pvalue.tsv", sep="\t", index=False)

    abl_df, imp_df = feature_ablation(expr, full_labels)
    subs_df = subsample_stability(expr, full_labels)
    feat_out = imp_df.merge(abl_df.groupby("panel_size", as_index=False)["ARI_vs_full"].mean().rename(columns={"ARI_vs_full": "mean_subset_ARI"}), how="cross")
    feat_out.to_csv(TAB / "feature_importance.tsv", sep="\t", index=False)
    subs_df.to_csv(TAB / "sample_subsample_stability.tsv", sep="\t", index=False)

    heat = metric_df.pivot(index="method_a", columns="method_b", values="ARI")
    px.imshow(heat, text_auto=".2f", aspect="auto", color_continuous_scale="Viridis", title="DM1/DM2 multi-method ARI").write_html(
        FIG / "robustness_method_consistency.html", include_plotlyjs="cdn"
    )
    go.Figure(data=go.Heatmap(z=consensus, colorscale="Viridis")).update_layout(title="Bootstrap consensus matrix").write_html(
        FIG / "bootstrap_consensus.html", include_plotlyjs="cdn"
    )
    fig = px.histogram(perm[perm["metric"] == "null_silhouette"], x="value", nbins=40, title="Permutation null distribution")
    obs = float(perm.loc[perm["metric"] == "observed_silhouette", "value"].iloc[0])
    fig.add_vline(x=obs, line_color="#ef4444")
    fig.write_html(FIG / "permutation_null_distribution.html", include_plotlyjs="cdn")
    ab_sum = pd.concat([
        abl_df.groupby("panel_size", as_index=False)["ARI_vs_full"].mean().assign(metric="feature_ablation"),
        subs_df.groupby("subsample_fraction", as_index=False)["ARI_vs_full"].mean().rename(columns={"subsample_fraction": "panel_size"}).assign(metric="sample_subsample"),
    ], ignore_index=True)
    px.line(ab_sum, x="panel_size", y="ARI_vs_full", color="metric", markers=True, title="Feature ablation and subsample stability").write_html(
        FIG / "feature_ablation_curve.html", include_plotlyjs="cdn"
    )

    summary = {
        "n_bootstrap": N_BOOT,
        "n_permutation": N_PERM,
        "mean_ari_offdiag": float(metric_df.loc[metric_df["method_a"] != metric_df["method_b"], "ARI"].mean()),
        "median_persistence": float(persistence["bootstrap_persistence"].median()),
        "permutation_p": float(perm.loc[perm["metric"] == "observed_silhouette", "empirical_p"].iloc[0]),
        "subsample_50_mean_ari": float(subs_df.loc[sub_df := subs_df["subsample_fraction"] == 0.5, "ARI_vs_full"].mean()),
    }
    write_json(TAB / "robustness_summary.json", summary)
    log_line(LOG, f"v17p2 layer1 done {summary}")


if __name__ == "__main__":
    main()

