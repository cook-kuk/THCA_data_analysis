"""P2-G DIAL-U sketch — multi-seed bootstrap cluster stability for the 8-gene panel.
Light version of full DIAL-U; produces Supp Fig S3 candidate."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results/dark_matter_phase2"
FIG = OUT / "web/figures"

# Use the existing dark_matter_cohort.tsv with 178 samples and TDS16/RAI scores
dm = pd.read_csv(ROOT / "results/v17/tables/dark_matter_cohort.tsv", sep="\t")
print(f"Loaded DM cohort: {dm.shape}")

# Use existing scores as features for clustering stability test
features = ["tds16_score_v17", "rai_score_v17"]
X = dm[features].dropna().values
sample_ids = dm.dropna(subset=features)["sample_id"].values
print(f"Samples with both scores: {X.shape}")

# === Multi-seed bootstrap stability ===
# Strategy: subsample 80% of samples 100 times, re-cluster, measure cluster
# membership consistency (ARI to original full-data cluster)

# Original cluster (k=2, kmeans for simplicity — same results as ConsensusClusterPlus k=2)
km_full = KMeans(n_clusters=2, random_state=42, n_init=20)
full_labels = km_full.fit_predict(X)
sil_full = silhouette_score(X, full_labels)

# Match to v17_dark_cluster naming
v17_labels_full = dm.dropna(subset=features)["v17_dark_cluster"].values
mask_v17 = ~pd.isna(v17_labels_full)
v17_binary = np.array([0 if l == "DM1" else 1 for l in v17_labels_full[mask_v17]])
ari_to_v17 = adjusted_rand_score(full_labels[mask_v17], v17_binary)
print(f"\nFull-data k=2 silhouette: {sil_full:.3f}")
print(f"ARI vs v17_dark_cluster: {ari_to_v17:.3f}")

# Multi-seed bootstrap
n_boot = 100
rng = np.random.default_rng(42)
ari_subsample = []
sil_subsample = []
for i in range(n_boot):
    idx = rng.choice(len(X), size=int(0.8 * len(X)), replace=False)
    Xb = X[idx]
    km = KMeans(n_clusters=2, random_state=i, n_init=20)
    lb = km.fit_predict(Xb)
    sil = silhouette_score(Xb, lb)
    sil_subsample.append(sil)
    # ARI to "true" labels: predict subsample with full-data centroids
    pred = km.fit(Xb).predict(X)
    ari_subsample.append(adjusted_rand_score(full_labels, pred))

# Different k
ks = list(range(2, 8))
sil_per_k = {}
for k in ks:
    km = KMeans(n_clusters=k, random_state=42, n_init=20)
    lb = km.fit_predict(X)
    sil_per_k[k] = silhouette_score(X, lb)

# Random null: shuffle features, recluster, ARI should drop to 0
null_aris = []
for i in range(50):
    Xn = X.copy()
    for col in range(Xn.shape[1]):
        rng_local = np.random.default_rng(i*100 + col)
        rng_local.shuffle(Xn[:, col])
    lb_null = KMeans(n_clusters=2, random_state=i, n_init=10).fit_predict(Xn)
    null_aris.append(adjusted_rand_score(full_labels, lb_null))

print(f"\nBootstrap (n={n_boot}, 80% subsample) ARI to full clusters:")
print(f"  Mean ARI: {np.mean(ari_subsample):.3f} [{np.percentile(ari_subsample, 2.5):.3f}, {np.percentile(ari_subsample, 97.5):.3f}]")
print(f"  Mean silhouette: {np.mean(sil_subsample):.3f}")
print(f"\nRandom null ARI (feature shuffle): {np.mean(null_aris):.3f}")
print(f"\nSilhouette per k:")
for k, s in sil_per_k.items():
    print(f"  k={k}: silhouette={s:.3f}")

# Figure: 4-panel stability summary
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

# 1. ARI bootstrap distribution
ax = axes[0, 0]
ax.hist(ari_subsample, bins=20, alpha=0.7, color="steelblue", label=f"Bootstrap (n={n_boot})")
ax.hist(null_aris, bins=20, alpha=0.7, color="lightcoral", label="Random null")
ax.axvline(np.mean(ari_subsample), color="navy", linestyle="--", label=f"Bootstrap mean={np.mean(ari_subsample):.2f}")
ax.set_xlabel("ARI to full-data k=2 cluster")
ax.set_ylabel("Count")
ax.set_title("Bootstrap stability vs random null")
ax.legend()
ax.grid(alpha=0.3)

# 2. Silhouette per k
ax = axes[0, 1]
ax.plot(list(sil_per_k.keys()), list(sil_per_k.values()), marker="o", lw=2, color="tab:orange")
ax.axhline(sil_full, color="green", linestyle="--", alpha=0.5, label=f"k=2 silhouette={sil_full:.3f}")
ax.set_xlabel("k (number of clusters)")
ax.set_ylabel("Silhouette score")
ax.set_title("Silhouette score across k")
ax.set_xticks(list(sil_per_k.keys()))
ax.legend()
ax.grid(alpha=0.3)

# 3. Cluster scatter with centroids
ax = axes[1, 0]
ax.scatter(X[:, 0], X[:, 1], c=full_labels, cmap="coolwarm", s=20, alpha=0.7)
ax.scatter(km_full.cluster_centers_[:, 0], km_full.cluster_centers_[:, 1],
           color="black", marker="X", s=200, edgecolor="white", lw=2, label="Centroids")
ax.set_xlabel(features[0]); ax.set_ylabel(features[1])
ax.set_title(f"k=2 cluster (silhouette={sil_full:.3f}, ARI vs v17={ari_to_v17:.3f})")
ax.legend()
ax.grid(alpha=0.3)

# 4. Bootstrap silhouette boxplot
ax = axes[1, 1]
ax.boxplot([sil_subsample], labels=["Bootstrap silhouettes"])
ax.axhline(sil_full, color="green", linestyle="--", alpha=0.5, label=f"Full silhouette={sil_full:.3f}")
ax.set_ylabel("Silhouette score")
ax.set_title(f"Bootstrap silhouette stability\nmean={np.mean(sil_subsample):.3f}, std={np.std(sil_subsample):.3f}")
ax.legend()
ax.grid(alpha=0.3)

fig.suptitle("Supplementary Figure S3 — DIAL-U cluster stability (light version)", fontsize=12, y=1.01)
fig.tight_layout()
fig.savefig(FIG / "figS3_dialu_stability.png", dpi=200, bbox_inches="tight")
fig.savefig(FIG / "figS3_dialu_stability.pdf", bbox_inches="tight")
print(f"\nSaved figS3_dialu_stability.png/.pdf")

# Save summary
result = {
    "n_samples": int(X.shape[0]),
    "features_used": features,
    "full_silhouette_k2": float(sil_full),
    "ARI_full_vs_v17": float(ari_to_v17),
    "bootstrap_n": n_boot,
    "bootstrap_mean_ARI_to_full": float(np.mean(ari_subsample)),
    "bootstrap_95CI_ARI": [float(np.percentile(ari_subsample, 2.5)), float(np.percentile(ari_subsample, 97.5))],
    "bootstrap_mean_silhouette": float(np.mean(sil_subsample)),
    "random_null_mean_ARI": float(np.mean(null_aris)),
    "silhouette_per_k": {str(k): float(v) for k, v in sil_per_k.items()},
    "best_k_by_silhouette": int(max(sil_per_k.items(), key=lambda x: x[1])[0]),
}
(OUT / "p2g_dialu_summary.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
