#!/usr/bin/env python3
"""P4 — Pan-genome top-5000 MAD vs TIERA67 candidate pool.

Reviewer Q3: "Does cluster definition depend on candidate pool restriction?"
- Apply UNRESTRICTED top-5000 MAD selection on TCGA-THCA
- Cluster with same algorithm (consensus k-means k=2)
- Measure ARI vs original DM1/DM2 labels
- Show overlap of top-discriminator gene set with TIERA67 / 8-gene
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/p4_pangenome_vs_tiera67"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Load
# ============================================================
print("=== 1. Load TCGA-THCA expression + DM labels ===")
expr = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                   sep="\t", index_col=0)
print(f"  expr: {expr.shape}")

labels = pd.read_csv(PROJ / "results/v17_realfix/R1A_cluster_labels.tsv", sep="\t")
labels["DM"] = labels["cluster"].str[:3]
samps = [s for s in labels["sample_id"] if s in expr.columns]
labels = labels[labels["sample_id"].isin(samps)].set_index("sample_id").loc[samps]
print(f"  matched samples: {len(samps)}, DM1={int((labels['DM']=='DM1').sum())} DM2={int((labels['DM']=='DM2').sum())}")

E = expr[samps]  # genes × samples
y_orig = (labels["DM"] == "DM1").astype(int).values

# Load TIERA67
tiera = [l.strip() for l in (PROJ / "metadata/tierA67_genes.txt").read_text().splitlines()
         if l.strip() and not l.strip().startswith("#") and not l.strip().startswith("[")]
GENE_8 = ['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1']

# ============================================================
# 2. Pan-genome top-5000 MAD selection
# ============================================================
print("\n=== 2. Pan-genome top-5000 MAD ===")
mad = (E.sub(E.median(axis=1), axis=0)).abs().median(axis=1)
top5000 = mad.sort_values(ascending=False).head(5000).index.tolist()
top1000 = mad.sort_values(ascending=False).head(1000).index.tolist()
top200 = mad.sort_values(ascending=False).head(200).index.tolist()
print(f"  top 5000 MAD selected from {len(mad)} genes")

# Overlap with TIERA67 / 8-gene
ovl_5000_tiera = sorted(set(top5000) & set(tiera))
ovl_5000_g8 = sorted(set(top5000) & set(GENE_8))
ovl_1000_tiera = sorted(set(top1000) & set(tiera))
ovl_1000_g8 = sorted(set(top1000) & set(GENE_8))
ovl_200_tiera = sorted(set(top200) & set(tiera))
ovl_200_g8 = sorted(set(top200) & set(GENE_8))
print(f"  TIERA67 captured by top-5000: {len(ovl_5000_tiera)}/{len(tiera)}")
print(f"  8-gene captured by top-5000: {len(ovl_5000_g8)}/8 — {ovl_5000_g8}")
print(f"  TIERA67 captured by top-1000: {len(ovl_1000_tiera)}/{len(tiera)}")
print(f"  8-gene captured by top-1000: {len(ovl_1000_g8)}/8 — {ovl_1000_g8}")
print(f"  TIERA67 captured by top-200: {len(ovl_200_tiera)}/{len(tiera)}")
print(f"  8-gene captured by top-200: {len(ovl_200_g8)}/8 — {ovl_200_g8}")

# ============================================================
# 3. Cluster on each candidate set, compare to original DM labels
# ============================================================
print("\n=== 3. Cluster k=2 on each candidate set, compare to DM ===")

def cluster_and_score(genes_in, label):
    g = [x for x in genes_in if x in E.index]
    X = E.loc[g].T.values  # samples × genes
    Xs = StandardScaler().fit_transform(X)
    km = KMeans(n_clusters=2, random_state=42, n_init=20)
    y_pred = km.fit_predict(Xs)
    ari = adjusted_rand_score(y_orig, y_pred)
    nmi = normalized_mutual_info_score(y_orig, y_pred)
    # Account for label flip
    ari = max(ari, adjusted_rand_score(y_orig, 1 - y_pred))
    return dict(label=label, n_genes=len(g), ARI=round(ari, 3), NMI=round(nmi, 3))

rows = []
rows.append(cluster_and_score(GENE_8, "8-gene panel"))
rows.append(cluster_and_score(tiera, "TIERA67 (full)"))
rows.append(cluster_and_score(top200, "Pan-genome top 200 MAD"))
rows.append(cluster_and_score(top1000, "Pan-genome top 1000 MAD"))
rows.append(cluster_and_score(top5000, "Pan-genome top 5000 MAD"))
# Also: TIERA67 minus 8-gene (to see if other categories matter)
rows.append(cluster_and_score([g for g in tiera if g not in GENE_8], "TIERA67 minus 8-gene"))
# 16-gene TDS_core
TDS_16 = ['DIO1', 'DIO2', 'DUOX1', 'DUOX2', 'FOXE1', 'GLIS3', 'NKX2-1', 'PAX8',
          'SLC26A4', 'SLC5A5', 'SLC5A8', 'TG', 'THRA', 'THRB', 'TPO', 'TSHR']
rows.append(cluster_and_score(TDS_16, "TDS_core 16"))
# Driver-only
DRIVER = ['BRAF', 'NRAS', 'HRAS', 'KRAS', 'RET', 'NTRK1', 'NTRK3', 'ALK', 'PAX8', 'PPARG', 'TERT', 'EIF1AX']
rows.append(cluster_and_score(DRIVER, "Driver_anchor 12"))
ari_df = pd.DataFrame(rows)
print(ari_df.to_string(index=False))
ari_df.to_csv(RES / "ari_comparison.tsv", sep="\t", index=False)

# ============================================================
# 4. Univariate Cohen d ranking — pan-genome vs TIERA67
# ============================================================
print("\n=== 4. Univariate Cohen's d ranking (DM1 vs DM2) on pan-genome ===")
# Compute Cohen's d for ALL genes
mu1 = E[labels[labels["DM"]=="DM1"].index].mean(axis=1)
mu0 = E[labels[labels["DM"]=="DM2"].index].mean(axis=1)
n1 = (labels["DM"] == "DM1").sum(); n0 = (labels["DM"] == "DM2").sum()
v1 = E[labels[labels["DM"]=="DM1"].index].var(axis=1, ddof=1)
v0 = E[labels[labels["DM"]=="DM2"].index].var(axis=1, ddof=1)
sp = np.sqrt(((n1-1)*v1 + (n0-1)*v0) / (n1+n0-2))
d = (mu1 - mu0) / sp.clip(lower=1e-9)
d_rank = pd.DataFrame({"gene": d.index, "cohen_d": d.values, "abs_d": d.abs().values})
d_rank = d_rank.sort_values("abs_d", ascending=False).reset_index(drop=True)
d_rank["rank"] = range(1, len(d_rank)+1)
d_rank["in_TIERA67"] = d_rank["gene"].isin(set(tiera))
d_rank["in_8gene"] = d_rank["gene"].isin(set(GENE_8))
d_rank["in_top5000_MAD"] = d_rank["gene"].isin(set(top5000))
d_rank.to_csv(RES / "pangenome_cohen_d_ranking.tsv", sep="\t", index=False)
print(f"  saved {len(d_rank)} gene ranking")
print("\n  Top 30 by |Cohen d| (pan-genome):")
print(d_rank.head(30).to_string(index=False))

# Where does the 8-gene rank?
print("\n  8-gene panel ranks within pan-genome |Cohen d|:")
print(d_rank[d_rank["in_8gene"]].to_string(index=False))

print("\n  TIERA67 in top 100 of pan-genome ranking:")
top100_tiera = d_rank.head(100)[d_rank.head(100)["in_TIERA67"]]
print(f"  {len(top100_tiera)}/100 top genes are in TIERA67")
print(f"  median rank of TIERA67 genes: {d_rank[d_rank['in_TIERA67']]['rank'].median():.0f}")
print(f"  median rank of 8-gene: {d_rank[d_rank['in_8gene']]['rank'].median():.0f}")
print(f"  median rank of all genes: {d_rank['rank'].median():.0f}")

# Hypergeometric test — is TIERA67 enriched in top 100 vs background?
N_total = len(d_rank); K_tiera = int(d_rank["in_TIERA67"].sum())
n_top = 100; k_obs = len(top100_tiera)
hg_p = stats.hypergeom.sf(k_obs - 1, N_total, K_tiera, n_top)
print(f"  Hypergeometric p (TIERA67 enrichment in top 100): {hg_p:.3g}")

# ============================================================
# 5. Top-N discriminator coverage ladder
# ============================================================
print("\n=== 5. Top-N discriminator coverage ladder ===")
ladder_rows = []
for N in [20, 50, 100, 200, 500, 1000, 2000, 5000]:
    topN = d_rank.head(N)
    overlap_t = topN["in_TIERA67"].sum()
    overlap_8 = topN["in_8gene"].sum()
    pct_t = round(100 * overlap_t / len(tiera), 1)
    pct_8 = round(100 * overlap_8 / 8, 1)
    ladder_rows.append(dict(top_N=N, TIERA67_in_topN=int(overlap_t), pct_of_TIERA=pct_t,
                             g8_in_topN=int(overlap_8), pct_of_g8=pct_8))
ladder_df = pd.DataFrame(ladder_rows)
print(ladder_df.to_string(index=False))
ladder_df.to_csv(RES / "topN_coverage_ladder.tsv", sep="\t", index=False)

# ============================================================
# 6. Summary
# ============================================================
summary = {
    "TIERA67_total": len(tiera),
    "g8_total": 8,
    "ari_table": rows,
    "g8_ranks_in_pangenome": d_rank[d_rank["in_8gene"]][["gene","rank","cohen_d"]].to_dict("records"),
    "top100_pangenome_TIERA_count": int(k_obs),
    "hypergeom_p_TIERA_top100": float(hg_p),
    "median_rank_TIERA67": int(d_rank[d_rank["in_TIERA67"]]["rank"].median()),
    "median_rank_g8": int(d_rank[d_rank["in_8gene"]]["rank"].median()),
    "median_rank_all": int(d_rank["rank"].median()),
    "topN_coverage_ladder": ladder_rows,
    "interpretation": ("Robustness of 8-gene cluster definition: ARI between cluster from 8-gene "
                       "vs unrestricted top-5000 MAD pan-genome cluster — high ARI implies the "
                       "8-gene panel captures the same axis as unbiased pan-genome clustering."),
}
(RES / "P4_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n✓ All outputs saved to {RES}")
