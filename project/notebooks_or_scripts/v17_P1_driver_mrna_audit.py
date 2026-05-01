#!/usr/bin/env python3
"""P1 — BRAF/RAS/TERT mRNA × mutation status (TCGA-THCA) + driver mRNA AUC
+ TIERA67 67-gene Cohen's d full ranking (sanity check for Q3/Q4 reviewer answer)."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/p1_driver_mrna_audit"
RES.mkdir(parents=True, exist_ok=True)

# Load TCGA expression (genes × samples)
print("=== Load TCGA-THCA expression ===")
expr = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                   sep="\t", index_col=0)
print(f"  expr: {expr.shape}")

# Load DM cluster labels
labels = pd.read_csv(PROJ / "results/v17_realfix/R1A_cluster_labels.tsv", sep="\t")
labels["DM"] = labels["cluster"].str[:3]  # DM1 or DM2
print(f"  R1A labels: {len(labels)}")

# Load mutation groups
muts = pd.read_csv(PROJ / "results/tables/tcga_thca_mutation_groups.tsv", sep="\t")
print(f"  mutations: {len(muts)}, columns: {muts.columns.tolist()}")

# Load TIERA67
tiera = [l.strip() for l in (PROJ / "metadata/tierA67_genes.txt").read_text().splitlines()
         if l.strip() and not l.strip().startswith("#") and not l.strip().startswith("[")]
print(f"  TIERA67: {len(tiera)} genes")

# Merge labels + muts (sample_id intersection)
m = labels.merge(muts, on="sample_id", how="inner")
print(f"  labels ∩ muts: {len(m)}")

# Get expression for the matched samples
samples = [s for s in m["sample_id"] if s in expr.columns]
expr_m = expr[samples]
m = m[m["sample_id"].isin(samples)].set_index("sample_id").loc[samples]
print(f"  matched samples: {len(samples)}")

# ============= Task 1-3: driver mRNA × mutation status =============
def quantify_driver(gene_name, mut_col, label):
    """driver mRNA expression by mutation status."""
    if gene_name not in expr_m.index:
        return None
    x = expr_m.loc[gene_name].values
    has_mut = m[mut_col].astype(bool).values
    n_mut = int(has_mut.sum()); n_wt = int((~has_mut).sum())
    if n_mut < 5 or n_wt < 5:
        return None
    mu_mut = x[has_mut].mean(); mu_wt = x[~has_mut].mean()
    sd_pool = np.sqrt(((n_mut-1)*x[has_mut].var(ddof=1) + (n_wt-1)*x[~has_mut].var(ddof=1)) / (n_mut+n_wt-2))
    cohen_d = (mu_mut - mu_wt) / max(sd_pool, 1e-9)
    fc = mu_mut - mu_wt  # log2 fold change
    u, p = stats.mannwhitneyu(x[has_mut], x[~has_mut], alternative="two-sided")
    return dict(gene=gene_name, label=label, n_mut=n_mut, n_wt=n_wt,
                mean_mut=round(mu_mut, 3), mean_wt=round(mu_wt, 3),
                log2_fc=round(fc, 3), cohen_d=round(cohen_d, 3),
                mw_p=p)

print("\n=== Task 1-3: Driver mRNA × mutation status ===")
results = []
# Task 1: BRAF V600E
results.append(quantify_driver("BRAF", "has_braf_v600e", "BRAF V600E carrier"))
# Task 2: RAS hotspot (any HRAS/NRAS/KRAS)
results.append(quantify_driver("HRAS", "has_ras_mut", "RAS hotspot any"))
results.append(quantify_driver("NRAS", "has_ras_mut", "RAS hotspot any"))
results.append(quantify_driver("KRAS", "has_ras_mut", "RAS hotspot any"))
# TERT — check column existence
for col in ["has_tert_promoter", "has_tert", "tert_status"]:
    if col in m.columns:
        results.append(quantify_driver("TERT", col, f"TERT promoter ({col})"))
        break
results = [r for r in results if r is not None]
df = pd.DataFrame(results)
print(df.to_string(index=False))
df.to_csv(RES / "driver_mrna_mutation_audit.tsv", sep="\t", index=False)

# ============= Task 4: driver mRNA AUC for DM1 vs DM2 =============
print("\n=== Task 4: Driver mRNA single-feature AUC for DM1 vs DM2 ===")
y = (m["DM"] == "DM1").astype(int).values  # DM1 = 1
auc_results = []
for gene in ["BRAF", "HRAS", "NRAS", "KRAS", "TERT"]:
    if gene in expr_m.index:
        x = expr_m.loc[gene].values
        # Try both directions
        auc1 = roc_auc_score(y, x)
        auc2 = roc_auc_score(y, -x)
        auc = max(auc1, auc2)
        cohen_d = abs((x[y==1].mean() - x[y==0].mean()) / max(np.std(x), 1e-9))
        auc_results.append(dict(gene=gene, n_dm1=int(y.sum()), n_dm2=int((1-y).sum()),
                                 auc=round(auc, 3), cohen_d=round(cohen_d, 3)))
auc_df = pd.DataFrame(auc_results)
print(auc_df.to_string(index=False))
auc_df.to_csv(RES / "driver_mrna_dm_auc.tsv", sep="\t", index=False)

# ============= Task 5: TIERA67 full Cohen's d ranking =============
print("\n=== Task 5: TIERA67 67-gene Cohen's d full ranking (DM1 vs DM2) ===")
tiera_in_expr = [g for g in tiera if g in expr_m.index]
print(f"  TIERA67 in expr: {len(tiera_in_expr)}/67")
ranking = []
for gene in tiera_in_expr:
    x = expr_m.loc[gene].values
    n1, n0 = int(y.sum()), int((1-y).sum())
    if n1 < 5 or n0 < 5:
        continue
    mu1 = x[y==1].mean(); mu0 = x[y==0].mean()
    sd_pool = np.sqrt(((n1-1)*x[y==1].var(ddof=1) + (n0-1)*x[y==0].var(ddof=1)) / (n1+n0-2))
    d = (mu1 - mu0) / max(sd_pool, 1e-9)
    auc = roc_auc_score(y, x); auc = max(auc, 1-auc)
    ranking.append(dict(gene=gene, cohen_d=round(d, 3), abs_d=round(abs(d), 3), auc=round(auc, 3)))
rank_df = pd.DataFrame(ranking).sort_values("abs_d", ascending=False).reset_index(drop=True)
rank_df["rank"] = range(1, len(rank_df)+1)
GENE_8 = ['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1']
rank_df["in_8gene"] = rank_df["gene"].isin(GENE_8)
DRIVERS = {"BRAF","HRAS","NRAS","KRAS","TERT","RET","NTRK1","NTRK3","ALK","PPARG","EIF1AX","TP53","CDKN2A","CDKN2B","PIK3CA","AKT1","PTEN","ATM"}
rank_df["is_driver"] = rank_df["gene"].isin(DRIVERS)

print("\n=== Top 20 by |Cohen's d| ===")
print(rank_df.head(20).to_string(index=False))
print(f"\n=== 8-gene panel ranks ===")
print(rank_df[rank_df["in_8gene"]].to_string(index=False))
print(f"\n=== Driver gene ranks (BRAF/RAS/TERT etc) ===")
print(rank_df[rank_df["is_driver"]].to_string(index=False))

rank_df.to_csv(RES / "top20_by_d_full_tiera67.tsv", sep="\t", index=False)

# ============= Summary JSON =============
summary = {
    "n_samples": len(samples),
    "DM1_n": int(y.sum()), "DM2_n": int((1-y).sum()),
    "driver_mrna_mutation_audit": df.to_dict("records"),
    "driver_mrna_dm_auc": auc_df.to_dict("records"),
    "tiera67_top20": rank_df.head(20).to_dict("records"),
    "8gene_ranks": rank_df[rank_df["in_8gene"]][["gene","rank","cohen_d","auc"]].to_dict("records"),
    "driver_ranks": rank_df[rank_df["is_driver"]][["gene","rank","cohen_d","auc"]].to_dict("records"),
}
(RES / "P1_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n✓ All outputs saved to {RES}")
