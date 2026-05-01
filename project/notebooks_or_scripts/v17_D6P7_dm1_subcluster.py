#!/usr/bin/env python3
"""D6-P7 — DM1 sub-A (n=72) vs sub-B (n=19) mechanism layer.

Re-derive sub-cluster membership from DM1-only KMeans k=2 on TIERA67/8-gene scores,
then characterize: DEG, GSEA, clinical, BCR/TLS, Hashimoto-like join, age cross.
"""
from __future__ import annotations
from pathlib import Path
import json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings("ignore")

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/d6p7_dm1_subcluster"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Load
# ============================================================
print("=== 1. Load TCGA + DM cluster + clinical ===")
expr = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                   sep="\t", index_col=0)
labels = pd.read_csv(PROJ / "results/v17_realfix/R1A_cluster_labels.tsv", sep="\t")
labels["DM"] = labels["cluster"].str[:3]
clin = pd.read_csv(PROJ / "results/tables/tcga_thca_clinical_extended.tsv", sep="\t")
muts = pd.read_csv(PROJ / "results/tables/tcga_thca_mutation_groups.tsv", sep="\t")

samps = [s for s in labels["sample_id"] if s in expr.columns]
labels = labels[labels["sample_id"].isin(samps)].set_index("sample_id").loc[samps]
clin = clin.set_index("sample_id")
muts = muts.set_index("sample_id")
print(f"  total: {len(samps)}, DM1={int((labels['DM']=='DM1').sum())}, DM2={int((labels['DM']=='DM2').sum())}")

# ============================================================
# 2. Re-derive DM1 sub-A vs sub-B (KMeans k=2 on TIERA67)
# ============================================================
print("\n=== 2. DM1 sub-cluster (KMeans k=2 on TIERA67) ===")
dm1_samps = labels[labels["DM"] == "DM1"].index.tolist()
print(f"  DM1 n={len(dm1_samps)}")

tiera = [l.strip() for l in (PROJ / "metadata/tierA67_genes.txt").read_text().splitlines()
         if l.strip() and not l.strip().startswith("#") and not l.strip().startswith("[")]
tiera_in = [g for g in tiera if g in expr.index]

X = expr.loc[tiera_in, dm1_samps].T.values
Xs = StandardScaler().fit_transform(X)
km = KMeans(n_clusters=2, random_state=42, n_init=20)
sub_labels = km.fit_predict(Xs)
counts = pd.Series(sub_labels).value_counts().sort_index()
print(f"  sub-cluster sizes (raw): {counts.to_dict()}")

# Larger cluster = sub-A, smaller = sub-B
larger = counts.idxmax()
sub_label_named = ["sub_A" if l == larger else "sub_B" for l in sub_labels]
sub = pd.Series(sub_label_named, index=dm1_samps, name="sub_cluster")
print(f"  sub-A n={(sub=='sub_A').sum()}, sub-B n={(sub=='sub_B').sum()}")

# Save
sub.to_csv(RES / "dm1_subcluster_labels.tsv", sep="\t")

# ============================================================
# 3. Differential expression (sub-B vs sub-A)
# ============================================================
print("\n=== 3. DEG sub-B vs sub-A (Welch t per gene + BH-FDR) ===")
A_samps = sub[sub == "sub_A"].index.tolist()
B_samps = sub[sub == "sub_B"].index.tolist()
mu_A = expr[A_samps].mean(axis=1)
mu_B = expr[B_samps].mean(axis=1)
v_A = expr[A_samps].var(axis=1, ddof=1)
v_B = expr[B_samps].var(axis=1, ddof=1)
nA, nB = len(A_samps), len(B_samps)
sp = np.sqrt(((nA-1)*v_A + (nB-1)*v_B) / (nA+nB-2)).clip(lower=1e-9)
d = (mu_B - mu_A) / sp
log2FC = mu_B - mu_A

# Welch t-test
from scipy.stats import ttest_ind
t_results = []
expr_A = expr[A_samps].values
expr_B = expr[B_samps].values
for i, gene in enumerate(expr.index):
    a = expr_A[i]; b = expr_B[i]
    if np.var(a) < 1e-12 and np.var(b) < 1e-12:
        p = 1.0
    else:
        try:
            _, p = ttest_ind(b, a, equal_var=False)
        except Exception:
            p = 1.0
    t_results.append(p)
deg = pd.DataFrame({
    "gene": expr.index,
    "log2FC_B_vs_A": log2FC.values,
    "cohen_d_B_vs_A": d.values,
    "p_value": t_results,
})
# BH-FDR
from statsmodels.stats.multitest import multipletests
deg["padj"] = multipletests(deg["p_value"], method="fdr_bh")[1]
deg = deg.sort_values("p_value").reset_index(drop=True)
n_sig = int((deg["padj"] < 0.05).sum())
n_up = int(((deg["padj"] < 0.05) & (deg["log2FC_B_vs_A"] > 0)).sum())
n_dn = int(((deg["padj"] < 0.05) & (deg["log2FC_B_vs_A"] < 0)).sum())
print(f"  DEGs (padj<0.05): {n_sig} (up={n_up}, dn={n_dn})")
deg.to_csv(RES / "dm1_subBvA_deg.tsv", sep="\t", index=False)

print("\n  Top 15 up in sub-B (low-RAI sub-cluster):")
print(deg[deg["log2FC_B_vs_A"] > 0].head(15)[["gene", "log2FC_B_vs_A", "cohen_d_B_vs_A", "padj"]].to_string(index=False))
print("\n  Top 15 down in sub-B:")
print(deg[deg["log2FC_B_vs_A"] < 0].head(15)[["gene", "log2FC_B_vs_A", "cohen_d_B_vs_A", "padj"]].to_string(index=False))

# ============================================================
# 4. Score profile
# ============================================================
print("\n=== 4. Score profile by sub-cluster ===")
GENE_8 = ['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1']
TDS_16 = ['DIO1', 'DIO2', 'DUOX1', 'DUOX2', 'FOXE1', 'GLIS3', 'NKX2-1', 'PAX8',
          'SLC26A4', 'SLC5A5', 'SLC5A8', 'TG', 'THRA', 'THRB', 'TPO', 'TSHR']
HLA_I = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"]
HLA_II = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
          "HLA-DMA", "HLA-DMB", "CIITA", "HLA-DOB"]
TLS = ["CCL19", "CCL21", "CXCL13", "CCR7", "CXCR5", "SELL", "LAMP3",
        "MS4A1", "CD79A", "CD79B", "PTGDS", "TRBC2"]
B_CELL = ["MS4A1", "CD19", "CD79A", "CD79B", "BANK1", "BLK"]
DEDIFF = ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "MMP9", "LOX"]

def zmean(genes, mat):
    keep = [g for g in genes if g in mat.index]
    if not keep:
        return None
    z = mat.loc[keep].apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9), axis=1)
    return z.mean(axis=0)

E_dm1 = expr[dm1_samps]
score_df = pd.DataFrame({
    "g8_RAI": zmean(GENE_8, E_dm1),
    "TDS16": zmean(TDS_16, E_dm1),
    "HLA_I": zmean(HLA_I, E_dm1),
    "HLA_II": zmean(HLA_II, E_dm1),
    "TLS": zmean(TLS, E_dm1),
    "B_cell": zmean(B_CELL, E_dm1),
    "Dediff_EMT": zmean(DEDIFF, E_dm1),
}, index=dm1_samps)
score_df["sub"] = sub.values

# Per-cluster median + comparison
cmp_rows = []
for col in ["g8_RAI", "TDS16", "HLA_I", "HLA_II", "TLS", "B_cell", "Dediff_EMT"]:
    a = score_df.loc[score_df["sub"] == "sub_A", col].values
    b = score_df.loc[score_df["sub"] == "sub_B", col].values
    sp = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) / (len(a)+len(b)-2))
    dval = (b.mean() - a.mean()) / max(sp, 1e-9)
    _, p = stats.mannwhitneyu(b, a, alternative="two-sided")
    cmp_rows.append(dict(score=col, median_A=round(np.median(a), 3), median_B=round(np.median(b), 3),
                          cohen_d=round(dval, 3), mw_p=p))
cmp_df = pd.DataFrame(cmp_rows)
print(cmp_df.to_string(index=False))
cmp_df.to_csv(RES / "subcluster_score_profile.tsv", sep="\t", index=False)
score_df.to_csv(RES / "subcluster_scores.tsv", sep="\t")

# ============================================================
# 5. Clinical phenotype × sub-cluster
# ============================================================
print("\n=== 5. Clinical phenotype × sub-cluster ===")
sd = score_df.join(clin[["age_at_diagnosis", "gender", "stage", "tumor_size_mm",
                          "vital_status", "os_days", "os_event"]], how="left")
clin_rows = []
for col in ["age_at_diagnosis", "tumor_size_mm", "os_days"]:
    a = sd.loc[sd["sub"] == "sub_A", col].dropna().values
    b = sd.loc[sd["sub"] == "sub_B", col].dropna().values
    if len(a) > 3 and len(b) > 3:
        sp = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) / (len(a)+len(b)-2))
        dval = (b.mean() - a.mean()) / max(sp, 1e-9)
        _, p = stats.mannwhitneyu(b, a, alternative="two-sided")
        clin_rows.append(dict(metric=col, mean_A=round(a.mean(), 2), mean_B=round(b.mean(), 2),
                                cohen_d=round(dval, 3), mw_p=p))
clin_df = pd.DataFrame(clin_rows)
print(clin_df.to_string(index=False))
clin_df.to_csv(RES / "subcluster_clinical.tsv", sep="\t", index=False)

# Young-onset proportion
sd["young"] = (sd["age_at_diagnosis"] < 45).astype(int)
sd["female"] = (sd["gender"] == "female").astype(int)
sd["adv_stage"] = sd["stage"].astype(str).isin(["Stage III", "Stage IVA", "Stage IVB", "Stage IVC"]).astype(int)

print("\n  Young-onset (<45 yr):")
print(pd.crosstab(sd["sub"], sd["young"]))
print("\n  Female:")
print(pd.crosstab(sd["sub"], sd["female"]))
print("\n  Advanced stage (III-IV):")
print(pd.crosstab(sd["sub"], sd["adv_stage"]))

# ============================================================
# 6. Mutation cross
# ============================================================
print("\n=== 6. BRAF / RAS / TERT × sub-cluster ===")
sd = sd.join(muts[["has_braf_v600e", "has_ras_mut"]], how="left")
print(pd.crosstab([sd["sub"]], sd["has_braf_v600e"]))
print(pd.crosstab([sd["sub"]], sd["has_ras_mut"]))

# ============================================================
# 7. Hashimoto-like join
# ============================================================
print("\n=== 7. Hashimoto-like × sub-cluster (P2 join) ===")
hashi_path = PROJ / "results/d4p2_tcga_hashimoto_signature/tcga_signature_scores.tsv"
if hashi_path.exists():
    hashi = pd.read_csv(hashi_path, sep="\t", index_col=0)
    join = sd.join(hashi[["sig_score", "hashi_otsu", "hashi_top20"]], how="left")
    for label in ["hashi_otsu", "hashi_top20"]:
        if label in join.columns:
            ct = pd.crosstab(join["sub"], join[label])
            print(f"\n  {label}:")
            print(ct)
            n11 = ct.loc["sub_A", 1] if ("sub_A" in ct.index and 1 in ct.columns) else 0
            n10 = ct.loc["sub_A", 0] if ("sub_A" in ct.index and 0 in ct.columns) else 0
            n21 = ct.loc["sub_B", 1] if ("sub_B" in ct.index and 1 in ct.columns) else 0
            n20 = ct.loc["sub_B", 0] if ("sub_B" in ct.index and 0 in ct.columns) else 0
            try:
                or_, p = stats.fisher_exact([[n11, n10], [n21, n20]])
                print(f"    Fisher OR={or_:.3f}, p={p:.3g}  (sub-A vs sub-B)")
            except Exception:
                pass
else:
    print(f"  Hashimoto-like signature scores not yet available at {hashi_path}")
    join = sd

# ============================================================
# 8. Save summary
# ============================================================
summary = {
    "n_DM1": len(dm1_samps), "n_subA": int((sub == "sub_A").sum()), "n_subB": int((sub == "sub_B").sum()),
    "deg_padj_05": n_sig, "deg_up_in_B": n_up, "deg_dn_in_B": n_dn,
    "score_profile": cmp_rows,
    "clinical": clin_rows,
}
(RES / "D6P7_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n✓ Outputs to {RES}")
