#!/usr/bin/env python3
"""P3 — GSE286332 PTC vs PTC+HT differential analysis (CRITICAL).

Gates Graves' paper 시나리오 2 viability:
  - PyDESeq2 DEG (NG=PTC, TH=PTC+HT)
  - GSEA Hallmark / KEGG / Reactome
  - 8-gene RAI panel score comparison
  - Apply TCGA-trained DM1/DM2 predictor → distribution
  - HLA-I / HLA-II module scores
  - Compare to TCGA Hashimoto-like signature (Q12 finding)
"""
from __future__ import annotations
from pathlib import Path
import gzip, json, warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

PROJ = Path("/opt/thyroid-dash/project")
DATA = Path("/data/thca/v17_korean/GSE286332/GSE286332_all_sample_rawdata.txt.gz")
RES = PROJ / "results/p3_gse286332"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Load count matrix + meta
# ============================================================
print("=== 1. Load GSE286332 ===")
df = pd.read_csv(DATA, sep="\t", low_memory=False)
print(f"  raw shape: {df.shape}")
print(f"  cols (first 30): {df.columns.tolist()[:30]}")

count_cols = [c for c in df.columns if c.endswith("_Read_Count")]
ng_cols = [c for c in count_cols if c.startswith("NG_")]
th_cols = [c for c in count_cols if c.startswith("TH_")]
fpkm_cols = [c for c in df.columns if c.endswith("_FPKM")]
print(f"  NG (PTC w/o HT): {len(ng_cols)}, TH (PTC+HT): {len(th_cols)}")

# Aggregate by Gene_Symbol (sum counts across transcripts)
counts = df.groupby("Gene_Symbol", as_index=True)[count_cols].sum()
counts = counts[counts.sum(axis=1) >= 10]
counts = counts.loc[~counts.index.isna()]
counts = counts.loc[counts.index.astype(str).str.len() > 0]
print(f"  count matrix (genes × samples): {counts.shape}")

# Sample metadata
meta = pd.DataFrame({
    "sample": count_cols,
    "group": ["PTC" if c.startswith("NG_") else "PTC_HT" for c in count_cols],
})
meta = meta.set_index("sample")
print(meta.groupby("group").size())

# ============================================================
# 2. PyDESeq2 DEG
# ============================================================
print("\n=== 2. PyDESeq2 DEG (PTC_HT vs PTC) ===")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

count_t = counts.T.astype(int)  # samples × genes
inference = DefaultInference(n_cpus=4)
dds = DeseqDataSet(counts=count_t, metadata=meta, design="~group", inference=inference, quiet=True)
dds.deseq2()
ds = DeseqStats(dds, contrast=["group", "PTC_HT", "PTC"], inference=inference, quiet=True)
ds.summary()
deg = ds.results_df.copy()
deg["gene"] = deg.index
deg = deg[["gene", "baseMean", "log2FoldChange", "stat", "pvalue", "padj"]]
deg = deg.sort_values("pvalue").reset_index(drop=True)

n_up_05 = int(((deg["padj"] < 0.05) & (deg["log2FoldChange"] > 0)).sum())
n_dn_05 = int(((deg["padj"] < 0.05) & (deg["log2FoldChange"] < 0)).sum())
n_sig = int((deg["padj"] < 0.05).sum())
print(f"  total DEGs (padj<0.05): {n_sig} (up={n_up_05}, dn={n_dn_05})")
deg.to_csv(RES / "deg_ptcht_vs_ptc.tsv", sep="\t", index=False)
print("  ✓ deg_ptcht_vs_ptc.tsv saved")

# Top hits
print("\n  Top 25 up in PTC+HT:")
print(deg[deg["log2FoldChange"] > 0].head(25)[["gene", "log2FoldChange", "padj"]].to_string(index=False))
print("\n  Top 25 down in PTC+HT:")
print(deg[deg["log2FoldChange"] < 0].head(25)[["gene", "log2FoldChange", "padj"]].to_string(index=False))

# ============================================================
# 3. GSEA / pre-ranked enrichment
# ============================================================
print("\n=== 3. GSEA pre-ranked (Hallmark / KEGG / Reactome) ===")
import gseapy as gp

# Rank by signed -log10(p) * sign(LFC), drop NaN
rank = deg.dropna(subset=["pvalue", "log2FoldChange"]).copy()
rank["score"] = -np.log10(rank["pvalue"].clip(lower=1e-300)) * np.sign(rank["log2FoldChange"])
rank_in = rank[["gene", "score"]].drop_duplicates(subset="gene")
rank_in = rank_in.sort_values("score", ascending=False).reset_index(drop=True)

gsea_results = {}
for gs_name in ["MSigDB_Hallmark_2020", "KEGG_2021_Human", "Reactome_2022"]:
    try:
        pre = gp.prerank(
            rnk=rank_in,
            gene_sets=gs_name,
            outdir=None,
            min_size=10, max_size=600,
            permutation_num=1000,
            seed=42,
            verbose=False,
        )
        res = pre.res2d.copy().sort_values("FDR q-val")
        res["geneset_db"] = gs_name
        gsea_results[gs_name] = res
        out = RES / f"gsea_{gs_name}.tsv"
        res.to_csv(out, sep="\t", index=False)
        print(f"  ✓ {gs_name}: {len(res)} terms, top FDR={res.iloc[0]['FDR q-val']:.3g}")
        # Show top up + down
        up = res[res["NES"] > 0].head(10)
        dn = res[res["NES"] < 0].head(10)
        print(f"    UP top 10 (Term | NES | FDR):")
        for _, r in up.iterrows():
            print(f"      {r['Term'][:60]:60s} | {r['NES']:6.2f} | {r['FDR q-val']:.3g}")
        print(f"    DN top 10:")
        for _, r in dn.iterrows():
            print(f"      {r['Term'][:60]:60s} | {r['NES']:6.2f} | {r['FDR q-val']:.3g}")
    except Exception as e:
        print(f"  ✗ {gs_name} failed: {e}")

# ============================================================
# 4. 8-gene RAI panel score
# ============================================================
print("\n=== 4. 8-gene RAI score ===")
GENE_8 = ['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1']

# Use FPKM table for biological score (matches our K2 8-gene scoring approach)
fpkm = df.groupby("Gene_Symbol", as_index=True)[fpkm_cols].sum()
fpkm.columns = [c.replace("_FPKM", "") for c in fpkm.columns]
log_fpkm = np.log2(fpkm + 1.0)

panel_in = [g for g in GENE_8 if g in log_fpkm.index]
print(f"  8-gene available: {panel_in}")

panel_expr = log_fpkm.loc[panel_in].T  # samples × genes
panel_expr["group"] = ["PTC" if s.startswith("NG_") else "PTC_HT" for s in panel_expr.index]

# Z-score each gene across samples then mean → RAI score per sample
panel_z = panel_expr[panel_in].apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9))
panel_expr["RAI_score_8gene"] = panel_z.mean(axis=1)

ng_score = panel_expr.loc[panel_expr["group"] == "PTC", "RAI_score_8gene"].values
th_score = panel_expr.loc[panel_expr["group"] == "PTC_HT", "RAI_score_8gene"].values
mw_u, mw_p = stats.mannwhitneyu(th_score, ng_score, alternative="two-sided")
t_t, t_p = stats.ttest_ind(th_score, ng_score, equal_var=False)
sd_pool = np.sqrt(((len(th_score)-1)*np.var(th_score, ddof=1) + (len(ng_score)-1)*np.var(ng_score, ddof=1)) / (len(th_score) + len(ng_score) - 2))
cohen_d = (th_score.mean() - ng_score.mean()) / max(sd_pool, 1e-9)

print(f"  PTC (n={len(ng_score)}): mean={ng_score.mean():.3f}, sd={ng_score.std(ddof=1):.3f}")
print(f"  PTC+HT (n={len(th_score)}): mean={th_score.mean():.3f}, sd={th_score.std(ddof=1):.3f}")
print(f"  Cohen's d = {cohen_d:.3f}, MW p = {mw_p:.3g}, Welch's t p = {t_p:.3g}")

panel_expr.to_csv(RES / "8gene_panel_per_sample.tsv", sep="\t")

# Per-gene comparison
gene_stats = []
for g in panel_in:
    a = panel_expr.loc[panel_expr["group"] == "PTC", g].values
    b = panel_expr.loc[panel_expr["group"] == "PTC_HT", g].values
    sp = np.sqrt(((len(b)-1)*np.var(b, ddof=1) + (len(a)-1)*np.var(a, ddof=1)) / (len(a)+len(b)-2))
    d = (b.mean() - a.mean()) / max(sp, 1e-9)
    _, p = stats.mannwhitneyu(b, a, alternative="two-sided")
    gene_stats.append(dict(gene=g, mean_PTC=round(a.mean(), 3), mean_PTC_HT=round(b.mean(), 3),
                           cohen_d=round(d, 3), mw_p=p))
gs_df = pd.DataFrame(gene_stats).sort_values("mw_p")
print("\n  per-gene 8-gene comparison (PTC+HT vs PTC):")
print(gs_df.to_string(index=False))
gs_df.to_csv(RES / "8gene_per_gene_compare.tsv", sep="\t", index=False)

# ============================================================
# 5. DM1/DM2 predictor — within-sample-centered profile (TCGA-validated)
# ============================================================
print("\n=== 5. DM1/DM2 prediction via TCGA-trained centered-profile classifier ===")
# Try the existing K2 predictor pickle if available
predictor_path = PROJ / "results/v17_korean/K2_dm_predictor.pkl"
labels_path = PROJ / "results/v17_korean/K2_korean_predictions_v4.tsv"
import pickle
dm_probs = None
try:
    if predictor_path.exists():
        with open(predictor_path, "rb") as f:
            pred_obj = pickle.load(f)
        print(f"  ✓ predictor loaded from {predictor_path}")
        # Apply centered-profile transformation
        sample_means = log_fpkm.mean(axis=0)  # per-sample mean over all genes
        centered = log_fpkm.subtract(sample_means, axis=1)
        feat = centered.loc[panel_in].T  # samples × 8 genes
        if hasattr(pred_obj, "predict_proba"):
            dm_probs = pred_obj.predict_proba(feat)
        elif isinstance(pred_obj, dict) and "model" in pred_obj:
            dm_probs = pred_obj["model"].predict_proba(feat)
except Exception as e:
    print(f"  ! predictor reload skipped: {e}")

if dm_probs is None:
    # Fallback: build a fresh TCGA-trained classifier on centered 8-gene profile
    print("  ! falling back: build fresh TCGA classifier on the fly")
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    tcga_expr = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                             sep="\t", index_col=0)
    cluster = pd.read_csv(PROJ / "results/v17_realfix/R1A_cluster_labels.tsv", sep="\t")
    cluster["DM"] = cluster["cluster"].str[:3]
    tcga_samps = [s for s in cluster["sample_id"] if s in tcga_expr.columns]
    cluster = cluster[cluster["sample_id"].isin(tcga_samps)].set_index("sample_id").loc[tcga_samps]
    tcga_g8 = [g for g in GENE_8 if g in tcga_expr.index]
    print(f"   TCGA 8-gene available: {tcga_g8}")
    tcga_means = tcga_expr.mean(axis=0)
    tcga_centered = tcga_expr.subtract(tcga_means, axis=1)
    Xtr = tcga_centered.loc[tcga_g8, tcga_samps].T.values
    ytr = (cluster["DM"] == "DM1").astype(int).values
    pipe = Pipeline([("sc", StandardScaler()), ("lr", LogisticRegression(C=1.0, max_iter=1000))])
    pipe.fit(Xtr, ytr)
    # Apply to GSE286332 — use SAME 8-gene order as TCGA training
    sample_means = log_fpkm.mean(axis=0)
    g286_centered = log_fpkm.subtract(sample_means, axis=1)
    Xte = g286_centered.loc[tcga_g8].T.values
    dm_probs = pipe.predict_proba(Xte)
    feat_index = log_fpkm.columns

dm_pred_df = pd.DataFrame(dm_probs, index=log_fpkm.columns, columns=["P_DM2", "P_DM1"])
dm_pred_df["group"] = ["PTC" if s.startswith("NG_") else "PTC_HT" for s in dm_pred_df.index]
dm_pred_df["DM_call"] = np.where(dm_pred_df["P_DM1"] >= 0.5, "DM1", "DM2")
print("\n  DM call table:")
print(pd.crosstab(dm_pred_df["group"], dm_pred_df["DM_call"]))
ng_p1 = dm_pred_df.loc[dm_pred_df["group"] == "PTC", "P_DM1"].values
th_p1 = dm_pred_df.loc[dm_pred_df["group"] == "PTC_HT", "P_DM1"].values
_, dm_mw_p = stats.mannwhitneyu(th_p1, ng_p1, alternative="two-sided")
print(f"  P(DM1) PTC mean={ng_p1.mean():.3f} | PTC+HT mean={th_p1.mean():.3f} | MW p={dm_mw_p:.3g}")
dm_pred_df.to_csv(RES / "dm12_predictions.tsv", sep="\t")

# ============================================================
# 6. HLA-I / HLA-II module scores
# ============================================================
print("\n=== 6. HLA-I / HLA-II module scores ===")
HLA_I = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"]
HLA_II = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
          "HLA-DMA", "HLA-DMB", "CIITA", "HLA-DOB"]

def module_score(genes, table):
    keep = [g for g in genes if g in table.index]
    if len(keep) < 3:
        return None, keep
    z = table.loc[keep].apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9), axis=1)
    return z.mean(axis=0), keep

hla1, k1 = module_score(HLA_I, log_fpkm)
hla2, k2 = module_score(HLA_II, log_fpkm)
print(f"  HLA-I genes used: {k1}")
print(f"  HLA-II genes used: {k2}")

hla_df = pd.DataFrame({"HLA_I": hla1, "HLA_II": hla2}, index=log_fpkm.columns)
hla_df["group"] = ["PTC" if s.startswith("NG_") else "PTC_HT" for s in hla_df.index]

for col in ["HLA_I", "HLA_II"]:
    a = hla_df.loc[hla_df["group"] == "PTC", col].values
    b = hla_df.loc[hla_df["group"] == "PTC_HT", col].values
    sp = np.sqrt(((len(b)-1)*np.var(b, ddof=1) + (len(a)-1)*np.var(a, ddof=1)) / (len(a)+len(b)-2))
    d = (b.mean() - a.mean()) / max(sp, 1e-9)
    _, mw_p = stats.mannwhitneyu(b, a, alternative="two-sided")
    print(f"  {col}: PTC mean={a.mean():.3f} | PTC+HT mean={b.mean():.3f} | Cohen d={d:.3f} | MW p={mw_p:.3g}")
hla_df.to_csv(RES / "hla_module_scores.tsv", sep="\t")

# ============================================================
# 7. Compare to TCGA Hashimoto-like signature (Q12 finding, optional)
# ============================================================
print("\n=== 7. TCGA Hashimoto-like signature concordance (optional) ===")
hashi_sig = PROJ / "results/q12_hashimoto/hashimoto_signature_genes.tsv"
if hashi_sig.exists():
    hsig = pd.read_csv(hashi_sig, sep="\t")
    hsig_genes = hsig["gene"].head(200).tolist() if "gene" in hsig.columns else []
    overlap = [g for g in hsig_genes if g in deg["gene"].values]
    if overlap:
        sub = deg[deg["gene"].isin(overlap)]
        n_concord = int((sub["log2FoldChange"] > 0).sum())
        print(f"  Hashimoto-sig overlap: {len(overlap)}/{len(hsig_genes)}, concordant up: {n_concord}")
else:
    print(f"  no Q12 hashimoto signature file at {hashi_sig} — skipping")

# ============================================================
# 8. Summary JSON
# ============================================================
print("\n=== 8. Summary JSON ===")
summary = {
    "n_PTC": len(ng_cols), "n_PTC_HT": len(th_cols),
    "deg_total_padj_05": n_sig, "deg_up": n_up_05, "deg_down": n_dn_05,
    "deg_top10_up": deg[deg["log2FoldChange"] > 0].head(10)[["gene", "log2FoldChange", "padj"]].to_dict("records"),
    "deg_top10_dn": deg[deg["log2FoldChange"] < 0].head(10)[["gene", "log2FoldChange", "padj"]].to_dict("records"),
    "panel_8gene_compare": {
        "mean_PTC": round(float(ng_score.mean()), 3),
        "mean_PTC_HT": round(float(th_score.mean()), 3),
        "cohen_d": round(float(cohen_d), 3),
        "mw_p": float(mw_p),
        "welch_p": float(t_p),
    },
    "DM_distribution": pd.crosstab(dm_pred_df["group"], dm_pred_df["DM_call"]).to_dict(),
    "P_DM1_mean_PTC": round(float(ng_p1.mean()), 3),
    "P_DM1_mean_PTC_HT": round(float(th_p1.mean()), 3),
    "DM1_mw_p": float(dm_mw_p),
    "hla_per_gene_compare": gs_df.to_dict("records"),
    "gsea_top": {
        k: v.head(10)[["Term", "NES", "FDR q-val"]].to_dict("records") for k, v in gsea_results.items()
    },
}

# HLA module summary
for col in ["HLA_I", "HLA_II"]:
    a = hla_df.loc[hla_df["group"] == "PTC", col].values
    b = hla_df.loc[hla_df["group"] == "PTC_HT", col].values
    summary[f"{col}_mean_PTC"] = round(float(a.mean()), 3)
    summary[f"{col}_mean_PTC_HT"] = round(float(b.mean()), 3)
    sp = np.sqrt(((len(b)-1)*np.var(b, ddof=1) + (len(a)-1)*np.var(a, ddof=1)) / (len(a)+len(b)-2))
    summary[f"{col}_cohen_d"] = round(float((b.mean() - a.mean()) / max(sp, 1e-9)), 3)
    _, p = stats.mannwhitneyu(b, a, alternative="two-sided")
    summary[f"{col}_mw_p"] = float(p)

(RES / "P3_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n✓ All outputs saved to {RES}")
print("\n=== DECISION GATE ===")
sig_strong = (n_sig > 50) and (gsea_results.get("MSigDB_Hallmark_2020") is not None and
                                (gsea_results["MSigDB_Hallmark_2020"]["FDR q-val"] < 0.1).any())
print(f"  N DEGs (padj<0.05): {n_sig} | strong: {sig_strong}")
print(f"  PTC vs PTC+HT 8-gene Cohen d: {cohen_d:.3f}")
print(f"  HLA-II Cohen d: {summary['HLA_II_cohen_d']:.3f}")
print(f"  → 시나리오 2 viability: {'✅ STRONG' if sig_strong and abs(summary['HLA_II_cohen_d']) > 0.5 else '🟡 MODERATE' if n_sig > 0 else '❌ FAIL'}")
