#!/usr/bin/env python3
"""D8-B — Korean GSE213647 (Lee 2024 n=632) replication of Hashimoto-like signature.

Fix D4-P2 zmean None error + run cross-cohort generalization on Korean PTC.
Confirms whether DM2-enriched Hashimoto-like axis reproduces in second
independent Korean cohort.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.mixture import GaussianMixture

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/d8b_korean_replication"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Load Korean GSE213647 (Lee 2024)
# ============================================================
print("=== 1. Load GSE213647 (Lee 2024) ===")
expr_kr = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_log2.tsv",
                       sep="\t", index_col=0)
print(f"  shape: {expr_kr.shape}")
# Map gene symbols → Ensembl
ensg_map = pd.read_csv(PROJ / "metadata/ensembl_to_symbol.tsv", sep="\t")
sym_to_ensg = dict(zip(ensg_map["gene_symbol"], ensg_map["ensembl_id"]))
def to_ensg(syms):
    out = []
    for s in syms:
        e = sym_to_ensg.get(s)
        if e and e in expr_kr.index:
            out.append(e)
    return out

# Load GSE286332 PTC+HT signature (D4-P2 method)
deg = pd.read_csv(PROJ / "results/p3_gse286332/deg_ptcht_vs_ptc.tsv", sep="\t")
sig_up = deg[(deg["padj"] < 0.01) & (deg["log2FoldChange"] > 1.5)].sort_values("padj").head(150)
sig_dn = deg[(deg["padj"] < 0.01) & (deg["log2FoldChange"] < -1.0)].sort_values("padj").head(50)

up_kr = to_ensg(sig_up["gene"].tolist())
dn_kr = to_ensg(sig_dn["gene"].tolist())
print(f"  Signature in Korean: up={len(up_kr)}/{len(sig_up)}, dn={len(dn_kr)}/{len(sig_dn)}")

if len(up_kr) < 10:
    print(f"  ! up gene overlap too small ({len(up_kr)}); using less strict signature")
    sig_up = deg[(deg["padj"] < 0.05) & (deg["log2FoldChange"] > 1.0)].sort_values("padj").head(200)
    up_kr = to_ensg(sig_up["gene"].tolist())
    print(f"  → expanded: up={len(up_kr)}")

# ============================================================
# 2. Per-sample signature score (handle empty dn case explicitly)
# ============================================================
print("\n=== 2. Korean per-sample signature score ===")
def zmean(genes, mat):
    keep = [g for g in genes if g in mat.index]
    if len(keep) < 3:
        return pd.Series(0.0, index=mat.columns)  # ★ FIX: don't return None
    z = mat.loc[keep].apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9), axis=1)
    return z.mean(axis=0)

sig_up_score = zmean(up_kr, expr_kr)
sig_dn_score = zmean(dn_kr, expr_kr)
sig_score = sig_up_score - sig_dn_score

# Modules — convert gene symbols to Ensembl
HLA_II = to_ensg(["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
                   "HLA-DMA", "HLA-DMB", "CIITA", "HLA-DOB"])
HLA_I = to_ensg(["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"])
B_CELL = to_ensg(["MS4A1", "CD19", "CD79A", "CD79B", "BANK1", "BLK"])
TLS = to_ensg(["CCL19", "CCL21", "CXCL13", "CCR7", "CXCR5", "SELL", "LAMP3",
                "MS4A1", "CD79A", "CD79B", "PTGDS", "TRBC2"])
GENE_8 = to_ensg(['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1'])

scores_kr = pd.DataFrame(index=expr_kr.columns)
scores_kr["sig_score"] = sig_score
scores_kr["HLA_II"] = zmean(HLA_II, expr_kr)
scores_kr["HLA_I"] = zmean(HLA_I, expr_kr)
scores_kr["B_cell"] = zmean(B_CELL, expr_kr)
scores_kr["TLS"] = zmean(TLS, expr_kr)
scores_kr["g8_RAI"] = zmean(GENE_8, expr_kr)
print(f"  Korean samples: {len(scores_kr)}")
print(f"  sig_score mean={scores_kr['sig_score'].mean():.3f}, sd={scores_kr['sig_score'].std():.3f}")
print(f"  range: [{scores_kr['sig_score'].min():.3f}, {scores_kr['sig_score'].max():.3f}]")

# ============================================================
# 3. Bimodality + Hashimoto-like calls
# ============================================================
print("\n=== 3. Bimodality + Hashimoto-like ===")
from scipy.stats import skew, kurtosis
sk = skew(scores_kr["sig_score"])
ku = kurtosis(scores_kr["sig_score"])
bimod = (sk**2 + 1) / (ku + 3 * (len(scores_kr)-1)**2 / ((len(scores_kr)-2)*(len(scores_kr)-3)))
print(f"  bimodality coef: {bimod:.3f}")

gmm = GaussianMixture(n_components=2, random_state=42).fit(scores_kr[["sig_score"]].values)
g_lab = gmm.predict(scores_kr[["sig_score"]].values)
hashi_comp = int(gmm.means_.flatten().argmax())
scores_kr["hashi_GMM"] = (g_lab == hashi_comp).astype(int)

def otsu(v, n_bins=64):
    h, e = np.histogram(v, bins=n_bins)
    p = h / h.sum()
    c = (e[:-1] + e[1:]) / 2
    best_var = -1; best_t = c[0]
    for i in range(1, n_bins):
        w0, w1 = p[:i].sum(), p[i:].sum()
        if w0 == 0 or w1 == 0:
            continue
        m0 = (p[:i] * c[:i]).sum() / w0
        m1 = (p[i:] * c[i:]).sum() / w1
        var = w0 * w1 * (m0 - m1)**2
        if var > best_var:
            best_var = var; best_t = c[i]
    return best_t

otsu_thr = otsu(scores_kr["sig_score"].values)
scores_kr["hashi_otsu"] = (scores_kr["sig_score"] >= otsu_thr).astype(int)
for pct in [10, 20, 30]:
    thr = scores_kr["sig_score"].quantile(1 - pct/100)
    scores_kr[f"hashi_top{pct}"] = (scores_kr["sig_score"] >= thr).astype(int)

print(f"  GMM hashi+: {scores_kr['hashi_GMM'].sum()}/{len(scores_kr)} ({100*scores_kr['hashi_GMM'].mean():.1f}%)")
print(f"  Otsu hashi+: {scores_kr['hashi_otsu'].sum()}/{len(scores_kr)} ({100*scores_kr['hashi_otsu'].mean():.1f}%)  [thr={otsu_thr:.3f}]")
print(f"  Top10/20/30: {scores_kr['hashi_top10'].sum()}/{scores_kr['hashi_top20'].sum()}/{scores_kr['hashi_top30'].sum()}")

# ============================================================
# 4. K2-derived DM call for Korean? — load Lee predictions if any
# ============================================================
print("\n=== 4. Korean DM cluster proxy ===")
# Lee 2024 doesn't have direct DM call in our pipeline. Approximate via 8-gene RAI score:
# - If we trained TCGA classifier on log2 FPKM 8-gene, score should map → P_DM2
# But we have a published p_DM2 calling? Let me check.
lee_pred_paths = [
    PROJ / "results/v17_korean/Lee2024_predictions.tsv",
    PROJ / "results/v17_korean/GSE213647_panel_score.tsv",
]
lee_pred = None
for p in lee_pred_paths:
    if p.exists():
        lee_pred = pd.read_csv(p, sep="\t")
        print(f"  Found Lee prediction: {p}")
        print(f"  cols: {lee_pred.columns.tolist()[:10]}")
        break
if lee_pred is None:
    print(f"  No Lee DM prediction file. Building proxy from 8-gene RAI score quartile:")
    # Q1 (lowest 25% RAI) ≈ DM2-like, Q4 (highest) ≈ DM1-like
    q25 = scores_kr["g8_RAI"].quantile(0.25)
    q75 = scores_kr["g8_RAI"].quantile(0.75)
    scores_kr["DM_proxy"] = pd.cut(scores_kr["g8_RAI"],
                                      bins=[-np.inf, q25, q75, np.inf],
                                      labels=["DM2-like", "Mid", "DM1-like"])
    print(f"  proxy distribution: {scores_kr['DM_proxy'].value_counts().to_dict()}")
else:
    # If we have actual predictions, merge
    if "g8_score" in lee_pred.columns or "DM_call" in lee_pred.columns or "p_DM2" in lee_pred.columns:
        if "DM_call" in lee_pred.columns:
            lee_pred_idx = lee_pred.set_index(lee_pred.columns[0])["DM_call"]
            scores_kr["DM_call"] = lee_pred_idx.reindex(scores_kr.index)
            print(scores_kr["DM_call"].value_counts())

# ============================================================
# 5. Hashimoto-like × DM proxy cross-tab
# ============================================================
print("\n=== 5. Hashimoto-like × DM proxy ===")
if "DM_proxy" in scores_kr.columns:
    for label in ["hashi_otsu", "hashi_top20", "hashi_GMM"]:
        ct = pd.crosstab(scores_kr["DM_proxy"], scores_kr[label])
        print(f"\n  {label}:")
        print(ct)
elif "DM_call" in scores_kr.columns:
    for label in ["hashi_otsu", "hashi_top20", "hashi_GMM"]:
        ct = pd.crosstab(scores_kr["DM_call"], scores_kr[label])
        if 1 in ct.columns and 0 in ct.columns:
            try:
                t = ct[[0,1]].values
                if t.shape == (2, 2):
                    odds, p = stats.fisher_exact(t)
                    print(f"  {label}: OR={odds:.3f}, p={p:.3g}")
                    print(ct)
            except Exception:
                pass

# ============================================================
# 6. Cross-cohort comparison summary
# ============================================================
print("\n=== 6. Cross-cohort Hashimoto-like prevalence summary ===")
print(f"  GSE286332 (n=18): PTC+HT 9/9 (100%)")
tcga_summary = json.loads((PROJ / "results/d4p2_tcga_hashimoto_signature/D4P2_summary.json").read_text())
print(f"  TCGA-THCA (n=500): GMM {tcga_summary['hashi_calls']['GMM']}/500 ({100*tcga_summary['hashi_calls']['GMM']/500:.1f}%)")
print(f"  TCGA-THCA (n=500): Otsu {tcga_summary['hashi_calls']['Otsu']}/500 ({100*tcga_summary['hashi_calls']['Otsu']/500:.1f}%)")
print(f"  Korean GSE213647 (n={len(scores_kr)}): GMM {scores_kr['hashi_GMM'].sum()}/{len(scores_kr)} ({100*scores_kr['hashi_GMM'].mean():.1f}%)")
print(f"  Korean GSE213647 (n={len(scores_kr)}): Otsu {scores_kr['hashi_otsu'].sum()}/{len(scores_kr)} ({100*scores_kr['hashi_otsu'].mean():.1f}%)")

# Save
scores_kr.to_csv(RES / "korean_GSE213647_hashimoto_scores.tsv", sep="\t")

summary = {
    "n_samples": len(scores_kr),
    "n_signature_up_in_korean": len(up_kr),
    "n_signature_dn_in_korean": len(dn_kr),
    "bimodality_coef": round(float(bimod), 3),
    "hashi_calls": {
        "GMM": int(scores_kr["hashi_GMM"].sum()),
        "Otsu": int(scores_kr["hashi_otsu"].sum()),
        "top10": int(scores_kr["hashi_top10"].sum()),
        "top20": int(scores_kr["hashi_top20"].sum()),
        "top30": int(scores_kr["hashi_top30"].sum()),
    },
    "hashi_pct": {
        "GMM": round(100 * scores_kr["hashi_GMM"].mean(), 1),
        "Otsu": round(100 * scores_kr["hashi_otsu"].mean(), 1),
    },
    "cross_cohort": {
        "TCGA_GMM_pct": round(100 * tcga_summary["hashi_calls"]["GMM"] / 500, 1),
        "TCGA_Otsu_pct": round(100 * tcga_summary["hashi_calls"]["Otsu"] / 500, 1),
        "Korean_GMM_pct": round(100 * scores_kr["hashi_GMM"].mean(), 1),
        "Korean_Otsu_pct": round(100 * scores_kr["hashi_otsu"].mean(), 1),
    },
}
(RES / "D8B_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n✓ Outputs to {RES}")
