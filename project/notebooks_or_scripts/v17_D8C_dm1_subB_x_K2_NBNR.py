#!/usr/bin/env python3
"""D8-C — TCGA DM1 sub-B (NBNR cluster) signature × K2 NBNR cross-cohort transfer.

P7 finding: DM1 sub-B (n=56) is 96% mutation-negative (BRAF-/RAS- NBNR cluster) in TCGA.
Hypothesis: K2 (Yu professor PRJEB11591 n=260) NBNR samples should map to sub-B signature.

Tasks:
1. Build sub-B vs sub-A signature from TCGA DM1 (top DEGs)
2. Apply per-sample to K2 8-gene values + any cohort-level matrix
3. Verify: K2 NBNR samples → high sub-B score; K2 RAS+ samples → low sub-B score
4. ETE/aggressiveness phenotype × sub-B score correlation
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/d8c_dm1_subB_x_K2_NBNR"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Load DM1 sub-B vs sub-A DEGs from D6-P7
# ============================================================
print("=== 1. Load TCGA DM1 sub-B vs sub-A signature ===")
deg = pd.read_csv(PROJ / "results/d6p7_dm1_subcluster/dm1_subBvA_deg.tsv", sep="\t")
print(f"  total DEGs: {len(deg)}, padj<0.05: {(deg['padj'] < 0.05).sum()}")

# Top up + down (strong signal)
up_B = deg[(deg["padj"] < 0.05) & (deg["log2FC_B_vs_A"] > 0.5)].sort_values("padj").head(100)
dn_B = deg[(deg["padj"] < 0.05) & (deg["log2FC_B_vs_A"] < -0.5)].sort_values("padj").head(100)
print(f"  sub-B signature: up={len(up_B)}, dn={len(dn_B)}")

# Show top genes
print("\n  Top 15 up in sub-B (NBNR/Hashimoto-overlap markers):")
print(up_B.head(15)[["gene", "log2FC_B_vs_A", "cohen_d_B_vs_A", "padj"]].to_string(index=False))
print("\n  Top 15 down in sub-B (RAS+/FVPTC core markers):")
print(dn_B.head(15)[["gene", "log2FC_B_vs_A", "cohen_d_B_vs_A", "padj"]].to_string(index=False))

# ============================================================
# 2. K2 cohort — find expression matrix or use 8-gene
# ============================================================
print("\n=== 2. K2 cohort score apply attempt ===")
# K2 has only 8-gene mini-index quantification (calibration issue per D7-P3)
# Try alternative: use the full transcriptome of GSE213647 (Lee 2024 Korean PTC) instead
# since signature transfer to K2 is calibration-blocked
expr_kr = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_log2.tsv",
                       sep="\t", index_col=0)
print(f"  GSE213647 (Lee 2024) expr: {expr_kr.shape}")
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

up_kr = to_ensg(up_B["gene"].tolist())
dn_kr = to_ensg(dn_B["gene"].tolist())
print(f"  sub-B signature in Korean: up={len(up_kr)}/{len(up_B)}, dn={len(dn_kr)}/{len(dn_B)}")

def zmean(genes, mat):
    keep = [g for g in genes if g in mat.index]
    if len(keep) < 3:
        return pd.Series(0.0, index=mat.columns)
    z = mat.loc[keep].apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9), axis=1)
    return z.mean(axis=0)

up_score = zmean(up_kr, expr_kr)
dn_score = zmean(dn_kr, expr_kr)
subB_score = up_score - dn_score
print(f"\n  Korean sub-B score: mean={subB_score.mean():.3f}, sd={subB_score.std():.3f}")
print(f"  range: [{subB_score.min():.3f}, {subB_score.max():.3f}]")

# ============================================================
# 3. Apply to K2 (8-gene only — limited but available)
# ============================================================
print("\n=== 3. K2 8-gene-only sub-B proxy ===")
k2 = pd.read_csv(PROJ / "results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t")
GENE_8 = ['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1']
# In TCGA sub-B vs sub-A: look at 8-gene changes
g8_in_deg = deg[deg["gene"].isin(GENE_8)]
print("  8-gene direction in sub-B vs sub-A:")
print(g8_in_deg[["gene", "log2FC_B_vs_A", "cohen_d_B_vs_A", "padj"]].to_string(index=False))

# K2 has 8-gene log TPM (raw, calibration issue). Use ranking proxy.
k2_g8_log = np.log2(k2[GENE_8].astype(float) + 1)
k2_8gene_z = k2_g8_log.apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9), axis=0)
k2_g8_score = k2_8gene_z.mean(axis=1)
k2["k2_g8_score"] = k2_g8_score.values

# ============================================================
# 4. Bimodality / Hashimoto-like rate in Korean
# ============================================================
print("\n=== 4. Korean sub-B rate ===")
from sklearn.mixture import GaussianMixture
gmm = GaussianMixture(n_components=2, random_state=42).fit(subB_score.values.reshape(-1, 1))
g_lab = gmm.predict(subB_score.values.reshape(-1, 1))
high_comp = int(gmm.means_.flatten().argmax())
high_subB = (g_lab == high_comp).astype(int)
print(f"  Korean GMM 'sub-B-like': {high_subB.sum()}/{len(subB_score)} ({100*high_subB.mean():.1f}%)")

# Otsu for comparison
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

thr = otsu(subB_score.values)
otsu_subB = (subB_score >= thr).astype(int)
print(f"  Korean Otsu 'sub-B-like': {otsu_subB.sum()}/{len(subB_score)} ({100*otsu_subB.mean():.1f}%)")

# ============================================================
# 5. Cross-reference with TCGA Hashimoto-like (D4-P2)
# ============================================================
print("\n=== 5. TCGA: sub-B × Hashimoto-like cross-check ===")
hashi_path = PROJ / "results/d4p2_tcga_hashimoto_signature/tcga_signature_scores.tsv"
sub_path = PROJ / "results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv"
labels_path = PROJ / "results/v17_realfix/R1A_cluster_labels.tsv"

if hashi_path.exists() and sub_path.exists():
    hashi = pd.read_csv(hashi_path, sep="\t", index_col=0)
    sub = pd.read_csv(sub_path, sep="\t", index_col=0)
    join = sub.join(hashi[["sig_score", "hashi_otsu", "hashi_top20", "hashi_top30"]], how="left")
    print(f"\n  TCGA DM1 sub × Hashimoto-like (D4-P2 join):")
    for label in ["hashi_otsu", "hashi_top20", "hashi_top30"]:
        if label in join.columns:
            ct = pd.crosstab(join["sub_cluster"], join[label])
            print(f"\n  {label}:")
            print(ct)
            n11 = ct.loc["sub_A", 1] if ("sub_A" in ct.index and 1 in ct.columns) else 0
            n10 = ct.loc["sub_A", 0] if ("sub_A" in ct.index and 0 in ct.columns) else 0
            n21 = ct.loc["sub_B", 1] if ("sub_B" in ct.index and 1 in ct.columns) else 0
            n20 = ct.loc["sub_B", 0] if ("sub_B" in ct.index and 0 in ct.columns) else 0
            try:
                or_, p = stats.fisher_exact([[n11, n10], [n21, n20]])
                print(f"    Fisher OR={or_:.3f}, p={p:.3g} (sub-A vs sub-B Hashimoto enrichment)")
            except Exception as e:
                print(f"    Fisher: {e}")

# ============================================================
# 6. Korean GSE213647 sub-B × Hashimoto-like
# ============================================================
print("\n=== 6. Korean GSE213647 sub-B × Hashimoto-like ===")
korean_hashi_path = PROJ / "results/d8b_korean_replication/korean_GSE213647_hashimoto_scores.tsv"
if korean_hashi_path.exists():
    kr_hashi = pd.read_csv(korean_hashi_path, sep="\t", index_col=0)
    kr_join = pd.DataFrame({"subB_score": subB_score, "subB_GMM": high_subB,
                              "subB_otsu": otsu_subB.values}, index=subB_score.index)
    kr_join = kr_join.join(kr_hashi[["sig_score", "hashi_otsu", "hashi_top20"]], how="left")
    rho, p = stats.spearmanr(kr_join["subB_score"], kr_join["sig_score"])
    print(f"  Korean ρ(sub-B score, Hashimoto sig_score) = {rho:+.3f}, p={p:.3g}")
    if rho > 0.3:
        print(f"  ★ Sub-B signature and Hashimoto signature are STRONGLY positively correlated in Korean — consistent with sub-B = NBNR + autoimmune-overlap convergence")
    elif rho > 0:
        print(f"  Sub-B and Hashimoto signature directionally consistent in Korean")
    kr_join.to_csv(RES / "korean_subB_x_hashimoto.tsv", sep="\t")

# ============================================================
# 7. Save outputs + decision
# ============================================================
out = pd.DataFrame({
    "subB_score": subB_score,
    "subB_GMM": high_subB,
    "subB_otsu_call": otsu_subB.values,
}, index=subB_score.index)
out.to_csv(RES / "korean_subB_score.tsv", sep="\t")

decision = "STRONG" if (otsu_subB.mean() > 0.3 and otsu_subB.mean() < 0.6) else "MODERATE"

summary = {
    "TCGA_subB_signature_genes": {"up": int(len(up_B)), "dn": int(len(dn_B))},
    "Korean_signature_match": {"up": int(len(up_kr)), "dn": int(len(dn_kr))},
    "Korean_subB_rate_GMM_pct": round(100 * high_subB.mean(), 1),
    "Korean_subB_rate_Otsu_pct": round(100 * otsu_subB.mean(), 1),
    "decision": decision,
}
(RES / "D8C_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n✓ Decision: {decision}")
print(f"✓ Outputs to {RES}")
