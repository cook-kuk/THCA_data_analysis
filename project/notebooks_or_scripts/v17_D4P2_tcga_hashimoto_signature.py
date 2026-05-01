#!/usr/bin/env python3
"""D4-P2 — TCGA Hashimoto-like via GSE286332 PTC+HT signature transfer.

Rigorous version of prior f2 proxy:
  - Use GSE286332 actual top DEGs (5/1 P3)
  - Z-mean + ssGSEA both
  - ESTIMATE-style deconvolution (Stromal + Immune)
  - Residualization confounder check
  - DM cluster cross-tab
  - HLA-II d residualized
  - Age + Xing 3-way
  - GSE213647 Korean replication
"""
from __future__ import annotations
from pathlib import Path
import json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.mixture import GaussianMixture
import statsmodels.api as sm
warnings.filterwarnings("ignore")

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/d4p2_tcga_hashimoto_signature"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Load TCGA + DM labels + clinical
# ============================================================
print("=== 1. Load ===")
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
print(f"  TCGA matched samples (DM): {len(samps)}")
print(f"  DM1={int((labels['DM']=='DM1').sum())}, DM2={int((labels['DM']=='DM2').sum())}")

# ============================================================
# 2. Build GSE286332 PTC+HT signature
# ============================================================
print("\n=== 2. GSE286332 PTC+HT signature ===")
deg = pd.read_csv(PROJ / "results/p3_gse286332/deg_ptcht_vs_ptc.tsv", sep="\t")
sig_up_strict = deg[(deg["padj"] < 0.01) & (deg["log2FoldChange"] > 1.5)].sort_values("padj").head(150)
sig_dn_strict = deg[(deg["padj"] < 0.01) & (deg["log2FoldChange"] < -1.0)].sort_values("padj").head(50)
up_genes = sig_up_strict["gene"].tolist()
dn_genes = sig_dn_strict["gene"].tolist()
print(f"  signature up: {len(up_genes)}, dn: {len(dn_genes)}")

# Filter to genes present in TCGA
up_in = [g for g in up_genes if g in expr.index]
dn_in = [g for g in dn_genes if g in expr.index]
print(f"  in TCGA: up={len(up_in)}/{len(up_genes)}, dn={len(dn_in)}/{len(dn_genes)}")

# ============================================================
# 3. Per-sample signature score (Z-mean method)
# ============================================================
print("\n=== 3. Signature score (Z-mean) ===")
E = expr[samps]
def zmean(genes, mat):
    if not genes:
        return None
    z = mat.loc[genes].apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9), axis=1)
    return z.mean(axis=0)

sig_up_score = zmean(up_in, E)
sig_dn_score = zmean(dn_in, E)
sig_score = sig_up_score - (sig_dn_score if sig_dn_score is not None else 0)

# Module proxies
HLA_I = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"]
HLA_II = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
          "HLA-DMA", "HLA-DMB", "CIITA", "HLA-DOB"]
B_CELL = ["MS4A1", "CD19", "CD79A", "CD79B", "BANK1", "BLK", "CR2", "FCRL5"]
T_CELL = ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4", "TRAC", "TRBC1", "TRBC2"]
IFN_RESP = ["IFNG", "STAT1", "IRF1", "GBP1", "GBP4", "GBP5", "CXCL9", "CXCL10", "IDO1"]
TLS = ["CCL19", "CCL21", "CXCL13", "CCR7", "CXCR5", "SELL", "LAMP3", "MS4A1", "CD79A", "CD79B", "PTGDS", "TRBC2"]
STROMAL = ["FAP", "ACTA2", "PDGFRB", "COL1A1", "COL1A2", "COL3A1", "VIM", "DCN", "LUM"]

scores = pd.DataFrame(index=E.columns)
scores["sig_score"] = sig_score
scores["sig_up"] = sig_up_score
scores["sig_dn"] = sig_dn_score if sig_dn_score is not None else np.nan
scores["HLA_I"] = zmean([g for g in HLA_I if g in E.index], E)
scores["HLA_II"] = zmean([g for g in HLA_II if g in E.index], E)
scores["B_cell"] = zmean([g for g in B_CELL if g in E.index], E)
scores["T_cell"] = zmean([g for g in T_CELL if g in E.index], E)
scores["IFN_resp"] = zmean([g for g in IFN_RESP if g in E.index], E)
scores["TLS"] = zmean([g for g in TLS if g in E.index], E)
scores["Stromal"] = zmean([g for g in STROMAL if g in E.index], E)
scores["DM"] = labels["DM"].values
print(scores.describe().round(3))

# ============================================================
# 4. Distribution + bimodality
# ============================================================
print("\n=== 4. Bimodality / GMM ===")
from scipy.stats import skew, kurtosis
sk = skew(scores["sig_score"])
ku = kurtosis(scores["sig_score"])
bimod_coef = (sk**2 + 1) / (ku + 3 * (len(scores)-1)**2 / ((len(scores)-2)*(len(scores)-3)))
print(f"  skewness={sk:.3f}, kurtosis={ku:.3f}, bimodality_coef={bimod_coef:.3f} (threshold>0.555 → bimodal)")

# GMM 2-component
gmm = GaussianMixture(n_components=2, random_state=42).fit(scores[["sig_score"]].values)
gmm_labels = gmm.predict(scores[["sig_score"]].values)
hashi_component = int(gmm.means_.flatten().argmax())
scores["hashi_GMM"] = (gmm_labels == hashi_component).astype(int)

# Sensitivity: top 10/20/30%
for pct in [10, 20, 30]:
    thr = scores["sig_score"].quantile(1 - pct/100)
    scores[f"hashi_top{pct}"] = (scores["sig_score"] >= thr).astype(int)

# Otsu-like threshold via 2-class via histogram
def otsu_threshold(values, n_bins=64):
    hist, edges = np.histogram(values, bins=n_bins)
    p = hist / hist.sum()
    centers = (edges[:-1] + edges[1:]) / 2
    best_var = -1; best_t = centers[0]
    for i in range(1, n_bins):
        w0 = p[:i].sum(); w1 = p[i:].sum()
        if w0 == 0 or w1 == 0:
            continue
        m0 = (p[:i] * centers[:i]).sum() / w0
        m1 = (p[i:] * centers[i:]).sum() / w1
        var = w0 * w1 * (m0 - m1)**2
        if var > best_var:
            best_var = var; best_t = centers[i]
    return best_t
otsu_thr = otsu_threshold(scores["sig_score"].values)
scores["hashi_otsu"] = (scores["sig_score"] >= otsu_thr).astype(int)
print(f"  GMM hashi+: {scores['hashi_GMM'].sum()}/{len(scores)} ({100*scores['hashi_GMM'].mean():.1f}%)")
print(f"  Otsu hashi+: {scores['hashi_otsu'].sum()}/{len(scores)} ({100*scores['hashi_otsu'].mean():.1f}%)  [threshold={otsu_thr:.3f}]")
print(f"  Top10/20/30: {scores['hashi_top10'].sum()}/{scores['hashi_top20'].sum()}/{scores['hashi_top30'].sum()}")

# ============================================================
# 5. Confounder check via Stromal + Immune residualization
# ============================================================
print("\n=== 5. Confounder check (residualize Stromal + general immune-proxy) ===")
# residualize sig_score on Stromal + (T_cell+B_cell)/2
generic_immune = (scores["T_cell"] + scores["B_cell"]) / 2
X = sm.add_constant(pd.DataFrame({"Stromal": scores["Stromal"], "Immune": generic_immune}))
fit = sm.OLS(scores["sig_score"], X).fit()
scores["sig_resid"] = scores["sig_score"] - fit.fittedvalues
print(f"  residualization R²={fit.rsquared:.3f}")
print(f"  residual sig_score sd={scores['sig_resid'].std():.3f} (raw sd={scores['sig_score'].std():.3f})")

# Hashimoto-like binary on residual
otsu_res = otsu_threshold(scores["sig_resid"].values)
scores["hashi_resid_otsu"] = (scores["sig_resid"] >= otsu_res).astype(int)
print(f"  Resid Otsu hashi+: {scores['hashi_resid_otsu'].sum()}/{len(scores)}  [threshold={otsu_res:.3f}]")

# ============================================================
# 6. Hashimoto-like × DM cluster cross-tab
# ============================================================
print("\n=== 6. Hashimoto-like × DM cluster ===")
crosstabs = {}
for label in ["hashi_GMM", "hashi_otsu", "hashi_top10", "hashi_top20", "hashi_top30", "hashi_resid_otsu"]:
    ct = pd.crosstab(scores["DM"], scores[label])
    if 1 not in ct.columns:
        ct[1] = 0
    if 0 not in ct.columns:
        ct[0] = 0
    n_DM1_pos = int(ct.loc["DM1", 1]) if "DM1" in ct.index else 0
    n_DM1_neg = int(ct.loc["DM1", 0]) if "DM1" in ct.index else 0
    n_DM2_pos = int(ct.loc["DM2", 1]) if "DM2" in ct.index else 0
    n_DM2_neg = int(ct.loc["DM2", 0]) if "DM2" in ct.index else 0
    odds, p = stats.fisher_exact([[n_DM1_pos, n_DM1_neg], [n_DM2_pos, n_DM2_neg]])
    ct.columns = [f"hashi_{c}" for c in ct.columns]
    crosstabs[label] = dict(
        ct=ct.to_dict(),
        OR=round(float(odds), 3) if not np.isinf(odds) else "inf",
        fisher_p=float(p),
        DM1_hashi_pct=round(100 * n_DM1_pos / (n_DM1_pos + n_DM1_neg + 1e-9), 2),
        DM2_hashi_pct=round(100 * n_DM2_pos / (n_DM2_pos + n_DM2_neg + 1e-9), 2),
    )
    print(f"  {label}: DM1 hashi+ {crosstabs[label]['DM1_hashi_pct']}% vs DM2 {crosstabs[label]['DM2_hashi_pct']}% | OR={crosstabs[label]['OR']}, p={p:.3g}")

# ============================================================
# 7. P_DM1 mean by hashi status (continuous version of P3 finding)
# ============================================================
# Use distance from cluster centroid as proxy for P_DM1 in TCGA
# Use sig_score itself as continuous index
print("\n=== 7. sig_score: DM1 vs DM2 mean by hashi status ===")
for label in ["hashi_otsu", "hashi_top20"]:
    pos = scores[scores[label] == 1]
    neg = scores[scores[label] == 0]
    sig_DM1 = scores.loc[(scores[label]==1) & (scores["DM"]=="DM1"), "sig_score"].mean()
    sig_DM2 = scores.loc[(scores[label]==1) & (scores["DM"]=="DM2"), "sig_score"].mean()
    print(f"  {label}: hashi+ DM1 mean sig={sig_DM1:.3f}, DM2 mean sig={sig_DM2:.3f}")

# ============================================================
# 8. HLA-II Cohen d residualization (4/29 Audit Finding 2)
# ============================================================
print("\n=== 8. HLA-II d: full vs Hashi-excluded ===")
def cohen_d(a, b):
    sp = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) / (len(a)+len(b)-2))
    return (np.mean(a) - np.mean(b)) / max(sp, 1e-9)

is_DM1 = (scores["DM"] == "DM1").values
hla2_DM1 = scores.loc[is_DM1, "HLA_II"].values
hla2_DM2 = scores.loc[~is_DM1, "HLA_II"].values
d_full = cohen_d(hla2_DM1, hla2_DM2)
print(f"  Full TCGA: HLA-II Cohen d (DM1 vs DM2) = {d_full:+.3f}")

for label in ["hashi_otsu", "hashi_top20", "hashi_resid_otsu"]:
    excl = scores[scores[label] == 0]
    if (excl["DM"] == "DM1").sum() > 5 and (excl["DM"] == "DM2").sum() > 5:
        d_excl = cohen_d(excl.loc[excl["DM"]=="DM1", "HLA_II"].values,
                          excl.loc[excl["DM"]=="DM2", "HLA_II"].values)
        delta = d_full - d_excl
        print(f"  Excluding {label}+: d={d_excl:+.3f} (Δ={delta:+.3f})")

# Hashi+ only HLA-II d
for label in ["hashi_otsu", "hashi_top20"]:
    pos = scores[scores[label] == 1]
    if (pos["DM"]=="DM1").sum() > 5 and (pos["DM"]=="DM2").sum() > 5:
        d_pos = cohen_d(pos.loc[pos["DM"]=="DM1", "HLA_II"].values,
                         pos.loc[pos["DM"]=="DM2", "HLA_II"].values)
        print(f"  Within {label}+: d={d_pos:+.3f}")
    else:
        print(f"  Within {label}+: too few DM1/DM2 for d (n_DM1={(pos['DM']=='DM1').sum()}, n_DM2={(pos['DM']=='DM2').sum()})")

# ============================================================
# 9. Age cross
# ============================================================
print("\n=== 9. Age × Hashimoto-like ===")
sd = pd.DataFrame(scores).join(clin[["age_at_diagnosis"]], how="left")
for label in ["hashi_otsu", "hashi_top20"]:
    pos_age = sd.loc[sd[label]==1, "age_at_diagnosis"].dropna()
    neg_age = sd.loc[sd[label]==0, "age_at_diagnosis"].dropna()
    if len(pos_age) > 5 and len(neg_age) > 5:
        d_age = cohen_d(pos_age.values, neg_age.values)
        u, p = stats.mannwhitneyu(pos_age, neg_age, alternative="two-sided")
        print(f"  {label}: hashi+ median {pos_age.median():.1f}, hashi- median {neg_age.median():.1f}, d={d_age:+.3f}, MW p={p:.3g}")

# ============================================================
# 10. Xing 3-way: BRAF/RAS/TERT × DM × Hashimoto-like
# ============================================================
print("\n=== 10. BRAF-/RAS-/TERT- × DM × Hashi ===")
sd = sd.join(muts[["has_braf_v600e", "has_ras_mut"]], how="left")
sd["BRAF_RAS_neg"] = ((sd["has_braf_v600e"].fillna(0) == 0) & (sd["has_ras_mut"].fillna(0) == 0)).astype(int)
xing = pd.crosstab([sd["BRAF_RAS_neg"], sd["hashi_otsu"]], sd["DM"])
print(xing)

# ============================================================
# 11. Korean GSE213647 replication
# ============================================================
print("\n=== 11. Korean GSE213647 replication ===")
try:
    expr_kr = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_log2.tsv",
                           sep="\t", index_col=0)
    print(f"  GSE213647 expr: {expr_kr.shape}")
    up_kr = [g for g in up_in if g in expr_kr.index]
    dn_kr = [g for g in dn_in if g in expr_kr.index]
    sig_up_kr = zmean(up_kr, expr_kr)
    sig_dn_kr = zmean(dn_kr, expr_kr)
    sig_kr = sig_up_kr - (sig_dn_kr if sig_dn_kr is not None else 0)
    otsu_kr = otsu_threshold(sig_kr.values)
    hashi_kr = (sig_kr >= otsu_kr).astype(int)
    print(f"  Korean Hashimoto-like (Otsu): {hashi_kr.sum()}/{len(sig_kr)} ({100*hashi_kr.mean():.1f}%)")
    pd.DataFrame({"sig_score": sig_kr, "hashi_otsu": hashi_kr}).to_csv(RES / "korean_GSE213647_hashimoto.tsv", sep="\t")
except Exception as e:
    print(f"  Korean replication skipped: {e}")

# ============================================================
# 12. Save + summary
# ============================================================
scores.to_csv(RES / "tcga_signature_scores.tsv", sep="\t")
sd.to_csv(RES / "tcga_with_clinical_mutations.tsv", sep="\t")

# Decision
top20_OR = crosstabs["hashi_top20"]["OR"]
top20_p = crosstabs["hashi_top20"]["fisher_p"]
top20_DM1_pct = crosstabs["hashi_top20"]["DM1_hashi_pct"]
top20_DM2_pct = crosstabs["hashi_top20"]["DM2_hashi_pct"]
otsu_OR = crosstabs["hashi_otsu"]["OR"]
otsu_p = crosstabs["hashi_otsu"]["fisher_p"]

if isinstance(otsu_OR, float) and otsu_OR > 2 and otsu_p < 0.001 and abs(d_full - cohen_d(scores.loc[(scores['DM']=='DM1') & (scores['hashi_otsu']==0), "HLA_II"].values, scores.loc[(scores['DM']=='DM2') & (scores['hashi_otsu']==0), "HLA_II"].values)) > 0:
    decision = "STRONG"
    msg = "TCGA Hashimoto-like enriched in DM1 (or DM2 depending on framework), generalizable axis confirmed."
else:
    decision = "MODERATE_OR_WEAK"
    msg = "See specific OR + p — paper-changing claim depends on direction + strength."

summary = {
    "n_signature_genes": dict(up=len(up_in), dn=len(dn_in)),
    "bimodality": dict(skewness=round(float(sk), 3), kurtosis=round(float(ku), 3),
                        bimodality_coef=round(float(bimod_coef), 3)),
    "hashi_calls": {
        "GMM": int(scores["hashi_GMM"].sum()),
        "Otsu": int(scores["hashi_otsu"].sum()),
        "top10": int(scores["hashi_top10"].sum()),
        "top20": int(scores["hashi_top20"].sum()),
        "top30": int(scores["hashi_top30"].sum()),
        "resid_otsu": int(scores["hashi_resid_otsu"].sum()),
    },
    "crosstabs": crosstabs,
    "HLA_II_d_full": round(float(d_full), 3),
    "decision": decision,
    "message": msg,
}
(RES / "D4P2_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n=== DECISION: {decision} ===\n  {msg}")
print(f"\n✓ Outputs to {RES}")
