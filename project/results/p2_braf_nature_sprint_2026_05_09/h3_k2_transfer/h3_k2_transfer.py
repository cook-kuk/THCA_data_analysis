"""H3: HT-13 RNA panel external cohort transfer to Korean data.

Decision after data search:
  - True K2 (PRJEB11591/Yoo2016 n=260) RNA-seq: only 8-gene mini-index kallisto
    quants exist; HT-13 genes (HLA-DR/DP/DQ A/B1, CD79A/B, MS4A1, AICDA, CXCL13,
    CCR6, IFNG) absent. Full-transcriptome re-quant of 33GB FASTQ would be
    multi-hour. Out of 30-min sprint budget.
  - GSE213647 (Lee 2024 Nat Commun, Korean SNUBH/CNUH/KRIBB, n=632): full
    bulk RNA-seq counts available. Use as Korean external cohort.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/h3_k2_transfer")
OUT.mkdir(parents=True, exist_ok=True)

HT13 = [
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "CD79A", "CD79B", "MS4A1", "AICDA", "CXCL13", "CCR6", "IFNG",
]

# --- 1. Korean cohort: GSE213647 expression -----------------------------------
print("[1/5] Loading GSE213647 counts ...")
counts = pd.read_csv(
    "/data/thca/v17_korean/GSE213647/expression_counts.tsv.gz",
    sep="\t",
    index_col=0,
)
print(f"   counts shape: {counts.shape}  (genes x samples)")

# Map ENSG (versioned) to gene symbol via gencode.v44 GTF
print("[2/5] Building ENSG -> symbol map from gencode v44 ...")
import re
gtf_path = "/data/thca/reference_kallisto/gencode.v44.basic.annotation.gtf"
ensg2sym: dict[str, str] = {}
needed = set(HT13)
with open(gtf_path) as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        if "\tgene\t" not in ln:
            continue
        m = re.search(r'gene_id "([^"]+)"', ln)
        s = re.search(r'gene_name "([^"]+)"', ln)
        if not m or not s:
            continue
        if s.group(1) in needed or s.group(1) in {
            "TG", "TPO", "SLC5A5", "PAX8", "NKX2-1", "FOXE1", "DIO1", "TSHR",
        }:
            ensg2sym[m.group(1)] = s.group(1)

sym2ensg = {v: k for k, v in ensg2sym.items()}
# strip version for matching (GSE213647 uses different ENSG versions than gencode v44)
sym2ensg_noversion = {v: k.split(".")[0] for k, v in ensg2sym.items()}
print(f"   resolved {len(sym2ensg)} symbols")

# build GSE213647 ENSG-noversion -> full index
counts_index_noversion = pd.Series(counts.index, index=[i.split(".")[0] for i in counts.index])

# subset to HT-13 genes
ht13_rows = []
missing = []
for g in HT13:
    ensg_nv = sym2ensg_noversion.get(g)
    if ensg_nv is None or ensg_nv not in counts_index_noversion.index:
        missing.append(g)
        continue
    full_ensg = counts_index_noversion.loc[ensg_nv]
    if isinstance(full_ensg, pd.Series):
        full_ensg = full_ensg.iloc[0]
    ht13_rows.append((g, counts.loc[full_ensg]))
print(f"   HT-13 found: {len(ht13_rows)} / 13   missing={missing}")

ht_df = pd.DataFrame({g: row for g, row in ht13_rows})
# normalize: log2(cpm + 1)
lib = counts.sum(axis=0)
cpm = (ht_df.T / lib * 1e6).T
log_cpm = np.log2(cpm + 1.0)

# attach metadata
meta = pd.read_csv(
    "/data/thca/v17_korean/GSE213647/sample_sheet_clinical.tsv", sep="\t"
).set_index("gsm")
log_cpm = log_cpm.loc[meta.index.intersection(log_cpm.index)]
meta = meta.loc[log_cpm.index]
print(f"   joined meta: {log_cpm.shape}")

# subB_GMM as DM1 proxy (DM1 sub-B = TCGA dedifferentiated NBNR-equivalent;
# memory v17_D6P7_dm1_subB_NBNR)
sub = pd.read_csv(
    "/data/thca/repo_results/d8c_dm1_subB_x_K2_NBNR/korean_subB_score.tsv",
    sep="\t",
    index_col=0,
)
sub.index = sub.index.astype(str)
hashi = pd.read_csv(
    "/data/thca/repo_results/d8b_korean_replication/korean_GSE213647_hashimoto_scores.tsv",
    sep="\t",
    index_col=0,
)
hashi.index = hashi.index.astype(str)

print("[3/5] Within-cohort z-score + HT score on Korean cohort ...")
# tumors only (PTC + PDTC + ATC, drop normals)
tumor_mask = meta["tissue_type"].isin(["PTC", "PDTC", "ATC"])
print(f"   tumor n = {int(tumor_mask.sum())}  normal n = {int((~tumor_mask).sum())}")

# within-cohort z (all samples, then restrict)
z_full = (log_cpm - log_cpm.mean()) / log_cpm.std(ddof=0)
ht_score_full = z_full.mean(axis=1)
results = pd.DataFrame({
    "ht_score_z": ht_score_full,
    "tissue_type": meta["tissue_type"],
    "histology": meta["histology"],
    "age": meta["age"],
    "sex": meta["sex"],
})
results = results.join(sub[["subB_score", "subB_GMM"]], how="left")
results = results.join(hashi[["hashi_GMM", "hashi_otsu", "sig_score"]], how="left")
results.to_csv(OUT / "h3_k2_results.tsv", sep="\t")

# --- 4. Korean unsupervised + label-cross-tab ---------------------------------
metrics: dict[str, object] = {"data_substrate": "GSE213647 (Lee 2024 Nat Commun, Korean PTC, n=632 full bulk RNA-seq)",
                              "true_K2_PRJEB11591_status": "8-gene mini-index quants only; HT-13 genes absent; would need full-transcriptome re-quant of 33GB FASTQ",
                              "n_HT13_genes_found": len(ht13_rows),
                              "missing_HT13": missing,
                              "n_total": int(log_cpm.shape[0]),
                              "n_tumor": int(tumor_mask.sum()),
                              "n_normal": int((~tumor_mask).sum())}

# Bimodality / Hartigan dip test
try:
    from scipy import stats as _sst
    # Gaussian mixture BIC: 1 vs 2 components on tumor HT scores
    from sklearn.mixture import GaussianMixture
    x_t = ht_score_full[tumor_mask].dropna().to_numpy().reshape(-1, 1)
    bic1 = GaussianMixture(n_components=1, random_state=0).fit(x_t).bic(x_t)
    bic2 = GaussianMixture(n_components=2, random_state=0).fit(x_t).bic(x_t)
    metrics["bic_1comp"] = float(bic1)
    metrics["bic_2comp"] = float(bic2)
    metrics["bic_delta_2_minus_1"] = float(bic2 - bic1)
    metrics["bic_supports_bimodal"] = bool(bic2 < bic1)
except Exception as e:
    metrics["bic_error"] = str(e)

# Cross-tab: HT score vs Hashimoto / vs DM1-sub-B / vs histology
def safe_auc(y, s):
    y = pd.Series(y).reset_index(drop=True)
    s = pd.Series(s).reset_index(drop=True)
    keep = y.notna() & s.notna()
    if keep.sum() < 10 or y[keep].nunique() < 2:
        return float("nan"), int(keep.sum())
    return float(roc_auc_score(y[keep], s[keep])), int(keep.sum())

# HT score discriminating Hashimoto (GMM)
auc_hashi, n_hashi = safe_auc(results["hashi_GMM"], results["ht_score_z"])
metrics["auc_HT_score_vs_hashi_GMM_all"] = auc_hashi
metrics["n_hashi_pair"] = n_hashi

# HT score discriminating DM1-sub-B (tumors only; sub-B = DM1-like)
res_t = results[tumor_mask]
auc_subB_t, n_subB_t = safe_auc(res_t["subB_GMM"], res_t["ht_score_z"])
metrics["auc_HT_score_vs_subB_GMM_tumor"] = auc_subB_t
metrics["n_subB_tumor_pair"] = n_subB_t

# HT score: tumor vs normal (sanity baseline)
auc_tvn, _ = safe_auc((meta["tissue_type"] != "Normal").astype(int), ht_score_full)
metrics["auc_HT_score_vs_tumor_vs_normal_all"] = auc_tvn

# HT score: ATC/PDTC vs PTC (dediff axis)
res_t2 = results[results["tissue_type"].isin(["PTC", "ATC", "PDTC"])].copy()
res_t2["dediff"] = res_t2["tissue_type"].isin(["ATC", "PDTC"]).astype(int)
auc_dediff, n_dediff = safe_auc(res_t2["dediff"], res_t2["ht_score_z"])
metrics["auc_HT_score_vs_ATCPDTC_vs_PTC"] = auc_dediff
metrics["n_dediff_pair"] = n_dediff

# Logistic regression CV on Korean tumors with subB_GMM as DM1 label
print("[4/5] LogReg(C=0.5) 5-fold CV on Korean tumors, label=subB_GMM ...")
mask_t = tumor_mask & results["subB_GMM"].notna()
X = z_full.loc[mask_t].to_numpy()  # within-cohort z (already)
y = results.loc[mask_t, "subB_GMM"].astype(int).to_numpy()
print(f"   n={len(y)}  pos(subB)={int(y.sum())}  neg={int(len(y)-y.sum())}")
oof = np.zeros(len(y))
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
for tr, te in skf.split(X, y):
    clf = LogisticRegression(C=0.5, max_iter=2000)
    clf.fit(X[tr], y[tr])
    oof[te] = clf.predict_proba(X[te])[:, 1]
auc_cv = float(roc_auc_score(y, oof))
metrics["k2_within_cohort_logreg_C0p5_5fold_AUC_subB"] = auc_cv

# bootstrap 95% CI
rng = np.random.default_rng(42)
boots = []
for _ in range(1000):
    idx = rng.integers(0, len(y), size=len(y))
    if len(np.unique(y[idx])) < 2:
        continue
    boots.append(roc_auc_score(y[idx], oof[idx]))
boots = np.asarray(boots)
metrics["k2_AUC_subB_95CI"] = [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]

# also try Hashimoto label (more direct HT-axis biology)
mask_h = tumor_mask & results["hashi_GMM"].notna()
Xh = z_full.loc[mask_h].to_numpy()
yh = results.loc[mask_h, "hashi_GMM"].astype(int).to_numpy()
oofh = np.zeros(len(yh))
for tr, te in skf.split(Xh, yh):
    clf = LogisticRegression(C=0.5, max_iter=2000)
    clf.fit(Xh[tr], yh[tr])
    oofh[te] = clf.predict_proba(Xh[te])[:, 1]
metrics["k2_within_cohort_logreg_C0p5_5fold_AUC_hashi"] = float(roc_auc_score(yh, oofh))
metrics["n_hashi_tumor"] = int(len(yh))

# --- 5. TCGA -> Korean transfer ----------------------------------------------
print("[5/5] TCGA -> GSE213647 transfer (within-each-cohort z-score) ...")
print("   loading TCGA z-score matrix ...")
# Read only HT-13 rows for memory
import csv
keep_genes = set(HT13)
header = None
rows = {}
with open("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv") as fh:
    rdr = csv.reader(fh, delimiter="\t")
    header = next(rdr)
    for row in rdr:
        if row[0] in keep_genes:
            rows[row[0]] = [float(v) if v not in ("", "NA") else np.nan for v in row[1:]]
        if len(rows) == len(HT13):
            break
print(f"   TCGA HT13 loaded: {len(rows)} / 13")
tcga_z = pd.DataFrame(rows, index=header[1:])  # samples x genes (TCGA z already gene-wise)
# TCGA z-score is per-gene across cohort -> already within-cohort z. good.
tcga_z = tcga_z[HT13]

# match TCGA samples to master (DM labels)
master = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv", sep="\t")
tcga_dm = master[["sample_id", "dm", "histology_subtype", "driver_anchor"]].copy()
# tcga_z index is "TCGA-XX-XXXX-01A" form; master sample_id needs check
print(f"   tcga_z first index: {tcga_z.index[0]}, master first sample_id: {tcga_dm['sample_id'].iloc[0]}")
# Align by short id (TCGA-XX-XXXX) prefix
tcga_z["short"] = tcga_z.index.str[:12]
tcga_dm["short"] = tcga_dm["sample_id"].str[:12]
joined = tcga_z.merge(tcga_dm[["short", "dm", "driver_anchor"]], on="short", how="left")
joined = joined[joined["dm"].isin(["DM1", "DM2"])].copy()
joined["y"] = (joined["dm"] == "DM1").astype(int)
print(f"   TCGA labelled: n={len(joined)} DM1={int(joined['y'].sum())} DM2={int(len(joined)-joined['y'].sum())}")

X_tcga = joined[HT13].to_numpy()
y_tcga = joined["y"].to_numpy()
clf_full = LogisticRegression(C=0.5, max_iter=2000)
clf_full.fit(X_tcga, y_tcga)
# train AUC sanity
metrics["tcga_train_AUC_HT13_DM1vDM2"] = float(roc_auc_score(y_tcga, clf_full.predict_proba(X_tcga)[:, 1]))

# 5-fold CV TCGA
oof_tcga = np.zeros(len(y_tcga))
for tr, te in skf.split(X_tcga, y_tcga):
    cc = LogisticRegression(C=0.5, max_iter=2000); cc.fit(X_tcga[tr], y_tcga[tr])
    oof_tcga[te] = cc.predict_proba(X_tcga[te])[:, 1]
metrics["tcga_5fold_CV_AUC_HT13_DM1vDM2"] = float(roc_auc_score(y_tcga, oof_tcga))

# BRAF subset CV
braf_mask = (joined["driver_anchor"] == "BRAF").to_numpy()
if braf_mask.sum() >= 20 and len(np.unique(y_tcga[braf_mask])) > 1:
    Xb = X_tcga[braf_mask]; yb = y_tcga[braf_mask]
    oofb = np.zeros(len(yb))
    skf_b = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    for tr, te in skf_b.split(Xb, yb):
        cc = LogisticRegression(C=0.5, max_iter=2000); cc.fit(Xb[tr], yb[tr])
        oofb[te] = cc.predict_proba(Xb[te])[:, 1]
    metrics["tcga_BRAF_5fold_CV_AUC_HT13"] = float(roc_auc_score(yb, oofb))
    metrics["tcga_BRAF_n"] = int(braf_mask.sum())
    metrics["tcga_BRAF_pos"] = int(yb.sum())

# Apply TCGA-trained model to Korean tumors (using Korean within-cohort z)
X_kor_tumor_full = z_full.loc[mask_t][HT13].to_numpy()
p_kor_subB = clf_full.predict_proba(X_kor_tumor_full)[:, 1]
auc_transfer_subB = float(roc_auc_score(y, p_kor_subB))
metrics["tcga_to_korean_transfer_AUC_subB"] = auc_transfer_subB
metrics["tcga_to_korean_transfer_direction_subB"] = (
    "correct" if auc_transfer_subB > 0.55 else ("flipped" if auc_transfer_subB < 0.45 else "null")
)
# also vs Hashimoto
X_kor_tumor_h = z_full.loc[mask_h][HT13].to_numpy()
p_kor_h = clf_full.predict_proba(X_kor_tumor_h)[:, 1]
auc_transfer_h = float(roc_auc_score(yh, p_kor_h))
metrics["tcga_to_korean_transfer_AUC_hashi"] = auc_transfer_h
metrics["tcga_to_korean_transfer_direction_hashi"] = (
    "correct" if auc_transfer_h > 0.55 else ("flipped" if auc_transfer_h < 0.45 else "null")
)

# correlation: TCGA-HT signature vs Hashimoto sig_score (independent module)
both = results.loc[mask_t.index].copy()
both["p_dm1_tcga_model"] = clf_full.predict_proba(z_full.loc[both.index][HT13].to_numpy())[:, 1]
keep = both[["p_dm1_tcga_model", "sig_score"]].dropna()
if len(keep) > 30:
    r, p = spearmanr(keep["p_dm1_tcga_model"], keep["sig_score"])
    metrics["spearman_pDM1tcga_vs_hashi_sigscore"] = {"rho": float(r), "p": float(p), "n": int(len(keep))}

# also store HT score vs hashi sig_score (cohort-internal sanity, expect strong)
keep2 = results[["ht_score_z", "sig_score"]].dropna()
r2, p2 = spearmanr(keep2["ht_score_z"], keep2["sig_score"])
metrics["spearman_HTscore_vs_hashi_sigscore"] = {"rho": float(r2), "p": float(p2), "n": int(len(keep2))}

# write JSON
with open(OUT / "h3_k2_metrics.json", "w") as fh:
    json.dump(metrics, fh, indent=2, default=str)

print("\n=== SUMMARY ===")
for k, v in metrics.items():
    print(f"  {k}: {v}")
print(f"\nResults: {OUT}")
