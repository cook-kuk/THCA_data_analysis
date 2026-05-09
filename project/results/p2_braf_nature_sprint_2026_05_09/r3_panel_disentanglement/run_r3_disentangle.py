"""R3 — panel disentanglement for HT-13 vs Ayers-TIS-18 vs Cabrita-TLS-12.

Reviewer Q to pre-empt: "is HT-13 just a rebranded TIS / TLS?"

Analyses
--------
1) Pairwise gene-level overlap (Jaccard, intersect lists).
2) Pairwise score correlation (gene-mean z) per stratum: TCGA all-tumor,
   BRAF-cPTC, RAS/FVPTC, full master DM1 vs DM2.
3) Predictive AUCs for DM1 vs DM2 (LogReg, 5-fold stratified CV mean ROC):
   - HT-13 full / TIS-18 / Cabrita-12
   - HT-13 minus TIS-overlap (HLA-DRA, IFNG removed)
   - HT-13 minus TLS-overlap (CD79B, MS4A1, CXCL13 removed)
   - HT-13 minus all overlap (HT-only unique = the 8 truly distinct genes)
   - Combined union (unique genes)
   - Residualized HT-13 score (after regressing out TIS); TIS residualized on HT-13.
4) B-cell sub-panel (CD79A/B + MS4A1 + AICDA, n=4) vs non-B HT
   (HLA-DR/DP/DQ + CXCL13 + CCR6 + IFNG, n=9) AUC.
5) Cross-cohort transfer to GSE213647 (Lee 2024 Korean cohort) using DM1 sub-B
   call as DM1-proxy: HT-13 vs TIS vs Cabrita Spearman with subB_score and AUC
   for subB_GMM=1 vs hashi_GMM controls.

Output files in this folder.
"""

from __future__ import annotations
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

OUT = Path(__file__).parent
TCGA_Z = "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv"
MASTER = "/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv"
GSE_COUNTS = "/data/thca/v17_korean/GSE213647/expression_counts.tsv.gz"
GSE_META = "/data/thca/v17_korean/GSE213647/sample_sheet_clinical.tsv"
GSE_HASHI = "/data/thca/repo_results/d8b_korean_replication/korean_GSE213647_hashimoto_scores.tsv"
GSE_SUBB = "/data/thca/repo_results/d8c_dm1_subB_x_K2_NBNR/korean_subB_score.tsv"
GTF = "/data/thca/reference_kallisto/gencode.v44.basic.annotation.gtf"

# ---- panels ----------------------------------------------------------------
HT13 = [
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "CD79A", "CD79B", "MS4A1", "AICDA", "CXCL13", "CCR6", "IFNG",
]
TIS18 = [
    "IFNG", "CXCL9", "CD8A", "GZMA", "GZMK", "HLA-DRA", "NKG7", "PSMB10",
    "IDO1", "STAT1", "CCL5", "TIGIT", "LAG3", "PDCD1LG2", "CD274", "CMKLR1",
    "CD276", "CXCR6", "HLA-DOB", "HLA-E",
]
TLS12 = [
    "CCL19", "CCL21", "CXCL13", "CCR7", "CXCR5", "SELL",
    "LAMP3", "CD79B", "CR2", "FCRL2", "MS4A1", "JCHAIN",
]
JCHAIN_ALIASES = ["JCHAIN", "IGJ"]  # GTF uses "JCHAIN" or "IGJ"

# Sub-panels
HT_BCELL4 = ["CD79A", "CD79B", "MS4A1", "AICDA"]
HT_NONB9 = [
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "CXCL13", "CCR6", "IFNG",
]


# ---- helpers ---------------------------------------------------------------
def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return float("nan")
    return len(a & b) / max(1, len(a | b))


def gene_mean_z(genes: list, expr: pd.DataFrame) -> pd.Series:
    """expr: gene x sample, already z-scored across samples per gene."""
    pres = [g for g in genes if g in expr.index]
    if not pres:
        return pd.Series(dtype=float)
    return expr.loc[pres].mean(axis=0)


def cv_auc(X: pd.DataFrame, y: np.ndarray, seed: int = 0) -> tuple[float, float]:
    """5-fold stratified CV pooled out-of-fold AUC, plus pooled std via subsamples."""
    if X.shape[0] < 10:
        return float("nan"), float("nan")
    X = X.copy().fillna(0.0)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof = np.zeros(len(y), dtype=float)
    fold_aucs = []
    for tr, te in skf.split(X, y):
        sc = StandardScaler().fit(X.iloc[tr])
        Xtr = sc.transform(X.iloc[tr])
        Xte = sc.transform(X.iloc[te])
        clf = LogisticRegression(max_iter=2000, C=1.0).fit(Xtr, y[tr])
        oof[te] = clf.predict_proba(Xte)[:, 1]
        if len(np.unique(y[te])) == 2:
            fold_aucs.append(roc_auc_score(y[te], oof[te]))
    pooled = roc_auc_score(y, oof) if len(np.unique(y)) == 2 else float("nan")
    return pooled, float(np.std(fold_aucs)) if fold_aucs else float("nan")


def auc_from_score(score: pd.Series, y: np.ndarray) -> float:
    s = score.values
    if len(np.unique(y)) < 2 or np.all(np.isnan(s)):
        return float("nan")
    mask = ~np.isnan(s)
    if mask.sum() < 5:
        return float("nan")
    return roc_auc_score(y[mask], s[mask])


def residualize(target: pd.Series, regressor: pd.Series) -> pd.Series:
    df = pd.concat({"y": target, "x": regressor}, axis=1).dropna()
    if len(df) < 5:
        return target * np.nan
    lr = LinearRegression().fit(df[["x"]].values, df["y"].values)
    pred = lr.predict(df[["x"]].values)
    res = df["y"].values - pred
    out = pd.Series(np.nan, index=target.index)
    out.loc[df.index] = res
    return out


# ---- load TCGA ----
print("[1] Loading TCGA z-score expression ...")
expr = pd.read_csv(TCGA_Z, sep="\t", index_col=0)
print(f"   shape: {expr.shape}")
master = pd.read_csv(MASTER, sep="\t")
master = master[master["normal_vs_tumor"] == "tumor"].copy()
master = master.set_index("sample_id")

# membership presence audit
all_panel_genes = set(HT13) | set(TIS18) | set(TLS12)
present = {g: g in expr.index for g in all_panel_genes}
# Cabrita JCHAIN may appear as IGJ
if "JCHAIN" not in expr.index:
    if "IGJ" in expr.index:
        # rename
        expr.loc["JCHAIN"] = expr.loc["IGJ"]
        present["JCHAIN"] = True
print("   panel-gene presence:")
for g, ok in sorted(present.items()):
    if not ok:
        print(f"     MISSING {g}")

# ---- 1) Overlap matrix ----
print("\n[2] Pairwise overlap ...")
sets = {"HT13": set(HT13), "TIS18": set(TIS18), "TLS12": set(TLS12)}
rows = []
for a in sets:
    for b in sets:
        ov = sets[a] & sets[b]
        rows.append({
            "panel_a": a, "panel_b": b,
            "n_a": len(sets[a]), "n_b": len(sets[b]),
            "n_overlap": len(ov),
            "jaccard": jaccard(sets[a], sets[b]),
            "overlap_genes": ",".join(sorted(ov)),
        })
ov_df = pd.DataFrame(rows)
ov_df.to_csv(OUT / "r3_overlap_jaccard.tsv", sep="\t", index=False)
print(ov_df[["panel_a", "panel_b", "n_overlap", "jaccard", "overlap_genes"]]
      .query("panel_a != panel_b").to_string(index=False))


# ---- 2) Stratum scores ----
def stratum_index(df, label):
    if label == "all_tumor":
        return df.index
    if label == "BRAF_cPTC":
        return df[(df["driver_anchor"] == "BRAF") & (df["histology_subtype"] == "cPTC")].index
    if label == "RAS_FVPTC":
        return df[(df["driver_anchor"] == "RAS") & (df["histology_subtype"] == "FVPTC")].index
    if label == "DM1_DM2_only":
        return df[df["dm"].isin(["DM1", "DM2"])].index
    raise KeyError(label)


strata = ["all_tumor", "BRAF_cPTC", "RAS_FVPTC", "DM1_DM2_only"]

print("\n[3] Per-stratum score correlations ...")
corr_rows = []
score_cache: dict[tuple[str, str], pd.Series] = {}
for stratum in strata:
    sids = master.index.intersection(expr.columns)
    sids = sids.intersection(stratum_index(master, stratum))
    sub = expr[sids]
    s_ht = gene_mean_z(HT13, sub)
    s_tis = gene_mean_z(TIS18, sub)
    s_tls = gene_mean_z(TLS12, sub)
    score_cache[(stratum, "HT13")] = s_ht
    score_cache[(stratum, "TIS18")] = s_tis
    score_cache[(stratum, "TLS12")] = s_tls
    for a, b in [("HT13", "TIS18"), ("HT13", "TLS12"), ("TIS18", "TLS12")]:
        sa = score_cache[(stratum, a)]
        sb = score_cache[(stratum, b)]
        m = pd.concat({"a": sa, "b": sb}, axis=1).dropna()
        if len(m) < 5:
            r_p = np.nan; r_s = np.nan
        else:
            r_p = stats.pearsonr(m["a"], m["b"])[0]
            r_s = stats.spearmanr(m["a"], m["b"])[0]
        corr_rows.append({
            "stratum": stratum, "n": len(sids),
            "pair": f"{a}__{b}",
            "pearson_r": r_p, "spearman_r": r_s,
        })

cor_df = pd.DataFrame(corr_rows)
cor_df.to_csv(OUT / "r3_score_correlations.tsv", sep="\t", index=False)
print(cor_df.to_string(index=False))

# ---- 3) AUCs (residualized + multiple panels) ----
print("\n[4] AUCs DM1 vs DM2 ...")

def panel_X(genes, expr, sids):
    pres = [g for g in genes if g in expr.index]
    return expr.loc[pres, sids].T

auc_rows = []

for stratum in strata:
    sids = master.index.intersection(expr.columns)
    sids = sids.intersection(stratum_index(master, stratum))
    # Restrict to DM1/DM2 (remove not_DM)
    keep = master.loc[sids][master.loc[sids, "dm"].isin(["DM1", "DM2"])].index
    if len(keep) < 10:
        continue
    y = (master.loc[keep, "dm"].values == "DM1").astype(int)
    n_dm1 = int(y.sum()); n_dm2 = int(len(y) - n_dm1)
    if n_dm1 < 3 or n_dm2 < 3:
        continue

    # 1) FULL panels — multivariate LogReg on raw genes
    for name, genes in [("HT13_full", HT13), ("TIS18_full", TIS18),
                        ("TLS12_full", TLS12)]:
        X = panel_X(genes, expr, keep)
        auc, sd = cv_auc(X, y)
        auc_rows.append({"stratum": stratum, "panel": name,
                          "n": len(keep), "n_dm1": n_dm1, "n_dm2": n_dm2,
                          "n_genes": X.shape[1], "auc_cv": auc, "auc_sd": sd})

    # 2) HT-13 minus overlapping genes
    ht_minus_tis = [g for g in HT13 if g not in set(TIS18)]
    ht_minus_tls = [g for g in HT13 if g not in set(TLS12)]
    ht_minus_all = [g for g in HT13 if g not in (set(TIS18) | set(TLS12))]
    union_unique = sorted(set(HT13) | set(TIS18) | set(TLS12))
    for name, genes in [("HT13_minus_TIS_overlap", ht_minus_tis),
                        ("HT13_minus_TLS_overlap", ht_minus_tls),
                        ("HT13_unique_only", ht_minus_all),
                        ("union_HT_TIS_TLS", union_unique)]:
        X = panel_X(genes, expr, keep)
        auc, sd = cv_auc(X, y)
        auc_rows.append({"stratum": stratum, "panel": name,
                          "n": len(keep), "n_dm1": n_dm1, "n_dm2": n_dm2,
                          "n_genes": X.shape[1], "auc_cv": auc, "auc_sd": sd})

    # 3) Score-based AUC (single mean-z per panel) — for residualization
    s_ht = gene_mean_z(HT13, expr[keep])
    s_tis = gene_mean_z(TIS18, expr[keep])
    s_tls = gene_mean_z(TLS12, expr[keep])
    auc_rows.append({"stratum": stratum, "panel": "HT13_score",
                      "n": len(keep), "n_dm1": n_dm1, "n_dm2": n_dm2,
                      "n_genes": len(HT13),
                      "auc_cv": auc_from_score(s_ht, y), "auc_sd": np.nan})
    auc_rows.append({"stratum": stratum, "panel": "TIS18_score",
                      "n": len(keep), "n_dm1": n_dm1, "n_dm2": n_dm2,
                      "n_genes": len(TIS18),
                      "auc_cv": auc_from_score(s_tis, y), "auc_sd": np.nan})
    auc_rows.append({"stratum": stratum, "panel": "TLS12_score",
                      "n": len(keep), "n_dm1": n_dm1, "n_dm2": n_dm2,
                      "n_genes": len(TLS12),
                      "auc_cv": auc_from_score(s_tls, y), "auc_sd": np.nan})

    # Residualized scores (the killer test)
    for tname, target_genes, tag_t in [("HT13", HT13, "HT13_score"),
                                       ("TIS18", TIS18, "TIS18_score"),
                                       ("TLS12", TLS12, "TLS12_score")]:
        s_t = gene_mean_z(target_genes, expr[keep])
        for rname, reg_genes in [("HT13", HT13), ("TIS18", TIS18),
                                 ("TLS12", TLS12)]:
            if rname == tname:
                continue
            s_r = gene_mean_z(reg_genes, expr[keep])
            res = residualize(s_t, s_r)
            auc_rows.append({"stratum": stratum,
                              "panel": f"{tname}_residualized_on_{rname}",
                              "n": len(keep), "n_dm1": n_dm1, "n_dm2": n_dm2,
                              "n_genes": np.nan,
                              "auc_cv": auc_from_score(res, y),
                              "auc_sd": np.nan})

auc_df = pd.DataFrame(auc_rows)
auc_df.to_csv(OUT / "r3_residualized_aucs.tsv", sep="\t", index=False)
print(auc_df.query("stratum=='DM1_DM2_only'").to_string(index=False))


# ---- 4) Sub-panel AUCs ----
print("\n[5] Sub-panel AUCs ...")
sub_rows = []
for stratum in strata:
    sids = master.index.intersection(expr.columns)
    sids = sids.intersection(stratum_index(master, stratum))
    keep = master.loc[sids][master.loc[sids, "dm"].isin(["DM1", "DM2"])].index
    if len(keep) < 10:
        continue
    y = (master.loc[keep, "dm"].values == "DM1").astype(int)
    n_dm1 = int(y.sum()); n_dm2 = int(len(y) - n_dm1)
    if n_dm1 < 3 or n_dm2 < 3:
        continue
    for name, genes in [("HT_Bcell4", HT_BCELL4),
                        ("HT_nonB9", HT_NONB9),
                        ("HT13_full", HT13)]:
        X = panel_X(genes, expr, keep)
        auc, sd = cv_auc(X, y)
        # Per-gene mean Cohen's d helper
        s = gene_mean_z(genes, expr[keep])
        a = s[y == 1].dropna(); b = s[y == 0].dropna()
        if len(a) > 1 and len(b) > 1:
            sd_p = np.sqrt(((len(a) - 1) * a.var() + (len(b) - 1) * b.var())
                            / (len(a) + len(b) - 2))
            d = (a.mean() - b.mean()) / sd_p if sd_p > 0 else np.nan
        else:
            d = np.nan
        sub_rows.append({"stratum": stratum, "panel": name,
                          "n": len(keep), "n_dm1": n_dm1, "n_dm2": n_dm2,
                          "n_genes": X.shape[1],
                          "auc_cv": auc, "auc_sd": sd,
                          "score_cohens_d": d})
sub_df = pd.DataFrame(sub_rows)
sub_df.to_csv(OUT / "r3_subpanel_aucs.tsv", sep="\t", index=False)
print(sub_df.to_string(index=False))


# ---- 5) Mutual information (decomposition) ----
print("\n[6] Mutual information (binned) ...")

def mi_binned(score: pd.Series, y: np.ndarray, bins: int = 5) -> float:
    """Estimate MI via equal-frequency binning."""
    s = score.dropna()
    yy = y[score.notna().values]
    if len(s) < 20 or len(np.unique(yy)) < 2:
        return float("nan")
    qs = np.quantile(s, np.linspace(0, 1, bins + 1))
    qs[0] -= 1e-9; qs[-1] += 1e-9
    bin_idx = np.digitize(s, qs[1:-1])
    # joint
    joint = pd.crosstab(bin_idx, yy).values.astype(float)
    n = joint.sum()
    if n == 0:
        return float("nan")
    p_joint = joint / n
    p_x = p_joint.sum(axis=1, keepdims=True)
    p_y = p_joint.sum(axis=0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(p_joint > 0, p_joint / (p_x * p_y), 1.0)
        mi = (p_joint * np.log(np.where(ratio > 0, ratio, 1.0))).sum()
    return float(mi)


def cmi_binned(target: pd.Series, y: np.ndarray, cond: pd.Series, bins: int = 4) -> float:
    """I(target ; y | cond) via binning of cond, weighted MI per bin."""
    df = pd.concat({"t": target, "c": cond}, axis=1).dropna()
    if df.empty:
        return float("nan")
    yy = y[target.notna().values & cond.notna().values]
    if len(df) != len(yy):
        # align
        yy = pd.Series(y, index=target.index).loc[df.index].values
    qs = np.quantile(df["c"], np.linspace(0, 1, bins + 1))
    qs[0] -= 1e-9; qs[-1] += 1e-9
    cb = np.digitize(df["c"], qs[1:-1])
    out = 0.0
    n = len(df)
    for b in np.unique(cb):
        mask = cb == b
        if mask.sum() < 10 or len(np.unique(yy[mask])) < 2:
            continue
        out += (mask.sum() / n) * mi_binned(df["t"][mask].reset_index(drop=True),
                                            yy[mask])
    return out


mi_rows = []
for stratum in ["DM1_DM2_only", "BRAF_cPTC"]:
    sids = master.index.intersection(expr.columns)
    sids = sids.intersection(stratum_index(master, stratum))
    keep = master.loc[sids][master.loc[sids, "dm"].isin(["DM1", "DM2"])].index
    if len(keep) < 20:
        continue
    y = (master.loc[keep, "dm"].values == "DM1").astype(int)
    s_ht = gene_mean_z(HT13, expr[keep])
    s_tis = gene_mean_z(TIS18, expr[keep])
    s_tls = gene_mean_z(TLS12, expr[keep])
    mi_ht = mi_binned(s_ht, y)
    mi_tis = mi_binned(s_tis, y)
    mi_tls = mi_binned(s_tls, y)
    cmi_ht_g_tis = cmi_binned(s_ht, y, s_tis)
    cmi_tis_g_ht = cmi_binned(s_tis, y, s_ht)
    cmi_ht_g_tls = cmi_binned(s_ht, y, s_tls)
    cmi_tls_g_ht = cmi_binned(s_tls, y, s_ht)
    mi_rows.append({
        "stratum": stratum, "n": len(keep),
        "MI_HT13": mi_ht, "MI_TIS18": mi_tis, "MI_TLS12": mi_tls,
        "CMI_HT13_given_TIS": cmi_ht_g_tis,
        "CMI_TIS_given_HT13": cmi_tis_g_ht,
        "CMI_HT13_given_TLS": cmi_ht_g_tls,
        "CMI_TLS_given_HT13": cmi_tls_g_ht,
    })

mi_df = pd.DataFrame(mi_rows)
mi_df.to_csv(OUT / "r3_mutual_information.tsv", sep="\t", index=False)
print(mi_df.to_string(index=False))


# ---- 6) Cross-cohort transfer to GSE213647 (Lee 2024 Korean) ----
print("\n[7] Cross-cohort transfer GSE213647 (Lee 2024) ...")

# Load all panel genes via ENSG mapping
import re
ALL_NEEDED = set(HT13) | set(TIS18) | set(TLS12) | {"IGJ"}
# CCR6 alias
ensg2sym: dict[str, str] = {}
needed = ALL_NEEDED
with open(GTF) as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        if "\tgene\t" not in ln:
            continue
        m = re.search(r'gene_id "([^"]+)"', ln)
        s = re.search(r'gene_name "([^"]+)"', ln)
        if not m or not s:
            continue
        if s.group(1) in needed:
            ensg2sym[m.group(1)] = s.group(1)

sym2ensg = {v: k for k, v in ensg2sym.items()}
sym2ensg_nv = {v: k.split(".")[0] for k, v in ensg2sym.items()}
print(f"   resolved {len(sym2ensg)} symbols from GTF")

print("   loading GSE213647 counts ...")
counts = pd.read_csv(GSE_COUNTS, sep="\t", index_col=0)
counts.index = counts.index.astype(str)
counts_index_nv = pd.Series(counts.index, index=[i.split(".")[0] for i in counts.index])

found_rows = {}
missing = []
for g in ALL_NEEDED:
    ensg_nv = sym2ensg_nv.get(g)
    if ensg_nv is None or ensg_nv not in counts_index_nv.index:
        missing.append(g)
        continue
    full = counts_index_nv.loc[ensg_nv]
    if isinstance(full, pd.Series):
        full = full.iloc[0]
    found_rows[g] = counts.loc[full]
print(f"   missing in GSE213647: {sorted(missing)}")

# Map JCHAIN <- IGJ if needed
if "JCHAIN" not in found_rows and "IGJ" in found_rows:
    found_rows["JCHAIN"] = found_rows["IGJ"]

panel_df = pd.DataFrame(found_rows)
# log2(cpm+1)
lib = counts.sum(axis=0)
cpm = panel_df.div(lib, axis=0) * 1e6
log_cpm = np.log2(cpm + 1.0)
# z-score per gene across all GSE213647 samples
gz = (log_cpm - log_cpm.mean(axis=0)) / log_cpm.std(axis=0)
print(f"   GSE213647 panel z shape: {gz.shape}")

# load metadata + subB labels
meta = pd.read_csv(GSE_META, sep="\t").set_index("gsm")
subb = pd.read_csv(GSE_SUBB, sep="\t", index_col=0)
subb.index.name = "gsm"
hashi = pd.read_csv(GSE_HASHI, sep="\t", index_col=0)
hashi.index.name = "gsm"
common = gz.index.intersection(meta.index).intersection(subb.index)
gz = gz.loc[common]
meta = meta.loc[common]
subb = subb.loc[common]
hashi = hashi.loc[common]

# Build scores per panel from log_cpm z
def gz_panel_score(genes):
    pres = [g for g in genes if g in gz.columns]
    return gz[pres].mean(axis=1) if pres else pd.Series(np.nan, index=gz.index)

s_ht = gz_panel_score(HT13)
s_tis = gz_panel_score(TIS18)
s_tls = gz_panel_score(TLS12)

# Tumor only
tumor_mask = meta["tissue_type"] != "Normal"
# DM1-proxy: subB_GMM == 1 (sub-B = TCGA DM1 dedifferentiated module)
subb_call = subb["subB_GMM"].astype(int)

x_rows = []
# 1) tumor-only AUC for subB=1 vs subB=0
m_t = tumor_mask
for sname, s in [("HT13", s_ht), ("TIS18", s_tis), ("TLS12", s_tls)]:
    df = pd.concat({"s": s, "y": subb_call, "tumor": tumor_mask}, axis=1).dropna()
    df = df[df["tumor"]]
    if df["y"].nunique() == 2 and len(df) >= 10:
        auc = roc_auc_score(df["y"], df["s"])
    else:
        auc = float("nan")
    x_rows.append({"cohort": "GSE213647_Lee2024",
                    "label": "subB_proxy_DM1_tumor_only",
                    "panel": sname, "n": int(df.shape[0]),
                    "auc": auc})

# 2) tumor vs normal — sanity check (panels should be ELEVATED in tumor, but
#    these immune panels reflect immune infiltrate, not disease-class)
for sname, s in [("HT13", s_ht), ("TIS18", s_tis), ("TLS12", s_tls)]:
    yy = (meta["tissue_type"] != "Normal").astype(int)
    df = pd.concat({"s": s, "y": yy}, axis=1).dropna()
    if df["y"].nunique() == 2 and len(df) >= 10:
        auc = roc_auc_score(df["y"], df["s"])
    else:
        auc = float("nan")
    x_rows.append({"cohort": "GSE213647_Lee2024", "label": "tumor_vs_normal",
                    "panel": sname, "n": int(df.shape[0]), "auc": auc})

# 3) hashimoto (HT-overlap) call
for sname, s in [("HT13", s_ht), ("TIS18", s_tis), ("TLS12", s_tls)]:
    df = pd.concat({"s": s, "y": hashi["hashi_GMM"].astype(int),
                    "tumor": tumor_mask}, axis=1).dropna()
    df = df[df["tumor"]]
    if df["y"].nunique() == 2 and len(df) >= 10:
        auc = roc_auc_score(df["y"], df["s"])
    else:
        auc = float("nan")
    x_rows.append({"cohort": "GSE213647_Lee2024",
                    "label": "hashi_GMM_tumor_only",
                    "panel": sname, "n": int(df.shape[0]), "auc": auc})

# 4) Spearman correlations between panel scores in Korean cohort
m_tumor = tumor_mask & s_ht.notna() & s_tis.notna() & s_tls.notna()
n_t = int(m_tumor.sum())
for a, sa in [("HT13", s_ht), ("TIS18", s_tis), ("TLS12", s_tls)]:
    for b, sb in [("HT13", s_ht), ("TIS18", s_tis), ("TLS12", s_tls)]:
        if a >= b:
            continue
        df = pd.concat({"a": sa, "b": sb}, axis=1)[tumor_mask].dropna()
        if len(df) < 10:
            continue
        r = stats.spearmanr(df["a"], df["b"])[0]
        x_rows.append({"cohort": "GSE213647_Lee2024",
                        "label": f"spearman_tumor_{a}_vs_{b}",
                        "panel": f"{a}__{b}", "n": len(df), "auc": r})

cross_df = pd.DataFrame(x_rows)
cross_df.to_csv(OUT / "r3_crosscohort_aucs.tsv", sep="\t", index=False)
print(cross_df.to_string(index=False))


# ---- summary JSON ----
summary = {
    "n_HT13": len(HT13),
    "n_TIS18": len(TIS18),
    "n_TLS12": len(TLS12),
    "overlap_HT13_TIS18": sorted(set(HT13) & set(TIS18)),
    "overlap_HT13_TLS12": sorted(set(HT13) & set(TLS12)),
    "overlap_TIS18_TLS12": sorted(set(TIS18) & set(TLS12)),
    "HT13_unique_only": sorted(set(HT13) - set(TIS18) - set(TLS12)),
    "jaccard_HT13_TIS18": jaccard(set(HT13), set(TIS18)),
    "jaccard_HT13_TLS12": jaccard(set(HT13), set(TLS12)),
    "jaccard_TIS18_TLS12": jaccard(set(TIS18), set(TLS12)),
}
with open(OUT / "r3_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2)

print("\n[done] Outputs in", OUT)
