"""H30: Comprehensive Lee 2024 (GSE213647 n=632) external replication.

Extends H3 (HT-13 within-cohort AUC=0.948) to full FA-12 + MAPK-9 panels,
HT x FA two-axis combo, DM1/DM2 surrogate via GMM, and driver / survival /
TERT replication where possible.

Driver / survival / TERT context:
  - GSE213647 SOFT (genotype: NA, treatment: NA for all 632) and
    Lee 2024 supplementary MOESM5 (only sample-id, GEO id, fixation,
    library kit, dx_year, age, sex; no driver, no TERT, no follow-up)
    do NOT carry per-sample driver / mutation / survival annotation.
  - Therefore H6 (PFI HR=5.91 / 0.21) and H8 (DM1xTERT) are NOT TESTABLE
    in Lee 2024 alone. Reported as "not testable, annotation absent in
    public data" (this is itself an honest finding for the manuscript).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture

OUT = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_braf_nature_sprint_2026_05_09/h30_lee2024_full"
)
OUT.mkdir(parents=True, exist_ok=True)

HT13 = [
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "CD79A", "CD79B", "MS4A1", "AICDA", "CXCL13", "CCR6", "IFNG",
]
FA12 = [
    "FASN", "ACACA", "ACLY", "SCD", "FADS1", "FADS2", "ELOVL6",
    "ACOX1", "CPT1A", "HMGCS2", "HADH", "ACADM",
]
MAPK9 = [
    "DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4",
    "ETV4", "ETV5", "PHLDA1", "CCND1",
]
PANEL_8GENE = ["TG", "TPO", "SLC5A5", "PAX8", "FOXE1", "DIO1", "TSHR", "NKX2-1"]

ALL_GENES = sorted(set(HT13 + FA12 + MAPK9 + PANEL_8GENE))

# ---------------------------------------------------------------------------
print("[1/8] Loading GSE213647 counts ...")
counts = pd.read_csv(
    "/data/thca/v17_korean/GSE213647/expression_counts.tsv.gz",
    sep="\t",
    index_col=0,
)
print(f"   counts shape: {counts.shape}")

print("[2/8] Building ENSG -> symbol map from gencode v44 ...")
gtf_path = "/data/thca/reference_kallisto/gencode.v44.basic.annotation.gtf"
ensg2sym: dict[str, str] = {}
needed = set(ALL_GENES)
with open(gtf_path) as fh:
    for ln in fh:
        if ln.startswith("#") or "\tgene\t" not in ln:
            continue
        m = re.search(r'gene_id "([^"]+)"', ln)
        s = re.search(r'gene_name "([^"]+)"', ln)
        if not m or not s:
            continue
        if s.group(1) in needed:
            ensg2sym[m.group(1)] = s.group(1)

sym2ensg_nv = {sym: eid.split(".")[0] for eid, sym in ensg2sym.items()}
counts_idx_nv = pd.Series(
    counts.index.values, index=[i.split(".")[0] for i in counts.index]
)


def fetch_panel(genes: list[str]) -> tuple[pd.DataFrame, list[str]]:
    rows: list[tuple[str, pd.Series]] = []
    missing: list[str] = []
    for g in genes:
        eid_nv = sym2ensg_nv.get(g)
        if eid_nv is None or eid_nv not in counts_idx_nv.index:
            missing.append(g)
            continue
        full = counts_idx_nv.loc[eid_nv]
        if isinstance(full, pd.Series):
            full = full.iloc[0]
        rows.append((g, counts.loc[full]))
    df = pd.DataFrame({g: r for g, r in rows})
    return df, missing


# Library size for log2(CPM+1)
lib = counts.sum(axis=0)


def to_log_cpm(df_counts: pd.DataFrame) -> pd.DataFrame:
    cpm = (df_counts.T / lib * 1e6).T
    return np.log2(cpm + 1.0)


print("[3/8] Building HT-13 / FA-12 / MAPK-9 / 8-gene panels ...")
ht_raw, miss_ht = fetch_panel(HT13)
fa_raw, miss_fa = fetch_panel(FA12)
mapk_raw, miss_mapk = fetch_panel(MAPK9)
panel8_raw, miss_p8 = fetch_panel(PANEL_8GENE)

ht_lc = to_log_cpm(ht_raw)
fa_lc = to_log_cpm(fa_raw)
mapk_lc = to_log_cpm(mapk_raw)
p8_lc = to_log_cpm(panel8_raw)
print(
    f"   HT-13: {ht_lc.shape[1]}/13 (missing {miss_ht})\n"
    f"   FA-12: {fa_lc.shape[1]}/12 (missing {miss_fa})\n"
    f"   MAPK-9: {mapk_lc.shape[1]}/9 (missing {miss_mapk})\n"
    f"   8-gene: {p8_lc.shape[1]}/8 (missing {miss_p8})"
)


# ---------------------------------------------------------------------------
print("[4/8] Within-cohort z-scoring + panel scores ...")


def panel_score(lc: pd.DataFrame) -> pd.Series:
    z = (lc - lc.mean()) / lc.std(ddof=0)
    return z.mean(axis=1)


ht_score = panel_score(ht_lc)
fa_score = panel_score(fa_lc)
mapk_score = panel_score(mapk_lc)
p8_score = panel_score(p8_lc)

# Attach metadata + prior calls
meta = pd.read_csv(
    "/data/thca/v17_korean/GSE213647/sample_sheet_clinical.tsv", sep="\t"
).set_index("gsm")
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

idx = (
    ht_score.index
    .intersection(fa_score.index)
    .intersection(mapk_score.index)
    .intersection(p8_score.index)
    .intersection(meta.index)
)

results = pd.DataFrame(
    {
        "ht_score": ht_score.loc[idx],
        "fa_score": fa_score.loc[idx],
        "mapk_score": mapk_score.loc[idx],
        "p8_score": p8_score.loc[idx],
        "tissue_type": meta.loc[idx, "tissue_type"],
        "histology": meta.loc[idx, "histology"],
        "age": meta.loc[idx, "age"],
        "sex": meta.loc[idx, "sex"],
    }
)
results = results.join(sub[["subB_score", "subB_GMM"]], how="left")
results = results.join(hashi[["hashi_GMM", "hashi_otsu", "sig_score"]], how="left")

# DM1/DM2 surrogate within Lee 2024:
# Definition: DM1 = high p8_score (TG/TPO/etc preserved differentiated module).
# Use 2-component GMM on tumor-only p8 distribution; "DM1" = upper component.
# Note: TCGA convention is DM1 = high differentiation panel (preserved 8-gene),
# DM2 = low. This matches r9_1 master labels.
tumor_mask = results["tissue_type"].isin(["PTC", "PDTC", "ATC"])
x_p8 = results.loc[tumor_mask, "p8_score"].dropna().to_numpy().reshape(-1, 1)
gmm = GaussianMixture(n_components=2, random_state=0).fit(x_p8)
means = gmm.means_.ravel()
hi_comp = int(np.argmax(means))
labels = gmm.predict(x_p8)
dm_call = pd.Series(
    np.where(labels == hi_comp, "DM1", "DM2"),
    index=results.loc[tumor_mask, "p8_score"].dropna().index,
)
results["dm_surrogate"] = dm_call
results.loc[~tumor_mask, "dm_surrogate"] = "Normal"

# Bimodality for Hashimoto cohort decision
bic1 = GaussianMixture(n_components=1, random_state=0).fit(x_p8).bic(x_p8)
bic2 = gmm.bic(x_p8)

results.to_csv(OUT / "h30_lee2024_panel_scores.tsv", sep="\t")
print(
    f"   tumor n={int(tumor_mask.sum())}  "
    f"DM1 surrogate n={int((dm_call=='DM1').sum())}  "
    f"DM2 n={int((dm_call=='DM2').sum())}  "
    f"BIC delta(2-1)={bic2-bic1:.1f}"
)


# ---------------------------------------------------------------------------
print("[5/8] Replication: panel effect sizes and AUCs ...")


def cohen_d(a, b):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    s = np.sqrt(((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1)) / (len(a)+len(b)-2))
    return float((a.mean() - b.mean()) / s) if s > 0 else float("nan")


def safe_auc(y, s):
    y = pd.Series(y).reset_index(drop=True)
    s = pd.Series(s).reset_index(drop=True)
    keep = y.notna() & s.notna()
    if keep.sum() < 10 or y[keep].nunique() < 2:
        return float("nan"), int(keep.sum())
    return float(roc_auc_score(y[keep], s[keep])), int(keep.sum())


tumor = results[tumor_mask].copy()
dm1_mask = tumor["dm_surrogate"] == "DM1"
dm2_mask = tumor["dm_surrogate"] == "DM2"

reps: list[dict] = []
# H2 effect sizes (DM1 vs DM2 in TCGA): HT d=+1.58, MAPK d=+1.81, FA d=-0.64
for panel_name, score_col, tcga_d in [
    ("HT-13", "ht_score", 1.583),
    ("MAPK-9", "mapk_score", 1.810),
    ("FA-12", "fa_score", -0.641),
]:
    a = tumor.loc[dm1_mask, score_col].dropna()
    b = tumor.loc[dm2_mask, score_col].dropna()
    d = cohen_d(a, b)
    try:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        u, p = float("nan"), float("nan")
    reps.append(
        {
            "test": f"{panel_name}_DM1vsDM2_cohen_d",
            "tcga_target": tcga_d,
            "lee2024_value": d,
            "n_DM1": int(dm1_mask.sum()),
            "n_DM2": int(dm2_mask.sum()),
            "p": float(p) if p == p else float("nan"),
            "verdict": (
                "REPLICATES_SIGN" if (np.sign(d) == np.sign(tcga_d) and abs(d) >= 0.3)
                else "WEAK_OR_FLIP"
            ),
        }
    )

# H3 within-cohort: HT score AUC vs Hashimoto (target 0.948)
auc_ht_hashi, n_ht_hashi = safe_auc(results["hashi_GMM"], results["ht_score"])
reps.append(
    {
        "test": "AUC_HT_vs_Hashimoto_GMM",
        "tcga_target": 0.948,
        "lee2024_value": auc_ht_hashi,
        "n": n_ht_hashi,
        "verdict": "REPLICATES" if (auc_ht_hashi >= 0.85) else "WEAK",
    }
)

# H1 echo: HT AUC discriminating dediff (ATC/PDTC vs PTC)
dediff_y = (tumor["tissue_type"].isin(["ATC", "PDTC"])).astype(int)
auc_ht_dediff, n_dediff = safe_auc(dediff_y, tumor["ht_score"])
auc_fa_dediff, _ = safe_auc(dediff_y, tumor["fa_score"])
auc_mapk_dediff, _ = safe_auc(dediff_y, tumor["mapk_score"])
reps.append(
    {"test": "AUC_HT_dediff_ATCPDTCvsPTC", "tcga_target": "n/a (dediff axis)",
     "lee2024_value": auc_ht_dediff, "n": n_dediff, "verdict": "INFO"}
)
reps.append(
    {"test": "AUC_FA_dediff_ATCPDTCvsPTC", "tcga_target": "n/a (dediff axis)",
     "lee2024_value": auc_fa_dediff, "n": n_dediff, "verdict": "INFO"}
)
reps.append(
    {"test": "AUC_MAPK_dediff_ATCPDTCvsPTC", "tcga_target": "n/a (dediff axis)",
     "lee2024_value": auc_mapk_dediff, "n": n_dediff, "verdict": "INFO"}
)

# DM-surrogate AUCs (within-cohort)
y_dm1 = (tumor["dm_surrogate"] == "DM1").astype(int)
auc_ht_dm, _ = safe_auc(y_dm1, tumor["ht_score"])
auc_mapk_dm, _ = safe_auc(y_dm1, tumor["mapk_score"])
auc_fa_dm, _ = safe_auc(y_dm1, -tumor["fa_score"])  # DM1 = low FA; flip
auc_p8_dm, _ = safe_auc(y_dm1, tumor["p8_score"])  # tautology, sanity
reps.append({"test": "AUC_HT_vs_DMsurrogate_DM1", "tcga_target": "n/a",
             "lee2024_value": auc_ht_dm, "n": int(tumor_mask.sum()), "verdict": "INFO"})
reps.append({"test": "AUC_MAPK_vs_DMsurrogate_DM1", "tcga_target": "n/a",
             "lee2024_value": auc_mapk_dm, "n": int(tumor_mask.sum()), "verdict": "INFO"})
reps.append({"test": "AUC_FAneg_vs_DMsurrogate_DM1", "tcga_target": "n/a",
             "lee2024_value": auc_fa_dm, "n": int(tumor_mask.sum()), "verdict": "INFO"})
reps.append({"test": "AUC_p8_vs_DMsurrogate_DM1_sanity", "tcga_target": ">=0.95 (tautology)",
             "lee2024_value": auc_p8_dm, "n": int(tumor_mask.sum()), "verdict": "SANITY"})

# Two-axis combo: HT + (-FA) z-score combo
combo = ((tumor["ht_score"] - tumor["ht_score"].mean()) / tumor["ht_score"].std(ddof=0)
         + (-(tumor["fa_score"] - tumor["fa_score"].mean()) / tumor["fa_score"].std(ddof=0)))
combo = combo / 2.0
auc_combo_dm, _ = safe_auc(y_dm1, combo)
auc_combo_hashi_t, n_ch = safe_auc(tumor["hashi_GMM"], combo)
reps.append(
    {
        "test": "AUC_HTplusFA_combo_vs_DMsurrogate",
        "tcga_target": 0.974,  # TCGA BRAF-cPTC target
        "lee2024_value": auc_combo_dm,
        "n": int(tumor_mask.sum()),
        "verdict": "REPLICATES" if auc_combo_dm >= 0.85 else "WEAK",
    }
)
reps.append(
    {
        "test": "AUC_HTplusFA_combo_vs_Hashimoto",
        "tcga_target": "n/a",
        "lee2024_value": auc_combo_hashi_t,
        "n": n_ch,
        "verdict": "INFO",
    }
)

# Pearson correlation HT x FA on tumors (TCGA target r = -0.16)
x = tumor["ht_score"].dropna()
y = tumor["fa_score"].dropna()
shared = x.index.intersection(y.index)
r_ht_fa, p_ht_fa = pearsonr(x.loc[shared], y.loc[shared])
r_ht_mapk, _ = pearsonr(tumor.loc[shared, "ht_score"], tumor.loc[shared, "mapk_score"])
r_fa_mapk, _ = pearsonr(tumor.loc[shared, "fa_score"], tumor.loc[shared, "mapk_score"])
reps.append(
    {
        "test": "pearson_HT_x_FA_tumor",
        "tcga_target": -0.156,
        "lee2024_value": float(r_ht_fa),
        "n": len(shared),
        "p": float(p_ht_fa),
        "verdict": (
            "REPLICATES_ORTHOGONAL" if abs(r_ht_fa) <= 0.30 else "AXES_NOT_ORTHOGONAL"
        ),
    }
)
reps.append(
    {"test": "pearson_HT_x_MAPK_tumor", "tcga_target": 0.233,
     "lee2024_value": float(r_ht_mapk), "n": len(shared), "verdict": "INFO"}
)
reps.append(
    {"test": "pearson_FA_x_MAPK_tumor", "tcga_target": 0.007,
     "lee2024_value": float(r_fa_mapk), "n": len(shared), "verdict": "INFO"}
)

# H6 BRAF-cPTC PFI HR=0.21 / DSS HR=8.78 / H8 DM1xTERT — NOT TESTABLE
for label, target in [
    ("PFI_HR_DM1vsDM2_BRAF_cPTC", 0.211),
    ("OS_HR_DM1vsDM2_BRAF_cPTC", 0.067),
    ("DSS_HR_DM1vsDM2_BRAF_cPTC", 8.78),
    ("DM1xTERT_PFI_HR", 6.94),
]:
    reps.append(
        {
            "test": label,
            "tcga_target": target,
            "lee2024_value": float("nan"),
            "n": 0,
            "verdict": (
                "NOT_TESTABLE_no_survival_in_GSE213647 "
                "(GEO SOFT genotype/treatment=NA; Lee 2024 MOESM5 has only "
                "fixation/library/dx_year/age/sex)"
            ),
        }
    )

rep_df = pd.DataFrame(reps)
rep_df.to_csv(OUT / "h30_lee2024_replication_table.tsv", sep="\t", index=False)


# ---------------------------------------------------------------------------
print("[6/8] Verdict counting ...")
testable = rep_df[~rep_df["verdict"].astype(str).str.startswith("NOT_TESTABLE")]
testable = testable[~testable["verdict"].isin(["INFO", "SANITY"])]
n_replicate = int(
    testable["verdict"].astype(str).str.contains(
        "REPLICATES", regex=False
    ).sum()
)
n_testable = int(len(testable))
n_blocked = int(
    rep_df["verdict"].astype(str).str.startswith("NOT_TESTABLE").sum()
)


# ---------------------------------------------------------------------------
print("[7/8] Save metrics JSON ...")
metrics = {
    "cohort": "GSE213647 Lee 2024 Nat Commun (Korean SNUBH/CNUH/KRIBB)",
    "n_total": int(len(results)),
    "n_tumor": int(tumor_mask.sum()),
    "n_normal": int((~tumor_mask).sum()),
    "panel_coverage": {
        "HT13": [ht_lc.shape[1], 13, miss_ht],
        "FA12": [fa_lc.shape[1], 12, miss_fa],
        "MAPK9": [mapk_lc.shape[1], 9, miss_mapk],
        "panel8_diff": [p8_lc.shape[1], 8, miss_p8],
    },
    "dm_surrogate": {
        "definition": "GMM(n=2) on tumor-only p8_score; DM1=upper, DM2=lower",
        "n_DM1": int(dm1_mask.sum()),
        "n_DM2": int(dm2_mask.sum()),
        "bic1": float(bic1),
        "bic2": float(bic2),
        "bic_supports_bimodal": bool(bic2 < bic1),
    },
    "panel_AUCs_DM1surrogate": {
        "HT": auc_ht_dm,
        "MAPK": auc_mapk_dm,
        "FA_neg": auc_fa_dm,
        "HTplusFA_combo": auc_combo_dm,
    },
    "axis_orthogonality_pearson": {
        "HT_x_FA": float(r_ht_fa),
        "HT_x_MAPK": float(r_ht_mapk),
        "FA_x_MAPK": float(r_fa_mapk),
    },
    "replication_summary": {
        "n_testable_findings": n_testable,
        "n_replicating": n_replicate,
        "n_NOT_testable_due_to_missing_annotation": n_blocked,
    },
    "annotation_gaps": {
        "BRAF_RAS_status": "NOT in GEO SOFT (genotype: NA), NOT in MOESM5",
        "TERT_promoter": "NOT in GEO SOFT, NOT in MOESM5",
        "follow_up_OS_PFI_DSS": "NOT in GEO SOFT, NOT in MOESM5",
        "implication": (
            "Driver-stratified and survival replications must come from a "
            "different external cohort (cBioPortal thca_tcga_pub already used "
            "in H6/H8; alternative: Korean follow-up audit via SNUBH IRB)"
        ),
    },
}
(OUT / "h30_metrics.json").write_text(json.dumps(metrics, indent=2))


# ---------------------------------------------------------------------------
print("[8/8] Done.")
print(json.dumps(metrics["replication_summary"], indent=2))
print(json.dumps(metrics["panel_AUCs_DM1surrogate"], indent=2))
print(json.dumps(metrics["axis_orthogonality_pearson"], indent=2))
