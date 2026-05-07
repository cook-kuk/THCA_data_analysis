#!/usr/bin/env python3
"""Paper 3 ICI Track B-lite — Phase 3 wrap-up.

(1) Singscore-style alternative scoring (rank-based, sample-independent) on
    GSE126698 + Hugo + Riaz_pre + Riaz_on → check ATC direction + DIAL-lite
    direction match the primary ssGSEA-like (mean) results.
(2) Riaz Pre→On delta split by responder (CR/PR) vs non-responder (PD).
(3) Combine 6 PNG figures into single 2x3 montage.
"""
from __future__ import annotations
import gzip, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

DATA_DIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_public_data")
REG_DIR  = Path("/home/seungho/personal/THCA_data_analysis/project/reports/paper3_ici")
OUT_DIR  = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_track_b_lite")
FIG_DIR  = OUT_DIR / "figures_png"

mods_df = pd.read_csv(REG_DIR / "paper3_ici_module_gene_list.tsv", sep="\t")
MODULES = mods_df.groupby("module")["gene_symbol"].apply(list).to_dict()
ALL_MODULE_GENES = sorted({g for genes in MODULES.values() for g in genes})

probe_map = pd.read_csv(OUT_DIR / "gpl570_probe_symbol_map_module_subset.tsv", sep="\t")
ENTREZ2SYM = {}
for _, r in probe_map.iterrows():
    if r.entrez_id and not pd.isna(r.entrez_id):
        ENTREZ2SYM[str(int(r.entrez_id))] = r.gene_symbol

def singscore_module(expr: pd.DataFrame, module_genes: list) -> pd.Series:
    """Simplified Singscore: rank genes per sample, take mean rank of module genes,
       normalize to [-1,1]."""
    present = [g for g in module_genes if g in expr.columns]
    if not present:
        return pd.Series([np.nan]*len(expr), index=expr.index)
    ranks = expr.rank(axis=1, method="average")  # rank across genes per sample
    n_genes = expr.shape[1]
    mod_rank = ranks[present].mean(axis=1)
    # normalize: (mean_rank - mid) / (n_genes/2)
    mid = (n_genes + 1) / 2.0
    score = (mod_rank - mid) / (n_genes / 2.0)
    return score

def score_singscore_all(expr: pd.DataFrame) -> pd.DataFrame:
    out = {}
    for mod, genes in MODULES.items():
        out[mod] = singscore_module(expr, genes)
    return pd.DataFrame(out)

# Load required cohort matrices (subset)
print("[loading expr matrices for Singscore]")

def load_gse126698():
    p = DATA_DIR/"GSE126698"/"GSE126698_DE_Thyroid_totalRNA_all.csv.gz"
    df = pd.read_csv(p)
    fpm_cols = [c for c in df.columns if c.endswith(".FPM")]
    sub = df[["names"]+fpm_cols].copy()
    sub = sub.dropna(subset=["names"])
    sub = sub.groupby("names").max(numeric_only=True)
    expr = sub.T
    expr.index = [c.replace(".FPM","") for c in expr.index]
    expr = np.log2(expr + 1.0)
    return expr

def load_hugo():
    import openpyxl
    p = DATA_DIR/"GSE78220"/"GSE78220_PatientFPKM.xlsx"
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    header = list(rows[0]); body = rows[1:]
    df = pd.DataFrame(body, columns=header).dropna(subset=["Gene"])
    df = df.groupby("Gene").max(numeric_only=True)
    expr = df.T
    expr.index = [str(c).replace(".baseline","") for c in expr.index]
    expr = expr.astype(float)
    return np.log2(expr + 1.0)

def load_riaz_all():
    p = DATA_DIR/"GSE91061"/"GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"
    df = pd.read_csv(p, index_col=0)
    df.index = df.index.astype(str)
    df = df.loc[df.index.intersection(ENTREZ2SYM.keys())]
    df["symbol"] = [ENTREZ2SYM[i] for i in df.index]
    df = df.groupby("symbol").max(numeric_only=True)
    expr = df.T.astype(float)
    return np.log2(expr + 1.0)

# (1) Singscore on GSE126698 + Hugo + Riaz_pre
print("[Singscore on GSE126698 + Hugo + Riaz]")
gse126_expr = load_gse126698()
gse126_ss = score_singscore_all(gse126_expr)
gse126_ss.index.name = "sample_id"
print(f"  GSE126698 Singscore shape: {gse126_ss.shape}")

hugo_expr = load_hugo()
hugo_ss = score_singscore_all(hugo_expr)

riaz_expr = load_riaz_all()
riaz_pre = riaz_expr.loc[[s for s in riaz_expr.index if "_Pre_" in s]]
riaz_on  = riaz_expr.loc[[s for s in riaz_expr.index if "_On_"  in s]]
riaz_pre_ss = score_singscore_all(riaz_pre)
riaz_on_ss  = score_singscore_all(riaz_on)

# Save Singscore tables
gse126_ss.to_csv(OUT_DIR / "singscore_GSE126698.tsv", sep="\t")
hugo_ss.to_csv(OUT_DIR / "singscore_Hugo.tsv", sep="\t")
riaz_pre_ss.to_csv(OUT_DIR / "singscore_Riaz_pre.tsv", sep="\t")
riaz_on_ss.to_csv(OUT_DIR / "singscore_Riaz_on.tsv", sep="\t")

# Sanity: ATC vs non-ATC in GSE126698 with Singscore
def gse126_subtype(s):
    return {"A":"ATC","F":"FTC","N":"NT","P":"PTC"}.get(s[0], "unknown")

gse126_ss_sub = gse126_ss.copy()
gse126_ss_sub["subtype"] = [gse126_subtype(s) for s in gse126_ss_sub.index]
ss_atc = gse126_ss_sub[gse126_ss_sub.subtype=="ATC"]
ss_other = gse126_ss_sub[gse126_ss_sub.subtype.isin(["PTC","FTC","NT"])]
ss_compare = []
for mod in MODULES:
    a = ss_atc[mod].dropna().values; b = ss_other[mod].dropna().values
    if len(a)<2 or len(b)<2: continue
    d = (a.mean()-b.mean()) / np.sqrt((a.std(ddof=1)**2 + b.std(ddof=1)**2)/2 + 1e-12)
    ss_compare.append({
        "module":mod, "n_atc":len(a), "n_other":len(b),
        "singscore_d":float(d), "ssgsea_d":None  # populated below
    })

# Compare with ssGSEA-like results from Phase 1
ssgsea_pooled = pd.read_csv(OUT_DIR / "pooled_module_scores_v2_with_subtype.tsv", sep="\t", index_col=0)
ssg_gse126 = ssgsea_pooled[ssgsea_pooled.cohort=="GSE126698"]
ssg_atc = ssg_gse126[ssg_gse126.subtype=="ATC"]
ssg_other = ssg_gse126[ssg_gse126.subtype.isin(["PTC","FTC","NT"])]
for row in ss_compare:
    mod = row["module"]
    a = ssg_atc[mod].dropna().values; b = ssg_other[mod].dropna().values
    if len(a)>=2 and len(b)>=2:
        d = (a.mean()-b.mean()) / np.sqrt((a.std(ddof=1)**2 + b.std(ddof=1)**2)/2 + 1e-12)
        row["ssgsea_d"] = float(d)
ss_cmp = pd.DataFrame(ss_compare)
ss_cmp["sign_match"] = ss_cmp.apply(
    lambda r: "yes" if (r.singscore_d>0)==(r.ssgsea_d>0) else "no",
    axis=1)
ss_cmp.to_csv(OUT_DIR / "singscore_vs_ssgsea_GSE126698_ATC_vs_other.tsv", sep="\t", index=False)
print("\nSingscore vs ssGSEA-like — GSE126698 ATC vs non-ATC d comparison:")
print(ss_cmp.round(2).to_string())

# Singscore DIAL-lite Hugo + Riaz_pre
print("\n[Singscore DIAL-lite]")
def parse_resp(acc):
    p = DATA_DIR / acc / f"{acc}_series_matrix.txt.gz"
    samples=[]; ch=[]
    with gzip.open(p,"rt",errors="replace") as fh:
        for line in fh:
            cols=line.rstrip("\n").split("\t")
            tag=cols[0]
            if tag=="!Sample_title": samples=[c.strip().strip('"') for c in cols[1:]]
            elif tag=="!Sample_characteristics_ch1": ch.append([c.strip().strip('"') for c in cols[1:]])
            if line.startswith("!series_matrix_table_begin"): break
    for row in ch:
        if any(("response" in c.lower() and "anti-pd-1 response" in c.lower()) or c.lower().startswith("response:") for c in row if c):
            return dict(zip(samples,row))
    return {}

def hugo_class(s):
    s=(s or "").lower()
    if "complete response" in s or "partial response" in s: return 1
    if "progressive disease" in s: return 0
    return np.nan
def riaz_class(s):
    s=(s or "").split(":")[-1].strip().lower()
    if s in ("prcr","cr","pr"): return 1
    if s in ("pd",): return 0
    return np.nan

hugo_resp = parse_resp("GSE78220")
riaz_resp = parse_resp("GSE91061")
hugo_ss_resp = hugo_ss.copy()
hugo_ss_resp["resp"] = [hugo_class(hugo_resp.get(s,"")) for s in hugo_ss_resp.index]
riaz_pre_ss_resp = riaz_pre_ss.copy()
riaz_pre_ss_resp["resp"] = [riaz_class(riaz_resp.get(s,"")) for s in riaz_pre_ss_resp.index]

dial_ss = []
for cname, scores in [("Hugo", hugo_ss_resp), ("Riaz_pre", riaz_pre_ss_resp)]:
    s = scores.dropna(subset=["resp"])
    pos = s[s.resp==1]; neg = s[s.resp==0]
    for mod in MODULES:
        a = pos[mod].dropna().values; b = neg[mod].dropna().values
        if len(a)<3 or len(b)<3: continue
        d = (a.mean()-b.mean()) / np.sqrt((a.std(ddof=1)**2 + b.std(ddof=1)**2)/2 + 1e-12)
        dial_ss.append({"cohort":cname,"module":mod,"singscore_d":float(d),"sign":"+" if d>0 else "-"})
dial_ss_df = pd.DataFrame(dial_ss)
piv_ss = dial_ss_df.pivot(index="module", columns="cohort", values="sign")
piv_ss["singscore_consistent"] = piv_ss.apply(lambda r: 1 if r.get("Hugo")==r.get("Riaz_pre") else 0, axis=1)
piv_ss.to_csv(OUT_DIR / "dial_lite_singscore_sign.tsv", sep="\t")
print(piv_ss.to_string())

# Combine vs ssGSEA-like sign consistency
piv_ssgsea = pd.read_csv(OUT_DIR / "dial_lite_sign_consistency_v2.tsv", sep="\t", index_col=0)
piv_ssgsea["ssgsea_consistent"] = piv_ssgsea["consistent"]
combined = piv_ssgsea[["ssgsea_consistent"]].join(piv_ss[["singscore_consistent"]])
combined["both_consistent"] = (combined.ssgsea_consistent.astype(int) & combined.singscore_consistent.astype(int))
combined.to_csv(OUT_DIR / "dial_lite_method_sensitivity.tsv", sep="\t")
print(combined.to_string())

# (2) Riaz Pre→On split by responder
print("\n[Riaz Pre→On delta — responder-stratified]")
delta_path = OUT_DIR / "riaz_pre_to_on_delta.tsv"
delta = pd.read_csv(delta_path, sep="\t")
# patient → responder class via Riaz response map
def patient_response(pt):
    # find any sample for this patient and use its response
    for sample, resp in riaz_resp.items():
        if sample.startswith(pt+"_"):
            return riaz_class(resp)
    return np.nan
delta["resp_class"] = [patient_response(pt) for pt in delta.patient]
delta.to_csv(OUT_DIR / "riaz_pre_to_on_delta_with_response.tsv", sep="\t", index=False)

# Per-module: responder vs non-responder mean delta + t-test
strat_rows = []
for col in [c for c in delta.columns if c.startswith("delta_")]:
    mod = col.replace("delta_","")
    r_pos = delta[delta.resp_class==1][col].dropna().values
    r_neg = delta[delta.resp_class==0][col].dropna().values
    if len(r_pos)>=2 and len(r_neg)>=2:
        d = (r_pos.mean()-r_neg.mean()) / np.sqrt((r_pos.std(ddof=1)**2 + r_neg.std(ddof=1)**2)/2 + 1e-12)
        t,p = stats.ttest_ind(r_pos, r_neg, equal_var=False)
    else:
        d=t=p=np.nan
    strat_rows.append({
        "module":mod,
        "n_resp":int(len(r_pos)), "mean_delta_resp":float(r_pos.mean()) if len(r_pos) else np.nan,
        "n_nonresp":int(len(r_neg)), "mean_delta_nonresp":float(r_neg.mean()) if len(r_neg) else np.nan,
        "delta_diff_d":float(d) if not np.isnan(d) else np.nan,
        "p":float(p) if not np.isnan(p) else np.nan,
    })
strat_df = pd.DataFrame(strat_rows)
strat_df.to_csv(OUT_DIR / "riaz_pre_to_on_responder_stratified.tsv", sep="\t", index=False)
print(strat_df.round(3).to_string())

# (3) Combined 6-figure montage
print("\n[6-figure montage]")
fig_files = [
    ("fig1_atc_vs_other_meta.png", "Fig 1 — ATC vs non-ATC pooled meta (4 thyroid cohorts)"),
    ("fig2_pca_scatter_by_cohort.png", "Fig 2 — Pooled PCA (n=415, 8 cohorts)"),
    ("fig3_subtype_heatmap.png", "Fig 3 — Module × cohort×subtype heatmap"),
    ("fig4_dial_lite_forest.png", "Fig 4 — DIAL-lite forest (Hugo + Riaz pre, 1000-bootstrap CI)"),
    ("fig5_dediff_vs_inflam_GSE126698.png", "Fig 5 — Dedifferentiation × inflammation (GSE126698)"),
    ("fig6_riaz_pre_to_on_delta.png", "Fig 6 — Riaz Pre→On per-patient delta (n=43 paired)"),
]
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
for ax, (fname, title) in zip(axes.flat, fig_files):
    img = mpimg.imread(FIG_DIR / fname)
    ax.imshow(img)
    ax.set_title(title, fontsize=10)
    ax.axis("off")
plt.suptitle("Paper 3 ICI Track B-lite — 6-figure overview (hypothesis-generating, not response prediction)",
             fontsize=12, y=0.995)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig_montage_6panel.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  saved fig_montage_6panel.png")

print("\n[done]")
