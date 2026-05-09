"""
Paper 1+2 BRAF Nature sprint H27 — immune checkpoint expression in DM1 vs DM2 BRAF-cPTC.

Question: Is DM1 (HT-overlap, immune-hot, B-cell/Tfh-rich) an ICI-candidate
phenotype within BRAF-mutant cPTC? We measure inhibitory + stimulatory
checkpoint, effector/cytotoxicity, antigen-presentation, Treg, and chemokine
panels, plus the Ayers 2017 Tumor Inflammation Signature (TIS, an FDA-validated
ICI-response surrogate).

Stratum
-------
master TSV: tissue_type=='Primary Tumor', dm in {DM1, DM2},
            molecular_subtype=='BRAF_like', histology_subtype=='cPTC' (n=110).

Outputs
-------
- h27_checkpoint_d.tsv     per-gene d, MW_p, BH_FDR, panel_class
- h27_composite_scores.tsv composite scores per sample + DM1 vs DM2 d / p
- h27_ici_likely_responder.tsv responder fraction + PFI association
- H27_REPORT.md            <400 word report
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

OUT = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_braf_nature_sprint_2026_05_09/h27_ici_panels"
)
OUT.mkdir(parents=True, exist_ok=True)

MASTER = "/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv"
EXPR_Z = "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv"
H6_CLIN = (
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_braf_nature_sprint_2026_05_09/h6_survival/h6_merged_clinical.tsv"
)
H10_PS = (
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_braf_nature_sprint_2026_05_09/h10_celltype_deconv/h10_per_sample_scores.tsv"
)


# ---------------------------------------------------------------------- panels
INHIBITORY = [
    "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "BTLA",
    "C10orf54",  # VISTA
    "VTCN1",     # B7-H4
    "IDO1", "IDO2",
]
STIMULATORY = [
    "CD40", "CD80", "CD86", "CD27", "CD28", "ICOS", "ICOSLG",
    "TNFRSF4",   # OX40
    "TNFRSF9",   # 4-1BB
    "TNFRSF18",  # GITR
]
EFFECTOR = ["GZMB", "PRF1", "GNLY", "IFNG", "TNF", "NKG7"]
ANTIGEN_PRES = ["TAP1", "TAP2", "B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA"]
TREG_SUPPR = ["FOXP3", "IL2RA", "IL10", "TGFB1"]
CYTOKINE_CHEMO = ["IFNG", "CXCL9", "CXCL10", "CXCL11", "CCL5"]

PANELS = {
    "inhibitory":   INHIBITORY,
    "stimulatory":  STIMULATORY,
    "effector":     EFFECTOR,
    "antigen_pres": ANTIGEN_PRES,
    "treg_suppr":   TREG_SUPPR,
    "cytok_chemo":  CYTOKINE_CHEMO,
}

# Ayers 2017 18-gene TIS (canonical: 18 genes; user listed 20 incl. extras)
TIS_AYERS = [
    "IFNG", "CXCL9", "CD8A", "GZMA", "GZMK", "HLA-DRA", "NKG7", "PSMB10",
    "IDO1", "STAT1", "CCL5", "TIGIT", "LAG3", "PDCD1LG2", "CD274", "CMKLR1",
    "CD276", "CXCR6", "HLA-DOB", "HLA-E",
]

# "Immune-hot" composite (as specified)
HOT5 = ["GZMB", "PRF1", "IFNG", "CD274", "PDCD1"]


# ---------------------------------------------------------------------- helpers
def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    sa, sb = a.var(ddof=1), b.var(ddof=1)
    pooled = np.sqrt(((len(a) - 1) * sa + (len(b) - 1) * sb) / (len(a) + len(b) - 2))
    if pooled == 0:
        return float("nan")
    return float((a.mean() - b.mean()) / pooled)


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    finite = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    if finite.sum() == 0:
        return out
    pf = p[finite]
    n = len(pf)
    order = np.argsort(pf)
    q = pf * n / np.arange(1, n + 1)[np.argsort(order)]  # placeholder, recompute below
    # proper BH
    sorted_p = pf[order]
    q_sorted = sorted_p * n / np.arange(1, n + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_unsorted = np.empty(n)
    q_unsorted[order] = np.clip(q_sorted, 0, 1)
    out[finite] = q_unsorted
    return out


def mw_p(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float); a = a[~np.isnan(a)]
    b = np.asarray(b, dtype=float); b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    try:
        _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        return float(p)
    except ValueError:
        return float("nan")


def gene_mean_z(genes: list, expr_df: pd.DataFrame) -> pd.Series:
    present = [g for g in genes if g in expr_df.index]
    if not present:
        return pd.Series(dtype=float)
    return expr_df.loc[present].mean(axis=0)


# ---------------------------------------------------------------------- load
print("[load] master")
master = pd.read_csv(MASTER, sep="\t", low_memory=False)
master = master[master["tissue_type"] == "Primary Tumor"].copy()
master = master[master["dm"].isin(["DM1", "DM2"])].copy()
braf_cptc = master[
    (master["molecular_subtype"] == "BRAF_like")
    & (master["histology_subtype"] == "cPTC")
].copy()
n_dm1 = (braf_cptc["dm"] == "DM1").sum()
n_dm2 = (braf_cptc["dm"] == "DM2").sum()
print(f"  BRAF_like+cPTC DM1/DM2 master n={len(braf_cptc)} (DM1={n_dm1}, DM2={n_dm2})")

print("[load] z-scored expression")
expr = pd.read_csv(EXPR_Z, sep="\t", index_col=0)
print(f"  expression shape: {expr.shape}")

# Restrict expr to BRAF_cPTC sample columns that exist
sids = [s for s in braf_cptc["sample_id"].tolist() if s in expr.columns]
print(f"  expr ∩ braf_cptc samples: {len(sids)}/{len(braf_cptc)}")
braf_cptc = braf_cptc[braf_cptc["sample_id"].isin(sids)].reset_index(drop=True)
dm1_ids = braf_cptc[braf_cptc["dm"] == "DM1"]["sample_id"].tolist()
dm2_ids = braf_cptc[braf_cptc["dm"] == "DM2"]["sample_id"].tolist()
print(f"  final stratum: DM1={len(dm1_ids)}  DM2={len(dm2_ids)}")

expr_b = expr[sids]


# ---------------------------------------------------------------------- 1) per-gene
print("[step 1] per-gene Cohen d (DM1-DM2) within BRAF-cPTC")
all_genes_panels = []
for panel, glist in PANELS.items():
    for g in glist:
        all_genes_panels.append((panel, g))

rows = []
for panel, g in all_genes_panels:
    if g not in expr_b.index:
        rows.append({
            "panel_class": panel, "gene": g,
            "n_DM1": 0, "n_DM2": 0,
            "mean_DM1": float("nan"), "mean_DM2": float("nan"),
            "cohen_d": float("nan"), "MW_p": float("nan"),
            "missing": True,
        })
        continue
    a = expr_b.loc[g, dm1_ids].astype(float).values
    b = expr_b.loc[g, dm2_ids].astype(float).values
    rows.append({
        "panel_class": panel, "gene": g,
        "n_DM1": int((~np.isnan(a)).sum()),
        "n_DM2": int((~np.isnan(b)).sum()),
        "mean_DM1": float(np.nanmean(a)) if len(a) else float("nan"),
        "mean_DM2": float(np.nanmean(b)) if len(b) else float("nan"),
        "cohen_d": cohens_d(a, b),
        "MW_p": mw_p(a, b),
        "missing": False,
    })
gdf = pd.DataFrame(rows)
gdf["BH_FDR"] = bh_fdr(gdf["MW_p"].values)
gdf = gdf.sort_values("cohen_d", ascending=False, na_position="last").reset_index(drop=True)
gdf.to_csv(OUT / "h27_checkpoint_d.tsv", sep="\t", index=False)
print(f"  wrote {OUT/'h27_checkpoint_d.tsv'}  rows={len(gdf)}")
print("\n=== TOP 20 ranked by Cohen d (DM1>DM2) ===")
print(gdf.head(20)[["panel_class", "gene", "n_DM1", "n_DM2",
                    "mean_DM1", "mean_DM2", "cohen_d", "MW_p", "BH_FDR"]]
        .to_string(index=False))


# ---------------------------------------------------------------------- 2) composites
print("\n[step 2] composite scores per sample")
per_sample = pd.DataFrame(index=sids)

# panel mean-z scores
for panel, glist in PANELS.items():
    s = gene_mean_z(glist, expr_b)
    per_sample[f"panel_{panel}"] = s.reindex(sids)

# ICI-favorable composite = inhibitory + effector
per_sample["ICI_favorable"] = (per_sample["panel_inhibitory"]
                                + per_sample["panel_effector"]) / 2.0
# Immune-hot 5-gene
per_sample["immune_hot5"] = gene_mean_z(HOT5, expr_b).reindex(sids)
# TIS
tis_present = [g for g in TIS_AYERS if g in expr_b.index]
print(f"  TIS Ayers genes present: {len(tis_present)}/{len(TIS_AYERS)} -> {tis_present}")
per_sample["TIS_Ayers"] = gene_mean_z(TIS_AYERS, expr_b).reindex(sids)

# DM annotation
per_sample["dm"] = braf_cptc.set_index("sample_id").reindex(sids)["dm"]

# composite-level d test
comp_rows = []
for col in ["panel_inhibitory", "panel_stimulatory", "panel_effector",
            "panel_antigen_pres", "panel_treg_suppr", "panel_cytok_chemo",
            "ICI_favorable", "immune_hot5", "TIS_Ayers"]:
    a = per_sample.loc[per_sample["dm"] == "DM1", col].astype(float).values
    b = per_sample.loc[per_sample["dm"] == "DM2", col].astype(float).values
    comp_rows.append({
        "composite": col,
        "n_DM1": int((~np.isnan(a)).sum()),
        "n_DM2": int((~np.isnan(b)).sum()),
        "mean_DM1": float(np.nanmean(a)),
        "mean_DM2": float(np.nanmean(b)),
        "cohen_d": cohens_d(a, b),
        "MW_p": mw_p(a, b),
    })
cdf = pd.DataFrame(comp_rows)
cdf["BH_FDR"] = bh_fdr(cdf["MW_p"].values)
per_sample.to_csv(OUT / "h27_composite_scores.tsv", sep="\t")
print(f"  wrote {OUT/'h27_composite_scores.tsv'}  shape={per_sample.shape}")
print("\n=== Composite-level d (DM1-DM2) within BRAF-cPTC ===")
print(cdf.to_string(index=False))


# ---------------------------------------------------------------------- 3) cross with H10 deconv
print("\n[step 3] cross-reference with H10 deconv (CD8 / Tfh)")
try:
    h10 = pd.read_csv(H10_PS, sep="\t", index_col=0)
    common = per_sample.index.intersection(h10.index)
    print(f"  H10 ∩ H27 samples: {len(common)}")
    merged = per_sample.loc[common].join(
        h10.loc[common, ["GS_CD8_T", "GS_CD4_Tfh", "TLS_Cabrita12"]],
        how="left",
    )
    cross = []
    for col_h27 in ["TIS_Ayers", "immune_hot5", "panel_inhibitory", "ICI_favorable"]:
        for col_h10 in ["GS_CD8_T", "GS_CD4_Tfh", "TLS_Cabrita12"]:
            sub = merged[[col_h27, col_h10]].dropna()
            if len(sub) < 10:
                continue
            r, p = stats.spearmanr(sub[col_h27], sub[col_h10])
            cross.append({"x": col_h27, "y": col_h10, "n": int(len(sub)),
                          "spearman_rho": float(r), "p": float(p)})
    cross_df = pd.DataFrame(cross)
    print(cross_df.to_string(index=False))
except FileNotFoundError:
    print("  H10 per-sample scores not found; skipping cross-ref")
    cross_df = pd.DataFrame()
    merged = per_sample.copy()
    merged["GS_CD8_T"] = np.nan
    merged["GS_CD4_Tfh"] = np.nan
    merged["TLS_Cabrita12"] = np.nan


# ---------------------------------------------------------------------- 4) ICI-likely-responder strata
print("\n[step 4] ICI-likely-responder strata + PFI association")

# Define cutoffs at within-stratum medians (BRAF-cPTC)
def quartile_high(s):
    return s >= s.median()

merged["TIS_high"] = quartile_high(merged["TIS_Ayers"])
merged["CD8_high"] = quartile_high(merged["GS_CD8_T"]) if "GS_CD8_T" in merged else False
merged["HLAI_high"] = quartile_high(merged["panel_antigen_pres"])
merged["dm"] = per_sample.reindex(merged.index)["dm"]

merged["ICI_likely"] = (
    (merged["dm"] == "DM1")
    & merged["TIS_high"].fillna(False)
    & merged["CD8_high"].fillna(False)
    & merged["HLAI_high"].fillna(False)
)

n_total = len(merged)
n_likely = int(merged["ICI_likely"].sum())
n_dm1 = int((merged["dm"] == "DM1").sum())
n_dm2 = int((merged["dm"] == "DM2").sum())
frac_total = n_likely / n_total if n_total else float("nan")
frac_dm1 = n_likely / n_dm1 if n_dm1 else float("nan")
print(f"  ICI-likely responders: {n_likely}/{n_total} ({frac_total:.1%}) of all BRAF-cPTC")
print(f"  ICI-likely / DM1 only: {n_likely}/{n_dm1} ({frac_dm1:.1%})")

# Merge PFI from H6
print("[load] H6 merged clinical")
h6 = pd.read_csv(H6_CLIN, sep="\t")
h6 = h6.set_index("sample_id")
merged_clin = merged.join(h6[["PFI", "PFI.time", "DSS", "DSS.time", "OS", "OS.time"]],
                          how="left")

# Wilcoxon: PFI events in ICI-likely vs others within DM1
def or_or_p(group_a, group_b, col="PFI"):
    a_ev = group_a[col].dropna().astype(float)
    b_ev = group_b[col].dropna().astype(float)
    a_pos = (a_ev == 1).sum(); a_n = len(a_ev)
    b_pos = (b_ev == 1).sum(); b_n = len(b_ev)
    table = np.array([[a_pos, a_n - a_pos], [b_pos, b_n - b_pos]])
    try:
        _, p = stats.fisher_exact(table)
    except ValueError:
        p = float("nan")
    return a_pos, a_n, b_pos, b_n, float(p)

# ICI-likely vs all-other-BRAF-cPTC
g_a = merged_clin[merged_clin["ICI_likely"]]
g_b = merged_clin[~merged_clin["ICI_likely"]]

resp_rows = []
for col in ["PFI", "DSS", "OS"]:
    if col not in merged_clin.columns:
        continue
    a_pos, a_n, b_pos, b_n, p = or_or_p(g_a, g_b, col=col)
    resp_rows.append({
        "endpoint": col,
        "n_ICIlikely": a_n,
        "events_ICIlikely": a_pos,
        "rate_ICIlikely": a_pos / a_n if a_n else float("nan"),
        "n_other": b_n,
        "events_other": b_pos,
        "rate_other": b_pos / b_n if b_n else float("nan"),
        "fisher_p": p,
    })
resp_df = pd.DataFrame(resp_rows)
print(resp_df.to_string(index=False))

# Logistic-y: PFI ~ ICI_likely (within DM1 only)
dm1_only = merged_clin[merged_clin["dm"] == "DM1"]
print(f"\n  Within DM1 only: ICI_likely={int(dm1_only['ICI_likely'].sum())}/{len(dm1_only)}")
dm1_resp = []
for col in ["PFI", "DSS", "OS"]:
    if col not in dm1_only.columns:
        continue
    g_a = dm1_only[dm1_only["ICI_likely"]]
    g_b = dm1_only[~dm1_only["ICI_likely"]]
    a_pos, a_n, b_pos, b_n, p = or_or_p(g_a, g_b, col=col)
    dm1_resp.append({
        "endpoint": col,
        "stratum": "DM1_only",
        "n_ICIlikely": a_n,
        "events_ICIlikely": a_pos,
        "rate_ICIlikely": a_pos / a_n if a_n else float("nan"),
        "n_other_DM1": b_n,
        "events_other_DM1": b_pos,
        "rate_other_DM1": b_pos / b_n if b_n else float("nan"),
        "fisher_p": p,
    })
dm1_resp_df = pd.DataFrame(dm1_resp)
if len(dm1_resp_df):
    print(dm1_resp_df.to_string(index=False))

# Save
out_resp = pd.DataFrame({
    "metric": [
        "n_total_BRAF_cPTC",
        "n_DM1",
        "n_DM2",
        "n_ICI_likely_responder",
        "fraction_of_total",
        "fraction_of_DM1",
    ],
    "value": [
        n_total, n_dm1, n_dm2, n_likely,
        round(frac_total, 4), round(frac_dm1, 4),
    ],
})
out_resp.to_csv(OUT / "h27_ici_likely_responder.tsv", sep="\t", index=False)
# Append PFI association rows as section 2
with open(OUT / "h27_ici_likely_responder.tsv", "a") as f:
    f.write("\n# clinical association — ICI-likely vs other (all BRAF-cPTC)\n")
    resp_df.to_csv(f, sep="\t", index=False)
    if len(dm1_resp_df):
        f.write("\n# clinical association — within DM1 only (ICI-likely vs other DM1)\n")
        dm1_resp_df.to_csv(f, sep="\t", index=False)
print(f"  wrote {OUT/'h27_ici_likely_responder.tsv'}")


# ---------------------------------------------------------------------- 5) headline json
tis_row = cdf[cdf["composite"] == "TIS_Ayers"].iloc[0]
hot_row = cdf[cdf["composite"] == "immune_hot5"].iloc[0]
inh_row = cdf[cdf["composite"] == "panel_inhibitory"].iloc[0]
ici_row = cdf[cdf["composite"] == "ICI_favorable"].iloc[0]

# Top 5 inhibitory
top_inhib = gdf[(gdf["panel_class"] == "inhibitory")
                & (~gdf["missing"])].sort_values("cohen_d", ascending=False).head(5)
top_effector = gdf[(gdf["panel_class"] == "effector")
                   & (~gdf["missing"])].sort_values("cohen_d", ascending=False).head(5)

headline = {
    "stratum": f"BRAF_like × cPTC × DM1/DM2 (n={n_total}; DM1={n_dm1}, DM2={n_dm2})",
    "TIS_Ayers": {
        "d_DM1_minus_DM2": float(tis_row["cohen_d"]),
        "MW_p": float(tis_row["MW_p"]),
        "BH_FDR": float(tis_row["BH_FDR"]),
        "mean_DM1": float(tis_row["mean_DM1"]),
        "mean_DM2": float(tis_row["mean_DM2"]),
        "n_genes_present": len(tis_present),
    },
    "immune_hot5": {
        "d_DM1_minus_DM2": float(hot_row["cohen_d"]),
        "MW_p": float(hot_row["MW_p"]),
    },
    "panel_inhibitory_mean_z": {
        "d_DM1_minus_DM2": float(inh_row["cohen_d"]),
        "MW_p": float(inh_row["MW_p"]),
    },
    "ICI_favorable_composite": {
        "d_DM1_minus_DM2": float(ici_row["cohen_d"]),
        "MW_p": float(ici_row["MW_p"]),
    },
    "top_inhibitory_DM1_up": top_inhib[["gene", "cohen_d", "BH_FDR"]]
        .to_dict(orient="records"),
    "top_effector_DM1_up": top_effector[["gene", "cohen_d", "BH_FDR"]]
        .to_dict(orient="records"),
    "ici_likely_responder_n": n_likely,
    "ici_likely_responder_frac_total": frac_total,
    "ici_likely_responder_frac_DM1": frac_dm1,
}
(OUT / "h27_headline.json").write_text(json.dumps(headline, indent=2))
print(f"\n  wrote {OUT/'h27_headline.json'}")
print(json.dumps(headline, indent=2))
