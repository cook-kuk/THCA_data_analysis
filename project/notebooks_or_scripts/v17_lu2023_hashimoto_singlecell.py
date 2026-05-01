#!/usr/bin/env python3
"""Q15 — Lu 2023 PTC single-cell Hashimoto-like signature × DM-axis cross-validation.

Q12 found 17% of Korean GSE213647 PTC (n=59/348) had HLA-II z>1 "Hashimoto-like" signature
with HLA-I p=4.7e-26 + DM1 panel_z p=9.6e-4. Test if same subset exists at single-cell scale.

Single-cell hypothesis (testable now):
  H_sc1: PTC cells with HLA-II z>1 (cell-level Hashimoto-like) have higher DM_score (less DM1)
  H_sc2: ~17% of PTC patient-level samples have enrichment for Hashimoto-like cells
  H_sc3: Hashimoto-like single-cells co-cluster on UMAP (immune-engaged tumor microenvironment)
"""
from __future__ import annotations
import json, sys, logging
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats, sparse
import plotly.graph_objects as go

ROOT = Path("/opt/thyroid-dash/project")
RES  = ROOT/"results/v17_lu2023_hashimoto"; RES.mkdir(parents=True, exist_ok=True)
FIG  = ROOT/"reports/html/figs_interactive/v17"
RPT  = ROOT/"reports/v17p35"
LOG  = ROOT/"logs/v17_lu2023_hashimoto.log"
ADATA_PATH = ROOT/"results/v17_lu2023/GSE193581_hvg_adata.h5ad"
RAW_DIR = Path("/data/thca/v17_lu2023_GSE193581/raw")
ANN_FILE = Path("/data/thca/v17_lu2023_GSE193581/celltype_annotation.txt.gz")

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.FileHandler(LOG, mode="w"), logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)
log.info("=== Lu 2023 single-cell Hashimoto signature × DM-axis ===")

BG, INK = "#0b0e12", "#F2F2F2"
DARK = dict(template="plotly_dark", paper_bgcolor=BG, plot_bgcolor=BG, font=dict(color=INK, size=14))
SEED = 42

PANEL = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
HLA_I = ["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","NLRC5"]
HLA_II = ["HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","CIITA"]

# ---------- step 1: rebuild full normalized expression for PANEL + HLA, all cells (memory-slim)
log.info("=== Step 1: rebuild log-normed expression, panel + HLA only (slim) ===")
ann = pd.read_csv(ANN_FILE, sep="\t")
ann.columns = ["sample_id","celltype"]
ann_by_sample = {s: g["celltype"].values for s,g in ann.groupby("sample_id")}

def histo(s):
    if s.startswith("PTC"): return "PTC"
    if s.startswith("ATC"): return "ATC"
    if s.startswith("NORM"): return "NORM"
    return "OTHER"

ALL_GENES = PANEL + HLA_I + HLA_II
log.info(f"target genes: {len(ALL_GENES)} ({len(PANEL)} panel + {len(HLA_I)} HLA-I + {len(HLA_II)} HLA-II)")

samples_data = []
for fp in sorted(RAW_DIR.glob("*_UMI.txt.gz")):
    sname = fp.name.split("_")[1]
    ann_key = sname if sname in ann_by_sample else (sname + "T")
    if ann_key not in ann_by_sample: continue
    df = pd.read_csv(fp, sep="\t", index_col=0)
    expected = ann_by_sample[ann_key]
    n = min(df.shape[1], len(expected))
    df = df.iloc[:, :n]
    cells = df.columns.tolist()[:n]
    ctypes = expected[:n]
    # QC
    cell_n_genes = (df > 0).sum(axis=0).values
    mt_mask = df.index.str.startswith("MT-")
    cell_pct_mt = 100 * df.loc[mt_mask].sum(axis=0).values / np.maximum(df.sum(axis=0).values, 1)
    keep = (cell_n_genes > 200) & (cell_n_genes < 6000) & (cell_pct_mt < 20)
    df = df.loc[:, keep]
    cells = [c for c,k in zip(cells, keep) if k]
    ctypes = [c for c,k in zip(ctypes, keep) if k]
    if df.shape[1] == 0: continue
    libsize = df.sum(axis=0).values
    df_norm = df.div(libsize, axis=1) * 1e4
    df_log = np.log1p(df_norm.values).astype(np.float32)
    present_genes = [g for g in ALL_GENES if g in df.index]
    gene_X = pd.DataFrame(df_log[df.index.get_indexer(present_genes), :],
                          index=present_genes, columns=cells)
    samples_data.append({
        "sample": sname, "histology": histo(sname),
        "cells": cells, "X": gene_X, "celltypes": np.array(ctypes),
    })
    log.info(f"  {sname} ({histo(sname)}): {df.shape[1]} cells, {len(present_genes)}/{len(ALL_GENES)} genes")

# Merge
all_X = pd.concat([sd["X"] for sd in samples_data], axis=1)
all_meta = pd.DataFrame({
    "cell": np.concatenate([sd["cells"] for sd in samples_data]),
    "sample": np.concatenate([[sd["sample"]]*len(sd["cells"]) for sd in samples_data]),
    "histology": np.concatenate([[sd["histology"]]*len(sd["cells"]) for sd in samples_data]),
    "author_celltype": np.concatenate([sd["celltypes"] for sd in samples_data]),
}).set_index("cell")
all_X = all_X.loc[:, all_meta.index]
log.info(f"merged: {all_X.shape}")
log.info(f"genes detected: panel={[g for g in PANEL if g in all_X.index]}")
log.info(f"  HLA-I detected: {[g for g in HLA_I if g in all_X.index]}")
log.info(f"  HLA-II detected: {[g for g in HLA_II if g in all_X.index]}")

# ---------- step 2: filter to malignant cells ----------
malignant_ct = ["Malignant cell","Epithelial cell","Thyroid follicular cell","Thyrocyte","Tumor cell","Cancer cell","Malignant"]
mal_mask = all_meta["author_celltype"].isin(malignant_ct)
if mal_mask.sum() == 0:
    log.warning("no cells matched malignant labels; using TG-high proxy")
    tg = all_X.loc["TG"] if "TG" in all_X.index else None
    if tg is not None:
        mal_mask = (tg >= tg.quantile(0.5))
    else:
        mal_mask = pd.Series(True, index=all_meta.index)
mal_meta = all_meta[mal_mask].copy()
mal_X = all_X.loc[:, mal_meta.index]
log.info(f"malignant cells: {mal_mask.sum()} ({mal_meta['histology'].value_counts().to_dict()})")

# ---------- step 3: compute z-scores per gene group (cohort-level z over malignant cells) ----------
def gene_score(X, gene_list):
    present = [g for g in gene_list if g in X.index]
    if not present: return None
    sub = X.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1) + 1e-9, axis=0)
    return z.mean(axis=0), present

dm_score, dm_genes = gene_score(mal_X, PANEL)
hlaI_score, hlaI_genes = gene_score(mal_X, HLA_I)
hlaII_score, hlaII_genes = gene_score(mal_X, HLA_II)

mal_meta["DM_score"] = dm_score.values
mal_meta["HLA_I_score"] = hlaI_score.values
mal_meta["HLA_II_score"] = hlaII_score.values
log.info(f"DM_score median={mal_meta['DM_score'].median():.3f}, HLA-I={mal_meta['HLA_I_score'].median():.3f}, HLA-II={mal_meta['HLA_II_score'].median():.3f}")

# DM class (median split)
mal_meta["DM_class"] = np.where(mal_meta["DM_score"] >= mal_meta["DM_score"].median(), "DM1_high", "DM2_low")
# Hashimoto-like (HLA-II z > 1)
mal_meta["hashi_like"] = mal_meta["HLA_II_score"] > 1.0

# ---------- step 4: Test H_sc1 — Hashimoto-like cells have higher DM_score? (PTC only)
ptc = mal_meta[mal_meta["histology"]=="PTC"].copy()
log.info(f"\n=== H_sc1 (PTC cells only): Hashimoto-like vs other ===")
log.info(f"PTC malignant cells: {len(ptc)}")
hashi_pct_cells = ptc["hashi_like"].mean()
log.info(f"PTC Hashimoto-like cell %: {100*hashi_pct_cells:.1f}%")
g_hashi = ptc[ptc["hashi_like"]]["DM_score"].values
g_other = ptc[~ptc["hashi_like"]]["DM_score"].values
u, p_h_dm = stats.mannwhitneyu(g_hashi, g_other, alternative="two-sided")
log.info(f"H_sc1 PTC: Hashimoto-like DM_score median={np.median(g_hashi):+.3f} vs other={np.median(g_other):+.3f}, MW p={p_h_dm:.3e}")

# Same: HLA-I score in Hashimoto-like vs other
g_hashi_I = ptc[ptc["hashi_like"]]["HLA_I_score"].values
g_other_I = ptc[~ptc["hashi_like"]]["HLA_I_score"].values
u_I, p_h_I = stats.mannwhitneyu(g_hashi_I, g_other_I, alternative="two-sided")
log.info(f"H_sc1 PTC: HLA-I: Hashimoto-like median={np.median(g_hashi_I):+.3f} vs other={np.median(g_other_I):+.3f}, MW p={p_h_I:.3e}")

# ---------- step 5: H_sc2 — patient-level Hashimoto enrichment ----------
log.info(f"\n=== H_sc2 (patient-level): % PTC samples with Hashimoto-cell enrichment ===")
per_sample = mal_meta[mal_meta["histology"]=="PTC"].groupby("sample").agg(
    n_cells=("DM_score","count"),
    median_DM=("DM_score","median"),
    pct_hashi_like=("hashi_like", "mean"),
    median_HLA_II=("HLA_II_score","median"),
).sort_values("pct_hashi_like", ascending=False)
log.info(f"Per-PTC-sample Hashimoto-like cell % distribution:\n{per_sample.round(3)}")
n_ptc_samples = len(per_sample)
n_hashi_enriched = int((per_sample["pct_hashi_like"] > 0.20).sum())  # >20% cells are Hashimoto-like
log.info(f"PTC samples with >20% Hashimoto-like cells: {n_hashi_enriched}/{n_ptc_samples} ({100*n_hashi_enriched/n_ptc_samples:.1f}%)")

# ---------- save TSVs ----------
mal_meta[["sample","histology","author_celltype","DM_score","HLA_I_score","HLA_II_score","DM_class","hashi_like"]].to_csv(
    RES/"GSE193581_cells_DM_HLA.tsv", sep="\t")
per_sample.to_csv(RES/"GSE193581_per_PTC_sample.tsv", sep="\t")

summary = {
    "cohort": "GSE193581 (Lu 2023, n=23 sample, 67k cells, 15.3k malignant, n=8.6k PTC malignant)",
    "panel_genes_present": dm_genes,
    "HLA_I_genes_present": hlaI_genes,
    "HLA_II_genes_present": hlaII_genes,
    "n_PTC_malignant_cells": int(len(ptc)),
    "n_PTC_samples": int(n_ptc_samples),
    "PTC_Hashimoto_like_cell_pct": round(100*hashi_pct_cells, 2),
    "PTC_samples_with_>20pct_Hashimoto_cells": int(n_hashi_enriched),
    "PTC_samples_with_>20pct_Hashimoto_pct": round(100*n_hashi_enriched/n_ptc_samples, 1),
    "Hsc1_DM_score_Hashimoto_vs_other": {
        "Hashimoto_median": round(float(np.median(g_hashi)), 3),
        "other_median": round(float(np.median(g_other)), 3),
        "MW_p": float(p_h_dm),
    },
    "Hsc1_HLA_I_Hashimoto_vs_other": {
        "Hashimoto_median": round(float(np.median(g_hashi_I)), 3),
        "other_median": round(float(np.median(g_other_I)), 3),
        "MW_p": float(p_h_I),
    },
    "comparison_to_Q12": "Q12 (Korean GSE213647 bulk PTC n=348): 17.0% Hashimoto-like, HLA-I MW p=4.7e-26, DM panel_z MW p=9.6e-4. Test if Lu 2023 single-cell PTC reproduces same subgroup.",
    "verdict": None,  # will fill below
}

# Compose verdict
verdict_parts = []
if hashi_pct_cells >= 0.10 and hashi_pct_cells <= 0.30:
    verdict_parts.append(f"Cell-level Hashimoto-like fraction ({100*hashi_pct_cells:.0f}%) in cross-validating Q12 Korean range (~17%)")
elif hashi_pct_cells > 0.30:
    verdict_parts.append(f"Cell-level Hashimoto-like fraction higher than Q12 ({100*hashi_pct_cells:.0f}% vs Q12 17%) — Lu cohort enriched")
else:
    verdict_parts.append(f"Cell-level Hashimoto-like fraction lower than Q12 ({100*hashi_pct_cells:.0f}% vs Q12 17%)")
if p_h_dm < 0.05:
    direction = "higher (less dedifferentiated)" if np.median(g_hashi) > np.median(g_other) else "lower (more dedifferentiated)"
    verdict_parts.append(f"Hashimoto-like cells DM_score {direction}, MW p={p_h_dm:.2e}")
else:
    verdict_parts.append(f"Hashimoto vs other DM_score not significant (MW p={p_h_dm:.2f})")
if p_h_I < 0.05 and np.median(g_hashi_I) > np.median(g_other_I):
    verdict_parts.append(f"Hashimoto-like cells HLA-I higher (consistent with immune-engaged), MW p={p_h_I:.2e}")
summary["verdict"] = "; ".join(verdict_parts)
log.info(f"\nVERDICT: {summary['verdict']}")

(RES/"GSE193581_hashimoto_summary.json").write_text(json.dumps(summary, indent=2))

# ---------- figures ----------
COLOR_HASHI = "#7ccfcd"; COLOR_OTHER = "#aaa"

# Fig A: PTC cells DM_score split by Hashimoto-like
fig = go.Figure()
fig.add_trace(go.Violin(y=g_hashi, name=f"Hashimoto-like (HLA-II z>1, n={len(g_hashi):,})",
                        marker=dict(color=COLOR_HASHI), box_visible=True, line_color="white", opacity=0.85))
fig.add_trace(go.Violin(y=g_other, name=f"Other PTC (n={len(g_other):,})",
                        marker=dict(color=COLOR_OTHER), box_visible=True, line_color="white", opacity=0.85))
fig.update_layout(title=dict(
    text=f"<b>Q15 — Lu 2023 PTC single-cell DM_score: Hashimoto-like vs other (n={len(ptc):,} cells)</b><br>"
         f"<sub style='color:#7ccfcd'>{100*hashi_pct_cells:.1f}% PTC cells Hashimoto-like (HLA-II z>1); MW p={p_h_dm:.2e}; "
         f"Hashimoto median={np.median(g_hashi):+.3f}, other={np.median(g_other):+.3f}</sub>",
    font=dict(size=14, color=INK)),
    yaxis=dict(title="DM_score (8-panel mean z)", gridcolor="rgba(255,255,255,0.08)"),
    showlegend=True, height=460, width=820, margin=dict(l=80,r=20,t=110,b=60), **DARK,
    legend=dict(bgcolor="rgba(0,0,0,0.3)"))
fig.write_html(FIG/"v17_lu2023_hashimoto_dm_split.html", include_plotlyjs="cdn", full_html=True)
log.info(f"wrote v17_lu2023_hashimoto_dm_split.html")

# Fig B: same for HLA-I
fig = go.Figure()
fig.add_trace(go.Violin(y=g_hashi_I, name=f"Hashimoto-like (n={len(g_hashi_I):,})",
                        marker=dict(color=COLOR_HASHI), box_visible=True, line_color="white", opacity=0.85))
fig.add_trace(go.Violin(y=g_other_I, name=f"Other PTC (n={len(g_other_I):,})",
                        marker=dict(color=COLOR_OTHER), box_visible=True, line_color="white", opacity=0.85))
fig.update_layout(title=dict(
    text=f"<b>Q15 — Lu 2023 PTC HLA-I expression: Hashimoto-like vs other</b><br>"
         f"<sub style='color:#7ccfcd'>MW p={p_h_I:.2e}; Hashimoto median={np.median(g_hashi_I):+.3f}, other={np.median(g_other_I):+.3f}</sub>",
    font=dict(size=14, color=INK)),
    yaxis=dict(title="HLA Class I score (mean z across 7 genes)", gridcolor="rgba(255,255,255,0.08)"),
    showlegend=True, height=460, width=820, margin=dict(l=80,r=20,t=110,b=60), **DARK,
    legend=dict(bgcolor="rgba(0,0,0,0.3)"))
fig.write_html(FIG/"v17_lu2023_hashimoto_hlaI.html", include_plotlyjs="cdn", full_html=True)
log.info(f"wrote v17_lu2023_hashimoto_hlaI.html")

# Fig C: per-PTC-sample bar — % Hashimoto-like cells
fig = go.Figure()
samples_sorted = per_sample.sort_values("pct_hashi_like", ascending=False)
colors = [COLOR_HASHI if v > 0.20 else COLOR_OTHER for v in samples_sorted["pct_hashi_like"]]
fig.add_trace(go.Bar(x=samples_sorted.index, y=samples_sorted["pct_hashi_like"]*100,
                     marker=dict(color=colors), name="% Hashimoto-like cells",
                     text=[f"{v*100:.0f}%" for v in samples_sorted["pct_hashi_like"]], textposition="outside"))
fig.add_hline(y=20, line=dict(color="white", dash="dash", width=2),
              annotation=dict(text="20% threshold", font=dict(color=INK, size=11), bgcolor="rgba(0,0,0,0.4)"))
fig.update_layout(title=dict(
    text=f"<b>Q15 — Per-PTC-sample Hashimoto-like cell % (n={n_ptc_samples} samples)</b><br>"
         f"<sub style='color:#7ccfcd'>Teal = sample with >20% Hashimoto-like cells (n={n_hashi_enriched}/{n_ptc_samples} = {100*n_hashi_enriched/n_ptc_samples:.0f}%); "
         f"Q12 Korean equivalent: 17% PTC samples Hashimoto-like</sub>",
    font=dict(size=14, color=INK)),
    yaxis=dict(title="% cells with HLA-II z>1 (Hashimoto-like)", range=[0, 100], gridcolor="rgba(255,255,255,0.08)"),
    xaxis=dict(tickfont=dict(size=12)),
    showlegend=False, height=440, width=820, margin=dict(l=80,r=20,t=110,b=70), **DARK)
fig.write_html(FIG/"v17_lu2023_hashimoto_per_sample.html", include_plotlyjs="cdn", full_html=True)
log.info(f"wrote v17_lu2023_hashimoto_per_sample.html")

# Fig D: 2D scatter — DM_score vs HLA-II score per PTC cell (color by Hashimoto-like)
fig = go.Figure()
fig.add_trace(go.Scattergl(x=g_other, y=ptc[~ptc["hashi_like"]]["HLA_II_score"].values, mode="markers",
                           marker=dict(color="rgba(170,170,170,0.4)", size=2.5),
                           name=f"Other PTC (n={len(g_other):,})", hoverinfo="skip"))
fig.add_trace(go.Scattergl(x=g_hashi, y=ptc[ptc["hashi_like"]]["HLA_II_score"].values, mode="markers",
                           marker=dict(color="rgba(124,207,205,0.7)", size=3),
                           name=f"Hashimoto-like (n={len(g_hashi):,})", hoverinfo="skip"))
fig.add_hline(y=1.0, line=dict(color="white", dash="dash", width=1))
fig.update_layout(title=dict(
    text=f"<b>Q15 — Lu 2023 PTC cells: DM_score vs HLA-II score (per cell)</b><br>"
         f"<sub style='color:#7ccfcd'>{len(ptc):,} PTC malignant cells; Hashimoto-like (HLA-II z>1, teal) {'shifted toward higher DM_score' if np.median(g_hashi)>np.median(g_other) else 'shifted toward lower DM_score'}</sub>",
    font=dict(size=14, color=INK)),
    xaxis=dict(title="DM_score (panel mean z)", showgrid=False, zeroline=True, zerolinecolor="rgba(255,255,255,0.2)"),
    yaxis=dict(title="HLA Class II score (mean z)", showgrid=False, zeroline=True, zerolinecolor="rgba(255,255,255,0.2)"),
    height=560, width=900, margin=dict(l=80,r=20,t=110,b=60), **DARK,
    legend=dict(itemsizing="constant", bgcolor="rgba(0,0,0,0.3)"))
fig.write_html(FIG/"v17_lu2023_hashimoto_2d_scatter.html", include_plotlyjs="cdn", full_html=True)
log.info(f"wrote v17_lu2023_hashimoto_2d_scatter.html")

log.info("=== DONE ===")
print(json.dumps(summary, indent=2, default=str))
