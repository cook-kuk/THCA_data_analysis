#!/usr/bin/env python3
"""GSE213647 (Lee 2024 SNUBH+CNUH+KRIBB) — v17p35 8-gene panel transfer."""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
import plotly.graph_objects as go

OUT = Path("/opt/thyroid-dash/project/results/v17_korean")
FIG = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17")
OUT.mkdir(parents=True, exist_ok=True)
BG = "#0b0e12"; INK = "#F2F2F2"
DARK = dict(template="plotly_dark", paper_bgcolor=BG, plot_bgcolor=BG, font=dict(color=INK, size=14))

GENES = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]

# Load
clin = pd.read_csv("/data/thca/v17_korean/GSE213647/sample_sheet_clinical.tsv", sep="\t")
expr = pd.read_csv("/data/thca/v17_korean/GSE213647/expression_counts.tsv.gz", sep="\t", index_col=0)
print(f"clin {clin.shape}, expr {expr.shape}")

# Map ENSG → symbol via biomart-style approach: use shared_genes.txt + gencode mapping not directly available
# Alternative: load TCGA shared_genes which already has symbols, find via a quick gene-symbol lookup
# Use a built-in mapping from a TSV in the project (look for any ENSG-symbol table)
import subprocess
candidates = subprocess.run("find /opt/thyroid-dash/project -iname '*ensg*' -o -iname '*gene_id*symbol*' -o -iname '*biomart*' 2>/dev/null | head",
                            shell=True, capture_output=True, text=True).stdout.strip().split("\n")
print(f"ENSG mapping candidates: {candidates[:5]}")

# Strategy: use a small inline gene mapping (these 8 genes are well-known)
# ENSG IDs from GENCODE v40 (matches GSE213647 STAR ReadsPerGene)
ENSG_MAP = {
    "ENSG00000105641": "SLC5A5",
    "ENSG00000115705": "TPO",
    "ENSG00000042832": "TG",
    "ENSG00000165409": "TSHR",
    "ENSG00000125618": "PAX8",
    "ENSG00000136352": "NKX2-1",
    "ENSG00000178919": "FOXE1",
    "ENSG00000211448": "DIO1",
}
# Strip version suffix from expr index
expr.index = expr.index.str.split(".").str[0]
print(f"\nCheck panel gene presence in GSE213647:")
for ensg, sym in ENSG_MAP.items():
    print(f"  {ensg} ({sym}): {'✓' if ensg in expr.index else '✗ MISSING'}")

# Subset to 8 panel
mat = expr.loc[[e for e in ENSG_MAP if e in expr.index]].copy()
mat.index = [ENSG_MAP[e] for e in mat.index]
mat = mat.reindex(GENES).dropna()
print(f"\npanel matrix: {mat.shape} (genes × samples)")

# Library size normalization (CPM) + log
lib = expr.sum(axis=0)
cpm = mat.div(lib, axis=1) * 1e6
log_cpm = np.log2(cpm + 1)
print(f"log_cpm range: {log_cpm.values.min():.2f} – {log_cpm.values.max():.2f}")

# Z-score per gene across all samples
z = log_cpm.sub(log_cpm.mean(axis=1), axis=0).div(log_cpm.std(axis=1)+1e-9, axis=0)
panel_score = z.mean(axis=0)  # mean of 8-gene z-score per sample = "DM-score proxy"

# Merge with clinical
clin = clin.set_index("gsm")
clin["panel_z"] = panel_score
clin["lib_size"] = lib

# Define dedifferentiation order
order = ["Normal","PTC","PDFP","UTC/ATC"]
clin_ord = clin[clin["histology"].isin(order)].copy()
clin_ord["histology"] = pd.Categorical(clin_ord["histology"], categories=order, ordered=True)

# Save
clin_ord.to_csv(OUT/"GSE213647_panel_score.tsv", sep="\t")
print(f"\nsaved {OUT}/GSE213647_panel_score.tsv")

# Statistics
print("\n=== Panel z-score by histology ===")
print(clin_ord.groupby("histology", observed=True)["panel_z"].describe()[["count","mean","std","50%"]])

# Trend test
from scipy import stats
groups = [clin_ord[clin_ord["histology"]==h]["panel_z"].dropna().values for h in order]
H, p_kw = stats.kruskal(*groups)
print(f"\nKruskal-Wallis across 4 groups: H={H:.2f}, p={p_kw:.2e}")
# Pairwise
for i, a in enumerate(order):
    for j, b in enumerate(order):
        if i >= j: continue
        u, p = stats.mannwhitneyu(groups[i], groups[j], alternative="two-sided")
        print(f"  {a} vs {b}: U={u:.0f}, p={p:.2e}")

# Library kit batch effect check
print("\n=== Panel z-score by library kit (batch check) ===")
print(clin_ord.groupby("library_kit", observed=True)["panel_z"].describe()[["count","mean","std","50%"]])
u, p = stats.mannwhitneyu(
    clin_ord[clin_ord["library_kit"]=="stranded mRNA LT sample prep kit"]["panel_z"].dropna(),
    clin_ord[clin_ord["library_kit"]=="TruSeq RNA access library"]["panel_z"].dropna(),
    alternative="two-sided")
print(f"  library kit difference: U={u:.0f}, p={p:.2e}")

# ============================================================
# Figures
# ============================================================
# 1. Box+strip per histology, colored by lib_kit
fig = go.Figure()
colors = {"Normal":"#2ECC71","PTC":"#F5A623","PDFP":"#E67E22","UTC/ATC":"#c24c4c"}
for h in order:
    sub = clin_ord[clin_ord["histology"]==h]
    fig.add_trace(go.Box(y=sub["panel_z"], x=[h]*len(sub),
                         name=f"{h} (n={len(sub)})",
                         marker=dict(color=colors[h]),
                         boxpoints="all", jitter=0.4, pointpos=0,
                         marker_size=4, line=dict(width=2)))
fig.update_layout(title=dict(text=f"<b>GSE213647 (Lee 2024, 분당서울대 SNUBH+CNUH+KRIBB) — 8-gene panel z-score by histology</b><br>" +
                                  f"<sub style='color:#7ccfcd'>Normal {len(clin_ord[clin_ord.histology=='Normal'])} → PTC {len(clin_ord[clin_ord.histology=='PTC'])} → PDFP {len(clin_ord[clin_ord.histology=='PDFP'])} → ATC {len(clin_ord[clin_ord.histology=='UTC/ATC'])}, " +
                                  f"Kruskal-Wallis p={p_kw:.2e}</sub>",
                             font=dict(size=14, color=INK)),
                  yaxis=dict(title=dict(text="8-gene panel mean z-score", font=dict(size=13, color=INK)),
                             gridcolor="rgba(255,255,255,0.08)",
                             tickfont=dict(size=12, color=INK), zeroline=True, zerolinecolor="rgba(255,255,255,0.2)"),
                  xaxis=dict(tickfont=dict(size=12, color=INK)),
                  showlegend=False, height=480, margin=dict(l=80,r=20,t=110,b=60), **DARK)
fig.write_html(FIG/"v17_korean_panel_score.html", include_plotlyjs="cdn", full_html=True)
print(f"\nwrote {FIG}/v17_korean_panel_score.html")

# 2. Per-gene heatmap (samples sorted by panel_z, grouped by histology)
sample_order = clin_ord.sort_values(["histology","panel_z"]).index.tolist()
zmat = z[sample_order]
hist_colors = [colors[clin_ord.loc[s,"histology"]] for s in sample_order]
# annotation row: histology bar
fig = go.Figure(go.Heatmap(z=zmat.values, y=zmat.index, x=list(range(len(sample_order))),
                           colorscale="RdBu_r", zmid=0, zmin=-3, zmax=3,
                           colorbar=dict(title=dict(text="z", font=dict(size=13, color=INK)),
                                         tickfont=dict(size=11, color=INK))))
# vertical separators between histology groups
sep_x = []
prev = None; cnt = 0
for s in sample_order:
    h = clin_ord.loc[s,"histology"]
    if prev is not None and h != prev:
        sep_x.append(cnt - 0.5)
    prev = h; cnt += 1
for x in sep_x:
    fig.add_vline(x=x, line=dict(color="white", width=2))
# group labels at bottom
prev = None; start = 0
for i, s in enumerate(sample_order):
    h = clin_ord.loc[s,"histology"]
    if prev is not None and h != prev:
        fig.add_annotation(x=(start+i-1)/2, y=-1.2, xref="x", yref="y", showarrow=False,
                           text=f"<b>{prev}</b><br>n={i-start}",
                           font=dict(color=colors[prev], size=12))
        start = i
    prev = h
fig.add_annotation(x=(start+len(sample_order)-1)/2, y=-1.2, xref="x", yref="y", showarrow=False,
                   text=f"<b>{prev}</b><br>n={len(sample_order)-start}",
                   font=dict(color=colors[prev], size=12))
fig.update_layout(title=dict(text=f"<b>GSE213647 × 8-gene panel z-score heatmap (n={len(sample_order)} samples sorted by histology, panel score)</b>",
                             font=dict(size=14, color=INK)),
                  yaxis=dict(autorange="reversed",
                             tickfont=dict(size=13, family="JetBrains Mono", color=INK)),
                  xaxis=dict(showticklabels=False,
                             title=dict(text="samples (Normal → PTC → PDFP → ATC)",
                                        font=dict(size=12, color=INK))),
                  height=460, margin=dict(l=120, r=40, t=60, b=70), **DARK)
fig.write_html(FIG/"v17_korean_panel_heatmap.html", include_plotlyjs="cdn", full_html=True)
print(f"wrote {FIG}/v17_korean_panel_heatmap.html")

# 3. Compare to TCGA-THCA panel z-score (cross-cohort consistency)
import json
summary = {
    "cohort": "GSE213647 (Lee 2024)",
    "origin": "Korea (CNUH + SNUBH + KRIBB)",
    "n_total": int(len(clin_ord)),
    "n_by_histology": {h: int(len(clin_ord[clin_ord.histology==h])) for h in order},
    "panel_z_median_by_histology": {h: float(clin_ord[clin_ord.histology==h]["panel_z"].median()) for h in order},
    "kruskal_p": float(p_kw),
    "library_kit_pvalue": float(p),
    "panel_genes_found": [ENSG_MAP[e] for e in ENSG_MAP if e in expr.index],
    "panel_genes_missing": [ENSG_MAP[e] for e in ENSG_MAP if e not in expr.index],
    "caveat": "BRAF/RAS/TERT per-sample calls NOT in GEO; full DM1/DM2 validation requires EGA DAC for EGAS00001003540 (Yoo 2019) or SNUBH partnership.",
}
(OUT/"GSE213647_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
print(f"\nwrote {OUT}/GSE213647_summary.json")
print(json.dumps(summary, indent=2, ensure_ascii=False))
