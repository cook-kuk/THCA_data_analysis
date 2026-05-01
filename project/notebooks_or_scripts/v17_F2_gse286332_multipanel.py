#!/usr/bin/env python3
"""F2 — GSE286332 PTC vs PTC+HT multi-panel (Pillar 2)."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/figures/F2_gse286332_multipanel"
RES.mkdir(parents=True, exist_ok=True)

deg = pd.read_csv(PROJ / "results/p3_gse286332/deg_ptcht_vs_ptc.tsv", sep="\t")
gsea_h = pd.read_csv(PROJ / "results/p3_gse286332/gsea_MSigDB_Hallmark_2020.tsv", sep="\t")
panel = pd.read_csv(PROJ / "results/p3_gse286332/8gene_panel_per_sample.tsv", sep="\t", index_col=0)
hla = pd.read_csv(PROJ / "results/p3_gse286332/hla_module_scores.tsv", sep="\t", index_col=0)

fig = plt.figure(figsize=(15, 11))
gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.30, left=0.07, right=0.96, top=0.93, bottom=0.07)

# ---------- Panel A — Volcano ----------
axA = fig.add_subplot(gs[0, 0])
deg_clean = deg.dropna(subset=["padj", "log2FoldChange"]).copy()
neg_log_p = -np.log10(deg_clean["padj"].clip(lower=1e-300))
sig = (deg_clean["padj"] < 0.05) & (deg_clean["log2FoldChange"].abs() > 1.0)
axA.scatter(deg_clean.loc[~sig, "log2FoldChange"], neg_log_p[~sig], s=2, c="#cccccc", alpha=0.5)
axA.scatter(deg_clean.loc[sig & (deg_clean["log2FoldChange"] > 0), "log2FoldChange"],
             neg_log_p[sig & (deg_clean["log2FoldChange"] > 0)], s=4, c="#cc4444", alpha=0.8)
axA.scatter(deg_clean.loc[sig & (deg_clean["log2FoldChange"] < 0), "log2FoldChange"],
             neg_log_p[sig & (deg_clean["log2FoldChange"] < 0)], s=4, c="#3a4ea0", alpha=0.8)

# Label top genes
top_up = deg_clean.nsmallest(8, "padj").query("log2FoldChange > 0")
top_dn = deg_clean[deg_clean["log2FoldChange"] < 0].nsmallest(5, "padj")
for _, r in top_up.iterrows():
    axA.annotate(r["gene"], (r["log2FoldChange"], -np.log10(max(r["padj"], 1e-300))),
                  fontsize=8, color="#cc4444", fontweight="bold",
                  xytext=(3, 3), textcoords="offset points")
for _, r in top_dn.iterrows():
    axA.annotate(r["gene"], (r["log2FoldChange"], -np.log10(max(r["padj"], 1e-300))),
                  fontsize=8, color="#3a4ea0", fontweight="bold",
                  xytext=(3, 3), textcoords="offset points")
axA.axhline(-np.log10(0.05), color="black", linestyle=":", lw=0.8)
axA.axvline(1.0, color="black", linestyle=":", lw=0.5)
axA.axvline(-1.0, color="black", linestyle=":", lw=0.5)
axA.set_xlabel("log2 FC (PTC+HT vs PTC)", fontsize=11)
axA.set_ylabel(r"$-\log_{10}$ padj", fontsize=11)
axA.set_title(f"A — Differential expression (n=29,672)\n10,380 DEGs (padj<0.05): {(deg_clean['padj']<0.05).sum()} = 6,004 up + 4,376 down",
               fontsize=11, loc="left", fontweight="bold")
axA.spines["top"].set_visible(False); axA.spines["right"].set_visible(False)

# ---------- Panel B — GSEA bar ----------
axB = fig.add_subplot(gs[0, 1])
gsea_clean = gsea_h.dropna(subset=["NES"]).copy()
gsea_clean["FDR q-val"] = gsea_clean["FDR q-val"].astype(float)
top_up_gsea = gsea_clean[gsea_clean["NES"] > 0].nsmallest(8, "FDR q-val")
top_dn_gsea = gsea_clean[gsea_clean["NES"] < 0].nsmallest(4, "FDR q-val")
top_combined = pd.concat([top_up_gsea, top_dn_gsea]).sort_values("NES")

colors_b = ["#3a4ea0" if v < 0 else "#cc4444" for v in top_combined["NES"]]
y_b = np.arange(len(top_combined))
axB.barh(y_b, top_combined["NES"], color=colors_b, edgecolor="black", lw=0.4)
axB.set_yticks(y_b)
axB.set_yticklabels([str(t)[:38] for t in top_combined["Term"]], fontsize=9)
axB.axvline(0, color="black", lw=0.5)
axB.set_xlabel("Normalized Enrichment Score (NES)", fontsize=11)
axB.set_title("B — Hallmark GSEA (top 8 up + 4 down)\nIFN-γ FDR=2e-4, Allograft Rejection NES=+2.12", fontsize=11, loc="left", fontweight="bold")
axB.spines["top"].set_visible(False); axB.spines["right"].set_visible(False)
# annotate FDR
for i, (_, r) in enumerate(top_combined.iterrows()):
    fdr = r["FDR q-val"]
    fdr_str = "0" if fdr == 0 else f"{fdr:.0e}"
    x_off = 0.05 if r["NES"] > 0 else -0.05
    axB.text(r["NES"] + x_off, i, fdr_str, fontsize=8, va="center",
              ha="left" if r["NES"] > 0 else "right")

# ---------- Panel C — 8-gene per-gene boxplot ----------
axC = fig.add_subplot(gs[1, 0])
panel["group"] = ["PTC" if s.startswith("NG_") else "PTC+HT" for s in panel.index]
g8_genes = ['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1']
data_pos = []; data_pht = []
for g in g8_genes:
    if g in panel.columns:
        data_pos.append(panel.loc[panel["group"] == "PTC", g].values)
        data_pht.append(panel.loc[panel["group"] == "PTC+HT", g].values)
positions_ptc = np.arange(len(g8_genes)) - 0.18
positions_pht = np.arange(len(g8_genes)) + 0.18
bp1 = axC.boxplot(data_pos, positions=positions_ptc, widths=0.32, patch_artist=True,
                    boxprops=dict(facecolor="#3a4ea0", alpha=0.7))
bp2 = axC.boxplot(data_pht, positions=positions_pht, widths=0.32, patch_artist=True,
                    boxprops=dict(facecolor="#cc4444", alpha=0.7))
axC.set_xticks(np.arange(len(g8_genes))); axC.set_xticklabels(g8_genes, rotation=20, fontsize=9)
axC.set_ylabel("log2 FPKM", fontsize=11)
axC.set_title("C — 8-gene panel per-gene (PTC blue vs PTC+HT red)\nPAX8 d=−2.32, NKX2-1 d=−1.92, FOXE1 d=−1.75; SLC5A5 preserved", fontsize=11, loc="left", fontweight="bold")
axC.spines["top"].set_visible(False); axC.spines["right"].set_visible(False)
# legend
axC.plot([], [], "s", color="#3a4ea0", label="PTC (n=9)", markersize=10)
axC.plot([], [], "s", color="#cc4444", label="PTC+HT (n=9)", markersize=10)
axC.legend(loc="upper right", fontsize=9.5)

# ---------- Panel D — HLA-II heatmap per sample ----------
axD = fig.add_subplot(gs[1, 1])
samples = list(hla.index)
ng_idx = [s for s in samples if s.startswith("NG_")]
th_idx = [s for s in samples if s.startswith("TH_")]
ordered = ng_idx + th_idx

# Re-load HLA-II raw expression to make heat
data_path = "/data/thca/v17_korean/GSE286332/GSE286332_all_sample_rawdata.txt.gz"
df_raw = pd.read_csv(data_path, sep="\t", low_memory=False)
fpkm_cols = [c for c in df_raw.columns if c.endswith("_FPKM")]
fpkm = df_raw.groupby("Gene_Symbol", as_index=True)[fpkm_cols].sum()
fpkm.columns = [c.replace("_FPKM", "") for c in fpkm.columns]
log_fpkm = np.log2(fpkm + 1)

hla2_genes = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
               "HLA-DMA", "HLA-DMB", "CIITA", "HLA-DOB"]
hla2_in = [g for g in hla2_genes if g in log_fpkm.index]
heat_d = log_fpkm.loc[hla2_in, ordered]
heat_z = heat_d.apply(lambda c: (c - c.mean()) / (c.std(ddof=1) + 1e-9), axis=1)
im = axD.imshow(heat_z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
axD.set_yticks(np.arange(len(hla2_in))); axD.set_yticklabels(hla2_in, fontsize=9)
axD.set_xticks(np.arange(len(ordered))); axD.set_xticklabels(ordered, rotation=90, fontsize=8)
axD.axvline(8.5, color="black", lw=2)
plt.colorbar(im, ax=axD, fraction=0.04, pad=0.04, label="Z-score")
axD.set_title("D — HLA-II module per sample (n=18)\nCohen d=+3.65, MW p=4e-4 (★ exceptional)", fontsize=11, loc="left", fontweight="bold")

fig.suptitle("Figure 2 — GSE286332 PTC vs PTC+Hashimoto's molecular dissection (Pillar 2)", fontsize=14, fontweight="bold", y=0.99)
plt.savefig(RES / "F2_gse286332_multipanel.pdf")
plt.savefig(RES / "F2_gse286332_multipanel.png", dpi=200)
plt.close()
print(f"✓ F2 saved to {RES}")
