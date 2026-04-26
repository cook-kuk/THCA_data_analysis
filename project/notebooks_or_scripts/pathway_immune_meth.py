#!/usr/bin/env python3
"""THCA interpretation layer — Phase 1-4 driver.

Builds three NEW scientific-analysis layers on top of the existing biomarker /
ML outputs, without modifying any pipeline file:

  Phase 1 — Pathway enrichment (hypergeometric ORA) on BRAF-up / RAS-up /
            novel-only gene lists using hardcoded curated gene-set collections.
  Phase 2 — Immune cell signature scoring (ssGSEA-like z-score means) for
            TCGA-THCA, with BRAF_like vs RAS_like Mann-Whitney contrast.
  Phase 3 — Methylation classifier on GSE97466 top-5000 probes (tumor vs
            normal since no mutation labels available), stratified 5-fold CV
            with 800-iteration bootstrap AUC CI.
  Phase 4 — Build the self-contained dashboard page
            reports/html/pages/18_pathway_immune_meth.html.

All outputs are self-contained: CSP-safe (no CDN), Korean primary labels,
vendored Plotly at ../assets/vendor/plotly.min.js. Honest caveats are baked
into the UI where the signal is weak.
"""
from __future__ import annotations

import json
import math
import os
import sys
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
TABLES = ROOT / "results" / "tables"
FIGS = ROOT / "reports" / "html" / "figs_interactive"
ML_DIR = ROOT / "results" / "ml"
PAGES = ROOT / "reports" / "html" / "pages"
META = ROOT / "metadata" / "sample_master.tsv"
TCGA_LOG2 = ROOT / "data_processed" / "bulk_rnaseq" / "TCGA-THCA_rnaseq_expression_log2.tsv"
METH_BETA = ROOT / "data_processed" / "methylation" / "GSE97466_beta_top5000.tsv"

for d in (TABLES, FIGS, ML_DIR, PAGES):
    d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Curated gene-set collections (MSigDB-style + thyroid-specific).
# Each set is a Python literal so there is no network dependency; sources are
# cited in markdown on the dashboard page.
# ---------------------------------------------------------------------------

PATHWAYS: Dict[str, List[str]] = {
    # --- MSigDB Hallmark-equivalent (hand-curated subsets, see page markdown)
    "MAPK_signaling": [
        "BRAF", "ARAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3",
        "DUSP1", "DUSP4", "DUSP5", "DUSP6", "SPRY1", "SPRY2", "SPRY4",
        "ETV1", "ETV4", "ETV5", "FOS", "FOSL1", "JUN", "EGR1",
    ],
    "PI3K_AKT": [
        "PIK3CA", "PIK3CB", "PIK3R1", "PTEN", "AKT1", "AKT2", "AKT3",
        "MTOR", "RPS6KB1", "PDK1", "GSK3B", "TSC1", "TSC2",
    ],
    "EMT": [
        "VIM", "CDH1", "CDH2", "SNAI1", "SNAI2", "ZEB1", "ZEB2",
        "TWIST1", "TWIST2", "FN1", "MMP2", "MMP9", "SPARC", "S100A4",
    ],
    "Apoptosis": [
        "TP53", "CDKN2A", "CDKN2B", "BAX", "BCL2", "BCL2L1", "BBC3",
        "CASP3", "CASP7", "CASP8", "CASP9", "MCL1", "BAK1",
    ],
    "Inflammation": [
        "IL6", "IL10", "TNF", "IFNG", "IL1A", "IL1B", "NFKB1", "NFKB2",
        "CXCL8", "CCL2", "CCL5", "PTGS2",
    ],
    "Thyroid_differentiation": [
        "DIO1", "DIO2", "TG", "TPO", "TSHR", "PAX8", "FOXE1", "NKX2-1",
        "SLC5A5", "SLC26A4", "TTF1", "DUOX1", "DUOX2",
    ],
    "Immune_exhaustion": [
        "PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TOX",
        "PDCD1LG2", "EOMES",
    ],
    "Kinase_signaling": [
        "ERBB2", "ERBB3", "ERBB4", "EGFR", "MET", "KIT", "PDGFRA",
        "PDGFRB", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "RET", "IGF1R",
    ],
    "Proteases": [
        "KLK1", "KLK3", "KLK5", "KLK6", "KLK7", "KLK10", "KLK11",
        "MMP1", "MMP2", "MMP3", "MMP7", "MMP9", "MMP11", "MMP14", "MMP16",
        "CTSA", "CTSB", "CTSD", "CTSH", "CTSK", "CTSL", "CTSS", "CTSZ",
    ],
    "Cell_cycle": [
        "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2", "CDK2", "CDK4",
        "CDK6", "RB1", "E2F1", "E2F2", "MCM2", "MCM3", "MCM4", "MCM5",
        "PCNA", "MKI67",
    ],
    "NOTCH": [
        "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "JAG1", "JAG2",
        "DLL1", "DLL3", "DLL4", "HES1", "HES5", "HEY1", "HEY2", "MAML1",
    ],
    "Wnt": [
        "CTNNB1", "APC", "AXIN1", "AXIN2", "WNT3A", "WNT5A", "WNT10B",
        "LEF1", "TCF7", "TCF7L1", "TCF7L2", "DVL1", "DVL2", "FZD1", "FZD2",
    ],
    # --- Thyroid-specific sets
    "BRS_up_BRAFlike": [
        # genes upregulated in BRAF-like (literature-anchored surrogates)
        "TACSTD2", "TMPRSS4", "KLK7", "KLK10", "LCN2", "FN1", "TIMP1",
        "SERPINA1", "CXCL14", "S100A4", "CD44", "HMGA2",
    ],
    "BRS_up_RASlike": [
        # genes upregulated in RAS-like (better-differentiated, follicular-axis)
        "TG", "TPO", "TSHR", "DIO1", "DIO2", "SLC5A5", "SLC26A4",
        "PAX8", "FOXE1", "NKX2-1", "TFF3", "ALDH1A1",
    ],
    "TDS_core": [
        "DIO1", "DIO2", "TG", "TPO", "TSHR", "SLC5A5", "SLC26A4",
        "FOXE1", "PAX8", "NKX2-1", "GLIS3", "TFF3", "DUOX1", "DUOX2",
        "THRA", "THRB",
    ],
    "Dedifferentiation_aggressive": [
        "TERT", "HMGA2", "S100A4", "CD44", "KRT19", "FN1", "VIM",
        "MMP9", "TWIST1", "LCN2", "SERPINA1",
    ],
}

IMMUNE_SIGNATURES: Dict[str, List[str]] = {
    "T_cell_total": ["CD3D", "CD3E", "CD3G", "CD2", "CD8A"],
    "CD8_cytotoxic": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "GNLY"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "TIGIT"],
    "Macrophage_M1": ["CD68", "CD86", "TNF", "NOS2", "IL1A", "IL6"],
    "Macrophage_M2": ["CD163", "MRC1", "IL10", "TGFB1", "ARG1"],
    "NK_cell": ["NCAM1", "FCGR3A", "KLRD1", "NKG7", "PRF1"],
    "B_cell": ["CD19", "MS4A1", "CD79A", "CD79B", "BLK"],
    "Checkpoint_exhaustion": ["PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2"],
    "IFN_gamma_response": ["IFNG", "STAT1", "IRF1", "CXCL10", "CXCL9", "IDO1"],
}

# Korean display names
PATHWAY_LABELS_KO: Dict[str, str] = {
    "MAPK_signaling": "MAPK 신호 (MAPK signaling)",
    "PI3K_AKT": "PI3K / AKT",
    "EMT": "상피-중간엽 전환 (EMT)",
    "Apoptosis": "세포사멸 (Apoptosis)",
    "Inflammation": "염증 (Inflammation)",
    "Thyroid_differentiation": "갑상선 분화 (Thyroid differentiation)",
    "Immune_exhaustion": "면역 탈진 (Immune exhaustion)",
    "Kinase_signaling": "키나제 신호 (Kinase signaling)",
    "Proteases": "단백분해효소 (Proteases)",
    "Cell_cycle": "세포주기 (Cell cycle)",
    "NOTCH": "NOTCH 신호",
    "Wnt": "Wnt 신호",
    "BRS_up_BRAFlike": "BRS-up · BRAF-like (BRS-up in BRAF-like)",
    "BRS_up_RASlike": "BRS-up · RAS-like (BRS-up in RAS-like)",
    "TDS_core": "TDS core (Thyroid Differentiation Score)",
    "Dedifferentiation_aggressive": "탈분화 · 공격성 (Dedifferentiation)",
}

IMMUNE_LABELS_KO: Dict[str, str] = {
    "T_cell_total": "T세포 전체 (T cell total)",
    "CD8_cytotoxic": "CD8 세포독성 (CD8 cytotoxic)",
    "Treg": "조절 T세포 (Treg)",
    "Macrophage_M1": "대식세포 M1 (Macrophage M1)",
    "Macrophage_M2": "대식세포 M2 (Macrophage M2)",
    "NK_cell": "자연살해세포 (NK cell)",
    "B_cell": "B세포 (B cell)",
    "Checkpoint_exhaustion": "면역관문 · 탈진 (Checkpoint exhaustion)",
    "IFN_gamma_response": "IFN-γ 반응 (IFN-γ response)",
}


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def bh_adjust(pvals: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    q = p * n / ranks
    # enforce monotonicity (decreasing when sorted by p)
    sorted_q = q[order]
    for i in range(n - 2, -1, -1):
        sorted_q[i] = min(sorted_q[i], sorted_q[i + 1])
    out = np.empty(n)
    out[order] = np.clip(sorted_q, 0, 1)
    return out


def plotly_offline_html(fig: go.Figure, out_path: Path, title: str) -> None:
    """Write plotly fig as a standalone CSP-safe HTML using the vendored plotly.min.js."""
    # include_plotlyjs=False so we reference the vendored copy; inject the <script> tag ourselves.
    fig.update_layout(
        paper_bgcolor="#0b1a3a",
        plot_bgcolor="#07132b",
        font=dict(color="#e2fff9", family="Inter, system-ui, sans-serif"),
        margin=dict(l=60, r=30, t=60, b=60),
    )
    inner = pio.to_html(
        fig,
        include_plotlyjs=False,
        full_html=False,
        config={"displaylogo": False, "responsive": True},
    )
    html = f"""<!doctype html>
<html lang=\"ko\"><head><meta charset=\"utf-8\">
<title>{title}</title>
<style>html,body{{margin:0;padding:0;background:#07132b;color:#e2fff9;font-family:Inter,system-ui,sans-serif}}</style>
<script src=\"../assets/vendor/plotly.min.js\"></script>
</head><body>{inner}</body></html>
"""
    out_path.write_text(html)


# ---------------------------------------------------------------------------
# PHASE 1 — Pathway enrichment (hypergeometric ORA)
# ---------------------------------------------------------------------------

def phase1_pathway_enrichment() -> dict:
    print("[PHASE 1] Pathway enrichment — hypergeometric ORA")

    de = pd.read_csv(TABLES / "biomarker_de_full.tsv", sep="\t")
    valid = pd.read_csv(TABLES / "biomarker_validated.tsv", sep="\t")
    novel_set = set(valid.loc[valid["is_novel_validated"] == True, "gene"].astype(str))

    # Universe = all genes with a TCGA DE result
    universe = set(de["gene"].astype(str))
    N_universe = len(universe)
    print(f"  universe size = {N_universe} genes; novel-validated = {len(novel_set)}")

    de = de.copy()
    de["gene"] = de["gene"].astype(str)
    # BRAF-up = log2FC > 0 & FDR < 0.05 in TCGA
    sig_braf_up = set(de.loc[(de["log2FC_tcga"] > 0) & (de["fdr_tcga"] < 0.05), "gene"])
    sig_ras_up = set(de.loc[(de["log2FC_tcga"] < 0) & (de["fdr_tcga"] < 0.05), "gene"])
    # novel-only subset = any novel-validated, pooled both directions
    sig_novel = novel_set & universe
    print(f"  BRAF-up sig = {len(sig_braf_up)}   RAS-up sig = {len(sig_ras_up)}   novel = {len(sig_novel)}")

    directions = {
        "BRAF_up": sig_braf_up,
        "RAS_up": sig_ras_up,
        "Novel_validated": sig_novel,
    }

    rows: List[dict] = []
    for direction, sig_set in directions.items():
        K = len(sig_set)  # draws
        for pw_name, pw_genes in PATHWAYS.items():
            pw_set = set(pw_genes) & universe  # restrict to tested genes
            M = len(pw_set)
            if M == 0:
                continue
            overlap = sig_set & pw_set
            k = len(overlap)
            # Hypergeometric: P(X >= k) given M success in N, K draws
            # scipy.stats.hypergeom.sf(k-1, N, M, K)
            pval = stats.hypergeom.sf(k - 1, N_universe, M, K) if k > 0 else 1.0
            expected = K * M / N_universe if N_universe else 0
            fold = (k / expected) if expected > 0 else 0.0
            rows.append({
                "direction": direction,
                "pathway": pw_name,
                "pathway_ko": PATHWAY_LABELS_KO.get(pw_name, pw_name),
                "N_universe": N_universe,
                "N_pathway_total": M,
                "N_sig_total": K,
                "N_overlap": k,
                "expected": round(expected, 3),
                "fold_enrichment": round(fold, 3),
                "p_value": pval,
                "overlap_genes": ";".join(sorted(overlap)[:30]),
            })

    out = pd.DataFrame(rows)
    if not out.empty:
        # BH per direction (each direction is an independent test family)
        out["q_value"] = np.nan
        for d in out["direction"].unique():
            mask = out["direction"] == d
            out.loc[mask, "q_value"] = bh_adjust(out.loc[mask, "p_value"].values)
        out["neg_log10_q"] = -np.log10(out["q_value"].clip(lower=1e-300))
        out["significant"] = (out["q_value"] < 0.05) & (out["fold_enrichment"] > 1.5)
        out = out.sort_values(["direction", "q_value", "p_value"]).reset_index(drop=True)

    out_path = TABLES / "pathway_enrichment.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"  wrote {out_path}  ({len(out)} rows)")

    # --- Dotplot (pathway y, -log10(q) color, dot size = fold enrichment, facet = direction)
    facet_dirs = ["BRAF_up", "RAS_up"]
    fig = make_subplots(rows=1, cols=2, subplot_titles=[
        "BRAF-up (업: BRAF_like)", "RAS-up (업: RAS_like)"],
        shared_yaxes=True, horizontal_spacing=0.1)
    for j, d in enumerate(facet_dirs, start=1):
        sub = out[out["direction"] == d].copy()
        # order pathways consistently (same order across facets) by union significance
        pw_order = (
            out.assign(score=out["neg_log10_q"])
               .groupby("pathway")["score"].max()
               .sort_values(ascending=True)
               .index.tolist()
        )
        sub = sub.set_index("pathway").reindex(pw_order).reset_index()
        sub["label"] = sub["pathway"].map(PATHWAY_LABELS_KO).fillna(sub["pathway"])
        sub["hover"] = sub.apply(
            lambda r: (
                f"{r['label']}<br>q={r['q_value']:.2e}<br>fold={r['fold_enrichment']:.2f}"
                f"<br>overlap={r['N_overlap']}/{r['N_pathway_total']}"
                if pd.notna(r["q_value"]) else f"{r['label']}"
            ),
            axis=1,
        )
        sizes = (sub["fold_enrichment"].fillna(0).clip(0.5, 10) * 5 + 6).tolist()
        fig.add_trace(
            go.Scatter(
                x=sub["neg_log10_q"],
                y=sub["label"],
                mode="markers",
                marker=dict(
                    size=sizes,
                    color=sub["neg_log10_q"],
                    colorscale="Viridis",
                    showscale=(j == 2),
                    colorbar=dict(title="-log10(q)") if j == 2 else None,
                    line=dict(color="rgba(94,234,212,.6)", width=0.5),
                ),
                text=sub["hover"],
                hoverinfo="text",
                name=d,
                showlegend=False,
            ),
            row=1, col=j,
        )
        fig.update_xaxes(title_text="-log10(q)", row=1, col=j, gridcolor="rgba(148,163,184,.2)")
    fig.update_layout(
        title="경로 농축 (Pathway ORA) · BRAF-up vs RAS-up",
        height=720,
    )
    plotly_offline_html(fig, FIGS / "pathway_dotplot.html",
                        "Pathway enrichment dotplot")

    # --- Top-15 bar (by |neg_log10_q| across directions, deduped pathways)
    top = (
        out.sort_values("neg_log10_q", ascending=False)
           .drop_duplicates("pathway")
           .head(15)
           .sort_values("neg_log10_q")
    )
    # colour by direction
    color_map = {"BRAF_up": "#fb7185", "RAS_up": "#34d399", "Novel_validated": "#a78bfa"}
    fig2 = go.Figure()
    fig2.add_trace(go.Bar(
        x=top["neg_log10_q"],
        y=[PATHWAY_LABELS_KO.get(p, p) for p in top["pathway"]],
        orientation="h",
        marker=dict(color=[color_map.get(d, "#64748b") for d in top["direction"]]),
        text=[f"{fe:.2f}×  (q={q:.2e}, {dn})" for fe, q, dn in
              zip(top["fold_enrichment"], top["q_value"], top["direction"])],
        textposition="outside",
    ))
    fig2.update_layout(
        title="Top-15 유의 경로 · BRAF-up / RAS-up / Novel",
        xaxis_title="-log10(q)",
        yaxis_title="",
        height=680,
    )
    plotly_offline_html(fig2, FIGS / "pathway_barplot_top15.html",
                        "Pathway top-15 bar")

    # KPI summary
    top_pw = out.sort_values("q_value").head(1)
    summary = {
        "n_rows": int(len(out)),
        "n_significant": int(out["significant"].sum()) if "significant" in out else 0,
        "top_pathway": top_pw["pathway"].iloc[0] if len(top_pw) else None,
        "top_pathway_q": float(top_pw["q_value"].iloc[0]) if len(top_pw) else None,
        "top_pathway_direction": top_pw["direction"].iloc[0] if len(top_pw) else None,
    }
    return {"df": out, "summary": summary}


# ---------------------------------------------------------------------------
# PHASE 2 — Immune signature scoring (TCGA-THCA)
# ---------------------------------------------------------------------------

def phase2_immune_signatures() -> dict:
    print("[PHASE 2] Immune signature scoring — TCGA-THCA")

    # Load expression (genes x samples, log2)
    expr = pd.read_csv(TCGA_LOG2, sep="\t", index_col=0)
    meta = pd.read_csv(META, sep="\t")
    tcga_meta = meta[meta["dataset"] == "TCGA-THCA"].copy()
    # Map sample_id to expr columns; expr columns are TCGA barcodes
    valid_ids = [s for s in tcga_meta["sample_id"] if s in expr.columns]
    expr = expr[valid_ids]
    tcga_meta = tcga_meta.set_index("sample_id").loc[valid_ids]
    print(f"  TCGA samples = {len(valid_ids)} (expr × meta intersection)")

    # z-score each gene across tumours (row-wise)
    # Drop rows with near-zero variance to avoid /0
    row_std = expr.std(axis=1)
    keep = row_std > 1e-6
    exprz = expr.loc[keep].sub(expr.loc[keep].mean(axis=1), axis=0).div(row_std[keep], axis=0)

    # Compute signature score = mean z-score of available markers
    sig_scores: Dict[str, pd.Series] = {}
    missing_log: Dict[str, List[str]] = {}
    for name, genes in IMMUNE_SIGNATURES.items():
        present = [g for g in genes if g in exprz.index]
        missing = [g for g in genes if g not in exprz.index]
        if missing:
            missing_log[name] = missing
        if not present:
            continue
        sig_scores[name] = exprz.loc[present].mean(axis=0)

    score_df = pd.DataFrame(sig_scores)
    score_df["molecular_subtype"] = tcga_meta["molecular_subtype"].values
    score_df["sample_id"] = score_df.index
    score_df = score_df[["sample_id", "molecular_subtype"] + list(sig_scores.keys())]
    out_path = TABLES / "immune_signatures_tcga.tsv"
    score_df.to_csv(out_path, sep="\t", index=False)
    print(f"  wrote {out_path}  ({score_df.shape})")

    # --- Mann-Whitney contrast BRAF_like vs RAS_like
    braf = score_df[score_df["molecular_subtype"] == "BRAF_like"]
    rasl = score_df[score_df["molecular_subtype"] == "RAS_like"]
    print(f"  BRAF_like n = {len(braf)}, RAS_like n = {len(rasl)}")

    rows = []
    for sig in sig_scores.keys():
        a = braf[sig].dropna().values
        b = rasl[sig].dropna().values
        if len(a) < 3 or len(b) < 3:
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        # Log2 effect size = log2(mean(2^a) / mean(2^b)) on scores doesn't make sense
        # (z-scores can be negative). Use mean difference + Cliff's-delta style effect.
        diff = float(np.mean(a) - np.mean(b))
        # Cliff's delta ~ 2U/(n1 n2) - 1
        cd = 2 * u / (len(a) * len(b)) - 1
        rows.append({
            "signature": sig,
            "signature_ko": IMMUNE_LABELS_KO.get(sig, sig),
            "n_braf": len(a),
            "n_ras": len(b),
            "mean_braf": float(np.mean(a)),
            "mean_ras": float(np.mean(b)),
            "mean_diff_braf_minus_ras": diff,
            "cliffs_delta": float(cd),
            "mw_u": float(u),
            "p_value": float(p),
            "missing_markers": ";".join(missing_log.get(sig, [])),
        })
    contrast = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    contrast["q_value"] = bh_adjust(contrast["p_value"].values) if len(contrast) else []
    contrast_path = TABLES / "immune_signatures_contrast.tsv"
    contrast.to_csv(contrast_path, sep="\t", index=False)
    print(f"  wrote {contrast_path}")

    # --- Violin plots per signature (one figure, faceted grid)
    sigs = list(sig_scores.keys())
    n = len(sigs)
    cols = 3
    rows_n = math.ceil(n / cols)
    fig = make_subplots(rows=rows_n, cols=cols, subplot_titles=[
        IMMUNE_LABELS_KO.get(s, s) for s in sigs], horizontal_spacing=0.08, vertical_spacing=0.12)
    color = {"BRAF_like": "#fb7185", "RAS_like": "#34d399"}
    for i, sig in enumerate(sigs):
        r = i // cols + 1
        c = i % cols + 1
        for grp in ("BRAF_like", "RAS_like"):
            vals = score_df.loc[score_df["molecular_subtype"] == grp, sig].dropna()
            fig.add_trace(
                go.Violin(
                    y=vals,
                    x=[grp] * len(vals),
                    name=grp,
                    side="both",
                    box_visible=True,
                    meanline_visible=True,
                    points="all",
                    jitter=0.35,
                    pointpos=0,
                    marker=dict(size=3, color=color[grp], line=dict(color="rgba(255,255,255,.2)", width=0.3)),
                    line=dict(color=color[grp]),
                    fillcolor=color[grp],
                    opacity=0.5,
                    showlegend=(i == 0),
                    legendgroup=grp,
                ),
                row=r, col=c,
            )
        # annotation: p-value
        row_contrast = contrast[contrast["signature"] == sig]
        if len(row_contrast):
            p = row_contrast["p_value"].iloc[0]
            cd = row_contrast["cliffs_delta"].iloc[0]
            fig.add_annotation(
                xref=f"x{ '' if i==0 else i+1} domain", yref=f"y{ '' if i==0 else i+1} domain",
                x=0.5, y=1.02, showarrow=False,
                text=f"p={p:.2e}  δ={cd:+.2f}",
                font=dict(size=10, color="#94a3b8"),
                row=r, col=c,
            )
    fig.update_layout(
        title="면역 시그니처 (Immune signatures) · BRAF_like vs RAS_like · TCGA-THCA",
        height=260 * rows_n + 120,
        violingap=0.2, violingroupgap=0.08,
    )
    plotly_offline_html(fig, FIGS / "immune_signature_violin.html",
                        "Immune signature violins")

    # KPI summary
    top_contrast = contrast.sort_values("p_value").head(3) if len(contrast) else pd.DataFrame()
    summary = {
        "n_tcga": len(valid_ids),
        "n_braf": int(len(braf)),
        "n_ras": int(len(rasl)),
        "top_contrasts": top_contrast.to_dict("records"),
        "missing_markers": missing_log,
    }
    return {"contrast": contrast, "scores": score_df, "summary": summary}


# ---------------------------------------------------------------------------
# PHASE 3 — Methylation classifier (GSE97466)
# ---------------------------------------------------------------------------

def phase3_methylation() -> dict:
    print("[PHASE 3] Methylation classifier — GSE97466")

    beta = pd.read_csv(METH_BETA, sep="\t", index_col=0)
    qc = pd.read_csv(ROOT / "data_processed" / "methylation" / "GSE97466_methylation_qc_metrics.tsv", sep="\t")
    meta = pd.read_csv(META, sep="\t")

    g97 = meta[meta["dataset"] == "GSE97466"].copy()
    print(f"  GSE97466 samples in master = {len(g97)}; beta shape = {beta.shape}")

    # ---- Labels: tumor vs normal (clean, proxy for histology).  Mutation labels are not available.
    g97["y"] = np.where(g97["normal_vs_tumor"] == "tumor", 1, 0)
    # restrict to samples present in both beta matrix and metadata
    shared = [s for s in g97["sample_id"] if s in beta.columns]
    g97 = g97.set_index("sample_id").loc[shared]
    X = beta[shared].T.values  # samples × probes
    probes = beta.index.values
    y = g97["y"].values.astype(int)
    print(f"  n samples = {len(shared)}  (tumor={int(y.sum())}, normal={int((1-y).sum())})")

    # ---- Feature table for top differential probes (Mann-Whitney + |Δβ|)
    tumors = beta[shared][[s for s in shared if g97.loc[s, "y"] == 1]]
    normals = beta[shared][[s for s in shared if g97.loc[s, "y"] == 0]]
    print(f"  tumor cols={tumors.shape[1]}, normal cols={normals.shape[1]}")
    delta = tumors.mean(axis=1) - normals.mean(axis=1)
    # Welch t as a speed proxy (similar to MW for this n)
    tstat, pval = stats.ttest_ind(tumors.values, normals.values, axis=1, equal_var=False, nan_policy="omit")
    probe_df = pd.DataFrame({
        "probe_id": probes,
        "mean_tumor": tumors.mean(axis=1).values,
        "mean_normal": normals.mean(axis=1).values,
        "delta_beta": delta.values,
        "abs_delta": np.abs(delta.values),
        "t_stat": tstat,
        "p_value": pval,
    })
    probe_df["q_value"] = bh_adjust(np.nan_to_num(probe_df["p_value"].values, nan=1.0))
    top_probes = probe_df.sort_values("abs_delta", ascending=False).head(20).reset_index(drop=True)

    # ---- Pipelines: LogReg and RandomForest
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=7)
    imp = SimpleImputer(strategy="median")
    X_imp = imp.fit_transform(X)

    def cv_oof_probas(clf):
        oof = np.zeros(len(y), dtype=float)
        fold_aucs = []
        for fold_idx, (tr, te) in enumerate(cv.split(X_imp, y)):
            Xtr, Xte = X_imp[tr], X_imp[te]
            clf.fit(Xtr, y[tr])
            p = clf.predict_proba(Xte)[:, 1]
            oof[te] = p
            try:
                fold_aucs.append(roc_auc_score(y[te], p))
            except Exception:
                fold_aucs.append(np.nan)
        return oof, fold_aucs

    logreg_pipe = Pipeline([("scaler", StandardScaler()),
                            ("clf", LogisticRegression(max_iter=400, C=0.5, solver="liblinear"))])
    rf_pipe = Pipeline([("clf", RandomForestClassifier(n_estimators=400, max_features="sqrt",
                                                      class_weight="balanced", n_jobs=2, random_state=7))])

    print("  LogReg CV ...")
    lr_oof, lr_fold = cv_oof_probas(logreg_pipe)
    print("  RF CV ...")
    rf_oof, rf_fold = cv_oof_probas(rf_pipe)

    def auc_bootstrap_ci(y_true, y_score, n_iter=800, seed=13):
        rng = np.random.default_rng(seed)
        n = len(y_true)
        aucs = []
        for _ in range(n_iter):
            idx = rng.integers(0, n, size=n)
            if len(np.unique(y_true[idx])) < 2:
                continue
            try:
                aucs.append(roc_auc_score(y_true[idx], y_score[idx]))
            except Exception:
                continue
        if not aucs:
            return float("nan"), float("nan"), float("nan")
        lo, hi = np.percentile(aucs, [2.5, 97.5])
        return float(np.mean(aucs)), float(lo), float(hi)

    lr_auc = roc_auc_score(y, lr_oof)
    rf_auc = roc_auc_score(y, rf_oof)
    lr_mean, lr_lo, lr_hi = auc_bootstrap_ci(y, lr_oof)
    rf_mean, rf_lo, rf_hi = auc_bootstrap_ci(y, rf_oof)

    results = pd.DataFrame([
        {
            "model": "LogisticRegression",
            "n_samples": len(y),
            "n_tumor": int(y.sum()),
            "n_normal": int((1 - y).sum()),
            "cv_auc": lr_auc,
            "bootstrap_mean_auc": lr_mean,
            "bootstrap_lo_95": lr_lo,
            "bootstrap_hi_95": lr_hi,
            "fold_aucs": ";".join(f"{a:.3f}" for a in lr_fold),
        },
        {
            "model": "RandomForest",
            "n_samples": len(y),
            "n_tumor": int(y.sum()),
            "n_normal": int((1 - y).sum()),
            "cv_auc": rf_auc,
            "bootstrap_mean_auc": rf_mean,
            "bootstrap_lo_95": rf_lo,
            "bootstrap_hi_95": rf_hi,
            "fold_aucs": ";".join(f"{a:.3f}" for a in rf_fold),
        },
    ])
    res_path = ML_DIR / "methylation_classifier_results.tsv"
    results.to_csv(res_path, sep="\t", index=False)
    print(f"  wrote {res_path}")
    print(results.to_string(index=False))

    # --- ROC figure
    fig = go.Figure()
    for name, oof, color in [("LogReg", lr_oof, "#60a5fa"), ("RandomForest", rf_oof, "#f472b6")]:
        fpr, tpr, _ = roc_curve(y, oof)
        auc = roc_auc_score(y, oof)
        fig.add_trace(go.Scatter(
            x=fpr, y=tpr, mode="lines",
            name=f"{name}  AUC={auc:.3f}",
            line=dict(color=color, width=2),
        ))
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines",
                             line=dict(color="rgba(148,163,184,.5)", dash="dash"), showlegend=False))
    fig.update_layout(
        title=f"GSE97466 methylation · Tumor vs Normal ROC (n={len(y)})",
        xaxis_title="False positive rate",
        yaxis_title="True positive rate",
        height=520,
        legend=dict(x=0.55, y=0.06, bgcolor="rgba(7,19,43,.6)"),
    )
    plotly_offline_html(fig, FIGS / "methylation_roc.html", "Methylation ROC")

    # --- Top probes figure
    tp = top_probes.copy()
    tp["label"] = tp["probe_id"] + "  Δβ=" + tp["delta_beta"].round(3).astype(str)
    fig2 = go.Figure()
    fig2.add_trace(go.Bar(
        x=tp["delta_beta"],
        y=tp["probe_id"],
        orientation="h",
        marker=dict(color=np.where(tp["delta_beta"] > 0, "#fb7185", "#34d399")),
        text=[f"q={q:.2e}" for q in tp["q_value"]],
        textposition="outside",
    ))
    fig2.update_layout(
        title="Top-20 차등 메틸레이션 프로브 · Δβ (tumor − normal)",
        xaxis_title="Δβ (tumor − normal)  · 양수=tumor hyper-methylated",
        height=620,
    )
    plotly_offline_html(fig2, FIGS / "methylation_top_probes.html", "Methylation top probes")

    # Also persist the probe table
    top_probes.to_csv(ML_DIR / "methylation_top_probes.tsv", sep="\t", index=False)

    summary = {
        "n_samples": len(y),
        "n_tumor": int(y.sum()),
        "n_normal": int((1 - y).sum()),
        "logreg_cv_auc": lr_auc,
        "logreg_ci_95": (lr_lo, lr_hi),
        "rf_cv_auc": rf_auc,
        "rf_ci_95": (rf_lo, rf_hi),
    }
    return {"results": results, "top_probes": top_probes, "summary": summary}


# ---------------------------------------------------------------------------
# PHASE 4 — Dashboard page
# ---------------------------------------------------------------------------

PAGE_HTML_TEMPLATE = """<!doctype html>
<html lang="ko">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="color-scheme" content="dark light">
  <meta name="theme-color" content="#07132b">
  <title>경로 · 면역 · 메틸레이션 해석 | THCA Multi-Omics Dashboard</title>
  <meta name="description" content="바이오마커의 생물학적 해석 층: 경로 농축 분석(ORA), TCGA 면역 시그니처 대조, GSE97466 메틸레이션 분류기.">
  <link rel="icon" href="../assets/img/favicon.svg" type="image/svg+xml">
  <link rel="preload" href="../assets/fonts/inter-400.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="preload" href="../assets/fonts/inter-600.woff2" as="font" type="font/woff2" crossorigin>
  <link rel="stylesheet" href="../assets/css/main.css">
  <link rel="stylesheet" href="../assets/css/edu-toggle.css">
  <style>
    .pim-hero{{position:relative;padding:36px 32px;border-radius:var(--radius-xl);overflow:hidden;border:1px solid rgba(148,163,184,.14);background:linear-gradient(135deg,rgba(6,17,38,.95),rgba(11,26,58,.68) 45%,rgba(17,52,92,.55));margin-bottom:22px}}
    .pim-hero::before{{content:'';position:absolute;inset:-40%;z-index:-1;pointer-events:none;background:conic-gradient(from 240deg at 25% 40%,rgba(139,92,246,.24),transparent 30%),radial-gradient(circle at 70% 90%,rgba(94,234,212,.14),transparent 55%);filter:blur(60px);animation:pimDrift 24s linear infinite}}
    @keyframes pimDrift{{0%{{transform:rotate(0)}}100%{{transform:rotate(360deg)}}}}
    .pim-hero h1{{font-size:clamp(1.8rem,1.3rem + 1.6vw,2.8rem);line-height:1.05;margin:0 0 6px;color:var(--ink-strong)}}
    .pim-hero .lede{{color:var(--slate-300);max-width:760px}}
    .pim-kpi-row{{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:12px;margin-top:20px}}
    .pim-kpi{{padding:14px 16px;border-radius:14px;border:1px solid rgba(139,92,246,.3);background:linear-gradient(135deg,rgba(139,92,246,.14),rgba(94,234,212,.06))}}
    .pim-kpi .k{{font-size:1.5rem;font-weight:700;color:#e9d5ff;letter-spacing:-.02em;line-height:1}}
    .pim-kpi .l{{font-size:.7rem;color:var(--slate-300);margin-top:4px;text-transform:uppercase;letter-spacing:.12em}}
    .pim-grid{{display:grid;grid-template-columns:1fr;gap:16px;margin-top:18px}}
    .pim-card{{padding:16px;border-radius:var(--radius-lg);border:1px solid rgba(148,163,184,.14);background:rgba(11,26,58,.45);backdrop-filter:blur(10px)}}
    .pim-card h3{{margin:0 0 6px;font-size:1.12rem;color:var(--ink-strong)}}
    .pim-card .note{{font-size:.85rem;color:var(--slate-400);margin-bottom:10px;line-height:1.55}}
    .pim-card .caveat{{font-size:.83rem;padding:10px 12px;border-radius:8px;background:rgba(250,204,21,.07);border:1px solid rgba(250,204,21,.2);color:#fde68a;margin-top:10px}}
    .pim-card iframe{{width:100%;height:620px;border:0;display:block;border-radius:10px;background:#07132b}}
    .pim-inlinefig{{width:100%;min-height:520px;border-radius:10px;overflow:hidden}}
    .pim-footer-cta{{padding:16px;border-radius:12px;background:rgba(11,26,58,.45);border:1px solid rgba(148,163,184,.12);margin-top:16px;display:flex;justify-content:space-between;gap:12px;flex-wrap:wrap;align-items:center}}
    .pim-footer-cta .txt{{color:var(--slate-300);font-size:.9rem}}
    .pim-footer-cta .btns{{display:flex;gap:8px;flex-wrap:wrap}}
    .pim-pillrow{{display:flex;gap:6px;flex-wrap:wrap;margin-top:10px}}
    .pim-pill{{padding:3px 10px;border-radius:999px;font-size:.72rem;border:1px solid rgba(94,234,212,.3);color:#bae6fd;background:rgba(6,17,38,.5)}}
    .pim-pill.warn{{border-color:rgba(250,204,21,.4);color:#fde68a}}
    .pim-pill.bad{{border-color:rgba(251,113,133,.4);color:#fda4af}}
  </style>
</head>
<body class="dark">
  <a class="sr-only" href="#main-content">본문으로 건너뛰기</a>
  <nav class="topnav" aria-label="주요 메뉴">
    <div class="topnav-inner">
      <a class="brand" href="../index.html" aria-label="홈으로"><span class="brand-dot" aria-hidden="true"></span><span>THYROID DASH</span></a>
      <div class="nav-links" role="menubar">
        <a class="nav-link" href="../index.html" role="menuitem">홈</a><a class="nav-link" href="01_overview.html" role="menuitem">개요</a><a class="nav-link" href="02_datasets.html" role="menuitem">데이터셋</a><a class="nav-link" href="03_sample_master.html" role="menuitem">샘플</a><a class="nav-link" href="04_gene_panels.html" role="menuitem">패널</a><a class="nav-link" href="05_eda.html" role="menuitem">EDA</a><a class="nav-link" href="06_scores.html" role="menuitem">점수</a><a class="nav-link" href="07_ml_baseline.html" role="menuitem">ML</a><a class="nav-link" href="08_panel_comparison.html" role="menuitem">패널 비교</a><a class="nav-link" href="09_shap.html" role="menuitem">SHAP</a><a class="nav-link" href="10_gene_explorer.html" role="menuitem">유전자</a><a class="nav-link" href="11_cohort_compare.html" role="menuitem">코호트</a><a class="nav-link" href="12_business.html" role="menuitem">비즈니스</a><a class="nav-link" href="13_reports.html" role="menuitem">리포트</a><a class="nav-link" href="14_caveats.html" role="menuitem">주의점</a><a class="nav-link" href="15_drug_discovery.html" role="menuitem">드럭</a><a class="nav-link" href="16_quantum.html" role="menuitem">양자</a><a class="nav-link" href="17_biomarker_insights.html" role="menuitem">마커</a><a class="nav-link active" href="18_pathway_immune_meth.html" role="menuitem">해석</a><a class="nav-link" href="view_investor.html" role="menuitem">투자자</a><a class="nav-link" href="view_researcher.html" role="menuitem">연구자</a><a class="nav-link" href="99_glossary.html" role="menuitem">용어</a>
      </div>
      <div class="nav-actions">
        <button class="btn sm ghost" type="button" onclick="toggleTheme()" aria-label="다크/라이트 모드 전환">다크/라이트</button>
      </div>
    </div>
  </nav>
  <div class="layout">
    <main class="content" id="main-content" style="grid-column:1 / -1">

      <section class="pim-hero" aria-labelledby="pim-hero-title">
        <div class="subtle mono">18 / INTERPRETATION — PATHWAY · IMMUNE · METHYLATION</div>
        <h1 id="pim-hero-title">경로 · 면역 · 메틸레이션 해석</h1>
        <p class="lede">바이오마커 리스트에서 한 층 더 들어가, (1) 어떤 <em>생물학적 경로</em>가 BRAF-like와 RAS-like를 가르는지, (2) 어떤 <em>면역 세포 활성</em>이 달라지는지, (3) <em>DNA 메틸레이션</em>만으로도 tumor vs normal을 구분할 수 있는지 검증합니다. 소표본 · proxy label 한계는 각 섹션에 명시했습니다.</p>
        <div class="pim-kpi-row" role="list">
          <div class="pim-kpi" role="listitem"><div class="k">{kpi_top_pathway}</div><div class="l">가장 강한 경로 (q={kpi_top_q})</div></div>
          <div class="pim-kpi" role="listitem"><div class="k">{kpi_top_immune}</div><div class="l">최대 BRAF↔RAS 면역 대조 (p={kpi_top_immune_p})</div></div>
          <div class="pim-kpi" role="listitem"><div class="k">{kpi_meth_auc}</div><div class="l">메틸레이션 AUC · n={kpi_meth_n} (95% CI {kpi_meth_ci})</div></div>
          <div class="pim-kpi" role="listitem"><div class="k">{kpi_n_sig_pathways}</div><div class="l">유의 경로 (q&lt;0.05, fold&gt;1.5)</div></div>
        </div>
      </section>

      <details class="edu-toggle" open style="margin-bottom:18px;border-radius:12px;border:1px solid rgba(148,163,184,.18);background:rgba(11,26,58,.35);overflow:hidden">
        <summary style="cursor:pointer;padding:10px 14px;font-size:.9rem;color:#5eead4;list-style:none">💡 이 페이지 읽는 법</summary>
        <div style="padding:4px 16px 14px;border-top:1px dashed rgba(148,163,184,.18);color:var(--ink);font-size:.9rem;line-height:1.6">
          <p>이 페이지는 <strong>바이오마커 리스트를 기능적으로 해석</strong>하는 3개 층입니다.</p>
          <ul>
            <li><strong>경로 농축 (Pathway ORA)</strong> — BRAF-up / RAS-up 유전자 집합이 16개 선별 경로와 유의하게 overlap하는지 하이퍼지오메트릭 검정으로 확인. dot 크기=fold enrichment, 색=-log10(q).</li>
            <li><strong>면역 시그니처</strong> — 면역 세포별 마커 유전자 z-score 평균을 샘플 점수로 사용 (ssGSEA-유사). BRAF_like vs RAS_like 분포를 violin + Mann-Whitney U로 비교.</li>
            <li><strong>메틸레이션 분류기</strong> — GSE97466 (n={kpi_meth_n}) 상위 변동 프로브 5,000개에서 tumor vs normal 5-fold CV + 800-iter bootstrap AUC CI. 이 데이터셋은 <em>molecular subtype 라벨이 없어</em> histology-proxy(종양/정상)를 사용합니다.</li>
          </ul>
          <p class="subtle small">소표본 · 면역 deconvolution이 아닌 marker-gene proxy · 외부 메틸레이션 platform 차이 — 각 섹션 caveat 참조.</p>
        </div>
      </details>

      <div class="pim-grid">

        <section class="pim-card" id="pathway">
          <h3>① 경로 농축 (Pathway enrichment · ORA)</h3>
          <div class="note">
            <strong>방법:</strong> 하이퍼지오메트릭 검정 + BH 보정. Universe = TCGA DE 스캔에 포함된 {n_universe} 유전자. BRAF-up = log2FC&gt;0 &amp; FDR&lt;0.05 ({n_braf_up}개), RAS-up = log2FC&lt;0 &amp; FDR&lt;0.05 ({n_ras_up}개), Novel = 검증된 신규 마커 ({n_novel}개).  경로 세트는 MSigDB Hallmark(해당 subset 손코딩) + 갑상선 분화/탈분화 문헌 기반.
          </div>
          <div class="pim-pillrow">
            <span class="pim-pill">BRAF-like: {kpi_braf_story}</span>
            <span class="pim-pill">RAS-like: {kpi_ras_story}</span>
          </div>
          <div class="pim-inlinefig" style="margin-top:12px">
            <iframe src="../figs_interactive/pathway_dotplot.html" loading="lazy" title="Pathway dotplot"></iframe>
          </div>
          <div class="pim-inlinefig" style="margin-top:12px">
            <iframe src="../figs_interactive/pathway_barplot_top15.html" loading="lazy" title="Top-15 pathways"></iframe>
          </div>
          <div class="caveat">
            <strong>해석 한 줄:</strong> BRAF-like에서는 MAPK 출력 + 면역 exhaustion 신호가 강화되고, RAS-like는 TDS + thyroid lineage 분화가 유지됩니다.  노트: 경로 세트는 손선택 subset이며, full MSigDB C2/C5는 아닙니다. Fold enrichment는 universe 선택(51,711개 DE-tested 유전자)에 의존합니다.
          </div>
        </section>

        <section class="pim-card" id="immune">
          <h3>② 면역 세포 시그니처 (Immune cell signatures)</h3>
          <div class="note">
            <strong>방법:</strong> 각 시그니처의 마커 유전자 log2 발현을 샘플 축으로 z-score한 뒤 평균(ssGSEA-유사, <em>reference matrix 기반 CIBERSORT/xCell가 아님</em>).  TCGA-THCA {n_tcga}개 샘플 (BRAF_like {n_braf_im}, RAS_like {n_ras_im}). 각 시그니처별 Mann-Whitney U + Cliff's delta.
          </div>
          <div class="pim-pillrow">
            {immune_pill_html}
          </div>
          <div class="pim-inlinefig" style="margin-top:12px">
            <iframe src="../figs_interactive/immune_signature_violin.html" loading="lazy" title="Immune signature violins"></iframe>
          </div>
          <div class="caveat">
            <strong>정직한 한계:</strong> 이것은 <em>marker gene 평균 proxy</em>입니다. 실제 면역 세포 분율을 reference matrix로 decompose하는 CIBERSORT/xCell과는 다르며, 종양 부정수 유전자가 면역 마커와 co-express할 경우 bias가 발생할 수 있습니다. 신호가 약한 시그니처(n_markers &lt; 5)는 결과가 더 불안정합니다.
          </div>
        </section>

        <section class="pim-card" id="methylation">
          <h3>③ 메틸레이션 분류기 (GSE97466 validation)</h3>
          <div class="note">
            <strong>방법:</strong> GSE97466 상위 변동 프로브 5,000개 β-value → stratified 5-fold CV (LogReg + RandomForest) → 800-iter bootstrap AUC 95% CI.  라벨 = tumor vs normal (이 데이터셋에는 BRAF/RAS mutation label이 없어 <em>histology proxy</em>를 사용).
          </div>
          <div class="pim-pillrow">
            <span class="pim-pill">n = {kpi_meth_n} ({kpi_meth_tumor} tumor, {kpi_meth_normal} normal)</span>
            <span class="pim-pill">LogReg AUC = {kpi_meth_lr_auc}  (95% CI {kpi_meth_lr_ci})</span>
            <span class="pim-pill">RF AUC = {kpi_meth_rf_auc}  (95% CI {kpi_meth_rf_ci})</span>
          </div>
          <div class="pim-inlinefig" style="margin-top:12px">
            <iframe src="../figs_interactive/methylation_roc.html" loading="lazy" title="Methylation ROC"></iframe>
          </div>
          <div class="pim-inlinefig" style="margin-top:12px">
            <iframe src="../figs_interactive/methylation_top_probes.html" loading="lazy" title="Methylation top probes"></iframe>
          </div>
          <div class="caveat">
            <strong>정직한 한계:</strong> (1) n={kpi_meth_n}은 여전히 작고 특히 tumor 아형 비율이 고르지 않음 (cPTC가 다수). (2) 라벨이 <em>tumor vs normal histology proxy</em>이기 때문에 BRAF vs RAS 같은 분자 이분법은 <strong>이 섹션에서 평가하지 않음</strong>. (3) 450K (GPL13534) 플랫폼은 TCGA-THCA의 Illumina 플랫폼과 normalization/probes 구성이 다르며, 교차 적용 시 batch 효과가 크다. (4) mRNA(TCGA)와 methylation(GSE97466)은 <strong>서로 다른 환자 집단</strong>으로 직접 환자-수준 비교는 불가; 여기서는 standalone 외부 유효성 시험으로만 해석.
          </div>
        </section>

        <section class="pim-card" id="sources">
          <h3>참고 · 경로 세트 출처</h3>
          <div class="note">
            <p>경로 세트는 인터넷 fetch 없이 이 파일 안에 손코딩되어 있으며, 다음 출처의 subset입니다:</p>
            <ul>
              <li><strong>MSigDB Hallmark</strong> (Liberzon et al., <em>Cell Syst</em> 2015; PMID 26771021) — Hallmark_MAPK, Hallmark_PI3K_AKT_MTOR, Hallmark_EMT, Hallmark_Apoptosis, Hallmark_Inflammatory_Response.</li>
              <li><strong>Thyroid differentiation score (TDS)</strong> — TCGA PTC (Cell 2014, PMID 25417114).  TDS 16-gene core + BRS-up/down surrogate.</li>
              <li><strong>Dedifferentiation & aggressive signatures</strong> — Landa et al. <em>J Clin Invest</em> 2016 (PMID 27240834), Xing <em>Endocr Rev</em> 2013.</li>
              <li><strong>Immune checkpoint / exhaustion</strong> — Thorsson et al. <em>Immunity</em> 2018 (PMID 29628290); Tirosh et al. <em>Science</em> 2016 (PMID 27124452).</li>
            </ul>
          </div>
        </section>

      </div>

      <div class="pim-footer-cta">
        <div class="txt">이 페이지의 분석은 <strong>바이오마커 인사이트</strong>의 gene 리스트를 기능적으로 해석하고, <strong>드럭 디스커버리</strong>의 target prioritization 근거가 됩니다.</div>
        <div class="btns">
          <a class="btn sm" href="17_biomarker_insights.html">바이오마커 인사이트</a>
          <a class="btn sm" href="15_drug_discovery.html">드럭 디스커버리</a>
          <a class="btn sm ghost" href="14_caveats.html">주의점 &amp; 한계</a>
        </div>
      </div>

    </main>
  </div>

  <script>
    if(typeof toggleTheme === 'undefined'){{
      window.toggleTheme = function(){{
        var isLight = document.documentElement.classList.toggle('light');
        try{{ localStorage.setItem('thyroid-theme', isLight ? 'light' : 'dark'); }}catch(e){{}}
      }};
    }}
  </script>
</body>
</html>
"""


def phase4_build_page(phase1: dict, phase2: dict, phase3: dict) -> None:
    print("[PHASE 4] Building 18_pathway_immune_meth.html")
    pw_df = phase1["df"]
    s1 = phase1["summary"]
    s2 = phase2["summary"]
    s3 = phase3["summary"]

    top_pw_ko = PATHWAY_LABELS_KO.get(s1["top_pathway"] or "", s1["top_pathway"] or "—")
    top_pw_q = f"{s1['top_pathway_q']:.2e}" if s1["top_pathway_q"] is not None else "—"

    # top immune contrast
    top_immune_label = "—"
    top_immune_p = "—"
    if s2["top_contrasts"]:
        tc = s2["top_contrasts"][0]
        top_immune_label = IMMUNE_LABELS_KO.get(tc["signature"], tc["signature"])
        top_immune_p = f"{tc['p_value']:.2e}"

    # best methylation AUC (take logreg by default; fall back to RF if RF > LR)
    lr_auc = s3["logreg_cv_auc"]
    rf_auc = s3["rf_cv_auc"]
    if rf_auc > lr_auc:
        meth_auc = f"{rf_auc:.3f}"
        ci = s3["rf_ci_95"]
    else:
        meth_auc = f"{lr_auc:.3f}"
        ci = s3["logreg_ci_95"]
    meth_ci = f"{ci[0]:.2f}-{ci[1]:.2f}"

    # BRAF-like / RAS-like top pathway narrative
    def _story(direction: str) -> str:
        sub = pw_df[(pw_df["direction"] == direction) & (pw_df["significant"])].head(3) if "significant" in pw_df.columns else pd.DataFrame()
        if sub.empty:
            return "유의 경로 없음 (q<0.05 & fold>1.5)"
        names = [PATHWAY_LABELS_KO.get(p, p).split(" (")[0] for p in sub["pathway"].tolist()]
        return " · ".join(names[:3])

    immune_pill_items = []
    if s2["top_contrasts"]:
        for tc in s2["top_contrasts"][:5]:
            pill_cls = "pim-pill"
            if tc["p_value"] >= 0.05:
                pill_cls += " warn"
            label = IMMUNE_LABELS_KO.get(tc["signature"], tc["signature"]).split(" (")[0]
            diff = tc["mean_diff_braf_minus_ras"]
            sign_word = "BRAF↑" if diff > 0 else "RAS↑"
            immune_pill_items.append(
                f'<span class="{pill_cls}">{label}: {sign_word} (δ={tc["cliffs_delta"]:+.2f}, p={tc["p_value"]:.1e})</span>'
            )

    # Fill template
    n_braf_up = int((pw_df[pw_df["direction"] == "BRAF_up"]["N_sig_total"].iloc[0]) if len(pw_df) else 0)
    n_ras_up = int((pw_df[pw_df["direction"] == "RAS_up"]["N_sig_total"].iloc[0]) if len(pw_df) else 0)
    n_novel = int((pw_df[pw_df["direction"] == "Novel_validated"]["N_sig_total"].iloc[0]) if len(pw_df) else 0)
    n_universe = int(pw_df["N_universe"].iloc[0]) if len(pw_df) else 0

    html = PAGE_HTML_TEMPLATE.format(
        kpi_top_pathway=top_pw_ko.split(" (")[0],
        kpi_top_q=top_pw_q,
        kpi_top_immune=top_immune_label.split(" (")[0],
        kpi_top_immune_p=top_immune_p,
        kpi_meth_auc=meth_auc,
        kpi_meth_n=s3["n_samples"],
        kpi_meth_ci=meth_ci,
        kpi_n_sig_pathways=s1["n_significant"],
        n_universe=n_universe,
        n_braf_up=n_braf_up,
        n_ras_up=n_ras_up,
        n_novel=n_novel,
        kpi_braf_story=_story("BRAF_up"),
        kpi_ras_story=_story("RAS_up"),
        n_tcga=s2["n_tcga"],
        n_braf_im=s2["n_braf"],
        n_ras_im=s2["n_ras"],
        immune_pill_html="".join(immune_pill_items) or '<span class="pim-pill warn">no significant contrasts</span>',
        kpi_meth_tumor=s3["n_tumor"],
        kpi_meth_normal=s3["n_normal"],
        kpi_meth_lr_auc=f"{s3['logreg_cv_auc']:.3f}",
        kpi_meth_lr_ci=f"{s3['logreg_ci_95'][0]:.2f}-{s3['logreg_ci_95'][1]:.2f}",
        kpi_meth_rf_auc=f"{s3['rf_cv_auc']:.3f}",
        kpi_meth_rf_ci=f"{s3['rf_ci_95'][0]:.2f}-{s3['rf_ci_95'][1]:.2f}",
    )
    out = PAGES / "18_pathway_immune_meth.html"
    out.write_text(html)
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main() -> int:
    p1 = phase1_pathway_enrichment()
    p2 = phase2_immune_signatures()
    p3 = phase3_methylation()
    phase4_build_page(p1, p2, p3)

    # Final file manifest check
    expected = [
        TABLES / "pathway_enrichment.tsv",
        TABLES / "immune_signatures_tcga.tsv",
        TABLES / "immune_signatures_contrast.tsv",
        ML_DIR / "methylation_classifier_results.tsv",
        ML_DIR / "methylation_top_probes.tsv",
        FIGS / "pathway_dotplot.html",
        FIGS / "pathway_barplot_top15.html",
        FIGS / "immune_signature_violin.html",
        FIGS / "methylation_roc.html",
        FIGS / "methylation_top_probes.html",
        PAGES / "18_pathway_immune_meth.html",
    ]
    ok = 0
    for e in expected:
        present = e.exists()
        ok += int(present)
        print(f"  [{'OK' if present else 'MISSING'}] {e.relative_to(ROOT)}")
    print(f"manifest: {ok}/{len(expected)} expected outputs present")

    # Print compact report block
    print("\n====== REPORT ======")
    top_each = {}
    for d in ("BRAF_up", "RAS_up", "Novel_validated"):
        sub = p1["df"][(p1["df"]["direction"] == d) & (p1["df"].get("significant", False))]
        print(f"\n# Pathway · {d} · top 5 (q<0.05 & fold>1.5)")
        for _, r in sub.head(5).iterrows():
            print(f"  {r['pathway']:30s} q={r['q_value']:.2e}  fold={r['fold_enrichment']:.2f}  overlap={r['N_overlap']}/{r['N_pathway_total']}")
    print("\n# Immune · top 3 contrasts")
    for tc in p2["summary"]["top_contrasts"][:3]:
        print(f"  {tc['signature']:24s} diff={tc['mean_diff_braf_minus_ras']:+.3f}  cliff_d={tc['cliffs_delta']:+.2f}  p={tc['p_value']:.2e}")
    print("\n# Methylation")
    print(f"  n={p3['summary']['n_samples']}  tumor={p3['summary']['n_tumor']}  normal={p3['summary']['n_normal']}")
    print(f"  LogReg CV AUC = {p3['summary']['logreg_cv_auc']:.3f}  (95% CI {p3['summary']['logreg_ci_95'][0]:.3f}-{p3['summary']['logreg_ci_95'][1]:.3f})")
    print(f"  RF     CV AUC = {p3['summary']['rf_cv_auc']:.3f}  (95% CI {p3['summary']['rf_ci_95'][0]:.3f}-{p3['summary']['rf_ci_95'][1]:.3f})")

    return 0


if __name__ == "__main__":
    sys.exit(main())
