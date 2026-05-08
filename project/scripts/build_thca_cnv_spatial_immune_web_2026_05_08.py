#!/usr/bin/env python3
"""Build professor-facing web report for the THCA CNV-spatial immune paper."""

from __future__ import annotations

import html
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
BASE = ROOT / "project/results/high_impact_topic_pilots_2026_05_08"
HUB = ROOT / "project/papers_hub_2026_05_04"
OUT = HUB / "thca_cnv_spatial_immune_decision_2026_05_08.html"
ASSET_DIR = HUB / "thca_cnv_spatial_immune_assets"
GSE_DIR = ROOT / "project/data/processed/GSE250521"
FIGURE_SAMPLES = ["GSM7980864_PTC-1", "GSM7980870_LPTC-3", "GSM7980872_ATC-1"]


PATHS = {
    "cbio": BASE / "data_search/cbio_class6_dm2_dm1_armdriver_validation.tsv",
    "external": BASE / "top_topic_validation_2026_05_08/external_evidence_summary.tsv",
    "territory_coh": BASE / "spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_coherence.tsv",
    "territory_stab": BASE / "spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_stability.tsv",
    "territory_ecosystem": BASE / "spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_ecosystem_enrichment.tsv",
    "infer_conc": BASE / "infercnv_lite_fast_2026_05_08/infercnv_lite_fast_concordance_to_arm_territory.tsv",
    "infer_coh": BASE / "infercnv_lite_fast_2026_05_08/infercnv_lite_fast_coherence.tsv",
    "infer_stab": BASE / "infercnv_lite_fast_2026_05_08/infercnv_lite_fast_stability.tsv",
    "infer_ecosystem": BASE / "infercnv_lite_fast_2026_05_08/infercnv_lite_fast_ecosystem_enrichment.tsv",
    "stress": BASE / "spatial_cnv_reviewer_stress_battery_2026_05_08/reviewer_stress_scorecard.tsv",
    "block": BASE / "spatial_cnv_reviewer_stress_battery_2026_05_08/spatial_block_permutation.tsv",
    "qc_pred": BASE / "spatial_cnv_reviewer_stress_battery_2026_05_08/qc_composition_predicts_territory.tsv",
    "loso": BASE / "spatial_cnv_reviewer_stress_battery_2026_05_08/loso_influence.tsv",
    "slide_summary": BASE / "spatial_cnv_reviewer_stress_battery_2026_05_08/slide_residual_summary.tsv",
    "condition_effects": BASE / "spatial_cnv_reviewer_stress_battery_2026_05_08/condition_pooled_residual_effects.tsv",
    "cov_ladder": BASE / "spatial_cnv_covariate_ladder_2026_05_08/spatial_cnv_covariate_ladder.tsv",
    "path_st": BASE / "pathology_to_st_cnv_territory_fast_2026_05_08/pathology_to_st_cnv_metrics.tsv",
    "path_inf": BASE / "infercnv_lite_fast_2026_05_08/pathology_to_infercnv_lite_fast_metrics.tsv",
    "haiku": BASE / "multimodal_haiku_lite/multimodal_haiku_lite_metrics.tsv",
    "territory_spots": BASE / "spatial_cnv_territory_gate_2026_05_08/spatial_cnv_territory_per_spot.tsv.gz",
    "infer_spots": BASE / "infercnv_lite_fast_2026_05_08/infercnv_lite_fast_per_spot.tsv.gz",
}


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(PATHS[name], sep="\t")


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def fmt(x, digits: int = 3) -> str:
    if x is None or pd.isna(x):
        return "NA"
    if isinstance(x, str):
        return html.escape(x)
    try:
        xf = float(x)
    except Exception:
        return html.escape(str(x))
    if math.isinf(xf):
        return "inf"
    if xf != 0 and abs(xf) < 1e-4:
        return f"{xf:.2e}"
    if abs(xf) >= 1000:
        return f"{xf:,.0f}"
    return f"{xf:.{digits}g}"


def pct(x) -> str:
    if pd.isna(x):
        return "NA"
    return f"{100 * float(x):.1f}%"


def df_to_table(df: pd.DataFrame, columns: list[str] | None = None, caption: str | None = None, critical: bool = False) -> str:
    if columns is not None:
        df = df[columns].copy()
    cls = "critical" if critical else ""
    cap = f"<caption>{html.escape(caption)}</caption>" if caption else ""
    out = [f'<table class="{cls}">{cap}<thead><tr>']
    for c in df.columns:
        out.append(f"<th>{html.escape(str(c))}</th>")
    out.append("</tr></thead><tbody>")
    for _, row in df.iterrows():
        row_cls = ""
        verdict = str(row.get("verdict", ""))
        if "fails" in verdict or "inconsistent" in verdict:
            row_cls = ' class="risk-row"'
        elif "strong" in verdict:
            row_cls = ' class="ok-row"'
        elif "moderate" in verdict:
            row_cls = ' class="mid-row"'
        out.append(f"<tr{row_cls}>")
        for c in df.columns:
            v = row[c]
            if isinstance(v, (float, np.floating, int, np.integer)):
                txt = fmt(v)
            else:
                txt = html.escape(str(v))
            out.append(f"<td>{txt}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def pill(text: str) -> str:
    cls = "pill"
    if "strong" in text:
        cls += " good"
    elif "moderate" in text:
        cls += " mid"
    elif "weak" in text:
        cls += " weak"
    elif "fails" in text or "inconsistent" in text:
        cls += " bad"
    return f'<span class="{cls}">{html.escape(text)}</span>'


def flow_svg() -> str:
    labels = [
        ("1", "Bulk genomic anchor", "TCGA/cBio CNV residual class"),
        ("2", "Spatial RNA territory", "Visium arm-level + gene-bin inferCNV-lite"),
        ("3", "Reviewer stress tests", "QC, slide, condition, block permutation"),
        ("4", "Surviving biology", "TAM + cytotoxic T + TLS immune ecology"),
        ("5", "Clinical upgrade", "SNUBH validation + outcome linkage"),
    ]
    x0, y, w, h, gap = 24, 28, 190, 92, 22
    parts = ['<svg viewBox="0 0 1120 160" role="img" aria-label="paper logic flow">']
    for i, (num, title, sub) in enumerate(labels):
        x = x0 + i * (w + gap)
        parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="#fff" stroke="#7B1F2A" stroke-width="2"/>')
        parts.append(f'<circle cx="{x+24}" cy="{y+26}" r="15" fill="#7B1F2A"/><text x="{x+24}" y="{y+31}" text-anchor="middle" font-size="15" font-weight="700" fill="#fff">{num}</text>')
        parts.append(f'<text x="{x+48}" y="{y+28}" font-size="15" font-weight="700" fill="#0F1A2E">{html.escape(title)}</text>')
        lines = sub.split(" + ")
        for j, line in enumerate(lines):
            parts.append(f'<text x="{x+18}" y="{y+56+j*17}" font-size="12" fill="#3A4658">{html.escape(line)}</text>')
        if i < len(labels) - 1:
            ax = x + w + 4
            parts.append(f'<path d="M {ax} {y+46} L {ax+gap-8} {y+46}" stroke="#7B1F2A" stroke-width="2" marker-end="url(#arrow)"/>')
    parts.insert(1, '<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 Z" fill="#7B1F2A"/></marker></defs>')
    parts.append("</svg>")
    return "".join(parts)


def stress_heatmap(score: pd.DataFrame) -> str:
    outcomes = ["Macrophage_TAM_z", "Tcell_cytotoxic_z", "TLS_B_z", "TROP2_raw", "TACSTD2_z"]
    label_sets = ["arm_level_expression_cnv", "gene_bin_infercnv_lite"]
    out = ['<div class="matrix">']
    out.append('<div class="matrix-head"></div>')
    for outcome in outcomes:
        out.append(f'<div class="matrix-head">{html.escape(outcome)}</div>')
    for label in label_sets:
        out.append(f'<div class="matrix-rowhead">{html.escape(label)}</div>')
        for outcome in outcomes:
            sub = score[(score["territory_label_set"] == label) & (score["outcome"] == outcome)]
            if sub.empty:
                out.append("<div></div>")
                continue
            r = sub.iloc[0]
            out.append(
                '<div class="matrix-cell">'
                f'{pill(str(r["verdict"]))}'
                f'<small>median beta {fmt(r["median_slide_beta"])}<br>+slides {int(r["n_positive_slides"])}/{int(r["n_slides"])}<br>block beta {fmt(r["block_beta"])}, p {fmt(r["block_perm_p"])}</small>'
                "</div>"
            )
    out.append("</div>")
    return "".join(out)


def bar_svg(title: str, rows: list[tuple[str, float, str]], xmin: float | None = None, xmax: float | None = None) -> str:
    vals = [v for _, v, _ in rows if not pd.isna(v)]
    if not vals:
        return ""
    xmin = min(vals + [0]) if xmin is None else xmin
    xmax = max(vals + [0]) if xmax is None else xmax
    if xmin == xmax:
        xmax = xmin + 1
    width, row_h = 900, 34
    height = 58 + row_h * len(rows)
    left, right = 240, 40
    plot_w = width - left - right
    zero_x = left + (0 - xmin) / (xmax - xmin) * plot_w
    parts = [f'<svg viewBox="0 0 {width} {height}" class="chart" role="img" aria-label="{html.escape(title)}">']
    parts.append(f'<text x="{left}" y="24" font-size="16" font-weight="700" fill="#0F1A2E">{html.escape(title)}</text>')
    parts.append(f'<line x1="{zero_x}" x2="{zero_x}" y1="40" y2="{height-14}" stroke="#333" stroke-width="1"/>')
    for i, (label, val, color) in enumerate(rows):
        y = 52 + i * row_h
        x = left + (min(val, 0) - xmin) / (xmax - xmin) * plot_w
        bw = abs(val) / (xmax - xmin) * plot_w
        parts.append(f'<text x="18" y="{y+17}" font-size="12" fill="#3A4658">{html.escape(label)}</text>')
        parts.append(f'<rect x="{x}" y="{y}" width="{max(bw, 1)}" height="18" rx="3" fill="{color}"/>')
        tx = x + bw + 6 if val >= 0 else x - 42
        parts.append(f'<text x="{tx}" y="{y+14}" font-size="12" fill="#0F1A2E">{fmt(val)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def short_metric(name: str) -> str:
    if name == "cbio_ARMDRIVER_signature_any":
        return "arm-driver signature"
    return (
        name.replace("cbio_ARMDRIVER_", "")
        .replace("ARM_SCNA_CLUSTER=", "cluster ")
        .replace("cnv_signature_any", "any CNV signature")
    )


def short_sample(sample_id: str) -> str:
    return sample_id.split("_")[-1]


def paper_figure(fig_id: str, title: str, body: str, caption: str, critical: bool = False) -> str:
    cls = "paper-figure critical-figure" if critical else "paper-figure"
    return (
        f'<figure id="{html.escape(fig_id.lower().replace(" ", "-"))}" class="{cls}">'
        f'<div class="fig-tag">{html.escape(fig_id)}</div>'
        f"<h3>{html.escape(title)}</h3>"
        f'<div class="figure-body">{body}</div>'
        f"<figcaption>{caption}</figcaption>"
        "</figure>"
    )


def make_tissue_thumbnails(sample_ids: list[str]) -> dict[str, str]:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    out: dict[str, str] = {}
    try:
        from PIL import Image
    except Exception:
        return out
    for sample_id in sample_ids:
        src = GSE_DIR / sample_id / "tissue_hires_image.png"
        dst = ASSET_DIR / f"{sample_id}_tissue_thumb.jpg"
        if not src.exists():
            continue
        if not dst.exists() or src.stat().st_mtime > dst.stat().st_mtime:
            with Image.open(src) as im:
                im = im.convert("RGB")
                im.thumbnail((1100, 900))
                im.save(dst, format="JPEG", quality=82, optimize=True)
        out[sample_id] = f"{ASSET_DIR.name}/{dst.name}"
    return out


def tissue_image_grid(thumbs: dict[str, str]) -> str:
    if not thumbs:
        return '<div class="notice warn">H&E thumbnails were not generated in this build.</div>'
    parts = ['<div class="image-grid">']
    for sample_id in FIGURE_SAMPLES:
        src = thumbs.get(sample_id)
        if not src:
            continue
        parts.append(
            '<div class="img-panel">'
            f'<img src="{html.escape(src)}" alt="{html.escape(sample_id)} H&amp;E thumbnail"/>'
            f'<div><strong>{html.escape(short_sample(sample_id))}</strong><br/>H&amp;E tissue thumbnail</div>'
            "</div>"
        )
    parts.append("</div>")
    return "".join(parts)


def decision_svg(cbio_any: pd.Series, qc_gene: float, loso_gene: float) -> str:
    boxes = [
        ("Bulk", f"DM2 {pct(cbio_any['DM2_rate'])}\\nDM1 {pct(cbio_any['DM1_rate'])}\\nOR {fmt(cbio_any['OR'])}", "#7B1F2A"),
        ("Spatial", "Arm + gene-bin\\nterritories\\ncoherent", "#34547A"),
        ("Stress", "TAM strong\\nT cell/TLS moderate\\nTROP2 rejected", "#3C6B4F"),
        ("QC Risk", f"Apparent AUROC\\n{fmt(qc_gene)}\\nLOSO {fmt(loso_gene)}", "#A57215"),
        ("Decision", "Immune ecosystem\\nmain claim\\nSNUBH validation", "#7B1F2A"),
    ]
    width, height = 1120, 330
    parts = [f'<svg viewBox="0 0 {width} {height}" class="figure-svg" role="img" aria-label="decision overview">']
    parts.append('<defs><marker id="fig-arrow" markerWidth="9" markerHeight="9" refX="8" refY="4.5" orient="auto"><path d="M0,0 L9,4.5 L0,9 Z" fill="#7B1F2A"/></marker></defs>')
    parts.append('<rect x="20" y="18" width="1080" height="292" rx="10" fill="#fffaf4" stroke="#d9d2bf"/>')
    x0, y0, bw, bh, gap = 48, 76, 176, 142, 34
    for i, (title, sub, color) in enumerate(boxes):
        x = x0 + i * (bw + gap)
        parts.append(f'<rect x="{x}" y="{y0}" width="{bw}" height="{bh}" rx="9" fill="#fff" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<rect x="{x}" y="{y0}" width="{bw}" height="34" rx="9" fill="{color}"/>')
        parts.append(f'<text x="{x+bw/2}" y="{y0+23}" text-anchor="middle" font-size="15" font-weight="700" fill="#fff">{html.escape(title)}</text>')
        for j, line in enumerate(sub.split("\\n")):
            parts.append(f'<text x="{x+bw/2}" y="{y0+62+j*22}" text-anchor="middle" font-size="18" font-weight="700" fill="#0F1A2E">{html.escape(line)}</text>')
        if i < len(boxes) - 1:
            parts.append(f'<path d="M{x+bw+8} {y0+71} L{x+bw+gap-8} {y0+71}" stroke="#7B1F2A" stroke-width="3" marker-end="url(#fig-arrow)"/>')
    parts.append('<rect x="60" y="245" width="475" height="42" rx="7" fill="#eaf4ec" stroke="#3C6B4F"/>')
    parts.append('<text x="76" y="271" font-size="16" font-weight="700" fill="#1d432d">KEEP: CNV-like territory associated with myeloid + lymphoid immune ecology</text>')
    parts.append('<rect x="575" y="245" width="465" height="42" rx="7" fill="#fff0ec" stroke="#962E2E"/>')
    parts.append('<text x="591" y="271" font-size="16" font-weight="700" fill="#7B1F2A">DROP: title-level TROP2-CNV niche claim until SNUBH/IHC validation</text>')
    parts.append("</svg>")
    return "".join(parts)


def bulk_cnv_svg(cbio: pd.DataFrame) -> str:
    rows = cbio.head(10).copy()
    width, row_h = 1120, 42
    height = 82 + row_h * len(rows)
    left, plot_w = 250, 520
    parts = [f'<svg viewBox="0 0 {width} {height}" class="figure-svg" role="img" aria-label="bulk CNV residual architecture">']
    parts.append('<text x="24" y="32" font-size="20" font-weight="700" fill="#0F1A2E">DM2 residual class is enriched for arm-driver CNV architecture</text>')
    parts.append('<text x="24" y="55" font-size="12" fill="#3A4658">Bars show fraction of samples positive for each feature; right labels show odds ratio and p value.</text>')
    parts.append(f'<line x1="{left}" y1="66" x2="{left+plot_w}" y2="66" stroke="#d9d2bf"/>')
    for i, (_, r) in enumerate(rows.iterrows()):
        y = 82 + i * row_h
        dm2 = float(r["DM2_rate"])
        dm1 = float(r["DM1_rate"])
        parts.append(f'<text x="24" y="{y+18}" font-size="13" fill="#0F1A2E">{html.escape(short_metric(str(r["metric"])))}</text>')
        parts.append(f'<rect x="{left}" y="{y}" width="{plot_w*dm2}" height="14" rx="3" fill="#7B1F2A"/>')
        parts.append(f'<rect x="{left}" y="{y+18}" width="{plot_w*dm1}" height="14" rx="3" fill="#A9B9C8"/>')
        parts.append(f'<text x="{left+plot_w*dm2+8}" y="{y+12}" font-size="11" fill="#7B1F2A">DM2 {pct(dm2)}</text>')
        parts.append(f'<text x="{left+plot_w*dm1+8}" y="{y+30}" font-size="11" fill="#3A4658">DM1 {pct(dm1)}</text>')
        parts.append(f'<text x="{left+plot_w+95}" y="{y+20}" font-size="12" fill="#0F1A2E">OR {fmt(r["OR"])} · p {fmt(r["p"], 2)}</text>')
    parts.append('<rect x="820" y="18" width="260" height="42" rx="7" fill="#fff3ef" stroke="#962E2E" stroke-width="2"/>')
    top = rows.iloc[0]
    parts.append(f'<text x="838" y="44" font-size="15" font-weight="700" fill="#7B1F2A">Anchor: OR {fmt(top["OR"])} · p {fmt(top["p"], 2)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def spot_map_svg(spots: pd.DataFrame, label_col: str, title: str, palette: dict[str, str]) -> str:
    width, height = 1120, 380
    panel_w, panel_h = 326, 275
    parts = [f'<svg viewBox="0 0 {width} {height}" class="figure-svg" role="img" aria-label="{html.escape(title)}">']
    parts.append(f'<text x="24" y="30" font-size="20" font-weight="700" fill="#0F1A2E">{html.escape(title)}</text>')
    for i, sample_id in enumerate(FIGURE_SAMPLES):
        sub = spots[spots["sample_id"].eq(sample_id)].copy()
        if sub.empty:
            continue
        if len(sub) > 3600:
            sub = sub.sample(3600, random_state=17)
        x_min, x_max = sub["pxl_col_in_fullres"].min(), sub["pxl_col_in_fullres"].max()
        y_min, y_max = sub["pxl_row_in_fullres"].min(), sub["pxl_row_in_fullres"].max()
        ox = 38 + i * 360
        oy = 64
        parts.append(f'<rect x="{ox}" y="{oy}" width="{panel_w}" height="{panel_h}" rx="8" fill="#fbf8f1" stroke="#d9d2bf"/>')
        for _, r in sub.sort_values(label_col).iterrows():
            x = ox + 18 + (float(r["pxl_col_in_fullres"]) - x_min) / max(x_max - x_min, 1) * (panel_w - 36)
            y = oy + 18 + (float(r["pxl_row_in_fullres"]) - y_min) / max(y_max - y_min, 1) * (panel_h - 42)
            val = str(r[label_col])
            color = palette.get(val, "#8a8a8a")
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.2" fill="{color}" opacity="0.78"/>')
        parts.append(f'<text x="{ox}" y="{oy+panel_h+24}" font-size="13" font-weight="700" fill="#0F1A2E">{html.escape(short_sample(sample_id))}</text>')
        parts.append(f'<text x="{ox+70}" y="{oy+panel_h+24}" font-size="12" fill="#3A4658">n={len(sub):,} displayed spots</text>')
    lx = 780
    ly = 32
    for j, (lab, color) in enumerate(palette.items()):
        parts.append(f'<circle cx="{lx+j*96}" cy="{ly}" r="6" fill="{color}"/>')
        parts.append(f'<text x="{lx+10+j*96}" y="{ly+4}" font-size="11" fill="#3A4658">{html.escape(lab)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def spatial_territory_figure(territory_spots: pd.DataFrame, infer_spots: pd.DataFrame, thumbs: dict[str, str]) -> str:
    arm_palette = {"T0_low": "#34547A", "T1_mid": "#C7BDA9", "T2_high": "#7B1F2A"}
    gene_palette = {"T0_low": "#34547A", "T1_mid": "#C7BDA9", "T2_high": "#7B1F2A"}
    return (
        tissue_image_grid(thumbs)
        + spot_map_svg(territory_spots, "cnv_territory", "Arm-level expression-CNV territory maps", arm_palette)
        + spot_map_svg(infer_spots, "infercnv_t00_territory", "Gene-bin inferCNV-lite territory maps", gene_palette)
    )


def coherence_svg(terr_coh: pd.DataFrame, infer_coh: pd.DataFrame, infer_conc: pd.DataFrame) -> str:
    arm = terr_coh[terr_coh["condition"].isin(["PTC", "LPTC", "ATC"])].copy()
    gene = infer_coh[
        infer_coh["condition"].isin(["PTC", "LPTC", "ATC"]) & infer_coh["label"].eq("infercnv_t00_territory")
    ].copy()
    merged = arm[["sample_id", "condition", "knn_same_territory_z"]].merge(
        gene[["sample_id", "same_knn_z"]], on="sample_id", how="inner"
    )
    merged = merged.sort_values(["condition", "sample_id"])
    width, height = 1120, 365
    left, top, plot_w, plot_h = 70, 52, 920, 230
    y_max = max(100, float(max(merged["knn_same_territory_z"].max(), merged["same_knn_z"].max())) * 1.08)
    parts = [f'<svg viewBox="0 0 {width} {height}" class="figure-svg" role="img" aria-label="territory coherence">']
    parts.append('<text x="24" y="30" font-size="20" font-weight="700" fill="#0F1A2E">Spatial territories are coherent in cancer slides, but remain expression-CNV calls</text>')
    parts.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#3A4658"/>')
    parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#3A4658"/>')
    for tick in [0, 25, 50, 75, 100]:
        y = top + plot_h - tick / y_max * plot_h
        parts.append(f'<line x1="{left-5}" x2="{left+plot_w}" y1="{y}" y2="{y}" stroke="#eee6d9"/>')
        parts.append(f'<text x="{left-12}" y="{y+4}" text-anchor="end" font-size="11" fill="#3A4658">{tick}</text>')
    xs = []
    for i, (_, r) in enumerate(merged.iterrows()):
        x = left + (i + 0.5) * plot_w / len(merged)
        xs.append(x)
        ya = top + plot_h - float(r["knn_same_territory_z"]) / y_max * plot_h
        yg = top + plot_h - float(r["same_knn_z"]) / y_max * plot_h
        parts.append(f'<line x1="{x}" x2="{x}" y1="{ya}" y2="{yg}" stroke="#b9ad98" stroke-width="1.5"/>')
        parts.append(f'<circle cx="{x}" cy="{ya}" r="5" fill="#34547A"/>')
        parts.append(f'<circle cx="{x}" cy="{yg}" r="5" fill="#7B1F2A"/>')
        parts.append(f'<text x="{x}" y="{top+plot_h+22}" text-anchor="middle" font-size="10" fill="#3A4658">{html.escape(short_sample(str(r["sample_id"])))}</text>')
    parts.append('<circle cx="780" cy="28" r="5" fill="#34547A"/><text x="790" y="32" font-size="12" fill="#3A4658">arm-level</text>')
    parts.append('<circle cx="870" cy="28" r="5" fill="#7B1F2A"/><text x="880" y="32" font-size="12" fill="#3A4658">gene-bin</text>')
    conc = infer_conc[infer_conc["scope"].eq("pooled_cancer")].iloc[0]
    parts.append('<rect x="765" y="300" width="318" height="42" rx="7" fill="#fff3ef" stroke="#7B1F2A"/>')
    parts.append(f'<text x="783" y="326" font-size="15" font-weight="700" fill="#7B1F2A">Arm/gene concordance ARI {fmt(conc["territory_ARI"])}</text>')
    parts.append("</svg>")
    return "".join(parts)


def stress_grid_svg(score: pd.DataFrame) -> str:
    outcomes = ["Macrophage_TAM_z", "Tcell_cytotoxic_z", "TLS_B_z", "TROP2_raw", "TACSTD2_z"]
    labels = ["arm_level_expression_cnv", "gene_bin_infercnv_lite"]
    verdict_colors = {
        "strong_survivor": "#3C6B4F",
        "moderate_survivor": "#A57215",
        "weak_or_contextual": "#B8893C",
        "fails_or_directionally_inconsistent": "#962E2E",
    }
    width, height = 1120, 325
    left, top, cell_w, cell_h = 250, 72, 160, 82
    parts = [f'<svg viewBox="0 0 {width} {height}" class="figure-svg" role="img" aria-label="stress scorecard">']
    parts.append('<text x="24" y="30" font-size="20" font-weight="700" fill="#0F1A2E">Reviewer stress scorecard: immune ecology survives, TROP2 does not</text>')
    for j, outcome in enumerate(outcomes):
        x = left + j * cell_w
        parts.append(f'<text x="{x+cell_w/2}" y="58" text-anchor="middle" font-size="12" font-weight="700" fill="#7B1F2A">{html.escape(outcome.replace("_z", "").replace("_raw", ""))}</text>')
    for i, label in enumerate(labels):
        y = top + i * (cell_h + 18)
        parts.append(f'<text x="24" y="{y+42}" font-size="13" font-weight="700" fill="#0F1A2E">{html.escape(label.replace("_", " "))}</text>')
        for j, outcome in enumerate(outcomes):
            x = left + j * cell_w
            r = score[(score["territory_label_set"].eq(label)) & (score["outcome"].eq(outcome))].iloc[0]
            color = verdict_colors.get(str(r["verdict"]), "#777")
            parts.append(f'<rect x="{x}" y="{y}" width="{cell_w-10}" height="{cell_h}" rx="7" fill="{color}" opacity="0.92"/>')
            parts.append(f'<text x="{x+10}" y="{y+23}" font-size="12" font-weight="700" fill="#fff">{html.escape(str(r["verdict"]).replace("_", " "))}</text>')
            parts.append(f'<text x="{x+10}" y="{y+47}" font-size="11" fill="#fff">slides + {int(r["n_positive_slides"])}/{int(r["n_slides"])}</text>')
            parts.append(f'<text x="{x+10}" y="{y+66}" font-size="11" fill="#fff">block beta {fmt(r["block_beta"])}</text>')
    parts.append('<rect x="250" y="270" width="790" height="34" rx="6" fill="#fff3ef" stroke="#962E2E"/>')
    parts.append('<text x="268" y="292" font-size="14" font-weight="700" fill="#7B1F2A">Title-level claim must track the green cells: macrophage/T-cell/TLS ecosystem, not TROP2/TACSTD2.</text>')
    parts.append("</svg>")
    return "".join(parts)


def qc_predictability_svg(qc_pred: pd.DataFrame) -> str:
    rows: list[tuple[str, float, str]] = []
    for label in ["arm_level_expression_cnv", "gene_bin_infercnv_lite"]:
        app = qc_pred[(qc_pred["territory_label_set"].eq(label)) & qc_pred["mode"].eq("within_slide_apparent")]["AUROC"].median()
        loso = qc_pred[(qc_pred["territory_label_set"].eq(label)) & qc_pred["mode"].eq("leave_one_slide_out")]["AUROC"].median()
        rows.append((f"{label.replace('_', ' ')} / within-slide apparent", float(app), "#962E2E"))
        rows.append((f"{label.replace('_', ' ')} / leave-one-slide-out", float(loso), "#34547A"))
    return bar_svg("QC/composition can predict territory within slide; cross-slide performance drops", rows, xmin=0, xmax=1.0)


def condition_heatmap_svg(cond: pd.DataFrame) -> str:
    outcomes = ["Macrophage_TAM_z", "Tcell_cytotoxic_z", "TLS_B_z", "TROP2_raw", "TACSTD2_z"]
    dat = cond[
        cond["covariate_set"].eq("qc_composition")
        & cond["condition"].isin(["PTC", "LPTC", "ATC"])
        & cond["outcome"].isin(outcomes)
    ].copy()
    labels = ["arm_level_expression_cnv", "gene_bin_infercnv_lite"]
    rows = [(label, condition) for label in labels for condition in ["PTC", "LPTC", "ATC"]]
    width, height = 1120, 420
    left, top, cell_w, cell_h = 270, 78, 150, 46
    vals = dat["beta_high"].dropna().abs()
    vmax = max(0.5, float(vals.quantile(0.95)) if not vals.empty else 0.5)
    parts = [f'<svg viewBox="0 0 {width} {height}" class="figure-svg" role="img" aria-label="condition residual heatmap">']
    parts.append('<text x="24" y="30" font-size="20" font-weight="700" fill="#0F1A2E">Condition-stratified residual effects after QC/composition adjustment</text>')
    for j, outcome in enumerate(outcomes):
        x = left + j * cell_w
        parts.append(f'<text x="{x+cell_w/2}" y="62" text-anchor="middle" font-size="12" font-weight="700" fill="#7B1F2A">{html.escape(outcome.replace("_z", "").replace("_raw", ""))}</text>')
    for i, (label, condition) in enumerate(rows):
        y = top + i * cell_h
        parts.append(f'<text x="24" y="{y+28}" font-size="12" font-weight="700" fill="#0F1A2E">{html.escape(label.replace("_", " "))} · {condition}</text>')
        for j, outcome in enumerate(outcomes):
            x = left + j * cell_w
            sub = dat[(dat["territory_label_set"].eq(label)) & dat["condition"].eq(condition) & dat["outcome"].eq(outcome)]
            val = float(sub.iloc[0]["beta_high"]) if not sub.empty else 0.0
            intensity = min(abs(val) / vmax, 1.0)
            color = "#7B1F2A" if val < 0 else "#3C6B4F"
            parts.append(f'<rect x="{x}" y="{y}" width="{cell_w-10}" height="{cell_h-8}" rx="5" fill="{color}" opacity="{0.18 + 0.72*intensity:.2f}"/>')
            parts.append(f'<text x="{x+cell_w/2-5}" y="{y+25}" text-anchor="middle" font-size="12" font-weight="700" fill="#0F1A2E">{fmt(val)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def covariate_flip_svg(cov_ladder: pd.DataFrame) -> str:
    keep = ["TROP2_raw", "TACSTD2_z", "Macrophage_TAM_z", "Tcell_cytotoxic_z", "TLS_B_z"]
    rows: list[tuple[str, float, str]] = []
    order = ["sample_only", "sample_qc", "sample_qc_composition"]
    colors = {"sample_only": "#A9B9C8", "sample_qc": "#A57215", "sample_qc_composition": "#7B1F2A"}
    for outcome in keep:
        for cov in order:
            sub = cov_ladder[
                cov_ladder["territory_label_set"].eq("arm_level_expression_cnv")
                & cov_ladder["outcome"].eq(outcome)
                & cov_ladder["covariate_set"].eq(cov)
            ]
            if not sub.empty:
                rows.append((f"{outcome} / {cov}", float(sub.iloc[0]["residual_beta_high"]), colors[cov]))
    return bar_svg("Covariate ladder exposes why TROP2 is unsafe as the title claim", rows, xmin=-0.7, xmax=2.1)


def snubh_ladder_svg() -> str:
    steps = [
        ("Now", "Public TCGA + Visium\\nstress-tested hypothesis", "#34547A"),
        ("SNUBH cohort", "FFPE thyroid cancer\\ndriver-negative enriched", "#7B1F2A"),
        ("CNV confirm", "low-pass WGS / SNP array\\nor targeted CNV panel", "#7B1F2A"),
        ("TME validate", "mIF/IHC: CD68 CD163\\nCD8 GZMB CD20 CXCL13", "#3C6B4F"),
        ("Clinical impact", "recurrence / RAI-refractory\\nDSS/OS decision model", "#A57215"),
    ]
    width, height = 1120, 260
    x0, y0, bw, bh, gap = 35, 78, 185, 112, 30
    parts = [f'<svg viewBox="0 0 {width} {height}" class="figure-svg" role="img" aria-label="SNUBH validation ladder">']
    parts.append('<defs><marker id="snubh-arrow" markerWidth="9" markerHeight="9" refX="8" refY="4.5" orient="auto"><path d="M0,0 L9,4.5 L0,9 Z" fill="#7B1F2A"/></marker></defs>')
    parts.append('<text x="24" y="32" font-size="20" font-weight="700" fill="#0F1A2E">Validation ladder needed for Nature Cancer / Nature Medicine ceiling</text>')
    for i, (title, sub, color) in enumerate(steps):
        x = x0 + i * (bw + gap)
        parts.append(f'<rect x="{x}" y="{y0}" width="{bw}" height="{bh}" rx="9" fill="#fff" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<text x="{x+16}" y="{y0+28}" font-size="16" font-weight="700" fill="{color}">{html.escape(title)}</text>')
        for j, line in enumerate(sub.split("\\n")):
            parts.append(f'<text x="{x+16}" y="{y0+58+j*20}" font-size="13" fill="#0F1A2E">{html.escape(line)}</text>')
        if i < len(steps) - 1:
            parts.append(f'<path d="M{x+bw+5} {y0+56} L{x+bw+gap-8} {y0+56}" stroke="#7B1F2A" stroke-width="3" marker-end="url(#snubh-arrow)"/>')
    parts.append('<rect x="225" y="210" width="650" height="32" rx="6" fill="#fff3ef" stroke="#962E2E"/>')
    parts.append('<text x="245" y="231" font-size="14" font-weight="700" fill="#7B1F2A">Critical access ask: outcome-linked CNV-confirmed FFPE cohort, not H&amp;E alone.</text>')
    parts.append("</svg>")
    return "".join(parts)


def metric_lookup(df: pd.DataFrame, mode: str, target: str, metric: str) -> float:
    sub = df[(df["mode"] == mode) & (df["target"] == target) & (df["metric"] == metric)]
    if sub.empty:
        return np.nan
    return float(sub.iloc[0]["value"])


def main() -> None:
    cbio = read("cbio")
    external = read("external")
    terr_coh = read("territory_coh")
    terr_stab = read("territory_stab")
    terr_eco = read("territory_ecosystem")
    infer_conc = read("infer_conc")
    infer_coh = read("infer_coh")
    infer_stab = read("infer_stab")
    infer_eco = read("infer_ecosystem")
    stress = read("stress")
    block = read("block")
    qc_pred = read("qc_pred")
    loso = read("loso")
    slide_summary = read("slide_summary")
    cond = read("condition_effects")
    cov_ladder = read("cov_ladder")
    path_st = read("path_st")
    path_inf = read("path_inf")
    haiku = read("haiku")
    territory_spots = pd.read_csv(
        PATHS["territory_spots"],
        sep="\t",
        usecols=[
            "sample_id",
            "condition",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
            "cnv_territory",
            "Macrophage_TAM_z",
            "Tcell_cytotoxic_z",
            "TLS_B_z",
        ],
    )
    infer_spots = pd.read_csv(
        PATHS["infer_spots"],
        sep="\t",
        usecols=[
            "sample_id",
            "condition",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
            "infercnv_t00_territory",
            "infercnv_aneuploidy_score",
        ],
    )

    cbio_any = cbio[cbio["metric"].eq("cbio_ARMDRIVER_signature_any")].iloc[0]
    conc_cancer = infer_conc[infer_conc["scope"].eq("pooled_cancer")].iloc[0]
    cancer_infer_coh = infer_coh[(infer_coh["condition"].isin(["PTC", "LPTC", "ATC"])) & infer_coh["label"].eq("infercnv_t00_territory")]
    cancer_infer_stab = infer_stab[infer_stab["condition"].isin(["PTC", "LPTC", "ATC"])]
    cancer_arm_coh = terr_coh[terr_coh["condition"].isin(["PTC", "LPTC", "ATC"])]
    cancer_arm_stab = terr_stab[terr_stab["condition"].isin(["PTC", "LPTC", "ATC"])]

    qc_loso_arm = qc_pred[(qc_pred["territory_label_set"].eq("arm_level_expression_cnv")) & qc_pred["mode"].eq("leave_one_slide_out")].iloc[0]
    qc_loso_gene = qc_pred[(qc_pred["territory_label_set"].eq("gene_bin_infercnv_lite")) & qc_pred["mode"].eq("leave_one_slide_out")].iloc[0]
    qc_app_arm = qc_pred[(qc_pred["territory_label_set"].eq("arm_level_expression_cnv")) & qc_pred["mode"].eq("within_slide_apparent")]["AUROC"].median()
    qc_app_gene = qc_pred[(qc_pred["territory_label_set"].eq("gene_bin_infercnv_lite")) & qc_pred["mode"].eq("within_slide_apparent")]["AUROC"].median()

    path_rows = [
        ("UNI image → ST-CNV high", metric_lookup(path_st, "UNI_image_only", "t00_territory_high", "AUROC"), "#A57215"),
        ("UNI image → inferCNV-lite high", metric_lookup(path_inf, "UNI_image_only", "infercnv_t00_high", "AUROC"), "#A57215"),
        ("UNI image → TAM", metric_lookup(path_st, "UNI_image_only", "Macrophage_TAM_z", "spearman"), "#3C6B4F"),
        ("UNI image → TACSTD2", metric_lookup(path_st, "UNI_image_only", "TACSTD2_z", "spearman"), "#7B1F2A"),
        ("UNI image → inferCNV aneuploidy", metric_lookup(path_inf, "UNI_image_only", "infercnv_aneuploidy_score", "spearman"), "#34547A"),
        ("UNI image → inferCNV bin SD", metric_lookup(path_inf, "UNI_image_only", "infercnv_bin_sd", "spearman"), "#34547A"),
    ]

    block_rows = []
    for _, r in stress.sort_values(["territory_label_set", "outcome"]).iterrows():
        if r["outcome"] in ["Macrophage_TAM_z", "Tcell_cytotoxic_z", "TLS_B_z", "TROP2_raw", "TACSTD2_z"]:
            color = "#3C6B4F" if r["block_beta"] > 0 else "#962E2E"
            block_rows.append((f'{r["territory_label_set"].replace("_", " ")} / {r["outcome"]}', float(r["block_beta"]), color))

    main_evidence = pd.DataFrame(
        [
            {
                "Step": "1. Bulk genomic anchor",
                "Data used": "cBioPortal/TCGA THCA Class6 DM1/DM2 + selected arm GISTIC/CNA calls",
                "Analysis": "DM2 vs DM1 enrichment for arm-driver CNV signature and selected arms",
                "Main result": f"DM2 {pct(cbio_any['DM2_rate'])} vs DM1 {pct(cbio_any['DM1_rate'])}; OR {fmt(cbio_any['OR'])}; p {fmt(cbio_any['p'])}",
                "Decision": "Strong enough to be the genomic spine.",
            },
            {
                "Step": "2. Spatial arm-level territory",
                "Data used": "GSE250521 Visium spots, GENCODE v48, hg38 cytoband arm mapping",
                "Analysis": "Arm-level expression-CNV proxy, per-slide 3-state territory clustering, KNN coherence, leave-arm-out stability",
                "Main result": f"Median KNN z {fmt(cancer_arm_coh['knn_same_territory_z'].median())}; median ARI {fmt(cancer_arm_stab['median_ARI_leave_arm_out'].median())}",
                "Decision": "Use as CNV-like RNA territory, not definitive clone phylogeny.",
            },
            {
                "Step": "3. Gene-order inferCNV-lite sensitivity",
                "Data used": "GSE250521 raw h5ad, normal thyroid Visium reference, 227 gene-order bins",
                "Analysis": "Normal-referenced bin deviations, smoothed gene-order bins, per-slide MiniBatchKMeans territories",
                "Main result": f"Median KNN z {fmt(cancer_infer_coh['same_knn_z'].median())}; median bin-bootstrap ARI {fmt(cancer_infer_stab['median_ARI_bin_bootstrap'].median())}; arm concordance ARI {fmt(conc_cancer['territory_ARI'])}",
                "Decision": "Supports territory sensitivity but remains expression-CNV.",
            },
            {
                "Step": "4. Reviewer stress battery",
                "Data used": "Territory labels + QC metrics + tissue-composition scores + spatial blocks",
                "Analysis": "Slide residuals, condition strata, LOSO influence, spatial block permutation, QC predictability",
                "Main result": "TAM survives strong; cytotoxic/TLS survive moderate; TROP2/TACSTD2 fail or become directionally inconsistent.",
                "Decision": "Final biology claim must be immune ecosystem, not TROP2.",
            },
            {
                "Step": "5. Pathology/multimodal boundary",
                "Data used": "UNI embeddings, local cell morphology features, stage text proxy, Visium ground truth",
                "Analysis": "Leave-one-slide-out prediction of territory, ecosystem programs, inferCNV-lite scores",
                "Main result": f"UNI→ST-CNV high AUROC {fmt(path_rows[0][1])}; UNI→inferCNV high AUROC {fmt(path_rows[1][1])}; UNI→TAM rho {fmt(path_rows[2][1])}",
                "Decision": "H&E cannot replace ST-CNV; use as boundary and phenotype reader.",
            },
        ]
    )

    score_table = stress.copy()
    score_table["verdict"] = score_table["verdict"].map(lambda x: x)
    score_table = score_table[
        [
            "territory_label_set",
            "outcome",
            "verdict",
            "median_slide_beta",
            "n_positive_slides",
            "n_slides",
            "condition_positive",
            "n_conditions",
            "block_beta",
            "block_perm_p",
            "loso_sign_changes",
        ]
    ]

    fig_map = pd.DataFrame(
        [
            ["Fig 1", "Decision overview and evidence cascade", "Manual synthesis from all stress reports", "Shows the final claim and why TROP2 was demoted", "Main text"],
            ["Fig 2", "Bulk CNV residual architecture", rel(PATHS["cbio"]), "Shows DM2 enrichment for recurrent arm-CNV signature", "Main text"],
            ["Fig 3", "Spatial expression-CNV territories", rel(PATHS["territory_coh"]) + "; " + rel(PATHS["infer_coh"]), "Shows territory spatial coherence and gene-bin sensitivity", "Main text"],
            ["Fig 4", "Reviewer stress scorecard", rel(PATHS["stress"]) + "; " + rel(PATHS["block"]), "Shows which biology survives slide/QC/block stress", "Main text key figure"],
            ["Fig 5", "Pathology/multimodal boundary", rel(PATHS["path_st"]) + "; " + rel(PATHS["path_inf"]), "Shows H&E/UNI cannot robustly recover CNV territory", "Main or supplement"],
            ["Fig 6", "SNUBH validation ladder", "Prospective/retrospective data collection plan", "Defines what data are needed to move to Nature Medicine/NEJM level", "Main final figure"],
            ["Supp Fig S1", "Full condition-stratified effects", rel(PATHS["condition_effects"]), "Shows ATC/LPTC/PTC consistency and failure points", "Supplement"],
            ["Supp Fig S2", "LOSO influence", rel(BASE / "spatial_cnv_reviewer_stress_battery_2026_05_08/loso_influence.tsv"), "Shows no single slide flips immune conclusion", "Supplement"],
            ["Supp Fig S3", "QC predictability audit", rel(PATHS["qc_pred"]), "Shows expression-CNV labels are QC-sensitive", "Supplement and limitation"],
            ["Supp Fig S4", "Covariate ladder", rel(PATHS["cov_ladder"]), "Shows sample-only → QC → QC+composition changes", "Supplement"],
            ["Supp Fig S5", "inferCNV-lite concordance", rel(PATHS["infer_conc"]), "Shows gene-bin territory concordance with arm-level territory", "Supplement"],
        ],
        columns=["Item", "Panel", "Data source", "What it tests", "Placement"],
    )

    data_needs = pd.DataFrame(
        [
            ["Minimum SNUBH cohort", "Retrospective FFPE thyroid cancer cohort, preferably 200-500+ cases", "Needed for Korean clinical validation; include PTC/FVPTC/PDTC/ATC and driver-negative cases", "Required"],
            ["Clinical labels", "Age, sex, stage, histology, ATA risk, recurrence, distant metastasis, disease-specific death, follow-up duration", "Converts molecular story into clinical utility", "Required"],
            ["Treatment/outcome labels", "RAI treatment, RAI avidity, structural incomplete response, Tg trend, recurrence-free survival, DSS/OS", "Needed for Nature Medicine/NEJM-level patient relevance", "Required"],
            ["Genomics", "BRAF/RAS/RET/NTRK/ALK/TERT/DICER1/DGCR8 + CNV profile by low-pass WGS, OncoScan/SNP array, or validated targeted-CNV panel", "Separates true driver-negative and confirms CNV residual class", "Required"],
            ["Pathology images", "Diagnostic H&E WSI with tumor ROI, stromal ROI, necrosis, lymphoid aggregate annotation", "Needed for pathology boundary and deployable classifier", "Required"],
            ["Spatial validation subset", "20-40 representative cases for Visium/CosMx/Xenium/GeoMx; enrich driver-negative/CNV-high and controls", "Moves from public ST hypothesis to institutional validation", "Required for Nature Cancer"],
            ["mIF/IHC panel", "PAX8/panCK, CD68, CD163/MRC1, CD8, GZMB, CD20, CXCL13, PD-L1; TROP2 optional", "Validates TAM/cytotoxic/TLS immune ecosystem in FFPE", "Required"],
            ["Independent pathology review", "Two pathologists, blinded annotation of tumor borders, immune aggregates, necrosis, fibrosis", "Reduces reviewer concern that signatures are tissue-quality artifacts", "Required"],
            ["Prospective utility layer", "Can this identify patients needing intensified follow-up or molecular testing?", "Needed for NEJM/Nature Medicine ceiling", "Optional now, required later"],
        ],
        columns=["Data to access", "Concrete requirement", "Why it matters", "Priority"],
    )

    journal_ladder = pd.DataFrame(
        [
            ["Current evidence", "Public TCGA/cBio + public Visium + UNI/H&E + stress tests", "Mechanistic-spatial hypothesis; no institutional clinical validation", "Strong preprint / specialty translational oncology / exploratory Nature Portfolio if packaged well"],
            ["Nature Cancer path", "Add SNUBH CNV-confirmed cohort + spatial/mIF validation + mechanistic immune ecology", "Fits cancer genomics, tumor evolution/heterogeneity, tumor-host interaction, systems biology scope", "Realistic high target if validation is clean"],
            ["Nature Medicine path", "Add clinically meaningful endpoint: recurrence, RAI-refractory course, DSS/OS, clinical decision model", "Needs human-health impact and unmet clinical need; biomarker/observational clinical study acceptable if strong", "Possible only with clinical utility"],
            ["NEJM / NEJM Evidence path", "Practice-changing clinical evidence, preferably prospective or multi-institutional clinical validation", "Must validate or challenge clinical practice; discovery omics alone is insufficient", "Not current target; build toward it"],
        ],
        columns=["Level", "What evidence is needed", "Journal logic", "Realistic target"],
    )

    public_data = pd.DataFrame(
        [
            ["cBioPortal/TCGA THCA", "CNA, mutation, clinical", "Bulk genomic anchor; Class6-DM2 vs DM1 arm-driver enrichment", rel(PATHS["cbio"])],
            ["GSE250521", "Visium ST + H&E", "Spatial expression-CNV territories; ecosystem maps; UNI tile bridge", "project/data/processed/GSE250521/*.h5ad"],
            ["GENCODE v48 + hg38 cytoband", "Gene coordinate and arm mapping", "Map Visium genes to chromosome arms and gene-order bins", "project/data/external_annotations/"],
            ["UNI foundation embeddings", "H&E tile embeddings", "Pathology-to-ST territory and ecosystem retrieval", "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/"],
            ["arXiv 2605.00925 Haiku", "Multimodal conceptual reference", "Supports boundary framing: multimodal retrieval, not overclaiming H&E→CNV", rel(BASE / "arxiv_2605_00925/haiku_2605_00925.txt")],
        ],
        columns=["Dataset/source", "Modality", "Use in this project", "Local file or access path"],
    )

    methods = pd.DataFrame(
        [
            ["Arm-level expression-CNV", "Map genes to chromosome arms; score arm-level expression deviations; build T00 signature from gain/loss arms", "Creates spot-level CNV-like RNA proxy", "Sensitive to RNA quality and cell composition; use cautious wording"],
            ["Spatial territory clustering", "Per-slide KMeans into low/mid/high territories ranked by T00 signature", "Tests whether CNV-like signals form spatial domains", "Not allele-specific clone phylogeny"],
            ["Gene-order inferCNV-lite", "Normal thyroid reference, gene-order bins, smoothed bin deviations, per-slide MiniBatchKMeans", "Sensitivity analysis closer to inferCNV", "Still expression-based; formal CalicoST/infercnvpy needed"],
            ["Reviewer stress battery", "Slide residuals, condition strata, LOSO, spatial block permutation, QC predictability", "Determines which claims survive reviewer confounds", "This is the main decision layer"],
            ["Pathology bridge", "Leave-one-slide-out UNI/morphology/stage models", "Tests whether H&E can recover CNV territory or ecosystem", "CNV territory not robustly visible from H&E"],
        ],
        columns=["Analysis", "How it was done", "What it asks", "Reviewer caveat"],
    )

    sources = """
    <ul>
      <li><a href="https://www.nature.com/natcancer/aims" target="_blank">Nature Cancer aims/scope</a>: cancer genomics, tumour evolution, tumour-host interactions, immunology, systems biology and multi-omics are in scope.</li>
      <li><a href="https://www.nature.com/nm/aims" target="_blank">Nature Medicine aims/scope</a>: translational and clinical research is judged by originality, timeliness, and impact on human health; clinical, observational, biomarker and hybrid studies are in scope.</li>
      <li><a href="https://store.nejm.org/signup/evidence/articles" target="_blank">NEJM Evidence article types</a>: emphasizes clinical trials and clinically grounded research that validates or challenges clinical findings.</li>
      <li><a href="https://arxiv.org/abs/2605.00925" target="_blank">arXiv 2605.00925</a>: used as multimodal inspiration; our data support a boundary claim rather than H&E replacement of CNV mapping.</li>
    </ul>
    """

    thumbs = make_tissue_thumbnails(FIGURE_SAMPLES)
    main_figures_html = [
        paper_figure(
            "Figure 1",
            "Decision overview: what survives reviewer stress",
            decision_svg(cbio_any, qc_app_gene, float(qc_loso_gene["AUROC"])),
            "Evidence cascade used to choose the paper direction. The red decision box is intentional: TROP2 is demoted because it fails direction/QC stress, while the immune-ecosystem claim survives.",
            critical=True,
        ),
        paper_figure(
            "Figure 2",
            "Bulk CNV residual architecture in TCGA/cBio",
            bulk_cnv_svg(cbio),
            f"DM2 residual thyroid cancers show strong enrichment for arm-driver CNV features. Anchor row: DM2 {pct(cbio_any['DM2_rate'])} vs DM1 {pct(cbio_any['DM1_rate'])}, OR {fmt(cbio_any['OR'])}, p {fmt(cbio_any['p'], 2)}.",
        ),
        paper_figure(
            "Figure 3",
            "Pathology image context plus spatial transcriptome CNV-like territories",
            spatial_territory_figure(territory_spots, infer_spots, thumbs),
            "Representative H&E thumbnails are paired with spot-level arm-expression and gene-bin inferCNV-lite territory maps from the same Visium slides. This is the visual bridge for the pathology + spatial-transcriptome story.",
            critical=True,
        ),
        paper_figure(
            "Figure 4",
            "Territory coherence and arm/gene-bin concordance",
            coherence_svg(terr_coh, infer_coh, infer_conc),
            "Cancer slides show strong spatial coherence for both arm-level and gene-bin territory definitions. This supports spatial organization, but the caption explicitly keeps the claim at CNV-like RNA territory level.",
        ),
        paper_figure(
            "Figure 5",
            "Reviewer stress result: immune ecology survives, TROP2 does not",
            stress_grid_svg(stress) + bar_svg("Spatial block permutation beta after QC/composition adjustment", block_rows, xmin=-1.8, xmax=1.0),
            "The main biological figure. TAM is strong in both territory definitions; cytotoxic T and TLS are moderate. TROP2/TACSTD2 are directionally unsafe after spatial/QC stress.",
            critical=True,
        ),
        paper_figure(
            "Figure 6",
            "Multimodal boundary and SNUBH validation ladder",
            bar_svg("Pathology/UNI prediction boundary", path_rows, xmin=0, xmax=0.8) + qc_predictability_svg(qc_pred) + snubh_ladder_svg(),
            "H&E/UNI does not robustly recover CNV territory, so the high-impact path is not H&E replacement. The upgrade path is outcome-linked SNUBH validation with CNV confirmation and immune mIF/IHC.",
        ),
    ]

    supp_figures_html = [
        paper_figure(
            "Supp Fig S1",
            "Condition-stratified residual effects",
            condition_heatmap_svg(cond),
            "QC/composition-adjusted high-vs-low territory effects across PTC, LPTC, and ATC. This shows where the immune claim is stable and where TROP2/TACSTD2 weaken.",
        ),
        paper_figure(
            "Supp Fig S2",
            "Covariate ladder explaining the TROP2 demotion",
            covariate_flip_svg(cov_ladder),
            "The sample-only model makes TROP2 look exciting, but adding QC and composition flips or shrinks the signal. This is why it should stay supplemental/optional.",
        ),
        paper_figure(
            "Supp Fig S3",
            "QC/composition territory predictability audit",
            qc_predictability_svg(qc_pred),
            "Within-slide apparent AUROC is very high, while leave-one-slide-out drops. This is a limitation figure, and it protects the manuscript from overclaiming CNV clone calls.",
        ),
    ]

    html_parts = [
        "<!DOCTYPE html>",
        '<html lang="ko"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width, initial-scale=1.0"/>',
        "<title>THCA CNV-Spatial Immune Paper Decision | 2026-05-08</title>",
        '<link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond:wght@500;600;700&family=Newsreader:wght@400;500;600&family=JetBrains+Mono:wght@400;500;700&family=Noto+Sans+KR:wght@400;500;700&display=swap" rel="stylesheet"/>',
        '<link rel="stylesheet" href="style.css"/>',
        """<style>
        .toc{position:sticky;top:0;z-index:50;background:#fff;border:1px solid var(--rule);padding:10px 14px;margin-bottom:20px;border-radius:6px;display:flex;gap:12px;flex-wrap:wrap;font-family:'JetBrains Mono',monospace;font-size:11px}
        .toc a{text-decoration:none;color:var(--accent)}
        .decision-title{font-family:'Cormorant Garamond',serif;font-size:30px;line-height:1.12;color:#7B1F2A;margin:8px 0 0}
        .claim{font-size:22px;font-weight:700;color:#0F1A2E;line-height:1.35}
        .subclaim{font-size:15px;color:#3A4658}
        .metric-strip{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px;margin:18px 0}
        .metric{background:#fff;border:1px solid var(--rule-soft);border-left:4px solid var(--accent);padding:14px;border-radius:6px}
        .metric .v{font-family:'Cormorant Garamond',serif;font-size:34px;font-weight:700;color:#7B1F2A;line-height:1}
        .metric .k{font-family:'JetBrains Mono',monospace;font-size:11px;color:#3A4658;margin-top:6px}
        .section-band{background:#fff;border:1px solid var(--rule-soft);border-radius:8px;padding:20px 24px;margin:22px 0}
        .matrix{display:grid;grid-template-columns:220px repeat(5,1fr);gap:8px;align-items:stretch;margin:14px 0}
        .matrix-head{font-family:'JetBrains Mono',monospace;font-size:11px;font-weight:700;color:#7B1F2A;background:#f3e8e3;padding:8px;border-radius:4px}
        .matrix-rowhead{font-family:'JetBrains Mono',monospace;font-size:11px;font-weight:700;color:#0F1A2E;background:#eee6d9;padding:12px;border-radius:4px}
        .matrix-cell{background:#fff;border:1px solid #d9d2bf;border-radius:6px;padding:10px;min-height:108px}
        .matrix-cell small{display:block;margin-top:7px;font-family:'JetBrains Mono',monospace;font-size:10px;line-height:1.45;color:#3A4658}
        .pill{display:inline-block;padding:3px 8px;border-radius:4px;background:#ddd;font-family:'JetBrains Mono',monospace;font-size:10px;font-weight:700}
        .pill.good{background:#3C6B4F;color:#fff}.pill.mid{background:#A57215;color:#fff}.pill.weak{background:#B8893C;color:#fff}.pill.bad{background:#962E2E;color:#fff}
        .ok-row td{background:#edf6ed}.mid-row td{background:#fff7e6}.risk-row td{background:#fff0ec}
        .chart{background:#fff;border:1px solid var(--rule-soft);border-radius:8px;margin:14px 0;padding:8px;width:100%;height:auto}
        .paper-figure{background:#fff;border:1px solid var(--rule-soft);border-radius:8px;margin:22px 0;padding:18px 20px}
        .paper-figure.critical-figure{border:3px solid #962E2E;box-shadow:0 5px 18px rgba(150,46,46,.12)}
        .paper-figure h3{font-size:23px;margin:6px 0 14px;color:#0F1A2E;border:none;padding:0}
        .fig-tag{display:inline-block;background:#7B1F2A;color:#fff;padding:4px 9px;border-radius:4px;font-family:'JetBrains Mono',monospace;font-size:11px;font-weight:700;letter-spacing:.06em;text-transform:uppercase}
        .figure-body{display:block}
        .figure-svg{display:block;width:100%;height:auto;background:#fff;border:1px solid var(--rule-soft);border-radius:8px;margin:12px 0}
        figcaption{font-size:13px;line-height:1.55;color:#3A4658;border-top:1px dotted var(--rule);padding-top:10px;margin-top:10px}
        .image-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px;margin:12px 0 16px}
        .img-panel{border:1px solid var(--rule-soft);border-radius:8px;background:#fbf8f1;padding:10px;font-size:12px;color:#3A4658}
        .img-panel img{width:100%;aspect-ratio:4/3;object-fit:cover;border-radius:5px;display:block;margin-bottom:8px}
        .compact-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:14px}
        .callout-red{border:4px solid #962E2E;background:#fff3ef;padding:18px 22px;border-radius:8px;box-shadow:0 4px 18px rgba(150,46,46,.18);margin:20px 0}
        .callout-red h2,.callout-red h3{margin-top:0;color:#962E2E;border:none}
        .danger-inline{color:#962E2E;font-weight:700}
        details.supp{background:#fff;border:1px solid var(--rule-soft);border-left:4px solid var(--blue);padding:14px 18px;margin:14px 0;border-radius:6px}
        details.supp summary{cursor:pointer;font-family:'JetBrains Mono',monospace;font-size:12px;color:#34547A;font-weight:700}
        .print-btn{position:fixed;right:16px;top:16px;z-index:99;background:#7B1F2A;color:#fff;border:0;border-radius:6px;padding:9px 13px;font-family:'JetBrains Mono',monospace;font-size:12px;cursor:pointer}
        @media(max-width:900px){.matrix{grid-template-columns:1fr}.matrix-head{display:none}.container{padding:24px}.hero{margin:-24px -24px 24px!important;padding:36px 24px 24px!important}}
        @media print{.print-btn,.toc{display:none}.section-band,details.supp{break-inside:avoid}body{background:#fff}.container{max-width:none}}
        </style>""",
        "</head><body>",
        '<button class="print-btn" onclick="window.print()">PDF 저장</button>',
        '<div class="container">',
        '<nav class="topnav"><a href="index.html">Hub</a><span class="center">THCA CNV-spatial immune decision board</span><a href="master_view.html">Master</a></nav>',
        '<header class="hero"><h1>THCA CNV-Spatial Immune Paper Decision Board</h1><div class="sub">유형원 교수님 판단용: 현재 증거, reviewer-risk, figure/table 계획, SNUBH 필요 데이터, NEJM/Nature Cancer/Nature Medicine 진입 조건을 한 페이지에 정리. 2026-05-08.</div><div class="meta"><span>Build: 2026-05-08</span><span>Data: TCGA/cBio + GSE250521 + UNI + stress battery</span><span>Status: decision-ready, validation-needed</span></div></header>',
        '<div class="toc"><a href="#decision">Decision</a><a href="#logic">Logic flow</a><a href="#figures">Figures/Tables</a><a href="#methods">Data & methods</a><a href="#stress">Stress tests</a><a href="#pathology">Pathology</a><a href="#snubh">SNUBH data</a><a href="#journal">Journal path</a><a href="#supp">Supplement</a></div>',
        '<section id="decision" class="callout-red"><h2>최종 결론</h2><div class="claim">지금 가장 방어 가능한 논문은 “TROP2-CNV niche”가 아니라 <br/>Copy-number residual thyroid cancers exhibit spatially organized myeloid and lymphoid immune ecosystems 입니다.</div><p class="subclaim">TROP2/TACSTD2는 raw high-low에서는 좋아 보였지만 QC/composition + spatial block stress에서 방향이 흔들립니다. 반대로 TAM, cytotoxic T, TLS는 arm-level과 gene-bin inferCNV-lite 양쪽에서 살아남습니다.</p></section>',
        '<div class="metric-strip">',
        f'<div class="metric"><div class="v">OR {fmt(cbio_any["OR"])}</div><div class="k">Bulk CNV residual anchor<br/>DM2 {pct(cbio_any["DM2_rate"])} vs DM1 {pct(cbio_any["DM1_rate"])}</div></div>',
        f'<div class="metric"><div class="v">{fmt(cancer_infer_coh["same_knn_z"].median())}</div><div class="k">Gene-bin inferCNV-lite<br/>cancer-slide median KNN z</div></div>',
        f'<div class="metric"><div class="v">{fmt(conc_cancer["territory_ARI"])}</div><div class="k">Gene-bin vs arm-level<br/>pooled cancer territory ARI</div></div>',
        f'<div class="metric"><div class="v">10/12</div><div class="k">TAM survives<br/>positive slides in both territory definitions</div></div>',
        f'<div class="metric"><div class="v">{fmt(qc_app_gene)}</div><div class="k">QC warning<br/>within-slide AUROC for gene-bin territory</div></div>',
        "</div>",
        '<section id="logic" class="section-band"><h2>논리 흐름</h2>',
        flow_svg(),
        '<p>흐름은 간단합니다. 먼저 bulk genomics에서 driver-negative residual group이 실제 CNV architecture를 가진다는 점을 잡습니다. 그다음 spatial transcriptome에서 CNV-like RNA territory를 만들고, reviewer가 물을 QC/composition/spatial autocorrelation 문제를 스트레스 테스트합니다. 살아남는 downstream biology는 TROP2가 아니라 TAM/cytotoxic/TLS immune interface입니다.</p>',
        df_to_table(main_evidence, critical=True, caption="MAIN TABLE 1. Evidence chain and decision at each gate"),
        "</section>",
        '<section id="figures" class="section-band"><h2>Main Figures</h2><p>이제 표가 아니라 실제 figure block을 먼저 보여줍니다. 메인 그림은 bulk anchor → pathology/spatial map → reviewer stress → multimodal boundary → SNUBH validation 순서로 배치했습니다.</p>',
        *main_figures_html,
        '<h3>Figure/table provenance map</h3>',
        df_to_table(fig_map, critical=True, caption="MAIN TABLE 2. Figure/table provenance map"),
        "</section>",
        '<section id="stress" class="section-band"><h2>Reviewer Stress Battery</h2><div class="callout-red"><h3>핵심 판정</h3><p><strong>Macrophage_TAM_z</strong>만 strong으로 살아남았습니다. <strong>Tcell_cytotoxic_z</strong>와 <strong>TLS_B_z</strong>는 moderate입니다. <span class="danger-inline">TROP2/TACSTD2는 메인 claim 금지</span>입니다.</p></div>',
        stress_heatmap(stress),
        df_to_table(score_table, critical=True, caption="MAIN TABLE 3. Strict reviewer stress scorecard"),
        bar_svg("Spatial block permutation beta after QC/composition adjustment", block_rows, xmin=-1.8, xmax=1.0),
        '<p>중요한 점은 block-level에서 TROP2/TACSTD2가 음수로 돌아선다는 것입니다. 따라서 raw high-low table에서는 TROP2가 강해 보여도, spatial block과 QC/composition을 넣은 reviewer 관점에서는 메인 생물학으로 쓸 수 없습니다.</p>',
        "</section>",
        '<section id="methods" class="section-band"><h2>Data & Analysis Methods</h2>',
        df_to_table(public_data, caption="DATA TABLE 1. Public/local data used"),
        df_to_table(methods, caption="DATA TABLE 2. What each analysis does and what it cannot claim"),
        '<div class="notice warn"><strong>용어 제한:</strong> 현재 spatial layer는 allele-specific CNV clone call이 아닙니다. 논문에서는 “CNV-like spatial RNA territory”, “expression-CNV territory”, “normal-referenced gene-bin inferCNV-lite sensitivity”라고 써야 합니다. CalicoST/infercnvpy 또는 DNA-CNV validation 전에는 “clone phylogeny”라고 쓰면 위험합니다.</div>',
        "</section>",
        '<section id="pathology" class="section-band"><h2>Pathology / Multimodal Boundary</h2><p>Haiku류 multimodal paper의 교훈은 “모든 molecular state를 H&E로 대체한다”가 아니라, 어떤 것은 morphology-visible이고 어떤 것은 molecular-hidden인지 경계를 보여주는 것입니다. 우리 데이터는 CNV territory가 hidden에 가깝고, 일부 phenotype/ecosystem이 visible하다는 쪽입니다.</p>',
        bar_svg("Pathology/UNI prediction boundary", path_rows, xmin=0, xmax=0.8),
        df_to_table(path_st.sort_values(["target", "mode", "metric"]).head(24), caption="SUPP TABLE. ST-CNV territory pathology bridge metrics"),
        '<div class="notice info">UNI image-only는 ST-CNV high territory AUROC 0.534, inferCNV-lite high AUROC 0.506입니다. 즉 H&E만으로 CNV territory를 robust하게 대체할 수 없습니다. 반면 TAM/TACSTD2 등 downstream phenotype은 더 잘 읽힙니다.</div>',
        "</section>",
        '<section id="snubh" class="section-band"><h2>분당서울대병원에서 필요한 데이터</h2><div class="callout-red"><h3>가장 중요한 요청</h3><p>이 논문을 Nature Cancer/Nature Medicine급으로 올리려면 SNUBH에서 <strong>driver-negative enriched FFPE cohort + CNV confirmation + mIF/IHC immune validation + clinical outcome</strong>이 필요합니다. 단순 H&E만 추가하면 ceiling이 크게 올라가지 않습니다.</p></div>',
        df_to_table(data_needs, critical=True, caption="MAIN TABLE 4. SNUBH data access request list"),
        '<h3>권장 SNUBH validation design</h3><ol><li>기존 FFPE thyroid cancer cohort에서 driver-negative 또는 mutation-negative/low-driver case를 enrichment합니다.</li><li>같은 case에서 BRAF/RAS/RET/NTRK/ALK/TERT/DICER1/DGCR8와 arm-level CNV를 확보합니다.</li><li>CNV residual high vs control을 정의한 뒤, mIF/IHC로 CD68/CD163/MRC1, CD8/GZMB, CD20/CXCL13, PD-L1을 봅니다.</li><li>종양 영역, stroma, lymphoid aggregate, necrosis/fibrosis를 pathologist-blinded ROI로 분리합니다.</li><li>최종 endpoint는 recurrence-free survival, structural incomplete response, RAI-refractory status, DSS/OS 중 가능한 것을 씁니다.</li></ol>',
        "</section>",
        '<section id="journal" class="section-band"><h2>NEJM / Nature Cancer / Nature Medicine로 가려면</h2>',
        df_to_table(journal_ladder, critical=True, caption="MAIN TABLE 5. Journal ladder and evidence gap"),
        '<div class="notice ok"><strong>Nature Cancer 쪽이 현재 가장 현실적입니다.</strong> Nature Cancer scope에는 cancer genetics/genomics, tumour evolution and heterogeneity, tumour-host interactions, immunology, systems biology/multi-omics가 모두 들어갑니다. 현재 결과는 이 방향과 맞습니다.</div>',
        '<div class="notice warn"><strong>Nature Medicine은 clinical utility가 필요합니다.</strong> biomarker/observational/hybrid study는 가능하지만, human health impact와 unmet clinical need가 전면에 있어야 합니다. SNUBH outcome linkage 없이는 어렵습니다.</div>',
        '<div class="notice danger"><strong>NEJM은 현재 단계의 목표가 아닙니다.</strong> discovery spatial omics로는 부족하고, 임상 의사결정을 바꾸는 prospective 또는 multi-institutional clinical evidence가 필요합니다.</div>',
        sources,
        "</section>",
        '<section id="supp" class="section-band"><h2>Supplement: full details</h2>',
        *supp_figures_html,
        '<details class="supp" open><summary>Supplement A. Bulk CNV residual evidence</summary>',
        df_to_table(cbio, caption="SUPP TABLE A1. cBio Class6 DM2 vs DM1 arm-driver validation"),
        df_to_table(external.head(28), caption="SUPP TABLE A2. External evidence snapshot"),
        "</details>",
        '<details class="supp" open><summary>Supplement B. Spatial territory coherence and sensitivity</summary>',
        df_to_table(terr_coh, caption="SUPP TABLE B1. Arm-level territory coherence by slide"),
        df_to_table(infer_conc, caption="SUPP TABLE B2. inferCNV-lite concordance to arm-level territory"),
        df_to_table(infer_stab, caption="SUPP TABLE B3. inferCNV-lite bin-bootstrap stability"),
        "</details>",
        '<details class="supp" open><summary>Supplement C. Ecosystem enrichment and covariate ladder</summary>',
        df_to_table(terr_eco[terr_eco["scope"].isin(["pooled_cancer_spots", "slide_median_effect"])], caption="SUPP TABLE C1. Arm-level high-vs-low ecosystem enrichment"),
        df_to_table(infer_eco[infer_eco["scope"].isin(["pooled_cancer_spots", "slide_median_effect"])], caption="SUPP TABLE C2. inferCNV-lite high-vs-low ecosystem enrichment"),
        df_to_table(cov_ladder, caption="SUPP TABLE C3. Covariate ladder"),
        "</details>",
        '<details class="supp" open><summary>Supplement D. Reviewer stress details</summary>',
        df_to_table(slide_summary[(slide_summary["scope"].eq("all_cancer_slides")) & slide_summary["outcome"].isin(["Macrophage_TAM_z", "Tcell_cytotoxic_z", "TLS_B_z", "TROP2_raw", "TACSTD2_z"])], caption="SUPP TABLE D1. Slide-level residual summary"),
        df_to_table(cond[(cond["covariate_set"].eq("qc_composition")) & cond["outcome"].isin(["Macrophage_TAM_z", "Tcell_cytotoxic_z", "TLS_B_z", "TROP2_raw", "TACSTD2_z"])], caption="SUPP TABLE D2. Condition-stratified residual effects"),
        df_to_table(block, caption="SUPP TABLE D3. Spatial block permutation"),
        df_to_table(qc_pred, caption="SUPP TABLE D4. QC/composition predicts territory"),
        "</details>",
        '<details class="supp" open><summary>Supplement E. Pathology and multimodal metrics</summary>',
        df_to_table(path_inf, caption="SUPP TABLE E1. Pathology to inferCNV-lite metrics"),
        df_to_table(haiku, caption="SUPP TABLE E2. Haiku-lite multimodal metrics"),
        "</details>",
        '<details class="supp" open><summary>Supplement F. File registry</summary>',
        '<ul>',
    ]
    registry = [
        ROOT / "project/reports/2026_05_08_post_qc_paper_decision.md",
        ROOT / "project/reports/2026_05_08_spatial_cnv_reviewer_stress_battery.md",
        ROOT / "project/reports/2026_05_08_infercnv_lite_fast.md",
        ROOT / "project/reports/2026_05_08_spatial_cnv_territory_gate.md",
        ROOT / "project/reports/2026_05_08_pathology_to_st_cnv_territory_fast.md",
        ROOT / "project/reports/2026_05_08_haiku_multimodal_review.md",
        ROOT / "project/notebooks_or_scripts/spatial_cnv_reviewer_stress_battery_2026_05_08.py",
        ROOT / "project/notebooks_or_scripts/infercnv_lite_fast_2026_05_08.py",
        ROOT / "project/notebooks_or_scripts/pathology_to_st_cnv_territory_fast_2026_05_08.py",
    ]
    for p in registry:
        html_parts.append(f'<li><code>{html.escape(rel(p))}</code></li>')
    html_parts.extend(
        [
            "</ul></details>",
            "</section>",
            "<footer><p>Generated from local TSV/MD outputs on 2026-05-08. This page is a decision board, not a final manuscript. Any claim marked exploratory should not appear as a title-level claim.</p></footer>",
            "</div></body></html>",
        ]
    )

    OUT.write_text("\n".join(html_parts), encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
