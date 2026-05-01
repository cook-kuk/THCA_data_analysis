#!/usr/bin/env python3
"""K2 Korean cohort histology × DM1/DM2 cross-tab figure (n=260 of 262 manifest)."""
from __future__ import annotations
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import plotly.graph_objects as go
from plotly.subplots import make_subplots

FIG_DIR = Path("/opt/thyroid-dash/project/submission/npj/figures")
HTML_DIR = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17")
FIG_DIR.mkdir(parents=True, exist_ok=True)
HTML_DIR.mkdir(parents=True, exist_ok=True)

DM1_C = "#3b82f6"
DM2_C = "#dc2626"
NORMAL_C = "#10b981"

CAT_ORDER = [
    "PTC (papillary)",
    "FVPTC (follicular variant)",
    "FTC (follicular carcinoma)",
    "FA (follicular adenoma)",
    "Matched normal",
]


def parse_cat(alias):
    if not isinstance(alias, str):
        return "unknown"
    if alias.endswith("-N"):
        return "Matched normal"
    m = re.match(r"SNU-GMI-([A-Z]+)\d+", alias)
    if not m:
        return "unknown"
    code = m.group(1)
    return {
        "PTC": "PTC (papillary)",
        "FV": "FVPTC (follicular variant)",
        "FTC": "FTC (follicular carcinoma)",
        "FA": "FA (follicular adenoma)",
    }.get(code, code)


def load():
    manifest = pd.read_csv("/opt/thyroid-dash/project/results/v17_korean/K1A_prjeb11591_runs.tsv", sep="\t")
    pred = pd.read_csv("/opt/thyroid-dash/project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t")
    manifest["category"] = manifest["sample_alias"].apply(parse_cat)
    df = pred.merge(manifest[["run_accession", "sample_alias", "category"]],
                    left_on="run", right_on="run_accession", how="left")
    return df


def build_html(df):
    cats = [c for c in CAT_ORDER if c in df["category"].unique()]
    counts = df.groupby(["category", "DM_call"]).size().unstack(fill_value=0)
    counts = counts.reindex(cats)
    n_dm1 = counts.get("DM1", pd.Series([0] * len(cats), index=cats)).values
    n_dm2 = counts.get("DM2", pd.Series([0] * len(cats), index=cats)).values
    totals = n_dm1 + n_dm2
    pct_dm1 = (n_dm1 / np.maximum(totals, 1) * 100).round(1)

    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{"type": "bar"}, {"type": "violin"}],
               [{"type": "table", "colspan": 2}, None]],
        column_widths=[0.45, 0.55],
        row_heights=[0.55, 0.45],
        horizontal_spacing=0.10,
        vertical_spacing=0.10,
        subplot_titles=(
            "<b>(a) 조직형 (histology) × DM1 / DM2 분류</b>",
            "<b>(b) 조직형별 p(DM2) 분포 (violin)</b>",
            "<b>(c) DM1 으로 분류된 14 명 — 조직형 + p(DM2) 정렬</b>",
        ),
    )

    fig.add_trace(go.Bar(
        x=cats, y=n_dm1, name="DM1 (탈분화 leaning)", marker_color=DM1_C,
        text=[f"{v}" if v > 0 else "" for v in n_dm1], textposition="inside",
        textfont=dict(color="white", size=12),
        hovertemplate="<b>%{x}</b><br>DM1 = %{y}<extra></extra>",
    ), row=1, col=1)
    fig.add_trace(go.Bar(
        x=cats, y=n_dm2, name="DM2 (분화 보존)", marker_color=DM2_C,
        text=[f"{v}" if v > 0 else "" for v in n_dm2], textposition="inside",
        textfont=dict(color="white", size=12),
        hovertemplate="<b>%{x}</b><br>DM2 = %{y}<extra></extra>",
    ), row=1, col=1)
    for i, (c, t, pct) in enumerate(zip(cats, totals, pct_dm1)):
        fig.add_annotation(x=c, y=t + 2, text=f"<b>{int(t)}</b><br>(DM1 {pct}%)",
                           showarrow=False, font=dict(size=11, color="#1f2937"),
                           row=1, col=1)

    for i, c in enumerate(cats):
        sub = df[df["category"] == c]
        fig.add_trace(go.Violin(
            y=sub["p_DM2"], x=[c] * len(sub), name=c,
            box_visible=True, meanline_visible=True, showlegend=False,
            line_color=NORMAL_C if c == "Matched normal" else (DM1_C if sub["DM_call"].eq("DM1").mean() > 0.1 else DM2_C),
            fillcolor="rgba(124,207,205,0.25)",
            hovertemplate="%{y:.3f}<extra>"+c+"</extra>",
        ), row=1, col=2)
    fig.add_hline(y=0.5, line_dash="dash", line_color="#444", row=1, col=2)
    fig.update_yaxes(title_text="p(DM2)", range=[-0.05, 1.05], row=1, col=2)

    dm1_df = df[df["DM_call"] == "DM1"].sort_values("p_DM2").reset_index(drop=True)
    fig.add_trace(go.Table(
        columnwidth=[60, 120, 80, 250],
        header=dict(
            values=["<b>#</b>", "<b>Run accession</b>", "<b>Sample alias</b>", "<b>p(DM2)</b>",
                    "<b>조직형 (histology category)</b>"][:4] + ["<b>조직형 (histology category)</b>"],
            fill_color="#1e3a8a", font=dict(color="white", size=12),
            align="left", height=32,
        ),
        cells=dict(
            values=[
                list(range(1, len(dm1_df) + 1)),
                dm1_df["run"].tolist(),
                dm1_df["sample_alias"].tolist(),
                [f"{v:.3f}" for v in dm1_df["p_DM2"]],
                dm1_df["category"].tolist(),
            ],
            fill_color=[["#dbeafe"] * len(dm1_df)],
            font=dict(color="#1f2937", size=11),
            align="left", height=24,
        ),
    ), row=2, col=1)

    fig.update_layout(
        barmode="stack",
        title=dict(
            text=f"<b>Fig Q13c</b> — 한국 코호트 (PRJEB11591 manifest n = 262, 처리 n = {len(df)}) 의 조직형 × DM1/DM2 분류",
            font=dict(size=14), x=0.5, xanchor="center",
        ),
        font=dict(family="-apple-system,Segoe UI,Noto Sans KR,sans-serif", size=11),
        height=900, width=1500,
        margin=dict(l=70, r=40, t=80, b=70),
        paper_bgcolor="white", plot_bgcolor="white",
        legend=dict(orientation="h", yanchor="bottom", y=1.04, xanchor="center", x=0.22),
    )
    fig.update_xaxes(title_text="조직형 (histology category)", row=1, col=1)
    fig.update_yaxes(title_text="환자 수 (n)", row=1, col=1)

    out_html = FIG_DIR / "Korean_K2_histology.html"
    fig.write_html(out_html, include_plotlyjs="cdn",
                   config={"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 2}})
    out_html2 = HTML_DIR / "Korean_K2_histology.html"
    fig.write_html(out_html2, include_plotlyjs="cdn",
                   config={"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 2}})
    print(f"  ✓ {out_html}")
    print(f"  ✓ {out_html2}")


def build_png(df):
    cats = [c for c in CAT_ORDER if c in df["category"].unique()]
    counts = df.groupby(["category", "DM_call"]).size().unstack(fill_value=0)
    counts = counts.reindex(cats)
    n_dm1 = counts.get("DM1", pd.Series([0] * len(cats), index=cats)).values
    n_dm2 = counts.get("DM2", pd.Series([0] * len(cats), index=cats)).values

    mpl.rcParams.update({"font.family": ["DejaVu Sans"], "font.size": 10})
    fig = plt.figure(figsize=(15, 7), dpi=120)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.2], wspace=0.25)

    ax = fig.add_subplot(gs[0])
    x = np.arange(len(cats))
    ax.bar(x, n_dm1, color=DM1_C, label="DM1 (dedifferentiation-leaning)", edgecolor="white")
    ax.bar(x, n_dm2, bottom=n_dm1, color=DM2_C, label="DM2 (differentiation-preserved)", edgecolor="white")
    for i, (a, b) in enumerate(zip(n_dm1, n_dm2)):
        if a > 0:
            ax.text(i, a / 2, str(int(a)), ha="center", va="center", color="white", fontweight="bold")
        if b > 0:
            ax.text(i, a + b / 2, str(int(b)), ha="center", va="center", color="white", fontweight="bold")
        ax.text(i, a + b + 1.5, f"n={int(a+b)}\nDM1 {a/(a+b)*100:.1f}%", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    short_labels = ["PTC", "FVPTC", "FTC", "FA", "Normal"]
    ax.set_xticklabels(short_labels[:len(cats)])
    ax.set_ylabel("sample count (n)")
    ax.set_title("a   Histology category × DM1/DM2", loc="left", fontweight="bold")
    ax.legend(loc="upper left", frameon=True, fontsize=9)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    ax2 = fig.add_subplot(gs[1])
    parts = ax2.violinplot(
        [df[df["category"] == c]["p_DM2"].values for c in cats],
        positions=range(len(cats)), showmeans=True, showmedians=True, widths=0.85,
    )
    for pc in parts["bodies"]:
        pc.set_facecolor("#7ccfcd"); pc.set_alpha(0.5); pc.set_edgecolor("black")
    ax2.axhline(0.5, color="#444", linestyle="--", linewidth=0.8)
    ax2.set_xticks(range(len(cats)))
    ax2.set_xticklabels(short_labels[:len(cats)])
    ax2.set_ylabel("p(DM2)  ←  TCGA-trained LogReg")
    ax2.set_ylim(-0.05, 1.05)
    ax2.set_title(f"b   Per-category p(DM2) distribution (Korean n={len(df)})", loc="left", fontweight="bold")
    ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)

    fig.suptitle(
        f"Fig Q13c. PRJEB11591 Korean cohort (manifest n=262, processed n={len(df)}) — histology category cross-tab.\n"
        f"DM1 calls enriched in benign FA (17.4%) and follicular tumors (FTC 16.7%, FVPTC 6.3%); depleted in PTC (2.6%) — biologically coherent.",
        fontsize=11.5, fontweight="bold", y=1.0,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out_png = FIG_DIR / "Korean_K2_histology.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(FIG_DIR / "Korean_K2_histology.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out_png}")


def main():
    df = load()
    build_html(df)
    build_png(df)
    # mirror to figs_interactive
    import shutil
    shutil.copy(FIG_DIR / "Korean_K2_histology.png", HTML_DIR / "Korean_K2_histology.png")
    shutil.copy(FIG_DIR / "Korean_K2_histology.pdf", HTML_DIR / "Korean_K2_histology.pdf")
    print("Done.")


if __name__ == "__main__":
    main()
