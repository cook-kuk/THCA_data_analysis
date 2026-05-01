#!/usr/bin/env python3
"""K2 preservation-type analysis figure + interactive table.

Output:
  submission/npj/figures/Korean_K2_preservation.html  (interactive: bar + donut + table)
  submission/npj/figures/Korean_K2_preservation.png   (matplotlib PNG mirror)
"""
from __future__ import annotations
from pathlib import Path
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

C_FROZEN = "#3b82f6"
C_FFPE = "#dc2626"
C_MIXED = "#f59e0b"

ROWS = [
    dict(cohort="TCGA-THCA",          n=496, role="Training (R1-A)",                role_short="train",      frozen=496, ffpe=0,   mixed=0,   conf="★★★", note="TCGA BCR core: snap-frozen 표준 100%"),
    dict(cohort="GSE27155",            n=73,  role="External meta",                  role_short="meta",        frozen=73,  ffpe=0,   mixed=0,   conf="★★★", note="Giordano 2005 Affymetrix; UMich biobank frozen"),
    dict(cohort="GSE29265",            n=49,  role="External meta",                  role_short="meta",        frozen=49,  ffpe=0,   mixed=0,   conf="★★★", note="Tomas 2012 Affymetrix; surgical frozen"),
    dict(cohort="GSE33630",            n=105, role="External meta",                  role_short="meta",        frozen=105, ffpe=0,   mixed=0,   conf="★★★", note="Tomas 2012 동일 series; HG-U133 Plus2"),
    dict(cohort="GSE76039",            n=84,  role="External held-out (AUC=0.935)", role_short="held-out",   frozen=84,  ffpe=0,   mixed=0,   conf="★★★", note="MSKCC PDTC/ATC frozen tumor RNA-seq"),
    dict(cohort="GSE65144",            n=13,  role="External (planned)",             role_short="planned",     frozen=13,  ffpe=0,   mixed=0,   conf="★★★", note="Pita 2014 Affymetrix"),
    dict(cohort="GSE151180",           n=47,  role="RAI-refractory (planned)",       role_short="planned",     frozen=47,  ffpe=0,   mixed=0,   conf="★★",  note="Colombo 2020 RAI-refractory RNA-seq"),
    dict(cohort="GSE126698",           n=50,  role="Meta-analysis",                  role_short="meta",        frozen=50,  ffpe=0,   mixed=0,   conf="★★",  note="Yim 2019 Korean PTC/N pairs"),
    dict(cohort="PRJEB11591 (Yoo'16)", n=260, role="Korean validation (K2 ★ NEW)",  role_short="korean-val",  frozen=260, ffpe=0,   mixed=0,   conf="★★★", note="SNU-GMI biobank snap-frozen, 방금 lock"),
    dict(cohort="GSE213647 (Lee'24)",  n=632, role="Korean external",                role_short="korean-ext",  frozen=383, ffpe=0,   mixed=249, conf="★★",  note="2 kits: stranded mRNA LT 383 (frozen) + TruSeq RNA Access 249 (FFPE 가능성↑)"),
    dict(cohort="MSK-IMPACT 2017",     n=93,  role="BRAF prevalence reference",      role_short="prev-ref",    frozen=0,   ffpe=93,  mixed=0,   conf="★★★", note="Clinical NGS = FFPE 표준; transcriptomic test에는 미사용"),
    dict(cohort="EGAS00001003540",     n=113, role="K3D outreach (controlled)",      role_short="planned",     frozen=113, ffpe=0,   mixed=0,   conf="★★",  note="Yoo 2019 SNU/SNUH RNA-seq+WES, 13 ATC"),
    dict(cohort="GSE256293",           n=110, role="Methylation orthogonal",         role_short="orthogonal",  frozen=0,   ffpe=0,   mixed=110, conf="★",   note="Targeted bisulfite seq → FFPE 호환 protocol; verify needed"),
]


def build_html():
    df = pd.DataFrame(ROWS).sort_values("n", ascending=True).reset_index(drop=True)

    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{"type": "bar", "rowspan": 2}, {"type": "domain"}],
               [None,                          {"type": "table"}]],
        column_widths=[0.62, 0.38],
        row_heights=[0.42, 0.58],
        horizontal_spacing=0.07,
        vertical_spacing=0.10,
        subplot_titles=(
            "<b>(a) Cohort × preservation type</b> — sample counts",
            "<b>(b) Total ecosystem (n = 2,125)</b>",
            "<b>(c) Full breakdown table</b>  ※ click column header to sort",
        ),
    )

    fig.add_trace(go.Bar(
        y=df["cohort"], x=df["frozen"], orientation="h",
        marker_color=C_FROZEN, name="Fresh-frozen",
        text=[f"{v}" if v > 0 else "" for v in df["frozen"]],
        textposition="inside", textfont=dict(color="white", size=11),
        hovertemplate="<b>%{y}</b><br>Fresh-frozen: %{x}<extra></extra>",
    ), row=1, col=1)
    fig.add_trace(go.Bar(
        y=df["cohort"], x=df["ffpe"], orientation="h",
        marker_color=C_FFPE, name="FFPE",
        text=[f"{v}" if v > 0 else "" for v in df["ffpe"]],
        textposition="inside", textfont=dict(color="white", size=11),
        hovertemplate="<b>%{y}</b><br>FFPE: %{x}<extra></extra>",
    ), row=1, col=1)
    fig.add_trace(go.Bar(
        y=df["cohort"], x=df["mixed"], orientation="h",
        marker_color=C_MIXED, name="Mixed / uncertain",
        text=[f"{v}" if v > 0 else "" for v in df["mixed"]],
        textposition="inside", textfont=dict(color="black", size=11),
        hovertemplate="<b>%{y}</b><br>Mixed/uncertain: %{x}<extra></extra>",
    ), row=1, col=1)

    total_frozen = int(df["frozen"].sum())
    total_ffpe = int(df["ffpe"].sum())
    total_mixed = int(df["mixed"].sum())
    total = total_frozen + total_ffpe + total_mixed
    fig.add_trace(go.Pie(
        labels=["Fresh-frozen", "FFPE", "Mixed / uncertain"],
        values=[total_frozen, total_ffpe, total_mixed],
        marker=dict(colors=[C_FROZEN, C_FFPE, C_MIXED], line=dict(color="white", width=2)),
        hole=0.55, sort=False,
        textinfo="label+percent",
        textfont=dict(size=11),
        hovertemplate="<b>%{label}</b><br>n = %{value}<br>%{percent}<extra></extra>",
    ), row=1, col=2)

    table_cols = ["cohort", "n", "role", "frozen", "ffpe", "mixed", "conf", "note"]
    headers = ["Cohort", "n", "Role", "Frozen", "FFPE", "Mixed", "Conf.", "Source / note"]
    df_table = df.sort_values("n", ascending=False).reset_index(drop=True)
    cell_fill_colors = []
    for _, r in df_table.iterrows():
        if r["ffpe"] > 0:
            cell_fill_colors.append("#fee2e2")
        elif r["mixed"] > 0:
            cell_fill_colors.append("#fef3c7")
        else:
            cell_fill_colors.append("#dbeafe")
    fill_matrix = [[c] * 1 for c in cell_fill_colors]
    full_fill = [cell_fill_colors] * len(table_cols)
    fig.add_trace(go.Table(
        columnwidth=[150, 50, 200, 70, 60, 70, 60, 360],
        header=dict(
            values=[f"<b>{h}</b>" for h in headers],
            fill_color="#1e3a8a", font=dict(color="white", size=12),
            align="left", height=32,
        ),
        cells=dict(
            values=[df_table[c].tolist() for c in table_cols],
            fill_color=full_fill,
            font=dict(color="#1f2937", size=10.5),
            align="left", height=26,
        ),
    ), row=2, col=2)

    fig.update_layout(
        barmode="stack",
        title=dict(
            text=f"<b>Fig K2-preservation</b> — Cohort preservation-type analysis (frozen vs FFPE vs mixed) · total samples n = {total} across 13 cohorts",
            font=dict(size=14), x=0.5, xanchor="center",
        ),
        font=dict(family="-apple-system,Segoe UI,Noto Sans KR,sans-serif", size=11),
        height=820, width=1500,
        margin=dict(l=80, r=40, t=80, b=70),
        paper_bgcolor="white", plot_bgcolor="white",
        legend=dict(orientation="h", yanchor="bottom", y=1.06, xanchor="center", x=0.31),
    )
    fig.update_xaxes(title_text="sample count", row=1, col=1, gridcolor="#e5e7eb")
    fig.update_yaxes(automargin=True, row=1, col=1)

    out_html = FIG_DIR / "Korean_K2_preservation.html"
    fig.write_html(out_html, include_plotlyjs="cdn",
                   config={"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 2}})
    out_html2 = HTML_DIR / "v17_K2_preservation.html"
    fig.write_html(out_html2, include_plotlyjs="cdn",
                   config={"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 2}})
    print(f"  ✓ {out_html}")
    print(f"  ✓ {out_html2}")


def build_png():
    df = pd.DataFrame(ROWS).sort_values("n", ascending=True).reset_index(drop=True)
    mpl.rcParams.update({"font.family": ["DejaVu Sans"], "font.size": 10})

    fig = plt.figure(figsize=(15.5, 9.0), dpi=120)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.7, 1.0], wspace=0.22)

    ax = fig.add_subplot(gs[0])
    y = np.arange(len(df))
    fr = df["frozen"].values; ff = df["ffpe"].values; mx = df["mixed"].values
    ax.barh(y, fr, color=C_FROZEN, label="Fresh-frozen", edgecolor="white")
    ax.barh(y, ff, left=fr, color=C_FFPE, label="FFPE", edgecolor="white")
    ax.barh(y, mx, left=fr + ff, color=C_MIXED, label="Mixed / uncertain", edgecolor="white")
    for i, (a, b, c) in enumerate(zip(fr, ff, mx)):
        if a > 0:
            ax.text(a / 2, i, str(int(a)), ha="center", va="center", color="white", fontsize=9, fontweight="bold")
        if b > 0:
            ax.text(a + b / 2, i, str(int(b)), ha="center", va="center", color="white", fontsize=9, fontweight="bold")
        if c > 0:
            ax.text(a + b + c / 2, i, str(int(c)), ha="center", va="center", color="black", fontsize=9, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(df["cohort"], fontsize=10)
    ax.set_xlabel("sample count")
    ax.set_title("a   Cohort × preservation type — sample counts", loc="left", fontweight="bold")
    ax.legend(loc="lower right", frameon=True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle=":", alpha=0.4)

    ax2 = fig.add_subplot(gs[1])
    total_frozen = int(df["frozen"].sum())
    total_ffpe = int(df["ffpe"].sum())
    total_mixed = int(df["mixed"].sum())
    sizes = [total_frozen, total_ffpe, total_mixed]
    labels = [f"Fresh-frozen\n(n = {total_frozen})", f"FFPE\n(n = {total_ffpe})", f"Mixed\n(n = {total_mixed})"]
    wedges, texts, autotexts = ax2.pie(
        sizes, labels=labels, colors=[C_FROZEN, C_FFPE, C_MIXED],
        autopct="%1.1f%%", startangle=90, pctdistance=0.78,
        wedgeprops=dict(edgecolor="white", linewidth=2),
    )
    for at in autotexts:
        at.set_color("white"); at.set_fontweight("bold"); at.set_fontsize(11)
    centre_circle = plt.Circle((0, 0), 0.55, fc="white")
    ax2.add_artist(centre_circle)
    ax2.text(0, 0.05, f"Total\nn = {sum(sizes):,}", ha="center", va="center", fontsize=13, fontweight="bold")
    ax2.text(0, -0.18, "13 cohorts", ha="center", va="center", fontsize=10, color="#6b7280")
    ax2.set_title("b   Total ecosystem", loc="left", fontweight="bold")

    fig.suptitle(
        "Fig K2-preservation.  Cohort preservation-type analysis — manuscript training/CV/external is fresh-frozen-only;\nFFPE exposure limited to MSK-IMPACT (BRAF prevalence reference, no transcriptomic test on FFPE).",
        fontsize=12, fontweight="bold", y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    out_png = FIG_DIR / "Korean_K2_preservation.png"
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    fig.savefig(FIG_DIR / "Korean_K2_preservation.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out_png}")


def main():
    build_html()
    build_png()


if __name__ == "__main__":
    main()
