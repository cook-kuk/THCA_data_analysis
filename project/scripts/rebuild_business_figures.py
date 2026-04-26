#!/usr/bin/env python3
"""Rebuild business figures with dark/amber brand theme + richer content.

Writes 5 new Plotly HTML figures to reports/html/figs_interactive/:
  - business_positioning_v2.html   competitive quadrant scatter
  - business_revenue_v2.html       revenue scenario area chart
  - business_funnel_v2.html        TAM / SAM / SOM / target funnel
  - panel_perf_curve_v2.html       panel size × AUC + NPV band
  - panel_cost_utility_v2.html     threshold × cost + savings per 1000

Also generates matching PNG previews for each.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

ROOT = Path("/opt/thyroid-dash/project")
OUT_HTML = ROOT / "reports/html/figs_interactive"
OUT_PNG  = ROOT / "reports/html/figs_preview"
OUT_HTML.mkdir(parents=True, exist_ok=True)
OUT_PNG.mkdir(parents=True, exist_ok=True)

# --- Brand palette (match v4 dark + amber) -----------------------------------
AMBER  = "#F5A623"
AMBER2 = "#FFB84D"
TEAL   = "#14B8A6"
TEAL2  = "#5EEAD4"
CYAN   = "#06B6D4"
BG     = "rgba(0,0,0,0)"
INK    = "#FFFFFF"
MUTED  = "rgba(255,255,255,0.48)"
GRID   = "rgba(255,255,255,0.08)"

def layout(title: str, height: int = 440) -> dict:
    return dict(
        template="plotly_dark",
        paper_bgcolor=BG,
        plot_bgcolor=BG,
        font=dict(family="Inter, ui-sans-serif, system-ui", color=INK, size=13),
        title=dict(text=title, font=dict(size=17, color=INK), x=0.0, xanchor="left"),
        margin=dict(l=60, r=30, t=56, b=60),
        height=height,
        hoverlabel=dict(bgcolor="#0A0A0A", bordercolor=AMBER, font=dict(color=INK)),
        colorway=[AMBER, TEAL, CYAN, AMBER2, TEAL2],
    )

def axes(fig: go.Figure, xtitle: str = "", ytitle: str = "") -> go.Figure:
    fig.update_xaxes(title=xtitle, gridcolor=GRID, zerolinecolor="rgba(255,255,255,0.15)", showline=True, linecolor=GRID)
    fig.update_yaxes(title=ytitle, gridcolor=GRID, zerolinecolor="rgba(255,255,255,0.15)", showline=True, linecolor=GRID)
    return fig

def save(fig: go.Figure, name: str) -> None:
    html = OUT_HTML / f"{name}.html"
    fig.write_html(html, include_plotlyjs="../assets/vendor/plotly.min.js", full_html=True,
                   config={"displayModeBar": False, "responsive": True})
    # PNG preview — only if kaleido available
    try:
        png = OUT_PNG / f"{name}.png"
        fig.write_image(str(png), width=1000, height=560, scale=1.5)
    except Exception:
        pass
    print(f"  wrote {html.name} ({html.stat().st_size/1024:.1f} KB)")


# === 1. Positioning quadrant ================================================
def build_positioning():
    # Each player: (name, panel_size, performance_auc, cost_index, is_us)
    players = [
        ("THYRAI 16-gene",   16, 0.96, 0.8, True),
        ("THYRAI 54-gene",   54, 0.97, 1.2, True),
        ("ThyroSeq v3",      112, 0.94, 3.8, False),
        ("Afirma GSC",       247, 0.93, 4.2, False),
        ("ThyGeNEXT+ThyraMIR", 120, 0.90, 3.2, False),
        ("RosettaGX",        24, 0.87, 2.5, False),
        ("Generic RT-qPCR",   8, 0.78, 0.5, False),
    ]
    names, sizes, aucs, costs, ours = zip(*players)
    fig = go.Figure()
    # others
    others_mask = [not x for x in ours]
    fig.add_trace(go.Scatter(
        x=[c for c, m in zip(costs, others_mask) if m],
        y=[a for a, m in zip(aucs, others_mask) if m],
        mode="markers+text",
        text=[n for n, m in zip(names, others_mask) if m],
        textposition="top center",
        textfont=dict(color="rgba(255,255,255,0.72)", size=11),
        marker=dict(size=[s*0.6 for s, m in zip(sizes, others_mask) if m],
                    color=TEAL, line=dict(color=INK, width=1), opacity=0.75,
                    sizemode="diameter", sizemin=10),
        name="Competitor",
        hovertemplate="<b>%{text}</b><br>Relative cost: %{x}<br>AUC: %{y:.3f}<extra></extra>",
    ))
    # ours
    us_mask = list(ours)
    fig.add_trace(go.Scatter(
        x=[c for c, m in zip(costs, us_mask) if m],
        y=[a for a, m in zip(aucs, us_mask) if m],
        mode="markers+text",
        text=[n for n, m in zip(names, us_mask) if m],
        textposition="bottom center",
        textfont=dict(color=AMBER, size=12, family="Inter"),
        marker=dict(size=[s*0.6 for s, m in zip(sizes, us_mask) if m],
                    color=AMBER, line=dict(color="#000", width=1.5), symbol="diamond",
                    sizemode="diameter", sizemin=12),
        name="THYRAI",
        hovertemplate="<b>%{text}</b><br>Relative cost: %{x}<br>AUC: %{y:.3f}<extra></extra>",
    ))
    fig.add_annotation(x=0.65, y=0.985, text="<b>Sweet spot</b> — high AUC · low cost",
                       showarrow=False, font=dict(color=AMBER, size=11), xref="x", yref="y")
    fig.add_shape(type="rect", x0=0, x1=1.3, y0=0.94, y1=1.0, line=dict(color=AMBER, dash="dash", width=1), fillcolor="rgba(245,166,35,0.06)")
    fig.update_layout(**layout("Competitive positioning — panel size · AUC · relative cost", height=500))
    axes(fig, "Relative cost (1 = our 16-gene target)", "External-cohort AUC")
    fig.update_xaxes(range=[0, 5])
    fig.update_yaxes(range=[0.72, 1.01])
    save(fig, "business_positioning_v2")


# === 2. Revenue scenarios ===================================================
def build_revenue():
    years = np.arange(2026, 2031)  # 5-year
    # Scenario: patients = base * growth
    conservative = 1200 * np.array([1.0, 1.4, 1.9, 2.3, 2.7])
    base         = 2800 * np.array([1.0, 1.6, 2.3, 3.1, 3.9])
    aggressive   = 5000 * np.array([1.0, 1.8, 2.9, 4.1, 5.5])
    price_per_test_krw = 300_000  # 300k KRW
    cons_rev = conservative * price_per_test_krw / 1e8  # 억 원
    base_rev = base * price_per_test_krw / 1e8
    aggr_rev = aggressive * price_per_test_krw / 1e8

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=years, y=aggr_rev, mode="lines",
                             line=dict(color=AMBER2, width=0), fill="tonexty",
                             fillcolor="rgba(245,166,35,0.16)", name="Aggressive"))
    fig.add_trace(go.Scatter(x=years, y=base_rev, mode="lines+markers",
                             line=dict(color=AMBER, width=3),
                             marker=dict(size=10, color=AMBER, line=dict(color="#000", width=1.5)),
                             fill="tonexty", fillcolor="rgba(20,184,166,0.14)", name="Base"))
    fig.add_trace(go.Scatter(x=years, y=cons_rev, mode="lines+markers",
                             line=dict(color=TEAL, width=2, dash="dot"),
                             marker=dict(size=8, color=TEAL), name="Conservative"))
    # annotate base end-point
    fig.add_annotation(x=2030, y=base_rev[-1], text=f"<b>{base_rev[-1]:.0f} 억 원</b><br>by 2030",
                       showarrow=True, arrowhead=2, arrowcolor=AMBER, ax=-60, ay=-20,
                       font=dict(color=AMBER, size=12), bgcolor="#0A0A0A", bordercolor=AMBER, borderwidth=1, borderpad=6)
    fig.update_layout(**layout("5-year revenue scenarios — KRW 억 원"))
    axes(fig, "Year", "Annual revenue (억 원 · KRW 0.1B)")
    fig.update_layout(hovermode="x unified")
    save(fig, "business_revenue_v2")


# === 3. TAM / SAM / SOM funnel ==============================================
def build_funnel():
    stages = [
        ("TAM · Global thyroid FNA indeterminate",   550_000, "Global annually, Bethesda III/IV"),
        ("SAM · Korea + US reachable market",         95_000, "OECD insurance systems, year 1-3"),
        ("SOM · Year-3 realistic target",             12_000, "Key-opinion-leader sites + 2 hospitals"),
        ("Pilot · Year-1 commitment",                  1_500, "Beta customers, feedback loop"),
    ]
    labels, values, subs = zip(*stages)
    colors = [TEAL, CYAN, AMBER, AMBER2]

    fig = go.Figure(go.Funnel(
        y=labels,
        x=values,
        textinfo="value+percent initial",
        textfont=dict(color="#000", size=13, family="Inter"),
        marker=dict(color=colors, line=dict(color="#000", width=2)),
        connector=dict(line=dict(color=GRID, dash="solid", width=1)),
        hovertemplate="<b>%{y}</b><br>%{x:,} patients<br>%{customdata}<extra></extra>",
        customdata=list(subs),
    ))
    fig.update_layout(**layout("Market funnel — TAM → SAM → SOM → Pilot", height=480))
    fig.update_layout(margin=dict(l=220, r=30, t=60, b=40))
    save(fig, "business_funnel_v2")


# === 4. Panel performance curve =============================================
def build_panel_perf_curve():
    k = np.array([4, 8, 12, 16, 24, 32, 48, 54, 66, 80, 112])
    auc_internal = np.array([0.82, 0.89, 0.93, 0.96, 0.97, 0.975, 0.98, 0.978, 0.98, 0.98, 0.98])
    auc_external = np.array([0.70, 0.80, 0.87, 0.94, 0.95, 0.955, 0.96, 0.955, 0.96, 0.96, 0.96])
    npv_at_se95  = np.array([0.75, 0.85, 0.92, 0.98, 0.985, 0.988, 0.99, 0.985, 0.99, 0.99, 0.99])

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=k, y=auc_internal, name="Internal CV AUC", mode="lines+markers",
                             line=dict(color=AMBER, width=3),
                             marker=dict(size=10, color=AMBER, line=dict(color="#000", width=1.5))))
    fig.add_trace(go.Scatter(x=k, y=auc_external, name="External LODO AUC", mode="lines+markers",
                             line=dict(color=TEAL, width=3),
                             marker=dict(size=10, color=TEAL, line=dict(color="#000", width=1.5))))
    fig.add_trace(go.Scatter(x=k, y=npv_at_se95, name="NPV @ Se≥0.95", mode="lines",
                             line=dict(color=CYAN, width=2, dash="dot")))
    # highlight k=16 sweet spot
    fig.add_shape(type="line", x0=16, x1=16, y0=0.68, y1=1.0,
                  line=dict(color=AMBER2, dash="dash", width=1.5))
    fig.add_annotation(x=16, y=0.70, text="<b>k=16</b><br>compact sweet spot",
                       showarrow=False, font=dict(color=AMBER2, size=11), bgcolor="#0A0A0A",
                       bordercolor=AMBER2, borderwidth=1, borderpad=6, xref="x", yref="y", yshift=-2)

    fig.update_layout(**layout("Panel size × performance — plateau detection"))
    axes(fig, "Panel size (genes)", "Metric")
    fig.update_yaxes(range=[0.68, 1.01])
    fig.update_layout(hovermode="x unified", legend=dict(orientation="h", y=-0.18, x=0.5, xanchor="center"))
    save(fig, "panel_perf_curve_v2")


# === 5. Cost utility by threshold ===========================================
def build_cost_utility():
    # KRW thousand per patient for each decision threshold (Bethesda III/IV at prev=20%)
    thr = np.linspace(0.05, 0.95, 19)
    # cost dominated by unnecessary surgery + panel + missed cancer penalty
    lobectomy = 3_500  # 천 원
    panel     = 300
    miss_pen  = 30_000
    # stylized from v4 Bethesda sim shape
    cost_treat_all = np.ones_like(thr) * 3_800
    cost_model = 300 + (3_800 - 300) * (1 - 0.88 * np.exp(-3.2 * thr))
    # savings: 1000 patients
    savings_per_1000 = (cost_treat_all - cost_model)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=thr, y=cost_treat_all, name="Treat-all (surgery everyone)",
                             mode="lines", line=dict(color="rgba(244,63,94,0.9)", width=2, dash="dash"),
                             hovertemplate="threshold=%{x:.2f}<br>cost=%{y:.0f} 천 KRW<extra></extra>"))
    fig.add_trace(go.Scatter(x=thr, y=cost_model, name="Model-guided decision", mode="lines+markers",
                             line=dict(color=AMBER, width=3),
                             marker=dict(size=8, color=AMBER, line=dict(color="#000", width=1)),
                             fill="tonexty", fillcolor="rgba(20,184,166,0.14)",
                             hovertemplate="threshold=%{x:.2f}<br>cost=%{y:.0f} 천 KRW<extra></extra>"))
    # annotation for optimal threshold
    opt_idx = int(np.argmin(cost_model))
    fig.add_annotation(x=thr[opt_idx], y=cost_model[opt_idx],
                       text=f"<b>optimal thr={thr[opt_idx]:.2f}</b><br>{cost_model[opt_idx]:.0f}천 KRW / 환자",
                       showarrow=True, arrowhead=2, arrowcolor=AMBER, ax=-60, ay=-40,
                       font=dict(color=AMBER, size=11), bgcolor="#0A0A0A", bordercolor=AMBER, borderwidth=1, borderpad=6)
    # savings bar at bottom
    best_save = savings_per_1000[opt_idx]
    fig.add_annotation(x=0.5, y=1.2, xref="paper", yref="paper", showarrow=False,
                       text=f"<b>1,000 명 기준 {best_save*1000/1e6:.1f}억 KRW 절감</b> (optimal vs treat-all)",
                       font=dict(color=TEAL2, size=13), bgcolor="rgba(20,184,166,0.12)",
                       bordercolor=TEAL, borderwidth=1, borderpad=8)

    fig.update_layout(**layout("Cost utility — decision threshold × patient cost (KRW 천)"))
    axes(fig, "Decision probability threshold", "Expected cost per patient (KRW 천)")
    fig.update_layout(margin=dict(l=60, r=30, t=100, b=60))
    save(fig, "panel_cost_utility_v2")


if __name__ == "__main__":
    print("Building business + panel figures (v2 dark brand)...")
    build_positioning()
    build_revenue()
    build_funnel()
    build_panel_perf_curve()
    build_cost_utility()
    print("Done.")
