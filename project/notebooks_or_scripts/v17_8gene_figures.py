#!/usr/bin/env python3
"""Generate the v17 8-gene network figures for the revision page."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

FIG = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17")
FIG.mkdir(parents=True, exist_ok=True)
DARK = dict(template="plotly_dark", paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)", font=dict(color="#EAEAEA"))

GENES = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]

# ============================================================
# A · Pathway enrichment bar chart (7 libraries top hit)
# ============================================================
pw_data = [
    ("Elsevier · Congenital Hypothyroidism", 7, 8, 19.93, 1.2e-19),
    ("Elsevier · NKX2-1 in Thyroid Dysgenesis", 6, 8, 16.57, 2.7e-17),
    ("Elsevier · PAX8 Targets in Thyroid Dysgenesis", 6, 8, 16.82, 1.5e-17),
    ("WikiPathway · Thyroid Hormones Production WP4746", 6, 8, 11.30, 5.0e-12),
    ("KEGG · Thyroid hormone synthesis", 5, 8, 9.45, 3.6e-10),
    ("GO BP · Thyroid Hormone Generation GO:0006590", 4, 8, 9.62, 2.4e-10),
    ("BioPlanet · TSH regulation of gene expression", 4, 8, 6.41, 3.9e-7),
    ("Reactome · Thyroxine Biosynthesis R-HSA-209968", 2, 8, 3.86, 1.4e-4),
]
df_pw = pd.DataFrame(pw_data, columns=["pathway","overlap","panel_n","neg_log10_p","p"])
df_pw = df_pw.sort_values("neg_log10_p")

fig = go.Figure()
fig.add_trace(go.Bar(
    x=df_pw["neg_log10_p"], y=df_pw["pathway"],
    orientation="h",
    text=[f"{o}/{n} · p={p:.1e}" for o,n,p in zip(df_pw["overlap"], df_pw["panel_n"], df_pw["p"])],
    textposition="outside",
    marker=dict(color=df_pw["neg_log10_p"], colorscale="Oranges", line=dict(color="#F5A623", width=1)),
))
fig.update_layout(title="8-gene panel — pathway enrichment across 7 libraries (sorted by −log10 adj P)",
                  xaxis=dict(title="−log10(adj P)", gridcolor="rgba(255,255,255,0.05)"),
                  yaxis=dict(gridcolor="rgba(255,255,255,0.0)"),
                  margin=dict(l=380, r=40, t=70, b=50), height=440, **DARK)
fig.write_html(FIG/"v17_8gene_pathway_bar.html", include_plotlyjs="cdn", full_html=True)
print("wrote pathway bar")

# ============================================================
# B · Panel vs Unsupervised top markers — Venn-like circle diagram
# ============================================================
panel = set(GENES)
unsup_top20 = ["DUSP5","SLC5A8","DUSP6","DIO1","TPO","KLK10","DUSP4","MET","DIO2","FOXP3",
               "HLA-DRA","IYD","CDKN2A","FOXE1","MMP9","FOSL1","SLC26A4","ETV4","ETV5","THRA"]
unsup = set(unsup_top20)
overlap = panel & unsup
panel_only = panel - unsup
unsup_only = unsup - panel

fig = go.Figure()
# two semi-transparent circles
fig.add_shape(type="circle", x0=0, y0=0, x1=4, y1=4,
              line=dict(color="#F5A623", width=3),
              fillcolor="rgba(245,166,35,0.15)")
fig.add_shape(type="circle", x0=2.5, y0=0, x1=6.5, y1=4,
              line=dict(color="#7ccfcd", width=3),
              fillcolor="rgba(124,207,205,0.15)")
# titles
fig.add_annotation(x=1.0, y=4.4, text="<b>8-gene panel</b><br>(canonical biology)",
                   showarrow=False, font=dict(color="#F5A623", size=14))
fig.add_annotation(x=5.5, y=4.4, text="<b>Unsupervised top 20</b><br>(data-driven discovery)",
                   showarrow=False, font=dict(color="#7ccfcd", size=14))
# panel only (left)
fig.add_annotation(x=1.0, y=2.3, text="<b>panel only · 5</b>",
                   showarrow=False, font=dict(color="#F5A623", size=12))
fig.add_annotation(x=1.0, y=1.5, text="SLC5A5 (NIS)<br>TG · TSHR<br>PAX8 · NKX2-1",
                   showarrow=False, font=dict(color="#EAEAEA", size=11))
# overlap (centre)
fig.add_annotation(x=3.5, y=2.3, text="<b>overlap · 3</b>",
                   showarrow=False, font=dict(color="#fff", size=12))
fig.add_annotation(x=3.5, y=1.5, text="TPO · DIO1<br>FOXE1",
                   showarrow=False, font=dict(color="#F5A623", size=11, family="JetBrains Mono"))
# unsup only (right)
fig.add_annotation(x=5.5, y=2.5, text="<b>unsup only · 17</b>",
                   showarrow=False, font=dict(color="#7ccfcd", size=12))
fig.add_annotation(x=5.5, y=1.4,
                   text="<b>MAPK feedback</b><br>DUSP5/6/4 · ETV4/5 · FOSL1<br><b>Inflammation</b><br>HLA-DRA · MMP9 · FOXP3<br><b>Differentiation alt.</b><br>SLC5A8 · DIO2 · IYD<br>SLC26A4 · THRA<br><b>Senescence/RTK</b><br>CDKN2A · MET · KLK10",
                   showarrow=False, font=dict(color="#EAEAEA", size=10), align="left")
fig.update_layout(title="Panel vs Unsupervised top markers — '원인' 8 panel vs '결과' unsupervised",
                  xaxis=dict(visible=False, range=[-1, 8]),
                  yaxis=dict(visible=False, range=[-0.5, 5], scaleanchor="x"),
                  height=520, margin=dict(l=20, r=20, t=70, b=20), **DARK)
fig.write_html(FIG/"v17_8gene_venn.html", include_plotlyjs="cdn", full_html=True)
print("wrote venn")

# ============================================================
# C · 8-gene functional network (manual node-edge layout)
# ============================================================
nodes = {
    # name : (x, y, role, color)
    "PAX8":   (1.0, 2.5, "TF",     "#7ccfcd"),
    "NKX2-1": (1.0, 1.7, "TF",     "#7ccfcd"),
    "FOXE1":  (1.0, 0.9, "TF",     "#7ccfcd"),
    "TSHR":   (3.0, 2.7, "Receptor","#F5A623"),
    "SLC5A5": (4.5, 2.0, "Uptake",  "#E74C3C"),
    "TPO":    (5.7, 2.7, "Synthesis","#F5A623"),
    "TG":     (5.7, 1.3, "Synthesis","#F5A623"),
    "DIO1":   (7.0, 2.0, "Activation","#9b59b6"),
}
edges = [
    ("TSHR","PAX8", "cAMP→PKA"),
    ("TSHR","NKX2-1", ""),
    ("PAX8","SLC5A5", "transcribe"),
    ("PAX8","TPO", ""),
    ("PAX8","TG", ""),
    ("NKX2-1","SLC5A5", ""),
    ("NKX2-1","TPO", ""),
    ("NKX2-1","TG", ""),
    ("FOXE1","TPO", ""),
    ("FOXE1","TG", ""),
    ("SLC5A5","TPO", "I⁻"),
    ("TPO","TG", "iodinate"),
    ("TG","DIO1", "T4 release"),
]
node_x, node_y, node_text, node_color, node_hover = [], [], [], [], []
for n,(x,y,role,col) in nodes.items():
    node_x.append(x); node_y.append(y); node_text.append(n); node_color.append(col)
    node_hover.append(f"<b>{n}</b><br>role: {role}")
edge_x, edge_y, edge_anno = [], [], []
for a,b,lab in edges:
    edge_x += [nodes[a][0], nodes[b][0], None]
    edge_y += [nodes[a][1], nodes[b][1], None]
    if lab:
        edge_anno.append(dict(x=(nodes[a][0]+nodes[b][0])/2, y=(nodes[a][1]+nodes[b][1])/2,
                              text=lab, showarrow=False, font=dict(color="#888", size=9)))
fig = go.Figure()
fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines",
                         line=dict(width=1.4, color="rgba(255,255,255,0.25)"),
                         hoverinfo="skip", showlegend=False))
fig.add_trace(go.Scatter(x=node_x, y=node_y, mode="markers+text",
                         text=node_text, textposition="middle center",
                         textfont=dict(color="#0a0d12", size=11, family="JetBrains Mono"),
                         marker=dict(size=58, color=node_color,
                                     line=dict(color="#fff", width=2)),
                         hovertext=node_hover, hoverinfo="text", showlegend=False))
# group labels
fig.add_annotation(x=1.0, y=3.4, text="<b>TF trio</b><br>(thyroid lineage masters)",
                   showarrow=False, font=dict(color="#7ccfcd", size=11))
fig.add_annotation(x=3.0, y=3.6, text="<b>TSH signal entry</b>",
                   showarrow=False, font=dict(color="#F5A623", size=11))
fig.add_annotation(x=4.5, y=2.8, text="<b>RAI uptake</b><br>(★ clinical key)",
                   showarrow=False, font=dict(color="#E74C3C", size=11))
fig.add_annotation(x=5.7, y=3.5, text="<b>Hormone synthesis</b>",
                   showarrow=False, font=dict(color="#F5A623", size=11))
fig.add_annotation(x=7.0, y=2.8, text="<b>T4→T3 activation</b>",
                   showarrow=False, font=dict(color="#9b59b6", size=11))
for a in edge_anno:
    fig.add_annotation(**a)
fig.update_layout(title="8-gene functional network — TSH→TF→NIS→TPO/TG→DIO1 chain",
                  xaxis=dict(visible=False, range=[0,8]),
                  yaxis=dict(visible=False, range=[0,4], scaleanchor="x"),
                  height=440, margin=dict(l=20, r=20, t=60, b=20), **DARK)
fig.write_html(FIG/"v17_8gene_network.html", include_plotlyjs="cdn", full_html=True)
print("wrote network")

# ============================================================
# D · Mechanism diagram — BRAF V600E → MAPK → NIS suppression → RAI loss
# ============================================================
fig = go.Figure()
# rectangles representing states
boxes = [
    (0.3, 1.7, 1.7, 2.5, "Normal thyrocyte<br>NIS↑ TPO↑ TG↑<br>RAI uptake", "#2ECC71"),
    (3.5, 1.7, 5.0, 2.5, "BRAF V600E<br>MAPK constitutively ON<br>CREB blocked", "#F5A623"),
    (6.5, 1.7, 8.0, 2.5, "PAX8 suppressed<br>NIS↓ TPO↓ TG↓<br>RAI uptake LOST", "#E74C3C"),
]
for x0,y0,x1,y1,lab,col in boxes:
    fig.add_shape(type="rect", x0=x0, y0=y0, x1=x1, y1=y1,
                  line=dict(color=col, width=2),
                  fillcolor=col.replace("#","rgba(") + "55" if False else None,
                  opacity=0.18)
    fig.add_annotation(x=(x0+x1)/2, y=(y0+y1)/2, text=lab,
                       showarrow=False, font=dict(color=col, size=11), align="center")
# arrows
fig.add_annotation(x=3.5, y=2.1, ax=2.0, ay=2.1, xref="x", yref="y", axref="x", ayref="y",
                   showarrow=True, arrowhead=2, arrowsize=1.2, arrowwidth=2, arrowcolor="#888")
fig.add_annotation(x=6.5, y=2.1, ax=5.0, ay=2.1, xref="x", yref="y", axref="x", ayref="y",
                   showarrow=True, arrowhead=2, arrowsize=1.2, arrowwidth=2, arrowcolor="#888")
# Therapeutic intervention layer
fig.add_shape(type="rect", x0=3.5, y0=0.4, x1=5.0, y1=1.1,
              line=dict(color="#7ccfcd", width=2, dash="dash"), opacity=0.20)
fig.add_annotation(x=4.25, y=0.75,
                   text="MAPK inhibition<br>(selumetinib, dabrafenib + trametinib)",
                   showarrow=False, font=dict(color="#7ccfcd", size=10), align="center")
fig.add_annotation(x=4.25, y=1.5, text="↑ block", showarrow=False, font=dict(color="#7ccfcd", size=10))
# Citations row
cit_y = -0.1
fig.add_annotation(x=4.0, y=cit_y,
                   text="Knauf 2003 (BRAF→NIS-)  ·  Chakravarty 2011 JCI (BRS52 MAPK rewiring)<br>" +
                        "Ho 2013 NEJM (selumetinib RAI redifferentiation)  ·  " +
                        "Iravani 2019 (dabrafenib+trametinib RAI restoration)",
                   showarrow=False, font=dict(color="#888", size=10), align="center")
fig.update_layout(title="MAPK constitutive activation → NIS suppression → RAI uptake loss",
                  xaxis=dict(visible=False, range=[0,8.5]),
                  yaxis=dict(visible=False, range=[-0.4, 3.0], scaleanchor="x"),
                  height=380, margin=dict(l=20,r=20,t=60,b=40), **DARK)
fig.write_html(FIG/"v17_mapk_nis_mechanism.html", include_plotlyjs="cdn", full_html=True)
print("wrote MAPK→NIS mechanism")

# ============================================================
# E · 8 gene heatmap on TCGA-THCA (sample × gene), colored by quad_group
# ============================================================
import sys
sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
try:
    sm = pd.read_csv("/opt/thyroid-dash/project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv", sep="\t")
    # Need expression matrix — try v17 cache
    expr_paths = [
        "/data/thca/data_processed/THCA/X_combined.npz",
        "/opt/thyroid-dash/project/data_processed/v5_cross_cancer/THCA/X_combined.npz",
    ]
    expr = None; gene_path = None
    for p in expr_paths:
        if Path(p).exists():
            expr = np.load(p)["X"]
            gp = Path(p).parent / "shared_genes.txt"
            if gp.exists():
                with open(gp) as fh: gene_path = [g.strip() for g in fh]
            break
    if expr is not None and gene_path is not None:
        gene_idx = {g:i for i,g in enumerate(gene_path)}
        cols = [g for g in GENES if g in gene_idx]
        if cols:
            mat = expr[:, [gene_idx[g] for g in cols]]
            # we don't have sample-id alignment guaranteed; just use first n samples
            n_show = min(200, mat.shape[0])
            mat = mat[:n_show]
            # z-score per gene
            mu = mat.mean(axis=0); sd = mat.std(axis=0)+1e-9
            z = (mat - mu) / sd
            # sort samples by mean expression (proxy DM2 -> DM1 axis)
            order = np.argsort(z.mean(axis=1))
            z = z[order]
            fig = go.Figure(go.Heatmap(z=z.T, x=[f"s{i}" for i in range(n_show)], y=cols,
                                       colorscale="RdBu_r", zmid=0, zmin=-3, zmax=3,
                                       colorbar=dict(title="z-score")))
            fig.update_layout(title=f"TCGA-THCA × 8-gene panel z-score heatmap (n={n_show}, sorted by mean panel expression)",
                              xaxis=dict(showticklabels=False, title="samples (sorted)"),
                              yaxis=dict(autorange="reversed"),
                              height=380, margin=dict(l=80,r=20,t=60,b=40), **DARK)
            fig.write_html(FIG/"v17_8gene_heatmap.html", include_plotlyjs="cdn", full_html=True)
            print(f"wrote heatmap (n={n_show}, genes={cols})")
        else:
            print("no panel genes found in shared_genes.txt; skipping heatmap")
    else:
        print("no expression matrix on disk; skipping heatmap")
except Exception as e:
    print(f"heatmap fail: {e}")

# ============================================================
# F · TERT × Age scatter with OS event overlay
# ============================================================
sm = pd.read_csv("/opt/thyroid-dash/project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv", sep="\t")
sm["tert"] = sm.get("tert_promoter_integrated","wildtype")
sm["tert_bin"] = (sm["tert"]=="mutated").astype(int)
sub = sm.dropna(subset=["age","os_event"]).copy()
sub["os_event_str"] = sub["os_event"].apply(lambda x: "event" if x==1 else "censored")
fig = px.scatter(sub, x="age", y="tds16_score_v17",
                 color="tert", symbol="os_event_str",
                 color_discrete_map={"mutated":"#E74C3C","wildtype":"#7ccfcd"},
                 symbol_map={"event":"x","censored":"circle"},
                 size_max=10,
                 title=f"Age × TDS16 colored by TERT (n={len(sub)}, OS events={int(sub['os_event'].sum())})")
fig.update_traces(marker=dict(size=8, opacity=0.7, line=dict(width=0.5, color="white")))
fig.update_layout(template="plotly_dark", paper_bgcolor="rgba(0,0,0,0)",
                  plot_bgcolor="rgba(0,0,0,0)", font=dict(color="#EAEAEA"),
                  height=420, margin=dict(l=70,r=20,t=60,b=50),
                  xaxis=dict(title="Age at diagnosis", gridcolor="rgba(255,255,255,0.05)"),
                  yaxis=dict(title="TDS16 score (v17)", gridcolor="rgba(255,255,255,0.05)"))
fig.write_html(FIG/"v17_age_tert_tds.html", include_plotlyjs="cdn", full_html=True)
print("wrote age×TDS×TERT scatter")

# ============================================================
# G · Quad group × TERT 4×2 stacked bar
# ============================================================
ct = pd.crosstab(sm["quad_group"].fillna("missing"), sm["tert"].fillna("wildtype"))
ct = ct.reindex([c for c in ["A_braf_only","B_ras_only","D_triple_negative"] if c in ct.index])
fig = go.Figure()
fig.add_trace(go.Bar(name="TERT wildtype", x=ct.index, y=ct.get("wildtype",0),
                     marker=dict(color="#7ccfcd"), text=ct.get("wildtype",0), textposition="inside"))
fig.add_trace(go.Bar(name="TERT mutated", x=ct.index, y=ct.get("mutated",0),
                     marker=dict(color="#E74C3C"), text=ct.get("mutated",0), textposition="inside"))
fig.update_layout(barmode="stack", title="Quad group × TERT promoter status (TCGA-THCA n=513)",
                  yaxis=dict(title="patients", gridcolor="rgba(255,255,255,0.05)"),
                  xaxis=dict(title=""),
                  height=380, margin=dict(l=70,r=20,t=60,b=60), **DARK)
fig.write_html(FIG/"v17_quad_tert_stack.html", include_plotlyjs="cdn", full_html=True)
print("wrote quad×TERT stack")

# ============================================================
# H · Cross-cohort TERT prevalence (TCGA → PDTC → ATC)
# ============================================================
prev = pd.DataFrame([
    ("TCGA-THCA primary PTC", 513, 7.0),
    ("MSK-thyroid 2016 PDTC", 84, 40.5),
    ("MSK-thyroid 2016 ATC",  33, 72.7),
], columns=["cohort","n","prevalence"])
fig = go.Figure(go.Bar(x=prev["cohort"], y=prev["prevalence"],
                       text=[f"{p:.1f}% (n={n})" for p,n in zip(prev["prevalence"], prev["n"])],
                       textposition="outside",
                       marker=dict(color=prev["prevalence"], colorscale="Reds",
                                   line=dict(color="#fff", width=1))))
fig.update_layout(title="TERT promoter mutation prevalence — PTC → PDTC → ATC monotonic ↑",
                  yaxis=dict(title="TERT+ %", range=[0, 90], gridcolor="rgba(255,255,255,0.05)"),
                  height=400, margin=dict(l=70,r=40,t=70,b=80), **DARK)
fig.add_annotation(x=1.0, y=85, xref="x", yref="y",
                   text="dedifferentiation 진행 단계와 monotonic하게 증가<br>v17 R3 trajectory 가설과 일치",
                   showarrow=False, font=dict(color="#7ccfcd", size=11))
fig.write_html(FIG/"v17_tert_prevalence.html", include_plotlyjs="cdn", full_html=True)
print("wrote tert prevalence")

print("\n=== ALL FIGURES WRITTEN ===")
for p in sorted(FIG.glob("v17_8gene_*.html")) + sorted(FIG.glob("v17_*.html")):
    print(f"  {p.name}  {p.stat().st_size} bytes")
