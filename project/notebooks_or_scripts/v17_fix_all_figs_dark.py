#!/usr/bin/env python3
"""Re-render every v17_*.html figure with solid dark background + readable text.
Goal: figures must look correct standalone (not just inside iframe)."""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

FIG = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17")

# SOLID dark — not transparent
BG = "#0b0e12"
INK = "#F2F2F2"
DARK = dict(template="plotly_dark", paper_bgcolor=BG,
            plot_bgcolor=BG, font=dict(color=INK, size=14))

GENES = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]

# ------------------------------------------------------------
# 1. Pathway enrichment bar
# ------------------------------------------------------------
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
fig = go.Figure(go.Bar(
    x=df_pw["neg_log10_p"], y=df_pw["pathway"], orientation="h",
    text=[f"{o}/{n} · p={p:.1e}" for o,n,p in zip(df_pw["overlap"], df_pw["panel_n"], df_pw["p"])],
    textposition="outside", textfont=dict(size=12, color=INK),
    marker=dict(color=df_pw["neg_log10_p"], colorscale="Oranges",
                line=dict(color="#F5A623", width=1)),
))
fig.update_layout(title=dict(text="<b>8-gene panel — pathway enrichment (7 libraries, sorted by −log10 adj P)</b>",
                             font=dict(size=15, color=INK)),
                  xaxis=dict(title=dict(text="−log10(adj P)", font=dict(size=13, color=INK)),
                             tickfont=dict(size=12, color=INK), gridcolor="rgba(255,255,255,0.08)"),
                  yaxis=dict(tickfont=dict(size=12, color=INK)),
                  margin=dict(l=380, r=60, t=70, b=50), height=460, **DARK)
fig.write_html(FIG/"v17_8gene_pathway_bar.html", include_plotlyjs="cdn", full_html=True)
print("✓ pathway_bar")

# ------------------------------------------------------------
# 2. Venn (cleaner, larger fonts)
# ------------------------------------------------------------
fig = go.Figure()
fig.add_shape(type="circle", x0=0, y0=0, x1=5, y1=5,
              line=dict(color="#F5A623", width=4),
              fillcolor="rgba(245,166,35,0.18)")
fig.add_shape(type="circle", x0=3.5, y0=0, x1=8.5, y1=5,
              line=dict(color="#7ccfcd", width=4),
              fillcolor="rgba(124,207,205,0.18)")
fig.add_annotation(x=1.5, y=5.5, text="<b>8-gene panel</b><br><i>(canonical biology)</i>",
                   showarrow=False, font=dict(color="#F5A623", size=18))
fig.add_annotation(x=7.0, y=5.5, text="<b>Unsupervised top 20</b><br><i>(data-driven)</i>",
                   showarrow=False, font=dict(color="#7ccfcd", size=18))
fig.add_annotation(x=1.6, y=2.7, text="<b>panel only<br>5 genes</b>",
                   showarrow=False, font=dict(color="#F5A623", size=15))
fig.add_annotation(x=1.6, y=1.5, text="SLC5A5 (NIS)<br>TG · TSHR<br>PAX8 · NKX2-1",
                   showarrow=False, font=dict(color=INK, size=14, family="JetBrains Mono"))
fig.add_annotation(x=4.25, y=2.7, text="<b>overlap<br>3 genes</b>",
                   showarrow=False, font=dict(color="#fff", size=15))
fig.add_annotation(x=4.25, y=1.6, text="TPO<br>DIO1<br>FOXE1",
                   showarrow=False, font=dict(color="#F5A623", size=14, family="JetBrains Mono"))
fig.add_annotation(x=6.9, y=2.7, text="<b>unsup only<br>17 genes</b>",
                   showarrow=False, font=dict(color="#7ccfcd", size=15))
fig.add_annotation(x=6.9, y=1.0,
                   text="<b>MAPK feedback</b>: DUSP5/6/4 · ETV4/5 · FOSL1<br>" +
                        "<b>Inflammation</b>: HLA-DRA · MMP9 · FOXP3<br>" +
                        "<b>Differentiation alt.</b>: SLC5A8 · DIO2 · IYD · SLC26A4 · THRA<br>" +
                        "<b>Senescence/RTK</b>: CDKN2A · MET · KLK10",
                   showarrow=False, font=dict(color=INK, size=12), align="left")
fig.update_layout(title=dict(text="<b>Panel vs Unsupervised top — overlap 3 genes</b>",
                             font=dict(size=16, color=INK)),
                  xaxis=dict(visible=False, range=[-0.5, 9]),
                  yaxis=dict(visible=False, range=[-0.5, 6.5], scaleanchor="x"),
                  height=560, margin=dict(l=20, r=20, t=70, b=20), **DARK)
fig.write_html(FIG/"v17_8gene_venn.html", include_plotlyjs="cdn", full_html=True)
print("✓ venn")

# ------------------------------------------------------------
# 3. Network
# ------------------------------------------------------------
nodes = {
    "PAX8":   (1.0, 3.0, "#7ccfcd"), "NKX2-1": (1.0, 2.0, "#7ccfcd"),
    "FOXE1":  (1.0, 1.0, "#7ccfcd"),
    "TSHR":   (3.2, 3.4, "#F5A623"), "SLC5A5": (4.8, 2.3, "#E74C3C"),
    "TPO":    (6.4, 3.0, "#F5A623"), "TG":     (6.4, 1.6, "#F5A623"),
    "DIO1":   (8.0, 2.3, "#9b59b6"),
}
edges = [("TSHR","PAX8"),("TSHR","NKX2-1"),("PAX8","SLC5A5"),("PAX8","TPO"),
         ("PAX8","TG"),("NKX2-1","SLC5A5"),("NKX2-1","TPO"),("NKX2-1","TG"),
         ("FOXE1","TPO"),("FOXE1","TG"),("SLC5A5","TPO"),("TPO","TG"),("TG","DIO1")]
edge_x, edge_y = [], []
for a,b in edges:
    edge_x += [nodes[a][0], nodes[b][0], None]
    edge_y += [nodes[a][1], nodes[b][1], None]
fig = go.Figure()
fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines",
                         line=dict(width=2.5, color="rgba(255,255,255,0.4)"),
                         hoverinfo="skip", showlegend=False))
nx_, ny_, nt_, nc_ = [], [], [], []
for n,(x,y,col) in nodes.items():
    nx_.append(x); ny_.append(y); nt_.append(n); nc_.append(col)
fig.add_trace(go.Scatter(x=nx_, y=ny_, mode="markers+text",
                         text=nt_, textposition="middle center",
                         textfont=dict(color="#0a0d12", size=14, family="JetBrains Mono"),
                         marker=dict(size=82, color=nc_,
                                     line=dict(color="#fff", width=3)),
                         hoverinfo="text", showlegend=False))
fig.add_annotation(x=1.0, y=4.0, text="<b>TF trio</b><br>(thyroid masters)",
                   showarrow=False, font=dict(color="#7ccfcd", size=14))
fig.add_annotation(x=3.2, y=4.4, text="<b>TSH signal entry</b>",
                   showarrow=False, font=dict(color="#F5A623", size=14))
fig.add_annotation(x=4.8, y=3.4, text="<b>RAI uptake ★</b><br>(임상 핵심)",
                   showarrow=False, font=dict(color="#E74C3C", size=14))
fig.add_annotation(x=6.4, y=4.1, text="<b>Hormone synthesis</b>",
                   showarrow=False, font=dict(color="#F5A623", size=14))
fig.add_annotation(x=8.0, y=3.4, text="<b>T4→T3 활성화</b>",
                   showarrow=False, font=dict(color="#9b59b6", size=14))
fig.update_layout(title=dict(text="<b>8-gene functional network</b>",
                             font=dict(size=16, color=INK)),
                  xaxis=dict(visible=False, range=[0, 9]),
                  yaxis=dict(visible=False, range=[0.2, 5], scaleanchor="x"),
                  height=500, margin=dict(l=20, r=20, t=60, b=20), **DARK)
fig.write_html(FIG/"v17_8gene_network.html", include_plotlyjs="cdn", full_html=True)
print("✓ network")

# ------------------------------------------------------------
# 4. MAPK→NIS mechanism
# ------------------------------------------------------------
fig = go.Figure()
boxes = [
    (0.2, 1.7, 2.6, 3.0, "Normal thyrocyte<br>NIS↑ TPO↑ TG↑<br><b>RAI uptake 정상</b>", "#2ECC71"),
    (3.6, 1.7, 6.0, 3.0, "BRAF V600E PTC<br>MAPK constitutive ON<br><b>CREB-PAX8 회로 차단</b>", "#F5A623"),
    (7.0, 1.7, 9.4, 3.0, "Dedifferentiated<br>NIS↓ TPO↓ TG↓<br><b>RAI uptake LOST</b>", "#E74C3C"),
]
for x0,y0,x1,y1,lab,col in boxes:
    fig.add_shape(type="rect", x0=x0, y0=y0, x1=x1, y1=y1,
                  line=dict(color=col, width=3), opacity=0.25)
    fig.add_annotation(x=(x0+x1)/2, y=(y0+y1)/2, text=lab,
                       showarrow=False, font=dict(color=col, size=14), align="center")
for ax, x in [(2.6, 3.6), (6.0, 7.0)]:
    fig.add_annotation(x=x, y=2.35, ax=ax, ay=2.35,
                       xref="x", yref="y", axref="x", ayref="y",
                       showarrow=True, arrowhead=3, arrowsize=1.6, arrowwidth=3,
                       arrowcolor="#ccc")
fig.add_shape(type="rect", x0=3.6, y0=0.3, x1=6.0, y1=1.3,
              line=dict(color="#7ccfcd", width=3, dash="dash"), opacity=0.25)
fig.add_annotation(x=4.8, y=0.8,
                   text="<b>MAPK inhibition</b><br>selumetinib · dabrafenib + trametinib",
                   showarrow=False, font=dict(color="#7ccfcd", size=13), align="center")
fig.add_annotation(x=4.8, y=1.55, text="↑ block", showarrow=False,
                   font=dict(color="#7ccfcd", size=13))
fig.add_annotation(x=4.7, y=-0.35,
                   text="<b>References</b>: Knauf 2003 · Chakravarty 2011 JCI · Ho 2013 NEJM · Iravani 2019",
                   showarrow=False, font=dict(color="#bbb", size=12), align="center")
fig.update_layout(title=dict(text="<b>BRAF V600E → MAPK → NIS suppression → RAI loss</b>",
                             font=dict(size=16, color=INK)),
                  xaxis=dict(visible=False, range=[0, 9.6]),
                  yaxis=dict(visible=False, range=[-0.7, 3.5], scaleanchor="x"),
                  height=460, margin=dict(l=20,r=20,t=60,b=40), **DARK)
fig.write_html(FIG/"v17_mapk_nis_mechanism.html", include_plotlyjs="cdn", full_html=True)
print("✓ mechanism")

# ------------------------------------------------------------
# 5. Heatmap
# ------------------------------------------------------------
try:
    expr = np.load("/opt/thyroid-dash/project/data_processed/v5_cross_cancer/THCA/X_combined.npz")["X"]
    with open("/opt/thyroid-dash/project/data_processed/v5_cross_cancer/THCA/shared_genes.txt") as fh:
        genes_all = [g.strip() for g in fh]
    idx = {g:i for i,g in enumerate(genes_all)}
    cols = [g for g in GENES if g in idx]
    mat = expr[:, [idx[g] for g in cols]]
    mu = mat.mean(0); sd = mat.std(0)+1e-9
    z = (mat - mu) / sd
    order = np.argsort(z.mean(axis=1))
    z = z[order][::3][:120]
    fig = go.Figure(go.Heatmap(z=z.T, y=cols, x=list(range(z.shape[0])),
                               colorscale="RdBu_r", zmid=0, zmin=-2.5, zmax=2.5,
                               colorbar=dict(title=dict(text="z-score", font=dict(size=14, color=INK)),
                                             tickfont=dict(size=12, color=INK))))
    fig.update_layout(title=dict(text=f"<b>TCGA-THCA × 8-gene panel z-score (n={z.shape[0]} samples, sorted)</b>",
                                 font=dict(size=15, color=INK)),
                      yaxis=dict(autorange="reversed",
                                 tickfont=dict(size=14, family="JetBrains Mono", color=INK)),
                      xaxis=dict(showticklabels=False,
                                 title=dict(text="← DM1 (low panel) ··· DM2 (high panel) →",
                                            font=dict(size=13, color=INK))),
                      height=440, margin=dict(l=120, r=20, t=80, b=60), **DARK)
    fig.write_html(FIG/"v17_8gene_heatmap.html", include_plotlyjs="cdn", full_html=True)
    print(f"✓ heatmap (n={z.shape[0]})")
except Exception as e:
    print(f"× heatmap: {e}")

# ------------------------------------------------------------
# 6. Cross-analysis matrix
# ------------------------------------------------------------
matrix = pd.read_csv("/opt/thyroid-dash/project/results/v17_tert_recovery/v3/v3_8panel_cross_analysis_matrix.tsv",
                     sep="\t").set_index("panel_gene")
ANALYSES = list(matrix.columns)
hits, labels = [], []
for g in GENES:
    row, lab = [], []
    for a in ANALYSES:
        v = str(matrix.loc[g, a]) if a in matrix.columns and g in matrix.index else ""
        if not v or v == "nan":
            row.append(0); lab.append("")
        elif "NS" in v or "lost" in v:
            row.append(1); lab.append("△")
        else:
            row.append(2); lab.append("✓")
    hits.append(row); labels.append(lab)
fig = go.Figure(go.Heatmap(z=np.array(hits), x=ANALYSES, y=GENES,
                           text=labels, texttemplate="<b>%{text}</b>",
                           textfont=dict(size=18, color=INK),
                           colorscale=[[0, "#1a1f2c"], [0.5, "#F5A623"], [1, "#2ECC71"]],
                           zmin=0, zmax=2, showscale=False))
fig.update_layout(title=dict(text="<b>8 panel × 11 analyses cross-check (✓ = sig/in, △ = NS, blank = absent)</b>",
                             font=dict(size=14, color=INK)),
                  xaxis=dict(tickangle=-30, side="top",
                             tickfont=dict(size=11, color=INK)),
                  yaxis=dict(autorange="reversed",
                             tickfont=dict(size=13, family="JetBrains Mono", color=INK)),
                  height=460, margin=dict(l=80, r=20, t=140, b=20), **DARK)
fig.write_html(FIG/"v17_8gene_cross_analysis.html", include_plotlyjs="cdn", full_html=True)
print("✓ cross_analysis")

# ------------------------------------------------------------
# 7. Panel size sensitivity
# ------------------------------------------------------------
df = pd.read_csv("/opt/thyroid-dash/project/results/v17_realfix/R5_alt_panel_table.tsv", sep="\t")
fig = go.Figure()
fig.add_trace(go.Bar(x=df["panel"], y=df["auc_logreg"],
                     error_y=dict(type="data", symmetric=False,
                                  array=df["ci_hi"]-df["auc_logreg"],
                                  arrayminus=df["auc_logreg"]-df["ci_lo"]),
                     name="LogReg", marker=dict(color="#F5A623"),
                     text=[f"{v:.3f}" for v in df["auc_logreg"]], textposition="outside",
                     textfont=dict(size=13, color=INK)))
fig.add_trace(go.Bar(x=df["panel"], y=df["auc_rf"], name="RandomForest",
                     marker=dict(color="#7ccfcd"),
                     text=[f"{v:.3f}" for v in df["auc_rf"]], textposition="outside",
                     textfont=dict(size=13, color=INK)))
fig.update_layout(title=dict(text="<b>Panel size sensitivity — 8 vs 10 vs 12 vs 16 gene</b>",
                             font=dict(size=15, color=INK)),
                  yaxis=dict(title=dict(text="5-fold CV AUC", font=dict(size=13, color=INK)),
                             range=[0.92, 1.00], gridcolor="rgba(255,255,255,0.08)",
                             tickfont=dict(size=12, color=INK)),
                  xaxis=dict(tickfont=dict(size=12, color=INK)),
                  legend=dict(font=dict(size=12, color=INK)),
                  barmode="group", height=440, margin=dict(l=70, r=20, t=70, b=70), **DARK)
fig.add_annotation(text="<b>8-gene이 16-gene 대비 ΔAUC < 0.013 — marginal cost</b>",
                   xref="paper", yref="paper", x=0.5, y=-0.18, showarrow=False,
                   font=dict(color="#7ccfcd", size=13))
fig.write_html(FIG/"v17_panel_size_sensitivity.html", include_plotlyjs="cdn", full_html=True)
print("✓ panel_size")

# ------------------------------------------------------------
# 8. BRAF positions
# ------------------------------------------------------------
braf_dist = pd.Series([275, 0, 1, 0, 0, 2],
    index=["V600E","V600K","K601E_RAS_like","V600_other","G_loop_RAS_like","other_BRAF"])
fig = go.Figure(go.Bar(x=braf_dist.index, y=braf_dist.values,
                       text=braf_dist.values, textposition="outside",
                       textfont=dict(size=14, color=INK),
                       marker=dict(color=["#F5A623","#c24c4c","#e8b86b","#7ccfcd","#9b59b6","#666"],
                                   line=dict(color="#fff", width=1))))
fig.update_layout(title=dict(text="<b>BRAF mutation position — TCGA-THCA n_patients=278</b>",
                             font=dict(size=15, color=INK)),
                  yaxis=dict(title=dict(text="patients", font=dict(size=13, color=INK)),
                             gridcolor="rgba(255,255,255,0.08)",
                             tickfont=dict(size=12, color=INK)),
                  xaxis=dict(tickfont=dict(size=11, color=INK), tickangle=-15),
                  height=440, margin=dict(l=70,r=20,t=70,b=80), **DARK)
fig.write_html(FIG/"v17_tert_v3_braf_positions.html", include_plotlyjs="cdn", full_html=True)
print("✓ braf_positions")

# ------------------------------------------------------------
# 9. Quad×TERT stack
# ------------------------------------------------------------
sm = pd.read_csv("/opt/thyroid-dash/project/results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv", sep="\t")
sm["tert"] = sm.get("tert_promoter_integrated","wildtype")
ct = pd.crosstab(sm["quad_group"].fillna("missing"), sm["tert"].fillna("wildtype"))
ct = ct.reindex([c for c in ["A_braf_only","B_ras_only","D_triple_negative"] if c in ct.index])
fig = go.Figure()
fig.add_trace(go.Bar(name="TERT wildtype", x=ct.index, y=ct.get("wildtype",0),
                     marker=dict(color="#7ccfcd"), text=ct.get("wildtype",0),
                     textposition="inside", textfont=dict(size=14, color="#0a0d12")))
fig.add_trace(go.Bar(name="TERT mutated", x=ct.index, y=ct.get("mutated",0),
                     marker=dict(color="#E74C3C"), text=ct.get("mutated",0),
                     textposition="inside", textfont=dict(size=14, color="#fff")))
fig.update_layout(barmode="stack",
                  title=dict(text="<b>Quad group × TERT promoter (TCGA-THCA n=513)</b>",
                             font=dict(size=15, color=INK)),
                  yaxis=dict(title=dict(text="patients", font=dict(size=13, color=INK)),
                             gridcolor="rgba(255,255,255,0.08)",
                             tickfont=dict(size=12, color=INK)),
                  xaxis=dict(tickfont=dict(size=12, color=INK)),
                  legend=dict(font=dict(size=12, color=INK)),
                  height=420, margin=dict(l=70,r=20,t=70,b=60), **DARK)
fig.write_html(FIG/"v17_quad_tert_stack.html", include_plotlyjs="cdn", full_html=True)
print("✓ quad_tert_stack")

# ------------------------------------------------------------
# 10. TERT prevalence
# ------------------------------------------------------------
prev = pd.DataFrame([
    ("TCGA-THCA primary PTC", 513, 7.0),
    ("MSK-thyroid 2016 PDTC", 84, 40.5),
    ("MSK-thyroid 2016 ATC",  33, 72.7),
], columns=["cohort","n","prevalence"])
fig = go.Figure(go.Bar(x=prev["cohort"], y=prev["prevalence"],
                       text=[f"{p:.1f}%<br>(n={n})" for p,n in zip(prev["prevalence"], prev["n"])],
                       textposition="outside", textfont=dict(size=14, color=INK),
                       marker=dict(color=prev["prevalence"], colorscale="Reds",
                                   line=dict(color="#fff", width=1))))
fig.update_layout(title=dict(text="<b>TERT promoter prevalence — PTC → PDTC → ATC monotonic ↑</b>",
                             font=dict(size=15, color=INK)),
                  yaxis=dict(title=dict(text="TERT+ %", font=dict(size=13, color=INK)),
                             range=[0, 95], gridcolor="rgba(255,255,255,0.08)",
                             tickfont=dict(size=12, color=INK)),
                  xaxis=dict(tickfont=dict(size=12, color=INK)),
                  height=440, margin=dict(l=70,r=40,t=70,b=80), **DARK)
fig.add_annotation(x=1.0, y=85, xref="x", yref="y",
                   text="<b>탈분화 진행과 monotonic ↑ — R3 trajectory 가설 일치</b>",
                   showarrow=False, font=dict(color="#7ccfcd", size=12))
fig.write_html(FIG/"v17_tert_prevalence.html", include_plotlyjs="cdn", full_html=True)
print("✓ tert_prevalence")

# ------------------------------------------------------------
# 11. Age × TDS × TERT scatter
# ------------------------------------------------------------
sm["tert_bin"] = (sm["tert"]=="mutated").astype(int)
sub = sm.dropna(subset=["age","os_event"]).copy()
sub["os_label"] = sub["os_event"].apply(lambda x: "event" if x==1 else "censored")
fig = px.scatter(sub, x="age", y="tds16_score_v17",
                 color="tert", symbol="os_label",
                 color_discrete_map={"mutated":"#E74C3C","wildtype":"#7ccfcd"},
                 symbol_map={"event":"x","censored":"circle"},
                 title=f"<b>Age × TDS16 × TERT (n={len(sub)}, OS events={int(sub['os_event'].sum())})</b>")
fig.update_traces(marker=dict(size=10, opacity=0.75, line=dict(width=0.6, color="#fff")))
fig.update_layout(template="plotly_dark", paper_bgcolor=BG, plot_bgcolor=BG,
                  font=dict(color=INK, size=13), height=440,
                  margin=dict(l=70,r=20,t=70,b=60),
                  xaxis=dict(title=dict(text="Age at diagnosis", font=dict(size=13, color=INK)),
                             gridcolor="rgba(255,255,255,0.08)",
                             tickfont=dict(size=12, color=INK)),
                  yaxis=dict(title=dict(text="TDS16 score", font=dict(size=13, color=INK)),
                             gridcolor="rgba(255,255,255,0.08)",
                             tickfont=dict(size=12, color=INK)),
                  legend=dict(font=dict(size=12, color=INK)))
fig.write_html(FIG/"v17_age_tert_tds.html", include_plotlyjs="cdn", full_html=True)
print("✓ age_tert_tds")

# ------------------------------------------------------------
# 12. KM curves (Quad + Quint)
# ------------------------------------------------------------
def km_curve(times, events):
    df = pd.DataFrame({"t":times, "e":events}).sort_values("t").reset_index(drop=True)
    n = len(df); s = 1.0
    ts, ss = [0.0], [1.0]; at_risk = n
    for _, r in df.iterrows():
        if r["e"]==1 and at_risk>0:
            s *= (1 - 1/at_risk)
        ts.append(r["t"]); ss.append(s); at_risk -= 1
    return ts, ss

m = pd.read_csv("/opt/thyroid-dash/project/results/v17_tert_recovery/v3/v3_merged_quint.tsv", sep="\t")
m = m.dropna(subset=["os_days","os_event"])
fig = make_subplots(rows=1, cols=2,
                    subplot_titles=("<b>Quad group KM</b>", "<b>BRAF×TERT split KM</b>"))
quad_pal = {"A_braf_only":"#F5A623","B_ras_only":"#7ccfcd","D_triple_negative":"#9b59b6"}
for g,col in quad_pal.items():
    sub = m[m["quad_group"]==g]
    if len(sub)<3: continue
    t,s = km_curve(sub["os_days"].values, sub["os_event"].values)
    fig.add_trace(go.Scatter(x=t,y=s,mode="lines",
        name=f"{g} (n={len(sub)})", line=dict(color=col,width=2.5),
        legendgroup="quad"), row=1, col=1)
quint_pal = {
    "BRAF_V600E·TERT-":"#F5A623","BRAF_V600E·TERT+":"#c24c4c",
    "RAS·TERT-":"#7ccfcd","RAS·TERT+":"#3498db",
    "TripleNeg·TERT-":"#9b59b6","TripleNeg·TERT+":"#6c2ba0",
}
for g,col in quint_pal.items():
    sub = m[m["quint_group"]==g]
    if len(sub)<3: continue
    t,s = km_curve(sub["os_days"].values, sub["os_event"].values)
    dash = "dash" if "TERT+" in g else "solid"
    fig.add_trace(go.Scatter(x=t,y=s,mode="lines",
        name=f"{g} (n={len(sub)},ev={int(sub['os_event'].sum())})",
        line=dict(color=col, width=2.5, dash=dash), legendgroup="quint"), row=1, col=2)
for c in [1,2]:
    fig.update_xaxes(title=dict(text="Days", font=dict(size=13, color=INK)),
                     tickfont=dict(size=11, color=INK), gridcolor="rgba(255,255,255,0.08)", row=1,col=c)
fig.update_yaxes(title=dict(text="OS probability", font=dict(size=13, color=INK)),
                 range=[0.5,1.02], tickfont=dict(size=11, color=INK),
                 gridcolor="rgba(255,255,255,0.08)", row=1,col=1)
fig.update_yaxes(range=[0.5,1.02], tickfont=dict(size=11, color=INK),
                 gridcolor="rgba(255,255,255,0.08)", row=1,col=2)
fig.update_layout(title=dict(text=f"<b>v17 TERT Survival — TCGA-THCA n={len(m)}, ev={int(m['os_event'].sum())}, BRAF×TERT logrank p=2.09×10⁻⁵</b>",
                             font=dict(size=14, color=INK)),
                  height=540, margin=dict(t=110,l=70,r=20,b=70),
                  legend=dict(font=dict(size=10, color=INK), tracegroupgap=8),
                  **DARK)
fig.update_annotations(font=dict(size=14, color=INK))  # subplot titles
fig.write_html(FIG/"v17_tert_v3_KM.html", include_plotlyjs="cdn", full_html=True)
print("✓ KM")

print("\nALL FIGURES RE-RENDERED WITH SOLID DARK BG + READABLE TEXT")
