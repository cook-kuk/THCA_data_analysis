#!/usr/bin/env python3
"""arcasHLA FINAL polish — 12 publication-quality figures + statistical tables.

Designed for full PRJEB11591 cohort (n≈260) but works on any subset.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats

RESULTS_DIR = Path("/opt/thyroid-dash/project/results/v17_korean/arcasHLA")
FIG_DIR = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17/hla_mega")
PRED = Path("/opt/thyroid-dash/project/results/v17_korean/K2_korean_predictions_v4.tsv")
MANIFEST = Path("/opt/thyroid-dash/project/results/v17_korean/K1A_prjeb11591_runs.tsv")

GENES = ["A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
RISK_ALLELES = {
    "DPB1*05:01": ("Korean Graves' risk ★", "#dc2626"),
    "B*46:01":    ("Asian Graves' risk ★", "#f59e0b"),
    "DRB1*04:01": ("Hashimoto/Graves' overlap", "#9b59b6"),
    "DRB1*15:01": ("Protective candidate", "#2ECC71"),
    "DRB1*03:01": ("Caucasian Graves' (rare in Asian)", "#7ccfcd"),
    "B*08:01":    ("Caucasian DR3 LD (rare in Asian)", "#94a3b8"),
}

# Published frequencies (Kim 2014 PLOS One n=413 Korean + 1000G EUR)
LIT_AF = {
    "DPB1*05:01": {"Korean": 0.180, "EUR": 0.060},
    "B*46:01":    {"Korean": 0.045, "EUR": 0.001},
    "DRB1*04:01": {"Korean": 0.018, "EUR": 0.075},
    "DRB1*15:01": {"Korean": 0.045, "EUR": 0.130},
    "DRB1*03:01": {"Korean": 0.003, "EUR": 0.090},
    "B*08:01":    {"Korean": 0.000, "EUR": 0.085},
}

COMMON = dict(include_plotlyjs="cdn",
              config={"displaylogo": False, "toImageButtonOptions": {"format":"png","scale":2}})


# ============= Data loading =============

def load_genotypes():
    rows = []
    for jp in sorted(RESULTS_DIR.glob("*.genotype.json")):
        run = jp.stem.replace(".genotype", "").rstrip("_1")
        try:
            d = json.loads(jp.read_text())
        except Exception:
            continue
        if not d or not any(g in d for g in GENES):
            continue
        row = {"run": run}
        for g in GENES:
            alleles = d.get(g, [])
            for i, a in enumerate(alleles[:2], 1):
                row[f"{g}_a{i}"] = a
                if a:
                    parts = a.split(":")
                    row[f"{g}_a{i}_4d"] = ":".join(parts[:2]) if len(parts) >= 2 else a
        rows.append(row)
    return pd.DataFrame(rows)


def merge_dm(geno):
    pred = pd.read_csv(PRED, sep="\t")[["run", "p_DM2", "DM_call"]]
    manifest = pd.read_csv(MANIFEST, sep="\t")[["run_accession", "sample_alias"]]
    df = geno.merge(pred, on="run").merge(manifest, left_on="run", right_on="run_accession", how="left")
    df["is_normal"] = df["sample_alias"].str.endswith("-N", na=False)

    def parse_cat(a):
        import re
        if not isinstance(a, str): return "unknown"
        if a.endswith("-N"): return "Matched normal"
        m = re.match(r"SNU-GMI-([A-Z]+)\d+", a)
        if not m: return "unknown"
        return {"PTC":"PTC", "FV":"FVPTC", "FTC":"FTC", "FA":"FA"}.get(m.group(1), m.group(1))
    df["histology"] = df["sample_alias"].apply(parse_cat)
    return df


def is_carrier(df, allele):
    gene = allele.split("*")[0]
    target = f"{gene}*{allele.split('*')[1]}"
    cols = [f"{gene}_a1_4d", f"{gene}_a2_4d"]
    cols = [c for c in cols if c in df.columns]
    return df.apply(lambda r: any(isinstance(r[c], str) and r[c] == target for c in cols), axis=1)


# ============= Statistical helper =============

def fisher_test(a, b, c, d):
    """Fisher's exact 2x2: [a, b] / [c, d]"""
    try:
        odds, p = stats.fisher_exact([[a, b], [c, d]])
    except Exception:
        return None, None
    return odds, p


def wilson_ci(k, n, alpha=0.05):
    """Wilson score 95% CI."""
    if n == 0: return (0, 0)
    z = stats.norm.ppf(1 - alpha/2)
    p = k / n
    den = 1 + z**2/n
    cen = (p + z**2/(2*n)) / den
    half = (z * np.sqrt(p*(1-p)/n + z**2/(4*n**2))) / den
    return (max(0, cen - half), min(1, cen + half))


def cohen_h(p1, p2):
    return 2 * (np.arcsin(np.sqrt(p1)) - np.arcsin(np.sqrt(p2)))


# ============= Figures (12) =============

C_DM1 = "#1f77b4"; C_DM2 = "#dc2626"; C_NORM = "#10b981"


def fig01_headline(df):
    """Top metrics tile board."""
    n = len(df)
    n_dm1 = int((df["DM_call"] == "DM1").sum())
    n_dm2 = int((df["DM_call"] == "DM2").sum())
    n_norm = int(df["is_normal"].sum())
    n_dpb1 = int(is_carrier(df, "DPB1*05:01").sum())
    n_b46 = int(is_carrier(df, "B*46:01").sum())
    n_drb1_15 = int(is_carrier(df, "DRB1*15:01").sum())

    metrics = [
        (f"{n}", "★ valid HLA genotypes"),
        (f"{n_dm1}/{n_dm2}", "DM1 / DM2 split"),
        (f"{n_norm}", "matched normal"),
        (f"{n_dpb1} ({100*n_dpb1/n:.1f}%)", "DPB1*05:01 carrier ★"),
        (f"{n_b46} ({100*n_b46/n:.1f}%)", "B*46:01 carrier ★"),
        (f"{n_drb1_15}", "DRB1*15:01 (protective)"),
    ]
    fig = go.Figure()
    for i, (val, lab) in enumerate(metrics):
        fig.add_annotation(x=i, y=0.6, text=f"<b style='font-size:32px;color:#dc2626'>{val}</b>",
                           showarrow=False, font=dict(size=28, color="#dc2626"))
        fig.add_annotation(x=i, y=0.15, text=lab, showarrow=False,
                           font=dict(size=11, color="#666"))
    fig.update_xaxes(visible=False, range=[-0.5, len(metrics)-0.5])
    fig.update_yaxes(visible=False, range=[0, 1])
    fig.update_layout(
        title=f"<b>arcasHLA Korean K2 (PRJEB11591) — Final results n={n}</b>",
        height=200, paper_bgcolor="white", plot_bgcolor="white",
        margin=dict(l=20, r=20, t=60, b=10),
    )
    out = FIG_DIR / "FINAL_01_headline.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig02_carrier_forest(df):
    """6 risk allele × overall + Korean ref + EUR ref forest plot with 95% CI."""
    n = len(df)
    rows = []
    for allele, (label, color) in RISK_ALLELES.items():
        car = is_carrier(df, allele)
        k = int(car.sum())
        af_ours = k / (2 * n)  # rough — assume mostly heterozygous
        car_pct = k / n
        ci_lo, ci_hi = wilson_ci(k, n)
        rows.append(dict(allele=allele, label=label, color=color,
                         k=k, n=n, car_pct=car_pct, ci_lo=ci_lo, ci_hi=ci_hi,
                         af_ours=af_ours,
                         af_kor=LIT_AF[allele]["Korean"],
                         af_eur=LIT_AF[allele]["EUR"]))
    df_ref = pd.DataFrame(rows)

    fig = go.Figure()
    y_pos = list(range(len(df_ref)))
    # Our cohort (red)
    for i, r in df_ref.iterrows():
        fig.add_trace(go.Scatter(
            x=[r["car_pct"]*100], y=[i],
            error_x=dict(type="data", symmetric=False,
                         array=[(r["ci_hi"]-r["car_pct"])*100],
                         arrayminus=[(r["car_pct"]-r["ci_lo"])*100]),
            mode="markers", marker=dict(size=14, color=r["color"], symbol="diamond",
                                        line=dict(width=2, color="white")),
            name=f"Ours ({r['allele']})", showlegend=False,
            hovertemplate=f"<b>{r['allele']}</b><br>Ours: {r['car_pct']*100:.1f}% [{r['ci_lo']*100:.1f}-{r['ci_hi']*100:.1f}]<br>k={r['k']}/{r['n']}<extra></extra>",
        ))
    # Kim 2014 Korean ref (orange)
    fig.add_trace(go.Scatter(
        x=[r["af_kor"]*200 for _, r in df_ref.iterrows()],  # AF×2 ≈ carrier rate for rare
        y=y_pos, mode="markers",
        marker=dict(size=12, color="#f59e0b", symbol="circle"),
        name="Kim 2014 Korean ref (n=413)",
        hovertemplate="Kim 2014 carrier: %{x:.1f}%<extra></extra>",
    ))
    # 1000G EUR (teal)
    fig.add_trace(go.Scatter(
        x=[r["af_eur"]*200 for _, r in df_ref.iterrows()],
        y=y_pos, mode="markers",
        marker=dict(size=12, color="#7ccfcd", symbol="circle"),
        name="1000G EUR ref",
        hovertemplate="1000G EUR carrier: %{x:.1f}%<extra></extra>",
    ))

    fig.update_layout(
        title=f"<b>Fig FINAL-2 — 6 Graves' risk allele carrier rate forest (Korean K2 n={n} + Kim 2014 + 1000G EUR)</b>",
        xaxis_title="carrier rate (%)",
        yaxis=dict(tickmode="array", tickvals=y_pos,
                   ticktext=[f"<b>{r['allele']}</b><br><sub>{r['label']}</sub>" for _,r in df_ref.iterrows()]),
        height=520, paper_bgcolor="white", plot_bgcolor="white",
        font=dict(family="-apple-system,Segoe UI,sans-serif", size=11),
        legend=dict(orientation="h", yanchor="top", y=1.10, x=0.5, xanchor="center"),
    )
    out = FIG_DIR / "FINAL_02_carrier_forest.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig03_dm_split_table(df):
    """Beautiful Plotly Table — risk allele × DM1/DM2 with Fisher test."""
    n = len(df)
    n_dm1 = int((df["DM_call"] == "DM1").sum())
    n_dm2 = int((df["DM_call"] == "DM2").sum())
    rows = []
    for allele, (label, _) in RISK_ALLELES.items():
        car = is_carrier(df, allele)
        in_dm1 = int(car[df["DM_call"]=="DM1"].sum())
        in_dm2 = int(car[df["DM_call"]=="DM2"].sum())
        non_dm1 = n_dm1 - in_dm1
        non_dm2 = n_dm2 - in_dm2
        odds, p = fisher_test(in_dm1, non_dm1, in_dm2, non_dm2)
        h = cohen_h(in_dm1/max(n_dm1,1), in_dm2/max(n_dm2,1))
        sig = "★" if p is not None and p < 0.05 else ("·" if p is not None and p < 0.1 else "")
        rows.append([
            f"<b>{allele}</b>",
            label,
            f"{int(car.sum())}/{n}",
            f"{in_dm1}/{n_dm1} ({100*in_dm1/max(n_dm1,1):.1f}%)",
            f"{in_dm2}/{n_dm2} ({100*in_dm2/max(n_dm2,1):.1f}%)",
            f"{odds:.2f}" if odds is not None else "—",
            f"{h:+.2f}",
            f"{p:.3f} {sig}" if p is not None else "—",
        ])
    cols = ["Allele", "Label", "Overall<br>carrier", "DM1<br>carrier", "DM2<br>carrier", "Fisher<br>OR", "Cohen's<br>h", "p-value"]
    fill = []
    for r in rows:
        if "★" in r[7]: fill.append("rgba(220,38,38,0.16)")
        elif "·" in r[7]: fill.append("rgba(245,158,11,0.10)")
        else: fill.append("rgba(124,207,205,0.06)")
    full_fill = [fill] * len(cols)
    fig = go.Figure(data=[go.Table(
        columnwidth=[80, 200, 80, 110, 110, 60, 60, 90],
        header=dict(values=[f"<b>{c}</b>" for c in cols], fill_color="#1e3a8a",
                    font=dict(color="white", size=12), align="left", height=44),
        cells=dict(values=list(zip(*rows)), fill_color=full_fill,
                   font=dict(color="#1f2937", size=11), align="left", height=36),
    )])
    fig.update_layout(
        title=f"<b>Fig FINAL-3 — DM1 vs DM2 risk allele cross-tab + Fisher exact + Cohen's h (n={n}, DM1={n_dm1}, DM2={n_dm2})</b>",
        height=460, margin=dict(l=20, r=20, t=80, b=20), paper_bgcolor="white",
        font=dict(family="-apple-system,Segoe UI,Noto Sans KR,sans-serif", size=11),
    )
    out = FIG_DIR / "FINAL_03_dm_table.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig04_histology_dpb1(df):
    """DPB1*05:01 × histology stacked bar."""
    target = "DPB1*05:01"
    car = is_carrier(df, target)
    df = df.copy()
    df["DPB1_carr"] = car
    cats = ["PTC", "FVPTC", "FTC", "FA", "Matched normal"]
    rows = []
    for c in cats:
        sub = df[df["histology"] == c]
        if len(sub) == 0: continue
        n = len(sub); k = int(sub["DPB1_carr"].sum())
        ci = wilson_ci(k, n)
        rows.append(dict(cat=c, n=n, k=k, pct=100*k/n, ci_lo=100*ci[0], ci_hi=100*ci[1]))
    rdf = pd.DataFrame(rows)
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=rdf["cat"], y=rdf["pct"],
        error_y=dict(type="data", symmetric=False,
                     array=rdf["ci_hi"]-rdf["pct"], arrayminus=rdf["pct"]-rdf["ci_lo"]),
        marker_color=["#3b82f6","#7ccfcd","#9b59b6","#f59e0b","#10b981"][:len(rdf)],
        text=[f"{r['k']}/{r['n']}<br><b>{r['pct']:.1f}%</b>" for _,r in rdf.iterrows()],
        textposition="outside", hovertemplate="<b>%{x}</b><br>%{y:.1f}% [%{error_y.arrayminus:.1f}, %{error_y.array:.1f}]<extra></extra>",
    ))
    fig.update_layout(
        title=f"<b>Fig FINAL-4 — DPB1*05:01 carrier × histology with Wilson 95% CI (Korean K2 n={len(df)})</b>",
        xaxis_title="histology category", yaxis_title="DPB1*05:01 carrier (%)",
        height=460, paper_bgcolor="white", plot_bgcolor="white",
        font=dict(family="-apple-system,Segoe UI,Noto Sans KR,sans-serif", size=11),
    )
    out = FIG_DIR / "FINAL_04_histology_dpb1.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig05_allele_frequency_panel(df):
    """6-panel grid: top alleles per HLA gene."""
    fig = make_subplots(rows=2, cols=3,
                        subplot_titles=[f"<b>HLA-{g}</b>" for g in ["A","B","C","DRB1","DQB1","DPB1"]],
                        vertical_spacing=0.18, horizontal_spacing=0.08)
    pos = [(1,1),(1,2),(1,3),(2,1),(2,2),(2,3)]
    for (gene, (r,c)) in zip(["A","B","C","DRB1","DQB1","DPB1"], pos):
        cols = [f"{gene}_a1_4d", f"{gene}_a2_4d"]
        cols = [col for col in cols if col in df.columns]
        if not cols: continue
        all_a = pd.concat([df[col].dropna() for col in cols])
        top = all_a.value_counts().head(8)
        # Highlight risk alleles
        colors = ["#dc2626" if a in RISK_ALLELES else "#3b82f6" for a in top.index]
        fig.add_trace(go.Bar(
            x=top.values, y=top.index, orientation="h",
            marker_color=colors, showlegend=False,
            text=[f"{v} ({100*v/(2*len(df)):.1f}%)" for v in top.values],
            textposition="outside",
            hovertemplate="<b>%{y}</b><br>n=%{x}<extra></extra>",
        ), row=r, col=c)
        fig.update_yaxes(autorange="reversed", row=r, col=c)
    fig.update_layout(
        title=f"<b>Fig FINAL-5 — Top 8 alleles per HLA gene (Korean K2 n={len(df)}, AF=count/(2×n)) — red = Graves' risk</b>",
        height=720, paper_bgcolor="white", plot_bgcolor="white",
        font=dict(family="-apple-system,Segoe UI,sans-serif", size=10),
    )
    out = FIG_DIR / "FINAL_05_allele_freq_grid.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig06_genotype_heatmap(df):
    """Sample × gene heatmap (DM-sorted, top alleles)."""
    df_sort = df.sort_values("p_DM2", ascending=False).reset_index(drop=True)
    genes = ["A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1"]
    # Top 6 alleles per gene as the cells
    cols, labels = [], []
    for g in genes:
        a1 = f"{g}_a1_4d"; a2 = f"{g}_a2_4d"
        if a1 not in df.columns: continue
        cols.append(a1); labels.append(f"{g} (a1)")
        cols.append(a2); labels.append(f"{g} (a2)")
    # Build matrix: encode each cell by top-allele rank within that gene
    Z = np.zeros((len(df_sort), len(cols)))
    for j, col in enumerate(cols):
        vc = df_sort[col].value_counts()
        rank = {a: i+1 for i, a in enumerate(vc.index[:30])}
        Z[:, j] = [rank.get(v, 31) if isinstance(v, str) else 0 for v in df_sort[col]]
    text = [[f"{df_sort.iloc[i][col]}" if isinstance(df_sort.iloc[i][col], str) else "—"
             for col in cols] for i in range(len(df_sort))]
    fig = go.Figure(go.Heatmap(
        z=Z, y=df_sort["run"], x=labels,
        colorscale=[[0,"#fef2f2"],[0.05,"#fee2e2"],[0.2,"#fca5a5"],[0.5,"#7ccfcd"],[1,"#1e3a8a"]],
        showscale=False, text=text, hovertemplate="<b>%{y}</b><br>%{x}: %{text}<extra></extra>",
    ))
    fig.update_layout(
        title=f"<b>Fig FINAL-6 — Per-sample HLA genotype matrix (Korean K2 n={len(df_sort)}, sorted by p(DM2) descending) — color = allele rarity rank</b>",
        height=max(600, 12*len(df_sort)),
        paper_bgcolor="white", plot_bgcolor="white",
        xaxis=dict(side="top"),
    )
    out = FIG_DIR / "FINAL_06_genotype_heatmap.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig07_dpb1_violin(df):
    """DPB1*05:01 carrier × p(DM2) distribution + Mann-Whitney p."""
    car = is_carrier(df, "DPB1*05:01")
    yes = df[car]["p_DM2"]; no = df[~car]["p_DM2"]
    if len(yes) > 0 and len(no) > 0:
        u, p_mw = stats.mannwhitneyu(yes, no, alternative="two-sided")
    else:
        p_mw = None
    fig = go.Figure()
    fig.add_trace(go.Violin(y=yes, name=f"Carrier (n={len(yes)})", marker_color="#dc2626",
                            box_visible=True, meanline_visible=True, points="all", pointpos=-1.1))
    fig.add_trace(go.Violin(y=no, name=f"Non-carrier (n={len(no)})", marker_color="#7ccfcd",
                            box_visible=True, meanline_visible=True, points="all", pointpos=-1.1))
    fig.add_hline(y=0.5, line_dash="dash", line_color="#444",
                  annotation_text="DM1/DM2 cut", annotation_position="right")
    p_str = f"p = {p_mw:.4f}" if p_mw is not None else "p = —"
    fig.update_layout(
        title=f"<b>Fig FINAL-7 — DPB1*05:01 carrier × p(DM2) distribution (Mann-Whitney {p_str})</b>",
        yaxis_title="predicted p(DM2)", yaxis=dict(range=[-0.05, 1.05]),
        height=480, paper_bgcolor="white", plot_bgcolor="white",
    )
    out = FIG_DIR / "FINAL_07_dpb1_violin.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig08_korean_vs_lit(df):
    """Group bar with our cohort (red) vs Kim 2014 Korean (orange) vs 1000G EUR (teal)."""
    n = len(df)
    rows = []
    for allele in RISK_ALLELES:
        k = int(is_carrier(df, allele).sum())
        ours_af = k / (2 * n)  # AF approximation
        rows.append(dict(allele=allele,
                         ours=ours_af*100,
                         kor=LIT_AF[allele]["Korean"]*100,
                         eur=LIT_AF[allele]["EUR"]*100))
    rdf = pd.DataFrame(rows)
    fig = go.Figure()
    fig.add_trace(go.Bar(x=rdf["allele"], y=rdf["kor"], name="Kim 2014 Korean (n=413)", marker_color="#f59e0b",
                         hovertemplate="Kim 2014 AF: %{y:.1f}%<extra></extra>"))
    fig.add_trace(go.Bar(x=rdf["allele"], y=rdf["eur"], name="1000G EUR ref", marker_color="#7ccfcd",
                         hovertemplate="1000G EUR AF: %{y:.1f}%<extra></extra>"))
    fig.add_trace(go.Bar(x=rdf["allele"], y=rdf["ours"], name=f"Ours K2 (n={n})", marker_color="#dc2626",
                         text=[f"{v:.1f}%" for v in rdf["ours"]], textposition="outside",
                         hovertemplate="Our AF: %{y:.1f}%<extra></extra>"))
    fig.update_layout(barmode="group",
        title=f"<b>Fig FINAL-8 — Allele frequency replication: Ours (K2 n={n}) vs Kim 2014 (Korean n=413) vs 1000G EUR</b>",
        yaxis_title="allele frequency (%)",
        height=480, paper_bgcolor="white", plot_bgcolor="white",
        font=dict(family="-apple-system,Segoe UI,Noto Sans KR,sans-serif", size=11),
    )
    out = FIG_DIR / "FINAL_08_lit_replication.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig09_zygosity(df):
    """Risk allele zygosity (homo/hetero/non-carrier) stacked."""
    n = len(df)
    rows = []
    for allele, (label, _) in RISK_ALLELES.items():
        gene = allele.split("*")[0]
        target = f"{gene}*{allele.split('*')[1]}"
        a1 = f"{gene}_a1_4d"; a2 = f"{gene}_a2_4d"
        a1_match = df[a1].apply(lambda x: isinstance(x, str) and x == target) if a1 in df else pd.Series([False]*n)
        a2_match = df[a2].apply(lambda x: isinstance(x, str) and x == target) if a2 in df else pd.Series([False]*n)
        homo = int((a1_match & a2_match).sum())
        het = int((a1_match ^ a2_match).sum())
        non = n - homo - het
        rows.append(dict(allele=allele, label=label, homo=homo, het=het, non=non,
                         pct_homo=100*homo/n, pct_het=100*het/n, pct_non=100*non/n))
    rdf = pd.DataFrame(rows)
    fig = go.Figure()
    fig.add_trace(go.Bar(y=rdf["allele"], x=rdf["pct_homo"], orientation="h", name="Homozygous",
                         marker_color="#dc2626",
                         text=[f"{v}" if v>0 else "" for v in rdf["homo"]], textposition="inside",
                         hovertemplate="%{y}<br>Homo: %{x:.1f}%<extra></extra>"))
    fig.add_trace(go.Bar(y=rdf["allele"], x=rdf["pct_het"], orientation="h", name="Heterozygous",
                         marker_color="#f59e0b",
                         text=[f"{v}" if v>0 else "" for v in rdf["het"]], textposition="inside",
                         hovertemplate="%{y}<br>Het: %{x:.1f}%<extra></extra>"))
    fig.add_trace(go.Bar(y=rdf["allele"], x=rdf["pct_non"], orientation="h", name="Non-carrier",
                         marker_color="#7ccfcd",
                         hovertemplate="%{y}<br>Non: %{x:.1f}%<extra></extra>"))
    fig.update_layout(barmode="stack",
        title=f"<b>Fig FINAL-9 — Risk allele zygosity breakdown (Korean K2 n={n})</b>",
        xaxis_title="proportion (%)",
        height=440, paper_bgcolor="white", plot_bgcolor="white",
        font=dict(family="-apple-system,Segoe UI,sans-serif", size=11),
    )
    out = FIG_DIR / "FINAL_09_zygosity.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig10_per_sample_table(df):
    """Beautiful per-sample HLA table — sortable, colored by risk allele carriers."""
    df_show = df.sort_values("p_DM2", ascending=False).copy()
    rows = []
    for _, r in df_show.iterrows():
        risk_flags = []
        for a in RISK_ALLELES:
            if is_carrier(pd.DataFrame([r]), a).iloc[0]:
                risk_flags.append(a)
        rows.append([
            r["run"],
            r.get("sample_alias", ""),
            r["histology"],
            r["DM_call"],
            f"{r['p_DM2']:.3f}",
            f"{r.get('A_a1_4d','')} / {r.get('A_a2_4d','')}",
            f"{r.get('B_a1_4d','')} / {r.get('B_a2_4d','')}",
            f"{r.get('C_a1_4d','')} / {r.get('C_a2_4d','')}",
            f"{r.get('DRB1_a1_4d','')} / {r.get('DRB1_a2_4d','')}",
            f"{r.get('DQB1_a1_4d','')} / {r.get('DQB1_a2_4d','')}",
            f"{r.get('DPB1_a1_4d','')} / {r.get('DPB1_a2_4d','')}",
            ", ".join(risk_flags) if risk_flags else "",
        ])
    cols = ["Run", "Sample", "Histology", "DM call", "p(DM2)",
            "HLA-A", "HLA-B", "HLA-C", "DRB1", "DQB1", "DPB1", "Risk alleles"]
    fill = ["rgba(220,38,38,0.10)" if "DPB1*05:01" in r[-1] or "B*46:01" in r[-1] else
            ("rgba(46,204,113,0.06)" if "DRB1*15:01" in r[-1] else "white")
            for r in rows]
    full_fill = [fill] * len(cols)
    fig = go.Figure(data=[go.Table(
        columnwidth=[100,120,80,60,60,140,140,140,140,140,140,180],
        header=dict(values=[f"<b>{c}</b>" for c in cols], fill_color="#1e3a8a",
                    font=dict(color="white", size=11), align="left", height=32),
        cells=dict(values=list(zip(*rows)), fill_color=full_fill,
                   font=dict(color="#1f2937", size=10), align="left", height=22),
    )])
    fig.update_layout(
        title=f"<b>Fig FINAL-10 — Per-sample HLA imputation table (n={len(df_show)}, sorted by p(DM2) descending)</b>",
        height=min(900, 50 + 24*len(df_show)),
        margin=dict(l=20, r=20, t=70, b=20), paper_bgcolor="white",
        font=dict(family="-apple-system,Segoe UI,sans-serif", size=10),
    )
    out = FIG_DIR / "FINAL_10_per_sample.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig11_dm1_vs_dm2_per_gene(df):
    """For each HLA gene, top alleles split by DM1/DM2."""
    fig = make_subplots(rows=2, cols=3,
                        subplot_titles=[f"<b>HLA-{g}</b>" for g in ["A","B","C","DRB1","DQB1","DPB1"]],
                        vertical_spacing=0.16, horizontal_spacing=0.08)
    pos = [(1,1),(1,2),(1,3),(2,1),(2,2),(2,3)]
    for (gene, (r,c)) in zip(["A","B","C","DRB1","DQB1","DPB1"], pos):
        a1 = f"{gene}_a1_4d"; a2 = f"{gene}_a2_4d"
        if a1 not in df: continue
        all_a = pd.concat([df[a1], df[a2]])
        top = all_a.value_counts().head(6)
        dm1_counts = []; dm2_counts = []
        for a in top.index:
            d_dm1 = df[df["DM_call"]=="DM1"]
            d_dm2 = df[df["DM_call"]=="DM2"]
            n1 = int((d_dm1[a1].eq(a) | d_dm1[a2].eq(a)).sum())
            n2 = int((d_dm2[a1].eq(a) | d_dm2[a2].eq(a)).sum())
            dm1_counts.append(n1); dm2_counts.append(n2)
        fig.add_trace(go.Bar(x=top.index, y=dm1_counts, name="DM1", marker_color=C_DM1,
                             showlegend=(r==1 and c==1),
                             hovertemplate="%{x}<br>DM1 carrier: %{y}<extra></extra>"), row=r, col=c)
        fig.add_trace(go.Bar(x=top.index, y=dm2_counts, name="DM2", marker_color=C_DM2,
                             showlegend=(r==1 and c==1),
                             hovertemplate="%{x}<br>DM2 carrier: %{y}<extra></extra>"), row=r, col=c)
    fig.update_layout(barmode="group",
        title=f"<b>Fig FINAL-11 — Top alleles × DM1/DM2 split per HLA gene (Korean K2 n={len(df)})</b>",
        height=720, paper_bgcolor="white", plot_bgcolor="white",
        font=dict(family="-apple-system,Segoe UI,sans-serif", size=10),
    )
    out = FIG_DIR / "FINAL_11_dm_per_gene.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig12_summary_dashboard(df):
    """One-page dashboard combining all key findings."""
    n = len(df)
    n_dm1 = int((df["DM_call"]=="DM1").sum())
    n_dm2 = int((df["DM_call"]=="DM2").sum())
    n_norm = int(df["is_normal"].sum())

    rows = []
    for allele, (label, color) in RISK_ALLELES.items():
        car = is_carrier(df, allele)
        k = int(car.sum())
        in_dm1 = int(car[df["DM_call"]=="DM1"].sum())
        in_dm2 = int(car[df["DM_call"]=="DM2"].sum())
        odds, p = fisher_test(in_dm1, n_dm1-in_dm1, in_dm2, n_dm2-in_dm2)
        rows.append([
            f"<b>{allele}</b>", label, f"{k}/{n} ({100*k/n:.1f}%)",
            f"{in_dm1}/{n_dm1} ({100*in_dm1/max(n_dm1,1):.1f}%)",
            f"{in_dm2}/{n_dm2} ({100*in_dm2/max(n_dm2,1):.1f}%)",
            f"{p:.3f}" if p is not None else "—",
        ])
    fig = make_subplots(rows=2, cols=1, row_heights=[0.4, 0.6],
                        specs=[[{"type":"table"}], [{"type":"bar"}]],
                        subplot_titles=("<b>Summary table — 6 risk allele × DM1/DM2 (Fisher)</b>",
                                        "<b>Carrier rate forest with Wilson 95% CI</b>"))
    fill = []
    for r in rows:
        try:
            p_v = float(r[5]) if r[5] != "—" else 1.0
        except: p_v = 1.0
        fill.append("rgba(220,38,38,0.16)" if p_v < 0.05 else "rgba(124,207,205,0.06)")
    full_fill = [fill] * 6
    fig.add_trace(go.Table(
        columnwidth=[80,180,100,110,110,80],
        header=dict(values=["<b>Allele</b>","<b>Label</b>","<b>Overall</b>","<b>DM1</b>","<b>DM2</b>","<b>Fisher p</b>"],
                    fill_color="#1e3a8a", font=dict(color="white", size=11), align="left", height=34),
        cells=dict(values=list(zip(*rows)), fill_color=full_fill,
                   font=dict(color="#1f2937", size=10.5), align="left", height=30),
    ), row=1, col=1)

    # Bottom: forest
    pos = list(range(len(RISK_ALLELES)))
    for i, (allele, (label, color)) in enumerate(RISK_ALLELES.items()):
        car = is_carrier(df, allele)
        k = int(car.sum())
        ci = wilson_ci(k, n)
        fig.add_trace(go.Scatter(x=[100*k/n], y=[i],
            error_x=dict(type="data", symmetric=False,
                         array=[100*ci[1] - 100*k/n], arrayminus=[100*k/n - 100*ci[0]]),
            mode="markers", marker=dict(size=14, color=color, symbol="diamond"),
            name=allele, showlegend=False,
            hovertemplate=f"<b>{allele}</b><br>%{{x:.1f}}% [CI {100*ci[0]:.1f}-{100*ci[1]:.1f}]<extra></extra>",
        ), row=2, col=1)

    fig.update_yaxes(tickmode="array", tickvals=pos,
                     ticktext=list(RISK_ALLELES.keys()), row=2, col=1)
    fig.update_xaxes(title="carrier rate (%)", row=2, col=1)
    fig.update_layout(
        title=f"<b>Fig FINAL-12 — Korean K2 (PRJEB11591) arcasHLA imputation summary dashboard · n={n}, DM1={n_dm1}, DM2={n_dm2}, Normal={n_norm}</b>",
        height=820, paper_bgcolor="white", plot_bgcolor="white",
        font=dict(family="-apple-system,Segoe UI,Noto Sans KR,sans-serif", size=11),
    )
    out = FIG_DIR / "FINAL_12_dashboard.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def main():
    print(f"=== Loading from {RESULTS_DIR} ===")
    geno = load_genotypes()
    print(f"  {len(geno)} valid genotypes")
    df = merge_dm(geno)
    print(f"  {len(df)} merged with DM1/DM2")
    if len(df) == 0:
        print("  No samples available")
        return
    print("\n=== Generating 12 polished figures ===")
    fig01_headline(df)
    fig02_carrier_forest(df)
    fig03_dm_split_table(df)
    fig04_histology_dpb1(df)
    fig05_allele_frequency_panel(df)
    fig06_genotype_heatmap(df)
    fig07_dpb1_violin(df)
    fig08_korean_vs_lit(df)
    fig09_zygosity(df)
    fig10_per_sample_table(df)
    fig11_dm1_vs_dm2_per_gene(df)
    fig12_summary_dashboard(df)
    df.to_csv(RESULTS_DIR / "K2_arcasHLA_FINAL.tsv", sep="\t", index=False)
    print(f"\n✓ All 12 figures saved to {FIG_DIR}")
    print(f"✓ Final TSV: K2_arcasHLA_FINAL.tsv (n={len(df)})")


if __name__ == "__main__":
    main()
