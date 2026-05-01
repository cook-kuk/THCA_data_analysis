#!/usr/bin/env python3
"""DEEP DM-independent Korean HLA landscape analysis.

8 additional figures focused on population genetics, allele combinations,
haplotype LD, supertype × histology — all WITHOUT requiring DM1/DM2 axis.
"""
from __future__ import annotations
from pathlib import Path
import json, re
from collections import Counter
from itertools import combinations
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats

RESULTS = Path("/opt/thyroid-dash/project/results/v17_korean/arcasHLA")
FIG_DIR = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17/hla_mega")
COMMON = dict(include_plotlyjs="cdn", config={"displaylogo": False, "toImageButtonOptions":{"format":"png","scale":2}})

# Pan-Asian reference frequencies (published)
PAN_ASIAN_AF = {
    "DPB1*05:01": {"Korean(Kim14)": 36.0, "Japanese": 33.0, "Han": 32.5, "EUR": 12.0, "AFR": 5.0},
    "DPB1*02:01": {"Korean(Kim14)": 27.0, "Japanese": 23.0, "Han": 28.0, "EUR": 10.0, "AFR": 8.0},
    "B*46:01":    {"Korean(Kim14)": 8.0,  "Japanese": 4.0,  "Han": 11.0, "EUR": 0.1,  "AFR": 0.05},
    "B*51:01":    {"Korean(Kim14)": 18.0, "Japanese": 16.0, "Han": 12.5, "EUR": 8.5,  "AFR": 4.0},
    "B*54:01":    {"Korean(Kim14)": 12.0, "Japanese": 13.0, "Han": 9.0,  "EUR": 0.0,  "AFR": 0.0},
    "DRB1*15:01": {"Korean(Kim14)": 12.5, "Japanese": 14.0, "Han": 12.0, "EUR": 26.0, "AFR": 9.0},
    "DRB1*09:01": {"Korean(Kim14)": 11.0, "Japanese": 13.5, "Han": 16.0, "EUR": 0.5,  "AFR": 4.0},
    "DRB1*04:05": {"Korean(Kim14)": 8.5,  "Japanese": 7.0,  "Han": 5.0,  "EUR": 0.5,  "AFR": 0.5},
    "A*24:02":    {"Korean(Kim14)": 36.0, "Japanese": 35.0, "Han": 30.0, "EUR": 18.0, "AFR": 5.0},
    "A*33:03":    {"Korean(Kim14)": 24.0, "Japanese": 18.0, "Han": 22.0, "EUR": 1.5,  "AFR": 1.5},
    "A*02:06":    {"Korean(Kim14)": 13.0, "Japanese": 12.0, "Han": 10.0, "EUR": 0.0,  "AFR": 0.0},
}


def load_full():
    rows = []
    for jp in sorted(RESULTS.glob("*.genotype.json")):
        run = jp.stem.replace(".genotype","").rstrip("_1")
        try: d = json.loads(jp.read_text())
        except: continue
        if not d or not any(g in d for g in ["A","B","C","DRB1","DQB1","DPB1"]): continue
        row = {"run": run}
        for g in ["A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1"]:
            for i, a in enumerate(d.get(g, [])[:2], 1):
                row[f"{g}_a{i}"] = a
                if a:
                    parts = a.split(":")
                    row[f"{g}_a{i}_4d"] = ":".join(parts[:2])
        rows.append(row)
    return pd.DataFrame(rows)


def merge_meta(df):
    pred = pd.read_csv("/opt/thyroid-dash/project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t")[["run","p_DM2","DM_call"]]
    manifest = pd.read_csv("/opt/thyroid-dash/project/results/v17_korean/K1A_prjeb11591_runs.tsv", sep="\t")[["run_accession","sample_alias"]]
    df = df.merge(pred, on="run", how="left").merge(manifest, left_on="run", right_on="run_accession", how="left")
    df["is_normal"] = df["sample_alias"].str.endswith("-N", na=False)
    def cat(a):
        if not isinstance(a, str): return "unknown"
        if a.endswith("-N"): return "Matched normal"
        m = re.match(r"SNU-GMI-([A-Z]+)\d+", a)
        return {"PTC":"PTC","FV":"FVPTC","FTC":"FTC","FA":"FA"}.get(m.group(1), m.group(1)) if m else "unknown"
    df["histology"] = df["sample_alias"].apply(cat)
    return df


def is_carrier(df, gene, allele):
    target = f"{gene}*{allele}"
    a1, a2 = f"{gene}_a1_4d", f"{gene}_a2_4d"
    return df.apply(lambda r: any(isinstance(r[c], str) and r[c]==target for c in [a1,a2] if c in df), axis=1)


# ============= 8 DEEP DM-INDEPENDENT FIGURES =============

def D1_pan_asian(df):
    """Korean K2 vs Pan-Asian + EUR + AFR allele frequencies."""
    n = len(df)
    rows = []
    for allele, refs in PAN_ASIAN_AF.items():
        g = allele.split("*")[0]; sub = allele.split("*")[1]
        a1, a2 = f"{g}_a1_4d", f"{g}_a2_4d"
        all_a = pd.concat([df[a1].dropna(), df[a2].dropna()])
        ours_count = (all_a == allele).sum()
        ours_af = 100 * ours_count / (2 * n)
        rows.append(dict(allele=allele, ours=ours_af, **refs))
    rdf = pd.DataFrame(rows)

    fig = go.Figure()
    pops = ["ours", "Korean(Kim14)", "Japanese", "Han", "EUR", "AFR"]
    colors = ["#dc2626", "#f59e0b", "#9b59b6", "#7ccfcd", "#3b82f6", "#94a3b8"]
    pop_labels = ["★ K2 ours", "Korean (Kim 2014)", "Japanese", "Han Chinese", "European 1000G", "African 1000G"]
    for pop, c, lab in zip(pops, colors, pop_labels):
        fig.add_trace(go.Bar(name=lab, x=rdf["allele"], y=rdf[pop], marker_color=c,
                             hovertemplate=f"<b>%{{x}}</b><br>{lab} AF: %{{y:.1f}}%<extra></extra>"))
    fig.update_layout(barmode="group",
        title=f"<b>Fig DEEP-1 — Pan-Asian + Caucasian + African allele frequency comparison (Korean K2 ours, n={n})</b>",
        xaxis_title="HLA allele", yaxis_title="allele frequency (%)",
        height=540, paper_bgcolor="white", plot_bgcolor="white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0.5, xanchor="center"),
    )
    out = FIG_DIR / "DEEP_01_pan_asian.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def D2_risk_burden(df):
    """Risk allele burden score per sample (DM-independent)."""
    risk_alleles_pos = [("DPB1","05:01",1.0), ("B","46:01",0.7), ("DRB1","04:05",0.5)]
    protective = [("DRB1","15:01",-0.5)]
    burden = np.zeros(len(df))
    for g, a, w in risk_alleles_pos + protective:
        c1 = df.get(f"{g}_a1_4d", pd.Series([None]*len(df))).apply(lambda x: 1 if isinstance(x,str) and x==f"{g}*{a}" else 0)
        c2 = df.get(f"{g}_a2_4d", pd.Series([None]*len(df))).apply(lambda x: 1 if isinstance(x,str) and x==f"{g}*{a}" else 0)
        burden += w * (c1 + c2)
    df = df.copy(); df["burden"] = burden

    fig = make_subplots(rows=1, cols=2, column_widths=[0.6, 0.4],
                        subplot_titles=("<b>(a) Risk burden distribution (n=" + str(len(df)) + ")</b>",
                                        "<b>(b) Risk burden × histology (DM-independent)</b>"))
    fig.add_trace(go.Histogram(x=df["burden"], nbinsx=20, marker_color="#dc2626",
                               hovertemplate="burden=%{x:.1f}<br>n=%{y}<extra></extra>"), row=1, col=1)
    pal = {"PTC":"#3b82f6","FVPTC":"#7ccfcd","FTC":"#9b59b6","FA":"#f59e0b","Matched normal":"#10b981","unknown":"#94a3b8"}
    for cat in ["Matched normal", "FA", "FTC", "FVPTC", "PTC"]:
        sub = df[df["histology"]==cat]["burden"]
        if len(sub)==0: continue
        fig.add_trace(go.Box(y=sub, name=cat, marker_color=pal[cat],
                             boxmean=True, showlegend=False), row=1, col=2)
    fig.update_layout(
        title=f"<b>Fig DEEP-2 — Composite risk allele burden score (DPB1*05:01 +1 dose × Asian Graves'+ B*46:01 +0.7 + DRB1*04:05 +0.5 − DRB1*15:01 0.5)</b>",
        height=480, paper_bgcolor="white", plot_bgcolor="white",
    )
    fig.update_xaxes(title="risk burden score", row=1, col=1)
    fig.update_yaxes(title="count", row=1, col=1)
    fig.update_yaxes(title="risk burden", row=1, col=2)
    out = FIG_DIR / "DEEP_02_burden.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def D3_histology_alleles(df):
    """Per-histology top allele frequency (5 cat × 6 gene heatmap)."""
    cats = ["Matched normal", "FA", "FTC", "FVPTC", "PTC"]
    n_cats = {c: len(df[df["histology"]==c]) for c in cats}
    risk_alleles = [
        ("DPB1*05:01","Korean Graves'"),
        ("DRB1*15:01","Protective"),
        ("B*46:01","Asian Graves'"),
        ("A*24:02","Asian"),
        ("A*02:06","Asian"),
        ("DRB1*04:05","Asian *04"),
        ("DRB1*09:01","Asian DR9"),
        ("DPB1*02:01","Asian DPB1"),
    ]
    Z = np.zeros((len(risk_alleles), len(cats)))
    for i, (allele, _) in enumerate(risk_alleles):
        g = allele.split("*")[0]
        sub = allele.split("*")[1]
        for j, cat in enumerate(cats):
            sub_df = df[df["histology"]==cat]
            if len(sub_df)==0: Z[i,j]=0; continue
            car = is_carrier(sub_df, g, sub)
            Z[i,j] = 100*car.sum()/len(sub_df) if len(sub_df)>0 else 0
    text = [[f"{v:.0f}%" for v in row] for row in Z]
    fig = go.Figure(go.Heatmap(
        z=Z, x=[f"{c}<br>(n={n_cats[c]})" for c in cats],
        y=[f"<b>{a}</b><br><sub>{l}</sub>" for a,l in risk_alleles],
        colorscale="Reds", zmin=0, zmax=70,
        text=text, texttemplate="%{text}", textfont=dict(size=12),
        colorbar=dict(title="carrier %", thickness=12),
    ))
    fig.update_layout(
        title=f"<b>Fig DEEP-3 — Allele carrier rate × histology heatmap (DM-independent, 8 alleles × 5 categories)</b>",
        height=540, paper_bgcolor="white", plot_bgcolor="white",
        margin=dict(l=180, r=20, t=80, b=80),
    )
    out = FIG_DIR / "DEEP_03_histology_heatmap.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def D4_co_occurrence(df):
    """Allele co-occurrence matrix (LD detection)."""
    target_alleles = [
        ("DPB1","05:01"), ("DPB1","02:01"), ("DRB1","15:01"), ("DRB1","09:01"),
        ("DRB1","04:05"), ("B","46:01"), ("B","51:01"), ("A","24:02"),
        ("A","02:06"), ("A","33:03"), ("DQB1","03:03"), ("DQB1","06:02"),
    ]
    labels = [f"{g}*{a}" for g,a in target_alleles]
    n = len(df)
    M = np.zeros((len(labels), len(labels)))
    for i, (g1, a1) in enumerate(target_alleles):
        c_i = is_carrier(df, g1, a1)
        for j, (g2, a2) in enumerate(target_alleles):
            c_j = is_carrier(df, g2, a2)
            both = (c_i & c_j).sum()
            obs = both / n
            exp = c_i.mean() * c_j.mean()
            # log2(observed/expected)
            if i == j:
                M[i,j] = 0
            elif exp > 0 and obs > 0:
                M[i,j] = np.log2(obs / exp)
            else:
                M[i,j] = 0
    fig = go.Figure(go.Heatmap(
        z=M, x=labels, y=labels,
        colorscale="RdBu_r", zmid=0, zmin=-2, zmax=3,
        text=[[f"{v:+.2f}" for v in row] for row in M], texttemplate="%{text}",
        textfont=dict(size=10),
        colorbar=dict(title="log2(O/E)<br>← negative<br>(positive→ LD)", thickness=12),
        hovertemplate="<b>%{y}</b> × <b>%{x}</b><br>log2(O/E) = %{z:.2f}<extra></extra>",
    ))
    fig.update_layout(
        title=f"<b>Fig DEEP-4 — Allele co-occurrence matrix · log2(observed/expected) (n={n}, DM-independent LD detection)</b>",
        height=620, paper_bgcolor="white", plot_bgcolor="white",
        xaxis=dict(side="top"), margin=dict(l=120, r=20, t=120, b=20),
    )
    out = FIG_DIR / "DEEP_04_co_occurrence.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def D5_pca(df):
    """PCA on allele dosage matrix → see if histology / DM clusters in HLA space."""
    # Build dosage matrix (sample × allele)
    target = [
        ("A","24:02"), ("A","33:03"), ("A","02:06"), ("A","02:01"),
        ("B","51:01"), ("B","44:03"), ("B","46:01"), ("B","58:01"),
        ("C","01:02"), ("C","03:03"), ("C","14:02"),
        ("DRB1","15:01"), ("DRB1","09:01"), ("DRB1","04:05"), ("DRB1","13:02"),
        ("DQB1","06:02"), ("DQB1","03:03"),
        ("DPB1","05:01"), ("DPB1","02:01"), ("DPB1","04:02"),
    ]
    X = np.zeros((len(df), len(target)))
    for j, (g, a) in enumerate(target):
        a1 = df.get(f"{g}_a1_4d", pd.Series([None]*len(df)))
        a2 = df.get(f"{g}_a2_4d", pd.Series([None]*len(df)))
        X[:, j] = (a1 == f"{g}*{a}").astype(int) + (a2 == f"{g}*{a}").astype(int)
    # PCA via SVD
    Xc = X - X.mean(axis=0)
    U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    pc1, pc2 = U[:,0]*s[0], U[:,1]*s[1]

    pal = {"PTC":"#3b82f6","FVPTC":"#7ccfcd","FTC":"#9b59b6","FA":"#f59e0b","Matched normal":"#10b981","unknown":"#94a3b8"}
    fig = go.Figure()
    for cat in pal:
        mask = df["histology"]==cat
        if mask.sum()==0: continue
        fig.add_trace(go.Scatter(
            x=pc1[mask], y=pc2[mask], mode="markers",
            marker=dict(size=10, color=pal[cat], line=dict(width=1, color="white"), opacity=0.7),
            name=f"{cat} (n={mask.sum()})",
            hovertemplate=f"<b>{cat}</b><br>"+"PC1=%{x:.2f}<br>PC2=%{y:.2f}<extra></extra>",
            text=df["run"][mask],
        ))
    var_explained = (s**2 / (s**2).sum())[:2] * 100
    fig.update_layout(
        title=f"<b>Fig DEEP-5 — HLA allele dosage PCA × histology (n={len(df)}, 20 alleles, DM-independent)</b>",
        xaxis_title=f"PC1 ({var_explained[0]:.1f}% variance)",
        yaxis_title=f"PC2 ({var_explained[1]:.1f}% variance)",
        height=560, paper_bgcolor="white", plot_bgcolor="white",
    )
    out = FIG_DIR / "DEEP_05_pca.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def D6_top_haplotypes(df):
    """Top 20 DRB1-DQB1 + DPB1-DPA1 haplotypes."""
    drb_dqb = []
    dpb_dpa = []
    for _, r in df.iterrows():
        for i in [1,2]:
            d1 = r.get(f"DRB1_a{i}_4d"); q1 = r.get(f"DQB1_a{i}_4d")
            if isinstance(d1,str) and isinstance(q1,str):
                drb_dqb.append(f"{d1} | {q1}")
            d2 = r.get(f"DPB1_a{i}_4d"); p1 = r.get(f"DPA1_a{i}_4d")
            if isinstance(d2,str) and isinstance(p1,str):
                dpb_dpa.append(f"{d2} | {p1}")
    top_drb = Counter(drb_dqb).most_common(15)
    top_dpb = Counter(dpb_dpa).most_common(15)
    fig = make_subplots(rows=1, cols=2, column_widths=[0.5, 0.5],
                        subplot_titles=("<b>(a) DRB1-DQB1 haplotype LD (top 15)</b>",
                                        "<b>(b) DPB1-DPA1 haplotype LD (top 15)</b>"))
    fig.add_trace(go.Bar(y=[h for h,_ in top_drb], x=[c for _,c in top_drb], orientation="h",
                         marker_color="#dc2626", showlegend=False,
                         text=[f"{c} ({100*c/(2*len(df)):.1f}%)" for _,c in top_drb], textposition="outside",
                         hovertemplate="%{y}<br>n=%{x}<extra></extra>"), row=1, col=1)
    fig.add_trace(go.Bar(y=[h for h,_ in top_dpb], x=[c for _,c in top_dpb], orientation="h",
                         marker_color="#9b59b6", showlegend=False,
                         text=[f"{c} ({100*c/(2*len(df)):.1f}%)" for _,c in top_dpb], textposition="outside",
                         hovertemplate="%{y}<br>n=%{x}<extra></extra>"), row=1, col=2)
    fig.update_yaxes(autorange="reversed", row=1, col=1)
    fig.update_yaxes(autorange="reversed", row=1, col=2)
    fig.update_layout(
        title=f"<b>Fig DEEP-6 — Class II haplotype LD blocks (DRB1-DQB1 + DPB1-DPA1, n={len(df)} samples × 2 chromosomes)</b>",
        height=600, paper_bgcolor="white", plot_bgcolor="white",
    )
    out = FIG_DIR / "DEEP_06_class2_haplotypes.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def D7_carriers_per_sample(df):
    """How many risk alleles per sample (carrier load distribution)."""
    risk = [("DPB1","05:01"),("B","46:01"),("DRB1","04:05"),("DRB1","09:01"),
            ("A","02:06"),("DRB1","15:01"),("DPB1","02:01"),("B","51:01"),("B","54:01")]
    counts = np.zeros(len(df))
    for g, a in risk:
        c = is_carrier(df, g, a)
        counts += c.astype(int)
    df = df.copy(); df["risk_count"] = counts.astype(int)
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("<b>(a) Distribution of risk allele count per sample</b>",
                                        "<b>(b) Risk allele count × histology</b>"))
    cnt = df["risk_count"].value_counts().sort_index()
    fig.add_trace(go.Bar(x=cnt.index, y=cnt.values, marker_color="#dc2626",
                         text=cnt.values, textposition="outside",
                         hovertemplate="count=%{x}<br>n=%{y}<extra></extra>"), row=1, col=1)
    pal = {"PTC":"#3b82f6","FVPTC":"#7ccfcd","FTC":"#9b59b6","FA":"#f59e0b","Matched normal":"#10b981"}
    for cat in ["Matched normal","FA","FTC","FVPTC","PTC"]:
        sub = df[df["histology"]==cat]["risk_count"]
        if len(sub)==0: continue
        fig.add_trace(go.Box(y=sub, name=cat, marker_color=pal[cat], boxmean=True, showlegend=False), row=1, col=2)
    fig.update_layout(
        title=f"<b>Fig DEEP-7 — Risk allele carrier count per sample (9 Asian-relevant alleles, DM-independent, n={len(df)})</b>",
        height=480, paper_bgcolor="white", plot_bgcolor="white",
    )
    fig.update_xaxes(title="number of risk alleles carried (0-9)", row=1, col=1)
    fig.update_yaxes(title="sample count", row=1, col=1)
    fig.update_yaxes(title="risk allele count", row=1, col=2)
    out = FIG_DIR / "DEEP_07_carrier_load.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def D8_dpb1_dr_combinations(df):
    """DPB1*05:01 × DRB1 combinations (which DRB1 pairs with risk DPB1)."""
    car_dpb = is_carrier(df, "DPB1", "05:01")
    dpb_carriers = df[car_dpb]
    drb_in_carriers = pd.concat([dpb_carriers["DRB1_a1_4d"].dropna(),
                                  dpb_carriers["DRB1_a2_4d"].dropna()])
    drb_in_non = pd.concat([df[~car_dpb]["DRB1_a1_4d"].dropna(),
                            df[~car_dpb]["DRB1_a2_4d"].dropna()])
    top = drb_in_carriers.value_counts().head(10).index
    rows = []
    for d in top:
        n_in_car = (drb_in_carriers == d).sum()
        n_in_non = (drb_in_non == d).sum()
        af_car = 100 * n_in_car / max(len(drb_in_carriers),1)
        af_non = 100 * n_in_non / max(len(drb_in_non),1)
        rows.append(dict(drb1=d, carrier_af=af_car, noncarrier_af=af_non, delta=af_car-af_non))
    rdf = pd.DataFrame(rows)
    fig = go.Figure()
    fig.add_trace(go.Bar(x=rdf["drb1"], y=rdf["carrier_af"], name="DPB1*05:01 carriers", marker_color="#dc2626",
                         hovertemplate="%{x}<br>carrier AF: %{y:.1f}%<extra></extra>"))
    fig.add_trace(go.Bar(x=rdf["drb1"], y=rdf["noncarrier_af"], name="DPB1*05:01 non-carriers", marker_color="#7ccfcd",
                         hovertemplate="%{x}<br>non-carrier AF: %{y:.1f}%<extra></extra>"))
    fig.update_layout(barmode="group",
        title=f"<b>Fig DEEP-8 — Top DRB1 alleles in DPB1*05:01 carriers vs non-carriers (LD test, n_carrier={car_dpb.sum()}, n_non={(~car_dpb).sum()})</b>",
        height=520, paper_bgcolor="white", plot_bgcolor="white",
        yaxis_title="DRB1 allele frequency (%)",
    )
    out = FIG_DIR / "DEEP_08_dpb1_drb1.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def main():
    df = merge_meta(load_full())
    print(f"=== n = {len(df)} samples ===")
    print("Generating 8 DEEP DM-independent figures...")
    D1_pan_asian(df)
    D2_risk_burden(df)
    D3_histology_alleles(df)
    D4_co_occurrence(df)
    D5_pca(df)
    D6_top_haplotypes(df)
    D7_carriers_per_sample(df)
    D8_dpb1_dr_combinations(df)
    print(f"✓ 8 DEEP figures saved to {FIG_DIR}")


if __name__ == "__main__":
    main()
