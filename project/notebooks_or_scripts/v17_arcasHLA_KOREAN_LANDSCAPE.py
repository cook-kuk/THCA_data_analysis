#!/usr/bin/env python3
"""DM-independent Korean HLA landscape figures (haplotype LD, homozygosity, supertype, HED).

Builds 6 additional figures that don't depend on DM1/DM2 axis — pure population genetics.
"""
from __future__ import annotations
from pathlib import Path
import json, re
from collections import Counter
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

RESULTS = Path("/opt/thyroid-dash/project/results/v17_korean/arcasHLA")
FIG_DIR = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17/hla_mega")
COMMON = dict(include_plotlyjs="cdn",
              config={"displaylogo": False, "toImageButtonOptions":{"format":"png","scale":2}})

# Class I supertypes (Sidney 2008, Lund-Johansen 2018)
SUPERTYPE = {
    "A*01": "A1", "A*36": "A1", "A*80": "A1",
    "A*02": "A2", "A*68": "A2", "A*69": "A2",
    "A*03": "A3", "A*11": "A3", "A*30": "A3", "A*31": "A3", "A*33": "A3", "A*68": "A3",
    "A*23": "A24", "A*24": "A24",
    "A*26": "A26",
    "B*07": "B7", "B*35": "B7", "B*51": "B7", "B*53": "B7", "B*54": "B7", "B*55": "B7", "B*56": "B7",
    "B*08": "B8",
    "B*27": "B27",
    "B*38": "B39", "B*39": "B39",
    "B*15": "B44", "B*40": "B44", "B*41": "B44", "B*44": "B44", "B*45": "B44", "B*47": "B44", "B*49": "B44", "B*50": "B44",
    "B*57": "B58", "B*58": "B58",
    "B*46": "B46-Asian",
}


def load():
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
    df = df.merge(pred, on="run").merge(manifest, left_on="run", right_on="run_accession", how="left")
    df["is_normal"] = df["sample_alias"].str.endswith("-N", na=False)
    def cat(a):
        if not isinstance(a, str): return "unknown"
        if a.endswith("-N"): return "Matched normal"
        m = re.match(r"SNU-GMI-([A-Z]+)\d+", a)
        return {"PTC":"PTC","FV":"FVPTC","FTC":"FTC","FA":"FA"}.get(m.group(1), m.group(1)) if m else "unknown"
    df["histology"] = df["sample_alias"].apply(cat)
    return df


# ============= 6 ADDITIONAL FIGURES =============

def fig_pop_top_alleles(df):
    """Top allele frequencies — Korean population characterization (DM-independent)."""
    fig = make_subplots(rows=2, cols=3,
                        subplot_titles=[f"<b>HLA-{g}</b><br><sub>top 8 alleles</sub>" for g in ["A","B","C","DRB1","DQB1","DPB1"]],
                        vertical_spacing=0.18, horizontal_spacing=0.10)
    pos = [(1,1),(1,2),(1,3),(2,1),(2,2),(2,3)]
    KOR_TOP_REF = {  # Kim 2014 Korean reference top frequencies (approximate)
        "A": {"A*24:02":21.0, "A*33:03":12.4, "A*02:01":11.0, "A*02:06":7.5},
        "B": {"B*51:01":9.4, "B*44:03":6.8, "B*40:01":4.8},
        "C": {"C*01:02":11.6, "C*03:03":6.5, "C*14:02":7.1},
        "DRB1": {"DRB1*15:01":12.3, "DRB1*04:05":8.5, "DRB1*09:01":7.0},
        "DQB1": {"DQB1*06:02":11.5, "DQB1*03:03":6.6},
        "DPB1": {"DPB1*05:01":35.5, "DPB1*02:01":23.8, "DPB1*04:02":8.7},
    }
    for (gene, (r,c)) in zip(["A","B","C","DRB1","DQB1","DPB1"], pos):
        a1, a2 = f"{gene}_a1_4d", f"{gene}_a2_4d"
        if a1 not in df: continue
        all_a = pd.concat([df[a1].dropna(), df[a2].dropna()])
        top = all_a.value_counts().head(8)
        af_ours = (top / (2 * len(df)) * 100).round(2)
        af_ref = [KOR_TOP_REF.get(gene, {}).get(a, 0) for a in top.index]
        fig.add_trace(go.Bar(x=top.index, y=af_ours.values, name="K2 (ours)", marker_color="#dc2626",
                             showlegend=(r==1 and c==1),
                             hovertemplate="<b>%{x}</b><br>K2 AF: %{y:.1f}%<extra></extra>"), row=r, col=c)
        fig.add_trace(go.Bar(x=top.index, y=af_ref, name="Kim 2014 Korean ref", marker_color="#f59e0b",
                             showlegend=(r==1 and c==1),
                             hovertemplate="<b>%{x}</b><br>Kim 2014 AF: %{y:.1f}%<extra></extra>"), row=r, col=c)
    fig.update_layout(barmode="group",
        title=f"<b>Fig KOR-1 — Korean HLA top allele frequency replication (K2 ours vs Kim 2014 reference n=413, DM-independent)</b>",
        height=720, paper_bgcolor="white", plot_bgcolor="white",
    )
    out = FIG_DIR / "KOR_01_top_alleles.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig_haplotype(df):
    """Top DRB1-DQA1-DQB1 haplotypes (Class II LD block)."""
    haps = []
    for _, r in df.iterrows():
        for i in [1, 2]:
            drb1 = r.get(f"DRB1_a{i}_4d")
            dqa1 = r.get(f"DQA1_a{i}_4d")
            dqb1 = r.get(f"DQB1_a{i}_4d")
            if all(isinstance(x, str) for x in [drb1, dqa1, dqb1]):
                haps.append(f"{drb1} | {dqa1} | {dqb1}")
    top = Counter(haps).most_common(15)
    labels = [h for h,_ in top]
    counts = [c for _,c in top]
    af = [c / (2 * len(df)) * 100 for c in counts]
    # Mark known haplotypes
    notable = {
        "DRB1*15:01 | DQA1*01:02 | DQB1*06:02": "★ DR15-DQ6 (Caucasian Graves' protective)",
        "DRB1*03:01 | DQA1*05:01 | DQB1*02:01": "★ DR3-DQ2 (Caucasian Graves' risk)",
        "DRB1*04:05 | DQA1*03:03 | DQB1*04:01": "★ DR4-DQ4 (Asian)",
        "DRB1*09:01 | DQA1*03:02 | DQB1*03:03": "★ DR9-DQ9 (Asian)",
        "DRB1*08:03 | DQA1*01:03 | DQB1*06:01": "★ DR8-DQ6 (Asian)",
    }
    annot = [notable.get(l, "") for l in labels]
    colors = ["#dc2626" if "★" in a else "#7ccfcd" for a in annot]
    fig = go.Figure(go.Bar(
        y=labels, x=af, orientation="h", marker_color=colors,
        text=[f"{c} obs ({a:.1f}%)<br>{ann}" for c,a,ann in zip(counts,af,annot)],
        textposition="outside",
        hovertemplate="%{y}<br>%{x:.1f}%<extra></extra>",
    ))
    fig.update_layout(
        title=f"<b>Fig KOR-2 — Korean K2 top DRB1-DQA1-DQB1 haplotypes (n={len(df)}, DM-independent) ★ literature-known marked</b>",
        xaxis_title="haplotype frequency (%)", yaxis=dict(autorange="reversed"),
        height=600, paper_bgcolor="white", plot_bgcolor="white",
    )
    out = FIG_DIR / "KOR_02_haplotype.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig_homozygosity(df):
    """Homozygosity rate per gene + comparison to HWE expected."""
    rows = []
    for g in ["A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1"]:
        a1, a2 = f"{g}_a1_4d", f"{g}_a2_4d"
        if a1 not in df: continue
        homo = ((df[a1] == df[a2]) & df[a1].notna()).sum()
        n = df[a1].notna().sum()
        # Compute HWE expected: Σ p_i² where p_i = AF for each allele
        all_a = pd.concat([df[a1].dropna(), df[a2].dropna()])
        afs = (all_a.value_counts() / len(all_a)).values
        hwe_expected = (afs ** 2).sum() * 100
        observed = 100 * homo / max(n, 1)
        rows.append(dict(gene=f"HLA-{g}", observed=observed, hwe=hwe_expected, n_homo=homo, n_total=n))
    rdf = pd.DataFrame(rows)
    fig = go.Figure()
    fig.add_trace(go.Bar(x=rdf["gene"], y=rdf["observed"], name="Observed (K2 ours)", marker_color="#dc2626",
                         text=[f"{int(n)}/{int(t)}<br>{v:.1f}%" for v,n,t in zip(rdf["observed"], rdf["n_homo"], rdf["n_total"])],
                         textposition="outside",
                         hovertemplate="%{x}<br>Observed: %{y:.1f}%<extra></extra>"))
    fig.add_trace(go.Bar(x=rdf["gene"], y=rdf["hwe"], name="HWE expected (Σp²)", marker_color="#7ccfcd",
                         hovertemplate="%{x}<br>HWE: %{y:.1f}%<extra></extra>"))
    fig.update_layout(barmode="group",
        title=f"<b>Fig KOR-3 — HLA homozygosity rate per gene · observed vs HWE expected (n={len(df)})</b>",
        yaxis_title="homozygosity (%)",
        height=460, paper_bgcolor="white", plot_bgcolor="white",
    )
    out = FIG_DIR / "KOR_03_homozygosity.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig_class1_supertype(df):
    """Class I supertype frequency distribution (HLA-A + B)."""
    supertypes = []
    for _, r in df.iterrows():
        for g in ["A","B"]:
            for i in [1, 2]:
                a = r.get(f"{g}_a{i}_4d")
                if isinstance(a, str):
                    two_d = a.split(":")[0]  # "A*02"
                    supertypes.append(SUPERTYPE.get(two_d, "Other"))
    cnt = Counter(supertypes)
    items = sorted(cnt.items(), key=lambda x: -x[1])
    labels = [k for k,_ in items]
    counts = [v for _,v in items]
    af = [v / sum(counts) * 100 for v in counts]
    colors = ["#dc2626" if "Asian" in l or "B27" in l else
              ("#f59e0b" if l.startswith("B") else "#3b82f6") for l in labels]
    fig = go.Figure(go.Bar(
        x=labels, y=af, marker_color=colors,
        text=[f"{c}<br>{a:.1f}%" for c,a in zip(counts, af)],
        textposition="outside",
        hovertemplate="<b>%{x}</b><br>count=%{customdata}<br>%{y:.1f}%<extra></extra>",
        customdata=counts,
    ))
    fig.update_layout(
        title=f"<b>Fig KOR-4 — Class I (HLA-A + B) supertype frequency distribution (n={len(df)} samples × 4 alleles each)</b>",
        xaxis_title="HLA Class I supertype", yaxis_title="frequency (%)",
        height=460, paper_bgcolor="white", plot_bgcolor="white",
    )
    out = FIG_DIR / "KOR_04_supertype.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig_hed_distribution(df):
    """HED (HLA Evolutionary Divergence) proxy — unique allele count per sample × histology."""
    rows = []
    for _, r in df.iterrows():
        c1 = set()
        c2 = set()
        for g in ["A","B","C"]:
            for i in [1,2]:
                a = r.get(f"{g}_a{i}_4d")
                if isinstance(a, str): c1.add(a)
        for g in ["DRB1","DQB1","DPB1"]:
            for i in [1,2]:
                a = r.get(f"{g}_a{i}_4d")
                if isinstance(a, str): c2.add(a)
        rows.append(dict(run=r["run"], hist=r.get("histology","unknown"),
                         dm=r.get("DM_call","NA"), c1_n=len(c1), c2_n=len(c2)))
    hed = pd.DataFrame(rows)
    fig = make_subplots(rows=1, cols=2, subplot_titles=(
        "<b>(a) Class I unique alleles (max 6)</b>",
        "<b>(b) Class II unique alleles (max 6)</b>"))
    pal = {"PTC":"#3b82f6","FVPTC":"#7ccfcd","FTC":"#9b59b6","FA":"#f59e0b","Matched normal":"#10b981","unknown":"#94a3b8"}
    for col, score in enumerate(["c1_n", "c2_n"]):
        for cat, color in pal.items():
            sub = hed[hed["hist"] == cat][score]
            if len(sub) == 0: continue
            fig.add_trace(go.Box(y=sub, name=cat, marker_color=color,
                                 boxmean=True, showlegend=(col==0)), row=1, col=col+1)
    fig.update_layout(
        title=f"<b>Fig KOR-5 — HED proxy (unique allele count) by histology · ICI response prediction layer (n={len(df)})</b>",
        height=460, paper_bgcolor="white", plot_bgcolor="white",
    )
    fig.update_yaxes(title="unique allele count", range=[2,7], row=1, col=1)
    fig.update_yaxes(title="unique allele count", range=[2,7], row=1, col=2)
    out = FIG_DIR / "KOR_05_hed.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def fig_a02_06(df):
    """Asian-specific HLA-A*02:06 — DM-independent population genetic check."""
    a02_06_carrier = df.apply(lambda r: any(
        isinstance(r.get(f"A_a{i}_4d"), str) and r.get(f"A_a{i}_4d") == "A*02:06"
        for i in [1,2]
    ), axis=1)
    n_car = int(a02_06_carrier.sum())
    n = len(df)
    fig = go.Figure(go.Pie(
        labels=[f"A*02:06 carrier (n={n_car})", f"Non-carrier (n={n-n_car})"],
        values=[n_car, n-n_car],
        marker=dict(colors=["#dc2626", "#7ccfcd"], line=dict(color="white", width=2)),
        hole=0.55, sort=False, textinfo="label+percent",
    ))
    pct = round(100*n_car/n, 1)
    fig.update_layout(
        title=f"<b>Fig KOR-6 — A*02:06 (Asian-specific) carrier rate · K2 {pct}% vs Korean ref ~12% (Kim 2014)</b>",
        annotations=[dict(text=f"{pct}%<br><sub>K2 ours</sub>", showarrow=False, x=0.5, y=0.5, font=dict(size=22))],
        height=440, paper_bgcolor="white",
    )
    out = FIG_DIR / "KOR_06_a02_06.html"
    fig.write_html(out, **COMMON); print(f"  ✓ {out.name}")


def main():
    df = merge_meta(load())
    print(f"=== n = {len(df)} samples ===")
    print("Generating 6 KOREAN-population figures (DM-independent)...")
    fig_pop_top_alleles(df)
    fig_haplotype(df)
    fig_homozygosity(df)
    fig_class1_supertype(df)
    fig_hed_distribution(df)
    fig_a02_06(df)
    print(f"✓ 6 figures saved to {FIG_DIR}")


if __name__ == "__main__":
    main()
