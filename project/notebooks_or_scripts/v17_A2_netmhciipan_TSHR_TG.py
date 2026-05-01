#!/usr/bin/env python3
"""A2 — NetMHCIIpan TSHR + Tg autoreactive epitope binding × Korean HLA Class II alleles.

Mechanism layer: WHY DPB1*05:01 is Korean Graves' risk allele.
Test: does DPB1*05:01 bind TSHR / Tg self-peptides STRONGER than other DPB1 alleles?
"""
from __future__ import annotations
import time, requests
from pathlib import Path
from io import StringIO
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

FIG_DIR = Path("/opt/thyroid-dash/project/reports/html/figs_interactive/v17/hla_mega")
RESULTS = Path("/opt/thyroid-dash/project/results/v17_korean/A2_netmhciipan")
RESULTS.mkdir(parents=True, exist_ok=True)

IEDB_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"

# Test HLA Class II alleles (Korean-relevant)
ALLELES = [
    ("HLA-DPA1*02:02/DPB1*05:01", "DPB1*05:01 ★ Korean Graves' risk", "risk"),
    ("HLA-DPA1*02:02/DPB1*02:01", "DPB1*02:01 (2nd common Korean)", "neutral"),
    ("HLA-DPA1*02:02/DPB1*04:02", "DPB1*04:02 (Korean common)", "neutral"),
    ("HLA-DRB1*15:01",            "DRB1*15:01 ★ protective", "protective"),
    ("HLA-DRB1*04:05",            "DRB1*04:05 (Korean *04 form)", "neutral"),
    ("HLA-DRB1*09:01",            "DRB1*09:01 (Asian common)", "neutral"),
    ("HLA-DRB1*04:01",            "DRB1*04:01 (EUR Hashimoto)", "EUR-risk"),
    ("HLA-DRB1*03:01",            "DRB1*03:01 (EUR Graves')", "EUR-risk"),
]

# TSHR A-subunit immunodominant epitopes (literature curated, Carayanniotis 2007 + Latrofa 2008)
# A-subunit residues 22-260 of TSHR (residue 1 = signal peptide cleavage)
TSHR_KEY_REGIONS = [
    ("TSHR A-subunit 22-49",    "MGCSSPPCECHQEEDFRVTCKDIQRIPSLPP"),
    ("TSHR A-subunit 56-83",    "TLKLIETHLRTIPSHAFSNLPNISRIYV"),
    ("TSHR A-subunit 132-159",  "DRDPFNNTTPVTGASRAVDLIKDLLLT"),  # immunodominant LRR3
    ("TSHR A-subunit 200-227",  "SQLQIKGLNFVNLDISTLPYNYKEKDF"),
    ("TSHR A-subunit 248-275",  "VEVLQSLLGAYLLGASGYSIKAFNQHK"),  # cleavage region
]

# Thyroglobulin (Tg) immunodominant epitopes (Latrofa 2017, Tomer 2003 review)
TG_KEY_REGIONS = [
    ("Tg 88-115",       "TGAQTLAILVPDLPGSSIPEDFKVTGAA"),
    ("Tg 247-275",      "FRRSGQYSYKRFFNGVLIPCQYGRVNF"),  # canonical Hash epitope
    ("Tg 1571-1599",    "DKPVVKLITNSIQKKGSLSRGTMEQQG"),  # Hash epitope #2
    ("Tg 2098-2125",    "VQVKQRLPCKSGTKLLPPKWSPEEVF"),
    ("Tg 2540-2560",    "QHVRFQCDEQCAYGLCISTV"),
]


def query_iedb(seq, allele, length=15):
    """Single IEDB MHCII API call."""
    data = {
        "method": "netmhciipan_ba",
        "sequence_text": seq,
        "allele": allele,
        "length": str(length),
    }
    try:
        r = requests.post(IEDB_URL, data=data, timeout=60)
        if r.status_code == 200:
            return r.text
    except Exception as e:
        print(f"  ! API error: {e}")
    return None


def parse_iedb(text):
    """Parse IEDB response into DataFrame."""
    if not text or text.startswith("Error"): return pd.DataFrame()
    lines = text.strip().split("\n")
    if len(lines) < 2: return pd.DataFrame()
    df = pd.read_csv(StringIO(text), sep="\t")
    return df


def run_predictions():
    """Run NetMHCIIpan for all allele × region combinations."""
    rows = []
    all_regions = [(name, seq, "TSHR") for name, seq in TSHR_KEY_REGIONS] + \
                  [(name, seq, "Tg") for name, seq in TG_KEY_REGIONS]
    total = len(ALLELES) * len(all_regions)
    done = 0
    for allele, label, group in ALLELES:
        for region_name, seq, protein in all_regions:
            done += 1
            print(f"  [{done}/{total}] {allele} × {region_name}", flush=True)
            text = query_iedb(seq, allele, length=15)
            df = parse_iedb(text)
            if len(df) > 0:
                # Best binder per allele × region
                best = df.iloc[df["ic50"].astype(float).idxmin()]
                rows.append(dict(
                    allele=allele, label=label, group=group,
                    region=region_name, protein=protein,
                    seq=seq, best_peptide=best["peptide"],
                    ic50=float(best["ic50"]), rank=float(best["rank"]),
                ))
            time.sleep(0.5)  # gentle on API
    return pd.DataFrame(rows)


def figure_heatmap(res):
    """Heatmap: TSHR/Tg region × HLA allele, IC50 (nM)."""
    pivot = res.pivot(index="region", columns="label", values="ic50")
    # Order regions: TSHR first then Tg
    order = [r[0] for r in TSHR_KEY_REGIONS] + [r[0] for r in TG_KEY_REGIONS]
    pivot = pivot.reindex(order)
    # Order alleles: risk → neutral → protective → EUR
    col_order = [l for _,l,_ in ALLELES]
    pivot = pivot[col_order]

    fig = go.Figure(data=go.Heatmap(
        z=pivot.values, x=pivot.columns, y=pivot.index,
        colorscale="RdBu_r", zmin=10, zmax=2000, reversescale=True,
        text=[[f"{v:.0f}" if not np.isnan(v) else "—" for v in row] for row in pivot.values],
        texttemplate="%{text}", textfont=dict(size=11),
        colorbar=dict(title="IC50 (nM)<br>← strong<br>weak →", thickness=12),
        hovertemplate="<b>%{y}</b> × %{x}<br>IC50 = %{z:.0f} nM<extra></extra>",
    ))
    # Add row separator between TSHR and Tg
    n_tshr = len(TSHR_KEY_REGIONS)
    fig.add_hline(y=n_tshr-0.5, line_color="black", line_width=2)
    fig.update_layout(
        title="<b>★ Fig A2 — NetMHCIIpan TSHR / Tg autoreactive epitope × Korean HLA Class II alleles binding (IC50 nM, lower = stronger binder)</b>",
        height=520, paper_bgcolor="white", plot_bgcolor="white",
        margin=dict(l=200, r=20, t=80, b=120),
        xaxis=dict(side="top", tickangle=-30),
    )
    out = FIG_DIR / "A2_netmhciipan_heatmap.html"
    fig.write_html(out, include_plotlyjs="cdn", config={"displaylogo": False, "toImageButtonOptions":{"format":"png","scale":2}})
    print(f"  ✓ {out.name}")


def figure_strong_binder_count(res):
    """Per-allele count of strong binders (IC50 < 500 nM)."""
    res["strong"] = (res["ic50"] < 500).astype(int)
    by_allele = res.groupby(["allele","label","group"]).agg(
        strong_n=("strong","sum"), total=("ic50","size"),
        median_ic50=("ic50","median")
    ).reset_index()
    by_allele["strong_pct"] = 100 * by_allele["strong_n"] / by_allele["total"]

    color_map = {"risk":"#dc2626","protective":"#10b981","EUR-risk":"#9b59b6","neutral":"#94a3b8"}
    colors = [color_map.get(g, "#888") for g in by_allele["group"]]

    fig = make_subplots(rows=1, cols=2, column_widths=[0.6, 0.4],
                        subplot_titles=(
                            "<b>(a) Strong binder count (IC50 &lt; 500 nM) per allele</b>",
                            "<b>(b) Median IC50 (lower = stronger overall)</b>"))
    fig.add_trace(go.Bar(y=by_allele["label"], x=by_allele["strong_n"], orientation="h",
                         marker_color=colors,
                         text=[f"{n}/{t} ({p:.0f}%)" for n,t,p in zip(by_allele["strong_n"], by_allele["total"], by_allele["strong_pct"])],
                         textposition="outside", showlegend=False,
                         hovertemplate="%{y}<br>strong n=%{x}<extra></extra>"), row=1, col=1)
    fig.add_trace(go.Bar(y=by_allele["label"], x=by_allele["median_ic50"], orientation="h",
                         marker_color=colors,
                         text=[f"{v:.0f}" for v in by_allele["median_ic50"]], textposition="outside", showlegend=False,
                         hovertemplate="%{y}<br>median IC50=%{x:.0f} nM<extra></extra>"), row=1, col=2)
    fig.update_yaxes(autorange="reversed", row=1, col=1)
    fig.update_yaxes(autorange="reversed", row=1, col=2)
    fig.update_layout(
        title="<b>Fig A2-summary — Korean Graves' risk DPB1*05:01 vs protective DRB1*15:01 vs EUR allele · TSHR + Tg self-peptide binding</b>",
        height=520, paper_bgcolor="white", plot_bgcolor="white",
        margin=dict(l=240, r=80, t=80, b=20),
    )
    out = FIG_DIR / "A2_strong_binder_summary.html"
    fig.write_html(out, include_plotlyjs="cdn", config={"displaylogo": False, "toImageButtonOptions":{"format":"png","scale":2}})
    print(f"  ✓ {out.name}")


def main():
    print(f"=== A2 NetMHCIIpan — {len(ALLELES)} alleles × {len(TSHR_KEY_REGIONS)+len(TG_KEY_REGIONS)} regions ===")
    res = run_predictions()
    print(f"\n=== {len(res)} predictions complete ===")
    res.to_csv(RESULTS / "A2_netmhciipan_results.tsv", sep="\t", index=False)
    print(f"  ✓ {RESULTS / 'A2_netmhciipan_results.tsv'}")
    print("\n=== Figures ===")
    figure_heatmap(res)
    figure_strong_binder_count(res)
    print(f"\n=== summary by allele ===")
    print(res.groupby(["label","group"]).agg(
        median_ic50=("ic50","median"),
        strong_n=("ic50", lambda x: (x<500).sum()),
        total=("ic50","size"),
    ).to_string())


if __name__ == "__main__":
    main()
