#!/usr/bin/env python3
from __future__ import annotations

import json

import gseapy as gp
import numpy as np
import pandas as pd
import plotly.express as px
import requests
from scipy import stats

from v17p3_common import FIG, LOG, MSIG_H, MSIG_R, SEED, TAB, THYROID_DIFF, bh_fdr, dm_labels, load_dark, log_line, parse_enrichr_text, parse_reactome_txt, score_expr, write_json


def fetch_kegg() -> dict[str, list[str]]:
    r = requests.get("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human", timeout=60)
    r.raise_for_status()
    return parse_enrichr_text(r.text)


def run_prerank(rank_series: pd.Series, gene_sets: dict[str, list[str]], source: str) -> pd.DataFrame:
    pre = gp.prerank(rnk=rank_series, gene_sets=gene_sets, min_size=5, max_size=500, seed=SEED, verbose=False, outdir=None)
    df = pre.res2d.copy().reset_index().rename(columns={"index": "Term"})
    df["source"] = source
    return df


def main() -> None:
    log_line(LOG, "v17p3 F2 start")
    expr, meta = load_dark()
    expr.columns = expr.columns.astype(str).str.upper()
    y = dm_labels(meta)
    deg_rows = []
    for gene in expr.columns:
        a = expr.loc[y == 0, gene].to_numpy(float)
        b = expr.loc[y == 1, gene].to_numpy(float)
        p = float(stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue)
        lfc = float(np.nanmean(b) - np.nanmean(a))
        deg_rows.append({"gene": gene, "log2fc_dm2_vs_dm1": lfc, "pvalue": p})
    deg = pd.DataFrame(deg_rows)
    deg["rank"] = deg["log2fc_dm2_vs_dm1"] * (-np.log10(np.clip(deg["pvalue"], 1e-300, 1.0)))
    rank = deg.set_index("gene")["rank"].dropna()

    hall = json.load(open(MSIG_H))
    react = parse_reactome_txt(MSIG_R)
    kegg = fetch_kegg()
    hall_df = run_prerank(rank, hall, "Hallmark").rename(columns={"Term": "pathway"})
    react_df = run_prerank(rank, react, "Reactome").rename(columns={"Term": "pathway"})
    kegg_df = run_prerank(rank, kegg, "KEGG").rename(columns={"Term": "pathway"})
    hall_df.to_csv(TAB / "F2_gsea_hallmark_proper.tsv", sep="\t", index=False)
    react_df.to_csv(TAB / "F2_gsea_reactome_proper.tsv", sep="\t", index=False)
    kegg_df.to_csv(TAB / "F2_gsea_kegg_proper.tsv", sep="\t", index=False)

    # Keep per-sample scoring intentionally compact. Hallmark is the main target;
    # Reactome/KEGG are used for pathway ranking tables rather than exhaustive per-sample matrices.
    ssgsea_df = pd.DataFrame(index=expr.index)
    for term, genes in hall.items():
        keep = [g.upper() for g in genes if g.upper() in expr.columns]
        if len(keep) >= 5:
            ssgsea_df[f"Hallmark:{term}"] = expr[keep].mean(axis=1)
    ssgsea_df["sample_id"] = expr.index
    ssgsea_df["cluster"] = meta["v17_dark_cluster"].values
    ssgsea_df.to_csv(TAB / "F2_ssgsea_per_sample.tsv", sep="\t", index=False)

    proxy = pd.read_csv("project/results/v17p2/tables/gsea_dm1_vs_dm2.tsv", sep="\t")
    comp = hall_df[["pathway", "NES", "FDR q-val"]].merge(proxy[["Term", "NES", "FDR q-val"]], left_on="pathway", right_on="Term", how="inner", suffixes=("_proper", "_proxy"))
    comp["sign_agree"] = np.sign(comp["NES_proper"]) == np.sign(comp["NES_proxy"])
    comp.to_csv(TAB / "F2_proxy_vs_proper_comparison.tsv", sep="\t", index=False)

    top = pd.concat([
        hall_df.sort_values("NES", ascending=False).head(12),
        hall_df.sort_values("NES", ascending=True).head(12),
    ], ignore_index=True)
    px.bar(top, x="NES", y="pathway", color="FDR q-val", orientation="h", title="Proper GSEA top pathways").write_html(
        FIG / "F2_gsea_top_pathways_dual_dir.html", include_plotlyjs="cdn"
    )
    hall_cols = [c for c in ssgsea_df.columns if c.startswith("Hallmark:")][:25]
    hm = ssgsea_df.sort_values("cluster")
    px.imshow(hm[hall_cols].T, aspect="auto", title="Hallmark ssGSEA heatmap").write_html(
        FIG / "F2_hallmark_heatmap_per_sample.html", include_plotlyjs="cdn"
    )
    hv = hall_df.copy()
    hv["mlog10fdr"] = -np.log10(np.clip(hv["FDR q-val"].astype(float), 1e-300, 1.0))
    px.scatter(hv, x="NES", y="mlog10fdr", hover_name="pathway", title="Hallmark pathway volcano").write_html(
        FIG / "F2_pathway_volcano.html", include_plotlyjs="cdn"
    )
    if not comp.empty:
        px.scatter(comp, x="NES_proxy", y="NES_proper", color="sign_agree", hover_name="pathway", title="Proxy vs proper GSEA").write_html(
            FIG / "F2_proxy_vs_proper_scatter.html", include_plotlyjs="cdn"
        )

    summary = {
        "hallmark_fdr_lt_0_01": int((hall_df["FDR q-val"].astype(float) < 0.01).sum()),
        "reactome_fdr_lt_0_05": int((react_df["FDR q-val"].astype(float) < 0.05).sum()),
        "proxy_proper_sign_agreement": float(comp["sign_agree"].mean()) if not comp.empty else None,
    }
    write_json(TAB / "F2_summary.json", summary)
    log_line(LOG, f"v17p3 F2 done {summary}")


if __name__ == "__main__":
    main()
