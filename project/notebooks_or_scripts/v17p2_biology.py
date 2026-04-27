#!/usr/bin/env python3
from __future__ import annotations

import warnings

import gseapy as gp
import numpy as np
import pandas as pd
import plotly.express as px
from scipy import stats

from v17p2_common import FIG, LOG, RAI_GENES, SEED, TAB, TDS16, bh_fdr, dm_binary_labels, load_dark_transcriptome, log_line, signature_score, write_json


MAPK_GENES = ["DUSP4", "DUSP5", "DUSP6", "ETV1", "ETV4", "ETV5", "FOS", "FOSB", "JUN", "JUNB", "SPRY2", "CCND1"]
THYROID_DIFF = ["PAX8", "NKX2-1", "FOXE1", "SLC5A5", "TPO", "TG", "TSHR", "DIO1", "DIO2", "DUOX2", "SLC26A4"]
IMMUNE_SETS = {
    "Treg_proxy": ["FOXP3", "IL2RA", "TIGIT", "CTLA4", "ICOS"],
    "M2_proxy": ["CD163", "MRC1", "CSF1R", "MSR1", "IL10"],
    "ExhaustedT_proxy": ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX"],
}
PATHWAY_PROXIES = {
    "Hallmark_EMT": ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "FN1", "COL1A1", "COL3A1", "LOX"],
    "Hallmark_Hypoxia": ["HIF1A", "VEGFA", "SLC2A1", "LDHA", "ENO1", "HK2", "CA9"],
    "Hallmark_Oxidative_Phosphorylation": ["NDUFA1", "NDUFB8", "COX4I1", "ATP5F1A", "SDHB", "UQCRC1"],
    "Hallmark_TGF_beta": ["TGFB1", "TGFBR1", "SMAD2", "SMAD3", "SMAD4", "SERPINE1", "COL1A1"],
    "Hallmark_Inflammatory_Response": ["IL6", "CXCL8", "TNF", "NFKB1", "RELA", "PTGS2", "CXCL2"],
    "Hallmark_E2F_Targets": ["E2F1", "MCM2", "MCM5", "PCNA", "MKI67", "CDC6", "TK1"],
    "Hallmark_G2M_Checkpoint": ["CCNB1", "CDK1", "AURKA", "TOP2A", "BUB1", "MAD2L1", "PLK1"],
    "Hallmark_Apoptosis": ["BAX", "BBC3", "CASP3", "CASP8", "FAS", "PMAIP1"],
    "Hallmark_MTORC1": ["RPTOR", "MTOR", "RHEB", "EIF4EBP1", "RPS6KB1", "MLST8"],
    "Hallmark_Interferon_Gamma": ["STAT1", "IRF1", "CXCL10", "HLA-A", "B2M", "GBP1", "IFI44L"],
}


def main() -> None:
    log_line(LOG, "v17p2 layer3 start")
    expr, meta = load_dark_transcriptome()
    y = dm_binary_labels(meta)
    rows = []
    for gene in expr.columns:
        a = expr.loc[y == 0, gene].to_numpy(float)
        b = expr.loc[y == 1, gene].to_numpy(float)
        lfc = float(np.nanmean(b) - np.nanmean(a))
        p = float(stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue)
        rows.append({"gene": gene, "log2fc_dm2_vs_dm1": lfc, "pvalue": p})
    deg = pd.DataFrame(rows)
    deg["fdr"] = bh_fdr(deg["pvalue"].to_numpy(float))
    deg["rank_score"] = deg["log2fc_dm2_vs_dm1"] * (-np.log10(np.clip(deg["pvalue"], 1e-300, 1.0)))
    deg = deg.sort_values("rank_score", ascending=False)
    deg.to_csv(TAB / "dm1_vs_dm2_deg_full.tsv", sep="\t", index=False)

    expr.columns = [str(c).upper() for c in expr.columns]
    deg["gene"] = deg["gene"].astype(str).str.upper()
    gsea_parts = []
    rank_series = deg.set_index("gene")["rank_score"]
    for gs in ["MSigDB_Hallmark_2020", "KEGG_2021_Human", "Reactome_2022"]:
        try:
            gene_sets = gp.get_library(name=gs, organism="Human")
            pre = gp.prerank(rnk=rank_series, gene_sets=gene_sets, min_size=10, max_size=500, seed=SEED, verbose=False, outdir=None)
            res = pre.res2d.copy()
            res["gene_set_source"] = gs
            gsea_parts.append(res)
        except Exception as e:
            log_line(LOG, f"v17p2 layer3 gsea skipped {gs}: {e}")
    if not gsea_parts:
        proxy_rows = []
        for term, genes in PATHWAY_PROXIES.items():
            keep = [g for g in genes if g in expr.columns]
            if len(keep) < 3:
                continue
            s = signature_score(expr, keep)
            a = s[y == 0]
            b = s[y == 1]
            diff = float(np.nanmean(b) - np.nanmean(a))
            p = float(stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue)
            pooled = float(np.nanstd(pd.concat([a, b]), ddof=1)) if len(pd.concat([a, b])) > 1 else np.nan
            nes = diff / pooled if pooled and np.isfinite(pooled) and pooled > 0 else diff
            proxy_rows.append({"Term": term, "NES": nes, "FDR q-val": p, "gene_set_source": "proxy"})
        proxy = pd.DataFrame(proxy_rows)
        if not proxy.empty:
            proxy["FDR q-val"] = bh_fdr(proxy["FDR q-val"].to_numpy(float))
            gsea_parts.append(proxy)
    gsea = pd.concat(gsea_parts, ignore_index=True) if gsea_parts else pd.DataFrame(columns=["Term", "NES", "FDR q-val", "gene_set_source"])
    gsea.to_csv(TAB / "gsea_dm1_vs_dm2.tsv", sep="\t", index=False)

    hallmark_scores = None
    try:
        hall = gp.get_library(name="MSigDB_Hallmark_2020", organism="Human")
        hall_scores = {k: signature_score(expr, v) for k, v in hall.items()}
        hallmark_scores = pd.DataFrame(hall_scores, index=expr.index)
        hallmark_scores["sample_id"] = expr.index
        hallmark_scores["cluster"] = meta["v17_dark_cluster"].values
        hallmark_scores.to_csv(TAB / "hallmark_score_per_sample.tsv", sep="\t", index=False)
    except Exception as e:
        log_line(LOG, f"v17p2 layer3 hallmark scoring fallback: {e}")
        hallmark_scores = pd.DataFrame({"sample_id": expr.index, "cluster": meta["v17_dark_cluster"].values})
        for name, genes in PATHWAY_PROXIES.items():
            hallmark_scores[name] = signature_score(expr, genes).values
        hallmark_scores.to_csv(TAB / "hallmark_score_per_sample.tsv", sep="\t", index=False)

    mapk = pd.DataFrame({
        "sample_id": expr.index,
        "cluster": meta["v17_dark_cluster"].values,
        "mapk_activity_score": signature_score(expr, MAPK_GENES).values,
        "thyroid_diff_score": signature_score(expr, THYROID_DIFF).values,
    })
    mapk.to_csv(TAB / "mapk_activity_score.tsv", sep="\t", index=False)

    immune = pd.DataFrame({"sample_id": expr.index, "cluster": meta["v17_dark_cluster"].values})
    for name, genes in IMMUNE_SETS.items():
        immune[name] = signature_score(expr, genes).values
    immune.to_csv(TAB / "cibersort_immune_fractions.tsv", sep="\t", index=False)

    top_gsea = gsea.sort_values("NES", ascending=False).head(15) if not gsea.empty else pd.DataFrame({"Term": [], "NES": []})
    px.bar(top_gsea, x="NES", y="Term", color="gene_set_source" if "gene_set_source" in top_gsea.columns else None, orientation="h", title="Top enriched pathways").write_html(
        FIG / "gsea_top_pathways_bar.html", include_plotlyjs="cdn"
    )
    if hallmark_scores.shape[1] > 4:
        keep = [c for c in hallmark_scores.columns if c not in {"sample_id", "cluster"}][:20]
        hm = hallmark_scores.sort_values("cluster")
        px.imshow(hm[keep].T, aspect="auto", title="Hallmark score heatmap").write_html(FIG / "hallmark_heatmap_per_sample.html", include_plotlyjs="cdn")
    else:
        px.scatter(title="Hallmark score heatmap unavailable").write_html(FIG / "hallmark_heatmap_per_sample.html", include_plotlyjs="cdn")
    immune_long = immune.melt(id_vars=["sample_id", "cluster"], var_name="cell_state", value_name="score")
    px.box(immune_long, x="cell_state", y="score", color="cluster", points="all", title="Immune composition proxy").write_html(
        FIG / "immune_composition_stacked.html", include_plotlyjs="cdn"
    )
    px.violin(mapk, x="cluster", y="mapk_activity_score", color="cluster", box=True, points="all", title="MAPK activity score").write_html(
        FIG / "mapk_activity_violin.html", include_plotlyjs="cdn"
    )
    px.violin(mapk, x="cluster", y="thyroid_diff_score", color="cluster", box=True, points="all", title="Thyroid differentiation score").write_html(
        FIG / "thyroid_diff_score_violin.html", include_plotlyjs="cdn"
    )

    hall_sig = 0
    if hallmark_scores is not None and hallmark_scores.shape[1] > 4:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pvals = []
            cols = [c for c in hallmark_scores.columns if c not in {"sample_id", "cluster"}]
            for c in cols:
                p = stats.ttest_ind(
                    hallmark_scores.loc[hallmark_scores["cluster"] == "DM1", c],
                    hallmark_scores.loc[hallmark_scores["cluster"] == "DM2", c],
                    equal_var=False,
                    nan_policy="omit",
                ).pvalue
                pvals.append(p)
            hall_sig = int((bh_fdr(np.asarray(pvals)) < 0.05).sum())
    summary = {
        "deg_genes": int(len(deg)),
        "gsea_terms_fdr_lt_0_05": int((gsea.get("FDR q-val", pd.Series(dtype=float)) < 0.05).sum()) if not gsea.empty else 0,
        "hallmark_fdr_lt_0_05": hall_sig,
        "mapk_dm2_minus_dm1": float(mapk.loc[mapk["cluster"] == "DM2", "mapk_activity_score"].mean() - mapk.loc[mapk["cluster"] == "DM1", "mapk_activity_score"].mean()),
    }
    write_json(TAB / "biology_summary.json", summary)
    log_line(LOG, f"v17p2 layer3 done {summary}")


if __name__ == "__main__":
    main()
