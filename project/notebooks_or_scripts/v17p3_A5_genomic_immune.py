#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px

from v17_common import load_local_tcga_mutation_summary
from v17p3_common import CYTOLYTIC, FIG, IMMUNE_EVASION, LOG, TAB, dm_labels, load_dark, log_line, score_expr, write_json


def main() -> None:
    log_line(LOG, "v17p3 A5 start")
    expr, meta = load_dark()
    expr.columns = expr.columns.astype(str).str.upper()
    mut = load_local_tcga_mutation_summary(LOG)
    mut["tcga12"] = mut["sample_id"].astype(str).str[:12]
    meta["tcga12"] = meta["sample_id"].astype(str).str[:12]
    mut["mutation_n"] = np.where(mut["mutation_genes"].fillna("").astype(str).str.len() > 0, mut["mutation_genes"].astype(str).str.split(",").str.len(), 0)
    burden = mut.groupby("tcga12").agg(total_mutations=("mutation_n", "sum"), unique_driver_strings=("mutation_genes", "nunique")).reset_index()
    df = meta.merge(burden, on="tcga12", how="left")
    df["cytolytic_score"] = score_expr(expr, CYTOLYTIC).values
    for g in IMMUNE_EVASION:
        if g in expr.columns:
            df[g] = expr[g].values
    df.to_csv(TAB / "A5_tmb_msi_per_cluster.tsv", sep="\t", index=False)
    imm_cols = [g for g in IMMUNE_EVASION if g in df.columns]
    imm = df[["sample_id", "v17_dark_cluster"] + imm_cols + ["cytolytic_score"]].copy()
    imm.to_csv(TAB / "A5_immune_evasion_genes.tsv", sep="\t", index=False)
    pd.DataFrame([{"metric": "CNV_burden", "status": "missing_local_source"}, {"metric": "aneuploidy", "status": "missing_local_source"}]).to_csv(
        TAB / "A5_cnv_aneuploidy_per_cluster.tsv", sep="\t", index=False
    )
    pd.DataFrame([{"immune_subtype": "missing", "status": "Thorsson subtype table not local"}]).to_csv(TAB / "A5_immune_subtype_distribution.tsv", sep="\t", index=False)

    px.box(df, x="v17_dark_cluster", y="total_mutations", color="v17_dark_cluster", points="all", title="Mutation burden by cluster").write_html(
        FIG / "A5_genomic_instability_panel.html", include_plotlyjs="cdn"
    )
    long = imm.melt(id_vars=["sample_id", "v17_dark_cluster"], var_name="feature", value_name="value")
    px.box(long, x="feature", y="value", color="v17_dark_cluster", points="all", title="Immune landscape").write_html(
        FIG / "A5_immune_landscape.html", include_plotlyjs="cdn"
    )
    if imm_cols:
        px.imshow(imm.sort_values("v17_dark_cluster")[imm_cols].T, aspect="auto", title="Immune evasion heatmap").write_html(
            FIG / "A5_immune_evasion_heatmap.html", include_plotlyjs="cdn"
        )
    summary = {"tmb_dm2_minus_dm1": float(df.loc[df['v17_dark_cluster']=='DM2','total_mutations'].mean() - df.loc[df['v17_dark_cluster']=='DM1','total_mutations'].mean())}
    write_json(TAB / "A5_summary.json", summary)
    log_line(LOG, f"v17p3 A5 done {summary}")


if __name__ == "__main__":
    main()
