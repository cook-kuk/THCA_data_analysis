#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from v17_common import FIG, TAB


def main() -> None:
    dark = pd.read_csv(TAB / "dark_matter_cluster_clinical.tsv", sep="\t")
    dark_bar = px.bar(
        dark,
        x="v17_dark_cluster",
        y="n_samples",
        color="v17_dark_cluster",
        text="n_samples",
        title="Dark Matter cluster size",
    )
    dark_bar.write_html(FIG / "dark_matter_cluster_size.html", include_plotlyjs="cdn")

    dark_prog = dark.melt(
        id_vars=["v17_dark_cluster"],
        value_vars=["tds16_mean", "rai_mean", "age_mean"],
        var_name="metric",
        value_name="value",
    )
    fig = px.bar(dark_prog, x="metric", y="value", color="v17_dark_cluster", barmode="group", title="DM cluster summary metrics")
    fig.write_html(FIG / "dark_matter_cluster_metrics.html", include_plotlyjs="cdn")

    quad = pd.read_csv(TAB / "quad_group_summary.tsv", sep="\t")
    qbar = px.bar(
        quad,
        x="quad_group",
        y="n_samples",
        color="quad_group",
        text="n_samples",
        title="Quad-group sample counts",
    )
    qbar.write_html(FIG / "quad_group_counts.html", include_plotlyjs="cdn")

    stage = pd.read_csv(TAB / "quad_group_stage_distribution.tsv", sep="\t")
    stage_long = stage.melt(id_vars=["quad_group"], var_name="stage", value_name="n")
    stage_long = stage_long[stage_long["n"] > 0].copy()
    sfig = px.bar(stage_long, x="quad_group", y="n", color="stage", title="Stage composition by v17 group")
    sfig.write_html(FIG / "quad_group_stage_stacked.html", include_plotlyjs="cdn")

    driver = pd.read_csv(TAB / "driver_landscape_v17_summary.tsv", sep="\t")
    dtop = driver.head(8).copy()
    dfig = px.treemap(dtop, path=[px.Constant("v17"), "driver_anchor_v17"], values="n_samples", title="Top v17 driver landscape treemap")
    dfig.write_html(FIG / "driver_landscape_treemap.html", include_plotlyjs="cdn")

    dial = pd.read_csv(TAB / "dial_audit_v17.tsv", sep="\t")
    dial_w = dial.pivot_table(index=["cohort", "cluster"], columns="ComBat_applied", values="identifiability").reset_index()
    if {"yes", "no"} <= set(dial_w.columns):
        dial_w["delta_identifiability"] = dial_w["yes"] - dial_w["no"]
        dcomp = px.bar(
            dial_w,
            x="cohort",
            y="delta_identifiability",
            color="cluster",
            barmode="group",
            title="Identifiability change after ComBat",
        )
        dcomp.write_html(FIG / "dial_identifiability_delta.html", include_plotlyjs="cdn")

    traj = pd.read_csv(TAB / "trajectory_pseudotime.tsv", sep="\t")
    thist = px.histogram(traj, x="pseudotime", color="dataset", nbins=40, marginal="box", title="Pseudotime distribution by cohort")
    thist.write_html(FIG / "trajectory_pseudotime_histogram.html", include_plotlyjs="cdn")

    top_dyn = pd.read_csv(TAB / "trajectory_dynamic_genes.tsv", sep="\t")
    top_pos = top_dyn.sort_values("spearman_rho", ascending=False).head(8)
    top_neg = top_dyn.sort_values("spearman_rho", ascending=True).head(8)
    dyn = pd.concat([top_pos, top_neg], axis=0).drop_duplicates("gene")
    gfig = px.bar(dyn.sort_values("spearman_rho"), x="spearman_rho", y="gene", orientation="h", color="spearman_rho", title="Top dynamic genes along pseudotime")
    gfig.write_html(FIG / "trajectory_dynamic_genes_bar.html", include_plotlyjs="cdn")

    if "v17_dark_cluster" in traj.columns:
        t2 = traj[pd.notna(traj["v17_dark_cluster"])].copy()
        if not t2.empty:
            box = px.box(t2, x="v17_dark_cluster", y="pseudotime", color="v17_dark_cluster", points="all", title="Pseudotime by dark-matter cluster")
            box.write_html(FIG / "trajectory_dark_cluster_pseudotime.html", include_plotlyjs="cdn")


if __name__ == "__main__":
    main()
