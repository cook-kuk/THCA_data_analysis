#!/usr/bin/env python3
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from v17_common import FIG, LOG, META, TAB, EXTRA_DRIVER_GENES, load_local_tcga_mutation_summary, load_sample_master, log_line, normalize_barcode


FUSION_KEYS = {
    "RET_fusion": ["RET"],
    "NTRK_fusion": ["NTRK1", "NTRK3"],
    "ALK_fusion": ["ALK"],
    "PAX8PPARG_fusion": ["PAX8", "PPARG"],
    "THADA_fusion": ["THADA"],
}


def parse_fusions() -> pd.DataFrame:
    proxy = pd.read_csv(META / "v3_fusion_anchor_tcga.tsv", sep="\t")
    proxy = proxy[proxy["dataset"] == "TCGA-THCA"].copy()
    proxy["tcga12"] = proxy["sample_id"].astype(str).map(normalize_barcode)
    proxy["fusion_pair"] = np.nan
    proxy["fusion_classes"] = "other"
    mapping = {
        "RET_fusion": ["RET"],
        "NTRK_fusion": ["NTRK1", "NTRK3"],
        "ALK_fusion": ["ALK"],
        "PAX8PPARG_fusion": ["PAX8", "PPARG"],
        "THADA_fusion": ["THADA"],
    }
    # v3 proxy is coarse. Keep it explicit: only "other" is observed locally here.
    proxy["fusion_source"] = "v3_fusion_anchor_proxy"
    return proxy[["sample_id", "tcga12", "v3_anchor_6class", "fusion_pair", "fusion_classes", "fusion_source"]]


def assign_driver_v17(row: pd.Series) -> str:
    genes = set(row.get("mutation_genes", []) or [])
    fusions = set(row.get("fusion_classes", []) or [])
    tert = row.get("tert_status", "wildtype")
    braf_v600e = pd.to_numeric(row.get("braf_v600e", 0), errors="coerce")
    ras_hotspot = pd.to_numeric(row.get("ras_hotspot", 0), errors="coerce")
    tert_promoter = pd.to_numeric(row.get("tert_promoter", 0), errors="coerce")
    if (pd.notna(braf_v600e) and int(braf_v600e) == 1) or "BRAF" in genes:
        return "BRAF"
    if (pd.notna(ras_hotspot) and int(ras_hotspot) == 1) or genes & {"HRAS", "KRAS", "NRAS"}:
        return "RAS"
    if tert == "mutated" or (pd.notna(tert_promoter) and int(tert_promoter) == 1) or "TERT" in genes:
        return "TERT"
    if "RET_fusion" in fusions or "RET" in genes:
        return "RET_fusion"
    if "NTRK_fusion" in fusions or genes & {"NTRK1", "NTRK3"}:
        return "NTRK_fusion"
    if "ALK_fusion" in fusions or "ALK" in genes:
        return "ALK_fusion"
    if "PAX8PPARG_fusion" in fusions or {"PAX8", "PPARG"} <= genes:
        return "PAX8PPARG"
    if genes & {"DICER1", "EIF1AX", "PPM1D"}:
        return "DICER1_EIF1AX_PPM1D"
    if "TP53" in genes:
        return "TP53"
    return row.get("driver_anchor", "unknown")


def main() -> None:
    log_line(LOG, "task3 start")
    sm = load_sample_master().copy()
    sm["tcga12"] = sm["sample_id"].map(normalize_barcode)
    tert_path = TAB / "sample_master_v17_tert.tsv"
    if tert_path.exists():
        tert_df = pd.read_csv(tert_path, sep="\t")[["sample_id", "tert_status"]].drop_duplicates()
        sm = sm.merge(tert_df, on="sample_id", how="left")
    else:
        sm["tert_status"] = np.where(sm["driver_anchor"] == "TERT", "mutated", "wildtype")

    fus = parse_fusions()
    mut = load_local_tcga_mutation_summary(LOG)
    mut["mutation_genes"] = mut["mutation_genes"].fillna("").astype(str).str.split(";").apply(lambda xs: [x for x in xs if x])
    mut = mut.drop(columns=["sample_id"], errors="ignore")
    fus_df = fus[["sample_id", "tcga12", "fusion_pair", "fusion_classes", "fusion_source", "v3_anchor_6class"]].copy()
    fus_df.to_csv(TAB / "fusion_calls_per_sample.tsv", sep="\t", index=False)

    tcga_fus = fus_df.groupby("tcga12").agg(
        fusion_pair=("fusion_pair", lambda x: ";".join(sorted({str(v) for v in x if pd.notna(v)}))),
        fusion_classes=("fusion_classes", lambda x: sorted({str(v) for v in x if pd.notna(v) and str(v) != ""})),
        v3_anchor_6class=("v3_anchor_6class", lambda x: next((str(v) for v in x if pd.notna(v)), "other")),
    ).reset_index()
    merged = sm.merge(mut, on="tcga12", how="left").merge(tcga_fus, on="tcga12", how="left")
    merged["fusion_classes"] = merged["fusion_classes"].apply(lambda x: x if isinstance(x, list) else [])
    merged["mutation_genes"] = merged["mutation_genes"].apply(lambda x: x if isinstance(x, list) else [])
    merged["v17_data_source"] = np.where(merged["dataset"].eq("TCGA-THCA"), "local_maf_plus_v3_proxy", "sample_master_v3_only")
    merged["driver_anchor_v17"] = merged.apply(assign_driver_v17, axis=1)
    # External cohorts currently do not have matched mutation/fusion calls here.
    ext_mask = merged["dataset"] != "TCGA-THCA"
    merged.loc[ext_mask, "driver_anchor_v17"] = merged.loc[ext_mask, "driver_anchor"]

    dark_path = TAB / "dark_matter_cohort.tsv"
    if dark_path.exists():
        dark = pd.read_csv(dark_path, sep="\t")[["sample_id", "v17_dark_cluster"]].drop_duplicates()
        merged = merged.merge(dark, on="sample_id", how="left")
    merged = merged.sort_values(["sample_id", "dataset"]).drop_duplicates("sample_id", keep="first").reset_index(drop=True)
    merged.to_csv(TAB / "sample_master_v17_full.tsv", sep="\t", index=False)

    summary = merged.groupby("driver_anchor_v17").agg(n_samples=("sample_id", "size")).reset_index().sort_values("n_samples", ascending=False)
    summary.to_csv(TAB / "driver_landscape_v17_summary.tsv", sep="\t", index=False)

    pie = px.pie(summary, values="n_samples", names="driver_anchor_v17", title="v17 driver landscape")
    pie.write_html(FIG / "driver_landscape_pie.html", include_plotlyjs="cdn")
    bar = px.bar(summary, x="driver_anchor_v17", y="n_samples", color="driver_anchor_v17", title="v17 driver landscape counts")
    bar.write_html(FIG / "driver_landscape_bar.html", include_plotlyjs="cdn")

    tcga_tumor = merged[(merged["dataset"] == "TCGA-THCA") & (merged["normal_vs_tumor"] != "normal")].copy()
    top_groups = summary.head(12)["driver_anchor_v17"].tolist()
    top_tcga = tcga_tumor[tcga_tumor["driver_anchor_v17"].isin(top_groups)].copy()
    mat = pd.crosstab(top_tcga["sample_id"], top_tcga["driver_anchor_v17"])
    mat = (mat > 0).astype(int).T
    onco = go.Figure(data=go.Heatmap(z=mat.values, x=mat.columns, y=mat.index, colorscale="Viridis"))
    onco.update_layout(title="v17 driver landscape oncoplot (TCGA tumors)")
    onco.write_html(FIG / "driver_landscape_oncoplot.html", include_plotlyjs="cdn")

    if dark_path.exists():
        if "v17_dark_cluster" in merged.columns:
            over = merged[pd.notna(merged["v17_dark_cluster"])].copy()
        else:
            dark = pd.read_csv(dark_path, sep="\t")[["sample_id", "v17_dark_cluster"]]
            over = merged.merge(dark, on="sample_id", how="inner")
        if not over.empty:
            ctab = pd.crosstab(over["v17_dark_cluster"], over["driver_anchor_v17"]).reset_index()
            ctab.to_csv(TAB / "dark_matter_fusion_overlay.tsv", sep="\t", index=False)
            overlay = px.imshow(ctab.set_index("v17_dark_cluster"), aspect="auto", title="Dark matter cluster by v17 driver")
            overlay.write_html(FIG / "dark_matter_fusion_overlay.html", include_plotlyjs="cdn")
            stack = px.bar(
                over.groupby(["v17_dark_cluster", "driver_anchor_v17"]).size().reset_index(name="n"),
                x="v17_dark_cluster",
                y="n",
                color="driver_anchor_v17",
                title="Dark matter cluster composition under v17 drivers",
            )
            stack.write_html(FIG / "dark_matter_fusion_overlay_bar.html", include_plotlyjs="cdn")

    confusion = pd.crosstab(merged["molecular_subtype"], merged["driver_anchor_v17"]).reset_index()
    confusion.to_csv(TAB / "v14_vs_v17_confusion.tsv", sep="\t", index=False)
    conf_fig = px.imshow(confusion.set_index("molecular_subtype"), aspect="auto", title="v14 molecular subtype vs v17 driver map")
    conf_fig.write_html(FIG / "v14_vs_v17_confusion.html", include_plotlyjs="cdn")

    log_line(LOG, f"task3 done drivers={json.dumps(summary.head(10).to_dict(orient='records'))}")


if __name__ == "__main__":
    main()
