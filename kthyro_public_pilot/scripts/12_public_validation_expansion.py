#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT.parent
OUT = ROOT / "results" / "public_validation_expansion"
TABLES = OUT / "tables"
FIGS = OUT / "figures"


def read(path: str) -> pd.DataFrame:
    p = PARENT / path
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t")


def cohens_d(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / pooled) if pooled and np.isfinite(pooled) else np.nan


def mwu(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def strength(effect: float, p: float, expected_match: bool | None = None) -> str:
    if expected_match is False:
        return "opposite"
    if not np.isfinite(effect):
        return "not_tested"
    ae = abs(effect)
    if ae >= 1.0 and (not np.isfinite(p) or p < 0.05):
        return "strong"
    if ae >= 0.5 and (not np.isfinite(p) or p < 0.20):
        return "moderate"
    if ae >= 0.3:
        return "weak"
    return "weak_or_none"


def add(rows: list[dict], dataset: str, validation_layer: str, axis: str, endpoint: str, n_a, n_b, effect, p, direction, caveat, source_file):
    rows.append(
        {
            "dataset": dataset,
            "validation_layer": validation_layer,
            "axis": axis,
            "endpoint": endpoint,
            "n_group_a": n_a,
            "n_group_b": n_b,
            "effect_size_cohens_d_or_reported": effect,
            "p_value": p,
            "direction": direction,
            "support_strength": strength(effect, p, None if direction != "opposite" else False),
            "caveat": caveat,
            "source_file": source_file,
        }
    )


def ingest_external_bulk(rows: list[dict]) -> None:
    f = "project/results/p_external_expression_validation/external_score_tests.tsv"
    df = read(f)
    if df.empty:
        return
    axis_map = {
        "RAI_8_score": "RAI differentiation / redifferentiation",
        "DM1_like_score": "Aggressive dedifferentiation",
        "THYROID_NONOVERLAP_score": "Thyroid lineage differentiation",
        "TDS_like_score": "Thyroid differentiation score",
        "STAT3_AP1_DNMT_score": "Stress/epigenetic dedifferentiation",
        "TACSTD2_z": "Tumor epithelial/aggressive marker",
    }
    for _, r in df.iterrows():
        for score, axis in axis_map.items():
            dc = f"{score}_d"
            pc = f"{score}_p"
            if dc not in df.columns:
                continue
            eff = pd.to_numeric(r[dc], errors="coerce")
            p = pd.to_numeric(r.get(pc), errors="coerce")
            add(
                rows,
                r["dataset"],
                "external_bulk_expression",
                axis,
                r["contrast"],
                r["n_group1"],
                r["n_group2"],
                eff,
                p,
                "reported",
                "Independent bulk/microarray validation; platform and metadata differ by cohort.",
                f,
            )


def ingest_gpl570(rows: list[dict]) -> None:
    f = "project/results/p_gpl570_validation/gpl570_score_tests.tsv"
    df = read(f)
    if df.empty:
        return
    for _, r in df.iterrows():
        test = str(r["test"])
        if "RAI_8_score" not in test and "TDS_like_score" not in test and "DM1_like_score" not in test:
            continue
        if "RAI_8_score" in test or "TDS_like_score" in test:
            axis = "RAI/thyroid differentiation"
        else:
            axis = "Aggressive dedifferentiation"
        add(
            rows,
            r["dataset"],
            "external_GPL570_bulk",
            axis,
            test.split(" :: ", 1)[-1],
            r["n_x"],
            r["n_y"],
            pd.to_numeric(r["cohens_d"], errors="coerce"),
            pd.to_numeric(r["MW_p_two_sided"], errors="coerce"),
            "reported",
            "Older array cohorts; good for directionality, not therapeutic response.",
            f,
        )


def ingest_rai_direct(rows: list[dict]) -> None:
    f = "project/results/aggressive_sprint_2026_05_06/gse151179_rai_scores.tsv"
    df = read(f)
    if df.empty:
        return
    tumor = df[df.get("is_tumor", False).astype(str).isin(["True", "true", "1"])] if "is_tumor" in df else df
    for score in ["rai8_z_mean", "rai6_z_mean", "SLC5A5", "TG", "TPO", "TSHR"]:
        if score not in tumor.columns:
            continue
        refr = tumor[tumor["response_refractory"].astype(str).isin(["1", "True", "true"])]
        avid = tumor[tumor["response_refractory"].astype(str).isin(["0", "False", "false"])]
        add(
            rows,
            "GSE151179",
            "direct_RAI_avid_refractory_bulk",
            "RAI differentiation / redifferentiation",
            f"refractory_vs_avid::{score}",
            len(refr),
            len(avid),
            cohens_d(refr[score], avid[score]),
            mwu(refr[score], avid[score]),
            "lower in refractory expected",
            "Direct public RAI label but small avid group and array platform; do not transfer as clinical predictor.",
            f,
        )

    f2 = "project/results/paper3_ici_track_b_lite/gse151179_rai_after_vs_before.tsv"
    d2 = read(f2)
    for _, r in d2.iterrows():
        if r["module"] not in {"HLA_class_I", "HLA_class_II", "IFNG_T_cell_inflamed", "myeloid_suppressive", "thyroid_differentiation"}:
            continue
        axis = "HLA/APM immune visibility" if "HLA" in r["module"] or "IFNG" in r["module"] else r["module"]
        add(
            rows,
            "GSE151179",
            "RAI_before_after_bulk",
            axis,
            f"after_vs_before::{r['module']}",
            r["n_after"],
            r["n_before"],
            r["cohens_d_after_minus_before"],
            r["p"],
            "after-minus-before",
            "Before/after RAI collection is not a controlled perturbation experiment.",
            f2,
        )


def ingest_landa_gse76039(rows: list[dict]) -> None:
    f = "project/results/paper3_ici_track_b_lite/scores_per_cohort/GSE76039_module_scores.tsv"
    predf = "project/results/v17_realfix/R3A_gse76039_predictions.tsv"
    df = read(f)
    pred = read(predf)
    if df.empty or pred.empty:
        return
    m = df.merge(pred, on="sample_id", how="inner")
    atc = m[m["y_ATC"].eq(1)]
    pdtc = m[m["y_ATC"].eq(0)]
    module_axis = {
        "thyroid_differentiation": "RAI/thyroid differentiation",
        "HLA_class_I": "HLA/APM immune visibility",
        "HLA_class_II": "APC/HLA-II niche",
        "IFNG_T_cell_inflamed": "IFN/cytotoxic inflammation",
        "myeloid_suppressive": "Myeloid suppression",
        "checkpoint_exhaustion": "Checkpoint/exhaustion",
    }
    for mod, axis in module_axis.items():
        add(
            rows,
            "GSE76039",
            "PDTC_ATC_bulk",
            axis,
            f"ATC_vs_PDTC::{mod}",
            len(atc),
            len(pdtc),
            cohens_d(atc[mod], pdtc[mod]),
            mwu(atc[mod], pdtc[mod]),
            "ATC-minus-PDTC",
            "PDTC/ATC cohort validates advanced-disease direction, not treatment response.",
            f,
        )


def ingest_geomx(rows: list[dict]) -> None:
    f = "project/results/spatial_full_2026_05_06/GSE301163_per_ROI_scores.tsv"
    df = read(f)
    if df.empty:
        return
    cancer = df[df["histology"].isin(["microPTC", "PDTC"])]
    normal = df[df["histology"].isin(["Normal thyroid", "Normal adjacent", "Non-neoplastic nodular area"])]
    for score, axis in {
        "RAI_8": "RAI differentiation / redifferentiation",
        "HLA_II": "APC/HLA-II niche",
        "Stromal": "CAF/stromal barrier",
        "Hypoxia": "Hypoxia/drug-delivery proxy",
        "Cell_cycle": "Proliferation/aggressive stress",
    }.items():
        add(
            rows,
            "GSE301163",
            "GeoMx_ROI_spatial",
            axis,
            f"cancer_ROI_vs_non_neoplastic_ROI::{score}",
            len(cancer),
            len(normal),
            cohens_d(cancer[score], normal[score]),
            mwu(cancer[score], normal[score]),
            "cancer-minus-nonneoplastic",
            "GeoMx ROI is spatial transcriptomics but not Visium; labels/ROI selection matter.",
            f,
        )
    panck = df[df["label"].eq("PanCK+")]
    vim = df[df["label"].eq("VIM+")]
    for score, axis in {"RAI_8": "RAI epithelial/stromal compartment", "Stromal": "CAF/stromal barrier", "HLA_II": "APC/HLA-II niche"}.items():
        add(
            rows,
            "GSE301163",
            "GeoMx_ROI_compartment",
            axis,
            f"PanCK_vs_VIM::{score}",
            len(panck),
            len(vim),
            cohens_d(panck[score], vim[score]),
            mwu(panck[score], vim[score]),
            "PanCK-minus-VIM",
            "Compartment contrast supports ROI logic, not patient outcome.",
            f,
        )


def ingest_scrna(rows: list[dict]) -> None:
    for f in [
        "project/results/nature_cancer_feasibility/gse232237_score_contrasts.tsv",
        "project/results/nature_cancer_feasibility/gse232237_epithelial_score_contrasts.tsv",
    ]:
        df = read(f)
        if df.empty:
            continue
        for _, r in df.iterrows():
            if r["score"] not in {"RAI_8_score", "DM1_like_score", "THYROID_NONOVERLAP_score", "TDS_like_score"}:
                continue
            axis = "RAI/thyroid differentiation" if "RAI" in r["score"] or "TDS" in r["score"] or "THYROID" in r["score"] else "Aggressive dedifferentiation"
            add(
                rows,
                "GSE232237",
                "Korean_scRNA_pseudobulk",
                axis,
                f"{r['contrast']}::{r['score']}",
                r.get("n_atc"),
                r.get("n_ptc"),
                r.get("cohen_d_atc_minus_ptc"),
                r.get("p_mwu"),
                "direction_match" if bool(r.get("direction_match")) else "opposite",
                str(r.get("note", "")) + " Public scRNA pseudo-bulk supports cell-state direction, not clinical response.",
                f,
            )

    f2 = "project/results/v17_lu2023/GSE193581_per_sample_summary.tsv"
    df2 = read(f2)
    if not df2.empty and {"histology", "mean_DM"}.issubset(df2.columns):
        atc = df2[df2["histology"].astype(str).eq("ATC")]
        ptc = df2[df2["histology"].astype(str).eq("PTC")]
        if len(atc) and len(ptc):
            add(
                rows,
                "GSE193581",
                "scRNA_per_sample_summary",
                "DM/dark-lineage prior axis",
                "ATC_vs_PTC::mean_DM",
                len(atc),
                len(ptc),
                cohens_d(atc["mean_DM"], ptc["mean_DM"]),
                mwu(atc["mean_DM"], ptc["mean_DM"]),
                "ATC-minus-PTC",
                "Legacy DM score sanity check; not the full K-Thyro axis set.",
                f2,
            )


def ingest_proteomics(rows: list[dict]) -> None:
    f = "project/results/proteogenomic_v1/paper3_mun2025_dediff_layer/module_dediff_trend.tsv"
    df = read(f)
    if not df.empty:
        for _, r in df.iterrows():
            add(
                rows,
                "Mun2025 thyroid proteogenomic layer",
                "proteomics_dedifferentiation_trend",
                r["module"],
                "PTC_to_PDTC_to_ATC_dedifferentiation_trend",
                np.nan,
                np.nan,
                r["spearman_dediff_axis"],
                r["p_trend"],
                "dediff-axis spearman",
                "Protein-layer support for module direction; publication-specific processed supplements.",
                f,
            )
    f2 = "project/results/proteogenomic_v1/processed/cellrepmed2026_CC_x_RAIsensitivity.tsv"
    df2 = read(f2)
    if not df2.empty:
        # Simple chi-square effect proxy: refractory fraction by subtype range.
        tmp = df2.copy()
        tmp["refractory_fraction"] = tmp["Refractory"] / (tmp["Avid"] + tmp["Refractory"])
        add(
            rows,
            "CellRepMed2026 thyroid multiomics",
            "clinical_RAIsensitivity_subtype_table",
            "RAI sensitivity / molecular subtype",
            "CC_subtype_vs_RAI_sensitivity_range",
            int(tmp["Avid"].sum()),
            int(tmp["Refractory"].sum()),
            float(tmp["refractory_fraction"].max() - tmp["refractory_fraction"].min()),
            np.nan,
            "range of refractory fraction",
            "Subtype-by-RAI table only; not reprocessed individual-level omics in this expansion.",
            f2,
        )


def make_figures(sig: pd.DataFrame) -> None:
    if sig.empty:
        return
    FIGS.mkdir(parents=True, exist_ok=True)
    plot = sig.copy()
    plot["signed_abs_effect"] = pd.to_numeric(plot["effect_size_cohens_d_or_reported"], errors="coerce")
    # Cap for readable heatmap
    plot["signed_abs_effect_cap"] = plot["signed_abs_effect"].clip(-4, 4)
    pivot = plot.pivot_table(
        index=["validation_layer", "dataset"],
        columns="axis",
        values="signed_abs_effect_cap",
        aggfunc=lambda x: np.nanmedian(x),
    )
    fig, ax = plt.subplots(figsize=(16, max(6, 0.35 * len(pivot))), facecolor="#111318")
    sns.heatmap(pivot, cmap="vlag", center=0, vmin=-4, vmax=4, linewidths=0.3, linecolor="#222831", ax=ax, cbar_kws={"label": "Effect size / signed statistic"})
    ax.set_title("External public-data validation expansion: axis directionality", color="white", fontsize=16, weight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(colors="white", labelsize=8)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_color("white")
    cbar = ax.collections[0].colorbar
    cbar.ax.yaxis.label.set_color("white")
    cbar.ax.tick_params(colors="white")
    fig.savefig(FIGS / "external_public_validation_axis_heatmap.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)

    counts = sig.groupby(["validation_layer", "support_strength"]).size().reset_index(name="n")
    fig, ax = plt.subplots(figsize=(13, 7), facecolor="#111318")
    sns.barplot(data=counts, y="validation_layer", x="n", hue="support_strength", ax=ax, palette={"strong": "#34d399", "moderate": "#60a5fa", "weak": "#fbbf24", "weak_or_none": "#9ca3af", "opposite": "#ef4444", "not_tested": "#6b7280"})
    ax.set_title("External validation support counts by public-data layer", color="white", fontsize=16, weight="bold")
    ax.set_xlabel("Number of tested contrasts/modules")
    ax.set_ylabel("")
    ax.set_facecolor("#111318")
    ax.tick_params(colors="white")
    ax.xaxis.label.set_color("white")
    leg = ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    for t in leg.get_texts():
        t.set_color("white")
    for spine in ax.spines.values():
        spine.set_color("#374151")
    fig.savefig(FIGS / "external_public_validation_support_counts.png", dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    ingest_external_bulk(rows)
    ingest_gpl570(rows)
    ingest_rai_direct(rows)
    ingest_landa_gse76039(rows)
    ingest_geomx(rows)
    ingest_scrna(rows)
    ingest_proteomics(rows)
    sig = pd.DataFrame(rows)
    sig.to_csv(TABLES / "external_public_validation_signal_matrix.tsv", sep="\t", index=False)
    if not sig.empty:
        summary = (
            sig.groupby(["validation_layer", "dataset", "support_strength"])
            .size()
            .reset_index(name="n_tests")
            .sort_values(["validation_layer", "dataset", "support_strength"])
        )
    else:
        summary = pd.DataFrame(columns=["validation_layer", "dataset", "support_strength", "n_tests"])
    summary.to_csv(TABLES / "external_public_validation_layer_summary.tsv", sep="\t", index=False)
    make_figures(sig)
    print(f"Wrote {len(sig)} validation tests to {TABLES / 'external_public_validation_signal_matrix.tsv'}")
    print(f"Wrote figures to {FIGS}")


if __name__ == "__main__":
    main()
