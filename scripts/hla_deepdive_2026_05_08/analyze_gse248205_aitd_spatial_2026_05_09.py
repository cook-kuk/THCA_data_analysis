#!/usr/bin/env python3
"""Analyze GSE248205 AITD Visium data for HLA/AP spatial mechanism support.

Paper 4 use: tissue mechanism support for AITD HLA genetics.
Paper 2 boundary: this is non-cancer AITD spatial expression evidence only;
it must not be framed as thyroid-cancer HLA allele/genotype evidence.
"""

from __future__ import annotations

import gzip
import json
import math
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.io import mmread
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "project/data/external/GSE248205/Processed files"
BASE = ROOT / "project/results/hla_two_paper_synthesis_2026_05_09"
OUT = BASE / "gse248205_aitd_spatial_validation"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
REPORT = ROOT / "project/reports/2026_05_09_GSE248205_AITD_SPATIAL_VALIDATION_KR.md"


MODULES = {
    "HLA_I": [
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "HLA-E",
        "HLA-F",
        "HLA-G",
        "B2M",
        "TAP1",
        "TAP2",
        "TAPBP",
        "PSMB8",
        "PSMB9",
        "PSMB10",
        "NLRC5",
    ],
    "HLA_II_AP": [
        "HLA-DRA",
        "HLA-DRB1",
        "HLA-DRB3",
        "HLA-DRB4",
        "HLA-DRB5",
        "HLA-DQA1",
        "HLA-DQA2",
        "HLA-DQB1",
        "HLA-DQB2",
        "HLA-DPA1",
        "HLA-DPB1",
        "HLA-DMA",
        "HLA-DMB",
        "HLA-DOA",
        "HLA-DOB",
        "CD74",
        "CIITA",
        "RFX5",
    ],
    "B_TLS": [
        "MS4A1",
        "CD79A",
        "CD79B",
        "CD19",
        "CD22",
        "CD27",
        "CD38",
        "MZB1",
        "JCHAIN",
        "IGHG1",
        "IGHG2",
        "IGKC",
        "CXCL13",
        "CCL19",
        "CCL21",
        "LTB",
        "LTBR",
        "CCR7",
        "POU2AF1",
        "AICDA",
    ],
    "T_IFNG": [
        "CD3D",
        "CD3E",
        "CD3G",
        "CD4",
        "CD8A",
        "CD8B",
        "IFNG",
        "GZMB",
        "PRF1",
        "CXCL9",
        "CXCL10",
        "CXCL11",
        "STAT1",
        "IRF1",
    ],
    "Myeloid_DC": [
        "LST1",
        "LYZ",
        "FCER1A",
        "CLEC10A",
        "CD1C",
        "LAMP3",
        "ITGAX",
        "CD68",
        "FCGR3A",
        "C1QA",
        "C1QB",
    ],
    "Thyrocyte": [
        "TG",
        "TPO",
        "EPCAM",
        "KRT8",
        "KRT18",
        "SLC5A5",
        "TSHR",
        "FOXE1",
        "NKX2-1",
    ],
    "CD74_MIF_axis": ["CD74", "MIF"],
}

MODULE_ORDER = ["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "Myeloid_DC", "Thyrocyte", "CD74_MIF_axis", "AP_TLS_composite"]
PRIMARY_MODULES = ["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "Myeloid_DC", "CD74_MIF_axis", "AP_TLS_composite"]
SAMPLE_ORDER = ["C1", "C2", "HT1", "HT2", "HT3", "GD1", "GD2", "GD3"]
GROUP_ORDER = ["Control", "HT", "GD"]
PALETTE = {"Control": "#72B7B2", "HT": "#54A24B", "GD": "#E45756"}


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)


def sample_group(sample: str) -> str:
    if sample.startswith("C"):
        return "Control"
    if sample.startswith("HT"):
        return "HT"
    if sample.startswith("GD"):
        return "GD"
    raise ValueError(f"Unknown sample group for {sample}")


def target_genes() -> set[str]:
    genes: set[str] = set()
    for gene_list in MODULES.values():
        genes.update(gene_list)
    return genes


def read_gzip_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n") for line in fh]


def read_sample(sample: str, genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_dir = DATA / sample
    features = pd.read_csv(
        sample_dir / f"{sample}_features.tsv.gz",
        sep="\t",
        header=None,
        names=["gene_id", "gene", "feature_type"],
        dtype=str,
    )
    barcodes = read_gzip_lines(sample_dir / f"{sample}_barcodes.tsv.gz")
    positions = pd.read_csv(
        sample_dir / f"{sample}_tissue_positions_list.csv",
        header=None,
        names=["barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"],
    )
    barcode_to_idx = {barcode: idx for idx, barcode in enumerate(barcodes)}
    tissue_positions = positions[positions["in_tissue"].eq(1)].copy()
    tissue_positions = tissue_positions[tissue_positions["barcode"].isin(barcode_to_idx)].copy()
    keep_cols = tissue_positions["barcode"].map(barcode_to_idx).to_numpy()

    matrix = mmread(str(sample_dir / f"{sample}_matrix.mtx.gz")).tocsr()
    matrix = matrix[:, keep_cols]
    tissue_positions = tissue_positions.assign(sample=sample, group=sample_group(sample)).reset_index(drop=True)

    totals = np.asarray(matrix.sum(axis=0)).ravel()
    valid = totals > 0
    gene_to_indices: dict[str, list[int]] = {}
    for idx, gene in enumerate(features["gene"].astype(str)):
        if gene in genes:
            gene_to_indices.setdefault(gene, []).append(idx)

    spot_rows = []
    spot_base = tissue_positions[["sample", "group", "barcode", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]].copy()
    for module, module_genes in MODULES.items():
        indices = [idx for gene in module_genes for idx in gene_to_indices.get(gene, [])]
        if indices:
            module_counts = np.asarray(matrix[indices, :].sum(axis=0)).ravel()
            score = np.full(matrix.shape[1], np.nan, dtype=float)
            score[valid] = np.log1p(module_counts[valid] / totals[valid] * 1e4)
            positive = module_counts > 0
        else:
            score = np.full(matrix.shape[1], np.nan, dtype=float)
            positive = np.zeros(matrix.shape[1], dtype=bool)
            module_counts = np.zeros(matrix.shape[1], dtype=float)
        tmp = spot_base.copy()
        tmp["module"] = module
        tmp["module_counts"] = module_counts
        tmp["spot_log1p_cp10k"] = score
        tmp["module_positive"] = positive
        tmp["n_module_genes_available"] = len({features.loc[idx, "gene"] for idx in indices}) if indices else 0
        spot_rows.append(tmp)

    spot_df = pd.concat(spot_rows, ignore_index=True)
    pivot = spot_df.pivot_table(
        index=["sample", "group", "barcode", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"],
        columns="module",
        values="spot_log1p_cp10k",
    ).reset_index()
    pivot["AP_TLS_composite"] = pivot[["HLA_II_AP", "B_TLS"]].mean(axis=1)
    composite_long = pivot[
        ["sample", "group", "barcode", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres", "AP_TLS_composite"]
    ].rename(columns={"AP_TLS_composite": "spot_log1p_cp10k"})
    composite_long["module"] = "AP_TLS_composite"
    composite_long["module_counts"] = np.nan
    composite_long["module_positive"] = composite_long["spot_log1p_cp10k"].gt(0)
    composite_long["n_module_genes_available"] = (
        spot_df[spot_df["module"].eq("HLA_II_AP")]["n_module_genes_available"].max()
        + spot_df[spot_df["module"].eq("B_TLS")]["n_module_genes_available"].max()
    )
    spot_df = pd.concat(
        [
            spot_df,
            composite_long[
                [
                    "sample",
                    "group",
                    "barcode",
                    "array_row",
                    "array_col",
                    "pxl_col_in_fullres",
                    "pxl_row_in_fullres",
                    "module",
                    "module_counts",
                    "spot_log1p_cp10k",
                    "module_positive",
                    "n_module_genes_available",
                ]
            ],
        ],
        ignore_index=True,
    )

    gene_rows = []
    for gene in sorted(genes):
        indices = gene_to_indices.get(gene, [])
        if indices:
            count = float(np.asarray(matrix[indices, :].sum(axis=0)).ravel().sum())
            expr = math.log1p(count / totals.sum() * 1e6) if totals.sum() > 0 else np.nan
            n_probes = len(indices)
        else:
            count = 0.0
            expr = np.nan
            n_probes = 0
        gene_rows.append(
            {
                "sample": sample,
                "group": sample_group(sample),
                "gene": gene,
                "pseudobulk_count": count,
                "log1p_cpm": expr,
                "n_feature_rows": n_probes,
            }
        )
    gene_df = pd.DataFrame(gene_rows)

    manifest = pd.DataFrame(
        [
            {
                "sample": sample,
                "group": sample_group(sample),
                "n_spots_in_tissue": int(matrix.shape[1]),
                "n_features": int(matrix.shape[0]),
                "total_umis_in_tissue": float(totals.sum()),
                "median_umis_per_spot": float(np.median(totals[valid])) if valid.any() else np.nan,
            }
        ]
    )
    return manifest, spot_df, gene_df


def score_sample_modules(gene_df: pd.DataFrame) -> pd.DataFrame:
    matrix = gene_df.pivot(index=["sample", "group"], columns="gene", values="log1p_cpm")
    z = matrix.copy()
    for gene in z.columns:
        sd = z[gene].std(ddof=0)
        z[gene] = (z[gene] - z[gene].mean()) / sd if sd and not pd.isna(sd) else 0.0

    out = matrix.reset_index()[["sample", "group"]].copy()
    for module, module_genes in MODULES.items():
        available = [gene for gene in module_genes if gene in z.columns and z[gene].notna().any()]
        out[f"{module}_score"] = z[available].mean(axis=1).to_numpy() if available else np.nan
        out[f"{module}_raw_log1p_cpm_mean"] = matrix[available].mean(axis=1).to_numpy() if available else np.nan
        out[f"{module}_n_genes"] = len(available)
    out["AP_TLS_composite_score"] = out[["HLA_II_AP_score", "B_TLS_score"]].mean(axis=1)
    out["AP_TLS_composite_raw_log1p_cpm_mean"] = out[["HLA_II_AP_raw_log1p_cpm_mean", "B_TLS_raw_log1p_cpm_mean"]].mean(axis=1)
    out["AP_TLS_composite_n_genes"] = out["HLA_II_AP_n_genes"] + out["B_TLS_n_genes"]
    return out


def build_spot_burden(spot_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    thresholds = (
        spot_df[spot_df["group"].eq("Control")]
        .groupby("module")["spot_log1p_cp10k"]
        .quantile(0.90)
        .rename("control_p90")
        .reset_index()
    )
    merged = spot_df.merge(thresholds, on="module", how="left")
    for (sample, group, module), d in merged.groupby(["sample", "group", "module"], sort=False):
        rows.append(
            {
                "sample": sample,
                "group": group,
                "module": module,
                "n_spots": int(d.shape[0]),
                "mean_spot_log1p_cp10k": float(d["spot_log1p_cp10k"].mean()),
                "median_spot_log1p_cp10k": float(d["spot_log1p_cp10k"].median()),
                "positive_fraction": float(d["module_positive"].mean()),
                "control_p90": float(d["control_p90"].iloc[0]) if d["control_p90"].notna().any() else np.nan,
                "fraction_above_control_p90": float((d["spot_log1p_cp10k"] > d["control_p90"]).mean()) if d["control_p90"].notna().any() else np.nan,
                "n_module_genes_available": int(d["n_module_genes_available"].max()),
            }
        )
    return pd.DataFrame(rows)


def cohen_d(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna().to_numpy()
    b = pd.to_numeric(b, errors="coerce").dropna().to_numpy()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = ((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) / (len(a) + len(b) - 2)
    return float((np.mean(a) - np.mean(b)) / math.sqrt(pooled)) if pooled > 0 else np.nan


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray([1.0 if pd.isna(x) else float(x) for x in pvals])
    order = np.argsort(p)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(p) + 1)
    q = p * len(p) / ranks
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    out = np.empty_like(q_sorted)
    out[order] = np.minimum(q_sorted, 1.0)
    return out.tolist()


def exact_permutation(values: np.ndarray, labels: np.ndarray) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels).astype(bool)
    ok = np.isfinite(values)
    values = values[ok]
    labels = labels[ok]
    n_pos = int(labels.sum())
    obs = float(values[labels].mean() - values[~labels].mean())
    extreme = 0
    total = 0
    for pos_idx in combinations(range(len(values)), n_pos):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(pos_idx)] = True
        diff = float(values[mask].mean() - values[~mask].mean())
        if abs(diff) >= abs(obs) - 1e-12:
            extreme += 1
        total += 1
    return obs, (extreme + 1) / (total + 1), total


def run_sample_contrasts(sample_scores: pd.DataFrame) -> pd.DataFrame:
    contrasts = [
        ("AITD_vs_Control", ["HT", "GD"], ["Control"]),
        ("HT_vs_Control", ["HT"], ["Control"]),
        ("GD_vs_Control", ["GD"], ["Control"]),
        ("HT_vs_GD", ["HT"], ["GD"]),
    ]
    rows = []
    for contrast_name, group_a_values, group_b_values in contrasts:
        subset = sample_scores[sample_scores["group"].isin(group_a_values + group_b_values)].copy()
        y = subset["group"].isin(group_a_values).to_numpy()
        for module in PRIMARY_MODULES:
            col = f"{module}_score"
            a = subset.loc[subset["group"].isin(group_a_values), col]
            b = subset.loc[subset["group"].isin(group_b_values), col]
            delta, p_exact, n_perm = exact_permutation(subset[col].to_numpy(), y)
            try:
                mwu = stats.mannwhitneyu(a, b, alternative="two-sided", method="auto")
                mwu_p = float(mwu.pvalue)
            except ValueError:
                mwu_p = np.nan
            try:
                auc = float(roc_auc_score(y.astype(int), subset[col].to_numpy()))
            except ValueError:
                auc = np.nan
            rows.append(
                {
                    "external_dataset": "GSE248205",
                    "contrast": contrast_name,
                    "group_a": "+".join(group_a_values),
                    "group_b": "+".join(group_b_values),
                    "module": module,
                    "n_group_a": int(a.notna().sum()),
                    "n_group_b": int(b.notna().sum()),
                    "mean_group_a": float(a.mean()),
                    "mean_group_b": float(b.mean()),
                    "delta_mean_score_group_a_minus_group_b": float(delta),
                    "cohen_d_group_a_minus_group_b": cohen_d(a, b),
                    "exact_two_sided_p": float(p_exact),
                    "n_exact_permutations": int(n_perm),
                    "mannwhitney_p": mwu_p,
                    "auc_group_a_vs_group_b": auc,
                    "boundary": "AITD spatial expression mechanism, not cancer HLA allele genotype",
                }
            )
    out = pd.DataFrame(rows)
    out["BH_FDR_exact_all_tests"] = bh_fdr(out["exact_two_sided_p"].tolist())
    return out


def sample_score_long(sample_scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for module in MODULE_ORDER:
        col = f"{module}_score"
        if col not in sample_scores:
            continue
        tmp = sample_scores[["sample", "group", col]].copy()
        tmp["module"] = module
        tmp["score"] = tmp[col]
        rows.append(tmp[["sample", "group", "module", "score"]])
    return pd.concat(rows, ignore_index=True)


def plot_figures(sample_scores: pd.DataFrame, score_long: pd.DataFrame, burden_df: pd.DataFrame, contrast_df: pd.DataFrame, spot_df: pd.DataFrame) -> None:
    plot_modules = ["HLA_I", "HLA_II_AP", "B_TLS", "CD74_MIF_axis", "AP_TLS_composite"]
    fig, ax = plt.subplots(figsize=(10.3, 4.9))
    plot_df = score_long[score_long["module"].isin(plot_modules)].copy()
    sns.pointplot(
        data=plot_df,
        x="module",
        y="score",
        hue="group",
        order=plot_modules,
        hue_order=GROUP_ORDER,
        palette=PALETTE,
        dodge=0.5,
        errorbar=None,
        markers="o",
        linestyles="",
        ax=ax,
    )
    sns.stripplot(
        data=plot_df,
        x="module",
        y="score",
        hue="group",
        order=plot_modules,
        hue_order=GROUP_ORDER,
        palette=PALETTE,
        dodge=True,
        jitter=0.12,
        linewidth=0.4,
        edgecolor="black",
        alpha=0.85,
        ax=ax,
    )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[: len(GROUP_ORDER)], labels[: len(GROUP_ORDER)], title="", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xlabel("")
    ax.set_ylabel("Sample module score (mean z of pseudobulk logCPM)")
    ax.set_title("GSE248205 AITD spatial pseudobulk module scores")
    fig.tight_layout()
    fig.savefig(FIGS / "F01_gse248205_sample_module_scores.png", dpi=220)
    plt.close(fig)

    primary = contrast_df[contrast_df["contrast"].isin(["AITD_vs_Control", "HT_vs_Control", "GD_vs_Control"])].copy()
    primary = primary[primary["module"].isin(["HLA_II_AP", "B_TLS", "CD74_MIF_axis", "AP_TLS_composite"])]
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    sns.barplot(
        data=primary,
        x="contrast",
        y="cohen_d_group_a_minus_group_b",
        hue="module",
        hue_order=["HLA_II_AP", "B_TLS", "CD74_MIF_axis", "AP_TLS_composite"],
        palette=["#E45756", "#54A24B", "#B279A2", "#F58518"],
        edgecolor="black",
        linewidth=0.5,
        ax=ax,
    )
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("")
    ax.set_ylabel("Cohen's d")
    ax.set_title("Sample-level exact-permutation contrasts")
    ax.tick_params(axis="x", rotation=20)
    ax.legend(title="", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(FIGS / "F02_gse248205_sample_effect_sizes.png", dpi=220)
    plt.close(fig)

    burden_plot = burden_df[burden_df["module"].isin(plot_modules)].copy()
    burden_plot["sample"] = pd.Categorical(burden_plot["sample"], SAMPLE_ORDER, ordered=True)
    heat = burden_plot.pivot(index="sample", columns="module", values="fraction_above_control_p90").loc[SAMPLE_ORDER, plot_modules]
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    sns.heatmap(heat, cmap="rocket_r", vmin=0, vmax=max(0.5, float(np.nanmax(heat.to_numpy()))), linewidths=0.4, linecolor="white", cbar_kws={"label": "Fraction above control p90"}, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Spatial burden above control 90th percentile")
    fig.tight_layout()
    fig.savefig(FIGS / "F03_gse248205_spot_burden_heatmap.png", dpi=220)
    plt.close(fig)

    spatial_samples = ["C1", "HT1", "GD1"]
    spatial = spot_df[spot_df["sample"].isin(spatial_samples) & spot_df["module"].eq("HLA_II_AP")].copy()
    vmax = float(spatial["spot_log1p_cp10k"].quantile(0.99))
    fig, axes = plt.subplots(1, len(spatial_samples), figsize=(11.0, 3.8), sharex=False, sharey=False)
    for ax, sample in zip(axes, spatial_samples, strict=True):
        d = spatial[spatial["sample"].eq(sample)]
        sc = ax.scatter(
            d["pxl_col_in_fullres"],
            -d["pxl_row_in_fullres"],
            c=d["spot_log1p_cp10k"],
            s=5,
            cmap="viridis",
            vmin=0,
            vmax=vmax,
            linewidths=0,
        )
        ax.set_title(sample)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
    cbar = fig.colorbar(sc, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
    cbar.set_label("HLA-II/AP spot score")
    fig.suptitle("Spatial localization of HLA-II/AP module", y=0.98)
    fig.savefig(FIGS / "F04_gse248205_hla2_spatial_snapshots.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def write_report(manifest_df: pd.DataFrame, contrast_df: pd.DataFrame, burden_df: pd.DataFrame) -> None:
    primary = contrast_df[
        contrast_df["contrast"].isin(["AITD_vs_Control", "HT_vs_Control", "GD_vs_Control"])
        & contrast_df["module"].isin(["HLA_II_AP", "B_TLS", "CD74_MIF_axis", "AP_TLS_composite"])
    ].sort_values(["contrast", "module"])
    primary_rows = []
    for r in primary.itertuples(index=False):
        primary_rows.append(
            f"| {r.contrast} | {r.module} | {r.cohen_d_group_a_minus_group_b:.2f} | "
            f"{r.delta_mean_score_group_a_minus_group_b:.2f} | {r.exact_two_sided_p:.3f} | "
            f"{r.auc_group_a_vs_group_b:.2f} |"
        )

    burden = burden_df[burden_df["module"].isin(["HLA_II_AP", "B_TLS", "CD74_MIF_axis", "AP_TLS_composite"])]
    burden_summary = (
        burden.groupby(["group", "module"], as_index=False)
        .agg(
            mean_fraction_above_control_p90=("fraction_above_control_p90", "mean"),
            mean_spot_score=("mean_spot_log1p_cp10k", "mean"),
        )
        .sort_values(["module", "group"])
    )
    burden_rows = []
    for r in burden_summary.itertuples(index=False):
        burden_rows.append(
            f"| {r.group} | {r.module} | {r.mean_fraction_above_control_p90:.3f} | {r.mean_spot_score:.2f} |"
        )

    n_spots = int(manifest_df["n_spots_in_tissue"].sum())
    group_counts = manifest_df["group"].value_counts().reindex(GROUP_ORDER).dropna().astype(int)
    count_text = ", ".join(f"{k} n={v}" for k, v in group_counts.items())
    report = f"""# GSE248205 AITD spatial HLA/AP 검증

## 결론

GSE248205는 2024 Nature Communications AITD spatial transcriptomics 자료로, 이번 재분석에서는 {count_text}, 총 {n_spots:,} in-tissue spots를 사용했다.

Paper 4에는 이 데이터를 **HLA allele genetics의 조직 기전 보강**으로 붙일 수 있다. 즉, 유전좌위 주장은 Paper 4의 HLA/GWAS 자료가 담당하고, GSE248205는 HT/GD 조직에서 HLA-II antigen-presentation, CD74/MIF, TLS/B-cell 축이 공간적으로 확장되는지를 보여준다.

Paper 2에는 암 allele 주장이 아니라 `HT background antigen-presentation ecosystem`을 설명하는 외부 tissue context로만 제한해서 넣어야 한다.

## Sample-level contrasts

| contrast | module | Cohen's d | delta score | exact p | AUC |
|---|---|---:|---:|---:|---:|
{chr(10).join(primary_rows)}

## Spatial burden

`fraction_above_control_p90`는 control spot 분포의 90th percentile보다 높은 spot 비율이다.

| group | module | fraction above control p90 | mean spot score |
|---|---|---:|---:|
{chr(10).join(burden_rows)}

## 논문 반영 포인트

1. Paper 4 Discussion/Mechanism figure에 `spatial AITD tissue validation` 패널을 추가한다.
2. GD와 HT를 분리해 보여주면 class-II/AP 공통축과 disease-specific tissue architecture를 동시에 주장할 수 있다.
3. 이 자료는 genotype 자료가 아니므로, allele association replication이라고 쓰면 안 된다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse248205_aitd_spatial_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse248205_aitd_spatial_validation/figures/`
- Source data: `project/data/external/GSE248205/`
"""
    REPORT.write_text(report, encoding="utf-8")
    (OUT / "GSE248205_AITD_SPATIAL_VALIDATION_KR.md").write_text(report, encoding="utf-8")


def main() -> None:
    ensure_dirs()
    genes = target_genes()
    manifests = []
    spots = []
    genes_long = []
    for sample in SAMPLE_ORDER:
        manifest, spot_df, gene_df = read_sample(sample, genes)
        manifests.append(manifest)
        spots.append(spot_df)
        genes_long.append(gene_df)

    manifest_df = pd.concat(manifests, ignore_index=True)
    spot_df = pd.concat(spots, ignore_index=True)
    gene_df = pd.concat(genes_long, ignore_index=True)
    burden_df = build_spot_burden(spot_df)
    sample_scores = score_sample_modules(gene_df).merge(manifest_df, on=["sample", "group"], how="left")
    score_long = sample_score_long(sample_scores)
    contrast_df = run_sample_contrasts(sample_scores)

    manifest_df.to_csv(TABLES / "T01_gse248205_sample_manifest.tsv", sep="\t", index=False)
    spot_df.to_csv(TABLES / "T02_gse248205_spot_module_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    gene_df.to_csv(TABLES / "T03_gse248205_target_gene_pseudobulk.tsv", sep="\t", index=False)
    burden_df.to_csv(TABLES / "T04_gse248205_spot_burden_by_sample.tsv", sep="\t", index=False)
    sample_scores.to_csv(TABLES / "T05_gse248205_sample_module_scores.tsv", sep="\t", index=False)
    score_long.to_csv(TABLES / "T06_gse248205_sample_module_scores_long.tsv", sep="\t", index=False)
    contrast_df.to_csv(TABLES / "T07_gse248205_sample_exact_contrasts.tsv", sep="\t", index=False)

    plot_figures(sample_scores, score_long, burden_df, contrast_df, spot_df)
    write_report(manifest_df, contrast_df, burden_df)

    manifest = {
        "dataset": "GSE248205",
        "source_dir": str(DATA.relative_to(ROOT)),
        "n_samples": int(manifest_df.shape[0]),
        "groups": manifest_df["group"].value_counts().to_dict(),
        "n_spots_in_tissue": int(manifest_df["n_spots_in_tissue"].sum()),
        "boundary": "AITD spatial expression mechanism only; no cancer HLA allele/genotype claims",
        "tables": sorted(str(p.relative_to(ROOT)) for p in TABLES.glob("*.tsv*")),
        "figures": sorted(str(p.relative_to(ROOT)) for p in FIGS.glob("*.png")),
        "report": str(REPORT.relative_to(ROOT)),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
