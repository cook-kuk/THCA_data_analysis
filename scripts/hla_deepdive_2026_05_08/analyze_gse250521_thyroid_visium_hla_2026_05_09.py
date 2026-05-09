#!/usr/bin/env python3
"""Analyze GSE250521 thyroid Visium data for HLA/AP/TLS spatial ecology.

Paper 2 use: thyroid cancer spatial-expression generalization. This is not
HT-specific and must not be framed as HLA allele/genotype evidence.
"""

from __future__ import annotations

import gzip
import io
import json
import math
import re
import tarfile
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
DATA = ROOT / "project/data/external/GSE250521"
RAW_TAR = DATA / "GSE250521_RAW.tar"
BASE = ROOT / "project/results/hla_two_paper_synthesis_2026_05_09"
OUT = BASE / "gse250521_thyroid_visium_hla_validation"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
REPORT = ROOT / "project/reports/2026_05_09_GSE250521_THYROID_VISIUM_HLA_VALIDATION_KR.md"


MODULES = {
    "HLA_I": ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "NLRC5"],
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
    "B_TLS": ["MS4A1", "CD79A", "CD79B", "CD19", "CD22", "CD27", "CD38", "MZB1", "JCHAIN", "IGHG1", "IGHG2", "IGKC", "CXCL13", "CCL19", "CCL21", "LTB", "LTBR", "CCR7", "POU2AF1", "AICDA"],
    "T_IFNG": ["CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "IFNG", "GZMB", "PRF1", "CXCL9", "CXCL10", "CXCL11", "STAT1", "IRF1"],
    "Myeloid_DC": ["LST1", "LYZ", "FCER1A", "CLEC10A", "CD1C", "LAMP3", "ITGAX", "CD68", "FCGR3A", "C1QA", "C1QB"],
    "Thyrocyte": ["TG", "TPO", "EPCAM", "KRT8", "KRT18", "SLC5A5", "TSHR", "FOXE1", "NKX2-1"],
    "CD74_MIF_axis": ["CD74", "MIF"],
}

MODULE_ORDER = ["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "Myeloid_DC", "Thyrocyte", "CD74_MIF_axis", "AP_TLS_composite"]
PRIMARY_MODULES = ["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "Myeloid_DC", "CD74_MIF_axis", "AP_TLS_composite"]
STAGE_ORDER = ["N", "PTC", "LPTC", "ATC"]
STAGE_RANK = {"N": 0, "PTC": 1, "LPTC": 2, "ATC": 3}
PALETTE = {"N": "#72B7B2", "PTC": "#F58518", "LPTC": "#B279A2", "ATC": "#E45756"}


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)


def target_genes() -> set[str]:
    genes: set[str] = set()
    for module_genes in MODULES.values():
        genes.update(module_genes)
    return genes


def parse_visium_sample(member: str) -> dict[str, str] | None:
    match = re.match(r"(GSM\d+)_([^_]+)_visium_matrix\.mtx\.gz$", member)
    if not match:
        return None
    gsm, label = match.groups()
    stage = label.split("-", 1)[0]
    return {
        "sample": f"{gsm}_{label}",
        "gsm": gsm,
        "label": label,
        "stage": stage,
        "stage_rank": STAGE_RANK[stage],
        "prefix": f"{gsm}_{label}_visium",
    }


def read_gzip_member_text(tar: tarfile.TarFile, name: str) -> io.TextIOWrapper:
    fileobj = tar.extractfile(name)
    if fileobj is None:
        raise FileNotFoundError(name)
    return io.TextIOWrapper(gzip.GzipFile(fileobj=fileobj, mode="rb"))


def read_gzip_member_binary(tar: tarfile.TarFile, name: str) -> gzip.GzipFile:
    fileobj = tar.extractfile(name)
    if fileobj is None:
        raise FileNotFoundError(name)
    return gzip.GzipFile(fileobj=fileobj, mode="rb")


def read_matrix_member(tar: tarfile.TarFile, name: str):
    with read_gzip_member_binary(tar, name) as gz:
        try:
            return mmread(gz).tocsr()
        except Exception:
            gz.seek(0)
            return mmread(io.BytesIO(gz.read())).tocsr()


def read_one_visium_sample(tar: tarfile.TarFile, sample: dict[str, str], genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prefix = sample["prefix"]
    features = pd.read_csv(
        read_gzip_member_text(tar, f"{prefix}_features.tsv.gz"),
        sep="\t",
        header=None,
        names=["gene_id", "gene", "feature_type"],
        dtype=str,
    )
    barcodes = pd.read_csv(
        read_gzip_member_text(tar, f"{prefix}_barcodes.tsv.gz"),
        sep="\t",
        header=None,
        names=["barcode"],
        dtype=str,
    )["barcode"].tolist()
    positions = pd.read_csv(
        read_gzip_member_text(tar, f"{prefix}_tissue_positions_list.csv.gz"),
        header=None,
        names=["barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"],
    )
    matrix = read_matrix_member(tar, f"{prefix}_matrix.mtx.gz")
    barcode_to_idx = {barcode: idx for idx, barcode in enumerate(barcodes)}
    positions = positions[positions["barcode"].isin(barcode_to_idx)].copy()
    positions["barcode_idx"] = positions["barcode"].map(barcode_to_idx).astype(int)

    total_counts_all = np.asarray(matrix.sum(axis=0)).ravel()
    n_genes_all = np.asarray((matrix > 0).sum(axis=0)).ravel()
    keep = positions["in_tissue"].eq(1).to_numpy()
    keep &= n_genes_all[positions["barcode_idx"].to_numpy()] >= 200
    keep &= total_counts_all[positions["barcode_idx"].to_numpy()] > 0
    positions = positions.loc[keep].copy().reset_index(drop=True)
    keep_cols = positions["barcode_idx"].to_numpy()
    matrix = matrix[:, keep_cols]

    total_counts = np.asarray(matrix.sum(axis=0)).ravel()
    n_genes_by_counts = np.asarray((matrix > 0).sum(axis=0)).ravel()
    gene_to_indices: dict[str, list[int]] = {}
    for idx, gene in enumerate(features["gene"].astype(str)):
        if gene in genes:
            gene_to_indices.setdefault(gene, []).append(idx)

    target_gene_counts = {}
    target_gene_expr = {}
    gene_rows = []
    for gene in sorted(genes):
        indices = gene_to_indices.get(gene, [])
        if indices:
            counts = np.asarray(matrix[indices, :].sum(axis=0)).ravel()
            expr = np.log1p(counts / total_counts * 1e4)
            pseudobulk_count = float(counts.sum())
            pseudobulk_expr = math.log1p(pseudobulk_count / total_counts.sum() * 1e6)
        else:
            counts = np.zeros(matrix.shape[1], dtype=float)
            expr = np.full(matrix.shape[1], np.nan, dtype=float)
            pseudobulk_count = 0.0
            pseudobulk_expr = np.nan
        target_gene_counts[gene] = counts
        target_gene_expr[gene] = expr
        gene_rows.append(
            {
                "sample": sample["sample"],
                "gsm": sample["gsm"],
                "label": sample["label"],
                "stage": sample["stage"],
                "stage_rank": sample["stage_rank"],
                "gene": gene,
                "pseudobulk_count": pseudobulk_count,
                "log1p_cpm": pseudobulk_expr,
                "n_feature_rows": len(indices),
            }
        )

    expr_df = pd.DataFrame(target_gene_expr)
    z_df = expr_df.copy()
    for gene in z_df.columns:
        sd = z_df[gene].std(ddof=0)
        z_df[gene] = (z_df[gene] - z_df[gene].mean()) / sd if sd and not pd.isna(sd) else 0.0

    spot_base = positions[["barcode", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]].copy()
    spot_base.insert(0, "sample", sample["sample"])
    spot_base.insert(1, "gsm", sample["gsm"])
    spot_base.insert(2, "label", sample["label"])
    spot_base.insert(3, "stage", sample["stage"])
    spot_base.insert(4, "stage_rank", sample["stage_rank"])
    spot_base["total_counts"] = total_counts
    spot_base["n_genes_by_counts"] = n_genes_by_counts

    spot_rows = []
    for module, module_genes in MODULES.items():
        available = [gene for gene in module_genes if gene in gene_to_indices]
        tmp = spot_base.copy()
        tmp["module"] = module
        tmp["n_module_genes_available"] = len(available)
        if available:
            module_counts = np.vstack([target_gene_counts[gene] for gene in available]).sum(axis=0)
            tmp["module_counts"] = module_counts
            tmp["spot_raw_mean_log1p_cp10k"] = expr_df[available].mean(axis=1).to_numpy()
            tmp["spot_module_score_within_sample"] = z_df[available].mean(axis=1).to_numpy()
            tmp["module_positive"] = module_counts > 0
        else:
            tmp["module_counts"] = 0.0
            tmp["spot_raw_mean_log1p_cp10k"] = np.nan
            tmp["spot_module_score_within_sample"] = np.nan
            tmp["module_positive"] = False
        spot_rows.append(tmp)

    composite = spot_base.copy()
    composite["module"] = "AP_TLS_composite"
    composite["n_module_genes_available"] = sum(
        len([gene for gene in MODULES[module] if gene in gene_to_indices]) for module in ["HLA_II_AP", "B_TLS"]
    )
    hla2 = spot_rows[list(MODULES).index("HLA_II_AP")]
    btls = spot_rows[list(MODULES).index("B_TLS")]
    composite["module_counts"] = hla2["module_counts"].to_numpy() + btls["module_counts"].to_numpy()
    composite["spot_raw_mean_log1p_cp10k"] = np.vstack([hla2["spot_raw_mean_log1p_cp10k"], btls["spot_raw_mean_log1p_cp10k"]]).mean(axis=0)
    composite["spot_module_score_within_sample"] = np.vstack([hla2["spot_module_score_within_sample"], btls["spot_module_score_within_sample"]]).mean(axis=0)
    composite["module_positive"] = composite["module_counts"].gt(0)
    spot_rows.append(composite)

    spot_df = pd.concat(spot_rows, ignore_index=True)
    gene_df = pd.DataFrame(gene_rows)
    manifest = pd.DataFrame(
        [
            {
                "sample": sample["sample"],
                "gsm": sample["gsm"],
                "label": sample["label"],
                "stage": sample["stage"],
                "stage_rank": sample["stage_rank"],
                "n_spots_post_qc": int(matrix.shape[1]),
                "n_features": int(matrix.shape[0]),
                "total_umis_post_qc": float(total_counts.sum()),
                "median_umis_per_spot": float(np.median(total_counts)),
                "median_genes_per_spot": float(np.median(n_genes_by_counts)),
                "n_target_genes_mapped": int(sum(1 for gene in genes if gene in gene_to_indices)),
            }
        ]
    )
    target_map = pd.DataFrame(
        [
            {
                "sample": sample["sample"],
                "gene": gene,
                "n_feature_rows": len(indices),
                "feature_ids": ",".join(features.loc[indices, "gene_id"].astype(str).tolist()) if indices else "",
            }
            for gene, indices in sorted(gene_to_indices.items())
        ]
    )
    return manifest, spot_df, gene_df, target_map


def load_visium_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genes = target_genes()
    with tarfile.open(RAW_TAR) as tar:
        members = tar.getnames()
        samples = [parse_visium_sample(member) for member in members]
        samples = [sample for sample in samples if sample is not None]
        samples = sorted(samples, key=lambda x: (x["stage_rank"], x["label"]))
        manifests = []
        spot_frames = []
        gene_frames = []
        map_frames = []
        for sample in samples:
            manifest, spot_df, gene_df, target_map = read_one_visium_sample(tar, sample, genes)
            manifests.append(manifest)
            spot_frames.append(spot_df)
            gene_frames.append(gene_df)
            map_frames.append(target_map)
    return (
        pd.concat(manifests, ignore_index=True),
        pd.concat(spot_frames, ignore_index=True),
        pd.concat(gene_frames, ignore_index=True),
        pd.concat(map_frames, ignore_index=True),
    )


def score_sample_modules(gene_df: pd.DataFrame) -> pd.DataFrame:
    matrix = gene_df.pivot(index=["sample", "gsm", "label", "stage", "stage_rank"], columns="gene", values="log1p_cpm")
    z = matrix.copy()
    for gene in z.columns:
        sd = z[gene].std(ddof=0)
        z[gene] = (z[gene] - z[gene].mean()) / sd if sd and not pd.isna(sd) else 0.0
    out = matrix.reset_index()[["sample", "gsm", "label", "stage", "stage_rank"]].copy()
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
    thresholds = (
        spot_df[spot_df["stage"].eq("N")]
        .groupby("module")["spot_raw_mean_log1p_cp10k"]
        .quantile(0.90)
        .rename("N_p90_raw_spot_score")
        .reset_index()
    )
    merged = spot_df.merge(thresholds, on="module", how="left")
    rows = []
    for (sample, label, stage, stage_rank, module), d in merged.groupby(["sample", "label", "stage", "stage_rank", "module"], sort=False):
        rows.append(
            {
                "sample": sample,
                "label": label,
                "stage": stage,
                "stage_rank": int(stage_rank),
                "module": module,
                "n_spots": int(d.shape[0]),
                "mean_raw_spot_score": float(d["spot_raw_mean_log1p_cp10k"].mean()),
                "median_raw_spot_score": float(d["spot_raw_mean_log1p_cp10k"].median()),
                "mean_within_sample_z_score": float(d["spot_module_score_within_sample"].mean()),
                "positive_fraction": float(d["module_positive"].mean()),
                "N_p90_raw_spot_score": float(d["N_p90_raw_spot_score"].iloc[0]) if d["N_p90_raw_spot_score"].notna().any() else np.nan,
                "fraction_above_N_p90": float((d["spot_raw_mean_log1p_cp10k"] > d["N_p90_raw_spot_score"]).mean())
                if d["N_p90_raw_spot_score"].notna().any()
                else np.nan,
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


def exact_permutation(values: np.ndarray, labels: np.ndarray) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels).astype(bool)
    ok = np.isfinite(values)
    values = values[ok]
    labels = labels[ok]
    n_pos = int(labels.sum())
    obs = float(values[labels].mean() - values[~labels].mean())
    extreme = total = 0
    for pos_idx in combinations(range(len(values)), n_pos):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(pos_idx)] = True
        diff = float(values[mask].mean() - values[~mask].mean())
        if abs(diff) >= abs(obs) - 1e-12:
            extreme += 1
        total += 1
    return obs, (extreme + 1) / (total + 1), total


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


def run_sample_contrasts(scores: pd.DataFrame) -> pd.DataFrame:
    contrast_defs = [
        ("PTC_vs_N", ["PTC"], ["N"]),
        ("LPTC_vs_N", ["LPTC"], ["N"]),
        ("ATC_vs_N", ["ATC"], ["N"]),
        ("Cancer_vs_N", ["PTC", "LPTC", "ATC"], ["N"]),
        ("LPTC_ATC_vs_PTC", ["LPTC", "ATC"], ["PTC"]),
        ("ATC_vs_PTC", ["ATC"], ["PTC"]),
    ]
    rows = []
    for cname, group_a, group_b in contrast_defs:
        subset = scores[scores["stage"].isin(group_a + group_b)].copy()
        y = subset["stage"].isin(group_a).to_numpy()
        for module in PRIMARY_MODULES:
            col = f"{module}_score"
            a = subset.loc[subset["stage"].isin(group_a), col]
            b = subset.loc[subset["stage"].isin(group_b), col]
            delta, p_exact, n_perm = exact_permutation(subset[col].to_numpy(), y)
            try:
                auc = float(roc_auc_score(y.astype(int), subset[col].to_numpy()))
            except ValueError:
                auc = np.nan
            rows.append(
                {
                    "external_dataset": "GSE250521",
                    "contrast": cname,
                    "group_a": "+".join(group_a),
                    "group_b": "+".join(group_b),
                    "module": module,
                    "n_group_a": int(a.notna().sum()),
                    "n_group_b": int(b.notna().sum()),
                    "mean_group_a": float(a.mean()),
                    "mean_group_b": float(b.mean()),
                    "delta_mean_score_group_a_minus_group_b": delta,
                    "cohen_d_group_a_minus_group_b": cohen_d(a, b),
                    "exact_two_sided_p": p_exact,
                    "n_exact_permutations": n_perm,
                    "mannwhitney_p": float(stats.mannwhitneyu(a, b, alternative="two-sided", method="auto").pvalue),
                    "auc_group_a_vs_group_b": auc,
                    "boundary": "Thyroid cancer spatial expression ecology only; no HLA allele/genotype claims",
                }
            )
    out = pd.DataFrame(rows)
    out["BH_FDR_exact_all_tests"] = bh_fdr(out["exact_two_sided_p"].tolist())
    return out


def run_stage_trends(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for module in PRIMARY_MODULES:
        col = f"{module}_score"
        rho, p_spearman = stats.spearmanr(scores["stage_rank"], scores[col], nan_policy="omit")
        slope, intercept, rvalue, p_lin, stderr = stats.linregress(scores["stage_rank"], scores[col])
        rows.append(
            {
                "external_dataset": "GSE250521",
                "module": module,
                "n_samples": int(scores[col].notna().sum()),
                "spearman_rho_stage_rank": float(rho),
                "spearman_p": float(p_spearman),
                "linear_slope_per_stage": float(slope),
                "linear_r": float(rvalue),
                "linear_p": float(p_lin),
                "boundary": "Thyroid cancer spatial expression stage trend only; no HLA allele/genotype claims",
            }
        )
    out = pd.DataFrame(rows)
    out["BH_FDR_spearman_all_modules"] = bh_fdr(out["spearman_p"].tolist())
    return out


def make_sample_long(scores: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for module in MODULE_ORDER:
        col = f"{module}_score"
        tmp = scores[["sample", "label", "stage", "stage_rank", col]].copy()
        tmp["module"] = module
        tmp["score"] = tmp[col]
        frames.append(tmp[["sample", "label", "stage", "stage_rank", "module", "score"]])
    return pd.concat(frames, ignore_index=True)


def plot_figures(scores: pd.DataFrame, long_scores: pd.DataFrame, contrasts: pd.DataFrame, trends: pd.DataFrame, burden: pd.DataFrame, spot_df: pd.DataFrame) -> None:
    plot_modules = ["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "AP_TLS_composite"]
    plot_df = long_scores[long_scores["module"].isin(plot_modules)].copy()
    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    sns.pointplot(data=plot_df, x="module", y="score", hue="stage", order=plot_modules, hue_order=STAGE_ORDER, palette=PALETTE, dodge=0.45, errorbar=None, ax=ax)
    sns.stripplot(data=plot_df, x="module", y="score", hue="stage", order=plot_modules, hue_order=STAGE_ORDER, palette=PALETTE, dodge=True, jitter=0.10, linewidth=0.4, edgecolor="black", alpha=0.85, ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[: len(STAGE_ORDER)], labels[: len(STAGE_ORDER)], title="", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xlabel("")
    ax.set_ylabel("Sample module score (mean z of pseudobulk logCPM)")
    ax.set_title("GSE250521 thyroid Visium pseudobulk HLA/AP ecology")
    fig.tight_layout()
    fig.savefig(FIGS / "F01_gse250521_stage_module_scores.png", dpi=220)
    plt.close(fig)

    trend_plot = trends.sort_values("spearman_rho_stage_rank", ascending=True)
    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    colors = ["#E45756" if x > 0 else "#72B7B2" for x in trend_plot["spearman_rho_stage_rank"]]
    ax.barh(trend_plot["module"], trend_plot["spearman_rho_stage_rank"], color=colors, edgecolor="black", linewidth=0.6)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Spearman rho vs stage rank")
    ax.set_ylabel("")
    ax.set_title("Stage trend across N/PTC/LPTC/ATC")
    for i, row in enumerate(trend_plot.itertuples(index=False)):
        val = row.spearman_rho_stage_rank
        ax.text(val + (0.03 if val >= 0 else -0.03), i, f"q={row.BH_FDR_spearman_all_modules:.2g}", va="center", ha="left" if val >= 0 else "right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(FIGS / "F02_gse250521_stage_trend_rho.png", dpi=220)
    plt.close(fig)

    burden_plot = burden[burden["module"].isin(plot_modules)].copy()
    label_order = sorted(burden_plot["label"].dropna().unique(), key=lambda x: (STAGE_RANK[x.split("-", 1)[0]], x))
    burden_plot["label"] = pd.Categorical(burden_plot["label"], label_order, ordered=True)
    heat = burden_plot.pivot(index="label", columns="module", values="fraction_above_N_p90").loc[:, plot_modules]
    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    sns.heatmap(heat, cmap="rocket_r", vmin=0, vmax=max(0.5, float(np.nanmax(heat.to_numpy()))), linewidths=0.4, linecolor="white", cbar_kws={"label": "Fraction above N p90"}, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Spot burden above normal-thyroid 90th percentile")
    fig.tight_layout()
    fig.savefig(FIGS / "F03_gse250521_spatial_burden_heatmap.png", dpi=220)
    plt.close(fig)

    spatial_labels = ["N-1", "PTC-1", "LPTC-1", "ATC-1"]
    spatial = spot_df[spot_df["label"].isin(spatial_labels) & spot_df["module"].eq("HLA_II_AP")].copy()
    vmax = float(spatial["spot_raw_mean_log1p_cp10k"].quantile(0.99))
    fig, axes = plt.subplots(1, len(spatial_labels), figsize=(12.2, 3.6), sharex=False, sharey=False)
    for ax, label in zip(axes, spatial_labels, strict=True):
        d = spatial[spatial["label"].eq(label)]
        sc = ax.scatter(
            d["pxl_col_in_fullres"],
            -d["pxl_row_in_fullres"],
            c=d["spot_raw_mean_log1p_cp10k"],
            s=4,
            cmap="viridis",
            vmin=0,
            vmax=vmax,
            linewidths=0,
        )
        ax.set_title(label)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
    cbar = fig.colorbar(sc, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
    cbar.set_label("HLA-II/AP raw spot score")
    fig.suptitle("HLA-II/AP spatial snapshots", y=0.98)
    fig.savefig(FIGS / "F04_gse250521_hla2_spatial_snapshots.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    primary = contrasts[contrasts["contrast"].isin(["PTC_vs_N", "LPTC_vs_N", "ATC_vs_N"])].copy()
    primary = primary[primary["module"].isin(["HLA_I", "HLA_II_AP", "T_IFNG", "AP_TLS_composite"])]
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    sns.barplot(data=primary, x="contrast", y="cohen_d_group_a_minus_group_b", hue="module", edgecolor="black", linewidth=0.5, ax=ax)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("")
    ax.set_ylabel("Cohen's d vs N")
    ax.set_title("Stage-specific sample-level effect sizes")
    ax.tick_params(axis="x", rotation=20)
    ax.legend(title="", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(FIGS / "F05_gse250521_stage_effect_sizes.png", dpi=220)
    plt.close(fig)


def write_report(manifest: pd.DataFrame, contrasts: pd.DataFrame, trends: pd.DataFrame, burden: pd.DataFrame) -> None:
    stage_counts = manifest.groupby("stage").agg(n_slides=("sample", "nunique"), n_spots=("n_spots_post_qc", "sum")).reindex(STAGE_ORDER).reset_index()
    count_text = ", ".join(f"{r.stage} slides={int(r.n_slides)}, spots={int(r.n_spots):,}" for r in stage_counts.itertuples(index=False))
    trend_rows = []
    for r in trends.sort_values("spearman_rho_stage_rank", ascending=False).itertuples(index=False):
        trend_rows.append(f"| {r.module} | {r.spearman_rho_stage_rank:.2f} | {r.spearman_p:.3g} | {r.BH_FDR_spearman_all_modules:.3g} | {r.linear_slope_per_stage:.2f} |")
    primary = contrasts[contrasts["contrast"].isin(["PTC_vs_N", "LPTC_vs_N", "ATC_vs_N"])].copy()
    primary = primary[primary["module"].isin(["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "Myeloid_DC", "AP_TLS_composite"])]
    contrast_rows = []
    for r in primary.sort_values(["contrast", "cohen_d_group_a_minus_group_b"], ascending=[True, False]).itertuples(index=False):
        contrast_rows.append(
            f"| {r.contrast} | {r.module} | {r.cohen_d_group_a_minus_group_b:.2f} | {r.delta_mean_score_group_a_minus_group_b:.2f} | {r.exact_two_sided_p:.3f} | {r.auc_group_a_vs_group_b:.2f} |"
        )
    burden_summary = (
        burden[burden["module"].isin(["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "AP_TLS_composite"])]
        .groupby(["stage", "module"], as_index=False)
        .agg(mean_fraction_above_N_p90=("fraction_above_N_p90", "mean"), mean_raw_spot_score=("mean_raw_spot_score", "mean"))
    )
    burden_summary["stage"] = pd.Categorical(burden_summary["stage"], STAGE_ORDER, ordered=True)
    burden_rows = []
    for r in burden_summary.sort_values(["stage", "module"]).itertuples(index=False):
        burden_rows.append(f"| {r.stage} | {r.module} | {r.mean_fraction_above_N_p90:.3f} | {r.mean_raw_spot_score:.2f} |")
    report = f"""# GSE250521 thyroid Visium HLA/AP spatial 외부검증

## 결론

GSE250521 Visium raw matrix를 직접 다운로드해 16개 slide를 재분석했다. QC 후 구성은 {count_text} 이다.

이 데이터는 HT-specific이 아니다. Paper 2에서는 **thyroid cancer progression/spatial ecology generalization**으로만 사용하고, HT-overlap의 핵심 증거는 GSE138198/GSE163203에 둔다. HLA allele/genotype/risk claim은 금지한다.

## Stage trend across N/PTC/LPTC/ATC

| module | Spearman rho | p | FDR | linear slope |
|---|---:|---:|---:|---:|
{chr(10).join(trend_rows)}

## Stage-specific contrasts vs normal thyroid

| contrast | module | Cohen's d | delta score | exact p | AUC |
|---|---|---:|---:|---:|---:|
{chr(10).join(contrast_rows)}

## Spatial burden above N p90

| stage | module | fraction above N p90 | mean raw spot score |
|---|---|---:|---:|
{chr(10).join(burden_rows)}

## 논문 반영 포인트

1. Positive 방향이면 Paper 2 supplement external validation panel에 `GSE250521 Visium progression spatial generalization`으로 넣는다.
2. HT-specific 결론에는 직접 쓰지 않는다.
3. `HLA allele`, `genotype`, `risk allele` 문구를 쓰지 않는다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse250521_thyroid_visium_hla_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse250521_thyroid_visium_hla_validation/figures/`
- Source data: `project/data/external/GSE250521/GSE250521_RAW.tar`
"""
    REPORT.write_text(report, encoding="utf-8")
    (OUT / "GSE250521_THYROID_VISIUM_HLA_VALIDATION_KR.md").write_text(report, encoding="utf-8")


def main() -> None:
    ensure_dirs()
    manifest, spot_df, gene_df, target_map = load_visium_data()
    sample_scores = score_sample_modules(gene_df)
    long_scores = make_sample_long(sample_scores)
    burden = build_spot_burden(spot_df)
    contrasts = run_sample_contrasts(sample_scores)
    trends = run_stage_trends(sample_scores)

    target_map.to_csv(TABLES / "T01_gse250521_target_gene_feature_map.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": sorted(target_genes() - set(target_map["gene"])), "status": "not_mapped_in_visium_features"}).to_csv(
        TABLES / "T02_gse250521_target_genes_missing.tsv", sep="\t", index=False
    )
    manifest.to_csv(TABLES / "T03_gse250521_sample_manifest_qc.tsv", sep="\t", index=False)
    gene_df.to_csv(TABLES / "T04_gse250521_pseudobulk_target_gene_expression.tsv", sep="\t", index=False)
    sample_scores.to_csv(TABLES / "T05_gse250521_sample_module_scores.tsv", sep="\t", index=False)
    long_scores.to_csv(TABLES / "T06_gse250521_sample_module_scores_long.tsv", sep="\t", index=False)
    burden.to_csv(TABLES / "T07_gse250521_spatial_module_burden.tsv", sep="\t", index=False)
    contrasts.to_csv(TABLES / "T08_gse250521_sample_module_contrasts.tsv", sep="\t", index=False)
    trends.to_csv(TABLES / "T09_gse250521_stage_trends.tsv", sep="\t", index=False)
    spot_df.to_csv(TABLES / "T10_gse250521_spot_module_scores.tsv.gz", sep="\t", index=False, compression="gzip")

    plot_figures(sample_scores, long_scores, contrasts, trends, burden, spot_df)
    write_report(manifest, contrasts, trends, burden)

    manifest_json = {
        "dataset": "GSE250521",
        "n_visium_samples": int(manifest.shape[0]),
        "stage_counts": {stage: int(n) for stage, n in manifest.groupby("stage")["sample"].nunique().items()},
        "n_spots_post_qc": int(manifest["n_spots_post_qc"].sum()),
        "n_target_genes_mapped": int(target_map["gene"].nunique()),
        "boundary": "Thyroid cancer spatial expression ecology only; no HLA allele/genotype claims",
        "tables": sorted(str(p.relative_to(ROOT)) for p in TABLES.glob("*.tsv*")),
        "figures": sorted(str(p.relative_to(ROOT)) for p in FIGS.glob("*.png")),
        "report": str(REPORT.relative_to(ROOT)),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest_json, indent=2), encoding="utf-8")
    print(json.dumps(manifest_json, indent=2))


if __name__ == "__main__":
    main()
