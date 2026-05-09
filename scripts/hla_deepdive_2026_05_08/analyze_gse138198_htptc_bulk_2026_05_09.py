#!/usr/bin/env python3
"""Analyze GSE138198 as an external HT/PTC bulk-expression validation dataset.

Paper 2 boundary: this uses HLA/AP gene-expression modules only. It does not
make HLA allele/genotype claims in cancer cohorts.
"""

from __future__ import annotations

import csv
import gzip
import io
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
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "project/data/external/GSE138198"
SERIES = DATA / "GSE138198_series_matrix.txt.gz"
ANNOT = DATA / "GPL6244.annot.gz"
BASE = ROOT / "project/results/hla_two_paper_synthesis_2026_05_09"
OUT = BASE / "gse138198_htptc_bulk_validation"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
REPORT = ROOT / "project/reports/2026_05_09_GSE138198_HTPTC_BULK_VALIDATION_KR.md"


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
}

MODULE_ORDER = ["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "Myeloid_DC", "Thyrocyte", "AP_TLS_composite"]
PRIMARY_MODULES = ["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "Myeloid_DC", "AP_TLS_composite"]
GROUP_ORDER = ["TN", "HT", "PTCwithoutHT", "PTCwithHT", "mPTC"]
PALETTE = {
    "TN": "#72B7B2",
    "HT": "#54A24B",
    "PTCwithoutHT": "#4C78A8",
    "PTCwithHT": "#F58518",
    "mPTC": "#B279A2",
}


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)


def target_genes() -> set[str]:
    genes: set[str] = set()
    for gene_list in MODULES.values():
        genes.update(gene_list)
    return genes


def split_geo_line(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\n")], delimiter="\t"))


def parse_sample_group(title: str, source: str) -> str:
    text = f"{title} {source}".lower()
    if "ptc w/ ht" in text or "with ht in background" in text:
        return "PTCwithHT"
    if "ptc w/o ht" in text or "without ht in background" in text:
        return "PTCwithoutHT"
    if title.lower().startswith("ht sample") or "hashimoto" in text:
        return "HT"
    if title.lower().startswith("mptc sample") or "micro ptc" in text:
        return "mPTC"
    if title.lower().startswith("tn sample") or "normal histology" in text:
        return "TN"
    raise ValueError(f"Cannot classify sample: title={title!r}, source={source!r}")


def read_series_matrix() -> tuple[pd.DataFrame, pd.DataFrame]:
    if not SERIES.exists():
        raise FileNotFoundError(f"Missing {SERIES}")

    meta_rows: dict[str, list[str]] = {}
    matrix_lines: list[str] = []
    in_table = False
    with gzip.open(SERIES, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                matrix_lines.append(line)
                continue
            for key in ["!Sample_title", "!Sample_geo_accession", "!Sample_source_name_ch1"]:
                if line.startswith(key):
                    parts = split_geo_line(line)
                    meta_rows[key] = parts[1:]
                    break

    expr = pd.read_csv(io.StringIO("".join(matrix_lines)), sep="\t", dtype={"ID_REF": str})
    expr["ID_REF"] = expr["ID_REF"].astype(str)

    titles = meta_rows["!Sample_title"]
    accessions = meta_rows["!Sample_geo_accession"]
    sources = meta_rows["!Sample_source_name_ch1"]
    sample_rows = []
    for idx, (gsm, title, source) in enumerate(zip(accessions, titles, sources, strict=True), start=1):
        group = parse_sample_group(title, source)
        sample_rows.append(
            {
                "sample": gsm,
                "gsm": gsm,
                "sample_index": idx,
                "title": title,
                "source": source,
                "group": group,
            }
        )
    sample_df = pd.DataFrame(sample_rows)
    return expr, sample_df


def read_annotation(genes: set[str]) -> pd.DataFrame:
    if not ANNOT.exists():
        raise FileNotFoundError(f"Missing {ANNOT}")

    rows: list[str] = []
    in_table = False
    with gzip.open(ANNOT, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                rows.append(line)

    annot = pd.read_csv(io.StringIO("".join(rows)), sep="\t", dtype=str, low_memory=False)
    annot = annot[["ID", "Gene symbol", "Gene title"]].fillna("")
    hits = []
    for row in annot.itertuples(index=False):
        probe_id = str(row.ID)
        symbols = [x.strip() for x in str(getattr(row, "_1")).split("///") if x.strip()]
        for symbol in symbols:
            if symbol in genes:
                hits.append(
                    {
                        "probe_id": probe_id,
                        "gene": symbol,
                        "all_symbols": getattr(row, "_1"),
                        "gene_title": getattr(row, "_2"),
                    }
                )
    out = pd.DataFrame(hits).drop_duplicates(["probe_id", "gene"]).sort_values(["gene", "probe_id"])
    missing = sorted(genes - set(out["gene"]))
    missing_df = pd.DataFrame({"gene": missing, "status": "not_mapped_on_GPL6244"})
    out.to_csv(TABLES / "T01_gse138198_probe_to_target_gene_map.tsv", sep="\t", index=False)
    missing_df.to_csv(TABLES / "T02_gse138198_target_genes_missing_from_gpl6244.tsv", sep="\t", index=False)
    return out


def build_gene_expression(expr: pd.DataFrame, sample_df: pd.DataFrame, probe_map: pd.DataFrame) -> pd.DataFrame:
    value_cols = sample_df["sample"].tolist()
    filtered = expr[expr["ID_REF"].isin(set(probe_map["probe_id"]))].copy()
    long = filtered.melt(id_vars="ID_REF", value_vars=value_cols, var_name="sample", value_name="expression")
    long["expression"] = pd.to_numeric(long["expression"], errors="coerce")
    long = long.merge(probe_map[["probe_id", "gene"]], left_on="ID_REF", right_on="probe_id", how="inner")
    gene_expr = (
        long.groupby(["sample", "gene"], as_index=False)
        .agg(expression=("expression", "mean"), n_probes=("probe_id", "nunique"))
        .merge(sample_df[["sample", "group", "title", "source"]], on="sample", how="left")
    )
    gene_expr = gene_expr.sort_values(["group", "sample", "gene"])
    return gene_expr


def score_modules(gene_df: pd.DataFrame) -> pd.DataFrame:
    matrix = gene_df.pivot(index=["sample", "group"], columns="gene", values="expression")
    z = matrix.copy()
    for gene in z.columns:
        sd = z[gene].std(ddof=0)
        z[gene] = (z[gene] - z[gene].mean()) / sd if sd and not pd.isna(sd) else 0.0

    out = matrix.reset_index()[["sample", "group"]].copy()
    for module, module_genes in MODULES.items():
        available = [gene for gene in module_genes if gene in z.columns]
        out[f"{module}_score"] = z[available].mean(axis=1).to_numpy() if available else np.nan
        out[f"{module}_raw_expression_mean"] = matrix[available].mean(axis=1).to_numpy() if available else np.nan
        out[f"{module}_n_genes"] = len(available)
    out["AP_TLS_composite_score"] = out[["HLA_II_AP_score", "B_TLS_score"]].mean(axis=1)
    out["AP_TLS_composite_raw_expression_mean"] = out[["HLA_II_AP_raw_expression_mean", "B_TLS_raw_expression_mean"]].mean(axis=1)
    out["AP_TLS_composite_n_genes"] = out["HLA_II_AP_n_genes"] + out["B_TLS_n_genes"]
    return out


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


def run_contrasts(sample_scores: pd.DataFrame) -> pd.DataFrame:
    contrasts = [
        ("PTCwithHT_vs_PTCwithoutHT", "PTCwithHT", "PTCwithoutHT"),
        ("HT_vs_TN", "HT", "TN"),
        ("PTCwithHT_vs_TN", "PTCwithHT", "TN"),
        ("PTCwithoutHT_vs_TN", "PTCwithoutHT", "TN"),
        ("PTCwithHT_vs_HT", "PTCwithHT", "HT"),
        ("mPTC_vs_TN", "mPTC", "TN"),
    ]
    rows = []
    for contrast_name, group_a, group_b in contrasts:
        subset = sample_scores[sample_scores["group"].isin([group_a, group_b])].copy()
        y = subset["group"].eq(group_a).to_numpy()
        for module in PRIMARY_MODULES:
            col = f"{module}_score"
            a = subset.loc[subset["group"].eq(group_a), col]
            b = subset.loc[subset["group"].eq(group_b), col]
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
                    "external_dataset": "GSE138198",
                    "contrast": contrast_name,
                    "group_a": group_a,
                    "group_b": group_b,
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
                    "boundary": "HLA/AP expression module, not HLA allele genotype",
                }
            )
    out = pd.DataFrame(rows)
    out["BH_FDR_exact_all_tests"] = bh_fdr(out["exact_two_sided_p"].tolist())
    return out


def make_score_long(sample_scores: pd.DataFrame) -> pd.DataFrame:
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


def plot_figures(sample_scores: pd.DataFrame, contrast_df: pd.DataFrame, score_long: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10.2, 5.1))
    plot_modules = ["HLA_I", "HLA_II_AP", "B_TLS", "AP_TLS_composite"]
    plot_df = score_long[score_long["module"].isin(plot_modules)].copy()
    sns.pointplot(
        data=plot_df,
        x="module",
        y="score",
        hue="group",
        order=plot_modules,
        hue_order=GROUP_ORDER,
        palette=PALETTE,
        dodge=0.55,
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
        alpha=0.82,
        linewidth=0.4,
        edgecolor="black",
        ax=ax,
    )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[: len(GROUP_ORDER)], labels[: len(GROUP_ORDER)], title="", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xlabel("")
    ax.set_ylabel("Module score (mean z of RMA expression)")
    ax.set_title("GSE138198 bulk validation: HT/PTC immune-antigen presentation axis")
    fig.tight_layout()
    fig.savefig(FIGS / "F01_gse138198_group_module_scores.png", dpi=220)
    plt.close(fig)

    primary = contrast_df[contrast_df["contrast"].eq("PTCwithHT_vs_PTCwithoutHT")].sort_values("cohen_d_group_a_minus_group_b", ascending=False)
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    colors = ["#F58518" if x > 0 else "#4C78A8" for x in primary["cohen_d_group_a_minus_group_b"]]
    ax.barh(primary["module"], primary["cohen_d_group_a_minus_group_b"], color=colors, edgecolor="black", linewidth=0.6)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Cohen's d (PTCwithHT - PTCwithoutHT)")
    ax.set_ylabel("")
    ax.set_title("Primary contrast exact-permutation effect sizes")
    for i, row in enumerate(primary.itertuples(index=False)):
        value = row.cohen_d_group_a_minus_group_b
        ax.text(
            value + (0.06 if value >= 0 else -0.06),
            i,
            f"p={row.exact_two_sided_p:.3g}\nAUC={row.auc_group_a_vs_group_b:.2f}",
            va="center",
            ha="left" if value >= 0 else "right",
            fontsize=8.5,
        )
    fig.tight_layout()
    fig.savefig(FIGS / "F02_gse138198_primary_effect_sizes.png", dpi=220)
    plt.close(fig)

    heat = sample_scores.set_index("sample")[[f"{m}_score" for m in MODULE_ORDER if f"{m}_score" in sample_scores]]
    ordered_samples = sample_scores.assign(group_order=sample_scores["group"].map({g: i for i, g in enumerate(GROUP_ORDER)})).sort_values(
        ["group_order", "sample"]
    )["sample"]
    heat = heat.loc[ordered_samples]
    heat.columns = [c.removesuffix("_score") for c in heat.columns]
    fig, ax = plt.subplots(figsize=(8.6, 7.2))
    sns.heatmap(heat, cmap="vlag", center=0, linewidths=0.4, linecolor="white", cbar_kws={"label": "Module score"}, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Sample-level module score heatmap")
    fig.tight_layout()
    fig.savefig(FIGS / "F03_gse138198_module_score_heatmap.png", dpi=220)
    plt.close(fig)

    context = contrast_df[contrast_df["module"].isin(["HLA_II_AP", "B_TLS", "AP_TLS_composite"])].copy()
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    sns.barplot(
        data=context,
        x="contrast",
        y="cohen_d_group_a_minus_group_b",
        hue="module",
        hue_order=["HLA_II_AP", "B_TLS", "AP_TLS_composite"],
        palette=["#E45756", "#54A24B", "#F58518"],
        edgecolor="black",
        linewidth=0.5,
        ax=ax,
    )
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("")
    ax.set_ylabel("Cohen's d")
    ax.set_title("Context contrasts for HLA-II/AP and TLS modules")
    ax.tick_params(axis="x", rotation=25)
    ax.legend(title="", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(FIGS / "F04_gse138198_context_contrast_effects.png", dpi=220)
    plt.close(fig)


def write_report(sample_df: pd.DataFrame, probe_map: pd.DataFrame, contrast_df: pd.DataFrame) -> None:
    primary = contrast_df[contrast_df["contrast"].eq("PTCwithHT_vs_PTCwithoutHT")].copy()
    primary = primary.sort_values("cohen_d_group_a_minus_group_b", ascending=False)
    primary_rows = []
    for r in primary.itertuples(index=False):
        primary_rows.append(
            f"| {r.module} | {r.cohen_d_group_a_minus_group_b:.2f} | "
            f"{r.delta_mean_score_group_a_minus_group_b:.2f} | {r.exact_two_sided_p:.4f} | "
            f"{r.BH_FDR_exact_all_tests:.4f} | {r.auc_group_a_vs_group_b:.2f} |"
        )

    context = contrast_df[
        contrast_df["contrast"].isin(["HT_vs_TN", "PTCwithHT_vs_TN", "PTCwithoutHT_vs_TN"])
        & contrast_df["module"].isin(["HLA_II_AP", "B_TLS", "AP_TLS_composite"])
    ].sort_values(["contrast", "module"])
    context_rows = []
    for r in context.itertuples(index=False):
        context_rows.append(
            f"| {r.contrast} | {r.module} | {r.cohen_d_group_a_minus_group_b:.2f} | "
            f"{r.exact_two_sided_p:.4f} | {r.auc_group_a_vs_group_b:.2f} |"
        )

    counts = sample_df["group"].value_counts().reindex(GROUP_ORDER).dropna().astype(int)
    count_text = ", ".join(f"{k} n={v}" for k, v in counts.items())
    n_probes = int(probe_map["probe_id"].nunique())
    n_genes = int(probe_map["gene"].nunique())
    report = f"""# GSE138198 HT/PTC bulk 외부검증

## 결론

GSE138198은 HT, PTC with HT background, PTC without HT background, mPTC, normal thyroid를 모두 포함하는 bulk expression 자료다. 이번 재분석에서는 {count_text} 구조로, {n_probes}개 probe가 {n_genes}개 target gene에 매핑되었다.

Paper 2에는 이 결과를 **GSE286332와 GSE163203 사이를 잇는 bulk validation**으로 추가할 수 있다. 핵심 주장은 암 코호트 HLA allele association이 아니라, `PTCwithHT`에서 HLA-II/AP 및 TLS/B-cell expression module이 `PTCwithoutHT`보다 높은지 검증하는 것이다.

## Primary contrast: PTCwithHT vs PTCwithoutHT

| module | Cohen's d | delta score | exact p | all-test FDR | AUC |
|---|---:|---:|---:|---:|---:|
{chr(10).join(primary_rows)}

## Disease-context contrasts

| contrast | module | Cohen's d | exact p | AUC |
|---|---|---:|---:|---:|
{chr(10).join(context_rows)}

## 논문 반영 포인트

1. Paper 2 Result에 `independent HT/PTC bulk validation in GSE138198`를 추가한다.
2. Figure ladder는 `GSE286332 bulk -> GSE138198 HT/PTC bulk -> GSE163203 scRNA -> TCGA-THCA network -> spatial TLS` 순서가 가장 강하다.
3. 문구는 반드시 `HLA-II antigen-presentation expression module`로 제한한다. 이 자료는 germline HLA allele을 직접 검정하지 않는다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse138198_htptc_bulk_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse138198_htptc_bulk_validation/figures/`
- Source data: `project/data/external/GSE138198/`
"""
    REPORT.write_text(report, encoding="utf-8")
    (OUT / "GSE138198_HTPTC_BULK_VALIDATION_KR.md").write_text(report, encoding="utf-8")


def main() -> None:
    ensure_dirs()
    genes = target_genes()
    expr, sample_df = read_series_matrix()
    probe_map = read_annotation(genes)
    gene_df = build_gene_expression(expr, sample_df, probe_map)
    sample_scores = score_modules(gene_df).merge(sample_df, on=["sample", "group"], how="left")
    contrast_df = run_contrasts(sample_scores)
    score_long = make_score_long(sample_scores)

    sample_df.to_csv(TABLES / "T03_gse138198_sample_metadata.tsv", sep="\t", index=False)
    gene_df.to_csv(TABLES / "T04_gse138198_target_gene_expression.tsv", sep="\t", index=False)
    sample_scores.to_csv(TABLES / "T05_gse138198_sample_module_scores.tsv", sep="\t", index=False)
    score_long.to_csv(TABLES / "T06_gse138198_module_scores_long.tsv", sep="\t", index=False)
    contrast_df.to_csv(TABLES / "T07_gse138198_exact_contrasts.tsv", sep="\t", index=False)

    plot_figures(sample_scores, contrast_df, score_long)
    write_report(sample_df, probe_map, contrast_df)

    manifest = {
        "dataset": "GSE138198",
        "series_matrix": str(SERIES.relative_to(ROOT)),
        "annotation": str(ANNOT.relative_to(ROOT)),
        "n_samples": int(sample_df.shape[0]),
        "groups": sample_df["group"].value_counts().to_dict(),
        "n_target_genes_mapped": int(probe_map["gene"].nunique()),
        "boundary": "HLA/AP expression module only; no cancer HLA allele/genotype claims",
        "tables": sorted(str(p.relative_to(ROOT)) for p in TABLES.glob("*.tsv")),
        "figures": sorted(str(p.relative_to(ROOT)) for p in FIGS.glob("*.png")),
        "report": str(REPORT.relative_to(ROOT)),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
