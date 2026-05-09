#!/usr/bin/env python3
"""Analyze GSE6004 PTC center/invasion/normal array data for HLA/AP ecology.

Paper 2 use: thyroid cancer expression geography support. This is not
HT-specific and must not be framed as HLA allele/genotype evidence.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import math
import re
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
DATA = ROOT / "project/data/external/GSE6004"
SERIES = DATA / "GSE6004_series_matrix.txt.gz"
ANNOT = DATA / "GPL570_from_geo.soft"
BASE = ROOT / "project/results/hla_two_paper_synthesis_2026_05_09"
OUT = BASE / "gse6004_ptc_invasion_array_validation"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
REPORT = ROOT / "project/reports/2026_05_09_GSE6004_PTC_INVASION_ARRAY_VALIDATION_KR.md"


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
GROUP_ORDER = ["Normal", "Center", "Invasion"]
PALETTE = {"Normal": "#72B7B2", "Center": "#F58518", "Invasion": "#E45756"}


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)


def target_genes() -> set[str]:
    genes: set[str] = set()
    for module_genes in MODULES.values():
        genes.update(module_genes)
    return genes


def split_geo_line(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\n")], delimiter="\t"))


def classify(title: str, description: str, characteristic: str) -> str:
    text = f"{title} {description} {characteristic}".lower()
    if "normal" in text:
        return "Normal"
    if "invasion" in text or "invasive" in text:
        return "Invasion"
    if "center" in text:
        return "Center"
    return "Other"


def patient_id(title: str, description: str) -> str:
    match = re.search(r"\b(T\d+)\b", f"{description} {title}")
    return match.group(1) if match else ""


def read_series() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta: dict[str, list[str]] = {}
    matrix_lines = []
    in_table = False
    wanted = {
        "!Sample_title",
        "!Sample_geo_accession",
        "!Sample_characteristics_ch1",
        "!Sample_description",
    }
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
            key = line.split("\t", 1)[0]
            if key in wanted:
                meta[key] = split_geo_line(line)[1:]
    expr = pd.read_csv(io.StringIO("".join(matrix_lines)), sep="\t", dtype={"ID_REF": str})
    expr["ID_REF"] = expr["ID_REF"].astype(str)
    rows = []
    for gsm, title, desc, char in zip(
        meta["!Sample_geo_accession"],
        meta["!Sample_title"],
        meta["!Sample_description"],
        meta["!Sample_characteristics_ch1"],
        strict=True,
    ):
        rows.append(
            {
                "sample": gsm,
                "gsm": gsm,
                "title": title,
                "description": desc,
                "characteristics": char,
                "patient": patient_id(title, desc),
                "group": classify(title, desc, char),
            }
        )
    return expr, pd.DataFrame(rows)


def read_annotation(genes: set[str]) -> pd.DataFrame:
    lines = []
    in_table = False
    with ANNOT.open("rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                lines.append(line)
    annot = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str, low_memory=False).fillna("")
    symbol_col = "Gene Symbol" if "Gene Symbol" in annot.columns else "Gene symbol"
    title_col = "Gene Title" if "Gene Title" in annot.columns else "Gene title"
    hits = []
    for probe_id, symbols, gene_title in annot[["ID", symbol_col, title_col]].itertuples(index=False, name=None):
        for symbol in [s.strip() for s in symbols.split("///") if s.strip()]:
            if symbol in genes:
                hits.append({"probe_id": probe_id, "gene": symbol, "all_symbols": symbols, "gene_title": gene_title})
    out = pd.DataFrame(hits).drop_duplicates(["probe_id", "gene"]).sort_values(["gene", "probe_id"])
    out.to_csv(TABLES / "T01_gse6004_probe_to_target_gene_map.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": sorted(genes - set(out["gene"])), "status": "not_mapped_on_GPL570"}).to_csv(
        TABLES / "T02_gse6004_target_genes_missing_from_gpl570.tsv", sep="\t", index=False
    )
    return out


def build_gene_expression(expr: pd.DataFrame, sample_df: pd.DataFrame, probe_map: pd.DataFrame) -> pd.DataFrame:
    value_cols = sample_df["sample"].tolist()
    filtered = expr[expr["ID_REF"].isin(set(probe_map["probe_id"]))].copy()
    long = filtered.melt(id_vars="ID_REF", value_vars=value_cols, var_name="sample", value_name="expression")
    long["expression"] = pd.to_numeric(long["expression"], errors="coerce")
    long = long.merge(probe_map[["probe_id", "gene"]], left_on="ID_REF", right_on="probe_id", how="inner")
    return (
        long.groupby(["sample", "gene"], as_index=False)
        .agg(expression=("expression", "mean"), n_probes=("probe_id", "nunique"))
        .merge(sample_df, on="sample", how="left")
        .sort_values(["group", "patient", "sample", "gene"])
    )


def score_modules(gene_df: pd.DataFrame) -> pd.DataFrame:
    matrix = gene_df.pivot(index=["sample", "group", "patient", "title", "description"], columns="gene", values="expression")
    z = matrix.copy()
    for gene in z.columns:
        sd = z[gene].std(ddof=0)
        z[gene] = (z[gene] - z[gene].mean()) / sd if sd and not pd.isna(sd) else 0.0
    out = matrix.reset_index()[["sample", "group", "patient", "title", "description"]].copy()
    for module, genes in MODULES.items():
        available = [gene for gene in genes if gene in z.columns]
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


def exact_permutation(values: np.ndarray, labels: np.ndarray) -> tuple[float, float, int]:
    values = np.asarray(values, float)
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


def run_contrasts(scores: pd.DataFrame) -> pd.DataFrame:
    contrast_defs = [
        ("Center_vs_Normal", "Center", "Normal"),
        ("Invasion_vs_Normal", "Invasion", "Normal"),
        ("Invasion_vs_Center", "Invasion", "Center"),
    ]
    rows = []
    for cname, group_a, group_b in contrast_defs:
        subset = scores[scores["group"].isin([group_a, group_b])].copy()
        y = subset["group"].eq(group_a).to_numpy()
        for module in PRIMARY_MODULES:
            col = f"{module}_score"
            a = subset.loc[subset["group"].eq(group_a), col]
            b = subset.loc[subset["group"].eq(group_b), col]
            delta, p_exact, n_perm = exact_permutation(subset[col].to_numpy(), y)
            try:
                auc = float(roc_auc_score(y.astype(int), subset[col].to_numpy()))
            except ValueError:
                auc = np.nan
            rows.append(
                {
                    "external_dataset": "GSE6004",
                    "contrast": cname,
                    "group_a": group_a,
                    "group_b": group_b,
                    "module": module,
                    "n_group_a": int(a.notna().sum()),
                    "n_group_b": int(b.notna().sum()),
                    "mean_group_a": float(a.mean()),
                    "mean_group_b": float(b.mean()),
                    "delta_mean_score_group_a_minus_group_b": delta,
                    "cohen_d_group_a_minus_group_b": cohen_d(a, b),
                    "exact_permutation_p": p_exact,
                    "n_exact_permutations": n_perm,
                    "mannwhitney_p": float(stats.mannwhitneyu(a, b, alternative="two-sided", method="auto").pvalue),
                    "welch_t_p": float(stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue),
                    "auc_group_a_vs_group_b": auc,
                    "boundary": "PTC expression geography only; no HLA allele/genotype claims",
                }
            )
    out = pd.DataFrame(rows)
    out["BH_FDR_exact_permutation_all_tests"] = bh_fdr(out["exact_permutation_p"].tolist())
    return out


def paired_delta_table(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    complete_patients = []
    for patient, sub in scores.groupby("patient"):
        groups = set(sub["group"])
        if {"Normal", "Center", "Invasion"}.issubset(groups):
            complete_patients.append(patient)
    for module in PRIMARY_MODULES:
        col = f"{module}_score"
        for contrast, group_a, group_b in [
            ("paired_Center_minus_Normal", "Center", "Normal"),
            ("paired_Invasion_minus_Normal", "Invasion", "Normal"),
            ("paired_Invasion_minus_Center", "Invasion", "Center"),
        ]:
            deltas = []
            for patient in complete_patients:
                sub = scores[scores["patient"].eq(patient)]
                a = sub.loc[sub["group"].eq(group_a), col].iloc[0]
                b = sub.loc[sub["group"].eq(group_b), col].iloc[0]
                deltas.append(a - b)
            if len(deltas) < 2:
                p_wilcoxon = np.nan
                p_t = np.nan
            else:
                try:
                    p_wilcoxon = float(stats.wilcoxon(deltas, alternative="two-sided", zero_method="wilcox").pvalue)
                except ValueError:
                    p_wilcoxon = 1.0
                p_t = float(stats.ttest_1samp(deltas, 0.0, nan_policy="omit").pvalue)
            rows.append(
                {
                    "external_dataset": "GSE6004",
                    "contrast": contrast,
                    "module": module,
                    "n_complete_pairs": len(deltas),
                    "mean_delta": float(np.mean(deltas)) if deltas else np.nan,
                    "median_delta": float(np.median(deltas)) if deltas else np.nan,
                    "wilcoxon_p": p_wilcoxon,
                    "paired_t_p": p_t,
                    "complete_patients": ",".join(complete_patients),
                    "boundary": "paired PTC expression geography only; no HLA allele/genotype claims",
                }
            )
    out = pd.DataFrame(rows)
    out["BH_FDR_wilcoxon_all_tests"] = bh_fdr(out["wilcoxon_p"].tolist())
    return out


def make_long(scores: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for module in MODULE_ORDER:
        col = f"{module}_score"
        tmp = scores[["sample", "group", "patient", "title", "description", col]].copy()
        tmp["module"] = module
        tmp["score"] = tmp[col]
        frames.append(tmp[["sample", "group", "patient", "title", "description", "module", "score"]])
    return pd.concat(frames, ignore_index=True)


def plot_figures(scores: pd.DataFrame, long_scores: pd.DataFrame, contrasts: pd.DataFrame, paired: pd.DataFrame) -> None:
    plot_modules = ["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "AP_TLS_composite"]
    plot_df = long_scores[long_scores["module"].isin(plot_modules)].copy()
    fig, ax = plt.subplots(figsize=(10.0, 4.8))
    sns.boxplot(data=plot_df, x="module", y="score", hue="group", order=plot_modules, hue_order=GROUP_ORDER, palette=PALETTE, showfliers=False, ax=ax)
    sns.stripplot(data=plot_df, x="module", y="score", hue="group", order=plot_modules, hue_order=GROUP_ORDER, palette=PALETTE, dodge=True, jitter=0.12, alpha=0.75, linewidth=0, ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:3], labels[:3], title="", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xlabel("")
    ax.set_ylabel("Module score (mean z of target genes)")
    ax.set_title("GSE6004 PTC center/invasion/normal HLA/AP ecology")
    fig.tight_layout()
    fig.savefig(FIGS / "F01_gse6004_group_module_scores.png", dpi=220)
    plt.close(fig)

    primary = contrasts[contrasts["contrast"].eq("Invasion_vs_Normal")].sort_values("cohen_d_group_a_minus_group_b", ascending=True)
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    colors = ["#E45756" if x > 0 else "#72B7B2" for x in primary["cohen_d_group_a_minus_group_b"]]
    ax.barh(primary["module"], primary["cohen_d_group_a_minus_group_b"], color=colors, edgecolor="black", linewidth=0.6)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Cohen's d (invasion - normal)")
    ax.set_ylabel("")
    ax.set_title("Invasive-region effect sizes")
    for i, row in enumerate(primary.itertuples(index=False)):
        val = row.cohen_d_group_a_minus_group_b
        ax.text(val + (0.05 if val >= 0 else -0.05), i, f"p={row.exact_permutation_p:.2g}\nAUC={row.auc_group_a_vs_group_b:.2f}", va="center", ha="left" if val >= 0 else "right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(FIGS / "F02_gse6004_invasion_effect_sizes.png", dpi=220)
    plt.close(fig)

    paired_plot = paired[paired["contrast"].eq("paired_Invasion_minus_Normal")].sort_values("mean_delta", ascending=True)
    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    ax.barh(paired_plot["module"], paired_plot["mean_delta"], color="#E45756", edgecolor="black", linewidth=0.6)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Mean paired delta (invasion - normal)")
    ax.set_ylabel("")
    ax.set_title("Within-patient paired trend (4 complete triplets)")
    fig.tight_layout()
    fig.savefig(FIGS / "F03_gse6004_paired_invasion_normal_delta.png", dpi=220)
    plt.close(fig)


def write_report(sample_df: pd.DataFrame, probe_map: pd.DataFrame, contrasts: pd.DataFrame, paired: pd.DataFrame) -> None:
    primary = contrasts[contrasts["contrast"].eq("Invasion_vs_Normal")].sort_values("cohen_d_group_a_minus_group_b", ascending=False)
    rows = []
    for r in primary.itertuples(index=False):
        rows.append(
            f"| {r.module} | {r.cohen_d_group_a_minus_group_b:.2f} | {r.delta_mean_score_group_a_minus_group_b:.2f} | "
            f"{r.exact_permutation_p:.3f} | {r.BH_FDR_exact_permutation_all_tests:.3f} | {r.auc_group_a_vs_group_b:.2f} |"
        )
    paired_primary = paired[paired["contrast"].eq("paired_Invasion_minus_Normal")].sort_values("mean_delta", ascending=False)
    paired_rows = []
    for r in paired_primary.itertuples(index=False):
        paired_rows.append(f"| {r.module} | {r.n_complete_pairs} | {r.mean_delta:.2f} | {r.wilcoxon_p:.3f} |")
    group_counts = sample_df.groupby("group").size().to_dict()
    report = f"""# GSE6004 PTC center/invasion array 외부검증

## 결론

GSE6004는 PTC 조직에서 normal thyroid, tumor center, invasive area를 비교한 GPL570 microarray 자료다. 이번 분석은 GEO series matrix와 GPL570 annotation을 직접 파싱해 HLA/AP/TLS target gene module을 재계산했다. 표본 구성은 Normal n={group_counts.get("Normal", 0)}, Center n={group_counts.get("Center", 0)}, Invasion n={group_counts.get("Invasion", 0)} 이다.

Paper 2에는 이 결과를 positive evidence로 쓰지 않는다. HT-specific cohort가 아니고, invasion-vs-normal 방향도 HLA/AP/TLS 상승이 아니므로 **negative/specificity stress-test**로만 보관한다. 핵심 외부검증은 GSE138198, GSE163203, GSE213647, GSE248205에 둔다.

## Primary contrast: invasion vs normal

| module | Cohen's d | delta score | exact permutation p | FDR | AUC |
|---|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

## Paired complete triplets

T2/T7/T8/T18 네 환자는 normal-center-invasion triplet이 모두 있어 paired delta를 보조로 계산했다.

| module | n pairs | mean invasion-normal delta | Wilcoxon p |
|---|---:|---:|---:|
{chr(10).join(paired_rows)}

## 논문 반영 포인트

1. Main text에는 넣지 않는다. Reviewer가 "모든 PTC geography에서 항상 HLA/AP/TLS가 올라가냐"라고 물을 때 보조 stress-test로만 사용한다.
2. Supplement에도 넣는다면 "not universal across tumor geography datasets"라는 한계/특이성 패널로 둔다.
3. `HLA allele`, `genotype`, `risk allele`, `HT-specific` 표현은 쓰지 않는다.

## 산출물

- Mapped target genes: {probe_map["gene"].nunique()}
- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse6004_ptc_invasion_array_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse6004_ptc_invasion_array_validation/figures/`
- Source data: `project/data/external/GSE6004/`
"""
    REPORT.write_text(report, encoding="utf-8")
    (OUT / "GSE6004_PTC_INVASION_ARRAY_VALIDATION_KR.md").write_text(report, encoding="utf-8")


def main() -> None:
    ensure_dirs()
    genes = target_genes()
    expr, sample_df = read_series()
    probe_map = read_annotation(genes)
    gene_expr = build_gene_expression(expr, sample_df, probe_map)
    scores = score_modules(gene_expr)
    long_scores = make_long(scores)
    contrasts = run_contrasts(scores)
    paired = paired_delta_table(scores)

    sample_df.to_csv(TABLES / "T03_gse6004_sample_metadata.tsv", sep="\t", index=False)
    gene_expr.to_csv(TABLES / "T04_gse6004_target_gene_expression.tsv", sep="\t", index=False)
    scores.to_csv(TABLES / "T05_gse6004_sample_module_scores.tsv", sep="\t", index=False)
    long_scores.to_csv(TABLES / "T06_gse6004_module_scores_long.tsv", sep="\t", index=False)
    contrasts.to_csv(TABLES / "T07_gse6004_module_contrasts.tsv", sep="\t", index=False)
    paired.to_csv(TABLES / "T08_gse6004_paired_complete_triplet_deltas.tsv", sep="\t", index=False)

    plot_figures(scores, long_scores, contrasts, paired)
    write_report(sample_df, probe_map, contrasts, paired)

    manifest = {
        "dataset": "GSE6004",
        "n_samples": int(sample_df.shape[0]),
        "groups": {str(k): int(v) for k, v in sample_df.groupby("group").size().items()},
        "n_complete_normal_center_invasion_triplets": int(
            sum({"Normal", "Center", "Invasion"}.issubset(set(sub["group"])) for _, sub in sample_df.groupby("patient"))
        ),
        "n_target_genes_mapped": int(probe_map["gene"].nunique()),
        "boundary": "PTC expression geography only; no HLA allele/genotype claims",
        "tables": sorted(str(p.relative_to(ROOT)) for p in TABLES.glob("*.tsv*")),
        "figures": sorted(str(p.relative_to(ROOT)) for p in FIGS.glob("*.png")),
        "report": str(REPORT.relative_to(ROOT)),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
