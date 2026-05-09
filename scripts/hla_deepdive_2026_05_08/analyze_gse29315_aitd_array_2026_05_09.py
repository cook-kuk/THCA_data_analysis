#!/usr/bin/env python3
"""Analyze GSE29315 HT vs thyroid hyperplasia array data.

Use: Paper 4 tissue-expression mechanism support, not HLA genotype evidence.
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
DATA = ROOT / "project/data/external/GSE29315"
SERIES = DATA / "GSE29315_series_matrix.txt.gz"
ANNOT = DATA / "GPL8300.annot.gz"
BASE = ROOT / "project/results/hla_two_paper_synthesis_2026_05_09"
OUT = BASE / "gse29315_aitd_array_validation"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
REPORT = ROOT / "project/reports/2026_05_09_GSE29315_AITD_ARRAY_VALIDATION_KR.md"


MODULES = {
    "HLA_I": ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5"],
    "HLA_II_AP": [
        "HLA-DRA",
        "HLA-DRB1",
        "HLA-DRB3",
        "HLA-DRB4",
        "HLA-DRB5",
        "HLA-DQA1",
        "HLA-DQB1",
        "HLA-DPA1",
        "HLA-DPB1",
        "HLA-DMA",
        "HLA-DMB",
        "CD74",
        "CIITA",
        "RFX5",
    ],
    "B_TLS": ["MS4A1", "CD79A", "CD79B", "CD19", "CD27", "MZB1", "JCHAIN", "CXCL13", "CCL19", "CCL21", "LTB", "CCR7"],
    "T_IFNG": ["CD3D", "CD3E", "CD4", "CD8A", "IFNG", "GZMB", "PRF1", "CXCL9", "CXCL10", "STAT1", "IRF1"],
    "Myeloid_DC": ["LST1", "LYZ", "FCER1A", "CD1C", "LAMP3", "ITGAX", "CD68", "C1QA", "C1QB"],
    "Thyrocyte": ["TG", "TPO", "EPCAM", "KRT8", "KRT18", "SLC5A5", "TSHR", "FOXE1", "NKX2-1"],
    "CD74_MIF_axis": ["CD74", "MIF"],
}

MODULE_ORDER = ["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "Myeloid_DC", "Thyrocyte", "CD74_MIF_axis", "AP_TLS_composite"]
PRIMARY_MODULES = ["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "Myeloid_DC", "CD74_MIF_axis", "AP_TLS_composite"]
GROUP_ORDER = ["Hyperplasia", "Hashimoto", "PTC", "FVPTC", "FTC", "FA", "Hurthle"]
PALETTE = {
    "Hyperplasia": "#72B7B2",
    "Hashimoto": "#54A24B",
    "PTC": "#F58518",
    "FVPTC": "#B279A2",
    "FTC": "#E45756",
    "FA": "#4C78A8",
    "Hurthle": "#9D755D",
}


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)


def target_genes() -> set[str]:
    out: set[str] = set()
    for genes in MODULES.values():
        out.update(genes)
    return out


def split_geo_line(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\n")], delimiter="\t"))


def classify(description: str, title: str) -> str:
    text = f"{description} {title}".lower()
    if "hashimoto" in text or title.upper().startswith("HS"):
        return "Hashimoto"
    if "hyperplasia" in text or title.upper().startswith("HYP"):
        return "Hyperplasia"
    if "follicular variant of papillary" in text or title.upper().startswith("FVPTC"):
        return "FVPTC"
    if "papillary thyroid carcinoma" in text or title.upper().startswith("PTC"):
        return "PTC"
    if "follicular thyroid carcinoma" in text or title.upper().startswith("FTC"):
        return "FTC"
    if "follicular adenoma" in text or title.upper().startswith("FA"):
        return "FA"
    if "hurthle" in text:
        return "Hurthle"
    return "Other"


def read_series() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta: dict[str, list[str]] = {}
    matrix_lines = []
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
            for key in ["!Sample_title", "!Sample_geo_accession", "!Sample_characteristics_ch1"]:
                if line.startswith(key):
                    meta[key] = split_geo_line(line)[1:]
                    break
    expr = pd.read_csv(io.StringIO("".join(matrix_lines)), sep="\t", dtype={"ID_REF": str})
    expr["ID_REF"] = expr["ID_REF"].astype(str)
    rows = []
    for gsm, title, desc in zip(meta["!Sample_geo_accession"], meta["!Sample_title"], meta["!Sample_characteristics_ch1"], strict=True):
        rows.append({"sample": gsm, "title": title, "clinical_description": desc, "group": classify(desc, title)})
    return expr, pd.DataFrame(rows)


def read_annotation(genes: set[str]) -> pd.DataFrame:
    rows = []
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
    annot = pd.read_csv(io.StringIO("".join(rows)), sep="\t", dtype=str, low_memory=False).fillna("")
    hits = []
    for row in annot[["ID", "Gene symbol", "Gene title"]].itertuples(index=False, name=None):
        probe_id, symbols, gene_title = row
        for symbol in [s.strip() for s in symbols.split("///") if s.strip()]:
            if symbol in genes:
                hits.append({"probe_id": probe_id, "gene": symbol, "all_symbols": symbols, "gene_title": gene_title})
    out = pd.DataFrame(hits).drop_duplicates(["probe_id", "gene"]).sort_values(["gene", "probe_id"])
    out.to_csv(TABLES / "T01_gse29315_probe_to_target_gene_map.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": sorted(genes - set(out["gene"])), "status": "not_mapped_on_GPL8300"}).to_csv(
        TABLES / "T02_gse29315_target_genes_missing_from_gpl8300.tsv", sep="\t", index=False
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
        .sort_values(["group", "sample", "gene"])
    )


def score_modules(gene_df: pd.DataFrame) -> pd.DataFrame:
    matrix = gene_df.pivot(index=["sample", "group", "title"], columns="gene", values="expression")
    z = matrix.copy()
    for gene in z.columns:
        sd = z[gene].std(ddof=0)
        z[gene] = (z[gene] - z[gene].mean()) / sd if sd and not pd.isna(sd) else 0.0
    out = matrix.reset_index()[["sample", "group", "title"]].copy()
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
    contrasts = [("Hashimoto_vs_Hyperplasia", "Hashimoto", "Hyperplasia"), ("PTC_vs_Hyperplasia", "PTC", "Hyperplasia")]
    rows = []
    for cname, group_a, group_b in contrasts:
        subset = scores[scores["group"].isin([group_a, group_b])].copy()
        y = subset["group"].eq(group_a).to_numpy()
        for module in PRIMARY_MODULES:
            col = f"{module}_score"
            a = subset.loc[subset["group"].eq(group_a), col]
            b = subset.loc[subset["group"].eq(group_b), col]
            delta, p_exact, n_perm = exact_permutation(subset[col].to_numpy(), y)
            mwu = stats.mannwhitneyu(a, b, alternative="two-sided", method="auto")
            rows.append(
                {
                    "external_dataset": "GSE29315",
                    "contrast": cname,
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
                    "mannwhitney_p": float(mwu.pvalue),
                    "auc_group_a_vs_group_b": float(roc_auc_score(y.astype(int), subset[col].to_numpy())),
                    "boundary": "AITD array expression mechanism, not HLA genotype",
                }
            )
    out = pd.DataFrame(rows)
    out["BH_FDR_exact_all_tests"] = bh_fdr(out["exact_two_sided_p"].tolist())
    return out


def make_long(scores: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for module in MODULE_ORDER:
        col = f"{module}_score"
        tmp = scores[["sample", "group", "title", col]].copy()
        tmp["module"] = module
        tmp["score"] = tmp[col]
        frames.append(tmp[["sample", "group", "title", "module", "score"]])
    return pd.concat(frames, ignore_index=True)


def plot(scores_long: pd.DataFrame, contrasts: pd.DataFrame) -> None:
    plot_modules = ["HLA_I", "HLA_II_AP", "B_TLS", "CD74_MIF_axis", "AP_TLS_composite"]
    plot_df = scores_long[scores_long["group"].isin(["Hyperplasia", "Hashimoto", "PTC"]) & scores_long["module"].isin(plot_modules)]
    fig, ax = plt.subplots(figsize=(9.8, 4.8))
    sns.stripplot(data=plot_df, x="module", y="score", hue="group", order=plot_modules, hue_order=["Hyperplasia", "Hashimoto", "PTC"], palette=PALETTE, dodge=True, jitter=0.12, linewidth=0.4, edgecolor="black", ax=ax)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xlabel("")
    ax.set_ylabel("Module score (mean z of array expression)")
    ax.set_title("GSE29315 AITD array module validation")
    ax.legend(title="", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(FIGS / "F01_gse29315_module_scores.png", dpi=220)
    plt.close(fig)

    primary = contrasts[contrasts["contrast"].eq("Hashimoto_vs_Hyperplasia")].sort_values("cohen_d_group_a_minus_group_b", ascending=False)
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    ax.barh(primary["module"], primary["cohen_d_group_a_minus_group_b"], color="#54A24B", edgecolor="black", linewidth=0.6)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Cohen's d (Hashimoto - Hyperplasia)")
    ax.set_ylabel("")
    ax.set_title("Primary HT vs hyperplasia effect sizes")
    for i, row in enumerate(primary.itertuples(index=False)):
        value = row.cohen_d_group_a_minus_group_b
        ax.text(value + (0.05 if value >= 0 else -0.05), i, f"p={row.exact_two_sided_p:.3g}\nAUC={row.auc_group_a_vs_group_b:.2f}", va="center", ha="left" if value >= 0 else "right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(FIGS / "F02_gse29315_ht_effect_sizes.png", dpi=220)
    plt.close(fig)


def write_report(sample_df: pd.DataFrame, probe_map: pd.DataFrame, contrasts: pd.DataFrame) -> None:
    primary = contrasts[contrasts["contrast"].eq("Hashimoto_vs_Hyperplasia")].sort_values("cohen_d_group_a_minus_group_b", ascending=False)
    rows = []
    for r in primary.itertuples(index=False):
        rows.append(
            f"| {r.module} | {r.cohen_d_group_a_minus_group_b:.2f} | {r.delta_mean_score_group_a_minus_group_b:.2f} | "
            f"{r.exact_two_sided_p:.4f} | {r.auc_group_a_vs_group_b:.2f} | {r.BH_FDR_exact_all_tests:.4f} |"
        )
    counts = sample_df["group"].value_counts()
    report = f"""# GSE29315 AITD array 외부검증

## 결론

GSE29315는 Affymetrix U95 array 기반 thyroid neoplasia/thyroiditis cohort이며, Hashimoto n={int(counts.get('Hashimoto', 0))}, hyperplasia n={int(counts.get('Hyperplasia', 0))}가 포함된다. GPL8300에서 target HLA/AP/TLS gene {probe_map['gene'].nunique()}개가 매핑되었다.

이 결과는 Paper 4의 allele genetics를 직접 반복검증하는 자료가 아니라, HT 조직에서 HLA/AP/TLS expression axis가 올라가는지 확인하는 보조 기전검증이다.

## Primary contrast: Hashimoto vs Hyperplasia

| module | Cohen's d | delta score | exact p | AUC | all-test FDR |
|---|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

## 논문 반영 포인트

1. GSE248205 spatial 결과의 독립 array corroboration으로 supplement에 배치한다.
2. 오래된 platform이라 gene coverage가 제한적이므로, main claim은 GSE248205와 genetics 자료에 둔다.
3. genotype/allele replication이라고 쓰지 않는다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse29315_aitd_array_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse29315_aitd_array_validation/figures/`
- Source data: `project/data/external/GSE29315/`
"""
    REPORT.write_text(report, encoding="utf-8")
    (OUT / "GSE29315_AITD_ARRAY_VALIDATION_KR.md").write_text(report, encoding="utf-8")


def main() -> None:
    ensure_dirs()
    genes = target_genes()
    expr, sample_df = read_series()
    probe_map = read_annotation(genes)
    gene_df = build_gene_expression(expr, sample_df, probe_map)
    scores = score_modules(gene_df)
    long = make_long(scores)
    contrasts = run_contrasts(scores)
    sample_df.to_csv(TABLES / "T03_gse29315_sample_metadata.tsv", sep="\t", index=False)
    gene_df.to_csv(TABLES / "T04_gse29315_target_gene_expression.tsv", sep="\t", index=False)
    scores.to_csv(TABLES / "T05_gse29315_sample_module_scores.tsv", sep="\t", index=False)
    long.to_csv(TABLES / "T06_gse29315_module_scores_long.tsv", sep="\t", index=False)
    contrasts.to_csv(TABLES / "T07_gse29315_exact_contrasts.tsv", sep="\t", index=False)
    plot(long, contrasts)
    write_report(sample_df, probe_map, contrasts)
    manifest = {
        "dataset": "GSE29315",
        "n_samples": int(sample_df.shape[0]),
        "groups": sample_df["group"].value_counts().to_dict(),
        "n_target_genes_mapped": int(probe_map["gene"].nunique()),
        "boundary": "AITD expression mechanism only; no HLA genotype claims",
        "tables": sorted(str(p.relative_to(ROOT)) for p in TABLES.glob("*.tsv")),
        "figures": sorted(str(p.relative_to(ROOT)) for p in FIGS.glob("*.png")),
        "report": str(REPORT.relative_to(ROOT)),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
