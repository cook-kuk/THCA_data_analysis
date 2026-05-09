#!/usr/bin/env python3
"""Analyze GSE213647 Korean thyroid bulk RNA-seq for HLA/AP ecology.

Paper 2 use: Korean thyroid expression generalization. This is not HT-specific
and must not be framed as HLA allele/genotype evidence.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import math
import re
import tarfile
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
DATA = ROOT / "project/data/external/GSE213647"
SERIES = DATA / "GSE213647_series_matrix.txt.gz"
RAW_TAR = DATA / "GSE213647_RAW.tar"
ENSEMBL_MAP = ROOT / "project/metadata/ensembl_to_symbol.tsv"
BASE = ROOT / "project/results/hla_two_paper_synthesis_2026_05_09"
OUT = BASE / "gse213647_korean_bulk_validation"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
REPORT = ROOT / "project/reports/2026_05_09_GSE213647_KOREAN_BULK_VALIDATION_KR.md"


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
PALETTE = {"Normal": "#72B7B2", "Tumor": "#F58518", "Benign": "#4C78A8", "Other": "#BAB0AC"}


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


def parse_characteristic(value: str) -> str:
    return value.split(":", 1)[1].strip() if ":" in value else value.strip()


def read_metadata() -> pd.DataFrame:
    meta: dict[str, list[list[str]]] = {}
    with gzip.open(SERIES, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession") or line.startswith("!Sample_characteristics_ch1"):
                parts = split_geo_line(line)
                meta.setdefault(parts[0], []).append(parts[1:])

    titles = meta["!Sample_title"][0]
    accessions = meta["!Sample_geo_accession"][0]
    characteristics = meta["!Sample_characteristics_ch1"]
    rows = []
    for idx, (gsm, title) in enumerate(zip(accessions, titles, strict=True)):
        values = [parse_characteristic(arr[idx]) for arr in characteristics]
        row = {
            "sample": gsm,
            "gsm": gsm,
            "title": title,
            "tissue": values[0] if len(values) > 0 else "",
            "cell_type": values[1] if len(values) > 1 else "",
            "cell_subtype": values[2] if len(values) > 2 else "",
        }
        row["group"] = row["cell_type"] if row["cell_type"] in {"Normal", "Tumor", "Benign"} else "Other"
        rows.append(row)
    return pd.DataFrame(rows)


def build_ensembl_target_map(genes: set[str]) -> pd.DataFrame:
    mapping = pd.read_csv(ENSEMBL_MAP, sep="\t", dtype=str)
    mapping["ensembl_id"] = mapping["ensembl_id"].str.replace(r"\.\d+$", "", regex=True)
    out = mapping[mapping["gene_symbol"].isin(genes)].drop_duplicates(["ensembl_id", "gene_symbol"]).rename(columns={"gene_symbol": "gene"})
    out.to_csv(TABLES / "T01_gse213647_ensembl_target_gene_map.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": sorted(genes - set(out["gene"])), "status": "not_mapped_in_local_ensembl_table"}).to_csv(
        TABLES / "T02_gse213647_target_genes_missing_from_map.tsv", sep="\t", index=False
    )
    return out


def read_counts_from_tar(sample_df: pd.DataFrame, target_map: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    id_to_gene = dict(zip(target_map["ensembl_id"], target_map["gene"], strict=False))
    target_ids = set(id_to_gene)
    member_to_meta = {}
    for member in tarfile.open(RAW_TAR).getnames():
        m = re.match(r"(GSM\d+)_", Path(member).name)
        if m:
            member_to_meta[m.group(1)] = member

    gene_rows = []
    qc_rows = []
    with tarfile.open(RAW_TAR) as tar:
        for row in sample_df.itertuples(index=False):
            member_name = member_to_meta.get(row.gsm)
            if not member_name:
                continue
            counts_by_gene = {gene: 0.0 for gene in target_map["gene"].unique()}
            special_counts = {}
            total_target = 0.0
            fileobj = tar.extractfile(member_name)
            if fileobj is None:
                continue
            with gzip.GzipFile(fileobj=fileobj, mode="rb") as gz:
                for raw in gz:
                    parts = raw.decode("utf-8").rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    gene_id = parts[0]
                    count = float(parts[1])
                    if gene_id.startswith("N_"):
                        special_counts[gene_id] = count
                        continue
                    stable_id = gene_id.split(".", 1)[0]
                    gene = id_to_gene.get(stable_id)
                    if gene:
                        counts_by_gene[gene] += count
                        total_target += count
            for gene, count in counts_by_gene.items():
                gene_rows.append(
                    {
                        "sample": row.sample,
                        "gsm": row.gsm,
                        "title": row.title,
                        "group": row.group,
                        "cell_type": row.cell_type,
                        "cell_subtype": row.cell_subtype,
                        "gene": gene,
                        "raw_count": count,
                    }
                )
            qc_rows.append(
                {
                    "sample": row.sample,
                    "gsm": row.gsm,
                    "title": row.title,
                    "group": row.group,
                    "cell_subtype": row.cell_subtype,
                    "target_gene_count_sum": total_target,
                    **special_counts,
                }
            )
    return pd.DataFrame(gene_rows), pd.DataFrame(qc_rows)


def normalize_gene_expression(gene_df: pd.DataFrame) -> pd.DataFrame:
    totals = gene_df.groupby("sample")["raw_count"].sum().rename("target_total_count").reset_index()
    out = gene_df.merge(totals, on="sample", how="left")
    out["log1p_target_cpm"] = np.log1p(out["raw_count"] / out["target_total_count"].replace(0, np.nan) * 1e6)
    return out


def score_modules(gene_df: pd.DataFrame) -> pd.DataFrame:
    matrix = gene_df.pivot(index=["sample", "group", "cell_subtype", "title"], columns="gene", values="log1p_target_cpm")
    z = matrix.copy()
    for gene in z.columns:
        sd = z[gene].std(ddof=0)
        z[gene] = (z[gene] - z[gene].mean()) / sd if sd and not pd.isna(sd) else 0.0
    out = matrix.reset_index()[["sample", "group", "cell_subtype", "title"]].copy()
    for module, genes in MODULES.items():
        available = [gene for gene in genes if gene in z.columns]
        out[f"{module}_score"] = z[available].mean(axis=1).to_numpy() if available else np.nan
        out[f"{module}_raw_log1p_target_cpm_mean"] = matrix[available].mean(axis=1).to_numpy() if available else np.nan
        out[f"{module}_n_genes"] = len(available)
    out["AP_TLS_composite_score"] = out[["HLA_II_AP_score", "B_TLS_score"]].mean(axis=1)
    out["AP_TLS_composite_raw_log1p_target_cpm_mean"] = out[["HLA_II_AP_raw_log1p_target_cpm_mean", "B_TLS_raw_log1p_target_cpm_mean"]].mean(axis=1)
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


def run_contrasts(scores: pd.DataFrame) -> pd.DataFrame:
    contrast_defs = [
        ("Tumor_vs_Normal_all", scores["group"].eq("Tumor"), scores["group"].eq("Normal")),
        ("PTC_Tumor_vs_PTC_Normal", scores["group"].eq("Tumor") & scores["cell_subtype"].eq("PTC"), scores["group"].eq("Normal") & scores["cell_subtype"].eq("PTC")),
        ("PD_Tumor_vs_PTC_Tumor", scores["group"].eq("Tumor") & scores["cell_subtype"].isin(["PD", "ATC"]), scores["group"].eq("Tumor") & scores["cell_subtype"].eq("PTC")),
    ]
    rows = []
    for cname, mask_a, mask_b in contrast_defs:
        subset = scores[mask_a | mask_b].copy()
        subset["is_group_a"] = mask_a[mask_a | mask_b].to_numpy()
        if subset["is_group_a"].sum() < 3 or (~subset["is_group_a"]).sum() < 3:
            continue
        for module in PRIMARY_MODULES:
            col = f"{module}_score"
            a = subset.loc[subset["is_group_a"], col]
            b = subset.loc[~subset["is_group_a"], col]
            try:
                auc = float(roc_auc_score(subset["is_group_a"].astype(int).to_numpy(), subset[col].to_numpy()))
            except ValueError:
                auc = np.nan
            rows.append(
                {
                    "external_dataset": "GSE213647",
                    "contrast": cname,
                    "module": module,
                    "n_group_a": int(a.notna().sum()),
                    "n_group_b": int(b.notna().sum()),
                    "mean_group_a": float(a.mean()),
                    "mean_group_b": float(b.mean()),
                    "delta_mean_score_group_a_minus_group_b": float(a.mean() - b.mean()),
                    "cohen_d_group_a_minus_group_b": cohen_d(a, b),
                    "mannwhitney_p": float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue),
                    "welch_t_p": float(stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue),
                    "auc_group_a_vs_group_b": auc,
                    "boundary": "Korean bulk HLA/AP expression ecology, not HLA allele genotype",
                }
            )
    out = pd.DataFrame(rows)
    out["BH_FDR_mannwhitney_all_tests"] = bh_fdr(out["mannwhitney_p"].tolist())
    return out


def make_long(scores: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for module in MODULE_ORDER:
        col = f"{module}_score"
        tmp = scores[["sample", "group", "cell_subtype", "title", col]].copy()
        tmp["module"] = module
        tmp["score"] = tmp[col]
        frames.append(tmp[["sample", "group", "cell_subtype", "title", "module", "score"]])
    return pd.concat(frames, ignore_index=True)


def plot_figures(scores: pd.DataFrame, score_long: pd.DataFrame, contrasts: pd.DataFrame) -> None:
    plot_modules = ["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "AP_TLS_composite"]
    plot_df = score_long[score_long["module"].isin(plot_modules) & score_long["group"].isin(["Normal", "Tumor"])].copy()
    fig, ax = plt.subplots(figsize=(10.0, 4.8))
    sns.boxplot(data=plot_df, x="module", y="score", hue="group", order=plot_modules, hue_order=["Normal", "Tumor"], palette=PALETTE, showfliers=False, ax=ax)
    sns.stripplot(data=plot_df.sample(min(len(plot_df), 2200), random_state=7), x="module", y="score", hue="group", order=plot_modules, hue_order=["Normal", "Tumor"], palette=PALETTE, dodge=True, jitter=0.16, alpha=0.25, linewidth=0, ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], title="", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xlabel("")
    ax.set_ylabel("Module score (mean z of target-gene CPM)")
    ax.set_title("GSE213647 Korean thyroid bulk RNA-seq")
    fig.tight_layout()
    fig.savefig(FIGS / "F01_gse213647_normal_tumor_module_scores.png", dpi=220)
    plt.close(fig)

    primary = contrasts[contrasts["contrast"].eq("PTC_Tumor_vs_PTC_Normal")].sort_values("cohen_d_group_a_minus_group_b", ascending=False)
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    colors = ["#F58518" if x > 0 else "#72B7B2" for x in primary["cohen_d_group_a_minus_group_b"]]
    ax.barh(primary["module"], primary["cohen_d_group_a_minus_group_b"], color=colors, edgecolor="black", linewidth=0.6)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Cohen's d (PTC tumor - PTC normal)")
    ax.set_ylabel("")
    ax.set_title("PTC-specific Korean bulk effect sizes")
    for i, row in enumerate(primary.itertuples(index=False)):
        val = row.cohen_d_group_a_minus_group_b
        ax.text(val + (0.05 if val >= 0 else -0.05), i, f"q={row.BH_FDR_mannwhitney_all_tests:.2g}\nAUC={row.auc_group_a_vs_group_b:.2f}", va="center", ha="left" if val >= 0 else "right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(FIGS / "F02_gse213647_ptc_effect_sizes.png", dpi=220)
    plt.close(fig)

    corr_df = scores[["HLA_II_AP_score", "B_TLS_score", "AP_TLS_composite_score", "Thyrocyte_score", "group", "cell_subtype"]].copy()
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    sns.scatterplot(data=corr_df, x="Thyrocyte_score", y="AP_TLS_composite_score", hue="group", palette=PALETTE, alpha=0.55, linewidth=0, ax=ax)
    rho, p = stats.spearmanr(corr_df["Thyrocyte_score"], corr_df["AP_TLS_composite_score"], nan_policy="omit")
    ax.set_title(f"AP/TLS vs thyrocyte differentiation (rho={rho:.2f}, p={p:.1e})")
    ax.set_xlabel("Thyrocyte module score")
    ax.set_ylabel("AP/TLS composite score")
    ax.legend(title="", frameon=False)
    fig.tight_layout()
    fig.savefig(FIGS / "F03_gse213647_ap_tls_vs_thyrocyte.png", dpi=220)
    plt.close(fig)


def write_report(sample_df: pd.DataFrame, target_map: pd.DataFrame, contrasts: pd.DataFrame, scores: pd.DataFrame) -> None:
    primary = contrasts[contrasts["contrast"].eq("PTC_Tumor_vs_PTC_Normal")].sort_values("cohen_d_group_a_minus_group_b", ascending=False)
    rows = []
    for r in primary.itertuples(index=False):
        rows.append(
            f"| {r.module} | {r.cohen_d_group_a_minus_group_b:.2f} | {r.delta_mean_score_group_a_minus_group_b:.2f} | "
            f"{r.mannwhitney_p:.2e} | {r.BH_FDR_mannwhitney_all_tests:.2e} | {r.auc_group_a_vs_group_b:.2f} |"
        )
    rho, p = stats.spearmanr(scores["Thyrocyte_score"], scores["AP_TLS_composite_score"], nan_policy="omit")
    counts = sample_df.groupby(["group", "cell_subtype"]).size().reset_index(name="n").sort_values(["group", "cell_subtype"])
    count_lines = "; ".join(f"{r.group}/{r.cell_subtype} n={r.n}" for r in counts.itertuples(index=False))
    report = f"""# GSE213647 Korean thyroid bulk RNA-seq 외부검증

## 결론

GSE213647은 한국 연구진이 생성한 thyroid bulk RNA-seq 632 samples 자료다. 이번 재분석은 STAR `ReadsPerGene.out.tab` 원자료에서 HLA/AP/TLS target gene만 직접 추출했다. 분석 sample 구성은 {count_lines} 이다.

Paper 2에는 이 결과를 **Korean thyroid expression generalization**으로 넣을 수 있다. 단, HT-specific 자료가 아니며 HLA allele/genotype 주장이 아니다.

## Primary contrast: PTC tumor vs PTC normal

| module | Cohen's d | delta score | Mann-Whitney p | FDR | AUC |
|---|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

## 추가 구조

Across all 632 samples, AP/TLS composite와 thyrocyte differentiation module의 Spearman rho={rho:.2f}, p={p:.2e}. 이 값은 Korean cohort에서 immune/AP axis와 thyroid differentiation state가 어떻게 함께 움직이는지 보여주는 보조 결과다.

## 논문 반영 포인트

1. Paper 2 supplement 또는 external validation panel에 Korean bulk RNA-seq generalization으로 배치한다.
2. HT-overlap claim의 핵심은 GSE138198/GSE163203에 두고, GSE213647은 ancestry/geography-independent thyroid expression robustness로 쓴다.
3. `HLA allele`, `genotype`, `risk allele` 문구를 쓰지 않는다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse213647_korean_bulk_validation/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse213647_korean_bulk_validation/figures/`
- Source data: `project/data/external/GSE213647/`
"""
    REPORT.write_text(report, encoding="utf-8")
    (OUT / "GSE213647_KOREAN_BULK_VALIDATION_KR.md").write_text(report, encoding="utf-8")


def main() -> None:
    ensure_dirs()
    genes = target_genes()
    sample_df = read_metadata()
    target_map = build_ensembl_target_map(genes)
    gene_counts, qc_df = read_counts_from_tar(sample_df, target_map)
    gene_expr = normalize_gene_expression(gene_counts)
    scores = score_modules(gene_expr)
    long = make_long(scores)
    contrasts = run_contrasts(scores)

    sample_df.to_csv(TABLES / "T03_gse213647_sample_metadata.tsv", sep="\t", index=False)
    qc_df.to_csv(TABLES / "T04_gse213647_sample_qc_target_counts.tsv", sep="\t", index=False)
    gene_expr.to_csv(TABLES / "T05_gse213647_target_gene_expression.tsv.gz", sep="\t", index=False, compression="gzip")
    scores.to_csv(TABLES / "T06_gse213647_sample_module_scores.tsv", sep="\t", index=False)
    long.to_csv(TABLES / "T07_gse213647_module_scores_long.tsv.gz", sep="\t", index=False, compression="gzip")
    contrasts.to_csv(TABLES / "T08_gse213647_module_contrasts.tsv", sep="\t", index=False)

    plot_figures(scores, long, contrasts)
    write_report(sample_df, target_map, contrasts, scores)

    manifest = {
        "dataset": "GSE213647",
        "n_samples": int(sample_df.shape[0]),
        "groups": {f"{group}|{subtype}": int(n) for (group, subtype), n in sample_df.groupby(["group", "cell_subtype"]).size().items()},
        "n_target_genes_mapped": int(target_map["gene"].nunique()),
        "boundary": "Korean bulk expression ecology only; no HLA allele/genotype claims",
        "tables": sorted(str(p.relative_to(ROOT)) for p in TABLES.glob("*.tsv*")),
        "figures": sorted(str(p.relative_to(ROOT)) for p in FIGS.glob("*.png")),
        "report": str(REPORT.relative_to(ROOT)),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
