#!/usr/bin/env python3
"""Analyze GSE184362 and GSE191288 scRNA-seq target modules.

Use: Paper 2 low-tier scRNA expression generalization only.
Do not use as HT-specific evidence because HT labels are not present in GEO
sample metadata or raw file names.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import math
import re
import tarfile
import tempfile
from itertools import combinations
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse, stats
from scipy.io import mmread
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
DATA184 = ROOT / "project/data/external/GSE184362"
DATA191 = ROOT / "project/data/external/GSE191288"
SERIES184 = DATA184 / "GSE184362_series_matrix.txt.gz"
SERIES191 = DATA191 / "GSE191288_series_matrix.txt.gz"
RAW184 = DATA184 / "GSE184362_RAW.tar"
RAW191 = DATA191 / "GSE191288_RAW.tar"
BASE = ROOT / "project/results/hla_two_paper_synthesis_2026_05_09"
OUT = BASE / "gse184362_gse191288_scrna_generalization"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
REPORT = ROOT / "project/reports/2026_05_09_GSE184362_GSE191288_SCRNA_GENERALIZATION_KR.md"


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
PALETTE184 = {"Paratumor": "#72B7B2", "Tumor": "#F58518", "Metastasis": "#E45756"}
PALETTE191 = {"NonTumor": "#72B7B2", "Tumor": "#F58518"}


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


def read_series_metadata(path: Path) -> pd.DataFrame:
    keys = {
        "!Sample_title",
        "!Sample_geo_accession",
        "!Sample_source_name_ch1",
        "!Sample_description",
        "!Sample_characteristics_ch1",
    }
    meta: dict[str, list[list[str]]] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            key = line.split("\t", 1)[0]
            if key in keys:
                meta.setdefault(key, []).append(split_geo_line(line)[1:])
    accessions = meta["!Sample_geo_accession"][0]
    rows = []
    for i, gsm in enumerate(accessions):
        row = {
            "gsm": gsm,
            "title": meta.get("!Sample_title", [[""] * len(accessions)])[0][i],
            "source": meta.get("!Sample_source_name_ch1", [[""] * len(accessions)])[0][i],
            "description": meta.get("!Sample_description", [[""] * len(accessions)])[0][i] if "!Sample_description" in meta else "",
        }
        for j, arr in enumerate(meta.get("!Sample_characteristics_ch1", []), start=1):
            row[f"characteristic_{j}"] = arr[i]
        rows.append(row)
    return pd.DataFrame(rows)


def decode_array(x) -> list[str]:
    return [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in x]


def classify_184(label: str, source: str) -> str:
    text = f"{label} {source}".lower()
    if re.search(r"(^|_)p($|_)", label.lower()) or "paratumor" in text:
        return "Paratumor"
    if "ln" in label.lower() or "sc" in label.lower() or "metastase" in text or "lymph node" in text or "subcutaneous" in text:
        return "Metastasis"
    return "Tumor"


def tissue_site_184(label: str) -> str:
    if label.endswith("_P"):
        return "Paratumor"
    if label.endswith("_T"):
        return "Tumor"
    if "LN" in label:
        return "LymphNode"
    if label.endswith("_SC"):
        return "Subcutaneous"
    return "Other"


def patient_184(label: str) -> str:
    match = re.search(r"PTC(\d+)", label)
    return f"PTC{match.group(1)}" if match else ""


def read_mtx_member(tar: tarfile.TarFile, name: str):
    fileobj = tar.extractfile(name)
    if fileobj is None:
        raise FileNotFoundError(name)
    with gzip.GzipFile(fileobj=fileobj, mode="rb") as gz:
        try:
            return mmread(gz).tocsr()
        except Exception:
            gz.seek(0)
            return mmread(io.BytesIO(gz.read())).tocsr()


def read_text_member(tar: tarfile.TarFile, name: str) -> io.TextIOWrapper:
    fileobj = tar.extractfile(name)
    if fileobj is None:
        raise FileNotFoundError(name)
    return io.TextIOWrapper(gzip.GzipFile(fileobj=fileobj, mode="rb"))


def module_gene_expression(matrix, genes_by_row: pd.Series, genes: set[str]) -> tuple[pd.DataFrame, dict[str, int]]:
    total_counts = float(matrix.sum())
    gene_to_indices: dict[str, list[int]] = {}
    for idx, gene in enumerate(genes_by_row.astype(str)):
        if gene in genes:
            gene_to_indices.setdefault(gene, []).append(idx)
    rows = []
    for gene in sorted(genes):
        idxs = gene_to_indices.get(gene, [])
        count = float(np.asarray(matrix[idxs, :].sum()).ravel()[0]) if idxs else 0.0
        rows.append(
            {
                "gene": gene,
                "pseudobulk_count": count,
                "log1p_cpm": math.log1p(count / total_counts * 1e6) if total_counts > 0 and idxs else np.nan,
                "n_feature_rows": len(idxs),
            }
        )
    return pd.DataFrame(rows), {gene: len(idxs) for gene, idxs in gene_to_indices.items()}


def analyze_gse184362() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = read_series_metadata(SERIES184)
    with tarfile.open(RAW184) as tar:
        members = tar.getnames()
        matrices = sorted([m for m in members if m.endswith("_matrix.mtx.gz")])
        rows = []
        gene_frames = []
        for matrix_member in matrices:
            match = re.match(r"(GSM\d+)_(.+)_matrix\.mtx\.gz$", matrix_member)
            if not match:
                continue
            gsm, label = match.groups()
            features_member = f"{gsm}_{label}_features.tsv.gz"
            barcodes_member = f"{gsm}_{label}_barcodes.tsv.gz"
            features = pd.read_csv(read_text_member(tar, features_member), sep="\t", header=None, names=["gene_id", "gene"], dtype=str)
            barcodes = pd.read_csv(read_text_member(tar, barcodes_member), sep="\t", header=None, names=["barcode"], dtype=str)
            matrix = read_mtx_member(tar, matrix_member)
            mrow = meta[meta["gsm"].eq(gsm)].iloc[0].to_dict()
            group = classify_184(label, mrow.get("source", ""))
            gene_df, gene_map = module_gene_expression(matrix, features["gene"], target_genes())
            gene_df = gene_df.assign(
                external_dataset="GSE184362",
                sample=f"{gsm}_{label}",
                gsm=gsm,
                label=label,
                patient=patient_184(label),
                compartment=group,
                tissue_site=tissue_site_184(label),
            )
            gene_frames.append(gene_df)
            rows.append(
                {
                    "external_dataset": "GSE184362",
                    "sample": f"{gsm}_{label}",
                    "gsm": gsm,
                    "label": label,
                    "patient": patient_184(label),
                    "compartment": group,
                    "tissue_site": tissue_site_184(label),
                    "title": mrow.get("title", ""),
                    "source": mrow.get("source", ""),
                    "n_cells": int(matrix.shape[1]),
                    "n_features": int(matrix.shape[0]),
                    "total_counts": float(matrix.sum()),
                    "n_target_genes_mapped": int(len(gene_map)),
                    "HT_label_status": "not_available_in_GEO_metadata_or_raw_filenames",
                }
            )
    manifest = pd.DataFrame(rows)
    gene_expr = pd.concat(gene_frames, ignore_index=True)
    return manifest, gene_expr


def parse_191_sample(member: str) -> dict[str, str] | None:
    match = re.match(r"(GSM\d+)_(.+)\.h5$", member)
    if not match:
        return None
    gsm, label = match.groups()
    return {"gsm": gsm, "label": label, "sample": f"{gsm}_{label}", "compartment": "NonTumor" if label == "NT" else "Tumor"}


def read_10x_h5_from_tar(tar: tarfile.TarFile, member: str):
    sample_file = tar.extractfile(member)
    if sample_file is None:
        raise FileNotFoundError(member)
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:
        tmp.write(sample_file.read())
        tmp.flush()
        with h5py.File(tmp.name, "r") as h5:
            g = h5["matrix"]
            shape = tuple(int(x) for x in g["shape"][()])
            mat = sparse.csc_matrix((g["data"][()], g["indices"][()], g["indptr"][()]), shape=shape).tocsr()
            genes = pd.Series(decode_array(g["features"]["name"][()]))
            barcodes = decode_array(g["barcodes"][()])
    return mat, genes, barcodes


def analyze_gse191288() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = read_series_metadata(SERIES191)
    with tarfile.open(RAW191) as tar:
        h5_members = sorted([m for m in tar.getnames() if m.endswith(".h5")])
        rows = []
        gene_frames = []
        for member in h5_members:
            sample = parse_191_sample(member)
            if sample is None:
                continue
            matrix, genes_by_row, barcodes = read_10x_h5_from_tar(tar, member)
            mrow = meta[meta["gsm"].eq(sample["gsm"])].iloc[0].to_dict()
            gene_df, gene_map = module_gene_expression(matrix, genes_by_row, target_genes())
            gene_df = gene_df.assign(
                external_dataset="GSE191288",
                sample=sample["sample"],
                gsm=sample["gsm"],
                label=sample["label"],
                patient=re.sub(r"[LR]$", "", sample["label"]) if sample["label"] != "NT" else "NT",
                compartment=sample["compartment"],
                tissue_site=sample["compartment"],
            )
            gene_frames.append(gene_df)
            rows.append(
                {
                    "external_dataset": "GSE191288",
                    "sample": sample["sample"],
                    "gsm": sample["gsm"],
                    "label": sample["label"],
                    "patient": re.sub(r"[LR]$", "", sample["label"]) if sample["label"] != "NT" else "NT",
                    "compartment": sample["compartment"],
                    "tissue_site": sample["compartment"],
                    "title": mrow.get("title", ""),
                    "source": mrow.get("source", ""),
                    "n_cells": int(matrix.shape[1]),
                    "n_features": int(matrix.shape[0]),
                    "total_counts": float(matrix.sum()),
                    "n_target_genes_mapped": int(len(gene_map)),
                    "HT_label_status": "not_available_in_GEO_metadata_or_raw_filenames",
                }
            )
    manifest = pd.DataFrame(rows)
    gene_expr = pd.concat(gene_frames, ignore_index=True)
    return manifest, gene_expr


def score_modules(gene_expr: pd.DataFrame) -> pd.DataFrame:
    matrix = gene_expr.pivot(index=["external_dataset", "sample", "gsm", "label", "patient", "compartment", "tissue_site"], columns="gene", values="log1p_cpm")
    z = matrix.copy()
    for gene in z.columns:
        sd = z[gene].std(ddof=0)
        z[gene] = (z[gene] - z[gene].mean()) / sd if sd and not pd.isna(sd) else 0.0
    out = matrix.reset_index()[["external_dataset", "sample", "gsm", "label", "patient", "compartment", "tissue_site"]].copy()
    for module, module_genes in MODULES.items():
        available = [gene for gene in module_genes if gene in z.columns and z[gene].notna().any()]
        out[f"{module}_score"] = z[available].mean(axis=1).to_numpy() if available else np.nan
        out[f"{module}_raw_log1p_cpm_mean"] = matrix[available].mean(axis=1).to_numpy() if available else np.nan
        out[f"{module}_n_genes"] = len(available)
    out["AP_TLS_composite_score"] = out[["HLA_II_AP_score", "B_TLS_score"]].mean(axis=1)
    out["AP_TLS_composite_raw_log1p_cpm_mean"] = out[["HLA_II_AP_raw_log1p_cpm_mean", "B_TLS_raw_log1p_cpm_mean"]].mean(axis=1)
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


def run_contrasts(scores: pd.DataFrame) -> pd.DataFrame:
    defs = [
        ("GSE184362", "Tumor_vs_Paratumor", ["Tumor"], ["Paratumor"]),
        ("GSE184362", "Metastasis_vs_Tumor", ["Metastasis"], ["Tumor"]),
        ("GSE191288", "Tumor_vs_NonTumor_qualitative", ["Tumor"], ["NonTumor"]),
    ]
    rows = []
    for dataset, cname, group_a, group_b in defs:
        subset = scores[scores["external_dataset"].eq(dataset) & scores["compartment"].isin(group_a + group_b)].copy()
        if subset.empty:
            continue
        y = subset["compartment"].isin(group_a).to_numpy()
        for module in PRIMARY_MODULES:
            col = f"{module}_score"
            a = subset.loc[subset["compartment"].isin(group_a), col]
            b = subset.loc[subset["compartment"].isin(group_b), col]
            if len(a) < 2 or len(b) < 2:
                p_exact = np.nan
                n_perm = 0
                delta = float(a.mean() - b.mean())
            else:
                delta, p_exact, n_perm = exact_permutation(subset[col].to_numpy(), y)
            try:
                auc = float(roc_auc_score(y.astype(int), subset[col].to_numpy()))
            except ValueError:
                auc = np.nan
            try:
                mwu_p = float(stats.mannwhitneyu(a, b, alternative="two-sided", method="auto").pvalue) if len(a) >= 2 and len(b) >= 2 else np.nan
            except ValueError:
                mwu_p = np.nan
            rows.append(
                {
                    "external_dataset": dataset,
                    "contrast": cname,
                    "group_a": "+".join(group_a),
                    "group_b": "+".join(group_b),
                    "module": module,
                    "n_group_a": int(a.notna().sum()),
                    "n_group_b": int(b.notna().sum()),
                    "mean_group_a": float(a.mean()),
                    "mean_group_b": float(b.mean()),
                    "delta_mean_score_group_a_minus_group_b": float(delta),
                    "cohen_d_group_a_minus_group_b": cohen_d(a, b),
                    "exact_two_sided_p": float(p_exact) if not pd.isna(p_exact) else np.nan,
                    "n_exact_permutations": int(n_perm),
                    "mannwhitney_p": mwu_p,
                    "auc_group_a_vs_group_b": auc,
                    "boundary": "scRNA expression generalization only; no HT-specific or HLA allele/genotype claims",
                }
            )
    out = pd.DataFrame(rows)
    out["BH_FDR_exact_all_tests"] = bh_fdr(out["exact_two_sided_p"].tolist())
    return out


def paired_184(scores: pd.DataFrame) -> pd.DataFrame:
    d = scores[scores["external_dataset"].eq("GSE184362") & scores["compartment"].isin(["Tumor", "Paratumor"])].copy()
    complete = []
    for patient, sub in d.groupby("patient"):
        if {"Tumor", "Paratumor"}.issubset(set(sub["compartment"])):
            complete.append(patient)
    rows = []
    for module in PRIMARY_MODULES:
        col = f"{module}_score"
        deltas = []
        for patient in complete:
            sub = d[d["patient"].eq(patient)]
            tumor = sub.loc[sub["compartment"].eq("Tumor"), col].iloc[0]
            para = sub.loc[sub["compartment"].eq("Paratumor"), col].iloc[0]
            deltas.append(tumor - para)
        try:
            p_w = float(stats.wilcoxon(deltas, zero_method="wilcox").pvalue) if len(deltas) >= 2 else np.nan
        except ValueError:
            p_w = 1.0
        rows.append(
            {
                "external_dataset": "GSE184362",
                "contrast": "paired_Tumor_minus_Paratumor",
                "module": module,
                "n_pairs": len(deltas),
                "mean_delta": float(np.mean(deltas)) if deltas else np.nan,
                "median_delta": float(np.median(deltas)) if deltas else np.nan,
                "wilcoxon_p": p_w,
                "patients": ",".join(complete),
                "boundary": "paired scRNA expression generalization only; no HT-specific or HLA allele/genotype claims",
            }
        )
    out = pd.DataFrame(rows)
    out["BH_FDR_wilcoxon_all_tests"] = bh_fdr(out["wilcoxon_p"].tolist())
    return out


def make_long(scores: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for module in MODULE_ORDER:
        col = f"{module}_score"
        tmp = scores[["external_dataset", "sample", "label", "patient", "compartment", "tissue_site", col]].copy()
        tmp["module"] = module
        tmp["score"] = tmp[col]
        frames.append(tmp[["external_dataset", "sample", "label", "patient", "compartment", "tissue_site", "module", "score"]])
    return pd.concat(frames, ignore_index=True)


def plot_figures(scores: pd.DataFrame, long_scores: pd.DataFrame, contrasts: pd.DataFrame, paired: pd.DataFrame) -> None:
    modules = ["HLA_I", "HLA_II_AP", "B_TLS", "T_IFNG", "AP_TLS_composite"]
    d184 = long_scores[long_scores["external_dataset"].eq("GSE184362") & long_scores["module"].isin(modules)].copy()
    fig, ax = plt.subplots(figsize=(10.0, 4.7))
    sns.boxplot(data=d184, x="module", y="score", hue="compartment", order=modules, hue_order=["Paratumor", "Tumor", "Metastasis"], palette=PALETTE184, showfliers=False, ax=ax)
    sns.stripplot(data=d184, x="module", y="score", hue="compartment", order=modules, hue_order=["Paratumor", "Tumor", "Metastasis"], palette=PALETTE184, dodge=True, jitter=0.12, alpha=0.75, linewidth=0.3, edgecolor="black", ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:3], labels[:3], title="", frameon=False, bbox_to_anchor=(1.01, 1.0), loc="upper left")
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xlabel("")
    ax.set_ylabel("Sample module score")
    ax.set_title("GSE184362 scRNA pseudobulk HLA/AP modules")
    fig.tight_layout()
    fig.savefig(FIGS / "F01_gse184362_compartment_module_scores.png", dpi=220)
    plt.close(fig)

    primary = contrasts[contrasts["external_dataset"].eq("GSE184362") & contrasts["contrast"].eq("Tumor_vs_Paratumor")].sort_values("cohen_d_group_a_minus_group_b")
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    colors = ["#F58518" if x > 0 else "#72B7B2" for x in primary["cohen_d_group_a_minus_group_b"]]
    ax.barh(primary["module"], primary["cohen_d_group_a_minus_group_b"], color=colors, edgecolor="black", linewidth=0.6)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Cohen's d (tumor - paratumor)")
    ax.set_title("GSE184362 tumor-vs-paratumor effect sizes")
    fig.tight_layout()
    fig.savefig(FIGS / "F02_gse184362_tumor_paratumor_effect_sizes.png", dpi=220)
    plt.close(fig)

    d191 = long_scores[long_scores["external_dataset"].eq("GSE191288") & long_scores["module"].isin(modules)].copy()
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    sns.stripplot(data=d191, x="module", y="score", hue="compartment", order=modules, hue_order=["NonTumor", "Tumor"], palette=PALETTE191, dodge=True, jitter=0.12, linewidth=0.4, edgecolor="black", ax=ax)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_xlabel("")
    ax.set_ylabel("Sample module score")
    ax.set_title("GSE191288 qualitative scRNA pseudobulk")
    ax.legend(title="", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(FIGS / "F03_gse191288_qualitative_module_scores.png", dpi=220)
    plt.close(fig)


def write_report(manifest: pd.DataFrame, contrasts: pd.DataFrame, paired: pd.DataFrame) -> None:
    c184 = manifest[manifest["external_dataset"].eq("GSE184362")].groupby("compartment").size().to_dict()
    c191 = manifest[manifest["external_dataset"].eq("GSE191288")].groupby("compartment").size().to_dict()
    primary = contrasts[contrasts["external_dataset"].eq("GSE184362") & contrasts["contrast"].isin(["Tumor_vs_Paratumor", "Metastasis_vs_Tumor"])]
    rows = []
    for r in primary.sort_values(["contrast", "cohen_d_group_a_minus_group_b"], ascending=[True, False]).itertuples(index=False):
        rows.append(
            f"| {r.contrast} | {r.module} | {r.cohen_d_group_a_minus_group_b:.2f} | {r.delta_mean_score_group_a_minus_group_b:.2f} | {r.exact_two_sided_p:.3f} | {r.auc_group_a_vs_group_b:.2f} |"
        )
    paired_rows = []
    for r in paired.sort_values("mean_delta", ascending=False).itertuples(index=False):
        paired_rows.append(f"| {r.module} | {r.n_pairs} | {r.mean_delta:.2f} | {r.wilcoxon_p:.3f} |")
    q191 = contrasts[contrasts["external_dataset"].eq("GSE191288")]
    q_rows = []
    for r in q191.sort_values("delta_mean_score_group_a_minus_group_b", ascending=False).itertuples(index=False):
        q_rows.append(f"| {r.module} | {r.n_group_a} | {r.n_group_b} | {r.delta_mean_score_group_a_minus_group_b:.2f} |")
    report = f"""# GSE184362/GSE191288 scRNA HLA/AP generalization

## 결론

두 dataset 모두 GEO sample metadata와 raw filename에서 HT/WHT/WOHT label을 확인할 수 없었다. 따라서 HT-specific evidence로 쓰면 안 된다.

대신 Paper 2 supplement에서 **low-tier scRNA expression generalization / boundary stress-test**로만 보관한다. GSE184362 구성은 Paratumor n={c184.get("Paratumor", 0)}, Tumor n={c184.get("Tumor", 0)}, Metastasis n={c184.get("Metastasis", 0)}이고, GSE191288은 Tumor n={c191.get("Tumor", 0)}, NonTumor n={c191.get("NonTumor", 0)}이다.

## GSE184362 sample-level contrasts

| contrast | module | Cohen's d | delta score | exact p | AUC |
|---|---|---:|---:|---:|---:|
{chr(10).join(rows)}

## GSE184362 paired tumor-paratumor

| module | n pairs | mean tumor-paratumor delta | Wilcoxon p |
|---|---:|---:|---:|
{chr(10).join(paired_rows)}

## GSE191288 qualitative only

GSE191288은 non-tumor가 1개뿐이어서 p-value를 해석하지 않는다.

| module | n tumor | n non-tumor | tumor - non-tumor delta |
|---|---:|---:|---:|
{chr(10).join(q_rows)}

## 논문 반영 포인트

1. Main claim에는 쓰지 않는다.
2. Supplement registry에는 `HT label unavailable; expression-only scRNA generalization`으로 둔다.
3. `HLA allele`, `genotype`, `risk allele`, `HT-specific` 표현을 쓰지 않는다.

## 산출물

- Tables: `project/results/hla_two_paper_synthesis_2026_05_09/gse184362_gse191288_scrna_generalization/tables/`
- Figures: `project/results/hla_two_paper_synthesis_2026_05_09/gse184362_gse191288_scrna_generalization/figures/`
- Source data: `project/data/external/GSE184362/`, `project/data/external/GSE191288/`
"""
    REPORT.write_text(report, encoding="utf-8")
    (OUT / "GSE184362_GSE191288_SCRNA_GENERALIZATION_KR.md").write_text(report, encoding="utf-8")


def main() -> None:
    ensure_dirs()
    manifest184, genes184 = analyze_gse184362()
    manifest191, genes191 = analyze_gse191288()
    manifest = pd.concat([manifest184, manifest191], ignore_index=True)
    gene_expr = pd.concat([genes184, genes191], ignore_index=True)
    scores = score_modules(gene_expr)
    long_scores = make_long(scores)
    contrasts = run_contrasts(scores)
    paired = paired_184(scores)

    manifest.to_csv(TABLES / "T01_scrna_sample_manifest.tsv", sep="\t", index=False)
    gene_expr.to_csv(TABLES / "T02_scrna_target_gene_pseudobulk_expression.tsv", sep="\t", index=False)
    scores.to_csv(TABLES / "T03_scrna_sample_module_scores.tsv", sep="\t", index=False)
    long_scores.to_csv(TABLES / "T04_scrna_sample_module_scores_long.tsv", sep="\t", index=False)
    contrasts.to_csv(TABLES / "T05_scrna_module_contrasts.tsv", sep="\t", index=False)
    paired.to_csv(TABLES / "T06_gse184362_paired_tumor_paratumor.tsv", sep="\t", index=False)

    plot_figures(scores, long_scores, contrasts, paired)
    write_report(manifest, contrasts, paired)

    manifest_json = {
        "datasets": ["GSE184362", "GSE191288"],
        "n_samples": {dataset: int(n) for dataset, n in manifest.groupby("external_dataset")["sample"].nunique().items()},
        "HT_label_status": "not_available_in_GEO_metadata_or_raw_filenames",
        "boundary": "scRNA expression generalization only; no HT-specific or HLA allele/genotype claims",
        "tables": sorted(str(p.relative_to(ROOT)) for p in TABLES.glob("*.tsv*")),
        "figures": sorted(str(p.relative_to(ROOT)) for p in FIGS.glob("*.png")),
        "report": str(REPORT.relative_to(ROOT)),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest_json, indent=2), encoding="utf-8")
    print(json.dumps(manifest_json, indent=2))


if __name__ == "__main__":
    main()
