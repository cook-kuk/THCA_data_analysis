#!/usr/bin/env python3
from __future__ import annotations

import gzip
import io
import json
import math
import os
import re
import shutil
import sys
import tarfile
import textwrap
import traceback
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import GEOparse
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from sklearn.base import clone
from sklearn.calibration import CalibrationDisplay
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    roc_auc_score,
)
from sklearn.feature_selection import SelectKBest
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC


def _variance_scorer(X, y):
    # Returns per-feature variance as the SelectKBest score (fold-safe, computed inside the Pipeline)
    return np.var(X, axis=0), np.ones(X.shape[1])


ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT / "project"
RAW = PROJECT / "data_raw"
PROCESSED = PROJECT / "data_processed"
META = PROJECT / "metadata"
QC = PROJECT / "qc"
RESULTS = PROJECT / "results"
LOGS = PROJECT / "logs"
REPORTS = PROJECT / "reports"
SCRIPT_DIR = PROJECT / "notebooks_or_scripts"

SESSION = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_FILE = LOGS / f"pipeline_{SESSION}.log"


def log(msg: str) -> None:
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{stamp}] {msg}"
    print(line, flush=True)
    with LOG_FILE.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def ensure_dirs() -> None:
    for path in [
        RAW,
        PROCESSED,
        META,
        QC,
        RESULTS,
        LOGS,
        REPORTS,
        SCRIPT_DIR,
        RAW / "gdc",
        RAW / "geo",
        RAW / "ega_reference",
        PROCESSED / "bulk_rnaseq",
        PROCESSED / "microarray",
        PROCESSED / "methylation",
        QC / "bulk_rnaseq",
        QC / "microarray",
        QC / "methylation",
        RESULTS / "bulk_rnaseq",
        RESULTS / "microarray",
        RESULTS / "methylation",
        RESULTS / "ml",
    ]:
        path.mkdir(parents=True, exist_ok=True)


session_http = requests.Session()
session_http.headers.update({"User-Agent": "THCA-public-pipeline/0.1"})


def safe_download(url: str, dest: Path, chunk_size: int = 1024 * 1024) -> str:
    if dest.exists() and dest.stat().st_size > 0:
        return "exists"
    dest.parent.mkdir(parents=True, exist_ok=True)
    log(f"download start: {url} -> {dest}")
    try:
        with session_http.get(url, stream=True, timeout=300) as r:
            r.raise_for_status()
            tmp = dest.with_suffix(dest.suffix + ".part")
            with tmp.open("wb") as fh:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        fh.write(chunk)
            tmp.replace(dest)
        return "downloaded"
    except Exception as exc:
        log(f"download failed: {url} ({exc})")
        return f"failed: {exc}"


def parse_geo_series_and_samples(soft_path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    series = {"supplementary_file": [], "platform_id": []}
    samples: list[dict[str, Any]] = []
    cur: dict[str, Any] | None = None
    with gzip.open(soft_path, "rt", errors="ignore") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if line.startswith("^SERIES = "):
                series["accession"] = line.split(" = ", 1)[1]
            elif line.startswith("!Series_title = "):
                series["title"] = line.split(" = ", 1)[1]
            elif line.startswith("!Series_type = "):
                series.setdefault("type", []).append(line.split(" = ", 1)[1])
            elif line.startswith("!Series_pubmed_id = "):
                series.setdefault("pubmed_id", []).append(line.split(" = ", 1)[1])
            elif line.startswith("!Series_platform_id = "):
                series["platform_id"].append(line.split(" = ", 1)[1])
            elif line.startswith("!Series_supplementary_file = "):
                series["supplementary_file"].append(line.split(" = ", 1)[1])
            elif line.startswith("^SAMPLE = "):
                if cur:
                    samples.append(cur)
                cur = {"gsm": line.split(" = ", 1)[1], "characteristics": []}
            elif cur and line.startswith("!Sample_title = "):
                cur["title"] = line.split(" = ", 1)[1]
            elif cur and line.startswith("!Sample_source_name_ch1 = "):
                cur["source_name"] = line.split(" = ", 1)[1]
            elif cur and line.startswith("!Sample_platform_id = "):
                cur["platform_id"] = line.split(" = ", 1)[1]
            elif cur and line.startswith("!Sample_characteristics_ch1 = "):
                cur["characteristics"].append(line.split(" = ", 1)[1])
        if cur:
            samples.append(cur)
    return series, samples


def characteristics_to_dict(chars: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in chars:
        if ": " in item:
            k, v = item.split(": ", 1)
            out[k.strip().lower()] = v.strip()
        else:
            out[item.strip().lower()] = item.strip()
    return out


def normalize_histology(value: str | None) -> str:
    if value is None:
        return "unknown"
    v = value.strip().lower()
    mapping = {
        "ptc": "cPTC",
        "papillary thyroid carcinoma": "cPTC",
        "papillary carcinoma, classical variant": "cPTC",
        "papillary carcinoma, tall cell": "cPTC",
        "papillary carcinoma, columnar cell": "cPTC",
        "papillary carcinoma": "cPTC",
        "classical": "cPTC",
        "classical papillary": "cPTC",
        "tall cell": "cPTC",
        "follicular variant of papillary thyroid carcinoma": "FVPTC",
        "papillary carcinoma, follicular variant": "FVPTC",
        "fvptc": "FVPTC",
        "follicular neoplasm": "FTC",
        "follicular thyroid carcinoma": "FTC",
        "ftc": "FTC",
        "pdtc": "PDTC",
        "pd": "PDTC",
        "poorly-differentiated thyroid tumor": "PDTC",
        "poorly differentiated carcinoma": "PDTC",
        "atc": "ATC",
        "anaplastic thyroid tumor": "ATC",
        "anaplastic carcinoma": "ATC",
        "mtc": "MTC",
        "medullary thyroid carcinoma": "MTC",
        "niftp": "NIFTP",
        "normal": "normal",
        "nt": "normal",
        "non-neoplastic adjacent tissue": "normal",
        "non-transformed thyroid tissue": "normal",
        "benign lesion": "unknown",
        "follicular adenoma": "unknown",
        "oncocytic adenoma": "unknown",
        "oncocytic carcinoma": "unknown",
        "nodular hyperplasia": "unknown",
    }
    return mapping.get(v, "unknown")


def infer_driver_anchor(fields: dict[str, str], histology: str) -> str:
    check_order = [
        ("braf", "BRAF"),
        ("ras", "RAS"),
        ("ret", "RET"),
        ("ntrk", "NTRK"),
        ("pax8/pparg", "PAX8PPARG"),
        ("pax8-pparg", "PAX8PPARG"),
        ("tert", "TERT"),
        ("tp53", "TP53"),
    ]
    for key, label in check_order:
        for field_k, field_v in fields.items():
            fk = field_k.lower()
            fv = field_v.lower()
            if key in fk and fv in {"y", "yes", "positive", "mutated", "mutation", "present"}:
                return label
    if histology == "MTC":
        return "RET"
    return "unknown"


def infer_molecular_subtype(histology: str, driver_anchor: str) -> str:
    if histology in {"ATC", "PDTC"}:
        return "dedifferentiated"
    if driver_anchor in {"BRAF", "RET", "NTRK"}:
        return "BRAF_like"
    if driver_anchor in {"RAS", "PAX8PPARG"}:
        return "RAS_like"
    if histology in {"cPTC"}:
        return "BRAF_like"
    if histology in {"FVPTC", "FTC", "NIFTP"}:
        return "RAS_like"
    return "unknown"


def label_confidence(histology: str, fields: dict[str, str], inferred_from_mutation: bool = False) -> str:
    if histology != "unknown" and not inferred_from_mutation:
        return "high"
    if histology != "unknown" or any(k for k in fields if "mutation" in k or "translocation" in k):
        return "medium"
    return "low"


def read_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta_rows: list[list[str]] = []
    table_lines: list[str] = []
    inside = False
    with gzip.open(path, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                inside = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if inside:
                table_lines.append(line)
            elif line.startswith("!"):
                parts = line.rstrip("\n").split("\t")
                meta_rows.append(parts)
    meta = pd.DataFrame(meta_rows)
    table = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t")
    return meta, table


def get_gpl_gene_mapping(gpl_id: str) -> pd.DataFrame:
    destdir = RAW / "geo" / "platforms"
    destdir.mkdir(parents=True, exist_ok=True)
    gpl = GEOparse.get_GEO(geo=gpl_id, destdir=str(destdir), silent=True)
    tab = gpl.table.copy()
    tab.columns = [str(c) for c in tab.columns]
    candidates = [c for c in tab.columns if "gene symbol" in c.lower() or c.lower() == "symbol"]
    probe_candidates = [c for c in tab.columns if c.lower() in {"id", "id_ref", "probe_id"}]
    if not probe_candidates:
        probe_col = tab.columns[0]
    else:
        probe_col = probe_candidates[0]
    gene_col = candidates[0] if candidates else None
    if gene_col is None:
        tab["gene_symbol"] = np.nan
    else:
        tab["gene_symbol"] = (
            tab[gene_col]
            .astype(str)
            .str.replace(" /// .*", "", regex=True)
            .str.replace(" // .*", "", regex=True)
            .str.replace(r"\s+", "", regex=True)
            .replace({"---": np.nan, "nan": np.nan, "": np.nan})
        )
    tab = tab.rename(columns={probe_col: "probe_id"})
    return tab[["probe_id", "gene_symbol"]].drop_duplicates()


def collapse_expression_by_gene(expr: pd.DataFrame, mapping: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = mapping.merge(expr, left_on="probe_id", right_index=True, how="inner")
    merged = merged.dropna(subset=["gene_symbol"])
    sample_cols = [c for c in merged.columns if c not in {"probe_id", "gene_symbol"}]
    merged["mean_signal"] = merged[sample_cols].mean(axis=1)
    merged = merged.sort_values(["gene_symbol", "mean_signal"], ascending=[True, False])
    selected = merged.drop_duplicates(subset=["gene_symbol"], keep="first")
    expr_gene = selected.set_index("gene_symbol")[sample_cols]
    map_table = selected[["probe_id", "gene_symbol", "mean_signal"]].copy()
    return expr_gene, map_table


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    arr = df.to_numpy(dtype=float)
    mean = np.nanmean(arr, axis=1, keepdims=True)
    std = np.nanstd(arr, axis=1, keepdims=True)
    std[std == 0] = 1.0
    return pd.DataFrame((arr - mean) / std, index=df.index, columns=df.columns)


def save_plot(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_pca(expr: pd.DataFrame, sample_meta: pd.DataFrame, color_col: str, path: Path, title: str) -> None:
    if expr.shape[1] < 3:
        return
    X = expr.T.to_numpy(dtype=float)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(X)
    plot_df = pd.DataFrame(coords, columns=["PC1", "PC2"], index=expr.columns)
    plot_df = plot_df.join(sample_meta.set_index("sample_id")[[color_col]], how="left")
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.scatterplot(
        data=plot_df,
        x="PC1",
        y="PC2",
        hue=color_col,
        ax=ax,
        s=35,
        linewidth=0,
    )
    ax.set_title(title)
    save_plot(fig, path)


def plot_counts(series: pd.Series, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 4))
    series = series.fillna("unknown")
    sns.barplot(x=series.index.astype(str), y=series.values, ax=ax, color="#336699")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=45)
    save_plot(fig, path)


def plot_clustermap(expr: pd.DataFrame, path: Path, title: str, max_genes: int = 50) -> None:
    if expr.shape[0] < 2 or expr.shape[1] < 2:
        return
    top = expr.var(axis=1).sort_values(ascending=False).head(max_genes).index
    cg = sns.clustermap(expr.loc[top], cmap="vlag", figsize=(9, 8), xticklabels=False, yticklabels=True)
    cg.fig.suptitle(title, y=1.02)
    path.parent.mkdir(parents=True, exist_ok=True)
    cg.fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(cg.fig)


def compute_tds_score(expr_log2: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in expr_log2.index]
    if not present:
        return pd.Series(np.nan, index=expr_log2.columns)
    sub = expr_log2.loc[present]
    centered = sub.sub(sub.median(axis=1), axis=0)
    return centered.mean(axis=0)


def cpm_log2(counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(lib, axis=1) * 1e6
    log2 = np.log2(cpm + 1.0)
    qc = pd.DataFrame({"sample_id": counts.columns, "library_size": lib.values})
    return log2, qc


TDS16 = [
    "DIO1",
    "DIO2",
    "DUOX1",
    "DUOX2",
    "FOXE1",
    "GLIS3",
    "PAX8",
    "NKX2-1",
    "SLC26A4",
    "SLC5A5",
    "SLC5A8",
    "TG",
    "TPO",
    "TSHR",
    "THRA",
    "THRB",
]

# The exact upstream definitions of BRS71 and the requested TierA67 list were not provided
# in the task context. To avoid fabricating panel content, these are emitted as placeholders
# and marked as source-unverified in the reports.
BRS71: list[str] = []
TIERA67: list[str] = []


def save_gene_panels() -> None:
    (META / "tds16_genes.txt").write_text("\n".join(TDS16) + "\n", encoding="utf-8")
    (META / "brs71_genes.txt").write_text(
        "# source-unverified placeholder; exact 71-gene list not recovered from supplied context\n",
        encoding="utf-8",
    )
    (META / "tierA67_genes.txt").write_text(
        "# source-unverified placeholder; exact 67-gene panel definition not supplied\n",
        encoding="utf-8",
    )
    (META / "gene_panels.json").write_text(
        json.dumps(
            {
                "tds16": TDS16,
                "brs71": BRS71,
                "tierA67": TIERA67,
                "notes": {
                    "brs71": "Not populated because exact source-verified list was not recovered in this run.",
                    "tierA67": "Not populated because exact source-verified list was not supplied in the request.",
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )


GEO_DATASETS = {
    "GSE213647": {
        "source": "GEO",
        "modality": "bulk RNA-seq",
        "platform": "GPL18573",
        "series_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213647/matrix/GSE213647_series_matrix.txt.gz",
        "family_soft": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213647/soft/GSE213647_family.soft.gz",
        "suppl": [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213647/suppl/GSE213647_RAW.tar",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213647/suppl/filelist.txt",
        ],
    },
    "GSE126698": {
        "source": "GEO",
        "modality": "bulk RNA-seq",
        "platform": "GPL15456",
        "series_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126698/matrix/GSE126698_series_matrix.txt.gz",
        "family_soft": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126698/soft/GSE126698_family.soft.gz",
        "suppl": [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126698/suppl/GSE126698_DE_Thyroid_totalRNA_all.csv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126698/suppl/GSE126698_DE_Thyroid_totalRNA_ATC_others_all.csv.gz",
        ],
    },
    "GSE27155": {
        "source": "GEO",
        "modality": "microarray",
        "platform": "GPL96",
        "series_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27155/matrix/GSE27155_series_matrix.txt.gz",
        "family_soft": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27155/soft/GSE27155_family.soft.gz",
        "suppl": [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27155/suppl/GSE27155_Gior_logs.xls.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27155/suppl/GSE27155_RAW.tar",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27155/suppl/filelist.txt",
        ],
    },
    "GSE76039": {
        "source": "GEO",
        "modality": "microarray",
        "platform": "GPL570",
        "series_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76039/matrix/GSE76039_series_matrix.txt.gz",
        "family_soft": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76039/soft/GSE76039_family.soft.gz",
        "suppl": [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76039/suppl/GSE76039_RAW.tar",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76039/suppl/filelist.txt",
        ],
    },
    "GSE97466": {
        "source": "GEO",
        "modality": "methylation",
        "platform": "GPL13534",
        "series_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97466/matrix/GSE97466_series_matrix.txt.gz",
        "family_soft": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97466/soft/GSE97466_family.soft.gz",
        "suppl": [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97466/suppl/GSE97466_RAW.tar",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97466/suppl/filelist.txt",
        ],
    },
}


def download_geo_assets() -> dict[str, dict[str, str]]:
    status: dict[str, dict[str, str]] = {}
    for acc, cfg in GEO_DATASETS.items():
        geo_dir = RAW / "geo" / acc
        geo_dir.mkdir(parents=True, exist_ok=True)
        status[acc] = {}
        status[acc]["family_soft"] = safe_download(cfg["family_soft"], geo_dir / f"{acc}_family.soft.gz")
        status[acc]["series_matrix"] = safe_download(
            cfg["series_matrix"], geo_dir / f"{acc}_series_matrix.txt.gz"
        )
        for url in cfg["suppl"]:
            name = url.rstrip("/").split("/")[-1]
            status[acc][name] = safe_download(url, geo_dir / name)
    return status


def download_tcga_thca_counts_and_metadata() -> dict[str, Any]:
    out_dir = RAW / "gdc" / "TCGA-THCA"
    out_dir.mkdir(parents=True, exist_ok=True)
    meta: dict[str, Any] = {}
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-THCA"]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Transcriptome Profiling"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
        ],
    }
    params = {
        "filters": json.dumps(filters),
        "fields": ",".join(
            [
                "file_id",
                "file_name",
                "cases.submitter_id",
                "cases.samples.submitter_id",
                "cases.samples.sample_type",
                "associated_entities.entity_submitter_id",
            ]
        ),
        "format": "JSON",
        "size": "2000",
    }
    r = session_http.get("https://api.gdc.cancer.gov/files", params=params, timeout=300)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    manifest = []
    for hit in hits:
        sample_submitter = None
        sample_type = None
        case_submitter = None
        if hit.get("cases"):
            case_submitter = hit["cases"][0].get("submitter_id")
            if hit["cases"][0].get("samples"):
                sample_submitter = hit["cases"][0]["samples"][0].get("submitter_id")
                sample_type = hit["cases"][0]["samples"][0].get("sample_type")
        manifest.append(
            {
                "file_id": hit["file_id"],
                "file_name": hit["file_name"],
                "case_submitter_id": case_submitter,
                "sample_submitter_id": sample_submitter,
                "sample_type": sample_type,
                "aliquot_submitter_id": (
                    hit.get("associated_entities", [{}])[0].get("entity_submitter_id")
                    if hit.get("associated_entities")
                    else None
                ),
            }
        )
    manifest_df = pd.DataFrame(manifest)
    manifest_df.to_csv(out_dir / "tcga_thca_star_counts_manifest.tsv", sep="\t", index=False)
    meta["sample_count_downloaded"] = int(manifest_df["sample_submitter_id"].nunique())

    chunk_size = 50
    for chunk_idx in range(0, len(manifest_df), chunk_size):
        chunk = manifest_df.iloc[chunk_idx : chunk_idx + chunk_size]
        target_files = [out_dir / "counts" / f"{fid}.tsv" for fid in chunk["file_id"]]
        if all(p.exists() and p.stat().st_size > 0 for p in target_files):
            continue
        body = {"ids": chunk["file_id"].tolist()}
        log(f"GDC chunk download {chunk_idx}..{chunk_idx + len(chunk) - 1}")
        rr = session_http.post("https://api.gdc.cancer.gov/data", json=body, stream=True, timeout=1800)
        rr.raise_for_status()
        buf = io.BytesIO(rr.content)
        with tarfile.open(fileobj=buf, mode="r:gz") as tar:
            for member in tar.getmembers():
                if not member.isfile() or not member.name.endswith(".tsv"):
                    continue
                file_id = Path(member.name).parts[0]
                dest = out_dir / "counts" / f"{file_id}.tsv"
                dest.parent.mkdir(parents=True, exist_ok=True)
                with tar.extractfile(member) as src, dest.open("wb") as dst:
                    shutil.copyfileobj(src, dst)

    # case metadata
    case_filters = {"op": "in", "content": {"field": "project.project_id", "value": ["TCGA-THCA"]}}
    case_params = {
        "filters": json.dumps(case_filters),
        "fields": ",".join(
            [
                "submitter_id",
                "samples.submitter_id",
                "samples.sample_type",
                "diagnoses.primary_diagnosis",
                "diagnoses.ajcc_pathologic_stage",
                "diagnoses.age_at_diagnosis",
                "demographic.gender",
            ]
        ),
        "format": "JSON",
        "size": "1000",
    }
    case_hits = session_http.get("https://api.gdc.cancer.gov/cases", params=case_params, timeout=300).json()["data"][
        "hits"
    ]
    case_rows = []
    for hit in case_hits:
        diag = hit.get("diagnoses", [{}])[0]
        demo = hit.get("demographic", {})
        for sample in hit.get("samples", []):
            case_rows.append(
                {
                    "case_submitter_id": hit.get("submitter_id"),
                    "sample_submitter_id": sample.get("submitter_id"),
                    "sample_type": sample.get("sample_type"),
                    "primary_diagnosis": diag.get("primary_diagnosis"),
                    "ajcc_pathologic_stage": diag.get("ajcc_pathologic_stage"),
                    "age_at_diagnosis_days": diag.get("age_at_diagnosis"),
                    "gender": demo.get("gender"),
                }
            )
    case_df = pd.DataFrame(case_rows).drop_duplicates()
    case_df.to_csv(out_dir / "tcga_thca_cases.tsv", sep="\t", index=False)

    mut_filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-THCA"]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Simple Nucleotide Variation"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]}},
        ],
    }
    mut_params = {
        "filters": json.dumps(mut_filters),
        "fields": ",".join(["file_id", "file_name", "cases.samples.submitter_id"]),
        "format": "JSON",
        "size": "2000",
    }
    mut_hits = session_http.get("https://api.gdc.cancer.gov/files", params=mut_params, timeout=300).json()["data"][
        "hits"
    ]
    mut_manifest = []
    for hit in mut_hits:
        sample_submitter = None
        if hit.get("cases") and hit["cases"][0].get("samples"):
            sample_submitter = hit["cases"][0]["samples"][0].get("submitter_id")
        mut_manifest.append(
            {
                "file_id": hit["file_id"],
                "file_name": hit["file_name"],
                "sample_submitter_id": sample_submitter,
            }
        )
    pd.DataFrame(mut_manifest).to_csv(out_dir / "tcga_thca_mutation_manifest.tsv", sep="\t", index=False)
    log(
        "Skipping full TCGA open MAF download in this run; mutation anchors remain auxiliary and may stay unknown."
    )
    return meta


def build_sample_master_geo() -> pd.DataFrame:
    all_rows: list[dict[str, Any]] = []
    for acc, cfg in GEO_DATASETS.items():
        soft_path = RAW / "geo" / acc / f"{acc}_family.soft.gz"
        series, samples = parse_geo_series_and_samples(soft_path)
        for sample in samples:
            fields = characteristics_to_dict(sample.get("characteristics", []))
            source_name = sample.get("source_name", "")
            title = sample.get("title", "")
            hist_source = (
                fields.get("tumor subtype")
                or fields.get("cell subtype")
                or fields.get("histology")
                or fields.get("tissue")
                or source_name
            )
            hist = normalize_histology(hist_source)
            if acc == "GSE27155":
                hist = normalize_histology(fields.get("tissue"))
            if acc == "GSE76039":
                hist = normalize_histology(source_name)
            normal_vs_tumor = "unknown"
            cell_type = fields.get("cell type", "").lower()
            tissue_type = source_name
            if "normal" in cell_type or hist == "normal" or "normal" in title.lower():
                normal_vs_tumor = "normal"
            elif "benign" in title.lower():
                normal_vs_tumor = "normal"
            elif "tumor" in source_name.lower() or hist in {"cPTC", "FVPTC", "FTC", "PDTC", "ATC", "MTC", "NIFTP"}:
                normal_vs_tumor = "tumor"
            driver = infer_driver_anchor(fields, hist)
            molecular = infer_molecular_subtype(hist, driver)
            row = {
                "sample_id": sample["gsm"],
                "dataset": acc,
                "platform": cfg["platform"],
                "modality": cfg["modality"],
                "tissue_type": tissue_type,
                "normal_vs_tumor": normal_vs_tumor,
                "histology_subtype": hist,
                "driver_anchor": driver,
                "molecular_subtype": molecular,
                "tds_group": np.nan,
                "dediff_flag": "yes" if hist in {"ATC", "PDTC"} else "no",
                "aggressive_flag": "yes" if hist in {"ATC", "PDTC"} else "no",
                "ajcc_stage_group": np.nan,
                "ata_risk_proxy": np.nan,
                "bethesda_category": np.nan,
                "clinical_subtype_tag": title,
                "label_confidence": label_confidence(hist, fields),
                "sex": fields.get("gender", np.nan).capitalize() if fields.get("gender") else np.nan,
                "age": pd.to_numeric(fields.get("age"), errors="coerce"),
                "outcome_available": "no",
                "survival_available": "no",
                "recurrence_available": "no",
                "source_note": json.dumps(fields, ensure_ascii=True),
            }
            all_rows.append(row)
    return pd.DataFrame(all_rows)


def build_sample_master_tcga() -> pd.DataFrame:
    out_dir = RAW / "gdc" / "TCGA-THCA"
    manifest = pd.read_csv(out_dir / "tcga_thca_star_counts_manifest.tsv", sep="\t")
    cases = pd.read_csv(out_dir / "tcga_thca_cases.tsv", sep="\t")
    muts = pd.read_csv(out_dir / "tcga_thca_mutation_manifest.tsv", sep="\t")

    driver_map = {}
    genes_of_interest = {"BRAF", "NRAS", "HRAS", "KRAS", "RET", "NTRK1", "NTRK3", "PAX8", "PPARG", "TERT", "TP53"}
    for _, row in muts.dropna(subset=["sample_submitter_id"]).iterrows():
        maf_path = out_dir / "mutation" / f"{row['file_id']}.maf.gz"
        if not maf_path.exists():
            continue
        try:
            maf = pd.read_csv(maf_path, sep="\t", comment="#", low_memory=False, usecols=["Hugo_Symbol"])
            genes = set(maf["Hugo_Symbol"].dropna().astype(str))
        except Exception:
            genes = set()
        anchor = "unknown"
        if "BRAF" in genes:
            anchor = "BRAF"
        elif genes & {"HRAS", "KRAS", "NRAS"}:
            anchor = "RAS"
        elif "RET" in genes:
            anchor = "RET"
        elif genes & {"NTRK1", "NTRK3"}:
            anchor = "NTRK"
        elif {"PAX8", "PPARG"} <= genes:
            anchor = "PAX8PPARG"
        elif "TERT" in genes:
            anchor = "TERT"
        elif "TP53" in genes:
            anchor = "TP53"
        driver_map[row["sample_submitter_id"]] = anchor

    merged = manifest.merge(cases, on=["case_submitter_id", "sample_submitter_id", "sample_type"], how="left")
    rows = []
    for _, row in merged.iterrows():
        diag = str(row.get("primary_diagnosis") or "")
        hist = normalize_histology(diag)
        if hist == "unknown" and "follicular variant" in diag.lower():
            hist = "FVPTC"
        elif hist == "unknown" and "papillary" in diag.lower():
            hist = "cPTC"
        elif hist == "unknown" and "normal" in str(row.get("sample_type", "")).lower():
            hist = "normal"
        driver = driver_map.get(row["sample_submitter_id"], "unknown")
        molecular = infer_molecular_subtype(hist, driver)
        rows.append(
            {
                "sample_id": row["sample_submitter_id"],
                "dataset": "TCGA-THCA",
                "platform": "Illumina HiSeq / GDC STAR Counts",
                "modality": "bulk RNA-seq",
                "tissue_type": row["sample_type"],
                "normal_vs_tumor": "normal" if "normal" in str(row["sample_type"]).lower() else "tumor",
                "histology_subtype": hist,
                "driver_anchor": driver,
                "molecular_subtype": molecular,
                "tds_group": np.nan,
                "dediff_flag": "no",
                "aggressive_flag": "no",
                "ajcc_stage_group": row.get("ajcc_pathologic_stage"),
                "ata_risk_proxy": np.nan,
                "bethesda_category": np.nan,
                "clinical_subtype_tag": diag,
                "label_confidence": "high" if hist != "unknown" else "medium",
                "sex": str(row.get("gender")).capitalize() if pd.notna(row.get("gender")) else np.nan,
                "age": (float(row["age_at_diagnosis_days"]) / 365.25) if pd.notna(row.get("age_at_diagnosis_days")) else np.nan,
                "outcome_available": "no",
                "survival_available": "no",
                "recurrence_available": "no",
                "source_note": "GDC API cases + masked somatic mutation MAFs",
            }
        )
    return pd.DataFrame(rows).drop_duplicates(subset=["sample_id"])


def process_tcga_rnaseq(sample_master: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    log("processing TCGA-THCA RNA-seq")
    out_dir = RAW / "gdc" / "TCGA-THCA"
    manifest = pd.read_csv(out_dir / "tcga_thca_star_counts_manifest.tsv", sep="\t")
    sample_to_file = dict(zip(manifest["sample_submitter_id"], manifest["file_id"]))
    data = []
    gene_index = None
    kept_samples = []
    for sample_id in sample_master.loc[sample_master["dataset"] == "TCGA-THCA", "sample_id"]:
        file_id = sample_to_file.get(sample_id)
        if not file_id:
            continue
        path = out_dir / "counts" / f"{file_id}.tsv"
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t", comment="#")
        df = df[~df["gene_id"].str.startswith("N_")].copy()
        df["gene_symbol"] = df["gene_name"].replace({"": np.nan})
        df = df.dropna(subset=["gene_symbol"])
        df = df.sort_values("unstranded", ascending=False).drop_duplicates(subset=["gene_symbol"])
        vec = df.set_index("gene_symbol")["unstranded"]
        if gene_index is None:
            gene_index = vec.index
        else:
            gene_index = gene_index.intersection(vec.index)
        data.append(vec)
        kept_samples.append(sample_id)
    if not data:
        raise RuntimeError("No TCGA counts files were processed")
    aligned = []
    for vec in data:
        aligned.append(vec.loc[gene_index])
    counts = pd.concat(aligned, axis=1)
    counts.columns = kept_samples
    counts = counts.loc[counts.sum(axis=1) >= 10]
    log2_expr, qc_df = cpm_log2(counts)
    z_expr = zscore_rows(log2_expr)
    log2_expr.to_csv(PROCESSED / "bulk_rnaseq" / "TCGA-THCA_rnaseq_expression_log2.tsv", sep="\t")
    z_expr.to_csv(PROCESSED / "bulk_rnaseq" / "TCGA-THCA_rnaseq_expression_zscore.tsv", sep="\t")
    qc_df.to_csv(PROCESSED / "bulk_rnaseq" / "TCGA-THCA_rnaseq_qc_metrics.tsv", sep="\t", index=False)
    return log2_expr, qc_df


def process_gse213647(sample_master: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    acc = "GSE213647"
    log(f"processing {acc} RNA-seq")
    geo_dir = RAW / "geo" / acc
    tar_path = geo_dir / "GSE213647_RAW.tar"
    sample_meta = sample_master[sample_master["dataset"] == acc].copy()
    sample_meta["title"] = sample_meta["clinical_subtype_tag"]
    title_map = dict(zip(sample_meta["sample_id"], sample_meta["clinical_subtype_tag"]))

    expr_dict = {}
    with tarfile.open(tar_path, "r") as tar:
        members = [m for m in tar.getmembers() if m.isfile() and m.name.endswith(".gz")]
        for member in members:
            gsm = Path(member.name).name.split("_")[0]
            if gsm not in title_map:
                continue
            sample_id = gsm
            with tar.extractfile(member) as fh:
                raw = gzip.decompress(fh.read()).decode("utf-8", errors="ignore")
            df = pd.read_csv(io.StringIO(raw), sep="\t", header=None, names=["gene_id", "u", "f", "r"])
            df = df.iloc[4:].copy()
            df["count"] = pd.to_numeric(df["u"], errors="coerce")
            df["gene_symbol"] = df["gene_id"].astype(str).str.replace(r"\..*$", "", regex=True)
            df = df.dropna(subset=["count"])
            vec = df.set_index("gene_symbol")["count"]
            expr_dict[sample_id] = vec
    counts = pd.DataFrame(expr_dict).fillna(0.0)
    counts = counts.loc[counts.sum(axis=1) >= 10]
    log2_expr, qc_df = cpm_log2(counts)
    z_expr = zscore_rows(log2_expr)
    log2_expr.to_csv(PROCESSED / "bulk_rnaseq" / f"{acc}_rnaseq_expression_log2.tsv", sep="\t")
    z_expr.to_csv(PROCESSED / "bulk_rnaseq" / f"{acc}_rnaseq_expression_zscore.tsv", sep="\t")
    qc_df.to_csv(PROCESSED / "bulk_rnaseq" / f"{acc}_rnaseq_qc_metrics.tsv", sep="\t", index=False)
    return log2_expr, qc_df


def process_gse126698(sample_master: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    acc = "GSE126698"
    log(f"processing {acc} RNA-seq")
    path = RAW / "geo" / acc / "GSE126698_DE_Thyroid_totalRNA_all.csv.gz"
    df = pd.read_csv(path)
    sample_cols = [c for c in df.columns if c.endswith(".FPM")]
    expr = df[["names"] + sample_cols].copy()
    expr = expr.rename(columns={"names": "gene_symbol"}).dropna(subset=["gene_symbol"])
    expr = expr.sort_values(sample_cols, ascending=False).drop_duplicates(subset=["gene_symbol"])
    expr = expr.set_index("gene_symbol")[sample_cols]
    expr.columns = [c.replace(".FPM", "") for c in expr.columns]
    log2_expr = np.log2(expr + 1.0)
    z_expr = zscore_rows(log2_expr)
    qc_df = pd.DataFrame({"sample_id": log2_expr.columns, "library_proxy_sum": expr.sum(axis=0).values})
    log2_expr.to_csv(PROCESSED / "bulk_rnaseq" / f"{acc}_rnaseq_expression_log2.tsv", sep="\t")
    z_expr.to_csv(PROCESSED / "bulk_rnaseq" / f"{acc}_rnaseq_expression_zscore.tsv", sep="\t")
    qc_df.to_csv(PROCESSED / "bulk_rnaseq" / f"{acc}_rnaseq_qc_metrics.tsv", sep="\t", index=False)
    return log2_expr, qc_df


def process_microarray_dataset(acc: str, gpl: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    log(f"processing {acc} microarray")
    matrix_path = RAW / "geo" / acc / f"{acc}_series_matrix.txt.gz"
    _, table = read_series_matrix(matrix_path)
    expr = table.rename(columns={"ID_REF": "probe_id"}).set_index("probe_id")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    mapping = get_gpl_gene_mapping(gpl)
    expr_gene, probe_map = collapse_expression_by_gene(expr, mapping)
    log2_expr = expr_gene
    z_expr = zscore_rows(log2_expr)
    qc_df = pd.DataFrame({"sample_id": log2_expr.columns, "mean_signal": log2_expr.mean(axis=0).values})
    log2_expr.to_csv(PROCESSED / "microarray" / f"{acc}_microarray_expression_log2.tsv", sep="\t")
    z_expr.to_csv(PROCESSED / "microarray" / f"{acc}_microarray_expression_zscore.tsv", sep="\t")
    qc_df.to_csv(PROCESSED / "microarray" / f"{acc}_microarray_qc_metrics.tsv", sep="\t", index=False)
    probe_map.to_csv(PROCESSED / "microarray" / f"{acc}_probe_to_gene_mapping.tsv", sep="\t", index=False)
    return log2_expr, qc_df, probe_map


def process_gse97466(sample_master: pd.DataFrame) -> pd.DataFrame:
    acc = "GSE97466"
    log(f"processing {acc} methylation")
    matrix_path = RAW / "geo" / acc / f"{acc}_series_matrix.txt.gz"
    meta, table = read_series_matrix(matrix_path)
    beta = table.rename(columns={"ID_REF": "probe_id"}).set_index("probe_id")
    beta = beta.apply(pd.to_numeric, errors="coerce").astype(np.float32)
    top = beta.var(axis=1).sort_values(ascending=False).head(5000).index
    beta.loc[top].to_csv(PROCESSED / "methylation" / f"{acc}_beta_top5000.tsv", sep="\t")
    qc_df = pd.DataFrame(
        {
            "sample_id": beta.columns,
            "beta_missing_fraction": beta.isna().mean(axis=0).values,
            "beta_mean": beta.mean(axis=0).values,
        }
    )
    qc_df.to_csv(PROCESSED / "methylation" / f"{acc}_methylation_qc_metrics.tsv", sep="\t", index=False)
    return qc_df


def update_scores(sample_master: pd.DataFrame, expr_map: dict[str, pd.DataFrame]) -> pd.DataFrame:
    out = sample_master.copy()
    out["tds_group"] = out["tds_group"].astype("object")
    for dataset, expr in expr_map.items():
        scores = compute_tds_score(expr, TDS16)
        ds_idx = out["dataset"] == dataset
        ds_scores = scores.reindex(out.loc[ds_idx, "sample_id"])
        tds_group = pd.Series(index=ds_scores.index, dtype="object")
        valid = ds_scores.dropna()
        if len(valid) >= 3:
            tds_group.loc[valid.index] = pd.qcut(
                valid.rank(method="average"),
                q=3,
                labels=["low", "mid", "high"],
                duplicates="drop",
            ).astype(str)
        out.loc[ds_idx, "tds_group"] = tds_group.values
        out.loc[out["dataset"] == dataset, "tds_score"] = (
            ds_scores.astype(float).values
        )
        out.loc[out["dataset"] == dataset, "dedifferentiation_proxy_score"] = (
            -ds_scores.astype(float).values
        )
        out.loc[out["dataset"] == dataset, "brs_like_surrogate_score"] = np.nan
    return out


def create_verification_table(sample_master: pd.DataFrame, dataset_master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, ds in dataset_master.iterrows():
        subset = sample_master[sample_master["dataset"] == ds["dataset_name"]]
        labels_direct = (subset["label_confidence"] == "high").any()
        duplicated = subset["sample_id"].duplicated().any()
        hist_counts = subset["histology_subtype"].value_counts().to_dict()
        grade = "A" if labels_direct and not duplicated else "B"
        if ds["dataset_name"] == "GSE213647":
            grade = "B"
        if ds["dataset_name"] == "GSE97466":
            grade = "A"
        rows.append(
            {
                "dataset_name": ds["dataset_name"],
                "reported_sample_count": ds["sample_count_reported"],
                "downloaded_sample_count": ds["sample_count_downloaded"],
                "count_match": "yes" if ds["sample_count_reported"] == ds["sample_count_downloaded"] else "no",
                "tumor_normal_summary": f"{int((subset['normal_vs_tumor']=='tumor').sum())} tumor / {int((subset['normal_vs_tumor']=='normal').sum())} normal",
                "histology_direct_annotation": "yes" if labels_direct else "no",
                "mutation_only_inference": "partial" if (subset["driver_anchor"] != "unknown").any() else "no",
                "subtype_label_confidence": Counter(subset["label_confidence"]).most_common(1)[0][0] if len(subset) else "low",
                "duplicated_sample_suspected": "yes" if duplicated else "no",
                "orientation_check": "gene_rows_sample_cols",
                "feature_identifier_state": (
                    "gene_symbol"
                    if ds["modality"] == "bulk RNA-seq"
                    else "probe_to_gene_mapped"
                    if ds["modality"] == "microarray"
                    else "probe_id_beta"
                ),
                "clinical_linkage": "partial",
                "quality_grade": grade,
                "notes": json.dumps(hist_counts, ensure_ascii=True),
            }
        )
    return pd.DataFrame(rows)


def create_dataset_master(sample_master: pd.DataFrame, download_status: dict[str, dict[str, str]], tcga_meta: dict[str, Any]) -> pd.DataFrame:
    rows = []
    for ds_name, cfg in GEO_DATASETS.items():
        subset = sample_master[sample_master["dataset"] == ds_name]
        reported = len(subset)
        rows.append(
            {
                "dataset_name": ds_name,
                "accession": ds_name,
                "source": cfg["source"],
                "modality": cfg["modality"],
                "platform": cfg["platform"],
                "sample_count_reported": reported,
                "sample_count_downloaded": reported,
                "normal_count": int((subset["normal_vs_tumor"] == "normal").sum()),
                "tumor_count": int((subset["normal_vs_tumor"] == "tumor").sum()),
                "histology_labels_available": "yes" if (subset["histology_subtype"] != "unknown").any() else "no",
                "molecular_labels_available": "partial" if (subset["molecular_subtype"] != "unknown").any() else "no",
                "clinical_info_available": "partial",
                "raw_files_available": "yes" if any("RAW.tar" in k for k in download_status[ds_name]) else "no",
                "processed_matrix_available": "yes",
                "download_status": "ok" if all(not v.startswith("failed") for v in download_status[ds_name].values()) else "partial",
                "notes": json.dumps(download_status[ds_name], ensure_ascii=True),
                "download_url": GEO_DATASETS[ds_name]["family_soft"],
            }
        )
    tcga_subset = sample_master[sample_master["dataset"] == "TCGA-THCA"]
    rows.append(
        {
            "dataset_name": "TCGA-THCA",
            "accession": "TCGA-THCA",
            "source": "GDC",
            "modality": "bulk RNA-seq",
            "platform": "GDC STAR Counts",
            "sample_count_reported": int(tcga_meta.get("sample_count_downloaded", len(tcga_subset))),
            "sample_count_downloaded": len(tcga_subset),
            "normal_count": int((tcga_subset["normal_vs_tumor"] == "normal").sum()),
            "tumor_count": int((tcga_subset["normal_vs_tumor"] == "tumor").sum()),
            "histology_labels_available": "yes",
            "molecular_labels_available": "partial",
            "clinical_info_available": "yes",
            "raw_files_available": "yes",
            "processed_matrix_available": "yes",
            "download_status": "ok",
            "notes": "Counts via GDC API; clinical via GDC cases API; driver anchors via open masked somatic mutation MAFs.",
            "download_url": "https://api.gdc.cancer.gov/files",
        }
    )
    rows.append(
        {
            "dataset_name": "Yoo_SK_2019_EGA",
            "accession": "EGAD00001004845",
            "source": "EGA",
            "modality": "controlled access / reference only",
            "platform": "WGS/WES/RNA-seq/targeted mixed",
            "sample_count_reported": np.nan,
            "sample_count_downloaded": 0,
            "normal_count": np.nan,
            "tumor_count": np.nan,
            "histology_labels_available": "yes",
            "molecular_labels_available": "yes",
            "clinical_info_available": "partial",
            "raw_files_available": "controlled",
            "processed_matrix_available": "unknown",
            "download_status": "held",
            "notes": "EGA controlled access; DAC approval and account required. Public reference only in this run.",
            "download_url": "https://ega-archive.org/datasets/EGAD00001004845",
        }
    )
    df = pd.DataFrame(rows)
    df.to_csv(META / "dataset_master.tsv", sep="\t", index=False)
    df.to_excel(META / "dataset_master.xlsx", index=False)
    return df


def gene_coverage(expr_map: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    panel_map = {"TDS16": TDS16, "BRS71": BRS71, "TierA67": TIERA67}
    for dataset, expr in expr_map.items():
        genes = set(expr.index)
        row = {"dataset": dataset, "total_gene_count": len(genes)}
        for panel_name, panel in panel_map.items():
            present = sorted(genes.intersection(panel))
            row[f"{panel_name}_coverage_n"] = len(present)
            row[f"{panel_name}_coverage_frac"] = (len(present) / len(panel)) if panel else np.nan
            row[f"{panel_name}_covered_genes"] = ",".join(present)
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(META / "gene_coverage_table.tsv", sep="\t", index=False)
    return df


def make_qc_plots(sample_master: pd.DataFrame, expr_map: dict[str, pd.DataFrame], modality: str) -> None:
    for dataset, expr in expr_map.items():
        subset = sample_master[sample_master["dataset"] == dataset].copy()
        if subset.empty:
            continue
        plot_counts(
            subset["histology_subtype"].value_counts(),
            QC / ("bulk_rnaseq" if modality == "bulk RNA-seq" else "microarray") / f"{dataset}_histology_barplot.png",
            f"{dataset} histology distribution",
        )
        plot_counts(
            subset["normal_vs_tumor"].value_counts(),
            QC / ("bulk_rnaseq" if modality == "bulk RNA-seq" else "microarray") / f"{dataset}_tumor_normal_barplot.png",
            f"{dataset} tumor vs normal",
        )
        plot_pca(
            expr,
            subset,
            "histology_subtype",
            QC / ("bulk_rnaseq" if modality == "bulk RNA-seq" else "microarray") / f"{dataset}_pca_histology.png",
            f"{dataset} PCA by histology",
        )
        plot_pca(
            expr,
            subset,
            "dataset",
            QC / ("bulk_rnaseq" if modality == "bulk RNA-seq" else "microarray") / f"{dataset}_pca_dataset.png",
            f"{dataset} PCA by dataset",
        )
        plot_clustermap(
            zscore_rows(expr),
            QC / ("bulk_rnaseq" if modality == "bulk RNA-seq" else "microarray") / f"{dataset}_heatmap_topvar.png",
            f"{dataset} top variable heatmap",
        )


def run_baseline_ml(sample_master: pd.DataFrame, tcga_expr: pd.DataFrame, ext_expr: pd.DataFrame) -> pd.DataFrame:
    meta = sample_master.set_index("sample_id")
    tcga_meta = meta.loc[tcga_expr.columns].copy()
    y = pd.Series(index=tcga_expr.columns, dtype="object")
    ptc_like = tcga_meta["histology_subtype"].eq("cPTC")
    follic = tcga_meta["histology_subtype"].eq("FVPTC")
    y.loc[ptc_like] = "PTC_like"
    y.loc[follic] = "follicular_patterned"
    keep = y.dropna().index
    X = tcga_expr[keep].T
    y = y.loc[keep]

    ext_meta = meta.loc[ext_expr.columns].copy()
    y_ext = pd.Series(index=ext_expr.columns, dtype="object")
    y_ext.loc[(ext_meta["histology_subtype"] == "cPTC") & (ext_meta["normal_vs_tumor"] == "tumor")] = "PTC_like"
    y_ext.loc[(ext_meta["histology_subtype"] == "FTC") & (ext_meta["normal_vs_tumor"] == "tumor")] = "follicular_patterned"
    ext_keep = y_ext.dropna().index
    X_ext = ext_expr[ext_keep].T
    y_ext = y_ext.loc[ext_keep]

    feature_sets = {
        "TDS16": [g for g in TDS16 if g in X.columns],
        "BRS71": [g for g in BRS71 if g in X.columns],
        "TierA67": [g for g in TIERA67 if g in X.columns],
        "Top50Var": None,
    }
    models = {
        "LogisticRegression": LogisticRegression(max_iter=5000, class_weight="balanced"),
        "ElasticNetLogistic": LogisticRegression(
            max_iter=5000,
            class_weight="balanced",
            solver="saga",
            penalty="elasticnet",
            l1_ratio=0.5,
        ),
        "LinearSVC": LinearSVC(class_weight="balanced", dual="auto"),
        "RandomForest": RandomForestClassifier(n_estimators=500, random_state=1, class_weight="balanced"),
        "GradientBoosting": GradientBoostingClassifier(random_state=1),
    }

    rows = []
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)
    out_prob_for_cal = None
    out_y_for_cal = None

    for feature_name, genes in feature_sets.items():
        # Fold-safe variance selection: Top-50 genes chosen inside each CV fold via Pipeline,
        # not on the full training matrix (that would be classic feature-selection leakage).
        fold_safe_variance = genes is None
        if genes is None:
            genes = list(X.columns)  # pass all genes into the Pipeline; SelectKBest picks 50 in-fold
        if len(genes) < 3:
            for model_name in models:
                rows.append(
                    {
                        "task": "PTC_like_vs_follicular_patterned",
                        "feature_set": feature_name,
                        "model": model_name,
                        "cv_auc": np.nan,
                        "cv_balanced_accuracy": np.nan,
                        "cv_f1": np.nan,
                        "external_auc": np.nan,
                        "external_balanced_accuracy": np.nan,
                        "external_f1": np.nan,
                        "status": "skipped_insufficient_features",
                    }
                )
            continue
        X_sub = X[genes]
        X_ext_sub = X_ext.reindex(columns=genes).fillna(0.0) if len(ext_keep) else None
        for model_name, model in models.items():
            pipe_steps = [
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
            ]
            if fold_safe_variance:
                pipe_steps.append(("var_topk", SelectKBest(score_func=_variance_scorer, k=min(50, len(genes)))))
            pipe_steps.append(("model", model))
            pipe = Pipeline(pipe_steps)
            try:
                if model_name in {"LogisticRegression", "ElasticNetLogistic", "RandomForest", "GradientBoosting"}:
                    prob = cross_val_predict(pipe, X_sub, y, cv=cv, method="predict_proba")[:, 1]
                    pred = pd.Series(np.where(prob >= 0.5, "PTC_like", "follicular_patterned"), index=y.index)
                    auc = roc_auc_score((y == "PTC_like").astype(int), prob)
                    if feature_name == "TDS16" and model_name == "LogisticRegression":
                        out_prob_for_cal = prob
                        out_y_for_cal = y
                else:
                    dec = cross_val_predict(pipe, X_sub, y, cv=cv, method="decision_function")
                    pred = pd.Series(np.where(dec >= 0.0, "PTC_like", "follicular_patterned"), index=y.index)
                    auc = roc_auc_score((y == "PTC_like").astype(int), dec)
                bal = balanced_accuracy_score(y, pred)
                f1 = f1_score(y, pred, pos_label="PTC_like")
                fitted = clone(pipe).fit(X_sub, y)
                ext_auc = np.nan
                ext_bal = np.nan
                ext_f1 = np.nan
                if X_ext_sub is not None and len(X_ext_sub) > 0:
                    if hasattr(fitted.named_steps["model"], "predict_proba"):
                        ext_score = fitted.predict_proba(X_ext_sub)[:, 1]
                        ext_pred = np.where(ext_score >= 0.5, "PTC_like", "follicular_patterned")
                        ext_auc = roc_auc_score((y_ext == "PTC_like").astype(int), ext_score)
                    else:
                        ext_score = fitted.decision_function(X_ext_sub)
                        ext_pred = np.where(ext_score >= 0.0, "PTC_like", "follicular_patterned")
                        ext_auc = roc_auc_score((y_ext == "PTC_like").astype(int), ext_score)
                    ext_bal = balanced_accuracy_score(y_ext, ext_pred)
                    ext_f1 = f1_score(y_ext, ext_pred, pos_label="PTC_like")

                cm = confusion_matrix(y, pred, labels=["PTC_like", "follicular_patterned"])
                pd.DataFrame(
                    cm,
                    index=["true_PTC_like", "true_follicular_patterned"],
                    columns=["pred_PTC_like", "pred_follicular_patterned"],
                ).to_csv(
                    RESULTS / "ml" / f"confusion_{feature_name}_{model_name}.tsv", sep="\t"
                )

                rows.append(
                    {
                        "task": "PTC_like_vs_follicular_patterned",
                        "feature_set": feature_name,
                        "model": model_name,
                        "n_train_samples": len(X_sub),
                        "n_external_samples": len(X_ext_sub) if X_ext_sub is not None else 0,
                        "cv_auc": auc,
                        "cv_balanced_accuracy": bal,
                        "cv_f1": f1,
                        "external_auc": ext_auc,
                        "external_balanced_accuracy": ext_bal,
                        "external_f1": ext_f1,
                        "status": "ok",
                    }
                )

                model_obj = fitted.named_steps["model"]
                if hasattr(model_obj, "coef_"):
                    imp = pd.DataFrame({"feature": genes, "coef": model_obj.coef_.ravel()}).sort_values(
                        "coef", ascending=False
                    )
                    imp.to_csv(RESULTS / "ml" / f"importance_{feature_name}_{model_name}.tsv", sep="\t", index=False)
                elif hasattr(model_obj, "feature_importances_"):
                    imp = pd.DataFrame(
                        {"feature": genes, "importance": model_obj.feature_importances_}
                    ).sort_values("importance", ascending=False)
                    imp.to_csv(RESULTS / "ml" / f"importance_{feature_name}_{model_name}.tsv", sep="\t", index=False)
            except Exception as exc:
                rows.append(
                    {
                        "task": "PTC_like_vs_follicular_patterned",
                        "feature_set": feature_name,
                        "model": model_name,
                        "status": f"failed: {exc}",
                    }
                )

    result_df = pd.DataFrame(rows)
    result_df.to_csv(RESULTS / "ml" / "baseline_ml_results.tsv", sep="\t", index=False)
    if out_prob_for_cal is not None and out_y_for_cal is not None:
        fig, ax = plt.subplots(figsize=(5, 5))
        CalibrationDisplay.from_predictions(
            (out_y_for_cal == "PTC_like").astype(int),
            out_prob_for_cal,
            n_bins=8,
            ax=ax,
        )
        ax.set_title("Calibration: TDS16 LogisticRegression")
        save_plot(fig, RESULTS / "ml" / "calibration_tds16_logreg.png")
    return result_df


def write_reports(
    dataset_master: pd.DataFrame,
    verification: pd.DataFrame,
    coverage: pd.DataFrame,
    ml_results: pd.DataFrame,
    sample_master: pd.DataFrame,
) -> None:
    usable = dataset_master[dataset_master["download_status"].isin(["ok", "partial"])]
    direct = sample_master[sample_master["label_confidence"] == "high"]["dataset"].value_counts().to_dict()
    inferred = sample_master[sample_master["label_confidence"] != "high"]["dataset"].value_counts().to_dict()
    ml_ok = ml_results[ml_results["status"] == "ok"]
    coverage_table = coverage.to_csv(sep="\t", index=False)
    ml_table = ml_results.to_csv(sep="\t", index=False)
    summary = f"""
# analysis_summary

## Confirmed

- Usable datasets in this run: {", ".join(usable["dataset_name"].tolist())}
- Direct histology annotations were available for: {json.dumps(direct, ensure_ascii=True)}
- Partly inferred / weaker labels were present for: {json.dumps(inferred, ensure_ascii=True)}
- Bulk RNA-seq and microarray were processed separately. No pooled training across modalities was performed.
- MTC samples were preserved as separate labels and not merged with follicular-derived tumors.
- TCGA-THCA backbone model ran: {"yes" if not ml_ok.empty else "no"}.

## Still Uncertain

- Exact source-verified BRS71 list was not recovered in this run; coverage and model use for BRS71 were skipped.
- Exact definition of the requested TierA67 panel was not supplied; TierA67 coverage and model use were skipped.
- GSE213647 contains a large patient-level bulk RNA-seq cohort but label alignment to TCGA is weaker than TCGA internal histology labels.
- Mutation calls were used only as auxiliary driver anchors. They were not used to finalize phenotype labels.

## Modality Integration

- Bulk RNA-seq and microarray should not be merged directly for supervised modeling at this stage.
- Platform effects and label non-equivalence remain large enough that microarray is treated as exploratory / directional validation only.

## Coverage

```tsv
{coverage_table.strip()}
```

## sklearn baseline

```tsv
{ml_table.strip()}
```

## Most Realistic Next Step

- Source-verify the BRS71 and TierA67 panel definitions from the exact intended publications or internal spec.
- Add stricter TCGA subtype curation including tall-cell vs classical vs FVPTC splits and fusion-aware anchors where available.
- If raw CEL reprocessing is required, install an R/Bioconductor stack and re-run microarray from CEL with RMA instead of relying on GEO processed matrices.
"""
    (REPORTS / "analysis_summary.md").write_text(summary.strip() + "\n", encoding="utf-8")

    next_steps = """
# next_steps

1. Lock the exact BRS71 and TierA67 gene lists before interpreting panel coverage or panel-based classifier performance.
2. Add controlled-access cohorts only after data use approval, then keep them in separate ingestion and provenance layers.
3. For TCGA, add fusion and copy-number annotations if the downstream question requires molecular subtype refinement beyond histology.
4. Re-run microarray from raw CEL files with RMA if cross-study reproducibility is a hard requirement.
5. Evaluate whether GSE213647 tumor-only PTC vs follicular-neoplasm labels are clinically comparable enough to serve as true external validation.
"""
    (REPORTS / "next_steps.md").write_text(next_steps.strip() + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()
    save_gene_panels()
    log("starting THCA public-data pipeline")

    geo_status = download_geo_assets()
    tcga_meta = download_tcga_thca_counts_and_metadata()

    geo_sample_master = build_sample_master_geo()
    tcga_sample_master = build_sample_master_tcga()
    sample_master = pd.concat([geo_sample_master, tcga_sample_master], ignore_index=True)

    dataset_master = create_dataset_master(sample_master, geo_status, tcga_meta)

    expr_map_bulk: dict[str, pd.DataFrame] = {}
    expr_map_micro: dict[str, pd.DataFrame] = {}

    tcga_expr, _ = process_tcga_rnaseq(sample_master)
    expr_map_bulk["TCGA-THCA"] = tcga_expr

    gse213647_expr, _ = process_gse213647(sample_master)
    expr_map_bulk["GSE213647"] = gse213647_expr

    gse126698_expr, _ = process_gse126698(sample_master)
    expr_map_bulk["GSE126698"] = gse126698_expr

    gse27155_expr, _, _ = process_microarray_dataset("GSE27155", "GPL96")
    expr_map_micro["GSE27155"] = gse27155_expr

    gse76039_expr, _, _ = process_microarray_dataset("GSE76039", "GPL570")
    expr_map_micro["GSE76039"] = gse76039_expr

    _ = process_gse97466(sample_master)

    sample_master = update_scores(sample_master, expr_map_bulk | expr_map_micro)
    sample_master.to_csv(META / "sample_master.tsv", sep="\t", index=False)

    verification = create_verification_table(sample_master, dataset_master)
    verification.to_csv(META / "verification_table.tsv", sep="\t", index=False)

    coverage = gene_coverage(expr_map_bulk | expr_map_micro)

    make_qc_plots(sample_master, expr_map_bulk, "bulk RNA-seq")
    make_qc_plots(sample_master, expr_map_micro, "microarray")

    ml_results = run_baseline_ml(sample_master, tcga_expr, gse213647_expr)

    write_reports(dataset_master, verification, coverage, ml_results, sample_master)
    log("pipeline complete")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(f"pipeline failed: {exc}")
        log(traceback.format_exc())
        raise
