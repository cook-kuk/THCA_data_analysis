from __future__ import annotations

import gzip
import json
import logging
import math
import os
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import yaml


def pilot_root() -> Path:
    return Path(__file__).resolve().parents[1]


def parent_root() -> Path:
    return pilot_root().parents[0]


def ensure_standard_dirs(root: Path | None = None) -> None:
    root = root or pilot_root()
    dirs = [
        "config",
        "scripts",
        "notebooks",
        "data_raw/tcga_thca",
        "data_raw/spatial",
        "data_raw/scrna",
        "data_raw/depmap_prism",
        "data_processed/tcga_thca",
        "data_processed/spatial",
        "data_processed/scrna",
        "data_processed/drug",
        "results/tables",
        "results/figures/spatial_maps",
        "results/reports",
        "results/logs",
    ]
    for d in dirs:
        (root / d).mkdir(parents=True, exist_ok=True)


def load_yaml(path: Path) -> dict:
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def config() -> dict:
    return load_yaml(pilot_root() / "config" / "config.yaml")


def gene_sets() -> dict:
    gs = load_yaml(pilot_root() / "config" / "gene_sets.yaml")
    return populate_dm1_gene_set(gs)


def setup_logging(name: str) -> logging.Logger:
    ensure_standard_dirs()
    logger = logging.getLogger(name)
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    stream = logging.StreamHandler(sys.stdout)
    stream.setFormatter(fmt)
    logger.addHandler(stream)
    log_path = pilot_root() / "results" / "logs" / f"{name}.log"
    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(fmt)
    logger.addHandler(file_handler)
    return logger


def warn(message: str) -> None:
    ensure_standard_dirs()
    path = pilot_root() / "results" / "logs" / "warnings.log"
    with path.open("a") as handle:
        handle.write(message.rstrip() + "\n")


def append_missing_data(section: str, message: str) -> None:
    path = pilot_root() / "MISSING_DATA.md"
    with path.open("a") as handle:
        handle.write(f"\n## {section}\n\n{message.strip()}\n")


def maybe_symlink(src: Path, dst: Path) -> Path:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        return dst
    try:
        dst.symlink_to(src)
    except OSError:
        shutil.copy2(src, dst)
    return dst


def clean_gene_symbol(value: object) -> str:
    s = str(value).strip()
    if not s:
        return s
    if " (" in s:
        s = s.split(" (", 1)[0]
    if "|" in s:
        parts = [p for p in s.split("|") if p and not p.upper().startswith("ENSG")]
        if parts:
            s = parts[-1]
    return s.upper()


def zscore_series(s: pd.Series) -> pd.Series:
    vals = pd.to_numeric(s, errors="coerce")
    mu = vals.mean()
    sd = vals.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(0.0, index=s.index)
    return (vals - mu) / sd


def zscore_df(df: pd.DataFrame) -> pd.DataFrame:
    return df.apply(zscore_series, axis=0)


def collapse_duplicate_genes(expr: pd.DataFrame) -> pd.DataFrame:
    """Collapse sample/cell/spot x gene expression by retaining highest mean duplicate."""
    expr = expr.copy()
    expr.columns = [clean_gene_symbol(c) for c in expr.columns]
    dupes = expr.columns[expr.columns.duplicated()].unique()
    if len(dupes) == 0:
        return expr.loc[:, [c for c in expr.columns if c]]
    pieces = []
    seen = set()
    for gene in expr.columns:
        if gene in seen or not gene:
            continue
        block = expr.loc[:, expr.columns == gene]
        if block.shape[1] == 1:
            pieces.append(block.iloc[:, 0].rename(gene))
        else:
            means = block.apply(pd.to_numeric, errors="coerce").mean(axis=0)
            best = int(np.nanargmax(means.to_numpy()))
            pieces.append(block.iloc[:, best].rename(gene))
        seen.add(gene)
    return pd.concat(pieces, axis=1)


def score_gene_set(
    expression_matrix: pd.DataFrame,
    genes: Iterable[str],
    method: str = "mean_z",
) -> tuple[pd.Series, dict]:
    genes = [clean_gene_symbol(g) for g in genes]
    available = [g for g in genes if g in expression_matrix.columns]
    meta = {"n_requested": len(genes), "n_found": len(available), "found_genes": ",".join(available)}
    if not available:
        return pd.Series(np.nan, index=expression_matrix.index), meta
    if method != "mean_z":
        warn(f"Unsupported method {method}; falling back to mean_z.")
    block = expression_matrix.loc[:, available].apply(pd.to_numeric, errors="coerce")
    scored = zscore_df(block).mean(axis=1)
    return scored, meta


def compute_module_scores(expr: pd.DataFrame, gs: dict, method: str = "mean_z") -> tuple[pd.DataFrame, pd.DataFrame]:
    score_df = pd.DataFrame(index=expr.index)
    coverage = []
    for name, payload in gs.items():
        short = payload.get("short_name", name)
        genes = payload.get("core", []) or []
        if short == "dm1_dark_lineage" and not genes:
            continue
        score, meta = score_gene_set(expr, genes, method)
        col = f"{short}_score"
        score_df[col] = score
        coverage.append({"gene_set": name, "score": col, **meta})
        if short == "vascular_delivery_proxy":
            neg = payload.get("barrier_negative", []) or []
            if neg:
                neg_score, neg_meta = score_gene_set(expr, neg, method)
                score_df["vascular_barrier_negative_score"] = neg_score
                coverage.append({"gene_set": f"{name}_barrier_negative", "score": "vascular_barrier_negative_score", **neg_meta})
    return score_df, pd.DataFrame(coverage)


def add_derived_axes(scores: pd.DataFrame, prefix: str = "") -> pd.DataFrame:
    out = scores.copy()
    def col(name: str) -> pd.Series:
        c = f"{name}_score"
        if c in out:
            return out[c]
        return pd.Series(0.0, index=out.index)

    rai = col("rai_differentiation")
    mapk = col("mapk")
    hla = col("hla_i_apm")
    ifn = col("ifn_activation")
    cyto = col("cytotoxic_t")
    epi = col("tumor_epithelial")
    caf = col("caf_ecm_tgfb")
    myeloid = col("myeloid_tam")
    hypoxia = col("hypoxia")
    vascular = col("vascular_delivery_proxy")
    prolif = col("proliferation")

    base = f"{prefix}_" if prefix else ""
    out[f"{base}immune_visibility_score"] = zscore_series(hla) + 0.5 * zscore_series(ifn) + 0.5 * zscore_series(cyto)
    out[f"{base}hla_low_invisible_score"] = -zscore_series(hla) + zscore_series(epi)
    out[f"{base}cd8_exclusion_proxy"] = zscore_series(caf) + zscore_series(myeloid) - zscore_series(cyto)
    out[f"{base}rai_restorable_score"] = zscore_series(rai) - 0.5 * zscore_series(mapk)
    out[f"{base}drug_delivery_failure_proxy"] = zscore_series(caf) + zscore_series(hypoxia) - zscore_series(vascular)
    out[f"{base}aggressive_dedifferentiation_score"] = -zscore_series(rai) + zscore_series(prolif) + zscore_series(mapk)
    out[f"{base}myeloid_caf_barrier_score"] = zscore_series(caf) + zscore_series(myeloid)
    return out


def patient_id_from_barcode(sample: object) -> str:
    s = str(sample)
    return "-".join(s.split("-")[:3]) if s.startswith("TCGA-") else s


def sample_type_code(sample: object) -> str:
    s = str(sample)
    parts = s.split("-")
    if len(parts) >= 4:
        return parts[3][:2]
    return ""


def stage_from_sample_id(sample_id: str) -> str:
    token = sample_id.split("_", 1)[-1]
    if "_N-" in sample_id or token.startswith("N-"):
        return "normal"
    if "_PTC-" in sample_id or token.startswith("PTC-"):
        return "PTC"
    if "_LPTC-" in sample_id or token.startswith("LPTC-"):
        return "locally_advanced_PTC"
    if "_ATC-" in sample_id or token.startswith("ATC-"):
        return "ATC"
    return "unknown"


def finite_spearman(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    from scipy.stats import spearmanr

    frame = pd.concat([x, y], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    n = frame.shape[0]
    if n < 3 or frame.iloc[:, 0].nunique() < 2 or frame.iloc[:, 1].nunique() < 2:
        return np.nan, np.nan, n
    r, p = spearmanr(frame.iloc[:, 0], frame.iloc[:, 1])
    return float(r), float(p), int(n)


def fdr_bh(pvals: Iterable[float]) -> list[float]:
    vals = np.array([np.nan if p is None else p for p in pvals], dtype=float)
    out = np.full(vals.shape, np.nan, dtype=float)
    mask = np.isfinite(vals)
    if not mask.any():
        return out.tolist()
    order = np.argsort(vals[mask])
    ranked = vals[mask][order]
    n = len(ranked)
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    tmp = np.empty(n)
    tmp[order] = np.minimum(q, 1.0)
    out[mask] = tmp
    return out.tolist()


def savefig(fig, path: Path, dpi: int = 300) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    if path.suffix.lower() != ".pdf":
        fig.savefig(path.with_suffix(".pdf"), dpi=dpi, bbox_inches="tight")


def classify_file(path: Path) -> str:
    s = str(path).lower()
    name = path.name.lower()
    if "gse250521" in s and (name.endswith(".h5ad") or "visium" in name or "matrix.mtx" in name):
        return "spatial_h5ad" if name.endswith(".h5ad") else "spatial_visium"
    if name.endswith(".h5ad") and ("scrna" in s or "single" in s or "lu2023" in s or "pu" in s):
        return "scrna_h5ad"
    if "tcga" in s and ("clinical" in name or "phenotype" in name or "survival" in name):
        return "tcga_clinical"
    if "tcga" in s and ("maf" in name or "mutation" in name or "driver" in name):
        return "tcga_mutation"
    if "tcga" in s and ("expression" in name or "geneexp" in name or "hiseq" in name):
        return "tcga_expression"
    if "depmap" in s or "ccle" in s or "achilles" in s:
        return "depmap"
    if "prism" in s or "drug" in s:
        return "prism"
    return "unknown"


def populate_dm1_gene_set(gs: dict) -> dict:
    """Populate DM1/dark-lineage list only from local files with explicit gene evidence."""
    root = Path("/home/seungho/personal/THCA_data_analysis")
    candidates = [
        root / "project/results/v17p35/tables/AMP4_8gene_model_coefficients.tsv",
        root / "project/src/02_qc_score/score_8gene_dm1.py",
        root / "PAPER1_DM1_FULL_2026_05_04.md",
    ]
    genes: list[str] = []
    for path in candidates:
        if not path.exists():
            continue
        try:
            text = path.read_text(errors="ignore")
        except OSError:
            continue
        for token in re.findall(r"\b[A-Z0-9]{2,12}(?:-[A-Z0-9]+)?\b", text):
            if token in {"DM1", "TCGA", "THCA", "PTC", "ATC", "RAI", "HLA", "AUC", "ROC", "PDF"}:
                continue
            if re.fullmatch(r"[A-Z][A-Z0-9-]{1,11}", token):
                genes.append(token)
    # Conservative whitelist from previous local 8-gene artifacts if present; otherwise empty.
    whitelist = {"PAX8", "TG", "TPO", "SLC5A5", "KRT19", "EPCAM", "FN1", "POSTN", "COL1A1", "SPP1", "TIMP1", "TROP2", "TACSTD2"}
    found = [g for g in genes if g in whitelist]
    seen = []
    for g in found:
        if g not in seen:
            seen.append(g)
    if seen:
        gs = dict(gs)
        payload = dict(gs.get("dm1_dark_lineage", {}))
        payload["core"] = seen
        payload["interpretation"] = payload.get("interpretation", "") + " Populated from local previous DM/dark-lineage artifacts; treat as exploratory."
        gs["dm1_dark_lineage"] = payload
    return gs


def write_not_run_table(path: Path, reason: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"status": "not_run", "reason": reason}]).to_csv(path, sep="\t", index=False)


def read_table_flexible(path: Path, **kwargs) -> pd.DataFrame:
    if str(path).endswith(".gz"):
        return pd.read_csv(path, sep=kwargs.pop("sep", "\t"), compression="gzip", **kwargs)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path, **kwargs)
    return pd.read_csv(path, sep=kwargs.pop("sep", "\t"), **kwargs)
