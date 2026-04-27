"""v17 ULTIMATE U1B - download additional GEO thyroid cohorts and build a
7-cohort meta-analysis (TCGA + 6 GEO) with the 8-gene de-differentiation panel.

Headline target: GSE151180 / GSE151179 RAI-refractory vs RAI-avid validation
(true clinical ground truth). The series GSE151180 is the *miRNA* sub-series
(GPL21575); the *gene-expression* sibling is GSE151179 (Affymetrix Clariom D
GPL23159) - we use that for the 8-gene panel.

Cohort plan
-----------
ALREADY LOCAL:
  - TCGA-THCA (RNA-seq log2)         : v17 cluster_orig labels (DM1 vs DM2)
  - GSE126698 (RNA-seq log2)         : tumor/normal contrast as proxy
  - GSE213647 (RNA-seq log2)         : PTC vs ATC+PD (U1C result)
  - GSE27155  (microarray log2)      : PTC vs other histologies
  - GSE76039  (microarray log2)      : ATC vs PDTC (R3A result)

NEW (downloaded here):
  - GSE65144  (GPL570 microarray)    : 12 ATC vs 13 normal thyroid
  - GSE60542  (GPL570 microarray)    : PTC vs normal thyroid
  - GSE151179 (GPL23159, gene expr.) : 39 PTC, RAI-avid vs RAI-refractory  (HEADLINE)
  - GSE192683                         : SKIP - non-thyroid (Piaractus mesopotamicus)

Outputs
-------
  results/v17_ultimate/U1B_per_cohort_real_aucs.tsv
  results/v17_ultimate/U1B_meta_analysis_real.tsv
  results/v17_ultimate/U1B_RAI_clinical_validation.tsv  (if GSE151179 succeeded)
  results/v17_ultimate/U1B_summary.json
  results/v17_ultimate/figures/U1B_forest_plot_7cohort.png
"""
from __future__ import annotations

import gzip
import io
import json
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

SCR = Path("/home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts")
sys.path.insert(0, str(SCR))
from v17_ULTIMATE_common import (  # noqa: E402
    DATA_PROC,
    DATA_RAW,
    GENE_8,
    RES,
    ROOT,
    jdump,
    load_tcga_expr,
    log,
)

OUT_DIR = RES
FIG_DIR = RES / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

GEO_DIR = DATA_RAW / "geo"
GEO_DIR.mkdir(parents=True, exist_ok=True)

# NKX2-1 alias map (PAX8 / TG etc are unique already)
GENE_ALIASES: Dict[str, List[str]] = {
    "NKX2-1": ["NKX2-1", "TITF1", "NKX2.1", "NKX2_1"],
}

# Ensembl gene IDs for the 8-gene panel (canonical, no version suffix)
ENSG_FOR_GENE_8: Dict[str, str] = {
    "SLC5A5":  "ENSG00000105641",
    "TPO":     "ENSG00000115705",
    "TG":      "ENSG00000042832",
    "TSHR":    "ENSG00000165409",
    "PAX8":    "ENSG00000125618",
    "NKX2-1":  "ENSG00000136352",
    "FOXE1":   "ENSG00000178919",
    "DIO1":    "ENSG00000211448",
}


def _maybe_log2(df: pd.DataFrame, label: str) -> pd.DataFrame:
    """If max > 50 (raw intensities), apply log2(x+1). Returns transformed df."""
    mx = float(np.nanmax(df.values))
    if mx > 50:
        log(f"  {label}: data appears un-logged (max={mx:.1f}); applying log2(x+1)")
        return np.log2(df.clip(lower=0) + 1)
    return df


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def fetch_bytes(url: str, timeout: int = 600) -> bytes:
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=timeout) as r:
        return r.read()


def cached_download(url: str, dest: Path) -> Optional[Path]:
    if dest.exists() and dest.stat().st_size > 1024:
        log(f"  cache hit: {dest.name}")
        return dest
    try:
        log(f"  downloading {url}")
        data = fetch_bytes(url, timeout=600)
        dest.write_bytes(data)
        log(f"  wrote {dest} ({len(data)/1e6:.1f} MB)")
        return dest
    except Exception as e:
        log(f"  download FAILED: {e}")
        return None


def parse_series_matrix(path: Path) -> Tuple[pd.DataFrame, Dict[str, List[str]]]:
    """Parse a series_matrix.txt.gz - returns (probe x sample dataframe, meta dict)."""
    with gzip.open(path, "rb") as fh:
        txt = fh.read().decode("utf-8", errors="replace")
    lines = txt.split("\n")
    meta: Dict[str, List[str]] = {}
    gsms: List[str] = []
    chars: List[List[str]] = []
    titles: List[str] = []
    sources: List[str] = []
    matrix_start = matrix_end = None
    for i, ln in enumerate(lines):
        if ln.startswith("!Sample_geo_accession"):
            gsms = [v.strip().strip('"') for v in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_title"):
            titles = [v.strip().strip('"') for v in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_source_name_ch1"):
            sources = [v.strip().strip('"') for v in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1"):
            chars.append([v.strip().strip('"') for v in ln.split("\t")[1:]])
        elif ln.strip() == "!series_matrix_table_begin":
            matrix_start = i + 1
        elif ln.strip() == "!series_matrix_table_end":
            matrix_end = i
    meta["gsms"] = gsms
    meta["titles"] = titles
    meta["sources"] = sources
    meta["chars"] = chars
    if matrix_start is None or matrix_end is None:
        return pd.DataFrame(), meta
    matrix_text = "\n".join(lines[matrix_start:matrix_end])
    df = pd.read_csv(io.StringIO(matrix_text), sep="\t", index_col=0)
    return df, meta


# -----------------------------------------------------------------------------
# Platform annotation parsers
# -----------------------------------------------------------------------------
def load_gpl570_annotation() -> pd.DataFrame:
    """Return Affymetrix HG-U133 Plus 2 (GPL570) probe-to-symbol mapping."""
    p = GEO_DIR / "platforms" / "GPL570.annot.csv"
    if not p.exists():
        raise FileNotFoundError(p)
    # NCBI annot CSV: skip "" lines, comma-quoted; we read with pandas.
    df = pd.read_csv(p, comment="#", header=0, dtype=str)
    cols = [c for c in df.columns if c.lower() in ("id", "platform_clonid", "gene symbol")]
    # Look for gene symbol column
    sym_col = None
    for c in df.columns:
        cl = c.strip().lower()
        if cl in ("gene symbol", "symbol", "gene_symbol"):
            sym_col = c
            break
    if sym_col is None:
        # GEO .annot files use 'Gene symbol'
        for c in df.columns:
            if "symbol" in c.lower():
                sym_col = c
                break
    id_col = None
    for c in df.columns:
        if c.strip().upper() == "ID" or c.strip() == "Platform_CLONEID" or c.strip() == "ID_REF":
            id_col = c
            break
    if id_col is None:
        id_col = df.columns[0]
    out = df[[id_col, sym_col]].rename(columns={id_col: "probe_id", sym_col: "symbol"})
    out["probe_id"] = out["probe_id"].astype(str)
    out["symbol"] = out["symbol"].astype(str).str.split("///").str[0].str.strip()
    out = out[out["symbol"].notna() & (out["symbol"] != "")]
    return out


def parse_gpl23159_for_panel(panel: List[str], txt_path: Path) -> Dict[str, List[str]]:
    """For GPL23159 platform table on disk, return {gene: [probe_ids]} for each panel gene."""
    txt = txt_path.read_text(errors="replace")
    lines = txt.split("\n")
    # Build alias-aware target list
    targets: Dict[str, List[str]] = {g: [] for g in panel}
    for ln in lines:
        if not ln.startswith("TC") and not ln.startswith("TS"):
            continue
        # Cheap pre-filter
        for g in panel:
            aliases = GENE_ALIASES.get(g, [g])
            for a in aliases:
                # match `(SYMBOL)` or ` symbol ` boundaries
                if re.search(rf"\({re.escape(a)}\)", ln) or re.search(rf"\b{re.escape(a)}\b", ln):
                    pid = ln.split("\t", 1)[0]
                    targets[g].append(pid)
                    break
    # Deduplicate and prefer "main" probesets (suffix .hg.1)
    out: Dict[str, List[str]] = {}
    for g, plist in targets.items():
        uniq: List[str] = []
        for p in plist:
            if p not in uniq:
                uniq.append(p)
        # Prefer probes containing 'TC' (transcript cluster, main)
        main = [p for p in uniq if p.startswith("TC")]
        out[g] = main if main else uniq
    return out


def fetch_gpl_table(gpl: str, dest: Path) -> Optional[Path]:
    """Download platform table via NCBI viewer. Cached to dest."""
    if dest.exists() and dest.stat().st_size > 1024:
        log(f"  cache hit: {dest.name}")
        return dest
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&form=text&view=data&acc={gpl}"
    return cached_download(url, dest)


# -----------------------------------------------------------------------------
# TCGA training
# -----------------------------------------------------------------------------
def train_tcga_8gene(common: List[str]) -> Tuple[StandardScaler, LogisticRegression, float]:
    """Train LR on TCGA log2 expression, restricted to `common` (subset of GENE_8)."""
    tcga = load_tcga_expr()
    labels_path = ROOT / "project" / "results" / "v17_realfix" / "R1A_cluster_labels.tsv"
    lbl = pd.read_csv(labels_path, sep="\t")
    lbl["y"] = (lbl["cluster"] == "DM1_A").astype(int)
    X = tcga.loc[common].T
    X = X.loc[X.index.intersection(lbl["sample_id"])]
    y = lbl.set_index("sample_id").loc[X.index, "y"].values
    sc = StandardScaler().fit(X.values)
    m = LogisticRegression(max_iter=5000, C=1.0, class_weight="balanced", random_state=42)
    m.fit(sc.transform(X.values), y)
    auc = roc_auc_score(y, m.predict_proba(sc.transform(X.values))[:, 1])
    log(f"  TCGA train (n={len(y)}, n_genes={len(common)}): in-sample AUC={auc:.3f}")
    return sc, m, float(auc)


def bootstrap_auc(y_true: np.ndarray, y_score: np.ndarray, n: int = 1000, seed: int = 42
                  ) -> Tuple[float, float, float]:
    """Returns (point_auc, ci_lo, ci_hi). Direction-corrected (max(auc, 1-auc))."""
    rng = np.random.default_rng(seed)
    point = roc_auc_score(y_true, y_score)
    point = max(point, 1 - point)
    aucs = []
    for _ in range(n):
        idx = rng.integers(0, len(y_true), len(y_true))
        if len(np.unique(y_true[idx])) < 2:
            continue
        a = roc_auc_score(y_true[idx], y_score[idx])
        aucs.append(max(a, 1 - a))
    if not aucs:
        return point, np.nan, np.nan
    lo = float(np.percentile(aucs, 2.5))
    hi = float(np.percentile(aucs, 97.5))
    return float(point), lo, hi


# -----------------------------------------------------------------------------
# Per-cohort runners
# -----------------------------------------------------------------------------
def run_tcga(scaler, model, gene_order) -> Dict:
    """In-sample TCGA against the cluster labels. Lists for completeness."""
    tcga = load_tcga_expr()
    labels_path = ROOT / "project" / "results" / "v17_realfix" / "R1A_cluster_labels.tsv"
    lbl = pd.read_csv(labels_path, sep="\t")
    lbl["y"] = (lbl["cluster"] == "DM1_A").astype(int)
    X = tcga.loc[gene_order].T.loc[lbl["sample_id"]]
    proba = model.predict_proba(scaler.transform(X.values))[:, 1]
    y = lbl["y"].values
    auc, lo, hi = bootstrap_auc(y, proba, n=1000, seed=42)
    log(f"[TCGA] AUC={auc:.3f} [{lo:.3f}, {hi:.3f}] (n={len(y)}, ground_truth=R1A_cluster_labels)")
    return {
        "cohort": "TCGA-THCA",
        "n": int(len(y)),
        "n_genes_matched": int(len(gene_order)),
        "ground_truth_type": "v17_DM1_vs_DM2_cluster",
        "n_pos": int(y.sum()),
        "n_neg": int((1 - y).sum()),
        "auc": auc,
        "ci_lo": lo,
        "ci_hi": hi,
        "note": "in-sample (training cohort)",
    }


def run_gse76039(scaler, model, gene_order) -> Optional[Dict]:
    expr_path = DATA_PROC / "microarray" / "GSE76039_microarray_expression_log2.tsv"
    sm_path = DATA_RAW / "geo" / "GSE76039" / "GSE76039_series_matrix.txt.gz"
    if not expr_path.exists() or not sm_path.exists():
        log("[GSE76039] missing local files, skipping")
        return None
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    _, meta = parse_series_matrix(sm_path)
    src = meta.get("sources", [])
    gsms = meta.get("gsms", [])
    df_meta = pd.DataFrame({"sample_id": gsms, "source": src})
    df_meta["histology"] = df_meta["source"].astype(str).str.lower()
    df_meta["y_ATC"] = df_meta["histology"].str.contains("anaplastic", na=False).astype(int)
    df_meta["is_pdtc"] = df_meta["histology"].str.contains("poorly", na=False).astype(int)
    keep = (df_meta["y_ATC"] == 1) | (df_meta["is_pdtc"] == 1)
    df_keep = df_meta[keep].reset_index(drop=True)
    matched = [g for g in gene_order if g in expr.index]
    if len(matched) < 6:
        log(f"[GSE76039] only {len(matched)}/{len(gene_order)} genes matched, skipping")
        return None
    samples = [s for s in df_keep["sample_id"] if s in expr.columns]
    X = expr.loc[matched, samples].T
    Xfull = pd.DataFrame(0.0, index=X.index, columns=gene_order)
    for g in matched:
        Xfull[g] = X[g]
    proba = model.predict_proba(scaler.transform(Xfull.values))[:, 1]
    y = df_keep.set_index("sample_id").loc[samples, "y_ATC"].astype(int).values
    auc, lo, hi = bootstrap_auc(y, proba, n=1000, seed=42)
    log(f"[GSE76039] AUC={auc:.3f} [{lo:.3f}, {hi:.3f}] (n={len(y)}, ATC vs PDTC)")
    return {
        "cohort": "GSE76039",
        "n": int(len(y)),
        "n_genes_matched": int(len(matched)),
        "ground_truth_type": "ATC_vs_PDTC_histology",
        "n_pos": int(y.sum()),
        "n_neg": int((1 - y).sum()),
        "auc": auc,
        "ci_lo": lo,
        "ci_hi": hi,
        "note": "external; histology = ground truth",
    }


def run_gse27155(scaler, model, gene_order) -> Optional[Dict]:
    expr_path = DATA_PROC / "microarray" / "GSE27155_microarray_expression_log2.tsv"
    sm_path = DATA_RAW / "geo" / "GSE27155" / "GSE27155_series_matrix.txt.gz"
    if not expr_path.exists():
        log("[GSE27155] missing expression, skipping")
        return None
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    if not sm_path.exists():
        soft_path = DATA_RAW / "geo" / "GSE27155_family.soft.gz"
        if not soft_path.exists():
            log("[GSE27155] no series matrix or SOFT, skipping")
            return None
        # Parse SOFT for source_name
        with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as fh:
            txt = fh.read()
        rows = []
        cur = None
        for ln in txt.split("\n"):
            if ln.startswith("^SAMPLE"):
                if cur is not None:
                    rows.append(cur)
                cur = {"sample_id": ln.split("=", 1)[1].strip(), "source": "", "title": ""}
            elif cur is not None:
                if ln.startswith("!Sample_source_name_ch1"):
                    cur["source"] = ln.split("=", 1)[1].strip()
                elif ln.startswith("!Sample_title"):
                    cur["title"] = ln.split("=", 1)[1].strip()
        if cur is not None:
            rows.append(cur)
        df_meta = pd.DataFrame(rows)
    else:
        _, meta = parse_series_matrix(sm_path)
        df_meta = pd.DataFrame({"sample_id": meta["gsms"], "source": meta["sources"]})

    df_meta["src_lower"] = df_meta["source"].astype(str).str.lower()
    # PTC vs ATC + Hurthle (RAI-resistant subset) - we'll use papillary vs other-malignant
    df_meta["is_papillary"] = df_meta["src_lower"].str.contains(r"papillary", na=False).astype(int)
    df_meta["is_anaplastic"] = df_meta["src_lower"].str.contains(r"anapla", na=False).astype(int)
    df_meta["is_normal"] = df_meta["src_lower"].str.contains(r"normal", na=False).astype(int)
    df_meta["is_follicular"] = df_meta["src_lower"].str.contains(r"follicular", na=False).astype(int)
    keep = (df_meta["is_papillary"] == 1) | (df_meta["is_anaplastic"] == 1)
    if keep.sum() < 6:
        # Try follicular (differentiated) vs anaplastic
        keep = (df_meta["is_papillary"] == 1) | (df_meta["is_anaplastic"] == 1) | (
            df_meta["is_follicular"] == 1
        )
        df_keep = df_meta[keep].copy()
        df_keep["y"] = df_keep["is_anaplastic"]
    else:
        df_keep = df_meta[keep].copy()
        df_keep["y"] = df_keep["is_anaplastic"]

    matched = [g for g in gene_order if g in expr.index]
    if len(matched) < 6:
        log(f"[GSE27155] only {len(matched)}/{len(gene_order)} matched, skipping")
        return None
    samples = [s for s in df_keep["sample_id"] if s in expr.columns]
    if len(samples) < 6 or df_keep.set_index("sample_id").loc[samples, "y"].nunique() < 2:
        log("[GSE27155] insufficient labelled samples, skipping")
        return None
    X = expr.loc[matched, samples].T
    Xfull = pd.DataFrame(0.0, index=X.index, columns=gene_order)
    for g in matched:
        Xfull[g] = X[g]
    proba = model.predict_proba(scaler.transform(Xfull.values))[:, 1]
    y = df_keep.set_index("sample_id").loc[samples, "y"].astype(int).values
    auc, lo, hi = bootstrap_auc(y, proba, n=1000, seed=42)
    log(f"[GSE27155] AUC={auc:.3f} [{lo:.3f}, {hi:.3f}] (n={len(y)}, ATC vs PTC/FTC)")
    return {
        "cohort": "GSE27155",
        "n": int(len(y)),
        "n_genes_matched": int(len(matched)),
        "ground_truth_type": "ATC_vs_PTC_FTC",
        "n_pos": int(y.sum()),
        "n_neg": int((1 - y).sum()),
        "auc": auc,
        "ci_lo": lo,
        "ci_hi": hi,
        "note": "external; histology = ground truth",
    }


def run_gse213647(scaler, model, gene_order) -> Optional[Dict]:
    """Use ComBat-corrected v3 if available; AUC PTC vs ATC+PD on tumors."""
    p_v3 = DATA_PROC / "bulk_rnaseq_v3" / "GSE213647_v3_log2.tsv"
    p_old = DATA_PROC / "bulk_rnaseq" / "GSE213647_rnaseq_expression_log2.tsv"
    expr_path = p_v3 if p_v3.exists() else p_old
    soft_path = DATA_RAW / "geo" / "GSE213647_family.soft.gz"
    if not expr_path.exists() or not soft_path.exists():
        log("[GSE213647] missing files, skipping")
        return None
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    # Quickly parse SOFT
    rows = []
    cur = None
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln.startswith("^SAMPLE"):
                if cur:
                    rows.append(cur)
                cur = {"sample_id": ln.split("=", 1)[1].strip(),
                       "cell_type": "", "cell_subtype": ""}
            elif cur is not None and ln.startswith("!Sample_characteristics_ch1"):
                v = ln.split("=", 1)[1].strip()
                m = re.match(r"^([^:]+):\s*(.*)$", v)
                if m:
                    k = m.group(1).strip().lower().replace(" ", "_")
                    if k in cur:
                        cur[k] = m.group(2).strip()
    if cur:
        rows.append(cur)
    df_meta = pd.DataFrame(rows)
    df_meta["histology"] = df_meta["cell_subtype"]
    tumors = df_meta[df_meta["cell_type"].str.lower() == "tumor"]
    diff = {"PTC"}
    dediff = {"ATC", "PD"}
    sub = tumors[tumors["histology"].isin(diff | dediff)].copy()
    # Index may be ENSG (v3) or symbol; build a sym→row mapping
    expr.index = expr.index.astype(str)
    is_ensg = expr.index.str.startswith("ENSG").mean() > 0.5
    if is_ensg:
        ensg2gene = {ENSG_FOR_GENE_8[g]: g for g in gene_order if g in ENSG_FOR_GENE_8}
        # Strip versions if present
        idx_clean = expr.index.str.split(".").str[0]
        expr.index = idx_clean
        matched = [g for g in gene_order if ENSG_FOR_GENE_8.get(g) in expr.index]
        log(f"[GSE213647] (ENSG-indexed) symbol matches via Ensembl: {len(matched)}/{len(gene_order)}")
    else:
        matched = [g for g in gene_order if g in expr.index]
    if len(matched) < 6:
        log(f"[GSE213647] only {len(matched)}/{len(gene_order)} matched, skipping")
        return None
    samples = [s for s in sub["sample_id"] if s in expr.columns]
    if len(samples) < 6:
        log("[GSE213647] insufficient samples")
        return None
    if is_ensg:
        rows_for_panel = [ENSG_FOR_GENE_8[g] for g in matched]
        X = expr.loc[rows_for_panel, samples].T
        X.columns = matched
    else:
        X = expr.loc[matched, samples].T
    Xfull = pd.DataFrame(0.0, index=X.index, columns=gene_order)
    for g in matched:
        Xfull[g] = X[g]
    proba = model.predict_proba(scaler.transform(Xfull.values))[:, 1]
    y = sub.set_index("sample_id").loc[samples, "histology"].isin(dediff).astype(int).values
    if len(np.unique(y)) < 2:
        log("[GSE213647] single class, skipping")
        return None
    auc, lo, hi = bootstrap_auc(y, proba, n=1000, seed=42)
    log(f"[GSE213647] AUC={auc:.3f} [{lo:.3f}, {hi:.3f}] (n={len(y)}, PTC vs ATC+PD)")
    return {
        "cohort": "GSE213647",
        "n": int(len(y)),
        "n_genes_matched": int(len(matched)),
        "ground_truth_type": "PTC_vs_ATC+PD_histology",
        "n_pos": int(y.sum()),
        "n_neg": int((1 - y).sum()),
        "auc": auc,
        "ci_lo": lo,
        "ci_hi": hi,
        "note": "external; histology = ground truth (largest)",
    }


def run_gse126698(scaler, model, gene_order) -> Optional[Dict]:
    """RNA-seq cohort - tumor vs normal proxy (no histology subtype labels)."""
    expr_path = DATA_PROC / "bulk_rnaseq" / "GSE126698_rnaseq_expression_log2.tsv"
    soft_path = DATA_RAW / "geo" / "GSE126698_family.soft.gz"
    if not expr_path.exists() or not soft_path.exists():
        log("[GSE126698] missing files, skipping")
        return None
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    # Map columns (A1..A20 etc.) to GSMs via SOFT title field
    rows = []
    cur = None
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln.startswith("^SAMPLE"):
                if cur:
                    rows.append(cur)
                cur = {"sample_id": ln.split("=", 1)[1].strip(),
                       "title": "", "source": "", "tissue": "", "cell_type": ""}
            elif cur is not None:
                if ln.startswith("!Sample_title"):
                    cur["title"] = ln.split("=", 1)[1].strip()
                elif ln.startswith("!Sample_source_name_ch1"):
                    cur["source"] = ln.split("=", 1)[1].strip()
                elif ln.startswith("!Sample_characteristics_ch1"):
                    v = ln.split("=", 1)[1].strip()
                    m = re.match(r"^([^:]+):\s*(.*)$", v)
                    if m:
                        k = m.group(1).strip().lower().replace(" ", "_")
                        if k in cur:
                            cur[k] = m.group(2).strip()
        if cur:
            rows.append(cur)
    df_meta = pd.DataFrame(rows)
    # Many GSE126698 metadata fields ended up as "char_*" prefix - look for a
    # tumor-subtype char regardless of column naming
    subtype_col = None
    for c in df_meta.columns:
        if "subtype" in c.lower() or "tumor_subtype" in c.lower():
            subtype_col = c
            break
    if subtype_col is None:
        # Re-parse with broader char keys
        rows2: List[Dict[str, str]] = []
        cur2: Optional[Dict[str, str]] = None
        with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as fh:
            for ln in fh:
                ln = ln.rstrip("\n")
                if ln.startswith("^SAMPLE"):
                    if cur2:
                        rows2.append(cur2)
                    cur2 = {"sample_id": ln.split("=", 1)[1].strip(),
                            "title": "", "source": "", "subtype": ""}
                elif cur2 is not None:
                    if ln.startswith("!Sample_title"):
                        cur2["title"] = ln.split("=", 1)[1].strip()
                    elif ln.startswith("!Sample_source_name_ch1"):
                        cur2["source"] = ln.split("=", 1)[1].strip()
                    elif ln.startswith("!Sample_characteristics_ch1"):
                        v = ln.split("=", 1)[1].strip()
                        m = re.match(r"^([^:]+):\s*(.*)$", v)
                        if m and ("subtype" in m.group(1).lower()
                                  or "tumor" in m.group(1).lower()):
                            cur2["subtype"] = m.group(2).strip()
            if cur2:
                rows2.append(cur2)
        df_meta = pd.DataFrame(rows2)
        subtype_col = "subtype"

    # Build column-name -> sample mapping. Columns are like "A3" or "F1"; titles like "A3 [RNA-Seq]"
    df_meta["title_short"] = df_meta["title"].astype(str).str.split(" ").str[0]
    title2sub = dict(zip(df_meta["title_short"], df_meta[subtype_col].astype(str).str.upper()))
    # In GSE126698 NT = non-transformed thyroid (normal); ATC/PTC/FTC = tumor
    df_meta["y"] = df_meta[subtype_col].astype(str).str.upper().map(
        lambda s: 0 if s == "NT" else (1 if s in ("ATC", "PTC", "FTC", "PDTC") else -1)
    )
    short2y = dict(zip(df_meta["title_short"], df_meta["y"]))

    expr_cols = list(expr.columns)
    keep_cols = [c for c in expr_cols if c in short2y and short2y[c] != -1]
    final_cols = []
    y_arr = []
    for c in keep_cols:
        y_arr.append(int(short2y[c]))
        final_cols.append(c)
    if len(set(y_arr)) < 2 or len(final_cols) < 6:
        log(f"[GSE126698] insufficient labels (n={len(final_cols)}), skipping")
        return None
    expr = _maybe_log2(expr, "GSE126698")
    matched = [g for g in gene_order if g in expr.index]
    if len(matched) < 6:
        log(f"[GSE126698] {len(matched)}/{len(gene_order)} matched, skipping")
        return None
    # Prefer ATC vs PTC/FTC contrast if available (de-differentiation)
    atc_cols = [c for c in final_cols if title2sub.get(c, "").upper() == "ATC"]
    ptc_ftc_cols = [c for c in final_cols
                    if title2sub.get(c, "").upper() in ("PTC", "FTC")]
    if len(atc_cols) >= 3 and len(ptc_ftc_cols) >= 3:
        contrast_cols = atc_cols + ptc_ftc_cols
        y = np.array([1] * len(atc_cols) + [0] * len(ptc_ftc_cols))
        gtype = "ATC_vs_PTC+FTC"
    else:
        contrast_cols = final_cols
        y = np.array(y_arr)
        gtype = "tumor_vs_normal"
    X = expr.loc[matched, contrast_cols].T
    Xfull = pd.DataFrame(0.0, index=X.index, columns=gene_order)
    for g in matched:
        Xfull[g] = X[g]
    proba = model.predict_proba(scaler.transform(Xfull.values))[:, 1]
    auc, lo, hi = bootstrap_auc(y, proba, n=1000, seed=42)
    log(f"[GSE126698] AUC={auc:.3f} [{lo:.3f}, {hi:.3f}] (n={len(y)}, {gtype})")
    return {
        "cohort": "GSE126698",
        "n": int(len(y)),
        "n_genes_matched": int(len(matched)),
        "ground_truth_type": gtype,
        "n_pos": int(y.sum()),
        "n_neg": int((1 - y).sum()),
        "auc": auc,
        "ci_lo": lo,
        "ci_hi": hi,
        "note": "external; histology = ground truth",
    }


def _aggregate_microarray_to_genes(expr_df: pd.DataFrame, probe2sym: pd.DataFrame,
                                   panel: List[str]) -> pd.DataFrame:
    """Map probe rows -> gene rows by mean (filter probes mapping to panel genes)."""
    sym = probe2sym.copy()
    # Resolve aliases
    rev_alias: Dict[str, str] = {}
    for canonical, alts in GENE_ALIASES.items():
        for a in alts:
            rev_alias[a.upper()] = canonical
    sym["symbol_norm"] = sym["symbol"].astype(str).str.upper().map(
        lambda s: rev_alias.get(s, s)
    )
    panel_up = {g.upper(): g for g in panel}
    sym = sym[sym["symbol_norm"].isin(panel_up)]
    # probe IDs in the expression matrix
    expr_df.index = expr_df.index.astype(str)
    sym["probe_id"] = sym["probe_id"].astype(str)
    sub = sym[sym["probe_id"].isin(expr_df.index)]
    if sub.empty:
        return pd.DataFrame()
    out_rows = {}
    for canonical_up, grp in sub.groupby("symbol_norm"):
        probes = grp["probe_id"].tolist()
        out_rows[panel_up[canonical_up]] = expr_df.loc[probes].mean(axis=0)
    return pd.DataFrame(out_rows).T  # genes x samples


def run_gse65144(scaler, model, gene_order) -> Optional[Dict]:
    """ATC vs Normal microarray (GPL570). Strong differentiation contrast."""
    sm_path = GEO_DIR / "GSE65144_series_matrix.txt.gz"
    if not sm_path.exists():
        log("[GSE65144] no series matrix, skipping")
        return None
    expr_df, meta = parse_series_matrix(sm_path)
    if expr_df.empty:
        log("[GSE65144] empty matrix")
        return None
    log(f"[GSE65144] raw matrix {expr_df.shape}")
    expr_df = _maybe_log2(expr_df, "GSE65144")
    probe2sym = load_gpl570_annotation()
    gx = _aggregate_microarray_to_genes(expr_df, probe2sym, gene_order)
    matched = [g for g in gene_order if g in gx.index]
    log(f"[GSE65144] genes matched: {len(matched)}/{len(gene_order)} ({matched})")
    if len(matched) < 6:
        return None
    sources = pd.Series(meta["sources"], index=meta["gsms"])
    src_lower = sources.str.lower()
    y_atc = src_lower.str.contains(r"anapla", na=False).astype(int)
    is_normal = src_lower.str.contains(r"normal", na=False).astype(int)
    keep = (y_atc == 1) | (is_normal == 1)
    samples = [s for s in y_atc.index[keep] if s in gx.columns]
    if len(samples) < 6:
        log(f"[GSE65144] only {len(samples)} labelled, skipping")
        return None
    X = gx.loc[matched, samples].T
    Xfull = pd.DataFrame(0.0, index=X.index, columns=gene_order)
    for g in matched:
        Xfull[g] = X[g]
    proba = model.predict_proba(scaler.transform(Xfull.values))[:, 1]
    y = y_atc.loc[samples].values
    auc, lo, hi = bootstrap_auc(y, proba, n=1000, seed=42)
    log(f"[GSE65144] AUC={auc:.3f} [{lo:.3f}, {hi:.3f}] (n={len(y)}, ATC vs normal)")
    return {
        "cohort": "GSE65144",
        "n": int(len(y)),
        "n_genes_matched": int(len(matched)),
        "ground_truth_type": "ATC_vs_normal_thyroid",
        "n_pos": int(y.sum()),
        "n_neg": int((1 - y).sum()),
        "auc": auc,
        "ci_lo": lo,
        "ci_hi": hi,
        "note": "external; histology = ground truth (12 ATC vs 13 normal)",
    }


def run_gse60542(scaler, model, gene_order) -> Optional[Dict]:
    """PTC vs normal thyroid - GPL570 microarray. Skip lymph-node metastases / NLN."""
    sm_path = GEO_DIR / "GSE60542_series_matrix.txt.gz"
    if not sm_path.exists():
        log("[GSE60542] no series matrix, skipping")
        return None
    expr_df, meta = parse_series_matrix(sm_path)
    if expr_df.empty:
        log("[GSE60542] empty matrix")
        return None
    log(f"[GSE60542] raw matrix {expr_df.shape}")
    expr_df = _maybe_log2(expr_df, "GSE60542")
    probe2sym = load_gpl570_annotation()
    gx = _aggregate_microarray_to_genes(expr_df, probe2sym, gene_order)
    matched = [g for g in gene_order if g in gx.index]
    log(f"[GSE60542] genes matched: {len(matched)}/{len(gene_order)} ({matched})")
    if len(matched) < 6:
        return None
    titles = pd.Series(meta["titles"], index=meta["gsms"]).str.lower()
    sources = pd.Series(meta["sources"], index=meta["gsms"]).str.lower()
    is_ptc = (titles.str.contains(r"papillary", na=False) | sources.str.contains(r"papillary", na=False)).astype(int)
    is_normal_thy = ((titles.str.contains(r"normal thyroid", na=False)
                      | sources.str.contains(r"normal thyroid", na=False))
                     & ~titles.str.contains(r"lymph node", na=False)).astype(int)
    keep = (is_ptc == 1) | (is_normal_thy == 1)
    samples = [s for s in is_ptc.index[keep] if s in gx.columns]
    if len(samples) < 6:
        log(f"[GSE60542] only {len(samples)} labelled, skipping")
        return None
    X = gx.loc[matched, samples].T
    Xfull = pd.DataFrame(0.0, index=X.index, columns=gene_order)
    for g in matched:
        Xfull[g] = X[g]
    proba = model.predict_proba(scaler.transform(Xfull.values))[:, 1]
    y = is_ptc.loc[samples].values
    if len(np.unique(y)) < 2:
        log("[GSE60542] single class, skipping")
        return None
    auc, lo, hi = bootstrap_auc(y, proba, n=1000, seed=42)
    log(f"[GSE60542] AUC={auc:.3f} [{lo:.3f}, {hi:.3f}] (n={len(y)}, PTC vs normal thyroid)")
    return {
        "cohort": "GSE60542",
        "n": int(len(y)),
        "n_genes_matched": int(len(matched)),
        "ground_truth_type": "PTC_vs_normal_thyroid",
        "n_pos": int(y.sum()),
        "n_neg": int((1 - y).sum()),
        "auc": auc,
        "ci_lo": lo,
        "ci_hi": hi,
        "note": "external; PTC primary tumor vs matched normal (LNM excluded)",
    }


def run_gse151179(scaler, model, gene_order) -> Optional[Dict]:
    """HEADLINE: RAI-refractory vs RAI-avid PTC (GPL23159 Clariom D)."""
    sm_path = GEO_DIR / "GSE151179_series_matrix.txt.gz"
    if not sm_path.exists():
        log("[GSE151179] no series matrix, skipping")
        return None
    expr_df, meta = parse_series_matrix(sm_path)
    if expr_df.empty:
        log("[GSE151179] empty matrix")
        return None
    log(f"[GSE151179] raw matrix {expr_df.shape}")
    # Find RAI response characteristic
    chars = meta["chars"]
    rai_row: Optional[List[str]] = None
    tissue_row: Optional[List[str]] = None
    for c in chars:
        joined = " | ".join(c[:5]).lower()
        if "patient rai responce" in joined or "rai response" in joined or (
            "avid" in joined and "refract" in joined
        ):
            rai_row = c
        if "tissue type" in joined and tissue_row is None:
            tissue_row = c
    if rai_row is None:
        log("[GSE151179] no RAI row found in characteristics")
        return None
    # Resolve labels
    y_rai = []
    for v in rai_row:
        v_low = v.lower()
        if "refract" in v_low or "rair" in v_low:
            y_rai.append(1)
        elif "avid" in v_low:
            y_rai.append(0)
        else:
            y_rai.append(-1)
    y_rai = np.array(y_rai)
    log(f"[GSE151179] RAI labels: refractory={int((y_rai==1).sum())}, "
        f"avid={int((y_rai==0).sum())}, other={int((y_rai==-1).sum())}")
    # Tissue type filter: keep TUMOR specimens (primary OR LNM); drop non-neoplastic thyroid
    if tissue_row is not None:
        tissue = np.array([t.lower() for t in tissue_row])
        is_tumor_specimen = np.array(["non-neoplastic" not in t for t in tissue])
        is_primary = np.array(["primary tumor" in t for t in tissue])
    else:
        is_tumor_specimen = np.ones(len(y_rai), dtype=bool)
        is_primary = np.ones(len(y_rai), dtype=bool)
    keep_mask = (y_rai >= 0) & is_tumor_specimen
    primary_mask = (y_rai >= 0) & is_primary
    log(f"[GSE151179] all tumor specimens with RAI label: {int(keep_mask.sum())}; "
        f"primary-tumor only: {int(primary_mask.sum())}")
    gsms = meta["gsms"]
    samples_keep = [g for i, g in enumerate(gsms) if keep_mask[i] and g in expr_df.columns]
    y = np.array([y_rai[i] for i, g in enumerate(gsms) if keep_mask[i] and g in expr_df.columns])
    # Build probe -> symbol
    panel_lookup_path = GEO_DIR / "platforms" / "GPL23159_panel_lookup.json"
    if panel_lookup_path.exists():
        log(f"[GSE151179] using cached panel lookup: {panel_lookup_path}")
        panel_lookup = json.loads(panel_lookup_path.read_text())
    else:
        gpl_table = GEO_DIR / "platforms" / "GPL23159.txt"
        gpl_table.parent.mkdir(parents=True, exist_ok=True)
        if not gpl_table.exists() or gpl_table.stat().st_size < 1000:
            url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&form=text&view=data&acc=GPL23159"
            log("[GSE151179] downloading GPL23159 platform table (~110 MB)")
            r = cached_download(url, gpl_table)
            if r is None:
                return None
        panel_lookup = parse_gpl23159_for_panel(gene_order, gpl_table)
        panel_lookup_path.write_text(json.dumps(panel_lookup, indent=2))
        log(f"[GSE151179] cached panel lookup -> {panel_lookup_path}")
    # Index normalisation: matrix probe IDs may include the integer ID for control
    expr_df.index = expr_df.index.astype(str)
    matched: List[str] = []
    out_rows: Dict[str, pd.Series] = {}
    for g in gene_order:
        probes = panel_lookup.get(g, [])
        avail = [p for p in probes if p in expr_df.index]
        if avail:
            out_rows[g] = expr_df.loc[avail].mean(axis=0)
            matched.append(g)
    log(f"[GSE151179] genes matched via probes: {len(matched)}/{len(gene_order)} ({matched})")
    if len(matched) < 6:
        log("[GSE151179] insufficient gene match, skipping")
        return None
    gx = pd.DataFrame(out_rows).T
    samples = [s for s in samples_keep if s in gx.columns]
    if len(np.unique(y)) < 2 or len(samples) < 6:
        log("[GSE151179] insufficient labelled samples after probe match")
        return None
    y = np.array([y[i] for i, g in enumerate(samples_keep) if g in samples])
    X = gx.loc[matched, samples].T
    Xfull = pd.DataFrame(0.0, index=X.index, columns=gene_order)
    for g in matched:
        Xfull[g] = X[g]
    proba = model.predict_proba(scaler.transform(Xfull.values))[:, 1]
    auc_all, lo_all, hi_all = bootstrap_auc(y, proba, n=1000, seed=42)
    log(f"[GSE151179] all-tumor-specimens AUC={auc_all:.3f} [{lo_all:.3f}, {hi_all:.3f}] "
        f"(n={len(y)}, RAI-refractory vs avid)")

    # Also compute primary-tumor-only AUC (cleaner headline; LNM exclude)
    primary_idx = []
    primary_y = []
    primary_p = []
    for i, gsm in enumerate(samples):
        # find original index in metadata
        try:
            j = gsms.index(gsm)
        except ValueError:
            continue
        if primary_mask[j]:
            primary_idx.append(i)
            primary_y.append(y[i])
            primary_p.append(proba[i])
    primary_y_arr = np.array(primary_y)
    primary_p_arr = np.array(primary_p)
    if len(np.unique(primary_y_arr)) >= 2 and len(primary_y_arr) >= 6:
        auc_p, lo_p, hi_p = bootstrap_auc(primary_y_arr, primary_p_arr, n=1000, seed=42)
        log(f"[GSE151179] PRIMARY-tumor-only AUC={auc_p:.3f} [{lo_p:.3f}, {hi_p:.3f}] "
            f"(n={len(primary_y_arr)}, refractory={int(primary_y_arr.sum())}, avid={int((1-primary_y_arr).sum())})")
    else:
        auc_p = lo_p = hi_p = float("nan")

    # Pick the cleaner / better-powered contrast as the headline:
    # use all-tumor-specimens (more power, n=22) but record both
    auc, lo, hi = auc_all, lo_all, hi_all

    # Save the patient-level table for the headline RAI clinical validation
    rai_df = pd.DataFrame({
        "sample_id": samples,
        "rai_label": ["refractory" if v == 1 else "avid" for v in y],
        "y": y,
        "p_DM1_score": proba,
        "n_genes_matched": int(len(matched)),
        "is_primary_tumor": [bool(primary_mask[gsms.index(s)]) if s in gsms else False
                             for s in samples],
    })
    rai_df.to_csv(OUT_DIR / "U1B_RAI_clinical_validation.tsv", sep="\t", index=False)
    log(f"[GSE151179] wrote U1B_RAI_clinical_validation.tsv (n={len(rai_df)})")
    return {
        "cohort": "GSE151179",
        "n": int(len(y)),
        "n_genes_matched": int(len(matched)),
        "ground_truth_type": "RAI_refractory_vs_avid",
        "n_pos": int(y.sum()),
        "n_neg": int((1 - y).sum()),
        "auc": auc,
        "ci_lo": lo,
        "ci_hi": hi,
        "auc_primary_tumor_only": float(auc_p),
        "ci_lo_primary": float(lo_p),
        "ci_hi_primary": float(hi_p),
        "n_primary": int(len(primary_y_arr)),
        "n_pos_primary": int(primary_y_arr.sum()) if len(primary_y_arr) else 0,
        "n_neg_primary": int((1 - primary_y_arr).sum()) if len(primary_y_arr) else 0,
        "note": "★ TRUE RAI clinical ground truth (papillary thyroid carcinoma); "
                "headline AUC uses all tumor specimens (primary + LNM); "
                "auc_primary_tumor_only is the conservative subset",
    }


# -----------------------------------------------------------------------------
# Meta-analysis (random effects on Fisher-Z transformed AUCs)
# -----------------------------------------------------------------------------
def meta_analysis(rows: List[Dict]) -> Dict:
    """Random-effects pooled AUC with DerSimonian-Laird I^2.

    AUC variance approximated via Hanley-McNeil; transform AUC to logit-space
    for pooling.
    """
    valid = [r for r in rows if r.get("auc") is not None and not pd.isna(r["auc"])]
    if len(valid) < 2:
        return {"error": "fewer than 2 valid cohorts", "k": len(valid)}
    # Hanley-McNeil variance of AUC
    def hm_var(auc, n_pos, n_neg):
        Q1 = auc / (2 - auc)
        Q2 = 2 * auc * auc / (1 + auc)
        if n_pos == 0 or n_neg == 0:
            return np.nan
        return (auc * (1 - auc) + (n_pos - 1) * (Q1 - auc * auc)
                + (n_neg - 1) * (Q2 - auc * auc)) / (n_pos * n_neg)

    aucs = np.array([r["auc"] for r in valid])
    vars_ = np.array([hm_var(r["auc"], r["n_pos"], r["n_neg"]) for r in valid])
    # Drop any with nan variance
    ok = ~np.isnan(vars_) & (vars_ > 0)
    aucs = aucs[ok]; vars_ = vars_[ok]
    valid_used = [v for v, k in zip(valid, ok) if k]
    if len(aucs) < 2:
        return {"error": "fewer than 2 cohorts after variance filtering"}
    # Fixed-effect weighted mean (for tau^2 estimation)
    w = 1 / vars_
    mu_fe = float(np.sum(w * aucs) / np.sum(w))
    Q = float(np.sum(w * (aucs - mu_fe) ** 2))
    df = len(aucs) - 1
    C = float(np.sum(w) - np.sum(w * w) / np.sum(w))
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    # Random-effects weights
    w_re = 1 / (vars_ + tau2)
    mu_re = float(np.sum(w_re * aucs) / np.sum(w_re))
    se_re = float(np.sqrt(1 / np.sum(w_re)))
    ci_lo = mu_re - 1.96 * se_re
    ci_hi = mu_re + 1.96 * se_re
    return {
        "k_cohorts": int(len(aucs)),
        "cohorts_used": [r["cohort"] for r in valid_used],
        "auc_pooled_random_effects": mu_re,
        "se_pooled": se_re,
        "ci_lo": float(ci_lo),
        "ci_hi": float(ci_hi),
        "Q": Q,
        "df": df,
        "tau2": tau2,
        "I2_percent": float(I2),
        "auc_pooled_fixed_effect": mu_fe,
    }


# -----------------------------------------------------------------------------
# Forest plot
# -----------------------------------------------------------------------------
def make_forest_plot(rows: List[Dict], pooled: Dict, out_png: Path):
    fig, ax = plt.subplots(figsize=(11, 6.5))
    rows_plot = [r for r in rows if r.get("auc") is not None and not pd.isna(r["auc"])]
    rows_plot = sorted(rows_plot, key=lambda r: r["cohort"] != "TCGA-THCA")  # TCGA first
    n = len(rows_plot)
    y_pos = np.arange(n + 2)[::-1]  # +1 spacer +1 pooled
    palette = {"TCGA-THCA": "#1f77b4", "GSE76039": "#ff7f0e", "GSE213647": "#2ca02c",
               "GSE27155": "#d62728", "GSE126698": "#9467bd", "GSE65144": "#8c564b",
               "GSE60542": "#e377c2", "GSE151179": "#000000"}
    for i, r in enumerate(rows_plot):
        c = palette.get(r["cohort"], "#7f7f7f")
        x = r["auc"]; lo = r["ci_lo"]; hi = r["ci_hi"]
        ax.errorbar([x], [y_pos[i]], xerr=[[x - lo], [hi - x]], fmt="o",
                    color=c, capsize=4, lw=1.5, markersize=10,
                    markeredgecolor="black", markerfacecolor=c)
        label = f"{r['cohort']} (n={r['n']}, {r['n_genes_matched']}/8 g)"
        if r["cohort"] == "GSE151179":
            label = f"★ {label}"
        ax.text(0.39, y_pos[i],
                f"{label}", ha="left", va="center", fontsize=9.5,
                fontweight="bold" if r["cohort"] == "GSE151179" else "normal")
        ax.text(1.02, y_pos[i],
                f"{x:.3f} [{lo:.3f}-{hi:.3f}]", ha="left", va="center", fontsize=9)
        # contrast tag
        ax.text(1.21, y_pos[i],
                r["ground_truth_type"].replace("_", " "),
                ha="left", va="center", fontsize=8, style="italic", color="dimgray")
    # Pooled diamond
    if "auc_pooled_random_effects" in pooled:
        m = pooled["auc_pooled_random_effects"]
        lo = pooled["ci_lo"]; hi = pooled["ci_hi"]
        yp = y_pos[-1]
        ax.fill([lo, m, hi, m, lo],
                [yp, yp + 0.25, yp, yp - 0.25, yp],
                color="#444", alpha=0.85)
        ax.text(0.39, yp,
                f"Pooled (RE): k={pooled['k_cohorts']}, I²={pooled['I2_percent']:.1f}%",
                ha="left", va="center", fontsize=10, fontweight="bold")
        ax.text(1.02, yp, f"{m:.3f} [{lo:.3f}-{hi:.3f}]",
                ha="left", va="center", fontsize=10, fontweight="bold")
    ax.axvline(0.5, ls="--", color="gray", lw=1)
    ax.set_xlim(0.39, 1.45)
    ax.set_xlabel("AUC (direction-corrected, 95% CI)")
    ax.set_yticks([])
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xticklabels(["0.5", "0.6", "0.7", "0.8", "0.9", "1.0"])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    title = ("U1B: 7-cohort meta-analysis of 8-gene de-differentiation panel\n"
             "(TCGA-THCA + 6 GEO; ★ = true RAI clinical ground truth)")
    ax.set_title(title, fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_png, dpi=160, bbox_inches="tight")
    fig.savefig(str(out_png).replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    log(f"  wrote {out_png}")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def ensure_downloads():
    log("=== ensuring downloads ===")
    targets = {
        "GSE65144_series_matrix.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65144/matrix/GSE65144_series_matrix.txt.gz",
        "GSE60542_series_matrix.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60542/matrix/GSE60542_series_matrix.txt.gz",
        "GSE151179_series_matrix.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151179/matrix/GSE151179_series_matrix.txt.gz",
        "GSE151180_series_matrix.txt.gz":
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151180/matrix/GSE151180_series_matrix.txt.gz",
    }
    for fname, url in targets.items():
        cached_download(url, GEO_DIR / fname)
    # GSE192683 = SKIP (non-thyroid Piaractus mesopotamicus)
    log("GSE192683 skipped: non-thyroid SuperSeries (Piaractus mesopotamicus muscle cells)")


def main() -> int:
    log("=== v17 ULTIMATE U1B: 7-cohort meta-analysis of 8-gene panel ===")
    ensure_downloads()

    # Train TCGA on full GENE_8 first; we'll re-train on cohort-specific intersections
    tcga = load_tcga_expr()
    tcga_have = [g for g in GENE_8 if g in tcga.index]
    log(f"GENE_8 in TCGA: {len(tcga_have)}/8 ({tcga_have})")
    scaler, model, train_auc = train_tcga_8gene(tcga_have)

    rows: List[Dict] = []

    # In-sample TCGA
    r = run_tcga(scaler, model, tcga_have)
    if r:
        rows.append(r)

    # Already-local cohorts
    for fn in (run_gse76039, run_gse27155, run_gse213647, run_gse126698):
        try:
            r = fn(scaler, model, tcga_have)
        except Exception as e:
            log(f"[{fn.__name__}] EXCEPTION: {e}")
            r = None
        if r:
            rows.append(r)

    # New cohorts
    for fn in (run_gse65144, run_gse60542, run_gse151179):
        try:
            r = fn(scaler, model, tcga_have)
        except Exception as e:
            log(f"[{fn.__name__}] EXCEPTION: {e}")
            r = None
        if r:
            rows.append(r)

    # Save per-cohort table
    df = pd.DataFrame(rows)
    main_cols = ["cohort", "n", "n_genes_matched", "ground_truth_type",
                 "n_pos", "n_neg", "auc", "ci_lo", "ci_hi", "note"]
    extra = [c for c in df.columns if c not in main_cols]
    df = df[main_cols + extra]
    for c in ("auc", "ci_lo", "ci_hi"):
        df[c] = df[c].round(4)
    out_tsv = OUT_DIR / "U1B_per_cohort_real_aucs.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    log(f"wrote {out_tsv}")
    log("Per-cohort AUCs:")
    print(df.to_string(index=False))

    # Meta-analysis
    pooled = meta_analysis(rows)
    log(f"Meta-analysis: {json.dumps(pooled, indent=2, default=str)}")
    pooled_df = pd.DataFrame([pooled])
    pooled_df.to_csv(OUT_DIR / "U1B_meta_analysis_real.tsv", sep="\t", index=False)

    # Identify GSE151179
    rai = next((r for r in rows if r["cohort"] == "GSE151179"), None)
    rai_status = "fail"
    if rai is not None and not pd.isna(rai["auc"]):
        rai_status = "success"
        log(f"★ RAI clinical validation succeeded: AUC={rai['auc']:.3f} "
            f"[{rai['ci_lo']:.3f},{rai['ci_hi']:.3f}] (n={rai['n']})")

    # Forest plot
    make_forest_plot(rows, pooled, FIG_DIR / "U1B_forest_plot_7cohort.png")

    # Summary JSON
    summary = {
        "sprint": "v17_ULTIMATE_U1B",
        "goal": "7-cohort meta-analysis of 8-gene de-differentiation panel",
        "gene_8_panel": GENE_8,
        "tcga_train_in_sample_auc": train_auc,
        "k_cohorts_attempted": len(rows),
        "cohorts": rows,
        "meta_analysis": pooled,
        "rai_clinical_validation": {
            "status": rai_status,
            "result": rai,
            "ground_truth_source": "GSE151179 patient_rai_responce field",
        },
        "skipped_cohorts": {
            "GSE192683": "non-thyroid SuperSeries (fish muscle cells)",
            "GSE151180": "miRNA sub-series (gene-expression sibling GSE151179 used instead)",
        },
        "outputs": {
            "per_cohort_tsv": str(OUT_DIR / "U1B_per_cohort_real_aucs.tsv"),
            "meta_tsv": str(OUT_DIR / "U1B_meta_analysis_real.tsv"),
            "rai_tsv": str(OUT_DIR / "U1B_RAI_clinical_validation.tsv"),
            "forest_png": str(FIG_DIR / "U1B_forest_plot_7cohort.png"),
        },
    }
    jdump(summary, OUT_DIR / "U1B_summary.json")
    log("=== U1B done ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
