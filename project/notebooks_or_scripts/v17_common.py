#!/usr/bin/env python3
from __future__ import annotations

import io
import json
import math
import tarfile
import time
import gzip
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
META = ROOT / "metadata"
RES = ROOT / "results" / "v17"
TAB = RES / "tables"
FIG = RES / "figs"
CACHE = RES / "cache"
RPT = ROOT / "reports" / "v17"
VENV_PY = ROOT / ".venv" / "bin" / "python"
LOG = RES / "orchestrator.log"

for d in (TAB, FIG, CACHE, RPT):
    d.mkdir(parents=True, exist_ok=True)

EXPR_PATHS = {
    "TCGA-THCA": Path("/data/thca/data_processed/bulk_rnaseq_v3/TCGA-THCA_v3_log2.tsv"),
    "GSE27155": Path("/data/thca/data_processed/microarray_v3/GSE27155_v3_log2.tsv"),
    "GSE126698": Path("/data/thca/data_processed/bulk_rnaseq_v3/GSE126698_v3_log2.tsv"),
    "GSE76039": Path("/data/thca/data_processed/microarray_v3/GSE76039_v3_log2.tsv"),
    "GSE213647": Path("/data/thca/data_processed/bulk_rnaseq_v3/GSE213647_v3_log2.tsv"),
}


def log_line(path: Path, msg: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{stamp}] {msg}"
    print(line, flush=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def clean_gene_list(path: Path) -> list[str]:
    return [x.strip() for x in path.read_text().splitlines() if x.strip() and not x.startswith("#")]


TDS16 = clean_gene_list(META / "tds16_genes.txt")
BRS71 = clean_gene_list(META / "brs71_genes.txt")
TIERA67 = clean_gene_list(META / "v3_tierA67_clean_genes.txt")
RAI_GENES = ["SLC5A5", "TPO", "TSHR", "TG", "DIO1", "DIO2", "PAX8", "NKX2-1"]
EXTRA_DRIVER_GENES = ["DICER1", "EIF1AX", "PPM1D", "EZH1", "EZH2", "ATM", "CHEK2", "TP53", "PIK3CA", "AKT1"]
TCGA_MUT_MANIFEST = Path("/data/thca/data_raw/gdc/TCGA-THCA/tcga_thca_mutation_manifest.tsv")
TCGA_MUT_DIR = Path("/data/thca/data_raw/gdc/TCGA-THCA/mutation")
TCGA_CLINICAL_EXT = ROOT / "results" / "tables" / "tcga_thca_clinical_extended.tsv"


def load_sample_master() -> pd.DataFrame:
    df = pd.read_csv(META / "sample_master_v3.tsv", sep="\t")
    return df


def normalize_barcode(x: str) -> str:
    x = str(x)
    if x.startswith("TCGA-"):
        return x[:12]
    return x


def read_expr(dataset: str, genes: Iterable[str] | None = None) -> pd.DataFrame:
    path = EXPR_PATHS[dataset]
    if not path.exists():
        raise FileNotFoundError(f"missing expression matrix: {path}")
    df = pd.read_csv(path, sep="\t")
    gene_col = df.columns[0]
    df = df.rename(columns={gene_col: "gene_symbol"})
    df["gene_symbol"] = df["gene_symbol"].astype(str)
    df = df.drop_duplicates("gene_symbol").set_index("gene_symbol")
    if genes is not None:
        keep = [g for g in genes if g in df.index]
        df = df.loc[keep]
    return df.T


def align_dataset_expr_meta(dataset: str, genes: Iterable[str] | None = None, tumor_only: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = read_expr(dataset, genes=genes)
    meta = load_sample_master()
    meta = meta[meta["dataset"] == dataset].copy()
    if tumor_only:
        meta = meta[meta["normal_vs_tumor"] != "normal"].copy()
    common = [s for s in meta["sample_id"] if s in expr.index]
    meta = meta.set_index("sample_id").loc[common].reset_index()
    expr = expr.loc[common].copy()
    return expr, meta


def signature_score(expr: pd.DataFrame, genes: Iterable[str]) -> pd.Series:
    keep = [g for g in genes if g in expr.columns]
    if not keep:
        return pd.Series(np.nan, index=expr.index)
    return expr[keep].mean(axis=1)


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    pvals = np.asarray(pvals, dtype=float)
    out = np.full_like(pvals, np.nan)
    mask = np.isfinite(pvals)
    if not mask.any():
        return out
    pv = pvals[mask]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adj = ranked * n / (np.arange(n) + 1.0)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    tmp = np.empty_like(adj)
    tmp[order] = adj
    out[mask] = tmp
    return out


def fetch_bytes(url: str, log_path: Path, retries: int = 3, timeout: int = 120) -> bytes:
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            log_line(log_path, f"fetch {url} attempt={attempt}")
            r = requests.get(url, timeout=timeout)
            r.raise_for_status()
            return r.content
        except Exception as e:
            last_err = e
            time.sleep(2 * attempt)
    raise RuntimeError(f"failed fetch after {retries} attempts: {url} ({last_err})")


def fetch_study_tar(log_path: Path) -> Path:
    out = CACHE / "thca_tcga_pan_can_atlas_2018.tar.gz"
    if not out.exists():
        urls = [
            "https://cbioportal-datahub.s3.amazonaws.com/thca_tcga_pan_can_atlas_2018.tar.gz",
            "https://github.com/cBioPortal/datahub-public/raw/master/public/thca_tcga_pan_can_atlas_2018.tar.gz",
        ]
        last_err = None
        for url in urls:
            try:
                data = fetch_bytes(url, log_path)
                out.write_bytes(data)
                break
            except Exception as e:
                last_err = e
        if not out.exists():
            raise RuntimeError(f"failed all study tar fetch routes ({last_err})")
    return out


def extract_from_tar(member_suffix: str, log_path: Path) -> pd.DataFrame:
    try:
        tar_path = fetch_study_tar(log_path)
        with tarfile.open(tar_path, "r:gz") as tf:
            names = [m.name for m in tf.getmembers() if m.isfile()]
            target = next((n for n in names if n.endswith(member_suffix)), None)
            if target is None:
                raise FileNotFoundError(f"{member_suffix} not found in {tar_path.name}")
            buf = tf.extractfile(target).read()
        return pd.read_csv(io.BytesIO(buf), sep="\t", comment="#", low_memory=False)
    except Exception as e:
        log_line(log_path, f"tar fallback for {member_suffix}: {e}")
        raw_url = f"https://github.com/cBioPortal/datahub-public/raw/master/public/thca_tcga_pan_can_atlas_2018/{member_suffix}"
        buf = fetch_bytes(raw_url, log_path, timeout=180)
        return pd.read_csv(io.BytesIO(buf), sep="\t", comment="#", low_memory=False)


def km_table(time_vals: np.ndarray, event_vals: np.ndarray) -> pd.DataFrame:
    order = np.argsort(time_vals)
    t = np.asarray(time_vals, dtype=float)[order]
    e = np.asarray(event_vals, dtype=int)[order]
    uniq = np.unique(t[e == 1])
    n_at_risk = len(t)
    surv = 1.0
    rows = [{"time": 0.0, "survival": 1.0}]
    for ut in uniq:
        d = int(((t == ut) & (e == 1)).sum())
        c = int(((t == ut) & (e == 0)).sum())
        if n_at_risk > 0:
            surv *= (1.0 - d / n_at_risk)
        rows.append({"time": float(ut), "survival": float(surv)})
        n_at_risk -= d + c
    return pd.DataFrame(rows)


def logrank_pair(time_a: np.ndarray, event_a: np.ndarray, time_b: np.ndarray, event_b: np.ndarray) -> float:
    times = np.unique(np.concatenate([time_a[event_a == 1], time_b[event_b == 1]]))
    O1 = E1 = V1 = 0.0
    for t in times:
        n1 = float((time_a >= t).sum())
        n2 = float((time_b >= t).sum())
        d1 = float(((time_a == t) & (event_a == 1)).sum())
        d2 = float(((time_b == t) & (event_b == 1)).sum())
        n = n1 + n2
        d = d1 + d2
        if n <= 1 or d == 0:
            continue
        e1 = d * (n1 / n)
        v1 = (n1 * n2 * d * (n - d)) / (n * n * (n - 1))
        O1 += d1
        E1 += e1
        V1 += v1
    if V1 <= 0:
        return np.nan
    z2 = (O1 - E1) ** 2 / V1
    return float(stats.chi2.sf(z2, 1))


def compact_df(df: pd.DataFrame, n: int = 10) -> str:
    if df.empty:
        return "(empty)"
    return df.head(n).to_string(index=False)


def load_tcga_clinical_extended() -> pd.DataFrame:
    df = pd.read_csv(TCGA_CLINICAL_EXT, sep="\t")
    df["tcga12"] = df["sample_id"].astype(str).map(normalize_barcode)
    return df


def load_local_tcga_mutation_summary(log_path: Path) -> pd.DataFrame:
    cache_path = CACHE / "tcga_thca_local_mutation_summary.tsv"
    if cache_path.exists():
        return pd.read_csv(cache_path, sep="\t")

    manifest = pd.read_csv(TCGA_MUT_MANIFEST, sep="\t")
    rows = []
    use_genes = set(["BRAF", "HRAS", "KRAS", "NRAS", "TERT", "RET", "NTRK1", "NTRK3", "ALK", "PAX8", "PPARG", "THADA"] + EXTRA_DRIVER_GENES)
    tert_positions = {1295228, 1295250, 1295113, 1295135}
    ras_hotspots = {"Q61", "G12", "G13"}
    for i, rec in manifest.iterrows():
        sample_id = str(rec["sample_submitter_id"])
        path = TCGA_MUT_DIR / f"{rec['file_id']}.maf.gz"
        if not path.exists():
            continue
        try:
            maf = pd.read_csv(
                path,
                sep="\t",
                comment="#",
                compression="gzip",
                low_memory=False,
                usecols=["Hugo_Symbol", "Variant_Classification", "Start_Position", "Protein_Change", "HGVSp_Short", "Tumor_Sample_Barcode"],
            )
        except ValueError:
            maf = pd.read_csv(
                path,
                sep="\t",
                comment="#",
                compression="gzip",
                low_memory=False,
            )
        maf.columns = [c.strip() for c in maf.columns]
        if "Hugo_Symbol" not in maf.columns:
            continue
        keep = maf["Hugo_Symbol"].astype(str).isin(use_genes)
        maf = maf.loc[keep].copy()
        if maf.empty:
            rows.append(
                {
                    "sample_id": sample_id,
                    "tcga12": normalize_barcode(sample_id),
                    "mutation_genes": "",
                    "braf_v600e": 0,
                    "ras_hotspot": 0,
                    "tert_promoter": 0,
                }
            )
            continue
        hgvs_col = "HGVSp_Short" if "HGVSp_Short" in maf.columns else None
        prot_col = "Protein_Change" if "Protein_Change" in maf.columns else None
        prot = maf[hgvs_col] if hgvs_col else maf[prot_col] if prot_col else pd.Series("", index=maf.index)
        prot = prot.astype(str)
        braf_v600e = ((maf["Hugo_Symbol"].astype(str) == "BRAF") & prot.str.contains("V600E", na=False)).any()
        ras_hotspot = (
            maf["Hugo_Symbol"].astype(str).isin(["HRAS", "KRAS", "NRAS"])
            & prot.str.contains("|".join(ras_hotspots), na=False)
        ).any()
        vc = maf["Variant_Classification"].astype(str) if "Variant_Classification" in maf.columns else pd.Series("", index=maf.index)
        pos = pd.to_numeric(maf["Start_Position"], errors="coerce") if "Start_Position" in maf.columns else pd.Series(np.nan, index=maf.index)
        tert_promoter = (
            (maf["Hugo_Symbol"].astype(str) == "TERT")
            & (vc.str.contains("5'Flank|5UTR|UTR|RNA", case=False, na=False) | pos.isin(tert_positions))
        ).any()
        genes = sorted(set(maf["Hugo_Symbol"].astype(str)))
        rows.append(
            {
                "sample_id": sample_id,
                "tcga12": normalize_barcode(sample_id),
                "mutation_genes": ";".join(genes),
                "braf_v600e": int(bool(braf_v600e)),
                "ras_hotspot": int(bool(ras_hotspot)),
                "tert_promoter": int(bool(tert_promoter)),
            }
        )
        if (i + 1) % 100 == 0:
            log_line(log_path, f"local mutation summary progress {i + 1}/{len(manifest)}")
    out = pd.DataFrame(rows).drop_duplicates("sample_id")
    out.to_csv(cache_path, sep="\t", index=False)
    return out
