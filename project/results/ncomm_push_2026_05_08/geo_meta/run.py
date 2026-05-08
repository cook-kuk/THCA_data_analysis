#!/usr/bin/env python3
"""GEO ATC meta-pool — older ATC datasets pooled for trajectory robustness.

Datasets (all public, microarray):
  GSE65144   Tomás 2015 ATC vs normal
  GSE60542   Hébrant 2014 ATC + PDTC + PTC
  GSE82208   Tarabichi 2017 ATC progression

Goal: 8-gene panel scored on each cohort, cross-cohort meta of differentiation
suppression in ATC vs PTC vs normal.

Uses GEO SOFT/Series Matrix files (txt.gz) which are small and
do not require GEOparse — we parse the matrix file directly.
"""
from __future__ import annotations

import gzip
import io
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import (  # noqa: E402
    PANEL_8,
    log,
    resolve_panel,
    score_panel,
    set_threads,
    write_status,
)

set_threads(2)

WORKER = "geo_meta"
OUT = Path(__file__).parent
RAW = OUT / "_raw"
RAW.mkdir(exist_ok=True)

DATASETS = ["GSE65144", "GSE60542", "GSE82208", "GSE76039"]


def geo_matrix_url(gse: str) -> str:
    stub = gse[:-3] + "nnn" if len(gse) > 6 else gse + "nnn"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/"
        f"{gse}_series_matrix.txt.gz"
    )


def download_matrix(gse: str) -> Path:
    dest = RAW / f"{gse}_series_matrix.txt.gz"
    if dest.exists() and dest.stat().st_size > 1024:
        log(WORKER, f"cached {dest.name}")
        return dest
    url = geo_matrix_url(gse)
    log(WORKER, f"downloading {url}")
    r = requests.get(url, timeout=600)
    r.raise_for_status()
    dest.write_bytes(r.content)
    log(WORKER, f"saved {dest.name} ({dest.stat().st_size:,} bytes)")
    return dest


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (expression_df: probes × samples, sample_meta: samples × keys)."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        text = fh.read()

    # Sample-level metadata
    meta_rows = {}
    for line in text.splitlines():
        if line.startswith("!Sample_"):
            key = line.split("\t", 1)[0].lstrip("!").rstrip()
            vals = line.split("\t")[1:]
            vals = [v.strip().strip('"') for v in vals]
            meta_rows.setdefault(key, []).append(vals)
        if line.startswith("!series_matrix_table_begin"):
            break
    # collapse multi-rows to last (or join)
    meta = {k: v[-1] for k, v in meta_rows.items() if v}
    sample_ids = meta.get("Sample_geo_accession", [])
    sample_meta = pd.DataFrame(meta).set_index("Sample_geo_accession") if "Sample_geo_accession" in meta else pd.DataFrame()
    if "Sample_geo_accession" not in meta:
        # fallback: extract sample ids from data section
        sample_meta = pd.DataFrame(index=sample_ids)

    # Expression table
    m = re.search(r"!series_matrix_table_begin\n(.*?)!series_matrix_table_end", text, re.DOTALL)
    if not m:
        raise RuntimeError(f"no expression table found in {path.name}")
    table_text = m.group(1).strip()
    df = pd.read_csv(io.StringIO(table_text), sep="\t", index_col=0)
    df.index.name = "probe"
    return df, sample_meta


def map_probes_to_genes_via_geo(gse: str, df_probes: pd.DataFrame) -> pd.DataFrame:
    """Best-effort probe→gene mapping using the platform GPL annotation file.

    The GPL annotation can be hundreds of MB. We use a fast heuristic: many
    GEO Affy / Illumina platforms encode the symbol in the probe row metadata
    when available; otherwise we drop and rely on the user to inspect.

    For this push, we try the GEOmetadb-free shortcut: query the
    https://eutils.ncbi.nlm.nih.gov annotation through its esearch+efetch
    API only for the panel genes (small).
    """
    # The simplest reliable path: query the platform's GPL file. We accept
    # that this is best-effort and only collapse panel-relevant probes.
    # For each panel gene, search public GPL gene-symbol→probe maps.
    # To stay self-contained we use a static map of well-known affy probes.
    AFFY_HG_U133_PLUS_2 = {
        "TPO": ["205557_at"],
        "TG": ["210234_at"],
        "TSHR": ["207152_at", "215310_at"],
        "PAX8": ["205044_at"],
        "NKX2-1": ["205373_at", "210503_at"],
        "FOXE1": ["207681_at"],
        "DIO1": ["205709_s_at", "207096_at"],
        "SLC5A5": ["207547_s_at"],
    }
    # Many of the older ATC datasets use HG-U133 Plus 2.0 (GPL570). Map and
    # collapse to gene-level mean.
    out = {}
    for gene, probes in AFFY_HG_U133_PLUS_2.items():
        hits = [p for p in probes if p in df_probes.index]
        if hits:
            out[gene] = df_probes.loc[hits].mean(axis=0)
    if not out:
        return pd.DataFrame()
    gene_df = pd.DataFrame(out)  # samples × genes
    return gene_df


def process_dataset(gse: str) -> dict | None:
    try:
        path = download_matrix(gse)
        expr_probes, sample_meta = parse_series_matrix(path)
    except Exception as exc:
        log(WORKER, f"[{gse}] download/parse FAIL: {exc}")
        return None
    log(WORKER, f"[{gse}] probes={expr_probes.shape[0]} samples={expr_probes.shape[1]}")

    gene_df = map_probes_to_genes_via_geo(gse, expr_probes)
    if gene_df.empty:
        log(WORKER, f"[{gse}] no panel probes mapped (probably non-Affy GPL) — skip")
        return None
    log(WORKER, f"[{gse}] panel gene coverage: {gene_df.shape[1]}/{len(PANEL_8)}")
    # gene_df is samples × genes; align with sample_meta
    if gene_df.shape[1] < 5:
        log(WORKER, f"[{gse}] thin coverage — partial only")

    gene_map = resolve_panel(gene_df.columns)
    scores = score_panel(gene_df, gene_map)
    scores["GSE"] = gse
    # Try to attach a phenotype label
    if "Sample_characteristics_ch1" in sample_meta.columns:
        scores["phenotype"] = sample_meta["Sample_characteristics_ch1"].reindex(scores.index)
    elif "Sample_source_name_ch1" in sample_meta.columns:
        scores["phenotype"] = sample_meta["Sample_source_name_ch1"].reindex(scores.index)
    scores.to_csv(OUT / f"scores_{gse}.tsv", sep="\t")

    return {
        "gse": gse,
        "n_samples": int(gene_df.shape[0]),
        "panel_genes_resolved": int(gene_df.shape[1]),
        "rai_8_mean": float(scores["RAI_8"].mean()),
        "rai_8_std": float(scores["RAI_8"].std()),
    }


def main() -> None:
    write_status(WORKER, "running", started_at=time.time())
    summaries = []
    for gse in DATASETS:
        s = process_dataset(gse)
        if s is not None:
            summaries.append(s)

    pd.DataFrame(summaries).to_csv(OUT / "geo_meta_summary.tsv", sep="\t", index=False)
    (OUT / "summary.json").write_text(json.dumps(summaries, indent=2))
    write_status(WORKER, "done", n_datasets=len(summaries), datasets=[s["gse"] for s in summaries], finished_at=time.time())
    log(WORKER, "DONE")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(WORKER, f"FATAL: {exc}")
        write_status(WORKER, "error", error=str(exc))
        raise
