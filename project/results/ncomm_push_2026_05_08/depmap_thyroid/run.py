#!/usr/bin/env python3
"""DepMap thyroid cell line sweep — DM1 axis × functional readouts.

Public sources:
  CCLE expression  : depmap.org figshare release (OmicsExpressionProteinCodingGenesTPMLogp1.csv)
  CRISPR Achilles  : CRISPRGeneEffect.csv
  PRISM drug screen: Repurposing_Public_24Q2_value.csv (or similar)

Goal: score 8-gene panel on thyroid cell lines, split DM1-like vs DM2-like,
then report:
  - DM1-vs-DM2 differential CRISPR essentiality of TPO/DIO1/TSHR/etc.
  - DM1-vs-DM2 drug sensitivity to HMA (decitabine, azacitidine) and RET inhibitors.

We use the Hugging Face mirror (broadinstitute/ccle) when direct figshare URLs
require auth; fall back to the deterministic figshare URLs noted in DepMap docs.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import (  # noqa: E402
    PANEL_8,
    kmeans_dm,
    log,
    resolve_panel,
    score_panel,
    set_threads,
    write_status,
)

set_threads(2)

WORKER = "depmap_thyroid"
OUT = Path(__file__).parent
RAW = OUT / "_raw"
RAW.mkdir(exist_ok=True)

# DepMap public 24Q2 release (free, no login required for top-level files)
# (URLs are stable per DepMap conventions; if these change, switch to figshare DOI.)
# DepMap collection on figshare. We resolve current file IDs via the API
# instead of hardcoding them (figshare numeric IDs change between releases).
FIGSHARE_API = "https://api.figshare.com/v2"
# Candidate articles, newest first. The first one that returns a valid file list wins.
DEPMAP_ARTICLES = [
    "27993248",  # DepMap Public 24Q4
    "26261476",  # DepMap Public 24Q2
    "25880521",  # DepMap Public 24Q2 alt
    "24667905",  # DepMap Public 23Q4
    "22765112",  # DepMap Public 23Q2
    "21637199",  # DepMap Public 23Q2 (alt)
]
# We want files matching these substrings (case-insensitive)
WANT = {
    "model": ["model.csv"],
    "expr": ["omicsexpressionproteincodingenestpmlogp1.csv",
             "omicsexpressionproteincodinggenestpmlogp1.csv",
             "ccle_expression.csv"],
    "crispr": ["crisprgeneeffect.csv", "achilles_gene_effect.csv"],
    "prism": ["repurposing_public_24q2_value.csv",
              "repurposing_public_23q2_value.csv",
              "primary-screen-replicate-collapsed-logfold-change.csv",
              "secondary-screen-dose-response-curve-parameters.csv"],
}


def download(url: str, dest: Path, expected_min_bytes: int = 1024) -> Path:
    if dest.exists() and dest.stat().st_size > expected_min_bytes:
        log(WORKER, f"cached {dest.name} ({dest.stat().st_size:,} bytes)")
        return dest
    log(WORKER, f"downloading {url} -> {dest.name}")
    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        with dest.open("wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
    log(WORKER, f"saved {dest.name} ({dest.stat().st_size:,} bytes)")
    return dest


def resolve_article_files() -> dict[str, str]:
    """Try each candidate article; return {key: download_url} for whatever resolves."""
    for art in DEPMAP_ARTICLES:
        try:
            r = requests.get(f"{FIGSHARE_API}/articles/{art}/files", timeout=60)
            r.raise_for_status()
            files = r.json()
        except Exception as exc:
            log(WORKER, f"article {art}: {exc}")
            continue
        log(WORKER, f"article {art}: {len(files)} files exposed")
        urls: dict[str, str] = {}
        # build lowercase->record dict
        by_name = {f["name"].lower(): f for f in files}
        for key, candidates in WANT.items():
            for cand in candidates:
                if cand in by_name:
                    urls[key] = by_name[cand]["download_url"]
                    break
        if "model" in urls and "expr" in urls:
            log(WORKER, f"article {art}: resolved keys = {sorted(urls)}")
            return urls
        log(WORKER, f"article {art}: missing essentials, trying next")
    return {}


def fetch_all() -> dict[str, Path]:
    urls = resolve_article_files()
    paths: dict[str, Path] = {}
    for key, url in urls.items():
        dest = RAW / f"{key}.csv"
        try:
            paths[key] = download(url, dest)
        except Exception as exc:
            log(WORKER, f"FAIL {key}: {exc}")
    return paths


def main() -> None:
    write_status(WORKER, "running", started_at=time.time())
    paths = fetch_all()
    if not all(k in paths for k in ("model", "expr")):
        raise RuntimeError("essential downloads (model/expr) failed")

    log(WORKER, "loading Model.csv")
    model = pd.read_csv(paths["model"])
    cols_lower = {c.lower(): c for c in model.columns}
    lineage_col = next(
        (cols_lower[k] for k in ("oncotreelineage", "oncotree_lineage", "lineage", "primarydisease")
         if k in cols_lower),
        None,
    )
    if lineage_col is None:
        raise RuntimeError(f"can't find lineage column in {model.columns.tolist()[:20]}")
    log(WORKER, f"using lineage column: {lineage_col}")

    thyroid_models = model[model[lineage_col].astype(str).str.contains("hyroid", na=False)].copy()
    log(WORKER, f"thyroid cell lines in DepMap: n={len(thyroid_models)}")
    thyroid_models.to_csv(OUT / "thyroid_models.tsv", sep="\t", index=False)

    log(WORKER, "loading expression (this is the big file ~1.5GB)")
    expr = pd.read_csv(paths["expr"], index_col=0)
    log(WORKER, f"expr shape={expr.shape}")
    # Columns are like 'TPO (7173)' — split off symbol
    new_cols = [c.split(" ", 1)[0] for c in expr.columns]
    expr.columns = new_cols
    thyroid_lines = [m for m in thyroid_models["ModelID"] if m in expr.index]
    log(WORKER, f"thyroid lines with expression: n={len(thyroid_lines)}")
    if not thyroid_lines:
        raise RuntimeError("no thyroid lines mapped to expression matrix")

    sub = expr.loc[thyroid_lines]
    gene_map = resolve_panel(sub.columns)
    missing = [g for g in PANEL_8 if g not in gene_map]
    log(WORKER, f"resolved {len(gene_map)}/{len(PANEL_8)} panel genes; missing={missing}")

    sub = sub.loc[:, list(gene_map.values())].dropna(how="any")
    scores = score_panel(sub, gene_map)
    dm = kmeans_dm(scores)
    out = pd.concat(
        [thyroid_models.set_index("ModelID").loc[sub.index, [lineage_col, "OncotreeSubtype" if "OncotreeSubtype" in thyroid_models.columns else lineage_col, "CellLineName" if "CellLineName" in thyroid_models.columns else lineage_col]], scores, dm],
        axis=1,
    )
    out.to_csv(OUT / "thyroid_cell_dm_calls.tsv", sep="\t")
    log(WORKER, f"thyroid DM calls: DM1={int((dm == 'DM1').sum())} DM2={int((dm == 'DM2').sum())}")

    # ---- CRISPR essentiality
    crispr_summary = None
    if "crispr" in paths:
        log(WORKER, "loading CRISPR essentiality (1.6GB)")
        crispr = pd.read_csv(paths["crispr"], index_col=0)
        crispr.columns = [c.split(" ", 1)[0] for c in crispr.columns]
        common = [m for m in sub.index if m in crispr.index]
        log(WORKER, f"crispr-overlap thyroid lines: {len(common)}")
        crispr_panel = crispr.loc[common, [g for g in gene_map.values() if g in crispr.columns]]
        crispr_panel.to_csv(OUT / "thyroid_crispr_panel.tsv", sep="\t")
        # DM1 vs DM2 differential essentiality (mean diff + Mann-Whitney)
        from scipy.stats import mannwhitneyu

        diffs = []
        for g in crispr_panel.columns:
            dm1_lines = [c for c in common if dm.get(c) == "DM1"]
            dm2_lines = [c for c in common if dm.get(c) == "DM2"]
            if len(dm1_lines) < 2 or len(dm2_lines) < 2:
                continue
            v1 = crispr_panel.loc[dm1_lines, g].dropna()
            v2 = crispr_panel.loc[dm2_lines, g].dropna()
            if len(v1) < 2 or len(v2) < 2:
                continue
            stat, p = mannwhitneyu(v1, v2, alternative="two-sided")
            diffs.append({
                "gene": g,
                "n_dm1": len(v1),
                "n_dm2": len(v2),
                "mean_dm1": float(v1.mean()),
                "mean_dm2": float(v2.mean()),
                "delta": float(v1.mean() - v2.mean()),
                "mw_p": float(p),
            })
        crispr_summary = pd.DataFrame(diffs)
        crispr_summary.to_csv(OUT / "thyroid_crispr_dm1_vs_dm2.tsv", sep="\t", index=False)
        log(WORKER, f"crispr diff table: {len(crispr_summary)} genes")

    # ---- PRISM drug screen
    prism_summary = None
    if "prism" in paths:
        try:
            log(WORKER, "loading PRISM drug screen")
            prism = pd.read_csv(paths["prism"], index_col=0)
            common = [m for m in sub.index if m in prism.index]
            log(WORKER, f"prism-overlap thyroid lines: {len(common)}")
            if common:
                # Find drugs of interest by name match
                target_substrings = (
                    "decitabine", "azacitidine", "5-aza", "selpercatinib",
                    "pralsetinib", "vandetanib", "lenvatinib", "sorafenib",
                    "cabozantinib", "vorinostat", "panobinostat",
                )
                hit_cols = [
                    c for c in prism.columns
                    if any(s in c.lower() for s in target_substrings)
                ]
                log(WORKER, f"PRISM target hits: {len(hit_cols)} columns matched")
                if hit_cols:
                    prism_sub = prism.loc[common, hit_cols]
                    prism_sub.to_csv(OUT / "thyroid_prism_panel.tsv", sep="\t")
                    rows = []
                    from scipy.stats import mannwhitneyu

                    dm1_lines = [c for c in common if dm.get(c) == "DM1"]
                    dm2_lines = [c for c in common if dm.get(c) == "DM2"]
                    for c in hit_cols:
                        v1 = prism_sub.loc[dm1_lines, c].dropna()
                        v2 = prism_sub.loc[dm2_lines, c].dropna()
                        if len(v1) < 2 or len(v2) < 2:
                            continue
                        stat, p = mannwhitneyu(v1, v2, alternative="two-sided")
                        rows.append({
                            "drug": c,
                            "n_dm1": len(v1),
                            "n_dm2": len(v2),
                            "mean_dm1": float(v1.mean()),
                            "mean_dm2": float(v2.mean()),
                            "delta": float(v1.mean() - v2.mean()),
                            "mw_p": float(p),
                        })
                    prism_summary = pd.DataFrame(rows)
                    prism_summary.to_csv(OUT / "thyroid_prism_dm1_vs_dm2.tsv", sep="\t", index=False)
                    log(WORKER, f"prism diff table: {len(prism_summary)} drugs")
        except Exception as exc:
            log(WORKER, f"PRISM block failed (non-fatal): {exc}")

    summary = {
        "n_thyroid_models": int(len(thyroid_models)),
        "n_thyroid_with_expr": int(len(sub)),
        "n_dm1": int((dm == "DM1").sum()),
        "n_dm2": int((dm == "DM2").sum()),
        "n_crispr_genes": int(len(crispr_summary)) if crispr_summary is not None else 0,
        "n_prism_drugs": int(len(prism_summary)) if prism_summary is not None else 0,
        "panel_gene_map": gene_map,
        "missing_panel_genes": missing,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    write_status(WORKER, "done", **{k: v for k, v in summary.items() if not isinstance(v, dict)}, finished_at=time.time())
    log(WORKER, "DONE")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(WORKER, f"FATAL: {exc}")
        write_status(WORKER, "error", error=str(exc))
        raise
