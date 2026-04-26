#!/usr/bin/env python3
"""v8 Task 3 - Pathway-level DIAL for THCA BRAF vs RAS.

Question: does the v5.1 post-ComBat label flip (DIAL up to 0.494 at 3000-gene
level) survive pathway-level aggregation? If pathway-level DIAL > 0.3 for any
classifier, the finding is robust and not a gene-level fluke.

Strategy:
  1. Load harmonized THCA X (392 samples x 11710 HGNC symbols, float32 log2-TPM).
  2. Load MSigDB Hallmark 2020 (50 pathways) via gseapy.get_library, with a
     plain-text Enrichr fallback.
  3. Compute sample x pathway activity matrix with ssGSEA (Barbie 2009).
  4. Apply pycombat_norm preserving Y (BRAF vs RAS) as covariate, B as batch
     (identical protocol to v5p1_common._combat_preserve).
  5. Reuse v5p1_common.compute_dial + LeaveOneGroupOut with all 5 classifiers.
  6. Write v8_pathway_dial.tsv (same schema as v5p1_dial_THCA.tsv + level col).
  7. Write v8_pathway_vs_gene.tsv comparing gene-level vs pathway-level DIAL.
"""
from __future__ import annotations

import os
import sys
import time
import json
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

PROJECT = Path("/opt/thyroid-dash/project")
sys.path.insert(0, str(PROJECT / "notebooks_or_scripts"))

from v5p1_common import (  # noqa: E402
    compute_dial,
    get_classifier_factories,
    log_line,
)

DATA = PROJECT / "data_processed" / "v5_cross_cancer" / "THCA"
RES = PROJECT / "results" / "v8_statgen"
RPT = PROJECT / "reports" / "v8"
LOGS = PROJECT / "logs"
CACHE = PROJECT / "data_raw" / "msigdb_cache"
for p in [RES, RPT, LOGS, CACHE]:
    p.mkdir(parents=True, exist_ok=True)

LOG = LOGS / "v8_pathway.log"
GENE_DIAL_TSV = PROJECT / "results" / "v5" / "v5p1_dial_THCA.tsv"


# ------------------------------------------------------------------
# 1. Hallmark loader with Enrichr fallback
# ------------------------------------------------------------------
def load_hallmark() -> dict:
    """Return dict pathway_name -> list of HGNC symbols.

    Tries gseapy.get_library first; on failure falls back to the Enrichr
    plain-text endpoint and caches to disk.
    """
    cache_json = CACHE / "MSigDB_Hallmark_2020.json"
    if cache_json.exists():
        log_line(LOG, f"[hallmark] loaded from cache {cache_json}")
        return json.loads(cache_json.read_text())

    try:
        import gseapy
        lib = gseapy.get_library("MSigDB_Hallmark_2020")
        if isinstance(lib, dict) and len(lib) >= 40:
            cache_json.write_text(json.dumps(lib))
            log_line(LOG, f"[hallmark] got {len(lib)} sets via gseapy.get_library")
            return lib
    except Exception as e:
        log_line(LOG, f"[hallmark] gseapy.get_library failed: {type(e).__name__}: {e}")

    url = (
        "https://maayanlab.cloud/Enrichr/geneSetLibrary"
        "?mode=text&libraryName=MSigDB_Hallmark_2020"
    )
    log_line(LOG, f"[hallmark] fallback: fetching {url}")
    with urllib.request.urlopen(url, timeout=120) as r:
        text = r.read().decode("utf-8", errors="ignore")

    lib: dict = {}
    for line in text.strip().split("\n"):
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        name = fields[0].strip()
        # Enrichr format is name \t blank \t gene \t gene ...; drop empties.
        genes = [g.strip() for g in fields[1:] if g.strip()]
        if name and genes:
            lib[name] = genes
    cache_json.write_text(json.dumps(lib))
    log_line(LOG, f"[hallmark] parsed {len(lib)} pathways from Enrichr text")
    return lib


# ------------------------------------------------------------------
# 2. ssGSEA with numpy-rank fallback
# ------------------------------------------------------------------
def compute_ssgsea(
    X: np.ndarray, genes: list, sample_names: list, gene_sets: dict
) -> pd.DataFrame:
    """Return samples x pathways DataFrame."""
    expr = pd.DataFrame(X.T, index=genes, columns=sample_names)
    expr = expr[~expr.index.duplicated(keep="first")]
    log_line(
        LOG,
        f"[ssgsea] expr matrix: {expr.shape[0]} genes x {expr.shape[1]} samples",
    )

    try:
        import gseapy

        t0 = time.time()
        res = gseapy.ssgsea(
            data=expr,
            gene_sets=gene_sets,
            outdir=None,
            sample_norm_method="rank",
            min_size=5,
            max_size=1000,
            permutation_num=0,
            threads=max(os.cpu_count() or 4, 4),
            no_plot=True,
            verbose=False,
        )
        dt = time.time() - t0
        # gseapy returns res.res2d with Term/ES/NES/sample columns. The wide
        # matrix is at res.resultsOnSamples (dict pathway -> series sample->es).
        # Build samples x pathways from res2d by pivoting.
        df = res.res2d.copy()
        log_line(LOG, f"[ssgsea] gseapy done in {dt:.1f}s; res2d cols={list(df.columns)[:6]} rows={len(df)}")
        # Expect cols: Name (sample), Term, ES, NES, ...
        score_col = "NES" if "NES" in df.columns else ("ES" if "ES" in df.columns else None)
        sample_col = "Name" if "Name" in df.columns else ("Sample" if "Sample" in df.columns else None)
        if score_col is None or sample_col is None or "Term" not in df.columns:
            raise RuntimeError(f"unexpected ssgsea columns: {list(df.columns)}")
        mat = df.pivot_table(
            index=sample_col, columns="Term", values=score_col, aggfunc="mean"
        )
        # Reindex to full sample list in original order
        mat = mat.reindex(sample_names)
        # Drop pathways that are all-NaN; forward-fill rare NaNs with col mean
        mat = mat.dropna(axis=1, how="all")
        mat = mat.fillna(mat.mean(axis=0))
        log_line(LOG, f"[ssgsea] activity matrix: {mat.shape}")
        return mat
    except Exception as e:
        log_line(
            LOG,
            f"[ssgsea] gseapy.ssgsea FAILED ({type(e).__name__}: {e}); "
            "degrading to numpy rank-mean fallback",
        )

    # Fallback: rank-normalize each sample, then per pathway take mean rank of
    # pathway genes. This is the simplest GSVA-like surrogate.
    from scipy.stats import rankdata

    ranks = np.zeros_like(expr.values, dtype=np.float32)
    for j in range(expr.shape[1]):
        ranks[:, j] = rankdata(expr.values[:, j], method="average")
    ranks = ranks / ranks.shape[0]  # normalize to (0,1]
    rank_df = pd.DataFrame(ranks, index=expr.index, columns=expr.columns)
    gene_set = set(expr.index)
    out_cols = {}
    for pw, gs in gene_sets.items():
        hits = [g for g in gs if g in gene_set]
        if len(hits) < 5:
            continue
        out_cols[pw] = rank_df.loc[hits].mean(axis=0).values
    mat = pd.DataFrame(out_cols, index=sample_names)
    log_line(LOG, f"[fallback] rank-mean matrix: {mat.shape}")
    return mat


# ------------------------------------------------------------------
# 3. Main DIAL loop
# ------------------------------------------------------------------
def main():
    log_line(LOG, "=" * 70)
    log_line(LOG, "v8 Task 3 pathway DIAL starting")

    # Load harmonized THCA
    X = np.load(DATA / "X_combined.npz")["X"].astype(np.float32)
    Y = np.loadtxt(DATA / "Y.tsv", dtype=str)
    B = np.loadtxt(DATA / "B.tsv", dtype=str)
    with open(DATA / "shared_genes.txt") as fh:
        genes = [l.strip() for l in fh if l.strip()]
    with open(DATA / "sample_names.txt") as fh:
        sample_names = [l.strip() for l in fh if l.strip()]
    assert X.shape == (len(sample_names), len(genes)), (
        f"X shape {X.shape} vs samples {len(sample_names)} genes {len(genes)}"
    )
    log_line(
        LOG,
        f"[load] X={X.shape} Y uniq={dict(zip(*np.unique(Y, return_counts=True)))} "
        f"B uniq={dict(zip(*np.unique(B, return_counts=True)))}",
    )

    # Load hallmark and compute pathway activity
    gene_sets = load_hallmark()
    log_line(LOG, f"[hallmark] {len(gene_sets)} pathways; median genes={int(np.median([len(v) for v in gene_sets.values()]))}")

    # Estimate ssGSEA runtime heuristically: 392 samples x 50 pathways with
    # 11710 gene background. Typical gseapy throughput ~1-2 samples/sec per
    # thread; with 8 threads we expect well under 30 minutes. Document anyway.
    log_line(LOG, "[plan] running gseapy.ssgsea with permutation_num=0 (scoring only)")

    t0 = time.time()
    pw_df = compute_ssgsea(X, genes, sample_names, gene_sets)
    log_line(LOG, f"[ssgsea] total wall: {time.time() - t0:.1f}s")

    # Persist activity matrix for audit
    act_tsv = RES / "v8_pathway_activity_THCA.tsv"
    pw_df.to_csv(act_tsv, sep="\t")
    log_line(LOG, f"[save] activity matrix -> {act_tsv}")

    X_pw = pw_df.values.astype(np.float32)
    n_samples, n_pw = X_pw.shape

    # CV splits: Leave-One-Dataset-Out on B
    from sklearn.model_selection import LeaveOneGroupOut

    classes_sorted = sorted(np.unique(Y).tolist())
    Y_bin = (np.asarray(Y) == classes_sorted[0]).astype(int)
    logo = LeaveOneGroupOut()
    splits = list(logo.split(np.arange(n_samples), Y_bin, groups=B))
    log_line(LOG, f"[cv] LODO splits: {len(splits)} folds")

    facs = get_classifier_factories()
    clf_order = ["LogReg_l2", "RandomForest", "XGBoost", "LogReg_elasticnet", "GradientBoosting"]
    if "XGBoost" not in facs:
        clf_order = [c if c != "XGBoost" else "HistGB" for c in clf_order]
    log_line(LOG, f"[clf] running {clf_order}")

    rows = []
    for clf_name in clf_order:
        if clf_name not in facs:
            log_line(LOG, f"[skip] {clf_name} not available")
            continue
        log_line(LOG, f"[run] {clf_name}")
        tt = time.time()
        res = compute_dial(X_pw, Y, B, facs[clf_name], splits)
        dt = time.time() - tt
        res.update(
            dict(
                cancer="THCA",
                classifier=clf_name,
                n_samples=n_samples,
                n_genes=n_pw,  # keep column name for schema-compat; value = # pathways
                semi_synthetic=False,
                seconds=round(dt, 2),
                is_label_flip=bool(res["auc_post"] < 0.5 and res["dial"] > 0.1),
                level="pathway",
            )
        )
        rows.append(res)
        log_line(
            LOG,
            f"[done] {clf_name} auc_pre={res['auc_pre']:.3f} auc_post={res['auc_post']:.3f} "
            f"dial={res['dial']:.3f} interp={res['interpretation']} ({dt:.1f}s)",
        )

    out_df = pd.DataFrame(rows)
    cols_order = [
        "auc_pre", "auc_flip_pre", "auc_post", "auc_flip_post", "dial",
        "batch_identifiability_post", "interpretation", "cancer", "classifier",
        "n_samples", "n_genes", "semi_synthetic", "seconds", "is_label_flip",
        "level",
    ]
    cols_order = [c for c in cols_order if c in out_df.columns]
    out_df = out_df[cols_order]
    out_tsv = RES / "v8_pathway_dial.tsv"
    out_df.to_csv(out_tsv, sep="\t", index=False)
    log_line(LOG, f"[save] pathway DIAL -> {out_tsv}")

    # ------------------------------------------------------------------
    # 4. Gene-vs-pathway comparison
    # ------------------------------------------------------------------
    gene_df = pd.read_csv(GENE_DIAL_TSV, sep="\t")
    gene_df = gene_df.rename(
        columns={"dial": "dial_gene", "auc_post": "auc_post_gene"}
    )[["classifier", "dial_gene", "auc_post_gene"]]
    pw_cmp = out_df.rename(
        columns={"dial": "dial_pathway", "auc_post": "auc_post_pathway"}
    )[["classifier", "dial_pathway", "auc_post_pathway"]]
    cmp_df = gene_df.merge(pw_cmp, on="classifier", how="outer")
    cmp_df = cmp_df[
        ["classifier", "dial_gene", "dial_pathway", "auc_post_gene", "auc_post_pathway"]
    ]
    cmp_tsv = RES / "v8_pathway_vs_gene.tsv"
    cmp_df.to_csv(cmp_tsv, sep="\t", index=False)
    log_line(LOG, f"[save] gene-vs-pathway -> {cmp_tsv}")
    log_line(LOG, "[summary]\n" + cmp_df.to_string(index=False))

    log_line(LOG, "v8 Task 3 pathway DIAL DONE")
    print("\n=== v8 Task 3 output files ===")
    print(out_tsv)
    print(cmp_tsv)
    print(act_tsv)


if __name__ == "__main__":
    main()
