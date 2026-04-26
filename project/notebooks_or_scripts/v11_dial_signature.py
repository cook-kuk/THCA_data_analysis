"""v11 Task 2 — DIAL-flip molecular signature.

Compute three gene rankings on THCA harmonized data:
  1. label-preserving direction (BRAF vs RAS, pre + post ComBat agreement)
  2. label-flip contributing direction (pre vs post disagreement)
  3. batch-aligned direction (TCGA vs GSE27155 mean difference)

Intersect the TOP 100 flip-contributing and TOP 100 batch-aligned lists to
produce the "confounded" signature, then enrich it via gseapy.enrichr.

Usage (logged):
    python -u notebooks_or_scripts/v11_dial_signature.py \
        2>&1 | tee logs/v11_dial_signature.log
"""

from __future__ import annotations

import os
import sys
import time
import traceback
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

PROJECT = Path("/opt/thyroid-dash/project")
DATA = PROJECT / "data_processed" / "v5_cross_cancer" / "THCA"
RES = PROJECT / "results" / "v11_novel_pathways" / "modules"
LOGS = PROJECT / "logs"
BIOMARKER_FULL = PROJECT / "results" / "tables" / "biomarker_de_full.tsv"

DRUGGABLE = [
    "TACSTD2", "TMPRSS4", "PLEKHA6", "CYP1B1",
    "LDLR", "GABRB2", "B3GNT3", "PTPRE",
]

TOP_K = 100
VARIANCE_FALLBACK_K = 5000


def log(msg: str) -> None:
    print(msg, flush=True)


def _to_unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if not np.isfinite(n) or n == 0:
        return v.astype(float)
    return (v / n).astype(float)


def _zscore_columns(X: np.ndarray) -> np.ndarray:
    mu = X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, ddof=0, keepdims=True)
    sd = np.where(sd <= 1e-12, 1.0, sd)
    return (X - mu) / sd


def load_inputs() -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str]]:
    log(f"[load] X from {DATA/'X_combined.npz'}")
    d = np.load(DATA / "X_combined.npz")
    X = np.asarray(d["X"], dtype=float)

    Y = pd.read_csv(DATA / "Y.tsv", sep="\t", header=None).iloc[:, -1].astype(str).values
    B = pd.read_csv(DATA / "B.tsv", sep="\t", header=None).iloc[:, -1].astype(str).values
    genes = [g.strip() for g in open(DATA / "shared_genes.txt").read().splitlines() if g.strip()]

    assert X.shape == (len(Y), len(genes)), f"mismatch: X{X.shape}, Y{len(Y)}, genes{len(genes)}"
    assert X.shape[0] == len(B), f"mismatch: X rows {X.shape[0]} vs B {len(B)}"

    log(f"[load] X={X.shape}, Y vc={dict(pd.Series(Y).value_counts())}, B vc={dict(pd.Series(B).value_counts())}")
    log(f"[load] genes: {len(genes)} (first={genes[:3]}, last={genes[-3:]})")
    return X, Y, B, genes


def combat_preserve(
    X: np.ndarray,
    Y: np.ndarray,
    B: np.ndarray,
    gene_list: List[str],
) -> Tuple[np.ndarray, List[str], bool]:
    """Return (X_post, genes_used, degraded).

    Mirrors v5p1_common._combat_preserve but with a variance-top-5000 fallback
    if the full-dimensional ComBat run blows up.
    """
    from inmoose.pycombat import pycombat_norm

    def _run(X_in: np.ndarray) -> np.ndarray:
        df = pd.DataFrame(X_in.T)
        covar = pd.get_dummies(pd.Series(Y).astype(str)).astype(float)
        force_no_mod = False
        for b in np.unique(B):
            sel = np.asarray(B) == b
            if len(np.unique(Y[sel])) < 2:
                force_no_mod = True
                break
        if force_no_mod:
            out = pycombat_norm(df.values, batch=list(B), covar_mod=None, par_prior=True)
        else:
            try:
                out = pycombat_norm(df.values, batch=list(B), covar_mod=covar.values, par_prior=True)
            except Exception as e:
                log(f"[combat] covar fit failed ({type(e).__name__}: {e}); retrying without covar")
                out = pycombat_norm(df.values, batch=list(B), covar_mod=None, par_prior=True)
        if hasattr(out, "values"):
            out = out.values
        return np.asarray(out).T

    log(f"[combat] attempting full-dimensional ComBat on X{X.shape}")
    t0 = time.time()
    try:
        X_post = _run(X)
        log(f"[combat] full-dim ComBat OK in {time.time()-t0:.1f}s")
        return X_post, gene_list, False
    except Exception as e:
        log(f"[combat] full-dim failed: {type(e).__name__}: {e}")
        traceback.print_exc()

    # Fallback: top-variance subset
    var = X.var(axis=0)
    keep = np.argsort(-var)[:VARIANCE_FALLBACK_K]
    keep = np.sort(keep)
    sub_genes = [gene_list[i] for i in keep]
    log(f"[combat] falling back to variance-top-{VARIANCE_FALLBACK_K} subset")
    t0 = time.time()
    X_post_sub = _run(X[:, keep])
    log(f"[combat] subset ComBat OK in {time.time()-t0:.1f}s, shape={X_post_sub.shape}")
    return X_post_sub, sub_genes, True


def mean_diff_direction(
    X: np.ndarray, group: np.ndarray, pos: str, neg: str
) -> np.ndarray:
    mpos = X[group == pos].mean(axis=0)
    mneg = X[group == neg].mean(axis=0)
    return _to_unit(mpos - mneg)


def save_ranked(
    path: Path,
    genes: List[str],
    w: np.ndarray,
    value_col: str,
    top_k: int,
) -> pd.DataFrame:
    order = np.argsort(-np.abs(w))
    df = pd.DataFrame({
        "gene": [genes[i] for i in order[:top_k]],
        value_col: w[order[:top_k]],
        "rank": np.arange(1, min(top_k, len(order)) + 1),
    })
    df.to_csv(path, sep="\t", index=False)
    log(f"[write] {path} ({len(df)} rows)")
    return df


def load_2773_biomarkers() -> set:
    try:
        df = pd.read_csv(BIOMARKER_FULL, sep="\t", low_memory=False)
        if "is_novel_validated" in df.columns:
            s = set(df[df["is_novel_validated"] == True]["gene"].dropna().astype(str).unique())
            log(f"[ref] loaded is_novel_validated set n={len(s)} from {BIOMARKER_FULL}")
            return s
    except Exception as e:
        log(f"[ref] could not load 2773 biomarker set: {e}")
    return set()


def run_enrichment(gene_list: List[str], out_path: Path) -> bool:
    """Run gseapy.enrichr; fall back to local hypergeometric. Return True if online."""
    libs = [
        "MSigDB_Hallmark_2020",
        "Reactome_2022",
        "GO_Biological_Process_2023",
        "KEGG_2021_Human",
    ]
    import gseapy
    records: List[dict] = []
    online = True
    for lib in libs:
        try:
            log(f"[enrichr] {lib} (online) on {len(gene_list)} genes")
            enr = gseapy.enrichr(
                gene_list=list(gene_list),
                gene_sets=lib,
                organism="Human",
                outdir=None,
                cutoff=1.0,
            )
            res = enr.results.copy()
            res["library"] = lib
            records.append(res)
        except Exception as e:
            log(f"[enrichr] {lib} online failed: {type(e).__name__}: {e}")
            online = False
            try:
                local = gseapy.get_library(name=lib, organism="Human")
                from scipy.stats import hypergeom
                # Universe size: union of all genes in the library (a pragmatic estimate).
                universe = set()
                for gs in local.values():
                    universe.update(gs)
                N = len(universe)
                sig = set(gene_list) & universe
                n = len(sig)
                rows = []
                for term, genes in local.items():
                    gs = set(genes)
                    K = len(gs)
                    overlap = sig & gs
                    k = len(overlap)
                    if k == 0:
                        continue
                    # P(X >= k) where X ~ Hypergeom(N, K, n)
                    p = hypergeom.sf(k - 1, N, K, n)
                    rows.append({
                        "Term": term,
                        "Overlap": f"{k}/{K}",
                        "P-value": p,
                        "Adjusted P-value": np.nan,  # BH below
                        "Genes": ";".join(sorted(overlap)),
                        "library": lib,
                    })
                rows_df = pd.DataFrame(rows).sort_values("P-value")
                # BH
                m = len(rows_df)
                if m:
                    rows_df = rows_df.reset_index(drop=True)
                    ranks = np.arange(1, m + 1)
                    bh = (rows_df["P-value"].values * m / ranks)
                    bh = np.minimum.accumulate(bh[::-1])[::-1]
                    rows_df["Adjusted P-value"] = np.clip(bh, 0, 1)
                records.append(rows_df)
                log(f"[enrichr] {lib} local fallback OK: {len(rows_df)} terms")
            except Exception as e2:
                log(f"[enrichr] {lib} local fallback failed: {type(e2).__name__}: {e2}")

    if not records:
        log(f"[enrichr] no results for any library")
        pd.DataFrame(columns=["term", "library", "overlap", "p_value", "adj_p", "genes"]).to_csv(
            out_path, sep="\t", index=False
        )
        return online

    all_df = pd.concat(records, ignore_index=True, sort=False)
    def _col(df, *names):
        for n in names:
            if n in df.columns:
                return df[n]
        return pd.Series([np.nan] * len(df))

    out = pd.DataFrame({
        "term": _col(all_df, "Term"),
        "library": _col(all_df, "library"),
        "overlap": _col(all_df, "Overlap"),
        "p_value": _col(all_df, "P-value"),
        "adj_p": _col(all_df, "Adjusted P-value"),
        "genes": _col(all_df, "Genes"),
    }).sort_values(["p_value"]).reset_index(drop=True)
    out.to_csv(out_path, sep="\t", index=False)
    log(f"[write] {out_path} ({len(out)} rows)")
    return online


def main() -> int:
    RES.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    log(f"[run] v11 Task 2 — DIAL signature; cwd={os.getcwd()}")

    X_raw, Y, B, genes = load_inputs()

    # Step 2: z-score per gene -> X_pre
    X_pre = _zscore_columns(X_raw)
    log(f"[pre] X_pre mean={X_pre.mean():.4e} std={X_pre.std():.4e}")

    # Step 3: covariate-preserving ComBat -> X_post (with fallback)
    X_post, genes_post, degraded = combat_preserve(X_pre, Y, B, genes)
    if degraded:
        log(f"[degrade] ComBat required variance-top-{VARIANCE_FALLBACK_K} fallback; "
            f"label-preserving + flip directions will be restricted to {len(genes_post)} genes")
        # Intersect to the reduced gene set for w_Y_pre/post comparison
        pos_map = {g: i for i, g in enumerate(genes)}
        idx = np.array([pos_map[g] for g in genes_post], dtype=int)
        X_pre_for_label = X_pre[:, idx]
        genes_label = genes_post
    else:
        X_pre_for_label = X_pre
        genes_label = genes

    # Step 4: label directions
    w_Y_pre = mean_diff_direction(X_pre_for_label, Y, pos="BRAF", neg="RAS")
    w_Y_post = mean_diff_direction(X_post, Y, pos="BRAF", neg="RAS")
    log(f"[dir] ||w_Y_pre||={np.linalg.norm(w_Y_pre):.4f} ||w_Y_post||={np.linalg.norm(w_Y_post):.4f}")

    # Step 5: preserve + flip
    w_preserve = w_Y_post + w_Y_pre
    w_flip = w_Y_pre - w_Y_post
    log(f"[dir] preserve stats: |mean|={np.abs(w_preserve).mean():.4e} max={np.abs(w_preserve).max():.4e}")
    log(f"[dir] flip stats: |mean|={np.abs(w_flip).mean():.4e} max={np.abs(w_flip).max():.4e}")

    # Step 6: batch direction (cohort mean difference over full X_pre)
    tcga_label = "TCGA-THCA" if "TCGA-THCA" in set(B) else sorted(set(B))[0]
    gse_label = [b for b in set(B) if b != tcga_label][0]
    log(f"[batch] TCGA label={tcga_label} ; other={gse_label}")
    w_B = mean_diff_direction(X_pre, B, pos=tcga_label, neg=gse_label)
    log(f"[batch] ||w_B||={np.linalg.norm(w_B):.4f} |mean|={np.abs(w_B).mean():.4e}")

    # Step 7: save top-100 lists
    df_preserve = save_ranked(
        RES / "dial_label_preserving_genes.tsv", genes_label, w_preserve, "w_preserve", TOP_K
    )
    df_flip = save_ranked(
        RES / "dial_flip_contributing_genes.tsv", genes_label, w_flip, "w_flip", TOP_K
    )
    df_batch = save_ranked(
        RES / "dial_batch_aligned_genes.tsv", genes, w_B, "w_B", TOP_K
    )

    # Step 8: intersection (flip & batch)
    inter = sorted(set(df_flip["gene"]) & set(df_batch["gene"]))
    log(f"[intersect] top-100 flip ∩ top-100 batch: n={len(inter)}")

    # Diagnostic: overlap at wider rank thresholds (helpful when top-100 is empty)
    abs_flip = np.abs(w_flip)
    abs_b = np.abs(w_B)
    # Need to index flip list on genes_label vs batch on full genes; align via name
    flip_order = np.argsort(-abs_flip)
    batch_order = np.argsort(-abs_b)
    for k in [200, 500, 1000, 2000]:
        s_flip = set(genes_label[i] for i in flip_order[:k])
        s_b = set(genes[i] for i in batch_order[:k])
        log(f"[intersect] top-{k} flip ∩ top-{k} batch: n={len(s_flip & s_b)}")

    # Build per-gene maps for annotation
    label_idx = {g: i for i, g in enumerate(genes_label)}
    full_idx = {g: i for i, g in enumerate(genes)}
    bio_2773 = load_2773_biomarkers()

    rows = []
    for g in inter:
        li = label_idx.get(g)
        fi = full_idx.get(g)
        rows.append({
            "gene": g,
            "w_flip": float(w_flip[li]) if li is not None else np.nan,
            "w_B": float(w_B[fi]) if fi is not None else np.nan,
            "w_preserve": float(w_preserve[li]) if li is not None else np.nan,
            "is_in_2773_biomarker": bool(g in bio_2773),
            "is_druggable_target": bool(g in set(DRUGGABLE)),
        })
    if rows:
        conf_df = pd.DataFrame(rows)
        conf_df = conf_df.reindex(conf_df["w_flip"].abs().sort_values(ascending=False).index).reset_index(drop=True)
    else:
        conf_df = pd.DataFrame(
            columns=["gene", "w_flip", "w_B", "w_preserve", "is_in_2773_biomarker", "is_druggable_target"]
        )
    conf_path = RES / "dial_confounded_signature.tsv"
    conf_df.to_csv(conf_path, sep="\t", index=False)
    log(f"[write] {conf_path} ({len(conf_df)} rows, "
        f"druggable_overlap={int(conf_df['is_druggable_target'].sum())}, "
        f"in_2773={int(conf_df['is_in_2773_biomarker'].sum())})")

    # Step 9: enrichment on intersection (on non-empty list)
    enr_path = RES / "confounded_enrichment.tsv"
    if len(conf_df) > 0:
        run_enrichment(conf_df["gene"].tolist(), enr_path)
    else:
        log(f"[enrichr] skipped (empty intersection)")
        pd.DataFrame(columns=["term", "library", "overlap", "p_value", "adj_p", "genes"]).to_csv(
            enr_path, sep="\t", index=False
        )

    # Short summary
    log("[summary] ---")
    log(f"[summary] preserve list:    {len(df_preserve)} rows -> {RES/'dial_label_preserving_genes.tsv'}")
    log(f"[summary] flip list:        {len(df_flip)} rows -> {RES/'dial_flip_contributing_genes.tsv'}")
    log(f"[summary] batch list:       {len(df_batch)} rows -> {RES/'dial_batch_aligned_genes.tsv'}")
    log(f"[summary] intersection:     n={len(conf_df)}")
    head10 = conf_df.head(10)[["gene", "w_flip", "w_B"]]
    log(f"[summary] top 10 confounded:\n{head10.to_string(index=False)}")
    try:
        enr_df = pd.read_csv(enr_path, sep="\t")
        if len(enr_df):
            top = enr_df.iloc[0]
            log(f"[summary] top enrichment term: '{top['term']}' ({top['library']}) "
                f"p={top['p_value']:.3e} adj={top['adj_p']:.3e}")
        else:
            log(f"[summary] no enrichment rows")
    except Exception as e:
        log(f"[summary] enrichment read failed: {e}")
    log(f"[summary] degraded (combat variance-top-{VARIANCE_FALLBACK_K} fallback): {degraded}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
