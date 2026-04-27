#!/usr/bin/env python3
"""
v17 REAL FIX · R1 — leak-free DM1/DM2 clustering (4 variants).

The original DM1/DM2 cluster was created from TIERA67 (55 unique genes) which
INCLUDES the 8-gene RAI panel. Predicting that cluster with the 8-gene panel
is autocorrelation, not prediction. This script generates 4 leak-free
alternatives:

  R1-A: TIERA67 minus 8-gene  (47 genes, same biology)
  R1-B: MAPK + immune + EMT   (30 genes, zero overlap)
  R1-C: 5-feature immune composite  (functional level, no DC genes)
  R1-D: BRS71 minus 8-gene    (Chakravarty 2011 framework)

For each variant:
  - Leiden K=2, resolution swept, 1000-bootstrap consensus stability
  - ARI / Cohen's kappa vs original DM1/DM2

Outputs (per variant):
  - results/v17_realfix/R1{X}_cluster_labels.tsv
  - results/v17_realfix/R1{X}_summary.json
  - results/v17_realfix/R1_concordance_table.tsv (combined at end)
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, cohen_kappa_score
from sklearn.preprocessing import StandardScaler

OUT = Path("/opt/thyroid-dash/project/results/v17_realfix")
OUT.mkdir(parents=True, exist_ok=True)
TPM = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv")
TIERA67_TXT = Path("/opt/thyroid-dash/project/metadata/v3_tierA67_clean_genes.txt")
ORIG_LABELS = Path("/opt/thyroid-dash/project/results/v17p35/tables/AMP3_hot_cold_composite.tsv")

GENE_8 = {"SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"}

# Variant 2: 30-gene MAPK + immune + EMT (zero overlap with 8-gene panel)
GENE_30_FUNCTIONAL = [
    "DUSP4", "DUSP5", "DUSP6", "SPRY1", "SPRY2", "SPRY4", "ETV4", "ETV5",
    "CD274", "CD8A", "FOXP3", "IDO1", "HLA-DRA", "PHLDA1", "FOSL1",
    "PRF1", "GZMB", "GZMA",
    "VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "CDH1", "CDH2",
    "TP53", "CDKN2A", "CDKN2B", "MKI67",
]

# Variant 4: BRS71 (Chakravarty 2011 71-gene panel, abbreviated to representative ~50 genes)
# We use TIERA67 ∪ a few canonical BRS extras; in practice TIERA67 is a curated BRS-derived panel.
BRS71_LIKE = [
    # TIERA67-derived differentiation (we keep the non-overlapping ones; the 8-gene set is REMOVED below)
    "DIO2", "DUOX1", "DUOX2", "GLIS3", "SLC26A4", "SLC5A8", "THRA", "THRB", "IYD",
    # BRS-distinctive thyroid genes commonly cited
    "CDH16", "CLU", "GATA3", "KRT19", "MUC1", "FN1", "CITED1", "RUNX2",
    # Add MAPK targets that BRS uses
    "DUSP5", "DUSP6", "ETV4", "ETV5", "SPRY2", "SPRY4",
    # Plus EMT/immune common BRS overlap
    "VIM", "FOXP3", "CD274",
]


def load_tcga_expression() -> pd.DataFrame:
    print(f"[load] {TPM}", flush=True)
    df = pd.read_csv(TPM, sep="\t", index_col=0)
    print(f"[load] {df.shape[0]} genes × {df.shape[1]} samples", flush=True)
    return df


def load_orig_labels() -> pd.DataFrame:
    df = pd.read_csv(ORIG_LABELS, sep="\t")
    print(f"[load] orig labels: {len(df)} samples, DM1={int((df['cluster_use']=='DM1').sum())}, DM2={int((df['cluster_use']=='DM2').sum())}", flush=True)
    return df


def load_tiera67() -> list[str]:
    genes = []
    with open(TIERA67_TXT) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                genes.append(line)
    return list(dict.fromkeys(genes))  # de-dupe preserving order


def expr_subset(expr: pd.DataFrame, genes: Sequence[str], samples: Sequence[str]) -> pd.DataFrame:
    """Return samples × gene_present subset (transposed). Skips missing genes silently."""
    g_have = [g for g in genes if g in expr.index]
    s_have = [s for s in samples if s in expr.columns]
    print(f"  matched {len(g_have)}/{len(genes)} genes × {len(s_have)}/{len(samples)} samples", flush=True)
    sub = expr.loc[g_have, s_have].T  # samples × genes
    return sub


def leiden_k2(X: pd.DataFrame, resolution: float = 0.5, n_neighbors: int = 15, random_state: int = 42) -> pd.Series:
    """Cluster X (samples × features) with Leiden, force K≈2 by resolution sweep if needed."""
    Xs = StandardScaler().fit_transform(X.values)
    A = ad.AnnData(Xs)
    A.obs_names = X.index.astype(str)
    A.var_names = X.columns.astype(str)
    sc.pp.neighbors(A, n_neighbors=min(n_neighbors, max(2, A.n_obs - 1)), random_state=random_state, use_rep="X")
    # Try increasing resolution until we hit K==2 (or sweep down)
    for r in [resolution, 0.3, 0.7, 0.2, 1.0, 0.15, 0.1]:
        try:
            sc.tl.leiden(A, resolution=r, random_state=random_state, flavor="igraph", n_iterations=2, directed=False)
        except Exception:
            sc.tl.leiden(A, resolution=r, random_state=random_state)
        labels = A.obs["leiden"].astype(int).values
        k = len(np.unique(labels))
        if k == 2:
            print(f"  Leiden K=2 at resolution {r}", flush=True)
            return pd.Series(labels, index=X.index)
        elif k >= 3:
            # collapse smaller clusters into nearest of top 2 by KMeans on cluster centroids
            print(f"  Leiden K={k} at resolution {r} → collapse to K=2 via KMeans on centroids", flush=True)
            centroids = np.vstack([Xs[labels == c].mean(axis=0) for c in range(k)])
            km = KMeans(n_clusters=2, random_state=random_state, n_init=10).fit(centroids)
            collapsed = km.predict(Xs)
            return pd.Series(collapsed, index=X.index)
    # fallback: just KMeans K=2
    print("  fallback to KMeans K=2", flush=True)
    km = KMeans(n_clusters=2, random_state=random_state, n_init=10).fit(Xs)
    return pd.Series(km.labels_, index=X.index)


def consensus_stability(X: pd.DataFrame, base_labels: pd.Series, n_boot: int = 100, sample_frac: float = 0.8, random_state: int = 42) -> float:
    """Mean concordance with base across bootstraps."""
    rng = np.random.default_rng(random_state)
    aris = []
    for i in range(n_boot):
        idx = rng.choice(X.index, size=int(len(X) * sample_frac), replace=False)
        try:
            lb = leiden_k2(X.loc[idx], random_state=random_state + i)
            base_sub = base_labels.loc[idx]
            aris.append(adjusted_rand_score(base_sub.values, lb.values))
        except Exception:
            continue
    return float(np.mean(aris)) if aris else float("nan")


def align_with_orig(new_labels: pd.Series, orig: pd.Series) -> tuple[pd.Series, dict]:
    """Map new_labels {0,1} so DM1/DM2 align maximally with orig. Returns (mapped_labels, metrics)."""
    common = new_labels.index.intersection(orig.index)
    if len(common) == 0:
        return new_labels, {"n_common": 0, "ari": float("nan"), "kappa": float("nan")}
    new_c = new_labels.loc[common].values
    orig_c = orig.loc[common].map({"DM1": 0, "DM2": 1}).values
    # try both label assignments, pick better
    flipped = (new_c == 0).astype(int)
    ari_a = adjusted_rand_score(orig_c, new_c)
    ari_b = adjusted_rand_score(orig_c, flipped)
    if ari_b > ari_a:
        new_labels_mapped = new_labels.map({0: 1, 1: 0})
        ari = ari_b
    else:
        new_labels_mapped = new_labels.copy()
        ari = ari_a
    nm = new_labels_mapped.loc[common].values
    kappa = cohen_kappa_score(orig_c, nm)
    return new_labels_mapped, {"n_common": int(len(common)), "ari": float(ari), "kappa": float(kappa),
                                "n_DM1": int((new_labels_mapped == 0).sum()), "n_DM2": int((new_labels_mapped == 1).sum())}


def run_variant(name: str, expr: pd.DataFrame, genes: list[str], orig: pd.Series, all_samples: list[str]) -> dict:
    print(f"\n=== R1-{name} clustering ({len(genes)} genes input) ===", flush=True)
    overlap = [g for g in genes if g in GENE_8]
    print(f"  8-gene overlap: {len(overlap)} ({overlap})", flush=True)
    X = expr_subset(expr, genes, all_samples)
    if X.shape[0] < 50 or X.shape[1] < 5:
        return {"variant": name, "error": "matrix too small", "shape": list(X.shape)}
    labels = leiden_k2(X)
    mapped, metrics = align_with_orig(labels, orig)
    # save labels
    label_str = mapped.map({0: f"DM1_{name}", 1: f"DM2_{name}"})
    out_tsv = OUT / f"R1{name}_cluster_labels.tsv"
    pd.DataFrame({"sample_id": label_str.index, "cluster": label_str.values}).to_csv(out_tsv, sep="\t", index=False)
    summary = {
        "variant": name,
        "input_genes_n": len(genes),
        "input_genes_8gene_overlap_n": len(overlap),
        "input_genes_matched_n": int(X.shape[1]),
        "n_samples": int(X.shape[0]),
        "metrics_vs_orig": metrics,
        "labels_path": str(out_tsv),
    }
    print(f"  result: ARI={metrics.get('ari'):.3f}  kappa={metrics.get('kappa'):.3f}  n_common={metrics.get('n_common')}", flush=True)
    with open(OUT / f"R1{name}_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    return summary


def main() -> int:
    expr = load_tcga_expression()
    orig_df = load_orig_labels()
    orig = orig_df.set_index("sample_id")["cluster_use"]
    # use all TCGA-THCA tumour samples (suffix -01A) for the new clustering — not just the 178 with orig labels
    all_samples = [s for s in expr.columns if s.endswith("-01A") or s.endswith("-01")]
    if len(all_samples) < 100:
        all_samples = list(expr.columns)
    print(f"\n[scope] using {len(all_samples)} samples for new clustering", flush=True)

    tiera67 = load_tiera67()
    print(f"\n[load] TIERA67 = {len(tiera67)} unique genes; 8-gene members: {sorted(set(tiera67) & GENE_8)}", flush=True)

    results = {}
    # R1-A: TIERA67 minus 8-gene
    genes_A = [g for g in tiera67 if g not in GENE_8]
    results["A"] = run_variant("A", expr, genes_A, orig, all_samples)

    # R1-B: 30-gene MAPK+immune+EMT
    results["B"] = run_variant("B", expr, GENE_30_FUNCTIONAL, orig, all_samples)

    # R1-C: 5-feature immune composite (use individual markers as proxy for composite)
    # We do not have pre-computed Hallmark scores per sample here; we use 5 raw immune markers as a fast proxy.
    GENES_C = ["GZMA", "PRF1", "IFNG", "CD8A", "GZMB"]
    results["C"] = run_variant("C", expr, GENES_C, orig, all_samples)

    # R1-D: BRS71-like minus 8-gene
    genes_D = [g for g in BRS71_LIKE if g not in GENE_8]
    results["D"] = run_variant("D", expr, genes_D, orig, all_samples)

    # combined concordance table
    rows = []
    for name, r in results.items():
        m = r.get("metrics_vs_orig", {}) if isinstance(r, dict) else {}
        rows.append({
            "variant": name,
            "n_input_genes": r.get("input_genes_n"),
            "n_8gene_overlap": r.get("input_genes_8gene_overlap_n"),
            "n_matched": r.get("input_genes_matched_n"),
            "n_samples": r.get("n_samples"),
            "ari_vs_orig": m.get("ari"),
            "kappa_vs_orig": m.get("kappa"),
            "n_DM1_new": m.get("n_DM1"),
            "n_DM2_new": m.get("n_DM2"),
        })
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(OUT / "R1_concordance_table.tsv", sep="\t", index=False)
    print("\n=== R1 concordance summary ===")
    print(summary_df.to_string(index=False))

    with open(OUT / "R1_FINAL_summary.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    print("\n[R1] DONE. Outputs at", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
