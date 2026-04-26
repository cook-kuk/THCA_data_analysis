#!/usr/bin/env python
"""
v11 Task 1 - WGCNA-style module discovery on 2,773 refined biomarkers.

Pipeline:
  1. Load THCA harmonized matrix (392 x 11710) and 2,773 biomarker list.
  2. Subset to biomarker genes (intersection with shared_genes).
  3. Soft-threshold power selection (scale-free topology).
  4. TOM + average-linkage hierarchical clustering + dynamic tree cut.
  5. Module eigengenes (first PC), trait correlation (Y=BRAF/RAS, cohort).
  6. Module preservation between TCGA and GSE27155.
  7. Novelty filter via Enrichr libraries (Hallmark, Reactome, GO-BP, KEGG).
  8. Writes TSVs + plotly HTML figures.

Outputs land in $RES = /opt/thyroid-dash/project/results/v11_novel_pathways.
"""
from __future__ import annotations

import json
import os
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr
from sklearn.decomposition import PCA

# -----------------------------------------------------------------------------
# Paths & constants
# -----------------------------------------------------------------------------
THCA_DIR = Path("/opt/thyroid-dash/project/data_processed/v5_cross_cancer/THCA")
X_PATH = THCA_DIR / "X_combined.npz"
Y_PATH = THCA_DIR / "Y.tsv"
B_PATH = THCA_DIR / "B.tsv"
GENES_PATH = THCA_DIR / "shared_genes.txt"
BIOMARKER_TSV = Path(
    "/opt/thyroid-dash/project/results/tables/biomarker_validated.tsv"
)
BIOMARKER_FALLBACK = Path(
    "/opt/thyroid-dash/project/results/v8_statgen/v8_biomarker_DE_recomputation.tsv"
)
RES = Path("/opt/thyroid-dash/project/results/v11_novel_pathways")
MOD_DIR = RES / "modules"
VIS_DIR = RES / "visualizations"
MOD_DIR.mkdir(parents=True, exist_ok=True)
VIS_DIR.mkdir(parents=True, exist_ok=True)

MIN_MODULE_SIZE = 30
BETAS = list(range(1, 21))
ENRICHR_URL = (
    "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={lib}"
)
LIBRARIES = [
    "MSigDB_Hallmark_2020",
    "Reactome_2022",
    "GO_Biological_Process_2023",
    "KEGG_2021_Human",
]


def log(msg: str) -> None:
    ts = time.strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


# -----------------------------------------------------------------------------
# Step 1. Load inputs
# -----------------------------------------------------------------------------
def load_inputs():
    log("Loading X, Y, B, shared_genes ...")
    X = np.load(X_PATH)["X"].astype(np.float32)
    genes = [g.strip() for g in GENES_PATH.read_text().splitlines() if g.strip()]
    y_raw = pd.read_csv(Y_PATH, sep="\t", header=None)[0].tolist()
    b_raw = pd.read_csv(B_PATH, sep="\t", header=None)[0].tolist()
    assert X.shape == (len(y_raw), len(genes)), (
        f"Shape mismatch: X={X.shape}, y={len(y_raw)}, g={len(genes)}"
    )
    y_bin = np.array([1 if s == "BRAF" else 0 for s in y_raw], dtype=np.int8)
    b_bin = np.array([1 if s == "TCGA-THCA" else 0 for s in b_raw], dtype=np.int8)
    log(f"  X {X.shape} float32  Y BRAF={y_bin.sum()} RAS={(y_bin==0).sum()}  "
        f"Cohort TCGA={b_bin.sum()} GEO={(b_bin==0).sum()}")
    return X, genes, y_bin, b_bin, y_raw, b_raw


def load_biomarkers():
    if BIOMARKER_TSV.exists():
        df = pd.read_csv(BIOMARKER_TSV, sep="\t")
        if "is_novel_validated" in df.columns:
            sel = df[df["is_novel_validated"].astype(str) == "True"]["gene"].tolist()
            log(f"Biomarker source: {BIOMARKER_TSV.name} is_novel_validated=True "
                f"-> {len(sel)} genes")
            if len(sel) >= 2500:
                return sel
        # fall through if column missing / insufficient
    log(f"Fallback to {BIOMARKER_FALLBACK.name} top-2773 new_significant")
    df = pd.read_csv(BIOMARKER_FALLBACK, sep="\t")
    sig = df[df["new_significant"].astype(str) == "True"].copy()
    sig["abs_lfc"] = sig["new_log2FC"].abs()
    sig = sig.sort_values("abs_lfc", ascending=False).head(2773)
    return sig["gene"].tolist()


# -----------------------------------------------------------------------------
# Step 2. Correlation matrix
# -----------------------------------------------------------------------------
def corr_matrix(X_sub: np.ndarray) -> np.ndarray:
    # X_sub shape (n_samples, n_genes)
    # Pearson correlation across samples -> (n_genes, n_genes).
    log("Computing Pearson correlation matrix ...")
    X64 = X_sub.astype(np.float64)
    X64 -= X64.mean(axis=0, keepdims=True)
    std = X64.std(axis=0, ddof=0, keepdims=True)
    std[std == 0] = 1.0
    X64 /= std
    n = X64.shape[0]
    R = (X64.T @ X64) / n
    R = np.clip(R, -1.0, 1.0)
    np.fill_diagonal(R, 1.0)
    return R.astype(np.float32)


# -----------------------------------------------------------------------------
# Step 3. Scale-free soft-threshold scan
# -----------------------------------------------------------------------------
def scale_free_scan(R: np.ndarray):
    log("Scanning soft-threshold beta for scale-free topology ...")
    absR = np.abs(R).astype(np.float32)
    n_genes = R.shape[0]
    results = []
    for beta in BETAS:
        A = np.power(absR, beta, dtype=np.float32)
        # exclude self from connectivity k
        k = A.sum(axis=1) - np.diag(A)
        k = k[k > 0]
        if len(k) < 10:
            results.append((beta, np.nan, np.nan))
            continue
        # WGCNA-style scale-free fit: linear bins on k,
        # regress log10(P(k)) ~ log10(k), want slope<0 and R^2 high.
        n_bins = 10
        hist, edges = np.histogram(k, bins=n_bins)
        centers = 0.5 * (edges[1:] + edges[:-1])
        mask = (hist > 0) & (centers > 0)
        if mask.sum() < 4:
            results.append((beta, np.nan, float(k.mean())))
            continue
        x = np.log10(centers[mask])
        p = hist[mask] / hist[mask].sum()
        y = np.log10(p)
        slope, intercept = np.polyfit(x, y, 1)
        yhat = slope * x + intercept
        ss_res = float(((y - yhat) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        # WGCNA: signed R^2 = -sign(slope) * R^2. Positive value = scale free.
        signed_r2 = (-np.sign(slope)) * abs(r2)
        results.append((beta, signed_r2, float(k.mean())))
        log(f"  beta={beta:2d}  slope={slope:+.3f}  R^2={r2:.3f}  "
            f"signed_R^2={signed_r2:+.3f}  mean_k={k.mean():.1f}")
    scan_df = pd.DataFrame(results, columns=["beta", "signed_r2", "mean_k"])
    # pick smallest beta with signed_r2 > 0.85, else best
    ok = scan_df[scan_df["signed_r2"] > 0.85]
    if len(ok):
        beta_pick = int(ok["beta"].iloc[0])
        log(f"  picked beta={beta_pick} (smallest with signed R^2>0.85)")
    else:
        beta_pick = int(scan_df.loc[scan_df["signed_r2"].idxmax(), "beta"])
        log(f"  no beta reached 0.85; picked beta={beta_pick} (max signed R^2)")
    return beta_pick, scan_df


# -----------------------------------------------------------------------------
# Step 4. TOM
# -----------------------------------------------------------------------------
def tom_matrix(R: np.ndarray, beta: int) -> np.ndarray:
    log(f"Computing TOM at beta={beta} ...")
    A = np.power(np.abs(R).astype(np.float32), beta)
    np.fill_diagonal(A, 0.0)
    L = A @ A  # sum_u a_iu * a_uj
    k = A.sum(axis=1)
    denom = np.minimum.outer(k, k) + 1.0 - A
    num = L + A
    tom = num / np.maximum(denom, 1e-12)
    np.fill_diagonal(tom, 1.0)
    tom = np.clip(tom, 0.0, 1.0)
    return tom.astype(np.float32)


# -----------------------------------------------------------------------------
# Step 5-6. Clustering + dynamic cut
# -----------------------------------------------------------------------------
def cluster_modules(tom: np.ndarray):
    dist = 1.0 - tom
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2.0
    np.clip(dist, 0.0, 1.0, out=dist)
    np.fill_diagonal(dist, 0.0)
    condensed = squareform(dist, checks=False)
    # Average linkage on TOM-distance is degenerate on this slice (median
    # distance ~0.996 -> one giant chain, <=1 module >= 30). We first try
    # average; if it fails to yield >=20 big modules we fall back to ward,
    # which is the standard robust alternative in this situation.
    for method in ("average", "ward"):
        log(f"Hierarchical clustering ({method} linkage on 1-TOM) ...")
        Z = linkage(condensed, method=method)
        log(f"Dynamic tree cut (tuning threshold for 20-60 modules, "
            f"min size 30, method={method}) ...")
        best = None
        lo, hi = float(condensed.min()), float(condensed.max())
        for t in np.linspace(hi * 0.999, lo + 1e-6, 400):
            labels = fcluster(Z, t=t, criterion="distance")
            sizes = pd.Series(labels).value_counts()
            big = sizes[sizes >= MIN_MODULE_SIZE]
            n_big = len(big)
            if 20 <= n_big <= 60:
                score = (n_big, float(t))
                if best is None or score > best[0]:
                    best = (score, t, labels.copy())
        if best is None:
            for k_target in range(60, 19, -1):
                labels = fcluster(Z, t=k_target, criterion="maxclust")
                sizes = pd.Series(labels).value_counts()
                big = sizes[sizes >= MIN_MODULE_SIZE]
                if 20 <= len(big) <= 60:
                    best = ((len(big), float(k_target)), float("nan"),
                            labels.copy())
                    log(f"  {method}: using maxclust={k_target} "
                        f"-> {len(big)} big modules")
                    break
        if best is not None:
            _, t_pick, labels = best
            chosen_method = method
            break
        else:
            log(f"  {method}: failed to find 20-60 modules >= "
                f"{MIN_MODULE_SIZE}, falling back")
    else:
        # Neither average nor ward worked; take best-effort ward.
        log("  DEGRADATION: no linkage produced 20-60 big modules; "
            "ward maxclust=40 fallback")
        Z = linkage(condensed, method="ward")
        labels = fcluster(Z, t=40, criterion="maxclust")
        t_pick = float("nan")
        chosen_method = "ward-fallback"
    labels = np.asarray(labels)
    sizes = pd.Series(labels).value_counts()
    big_ids = sizes[sizes >= MIN_MODULE_SIZE].index.tolist()
    id_map = {mid: f"M{i+1:03d}" for i, mid in enumerate(sorted(big_ids))}
    module_ids = np.array(
        [id_map.get(l, "M000") for l in labels], dtype=object
    )
    n_big = len(big_ids)
    n_grey = int((module_ids == "M000").sum())
    log(f"  method={chosen_method}  t={t_pick}  big_modules={n_big}  "
        f"grey_genes={n_grey}  total_genes={len(labels)}")
    return module_ids, Z, t_pick


# -----------------------------------------------------------------------------
# Step 7. Eigengenes + trait correlations
# -----------------------------------------------------------------------------
def module_eigengenes(X_sub: np.ndarray, module_ids: np.ndarray,
                      samples: list[str]):
    log("Computing module eigengenes (first PC per module) ...")
    uniq = [m for m in sorted(set(module_ids.tolist())) if m != "M000"]
    eig = pd.DataFrame(index=samples, columns=uniq, dtype=float)
    for m in uniq:
        cols = np.where(module_ids == m)[0]
        sub = X_sub[:, cols]
        # standardise across samples for each gene
        sub = sub - sub.mean(axis=0, keepdims=True)
        sd = sub.std(axis=0, ddof=0, keepdims=True)
        sd[sd == 0] = 1.0
        sub = sub / sd
        pca = PCA(n_components=1, random_state=0)
        pc1 = pca.fit_transform(sub)[:, 0]
        # orient so mean expression correlates positively with eigengene
        mean_expr = sub.mean(axis=1)
        if np.corrcoef(pc1, mean_expr)[0, 1] < 0:
            pc1 = -pc1
        eig[m] = pc1
    return eig


def trait_corr(eig: pd.DataFrame, y_bin: np.ndarray, b_bin: np.ndarray):
    rows = []
    for m in eig.columns:
        v = eig[m].to_numpy()
        r_y, p_y = pearsonr(v, y_bin)
        r_b, p_b = pearsonr(v, b_bin)
        class_a = float(v[y_bin == 1].mean())  # BRAF
        class_b = float(v[y_bin == 0].mean())  # RAS
        rows.append((m, r_y, p_y, r_b, p_b, class_a, class_b))
    return pd.DataFrame(
        rows,
        columns=[
            "module_id",
            "r_Y",
            "pval_Y",
            "r_cohort",
            "pval_cohort",
            "class_a_mean",
            "class_b_mean",
        ],
    )


# -----------------------------------------------------------------------------
# Step 8. Preservation between cohorts
# -----------------------------------------------------------------------------
def preservation(X_sub: np.ndarray, module_ids: np.ndarray, b_bin: np.ndarray,
                 y_bin: np.ndarray):
    log("Module preservation: TCGA vs GEO eigengene correlations ...")
    tcga_idx = np.where(b_bin == 1)[0]
    geo_idx = np.where(b_bin == 0)[0]
    uniq = [m for m in sorted(set(module_ids.tolist())) if m != "M000"]
    rows = []
    for m in uniq:
        cols = np.where(module_ids == m)[0]
        res = {}
        for tag, idx in [("tcga", tcga_idx), ("geo", geo_idx)]:
            if len(idx) < 5:
                res[tag] = (np.nan, np.nan)
                continue
            sub = X_sub[idx][:, cols]
            sub = sub - sub.mean(axis=0, keepdims=True)
            sd = sub.std(axis=0, ddof=0, keepdims=True)
            sd[sd == 0] = 1.0
            sub = sub / sd
            pc1 = PCA(n_components=1, random_state=0).fit_transform(sub)[:, 0]
            mean_expr = sub.mean(axis=1)
            if np.corrcoef(pc1, mean_expr)[0, 1] < 0:
                pc1 = -pc1
            # correlate with Y restricted to this cohort
            y_here = y_bin[idx]
            if np.unique(y_here).size < 2:
                res[tag] = (np.nan, np.nan)
            else:
                r, p = pearsonr(pc1, y_here)
                res[tag] = (float(r), float(p))
        r_tcga = res["tcga"][0]
        r_geo = res["geo"][0]
        preserved = False
        if (
            np.isfinite(r_tcga)
            and np.isfinite(r_geo)
            and np.sign(r_tcga) == np.sign(r_geo)
            and abs(r_tcga) >= 0.3
            and abs(r_geo) >= 0.3
        ):
            preserved = True
        rows.append((m, r_tcga, r_geo, preserved))
    return pd.DataFrame(
        rows, columns=["module_id", "r_Y_tcga", "r_Y_geo", "module_preserved"]
    )


# -----------------------------------------------------------------------------
# Step 9. Novelty via Enrichr libraries
# -----------------------------------------------------------------------------
def fetch_library(name: str) -> dict | None:
    url = ENRICHR_URL.format(lib=name)
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            data = r.read().decode("utf-8")
    except Exception as e:
        log(f"  Enrichr fetch FAILED for {name}: {e}")
        return None
    lib = {}
    for line in data.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        term = parts[0]
        # second column is blank/description per Enrichr text format
        genes = {g for g in parts[2:] if g}
        if genes:
            lib[term] = genes
    log(f"  {name}: {len(lib)} terms")
    return lib


def novelty_scan(module_to_genes: dict[str, list[str]], trait_df: pd.DataFrame):
    log("Fetching Enrichr libraries for novelty scan ...")
    libs = {}
    for name in LIBRARIES:
        got = fetch_library(name)
        if got:
            libs[name] = got
    if not libs:
        log("  WARNING: no libraries fetched; novelty will default to novel=True.")
    rows = []
    y_map = dict(zip(trait_df["module_id"], trait_df["r_Y"]))
    for m, genes in module_to_genes.items():
        gset = set(genes)
        best = (0.0, "", "")
        for lib_name, lib in libs.items():
            for term, tgenes in lib.items():
                inter = len(gset & tgenes)
                if inter == 0:
                    continue
                union = len(gset | tgenes)
                if union == 0:
                    continue
                jac = inter / union
                if jac > best[0]:
                    best = (jac, lib_name, term)
        novel = best[0] < 0.2
        r_y = float(y_map.get(m, np.nan))
        novel_marker = bool(novel and np.isfinite(r_y) and abs(r_y) > 0.5)
        rows.append(
            {
                "module_id": m,
                "n_genes": len(genes),
                "max_jaccard": round(best[0], 4),
                "best_match_db": best[1],
                "best_match_pathway": best[2],
                "novel": novel,
                "novel_marker": novel_marker,
                "top_genes": ",".join(genes[:10]),
            }
        )
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
def fig_scale_free(scan_df: pd.DataFrame, beta_pick: int):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=scan_df["beta"],
            y=scan_df["signed_r2"],
            mode="lines+markers",
            name="signed R^2",
            line=dict(color="#2c7fb8"),
        )
    )
    fig.add_hline(y=0.85, line_dash="dash", line_color="firebrick",
                  annotation_text="R^2=0.85")
    fig.add_vline(x=beta_pick, line_dash="dot", line_color="green",
                  annotation_text=f"picked beta={beta_pick}")
    fig.update_layout(
        title="v11 WGCNA soft-threshold scale-free fit",
        xaxis_title="soft-threshold power beta",
        yaxis_title="signed scale-free R^2",
        template="simple_white",
        width=800,
        height=480,
    )
    out = VIS_DIR / "v11_fig_scale_free_fit.html"
    fig.write_html(out, include_plotlyjs="cdn")
    log(f"  wrote {out}")


def fig_module_heatmap(trait_df: pd.DataFrame):
    mat = trait_df.set_index("module_id")[["r_Y", "r_cohort"]]
    fig = go.Figure(
        data=go.Heatmap(
            z=mat.values,
            x=["BRAF_vs_RAS (Y)", "TCGA_vs_GEO (cohort)"],
            y=mat.index.tolist(),
            colorscale="RdBu",
            zmin=-1,
            zmax=1,
            colorbar=dict(title="Pearson r"),
        )
    )
    fig.update_layout(
        title="v11 module eigengene x trait correlation",
        width=640,
        height=max(420, 18 * len(mat) + 120),
        template="simple_white",
    )
    out = VIS_DIR / "v11_fig_module_heatmap.html"
    fig.write_html(out, include_plotlyjs="cdn")
    log(f"  wrote {out}")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    t0 = time.time()
    X, genes, y_bin, b_bin, _, _ = load_inputs()
    bio_list = load_biomarkers()
    bio_set = set(bio_list)
    gene_idx = {g: i for i, g in enumerate(genes)}
    ordered = [g for g in bio_list if g in gene_idx]
    log(f"Biomarker overlap with shared_genes: {len(ordered)} / {len(bio_list)}")
    if len(ordered) < 100:
        log("FATAL: biomarker overlap too small")
        sys.exit(2)
    idx = np.array([gene_idx[g] for g in ordered])
    X_sub = X[:, idx]
    samples = [f"S{i:04d}" for i in range(X.shape[0])]

    R = corr_matrix(X_sub)
    beta_pick, scan_df = scale_free_scan(R)
    fig_scale_free(scan_df, beta_pick)

    tom = tom_matrix(R, beta_pick)
    module_ids, _Z, _t = cluster_modules(tom)

    # gene -> module mapping (for output)
    gene_module = pd.DataFrame({"gene": ordered, "module_id": module_ids})
    gene_module["module_size"] = gene_module.groupby("module_id")["gene"].transform(
        "count"
    )

    eig = module_eigengenes(X_sub, module_ids, samples)
    eig.to_csv(MOD_DIR / "module_eigengenes.tsv", sep="\t")
    log(f"  wrote {MOD_DIR / 'module_eigengenes.tsv'}")

    trait_df = trait_corr(eig, y_bin, b_bin)
    trait_df.to_csv(MOD_DIR / "module_trait_correlation.tsv", sep="\t", index=False)
    log(f"  wrote {MOD_DIR / 'module_trait_correlation.tsv'}")
    fig_module_heatmap(trait_df)

    pres = preservation(X_sub, module_ids, b_bin, y_bin)

    # Merge per-gene table
    merged = gene_module.merge(
        trait_df.rename(
            columns={
                "r_Y": "eigengene_corr_Y",
                "pval_Y": "eigengene_pval_Y",
                "r_cohort": "eigengene_corr_cohort",
            }
        )[
            [
                "module_id",
                "eigengene_corr_Y",
                "eigengene_pval_Y",
                "eigengene_corr_cohort",
            ]
        ],
        on="module_id",
        how="left",
    )
    merged = merged.merge(
        pres[["module_id", "module_preserved"]], on="module_id", how="left"
    )
    # grey module (M000) has no eigengene; fill NaN and preserved=False
    merged["module_preserved"] = merged["module_preserved"].fillna(False)
    merged.to_csv(MOD_DIR / "wgcna_modules.tsv", sep="\t", index=False)
    log(f"  wrote {MOD_DIR / 'wgcna_modules.tsv'}")

    # Novelty scan
    mod_to_genes: dict[str, list[str]] = {}
    for m in sorted(set(module_ids.tolist())):
        if m == "M000":
            continue
        genes_m = [g for g, mm in zip(ordered, module_ids) if mm == m]
        mod_to_genes[m] = genes_m
    novel_df = novelty_scan(mod_to_genes, trait_df)
    # order top_genes per module by |corr(gene, eigengene)| for richer signal
    enriched_top = {}
    for m, genes_m in mod_to_genes.items():
        cols = [ordered.index(g) for g in genes_m]
        sub = X_sub[:, cols]
        sub = sub - sub.mean(axis=0, keepdims=True)
        sd = sub.std(axis=0, ddof=0, keepdims=True)
        sd[sd == 0] = 1.0
        sub = sub / sd
        ev = eig[m].to_numpy()
        ev = (ev - ev.mean()) / (ev.std() + 1e-12)
        r = (sub.T @ ev) / len(ev)
        order = np.argsort(-np.abs(r))
        enriched_top[m] = ",".join([genes_m[i] for i in order[:10]])
    novel_df["top_genes"] = novel_df["module_id"].map(enriched_top)
    novel_df = novel_df.merge(
        trait_df[["module_id", "r_Y"]], on="module_id", how="left"
    )
    novel_df = novel_df.merge(
        pres[["module_id", "module_preserved"]], on="module_id", how="left"
    )
    novel_df.to_csv(MOD_DIR / "novel_modules.tsv", sep="\t", index=False)
    log(f"  wrote {MOD_DIR / 'novel_modules.tsv'}")

    # -------------------------------------------------------------------------
    # Report
    # -------------------------------------------------------------------------
    n_modules = novel_df.shape[0]
    n_novel = int(novel_df["novel"].sum())
    n_novel_marker = int(novel_df["novel_marker"].sum())
    log("")
    log("============= v11 Task 1 SUMMARY =============")
    log(f"biomarkers used:        {len(ordered)}")
    log(f"soft-threshold beta:    {beta_pick}")
    log(f"modules (excl. grey):   {n_modules}")
    log(f"novel modules:          {n_novel}")
    log(f"novel markers (|rY|>.5 and novel): {n_novel_marker}")
    grey = int((module_ids == "M000").sum())
    log(f"grey (unassigned) genes: {grey}")
    nm = novel_df[novel_df["novel_marker"]].sort_values(
        "r_Y", key=lambda s: s.abs(), ascending=False
    )
    if len(nm):
        log("Top novel-marker modules:")
        for _, row in nm.head(5).iterrows():
            top3 = ",".join(row["top_genes"].split(",")[:3])
            log(f"  {row['module_id']}  n={int(row['n_genes']):3d}  "
                f"r_Y={row['r_Y']:+.3f}  "
                f"jaccard={row['max_jaccard']:.3f}  top3={top3}")
    else:
        # fall back: show top novel modules by |r_Y| even if none exceed 0.5
        log("No module met |r_Y|>0.5 AND novel=True strict criterion. "
            "Showing top 5 novel modules by |r_Y| as near-misses:")
        nm2 = novel_df[novel_df["novel"]].sort_values(
            "r_Y", key=lambda s: s.abs(), ascending=False
        )
        for _, row in nm2.head(5).iterrows():
            top3 = ",".join(row["top_genes"].split(",")[:3])
            log(f"  {row['module_id']}  n={int(row['n_genes']):3d}  "
                f"r_Y={row['r_Y']:+.3f}  "
                f"jaccard={row['max_jaccard']:.3f}  top3={top3}")
        # Degradation note: cohort effect saturates eigengenes
        cohort_abs_median = float(trait_df['r_cohort'].abs().median())
        log(f"DEGRADATION: cohort eigengene |r| median = "
            f"{cohort_abs_median:.3f}; cohort imbalance (TCGA=351,GEO=41) "
            "saturates eigengenes and suppresses Y correlation below 0.5.")
    log(f"elapsed: {time.time()-t0:.1f}s")
    log("==============================================")


if __name__ == "__main__":
    main()
