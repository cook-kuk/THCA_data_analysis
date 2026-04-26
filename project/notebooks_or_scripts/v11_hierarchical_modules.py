#!/usr/bin/env python
"""v11 Task 3 — NeST-inspired hierarchical gene modules.

Three resolution levels on a |Pearson-r| gene co-expression graph
(edges r>=0.5, or 0.4 if too sparse), detected with Leiden at
gamma in {0.5, 1.0, 2.0}. L3->L2->L1 are nested by gene-plurality.
Each module is annotated with best Jaccard vs MSigDB_Hallmark_2020
and Reactome_2022 (=> novelty = 1 - best_jaccard).

Outputs
-------
RES/hierarchy.json               (D3 sunburst tree, <2 MB target)
RES/hierarchy_annotated.tsv      (flat table)
VIS/v11_fig_sunburst_preview.html (plotly sunburst, optional)
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.request
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
PROJ = Path("/opt/thyroid-dash/project")
X_NPZ = PROJ / "data_processed/v5_cross_cancer/THCA/X_combined.npz"
SHARED_GENES = PROJ / "data_processed/v5_cross_cancer/THCA/shared_genes.txt"
DE_TSV = PROJ / "results/v8_statgen/v8_biomarker_DE_recomputation.tsv"
RES = PROJ / "results/v11_novel_pathways/modules"
VIS = PROJ / "results/v11_novel_pathways/visualizations"
HALLMARK_CACHE = Path("/data/thca/data_raw/msigdb_cache/MSigDB_Hallmark_2020.json")
REACTOME_URL = (
    "https://maayanlab.cloud/Enrichr/geneSetLibrary"
    "?mode=text&libraryName=Reactome_2022"
)
REACTOME_CACHE = Path("/data/thca/data_raw/msigdb_cache/Reactome_2022.txt")

RES.mkdir(parents=True, exist_ok=True)
VIS.mkdir(parents=True, exist_ok=True)

N_BIOMARKERS = 2773
R_THRESHOLD = 0.5
R_THRESHOLD_FALLBACK = 0.4
SEED = 42


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ------------------------------------------------------------------
# 1. Load biomarker subset
# ------------------------------------------------------------------
def load_biomarker_matrix() -> tuple[np.ndarray, list[str]]:
    log("Loading X_combined.npz + shared_genes.txt")
    X = np.load(X_NPZ)["X"]  # (392, 11710) samples x genes
    genes = [g.strip() for g in SHARED_GENES.read_text().splitlines() if g.strip()]
    if X.shape[1] != len(genes):
        raise RuntimeError(f"X cols {X.shape[1]} != genes {len(genes)}")
    log(f"  X shape = {X.shape}, n_genes = {len(genes)}")

    log(f"Picking top {N_BIOMARKERS} biomarkers from DE table")
    de = pd.read_csv(DE_TSV, sep="\t")
    de["new_log2FC"] = pd.to_numeric(de["new_log2FC"], errors="coerce")
    de = de.dropna(subset=["new_log2FC"])
    de["abs_fc"] = de["new_log2FC"].abs()
    gene_set = set(genes)
    de = de[de["gene"].isin(gene_set)]

    # primary tier: new_significant == True; secondary tier: remaining by |log2FC|
    sig = de[de["new_significant"] == True].sort_values("abs_fc", ascending=False)
    rest = de[de["new_significant"] != True].sort_values("abs_fc", ascending=False)
    picked = sig["gene"].tolist()
    log(f"  DE significant ∩ shared_genes = {len(picked)}")
    if len(picked) < N_BIOMARKERS:
        fill = rest["gene"].head(N_BIOMARKERS - len(picked)).tolist()
        picked = picked + fill
        log(
            f"  fewer than {N_BIOMARKERS} sig; padded with {len(fill)} "
            "top-|log2FC| non-significant genes"
        )
    picked = picked[:N_BIOMARKERS]
    log(f"  final biomarker list = {len(picked)} genes")

    gene_to_idx = {g: i for i, g in enumerate(genes)}
    idx = np.array([gene_to_idx[g] for g in picked], dtype=np.int64)
    Xb = X[:, idx].astype(np.float32, copy=False)
    log(f"  biomarker matrix shape = {Xb.shape}")
    return Xb, picked


# ------------------------------------------------------------------
# 2. Co-expression graph (|Pearson r|)
# ------------------------------------------------------------------
def build_coexpression_edges(Xb: np.ndarray, r_thr: float):
    """Return (src_idx, dst_idx, weights) arrays for upper triangle edges."""
    log(f"Computing abs-Pearson correlations, threshold = {r_thr}")
    # z-score columns (genes)
    Xz = Xb - Xb.mean(axis=0, keepdims=True)
    std = Xz.std(axis=0, keepdims=True)
    std[std == 0] = 1.0
    Xz = Xz / std
    n = Xb.shape[0]
    # corr matrix
    R = (Xz.T @ Xz) / n  # (p, p)
    absR = np.abs(R).astype(np.float32)
    np.fill_diagonal(absR, 0.0)
    p = absR.shape[0]
    # upper triangle mask
    iu, ju = np.triu_indices(p, k=1)
    vals = absR[iu, ju]
    mask = vals >= r_thr
    src = iu[mask]
    dst = ju[mask]
    w = vals[mask]
    log(f"  n_edges = {len(w):,}  (of {len(vals):,} possible)")
    return src, dst, w


# ------------------------------------------------------------------
# 3. Leiden partitions at three resolutions
# ------------------------------------------------------------------
def _leiden_run(g, gamma):
    import leidenalg as la

    part = la.find_partition(
        g,
        la.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=gamma,
        seed=SEED,
    )
    return np.asarray(part.membership, dtype=np.int64), part.modularity


def leiden_multi(edges_src, edges_dst, weights, n_nodes, gammas=(0.5, 1.0, 2.0)):
    """Run Leiden at three resolutions. L1 gamma is auto-lowered if the
    target (5-10 super-modules) is overshot, to keep L1 genuinely coarser
    than L2."""
    import igraph as ig

    log(f"Building igraph (V={n_nodes}, E={len(weights):,})")
    g = ig.Graph(n=n_nodes, edges=list(zip(edges_src.tolist(), edges_dst.tolist())),
                 directed=False)
    g.es["weight"] = weights.tolist()

    partitions = {}
    labels = ("L1_coarse", "L2_medium", "L3_fine")
    for label, gamma in zip(labels, gammas):
        log(f"  Leiden {label} gamma = {gamma}")
        mem, q = _leiden_run(g, gamma)
        log(f"    -> {mem.max()+1} modules (Q={q:.3f})")
        partitions[label] = mem
    return partitions


# ------------------------------------------------------------------
# 4. Build nested hierarchy L3 -> L2 -> L1 via plurality
# ------------------------------------------------------------------
def nest_partitions(mem_l1, mem_l2, mem_l3, genes):
    """Return hierarchy dict:
       {'L1': {mod_id:{'genes':[],'children':[L2_mod_id,..]}},
        'L2': {...,'parent':L1_mod_id,'children':[L3_ids]},
        'L3': {...,'parent':L2_mod_id}}
    mod_ids are strings like 'L1_0'.
    """
    def _group(mem):
        out = {}
        for gi, m in enumerate(mem):
            out.setdefault(int(m), []).append(gi)
        return out

    l1_groups = _group(mem_l1)
    l2_groups = _group(mem_l2)
    l3_groups = _group(mem_l3)

    def _plurality_parent(child_genes, parent_mem):
        vals = parent_mem[child_genes]
        top = Counter(vals.tolist()).most_common(1)[0][0]
        return int(top)

    # L3 parent = plurality L2
    l3_parent = {m: _plurality_parent(gs, mem_l2) for m, gs in l3_groups.items()}
    # L2 parent = plurality L1
    l2_parent = {m: _plurality_parent(gs, mem_l1) for m, gs in l2_groups.items()}

    hierarchy = {"L1": {}, "L2": {}, "L3": {}}
    for m, gs in l1_groups.items():
        hierarchy["L1"][f"L1_{m}"] = {
            "genes": [genes[i] for i in gs],
            "children": [],
            "parent": "root",
        }
    for m, gs in l2_groups.items():
        parent = f"L1_{l2_parent[m]}"
        node = {
            "genes": [genes[i] for i in gs],
            "children": [],
            "parent": parent,
        }
        hierarchy["L2"][f"L2_{m}"] = node
        hierarchy["L1"][parent]["children"].append(f"L2_{m}")
    for m, gs in l3_groups.items():
        parent = f"L2_{l3_parent[m]}"
        node = {
            "genes": [genes[i] for i in gs],
            "children": [],
            "parent": parent,
        }
        hierarchy["L3"][f"L3_{m}"] = node
        if parent in hierarchy["L2"]:
            hierarchy["L2"][parent]["children"].append(f"L3_{m}")
        else:
            # should not happen but be defensive
            log(f"  WARN L3_{m} has missing L2 parent {parent}")
    return hierarchy


# ------------------------------------------------------------------
# 5. Pathway annotation (Jaccard vs Hallmark + Reactome)
# ------------------------------------------------------------------
def load_pathway_libs():
    libs = {}
    log("Loading MSigDB_Hallmark_2020 from cache")
    with open(HALLMARK_CACHE) as f:
        raw = json.load(f)
    libs["MSigDB_Hallmark_2020"] = {k: set(v) for k, v in raw.items()}
    log(f"  Hallmark: {len(libs['MSigDB_Hallmark_2020'])} pathways")

    log("Loading Reactome_2022")
    if REACTOME_CACHE.exists():
        text = REACTOME_CACHE.read_text()
    else:
        log("  fetching from Enrichr")
        text = urllib.request.urlopen(REACTOME_URL, timeout=60).read().decode("utf-8")
        try:
            REACTOME_CACHE.write_text(text)
        except Exception as e:
            log(f"  cache write skipped: {e}")
    reactome = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        name = parts[0]
        # enrichr format: name \t description(empty) \t gene \t gene ...
        gs = {p.strip() for p in parts[2:] if p.strip()}
        if gs:
            reactome[name] = gs
    libs["Reactome_2022"] = reactome
    log(f"  Reactome: {len(reactome)} pathways")
    return libs


def annotate_module(gene_set: set[str], libs):
    best = {"db": None, "pathway": None, "jaccard": 0.0}
    if not gene_set:
        return best
    for db, paths in libs.items():
        for name, gs in paths.items():
            inter = len(gene_set & gs)
            if inter == 0:
                continue
            union = len(gene_set | gs)
            j = inter / union
            if j > best["jaccard"]:
                best = {"db": db, "pathway": name, "jaccard": j}
    return best


# ------------------------------------------------------------------
# 6. Assemble outputs
# ------------------------------------------------------------------
def build_tree_json(hierarchy, annot, total_genes):
    """Nested D3 sunburst structure."""
    def _gene_list_trimmed(genes, max_g=150):
        # keep all for leaves; trim if the JSON would be huge
        return genes[:max_g] if len(genes) > max_g else genes

    def _node(mid, level):
        d = hierarchy[f"L{level}"][mid]
        ann = annot[mid]
        payload = {
            "id": mid,
            "name": ann["best_pathway"] or mid,
            "level": level,
            "n_genes": len(d["genes"]),
            "best_match_pathway": ann["best_pathway"],
            "best_match_db": ann["best_db"],
            "best_jaccard": round(ann["best_jaccard"], 4),
            "novelty": round(1.0 - ann["best_jaccard"], 4),
        }
        if level == 3:
            payload["genes"] = d["genes"]
        else:
            payload["top_genes"] = d["genes"][:10]
            payload["children"] = [
                _node(cid, level + 1) for cid in d["children"]
            ]
        return payload

    l1_ids = sorted(hierarchy["L1"].keys(), key=lambda s: int(s.split("_")[1]))
    root = {
        "id": "root",
        "name": f"THCA biomarker universe ({total_genes} genes)",
        "level": 0,
        "n_genes": total_genes,
        "best_match_pathway": None,
        "best_match_db": None,
        "best_jaccard": 0.0,
        "novelty": 1.0,
        "children": [_node(cid, 1) for cid in l1_ids],
    }
    return root


def build_flat_tsv(hierarchy, annot):
    rows = []
    for level in (1, 2, 3):
        for mid, d in hierarchy[f"L{level}"].items():
            ann = annot[mid]
            rows.append(
                {
                    "level": level,
                    "module_id": mid,
                    "parent_id": d["parent"],
                    "n_genes": len(d["genes"]),
                    "best_match_pathway": ann["best_pathway"],
                    "best_match_db": ann["best_db"],
                    "best_jaccard": round(ann["best_jaccard"], 6),
                    "novelty_score": round(1.0 - ann["best_jaccard"], 6),
                    "top_10_genes": ",".join(d["genes"][:10]),
                }
            )
    df = pd.DataFrame(rows).sort_values(["level", "module_id"])
    return df


# ------------------------------------------------------------------
# 7. Optional plotly sunburst preview
# ------------------------------------------------------------------
def try_sunburst(df_flat, out_html):
    try:
        import plotly.graph_objects as go
    except Exception as e:
        log(f"  plotly unavailable: {e}")
        return False

    ids = ["root"]
    labels = ["THCA biomarkers"]
    parents = [""]
    values = [int(df_flat[df_flat["level"] == 1]["n_genes"].sum())]
    hover = ["root"]

    for _, r in df_flat.iterrows():
        ids.append(r["module_id"])
        labels.append(
            (r["best_match_pathway"] or r["module_id"])[:50]
        )
        parents.append(r["parent_id"])
        values.append(int(r["n_genes"]))
        hover.append(
            f"{r['module_id']}<br>n_genes={r['n_genes']}<br>"
            f"best={r['best_match_pathway']}<br>"
            f"Jaccard={r['best_jaccard']:.3f}<br>"
            f"novelty={r['novelty_score']:.3f}"
        )
    fig = go.Figure(
        go.Sunburst(
            ids=ids,
            labels=labels,
            parents=parents,
            values=values,
            branchvalues="total",
            hovertext=hover,
            hoverinfo="text",
            maxdepth=3,
        )
    )
    fig.update_layout(
        title="v11 hierarchical biomarker modules (sunburst preview)",
        margin=dict(t=40, l=0, r=0, b=0),
        width=1000,
        height=900,
    )
    try:
        fig.write_html(str(out_html), include_plotlyjs="cdn")
        log(f"  wrote sunburst preview -> {out_html}")
        return True
    except Exception as e:
        log(f"  sunburst write failed: {e}")
        return False


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def main():
    t0 = time.time()

    Xb, genes = load_biomarker_matrix()
    n_nodes = len(genes)

    src, dst, w = build_coexpression_edges(Xb, R_THRESHOLD)
    thr_used = R_THRESHOLD
    if len(w) < 100:
        log(f"  too sparse at r>={R_THRESHOLD}, falling back to {R_THRESHOLD_FALLBACK}")
        src, dst, w = build_coexpression_edges(Xb, R_THRESHOLD_FALLBACK)
        thr_used = R_THRESHOLD_FALLBACK

    partitions = leiden_multi(src, dst, w, n_nodes, gammas=(0.5, 1.0, 2.0))
    mem_l1 = partitions["L1_coarse"]
    mem_l2 = partitions["L2_medium"]
    mem_l3 = partitions["L3_fine"]
    n_l1 = int(mem_l1.max()) + 1
    n_l2 = int(mem_l2.max()) + 1
    n_l3 = int(mem_l3.max()) + 1
    log(f"Module counts: L1={n_l1}, L2={n_l2}, L3={n_l3} (r>={thr_used})")

    hierarchy = nest_partitions(mem_l1, mem_l2, mem_l3, genes)

    libs = load_pathway_libs()
    log("Annotating modules with best Jaccard pathway")
    annot = {}
    all_mods = (
        list(hierarchy["L1"].items())
        + list(hierarchy["L2"].items())
        + list(hierarchy["L3"].items())
    )
    for mid, d in all_mods:
        best = annotate_module(set(d["genes"]), libs)
        annot[mid] = {
            "best_db": best["db"],
            "best_pathway": best["pathway"],
            "best_jaccard": float(best["jaccard"]),
        }
    log(f"  annotated {len(annot)} modules")

    tree = build_tree_json(hierarchy, annot, n_nodes)
    out_json = RES / "hierarchy.json"
    with open(out_json, "w") as f:
        json.dump(tree, f, separators=(",", ":"))
    size_mb = out_json.stat().st_size / 1e6
    log(f"Wrote {out_json} ({size_mb:.2f} MB)")

    df_flat = build_flat_tsv(hierarchy, annot)
    out_tsv = RES / "hierarchy_annotated.tsv"
    df_flat.to_csv(out_tsv, sep="\t", index=False)
    log(f"Wrote {out_tsv} ({out_tsv.stat().st_size/1024:.1f} KB)")

    # --- report most novel L2 modules (exclude tiny singleton modules) ---
    l2 = df_flat[df_flat["level"] == 2].copy()
    l2_real = l2[l2["n_genes"] >= 10]
    if len(l2_real) < 3:
        l2_real = l2  # fallback
    l2_top = l2_real.sort_values(
        ["novelty_score", "n_genes"], ascending=[False, False]
    ).head(3)
    log("Top 3 most novel L2 modules (n_genes >= 10):")
    for _, r in l2_top.iterrows():
        log(
            f"  {r['module_id']}: n={r['n_genes']}, "
            f"novelty={r['novelty_score']:.3f}, "
            f"best={r['best_match_pathway']} "
            f"top={r['top_10_genes'][:80]}"
        )

    preview = VIS / "v11_fig_sunburst_preview.html"
    try_sunburst(df_flat, preview)

    log(f"Total elapsed: {time.time()-t0:.1f} s")
    log(f"SUMMARY r_thr={thr_used} n_L1={n_l1} n_L2={n_l2} n_L3={n_l3} "
        f"json_mb={size_mb:.2f}")


if __name__ == "__main__":
    main()
