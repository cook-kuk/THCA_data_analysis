"""
Track 33 — Spatial transcriptomic TLS niche × HLA-II niche overlap.

Boundary: HLA gene-expression module ONLY, not allele genotype. Bridge zone
(HT+PTC) is hypothesis-generation only — no allele claim, no causal claim.
See project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md (Section 1.1).

Inputs:
  GSE250521 Visium spatial atlas — 16 samples (N/PTC/LPTC/ATC × 4 each):
    full-gene scored h5ad : project/data/processed/GSE250521/<gsm>/<gsm>.scored.h5ad
    MOSCOT mapped celltype : project/results/pantheonos_demo/<sample>/visium_with_mapped.h5ad

Goal:
  Niche-level (NOT spot-level) TLS × HLA-II co-localization. Track 12 reported
  spot-level rho ~ -0.18 due to multicell-spot composition artifacts. This
  Track 33 defines TLS niches as connected components of TLS-marker-high spots
  on the hex grid, then aggregates HLA-II module score per niche vs background.

Outputs:
  project/results/hla_deepdive_2026_05_08/track33_spatial_tls_hla/
"""

from __future__ import annotations

import json
import os
import warnings
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats, sparse
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

# ---------- paths ----------
ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track33_spatial_tls_hla"
FIG = OUT / "figs"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

VISIUM_DIR = ROOT / "project/data/processed/GSE250521"
PANTHEON_DIR = ROOT / "project/results/pantheonos_demo"

# ---------- gene lists ----------
HLA_II = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
          "HLA-DQB1", "HLA-DMA", "HLA-DMB", "CIITA", "CD74"]

# core TLS scoring markers (high specificity / present in atlas)
TLS_CORE = ["CXCL13", "MS4A1", "CD79A"]
# extended TLS markers for zone analysis
TLS_FULL = ["CXCL13", "CXCR5", "CCR6", "CCR7", "MS4A1", "CD79A", "CD79B",
            "CD3D", "IL21", "BCL6", "AICDA", "IGHV1-69"]

# zone markers
B_ZONE = ["MS4A1", "CD79A", "CD79B"]
T_ZONE = ["CD3D"]
FDC_ZONE = ["CXCL13"]

DM1_GENES = ["TG", "TPO", "TSHR", "SLC5A5", "DUOX1", "DUOX2", "DIO1", "DIO2"]


# ---------- sample registry ----------
def list_samples():
    """Return list of (gsm_dir_name, sample_short, stage)."""
    out = []
    for d in sorted(VISIUM_DIR.iterdir()):
        if not d.is_dir() or not d.name.startswith("GSM"):
            continue
        gsm = d.name
        # GSM7980860_N-1 -> short "N-1", stage "N"
        short = gsm.split("_", 1)[1]
        stage = short.split("-")[0]
        out.append((gsm, short, stage))
    return out


# ---------- helpers ----------
def _genes_in(a, glist):
    s = set(a.var_names.tolist())
    return [g for g in glist if g in s]


def _zscore(v):
    v = np.asarray(v, dtype=float)
    s = v.std()
    if s < 1e-9:
        return v - v.mean()
    return (v - v.mean()) / s


def hex_neighbors(array_row, array_col):
    """Return list of neighbor indices for each spot using Visium hex offsets.

    Visium uses an oddr layout: rows alternate offsets. Neighbor offsets in
    (row, col) for a hex grid where col steps by 2 within row, by 1 between
    rows: (-2, 0), (+2, 0), (-1, -1), (-1, +1), (+1, -1), (+1, +1).
    """
    coords = np.stack([array_row, array_col], axis=1).astype(int)
    n = coords.shape[0]
    coord_to_idx = {tuple(c): i for i, c in enumerate(coords)}
    offsets = [(-2, 0), (2, 0), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    neigh = [[] for _ in range(n)]
    for i, (r, c) in enumerate(coords):
        for dr, dc in offsets:
            j = coord_to_idx.get((r + dr, c + dc))
            if j is not None:
                neigh[i].append(j)
    return neigh


def smooth_score(score, neigh, alpha=0.5):
    """k-NN spatial smoothing: new = (1-alpha)*self + alpha*mean(neighbors)."""
    score = np.asarray(score, dtype=float)
    out = score.copy()
    for i, ns in enumerate(neigh):
        if not ns:
            continue
        out[i] = (1 - alpha) * score[i] + alpha * np.mean(score[ns])
    return out


def connected_components(node_set, neigh):
    """Connected components on a subset of nodes given adjacency `neigh`."""
    comp = {}
    cid = 0
    for s in node_set:
        if s in comp:
            continue
        # BFS
        cid += 1
        stack = [s]
        while stack:
            u = stack.pop()
            if u in comp:
                continue
            comp[u] = cid
            for v in neigh[u]:
                if v in node_set and v not in comp:
                    stack.append(v)
    return comp


# ---------- per-sample analysis ----------
def analyze_sample(gsm: str, short: str, stage: str) -> dict:
    p_full = VISIUM_DIR / gsm / f"{gsm}.scored.h5ad"
    p_pan = PANTHEON_DIR / short / "visium_with_mapped.h5ad"
    if not p_full.exists():
        return {"sample": short, "stage": stage, "status": "missing_full"}
    a = sc.read_h5ad(p_full)

    # bring mapped_celltype across (index match)
    if p_pan.exists():
        b = sc.read_h5ad(p_pan)
        b_obs = b.obs[["mapped_celltype", "mapped_celltype_conf"]].copy()
        a.obs = a.obs.join(b_obs, how="left")
    else:
        a.obs["mapped_celltype"] = np.nan
        a.obs["mapped_celltype_conf"] = np.nan

    # --- score modules (HLA-II, TLS_CORE, B/T/FDC zones) ---
    h2 = _genes_in(a, HLA_II)
    tls_core = _genes_in(a, TLS_CORE)
    bz = _genes_in(a, B_ZONE)
    tz = _genes_in(a, T_ZONE)
    fdc = _genes_in(a, FDC_ZONE)
    dm1g = _genes_in(a, DM1_GENES)
    sc.tl.score_genes(a, h2, score_name="HLA_II_score", use_raw=False)
    sc.tl.score_genes(a, tls_core, score_name="TLS_score_raw", use_raw=False)
    if bz:
        sc.tl.score_genes(a, bz, score_name="Bzone_score", use_raw=False)
    else:
        a.obs["Bzone_score"] = np.nan
    if tz:
        sc.tl.score_genes(a, tz, score_name="Tzone_score", use_raw=False)
    else:
        a.obs["Tzone_score"] = np.nan
    if fdc:
        sc.tl.score_genes(a, fdc, score_name="FDC_score", use_raw=False)
    else:
        a.obs["FDC_score"] = np.nan
    if dm1g:
        sc.tl.score_genes(a, dm1g, score_name="DM1_score_check", use_raw=False)

    # --- spatial smoothing on hex grid ---
    rows = a.obs["array_row"].astype(int).values
    cols = a.obs["array_col"].astype(int).values
    neigh = hex_neighbors(rows, cols)

    a.obs["HLA_II_smooth"] = smooth_score(a.obs["HLA_II_score"].values, neigh, alpha=0.5)
    a.obs["TLS_score_smooth"] = smooth_score(a.obs["TLS_score_raw"].values, neigh, alpha=0.5)

    # --- TLS niche definition: top 10% TLS smoothed ≥ z=0 ---
    tls_v = a.obs["TLS_score_smooth"].values
    thr = np.quantile(tls_v, 0.90)
    tls_hi = set(int(i) for i in np.where(tls_v >= thr)[0])
    # require additional positivity: spot-level core marker > 0 in raw mean (sanity)
    comp = connected_components(tls_hi, neigh)

    # filter niches: drop singleton components (size==1) and require at least 3 spots
    sizes = pd.Series(list(comp.values())).value_counts()
    valid_cids = set(sizes[sizes >= 3].index.tolist())
    spot_niche = np.full(a.n_obs, 0, dtype=int)  # 0 = not in niche
    for spot_i, cid in comp.items():
        if cid in valid_cids:
            spot_niche[spot_i] = cid
    a.obs["TLS_niche_id"] = spot_niche
    a.obs["in_TLS_niche"] = (spot_niche > 0).astype(int)

    # --- niche-level stats ---
    niche_rows = []
    pxl_r = a.obs["pxl_row_in_fullres"].astype(float).values
    pxl_c = a.obs["pxl_col_in_fullres"].astype(float).values
    for cid in sorted(valid_cids):
        idx = np.where(spot_niche == cid)[0]
        if len(idx) == 0:
            continue
        bg_idx = np.where(spot_niche == 0)[0]
        d = {
            "sample": short,
            "stage": stage,
            "niche_id": int(cid),
            "n_spots": int(len(idx)),
            "mean_TLS_smooth": float(tls_v[idx].mean()),
            "mean_HLA_II_smooth": float(a.obs["HLA_II_smooth"].values[idx].mean()),
            "bg_mean_HLA_II_smooth": float(a.obs["HLA_II_smooth"].values[bg_idx].mean()) if len(bg_idx) > 0 else np.nan,
            "mean_Bzone": float(a.obs["Bzone_score"].values[idx].mean()) if "Bzone_score" in a.obs else np.nan,
            "mean_Tzone": float(a.obs["Tzone_score"].values[idx].mean()) if "Tzone_score" in a.obs else np.nan,
            "mean_FDC": float(a.obs["FDC_score"].values[idx].mean()) if "FDC_score" in a.obs else np.nan,
            "centroid_row": float(np.median(pxl_r[idx])),
            "centroid_col": float(np.median(pxl_c[idx])),
        }
        niche_rows.append(d)
    niche_df = pd.DataFrame(niche_rows)

    # --- niche vs background HLA-II two-sample test (per sample) ---
    in_niche = np.where(spot_niche > 0)[0]
    out_niche = np.where(spot_niche == 0)[0]
    if len(in_niche) >= 5 and len(out_niche) >= 5:
        h_in = a.obs["HLA_II_smooth"].values[in_niche]
        h_out = a.obs["HLA_II_smooth"].values[out_niche]
        t, p = stats.ttest_ind(h_in, h_out, equal_var=False)
        # cohen d
        pooled = np.sqrt((h_in.var(ddof=1) + h_out.var(ddof=1)) / 2.0) + 1e-9
        d_eff = (h_in.mean() - h_out.mean()) / pooled
    else:
        t, p, d_eff = np.nan, np.nan, np.nan

    # --- distance gradient: distance from each spot to nearest niche centroid ---
    dist_rows = []
    if len(niche_df) > 0:
        ncent = niche_df[["centroid_row", "centroid_col"]].values
        tree = cKDTree(ncent)
        spot_xy = np.stack([pxl_r, pxl_c], axis=1)
        dist, _ = tree.query(spot_xy, k=1)
        # bin distances (microns ≈ pixel at full-res; for relative use within sample)
        # bin into deciles by distance among non-niche spots
        non_niche_dist = dist[spot_niche == 0]
        if len(non_niche_dist) > 10:
            edges = np.quantile(non_niche_dist, np.linspace(0, 1, 11))
            edges[0] = -1e-9
            dbin = np.digitize(dist, edges) - 1
            dbin = np.clip(dbin, 0, 9)
            for k in range(10):
                m = (dbin == k) & (spot_niche == 0)
                if m.sum() > 0:
                    dist_rows.append({
                        "sample": short,
                        "stage": stage,
                        "decile": int(k),
                        "n_spots": int(m.sum()),
                        "mean_dist": float(dist[m].mean()),
                        "mean_HLA_II_smooth": float(a.obs["HLA_II_smooth"].values[m].mean()),
                    })
        # IN niche row (decile = -1)
        m = spot_niche > 0
        if m.sum() > 0:
            dist_rows.append({
                "sample": short, "stage": stage, "decile": -1,
                "n_spots": int(m.sum()),
                "mean_dist": 0.0,
                "mean_HLA_II_smooth": float(a.obs["HLA_II_smooth"].values[m].mean()),
            })
    dist_df = pd.DataFrame(dist_rows)

    # --- HLA-II producers within TLS niche (per cell-type) ---
    prod_rows = []
    if "mapped_celltype" in a.obs.columns and a.obs["mapped_celltype"].notna().any():
        in_idx = np.where(spot_niche > 0)[0]
        if len(in_idx) >= 5:
            ct = a.obs["mapped_celltype"].values
            hla = a.obs["HLA_II_smooth"].values
            # take both inside-niche and overall
            for label in pd.unique(ct[in_idx]):
                if pd.isna(label):
                    continue
                m_in = (np.array([str(x) for x in ct]) == str(label)) & (spot_niche > 0)
                m_out = (np.array([str(x) for x in ct]) == str(label)) & (spot_niche == 0)
                prod_rows.append({
                    "sample": short, "stage": stage, "celltype": str(label),
                    "n_in_niche": int(m_in.sum()),
                    "mean_HLA_II_in_niche": float(hla[m_in].mean()) if m_in.sum() > 0 else np.nan,
                    "n_out_niche": int(m_out.sum()),
                    "mean_HLA_II_out_niche": float(hla[m_out].mean()) if m_out.sum() > 0 else np.nan,
                })
    prod_df = pd.DataFrame(prod_rows)

    # --- DM1 × niche presence ---
    dm1_score = a.obs["DM1_like_score"].values if "DM1_like_score" in a.obs else a.obs.get("DM1_score_check", pd.Series(np.nan, index=a.obs.index)).values
    n_niche = len(valid_cids)
    n_niche_spots = int((spot_niche > 0).sum())
    sample_summary = {
        "sample": short,
        "stage": stage,
        "n_spots": int(a.n_obs),
        "n_TLS_niches": int(n_niche),
        "n_TLS_niche_spots": int(n_niche_spots),
        "frac_niche_spots": float(n_niche_spots / max(a.n_obs, 1)),
        "median_niche_size": float(np.median([s for s in sizes.tolist() if s >= 3])) if any(s >= 3 for s in sizes.tolist()) else np.nan,
        "mean_HLA_II_in_niche": float(a.obs["HLA_II_smooth"].values[in_niche].mean()) if len(in_niche) > 0 else np.nan,
        "mean_HLA_II_out_niche": float(a.obs["HLA_II_smooth"].values[out_niche].mean()) if len(out_niche) > 0 else np.nan,
        "diff_HLA_II_in_minus_out": (float(a.obs["HLA_II_smooth"].values[in_niche].mean()) - float(a.obs["HLA_II_smooth"].values[out_niche].mean())) if (len(in_niche) > 0 and len(out_niche) > 0) else np.nan,
        "ttest_t": float(t) if not np.isnan(t) else np.nan,
        "ttest_p": float(p) if not np.isnan(p) else np.nan,
        "cohen_d": float(d_eff) if not np.isnan(d_eff) else np.nan,
        "spot_level_rho_TLS_HLA_II": float(stats.spearmanr(tls_v, a.obs["HLA_II_smooth"].values).statistic),
        "mean_DM1_score": float(np.nanmean(dm1_score)) if dm1_score is not None else np.nan,
    }

    # --- spatial figure per sample ---
    if a.n_obs > 0:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        for ax, (col, title, cmap) in zip(axes, [
            ("TLS_score_smooth", f"{short} — TLS score (smooth)", "Reds"),
            ("HLA_II_smooth", f"{short} — HLA-II score (smooth)", "Blues"),
            ("TLS_niche_id", f"{short} — TLS niches (n={n_niche})", "tab20"),
        ]):
            v = a.obs[col].values
            sc_h = ax.scatter(pxl_c, -pxl_r, c=v, s=8, cmap=cmap)
            ax.set_aspect("equal")
            ax.set_title(title, fontsize=10)
            ax.set_xticks([])
            ax.set_yticks([])
            plt.colorbar(sc_h, ax=ax, fraction=0.046, pad=0.04)
        plt.tight_layout()
        plt.savefig(FIG / f"sample_{short}_TLS_HLA_II.png", dpi=130, bbox_inches="tight")
        plt.close(fig)

    return {
        "summary": sample_summary,
        "niche_df": niche_df,
        "dist_df": dist_df,
        "prod_df": prod_df,
    }


# ---------- main ----------
def main():
    samples = list_samples()
    all_summaries = []
    all_niche = []
    all_dist = []
    all_prod = []
    for gsm, short, stage in samples:
        print(f"[*] {short} ({stage})")
        try:
            res = analyze_sample(gsm, short, stage)
        except Exception as e:
            print(f"    !! {e}")
            continue
        if isinstance(res, dict) and "summary" in res:
            all_summaries.append(res["summary"])
            if len(res["niche_df"]) > 0:
                all_niche.append(res["niche_df"])
            if len(res["dist_df"]) > 0:
                all_dist.append(res["dist_df"])
            if len(res["prod_df"]) > 0:
                all_prod.append(res["prod_df"])

    summary_df = pd.DataFrame(all_summaries)
    niche_df = pd.concat(all_niche, ignore_index=True) if all_niche else pd.DataFrame()
    dist_df = pd.concat(all_dist, ignore_index=True) if all_dist else pd.DataFrame()
    prod_df = pd.concat(all_prod, ignore_index=True) if all_prod else pd.DataFrame()

    summary_df.to_csv(TAB / "T1_per_sample_summary.tsv", sep="\t", index=False)
    niche_df.to_csv(TAB / "T2_per_niche_stats.tsv", sep="\t", index=False)
    dist_df.to_csv(TAB / "T3_distance_gradient.tsv", sep="\t", index=False)
    prod_df.to_csv(TAB / "T4_HLA_II_producers_per_celltype.tsv", sep="\t", index=False)

    # ----- aggregate plots / tables -----
    # 1. cross-sample niche-level paired diff (in-niche vs background)
    fig, ax = plt.subplots(figsize=(9, 5))
    if len(summary_df):
        order = ["N", "PTC", "LPTC", "ATC"]
        df_p = summary_df.copy()
        df_p["stage"] = pd.Categorical(df_p["stage"], categories=order, ordered=True)
        df_p = df_p.sort_values(["stage", "sample"])
        x = np.arange(len(df_p))
        ax.bar(x - 0.2, df_p["mean_HLA_II_in_niche"], 0.4, label="In-TLS-niche", color="#c0392b")
        ax.bar(x + 0.2, df_p["mean_HLA_II_out_niche"], 0.4, label="Background", color="#7f8c8d")
        ax.set_xticks(x)
        ax.set_xticklabels(df_p["sample"], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("HLA-II smoothed score")
        ax.set_title("F1 — Per-sample HLA-II: TLS niche vs background")
        ax.legend()
        ax.axhline(0, color="k", lw=0.5)
        plt.tight_layout()
        plt.savefig(FIG / "F1_per_sample_HLA_II_in_vs_out_niche.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # 2. effect-size forest per sample
    if "cohen_d" in summary_df:
        fig, ax = plt.subplots(figsize=(7, 7))
        df_f = summary_df.dropna(subset=["cohen_d"]).copy()
        df_f["stage"] = pd.Categorical(df_f["stage"], categories=["N", "PTC", "LPTC", "ATC"], ordered=True)
        df_f = df_f.sort_values(["stage", "sample"])
        y = np.arange(len(df_f))
        cmap = {"N": "#2980b9", "PTC": "#27ae60", "LPTC": "#f39c12", "ATC": "#c0392b"}
        ax.scatter(df_f["cohen_d"], y, c=[cmap.get(s, "k") for s in df_f["stage"]], s=60)
        for yi, (_, row) in zip(y, df_f.iterrows()):
            ax.text(row["cohen_d"] + 0.04, yi, f"p={row['ttest_p']:.1e}", fontsize=7, va="center")
        ax.set_yticks(y)
        ax.set_yticklabels(df_f["sample"])
        ax.axvline(0, color="k", lw=0.5)
        ax.set_xlabel("Cohen d (HLA-II in-niche − out-niche)")
        ax.set_title("F2 — Forest: TLS niche × HLA-II Cohen d, per sample")
        plt.tight_layout()
        plt.savefig(FIG / "F2_forest_cohen_d_per_sample.png", dpi=140, bbox_inches="tight")
        plt.close(fig)

    # 3. stage-aggregated boxplot
    fig, ax = plt.subplots(figsize=(7, 5))
    if len(summary_df):
        order = ["N", "PTC", "LPTC", "ATC"]
        sns.boxplot(data=summary_df, x="stage", y="diff_HLA_II_in_minus_out", order=order,
                    palette={"N": "#2980b9", "PTC": "#27ae60", "LPTC": "#f39c12", "ATC": "#c0392b"}, ax=ax)
        sns.stripplot(data=summary_df, x="stage", y="diff_HLA_II_in_minus_out",
                      order=order, color="k", size=4, alpha=0.7, ax=ax)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylabel("HLA-II diff (in-niche − out)")
        ax.set_title("F3 — Stage stratification: HLA-II at TLS niches")
        plt.tight_layout()
        plt.savefig(FIG / "F3_stage_stratification_diff.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # 4. distance gradient pooled
    if len(dist_df):
        fig, ax = plt.subplots(figsize=(8, 5))
        order = ["N", "PTC", "LPTC", "ATC"]
        agg = dist_df.groupby(["stage", "decile"]).agg(
            mean_HLA_II=("mean_HLA_II_smooth", "mean"),
            sem_HLA_II=("mean_HLA_II_smooth", lambda x: x.std() / max(np.sqrt(len(x)), 1)),
            n=("mean_HLA_II_smooth", "count"),
        ).reset_index()
        cmap = {"N": "#2980b9", "PTC": "#27ae60", "LPTC": "#f39c12", "ATC": "#c0392b"}
        for stg in order:
            sub = agg[agg["stage"] == stg].sort_values("decile")
            if len(sub) == 0:
                continue
            ax.errorbar(sub["decile"], sub["mean_HLA_II"], yerr=sub["sem_HLA_II"],
                        marker="o", label=stg, color=cmap.get(stg, "k"), capsize=3)
        ax.axvline(-0.5, color="k", lw=0.5, ls="--")
        ax.set_xticks([-1, 0, 2, 4, 6, 8, 9])
        ax.set_xticklabels(["IN niche", "D0\n(near)", "D2", "D4", "D6", "D8", "D9\n(far)"], fontsize=8)
        ax.set_ylabel("HLA-II smoothed score")
        ax.set_xlabel("Distance bin from nearest TLS niche centroid")
        ax.set_title("F4 — Distance gradient HLA-II vs TLS niche")
        ax.legend()
        plt.tight_layout()
        plt.savefig(FIG / "F4_distance_gradient.png", dpi=140, bbox_inches="tight")
        plt.close(fig)

    # 5. niche size vs HLA-II intensity
    if len(niche_df):
        fig, ax = plt.subplots(figsize=(7, 5))
        cmap = {"N": "#2980b9", "PTC": "#27ae60", "LPTC": "#f39c12", "ATC": "#c0392b"}
        for stg, sub in niche_df.groupby("stage"):
            ax.scatter(sub["n_spots"], sub["mean_HLA_II_smooth"], label=stg,
                       color=cmap.get(stg, "k"), s=40, alpha=0.7)
        # global linear fit
        x = niche_df["n_spots"].values
        y = niche_df["mean_HLA_II_smooth"].values
        if len(x) > 3:
            r = stats.spearmanr(x, y)
            ax.set_title(f"F5 — Niche size × HLA-II intensity (Spearman ρ={r.statistic:.2f}, p={r.pvalue:.2g}, n={len(x)})")
        ax.set_xlabel("Niche size (n spots)")
        ax.set_ylabel("Mean HLA-II in niche")
        ax.legend()
        plt.tight_layout()
        plt.savefig(FIG / "F5_niche_size_vs_HLA_II.png", dpi=140, bbox_inches="tight")
        plt.close(fig)

    # 6. zone analysis (B / T / FDC means within niche per stage)
    if len(niche_df):
        fig, axes = plt.subplots(1, 3, figsize=(13, 4))
        order = ["N", "PTC", "LPTC", "ATC"]
        for ax, col, title in zip(axes,
                                   ["mean_Bzone", "mean_Tzone", "mean_FDC"],
                                   ["B-zone (MS4A1/CD79A/B)", "T-zone (CD3D)", "FDC (CXCL13)"]):
            sns.boxplot(data=niche_df, x="stage", y=col, order=order,
                        palette={"N": "#2980b9", "PTC": "#27ae60", "LPTC": "#f39c12", "ATC": "#c0392b"},
                        ax=ax)
            sns.stripplot(data=niche_df, x="stage", y=col, order=order,
                          color="k", size=3, alpha=0.6, ax=ax)
            ax.axhline(0, color="k", lw=0.5)
            ax.set_title(title, fontsize=10)
        plt.suptitle("F6 — Per-niche zone composition by stage")
        plt.tight_layout()
        plt.savefig(FIG / "F6_zone_composition.png", dpi=140, bbox_inches="tight")
        plt.close(fig)

    # 7. HLA-II producers ranking (across niches, per cell type)
    if len(prod_df):
        agg_prod = prod_df.groupby("celltype").agg(
            mean_HLA_II_in=("mean_HLA_II_in_niche", "mean"),
            mean_HLA_II_out=("mean_HLA_II_out_niche", "mean"),
            n_in_total=("n_in_niche", "sum"),
            n_samples=("sample", "nunique"),
        ).reset_index()
        agg_prod["delta_in_minus_out"] = agg_prod["mean_HLA_II_in"] - agg_prod["mean_HLA_II_out"]
        agg_prod = agg_prod.sort_values("mean_HLA_II_in", ascending=False)
        agg_prod.to_csv(TAB / "T5_HLA_II_producers_aggregate.tsv", sep="\t", index=False)

        fig, ax = plt.subplots(figsize=(7, 4))
        x = np.arange(len(agg_prod))
        ax.bar(x - 0.2, agg_prod["mean_HLA_II_in"], 0.4, color="#c0392b", label="In niche")
        ax.bar(x + 0.2, agg_prod["mean_HLA_II_out"], 0.4, color="#7f8c8d", label="Background")
        ax.set_xticks(x)
        ax.set_xticklabels(agg_prod["celltype"], rotation=20, ha="right")
        ax.set_ylabel("Mean HLA-II smoothed score")
        ax.set_title("F7 — HLA-II producers in TLS niche (per cell type)")
        ax.axhline(0, color="k", lw=0.5)
        ax.legend()
        plt.tight_layout()
        plt.savefig(FIG / "F7_HLA_II_producers_by_celltype.png", dpi=140, bbox_inches="tight")
        plt.close(fig)
    else:
        agg_prod = pd.DataFrame()

    # 8. DM1 score × niche count per sample (extra)
    if len(summary_df):
        fig, ax = plt.subplots(figsize=(7, 5))
        order = ["N", "PTC", "LPTC", "ATC"]
        cmap = {"N": "#2980b9", "PTC": "#27ae60", "LPTC": "#f39c12", "ATC": "#c0392b"}
        # robust DM1 plotting: clip to 1-99 percentile to avoid extreme outliers
        dm1_vals = summary_df["mean_DM1_score"].astype(float).values
        finite = dm1_vals[np.isfinite(dm1_vals)]
        if len(finite) > 2:
            lo, hi = np.percentile(finite, [1, 99])
            dm1_clip = np.clip(dm1_vals, lo, hi)
        else:
            dm1_clip = dm1_vals
        df_p = summary_df.copy()
        df_p["dm1_plot"] = dm1_clip
        for stg in order:
            sub = df_p[df_p["stage"] == stg]
            if len(sub) == 0:
                continue
            ax.scatter(sub["dm1_plot"], sub["n_TLS_niches"], color=cmap[stg], s=70, label=stg)
            for _, r in sub.iterrows():
                ax.text(r["dm1_plot"], r["n_TLS_niches"] + 0.2, r["sample"], fontsize=7)
        if df_p["dm1_plot"].notna().sum() >= 4:
            r = stats.spearmanr(df_p["dm1_plot"], df_p["n_TLS_niches"], nan_policy="omit")
            ax.set_title(f"F8 — DM1 score × number of TLS niches (Spearman ρ={r.statistic:.2f}, p={r.pvalue:.2g})")
        ax.set_xlabel("Mean DM1 score (sample, clipped 1-99%)")
        ax.set_ylabel("# TLS niches")
        ax.legend()
        plt.tight_layout()
        plt.savefig(FIG / "F8_DM1_x_TLS_niche_count.png", dpi=140)
        plt.close(fig)

    # ----- aggregate stage-level Cohen d -----
    stage_d_rows = []
    for stg in ["N", "PTC", "LPTC", "ATC"]:
        sub = summary_df[summary_df["stage"] == stg]
        ds = sub["cohen_d"].dropna().values
        if len(ds) > 0:
            stage_d_rows.append({
                "stage": stg,
                "n_samples": int(len(ds)),
                "mean_cohen_d": float(np.mean(ds)),
                "sem_cohen_d": float(np.std(ds, ddof=1) / max(np.sqrt(len(ds)), 1)) if len(ds) > 1 else np.nan,
                "mean_diff_in_minus_out": float(sub["diff_HLA_II_in_minus_out"].mean()),
                "median_n_TLS_niches": float(sub["n_TLS_niches"].median()),
            })
    stage_d_df = pd.DataFrame(stage_d_rows)
    stage_d_df.to_csv(TAB / "T6_stage_aggregate_effect.tsv", sep="\t", index=False)

    # ----- summary JSON -----
    summary_json = {
        "track": "Track 33 — Spatial TLS niche × HLA-II co-localization",
        "boundary": "HLA gene-expression module ONLY, not allele genotype. Bridge zone hypothesis-only. Per HLA_CANCER_SEPARATION_RULES Section 1.1.",
        "n_samples_processed": int(len(summary_df)),
        "n_total_niches": int(len(niche_df)),
        "spot_level_TLS_HLA_II_rho_pooled": float(np.nanmean(summary_df["spot_level_rho_TLS_HLA_II"].values)) if "spot_level_rho_TLS_HLA_II" in summary_df else np.nan,
        "niche_level_mean_diff_in_minus_out_global": float(summary_df["diff_HLA_II_in_minus_out"].mean()) if "diff_HLA_II_in_minus_out" in summary_df else np.nan,
        "niche_level_mean_cohen_d_global": float(summary_df["cohen_d"].mean()) if "cohen_d" in summary_df else np.nan,
        "stage_aggregate": stage_d_rows,
        "outputs": {
            "tables": [str(p) for p in sorted(TAB.glob("T*.tsv"))],
            "figures": [str(p) for p in sorted(FIG.glob("F*.png")) + sorted(FIG.glob("sample_*.png"))],
        },
    }
    if len(agg_prod):
        summary_json["HLA_II_producers_top"] = agg_prod.head(10).to_dict(orient="records")

    (OUT / "track33_summary.json").write_text(json.dumps(summary_json, indent=2, default=str))

    print("[+] Done.")
    print(f"    n_samples = {len(summary_df)}")
    print(f"    n_niches  = {len(niche_df)}")
    if len(summary_df):
        print(f"    global cohen_d = {summary_df['cohen_d'].mean():.3f}")
        print(f"    spot-level rho TLS×HLA-II = {summary_df['spot_level_rho_TLS_HLA_II'].mean():.3f}")
    print(f"    out: {OUT}")


if __name__ == "__main__":
    main()
