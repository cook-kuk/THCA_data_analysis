#!/usr/bin/env python3
"""v15: K2 score overlay + GSE250521 spatial MAPK-panel + DepMap/PRISM overlay.

This is a reviewer-reserve extension of the Paper 1 Fig 8 mechanism layer.
It deliberately avoids voice-protected manuscript prose and writes only
analysis tables/figures.
"""
from __future__ import annotations

import gzip
import json
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

warnings.filterwarnings("ignore")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project" / "results" / "p_deconv_2026_05_08"
OUT.mkdir(parents=True, exist_ok=True)

PANEL_8 = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
MAPK_OUT = ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1", "CCND1"]

TCGA_LOG2 = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv")
TCGA_LABELS = ROOT / "project" / "results" / "v17_realfix" / "R1A_cluster_labels.tsv"
LEE_LOG2 = Path("/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_log2.tsv")
LEE_META = ROOT / "project" / "results" / "v17_korean" / "GSE213647_panel_score.tsv"
GENE_MAP = Path("/data/thca/repo_results/v17p3/tables/F1_gene_recovery_mapping.tsv")
K2_TPM = ROOT / "project" / "results" / "v17_korean" / "K2_8gene_tpm_matrix_v4.tsv"
K2_PREV = ROOT / "project" / "results" / "v17_korean" / "K2_korean_predictions_v4.tsv"

GSE250521_META = ROOT / "project" / "data" / "processed" / "GSE250521" / "sample_metadata.tsv"

PRISM = Path("/data/thca/repo_results/paper11_pancancer/phase_H_prism/dm1_drug_dependency.tsv")
DEPMAP = Path("/data/thca/repo_results/paper11_pancancer/phase_D_depmap/dm1_high_low_dependency.tsv")


def center_rows(x: np.ndarray) -> np.ndarray:
    return x - x.mean(axis=1, keepdims=True)


def read_symbol_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    first = df.columns[0]
    return df.set_index(first)


def collapse_ensembl_to_symbol(df: pd.DataFrame, gene_map: pd.DataFrame) -> pd.DataFrame:
    ens_to_sym = dict(zip(gene_map["ensembl"], gene_map["symbol"]))
    keep = [idx for idx in df.index if idx in ens_to_sym]
    out = df.loc[keep].copy()
    out.index = [ens_to_sym[idx] for idx in out.index]
    return out.groupby(out.index).mean()


def fit_tcga_centered_panel_model() -> tuple[StandardScaler, LogisticRegression, pd.DataFrame]:
    tcga = read_symbol_matrix(TCGA_LOG2)
    labels = pd.read_csv(TCGA_LABELS, sep="\t")
    labels["is_DM2"] = labels["cluster"].astype(str).str.startswith("DM2").astype(int)
    samples = [s for s in labels["sample_id"] if s in tcga.columns]
    genes = [g for g in PANEL_8 if g in tcga.index]
    x = tcga.loc[genes, samples].T[genes].to_numpy(dtype=float)
    y = labels.set_index("sample_id").loc[samples, "is_DM2"].to_numpy(dtype=int)
    x_centered = center_rows(x)
    scaler = StandardScaler().fit(x_centered)
    x_scaled = scaler.transform(x_centered)
    model = LogisticRegression(C=1.0, max_iter=3000, random_state=42)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    oof = cross_val_predict(model, x_scaled, y, cv=cv, method="predict_proba")[:, 1]
    auc_oof = roc_auc_score(y, oof)
    model.fit(x_scaled, y)
    train_p = model.predict_proba(x_scaled)[:, 1]
    tcga_scores = pd.DataFrame(
        {
            "sample_id": samples,
            "source": "TCGA-THCA",
            "group": np.where(y == 1, "TCGA_DM2", "TCGA_DM1"),
            "p_DM2": train_p,
            "p_DM2_oof": oof,
            "DM_call": np.where(oof >= 0.5, "DM2", "DM1"),
            "score_basis": "centered_panel8_log2",
            "n_panel_genes": len(genes),
        }
    )
    tcga_scores.attrs["auc_oof"] = float(auc_oof)
    tcga_scores.attrs["auc_train"] = float(roc_auc_score(y, train_p))
    return scaler, model, tcga_scores


def score_external_panel(
    expr: pd.DataFrame,
    sample_ids: list[str],
    scaler: StandardScaler,
    model: LogisticRegression,
    source: str,
    group_values: pd.Series | None = None,
    score_basis: str = "centered_panel8_log2",
) -> pd.DataFrame:
    genes = [g for g in PANEL_8 if g in expr.index]
    x = expr.loc[genes, sample_ids].T[genes].to_numpy(dtype=float)
    p = model.predict_proba(scaler.transform(center_rows(x)))[:, 1]
    groups = group_values.reindex(sample_ids).fillna(source).astype(str).to_numpy() if group_values is not None else np.repeat(source, len(sample_ids))
    return pd.DataFrame(
        {
            "sample_id": sample_ids,
            "source": source,
            "group": groups,
            "p_DM2": p,
            "p_DM2_oof": np.nan,
            "DM_call": np.where(p >= 0.5, "DM2", "DM1"),
            "score_basis": score_basis,
            "n_panel_genes": len(genes),
        }
    )


def run_k2_overlay() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    print("\n[v15A] K2 centered-panel score overlay")
    scaler, model, tcga_scores = fit_tcga_centered_panel_model()
    auc_oof = tcga_scores.attrs["auc_oof"]
    print(f"  TCGA centered-panel 5-fold OOF AUC: {auc_oof:.4f}")

    gene_map = pd.read_csv(GENE_MAP, sep="\t")
    lee = read_symbol_matrix(LEE_LOG2)
    lee = collapse_ensembl_to_symbol(lee, gene_map)
    lee_meta = pd.read_csv(LEE_META, sep="\t").set_index("gsm")
    lee_samples = [s for s in lee.columns if s in lee_meta.index]
    lee_group = lee_meta["histology"].replace({"UTC/ATC": "ATC"})
    lee_scores = score_external_panel(
        lee,
        lee_samples,
        scaler,
        model,
        "Lee_GSE213647",
        group_values=lee_group.map(lambda x: f"Lee_{x}"),
    )

    k2 = pd.read_csv(K2_TPM, sep="\t").set_index("run")
    k2_expr = np.log2(k2[PANEL_8] + 1.0).T
    k2_scores = score_external_panel(
        k2_expr,
        list(k2_expr.columns),
        scaler,
        model,
        "K2_PRJEB11591",
        group_values=pd.Series("K2_PRJEB11591", index=k2_expr.columns),
        score_basis="centered_panel8_log2_mini_index_tpm",
    )
    if K2_PREV.exists():
        prev = pd.read_csv(K2_PREV, sep="\t").set_index("run")
        k2_scores["p_DM2_prev_v4"] = prev["p_DM2"].reindex(k2_scores["sample_id"]).to_numpy()
        k2_scores["abs_delta_prev_v4"] = (k2_scores["p_DM2"] - k2_scores["p_DM2_prev_v4"]).abs()
    else:
        k2_scores["p_DM2_prev_v4"] = np.nan
        k2_scores["abs_delta_prev_v4"] = np.nan

    all_scores = pd.concat([tcga_scores, lee_scores, k2_scores], ignore_index=True)
    all_scores.to_csv(OUT / "v15_k2_panel_profile_scores.tsv", sep="\t", index=False)

    rows = []
    for (source, group), sub in all_scores.groupby(["source", "group"], sort=False):
        vals = sub["p_DM2"].dropna()
        rows.append(
            {
                "source": source,
                "group": group,
                "n": int(len(vals)),
                "median_p_DM2": float(vals.median()),
                "mean_p_DM2": float(vals.mean()),
                "q25_p_DM2": float(vals.quantile(0.25)),
                "q75_p_DM2": float(vals.quantile(0.75)),
                "frac_DM2_call": float((vals >= 0.5).mean()),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "v15_k2_score_distribution_summary.tsv", sep="\t", index=False)
    k2_delta_max = float(k2_scores["abs_delta_prev_v4"].max(skipna=True)) if "abs_delta_prev_v4" in k2_scores else np.nan
    metrics = {
        "tcga_centered_panel_oof_auc": float(auc_oof),
        "tcga_centered_panel_train_auc": float(tcga_scores.attrs["auc_train"]),
        "k2_n": int(len(k2_scores)),
        "k2_dm2_calls": int((k2_scores["p_DM2"] >= 0.5).sum()),
        "k2_dm1_calls": int((k2_scores["p_DM2"] < 0.5).sum()),
        "k2_median_p_DM2": float(k2_scores["p_DM2"].median()),
        "k2_max_abs_delta_vs_previous_v4": k2_delta_max,
        "lee_n": int(len(lee_scores)),
    }
    print(f"  K2 DM2 calls: {metrics['k2_dm2_calls']}/{metrics['k2_n']}; median p_DM2={metrics['k2_median_p_DM2']:.3f}")
    return all_scores, summary, metrics


def get_expr(adata, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan)
    j = adata.var_names.get_loc(gene)
    col = adata.X[:, j]
    if sp.issparse(col):
        return col.toarray().ravel()
    return np.asarray(col).ravel()


def module_z(adata, genes: list[str]) -> tuple[np.ndarray, int]:
    vals = []
    for g in genes:
        v = get_expr(adata, g)
        if np.all(~np.isfinite(v)):
            continue
        sd = np.nanstd(v)
        if not np.isfinite(sd) or sd == 0:
            vals.append(np.zeros_like(v, dtype=float))
        else:
            vals.append((v - np.nanmean(v)) / (sd + 1e-9))
    if not vals:
        return np.full(adata.n_obs, np.nan), 0
    return np.nanmean(np.vstack(vals), axis=0), len(vals)


def spearman_block(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5:
        return np.nan, np.nan, int(keep.sum())
    rho, p = stats.spearmanr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def residual_spearman_block(x: np.ndarray, y: np.ndarray, covar: np.ndarray) -> tuple[float, float, int]:
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(covar)
    if keep.sum() < 10 or np.nanstd(covar[keep]) == 0:
        return np.nan, np.nan, int(keep.sum())
    xx = covar[keep].reshape(-1, 1)
    x_res = x[keep] - LinearRegression().fit(xx, x[keep]).predict(xx)
    y_res = y[keep] - LinearRegression().fit(xx, y[keep]).predict(xx)
    rho, p = stats.spearmanr(x_res, y_res)
    return float(rho), float(p), int(keep.sum())


def run_spatial_mapk_panel() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    print("\n[v15B] GSE250521 spatial MAPK output x Panel-8")
    import scanpy as sc

    meta = pd.read_csv(GSE250521_META, sep="\t")
    per_sample = []
    spot_chunks = []
    for _, r in meta.iterrows():
        sid = r["sample_id"]
        stage = r["stage_inferred"]
        raw_h5 = ROOT / str(r["h5ad"]).replace("project/", "project/")
        scored_h5 = raw_h5.parent / (raw_h5.stem.replace(".raw", "") + ".scored.h5ad")
        use_h5 = scored_h5 if scored_h5.exists() else raw_h5
        print(f"  loading {sid} ({stage})")
        adata = sc.read_h5ad(use_h5)
        if "log1p" not in (adata.uns or {}):
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

        mapk, n_mapk = module_z(adata, MAPK_OUT)
        panel, n_panel = module_z(adata, PANEL_8)
        epithelial = adata.obs["Epithelial_score"].astype(float).to_numpy() if "Epithelial_score" in adata.obs else np.full(adata.n_obs, np.nan)
        epi_cut = np.nanmedian(epithelial) if np.isfinite(epithelial).any() else np.nan
        epi_mask = epithelial >= epi_cut if np.isfinite(epi_cut) else np.ones(adata.n_obs, dtype=bool)
        all_rho, all_p, n_all = spearman_block(mapk, panel)
        epi_rho, epi_p, n_epi = spearman_block(mapk[epi_mask], panel[epi_mask])
        all_resid_rho, all_resid_p, n_all_resid = residual_spearman_block(mapk, panel, epithelial)
        epi_resid_rho, epi_resid_p, n_epi_resid = residual_spearman_block(
            mapk[epi_mask], panel[epi_mask], epithelial[epi_mask]
        )
        per_sample.append(
            {
                "sample_id": sid,
                "stage": stage,
                "n_spots": int(adata.n_obs),
                "n_MAPK_genes": int(n_mapk),
                "n_panel_genes": int(n_panel),
                "rho_mapk_panel_all_spots": all_rho,
                "p_mapk_panel_all_spots": all_p,
                "n_all_corr": n_all,
                "rho_mapk_panel_epithelial_enriched": epi_rho,
                "p_mapk_panel_epithelial_enriched": epi_p,
                "n_epi_corr": n_epi,
                "rho_mapk_panel_all_spots_resid_epithelial": all_resid_rho,
                "p_mapk_panel_all_spots_resid_epithelial": all_resid_p,
                "n_all_resid_corr": n_all_resid,
                "rho_mapk_panel_epithelial_enriched_resid_epithelial": epi_resid_rho,
                "p_mapk_panel_epithelial_enriched_resid_epithelial": epi_resid_p,
                "n_epi_resid_corr": n_epi_resid,
                "epithelial_score_median": float(epi_cut) if np.isfinite(epi_cut) else np.nan,
                "mean_MAPK_z": float(np.nanmean(mapk)),
                "mean_panel8_z": float(np.nanmean(panel)),
            }
        )
        obs = adata.obs
        spot_chunks.append(
            pd.DataFrame(
                {
                    "sample_id": sid,
                    "stage": stage,
                    "spot_id": obs.index.astype(str).to_numpy(),
                    "array_row": obs["array_row"].to_numpy() if "array_row" in obs else np.nan,
                    "array_col": obs["array_col"].to_numpy() if "array_col" in obs else np.nan,
                    "MAPK_z": mapk,
                    "panel8_z": panel,
                    "DM1_like_z": -panel,
                    "Epithelial_score": epithelial,
                    "epithelial_enriched": epi_mask.astype(int),
                }
            )
        )

    sample_df = pd.DataFrame(per_sample)
    spot_df = pd.concat(spot_chunks, ignore_index=True)
    sample_df.to_csv(OUT / "v15_spatial_mapk_panel_per_sample.tsv", sep="\t", index=False)
    with gzip.open(OUT / "v15_spatial_mapk_panel_spot_scores.tsv.gz", "wt") as f:
        spot_df.to_csv(f, sep="\t", index=False)

    stage_rows = []
    for stage, sub in spot_df.groupby("stage", sort=False):
        for subset, s2 in [("all_spots", sub), ("epithelial_enriched", sub[sub["epithelial_enriched"] == 1])]:
            rho, p, n = spearman_block(s2["MAPK_z"].to_numpy(), s2["panel8_z"].to_numpy())
            stage_rows.append({"stage": stage, "subset": subset, "n_spots": n, "rho_mapk_panel": rho, "p_mapk_panel": p})
            rrho, rp, rn = residual_spearman_block(
                s2["MAPK_z"].to_numpy(),
                s2["panel8_z"].to_numpy(),
                s2["Epithelial_score"].to_numpy(),
            )
            stage_rows.append(
                {
                    "stage": stage,
                    "subset": f"{subset}_resid_epithelial",
                    "n_spots": rn,
                    "rho_mapk_panel": rrho,
                    "p_mapk_panel": rp,
                }
            )
    for subset, s2 in [
        ("tumor_all_spots", spot_df[spot_df["stage"] != "PT"]),
        ("tumor_epithelial_enriched", spot_df[(spot_df["stage"] != "PT") & (spot_df["epithelial_enriched"] == 1)]),
    ]:
        rho, p, n = spearman_block(s2["MAPK_z"].to_numpy(), s2["panel8_z"].to_numpy())
        stage_rows.append({"stage": "TUMOR_ONLY", "subset": subset, "n_spots": n, "rho_mapk_panel": rho, "p_mapk_panel": p})
        rrho, rp, rn = residual_spearman_block(
            s2["MAPK_z"].to_numpy(),
            s2["panel8_z"].to_numpy(),
            s2["Epithelial_score"].to_numpy(),
        )
        stage_rows.append(
            {
                "stage": "TUMOR_ONLY",
                "subset": f"{subset}_resid_epithelial",
                "n_spots": rn,
                "rho_mapk_panel": rrho,
                "p_mapk_panel": rp,
            }
        )
    stage_df = pd.DataFrame(stage_rows)
    stage_df.to_csv(OUT / "v15_spatial_mapk_panel_stage_summary.tsv", sep="\t", index=False)
    tumor_epi = stage_df[(stage_df["stage"] == "TUMOR_ONLY") & (stage_df["subset"] == "tumor_epithelial_enriched")]
    tumor_epi_resid = stage_df[
        (stage_df["stage"] == "TUMOR_ONLY")
        & (stage_df["subset"] == "tumor_epithelial_enriched_resid_epithelial")
    ]
    metrics = {
        "n_samples": int(sample_df.shape[0]),
        "n_spots_total": int(spot_df.shape[0]),
        "tumor_epithelial_enriched_rho": float(tumor_epi["rho_mapk_panel"].iloc[0]) if len(tumor_epi) else np.nan,
        "tumor_epithelial_enriched_p": float(tumor_epi["p_mapk_panel"].iloc[0]) if len(tumor_epi) else np.nan,
        "tumor_epithelial_enriched_resid_epithelial_rho": float(tumor_epi_resid["rho_mapk_panel"].iloc[0]) if len(tumor_epi_resid) else np.nan,
        "tumor_epithelial_enriched_resid_epithelial_p": float(tumor_epi_resid["p_mapk_panel"].iloc[0]) if len(tumor_epi_resid) else np.nan,
        "median_per_sample_rho_all_tumor": float(sample_df.loc[sample_df["stage"] != "PT", "rho_mapk_panel_all_spots"].median()),
        "median_per_sample_rho_epi_tumor": float(sample_df.loc[sample_df["stage"] != "PT", "rho_mapk_panel_epithelial_enriched"].median()),
        "median_per_sample_rho_epi_tumor_resid_epithelial": float(
            sample_df.loc[sample_df["stage"] != "PT", "rho_mapk_panel_epithelial_enriched_resid_epithelial"].median()
        ),
    }
    print(
        "  tumor epithelial-enriched rho="
        f"{metrics['tumor_epithelial_enriched_rho']:.3f}; "
        f"median sample rho={metrics['median_per_sample_rho_epi_tumor']:.3f}; "
        f"resid rho={metrics['tumor_epithelial_enriched_resid_epithelial_rho']:.3f}"
    )
    return sample_df, stage_df, metrics


MAPK_NAME_PAT = re.compile(
    r"(MEK|RAF|ERK|MAPK|BRAF|trametinib|selumetinib|cobimetinib|binimetinib|ulixertinib|lifirafenib|AZD-0364|AZ 628|RAF709|LY3214996|CC-90003|ML786)",
    flags=re.I,
)
KNOWN_MAPK = {"AZD-0364", "RAF709", "AZ 628", "CC-90003", "ML786", "LY3214996", "lifirafenib"}
KNOWN_OTHER = {
    "elacestrant": "ER/SERD",
    "SR-1078": "ROR agonist",
    "lylamine": "amine/lipid",
    "GDC-0077": "PI3K",
}


def classify_drug(name: str) -> tuple[str, bool]:
    s = str(name)
    if s in KNOWN_MAPK or MAPK_NAME_PAT.search(s):
        return "MAPK_pathway", True
    if s in KNOWN_OTHER:
        return KNOWN_OTHER[s], False
    return "other_or_unannotated", False


def run_depmap_prism_overlay() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    print("\n[v15C] DepMap/PRISM MAPK-inhibitor overlay")
    pr = pd.read_csv(PRISM, sep="\t").sort_values("p").reset_index(drop=True)
    classes = pr["name"].map(classify_drug)
    pr["pathway_class"] = [x[0] for x in classes]
    pr["is_mapk_pathway"] = [x[1] for x in classes]
    pr["rank_by_p"] = np.arange(1, len(pr) + 1)
    pr["rank_by_abs_d"] = pr["cohens_d"].abs().rank(method="first", ascending=False).astype(int)
    pr.to_csv(OUT / "v15_prism_mapk_overlay.tsv", sep="\t", index=False)

    enrich_rows = []
    n_total = len(pr)
    n_mapk = int(pr["is_mapk_pathway"].sum())
    for label, mask in {
        "FDR_lt_0.05": pr["fdr"] < 0.05,
        "p_lt_1e-4": pr["p"] < 1e-4,
        "top15_by_p": pr["rank_by_p"] <= 15,
        "top15_by_abs_d": pr["rank_by_abs_d"] <= 15,
    }.items():
        n_hit = int(mask.sum())
        n_mapk_hit = int((mask & pr["is_mapk_pathway"]).sum())
        p_hyper = float(stats.hypergeom.sf(n_mapk_hit - 1, n_total, n_mapk, n_hit)) if n_hit and n_mapk else np.nan
        enrich_rows.append(
            {
                "set": label,
                "n_total_drugs": n_total,
                "n_mapk_universe": n_mapk,
                "n_hit": n_hit,
                "n_mapk_hit": n_mapk_hit,
                "mapk_fraction_hit": n_mapk_hit / n_hit if n_hit else np.nan,
                "fold_enrichment": (n_mapk_hit / n_hit) / (n_mapk / n_total) if n_hit and n_mapk else np.nan,
                "hypergeom_p": p_hyper,
            }
        )
    enrich = pd.DataFrame(enrich_rows)
    enrich.to_csv(OUT / "v15_prism_mapk_enrichment.tsv", sep="\t", index=False)

    dep = pd.read_csv(DEPMAP, sep="\t").drop_duplicates("gene").sort_values("p")
    dep["dependency_interpretation"] = np.where(
        dep["cohens_d"] < 0,
        "more_essential_in_DM1_high",
        "more_essential_in_DM1_low",
    )
    dep.to_csv(OUT / "v15_depmap_dependency_overlay.tsv", sep="\t", index=False)
    fdr = enrich[enrich["set"] == "FDR_lt_0.05"].iloc[0].to_dict()
    metrics = {
        "prism_n_drugs": int(n_total),
        "prism_n_mapk_universe_regex": int(n_mapk),
        "prism_fdr05_n": int(fdr["n_hit"]),
        "prism_fdr05_n_mapk": int(fdr["n_mapk_hit"]),
        "prism_fdr05_hypergeom_p": float(fdr["hypergeom_p"]),
        "prism_top_drug": str(pr.iloc[0]["name"]),
        "prism_top_drug_cohens_d": float(pr.iloc[0]["cohens_d"]),
        "depmap_top_dependency": str(dep.iloc[0]["gene"]),
        "depmap_top_dependency_cohens_d": float(dep.iloc[0]["cohens_d"]),
        "depmap_top_dependency_p": float(dep.iloc[0]["p"]),
    }
    print(
        f"  PRISM FDR<0.05 MAPK hits: {metrics['prism_fdr05_n_mapk']}/"
        f"{metrics['prism_fdr05_n']} (p={metrics['prism_fdr05_hypergeom_p']:.2e})"
    )
    return pr, enrich, dep, metrics


def plot_v15(
    k2_scores: pd.DataFrame,
    k2_summary: pd.DataFrame,
    spatial_sample: pd.DataFrame,
    spatial_stage: pd.DataFrame,
    prism: pd.DataFrame,
    depmap: pd.DataFrame,
    metrics: dict,
) -> None:
    print("\n[v15] plotting composite figure")
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    ax1, ax2, ax3, ax4, ax5, ax6 = axes.ravel()

    # A. Cross-cohort centered-panel probability overlay.
    group_order = ["TCGA_DM1", "TCGA_DM2", "Lee_Normal", "Lee_PTC", "Lee_PDFP", "Lee_ATC", "K2_PRJEB11591"]
    labels = ["TCGA\nDM1", "TCGA\nDM2", "Lee\nNormal", "Lee\nPTC", "Lee\nPDFP", "Lee\nATC", "K2\nPRJEB11591"]
    data = [k2_scores.loc[k2_scores["group"] == g, "p_DM2"].dropna().to_numpy() for g in group_order]
    parts = ax1.violinplot(data, showmeans=False, showextrema=False, showmedians=True)
    palette = ["#d95f02", "#1b9e77", "#66a61e", "#7570b3", "#e7298a", "#a6761d", "#2c7fb8"]
    for body, color in zip(parts["bodies"], palette):
        body.set_facecolor(color)
        body.set_alpha(0.65)
        body.set_edgecolor("#222")
    parts["cmedians"].set_color("#111")
    parts["cmedians"].set_linewidth(1.5)
    rng = np.random.default_rng(42)
    for i, vals in enumerate(data, start=1):
        if len(vals) == 0:
            continue
        plot_vals = vals if len(vals) <= 220 else rng.choice(vals, 220, replace=False)
        ax1.scatter(rng.normal(i, 0.045, len(plot_vals)), plot_vals, s=9, color="#111", alpha=0.35)
    ax1.axhline(0.5, color="#555", lw=0.8, ls=":")
    ax1.set_xticks(np.arange(1, len(labels) + 1))
    ax1.set_xticklabels(labels, fontsize=9)
    ax1.set_ylabel("Centered-panel classifier p_DM2")
    ax1.set_title("A. K2 score overlay: scale-invariant 8-gene profile", loc="left", weight="bold")
    ax1.grid(axis="y", alpha=0.25)

    # B. K2 distribution.
    k2 = k2_scores[k2_scores["source"] == "K2_PRJEB11591"]["p_DM2"].dropna()
    ax2.hist(k2, bins=np.linspace(0, 1, 21), color="#2c7fb8", edgecolor="white", alpha=0.9)
    ax2.axvline(0.5, color="#222", lw=1, ls=":")
    ax2.axvline(k2.median(), color="#f6c85f", lw=2)
    ax2.set_xlabel("p_DM2")
    ax2.set_ylabel("K2 samples")
    ax2.set_title(
        f"B. K2 PRJEB11591: {int((k2 >= 0.5).sum())}/{len(k2)} DM2 calls; median={k2.median():.3f}",
        loc="left",
        weight="bold",
    )
    ax2.grid(axis="y", alpha=0.25)

    # C. Spatial per-sample Spearman.
    stage_order = ["PT", "PTC", "LPTC", "ATC"]
    stage_color = {"PT": "#7f7f7f", "PTC": "#2c7fb8", "LPTC": "#e08214", "ATC": "#b2182b"}
    for i, st in enumerate(stage_order):
        sub = spatial_sample[spatial_sample["stage"] == st]
        y = sub["rho_mapk_panel_epithelial_enriched"].to_numpy()
        ax3.scatter(np.full(len(y), i) + rng.normal(0, 0.04, len(y)), y, s=70, color=stage_color[st], edgecolor="#111", alpha=0.9)
        if len(y):
            ax3.plot([i - 0.22, i + 0.22], [np.nanmedian(y), np.nanmedian(y)], color="#111", lw=2)
    ax3.axhline(0, color="#555", ls=":", lw=0.8)
    ax3.set_xticks(np.arange(len(stage_order)))
    ax3.set_xticklabels(stage_order)
    ax3.set_ylabel("Spearman rho (MAPK output vs Panel-8)\nepithelial-enriched spots")
    ax3.set_title("C. Spatial within-sample anti-correlation", loc="left", weight="bold")
    ax3.set_ylim(-1, 0.55)
    ax3.grid(axis="y", alpha=0.25)

    # D. Stage pooled spatial summary.
    stage_sub = spatial_stage[
        spatial_stage["subset"].isin(["epithelial_enriched", "epithelial_enriched_resid_epithelial"])
    ].copy()
    x = np.arange(len(stage_order))
    width = 0.36
    for offset, subset, color, label in [
        (-width / 2, "epithelial_enriched", "#08519c", "epithelial enriched"),
        (width / 2, "epithelial_enriched_resid_epithelial", "#f6c85f", "residualized on epithelial score"),
    ]:
        vals = [
            float(stage_sub[(stage_sub["stage"] == st) & (stage_sub["subset"] == subset)]["rho_mapk_panel"].iloc[0])
            if len(stage_sub[(stage_sub["stage"] == st) & (stage_sub["subset"] == subset)])
            else np.nan
            for st in stage_order
        ]
        ax4.bar(x + offset, vals, width=width, color=color, edgecolor="#111", label=label)
    ax4.axhline(0, color="#555", ls=":", lw=0.8)
    ax4.set_xticks(x)
    ax4.set_xticklabels(stage_order)
    ax4.set_ylabel("Pooled spot-level Spearman rho")
    ax4.set_title("D. Stage-pooled spatial signal: expected anti-correlation not observed", loc="left", weight="bold")
    ax4.legend(frameon=False, fontsize=9)
    ax4.set_ylim(-1, 0.55)
    ax4.grid(axis="y", alpha=0.25)

    # E. PRISM top drugs.
    top_pr = prism.sort_values("p").head(15).iloc[::-1]
    colors = np.where(top_pr["is_mapk_pathway"], "#d95f02", "#7570b3")
    ax5.barh(np.arange(len(top_pr)), top_pr["cohens_d"], color=colors, edgecolor="#111", alpha=0.88)
    ax5.axvline(0, color="#555", ls=":", lw=0.8)
    ax5.set_yticks(np.arange(len(top_pr)))
    ax5.set_yticklabels([f"{r['name']}  FDR={r['fdr']:.1e}" for _, r in top_pr.iterrows()], fontsize=8)
    ax5.set_xlabel("Cohen's d (PRISM LFC: DM1-high minus DM1-low)\nnegative = more growth suppression in DM1-high lines")
    ax5.set_title("E. PRISM DM1-high selective drugs", loc="left", weight="bold")
    ax5.legend(handles=[Patch(facecolor="#d95f02", label="MAPK pathway"), Patch(facecolor="#7570b3", label="other")], frameon=False, fontsize=9)
    ax5.grid(axis="x", alpha=0.25)

    # F. DepMap CRISPR dependencies.
    top_dep = depmap.reindex(depmap["cohens_d"].abs().sort_values(ascending=False).head(12).index).iloc[::-1]
    dep_colors = np.where(top_dep["cohens_d"] < 0, "#b2182b", "#2166ac")
    ax6.barh(np.arange(len(top_dep)), top_dep["cohens_d"], color=dep_colors, edgecolor="#111", alpha=0.88)
    ax6.axvline(0, color="#555", ls=":", lw=0.8)
    ax6.set_yticks(np.arange(len(top_dep)))
    ax6.set_yticklabels([f"{r['gene']}  p={r['p']:.1e}" for _, r in top_dep.iterrows()], fontsize=8)
    ax6.set_xlabel("Cohen's d (CRISPR effect: DM1-high minus DM1-low)\nnegative = more essential in DM1-high lines")
    ax6.set_title("F. DepMap CRISPR dependency overlay", loc="left", weight="bold")
    ax6.grid(axis="x", alpha=0.25)

    fig.suptitle(
        "Figure SX_v15. Reviewer-reserve extensions: K2 profile transfer, spatial MAPK-panel anti-correlation, and PRISM/DepMap actionability.",
        fontsize=13,
        weight="bold",
        y=1.01,
    )
    fig.tight_layout()
    fig.savefig(OUT / "Fig_SX_v15_K2_spatial_DEPMap.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "Fig_SX_v15_K2_spatial_DEPMap.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {OUT / 'Fig_SX_v15_K2_spatial_DEPMap.png'}")


def main() -> None:
    k2_scores, k2_summary, k2_metrics = run_k2_overlay()
    spatial_sample, spatial_stage, spatial_metrics = run_spatial_mapk_panel()
    prism, enrich, depmap, depmap_metrics = run_depmap_prism_overlay()

    all_metrics = {
        "v15A_k2_overlay": k2_metrics,
        "v15B_spatial_mapk_panel": spatial_metrics,
        "v15C_depmap_prism": depmap_metrics,
        "generated_at": "2026-05-09",
        "output_dir": str(OUT),
    }
    plot_v15(k2_scores, k2_summary, spatial_sample, spatial_stage, prism, depmap, all_metrics)
    with open(OUT / "v15_summary.json", "w") as f:
        json.dump(all_metrics, f, indent=2, default=str)

    print("\n[v15] done")
    print(json.dumps(all_metrics, indent=2, default=str))


if __name__ == "__main__":
    main()
