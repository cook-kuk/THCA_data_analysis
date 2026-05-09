"""v15: final mechanism-extension bundle for Paper 1 Fig 8.

Authorized extensions:
  A. K2 mini-index Panel-only overlay using the scale-invariant TCGA classifier.
  B. GSE250521 spatial MAPK output x Panel-8 per tumor slide.
  C. DepMap/PRISM MAPK-inhibitor sensitivity overlay for DM1-high cell lines.

Outputs live beside v13/v14 so the Fig 8 mechanism dossier has one audit trail.
"""
from __future__ import annotations

import json
import math
import re
import warnings
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/p_deconv_2026_05_08"
OUT.mkdir(parents=True, exist_ok=True)

PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
PANEL_8_CANON = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]
TDS_16 = [
    "DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3", "NKX2-1", "PAX8",
    "SLC26A4", "SLC5A5", "SLC5A8", "TG", "THRA", "THRB", "TPO", "TSHR",
]
MAPK_9 = ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "PHLDA1", "CCND1"]

TCGA_LOG2 = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv")
TCGA_LABELS = ROOT / "project/results/v17_realfix/R1A_cluster_labels.tsv"
LEE_LOG2 = Path("/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_log2.tsv")
LEE_META_PANEL = ROOT / "project/results/v17_korean/GSE213647_panel_score.tsv"
GENE_MAP = Path("/data/thca/repo_results/v17p3/tables/F1_gene_recovery_mapping.tsv")
K2_TPM = Path("/data/thca/repo_results/v17_korean/K2_8gene_tpm_matrix_v4.tsv")
K2_PRED = ROOT / "project/results/v17_korean/K2_korean_predictions_v4.tsv"
SPATIAL_DIR = ROOT / "project/data/processed/GSE250521"
SPATIAL_SCORES = ROOT / "project/results/01_spatial_score/all_spots_scored.tsv.gz"
PRISM_TABLE = ROOT / "project/results/paper11_pancancer/phase_H_prism/dm1_drug_dependency.tsv"
PRISM_RAW = Path("/data/thca/repo_results/p3_p9_full_execution/paper9/raw/Repurposing_Public_24Q2_LFC_COLLAPSED.csv")
H9_RANKING = ROOT / "project/results/p2_braf_nature_sprint_2026_05_09/h9_mapk_therapy/h9_prism_drug_ranking.tsv"
CELLS_DM1 = ROOT / "project/results/paper11_pancancer/phase_D_depmap/celllines_dm1_crispr.tsv"
DEPMAP_DEP = ROOT / "project/results/paper11_pancancer/phase_D_depmap/dm1_high_low_dependency.tsv"


def center_rows(x: np.ndarray) -> np.ndarray:
    return x - np.nanmean(x, axis=1, keepdims=True)


def bh_fdr(pvals: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvals, errors="coerce").to_numpy(float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return pd.Series(out, index=pvals.index)
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    m = len(ranked)
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q_full = np.empty(m, dtype=float)
    q_full[order] = np.clip(q, 0, 1)
    out[np.where(ok)[0]] = q_full
    return pd.Series(out, index=pvals.index)


def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    sp = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / sp) if sp > 0 else float("nan")


def fisher_ci(r: float, n: int) -> tuple[float, float]:
    if n < 4 or not np.isfinite(r) or abs(r) >= 1:
        return (float("nan"), float("nan"))
    z = np.arctanh(r)
    se = 1.0 / math.sqrt(n - 3)
    return tuple(np.tanh([z - 1.96 * se, z + 1.96 * se]))


def fisher_pool(rows: pd.DataFrame, col: str = "spearman_rho") -> dict[str, float]:
    sub = rows.dropna(subset=[col, "n"]).copy()
    sub = sub[(sub["n"] > 3) & (sub[col].abs() < 1)]
    if sub.empty:
        return {"rho": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "n": 0, "k": 0}
    z = np.arctanh(sub[col].to_numpy(float))
    n = sub["n"].to_numpy(float)
    w = n - 3
    mz = float((w * z).sum() / w.sum())
    se = 1.0 / math.sqrt(float(w.sum()))
    return {
        "rho": float(np.tanh(mz)),
        "ci_lo": float(np.tanh(mz - 1.96 * se)),
        "ci_hi": float(np.tanh(mz + 1.96 * se)),
        "n": int(n.sum()),
        "k": int(len(sub)),
    }


def zscore_cols(df: pd.DataFrame) -> pd.DataFrame:
    arr = df.to_numpy(float)
    mu = np.nanmean(arr, axis=1, keepdims=True)
    sd = np.nanstd(arr, axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = np.nan
    return pd.DataFrame((arr - mu) / sd, index=df.index, columns=df.columns)


def map_lee_to_symbols(lee: pd.DataFrame, target_genes: list[str]) -> pd.DataFrame:
    gmap = pd.read_csv(GENE_MAP, sep="\t")
    ens_to_sym = dict(zip(gmap["ensembl"], gmap["symbol"]))
    ens_target = [e for e, s in ens_to_sym.items() if s in set(target_genes)]
    sub = lee.loc[lee.index.intersection(ens_target)].copy()
    sub.index = [ens_to_sym[e] for e in sub.index]
    return sub.groupby(sub.index).mean()


def run_v15a_k2_overlay() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    print("\n=== v15A: K2 Panel-only calibrated overlay ===")
    tcga = pd.read_csv(TCGA_LOG2, sep="\t", index_col=0)
    labels = pd.read_csv(TCGA_LABELS, sep="\t")
    labels["DM"] = labels["cluster"].str.replace(r"_[A-Z]$", "", regex=True)
    labels["y_dm2"] = labels["cluster"].str.startswith("DM2").astype(int)
    common_samples = [s for s in labels["sample_id"] if s in tcga.columns]
    x_train = tcga.loc[PANEL_8, common_samples].T
    y = labels.set_index("sample_id").loc[common_samples, "y_dm2"].to_numpy(int)
    dm = labels.set_index("sample_id").loc[common_samples, "DM"]
    x_train_c = center_rows(x_train.to_numpy(float))

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    pipe = make_pipeline(
        StandardScaler(),
        LogisticRegression(C=1.0, max_iter=2000, random_state=42),
    )
    p_cv = cross_val_predict(pipe, x_train_c, y, cv=cv, method="predict_proba")[:, 1]
    cv_auc = float(roc_auc_score(y, p_cv))

    scaler = StandardScaler().fit(x_train_c)
    model = LogisticRegression(C=1.0, max_iter=2000, random_state=42).fit(scaler.transform(x_train_c), y)
    p_tcga_fit = model.predict_proba(scaler.transform(x_train_c))[:, 1]
    train_auc = float(roc_auc_score(y, p_tcga_fit))

    rows = []
    for sample, p, p_fit, dm_call in zip(common_samples, p_cv, p_tcga_fit, dm):
        rows.append({
            "cohort": "TCGA_THCA",
            "sample_id": sample,
            "group": dm_call,
            "p_DM2_cv": float(p),
            "p_DM2_full_model": float(p_fit),
            "score_for_overlay": float(p),
            "source": "TCGA 5-fold out-of-fold prediction",
        })

    lee = pd.read_csv(LEE_LOG2, sep="\t", index_col=0)
    lee_sub = map_lee_to_symbols(lee, PANEL_8)
    lee_x = lee_sub.reindex(PANEL_8).T
    lee_x_c = center_rows(lee_x.to_numpy(float))
    p_lee = model.predict_proba(scaler.transform(lee_x_c))[:, 1]
    lee_meta = pd.read_csv(LEE_META_PANEL, sep="\t")
    lee_meta = lee_meta.drop_duplicates("gsm").set_index("gsm")
    for sample, p in zip(lee_x.index, p_lee):
        tissue = str(lee_meta["tissue_type"].get(sample, "unknown")) if "tissue_type" in lee_meta else "unknown"
        histology = str(lee_meta["histology"].get(sample, "")) if "histology" in lee_meta else ""
        rows.append({
            "cohort": "Lee_GSE213647",
            "sample_id": sample,
            "group": tissue,
            "histology": histology,
            "p_DM2_cv": np.nan,
            "p_DM2_full_model": float(p),
            "score_for_overlay": float(p),
            "source": "TCGA-trained centered-profile model",
        })

    k2 = pd.read_csv(K2_TPM, sep="\t", index_col=0)
    xk = np.log2(k2[PANEL_8].to_numpy(float) + 1.0)
    p_k2 = model.predict_proba(scaler.transform(center_rows(xk)))[:, 1]
    k2_pred = pd.read_csv(K2_PRED, sep="\t") if K2_PRED.exists() else pd.DataFrame()
    k2_pred = k2_pred.set_index("run") if not k2_pred.empty else pd.DataFrame(index=k2.index)
    for sample, p in zip(k2.index, p_k2):
        prior_call = str(k2_pred["DM_call"].get(sample, "DM2" if p >= 0.5 else "DM1")) if "DM_call" in k2_pred else ("DM2" if p >= 0.5 else "DM1")
        rows.append({
            "cohort": "K2_PRJEB11591",
            "sample_id": sample,
            "group": prior_call,
            "p_DM2_cv": np.nan,
            "p_DM2_full_model": float(p),
            "score_for_overlay": float(p),
            "source": "K2 8-gene mini-index; centered-profile correction",
        })

    per_sample = pd.DataFrame(rows)
    per_sample.to_csv(OUT / "v15A_k2_panel_overlay_per_sample.tsv", sep="\t", index=False)

    summary_rows = []
    for (cohort, group), sub in per_sample.groupby(["cohort", "group"]):
        v = sub["score_for_overlay"].dropna()
        summary_rows.append({
            "cohort": cohort,
            "group": group,
            "n": int(len(v)),
            "median_p_DM2": float(v.median()) if len(v) else np.nan,
            "mean_p_DM2": float(v.mean()) if len(v) else np.nan,
            "q25_p_DM2": float(v.quantile(0.25)) if len(v) else np.nan,
            "q75_p_DM2": float(v.quantile(0.75)) if len(v) else np.nan,
            "pct_p_DM2_ge_0p5": float((v >= 0.5).mean() * 100) if len(v) else np.nan,
            "pct_p_DM2_ge_0p9": float((v >= 0.9).mean() * 100) if len(v) else np.nan,
        })
    by_group = pd.DataFrame(summary_rows).sort_values(["cohort", "group"])
    by_group.to_csv(OUT / "v15A_k2_panel_overlay_by_group.tsv", sep="\t", index=False)

    summary = {
        "tcga_centered_profile_cv_auc": cv_auc,
        "tcga_centered_profile_train_auc": train_auc,
        "tcga_training_n": int(len(y)),
        "tcga_training_dm1_dm2": f"{int((y == 0).sum())}:{int((y == 1).sum())}",
        "k2_n": int((per_sample["cohort"] == "K2_PRJEB11591").sum()),
        "k2_n_dm2_like_p_ge_0p5": int(((per_sample["cohort"] == "K2_PRJEB11591") & (per_sample["score_for_overlay"] >= 0.5)).sum()),
        "k2_pct_dm2_like_p_ge_0p5": float((per_sample.loc[per_sample["cohort"] == "K2_PRJEB11591", "score_for_overlay"] >= 0.5).mean() * 100),
        "k2_median_p_dm2": float(per_sample.loc[per_sample["cohort"] == "K2_PRJEB11591", "score_for_overlay"].median()),
        "calibration_boundary": "K2 uses only 8 genes and no MAPK genes; it is score-distribution evidence only, not a MAPK x Panel cohort.",
    }
    (OUT / "v15A_k2_panel_overlay_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2)[:1200])
    return per_sample, by_group, summary


def extract_log_cp10k(adata, genes: list[str]) -> tuple[pd.DataFrame, list[str]]:
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return pd.DataFrame(index=adata.obs_names), []
    idx = [int(np.where(adata.var_names == g)[0][0]) for g in present]
    x_sub = adata.X[:, idx]
    if sparse.issparse(x_sub):
        x_sub = x_sub.toarray()
    else:
        x_sub = np.asarray(x_sub)
    lib = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)
    lib[~np.isfinite(lib) | (lib <= 0)] = np.nan
    log = np.log1p(x_sub / lib[:, None] * 1e4)
    return pd.DataFrame(log, index=adata.obs_names, columns=present), present


def partial_spearman(x: pd.Series, y: pd.Series, cov: pd.DataFrame) -> tuple[float, float, int]:
    df = pd.concat([x.rename("x"), y.rename("y"), cov], axis=1).dropna()
    if len(df) < 8:
        return (np.nan, np.nan, len(df))
    rx = stats.rankdata(df["x"].to_numpy(float))
    ry = stats.rankdata(df["y"].to_numpy(float))
    c_rank = np.column_stack([np.ones(len(df))] + [stats.rankdata(df[c].to_numpy(float)) for c in cov.columns])
    bx, *_ = np.linalg.lstsq(c_rank, rx, rcond=None)
    by, *_ = np.linalg.lstsq(c_rank, ry, rcond=None)
    ex = rx - c_rank @ bx
    ey = ry - c_rank @ by
    r, p = stats.pearsonr(ex, ey)
    return (float(r), float(p), int(len(df)))


def run_v15b_spatial_mapk_panel() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    print("\n=== v15B: GSE250521 spatial MAPK x Panel ===")
    score_meta = pd.read_csv(SPATIAL_SCORES, sep="\t")
    score_meta["sample_spot"] = score_meta["sample_id"].astype(str) + "::" + score_meta["spot_id"].astype(str)
    spot_rows = []
    slide_rows = []

    for h5 in sorted(SPATIAL_DIR.glob("GSM*/GSM*.raw.h5ad")):
        sample = h5.parent.name
        a = ad.read_h5ad(h5)
        obs = a.obs.copy()
        obs["spot_id"] = obs.index.astype(str)
        stage = str(obs["stage"].iloc[0]) if "stage" in obs.columns else sample.split("_", 1)[1].split("-")[0]
        mat, present = extract_log_cp10k(a, sorted(set(PANEL_8_CANON + MAPK_9 + TDS_16)))
        if mat.empty:
            continue
        z = zscore_cols(mat.T).T
        panel_present = [g for g in PANEL_8_CANON if g in z.columns]
        mapk_present = [g for g in MAPK_9 if g in z.columns]
        tds_present = [g for g in TDS_16 if g in z.columns]
        df = pd.DataFrame({
            "sample_id": sample,
            "stage": stage,
            "spot_id": z.index.astype(str),
            "panel8_score": z[panel_present].mean(axis=1) if panel_present else np.nan,
            "mapk9_score": z[mapk_present].mean(axis=1) if mapk_present else np.nan,
            "tds16_score": z[tds_present].mean(axis=1) if tds_present else np.nan,
            "n_panel_genes": len(panel_present),
            "n_mapk_genes": len(mapk_present),
            "n_tds16_genes": len(tds_present),
        })
        df["sample_spot"] = df["sample_id"].astype(str) + "::" + df["spot_id"].astype(str)
        meta_cols = [
            "sample_spot", "total_counts", "n_genes_by_counts", "pct_counts_mt",
            "array_row", "array_col", "in_tissue", "DM1_like_score", "RAI_8_score", "TDS_like_score",
        ]
        df = df.merge(score_meta[[c for c in meta_cols if c in score_meta.columns]], on="sample_spot", how="left")
        df = df[df.get("in_tissue", 1).fillna(1).astype(int) == 1].copy()
        df["log_total_counts"] = np.log1p(df["total_counts"].astype(float))
        df["log_n_genes"] = np.log1p(df["n_genes_by_counts"].astype(float))
        spot_rows.append(df)

        sub = df.dropna(subset=["mapk9_score", "panel8_score"])
        rho, p = stats.spearmanr(sub["mapk9_score"], sub["panel8_score"])
        rr, rp, rn = partial_spearman(
            sub["mapk9_score"],
            sub["panel8_score"],
            sub[["log_total_counts", "log_n_genes", "pct_counts_mt"]],
        )
        lo, hi = fisher_ci(float(rho), len(sub))
        slide_rows.append({
            "sample_id": sample,
            "stage": stage,
            "n": int(len(sub)),
            "spearman_rho": float(rho),
            "spearman_p": float(p),
            "ci_lo": lo,
            "ci_hi": hi,
            "partial_spearman_rho_qc": rr,
            "partial_spearman_p_qc": rp,
            "partial_n_qc": rn,
            "panel_genes_used": len(panel_present),
            "mapk_genes_used": len(mapk_present),
            "tds16_genes_used": len(tds_present),
        })
        print(f"  {sample:18s} {stage:4s} n={len(sub):5d} rho={rho:+.3f} p={p:.1e} partial={rr:+.3f}")

    per_spot = pd.concat(spot_rows, ignore_index=True)
    per_slide = pd.DataFrame(slide_rows)
    per_slide["spearman_fdr"] = bh_fdr(per_slide["spearman_p"])
    per_slide.to_csv(OUT / "v15B_spatial_mapk_panel_per_slide.tsv", sep="\t", index=False)
    per_spot.to_csv(OUT / "v15B_spatial_mapk_panel_per_spot.tsv.gz", sep="\t", index=False, compression="gzip")

    stage_rows = []
    for stage, sub in per_spot.dropna(subset=["mapk9_score", "panel8_score"]).groupby("stage"):
        rho, p = stats.spearmanr(sub["mapk9_score"], sub["panel8_score"])
        rr, rp, rn = partial_spearman(
            sub["mapk9_score"],
            sub["panel8_score"],
            sub[["log_total_counts", "log_n_genes", "pct_counts_mt"]],
        )
        lo, hi = fisher_ci(float(rho), len(sub))
        stage_rows.append({
            "stage": stage,
            "n": int(len(sub)),
            "n_slides": int(sub["sample_id"].nunique()),
            "spearman_rho": float(rho),
            "spearman_p": float(p),
            "ci_lo": lo,
            "ci_hi": hi,
            "partial_spearman_rho_qc": rr,
            "partial_spearman_p_qc": rp,
            "partial_n_qc": rn,
        })
    by_stage = pd.DataFrame(stage_rows)
    by_stage.to_csv(OUT / "v15B_spatial_mapk_panel_by_stage.tsv", sep="\t", index=False)

    tumor_slides = per_slide[per_slide["stage"].isin(["PTC", "LPTC", "ATC"])].copy()
    pool_raw = fisher_pool(tumor_slides, "spearman_rho")
    pool_partial = fisher_pool(tumor_slides.rename(columns={"partial_spearman_rho_qc": "rho"}), "rho")
    summary = {
        "n_spots_total": int(len(per_spot)),
        "n_spots_tumor": int(per_spot[per_spot["stage"].isin(["PTC", "LPTC", "ATC"])].shape[0]),
        "n_slides_total": int(per_slide.shape[0]),
        "n_tumor_slides": int(tumor_slides.shape[0]),
        "tumor_slide_fixed_effect_rho": pool_raw,
        "tumor_slide_qc_partial_fixed_effect_rho": pool_partial,
        "stage_summary": by_stage.to_dict(orient="records"),
        "direction_check": "Observed per-spot MAPK x Panel association is positive, not the expected anti-correlation; QC-partial association is much smaller.",
        "boundary": "Spatial result is a tissue-localization stress test, not a mechanism lock. It does not support using GSE250521 as a spatial MAPK-to-panel anti-correlation confirmation.",
    }
    (OUT / "v15B_spatial_mapk_panel_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2)[:1600])
    return per_spot, per_slide, by_stage, summary


MEK = {"selumetinib", "trametinib", "azd-0364", "azd0364", "cobimetinib", "binimetinib", "pimasertib", "refametinib", "pd0325901", "pd-0325901", "mirdametinib", "tak-733"}
BRAF_RAF = {"dabrafenib", "vemurafenib", "encorafenib", "raf265", "raf709", "az-628", "az628", "lifirafenib", "tovorafenib", "sorafenib", "lxh254"}
ERK = {"sch772984", "sch-772984", "ulixertinib", "bvd-523", "bvd523", "ly3214996", "ravoxertinib", "gdc-0994", "gdc0994", "vx-11e", "vx11e", "ml786", "cc-90003", "cc90003"}
PI3K = {"gdc-0077", "gdc0077", "inavolisib"}


def classify_drug(name: str, h9_class: str = "") -> str:
    h9 = str(h9_class).strip()
    if h9:
        return h9
    n = str(name).lower().replace(" ", "-")
    if n in MEK:
        return "MEK"
    if n in BRAF_RAF:
        return "BRAF/RAF"
    if n in ERK:
        return "ERK"
    if n in PI3K:
        return "PI3K"
    return "non-MAPK/unknown"


def run_v15c_prism_depmap_overlay() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    print("\n=== v15C: PRISM/DepMap MAPK inhibitor overlay ===")
    prism = pd.read_csv(PRISM_TABLE, sep="\t")
    prism = prism.sort_values(["fdr", "cohens_d"]).reset_index(drop=True)
    prism["rank_fdr"] = np.arange(1, len(prism) + 1)
    if H9_RANKING.exists():
        h9 = pd.read_csv(H9_RANKING, sep="\t", usecols=["broad_id", "name", "mapk_class"])
        h9 = h9.drop_duplicates("broad_id")
        prism = prism.merge(h9[["broad_id", "mapk_class"]], on="broad_id", how="left")
    else:
        prism["mapk_class"] = ""
    prism["mapk_class"] = [classify_drug(n, c) for n, c in zip(prism["name"], prism["mapk_class"].fillna(""))]
    prism["is_canonical_mapk"] = prism["mapk_class"].isin(["MEK", "BRAF/RAF", "ERK"])
    prism["is_mapk_pi3k_axis"] = prism["mapk_class"].isin(["MEK", "BRAF/RAF", "ERK", "PI3K"])
    prism.to_csv(OUT / "v15C_prism_annotated_drugs.tsv", sep="\t", index=False)

    top15 = prism.head(15).copy()
    top15.to_csv(OUT / "v15C_prism_top15_annotated.tsv", sep="\t", index=False)
    sig = prism[prism["fdr"] < 0.05].copy()
    focus_ids = sig.loc[sig["is_canonical_mapk"], "broad_id"].tolist()

    cell_scores = pd.DataFrame()
    if focus_ids:
        print(f"  loading PRISM raw LFC for {len(focus_ids)} canonical MAPK FDR<0.05 drugs")
        raw_iter = pd.read_csv(PRISM_RAW, usecols=["row_id", "broad_id", "LFC"], chunksize=250000)
        chunks = []
        focus = set(focus_ids)
        for chunk in raw_iter:
            sub = chunk[chunk["broad_id"].isin(focus)].copy()
            if len(sub):
                sub["ModelID"] = sub["row_id"].str.split("::").str[0]
                chunks.append(sub[["ModelID", "broad_id", "LFC"]])
        raw = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=["ModelID", "broad_id", "LFC"])
        cell_drug = raw.groupby(["ModelID", "broad_id"])["LFC"].mean().reset_index()
        wide = cell_drug.pivot(index="ModelID", columns="broad_id", values="LFC")
        wide["mapk_inhibitor_mean_LFC"] = wide.mean(axis=1)
        cells = pd.read_csv(CELLS_DM1, sep="\t", usecols=["ModelID", "StrippedCellLineName", "OncotreeLineage", "DM1_like_score"])
        cell_scores = cells.merge(wide[["mapk_inhibitor_mean_LFC"]].reset_index(), on="ModelID", how="inner")
        q1, q2 = cell_scores["DM1_like_score"].quantile([1 / 3, 2 / 3])
        cell_scores["dm1_bin"] = np.where(
            cell_scores["DM1_like_score"] <= q1,
            "low",
            np.where(cell_scores["DM1_like_score"] >= q2, "high", "mid"),
        )
        valid = cell_scores[["DM1_like_score", "mapk_inhibitor_mean_LFC"]].dropna()
        if len(valid) >= 3 and valid["DM1_like_score"].std(ddof=1) > 0 and valid["mapk_inhibitor_mean_LFC"].std(ddof=1) > 0:
            rho, p = stats.spearmanr(valid["DM1_like_score"], valid["mapk_inhibitor_mean_LFC"])
        else:
            rho = p = np.nan
        high = cell_scores.loc[cell_scores["dm1_bin"] == "high", "mapk_inhibitor_mean_LFC"].dropna().to_numpy()
        low = cell_scores.loc[cell_scores["dm1_bin"] == "low", "mapk_inhibitor_mean_LFC"].dropna().to_numpy()
        cell_scores.to_csv(OUT / "v15C_prism_mapk_cellline_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    else:
        rho = p = d = np.nan
        high = low = np.array([])

    d = cohen_d(high, low) if len(high) and len(low) else np.nan
    dep = pd.read_csv(DEPMAP_DEP, sep="\t")
    dep = dep.drop_duplicates("gene").sort_values("cohens_d")
    dep.to_csv(OUT / "v15C_depmap_dependency_overlay.tsv", sep="\t", index=False)

    summary = {
        "n_prism_drugs_tested": int(len(prism)),
        "n_prism_fdr05": int((prism["fdr"] < 0.05).sum()),
        "n_prism_fdr05_dm1_high_selective": int(((prism["fdr"] < 0.05) & (prism["cohens_d"] < 0)).sum()),
        "n_fdr05_canonical_mapk": int(sig["is_canonical_mapk"].sum()),
        "n_fdr05_mapk_pi3k_axis": int(sig["is_mapk_pi3k_axis"].sum()),
        "top_drug": prism.iloc[0][["name", "cohens_d", "fdr", "mapk_class"]].to_dict(),
        "mapk_focus_drug_names": sig.loc[sig["is_canonical_mapk"], "name"].tolist(),
        "mapk_focus_celllines_n": int(len(cell_scores)) if len(cell_scores) else 0,
        "mapk_focus_celllines_nonmissing_n": int(cell_scores["mapk_inhibitor_mean_LFC"].notna().sum()) if len(cell_scores) else 0,
        "mapk_focus_spearman_dm1_vs_lfc": float(rho) if np.isfinite(rho) else np.nan,
        "mapk_focus_spearman_p": float(p) if np.isfinite(p) else np.nan,
        "mapk_focus_high_vs_low_cohens_d_lfc": float(d) if np.isfinite(d) else np.nan,
        "depmap_top_dependencies": dep.head(5).to_dict(orient="records"),
        "boundary": "PRISM/DepMap supports MAPK-pathway vulnerability in DM1-high models, but it does not measure expression restoration after MEK/RAF/ERK treatment.",
    }
    (OUT / "v15C_prism_depmap_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2)[:1800])
    return prism, cell_scores, dep, summary


def plot_v15(k2: pd.DataFrame, spatial_slide: pd.DataFrame, prism: pd.DataFrame, dep: pd.DataFrame, summaries: dict) -> None:
    print("\n=== v15 figure ===")
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.26)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # A. K2 overlay, score-distribution only
    groups = [
        ("TCGA_THCA", "DM1", "#d8453a", "TCGA DM1"),
        ("TCGA_THCA", "DM2", "#2c7fb8", "TCGA DM2"),
        ("Lee_GSE213647", "PTC", "#7a4ca8", "Lee PTC"),
        ("K2_PRJEB11591", "DM2", "#1c8e6d", "K2 predicted DM2"),
        ("K2_PRJEB11591", "DM1", "#ff9b3d", "K2 predicted DM1"),
    ]
    bins = np.linspace(0, 1, 32)
    for cohort, group, color, label in groups:
        vals = k2[(k2["cohort"] == cohort) & (k2["group"] == group)]["score_for_overlay"].dropna()
        if len(vals):
            ax_a.hist(vals, bins=bins, histtype="step", lw=1.8, color=color, density=True, label=f"{label} n={len(vals)}")
    ax_a.axvline(0.5, color="#555", ls="--", lw=0.8)
    ax_a.set_xlabel("TCGA-centered profile classifier P(DM2)")
    ax_a.set_ylabel("Density")
    ax_a.set_title("A. K2 Panel-only overlay: score-distribution evidence, not MAPK correlation", loc="left", fontsize=10, weight="bold")
    ax_a.legend(fontsize=7.5, frameon=False)
    ax_a.grid(axis="y", alpha=0.25)

    # B. Spatial per-slide forest
    stage_order = ["PT", "PTC", "LPTC", "ATC"]
    colors = {"PT": "#8796a8", "PTC": "#2c7fb8", "LPTC": "#ff9b3d", "ATC": "#d8453a"}
    s = spatial_slide.copy()
    s["stage_rank"] = s["stage"].map({k: i for i, k in enumerate(stage_order)})
    s = s.sort_values(["stage_rank", "sample_id"])
    y = np.arange(len(s))
    for i, (_, r) in enumerate(s.iterrows()):
        c = colors.get(r["stage"], "#888")
        ax_b.plot([r["ci_lo"], r["ci_hi"]], [i, i], color=c, lw=1.4, alpha=0.8)
        ax_b.plot(r["spearman_rho"], i, "s", color=c, mec="black", mew=0.4, ms=5 + 5 * math.sqrt(min(r["n"], 5000) / 5000))
    pool = summaries["spatial"]["tumor_slide_fixed_effect_rho"]
    py = len(s) + 0.4
    ax_b.plot([pool["ci_lo"], pool["ci_hi"]], [py, py], color="black", lw=2)
    ax_b.plot(pool["rho"], py, "D", color="gold", mec="black", mew=0.6, ms=9)
    ax_b.axvline(0, color="#777", ls="--", lw=0.8)
    ax_b.set_yticks(list(y) + [py])
    ax_b.set_yticklabels([f"{r.sample_id.replace('GSM', '')} ({r.stage}, n={int(r.n)})" for _, r in s.iterrows()] + [f"Tumor pooled (k={pool['k']}, n={pool['n']})"], fontsize=7)
    ax_b.invert_yaxis()
    xmin = min(-0.12, float(np.nanmin(s["ci_lo"])) - 0.04, pool["ci_lo"] - 0.04)
    xmax = max(0.70, float(np.nanmax(s["ci_hi"])) + 0.04, pool["ci_hi"] + 0.04)
    ax_b.set_xlim(xmin, xmax)
    ax_b.set_xlabel("Per-spot Spearman rho: MAPK-9 x Panel-8")
    ax_b.set_title("B. Spatial Lu 2023 GSE250521: positive MAPK x Panel co-localization, QC-sensitive", loc="left", fontsize=10, weight="bold")
    ax_b.grid(axis="x", alpha=0.25)

    # C. PRISM top hits
    top = prism.head(15).sort_values("cohens_d", ascending=True)
    class_color = {
        "MEK": "#d8453a",
        "BRAF/RAF": "#ff9b3d",
        "ERK": "#7a4ca8",
        "PI3K": "#2c7fb8",
        "non-MAPK/unknown": "#8796a8",
    }
    y = np.arange(len(top))
    ax_c.barh(y, top["cohens_d"], color=[class_color.get(c, "#8796a8") for c in top["mapk_class"]], edgecolor="black", lw=0.4)
    ax_c.set_yticks(y)
    ax_c.set_yticklabels([f"{r['name']}  ({r['mapk_class']})" for _, r in top.iterrows()], fontsize=7.5)
    ax_c.axvline(0, color="#777", ls="--", lw=0.8)
    ax_c.set_xlabel("Cohen's d on PRISM LFC (DM1-high - DM1-low)")
    ax_c.set_title("C. PRISM: DM1-high cell lines are selectively sensitive to MAPK-axis inhibitors", loc="left", fontsize=10, weight="bold")
    ax_c.grid(axis="x", alpha=0.25)

    # D. DepMap dependency overlay
    dep_show = dep.head(12).sort_values("cohens_d")
    y = np.arange(len(dep_show))
    ax_d.barh(y, dep_show["cohens_d"], color=["#d8453a" if x < 0 else "#2c7fb8" for x in dep_show["cohens_d"]], edgecolor="black", lw=0.4)
    ax_d.set_yticks(y)
    ax_d.set_yticklabels([f"{r['gene']}  p={r['p']:.1e}" for _, r in dep_show.iterrows()], fontsize=8)
    ax_d.axvline(0, color="#777", ls="--", lw=0.8)
    ax_d.set_xlabel("CRISPR GeneEffect Cohen's d (DM1-high - DM1-low)")
    ax_d.set_title("D. DepMap: MYC/NAMPT vulnerabilities; thyroid TFs are not the therapeutic lever", loc="left", fontsize=10, weight="bold")
    ax_d.grid(axis="x", alpha=0.25)

    fig.suptitle(
        "Supplementary Figure SX_v15. Paper 1 Fig 8 mechanism extensions: K2 score-only overlay, spatial MAPK x Panel, and PRISM/DepMap therapy",
        fontsize=13,
        weight="bold",
        y=0.995,
    )
    fig.savefig(OUT / "Fig_SX_v15_mechanism_extensions.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "Fig_SX_v15_mechanism_extensions.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUT}/Fig_SX_v15_mechanism_extensions.{{png,pdf}}")


def main() -> None:
    k2, k2_group, s_a = run_v15a_k2_overlay()
    spatial_spot, spatial_slide, spatial_stage, s_b = run_v15b_spatial_mapk_panel()
    prism, cell_scores, dep, s_c = run_v15c_prism_depmap_overlay()
    summaries = {"k2": s_a, "spatial": s_b, "prism": s_c}
    (OUT / "v15_mechanism_extensions_summary.json").write_text(json.dumps(summaries, indent=2, default=str))
    plot_v15(k2, spatial_slide, prism, dep, summaries)
    print("\n=== DONE v15 ===")


if __name__ == "__main__":
    main()
