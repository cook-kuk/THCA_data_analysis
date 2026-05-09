"""
Pathology finding scout: spatial-DM1 morphology transfer.

Question
--------
Can a morphology axis trained from GSE250521 spatial H&E UNI tiles to predict
spot-level DM1_like_score transfer to TCGA-THCA WSI tile embeddings?

This is intentionally a scout, not a manuscript claim generator. It reports
both signal and confound checks so weak/unsafe findings can be killed quickly.

Outputs
-------
analysis_supp/pathology_scout_2026_05_09/
  PATHOLOGY_SCOUT_REPORT.md
  PATHOLOGY_SCOUT_SUMMARY.json
  spatial_spot_loso_predictions.tsv
  spatial_slide_summary.tsv
  tcga_slide_pathology_scores.tsv
  tcga_association_tests.tsv
  tcga_confound_residual_tests.tsv
  fig_pathology_scout_spatial_loso.{png,pdf}
  fig_pathology_scout_tcga_projection.{png,pdf}
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler


REPO = Path(os.environ.get("THCA_REPO_ROOT", Path.cwd()))
P2_ROOT = Path(
    os.environ.get(
        "THCA_P2_ROOT",
        REPO / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07",
    )
)
SPATIAL_DIR = P2_ROOT / "phase1_gse250521"
SPATIAL_SCORE_PATH = REPO / "project/results/01_spatial_score/all_spots_scored.tsv.gz"
TCGA_MANIFEST = P2_ROOT / "phase2_tcga_clam/slide_manifest.tsv"
TCGA_FEATURE_DIR = P2_ROOT / "phase2_tcga_clam_UNI/features"
TCGA_UNI_PREDS = P2_ROOT / "phase2_tcga_clam_UNI/clam_per_slide_predictions.tsv"
TCGA_UNI_LOTO = P2_ROOT / "analysis_supp/audit_uni_loto/uni_loto_oof_predictions.tsv"
TCGA_META = REPO / "project/metadata/sample_master_v3.tsv"
OUT_DIR = P2_ROOT / "analysis_supp/pathology_scout_2026_05_09"

N_PCS = 32
RIDGE_ALPHA = 100.0
RANDOM_STATE = 17


def load_pt(path: Path) -> np.ndarray:
    try:
        obj = torch.load(path, map_location="cpu", weights_only=True)
    except TypeError:
        obj = torch.load(path, map_location="cpu")
    if isinstance(obj, torch.Tensor):
        arr = obj.detach().cpu().numpy()
    elif isinstance(obj, dict):
        for key in ("features", "embeddings", "x"):
            if key in obj:
                val = obj[key]
                arr = val.detach().cpu().numpy() if isinstance(val, torch.Tensor) else np.asarray(val)
                break
        else:
            raise ValueError(f"No feature tensor key found in {path}")
    else:
        arr = np.asarray(obj)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D feature array in {path}, got {arr.shape}")
    return arr.astype(np.float32, copy=False)


def infer_stage(slide: str) -> str:
    if "_N-" in slide:
        return "PT"
    if "_PTC-" in slide:
        return "PTC"
    if "_LPTC-" in slide:
        return "LPTC"
    if "_ATC-" in slide:
        return "ATC"
    return "unknown"


def safe_spearman(x: Iterable[float], y: Iterable[float]) -> tuple[float, float]:
    xx = np.asarray(list(x), dtype=float)
    yy = np.asarray(list(y), dtype=float)
    mask = np.isfinite(xx) & np.isfinite(yy)
    if mask.sum() < 3 or np.nanstd(xx[mask]) == 0 or np.nanstd(yy[mask]) == 0:
        return float("nan"), float("nan")
    rho, p = spearmanr(xx[mask], yy[mask])
    return float(rho), float(p)


def auc_or_nan(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y).astype(int)
    score = np.asarray(score).astype(float)
    mask = np.isfinite(score)
    y = y[mask]
    score = score[mask]
    if len(np.unique(y)) != 2:
        return float("nan")
    return float(roc_auc_score(y, score))


def mwu_test(pos: np.ndarray, neg: np.ndarray) -> float:
    pos = np.asarray(pos, dtype=float)
    neg = np.asarray(neg, dtype=float)
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    if len(pos) < 2 or len(neg) < 2:
        return float("nan")
    return float(mannwhitneyu(pos, neg, alternative="two-sided").pvalue)


def train_predict_model(x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray) -> np.ndarray:
    scaler = StandardScaler()
    x_train_z = scaler.fit_transform(x_train)
    x_test_z = scaler.transform(x_test)
    n_comp = min(N_PCS, x_train_z.shape[0] - 1, x_train_z.shape[1])
    pca = PCA(n_components=n_comp, random_state=RANDOM_STATE)
    x_train_pc = pca.fit_transform(x_train_z)
    x_test_pc = pca.transform(x_test_z)
    y_mean = float(np.mean(y_train))
    y_sd = float(np.std(y_train, ddof=0))
    if y_sd == 0:
        y_sd = 1.0
    ridge = Ridge(alpha=RIDGE_ALPHA)
    ridge.fit(x_train_pc, (y_train - y_mean) / y_sd)
    return ridge.predict(x_test_pc)


def fit_final_spatial_model(x: np.ndarray, y: np.ndarray):
    scaler = StandardScaler()
    x_z = scaler.fit_transform(x)
    n_comp = min(N_PCS, x_z.shape[0] - 1, x_z.shape[1])
    pca = PCA(n_components=n_comp, random_state=RANDOM_STATE)
    x_pc = pca.fit_transform(x_z)
    y_mean = float(np.mean(y))
    y_sd = float(np.std(y, ddof=0))
    if y_sd == 0:
        y_sd = 1.0
    ridge = Ridge(alpha=RIDGE_ALPHA)
    ridge.fit(x_pc, (y - y_mean) / y_sd)
    pred_train = ridge.predict(x_pc)
    return scaler, pca, ridge, pred_train


def project_model(x: np.ndarray, scaler: StandardScaler, pca: PCA, ridge: Ridge) -> np.ndarray:
    return ridge.predict(pca.transform(scaler.transform(x)))


def load_spatial() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    embeds = np.load(SPATIAL_DIR / "uni_embeddings_size224.npz")["embeddings"].astype(np.float32)
    meta = pd.read_csv(SPATIAL_DIR / "uni_embed_metadata_size224.tsv", sep="\t")
    scores = pd.read_csv(SPATIAL_SCORE_PATH, sep="\t")
    meta = meta.reset_index().rename(columns={"index": "embed_idx"})
    meta["key"] = meta["slide"].astype(str) + "_" + meta["spot_id"].astype(str)
    scores["key"] = scores["sample_id"].astype(str) + "_" + scores["spot_id"].astype(str)
    keep_cols = [
        "key",
        "sample_id",
        "stage",
        "in_tissue",
        "DM1_like_score",
        "RAI_8_score",
        "TDS_like_score",
        "CAF_ECM_score",
        "EMT_score",
        "Hypoxia_score",
        "Proliferation_score",
        "Epithelial_score",
    ]
    merged = meta.merge(scores[keep_cols], on="key", how="inner")
    merged["stage"] = merged["stage"].fillna(merged["slide"].map(infer_stage))
    mask = (merged["in_tissue"] == 1) & merged["DM1_like_score"].notna()
    merged = merged.loc[mask].copy()
    x = embeds[merged["embed_idx"].to_numpy()]
    y = merged["DM1_like_score"].astype(float).to_numpy()
    y = (y - np.nanmean(y)) / (np.nanstd(y) if np.nanstd(y) > 0 else 1.0)
    return merged.reset_index(drop=True), x, y


def spatial_leave_one_slide_out(meta: pd.DataFrame, x: np.ndarray, y: np.ndarray) -> pd.DataFrame:
    rows = []
    for held in sorted(meta["slide"].unique()):
        train = meta["slide"].to_numpy() != held
        test = ~train
        pred = train_predict_model(x[train], y[train], x[test])
        part = meta.loc[test, ["slide", "spot_id", "stage", "DM1_like_score"]].copy()
        part["target_DM1_like_z"] = y[test]
        part["pred_morph_DM1_z"] = pred
        rows.append(part)
    return pd.concat(rows, ignore_index=True)


def summarize_spatial(pred: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for slide, sub in pred.groupby("slide"):
        rho, p = safe_spearman(sub["target_DM1_like_z"], sub["pred_morph_DM1_z"])
        rows.append(
            {
                "slide": slide,
                "stage": sub["stage"].iloc[0],
                "n_spots": int(len(sub)),
                "mean_target_DM1_like_z": float(sub["target_DM1_like_z"].mean()),
                "mean_pred_morph_DM1_z": float(sub["pred_morph_DM1_z"].mean()),
                "rho_spot_target_vs_pred": rho,
                "p_spot_target_vs_pred": p,
            }
        )
    return pd.DataFrame(rows)


def load_tcga_metadata() -> pd.DataFrame:
    man = pd.read_csv(TCGA_MANIFEST, sep="\t")
    man["tss"] = man["submitter_id"].str.split("-").str[1]
    man["case_barcode"] = man["submitter_id"].str[:12]
    pred = pd.read_csv(TCGA_UNI_PREDS, sep="\t")
    pred = pred.rename(columns={"slide": "file_id", "prob_DM1": "prob_DM1_uni_oof"})
    out = man.merge(pred[["file_id", "fold", "label", "prob_DM1_uni_oof"]], on="file_id", how="left")
    if TCGA_UNI_LOTO.exists():
        loto = pd.read_csv(TCGA_UNI_LOTO, sep="\t")
        loto = loto.rename(columns={"slide": "file_id"})
        cols = [c for c in ["file_id", "prob_DM1_uni_loto"] if c in loto.columns]
        out = out.merge(loto[cols], on="file_id", how="left")
    meta = pd.read_csv(TCGA_META, sep="\t")
    meta["case_barcode"] = meta["sample_id"].astype(str).str[:12]
    meta = meta[(meta["dataset"] == "TCGA-THCA") & (meta["normal_vs_tumor"] == "tumor")].copy()
    meta = meta.drop_duplicates("case_barcode", keep="first")
    cols = [
        "case_barcode",
        "sample_id",
        "histology_subtype",
        "driver_anchor",
        "molecular_subtype",
        "tds_group",
        "dediff_flag",
        "aggressive_flag",
        "ajcc_stage_group",
        "sex",
        "age",
        "tds_score",
        "dedifferentiation_proxy_score",
        "brs_like_surrogate_score",
        "v3_anchor",
    ]
    out = out.merge(meta[cols], on="case_barcode", how="left")
    return out


def project_tcga(
    tcga: pd.DataFrame, scaler: StandardScaler, pca: PCA, ridge: Ridge, spatial_pred: np.ndarray
) -> pd.DataFrame:
    q75 = float(np.quantile(spatial_pred, 0.75))
    q90 = float(np.quantile(spatial_pred, 0.90))
    rows = []
    for row in tcga.itertuples(index=False):
        fpath = TCGA_FEATURE_DIR / f"{row.file_id}.pt"
        if not fpath.exists():
            continue
        feats = load_pt(fpath)
        pred = project_model(feats, scaler, pca, ridge)
        rec = row._asdict()
        rec.update(
            {
                "n_tiles": int(feats.shape[0]),
                "morph_dm1_mean": float(np.mean(pred)),
                "morph_dm1_median": float(np.median(pred)),
                "morph_dm1_sd": float(np.std(pred)),
                "morph_dm1_q90": float(np.quantile(pred, 0.90)),
                "morph_dm1_frac_spatial_q75": float(np.mean(pred >= q75)),
                "morph_dm1_frac_spatial_q90": float(np.mean(pred >= q90)),
            }
        )
        rows.append(rec)
    return pd.DataFrame(rows)


def association_tests(tcga: pd.DataFrame) -> pd.DataFrame:
    score_cols = [
        "morph_dm1_mean",
        "morph_dm1_q90",
        "morph_dm1_frac_spatial_q75",
        "morph_dm1_frac_spatial_q90",
    ]
    tests = [
        ("DM1_vs_DM2", "dm", ["DM1"], ["DM2"]),
        ("DM1_vs_nonDM1", "dm", ["DM1"], ["DM2", "not_DM"]),
        ("RAS_like_vs_BRAF_like", "molecular_subtype", ["RAS_like"], ["BRAF_like"]),
        ("FVPTC_vs_cPTC", "histology_subtype", ["FVPTC"], ["cPTC"]),
        ("TDS_high_vs_low", "tds_group", ["high"], ["low"]),
        ("Stage_IIIplus_vs_I", "ajcc_stage_group", ["Stage III", "Stage IVA", "Stage IVB", "Stage IVC"], ["Stage I"]),
    ]
    rows = []
    for score in score_cols:
        for test_name, col, pos_vals, neg_vals in tests:
            if col not in tcga.columns:
                continue
            sub = tcga[tcga[col].isin(pos_vals + neg_vals) & tcga[score].notna()].copy()
            if sub.empty:
                continue
            y = sub[col].isin(pos_vals).astype(int).to_numpy()
            s = sub[score].astype(float).to_numpy()
            pos = s[y == 1]
            neg = s[y == 0]
            rows.append(
                {
                    "score": score,
                    "test": test_name,
                    "n": int(len(sub)),
                    "n_pos": int(y.sum()),
                    "n_neg": int(len(y) - y.sum()),
                    "median_pos": float(np.nanmedian(pos)) if len(pos) else float("nan"),
                    "median_neg": float(np.nanmedian(neg)) if len(neg) else float("nan"),
                    "delta_median_pos_minus_neg": float(np.nanmedian(pos) - np.nanmedian(neg))
                    if len(pos) and len(neg)
                    else float("nan"),
                    "auc_pos_high": auc_or_nan(y, s),
                    "mannwhitney_p": mwu_test(pos, neg),
                }
            )
    for score in score_cols + ["prob_DM1_uni_oof", "prob_DM1_uni_loto"]:
        if score not in tcga.columns:
            continue
        for target in ["tds_score", "dedifferentiation_proxy_score", "prob_DM1_uni_oof", "prob_DM1_uni_loto"]:
            if target == score or target not in tcga.columns:
                continue
            sub = tcga[[score, target]].dropna()
            if len(sub) < 8:
                continue
            rho, p = safe_spearman(sub[score], sub[target])
            rows.append(
                {
                    "score": score,
                    "test": f"spearman_vs_{target}",
                    "n": int(len(sub)),
                    "n_pos": "",
                    "n_neg": "",
                    "median_pos": "",
                    "median_neg": "",
                    "delta_median_pos_minus_neg": "",
                    "auc_pos_high": "",
                    "mannwhitney_p": "",
                    "spearman_rho": rho,
                    "spearman_p": p,
                }
            )
    return pd.DataFrame(rows)


def residualize_score(df: pd.DataFrame, score: str, covars: list[str]) -> np.ndarray:
    sub = df.copy()
    x_parts = []
    for cov in covars:
        if cov not in sub.columns:
            continue
        d = pd.get_dummies(sub[cov].fillna("missing").astype(str), prefix=cov, drop_first=True)
        if d.shape[1]:
            x_parts.append(d)
    if not x_parts:
        return sub[score].astype(float).to_numpy()
    x = pd.concat(x_parts, axis=1)
    y = sub[score].astype(float).to_numpy()
    mask = np.isfinite(y)
    resid = np.full(len(sub), np.nan)
    if mask.sum() <= x.shape[1] + 3:
        return resid
    lr = LinearRegression()
    lr.fit(x.loc[mask].to_numpy(dtype=float), y[mask])
    resid[mask] = y[mask] - lr.predict(x.loc[mask].to_numpy(dtype=float))
    return resid


def confound_residual_tests(tcga: pd.DataFrame) -> pd.DataFrame:
    score_cols = ["morph_dm1_mean", "morph_dm1_q90", "morph_dm1_frac_spatial_q75"]
    covar_sets = {
        "tss_only": ["tss"],
        "tss_sex": ["tss", "sex"],
        "tss_histology_molecular_sex": ["tss", "histology_subtype", "molecular_subtype", "sex"],
    }
    rows = []
    base = tcga[tcga["dm"].isin(["DM1", "DM2"])].copy()
    for score in score_cols:
        for label, covars in covar_sets.items():
            sub = base[base[score].notna()].copy()
            if len(sub) < 10:
                continue
            resid = residualize_score(sub, score, covars)
            y = (sub["dm"] == "DM1").astype(int).to_numpy()
            rows.append(
                {
                    "score": score,
                    "target": "DM1_vs_DM2",
                    "residualized_against": label,
                    "n": int(np.isfinite(resid).sum()),
                    "n_pos": int(y[np.isfinite(resid)].sum()),
                    "n_neg": int(np.isfinite(resid).sum() - y[np.isfinite(resid)].sum()),
                    "auc_pos_high_after_resid": auc_or_nan(y[np.isfinite(resid)], resid[np.isfinite(resid)]),
                    "mannwhitney_p_after_resid": mwu_test(
                        resid[(y == 1) & np.isfinite(resid)], resid[(y == 0) & np.isfinite(resid)]
                    ),
                    "median_pos_resid": float(np.nanmedian(resid[y == 1])),
                    "median_neg_resid": float(np.nanmedian(resid[y == 0])),
                }
            )
    for score in score_cols:
        sub = tcga[tcga[score].notna()].copy()
        if len(sub) < 10:
            continue
        for cov in ["tss", "histology_subtype", "molecular_subtype", "sex"]:
            if cov not in sub.columns:
                continue
            x = pd.get_dummies(sub[cov].fillna("missing").astype(str), drop_first=True)
            if x.shape[1] == 0:
                continue
            y = sub[score].astype(float).to_numpy()
            lr = LinearRegression().fit(x.to_numpy(dtype=float), y)
            r2 = float(lr.score(x.to_numpy(dtype=float), y))
            rows.append(
                {
                    "score": score,
                    "target": f"{score}_explained_by_{cov}",
                    "residualized_against": cov,
                    "n": int(len(sub)),
                    "n_pos": "",
                    "n_neg": "",
                    "auc_pos_high_after_resid": "",
                    "mannwhitney_p_after_resid": "",
                    "median_pos_resid": "",
                    "median_neg_resid": "",
                    "r2": r2,
                }
            )
    return pd.DataFrame(rows)


def plot_spatial(pred: pd.DataFrame, slide_summary: pd.DataFrame, out_dir: Path) -> None:
    stage_order = ["PT", "PTC", "LPTC", "ATC"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axes[0]
    ax.scatter(pred["target_DM1_like_z"], pred["pred_morph_DM1_z"], s=8, alpha=0.35, color="#2f6f8f")
    rho, p = safe_spearman(pred["target_DM1_like_z"], pred["pred_morph_DM1_z"])
    ax.set_title(f"Leave-slide-out spots: rho={rho:.2f}, p={p:.1e}")
    ax.set_xlabel("Spatial transcriptomic DM1_like score (z)")
    ax.set_ylabel("H&E UNI-predicted morphology score (z)")
    ax.axhline(0, color="#999", lw=0.8)
    ax.axvline(0, color="#999", lw=0.8)

    ax = axes[1]
    data = [slide_summary.loc[slide_summary["stage"] == s, "mean_pred_morph_DM1_z"].dropna() for s in stage_order]
    ax.boxplot(data, labels=stage_order, patch_artist=True)
    for i, vals in enumerate(data, start=1):
        x = np.full(len(vals), i) + np.linspace(-0.08, 0.08, len(vals))
        ax.scatter(x, vals, s=36, color="#7b1f2a", alpha=0.85, zorder=3)
    ax.set_title("Slide mean morphology-DM1 score by spatial stage")
    ax.set_ylabel("Mean held-out predicted score")
    ax.axhline(0, color="#999", lw=0.8)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_pathology_scout_spatial_loso.png", dpi=180)
    fig.savefig(out_dir / "fig_pathology_scout_spatial_loso.pdf")
    plt.close(fig)


def plot_tcga(tcga: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
    panels = [
        ("dm", ["DM1", "DM2", "not_DM"], "morph_dm1_mean", "Morphology-DM1 by TCGA DM group"),
        ("molecular_subtype", ["BRAF_like", "RAS_like", "unknown"], "morph_dm1_mean", "By molecular subtype"),
        ("tds_group", ["low", "mid", "high"], "morph_dm1_mean", "By TDS group"),
    ]
    for ax, (col, order, score, title) in zip(axes.flat[:3], panels):
        data = [tcga.loc[tcga[col] == g, score].dropna() for g in order]
        ax.boxplot(data, labels=order, patch_artist=True)
        for i, vals in enumerate(data, start=1):
            if len(vals) == 0:
                continue
            x = np.full(len(vals), i) + np.linspace(-0.12, 0.12, len(vals))
            ax.scatter(x, vals, s=26, alpha=0.8, color="#2f6f8f", zorder=3)
        ax.set_title(title)
        ax.set_ylabel(score)
        ax.axhline(0, color="#999", lw=0.8)
        ax.tick_params(axis="x", rotation=20)
    ax = axes.flat[3]
    sub = tcga[["morph_dm1_mean", "prob_DM1_uni_loto"]].dropna()
    if len(sub) >= 3:
        rho, p = safe_spearman(sub["morph_dm1_mean"], sub["prob_DM1_uni_loto"])
        ax.scatter(sub["morph_dm1_mean"], sub["prob_DM1_uni_loto"], s=35, alpha=0.8, color="#7b1f2a")
        ax.set_title(f"Versus UNI LOTO DM1 probability: rho={rho:.2f}, p={p:.1e}")
        ax.set_xlabel("Spatial-trained morphology-DM1 score")
        ax.set_ylabel("UNI LOTO prob_DM1")
    else:
        ax.text(0.5, 0.5, "UNI LOTO predictions unavailable", ha="center", va="center")
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_pathology_scout_tcga_projection.png", dpi=180)
    fig.savefig(out_dir / "fig_pathology_scout_tcga_projection.pdf")
    plt.close(fig)


def get_assoc_row(assoc: pd.DataFrame, score: str, test: str) -> dict:
    sub = assoc[(assoc["score"] == score) & (assoc["test"] == test)]
    if sub.empty:
        return {}
    return sub.iloc[0].to_dict()


def write_report(summary: dict, out_dir: Path) -> None:
    lines = []
    lines.append("# Pathology finding scout - spatial-DM1 morphology transfer")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append(f"- Verdict: **{summary['verdict']}**")
    lines.append(f"- Spatial leave-slide-out Spearman rho: {summary['spatial_loso_overall_rho']:.3f}")
    lines.append(f"- Spatial median per-slide Spearman rho: {summary['spatial_loso_median_slide_rho']:.3f}")
    lines.append(
        f"- TCGA DM1 vs DM2 AUC from spatial-trained morphology score: "
        f"{summary['tcga_dm1_vs_dm2_auc']:.3f}"
    )
    lines.append(
        f"- TSS-residualized TCGA DM1 vs DM2 AUC: "
        f"{summary['tcga_dm1_vs_dm2_auc_tss_resid']:.3f}"
    )
    lines.append(
        f"- Spearman versus UNI LOTO DM1 probability: "
        f"{summary['rho_morph_vs_uni_loto']:.3f}"
    )
    lines.append("")
    lines.append("## Candidate finding")
    lines.append("")
    if summary["verdict"] == "CANDIDATE_PASS":
        lines.append(
            "A spatial-transcriptomic DM1 morphology axis learned from GSE250521 H&E UNI "
            "tiles transfers to TCGA-THCA WSI embeddings and ranks TCGA DM1 above DM2. "
            "This is a candidate pathology finding, not yet a main-text claim."
        )
    elif summary["verdict"] == "SPATIAL_ONLY":
        lines.append(
            "H&E tiles carry within-spatial-cohort DM1 transcriptomic information, but "
            "the transfer to TCGA DM1 labels is too weak or confounded for a new finding."
        )
    else:
        lines.append(
            "No defensible new pathology finding was found in this scout. The result should "
            "be treated as a negative control unless later validation improves it."
        )
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append(
        "- Fit a PCA(32)+ridge model from GSE250521 UNI tile embeddings to spot-level "
        "DM1_like_score, evaluated by leave-one-slide-out prediction."
    )
    lines.append(
        "- Refit the same model on all spatial spots, projected it to TCGA UNI tile "
        "features, and summarized each WSI by mean, q90, and high-tile fractions."
    )
    lines.append(
        "- Tested TCGA associations with DM group, molecular subtype, histology, TDS, "
        "and conservative residualization against TSS/sex/histology/molecular subtype."
    )
    lines.append("")
    lines.append("## Key files")
    lines.append("")
    for name in [
        "PATHOLOGY_SCOUT_SUMMARY.json",
        "spatial_spot_loso_predictions.tsv",
        "spatial_slide_summary.tsv",
        "tcga_slide_pathology_scores.tsv",
        "tcga_association_tests.tsv",
        "tcga_confound_residual_tests.tsv",
        "fig_pathology_scout_spatial_loso.png",
        "fig_pathology_scout_tcga_projection.png",
    ]:
        lines.append(f"- `{name}`")
    lines.append("")
    (out_dir / "PATHOLOGY_SCOUT_REPORT.md").write_text("\n".join(lines).rstrip() + "\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    spatial_meta, x_spatial, y_spatial = load_spatial()
    spatial_pred = spatial_leave_one_slide_out(spatial_meta, x_spatial, y_spatial)
    slide_summary = summarize_spatial(spatial_pred)
    scaler, pca, ridge, spatial_pred_final = fit_final_spatial_model(x_spatial, y_spatial)
    tcga_meta = load_tcga_metadata()
    tcga_scores = project_tcga(tcga_meta, scaler, pca, ridge, spatial_pred_final)
    assoc = association_tests(tcga_scores)
    resid = confound_residual_tests(tcga_scores)

    spatial_pred.to_csv(OUT_DIR / "spatial_spot_loso_predictions.tsv", sep="\t", index=False, na_rep="NA")
    slide_summary.to_csv(OUT_DIR / "spatial_slide_summary.tsv", sep="\t", index=False, na_rep="NA")
    tcga_scores.to_csv(OUT_DIR / "tcga_slide_pathology_scores.tsv", sep="\t", index=False, na_rep="NA")
    assoc.to_csv(OUT_DIR / "tcga_association_tests.tsv", sep="\t", index=False, na_rep="NA")
    resid.to_csv(OUT_DIR / "tcga_confound_residual_tests.tsv", sep="\t", index=False, na_rep="NA")

    plot_spatial(spatial_pred, slide_summary, OUT_DIR)
    plot_tcga(tcga_scores, OUT_DIR)

    spatial_rho, spatial_p = safe_spearman(
        spatial_pred["target_DM1_like_z"], spatial_pred["pred_morph_DM1_z"]
    )
    median_slide_rho = float(np.nanmedian(slide_summary["rho_spot_target_vs_pred"]))
    dm_row = get_assoc_row(assoc, "morph_dm1_mean", "DM1_vs_DM2")
    dm_auc = float(dm_row.get("auc_pos_high", float("nan")))
    dm_p = float(dm_row.get("mannwhitney_p", float("nan")))
    resid_row = resid[
        (resid["score"] == "morph_dm1_mean")
        & (resid["target"] == "DM1_vs_DM2")
        & (resid["residualized_against"] == "tss_only")
    ]
    resid_auc = float(resid_row["auc_pos_high_after_resid"].iloc[0]) if not resid_row.empty else float("nan")
    morph_uni = assoc[
        (assoc["score"] == "morph_dm1_mean") & (assoc["test"] == "spearman_vs_prob_DM1_uni_loto")
    ]
    rho_uni = float(morph_uni["spearman_rho"].iloc[0]) if not morph_uni.empty else float("nan")
    p_uni = float(morph_uni["spearman_p"].iloc[0]) if not morph_uni.empty else float("nan")
    morph_tds = assoc[(assoc["score"] == "morph_dm1_mean") & (assoc["test"] == "spearman_vs_tds_score")]
    rho_tds = float(morph_tds["spearman_rho"].iloc[0]) if not morph_tds.empty else float("nan")
    p_tds = float(morph_tds["spearman_p"].iloc[0]) if not morph_tds.empty else float("nan")

    if spatial_rho >= 0.20 and dm_auc >= 0.65 and resid_auc >= 0.60:
        verdict = "CANDIDATE_PASS"
    elif spatial_rho >= 0.20:
        verdict = "SPATIAL_ONLY"
    else:
        verdict = "NO_GO"

    summary = {
        "verdict": verdict,
        "n_spatial_spots": int(len(spatial_pred)),
        "n_spatial_slides": int(spatial_pred["slide"].nunique()),
        "spatial_loso_overall_rho": spatial_rho,
        "spatial_loso_overall_p": spatial_p,
        "spatial_loso_median_slide_rho": median_slide_rho,
        "n_tcga_slides_scored": int(len(tcga_scores)),
        "tcga_dm1_vs_dm2_auc": dm_auc,
        "tcga_dm1_vs_dm2_mannwhitney_p": dm_p,
        "tcga_dm1_vs_dm2_auc_tss_resid": resid_auc,
        "rho_morph_vs_uni_loto": rho_uni,
        "p_morph_vs_uni_loto": p_uni,
        "rho_morph_vs_tds_score": rho_tds,
        "p_morph_vs_tds_score": p_tds,
        "model": {"pca_components": N_PCS, "ridge_alpha": RIDGE_ALPHA},
    }
    (OUT_DIR / "PATHOLOGY_SCOUT_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_report(summary, OUT_DIR)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
