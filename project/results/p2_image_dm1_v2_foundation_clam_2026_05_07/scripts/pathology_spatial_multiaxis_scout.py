"""
Spatial pathology multi-axis scout.

Uses GSE250521 spatial H&E UNI tile embeddings and spot-level RNA program
scores to find which pathology programs are most visibly encoded in local
H&E morphology under leave-one-slide-out validation.

This is a discovery scan. It does not claim external generalization.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler


REPO = Path(os.environ.get("THCA_REPO_ROOT", Path.cwd()))
P2_ROOT = Path(
    os.environ.get(
        "THCA_P2_ROOT",
        REPO / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07",
    )
)
SPATIAL_DIR = P2_ROOT / "phase1_gse250521"
SCORES = REPO / "project/results/01_spatial_score/all_spots_scored.tsv.gz"
OUT_DIR = P2_ROOT / "analysis_supp/pathology_spatial_multiaxis_2026_05_09"

TARGETS = [
    "DM1_like_score",
    "RAI_8_score",
    "TDS_like_score",
    "CAF_ECM_score",
    "EMT_score",
    "Hypoxia_score",
    "Proliferation_score",
    "Epithelial_score",
]
N_PCS = 32
RIDGE_ALPHA = 100.0


def safe_spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4 or np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return float("nan"), float("nan")
    rho, p = spearmanr(x[m], y[m])
    return float(rho), float(p)


def load_data():
    emb = np.load(SPATIAL_DIR / "uni_embeddings_size224.npz")["embeddings"].astype(np.float32)
    meta = pd.read_csv(SPATIAL_DIR / "uni_embed_metadata_size224.tsv", sep="\t")
    scores = pd.read_csv(SCORES, sep="\t")
    meta = meta.reset_index().rename(columns={"index": "embed_idx"})
    meta["key"] = meta["slide"].astype(str) + "_" + meta["spot_id"].astype(str)
    scores["key"] = scores["sample_id"].astype(str) + "_" + scores["spot_id"].astype(str)
    cols = ["key", "sample_id", "stage", "in_tissue"] + TARGETS
    df = meta.merge(scores[cols], on="key", how="inner")
    df = df[(df["in_tissue"] == 1)].copy()
    x = emb[df["embed_idx"].to_numpy()]
    return df.reset_index(drop=True), x


def predict_loso(x, y, slides):
    out = np.full(len(y), np.nan)
    for held in sorted(pd.unique(slides)):
        train = slides != held
        test = ~train
        yt = y[train]
        yt_mean = float(np.nanmean(yt))
        yt_sd = float(np.nanstd(yt))
        if yt_sd == 0:
            yt_sd = 1.0
        scaler = StandardScaler()
        xtr = scaler.fit_transform(x[train])
        xte = scaler.transform(x[test])
        n_comp = min(N_PCS, xtr.shape[0] - 1, xtr.shape[1])
        pca = PCA(n_components=n_comp, random_state=7)
        xtr_pc = pca.fit_transform(xtr)
        xte_pc = pca.transform(xte)
        ridge = Ridge(alpha=RIDGE_ALPHA)
        ridge.fit(xtr_pc, (yt - yt_mean) / yt_sd)
        out[test] = ridge.predict(xte_pc)
    return out


def run_scan(df, x):
    pred_rows = []
    summary_rows = []
    slide_rows = []
    slides = df["slide"].to_numpy()
    for target in TARGETS:
        y_raw = pd.to_numeric(df[target], errors="coerce").to_numpy(dtype=float)
        mask = np.isfinite(y_raw)
        y = np.full(len(y_raw), np.nan)
        y[mask] = (y_raw[mask] - np.nanmean(y_raw[mask])) / (np.nanstd(y_raw[mask]) or 1.0)
        pred = predict_loso(x[mask], y[mask], slides[mask])
        full_pred = np.full(len(y), np.nan)
        full_pred[mask] = pred
        rho, p = safe_spearman(y, full_pred)
        per_slide_rhos = []
        for slide, idx in df.groupby("slide").groups.items():
            idx = np.asarray(list(idx))
            srho, sp = safe_spearman(y[idx], full_pred[idx])
            per_slide_rhos.append(srho)
            slide_rows.append(
                {
                    "target": target,
                    "slide": slide,
                    "stage": df.loc[idx[0], "stage"],
                    "n_spots": int(np.isfinite(y[idx]).sum()),
                    "rho": srho,
                    "p": sp,
                    "mean_target_z": float(np.nanmean(y[idx])),
                    "mean_pred_z": float(np.nanmean(full_pred[idx])),
                }
            )
        summary_rows.append(
            {
                "target": target,
                "n_spots": int(mask.sum()),
                "overall_spearman_rho": rho,
                "overall_spearman_p": p,
                "median_slide_rho": float(np.nanmedian(per_slide_rhos)),
                "n_slides_rho_gt_0_2": int(np.nansum(np.asarray(per_slide_rhos) > 0.2)),
                "n_slides_rho_gt_0_3": int(np.nansum(np.asarray(per_slide_rhos) > 0.3)),
            }
        )
        tmp = df[["slide", "spot_id", "stage"]].copy()
        tmp["target"] = target
        tmp["target_z"] = y
        tmp["pred_z"] = full_pred
        pred_rows.append(tmp)
    return (
        pd.DataFrame(summary_rows).sort_values("overall_spearman_rho", ascending=False),
        pd.DataFrame(slide_rows),
        pd.concat(pred_rows, ignore_index=True),
    )


def plot_summary(summary, slide_summary, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ordered = summary.sort_values("overall_spearman_rho")
    axes[0].barh(ordered["target"], ordered["overall_spearman_rho"], color="#2f6f8f")
    axes[0].axvline(0, color="#999", lw=0.8)
    axes[0].set_xlabel("Overall leave-slide-out Spearman rho")
    axes[0].set_title("Which RNA programs are visible in H&E UNI morphology?")

    top_targets = summary.head(4)["target"].tolist()
    data = [
        slide_summary.loc[slide_summary["target"] == t, "rho"].dropna().to_numpy()
        for t in top_targets
    ]
    axes[1].boxplot(data, labels=top_targets, patch_artist=True)
    for i, vals in enumerate(data, 1):
        if len(vals) == 0:
            continue
        axes[1].scatter(np.full(len(vals), i) + np.linspace(-0.1, 0.1, len(vals)), vals, s=28, color="#7b1f2a")
    axes[1].axhline(0, color="#999", lw=0.8)
    axes[1].set_ylabel("Per-slide Spearman rho")
    axes[1].set_title("Slide-level robustness for top programs")
    axes[1].tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_spatial_multiaxis_summary.png", dpi=180)
    fig.savefig(out_dir / "fig_spatial_multiaxis_summary.pdf")
    plt.close(fig)


def write_report(summary, out_dir):
    top = summary.iloc[0]
    verdict = "CANDIDATE_PASS" if top.overall_spearman_rho >= 0.35 and top.n_slides_rho_gt_0_2 >= 10 else "WEAK"
    lines = [
        "# Spatial pathology multi-axis scout",
        "",
        "## Verdict",
        "",
        f"- Verdict: **{verdict}**",
        f"- Top morphology-visible RNA program: **{top.target}**",
        f"- Overall leave-slide-out Spearman rho: {top.overall_spearman_rho:.3f}",
        f"- Slides with rho > 0.2: {int(top.n_slides_rho_gt_0_2)}/16",
        "",
        "## Ranked programs",
        "",
        "| rank | target | overall rho | p | median slide rho | slides rho>0.2 |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for i, row in enumerate(summary.itertuples(index=False), 1):
        lines.append(
            f"| {i} | {row.target} | {row.overall_spearman_rho:.3f} | "
            f"{row.overall_spearman_p:.2e} | {row.median_slide_rho:.3f} | "
            f"{int(row.n_slides_rho_gt_0_2)}/16 |"
        )
    lines += [
        "",
        "## Interpretation boundary",
        "",
        "This is a spatial-cohort pathology finding scout. It shows what local H&E "
        "morphology can reconstruct from RNA spots under leave-one-slide-out "
        "validation. It is not yet external WSI generalization.",
    ]
    (out_dir / "SPATIAL_MULTIAXIS_REPORT.md").write_text("\n".join(lines) + "\n")
    return verdict


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df, x = load_data()
    summary, slide_summary, preds = run_scan(df, x)
    summary.to_csv(OUT_DIR / "spatial_multiaxis_summary.tsv", sep="\t", index=False, na_rep="NA")
    slide_summary.to_csv(OUT_DIR / "spatial_multiaxis_slide_summary.tsv", sep="\t", index=False, na_rep="NA")
    preds.to_csv(OUT_DIR / "spatial_multiaxis_spot_predictions.tsv", sep="\t", index=False, na_rep="NA")
    plot_summary(summary, slide_summary, OUT_DIR)
    verdict = write_report(summary, OUT_DIR)
    out = {
        "verdict": verdict,
        "top_target": str(summary.iloc[0]["target"]),
        "top_overall_rho": float(summary.iloc[0]["overall_spearman_rho"]),
        "top_p": float(summary.iloc[0]["overall_spearman_p"]),
        "top_slides_rho_gt_0_2": int(summary.iloc[0]["n_slides_rho_gt_0_2"]),
        "n_targets": int(len(summary)),
    }
    (OUT_DIR / "SPATIAL_MULTIAXIS_SUMMARY.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
