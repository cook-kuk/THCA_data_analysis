#!/usr/bin/env python3
"""GSE250521 stage-stratified and leave-one-stage-out controls.

This control asks whether the Path2Space-inspired UNI signal is only a stage
surrogate. It has two pieces:

1. Re-score the existing leave-one-slide-out predictions after centering within
   stage and summarize per-stage slide-level performance.
2. Train leave-one-stage-out models where all slides from one stage are held
   out together.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.preprocessing import StandardScaler


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
P2_ROOT = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
BASE = P2_ROOT / "analysis_supp/path2space_inspired_reanalysis_2026_05_09"
PHASE1 = P2_ROOT / "phase1_gse250521"
OUT = P2_ROOT / "analysis_supp/path2space_stage_generalization_controls_2026_05_09"

TARGETS = ["DM1_like_score", "TDS_like_score", "MAPK_output_score"]
STAGE_ORDER = ["PT", "PTC", "LPTC", "ATC"]
N_PCS = 64
RIDGE_ALPHA = 100.0
RANDOM_SEED = 9259


def safe_spearman(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 5 or np.nanstd(x[keep]) == 0 or np.nanstd(y[keep]) == 0:
        return np.nan, np.nan, int(keep.sum())
    rho, p = stats.spearmanr(x[keep], y[keep])
    return float(rho), float(p), int(keep.sum())


def center_by_group(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups).astype(str)
    out = np.full(len(values), np.nan)
    for group in pd.unique(groups):
        idx = groups == group
        out[idx] = values[idx] - np.nanmean(values[idx])
    return out


def residualize_within_slide(values: np.ndarray, covars: np.ndarray, slides: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    covars = np.asarray(covars, dtype=float)
    slides = np.asarray(slides).astype(str)
    out = np.full(len(values), np.nan)
    for slide in pd.unique(slides):
        idx = np.where(slides == slide)[0]
        y = values[idx]
        x = covars[idx]
        keep = np.isfinite(y) & np.isfinite(x).all(axis=1)
        if keep.sum() < 30:
            continue
        xs = StandardScaler().fit_transform(x[keep])
        fit = LinearRegression().fit(xs, y[keep])
        out[idx[keep]] = y[keep] - fit.predict(xs)
    return out


def load_inputs() -> tuple[pd.DataFrame, np.ndarray]:
    pred = pd.read_csv(BASE / "path2space_spot_predictions.tsv.gz", sep="\t")
    input_df = pd.read_csv(BASE / "path2space_input_spots.tsv.gz", sep="\t")
    emb = np.load(PHASE1 / "uni_embeddings_size224.npz")["embeddings"].astype(np.float32)
    input_df = input_df.sort_values("embed_idx").reset_index(drop=True)
    key_cols = ["sample_id", "spot_id", "stage", "array_row", "array_col", "total_counts", "n_genes_by_counts", "pct_counts_mt"]
    meta = input_df[key_cols].copy()
    df = meta.merge(pred, on=["sample_id", "spot_id", "stage", "array_row", "array_col"], how="inner")
    if len(df) != emb.shape[0]:
        raise RuntimeError(f"Expected {emb.shape[0]} embedded spots, merged {len(df)}")
    return df, emb


def loso_stage_predict(x: np.ndarray, y: np.ndarray, stages: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    stages = np.asarray(stages).astype(str)
    pred = np.full(len(y), np.nan)
    finite_x = np.isfinite(x).all(axis=1)
    for held in STAGE_ORDER:
        train = (stages != held) & finite_x & np.isfinite(y)
        test = (stages == held) & finite_x
        if train.sum() < 50 or test.sum() == 0:
            continue
        scaler = StandardScaler()
        x_train = scaler.fit_transform(x[train])
        x_test = scaler.transform(x[test])
        n_comp = min(N_PCS, x_train.shape[0] - 1, x_train.shape[1])
        pca = PCA(n_components=n_comp, random_state=RANDOM_SEED)
        x_train = pca.fit_transform(x_train)
        x_test = pca.transform(x_test)
        mean = float(np.nanmean(y[train]))
        sd = float(np.nanstd(y[train]))
        if not np.isfinite(sd) or sd == 0:
            sd = 1.0
        model = Ridge(alpha=RIDGE_ALPHA)
        model.fit(x_train, (y[train] - mean) / sd)
        pred[test] = model.predict(x_test)
    return pred


def summarize_existing_loso(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    slides = df["sample_id"].astype(str).to_numpy()
    stages = df["stage"].astype(str).to_numpy()
    rows = []
    per_stage_rows = []
    for target in TARGETS:
        obs = df[f"obs_{target}_smooth8"].to_numpy(dtype=float)
        pred = df[f"pred_{target}_smooth8"].to_numpy(dtype=float)
        pooled, pp, pn = safe_spearman(obs, pred)
        slide_centered, sp, sn = safe_spearman(center_by_group(obs, slides), center_by_group(pred, slides))
        stage_centered, gp, gn = safe_spearman(center_by_group(obs, stages), center_by_group(pred, stages))
        rows.append(
            {
                "target": target,
                "evaluation": "existing_loso",
                "pooled_rho": pooled,
                "pooled_p": pp,
                "slide_centered_rho": slide_centered,
                "slide_centered_p": sp,
                "stage_centered_rho": stage_centered,
                "stage_centered_p": gp,
                "n": pn,
            }
        )
        for stage in STAGE_ORDER:
            idx = np.where(stages == stage)[0]
            rho, p, n = safe_spearman(obs[idx], pred[idx])
            sc_rho, sc_p, _ = safe_spearman(center_by_group(obs[idx], slides[idx]), center_by_group(pred[idx], slides[idx]))
            per_slide = []
            for slide in pd.unique(slides[idx]):
                sidx = idx[slides[idx] == slide]
                srho, _, _ = safe_spearman(obs[sidx], pred[sidx])
                per_slide.append(srho)
            per_stage_rows.append(
                {
                    "target": target,
                    "evaluation": "existing_loso",
                    "stage": stage,
                    "n_spots": n,
                    "pooled_within_stage_rho": rho,
                    "pooled_within_stage_p": p,
                    "slide_centered_within_stage_rho": sc_rho,
                    "slide_centered_within_stage_p": sc_p,
                    "median_slide_rho": float(np.nanmedian(per_slide)),
                    "slides_rho_gt_0_2": int(np.nansum(np.asarray(per_slide) > 0.2)),
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(per_stage_rows)


def summarize_leave_stage_out(df: pd.DataFrame, emb: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    stages = df["stage"].astype(str).to_numpy()
    slides = df["sample_id"].astype(str).to_numpy()
    coord_qc = np.column_stack(
        [
            df[["array_row", "array_col"]].to_numpy(dtype=float),
            np.log1p(df["total_counts"].to_numpy(dtype=float)),
            np.log1p(df["n_genes_by_counts"].to_numpy(dtype=float)),
            df["pct_counts_mt"].to_numpy(dtype=float),
        ]
    )

    pred_out = df[["sample_id", "spot_id", "stage", "array_row", "array_col"]].copy()
    rows = []
    per_stage_rows = []
    for target in TARGETS:
        raw = df[f"obs_{target}_smooth8"].to_numpy(dtype=float)
        resid = residualize_within_slide(raw, coord_qc, slides)
        for mode, obs in [("raw_smoothed", raw), ("coord_qc_residual_target", resid)]:
            pred = loso_stage_predict(emb, obs, stages)
            pred_out[f"obs_{target}_{mode}"] = obs
            pred_out[f"pred_{target}_{mode}"] = pred
            pooled, pp, pn = safe_spearman(obs, pred)
            stage_centered, gp, gn = safe_spearman(center_by_group(obs, stages), center_by_group(pred, stages))
            rows.append(
                {
                    "target": target,
                    "mode": mode,
                    "evaluation": "leave_one_stage_out",
                    "pooled_rho": pooled,
                    "pooled_p": pp,
                    "stage_centered_rho": stage_centered,
                    "stage_centered_p": gp,
                    "n": pn,
                }
            )
            for stage in STAGE_ORDER:
                idx = np.where(stages == stage)[0]
                rho, p, n = safe_spearman(obs[idx], pred[idx])
                sc_rho, sc_p, _ = safe_spearman(center_by_group(obs[idx], slides[idx]), center_by_group(pred[idx], slides[idx]))
                per_stage_rows.append(
                    {
                        "target": target,
                        "mode": mode,
                        "evaluation": "leave_one_stage_out",
                        "held_out_stage": stage,
                        "n_spots": n,
                        "rho": rho,
                        "p": p,
                        "slide_centered_rho": sc_rho,
                        "slide_centered_p": sc_p,
                    }
                )
    return pd.DataFrame(rows), pd.DataFrame(per_stage_rows), pred_out


def write_figure(existing_stage: pd.DataFrame, lostage_stage: pd.DataFrame, summary: dict) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.2))

    dm1_existing = existing_stage[existing_stage["target"] == "DM1_like_score"].set_index("stage").reindex(STAGE_ORDER)
    x = np.arange(len(STAGE_ORDER))
    axes[0, 0].bar(x, dm1_existing["pooled_within_stage_rho"], color="#2b6f73", label="pooled")
    axes[0, 0].scatter(x, dm1_existing["median_slide_rho"], color="#d6b25e", s=55, label="median slide")
    axes[0, 0].axhline(0, color="#333", lw=0.8)
    axes[0, 0].set_xticks(x)
    axes[0, 0].set_xticklabels(STAGE_ORDER)
    axes[0, 0].set_ylim(-0.15, 0.75)
    axes[0, 0].set_ylabel("Spearman rho")
    axes[0, 0].set_title("A. Existing LOSO DM1/RAI by stage")
    axes[0, 0].legend(frameon=False)

    heat = existing_stage.pivot(index="target", columns="stage", values="slide_centered_within_stage_rho").reindex(TARGETS)[STAGE_ORDER]
    im = axes[0, 1].imshow(heat.to_numpy(dtype=float), vmin=-0.2, vmax=0.7, cmap="RdBu_r", aspect="auto")
    axes[0, 1].set_xticks(np.arange(len(STAGE_ORDER)))
    axes[0, 1].set_xticklabels(STAGE_ORDER)
    axes[0, 1].set_yticks(np.arange(len(TARGETS)))
    axes[0, 1].set_yticklabels([t.replace("_score", "") for t in TARGETS])
    axes[0, 1].set_title("B. Existing LOSO slide-centered rho")
    for i in range(len(TARGETS)):
        for j in range(len(STAGE_ORDER)):
            val = heat.iloc[i, j]
            axes[0, 1].text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=10)
    fig.colorbar(im, ax=axes[0, 1], fraction=0.046, pad=0.04)

    dm1_lostage = lostage_stage[
        (lostage_stage["target"] == "DM1_like_score") & (lostage_stage["mode"] == "raw_smoothed")
    ].set_index("held_out_stage").reindex(STAGE_ORDER)
    dm1_lostage_res = lostage_stage[
        (lostage_stage["target"] == "DM1_like_score") & (lostage_stage["mode"] == "coord_qc_residual_target")
    ].set_index("held_out_stage").reindex(STAGE_ORDER)
    axes[1, 0].bar(x - 0.18, dm1_lostage["rho"], width=0.35, color="#2b6f73", label="raw")
    axes[1, 0].bar(x + 0.18, dm1_lostage_res["rho"], width=0.35, color="#d6b25e", label="coord+QC residual")
    axes[1, 0].axhline(0, color="#333", lw=0.8)
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels(STAGE_ORDER)
    axes[1, 0].set_ylim(-0.25, 0.58)
    axes[1, 0].set_ylabel("Held-stage rho")
    axes[1, 0].set_title("C. Leave-one-stage-out DM1/RAI")
    axes[1, 0].legend(frameon=False)

    heat2 = (
        lostage_stage[lostage_stage["mode"] == "raw_smoothed"]
        .pivot(index="target", columns="held_out_stage", values="rho")
        .reindex(TARGETS)[STAGE_ORDER]
    )
    im2 = axes[1, 1].imshow(heat2.to_numpy(dtype=float), vmin=-0.3, vmax=0.45, cmap="RdBu_r", aspect="auto")
    axes[1, 1].set_xticks(np.arange(len(STAGE_ORDER)))
    axes[1, 1].set_xticklabels(STAGE_ORDER)
    axes[1, 1].set_yticks(np.arange(len(TARGETS)))
    axes[1, 1].set_yticklabels([t.replace("_score", "") for t in TARGETS])
    axes[1, 1].set_title("D. Leave-one-stage-out raw rho")
    for i in range(len(TARGETS)):
        for j in range(len(STAGE_ORDER)):
            val = heat2.iloc[i, j]
            axes[1, 1].text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=10)
    fig.colorbar(im2, ax=axes[1, 1], fraction=0.046, pad=0.04)

    fig.suptitle("GSE250521 stage-stratified controls for UNI-to-DM1/RAI", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "fig_gse250521_stage_generalization_controls.png", dpi=220)
    fig.savefig(OUT / "fig_gse250521_stage_generalization_controls.pdf")
    plt.close(fig)


def write_report(existing_summary: pd.DataFrame, existing_stage: pd.DataFrame, lostage_summary: pd.DataFrame, lostage_stage: pd.DataFrame, summary: dict) -> None:
    lines = [
        "# GSE250521 stage-generalization controls",
        "",
        "## Verdict",
        "",
        f"- Existing leave-one-slide-out DM1/RAI stage-centered rho: **{summary['dm1_existing_stage_centered_rho']:.3f}**.",
        f"- Existing leave-one-slide-out DM1/RAI by-stage pooled rho: PT **{summary['dm1_existing_pt_rho']:.3f}**, PTC **{summary['dm1_existing_ptc_rho']:.3f}**, LPTC **{summary['dm1_existing_lptc_rho']:.3f}**, ATC **{summary['dm1_existing_atc_rho']:.3f}**.",
        f"- Leave-one-stage-out raw DM1/RAI pooled rho: **{summary['dm1_lostage_raw_pooled_rho']:.3f}**; stage-centered rho **{summary['dm1_lostage_raw_stage_centered_rho']:.3f}**.",
        f"- Leave-one-stage-out coord+QC residual-target DM1/RAI pooled rho: **{summary['dm1_lostage_resid_pooled_rho']:.3f}**; stage-centered rho **{summary['dm1_lostage_resid_stage_centered_rho']:.3f}**.",
        "",
        "## Interpretation",
        "",
        "The existing leave-one-slide-out signal is not only a global stage mean because stage-centered evaluation remains positive. The raw signal also transfers in leave-one-stage-out evaluation, especially into PTC and LPTC. The caveat is the stricter coord+QC residual target: it stays positive but is weaker, and ATC is the least stable stage.",
        "",
        "Safe wording: use GSE250521 as a local image-to-DM1/RAI spatial result that is not just stage mean structure, with stage heterogeneity and weaker residual-stage transfer disclosed.",
        "",
        "## Existing Leave-One-Slide-Out Summary",
        "",
        existing_summary.round(4).to_markdown(index=False),
        "",
        "## Existing Leave-One-Slide-Out By Stage",
        "",
        existing_stage.round(4).to_markdown(index=False),
        "",
        "## Leave-One-Stage-Out Summary",
        "",
        lostage_summary.round(4).to_markdown(index=False),
        "",
        "## Leave-One-Stage-Out By Held Stage",
        "",
        lostage_stage.round(4).to_markdown(index=False),
    ]
    (OUT / "GSE250521_STAGE_GENERALIZATION_REPORT.md").write_text("\n".join(lines) + "\n")

    summary_lines = [
        "# GSE250521 Stage-Generalization Controls",
        "",
        "## Verdict",
        "",
        "This is a caveat-bearing control.",
        "",
        "| Layer | Result | Interpretation |",
        "|---|---:|---|",
        f"| Existing LOSO DM1/RAI stage-centered rho | {summary['dm1_existing_stage_centered_rho']:.3f} | positive after stage centering |",
        f"| Existing LOSO DM1/RAI PT/PTC/LPTC/ATC rho | {summary['dm1_existing_pt_rho']:.3f} / {summary['dm1_existing_ptc_rho']:.3f} / {summary['dm1_existing_lptc_rho']:.3f} / {summary['dm1_existing_atc_rho']:.3f} | ATC is the weak stage |",
        f"| Leave-one-stage-out raw DM1/RAI pooled rho | {summary['dm1_lostage_raw_pooled_rho']:.3f} | raw cross-stage transfer support |",
        f"| Leave-one-stage-out residual DM1/RAI pooled rho | {summary['dm1_lostage_resid_pooled_rho']:.3f} | residual transfer is weaker caveat |",
        "",
        "Safe wording: stage-centered local signal remains and raw leave-one-stage-out transfer is positive, but coord+QC residual transfer and ATC stability are weaker.",
    ]
    (OUT / "SUMMARY.md").write_text("\n".join(summary_lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df, emb = load_inputs()
    existing_summary, existing_stage = summarize_existing_loso(df)
    lostage_summary, lostage_stage, pred_out = summarize_leave_stage_out(df, emb)

    existing_summary.to_csv(OUT / "gse250521_existing_loso_stage_centering.tsv", sep="\t", index=False, na_rep="NA")
    existing_stage.to_csv(OUT / "gse250521_existing_loso_by_stage.tsv", sep="\t", index=False, na_rep="NA")
    lostage_summary.to_csv(OUT / "gse250521_leave_one_stage_out_summary.tsv", sep="\t", index=False, na_rep="NA")
    lostage_stage.to_csv(OUT / "gse250521_leave_one_stage_out_by_stage.tsv", sep="\t", index=False, na_rep="NA")
    pred_out.to_csv(OUT / "gse250521_leave_one_stage_out_predictions.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")

    dm1_existing = existing_summary[existing_summary["target"] == "DM1_like_score"].iloc[0]
    dm1_stage = existing_stage[existing_stage["target"] == "DM1_like_score"].set_index("stage")
    dm1_lostage_raw = lostage_summary[
        (lostage_summary["target"] == "DM1_like_score") & (lostage_summary["mode"] == "raw_smoothed")
    ].iloc[0]
    dm1_lostage_resid = lostage_summary[
        (lostage_summary["target"] == "DM1_like_score") & (lostage_summary["mode"] == "coord_qc_residual_target")
    ].iloc[0]
    summary = {
        "n_spots": int(len(df)),
        "n_slides": int(df["sample_id"].nunique()),
        "dm1_existing_pooled_rho": float(dm1_existing["pooled_rho"]),
        "dm1_existing_slide_centered_rho": float(dm1_existing["slide_centered_rho"]),
        "dm1_existing_stage_centered_rho": float(dm1_existing["stage_centered_rho"]),
        "dm1_existing_pt_rho": float(dm1_stage.loc["PT", "pooled_within_stage_rho"]),
        "dm1_existing_ptc_rho": float(dm1_stage.loc["PTC", "pooled_within_stage_rho"]),
        "dm1_existing_lptc_rho": float(dm1_stage.loc["LPTC", "pooled_within_stage_rho"]),
        "dm1_existing_atc_rho": float(dm1_stage.loc["ATC", "pooled_within_stage_rho"]),
        "dm1_lostage_raw_pooled_rho": float(dm1_lostage_raw["pooled_rho"]),
        "dm1_lostage_raw_stage_centered_rho": float(dm1_lostage_raw["stage_centered_rho"]),
        "dm1_lostage_resid_pooled_rho": float(dm1_lostage_resid["pooled_rho"]),
        "dm1_lostage_resid_stage_centered_rho": float(dm1_lostage_resid["stage_centered_rho"]),
    }
    (OUT / "GSE250521_STAGE_GENERALIZATION_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_figure(existing_stage, lostage_stage, summary)
    write_report(existing_summary, existing_stage, lostage_summary, lostage_stage, summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
