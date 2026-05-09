#!/usr/bin/env python3
"""v16: diagnostic stress tests for v15 spatial and PRISM extensions.

Purpose:
  - determine whether any reasonable GSE250521 spatial stratum/adjustment
    recovers the expected MAPK x Panel anti-correlation;
  - quantify how strong and how narrow the PRISM MAPK-axis enrichment is.

This is reviewer-surface diagnostics, not protected manuscript prose.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/p_deconv_2026_05_08"
SPOT = OUT / "v15B_spatial_mapk_panel_per_spot.tsv.gz"
SPATIAL_META = ROOT / "project/results/01_spatial_score/all_spots_scored.tsv.gz"
PRISM = OUT / "v15C_prism_annotated_drugs.tsv"
PRISM_CELLS = OUT / "v15C_prism_mapk_cellline_scores.tsv.gz"

TUMOR_STAGES = {"PTC", "LPTC", "ATC"}


def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    sp = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / sp) if sp > 0 else float("nan")


def spearman(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    df = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    if len(df) < 5 or df["x"].nunique() < 2 or df["y"].nunique() < 2:
        return (np.nan, np.nan, int(len(df)))
    r, p = stats.spearmanr(df["x"], df["y"])
    return (float(r), float(p), int(len(df)))


def partial_spearman(x: pd.Series, y: pd.Series, cov: pd.DataFrame) -> tuple[float, float, int]:
    df = pd.concat([x.rename("x"), y.rename("y"), cov], axis=1).dropna()
    cov_cols = list(cov.columns)
    cov_cols = [c for c in cov_cols if df[c].nunique(dropna=True) > 1]
    if len(df) < 20 or df["x"].nunique() < 2 or df["y"].nunique() < 2:
        return (np.nan, np.nan, int(len(df)))
    if not cov_cols:
        return spearman(df["x"], df["y"])
    rx = stats.rankdata(df["x"].to_numpy(float))
    ry = stats.rankdata(df["y"].to_numpy(float))
    design = [np.ones(len(df))]
    for c in cov_cols:
        design.append(stats.rankdata(df[c].to_numpy(float)))
    c_rank = np.column_stack(design)
    bx, *_ = np.linalg.lstsq(c_rank, rx, rcond=None)
    by, *_ = np.linalg.lstsq(c_rank, ry, rcond=None)
    ex = rx - c_rank @ bx
    ey = ry - c_rank @ by
    r, p = stats.pearsonr(ex, ey)
    return (float(r), float(p), int(len(df)))


def fisher_pool(rows: pd.DataFrame, rho_col: str, n_col: str = "n") -> dict[str, float]:
    sub = rows.dropna(subset=[rho_col, n_col]).copy()
    sub = sub[(sub[n_col] > 3) & (sub[rho_col].abs() < 1)]
    if sub.empty:
        return {"rho": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "n": 0, "k": 0}
    z = np.arctanh(sub[rho_col].to_numpy(float))
    w = sub[n_col].to_numpy(float) - 3
    mz = float((z * w).sum() / w.sum())
    se = 1 / math.sqrt(float(w.sum()))
    return {
        "rho": float(np.tanh(mz)),
        "ci_lo": float(np.tanh(mz - 1.96 * se)),
        "ci_hi": float(np.tanh(mz + 1.96 * se)),
        "n": int(sub[n_col].sum()),
        "k": int(len(sub)),
    }


def hypergeom_p(n_total: int, n_success: int, n_draw: int, n_success_draw: int) -> float:
    if n_total <= 0 or n_success <= 0 or n_draw <= 0:
        return float("nan")
    return float(stats.hypergeom.sf(n_success_draw - 1, n_total, n_success, n_draw))


def load_spatial() -> pd.DataFrame:
    spot = pd.read_csv(SPOT, sep="\t")
    meta_cols = [
        "sample_id", "spot_id", "CAF_ECM_score", "EMT_score", "Hypoxia_score",
        "Proliferation_score", "Epithelial_score", "pxl_row_in_fullres", "pxl_col_in_fullres",
    ]
    meta = pd.read_csv(SPATIAL_META, sep="\t", usecols=meta_cols)
    df = spot.merge(meta, on=["sample_id", "spot_id"], how="left")
    df["is_tumor_stage"] = df["stage"].isin(TUMOR_STAGES)
    df["high_epithelial_median"] = df["Epithelial_score"] >= df.groupby("sample_id")["Epithelial_score"].transform("median")
    df["high_epithelial_q75"] = df["Epithelial_score"] >= df.groupby("sample_id")["Epithelial_score"].transform(lambda s: s.quantile(0.75))
    df["low_epithelial_q25"] = df["Epithelial_score"] <= df.groupby("sample_id")["Epithelial_score"].transform(lambda s: s.quantile(0.25))
    df["epithelial_quartile"] = (
        df.groupby("sample_id")["Epithelial_score"]
        .transform(lambda s: pd.qcut(s.rank(method="first"), 4, labels=["Q1_low", "Q2", "Q3", "Q4_high"]))
        .astype(str)
    )
    return df


def run_spatial_diagnostics(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    print("\n[v16A] spatial diagnostics")
    cov_models = {
        "raw": [],
        "QC": ["log_total_counts", "log_n_genes", "pct_counts_mt"],
        "epithelial": ["Epithelial_score"],
        "cell_state": ["Epithelial_score", "CAF_ECM_score", "EMT_score", "Hypoxia_score", "Proliferation_score"],
        "QC_plus_cell_state": [
            "log_total_counts", "log_n_genes", "pct_counts_mt", "Epithelial_score",
            "CAF_ECM_score", "EMT_score", "Hypoxia_score", "Proliferation_score",
        ],
    }
    subsets: list[tuple[str, pd.Series]] = [
        ("all_spots", pd.Series(True, index=df.index)),
        ("tumor_stages", df["is_tumor_stage"]),
        ("PT_only", df["stage"].eq("PT")),
        ("PTC_only", df["stage"].eq("PTC")),
        ("LPTC_only", df["stage"].eq("LPTC")),
        ("ATC_only", df["stage"].eq("ATC")),
        ("tumor_high_epithelial_median", df["is_tumor_stage"] & df["high_epithelial_median"]),
        ("tumor_high_epithelial_q75", df["is_tumor_stage"] & df["high_epithelial_q75"]),
        ("tumor_low_epithelial_q25", df["is_tumor_stage"] & df["low_epithelial_q25"]),
    ]
    pairs = [
        ("MAPK9_vs_Panel8", "mapk9_score", "panel8_score"),
        ("MAPK9_vs_TDS16", "mapk9_score", "tds16_score"),
        ("MAPK9_vs_RAI8_existing", "mapk9_score", "RAI_8_score"),
        ("MAPK9_vs_DM1_existing", "mapk9_score", "DM1_like_score"),
        ("MAPK9_vs_TDS_existing", "mapk9_score", "TDS_like_score"),
    ]

    rows = []
    for subset_name, mask in subsets:
        sub = df.loc[mask].copy()
        for pair_name, xcol, ycol in pairs:
            for model, covs in cov_models.items():
                if model == "raw":
                    rho, p, n = spearman(sub[xcol], sub[ycol])
                else:
                    rho, p, n = partial_spearman(sub[xcol], sub[ycol], sub[covs])
                rows.append({
                    "subset": subset_name,
                    "pair": pair_name,
                    "adjustment": model,
                    "n": n,
                    "rho": rho,
                    "p": p,
                    "negative_direction": bool(np.isfinite(rho) and rho < 0),
                    "covariates": ",".join(covs) if covs else "none",
                })
    by_subset = pd.DataFrame(rows)
    by_subset.to_csv(OUT / "v16_spatial_adjustment_grid.tsv", sep="\t", index=False)

    slide_rows = []
    for sample_id, sub in df.groupby("sample_id", sort=False):
        for model, covs in cov_models.items():
            if model == "raw":
                rho, p, n = spearman(sub["mapk9_score"], sub["panel8_score"])
            else:
                rho, p, n = partial_spearman(sub["mapk9_score"], sub["panel8_score"], sub[covs])
            slide_rows.append({
                "sample_id": sample_id,
                "stage": sub["stage"].iloc[0],
                "adjustment": model,
                "n": n,
                "rho_mapk_panel": rho,
                "p": p,
                "negative_direction": bool(np.isfinite(rho) and rho < 0),
            })
    per_slide = pd.DataFrame(slide_rows)
    per_slide.to_csv(OUT / "v16_spatial_per_slide_adjusted.tsv", sep="\t", index=False)

    quartile_rows = []
    for (stage, q), sub in df.groupby(["stage", "epithelial_quartile"], sort=False):
        rho, p, n = spearman(sub["mapk9_score"], sub["panel8_score"])
        full_rho, full_p, full_n = partial_spearman(
            sub["mapk9_score"],
            sub["panel8_score"],
            sub[cov_models["QC_plus_cell_state"]],
        )
        quartile_rows.append({
            "stage": stage,
            "epithelial_quartile": q,
            "n": n,
            "rho_raw": rho,
            "p_raw": p,
            "rho_QC_plus_cell_state": full_rho,
            "p_QC_plus_cell_state": full_p,
            "n_QC_plus_cell_state": full_n,
        })
    quartiles = pd.DataFrame(quartile_rows)
    quartiles.to_csv(OUT / "v16_spatial_epithelial_quartile_grid.tsv", sep="\t", index=False)

    tumor_slides = per_slide[per_slide["stage"].isin(TUMOR_STAGES)]
    pooled = {
        model: fisher_pool(tumor_slides[tumor_slides["adjustment"] == model], "rho_mapk_panel")
        for model in cov_models
    }
    panel_rows = by_subset[(by_subset["pair"] == "MAPK9_vs_Panel8") & (by_subset["adjustment"].isin(cov_models))]
    summary = {
        "n_spots": int(len(df)),
        "n_tumor_spots": int(df["is_tumor_stage"].sum()),
        "n_slides": int(df["sample_id"].nunique()),
        "tumor_slide_pooled_by_adjustment": pooled,
        "negative_raw_slide_count": int(((per_slide["adjustment"] == "raw") & (per_slide["rho_mapk_panel"] < 0)).sum()),
        "negative_full_adjusted_slide_count": int(((per_slide["adjustment"] == "QC_plus_cell_state") & (per_slide["rho_mapk_panel"] < 0)).sum()),
        "grid_negative_count": int(panel_rows["negative_direction"].sum()),
        "grid_test_count": int(panel_rows.shape[0]),
        "bottom_line": "No spatial stratum or adjustment model recovers a robust MAPK x Panel anti-correlation; the positive association is strongly attenuated by QC/cell-state adjustment.",
    }
    print(json.dumps(summary, indent=2)[:1600])
    return by_subset, per_slide, quartiles, summary


def run_prism_diagnostics() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    print("\n[v16B] PRISM diagnostics")
    pr = pd.read_csv(PRISM, sep="\t").sort_values(["fdr", "cohens_d"]).reset_index(drop=True)
    pr["rank"] = np.arange(1, len(pr) + 1)
    n_total = len(pr)
    n_canon = int(pr["is_canonical_mapk"].sum())
    n_axis = int(pr["is_mapk_pi3k_axis"].sum())
    rows = []
    for kind, values in {
        "top_k": [5, 7, 10, 11, 15, 25, 50, 100],
        "fdr_threshold": [0.001, 0.01, 0.05, 0.1],
    }.items():
        for val in values:
            if kind == "top_k":
                sub = pr.head(int(val))
                label = f"top_{int(val)}"
            else:
                sub = pr[pr["fdr"] < float(val)]
                label = f"fdr_lt_{val:g}"
            n = len(sub)
            c = int(sub["is_canonical_mapk"].sum())
            a = int(sub["is_mapk_pi3k_axis"].sum())
            rows.append({
                "set": label,
                "criterion": kind,
                "n_hit": n,
                "n_canonical_mapk": c,
                "n_mapk_pi3k_axis": a,
                "canonical_mapk_fraction": c / n if n else np.nan,
                "mapk_pi3k_fraction": a / n if n else np.nan,
                "canonical_fold_enrichment": (c / n) / (n_canon / n_total) if n and n_canon else np.nan,
                "mapk_pi3k_fold_enrichment": (a / n) / (n_axis / n_total) if n and n_axis else np.nan,
                "canonical_hypergeom_p": hypergeom_p(n_total, n_canon, n, c),
                "mapk_pi3k_hypergeom_p": hypergeom_p(n_total, n_axis, n, a),
            })
    enrich = pd.DataFrame(rows)
    enrich.to_csv(OUT / "v16_prism_enrichment_sensitivity.tsv", sep="\t", index=False)

    cells = pd.read_csv(PRISM_CELLS, sep="\t").dropna(subset=["DM1_like_score", "mapk_inhibitor_mean_LFC"]).copy()
    line_rows = []
    groups = [("all_nonmissing", cells), ("exclude_thyroid", cells[cells["OncotreeLineage"] != "Thyroid"])]
    groups.extend((f"lineage_{lineage}", sub) for lineage, sub in cells.groupby("OncotreeLineage") if len(sub) >= 10)
    for label, sub in groups:
        rho, p, n = spearman(sub["DM1_like_score"], sub["mapk_inhibitor_mean_LFC"])
        q1, q3 = sub["DM1_like_score"].quantile([1 / 3, 2 / 3])
        high = sub.loc[sub["DM1_like_score"] >= q3, "mapk_inhibitor_mean_LFC"].to_numpy()
        low = sub.loc[sub["DM1_like_score"] <= q1, "mapk_inhibitor_mean_LFC"].to_numpy()
        line_rows.append({
            "group": label,
            "n": int(len(sub)),
            "n_unique_dm1_scores": int(sub["DM1_like_score"].nunique()),
            "rho_dm1_vs_mapk_lfc": rho,
            "p": p,
            "high_vs_low_cohens_d_lfc": cohen_d(high, low),
            "median_lfc": float(sub["mapk_inhibitor_mean_LFC"].median()),
            "lineage": label.replace("lineage_", "") if label.startswith("lineage_") else label,
        })
    lineage = pd.DataFrame(line_rows).sort_values(["group"])
    lineage.to_csv(OUT / "v16_prism_lineage_sensitivity.tsv", sep="\t", index=False)

    top7 = enrich[enrich["set"] == "top_7"].iloc[0]
    top11 = enrich[enrich["set"] == "top_11"].iloc[0]
    top15 = enrich[enrich["set"] == "top_15"].iloc[0]
    fdr05 = enrich[enrich["set"] == "fdr_lt_0.05"].iloc[0]
    all_row = lineage[lineage["group"] == "all_nonmissing"].iloc[0]
    ex_thy = lineage[lineage["group"] == "exclude_thyroid"].iloc[0]
    thyroid = lineage[lineage["group"] == "lineage_Thyroid"].iloc[0] if (lineage["group"] == "lineage_Thyroid").any() else None
    summary = {
        "n_total_prism_drugs": int(n_total),
        "n_canonical_mapk_universe": int(n_canon),
        "n_mapk_pi3k_universe": int(n_axis),
        "top7_canonical_mapk": f"{int(top7['n_canonical_mapk'])}/{int(top7['n_hit'])}",
        "top11_canonical_mapk": f"{int(top11['n_canonical_mapk'])}/{int(top11['n_hit'])}",
        "top15_canonical_mapk": f"{int(top15['n_canonical_mapk'])}/{int(top15['n_hit'])}",
        "fdr05_canonical_mapk": f"{int(fdr05['n_canonical_mapk'])}/{int(fdr05['n_hit'])}",
        "fdr05_canonical_hypergeom_p": float(fdr05["canonical_hypergeom_p"]),
        "all_celllines_rho": float(all_row["rho_dm1_vs_mapk_lfc"]),
        "all_celllines_p": float(all_row["p"]),
        "exclude_thyroid_rho": float(ex_thy["rho_dm1_vs_mapk_lfc"]),
        "exclude_thyroid_p": float(ex_thy["p"]),
        "thyroid_line_n": int(thyroid["n"]) if thyroid is not None else 0,
        "thyroid_unique_dm1_scores": int(thyroid["n_unique_dm1_scores"]) if thyroid is not None else 0,
        "bottom_line": "PRISM support is real but narrower than 'top 15 all MAPK': top 7 are canonical MAPK; 7/11 FDR<0.05 are canonical MAPK; top 15 includes non-MAPK hits.",
    }
    print(json.dumps(summary, indent=2))
    return enrich, lineage, pr, summary


def plot(by_subset: pd.DataFrame, quartiles: pd.DataFrame, enrich: pd.DataFrame, lineage: pd.DataFrame) -> None:
    print("\n[v16] plotting")
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    ax1, ax2, ax3, ax4 = axes.ravel()

    show_subsets = ["all_spots", "tumor_stages", "PTC_only", "LPTC_only", "ATC_only", "tumor_high_epithelial_q75"]
    show_adj = ["raw", "QC", "QC_plus_cell_state"]
    p = by_subset[(by_subset["pair"] == "MAPK9_vs_Panel8") & (by_subset["subset"].isin(show_subsets)) & (by_subset["adjustment"].isin(show_adj))].copy()
    p["subset"] = pd.Categorical(p["subset"], show_subsets, ordered=True)
    p["adjustment"] = pd.Categorical(p["adjustment"], show_adj, ordered=True)
    width = 0.24
    x = np.arange(len(show_subsets))
    colors = {"raw": "#d8453a", "QC": "#ff9b3d", "QC_plus_cell_state": "#2c7fb8"}
    for i, adj in enumerate(show_adj):
        vals = p[p["adjustment"] == adj].sort_values("subset")["rho"].to_numpy()
        ax1.bar(x + (i - 1) * width, vals, width=width, color=colors[adj], label=adj, edgecolor="black", lw=0.4)
    ax1.axhline(0, color="#555", lw=0.9)
    ax1.set_xticks(x)
    ax1.set_xticklabels(["all", "tumor", "PTC", "LPTC", "ATC", "tumor epi Q4"], rotation=25, ha="right")
    ax1.set_ylabel("Spearman / partial Spearman rho")
    ax1.set_title("A. Spatial MAPK x Panel: adjustment attenuates but does not flip", loc="left", fontsize=10, weight="bold")
    ax1.legend(frameon=False, fontsize=8)
    ax1.grid(axis="y", alpha=0.25)

    q = quartiles.pivot(index="stage", columns="epithelial_quartile", values="rho_raw").reindex(["PT", "PTC", "LPTC", "ATC"])
    q = q[["Q1_low", "Q2", "Q3", "Q4_high"]]
    im = ax2.imshow(q.to_numpy(float), vmin=-0.1, vmax=0.65, cmap="RdYlBu_r", aspect="auto")
    ax2.set_xticks(np.arange(q.shape[1]))
    ax2.set_xticklabels(q.columns, rotation=25, ha="right")
    ax2.set_yticks(np.arange(q.shape[0]))
    ax2.set_yticklabels(q.index)
    for i in range(q.shape[0]):
        for j in range(q.shape[1]):
            ax2.text(j, i, f"{q.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8, color="black")
    ax2.set_title("B. Spatial rho by within-slide epithelial-score quartile", loc="left", fontsize=10, weight="bold")
    fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)

    top = enrich[enrich["criterion"] == "top_k"].copy()
    ax3.plot(top["n_hit"], top["canonical_mapk_fraction"], marker="o", color="#d8453a", label="canonical MAPK")
    ax3.plot(top["n_hit"], top["mapk_pi3k_fraction"], marker="s", color="#2c7fb8", label="MAPK + PI3K")
    baseline_c = (
        top["canonical_mapk_fraction"].iloc[0] / top["canonical_fold_enrichment"].iloc[0]
        if len(top) and np.isfinite(top["canonical_fold_enrichment"].iloc[0])
        else np.nan
    )
    ax3.axhline(baseline_c, color="#d8453a", ls=":", lw=0.8)
    ax3.set_xlabel("Top-k drugs by FDR rank")
    ax3.set_ylabel("Fraction of hits")
    ax3.set_ylim(0, 1.05)
    ax3.set_title("C. PRISM enrichment: top 7 all MAPK, top 15 not all MAPK", loc="left", fontsize=10, weight="bold")
    ax3.legend(frameon=False, fontsize=8)
    ax3.grid(alpha=0.25)

    lin = lineage[(lineage["group"].isin(["all_nonmissing", "exclude_thyroid"])) | (lineage["group"].str.startswith("lineage_"))].copy()
    lin = lin[np.isfinite(lin["rho_dm1_vs_mapk_lfc"])]
    lin = lin.sort_values("rho_dm1_vs_mapk_lfc").head(12)
    y = np.arange(len(lin))
    ax4.barh(y, lin["rho_dm1_vs_mapk_lfc"], color=["#d8453a" if v < 0 else "#2c7fb8" for v in lin["rho_dm1_vs_mapk_lfc"]], edgecolor="black", lw=0.4)
    ax4.set_yticks(y)
    ax4.set_yticklabels([f"{r.lineage} n={int(r.n)}" for _, r in lin.iterrows()], fontsize=8)
    ax4.axvline(0, color="#555", lw=0.8)
    ax4.set_xlabel("Spearman rho: DM1 score vs MAPK-inhibitor mean LFC")
    ax4.set_title("D. PRISM lineage sensitivity", loc="left", fontsize=10, weight="bold")
    ax4.grid(axis="x", alpha=0.25)

    fig.suptitle("Supplementary Figure SX_v16. Spatial and PRISM diagnostic stress tests", fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT / "Fig_SX_v16_spatial_prism_diagnostics.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "Fig_SX_v16_spatial_prism_diagnostics.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    spatial = load_spatial()
    by_subset, per_slide, quartiles, spatial_summary = run_spatial_diagnostics(spatial)
    enrich, lineage, prism, prism_summary = run_prism_diagnostics()
    summary = {"spatial": spatial_summary, "prism": prism_summary}
    (OUT / "v16_spatial_prism_diagnostics_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    plot(by_subset, quartiles, enrich, lineage)
    print("\n=== DONE v16 diagnostics ===")


if __name__ == "__main__":
    main()
