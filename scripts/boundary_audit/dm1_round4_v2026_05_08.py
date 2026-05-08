#!/usr/bin/env python3
"""DM1 Round 4 — methylation 4-pillar + Landa + sub-A/B + spatial — Paper 1 paper-blocking.

Layers:
  1. Methylation 4th pillar (TCGA n=504, 8-gene beta values + DM-stratified)
  2. Landa GSE76039 add to robustness (n=38 PDTC/ATC)
  3. DM1 sub-A vs sub-B detail (age/stage/OS/sub-DEG)
  4. Spatial niche cross-stage (GSE250521 ATC/LPTC/PTC/normal)

Outputs in project/results/dm1_robustness_v2026_05_08/round4/.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "project" / "results" / "dm1_robustness_v2026_05_08" / "round4"
OUT.mkdir(parents=True, exist_ok=True)

METHYL_SAMPLES = Path("/data/thca/repo_results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv")
METHYL_PERGENE = Path("/data/thca/repo_results/audit_2026_04_30/round5/r5_2_per_gene_methylation_DM.tsv")
METHYL_EXPR_CORR = Path("/data/thca/repo_results/aggressive_sprint_2026_05_06/methylation_expression_correlations.tsv")
LANDA_EXPR = Path("/data/thca/data_processed/microarray_v3/GSE76039_v3_log2.tsv")
LANDA_META = Path("/data/thca/repo_results/p_landa_2016/sample_metadata.tsv")
LANDA_SCORES = Path("/data/thca/repo_results/p_landa_2016/scores.tsv")
DM_SIG = REPO / "project" / "results" / "d4p2_tcga_hashimoto_signature" / "tcga_signature_scores.tsv"
SUB_CLIN = REPO / "project" / "results" / "d6p7_dm1_subcluster" / "subcluster_clinical.tsv"
SUB_DEG = REPO / "project" / "results" / "d6p7_dm1_subcluster" / "dm1_subBvA_deg.tsv"
SUB_LABELS = REPO / "project" / "results" / "d6p7_dm1_subcluster" / "dm1_subcluster_labels.tsv"
SPATIAL_NICHE = Path("/data/thca/repo_results/spatial_full_2026_05_06/spatial_A_DM1_niche_density.tsv")
SPATIAL_INTEG = Path("/data/thca/repo_results/spatial_full_2026_05_06/spatial_F_cross_cohort_integration.tsv")

PANEL = ["TG", "TPO", "TSHR", "PAX8", "FOXE1", "NKX2-1", "DIO1", "SLC5A5"]


def cohens_d(x: np.ndarray, y: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    sp = math.sqrt(((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / (len(x) + len(y) - 2))
    return (x.mean() - y.mean()) / sp if sp > 0 else float("nan")


# ============== 1. Methylation 4th pillar ==============
def methylation_pillar():
    print("\n[1] Methylation 4th pillar (TCGA n=504)")
    msam = pd.read_csv(METHYL_SAMPLES, sep="\t")
    mgene = pd.read_csv(METHYL_PERGENE, sep="\t")
    print(f"  methyl samples: {len(msam)}, per-gene rows: {len(mgene)}")

    # Compute mean_8g_beta DM1 vs DM2 from per-sample joining DM call
    sig = pd.read_csv(DM_SIG, sep="\t", index_col=0)
    sig["sample_short"] = sig.index.to_series().apply(lambda s: "-".join(str(s).split("-")[:3]))
    msam["sample_short"] = msam["sample_short"].astype(str)
    merged = msam.merge(sig.reset_index()[["sample_short", "DM"]], on="sample_short", how="inner")
    print(f"  merged with DM call: {len(merged)}; DM1: {(merged['DM']=='DM1').sum()}, DM2: {(merged['DM']=='DM2').sum()}")

    d_panel = float("nan")
    p_panel = float("nan")
    if len(merged) > 0 and merged["DM"].nunique() > 1:
        a = merged.loc[merged["DM"] == "DM1", "mean_8g_beta"].to_numpy()
        b = merged.loc[merged["DM"] == "DM2", "mean_8g_beta"].to_numpy()
        d_panel = cohens_d(a, b)
        try:
            p_panel = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except Exception:
            p_panel = float("nan")
        print(f"  mean_8g_beta DM1 vs DM2: d={d_panel:+.3f}, p={p_panel:.2e}")

    # Recompute per-gene Cohen's d using CURRENT DM convention
    # (per-gene file uses an older DM convention; we want the modern call).
    per_gene_rows = []
    for gene in PANEL:
        if gene not in merged.columns:
            continue
        a = merged.loc[merged["DM"] == "DM1", gene].dropna().to_numpy()
        b = merged.loc[merged["DM"] == "DM2", gene].dropna().to_numpy()
        if len(a) < 2 or len(b) < 2:
            continue
        dg = cohens_d(a, b)
        try:
            pg = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except Exception:
            pg = float("nan")
        per_gene_rows.append({
            "gene": gene, "n_DM1": int(len(a)), "n_DM2": int(len(b)),
            "mean_β_DM1": float(a.mean()), "mean_β_DM2": float(b.mean()),
            "Δβ": float(a.mean() - b.mean()), "cohens_d": float(dg), "mw_p": pg,
        })
    pg_current = pd.DataFrame(per_gene_rows)
    pg_current.to_csv(OUT / "methylation_per_gene_current_DM_convention.tsv", sep="\t", index=False)
    # Use the recomputed table for the figure (consistent with panel-level result)
    mgene_sorted = pg_current.sort_values("cohens_d")
    fig, ax = plt.subplots(figsize=(8, 5))
    y = np.arange(len(mgene_sorted))
    bars = ax.barh(y, mgene_sorted["cohens_d"],
                   color=["#d62728" if d > 0 else "#1f77b4" for d in mgene_sorted["cohens_d"]],
                   edgecolor="black", alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['gene']}  Δβ={r['Δβ']:+.2f}, p={r['mw_p']:.1e}" for _, r in mgene_sorted.iterrows()], fontsize=9)
    ax.axvline(0, color="#888", linestyle=":", linewidth=0.8)
    ax.set_xlabel("Cohen's d (DM1 vs DM2 promoter β-value)")
    ax.set_title(
        "Methylation 4th pillar — TCGA-THCA HM450 promoter β per panel gene by DM call\n"
        "(positive d = DM1 has HIGHER β = MORE methylation = expression silencing pattern)",
        fontsize=10)
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "methylation_per_gene.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "methylation_per_gene.pdf", bbox_inches="tight")
    plt.close(fig)

    # Mean β panel score boxplot DM1 vs DM2
    fig2, ax2 = plt.subplots(figsize=(6, 4.5))
    if len(merged) > 0:
        groups_data = [merged.loc[merged["DM"] == g, "mean_8g_beta"].dropna().to_numpy() for g in ("DM1", "DM2")]
        bp = ax2.boxplot(groups_data, tick_labels=[f"DM1\nn={len(groups_data[0])}", f"DM2\nn={len(groups_data[1])}"],
                         patch_artist=True, showfliers=False)
        for patch, color in zip(bp["boxes"], ["#fdae61", "#abd9e9"]):
            patch.set_facecolor(color)
            patch.set_edgecolor("#222")
        rng = np.random.default_rng(42)
        for i, gd in enumerate(groups_data, start=1):
            ax2.scatter(rng.normal(i, 0.05, len(gd)), gd, s=12, color="#222", alpha=0.55, zorder=3)
        ax2.set_ylabel("mean_8g_beta (8-gene promoter β-value)")
        ax2.set_title(f"4th pillar — TCGA HM450 mean 8-gene β by DM (d={d_panel:+.2f}, p={p_panel:.2g})", fontsize=10)
        ax2.grid(axis="y", linestyle=":", alpha=0.3)
        plt.tight_layout()
        fig2.savefig(OUT / "methylation_mean_panel.png", dpi=160, bbox_inches="tight")
        fig2.savefig(OUT / "methylation_mean_panel.pdf", bbox_inches="tight")
        plt.close(fig2)

    # Methylation × expression correlation table
    if METHYL_EXPR_CORR.exists():
        corr = pd.read_csv(METHYL_EXPR_CORR, sep="\t")
        corr.to_csv(OUT / "methylation_expression_correlations.tsv", sep="\t", index=False)

    return mgene_sorted, merged


# ============== 2. Landa GSE76039 add to robustness ==============
def landa_add():
    print("\n[2] Landa GSE76039 (n=38 PDTC/ATC)")
    if not LANDA_EXPR.exists():
        print("  Landa expression matrix not found — skipping")
        return None
    expr = pd.read_csv(LANDA_EXPR, sep="\t", index_col=0)
    print(f"  Landa matrix: {expr.shape[0]} genes × {expr.shape[1]} samples")
    panel_in = [g for g in PANEL if g in expr.index]
    print(f"  panel genes present: {panel_in}")

    sub = expr.loc[panel_in]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    score = z.mean(axis=0)

    score_df = score.to_frame("g8_RAI")
    if LANDA_META.exists():
        meta = pd.read_csv(LANDA_META, sep="\t").set_index("sample_id")
        score_df["group"] = meta["histology"].reindex(score_df.index).values
    else:
        score_df["group"] = "unknown"
    score_df.to_csv(OUT / "landa_per_sample.tsv", sep="\t")
    counts = score_df["group"].value_counts().to_dict()
    print(f"  Landa group counts: {counts}")

    # Effect size PDTC vs ATC
    rows = []
    for la, lb in [("PDTC", "ATC")]:
        a = score_df.loc[score_df["group"] == la, "g8_RAI"].dropna().to_numpy()
        b = score_df.loc[score_df["group"] == lb, "g8_RAI"].dropna().to_numpy()
        if len(a) < 2 or len(b) < 2:
            continue
        dval = cohens_d(a, b)
        pval = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue) if len(a) > 1 and len(b) > 1 else float("nan")
        rows.append({"contrast": f"{la} vs {lb}", "n_a": int(len(a)), "n_b": int(len(b)),
                     "mean_a": float(a.mean()), "mean_b": float(b.mean()),
                     "cohens_d": float(dval), "mwu_p": pval})
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "landa_contrasts.tsv", sep="\t", index=False)
    print(res.to_string(index=False))

    # boxplot
    if len(score_df.dropna()) > 5:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        order = ["PDTC", "ATC"]
        gs = [g for g in order if (score_df["group"] == g).any()]
        data = [score_df.loc[score_df["group"] == g, "g8_RAI"].dropna().to_numpy() for g in gs]
        bp = ax.boxplot(data, tick_labels=[f"{g}\nn={len(d)}" for g, d in zip(gs, data)],
                        patch_artist=True, showfliers=False)
        palette = {"PDTC": "#fdae61", "ATC": "#d7191c"}
        for patch, g in zip(bp["boxes"], gs):
            patch.set_facecolor(palette[g])
            patch.set_edgecolor("#222")
        rng = np.random.default_rng(42)
        for i, gd in enumerate(data, start=1):
            ax.scatter(rng.normal(i, 0.05, len(gd)), gd, s=12, color="#222", alpha=0.55, zorder=3)
        ax.set_ylabel("8-gene panel score (within-sample z mean)")
        title_d = res["cohens_d"].iloc[0] if len(res) else float("nan")
        title_p = res["mwu_p"].iloc[0] if len(res) else float("nan")
        ax.set_title(f"Landa 2016 GSE76039 (n=37) — 8-gene panel PDTC vs ATC  (d={title_d:+.2f}, p={title_p:.2g})", fontsize=10)
        ax.grid(axis="y", linestyle=":", alpha=0.3)
        plt.tight_layout()
        fig.savefig(OUT / "landa_box.png", dpi=160, bbox_inches="tight")
        fig.savefig(OUT / "landa_box.pdf", bbox_inches="tight")
        plt.close(fig)
    return res


# ============== 3. DM1 sub-A vs sub-B detail figure ==============
def sub_AB_detail():
    print("\n[3] DM1 sub-A vs sub-B detail")
    clin = pd.read_csv(SUB_CLIN, sep="\t")
    deg = pd.read_csv(SUB_DEG, sep="\t")
    labels = pd.read_csv(SUB_LABELS, sep="\t", index_col=0) if SUB_LABELS.exists() else None
    print(f"  sub clinical metrics: {len(clin)}, DEGs: {len(deg)}")

    # Top DEGs (by |cohen_d|, padj < 0.05)
    sig_deg = deg[(deg["padj"].astype(float) < 0.05)].copy()
    sig_deg["abs_d"] = sig_deg["cohen_d_B_vs_A"].abs()
    top_up = sig_deg[sig_deg["cohen_d_B_vs_A"] > 0].sort_values("abs_d", ascending=False).head(15)
    top_dn = sig_deg[sig_deg["cohen_d_B_vs_A"] < 0].sort_values("abs_d", ascending=False).head(15)
    print(f"  significant DEGs: {len(sig_deg)}; up sub-B: {(sig_deg['cohen_d_B_vs_A'] > 0).sum()}, down: {(sig_deg['cohen_d_B_vs_A'] < 0).sum()}")

    # 2-panel figure: clinical metrics + top DEG bar
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(5, 0.35 * (len(top_up) + len(top_dn)) + 2)))

    # Panel 1: clinical metrics with cohen d
    clin = clin.dropna(subset=["cohen_d"])
    y1 = np.arange(len(clin))
    colors1 = ["#d62728" if r["cohen_d"] > 0 else "#1f77b4" for _, r in clin.iterrows()]
    ax1.barh(y1, clin["cohen_d"], color=colors1, edgecolor="black", alpha=0.85)
    ax1.set_yticks(y1)
    ax1.set_yticklabels([f"{r['metric']}  Δ={r['mean_B'] - r['mean_A']:+.2f}, p={r['mw_p']:.2g}" for _, r in clin.iterrows()], fontsize=9)
    ax1.axvline(0, color="#888", linestyle=":", linewidth=0.8)
    ax1.set_xlabel("Cohen's d (sub_B vs sub_A)")
    ax1.set_title("DM1 sub-B vs sub-A — clinical features", fontsize=10)
    ax1.grid(axis="x", linestyle=":", alpha=0.3)
    ax1.invert_yaxis()

    # Panel 2: top up + down DEGs
    deg_show = pd.concat([top_up, top_dn]).sort_values("cohen_d_B_vs_A")
    y2 = np.arange(len(deg_show))
    colors2 = ["#d62728" if d > 0 else "#1f77b4" for d in deg_show["cohen_d_B_vs_A"]]
    ax2.barh(y2, deg_show["cohen_d_B_vs_A"], color=colors2, edgecolor="black", alpha=0.85)
    ax2.set_yticks(y2)
    ax2.set_yticklabels([f"{r['gene']}  log2FC={r['log2FC_B_vs_A']:+.2f}, padj={r['padj']:.1e}" for _, r in deg_show.iterrows()], fontsize=7)
    ax2.axvline(0, color="#888", linestyle=":", linewidth=0.8)
    ax2.set_xlabel("Cohen's d (sub_B vs sub_A)")
    ax2.set_title("DM1 sub-B vs sub-A — top 15 up + top 15 down DEGs", fontsize=10)
    ax2.grid(axis="x", linestyle=":", alpha=0.3)
    ax2.invert_yaxis()

    fig.suptitle("DM1 sub-A vs sub-B characterization (n=84 vs 56) — Round 4 detail (per memory v17_D6P7)", fontsize=11, y=1.01)
    plt.tight_layout()
    fig.savefig(OUT / "sub_AB_detail.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "sub_AB_detail.pdf", bbox_inches="tight")
    plt.close(fig)

    # also save concise summary
    summary = {
        "n_clinical_metrics": int(len(clin)),
        "n_significant_DEGs": int(len(sig_deg)),
        "n_DEG_up_sub_B": int((sig_deg["cohen_d_B_vs_A"] > 0).sum()),
        "n_DEG_down_sub_B": int((sig_deg["cohen_d_B_vs_A"] < 0).sum()),
        "top_up_genes": top_up["gene"].head(10).tolist(),
        "top_down_genes": top_dn["gene"].head(10).tolist(),
    }
    (OUT / "sub_AB_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


# ============== 4. Spatial niche cross-stage ==============
def spatial_niche():
    print("\n[4] Spatial niche cross-stage (GSE250521)")
    if not SPATIAL_NICHE.exists():
        return None
    n = pd.read_csv(SPATIAL_NICHE, sep="\t")
    integ = pd.read_csv(SPATIAL_INTEG, sep="\t") if SPATIAL_INTEG.exists() else pd.DataFrame()

    # Plot DM1 vs TDS Spearman ρ across stages
    fig, ax = plt.subplots(figsize=(9, 5))
    if "stage" in n.columns and "rho_DM1_vs_TDS" in n.columns:
        order = ["normal", "PTC", "LPTC", "ATC"]
        n["stage"] = n["stage"].astype(str)
        present_stages = [s for s in order if (n["stage"] == s).any()]
        x = np.arange(len(present_stages))
        rhos = [float(n.loc[n["stage"] == s, "rho_DM1_vs_TDS"].iloc[0]) for s in present_stages]
        ns = [int(n.loc[n["stage"] == s, "n_spots"].iloc[0]) for s in present_stages]
        bars = ax.bar(x, rhos, color=["#2c7fb8", "#fdae61", "#d7191c", "#9e006e"][:len(present_stages)],
                      edgecolor="black", alpha=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{s}\nn={n_spots:,} spots" for s, n_spots in zip(present_stages, ns)], fontsize=10)
        ax.set_ylabel("Spearman ρ (DM1 score vs TDS score)\nper-spot, within stage")
        ax.set_title("Spatial DM1 ↔ TDS niche correlation by stage (GSE250521 Visium)\nDM1 score is consistently anti-correlated with TDS (thyroid differentiation), strongest in LPTC/PTC", fontsize=10)
        ax.axhline(0, color="#888", linestyle=":", linewidth=0.8)
        ax.set_ylim(-1, 0.2)
        ax.grid(axis="y", linestyle=":", alpha=0.3)
        for i, (bar, rho) in enumerate(zip(bars, rhos)):
            ax.text(bar.get_x() + bar.get_width() / 2, rho - 0.04, f"ρ={rho:+.3f}",
                    ha="center", fontsize=10, color="white" if abs(rho) > 0.5 else "black")
    plt.tight_layout()
    fig.savefig(OUT / "spatial_niche_cross_stage.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "spatial_niche_cross_stage.pdf", bbox_inches="tight")
    plt.close(fig)

    n.to_csv(OUT / "spatial_niche_per_stage.tsv", sep="\t", index=False)
    return n


def main():
    res = {}
    mg, merged = methylation_pillar()
    res["methylation_per_gene"] = mg.to_dict(orient="records")
    res["methylation_n_merged_with_DM"] = int(len(merged))
    landa = landa_add()
    res["landa_contrasts"] = landa.to_dict(orient="records") if landa is not None else []
    res["sub_AB_detail"] = sub_AB_detail()
    sp = spatial_niche()
    if sp is not None:
        res["spatial_niche"] = sp.to_dict(orient="records")
    res["generated_at"] = "2026-05-08 round4"

    with open(OUT / "round4_summary.json", "w") as f:
        json.dump(res, f, indent=2, default=str)

    print(f"\nwrote {OUT.relative_to(REPO)}/")
    print("Outputs:")
    for f in sorted(OUT.iterdir()):
        print(f"  {f.relative_to(OUT)}")


if __name__ == "__main__":
    main()
