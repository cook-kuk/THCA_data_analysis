#!/usr/bin/env python3
"""DM1 Round 5 — pan-cancer prognostic + fusion + K2 + GSEA + driver × DM.

Layers:
  1. Pan-cancer DM1 Cox forest (per memory paper11_pancancer)
  2. Fusion enrichment in DM1 (Paper 1 claim P1_C2)
  3. K2 (PRJEB11591) Korean cohort 8-gene distribution
  4. GSEA top pathway enrichment in GSE286332 (PTC+HT vs PTC)
  5. driver_anchor × DM stratified table

Outputs in project/results/dm1_robustness_v2026_05_08/round5/.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "project" / "results" / "dm1_robustness_v2026_05_08" / "round5"
OUT.mkdir(parents=True, exist_ok=True)

PANCAN_COX = Path("/data/thca/repo_results/paper11_pancancer/phase_B_survival/cox_per_lineage.tsv")
PANCAN_STATS = Path("/data/thca/repo_results/paper11_pancancer/pancan_lineage_stats.tsv")
FUSION = Path("/data/thca/repo_results/v17/tables/fusion_calls_per_sample.tsv")
DM_SIG = REPO / "project" / "results" / "d4p2_tcga_hashimoto_signature" / "tcga_signature_scores.tsv"
SAMPLE_MASTER = Path("/data/thca/repo_results/v17/tables/sample_master_v17_tert.tsv")
K2_TPM = Path("/data/thca/repo_results/v17_korean/K2_8gene_tpm_matrix_v4.tsv")
K2_SUMMARY = Path("/data/thca/repo_results/v17_korean/K2_PRJEB11591_q13_summary.json")
GSEA_HALLMARK = REPO / "project" / "results" / "p3_gse286332" / "gsea_MSigDB_Hallmark_2020.tsv"


# ============== 1. Pan-cancer Cox forest ==============
def pancan_prognostic():
    print("\n[1] Pan-cancer DM1 Cox forest")
    cox = pd.read_csv(PANCAN_COX, sep="\t")
    print(f"  cox rows: {len(cox)}; lineages: {cox['lineage'].nunique()}")

    # Significant FDR<0.1 OS or DSS
    sig = cox[(cox["fdr"] < 0.1) & (cox["outcome"].isin(["OS", "DSS"]))].copy()
    print(f"  significant (FDR<0.1, OS|DSS): {len(sig)}")
    sig.to_csv(OUT / "pancan_significant_lineages.tsv", sep="\t", index=False)

    # Forest
    plot_df = cox[(cox["fdr"] < 0.25) & (cox["HR"] > 0) & (cox["HR_upper95"] < 1e3)].copy()
    plot_df = plot_df.sort_values(["outcome", "HR"], ascending=[True, False])
    if len(plot_df):
        fig, ax = plt.subplots(figsize=(11, max(5, 0.32 * len(plot_df) + 1.5)))
        y = np.arange(len(plot_df))
        colors = []
        for _, r in plot_df.iterrows():
            if r["fdr"] < 0.05 and r["HR"] > 1:
                colors.append("#d62728")
            elif r["fdr"] < 0.1 and r["HR"] > 1:
                colors.append("#fdae61")
            else:
                colors.append("#7f7f7f")
        for i, (_, r) in enumerate(plot_df.iterrows()):
            ax.errorbar(r["HR"], i,
                        xerr=[[max(r["HR"] - r["HR_lower95"], 0)], [max(r["HR_upper95"] - r["HR"], 0)]],
                        fmt="s", color=colors[i], ecolor=colors[i], capsize=3, markersize=6)
        ax.axvline(1.0, color="#888", linestyle=":", linewidth=0.8)
        ax.set_xscale("log")
        ax.set_yticks(y)
        labels = [f"{r['outcome']} | {r['lineage']}  n={r['n']}, ev={r['events']}, FDR={r['fdr']:.2g}" for _, r in plot_df.iterrows()]
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("Hazard ratio (DM1_like + age covariates; log scale)", fontsize=9)
        ax.set_title(
            f"Pan-cancer DM1 prognostic Cox HR forest (per memory paper11_pancancer)\n"
            f"FDR<0.05 red ({(plot_df['fdr']<0.05).sum()}), FDR<0.10 orange ({((plot_df['fdr']<0.1) & (plot_df['fdr']>=0.05)).sum()}), grey rest", fontsize=10)
        ax.grid(axis="x", linestyle=":", alpha=0.3)
        plt.tight_layout()
        fig.savefig(OUT / "pancan_cox_forest.png", dpi=160, bbox_inches="tight")
        fig.savefig(OUT / "pancan_cox_forest.pdf", bbox_inches="tight")
        plt.close(fig)

    return sig


# ============== 2. Fusion enrichment in DM1 ==============
def fusion_dm():
    print("\n[2] Fusion enrichment in DM1")
    fus = pd.read_csv(FUSION, sep="\t")
    sig = pd.read_csv(DM_SIG, sep="\t", index_col=0)
    fus["sample_short"] = fus["sample_id"].astype(str)
    sig["sample_short"] = sig.index.to_series().astype(str)
    merged = sig.merge(fus[["sample_short", "fusion_classes", "v3_anchor_6class"]], on="sample_short", how="left")
    print(f"  merged: {len(merged)}; fusion_classes distribution:")
    print(merged["fusion_classes"].value_counts())

    # Build binary fusion flag (anything other than 'other'/null = fusion-positive)
    merged["fusion_pos"] = (~merged["fusion_classes"].isin(["other", "", np.nan])) & merged["fusion_classes"].notna()
    merged["fusion_pos"] = merged["fusion_pos"].astype(int)

    ct = pd.crosstab(merged["DM"], merged["fusion_pos"])
    ct.to_csv(OUT / "fusion_dm_crosstab.tsv", sep="\t")
    print(f"\n  DM × fusion_pos crosstab:\n{ct}")

    # Fisher
    if ct.shape == (2, 2):
        try:
            or_, p = stats.fisher_exact(ct.values)
            print(f"  Fisher OR={or_:.3f}, p={p:.2e}")
            res = {"DM_fusion_OR": float(or_), "DM_fusion_p": float(p),
                   "n_DM1_fusion_pos": int(ct.loc["DM1", 1]) if "DM1" in ct.index and 1 in ct.columns else 0,
                   "n_DM1_fusion_neg": int(ct.loc["DM1", 0]) if "DM1" in ct.index and 0 in ct.columns else 0,
                   "n_DM2_fusion_pos": int(ct.loc["DM2", 1]) if "DM2" in ct.index and 1 in ct.columns else 0,
                   "n_DM2_fusion_neg": int(ct.loc["DM2", 0]) if "DM2" in ct.index and 0 in ct.columns else 0}
        except Exception:
            res = {"err": "fisher failed"}
    else:
        res = {"shape": str(ct.shape)}

    # Per fusion-class breakdown (informative — RET, NTRK, BRAF-fusion etc.)
    by_class = pd.crosstab(merged["fusion_classes"], merged["DM"])
    by_class["DM1_pct"] = by_class.get("DM1", 0) / (by_class.get("DM1", 0) + by_class.get("DM2", 0)) * 100
    by_class.to_csv(OUT / "fusion_class_by_DM.tsv", sep="\t")
    print(f"\n  Per-fusion-class DM distribution:\n{by_class}")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 4))
    show_classes = by_class[(by_class.get("DM1", 0) + by_class.get("DM2", 0)) >= 5].copy().sort_values("DM1_pct", ascending=False)
    if len(show_classes):
        x = np.arange(len(show_classes))
        ax.bar(x, show_classes["DM1_pct"], color="#d62728", alpha=0.85, edgecolor="black")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{idx}\n(n={int(show_classes.loc[idx].get('DM1', 0) + show_classes.loc[idx].get('DM2', 0))})" for idx in show_classes.index],
                            rotation=30, ha="right", fontsize=9)
        ax.axhline(50, color="#888", linestyle=":", linewidth=0.8)
        ax.set_ylabel("% DM1 within fusion class")
        ax.set_title(f"DM1 enrichment by fusion class (TCGA-THCA, n_total={len(merged)})\n(50% line = baseline)", fontsize=10)
        for i, idx in enumerate(show_classes.index):
            ax.text(i, show_classes.loc[idx, "DM1_pct"] + 2, f"{show_classes.loc[idx, 'DM1_pct']:.0f}%", ha="center", fontsize=9)
        ax.set_ylim(0, 110)
        ax.grid(axis="y", linestyle=":", alpha=0.3)
        plt.tight_layout()
        fig.savefig(OUT / "fusion_dm_enrichment.png", dpi=160, bbox_inches="tight")
        fig.savefig(OUT / "fusion_dm_enrichment.pdf", bbox_inches="tight")
        plt.close(fig)

    return res


# ============== 3. K2 cohort distribution ==============
def k2_layer():
    print("\n[3] K2 (PRJEB11591) Korean cohort")
    summary = json.loads(K2_SUMMARY.read_text())
    print(f"  K2 summary: {summary['n_samples_processed_so_far']}/262 samples; DM1={summary['n_DM1']}, DM2={summary['n_DM2']}, %DM2={summary['pct_DM2']}%")
    print(f"  TCGA baseline: {summary['tcga_dm1_dm2_baseline']}")

    if K2_TPM.exists():
        k2 = pd.read_csv(K2_TPM, sep="\t", index_col=0)
        print(f"  K2 TPM matrix: {k2.shape}")
        # Within-sample z then panel mean
        z = k2.sub(k2.mean(axis=1), axis=0).div(k2.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
        score = z.mean(axis=0).rename("g8_RAI")
        score.to_csv(OUT / "k2_per_sample_score.tsv", sep="\t")

        # Distribution plot
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(score.dropna().to_numpy(), bins=30, color="#abd9e9", edgecolor="black", alpha=0.85)
        ax.axvline(score.dropna().median(), color="#d62728", linestyle="--", linewidth=1, label=f"K2 median={score.dropna().median():+.2f}")
        ax.axvline(0, color="#888", linestyle=":", linewidth=0.8, label="mean reference")
        ax.set_xlabel("8-gene RAI score (within-sample z mean)")
        ax.set_ylabel("# K2 samples")
        ax.set_title(
            f"K2 (PRJEB11591 Yoo SK Korean) 8-gene panel distribution\n"
            f"n_processed={k2.shape[1]}; %DM2={summary.get('pct_DM2', '?')}% per memory v17_korean_k2_calibration",
            fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(axis="y", linestyle=":", alpha=0.3)
        plt.tight_layout()
        fig.savefig(OUT / "k2_distribution.png", dpi=160, bbox_inches="tight")
        fig.savefig(OUT / "k2_distribution.pdf", bbox_inches="tight")
        plt.close(fig)

    return summary


# ============== 4. GSEA top pathways ==============
def gsea_layer():
    print("\n[4] GSEA top pathways (GSE286332 PTC+HT vs PTC)")
    g = pd.read_csv(GSEA_HALLMARK, sep="\t")
    g = g.dropna(subset=["NES", "FDR q-val"]).copy()
    g["FDR q-val"] = pd.to_numeric(g["FDR q-val"], errors="coerce")
    g_sig = g[g["FDR q-val"] < 0.05].copy()
    g_sig["abs_NES"] = g_sig["NES"].abs()
    g_top = g_sig.sort_values("abs_NES", ascending=False).head(20)

    g_top.to_csv(OUT / "gsea_top_hallmarks.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(10, max(4, 0.32 * len(g_top) + 1.5)))
    g_top_sorted = g_top.sort_values("NES")
    y = np.arange(len(g_top_sorted))
    colors = ["#d62728" if n > 0 else "#1f77b4" for n in g_top_sorted["NES"]]
    ax.barh(y, g_top_sorted["NES"], color=colors, edgecolor="black", alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['Term']}  FDR={r['FDR q-val']:.2g}" for _, r in g_top_sorted.iterrows()], fontsize=9)
    ax.axvline(0, color="#888", linestyle=":", linewidth=0.8)
    ax.set_xlabel("NES (positive = up in PTC+HT, negative = down)")
    ax.set_title("GSEA Hallmark — GSE286332 PTC+HT vs PTC (top 20 by |NES|, FDR<0.05)\nUp = immune/inflammation; down = fatty acid metabolism (thyroid identity)", fontsize=10)
    ax.grid(axis="x", linestyle=":", alpha=0.3)
    plt.tight_layout()
    fig.savefig(OUT / "gsea_top_hallmarks.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "gsea_top_hallmarks.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  top 20 GSEA pathways saved")
    return g_top


# ============== 5. Driver × DM stratified ==============
def driver_dm():
    print("\n[5] driver_anchor × DM")
    sm = pd.read_csv(SAMPLE_MASTER, sep="\t")
    sig = pd.read_csv(DM_SIG, sep="\t", index_col=0)
    sm["sample_short"] = sm["sample_id"].astype(str)
    sig["sample_short"] = sig.index.to_series().astype(str)
    merged = sig.merge(sm[["sample_short", "driver_anchor", "molecular_subtype", "ata_risk_proxy"]], on="sample_short", how="left")
    merged = merged.dropna(subset=["DM", "driver_anchor"])

    ct = pd.crosstab(merged["driver_anchor"], merged["DM"])
    ct["total"] = ct.sum(axis=1)
    ct["DM1_pct"] = ct.get("DM1", 0) / ct["total"] * 100
    ct.to_csv(OUT / "driver_dm.tsv", sep="\t")
    print(f"  driver × DM:\n{ct}")

    # Plot — exclude tiny groups (n<5)
    show = ct[ct["total"] >= 5].sort_values("DM1_pct", ascending=False)
    if len(show):
        fig, ax = plt.subplots(figsize=(8, 4.5))
        x = np.arange(len(show))
        ax.bar(x, show["DM1_pct"], color="#fdae61", edgecolor="black", alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{idx}\n(n={int(show.loc[idx, 'total'])})" for idx in show.index], fontsize=9)
        ax.axhline(28, color="#d62728", linestyle="--", linewidth=1, label="TCGA baseline DM1=28%")
        ax.axhline(50, color="#888", linestyle=":", linewidth=0.8, alpha=0.5)
        ax.set_ylabel("% DM1 within driver class")
        ax.set_title(f"DM1 enrichment by driver_anchor (TCGA-THCA, n_total={len(merged)})", fontsize=10)
        for i, idx in enumerate(show.index):
            ax.text(i, show.loc[idx, "DM1_pct"] + 2, f"{show.loc[idx, 'DM1_pct']:.0f}%", ha="center", fontsize=9)
        ax.set_ylim(0, 110)
        ax.legend(fontsize=8)
        ax.grid(axis="y", linestyle=":", alpha=0.3)
        plt.tight_layout()
        fig.savefig(OUT / "driver_dm.png", dpi=160, bbox_inches="tight")
        fig.savefig(OUT / "driver_dm.pdf", bbox_inches="tight")
        plt.close(fig)

    return ct


def main():
    res = {}
    res["pancan_significant_lineages"] = pancan_prognostic().to_dict(orient="records")
    res["fusion_enrichment"] = fusion_dm()
    res["k2_summary"] = k2_layer()
    res["gsea_top"] = gsea_layer().to_dict(orient="records")
    res["driver_dm"] = driver_dm().to_dict()
    res["generated_at"] = "2026-05-08 round5"

    with open(OUT / "round5_summary.json", "w") as f:
        json.dump(res, f, indent=2, default=str)

    print(f"\nwrote {OUT.relative_to(REPO)}/")
    print("Outputs:")
    for f in sorted(OUT.iterdir()):
        print(f"  {f.relative_to(OUT)}")


if __name__ == "__main__":
    main()
