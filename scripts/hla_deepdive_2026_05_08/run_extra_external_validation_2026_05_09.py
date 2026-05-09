#!/usr/bin/env python3
"""Build extra external-validation assets for the two HLA manuscripts.

The script keeps the Paper 2 boundary explicit: cancer-cohort HLA results are
gene-expression modules, not HLA allele associations.
"""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation"
TABLES = RESULTS / "tables"
FIGS = RESULTS / "figures"
REPORTS = ROOT / "project/reports"


TRACK = ROOT / "project/results/hla_deepdive_2026_05_08"


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)


def read_tsv(rel: str, **kwargs) -> pd.DataFrame:
    path = ROOT / rel
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", **kwargs)


def read_csv(rel: str, **kwargs) -> pd.DataFrame:
    path = ROOT / rel
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, **kwargs)


def cohen_d(group_a: pd.Series, group_b: pd.Series) -> float:
    a = pd.to_numeric(group_a, errors="coerce").dropna().to_numpy()
    b = pd.to_numeric(group_b, errors="coerce").dropna().to_numpy()
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    va = np.var(a, ddof=1)
    vb = np.var(b, ddof=1)
    pooled = ((len(a) - 1) * va + (len(b) - 1) * vb) / (len(a) + len(b) - 2)
    if pooled <= 0:
        return float("nan")
    return float((np.mean(a) - np.mean(b)) / math.sqrt(pooled))


def bh_fdr(pvals: list[float]) -> list[float]:
    arr = np.array([1.0 if pd.isna(p) else float(p) for p in pvals])
    order = np.argsort(arr)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(arr) + 1)
    q = arr * len(arr) / ranks
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    out = np.empty_like(q_sorted)
    out[order] = np.minimum(q_sorted, 1.0)
    return out.tolist()


def spearman_pair(df: pd.DataFrame, x: str, y: str) -> dict:
    sub = df[[x, y]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(sub) < 3:
        return {"var_x": x, "var_y": y, "n": len(sub), "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(sub[x], sub[y])
    return {"var_x": x, "var_y": y, "n": len(sub), "rho": float(rho), "p": float(p)}


def fmt_p(p: float) -> str:
    if pd.isna(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def run_gse286332_ht_overlap() -> pd.DataFrame:
    hla = read_tsv("project/results/p3_gse286332/hla_module_scores.tsv")
    tls = read_tsv("project/results/d5p6_bcr_repertoire/tls_score_per_sample.tsv")
    hla = hla.rename(columns={"Unnamed: 0": "sample"})
    tls = tls.rename(columns={"Unnamed: 0": "sample"})
    merged = hla.merge(tls[["sample", "TLS_score"]], on="sample", how="inner")
    merged.to_csv(TABLES / "T01a_gse286332_hla_tls_per_sample.tsv", sep="\t", index=False)

    rows = []
    for score in ["HLA_I", "HLA_II", "TLS_score"]:
        a = merged.loc[merged["group"] == "PTC_HT", score]
        b = merged.loc[merged["group"] == "PTC", score]
        mwu = stats.mannwhitneyu(a, b, alternative="two-sided", method="auto")
        rows.append(
            {
                "score": score,
                "external_dataset": "GSE286332",
                "contrast": "PTC_HT_vs_PTC",
                "n_PTC_HT": int(a.notna().sum()),
                "n_PTC": int(b.notna().sum()),
                "mean_PTC_HT": float(a.mean()),
                "mean_PTC": float(b.mean()),
                "median_PTC_HT": float(a.median()),
                "median_PTC": float(b.median()),
                "delta_mean": float(a.mean() - b.mean()),
                "cohen_d_PTC_HT_minus_PTC": cohen_d(a, b),
                "mannwhitney_p": float(mwu.pvalue),
                "boundary": "HLA expression module, not allele genotype",
            }
        )
    effects = pd.DataFrame(rows)
    effects["BH_FDR_3_tests"] = bh_fdr(effects["mannwhitney_p"].tolist())
    effects.to_csv(TABLES / "T01_gse286332_ht_overlap_hla_tls_effects.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    colors = ["#5B8FF9", "#61DDAA", "#F6BD16"]
    ax.bar(effects["score"], effects["cohen_d_PTC_HT_minus_PTC"], color=colors, edgecolor="black", linewidth=0.6)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Cohen's d (PTC_HT - PTC)")
    ax.set_title("GSE286332 HT-overlap validation")
    for i, row in effects.iterrows():
        ax.text(
            i,
            row["cohen_d_PTC_HT_minus_PTC"] + 0.05,
            f"p={fmt_p(row['mannwhitney_p'])}\nFDR={fmt_p(row['BH_FDR_3_tests'])}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    fig.tight_layout()
    fig.savefig(FIGS / "F01_gse286332_ht_overlap_hla_tls_effects.png", dpi=220)
    plt.close(fig)
    return effects


def run_tcga_network() -> tuple[pd.DataFrame, pd.DataFrame]:
    tcga = read_tsv(
        "project/results/hla_deepdive_2026_05_08/track26_ifng_hla_dm1/tables/T10_per_sample_combined.tsv"
    )
    tls = read_tsv("project/results/audit_2026_04_30/round10/r10_5_tcga_tls_per_sample.tsv")
    rai = read_tsv(
        "project/results/hla_deepdive_2026_05_08/track34_rai_hla_decouple/tables/T01_per_sample_scores.tsv"
    )
    tls = tls.rename(columns={"sample_short": "sample"})
    duplicate_audit = pd.DataFrame(
        [
            {
                "table": "track26_T10_per_sample_combined",
                "n_rows_before_dedupe": len(tcga),
                "n_unique_samples": tcga["sample"].nunique(),
                "n_duplicate_rows": len(tcga) - tcga["sample"].nunique(),
            },
            {
                "table": "audit_round10_TLS_per_sample",
                "n_rows_before_dedupe": len(tls),
                "n_unique_samples": tls["sample"].nunique(),
                "n_duplicate_rows": len(tls) - tls["sample"].nunique(),
            },
            {
                "table": "track34_RAI_per_sample_scores",
                "n_rows_before_dedupe": len(rai),
                "n_unique_samples": rai["sample"].nunique(),
                "n_duplicate_rows": len(rai) - rai["sample"].nunique(),
            },
        ]
    )
    duplicate_audit.to_csv(TABLES / "T02d_tcga_sample_duplicate_audit.tsv", sep="\t", index=False)
    tcga = tcga.drop_duplicates("sample", keep="first").copy()
    tls = tls.drop_duplicates("sample", keep="first").copy()
    rai = rai.drop_duplicates("sample", keep="first").copy()
    tcga = tcga.merge(tls[["sample", "TLS_score"]], on="sample", how="left")
    tcga = tcga.merge(rai[["sample", "RAI_module"]], on="sample", how="left")
    tcga.to_csv(TABLES / "T02a_tcga_hla_ifng_tls_rai_per_sample.tsv", sep="\t", index=False)

    variables = [
        "DM1_use",
        "HLA1_score",
        "HLA2_score",
        "IFNG_hallmark",
        "IFNG_ayers",
        "TLS_score",
        "RAI_module",
        "Immune_proxy",
        "Purity_proxy",
    ]
    labels = {
        "DM1_use": "DM1",
        "HLA1_score": "HLA-I",
        "HLA2_score": "HLA-II",
        "IFNG_hallmark": "IFNG Hallmark",
        "IFNG_ayers": "IFNG Ayers",
        "TLS_score": "TLS",
        "RAI_module": "RAI",
        "Immune_proxy": "Immune proxy",
        "Purity_proxy": "Purity proxy",
    }

    rows = []
    for i, x in enumerate(variables):
        for y in variables[i + 1 :]:
            row = spearman_pair(tcga, x, y)
            row["var_x_label"] = labels[x]
            row["var_y_label"] = labels[y]
            rows.append(row)
    network = pd.DataFrame(rows)
    network.to_csv(TABLES / "T02_tcga_hla_ifng_tls_rai_network.tsv", sep="\t", index=False)
    network.loc[(network["p"] < 1e-3) & (network["rho"].abs() >= 0.3)].sort_values(
        "rho", key=lambda s: s.abs(), ascending=False
    ).to_csv(TABLES / "T02b_tcga_strong_network_edges.tsv", sep="\t", index=False)

    mat = pd.DataFrame(np.eye(len(variables)), index=[labels[v] for v in variables], columns=[labels[v] for v in variables])
    for _, row in network.iterrows():
        mat.loc[row["var_x_label"], row["var_y_label"]] = row["rho"]
        mat.loc[row["var_y_label"], row["var_x_label"]] = row["rho"]
    mat.to_csv(TABLES / "T02c_tcga_hla_ifng_tls_rai_rho_matrix.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    im = ax.imshow(mat.values.astype(float), vmin=-1, vmax=1, cmap="coolwarm")
    ax.set_xticks(range(len(mat.columns)))
    ax.set_yticks(range(len(mat.index)))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right")
    ax.set_yticklabels(mat.index)
    for r in range(mat.shape[0]):
        for c in range(mat.shape[1]):
            ax.text(c, r, f"{mat.iloc[r, c]:.2f}", ha="center", va="center", fontsize=8)
    ax.set_title("TCGA-THCA HLA/IFNG/TLS/RAI module network")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Spearman rho")
    fig.tight_layout()
    fig.savefig(FIGS / "F02_tcga_hla_ifng_tls_rai_network.png", dpi=220)
    plt.close(fig)
    return network, mat


def run_scrna_concordance() -> tuple[pd.DataFrame, pd.DataFrame]:
    pu = read_csv(
        "project/results/hla_deepdive_2026_05_08/track12_scrna_hla_celltype/tables/T3_pu_HLA_module_per_celltype.csv"
    )
    lu = read_csv(
        "project/results/hla_deepdive_2026_05_08/track12_scrna_hla_celltype/tables/T3_lu_HLA_module_per_celltype.csv"
    )
    shared = pu.merge(lu, on="celltype", suffixes=("_pu", "_lu"))
    shared = shared[(shared["HLA_II_score_count_pu"] > 0) & (shared["HLA_II_score_count_lu"] > 0)].copy()
    shared = shared.dropna(subset=["HLA_II_score_mean_pu", "HLA_II_score_mean_lu"])
    shared["rank_HLA_II_pu"] = shared["HLA_II_score_mean_pu"].rank(ascending=False, method="min")
    shared["rank_HLA_II_lu"] = shared["HLA_II_score_mean_lu"].rank(ascending=False, method="min")
    shared["mean_HLA_II_across_scRNA"] = shared[["HLA_II_score_mean_pu", "HLA_II_score_mean_lu"]].mean(axis=1)
    shared.sort_values("mean_HLA_II_across_scRNA", ascending=False).to_csv(
        TABLES / "T03_scrna_celltype_hla_concordance.tsv", sep="\t", index=False
    )

    rho2, p2 = stats.spearmanr(shared["HLA_II_score_mean_pu"], shared["HLA_II_score_mean_lu"])
    rho1, p1 = stats.spearmanr(shared["HLA_I_score_mean_pu"], shared["HLA_I_score_mean_lu"])
    top3_pu = set(shared.sort_values("HLA_II_score_mean_pu", ascending=False).head(3)["celltype"])
    top3_lu = set(shared.sort_values("HLA_II_score_mean_lu", ascending=False).head(3)["celltype"])
    apc = ["DC", "Myeloid", "B_cell", "Plasma"]
    thy = shared.loc[shared["celltype"] == "Thyrocyte", "mean_HLA_II_across_scRNA"].iloc[0]
    apc_mean = shared.loc[shared["celltype"].isin(apc), "mean_HLA_II_across_scRNA"].mean()
    summary = pd.DataFrame(
        [
            {
                "metric": "HLA-II cell-type mean concordance, Pu vs Lu",
                "n_celltypes": len(shared),
                "spearman_rho": float(rho2),
                "p": float(p2),
                "interpretation": "cell-type source is externally reproducible",
            },
            {
                "metric": "HLA-I cell-type mean concordance, Pu vs Lu",
                "n_celltypes": len(shared),
                "spearman_rho": float(rho1),
                "p": float(p1),
                "interpretation": "HLA-I source ranking is also reproducible",
            },
            {
                "metric": "HLA-II top-3 overlap",
                "n_celltypes": 3,
                "spearman_rho": np.nan,
                "p": np.nan,
                "interpretation": ",".join(sorted(top3_pu & top3_lu)),
            },
            {
                "metric": "APC mean HLA-II vs thyrocyte HLA-II",
                "n_celltypes": len(apc) + 1,
                "spearman_rho": np.nan,
                "p": np.nan,
                "interpretation": f"APC_mean={apc_mean:.3f}; thyrocyte={thy:.3f}; ratio={(apc_mean / thy if thy else np.nan):.2f}",
            },
        ]
    )
    summary.to_csv(TABLES / "T03b_scrna_concordance_summary.tsv", sep="\t", index=False)

    plot_df = shared.sort_values("mean_HLA_II_across_scRNA", ascending=False)
    x = np.arange(len(plot_df))
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    ax.bar(x - 0.18, plot_df["HLA_II_score_mean_pu"], width=0.36, label="GSE184362 Pu", color="#5B8FF9")
    ax.bar(x + 0.18, plot_df["HLA_II_score_mean_lu"], width=0.36, label="GSE193581 Lu", color="#61DDAA")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df["celltype"], rotation=45, ha="right")
    ax.set_ylabel("Mean HLA-II module")
    ax.set_title(f"scRNA cell-type HLA-II concordance (rho={rho2:.2f})")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGS / "F03_scrna_hla2_celltype_concordance.png", dpi=220)
    plt.close(fig)
    return shared, summary


def run_spatial_tls() -> tuple[pd.DataFrame, pd.DataFrame]:
    sample = read_tsv("project/results/hla_deepdive_2026_05_08/track33_spatial_tls_hla/tables/T1_per_sample_summary.tsv")
    stage = read_tsv("project/results/hla_deepdive_2026_05_08/track33_spatial_tls_hla/tables/T6_stage_aggregate_effect.tsv")
    sample = sample.copy()
    sample["positive_HLA_II_in_TLS_niche"] = sample["diff_HLA_II_in_minus_out"] > 0
    sample["nominal_positive"] = (sample["diff_HLA_II_in_minus_out"] > 0) & (sample["ttest_p"] < 0.05)
    sample.to_csv(TABLES / "T04_spatial_tls_hla_per_sample_sign_consistency.tsv", sep="\t", index=False)
    stage.to_csv(TABLES / "T04b_spatial_tls_hla_stage_effect.tsv", sep="\t", index=False)

    summaries = []
    for label, sub in [("all_samples", sample), ("cancer_samples", sample[sample["stage"] != "N"])]:
        n = len(sub)
        n_pos = int(sub["positive_HLA_II_in_TLS_niche"].sum())
        n_sig = int(sub["nominal_positive"].sum())
        summaries.append(
            {
                "subset": label,
                "n_samples": n,
                "n_positive_diff": n_pos,
                "n_nominal_positive_p_lt_0_05": n_sig,
                "binom_p_positive_vs_0_5": float(stats.binomtest(n_pos, n, p=0.5, alternative="greater").pvalue),
                "mean_cohen_d": float(sub["cohen_d"].mean()),
                "mean_diff_HLA_II_in_minus_out": float(sub["diff_HLA_II_in_minus_out"].mean()),
            }
        )
    summary = pd.DataFrame(summaries)
    summary.to_csv(TABLES / "T04c_spatial_tls_hla_summary.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    ax.bar(stage["stage"], stage["mean_cohen_d"], yerr=stage["sem_cohen_d"], color="#E8684A", alpha=0.85)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Mean Cohen's d")
    ax.set_title("Spatial TLS niche HLA-II enrichment")
    fig.tight_layout()
    fig.savefig(FIGS / "F04_spatial_tls_hla_stage_effect.png", dpi=220)
    plt.close(fig)
    return sample, summary


def run_ici_context() -> pd.DataFrame:
    meta = read_tsv("project/results/hla_deepdive_2026_05_08/track18_ici_hla/tables/T3_meta_OR_response.tsv")
    out = meta.copy()
    out["external_role"] = np.where(
        out["score"].eq("HLA_class_I"),
        "external ICI response validation of antigen-presentation context",
        np.where(out["score"].eq("HLA_class_II"), "negative/control translational context", "immune-context comparator"),
    )
    out["boundary"] = "HLA expression module, not allele genotype"
    out.to_csv(TABLES / "T05_ici_hla_response_external_context.tsv", sep="\t", index=False)
    return out


def classify_thyroid_branch(tag: str) -> str:
    tag_u = tag.upper()
    if tag_u in {"GD", "GD-OPHT", "GD-OPHT-BROAD", "GD-EUR-EUR"}:
        return "strict_GD"
    if tag_u == "AI-HYPER":
        return "autoimmune_hyperthyroid"
    if any(k in tag_u for k in ["SELF-HYPER", "E05", "CARBIM", "THYROTOX"]):
        return "broad_hyperthyroid"
    if any(k in tag_u for k in ["HT", "HYPOTHY", "LEVO", "E03", "CHRONIC", "THYROIDITIS", "E06"]):
        return "HT_or_hypothyroid"
    if any(k in tag_u for k in ["PSOR", "RA", "T1D"]):
        return "non_thyroid_autoimmune"
    return "broad_thyroid"


def run_paper4_external() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gwas = read_tsv("project/results/hla_deepdive_2026_05_08/track7_gwas_catalog_mhc/tables/T08_mhc_share_per_trait.tsv")
    gwas = gwas.sort_values("MHC_GWS_share", ascending=False)
    gwas.to_csv(TABLES / "T06_gwas_catalog_mhc_share_external.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    ax.bar(gwas["trait_label"], gwas["MHC_GWS_share"], color="#5D7092")
    ax.set_ylabel("Genome-wide significant share in MHC")
    ax.set_title("GWAS Catalog thyroid autoimmunity MHC concentration")
    ax.set_xticks(range(len(gwas)))
    ax.set_xticklabels(gwas["trait_label"], rotation=35, ha="right")
    for i, row in enumerate(gwas.itertuples(index=False)):
        ax.text(i, row.MHC_GWS_share + 0.006, f"{row.n_GWS_in_MHC}/{row.n_GWS_total}", ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "F05_paper4_gwas_mhc_share.png", dpi=220)
    plt.close(fig)

    tags = read_tsv("project/results/hla_deepdive_2026_05_08/track11_finngen_panukbb/tables/T05_track1_tagsnp_transancestry.tsv")
    tags["thyroid_branch"] = tags["tag"].map(classify_thyroid_branch)
    tags["direction"] = np.where(tags["beta"] > 0, "positive", "negative")
    top_tags = (
        tags.sort_values("pval")
        .groupby(["allele", "thyroid_branch"], as_index=False)
        .head(3)
        .sort_values(["allele", "thyroid_branch", "pval"])
    )
    top_tags.to_csv(TABLES / "T07_track1_tag_snp_top_external_signals.tsv", sep="\t", index=False)

    sig = tags[tags["pval"] < 5e-8].copy()
    branch_summary = (
        sig.groupby(["allele", "thyroid_branch", "direction"], dropna=False)
        .agg(n_gws=("pval", "size"), min_p=("pval", "min"), median_beta=("beta", "median"))
        .reset_index()
        .sort_values(["allele", "thyroid_branch", "direction"])
    )
    branch_summary.to_csv(TABLES / "T07b_track1_tag_snp_branch_direction_summary.tsv", sep="\t", index=False)

    threshold_rows = []
    for threshold_label, threshold in [
        ("gws_5e-8", 5e-8),
        ("strong_1e-6", 1e-6),
        ("suggestive_1e-5", 1e-5),
        ("nominal_1e-4", 1e-4),
        ("nominal_0.05", 0.05),
    ]:
        sub = tags[tags["pval"] < threshold]
        grouped = (
            sub.groupby(["allele", "thyroid_branch", "direction"], dropna=False)
            .agg(n_signals=("pval", "size"), min_p=("pval", "min"), median_beta=("beta", "median"))
            .reset_index()
        )
        grouped["threshold_label"] = threshold_label
        grouped["p_threshold"] = threshold
        threshold_rows.append(grouped)
    threshold_summary = pd.concat(threshold_rows, ignore_index=True).sort_values(
        ["allele", "thyroid_branch", "direction", "p_threshold"]
    )
    threshold_summary.to_csv(TABLES / "T07c_track1_tag_snp_branch_threshold_summary.tsv", sep="\t", index=False)

    afnd = read_tsv("project/results/hla_deepdive_2026_05_08/track8_afnd_extended/tables/track1_korean_baseline.tsv")
    tri = read_tsv("project/results/hla_deepdive_2026_05_08/track8_afnd_extended/tables/korean_triangulation.tsv")
    focus = ["DPB1*05:01", "B*46:01", "A*02:07", "C*01:02", "DRB1*07:01", "DQB1*02:01"]
    afnd_rows = []
    for allele in focus:
        x = afnd[afnd["allele"] == allele].copy()
        tri_x = tri[tri["allele"] == allele]
        row = {
            "allele": allele,
            "AFND_available": bool(len(x)),
            "n_countries": int(x["country"].nunique()) if len(x) else 0,
            "korea_freq": np.nan,
            "korea_n": np.nan,
            "country_min_freq": np.nan,
            "country_max_freq": np.nan,
            "country_of_max": "",
            "korea_rank_high_to_low": np.nan,
            "eastasia_pool_freq": np.nan,
            "interpretation": "not available in Track 8 AFND focus extraction",
        }
        if len(x):
            sx = x.sort_values("weighted_mean_freq", ascending=False).reset_index(drop=True)
            korea = sx[sx["country"] == "South Korea"]
            row.update(
                {
                    "korea_freq": float(korea["weighted_mean_freq"].iloc[0]) if len(korea) else np.nan,
                    "korea_n": int(korea["n_pooled_individuals"].iloc[0]) if len(korea) else np.nan,
                    "country_min_freq": float(sx["weighted_mean_freq"].min()),
                    "country_max_freq": float(sx["weighted_mean_freq"].max()),
                    "country_of_max": str(sx.loc[sx["weighted_mean_freq"].idxmax(), "country"]),
                    "korea_rank_high_to_low": int(korea.index[0] + 1) if len(korea) else np.nan,
                    "interpretation": "population-frequency context available; not a disease association",
                }
            )
        if len(tri_x):
            row["eastasia_pool_freq"] = float(tri_x["EastAsia_pool_freq"].iloc[0])
        afnd_rows.append(row)
    afnd_context = pd.DataFrame(afnd_rows)
    afnd_context.to_csv(TABLES / "T08_afnd_korean_baseline_context.tsv", sep="\t", index=False)

    return gwas, top_tags, branch_summary, threshold_summary, afnd_context


def run_external_st_secondary() -> pd.DataFrame:
    path = ROOT / "project_external_st/results/meta/sample_level_score_summary.tsv"
    if not path.exists():
        out = pd.DataFrame(
            [
                {
                    "status": "missing",
                    "note": "project_external_st sample summary not present",
                }
            ]
        )
        out.to_csv(TABLES / "T09_external_st_dm1_rai_secondary.tsv", sep="\t", index=False)
        return out

    st = pd.read_csv(path, sep="\t")
    rows = []
    for dataset, sub in st.groupby("dataset"):
        for score in [
            "mean_RAI_8_score_raw_epi50",
            "mean_DM1_like_score_raw_epi50",
            "mean_TDS_overlap_score_raw_epi50",
            "mean_CAF_ECM_score_raw_epi50",
        ]:
            for condition, csub in sub.groupby("condition"):
                rows.append(
                    {
                        "dataset": dataset,
                        "condition": condition,
                        "score": score,
                        "n_slides": int(csub["sample_id"].nunique()),
                        "mean": float(csub[score].mean()),
                        "median": float(csub[score].median()),
                        "sd": float(csub[score].std(ddof=1)) if len(csub) > 1 else np.nan,
                        "boundary": "secondary DM1/RAI spatial validation; no HLA module in this score table",
                    }
                )
    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "T09_external_st_dm1_rai_secondary.tsv", sep="\t", index=False)
    return out


def build_validation_ladder(
    gse: pd.DataFrame,
    network: pd.DataFrame,
    scrna_summary: pd.DataFrame,
    spatial_summary: pd.DataFrame,
    ici: pd.DataFrame,
    gwas: pd.DataFrame,
    branch_summary: pd.DataFrame,
    threshold_summary: pd.DataFrame,
    afnd_context: pd.DataFrame,
) -> pd.DataFrame:
    def get_rho(x: str, y: str) -> float:
        hit = network[
            ((network["var_x"] == x) & (network["var_y"] == y)) | ((network["var_x"] == y) & (network["var_y"] == x))
        ]
        return float(hit["rho"].iloc[0]) if len(hit) else np.nan

    hla2_effect = gse.loc[gse["score"] == "HLA_II"].iloc[0]
    tls_effect = gse.loc[gse["score"] == "TLS_score"].iloc[0]
    hla_i_ici = ici.loc[ici["score"] == "HLA_class_I"].iloc[0]
    hla_ii_ici = ici.loc[ici["score"] == "HLA_class_II"].iloc[0]
    gd = gwas.loc[gwas["trait_label"] == "Graves disease"].iloc[0]
    hypo = gwas.loc[gwas["trait_label"] == "Hypothyroidism"].iloc[0]
    dpb1_hyper_strong = threshold_summary[
        (threshold_summary["allele"] == "HLA-DPB1*05:01")
        & (threshold_summary["thyroid_branch"] == "strict_GD")
        & (threshold_summary["direction"] == "positive")
        & (threshold_summary["threshold_label"] == "strong_1e-6")
    ]["n_signals"].sum()
    dpb1_hyper_nominal_1e4 = threshold_summary[
        (threshold_summary["allele"] == "HLA-DPB1*05:01")
        & (threshold_summary["thyroid_branch"] == "strict_GD")
        & (threshold_summary["direction"] == "positive")
        & (threshold_summary["threshold_label"] == "nominal_1e-4")
    ]["n_signals"].sum()
    dpb1_broad_hyper_negative_gws = branch_summary[
        (branch_summary["allele"] == "HLA-DPB1*05:01")
        & (branch_summary["thyroid_branch"] == "broad_hyperthyroid")
        & (branch_summary["direction"] == "negative")
    ]["n_gws"].sum()
    dpb1_hypo_neg = branch_summary[
        (branch_summary["allele"] == "HLA-DPB1*05:01")
        & (branch_summary["thyroid_branch"] == "HT_or_hypothyroid")
        & (branch_summary["direction"] == "negative")
    ]["n_gws"].sum()
    dpb1_afnd = afnd_context.loc[afnd_context["allele"] == "DPB1*05:01"].iloc[0]

    rows = [
        {
            "paper": "Paper 4",
            "hypothesis": "DPB1*05:01 is a Pan-Asian GD susceptibility anchor",
            "new_external_validation": f"FinnGen/PanUKBB tag SNP: {int(dpb1_hyper_strong)} strict-GD positive signal at p<1e-6 and {int(dpb1_hyper_nominal_1e4)} at p<1e-4; AFND Korea freq {dpb1_afnd['korea_freq']:.3f}",
            "support_score_1_to_5": 5,
            "manuscript_action": "Main genetic architecture figure + external validation panel",
            "boundary": "allele-level claim allowed only in non-cancer GD/AITD cohorts",
        },
        {
            "paper": "Paper 4",
            "hypothesis": "GD and HT/hypothyroid branches have divergent MHC directionality",
            "new_external_validation": f"DPB1*05:01 tag has {int(dpb1_hypo_neg)} genome-wide negative HT/hypothyroid signals and {int(dpb1_broad_hyper_negative_gws)} broad-hyperthyroid negative signal, indicating endpoint sensitivity",
            "support_score_1_to_5": 4,
            "manuscript_action": "Add disease-branch divergence hypothesis as Discussion + Supplementary forest",
            "boundary": "tag-SNP proxy, not direct HLA imputation",
        },
        {
            "paper": "Paper 4",
            "hypothesis": "Thyroid autoimmunity has an MHC-concentrated architecture",
            "new_external_validation": f"GWAS Catalog: GD MHC GWS share {gd['MHC_GWS_share']:.3f}; hypothyroidism {hypo['MHC_GWS_share']:.3f}",
            "support_score_1_to_5": 4,
            "manuscript_action": "Use as independent public-genetics support for MHC-centered framing",
            "boundary": "locus-level external support, not allele-specific proof",
        },
        {
            "paper": "Paper 2",
            "hypothesis": "HT-overlap PTC has a higher antigen-presentation/TLS state",
            "new_external_validation": f"GSE286332: HLA-II d={hla2_effect['cohen_d_PTC_HT_minus_PTC']:.2f}, FDR={hla2_effect['BH_FDR_3_tests']:.3g}; TLS d={tls_effect['cohen_d_PTC_HT_minus_PTC']:.2f}",
            "support_score_1_to_5": 4,
            "manuscript_action": "Move from exploratory aside to main external HT-overlap validation panel",
            "boundary": "HLA expression module only",
        },
        {
            "paper": "Paper 2",
            "hypothesis": "IFN-gamma is the gateway linking DM1 to HLA induction",
            "new_external_validation": f"TCGA network: IFNG-HLA-I rho={get_rho('IFNG_hallmark','HLA1_score'):.2f}; IFNG-HLA-II rho={get_rho('IFNG_hallmark','HLA2_score'):.2f}",
            "support_score_1_to_5": 5,
            "manuscript_action": "Main mechanism schematic; observational mediation caveat retained",
            "boundary": "transcriptomic module; not interventional causality",
        },
        {
            "paper": "Paper 2",
            "hypothesis": "HLA-II signal is APC/TLS-niche dominated, not simply tumor-cell intrinsic",
            "new_external_validation": f"scRNA Pu/Lu HLA-II cell-type concordance rho={scrna_summary.iloc[0]['spearman_rho']:.2f}; spatial TLS all-sample positive {int(spatial_summary.iloc[0]['n_positive_diff'])}/{int(spatial_summary.iloc[0]['n_samples'])}",
            "support_score_1_to_5": 5,
            "manuscript_action": "Add single-cell and spatial localization as a main claim",
            "boundary": "cell-state localization only",
        },
        {
            "paper": "Paper 2",
            "hypothesis": "HLA-I, but not HLA-II, carries pan-cancer ICI response signal",
            "new_external_validation": f"ICI meta: HLA-I OR={hla_i_ici['OR']:.2f}, p={fmt_p(hla_i_ici['p'])}; HLA-II OR={hla_ii_ici['OR']:.2f}, p={fmt_p(hla_ii_ici['p'])}",
            "support_score_1_to_5": 3,
            "manuscript_action": "Use as translational context, not as thyroid-specific efficacy claim",
            "boundary": "pan-cancer ICI context; not THCA treatment prediction",
        },
    ]
    ladder = pd.DataFrame(rows)
    ladder.to_csv(TABLES / "T10_external_validation_ladder.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    colors = ladder["paper"].map({"Paper 4": "#5B8FF9", "Paper 2": "#E8684A"})
    y = np.arange(len(ladder))
    ax.barh(y, ladder["support_score_1_to_5"], color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels([h[:58] + ("..." if len(h) > 58 else "") for h in ladder["hypothesis"]], fontsize=8)
    ax.set_xlim(0, 5)
    ax.set_xlabel("External-validation support score")
    ax.set_title("Two-paper external validation ladder")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIGS / "F06_external_validation_ladder.png", dpi=220)
    plt.close(fig)
    return ladder


def write_report(
    gse: pd.DataFrame,
    network: pd.DataFrame,
    scrna_summary: pd.DataFrame,
    spatial_summary: pd.DataFrame,
    ici: pd.DataFrame,
    gwas: pd.DataFrame,
    branch_summary: pd.DataFrame,
    threshold_summary: pd.DataFrame,
    afnd_context: pd.DataFrame,
    ladder: pd.DataFrame,
    st_secondary: pd.DataFrame,
) -> Path:
    hla2 = gse.loc[gse["score"] == "HLA_II"].iloc[0]
    hlai = gse.loc[gse["score"] == "HLA_I"].iloc[0]
    tls = gse.loc[gse["score"] == "TLS_score"].iloc[0]

    def net_line(a: str, b: str) -> str:
        hit = network[
            ((network["var_x"] == a) & (network["var_y"] == b)) | ((network["var_x"] == b) & (network["var_y"] == a))
        ].iloc[0]
        return f"rho={hit['rho']:.3f}, p={fmt_p(hit['p'])}, n={int(hit['n'])}"

    hla_i_ici = ici.loc[ici["score"] == "HLA_class_I"].iloc[0]
    hla_ii_ici = ici.loc[ici["score"] == "HLA_class_II"].iloc[0]
    gd = gwas.loc[gwas["trait_label"] == "Graves disease"].iloc[0]
    hyp = gwas.loc[gwas["trait_label"] == "Hypothyroidism"].iloc[0]
    dpb1 = afnd_context.loc[afnd_context["allele"] == "DPB1*05:01"].iloc[0]
    dpb1_branch = branch_summary[branch_summary["allele"] == "HLA-DPB1*05:01"]
    dpb1_hyper_strong = int(
        threshold_summary[
            (threshold_summary["allele"] == "HLA-DPB1*05:01")
            & (threshold_summary["thyroid_branch"] == "strict_GD")
            & (threshold_summary["direction"] == "positive")
            & (threshold_summary["threshold_label"] == "strong_1e-6")
        ]["n_signals"].sum()
    )
    dpb1_hyper_nominal_1e4 = int(
        threshold_summary[
            (threshold_summary["allele"] == "HLA-DPB1*05:01")
            & (threshold_summary["thyroid_branch"] == "strict_GD")
            & (threshold_summary["direction"] == "positive")
            & (threshold_summary["threshold_label"] == "nominal_1e-4")
        ]["n_signals"].sum()
    )
    dpb1_broad_hyper_negative_gws = int(
        dpb1_branch[
            (dpb1_branch["thyroid_branch"] == "broad_hyperthyroid") & (dpb1_branch["direction"] == "negative")
        ]["n_gws"].sum()
    )
    dpb1_hypo_neg = int(
        dpb1_branch[
            (dpb1_branch["thyroid_branch"] == "HT_or_hypothyroid") & (dpb1_branch["direction"] == "negative")
        ]["n_gws"].sum()
    )

    lines = [
        "# HLA Two-Paper Extra External Validation Pack (2026-05-09)",
        "",
        "## Executive Upgrade",
        "",
        "이번 추가 분석은 두 논문을 같은 HLA라는 단어로 억지 연결하지 않고, 서로 다른 증거축으로 수준을 올리는 구조입니다.",
        "Paper 4는 비암 GD/AITD germline HLA allele 논문이고, Paper 2는 암 코호트의 HLA gene-expression module/immune-context 논문입니다.",
        "",
        "## Paper 2: HT-PTC Immune/HLA Module",
        "",
        f"- GSE286332 HT-overlap 재검정: HLA-I d={hlai['cohen_d_PTC_HT_minus_PTC']:.2f}, p={fmt_p(hlai['mannwhitney_p'])}; HLA-II d={hla2['cohen_d_PTC_HT_minus_PTC']:.2f}, p={fmt_p(hla2['mannwhitney_p'])}, FDR={fmt_p(hla2['BH_FDR_3_tests'])}; TLS d={tls['cohen_d_PTC_HT_minus_PTC']:.2f}, p={fmt_p(tls['mannwhitney_p'])}.",
        f"- TCGA-THCA module network: DM1-HLA-II {net_line('DM1_use','HLA2_score')}; IFNG-HLA-II {net_line('IFNG_hallmark','HLA2_score')}; TLS-HLA-II {net_line('TLS_score','HLA2_score')}; RAI-HLA-II {net_line('RAI_module','HLA2_score')}.",
        f"- scRNA 외부 검증: Pu/Lu 두 scRNA 코호트의 HLA-II cell-type rank concordance rho={scrna_summary.iloc[0]['spearman_rho']:.2f}, p={fmt_p(scrna_summary.iloc[0]['p'])}; top producers는 {scrna_summary.iloc[2]['interpretation']}.",
        f"- Spatial validation: GSE250521 TLS niche 안쪽 HLA-II가 바깥보다 높은 샘플 {int(spatial_summary.iloc[0]['n_positive_diff'])}/{int(spatial_summary.iloc[0]['n_samples'])}, binomial p={fmt_p(spatial_summary.iloc[0]['binom_p_positive_vs_0_5'])}; cancer subset mean d={spatial_summary.iloc[1]['mean_cohen_d']:.2f}.",
        f"- ICI pan-cancer context: HLA-I response OR={hla_i_ici['OR']:.2f}, p={fmt_p(hla_i_ici['p'])}; HLA-II OR={hla_ii_ici['OR']:.2f}, p={fmt_p(hla_ii_ici['p'])}. 이 결과는 HLA-I translational context는 살리고, HLA-II는 thyroid/TLS biology 중심으로 제한하는 데 쓰는 것이 안전합니다.",
        "",
        "## Paper 4: GD/AITD HLA Allele Genetics",
        "",
        f"- GWAS Catalog 외부 검증: GD genome-wide significant signals 중 MHC share={gd['MHC_GWS_share']:.3f} ({int(gd['n_GWS_in_MHC'])}/{int(gd['n_GWS_total'])}); hypothyroidism MHC share={hyp['MHC_GWS_share']:.3f} ({int(hyp['n_GWS_in_MHC'])}/{int(hyp['n_GWS_total'])}).",
        f"- FinnGen/PanUKBB tag-SNP 검증: DPB1*05:01 proxy는 strict Graves branch에서 positive signal p<1e-6 {dpb1_hyper_strong}개, p<1e-4 {dpb1_hyper_nominal_1e4}개이고, HT/hypothyroid branch에서 genome-wide negative signal {dpb1_hypo_neg}개입니다. Broad hyperthyroid endpoint에서는 negative genome-wide signal {dpb1_broad_hyper_negative_gws}개가 있어 endpoint sensitivity panel로 분리해야 합니다. 직접 HLA imputation이 아니라 tag-SNP proxy로 표기해야 합니다.",
        f"- AFND population context: DPB1*05:01 South Korea weighted frequency={dpb1['korea_freq']:.3f}, East-Asia pool={dpb1['eastasia_pool_freq']:.3f}; disease association이 아니라 background-frequency calibration으로만 사용합니다.",
        "",
        "## Manuscript-Level Actions",
        "",
        "- Paper 4 main Figure: Pan-Asian forest + GWAS MHC concentration + FinnGen/PanUKBB DPB1 branch divergence + AFND baseline context.",
        "- Paper 2 main Figure: HT-overlap HLA/TLS effect + TCGA IFNG-HLA network + scRNA APC source + spatial TLS-HLA niche + ICI context panel.",
        "- Paper 2 allele 후보(DRB1*04:05 등)는 n=9 vs 9에서 FDR-negative이므로 main claim이 아니라 power/future cohort box로 내려야 합니다.",
        "",
        "## Output Index",
        "",
        f"- Tables: `{TABLES.relative_to(ROOT)}`",
        f"- Figures: `{FIGS.relative_to(ROOT)}`",
        f"- Validation ladder: `{(TABLES / 'T10_external_validation_ladder.tsv').relative_to(ROOT)}`",
        "",
        "## Boundary",
        "",
        "Cancer-cohort HLA는 gene-expression module입니다. HLA allele association은 Paper 4의 비암 GD/AITD 유전학에서만 주장합니다.",
    ]

    report = REPORTS / "2026_05_09_HLA_EXTRA_EXTERNAL_VALIDATION_KR.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    shutil.copy2(report, RESULTS / "HLA_EXTRA_EXTERNAL_VALIDATION_KR.md")
    return report


def main() -> None:
    ensure_dirs()
    gse = run_gse286332_ht_overlap()
    network, _ = run_tcga_network()
    _, scrna_summary = run_scrna_concordance()
    _, spatial_summary = run_spatial_tls()
    ici = run_ici_context()
    gwas, top_tags, branch_summary, threshold_summary, afnd_context = run_paper4_external()
    st_secondary = run_external_st_secondary()
    ladder = build_validation_ladder(
        gse, network, scrna_summary, spatial_summary, ici, gwas, branch_summary, threshold_summary, afnd_context
    )
    report = write_report(
        gse,
        network,
        scrna_summary,
        spatial_summary,
        ici,
        gwas,
        branch_summary,
        threshold_summary,
        afnd_context,
        ladder,
        st_secondary,
    )

    manifest = {
        "results_dir": str(RESULTS),
        "tables_dir": str(TABLES),
        "figures_dir": str(FIGS),
        "report": str(report),
        "boundary": "Paper 2 uses HLA expression modules only; no cancer-cohort allele association claims.",
        "n_tables": len(list(TABLES.glob("*.tsv"))),
        "n_figures": len(list(FIGS.glob("*.png"))),
    }
    (RESULTS / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
