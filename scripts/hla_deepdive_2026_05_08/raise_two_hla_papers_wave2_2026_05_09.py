#!/usr/bin/env python3
"""Second level-up wave for the two HLA manuscripts.

Adds reviewer-facing robustness layers:
1) Paper 4 allele triangulation and endpoint robustness.
2) Paper 2 exact small-n permutation/bootstrap and TCGA HLA-II dominance.

Paper 2 remains HLA expression-module only in cancer cohorts.
"""

from __future__ import annotations

import json
import math
import shutil
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "project/results/hla_two_paper_synthesis_2026_05_09"
OUT = BASE / "level_up_wave2"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
REPORT = ROOT / "project/reports/2026_05_09_HLA_LEVEL_UP_WAVE2_KR.md"


def ensure_dirs() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)


def read_tsv(rel: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / rel, sep="\t")


def read_csv(rel: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / rel)


def cohen_d(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna().to_numpy()
    b = pd.to_numeric(b, errors="coerce").dropna().to_numpy()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = ((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) / (len(a) + len(b) - 2)
    return float((np.mean(a) - np.mean(b)) / math.sqrt(pooled)) if pooled > 0 else np.nan


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray([1.0 if pd.isna(x) else float(x) for x in pvals])
    order = np.argsort(p)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(p) + 1)
    q = p * len(p) / ranks
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    out = np.empty_like(q_sorted)
    out[order] = np.minimum(q_sorted, 1.0)
    return out.tolist()


def classify_branch(tag: str) -> str:
    tag_u = str(tag).upper()
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


def signed_logp(beta: float, pval: float, cap: float = 50.0) -> float:
    if pd.isna(beta) or pd.isna(pval) or pval <= 0:
        return np.sign(beta) * cap if not pd.isna(beta) else np.nan
    return float(np.sign(beta) * min(-math.log10(pval), cap))


def minmax(series: pd.Series, fill: float = 0.0) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    if s.notna().sum() == 0:
        return pd.Series(fill, index=s.index)
    lo, hi = s.min(), s.max()
    if hi == lo:
        return pd.Series(1.0, index=s.index)
    return (s - lo) / (hi - lo)


def exact_label_permutation(values: np.ndarray, labels: np.ndarray) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels).astype(bool)
    n = len(values)
    n_pos = int(labels.sum())
    obs = float(values[labels].mean() - values[~labels].mean())
    abs_obs = abs(obs)
    extreme = 0
    total = 0
    all_idx = range(n)
    for pos_idx in combinations(all_idx, n_pos):
        mask = np.zeros(n, dtype=bool)
        mask[list(pos_idx)] = True
        diff = float(values[mask].mean() - values[~mask].mean())
        if abs(diff) >= abs_obs - 1e-12:
            extreme += 1
        total += 1
    return obs, (extreme + 1) / (total + 1), total


def bootstrap_d(a: np.ndarray, b: np.ndarray, n_boot: int = 10000, seed: int = 20260509) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        aa = rng.choice(a, size=len(a), replace=True)
        bb = rng.choice(b, size=len(b), replace=True)
        vals.append(cohen_d(pd.Series(aa), pd.Series(bb)))
    return float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def paper4_triangulation() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = read_tsv("project/results/hla_deepdive_2026_05_08/track10_korean_lit/tables/T05_panasian_GD_DL_random_effects_v3.tsv")
    spec = read_tsv("project/results/hla_deepdive_2026_05_08/track31_cross_autoimmune/tables/T05_specificity_scores.tsv")
    tags = read_tsv("project/results/hla_deepdive_2026_05_08/track11_finngen_panukbb/tables/T05_track1_tagsnp_transancestry.tsv")
    afnd = read_tsv("project/results/hla_deepdive_2026_05_08/track8_afnd_extended/tables/track1_korean_baseline.tsv")

    tags = tags.copy()
    tags["branch"] = tags["tag"].map(classify_branch)
    tags["signed_log10p"] = tags.apply(lambda r: signed_logp(r["beta"], r["pval"], cap=50.0), axis=1)
    tags["allele_clean"] = tags["allele"].str.replace("HLA-", "", regex=False)
    branch = (
        tags.sort_values("pval")
        .groupby(["allele_clean", "branch"], as_index=False)
        .first()[["allele_clean", "branch", "tag", "beta", "pval", "signed_log10p", "ancestry"]]
    )
    branch.to_csv(TABLES / "T02_paper4_endpoint_best_tag_signal.tsv", sep="\t", index=False)

    rows = []
    for _, r in meta.iterrows():
        allele = r["allele"]
        sp = spec[spec["allele"] == allele]
        af = afnd[afnd["allele"] == allele]
        br = branch[branch["allele_clean"] == allele]
        strict = br[br["branch"] == "strict_GD"].sort_values("pval").head(1)
        ht = br[br["branch"] == "HT_or_hypothyroid"].sort_values("pval").head(1)
        broad_h = br[br["branch"] == "broad_hyperthyroid"].sort_values("pval").head(1)
        rows.append(
            {
                "allele": allele,
                "k_meta": int(r["k"]),
                "pooled_or": float(r["pooled_or"]),
                "ci_lo": float(r["ci_lo"]),
                "ci_hi": float(r["ci_hi"]),
                "p_random": float(r["p_random"]),
                "I2_pct": float(r["I2_pct"]),
                "meta_abs_logOR": abs(math.log(float(r["pooled_or"]))),
                "meta_strength_neglog10p": min(-math.log10(max(float(r["p_random"]), 1e-300)), 50.0),
                "heterogeneity_penalty": min(float(r["I2_pct"]) / 100.0, 1.0),
                "GD_specificity_score": float(sp["GD_specificity_score"].iloc[0]) if len(sp) else np.nan,
                "strict_GD_best_beta": float(strict["beta"].iloc[0]) if len(strict) else np.nan,
                "strict_GD_best_p": float(strict["pval"].iloc[0]) if len(strict) else np.nan,
                "strict_GD_signed_log10p": float(strict["signed_log10p"].iloc[0]) if len(strict) else np.nan,
                "HT_hypothyroid_best_beta": float(ht["beta"].iloc[0]) if len(ht) else np.nan,
                "HT_hypothyroid_best_p": float(ht["pval"].iloc[0]) if len(ht) else np.nan,
                "HT_hypothyroid_signed_log10p": float(ht["signed_log10p"].iloc[0]) if len(ht) else np.nan,
                "broad_hyperthyroid_best_beta": float(broad_h["beta"].iloc[0]) if len(broad_h) else np.nan,
                "broad_hyperthyroid_best_p": float(broad_h["pval"].iloc[0]) if len(broad_h) else np.nan,
                "afnd_korea_freq": float(af.loc[af["country"] == "South Korea", "weighted_mean_freq"].iloc[0])
                if len(af[af["country"] == "South Korea"])
                else np.nan,
                "afnd_n_countries": int(af["country"].nunique()) if len(af) else 0,
            }
        )
    tri = pd.DataFrame(rows)
    tri["score_meta"] = minmax(tri["meta_strength_neglog10p"])
    tri["score_specificity"] = minmax(tri["GD_specificity_score"]).fillna(0.0)
    tri["score_external_strict_GD"] = minmax(tri["strict_GD_signed_log10p"].clip(lower=0)).fillna(0.0)
    tri["score_population_context"] = (tri["afnd_n_countries"] > 0).astype(float)
    tri["score_low_heterogeneity"] = 1.0 - tri["heterogeneity_penalty"].fillna(1.0)
    tri["triangulation_score_100"] = (
        30 * tri["score_meta"]
        + 25 * tri["score_specificity"]
        + 20 * tri["score_external_strict_GD"]
        + 15 * tri["score_low_heterogeneity"]
        + 10 * tri["score_population_context"]
    ).round(1)
    tri["paper4_role"] = np.where(
        tri["allele"].eq("C*01:02"),
        "GD-specificity lead",
        np.where(tri["allele"].eq("DPB1*05:01"), "AITD-broad class-II anchor", "supporting architecture allele"),
    )
    tri = tri.sort_values("triangulation_score_100", ascending=False)
    tri.to_csv(TABLES / "T01_paper4_allele_triangulation_score.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    x = np.log(tri["pooled_or"].astype(float))
    y = tri["GD_specificity_score"].fillna(0)
    sizes = 80 + 4.0 * tri["triangulation_score_100"]
    colors = tri["I2_pct"]
    sc = ax.scatter(x, y, s=sizes, c=colors, cmap="viridis_r", edgecolor="black", linewidth=0.7)
    ax.axvline(0, color="black", linewidth=0.8)
    for _, row in tri.iterrows():
        ax.text(math.log(row["pooled_or"]), 0 if pd.isna(row["GD_specificity_score"]) else row["GD_specificity_score"] + 0.25, row["allele"], fontsize=9, ha="center")
    ax.set_xlabel("log pooled OR in Pan-Asian GD meta-analysis")
    ax.set_ylabel("GD specificity score")
    ax.set_title("Paper 4 allele triangulation map")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("I2 heterogeneity (%)")
    fig.tight_layout()
    fig.savefig(FIGS / "F01_paper4_allele_triangulation_map.png", dpi=220)
    plt.close(fig)

    heat = branch[branch["allele_clean"].isin(tri["allele"])].pivot_table(
        index="allele_clean", columns="branch", values="signed_log10p", aggfunc="first"
    )
    wanted = ["strict_GD", "autoimmune_hyperthyroid", "broad_hyperthyroid", "HT_or_hypothyroid", "broad_thyroid", "non_thyroid_autoimmune"]
    heat = heat.reindex(index=tri["allele"], columns=[c for c in wanted if c in heat.columns])
    heat.to_csv(TABLES / "T02b_paper4_endpoint_signed_log10p_matrix.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(8.8, 4.6))
    vals = heat.fillna(0).to_numpy(dtype=float)
    im = ax.imshow(vals, cmap="coolwarm", vmin=-12, vmax=12, aspect="auto")
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=35, ha="right")
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            text = "" if pd.isna(heat.iloc[i, j]) else f"{heat.iloc[i, j]:.1f}"
            ax.text(j, i, text, ha="center", va="center", fontsize=8)
    ax.set_title("Endpoint robustness and phenotype sensitivity")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("signed -log10(p) at best tag SNP")
    fig.tight_layout()
    fig.savefig(FIGS / "F02_paper4_endpoint_robustness_heatmap.png", dpi=220)
    plt.close(fig)
    return tri, branch


def gse286332_exact_robustness() -> tuple[pd.DataFrame, pd.DataFrame]:
    hla = read_tsv("project/results/p3_gse286332/hla_module_scores.tsv").rename(columns={"Unnamed: 0": "sample"})
    tls = read_tsv("project/results/d5p6_bcr_repertoire/tls_score_per_sample.tsv").rename(columns={"Unnamed: 0": "sample"})
    df = hla.merge(tls[["sample", "TLS_score"]], on="sample", how="inner")
    df["is_PTC_HT"] = df["group"].eq("PTC_HT")

    # A no-model composite. Z-scores are descriptive and explicitly not a fitted classifier.
    for score in ["HLA_I", "HLA_II", "TLS_score"]:
        df[f"{score}_z"] = (df[score] - df[score].mean()) / df[score].std(ddof=1)
    df["AP_TLS_composite"] = df[["HLA_I_z", "HLA_II_z", "TLS_score_z"]].mean(axis=1)
    df.to_csv(TABLES / "T03a_gse286332_expression_tls_composite_per_sample.tsv", sep="\t", index=False)

    rows = []
    for score in ["HLA_I", "HLA_II", "TLS_score", "AP_TLS_composite"]:
        vals = df[score].to_numpy(dtype=float)
        labels = df["is_PTC_HT"].to_numpy()
        obs_diff, p_exact, n_perm = exact_label_permutation(vals, labels)
        a = vals[labels]
        b = vals[~labels]
        d = cohen_d(pd.Series(a), pd.Series(b))
        lo, hi = bootstrap_d(a, b)
        auc = roc_auc_score(labels.astype(int), vals)
        rows.append(
            {
                "score": score,
                "n_PTC_HT": int(labels.sum()),
                "n_PTC": int((~labels).sum()),
                "mean_diff_PTC_HT_minus_PTC": obs_diff,
                "cohen_d": d,
                "bootstrap_d_ci95_lo": lo,
                "bootstrap_d_ci95_hi": hi,
                "exact_label_permutation_p": p_exact,
                "n_exact_label_assignments": n_perm,
                "rank_auc_descriptive": auc,
                "boundary": "HLA expression/TLS scores only; no cancer-cohort HLA allele association",
            }
        )
    out = pd.DataFrame(rows)
    out["BH_FDR_exact_p"] = bh_fdr(out["exact_label_permutation_p"].tolist())
    out.to_csv(TABLES / "T03_gse286332_exact_permutation_bootstrap.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(7.4, 4.5))
    xs = np.arange(len(out))
    ax.errorbar(
        xs,
        out["cohen_d"],
        yerr=[out["cohen_d"] - out["bootstrap_d_ci95_lo"], out["bootstrap_d_ci95_hi"] - out["cohen_d"]],
        fmt="o",
        color="#E8684A",
        ecolor="black",
        capsize=4,
        markersize=8,
    )
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels(out["score"], rotation=20, ha="right")
    ax.set_ylabel("Cohen's d with bootstrap 95% CI")
    ax.set_title("GSE286332 exact small-n robustness")
    for i, r in out.iterrows():
        ax.text(i, r["bootstrap_d_ci95_hi"] + 0.1, f"exact p={r['exact_label_permutation_p']:.4f}", ha="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "F03_gse286332_exact_permutation_bootstrap.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    plot_df = df.sort_values("AP_TLS_composite")
    ax.scatter(plot_df["AP_TLS_composite"], np.arange(len(plot_df)), c=plot_df["is_PTC_HT"].map({True: "#E8684A", False: "#5B8FF9"}), s=60, edgecolor="black")
    ax.axvline(plot_df.loc[plot_df["is_PTC_HT"], "AP_TLS_composite"].median(), color="#E8684A", linestyle="--", linewidth=1)
    ax.axvline(plot_df.loc[~plot_df["is_PTC_HT"], "AP_TLS_composite"].median(), color="#5B8FF9", linestyle="--", linewidth=1)
    ax.set_xlabel("Antigen-presentation/TLS composite")
    ax.set_ylabel("Samples sorted by composite")
    ax.set_title("PTC+HT separation by expression/TLS composite")
    fig.tight_layout()
    fig.savefig(FIGS / "F04_gse286332_ap_tls_composite_separation.png", dpi=220)
    plt.close(fig)
    return out, df


def tcga_hla2_dominance() -> pd.DataFrame:
    tcga = read_tsv(
        "project/results/hla_two_paper_synthesis_2026_05_09/extra_external_validation/tables/T02a_tcga_hla_ifng_tls_rai_per_sample.tsv"
    )
    cols = ["HLA2_score", "IFNG_hallmark", "TLS_score", "DM1_use", "RAI_module", "Immune_proxy", "Purity_proxy"]
    df = tcga[cols].apply(pd.to_numeric, errors="coerce").dropna().copy()
    y = df["HLA2_score"].to_numpy(dtype=float)
    predictors = ["IFNG_hallmark", "TLS_score", "DM1_use", "RAI_module", "Immune_proxy", "Purity_proxy"]
    X = df[predictors].copy()
    X = (X - X.mean()) / X.std(ddof=1)
    y_z = (y - y.mean()) / y.std(ddof=1)

    X_full = np.column_stack([np.ones(len(X)), X.to_numpy()])
    beta = np.linalg.lstsq(X_full, y_z, rcond=None)[0]
    pred = X_full @ beta
    ss_res = float(np.sum((y_z - pred) ** 2))
    ss_tot = float(np.sum((y_z - y_z.mean()) ** 2))
    full_r2 = 1.0 - ss_res / ss_tot

    rows = []
    for i, p in enumerate(predictors, start=1):
        X_drop = X.drop(columns=[p])
        Xd = np.column_stack([np.ones(len(X_drop)), X_drop.to_numpy()])
        bd = np.linalg.lstsq(Xd, y_z, rcond=None)[0]
        pred_d = Xd @ bd
        r2_drop = 1.0 - float(np.sum((y_z - pred_d) ** 2)) / ss_tot
        rho, pval = stats.spearmanr(df[p], df["HLA2_score"])
        rows.append(
            {
                "predictor": p,
                "n_samples": len(df),
                "univariate_spearman_rho_with_HLA2": float(rho),
                "univariate_spearman_p": float(pval),
                "multivariable_standardized_beta": float(beta[i]),
                "full_model_R2": full_r2,
                "drop_one_R2": r2_drop,
                "partial_R2_loss_if_removed": full_r2 - r2_drop,
                "boundary": "TCGA HLA-II expression-module dominance; not allele genotype",
            }
        )
    out = pd.DataFrame(rows).sort_values("partial_R2_loss_if_removed", ascending=False)
    out.to_csv(TABLES / "T04_tcga_hla2_driver_dominance.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    ax.bar(out["predictor"], out["partial_R2_loss_if_removed"], color="#61DDAA", edgecolor="black", linewidth=0.6)
    ax.set_ylabel("Drop-one R2 loss for HLA-II model")
    ax.set_title(f"TCGA HLA-II expression-module dominance (full R2={full_r2:.2f})")
    ax.set_xticks(range(len(out)))
    ax.set_xticklabels(out["predictor"], rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(FIGS / "F05_tcga_hla2_driver_dominance.png", dpi=220)
    plt.close(fig)
    return out


def build_wave2_ladder(tri: pd.DataFrame, gse: pd.DataFrame, dominance: pd.DataFrame) -> pd.DataFrame:
    c0102 = tri[tri["allele"] == "C*01:02"].iloc[0]
    dpb1 = tri[tri["allele"] == "DPB1*05:01"].iloc[0]
    composite = gse[gse["score"] == "AP_TLS_composite"].iloc[0]
    top_driver = dominance.iloc[0]
    rows = [
        {
            "paper": "Paper 4",
            "upgrade": "Allele triangulation figure",
            "why_it_raises_level": f"C*01:02 becomes a specificity-led candidate (score {c0102['triangulation_score_100']:.1f}); DPB1*05:01 remains class-II AITD anchor (score {dpb1['triangulation_score_100']:.1f}).",
            "reviewer_attack_closed": "This is just a list of old alleles.",
            "new_asset": "F01_paper4_allele_triangulation_map.png",
        },
        {
            "paper": "Paper 4",
            "upgrade": "Endpoint-sensitivity heatmap",
            "why_it_raises_level": "Strict GD, autoimmune hyperthyroid, broad hyperthyroid and HT/hypothyroid are separated instead of being mixed.",
            "reviewer_attack_closed": "External GWAS direction is inconsistent.",
            "new_asset": "F02_paper4_endpoint_robustness_heatmap.png",
        },
        {
            "paper": "Paper 2",
            "upgrade": "Exact small-n GSE286332 robustness",
            "why_it_raises_level": f"Composite AP/TLS score separates PTC+HT with d={composite['cohen_d']:.2f}, exact p={composite['exact_label_permutation_p']:.4f}, AUC={composite['rank_auc_descriptive']:.2f}.",
            "reviewer_attack_closed": "n=9 vs 9 is too small to trust asymptotic p-values.",
            "new_asset": "F03_gse286332_exact_permutation_bootstrap.png",
        },
        {
            "paper": "Paper 2",
            "upgrade": "TCGA HLA-II dominance model",
            "why_it_raises_level": f"{top_driver['predictor']} is the top drop-one contributor to HLA-II expression (R2 loss {top_driver['partial_R2_loss_if_removed']:.3f}).",
            "reviewer_attack_closed": "DM1-HLA-II is a loose correlation without mechanism prioritization.",
            "new_asset": "F05_tcga_hla2_driver_dominance.png",
        },
    ]
    ladder = pd.DataFrame(rows)
    ladder.to_csv(TABLES / "T05_wave2_level_up_ladder.tsv", sep="\t", index=False)
    return ladder


def write_report(tri: pd.DataFrame, gse: pd.DataFrame, dominance: pd.DataFrame, ladder: pd.DataFrame) -> None:
    c0102 = tri[tri["allele"] == "C*01:02"].iloc[0]
    dpb1 = tri[tri["allele"] == "DPB1*05:01"].iloc[0]
    hla2 = gse[gse["score"] == "HLA_II"].iloc[0]
    comp = gse[gse["score"] == "AP_TLS_composite"].iloc[0]
    top = dominance.iloc[0]
    lines = [
        "# HLA Two-Paper Level-Up Wave 2",
        "",
        "## 결론",
        "",
        "두 논문은 계속 2개로 유지하는 것이 맞습니다. 이번 wave는 새 논문을 쪼개기보다, 각 논문 안에 reviewer-proof 정량 방어층을 추가했습니다.",
        "",
        "## Paper 4 강화",
        "",
        f"- Allele triangulation: C*01:02는 GD-specificity lead로 triangulation score {c0102['triangulation_score_100']:.1f}; DPB1*05:01은 AITD-broad class-II anchor로 score {dpb1['triangulation_score_100']:.1f}.",
        "- Endpoint robustness: strict Graves, autoimmune hyperthyroid, broad hyperthyroid, HT/hypothyroid endpoint를 분리해 tag-SNP direction conflict를 phenotype sensitivity로 처리합니다.",
        "- 이 레이어가 Paper 4를 단순 meta-analysis에서 population-resolved HLA architecture/resource 논문으로 올립니다.",
        "",
        "## Paper 2 강화",
        "",
        f"- GSE286332 exact permutation: HLA-II d={hla2['cohen_d']:.2f}, exact p={hla2['exact_label_permutation_p']:.4f}; AP/TLS composite d={comp['cohen_d']:.2f}, exact p={comp['exact_label_permutation_p']:.4f}, AUC={comp['rank_auc_descriptive']:.2f}.",
        f"- TCGA HLA-II dominance: full model R2={top['full_model_R2']:.2f}; top drop-one contributor는 {top['predictor']} (R2 loss={top['partial_R2_loss_if_removed']:.3f}).",
        "- 이 레이어가 Paper 2를 underpowered allele story가 아니라 HT-overlap antigen-presentation/TLS immune-state 논문으로 고정합니다.",
        "",
        "## 새 산출물",
        "",
        f"- Tables: `{TABLES.relative_to(ROOT)}`",
        f"- Figures: `{FIGS.relative_to(ROOT)}`",
        "- `T01_paper4_allele_triangulation_score.tsv`",
        "- `T03_gse286332_exact_permutation_bootstrap.tsv`",
        "- `T04_tcga_hla2_driver_dominance.tsv`",
        "- `T05_wave2_level_up_ladder.tsv`",
        "",
        "## Boundary",
        "",
        "Paper 2의 HLA는 cancer cohort에서 expression module입니다. Allele association은 prospective validation roadmap에만 둡니다.",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    shutil.copy2(REPORT, OUT / "HLA_LEVEL_UP_WAVE2_KR.md")


def main() -> None:
    ensure_dirs()
    tri, branch = paper4_triangulation()
    gse, gse_per_sample = gse286332_exact_robustness()
    dominance = tcga_hla2_dominance()
    ladder = build_wave2_ladder(tri, gse, dominance)
    write_report(tri, gse, dominance, ladder)
    manifest = {
        "results_dir": str(OUT),
        "tables_dir": str(TABLES),
        "figures_dir": str(FIGS),
        "report": str(REPORT),
        "boundary": "Paper 2 HLA is expression-module only in cancer cohorts.",
        "n_tables": len(list(TABLES.glob("*.tsv"))),
        "n_figures": len(list(FIGS.glob("*.png"))),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
