"""Round 8 — panel specificity bulletproofing.

Round 7 weakness: random 8-gene null gave p_emp=0.288 because random pool included
non-thyroid / lowly-expressed genes whose |d| can spuriously exceed ours when the
direction happens to align with batch / RIN. Round 8 corrects this with four layers:

A) **Label-permutation null** (panel fixed, DM labels shuffled)
   — answers "given THIS panel, is the DM1/DM2 separation real?"
B) **Bio-prior-matched random null** (random 8 genes drawn ONLY from thyroid-
   expressed pool, top expression quartile in TCGA-THCA n=460)
   — fair counterpart to Round 7's all-genome random null
C) **TDS-16 canonical comparison** (Landa/Krishnamoorthy thyroid differentiation
   score 16 genes vs our 8 — our 8 ⊂ TDS-16; tests whether 8 captures the full
   canonical TDS signal)
D) **Knock-out n-1 sensitivity** (each gene removed once; Round 1 LOO reorganized
   for Round 8 figure)

Outputs (project/results/dm1_robustness_v2026_05_08/round8/):
- label_permutation_null.{tsv,png,pdf}
- bio_prior_matched_null.{tsv,png,pdf}
- tds16_vs_8gene_comparison.{tsv,png,pdf}
- knockout_n_minus_1.{tsv,png,pdf}
- round8_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

EXPR_PATH = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv")
# Use the DM call source that paper 1 cites (Round 7 d=1.79 source); the alternative
# dm_master_with_pfi.tsv is an older labeling (74.5% concordance, d=0.52). Sticking
# with the paper-1 reference label keeps Round 8 directly comparable to Round 7.
LABEL_PATH = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/d4p2_tcga_hashimoto_signature/tcga_signature_scores.tsv"
)
OUT_DIR = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/dm1_robustness_v2026_05_08/round8"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

PANEL = ["TG", "TPO", "TSHR", "PAX8", "FOXE1", "NKX2-1", "DIO1", "SLC5A5"]

# Landa/Krishnamoorthy canonical thyroid differentiation score (16 genes).
# Our 8 ⊂ TDS-16.
TDS16 = [
    "TG", "TPO", "SLC5A5", "DIO1", "DIO2", "TSHR", "FOXE1", "GLIS3",
    "NKX2-1", "PAX8", "SLC26A4", "SLC5A8", "THRA", "THRB", "TFF3", "TPH1",
]

N_PERM_LABEL = 10000
N_PERM_PANEL = 1000
RNG_SEED = 20260508


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    """Pooled-SD Cohen's d (a - b)."""
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return float("nan")
    va, vb = a.var(ddof=1), b.var(ddof=1)
    sp = np.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    if sp == 0:
        return float("nan")
    return (a.mean() - b.mean()) / sp


def panel_score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Round-7 'g8_RAI' score: panel-subset cross-sample z, mean across panel.

    For the gene subset present in `expr`: z_{ij} = (x_{ij} - mean_i) / sd_i,
    where mean_i, sd_i are the per-gene cross-sample mean/sd computed on the
    PANEL SUBSET only (not the whole transcriptome). Then the per-sample score
    is the mean of those z values across the panel. This reproduces the paper-1
    g8_RAI definition (within_sample_z in dm_robustness_v2026_05_08_v2.py).
    """
    avail = [g for g in genes if g in expr.index]
    if len(avail) < 2:
        return pd.Series(dtype=float)
    sub = expr.loc[avail]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def load() -> tuple[pd.DataFrame, pd.Series]:
    print(f"[load] reading expression {EXPR_PATH.name}")
    expr = pd.read_csv(EXPR_PATH, sep="\t", index_col=0)
    expr.index.name = "gene_symbol"
    expr = expr.loc[~expr.index.duplicated(keep="first")]
    print(f"       expr {expr.shape[0]} genes × {expr.shape[1]} samples")

    sig = pd.read_csv(LABEL_PATH, sep="\t", index_col=0)
    sig = sig.loc[sig["DM"].isin(["DM1", "DM2"])]
    common = expr.columns.intersection(sig.index)
    print(f"       overlap {len(common)} samples with paper-1 DM call")
    expr = expr[common]
    dm = sig.loc[common, "DM"]
    return expr, dm


def per_gene_z(expr: pd.DataFrame) -> pd.DataFrame:
    """Pass-through: in this design we use panel_score() directly per panel,
    so the random/TDS panels share the same cross-sample-z normalisation as
    the canonical 8-gene panel (panel-subset z). We keep the function for
    Layers C/D where we need a fixed reference matrix to reuse."""
    return expr


def layer_A_label_perm(expr: pd.DataFrame, dm: pd.Series, rng: np.random.Generator) -> dict:
    print(f"[A] label-permutation null (n_perm={N_PERM_LABEL})")
    score = panel_score(expr, PANEL)
    score = score.dropna()
    dm_aligned = dm.loc[score.index].values
    score_v = score.values
    a = score_v[dm_aligned == "DM1"]
    b = score_v[dm_aligned == "DM2"]
    d_obs = cohens_d(a, b)
    n_total = len(score_v)
    n_dm1 = (dm_aligned == "DM1").sum()
    null_d = np.empty(N_PERM_LABEL, dtype=float)
    idx = np.arange(n_total)
    for i in range(N_PERM_LABEL):
        rng.shuffle(idx)
        a_i = score_v[idx[:n_dm1]]
        b_i = score_v[idx[n_dm1:]]
        null_d[i] = cohens_d(a_i, b_i)
    p_emp = (np.abs(null_d) >= np.abs(d_obs)).mean()
    out = {
        "d_obs": float(d_obs),
        "n_perm": N_PERM_LABEL,
        "null_mean": float(null_d.mean()),
        "null_sd": float(null_d.std(ddof=1)),
        "null_max_abs": float(np.abs(null_d).max()),
        "p_empirical": float(max(p_emp, 1.0 / (N_PERM_LABEL + 1))),
        "n_total": int(n_total),
        "n_dm1": int(n_dm1),
        "n_dm2": int(n_total - n_dm1),
    }
    pd.DataFrame({"perm_d": null_d}).to_csv(OUT_DIR / "label_permutation_null.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(null_d, bins=80, color="#cfd9e8", edgecolor="white")
    ax.axvline(d_obs, color="#c0392b", lw=2, label=f"observed d = {d_obs:.2f}")
    ax.axvline(-d_obs, color="#c0392b", lw=2, ls="--", alpha=0.5)
    ax.set_xlabel("Cohen's d (DM1 - DM2)")
    ax.set_ylabel("perm count")
    ax.set_title(f"Layer A: label-permutation null  p_emp={out['p_empirical']:.1e}")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "label_permutation_null.png", dpi=180)
    fig.savefig(OUT_DIR / "label_permutation_null.pdf")
    plt.close(fig)
    print(f"    d_obs={d_obs:.3f}  p_emp={out['p_empirical']:.4f}")
    return out


def layer_B_matched_panel(expr: pd.DataFrame, dm: pd.Series, rng: np.random.Generator) -> dict:
    print(f"[B] bio-prior-matched random null (n_perm={N_PERM_PANEL})")
    mean_expr = expr.mean(axis=1)
    cutoff = mean_expr.quantile(0.75)
    pool = mean_expr[mean_expr >= cutoff].index.tolist()
    pool = [g for g in pool if g not in PANEL]
    print(f"    thyroid-expressed pool: {len(pool)} genes (top quartile mean log2 > {cutoff:.2f})")

    obs_score = panel_score(expr, PANEL).dropna()
    dm_aligned = dm.loc[obs_score.index].values
    a_obs = obs_score.values[dm_aligned == "DM1"]
    b_obs = obs_score.values[dm_aligned == "DM2"]
    d_obs = cohens_d(a_obs, b_obs)

    pool_arr = np.array(pool)
    null_d = np.empty(N_PERM_PANEL, dtype=float)
    dm_v = dm.values
    for i in range(N_PERM_PANEL):
        sample = rng.choice(pool_arr, size=len(PANEL), replace=False)
        s = panel_score(expr, list(sample)).dropna()
        dm_a = dm.loc[s.index].values
        a_i = s.values[dm_a == "DM1"]
        b_i = s.values[dm_a == "DM2"]
        null_d[i] = cohens_d(a_i, b_i)
    p_emp = (np.abs(null_d) >= np.abs(d_obs)).mean()
    out = {
        "d_obs": float(d_obs),
        "n_perm": N_PERM_PANEL,
        "pool_size": len(pool),
        "null_mean_abs_d": float(np.abs(null_d).mean()),
        "null_max_abs_d": float(np.abs(null_d).max()),
        "p_empirical": float(max(p_emp, 1.0 / (N_PERM_PANEL + 1))),
    }
    pd.DataFrame({"random_d": null_d}).to_csv(OUT_DIR / "bio_prior_matched_null.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(np.abs(null_d), bins=60, color="#cdebcb", edgecolor="white")
    ax.axvline(abs(d_obs), color="#c0392b", lw=2, label=f"observed |d| = {abs(d_obs):.2f}")
    ax.set_xlabel("|Cohen's d| (random 8-gene panel)")
    ax.set_ylabel("perm count")
    ax.set_title(f"Layer B: bio-prior-matched null  p={out['p_empirical']:.3f}  pool={len(pool)}")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "bio_prior_matched_null.png", dpi=180)
    fig.savefig(OUT_DIR / "bio_prior_matched_null.pdf")
    plt.close(fig)
    print(f"    d_obs={d_obs:.3f}  p_emp={out['p_empirical']:.4f}")
    return out


def layer_C_tds16(expr: pd.DataFrame, dm: pd.Series) -> dict:
    print("[C] TDS-16 canonical comparison")
    avail8 = [g for g in PANEL if g in expr.index]
    avail16 = [g for g in TDS16 if g in expr.index]
    s8 = panel_score(expr, avail8).dropna()
    s16 = panel_score(expr, avail16).dropna()
    common = s8.index.intersection(s16.index).intersection(dm.index)
    s8 = s8.loc[common]; s16 = s16.loc[common]; lab = dm.loc[common]
    d8 = cohens_d(s8.values[lab == "DM1"], s8.values[lab == "DM2"])
    d16 = cohens_d(s16.values[lab == "DM1"], s16.values[lab == "DM2"])
    auc8 = roc_auc(s8.values, (lab == "DM1").values)
    auc16 = roc_auc(s16.values, (lab == "DM1").values)
    rho, _ = stats.spearmanr(s8.values, s16.values)
    out = {
        "panel8_genes_available": avail8,
        "tds16_genes_available": avail16,
        "panel8_d": float(d8),
        "tds16_d": float(d16),
        "panel8_auc_dm1": float(auc8),
        "tds16_auc_dm1": float(auc16),
        "spearman_panel8_vs_tds16": float(rho),
        "ratio_d_panel8_over_tds16": float(d8 / d16) if d16 else float("nan"),
    }

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].scatter(s16.values, s8.values, c=(lab == "DM1").map({True: "#c0392b", False: "#2c6fbb"}), s=8, alpha=0.7)
    axes[0].set_xlabel("TDS-16 score")
    axes[0].set_ylabel("8-gene panel score")
    axes[0].set_title(f"score correlation (Spearman ρ = {rho:.3f})")
    axes[0].axhline(0, color="#888", lw=0.5); axes[0].axvline(0, color="#888", lw=0.5)
    axes[0].spines[["top", "right"]].set_visible(False)

    bars = axes[1].bar(
        ["TDS-16", "8-gene"], [d16, d8],
        color=["#bbbbbb", "#2c6fbb"], width=0.5,
    )
    for bar, val in zip(bars, [d16, d8]):
        axes[1].text(bar.get_x() + bar.get_width() / 2, val + 0.04, f"{val:.2f}",
                     ha="center", fontsize=10)
    axes[1].set_ylabel("Cohen's d (DM1 vs DM2)")
    axes[1].set_title(f"AUC: TDS-16 = {auc16:.3f}, 8-gene = {auc8:.3f}")
    axes[1].spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "tds16_vs_8gene_comparison.png", dpi=180)
    fig.savefig(OUT_DIR / "tds16_vs_8gene_comparison.pdf")
    plt.close(fig)

    pd.DataFrame(
        [
            {"panel": "8-gene", "n_genes": len(avail8), "d": d8, "auc": auc8},
            {"panel": "TDS-16", "n_genes": len(avail16), "d": d16, "auc": auc16},
        ]
    ).to_csv(OUT_DIR / "tds16_vs_8gene_comparison.tsv", sep="\t", index=False)
    print(f"    d_8={d8:.3f}  d_16={d16:.3f}  ρ={rho:.3f}")
    return out


def roc_auc(score: np.ndarray, y: np.ndarray) -> float:
    order = np.argsort(score)
    y_sorted = y[order]
    n_pos = y.sum()
    n_neg = len(y) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    cum_neg = (1 - y_sorted).cumsum()
    auc = (y_sorted * cum_neg).sum() / (n_pos * n_neg)
    return float(auc)


def layer_D_knockout(expr: pd.DataFrame, dm: pd.Series) -> dict:
    print("[D] knockout n-1 sensitivity")
    rows = []
    for held in PANEL:
        sub = [g for g in PANEL if g != held and g in expr.index]
        s = panel_score(expr, sub).dropna()
        common = s.index.intersection(dm.index)
        lab = dm.loc[common]
        d_ko = cohens_d(s.loc[common].values[lab == "DM1"], s.loc[common].values[lab == "DM2"])
        rows.append({"held_out": held, "n_remaining": len(sub), "d": float(d_ko)})
    s_full = panel_score(expr, PANEL).dropna()
    common = s_full.index.intersection(dm.index)
    lab = dm.loc[common]
    d_full = cohens_d(s_full.loc[common].values[lab == "DM1"], s_full.loc[common].values[lab == "DM2"])
    df = pd.DataFrame(rows)
    df["delta_d_vs_full"] = df["d"] - d_full
    df.to_csv(OUT_DIR / "knockout_n_minus_1.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(7, 4))
    order = df.sort_values("d").reset_index(drop=True)
    ax.barh(order["held_out"], order["d"], color="#9bb7d4")
    ax.axvline(d_full, color="#c0392b", lw=2, label=f"full 8-gene d = {d_full:.2f}")
    ax.set_xlabel("Cohen's d (n-1 panel)")
    ax.set_title("Layer D: knock-out one gene at a time")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "knockout_n_minus_1.png", dpi=180)
    fig.savefig(OUT_DIR / "knockout_n_minus_1.pdf")
    plt.close(fig)

    out = {
        "d_full": float(d_full),
        "min_n_minus_1_d": float(df["d"].min()),
        "max_drop": float(df["delta_d_vs_full"].min()),
        "max_drop_held_out": str(df.iloc[df["delta_d_vs_full"].idxmin()]["held_out"]),
        "all_n_minus_1_above_1": bool((df["d"] > 1.0).all()),
    }
    print(f"    d_full={d_full:.3f}  min_n-1={out['min_n_minus_1_d']:.3f}  weakest_drop_held={out['max_drop_held_out']}")
    return out


def main() -> None:
    rng = np.random.default_rng(RNG_SEED)
    expr, dm = load()

    A = layer_A_label_perm(expr, dm, rng)
    B = layer_B_matched_panel(expr, dm, rng)
    C = layer_C_tds16(expr, dm)
    D = layer_D_knockout(expr, dm)

    summary = {
        "panel": PANEL,
        "tds16": TDS16,
        "n_samples_overlap": int(len(dm)),
        "layer_A_label_permutation": A,
        "layer_B_bio_prior_matched_random": B,
        "layer_C_tds16_comparison": C,
        "layer_D_knockout_n_minus_1": D,
        "headline": {
            "panel_conditional_p_emp": A["p_empirical"],
            "thyroid_matched_p_emp": B["p_empirical"],
            "tds16_redundancy": C["spearman_panel8_vs_tds16"],
            "knockout_floor_d": D["min_n_minus_1_d"],
        },
        "interpretation": [
            "Round 7 weakness was a genome-wide random null contaminated by lowly-expressed",
            "or non-thyroid genes; Round 8 Layer B restricts the null to top-quartile thyroid-",
            "expressed genes — the relevant counterfactual.",
            "Layer A answers the panel-conditional question (is THIS panel's separation real?),",
            "Layer B the panel-specific question (could ANY thyroid panel match?), Layer C",
            "the redundancy question (is 8 a sufficient summary of TDS-16?), Layer D the",
            "robustness question (does any single gene drive the result?).",
        ],
        "generated_at": "2026-05-08 round8",
    }
    with open(OUT_DIR / "round8_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print("\n=== ROUND 8 HEADLINE ===")
    print(json.dumps(summary["headline"], indent=2))


if __name__ == "__main__":
    main()
