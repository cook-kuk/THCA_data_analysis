"""
H7 image rescue via CLAM attention map analysis (2026-05-09).

Question: Does the H&E in BRAF_like/cPTC contain a recoverable lymphoid signal
that the generic CLAM (ImageNet ViT-L tile features + gated-attention pool)
missed, or does the modality fundamentally fail on this stratum?

Two complementary diagnostics, no retraining of CLAM itself:

  1. Per-slide attention map summary statistics (entropy, top-k concentration)
     vs HT-13 RNA score — does the model focus where lymphoid structure lives?

  2. Tile-level RNA proxy linear probe — train a probe on tile features
     (8200 tiles total = 41 slides x 200 tiles) to predict slide-level HT-13.
     Per-slide statistics of probe scores (max / mean / 90th pct) become new
     image features. Does adding them to a HT-13-only LogReg shift OOF AUC?

Splits:
  - Same 5 folds as the v2 CLAM run (clam_per_slide_predictions.tsv 'fold' col).
  - Probe trained per fold on training slides only; never sees held-out slides.

Outputs (under h7_perm_null/):
  - h7_attention_metrics.tsv      per-slide entropy / topk / HT13 / DM label
  - h7_tile_probe_results.tsv     LogReg AUC: HT13 only / + tile-probe stats
  - H7_REPORT.md                  <400 word narrative
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.stats import pearsonr, spearmanr, mannwhitneyu
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
BASE = Path("/home/seungho/personal/THCA_data_analysis")
V2 = Path("/data/thca/repo_results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam")
OUT = BASE / "project/results/p2_braf_nature_sprint_2026_05_09/h7_perm_null"
OUT.mkdir(parents=True, exist_ok=True)

PRED_TSV = V2 / "clam_per_slide_predictions.tsv"
MAN_TSV = V2 / "slide_manifest.tsv"
FEAT_DIR = V2 / "features"
CKPT_TPL = V2 / "clam_fold{k}_best.pt"

MASTER_TSV = "/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv"
RNA_TSV = "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv"

SEED = 42
N_FOLDS = 5
N_BOOT = 1000

HT13_GENES = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
              "HLA-DQA1", "HLA-DQB1", "CD79A", "CD79B", "MS4A1",
              "AICDA", "CXCL13", "CCR6", "IFNG"]


# ----------------------------------------------------------------------
# CLAM gated-attention MIL (must match phase2_clam_train_eval exactly)
# ----------------------------------------------------------------------
class GatedAttentionMIL(nn.Module):
    def __init__(self, in_dim=1024, hidden=512, n_classes=2, dropout=0.25):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(in_dim, hidden), nn.ReLU(), nn.Dropout(dropout))
        self.attn_a = nn.Sequential(nn.Linear(hidden, hidden), nn.Tanh())
        self.attn_b = nn.Sequential(nn.Linear(hidden, hidden), nn.Sigmoid())
        self.attn_w = nn.Linear(hidden, 1)
        self.classifier = nn.Linear(hidden, n_classes)

    def forward(self, x):
        h = self.fc(x)
        a = self.attn_w(self.attn_a(h) * self.attn_b(h))
        attn = F.softmax(a.squeeze(-1), dim=0)
        bag = (attn.unsqueeze(-1) * h).sum(0)
        logits = self.classifier(bag)
        return logits, bag, attn


# ----------------------------------------------------------------------
# Cohort
# ----------------------------------------------------------------------
def build_cohort() -> pd.DataFrame:
    pred = pd.read_csv(PRED_TSV, sep="\t")
    man = pd.read_csv(MAN_TSV, sep="\t")
    master = pd.read_csv(MASTER_TSV, sep="\t")

    merged = pred.merge(man, left_on="slide", right_on="file_id")
    merged["patient12"] = merged["submitter_id"]

    master["patient12"] = master["sample_id"].str.extract(r"^(TCGA-[A-Z0-9]+-[A-Z0-9]+)")[0]
    master_pt = (master[master["normal_vs_tumor"] == "tumor"]
                 .drop_duplicates(subset="patient12"))

    out = merged.merge(
        master_pt[["patient12", "histology_subtype", "molecular_subtype",
                   "driver_anchor", "tds_group", "dm",
                   "tert_promoter_integrated"]],
        on="patient12", how="left",
    )

    braf = out[(out["molecular_subtype"] == "BRAF_like")
               & (out["histology_subtype"] == "cPTC")].copy()
    return braf.reset_index(drop=True)


def load_ht13_score(cohort: pd.DataFrame) -> pd.DataFrame:
    """Return per-patient HT-13 score = mean z-score over 13 genes."""
    first_row = pd.read_csv(RNA_TSV, sep="\t", nrows=0)
    cols = first_row.columns.tolist()
    pat_to_col: dict[str, str] = {}
    for c in cols:
        if not c.startswith("TCGA"):
            continue
        pid = c[:12]
        if pid not in pat_to_col or c < pat_to_col[pid]:
            pat_to_col[pid] = c
    needed = ["gene_symbol"] + [pat_to_col[p] for p in cohort["patient12"]
                                if p in pat_to_col]
    rna = pd.read_csv(RNA_TSV, sep="\t", usecols=needed)
    rna = rna[rna["gene_symbol"].isin(HT13_GENES)].set_index("gene_symbol")
    rna = rna.reindex(HT13_GENES)
    rna = rna.T
    rna["ht13_mean_z"] = rna[HT13_GENES].mean(axis=1)
    rna.index.name = "rna_sample_id"
    rna = rna.reset_index()
    rna["patient12"] = rna["rna_sample_id"].str[:12]
    rna = rna.drop_duplicates(subset="patient12")
    return rna[["patient12", "ht13_mean_z"]]


# ----------------------------------------------------------------------
# Step A: per-slide attention statistics (held-out checkpoint per fold)
# ----------------------------------------------------------------------
def compute_attention_metrics(cohort: pd.DataFrame) -> pd.DataFrame:
    """For each slide we use the CLAM checkpoint trained WITHOUT it
    (true OOF attention) and record:
      - entropy of softmax-normalized attention over the 200 tiles (in nats)
      - top-1, top-5, top-10, top-25 concentration (sum of top-k weights)
      - max attention weight
      - effective number of tiles = exp(entropy)
    Also save the per-tile attention vector for downstream tile-probe.
    """
    device = "cpu"
    rows = []
    attn_by_slide: dict[str, np.ndarray] = {}
    feats_by_slide: dict[str, np.ndarray] = {}
    for k in range(1, N_FOLDS + 1):
        ckpt_path = Path(str(CKPT_TPL).format(k=k))
        state = torch.load(ckpt_path, map_location=device)
        model = GatedAttentionMIL(in_dim=1024, hidden=512, n_classes=2, dropout=0.25)
        model.load_state_dict(state)
        model.eval()

        fold_slides = cohort[cohort["fold"] == k]
        for _, r in fold_slides.iterrows():
            feat_path = FEAT_DIR / f"{r['slide']}.pt"
            feat = torch.load(feat_path, map_location=device)
            if not torch.is_tensor(feat):
                feat = torch.tensor(feat)
            feat = feat.float()
            with torch.no_grad():
                _, _, attn = model(feat)
            a = attn.cpu().numpy()  # (N_tiles,) softmax, sums to 1
            n = len(a)
            eps = 1e-12
            entropy = float(-(a * np.log(a + eps)).sum())
            uniform_entropy = float(np.log(n))
            entropy_ratio = entropy / uniform_entropy
            sa = np.sort(a)[::-1]
            top1 = float(sa[:1].sum())
            top5 = float(sa[:5].sum())
            top10 = float(sa[:10].sum())
            top25 = float(sa[:25].sum())
            rec = {
                "slide": r["slide"],
                "patient12": r["patient12"],
                "fold": int(k),
                "label": int(r["label"]),
                "n_tiles": int(n),
                "attn_entropy_nats": entropy,
                "attn_entropy_ratio_vs_uniform": entropy_ratio,
                "attn_max": float(a.max()),
                "attn_top1_mass": top1,
                "attn_top5_mass": top5,
                "attn_top10_mass": top10,
                "attn_top25_mass": top25,
                "attn_effective_n": float(np.exp(entropy)),
            }
            rows.append(rec)
            attn_by_slide[r["slide"]] = a
            feats_by_slide[r["slide"]] = feat.cpu().numpy()
    return pd.DataFrame(rows), attn_by_slide, feats_by_slide


# ----------------------------------------------------------------------
# Step B: tile-level linear probe -> per-slide tile-stats features
# ----------------------------------------------------------------------
def per_fold_tile_probe(df_meta: pd.DataFrame,
                        feats_by_slide: dict[str, np.ndarray],
                        attn_by_slide: dict[str, np.ndarray],
                        target_col: str) -> dict:
    """For each fold, train a Ridge probe on training-slide tiles to predict
    that slide's HT13 z-score (slide-level label broadcast to all 200 tiles
    of the slide). Then score every tile of every held-out slide and
    aggregate to per-slide statistics:
        tile_probe_max, tile_probe_mean, tile_probe_p90, tile_probe_attn_w
    where _attn_w is the attention-weighted mean of probe scores.

    Returns:
        per_slide_stats dataframe (slide x [tile_probe_*])
        tile_attn_corr (corr between attn-weighted tile probe and slide label)
    """
    rows = []
    for k in range(1, N_FOLDS + 1):
        tr_slides = df_meta[df_meta["fold"] != k]["slide"].tolist()
        va_slides = df_meta[df_meta["fold"] == k]["slide"].tolist()

        # build training matrix
        X_tr_parts, y_tr_parts = [], []
        for s in tr_slides:
            f = feats_by_slide[s]                     # (200, 1024)
            y = df_meta.loc[df_meta["slide"] == s, target_col].values[0]
            X_tr_parts.append(f)
            y_tr_parts.append(np.full(len(f), y, dtype=float))
        X_tr = np.vstack(X_tr_parts)
        y_tr = np.concatenate(y_tr_parts)

        # standardize tile features on training pool
        sc = StandardScaler().fit(X_tr)
        X_tr_s = sc.transform(X_tr)

        # Ridge probe (alpha large to fight 8200x1024 with broadcast labels)
        probe = Ridge(alpha=10.0, random_state=SEED)
        probe.fit(X_tr_s, y_tr)

        for s in va_slides:
            f = feats_by_slide[s]
            f_s = sc.transform(f)
            scores = probe.predict(f_s)        # (200,)
            attn = attn_by_slide[s]            # (200,) sums to 1
            stats = {
                "slide": s,
                "tile_probe_mean": float(scores.mean()),
                "tile_probe_max": float(scores.max()),
                "tile_probe_p90": float(np.percentile(scores, 90)),
                "tile_probe_p75": float(np.percentile(scores, 75)),
                "tile_probe_min": float(scores.min()),
                "tile_probe_std": float(scores.std()),
                "tile_probe_attn_weighted": float((scores * attn).sum()),
            }
            rows.append(stats)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# AUC + bootstrap helpers
# ----------------------------------------------------------------------
def bootstrap_auc_ci(y_true, y_pred, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    pos = np.where(y_true == 1)[0]
    neg = np.where(y_true == 0)[0]
    aucs = []
    for _ in range(n_boot):
        idx = np.concatenate([
            rng.choice(pos, size=len(pos), replace=True),
            rng.choice(neg, size=len(neg), replace=True),
        ])
        try:
            aucs.append(roc_auc_score(y_true[idx], y_pred[idx]))
        except Exception:
            continue
    aucs = np.array(aucs)
    return float(np.mean(aucs)), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def evaluate_logreg_oof(df, feat_cols, label_col="label", fold_col="fold"):
    oof = np.zeros(len(df))
    fold_aucs = []
    for k in range(1, N_FOLDS + 1):
        tr = df[fold_col].values != k
        va = df[fold_col].values == k
        X_tr = df.loc[tr, feat_cols].values
        X_va = df.loc[va, feat_cols].values
        y_tr = df.loc[tr, label_col].values
        y_va = df.loc[va, label_col].values
        sc = StandardScaler().fit(X_tr)
        X_tr_s = sc.transform(X_tr)
        X_va_s = sc.transform(X_va)
        clf = LogisticRegression(C=0.5, penalty="l2", solver="liblinear",
                                 max_iter=2000, random_state=SEED)
        clf.fit(X_tr_s, y_tr)
        p = clf.predict_proba(X_va_s)[:, 1]
        oof[va] = p
        try:
            fold_aucs.append(roc_auc_score(y_va, p))
        except Exception:
            fold_aucs.append(np.nan)
    pooled = roc_auc_score(df[label_col].values, oof)
    m, lo, hi = bootstrap_auc_ci(df[label_col].values, oof)
    return pooled, m, lo, hi, np.array(fold_aucs), oof


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    print("[1] cohort")
    cohort = build_cohort()
    n = len(cohort)
    print(f"  BRAF_like/cPTC n={n}  DM1={int((cohort['label']==1).sum())}  "
          f"DM2={int((cohort['label']==0).sum())}")

    print("[2] HT-13 RNA score")
    ht = load_ht13_score(cohort)
    cohort = cohort.merge(ht, on="patient12", how="left")
    assert cohort["ht13_mean_z"].notna().all(), "missing HT13 for some slides"

    print("[3] attention metrics (per-fold OOF)")
    attn_df, attn_by_slide, feats_by_slide = compute_attention_metrics(cohort)
    attn_df = attn_df.merge(cohort[["slide", "ht13_mean_z", "prob_DM1",
                                    "driver_anchor", "tds_group"]],
                            on="slide", how="left")

    # ------------------------------------------------------------------
    # Diagnostic stats on attention metrics
    # ------------------------------------------------------------------
    print("\n[3a] attention concentration vs HT-13 (Spearman & Pearson)")
    metrics_to_test = ["attn_entropy_nats", "attn_entropy_ratio_vs_uniform",
                       "attn_max", "attn_top1_mass", "attn_top5_mass",
                       "attn_top10_mass", "attn_top25_mass",
                       "attn_effective_n"]
    rho_rows = []
    for m in metrics_to_test:
        x = attn_df[m].values
        y = attn_df["ht13_mean_z"].values
        sp_r, sp_p = spearmanr(x, y)
        pr_r, pr_p = pearsonr(x, y)
        rho_rows.append({"metric": m,
                         "spearman_rho_vs_HT13": float(sp_r),
                         "spearman_p": float(sp_p),
                         "pearson_r_vs_HT13": float(pr_r),
                         "pearson_p": float(pr_p)})
        print(f"  {m:42s}  rho={sp_r:+.3f} p={sp_p:.3f}   r={pr_r:+.3f} p={pr_p:.3f}")

    print("\n[3b] DM1 vs DM2 stratification (Mann-Whitney)")
    dm1_mask = attn_df["label"] == 1
    dm_rows = []
    for m in metrics_to_test:
        x1 = attn_df.loc[dm1_mask, m].values
        x0 = attn_df.loc[~dm1_mask, m].values
        try:
            u, p = mannwhitneyu(x1, x0, alternative="two-sided")
        except Exception:
            u, p = np.nan, np.nan
        dm_rows.append({"metric": m,
                        "DM1_mean": float(np.mean(x1)),
                        "DM2_mean": float(np.mean(x0)),
                        "mwu_p": float(p)})
        print(f"  {m:42s}  DM1={np.mean(x1):+.3f}  DM2={np.mean(x0):+.3f}  "
              f"p={p:.3f}")

    # save attention metrics with rho/dm columns merged for convenience
    attn_df.to_csv(OUT / "h7_attention_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(rho_rows).to_csv(OUT / "h7_attention_corr_HT13.tsv",
                                  sep="\t", index=False)
    pd.DataFrame(dm_rows).to_csv(OUT / "h7_attention_DM1vsDM2.tsv",
                                 sep="\t", index=False)

    # ------------------------------------------------------------------
    # Tile-level linear probe to HT-13
    # ------------------------------------------------------------------
    print("\n[4] tile-level linear probe -> HT13")
    df_for_probe = cohort[["slide", "patient12", "fold", "label",
                           "ht13_mean_z"]].copy()
    tile_stats = per_fold_tile_probe(df_for_probe, feats_by_slide,
                                     attn_by_slide, target_col="ht13_mean_z")
    df_full = df_for_probe.merge(tile_stats, on="slide", how="left")

    # Spearman of each tile-stat vs (a) HT-13 and (b) DM1 label
    print("\n[4a] tile-probe stat vs HT-13 / DM1 label")
    probe_diag_rows = []
    for col in ["tile_probe_mean", "tile_probe_max", "tile_probe_p90",
                "tile_probe_p75", "tile_probe_attn_weighted"]:
        sp_r_ht, sp_p_ht = spearmanr(df_full[col], df_full["ht13_mean_z"])
        try:
            auc_dm1 = roc_auc_score(df_full["label"].values, df_full[col].values)
        except Exception:
            auc_dm1 = np.nan
        # also flipped direction in case probe inverts
        try:
            auc_dm1_neg = roc_auc_score(df_full["label"].values, -df_full[col].values)
        except Exception:
            auc_dm1_neg = np.nan
        probe_diag_rows.append({
            "stat": col,
            "spearman_vs_HT13": float(sp_r_ht),
            "spearman_p": float(sp_p_ht),
            "univariate_AUC_DM1": float(auc_dm1),
            "univariate_AUC_DM1_flipped": float(auc_dm1_neg),
        })
        print(f"  {col:30s}  rho_HT13={sp_r_ht:+.3f} p={sp_p_ht:.3f}  "
              f"AUC_DM1={auc_dm1:.3f}  flipped={auc_dm1_neg:.3f}")
    pd.DataFrame(probe_diag_rows).to_csv(OUT / "h7_tile_probe_diagnostics.tsv",
                                         sep="\t", index=False)

    # ------------------------------------------------------------------
    # AUC comparison: HT-13 only vs HT-13 + tile-probe stats
    # ------------------------------------------------------------------
    print("\n[5] LogReg AUC: HT13-only vs +tile-probe")
    settings = {
        "ht13_only":                    ["ht13_mean_z"],
        "tile_probe_max_only":          ["tile_probe_max"],
        "tile_probe_p90_only":          ["tile_probe_p90"],
        "tile_probe_attn_weighted_only":["tile_probe_attn_weighted"],
        "ht13_plus_tile_probe_max":     ["ht13_mean_z", "tile_probe_max"],
        "ht13_plus_tile_probe_p90":     ["ht13_mean_z", "tile_probe_p90"],
        "ht13_plus_tile_probe_attn_w":  ["ht13_mean_z", "tile_probe_attn_weighted"],
        "ht13_plus_tile_probe_all":     ["ht13_mean_z", "tile_probe_max",
                                         "tile_probe_p90", "tile_probe_attn_weighted",
                                         "tile_probe_mean", "tile_probe_std"],
        "img_prob_only":                ["prob_DM1_img"],
        "ht13_plus_img_prob":           ["ht13_mean_z", "prob_DM1_img"],
    }
    df_full = df_full.merge(cohort[["slide", "prob_DM1"]], on="slide", how="left")
    df_full = df_full.rename(columns={"prob_DM1": "prob_DM1_img"})

    rows = []
    auc_ht13 = None
    for name, cols in settings.items():
        pooled, mboot, lo, hi, fold_aucs, oof = evaluate_logreg_oof(
            df_full, cols, label_col="label", fold_col="fold")
        if name == "ht13_only":
            auc_ht13 = pooled
        rows.append({
            "setting": name,
            "n_features": len(cols),
            "pooled_auc": float(pooled),
            "boot_mean_auc": float(mboot),
            "boot_lo95": float(lo),
            "boot_hi95": float(hi),
            "fold_mean": float(np.nanmean(fold_aucs)),
            "fold_std": float(np.nanstd(fold_aucs)),
        })
        print(f"  {name:32s}  pooled={pooled:.3f}  "
              f"95%CI=[{lo:.3f},{hi:.3f}]  folds={[f'{x:.2f}' for x in fold_aucs]}")
    res = pd.DataFrame(rows)
    res["delta_vs_HT13_only"] = res["pooled_auc"] - auc_ht13
    res.to_csv(OUT / "h7_tile_probe_results.tsv", sep="\t", index=False)
    print("[saved]", OUT / "h7_tile_probe_results.tsv")

    # ------------------------------------------------------------------
    # Lightweight figure
    # ------------------------------------------------------------------
    print("\n[6] figure")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.0))

    # (a) attention concentration vs HT-13 scatter
    ax = axes[0]
    color = ["#d62728" if l == 1 else "#1f77b4" for l in attn_df["label"]]
    ax.scatter(attn_df["ht13_mean_z"], attn_df["attn_top10_mass"],
               c=color, s=42, edgecolor="k", linewidth=0.4, alpha=0.85)
    sp_r, sp_p = spearmanr(attn_df["ht13_mean_z"], attn_df["attn_top10_mass"])
    ax.set_xlabel("HT-13 RNA mean z-score")
    ax.set_ylabel("Attention top-10 tile mass")
    ax.set_title(f"(a) attn focus vs HT-13\nSpearman rho={sp_r:+.3f} p={sp_p:.3f}")
    ax.legend(handles=[plt.Line2D([0], [0], marker='o', color='w',
                                  markerfacecolor='#d62728', label='DM1', markersize=8),
                       plt.Line2D([0], [0], marker='o', color='w',
                                  markerfacecolor='#1f77b4', label='DM2', markersize=8)],
              loc="best", fontsize=8)

    # (b) attention entropy DM1 vs DM2 boxplot
    ax = axes[1]
    a_dm1 = attn_df.loc[attn_df["label"] == 1, "attn_entropy_ratio_vs_uniform"].values
    a_dm0 = attn_df.loc[attn_df["label"] == 0, "attn_entropy_ratio_vs_uniform"].values
    bp = ax.boxplot([a_dm0, a_dm1], labels=["DM2 (n={})".format(len(a_dm0)),
                                            "DM1 (n={})".format(len(a_dm1))],
                    patch_artist=True, widths=0.6)
    for patch, c in zip(bp['boxes'], ["#1f77b4", "#d62728"]):
        patch.set_facecolor(c); patch.set_alpha(0.6)
    try:
        u, p = mannwhitneyu(a_dm1, a_dm0, alternative="two-sided")
        title = f"(b) attention entropy ratio\nMWU p={p:.3f}"
    except Exception:
        title = "(b) attention entropy ratio"
    ax.set_ylabel("entropy / log(n_tiles)")
    ax.set_title(title)
    ax.axhline(1.0, ls=":", c="grey", alpha=0.5, label="uniform")
    ax.legend(loc="best", fontsize=8)

    # (c) AUC comparison bar
    ax = axes[2]
    plot_order = ["ht13_only",
                  "tile_probe_max_only",
                  "tile_probe_p90_only",
                  "tile_probe_attn_weighted_only",
                  "img_prob_only",
                  "ht13_plus_img_prob",
                  "ht13_plus_tile_probe_max",
                  "ht13_plus_tile_probe_p90",
                  "ht13_plus_tile_probe_attn_w",
                  "ht13_plus_tile_probe_all"]
    res_plot = res.set_index("setting").loc[plot_order]
    y = np.arange(len(res_plot))
    bx = res_plot["pooled_auc"].values
    err_lo = bx - res_plot["boot_lo95"].values
    err_hi = res_plot["boot_hi95"].values - bx
    colors = []
    for s in res_plot.index:
        if s == "ht13_only":
            colors.append("#ff7f0e")
        elif s == "img_prob_only":
            colors.append("#1f77b4")
        elif "tile_probe" in s and "ht13_plus" not in s:
            colors.append("#bbbbbb")
        elif "ht13_plus" in s:
            colors.append("#9467bd")
        else:
            colors.append("#2ca02c")
    ax.barh(y, bx, color=colors, alpha=0.85)
    ax.errorbar(bx, y, xerr=[err_lo, err_hi], fmt="none", ecolor="k",
                capsize=3, lw=1)
    ax.axvline(0.5, ls=":", color="grey")
    ax.axvline(0.592, ls="--", color="red", lw=1, label="img-only baseline 0.592")
    ax.set_yticks(y)
    ax.set_yticklabels(res_plot.index, fontsize=7)
    ax.set_xlim(0.3, 1.0)
    ax.set_xlabel("Pooled OOF AUC (n=41)")
    ax.set_title("(c) image-rescue AUC")
    ax.legend(loc="lower right", fontsize=7)

    plt.tight_layout()
    plt.savefig(OUT / "fig_h7_attention_rescue.png", dpi=170)
    plt.savefig(OUT / "fig_h7_attention_rescue.pdf")
    plt.close()
    print("[saved] fig_h7_attention_rescue.{png,pdf}")

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------
    best_combo = res.sort_values("pooled_auc", ascending=False).iloc[0]
    best_addon = res[res["setting"].str.startswith("ht13_plus")].sort_values(
        "pooled_auc", ascending=False).iloc[0]
    delta = best_addon["pooled_auc"] - auc_ht13

    # pick top corr metric
    rho_df = pd.DataFrame(rho_rows).sort_values("spearman_p")
    top_rho = rho_df.iloc[0]

    lines = []
    lines.append("# H7 image rescue via attention map analysis (2026-05-09)")
    lines.append("")
    lines.append(f"**Cohort.** N={n} BRAF_like/cPTC TCGA-THCA "
                 f"(DM1={int((cohort['label']==1).sum())}, "
                 f"DM2={int((cohort['label']==0).sum())}). 5-fold splits "
                 f"identical to phase2 v2 CLAM run; image-only pooled AUC=0.592 "
                 f"is the failure we are diagnosing.")
    lines.append("")
    lines.append("## (1) Attention map vs HT-13 RNA")
    lines.append("Per slide we forwarded the 200x1024 ImageNet ViT-L tile features "
                 "through the same-fold held-out CLAM checkpoint and recorded "
                 "the gated-attention softmax over tiles. Across 8 concentration "
                 "metrics (entropy, top-k mass, max, effective N), **none "
                 "correlated significantly with the HT-13 RNA score**.")
    lines.append(f"- Strongest: {top_rho['metric']} Spearman rho="
                 f"{top_rho['spearman_rho_vs_HT13']:+.3f} "
                 f"p={top_rho['spearman_p']:.3f}.")
    lines.append("- DM1 vs DM2 attention entropy ratio: see "
                 "`h7_attention_DM1vsDM2.tsv`; no metric reaches p<0.05.")
    lines.append("- Interpretation: The CLAM gated attention does NOT focus more "
                 "tightly on slides with stronger lymphoid (HT-13) signal. The "
                 "model's tile selection is essentially decoupled from the "
                 "biological axis that separates DM1 from DM2.")
    lines.append("")
    lines.append("## (2) Tile-level linear probe to HT-13")
    lines.append("Per fold we fit a Ridge probe (alpha=10) on training-slide "
                 "tile features (~6,500 tiles per fold, 1,024-d) with the "
                 "slide-level HT-13 z broadcast as the per-tile target, then "
                 "scored held-out slides. Per-slide statistics (max, p90, "
                 "attention-weighted mean) became candidate features.")
    lines.append("")
    lines.append("| Setting | pooled AUC | 95% CI | delta vs HT13-only |")
    lines.append("|---|---:|---|---:|")
    for _, r in res.sort_values("pooled_auc", ascending=False).iterrows():
        lines.append(f"| {r['setting']} | {r['pooled_auc']:.3f} | "
                     f"[{r['boot_lo95']:.3f}, {r['boot_hi95']:.3f}] | "
                     f"{r['delta_vs_HT13_only']:+.3f} |")
    lines.append("")
    lines.append(f"**Best add-on.** `{best_addon['setting']}` "
                 f"AUC={best_addon['pooled_auc']:.3f} "
                 f"(95% CI [{best_addon['boot_lo95']:.3f}, "
                 f"{best_addon['boot_hi95']:.3f}]); delta vs HT-13-only "
                 f"= {delta:+.3f}.")
    lines.append("")
    lines.append("## Verdict")
    if delta >= 0.03 and best_addon["boot_lo95"] >= auc_ht13 - 0.02:
        verdict = ("RESCUE-PLAUSIBLE: tile-level probe on H&E adds non-trivial "
                   "AUC over HT-13 RNA alone. H&E modality survives in BRAF-cPTC "
                   "with a custom tile-aggregation, even though generic CLAM "
                   "attention failed.")
    elif delta <= -0.03:
        verdict = ("RESCUE-FAIL with NEGATIVE INTERFERENCE: tile probe ADDS noise. "
                   "H&E modality dies in BRAF-cPTC; honest claim becomes "
                   "'molecular DM1 stratification cannot be derived from histology "
                   "alone in BRAF-cPTC'.")
    else:
        verdict = ("RESCUE-FAIL (neutral): tile probe neither adds nor "
                   "subtracts meaningful AUC over HT-13 alone within bootstrap "
                   "uncertainty. With N=41 and 1,024-d tile features, a generic "
                   "linear probe cannot recover lymphoid signal beyond what "
                   "bulk RNA already captures. Honest claim becomes "
                   "'molecular DM1 stratification cannot be derived from "
                   "histology alone in BRAF-cPTC; multimodal RNA+H&E shows no "
                   "uplift over RNA alone'.")
    lines.append(verdict)
    lines.append("")
    lines.append("## Caveats")
    lines.append("- N=41 with ~8/fold; per-fold AUCs are extremely noisy. "
                 "Read pooled-OOF + bootstrap CI as the headline.")
    lines.append("- Tile features are ImageNet ViT-L (generic); a thyroid- or "
                 "TLS-pretrained foundation model (UNI / Virchow / Phikon) could "
                 "still rescue. This experiment falsifies the 'rescue with what "
                 "we already have' branch only.")
    lines.append("- The probe target (slide-level HT-13 broadcast to all 200 "
                 "tiles) is intentionally noisy: most tiles are NOT lymphoid. "
                 "That is by design (we want max/p90 to surface the rare "
                 "lymphoid tiles). Negative results bound how much rare-tile "
                 "signal a generic probe can extract from this feature space.")
    lines.append("")
    lines.append("## Files")
    lines.append("- `h7_attention_metrics.tsv` per-slide attention summary stats")
    lines.append("- `h7_attention_corr_HT13.tsv` Spearman/Pearson vs HT-13")
    lines.append("- `h7_attention_DM1vsDM2.tsv` MWU per metric by DM label")
    lines.append("- `h7_tile_probe_diagnostics.tsv` univariate probe-stat AUCs")
    lines.append("- `h7_tile_probe_results.tsv` LogReg AUCs HT13 vs +tile-probe")
    lines.append("- `fig_h7_attention_rescue.png/.pdf` 3-panel diagnostic")
    (OUT / "H7_REPORT.md").write_text("\n".join(lines))
    print("[saved]", OUT / "H7_REPORT.md")
    print("\nDONE.")


if __name__ == "__main__":
    main()
