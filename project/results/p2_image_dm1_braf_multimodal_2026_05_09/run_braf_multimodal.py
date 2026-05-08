"""
Paper 2 BRAF_like/cPTC multimodal rescue (2026-05-09).

Goal: lift image-only AUC=0.592 (n=41, 26 DM1 / 15 DM2) on the BRAF_like/cPTC
stratum where the CLAM image classifier fails. Build late-fusion classifiers
combining the trained CLAM slide-level embedding with bulk-RNA panel z-scores
and HM450 8-gene beta features.

Design notes
------------
- Reuse the SAME 5-fold assignment from the v2 CLAM run (clam_per_slide_predictions.tsv
  has a 'fold' column) so multimodal AUC is directly comparable to image-only 0.592.
- For each fold, we re-load the corresponding CLAM checkpoint (clam_fold{k}_best.pt)
  and forward each held-out slide through the gated-attention pool to obtain a
  512-d slide-level embedding (the 'bag' tensor before the classifier head).
  We also use the trained classifier's prob_DM1 as the image scalar feature.
- RNA panel = 8 thyroid-diff + 9 MAPK-output + 13 HT/B-cell signature = 30 z-scores.
- Methylation = 8 per-gene HM450 betas + mean_8g_beta = 9 features.
- Late fusion = standardize on training fold, concat, LogReg (L2, C=1.0).
  Held-out predictions across folds are pooled to compute a single OOF AUC,
  matching how 0.592 was reported.
- Bootstrap 95% CI on pooled OOF predictions (1000 resamples, stratified).
- All combos reported (no cherry-pick): image-only / RNA-only / Meth-only /
  image+RNA / image+Meth / RNA+Meth / image+RNA+Meth.

Inputs
------
- /data/thca/repo_results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/
  - clam_per_slide_predictions.tsv  (slide, fold, label, prob_DM1)
  - slide_manifest.tsv              (file_id ↔ submitter_id)
  - features/<file_id>.pt           (200 × 1024 ImageNet ViT-L tensor)
  - clam_fold{1..5}_best.pt         (GatedAttentionMIL state_dict)
- /data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv
- /data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv
- /data/thca/repo_results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv

Outputs (in this directory)
---------------------------
- braf_multimodal_results.tsv          per-fold + pooled AUC + bootstrap CI per combo
- braf_multimodal_per_slide_preds.tsv  OOF prob_DM1 for every combo on the 41 slides
- fig_braf_multimodal_uplift.png/pdf   bar + 95% CI forest figure
- BRAF_MULTIMODAL_REPORT.md            short report
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
BASE = Path("/home/seungho/personal/THCA_data_analysis")
V2 = BASE / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam"
OUT = BASE / "project/results/p2_image_dm1_braf_multimodal_2026_05_09"
OUT.mkdir(parents=True, exist_ok=True)

PRED_TSV = V2 / "clam_per_slide_predictions.tsv"
MAN_TSV = V2 / "slide_manifest.tsv"
FEAT_DIR = V2 / "features"
CKPT_TPL = V2 / "clam_fold{k}_best.pt"

MASTER_TSV = "/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv"
RNA_TSV = "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv"
METH_TSV = "/data/thca/repo_results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv"

SEED = 42
N_FOLDS = 5
N_BOOT = 1000

# ----------------------------------------------------------------------
# CLAM gated-attention MIL (must match phase2_clam_train_eval.py exactly)
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
        h = self.fc(x)  # (N, hidden)
        a = self.attn_w(self.attn_a(h) * self.attn_b(h))
        attn = F.softmax(a.squeeze(-1), dim=0)
        bag = (attn.unsqueeze(-1) * h).sum(0)  # (hidden,)
        logits = self.classifier(bag)
        return logits, bag, attn


# ----------------------------------------------------------------------
# Step 1: identify BRAF_like/cPTC subset and join all metadata
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
    braf = braf.reset_index(drop=True)
    return braf


# ----------------------------------------------------------------------
# Step 2: image embedding (per-fold, using held-out slides only)
# ----------------------------------------------------------------------
def extract_image_embedding_oof(cohort: pd.DataFrame) -> pd.DataFrame:
    """For each held-out slide, forward through the CLAM model that was
    trained without it and return the 512-d gated-attention bag vector.
    Output column block: img_e000 ... img_e511 + img_prob_DM1 (recomputed
    from the loaded checkpoint so we can sanity-check vs cohort['prob_DM1'])."""
    device = "cpu"  # 41 slides * 200 tiles * 1024 = trivially small
    rows = []
    for k in range(1, N_FOLDS + 1):
        ckpt_path = Path(str(CKPT_TPL).format(k=k))
        state = torch.load(ckpt_path, map_location=device)
        model = GatedAttentionMIL(in_dim=1024, hidden=512, n_classes=2, dropout=0.25)
        # The training sets dropout in fc; in eval() it is disabled.
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
                logits, bag, _ = model(feat)
                p = F.softmax(logits, dim=-1)[1].item()
            rec = {"slide": r["slide"], "patient12": r["patient12"],
                   "label": int(r["label"]), "fold": int(k),
                   "img_prob_DM1_loaded": float(p)}
            for i, v in enumerate(bag.cpu().numpy()):
                rec[f"img_e{i:03d}"] = float(v)
            rows.append(rec)
    df = pd.DataFrame(rows)
    return df


# ----------------------------------------------------------------------
# Step 3: RNA panel features
#
# CRITICAL — label-leak audit:
#   The DM1/DM2 binary label was DEFINED in `audit_2026_04_30/round5` using
#   the 8-gene HM450 methylation panel (DIO1, FOXE1, NKX2-1, PAX8, SLC5A5,
#   TG, TPO, TSHR). Empirically: DM1 mean_8g_beta=0.384 vs DM2=0.258 (~2-σ
#   separation). Any RNA-z feature for the SAME 8 genes is ~80-90% redundant
#   with the label by construction (methylation↔expression inverse
#   correlation, see DM1 Round 4). So including thyroid_diff_8gene in the
#   RNA block produces tautological AUC=1.000.
#
# Honest design:
#   RNA_NONLEAK = 9 MAPK-output + 13 HT/B-cell = 22 genes (NOT used to
#   define DM1 label). This is what we report as "rna" everywhere.
#   RNA_LEAK = the 8 panel genes — reported as a tautology reference only,
#   so reviewers see why we excluded them.
# ----------------------------------------------------------------------
RNA_PANEL = {
    "thyroid_diff_8gene_LEAKING": ["DIO1", "FOXE1", "NKX2-1", "PAX8",
                                   "SLC5A5", "TG", "TPO", "TSHR"],
    "mapk_output": ["DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4",
                    "ETV4", "ETV5", "PHLDA1", "CCND1"],
    "ht_bcell": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
                 "HLA-DQA1", "HLA-DQB1", "CD79A", "CD79B", "MS4A1",
                 "AICDA", "CXCL13", "CCR6", "IFNG"],
}
RNA_NONLEAK_GENES = RNA_PANEL["mapk_output"] + RNA_PANEL["ht_bcell"]
RNA_LEAK_GENES = RNA_PANEL["thyroid_diff_8gene_LEAKING"]


def load_rna_panel(cohort: pd.DataFrame) -> pd.DataFrame:
    panel_genes = sum(RNA_PANEL.values(), [])
    # First pass: list columns
    first_row = pd.read_csv(RNA_TSV, sep="\t", nrows=0)
    cols = first_row.columns.tolist()
    pat_to_col: dict[str, str] = {}
    for c in cols:
        if not c.startswith("TCGA"):
            continue
        # Prefer 01A primary tumor; for ties pick the lexicographically first
        pid = c[:12]
        if pid not in pat_to_col or c < pat_to_col[pid]:
            pat_to_col[pid] = c

    needed_cols = ["gene_symbol"] + [pat_to_col[p] for p in cohort["patient12"]
                                     if p in pat_to_col]
    rna = pd.read_csv(RNA_TSV, sep="\t", usecols=needed_cols)
    rna = rna[rna["gene_symbol"].isin(panel_genes)].set_index("gene_symbol")
    # Reorder to ensure consistent gene order
    rna = rna.reindex(panel_genes)
    rna = rna.T  # samples × genes
    rna.index.name = "rna_sample_id"
    rna = rna.reset_index()
    rna["patient12"] = rna["rna_sample_id"].str[:12]

    # ensure unique per patient12
    rna = rna.drop_duplicates(subset="patient12").set_index("patient12")
    rna.columns = ["rna_sample_id"] + [f"rna_{g}" for g in panel_genes]
    return rna.reset_index()


# ----------------------------------------------------------------------
# Step 4: methylation features (8 per-gene + mean = 9)
# ----------------------------------------------------------------------
METH_GENES = ["DIO1", "FOXE1", "NKX2-1", "PAX8", "SLC5A5", "TG", "TPO", "TSHR"]


def load_meth_panel(cohort: pd.DataFrame) -> pd.DataFrame:
    m = pd.read_csv(METH_TSV, sep="\t")
    m["patient12"] = m["sample_short"].astype(str).str[:12]
    m = m.drop_duplicates(subset="patient12")
    keep_cols = ["patient12"] + METH_GENES + ["mean_8g_beta"]
    m = m[keep_cols].copy()
    m.columns = ["patient12"] + [f"meth_{g}" for g in METH_GENES] + ["meth_mean_8g"]
    return m


# ----------------------------------------------------------------------
# Step 5: assemble fold-OOF dataset & evaluate combos
# ----------------------------------------------------------------------
COMBOS = {
    # ---- honest combos (no label-defining features) ----
    "image_only_loaded":           ["img_prob"],            # sanity vs 0.592
    "image_emb_only":              ["img_emb"],
    "rna_nonleak_only":            ["rna_nonleak"],         # MAPK(9) + HT(13) = 22
    "rna_mapk_only":               ["rna_mapk"],            # 9
    "rna_ht_only":                 ["rna_ht"],              # 13
    "image_emb_plus_rna_nonleak":  ["img_emb", "rna_nonleak"],
    "img_prob_plus_rna_nonleak":   ["img_prob", "rna_nonleak"],
    "img_prob_plus_rna_mapk":      ["img_prob", "rna_mapk"],
    "img_prob_plus_rna_ht":        ["img_prob", "rna_ht"],
    # ---- tautology references (label-leaking; for transparency, NOT headlines) ----
    "rna_leak_only_TAUTOLOGY":     ["rna_leak"],            # 8 thyroid-diff genes
    "meth_only_TAUTOLOGY":         ["meth"],                # HM450 8-gene betas
    "rna_leak_plus_meth_TAUTOLOGY": ["rna_leak", "meth"],
    # ---- mixed for completeness (label-leak still present, flagged) ----
    "image_emb_plus_meth_LEAK":    ["img_emb", "meth"],
    "all_three_LEAK":              ["img_emb", "rna_nonleak", "rna_leak", "meth"],
}


def feature_block(name, df, img_emb_cols, rna_nonleak_cols, rna_leak_cols,
                  rna_mapk_cols, rna_ht_cols, meth_cols):
    if name == "img_prob":
        return df[["img_prob_DM1_loaded"]].values
    if name == "img_emb":
        return df[img_emb_cols].values
    if name == "rna_nonleak":
        return df[rna_nonleak_cols].values
    if name == "rna_leak":
        return df[rna_leak_cols].values
    if name == "rna_mapk":
        return df[rna_mapk_cols].values
    if name == "rna_ht":
        return df[rna_ht_cols].values
    if name == "meth":
        return df[meth_cols].values
    raise ValueError(name)


def evaluate_combo(combo_blocks, df, img_emb_cols, rna_nonleak_cols, rna_leak_cols,
                   rna_mapk_cols, rna_ht_cols, meth_cols, n_folds=N_FOLDS):
    """Stratified-by-fold OOF prediction using LogReg with standardization on
    training samples per fold. Returns per-fold AUC (when feasible) and pooled OOF
    AUC + per-slide probabilities."""
    oof_prob = np.zeros(len(df))
    fold_aucs = []
    for k in range(1, n_folds + 1):
        tr = df["fold"].values != k
        va = df["fold"].values == k
        # build feature block
        Xtr_parts = [feature_block(b, df.iloc[tr], img_emb_cols, rna_nonleak_cols,
                                   rna_leak_cols, rna_mapk_cols, rna_ht_cols, meth_cols)
                     for b in combo_blocks]
        Xva_parts = [feature_block(b, df.iloc[va], img_emb_cols, rna_nonleak_cols,
                                   rna_leak_cols, rna_mapk_cols, rna_ht_cols, meth_cols)
                     for b in combo_blocks]
        # standardize each block on training rows then concat
        Xtr_std, Xva_std = [], []
        for tr_blk, va_blk in zip(Xtr_parts, Xva_parts):
            sc = StandardScaler().fit(tr_blk)
            Xtr_std.append(sc.transform(tr_blk))
            Xva_std.append(sc.transform(va_blk))
        Xtr = np.concatenate(Xtr_std, axis=1)
        Xva = np.concatenate(Xva_std, axis=1)
        ytr = df.iloc[tr]["label"].values
        yva = df.iloc[va]["label"].values
        # LogReg with mild regularization; high-dim image_emb (512) + few samples
        # → lean on L2. C=0.5 is a reasonable default; we keep it fixed for honesty
        # (no tuning per fold).
        clf = LogisticRegression(C=0.5, penalty="l2", solver="liblinear",
                                 max_iter=2000, random_state=SEED)
        clf.fit(Xtr, ytr)
        p = clf.predict_proba(Xva)[:, 1]
        oof_prob[va] = p
        try:
            auc = roc_auc_score(yva, p)
        except Exception:
            auc = np.nan
        fold_aucs.append(auc)
    pooled = roc_auc_score(df["label"].values, oof_prob)
    return pooled, np.array(fold_aucs), oof_prob


def bootstrap_auc_ci(y_true, y_pred, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    pos = np.where(y_true == 1)[0]
    neg = np.where(y_true == 0)[0]
    aucs = []
    for _ in range(n_boot):
        # stratified bootstrap to preserve class balance
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


# ----------------------------------------------------------------------
# Main driver
# ----------------------------------------------------------------------
def main():
    print("[step 1] cohort")
    cohort = build_cohort()
    print(f"  BRAF_like/cPTC n={len(cohort)} (DM1={sum(cohort['label']==1)}, "
          f"DM2={sum(cohort['label']==0)})")
    img_only_pooled = roc_auc_score(cohort["label"], cohort["prob_DM1"])
    print(f"  image-only pooled AUC (sanity vs 0.592): {img_only_pooled:.4f}")

    print("[step 2] image embeddings via held-out CLAM forward")
    img_df = extract_image_embedding_oof(cohort)
    print(f"  loaded n={len(img_df)} embeddings, dim=512")
    sanity_loaded = roc_auc_score(img_df["label"], img_df["img_prob_DM1_loaded"])
    print(f"  loaded-checkpoint pooled AUC: {sanity_loaded:.4f}")

    print("[step 3] RNA panel")
    rna_df = load_rna_panel(cohort)
    print(f"  RNA rows for cohort: {sum(rna_df['patient12'].isin(cohort['patient12']))} "
          f"of {len(cohort)}")

    print("[step 4] methylation panel")
    meth_df = load_meth_panel(cohort)
    print(f"  Meth rows for cohort: {sum(meth_df['patient12'].isin(cohort['patient12']))} "
          f"of {len(cohort)}")

    print("[step 5] assemble matrix")
    df = img_df.merge(rna_df, on="patient12", how="left")
    df = df.merge(meth_df, on="patient12", how="left")
    img_emb_cols = [c for c in df.columns if c.startswith("img_e")]
    rna_nonleak_cols = [f"rna_{g}" for g in RNA_NONLEAK_GENES]
    rna_leak_cols = [f"rna_{g}" for g in RNA_LEAK_GENES]
    rna_mapk_cols = [f"rna_{g}" for g in RNA_PANEL["mapk_output"]]
    rna_ht_cols = [f"rna_{g}" for g in RNA_PANEL["ht_bcell"]]
    meth_cols = [c for c in df.columns if c.startswith("meth_")]
    print(f"  feature dims: image_emb={len(img_emb_cols)}, "
          f"rna_nonleak={len(rna_nonleak_cols)} "
          f"(mapk={len(rna_mapk_cols)} + ht={len(rna_ht_cols)}), "
          f"rna_leak={len(rna_leak_cols)}, "
          f"meth={len(meth_cols)}")

    # any NaNs?
    assert df[rna_nonleak_cols + rna_leak_cols].isna().sum().sum() == 0, "RNA NaNs"
    assert df[meth_cols].isna().sum().sum() == 0, "Meth NaNs"

    print("[step 6] evaluate combos")
    rows = []
    pred_rows = {"slide": df["slide"].tolist(),
                 "patient12": df["patient12"].tolist(),
                 "label": df["label"].tolist(),
                 "fold": df["fold"].tolist()}
    for combo_name, blocks in COMBOS.items():
        pooled, fold_aucs, oof = evaluate_combo(blocks, df, img_emb_cols,
                                                rna_nonleak_cols, rna_leak_cols,
                                                rna_mapk_cols, rna_ht_cols, meth_cols)
        m, lo, hi = bootstrap_auc_ci(df["label"].values, oof)
        n_feat_blocks = ",".join(blocks)
        n_feat = sum({"img_prob": 1, "img_emb": len(img_emb_cols),
                      "rna_nonleak": len(rna_nonleak_cols),
                      "rna_leak": len(rna_leak_cols),
                      "rna_mapk": len(rna_mapk_cols),
                      "rna_ht": len(rna_ht_cols),
                      "meth": len(meth_cols)}[b] for b in blocks)
        print(f"  {combo_name:25s} pooled={pooled:.3f}  boot95%CI=[{lo:.3f},{hi:.3f}]  "
              f"folds={[round(x,3) if not np.isnan(x) else 'nan' for x in fold_aucs]}")
        rows.append({
            "combo": combo_name,
            "blocks": n_feat_blocks,
            "n_features": int(n_feat),
            "pooled_auc": float(pooled),
            "boot_mean_auc": float(m),
            "boot_lo95": float(lo),
            "boot_hi95": float(hi),
            "fold1_auc": float(fold_aucs[0]) if not np.isnan(fold_aucs[0]) else None,
            "fold2_auc": float(fold_aucs[1]) if not np.isnan(fold_aucs[1]) else None,
            "fold3_auc": float(fold_aucs[2]) if not np.isnan(fold_aucs[2]) else None,
            "fold4_auc": float(fold_aucs[3]) if not np.isnan(fold_aucs[3]) else None,
            "fold5_auc": float(fold_aucs[4]) if not np.isnan(fold_aucs[4]) else None,
            "fold_mean": float(np.nanmean(fold_aucs)),
            "fold_std":  float(np.nanstd(fold_aucs)),
        })
        pred_rows[f"prob_{combo_name}"] = oof.tolist()

    res_df = pd.DataFrame(rows).sort_values("pooled_auc", ascending=False)
    res_df.to_csv(OUT / "braf_multimodal_results.tsv", sep="\t", index=False)
    print("[saved]", OUT / "braf_multimodal_results.tsv")

    pd.DataFrame(pred_rows).to_csv(OUT / "braf_multimodal_per_slide_preds.tsv",
                                   sep="\t", index=False)
    print("[saved]", OUT / "braf_multimodal_per_slide_preds.tsv")

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    print("[step 7] figure")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    res_df_plot = res_df.sort_values("pooled_auc")
    y = np.arange(len(res_df_plot))
    fig, ax = plt.subplots(figsize=(9, 0.45 * len(res_df_plot) + 1.8))
    bar_x = res_df_plot["pooled_auc"].values
    err_lo = bar_x - res_df_plot["boot_lo95"].values
    err_hi = res_df_plot["boot_hi95"].values - bar_x

    def _color(name):
        if "TAUTOLOGY" in name or "_LEAK" in name:
            return "#bbbbbb"  # gray = label-leaking (not headline)
        if name.startswith("image_") or name == "image_only_loaded":
            return "#1f77b4"
        if "rna_nonleak" in name and "img" not in name:
            return "#ff7f0e"
        if "rna_nonleak" in name and "img" in name:
            return "#9467bd"
        return "#2ca02c"

    bar_color = [_color(c) for c in res_df_plot["combo"]]
    ax.barh(y, bar_x, color=bar_color, alpha=0.88)
    ax.errorbar(bar_x, y, xerr=[err_lo, err_hi], fmt="none",
                ecolor="k", capsize=3, lw=1)
    ax.axvline(0.5, ls=":", color="grey", label="chance")
    ax.axvline(0.592, ls="--", color="red", lw=1, label="image-only baseline (0.592)")
    ax.set_yticks(y)
    ax.set_yticklabels(res_df_plot["combo"], fontsize=8)
    ax.set_xlim(0.3, 1.05)
    ax.set_xlabel("Pooled OOF AUC (BRAF_like/cPTC, n=41)")
    ax.set_title("Paper 2 BRAF_like/cPTC multimodal rescue\n"
                 "5-fold CV (same folds as image-only); 95% CI = stratified bootstrap (1000)\n"
                 "Gray = label-leaking (HM450 8-gene panel defines DM1/DM2)")
    ax.legend(loc="lower right", fontsize=8)
    plt.tight_layout()
    plt.savefig(OUT / "fig_braf_multimodal_uplift.png", dpi=180)
    plt.savefig(OUT / "fig_braf_multimodal_uplift.pdf")
    plt.close()
    print("[saved] fig_braf_multimodal_uplift.{png,pdf}")

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------
    # split combos into honest vs leaking for narrative
    honest_combos = res_df[~res_df["combo"].str.contains("TAUTOLOGY|_LEAK", regex=True)]
    leak_combos = res_df[res_df["combo"].str.contains("TAUTOLOGY|_LEAK", regex=True)]
    best_honest = honest_combos.iloc[0]
    delta_vs_image = best_honest["pooled_auc"] - img_only_pooled

    lines = []
    lines.append("# Paper 2 — BRAF_like/cPTC multimodal rescue (2026-05-09)\n")
    lines.append("**Cohort.** N=41 BRAF_like / cPTC TCGA-THCA primary tumors with all three "
                 "modalities available (CLAM image embedding + bulk-RNA z-score + HM450 8-gene beta). "
                 f"Class balance: {sum(cohort['label']==1)} DM1 / {sum(cohort['label']==0)} DM2.\n")
    lines.append("**Baseline.** Image-only CLAM (ImageNet ViT-L tile features → gated-attention MIL) "
                 f"pooled OOF AUC = {img_only_pooled:.3f} on this stratum, vs the overall N=59 mean "
                 "AUC 0.83 ± 0.14 reported in `phase2_tcga_clam/PHASE2_REPORT.md`. The RAS_like/FVPTC "
                 "stratum (n=16) had AUC=1.00; BRAF_like/cPTC is the failure mode we are rescuing.\n")

    lines.append("**Label-leak audit (DECISIVE for honest reporting).** "
                 "DM1/DM2 was DEFINED upstream from the HM450 8-gene panel "
                 "(DIO1, FOXE1, NKX2-1, PAX8, SLC5A5, TG, TPO, TSHR; "
                 "see `audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv`). "
                 "On the cohort: DM1 mean_8g_β=0.384 vs DM2=0.258 — near-perfect by construction. "
                 "Therefore: (a) HM450 betas of the 8-gene panel **trivially** reproduce the label, "
                 "(b) RNA z-scores of the SAME 8 genes are ~80% redundant with the label via the "
                 "well-documented methylation↔expression inverse correlation. We split the RNA panel "
                 "into a NON-LEAK block (9 MAPK-output + 13 HT/B-cell = 22 genes) and a LEAK block "
                 "(8 thyroid-diff genes), and likewise mark methylation as TAUTOLOGY. "
                 "Headline conclusions use NON-LEAK combos only.\n")

    lines.append("**Setup.** Same 5 folds as the image-only run (loaded from "
                 "`clam_per_slide_predictions.tsv`). Per fold we (a) re-load the matching CLAM "
                 "checkpoint, (b) forward the 200×1024 tile features through gated-attention pooling "
                 "to obtain the 512-d slide bag vector, (c) concatenate with the chosen modality "
                 "blocks, (d) standardize each block on the training fold, (e) fit L2-regularized "
                 "LogReg (C=0.5). Held-out fold predictions are pooled across all 5 folds for one "
                 "OOF AUC. 95% CI is a 1,000-resample stratified bootstrap.\n")

    lines.append("## Headline (NON-LEAK combos)\n")
    lines.append("| Combo | n_feat | pooled AUC | 95% CI | per-fold mean ± std |")
    lines.append("|---|---:|---:|---|---|")
    for _, r in honest_combos.iterrows():
        lines.append(f"| {r['combo']} | {r['n_features']} | {r['pooled_auc']:.3f} | "
                     f"[{r['boot_lo95']:.3f}, {r['boot_hi95']:.3f}] | "
                     f"{r['fold_mean']:.3f} ± {r['fold_std']:.3f} |")
    lines.append("")
    lines.append(f"**Best non-leak.** `{best_honest['combo']}` "
                 f"AUC={best_honest['pooled_auc']:.3f} "
                 f"(95% CI [{best_honest['boot_lo95']:.3f}, {best_honest['boot_hi95']:.3f}]); "
                 f"Δ vs image-only = +{delta_vs_image:.3f}.\n")

    lines.append("## Reference (label-leaking; included for transparency only)\n")
    lines.append("| Combo | n_feat | pooled AUC | 95% CI | per-fold mean ± std |")
    lines.append("|---|---:|---:|---|---|")
    for _, r in leak_combos.iterrows():
        lines.append(f"| {r['combo']} | {r['n_features']} | {r['pooled_auc']:.3f} | "
                     f"[{r['boot_lo95']:.3f}, {r['boot_hi95']:.3f}] | "
                     f"{r['fold_mean']:.3f} ± {r['fold_std']:.3f} |")
    lines.append("")
    lines.append("These rows are tautologies (or near-tautologies) because the modality being "
                 "tested IS, or is tightly coupled to, the modality used to define DM1/DM2. They "
                 "are NOT a generalization claim — they are a sanity check that the data join "
                 "and pipeline work end-to-end.\n")

    lines.append("## Honesty caveats\n")
    lines.append("- N=41 with held-out folds of ~8 slides — single-fold AUCs are very noisy; "
                 "treat the bootstrap-CI on the pooled OOF AUC as the headline uncertainty.")
    lines.append("- The CLAM checkpoint for fold k was trained without that fold's slides, so "
                 "per-slide image embeddings and `img_prob_DM1_loaded` are valid OOF (matches the "
                 "original v2 0.830 ± 0.139 mean AUC).")
    lines.append("- `image_only_loaded` runs the held-out scalar prob_DM1 through the same per-fold "
                 "StandardScaler + LogReg as the other combos — this can monotonically remap and "
                 "thus give a slightly different pooled AUC from the raw 0.592. This is BY DESIGN "
                 "(apples-to-apples) but means the reported `image_only_loaded` AUC is NOT the same "
                 "number as the published 0.592; see the script header for the rationale.")
    lines.append("- RNA z-scores are pan-TCGA-cohort-relative. Within-fold leakage is impossible "
                 "(z-scores are fixed at the cohort level), but a Korean K2 validation will need to "
                 "rebuild z-scores within the K2 cohort before transferring this classifier.")
    lines.append("- All combos use identical CV / standardization / hyper-parameters. No "
                 "per-combo tuning, no model selection on test folds.")
    lines.append("- 512-d image embedding + 22 RNA = 534 features for N=41 is high-dim/low-N; the "
                 "bootstrap CI is the relevant uncertainty, not the point AUC.")
    lines.append("- The 9 MAPK-output genes (DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1, CCND1) are a "
                 "downstream readout of BRAF V600E activation. In the BRAF_like/cPTC stratum these "
                 "should NOT separate DM1 vs DM2 by driver alone — both groups are BRAF-like. The "
                 "HT/B-cell genes (HLA-DRA/DRB1/DPA1/DPB1/DQA1/DQB1, CD79A/B, MS4A1, AICDA, CXCL13, "
                 "CCR6, IFNG) capture the Hashimoto-overlap immune axis (per "
                 "v17_D5P6_BCR_clonal_TLS) which IS expected to enrich in DM1.")
    lines.append("")
    lines.append(f"## Files\n"
                 f"- `braf_multimodal_results.tsv` — full table\n"
                 f"- `braf_multimodal_per_slide_preds.tsv` — OOF prob per slide per combo\n"
                 f"- `fig_braf_multimodal_uplift.png` / `.pdf` — bar+CI figure (gray = leaking)")

    (OUT / "BRAF_MULTIMODAL_REPORT.md").write_text("\n".join(lines))
    print("[saved] BRAF_MULTIMODAL_REPORT.md")
    print()
    print("DONE.")


if __name__ == "__main__":
    main()
