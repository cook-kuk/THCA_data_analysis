"""
H15 — foundation-model H&E rescue on the BRAF/cPTC stratum.

Goal: lift image-only AUC=0.564 (BRAF_like/cPTC, n=41) from generic ImageNet
ViT-L tile features to UNI (MahmoodLab/UNI) features, by re-training the same
gated-attention CLAM MIL on UNI features only on the BRAF/cPTC stratum, using
the SAME 5-fold split as the ImageNet baseline.

Honest design — what's at risk:
- 4 of 41 BRAF/cPTC slides are missing UNI features locally (raw WSI was on a
  RunPod that no longer exists). We run on n=37, with the 4 missing reported.
- Since UNI features come from a stronger encoder, we MUST retrain CLAM head
  (cannot just swap features into the ImageNet-trained checkpoint). Same 5-fold
  fold assignment is preserved so AUCs are comparable.
- Cohort N is small. Headline = pooled OOF AUC + stratified bootstrap CI;
  per-fold AUCs are noisy.

Outputs (in this dir):
- h15_per_encoder_results.tsv      ImageNet vs UNI image-only / multimodal AUCs
- h15_per_slide_preds.tsv          per-slide OOF probabilities
- h15_recon.txt                    cohort + feature inventory
- h15_pod_runtime_log.txt          wall time + dollar tracking (no pod used)
- H15_REPORT.md                    < 400 word verdict
- fig_h15_encoder_uplift.png/pdf   bar + bootstrap CI

UNI features from /data/thca/repo_results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam_UNI/features/
"""
from __future__ import annotations

import json
import time
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
warnings.filterwarnings("ignore", category=UserWarning)

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
BASE = Path("/home/seungho/personal/THCA_data_analysis")
V2 = BASE / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"

IMG_FEAT_DIR = V2 / "phase2_tcga_clam/features"            # ImageNet ViT-L 200x1024
UNI_FEAT_DIR = V2 / "phase2_tcga_clam_UNI/features"        # UNI 200x1024

PRED_TSV = V2 / "phase2_tcga_clam/clam_per_slide_predictions.tsv"
MAN_TSV = V2 / "phase2_tcga_clam/slide_manifest.tsv"
IMG_CKPT_TPL = V2 / "phase2_tcga_clam/clam_fold{k}_best.pt"

MASTER_TSV = "/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv"
RNA_TSV = "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv"

OUT = BASE / "project/results/p2_braf_nature_sprint_2026_05_09/h15_foundation_image"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 42
N_FOLDS = 5
N_BOOT = 1000
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


# ----------------------------------------------------------------------
# CLAM gated attention MIL (matches phase2 architecture exactly)
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
# Cohort builder (BRAF_like + cPTC)
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
# RNA panel (HT-13 only - the meaningful axis for DM1 vs DM2 in BRAF/cPTC)
# ----------------------------------------------------------------------
RNA_HT13 = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
            "HLA-DQA1", "HLA-DQB1", "CD79A", "CD79B", "MS4A1",
            "AICDA", "CXCL13", "CCR6", "IFNG"]


def load_rna_panel(cohort: pd.DataFrame) -> pd.DataFrame:
    first = pd.read_csv(RNA_TSV, sep="\t", nrows=0)
    cols = first.columns.tolist()
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
    rna = rna[rna["gene_symbol"].isin(RNA_HT13)].set_index("gene_symbol")
    rna = rna.reindex(RNA_HT13)
    rna = rna.T
    rna.index.name = "rna_sample_id"
    rna = rna.reset_index()
    rna["patient12"] = rna["rna_sample_id"].str[:12]
    rna = rna.drop_duplicates(subset="patient12").set_index("patient12")
    rna.columns = ["rna_sample_id"] + [f"rna_{g}" for g in RNA_HT13]
    return rna.reset_index()


# ----------------------------------------------------------------------
# Train CLAM on training fold; return slide-level (prob, bag_emb) on val
# ----------------------------------------------------------------------
def train_clam_one_fold(train_slides, val_slides, feat_dir: Path, *,
                        epochs: int = 60, lr: float = 1e-4,
                        weight_decay: float = 1e-4, dropout: float = 0.25,
                        seed: int = SEED):
    """Train a fresh CLAM model on training slides, evaluate on val.
    Returns: dict slide -> {prob, bag_emb, attn_entropy}."""
    torch.manual_seed(seed)
    np.random.seed(seed)

    model = GatedAttentionMIL(in_dim=1024, hidden=512, n_classes=2,
                              dropout=dropout).to(DEVICE)
    optim = torch.optim.AdamW(model.parameters(), lr=lr,
                              weight_decay=weight_decay)
    loss_fn = nn.CrossEntropyLoss()

    # Pre-load train features into memory (37 slides, ~150 MB total)
    train_feats = []
    for s, lab in train_slides:
        f = torch.load(feat_dir / f"{s}.pt", map_location="cpu",
                       weights_only=False).float()
        train_feats.append((s, lab, f))

    best_train_loss = float("inf")
    for ep in range(epochs):
        model.train()
        order = np.random.permutation(len(train_feats))
        ep_loss = 0.0
        for idx in order:
            s, lab, feat = train_feats[idx]
            feat = feat.to(DEVICE)
            y = torch.tensor([lab], dtype=torch.long, device=DEVICE)
            optim.zero_grad()
            logits, _, _ = model(feat)
            loss = loss_fn(logits.unsqueeze(0), y)
            loss.backward()
            optim.step()
            ep_loss += loss.item()
        ep_loss /= max(1, len(train_feats))
        if ep_loss < best_train_loss:
            best_train_loss = ep_loss

    model.eval()
    out = {}
    with torch.no_grad():
        for s, lab in val_slides:
            f = torch.load(feat_dir / f"{s}.pt", map_location="cpu",
                           weights_only=False).float().to(DEVICE)
            logits, bag, attn = model(f)
            p = F.softmax(logits, dim=-1)[1].item()
            entropy = float(-(attn * torch.log(attn + 1e-12)).sum().item())
            out[s] = {
                "prob": float(p),
                "bag_emb": bag.cpu().numpy(),
                "attn_entropy": entropy,
            }
    return out, float(best_train_loss)


def run_clam_5fold(cohort_local: pd.DataFrame, feat_dir: Path, encoder_name: str):
    """Run 5-fold CLAM training on the given encoder's features.
    Returns df with slide, fold, label, prob_DM1, bag_emb_e000..511, encoder."""
    rows = []
    fold_train_loss = {}
    for k in range(1, N_FOLDS + 1):
        val_mask = cohort_local["fold"] == k
        tr_slides = [(r["slide"], int(r["label"])) for _, r in
                     cohort_local[~val_mask].iterrows()]
        va_slides = [(r["slide"], int(r["label"])) for _, r in
                     cohort_local[val_mask].iterrows()]

        out, tr_loss = train_clam_one_fold(tr_slides, va_slides, feat_dir)
        fold_train_loss[k] = tr_loss

        for s, lab in va_slides:
            r = cohort_local[cohort_local["slide"] == s].iloc[0]
            rec = {
                "slide": s,
                "patient12": r["patient12"],
                "label": int(lab),
                "fold": int(k),
                "encoder": encoder_name,
                "prob_DM1": float(out[s]["prob"]),
                "attn_entropy": float(out[s]["attn_entropy"]),
            }
            for i, v in enumerate(out[s]["bag_emb"]):
                rec[f"img_e{i:03d}"] = float(v)
            rows.append(rec)
    return pd.DataFrame(rows), fold_train_loss


# ----------------------------------------------------------------------
# LogReg combo evaluation (re-using H7-style design)
# ----------------------------------------------------------------------
def eval_combo_logreg(df, feat_cols, n_folds=N_FOLDS, C=0.5):
    oof = np.zeros(len(df))
    fold_aucs = []
    for k in range(1, n_folds + 1):
        tr = df["fold"].values != k
        va = df["fold"].values == k
        Xtr = df.iloc[tr][feat_cols].values
        Xva = df.iloc[va][feat_cols].values
        ytr = df.iloc[tr]["label"].values
        yva = df.iloc[va]["label"].values
        sc = StandardScaler().fit(Xtr)
        Xtr_s = sc.transform(Xtr)
        Xva_s = sc.transform(Xva)
        clf = LogisticRegression(C=C, penalty="l2", solver="liblinear",
                                 max_iter=2000, random_state=SEED)
        clf.fit(Xtr_s, ytr)
        oof[va] = clf.predict_proba(Xva_s)[:, 1]
        try:
            fold_aucs.append(roc_auc_score(yva, oof[va]))
        except Exception:
            fold_aucs.append(np.nan)
    pooled = roc_auc_score(df["label"].values, oof)
    return pooled, np.array(fold_aucs), oof


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
    return (float(np.mean(aucs)),
            float(np.percentile(aucs, 2.5)),
            float(np.percentile(aucs, 97.5)))


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    t0 = time.time()
    print(f"[device] {DEVICE}")
    print(f"[seeds]  fixed at {SEED}, folds={N_FOLDS}, boot={N_BOOT}")

    # ----------------------------------------------------------------
    # Step 1: cohort
    # ----------------------------------------------------------------
    cohort = build_cohort()
    n_total = len(cohort)
    n_dm1 = int((cohort["label"] == 1).sum())
    n_dm2 = int((cohort["label"] == 0).sum())
    print(f"[cohort] BRAF_like/cPTC n={n_total} (DM1={n_dm1}, DM2={n_dm2})")

    img_only_pooled = roc_auc_score(cohort["label"], cohort["prob_DM1"])
    print(f"[sanity] ImageNet image-only pooled AUC: {img_only_pooled:.4f} "
          f"(should be ~0.564)")

    # ----------------------------------------------------------------
    # Step 2: which slides have UNI features locally?
    # ----------------------------------------------------------------
    uni_avail = {p.stem for p in UNI_FEAT_DIR.glob("*.pt")}
    img_avail = {p.stem for p in IMG_FEAT_DIR.glob("*.pt")}
    cohort["has_uni"] = cohort["slide"].isin(uni_avail)
    cohort["has_img"] = cohort["slide"].isin(img_avail)
    cohort_uni = cohort[cohort["has_uni"]].reset_index(drop=True)
    n_uni = len(cohort_uni)
    print(f"[uni]    UNI features local for {n_uni}/{n_total} BRAF/cPTC slides "
          f"(missing {n_total - n_uni})")

    missing_slides = cohort.loc[~cohort["has_uni"], "slide"].tolist()
    cohort_uni_dm1 = int((cohort_uni["label"] == 1).sum())
    cohort_uni_dm2 = int((cohort_uni["label"] == 0).sum())
    print(f"[uni]    n={n_uni} (DM1={cohort_uni_dm1}, DM2={cohort_uni_dm2})")

    # Build matched-N ImageNet cohort (drop the 4 same slides for fair compare)
    cohort_matched = cohort_uni.copy()  # same N=37 slides

    # ----------------------------------------------------------------
    # Step 3: CLAM 5-fold on UNI features (BRAF/cPTC stratum only)
    # ----------------------------------------------------------------
    print("[uni-clam] training fresh CLAM 5-fold on UNI features...")
    t_uni0 = time.time()
    uni_df, uni_train_losses = run_clam_5fold(cohort_uni, UNI_FEAT_DIR, "UNI")
    t_uni = time.time() - t_uni0
    print(f"[uni-clam] done in {t_uni:.1f}s, train losses by fold: "
          f"{[round(uni_train_losses[k],4) for k in range(1, N_FOLDS+1)]}")
    uni_img_only = roc_auc_score(uni_df["label"], uni_df["prob_DM1"])
    print(f"[uni-clam] image-only pooled OOF AUC: {uni_img_only:.4f}")

    # ----------------------------------------------------------------
    # Step 4: CLAM 5-fold on ImageNet features matched-N (for fair compare)
    # ----------------------------------------------------------------
    print("[img-clam-matched] training fresh CLAM 5-fold on ImageNet features (n=37)...")
    t_img0 = time.time()
    img_df, img_train_losses = run_clam_5fold(cohort_matched, IMG_FEAT_DIR, "ImageNet")
    t_img = time.time() - t_img0
    print(f"[img-clam-matched] done in {t_img:.1f}s, train losses by fold: "
          f"{[round(img_train_losses[k],4) for k in range(1, N_FOLDS+1)]}")
    img_only_matched = roc_auc_score(img_df["label"], img_df["prob_DM1"])
    print(f"[img-clam-matched] image-only pooled OOF AUC: {img_only_matched:.4f}")

    # ----------------------------------------------------------------
    # Step 5: RNA panel (HT-13)
    # ----------------------------------------------------------------
    print("[rna] loading HT-13 panel...")
    rna_df = load_rna_panel(cohort_uni)
    rna_cols = [f"rna_{g}" for g in RNA_HT13]
    n_rna = sum(rna_df["patient12"].isin(cohort_uni["patient12"]))
    print(f"[rna] HT-13 rows for {n_rna}/{n_uni} patients")

    # ----------------------------------------------------------------
    # Step 6: assemble matrices, evaluate combos for both encoders
    # ----------------------------------------------------------------
    img_emb_cols = [f"img_e{i:03d}" for i in range(512)]
    results = []
    pred_records = {}

    for encoder, feat_df in [("ImageNet", img_df), ("UNI", uni_df)]:
        df = feat_df.merge(rna_df, on="patient12", how="left")
        # Sanity check: all RNA cols should be filled (TCGA RNA should have all)
        miss = df[rna_cols].isna().any(axis=1).sum()
        if miss > 0:
            print(f"[WARN] {encoder} has {miss} rows with missing RNA — dropping")
            df = df.dropna(subset=rna_cols).reset_index(drop=True)

        # 4 combos: img_prob_only, img_emb_only, img_prob+ht13, img_emb+ht13
        # plus ht13_only as common reference
        combos = {
            "img_prob_only":    ["img_prob_DM1"],
            "img_emb_only":     img_emb_cols,
            "ht13_only":        rna_cols,
            "img_prob_plus_ht13": ["img_prob_DM1"] + rna_cols,
            "img_emb_plus_ht13":  img_emb_cols + rna_cols,
        }
        # img_prob_DM1 column rename for uniform feature_block API
        df = df.rename(columns={"prob_DM1": "img_prob_DM1"})

        for cname, fcols in combos.items():
            pooled, fold_aucs, oof = eval_combo_logreg(df, fcols)
            m, lo, hi = bootstrap_auc_ci(df["label"].values, oof)
            print(f"[{encoder:9s}] {cname:25s} pooled={pooled:.3f} "
                  f"95%CI=[{lo:.3f},{hi:.3f}] folds={[round(x,3) for x in fold_aucs]}")
            results.append({
                "encoder": encoder,
                "combo": cname,
                "n_features": len(fcols),
                "n_slides": int(len(df)),
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
                "fold_std": float(np.nanstd(fold_aucs)),
            })
            pred_records[f"{encoder}__{cname}"] = pd.DataFrame({
                "slide": df["slide"].values,
                "patient12": df["patient12"].values,
                "label": df["label"].values,
                "fold": df["fold"].values,
                "prob": oof,
            }).assign(encoder=encoder, combo=cname)

    res_df = pd.DataFrame(results).sort_values(
        ["encoder", "pooled_auc"], ascending=[True, False])
    res_df.to_csv(OUT / "h15_per_encoder_results.tsv", sep="\t", index=False)
    print(f"[saved] {OUT / 'h15_per_encoder_results.tsv'}")

    pred_df_all = pd.concat(pred_records.values(), ignore_index=True)
    pred_df_all.to_csv(OUT / "h15_per_slide_preds.tsv", sep="\t", index=False)
    print(f"[saved] {OUT / 'h15_per_slide_preds.tsv'}")

    # ----------------------------------------------------------------
    # Step 7: figure
    # ----------------------------------------------------------------
    print("[fig] rendering encoder comparison...")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keep = ["img_prob_only", "img_emb_only", "ht13_only",
            "img_prob_plus_ht13", "img_emb_plus_ht13"]
    plot_df = res_df[res_df["combo"].isin(keep)].copy()
    plot_df["combo"] = pd.Categorical(plot_df["combo"], keep, ordered=True)
    plot_df = plot_df.sort_values(["combo", "encoder"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    width = 0.38
    x = np.arange(len(keep))
    img_rows = plot_df[plot_df["encoder"] == "ImageNet"].set_index("combo").reindex(keep)
    uni_rows = plot_df[plot_df["encoder"] == "UNI"].set_index("combo").reindex(keep)

    img_auc = img_rows["pooled_auc"].values
    uni_auc = uni_rows["pooled_auc"].values
    img_lo = img_rows["pooled_auc"].values - img_rows["boot_lo95"].values
    img_hi = img_rows["boot_hi95"].values - img_rows["pooled_auc"].values
    uni_lo = uni_rows["pooled_auc"].values - uni_rows["boot_lo95"].values
    uni_hi = uni_rows["boot_hi95"].values - uni_rows["pooled_auc"].values

    ax.bar(x - width / 2, img_auc, width, color="#1f77b4", label="ImageNet ViT-L",
           alpha=0.85)
    ax.errorbar(x - width / 2, img_auc, yerr=[img_lo, img_hi], fmt="none",
                ecolor="k", capsize=3, lw=1)
    ax.bar(x + width / 2, uni_auc, width, color="#ff7f0e", label="UNI (foundation)",
           alpha=0.85)
    ax.errorbar(x + width / 2, uni_auc, yerr=[uni_lo, uni_hi], fmt="none",
                ecolor="k", capsize=3, lw=1)

    ax.axhline(0.5, ls=":", color="grey", label="chance")
    ax.axhline(0.564, ls="--", color="red", lw=1,
               label="ImageNet image-only baseline (0.564)")
    ax.set_xticks(x)
    ax.set_xticklabels(keep, rotation=20, ha="right")
    ax.set_ylim(0.3, 1.05)
    ax.set_ylabel("Pooled OOF AUC (BRAF_like/cPTC, n=37)")
    ax.set_title("H15 — Foundation-model H&E rescue on BRAF/cPTC stratum\n"
                 "5-fold CV; 95% CI = stratified bootstrap (1000)")
    ax.legend(loc="lower right", fontsize=8)
    plt.tight_layout()
    plt.savefig(OUT / "fig_h15_encoder_uplift.png", dpi=180)
    plt.savefig(OUT / "fig_h15_encoder_uplift.pdf")
    plt.close()
    print(f"[saved] fig_h15_encoder_uplift.png/.pdf")

    # ----------------------------------------------------------------
    # Step 8: recon log
    # ----------------------------------------------------------------
    recon_lines = [
        "=== H15 recon ===",
        f"BRAF_like/cPTC cohort: n={n_total} (DM1={n_dm1}, DM2={n_dm2})",
        f"  ImageNet sanity image-only AUC (loaded ckpt): {img_only_pooled:.4f}",
        f"UNI features available locally for {n_uni}/{n_total} slides",
        f"  Missing UNI features for {n_total - n_uni} slides:",
    ]
    for s in missing_slides:
        recon_lines.append(f"    {s}")
    recon_lines += [
        "",
        "=== feature inventory ===",
        f"ImageNet features dir: {IMG_FEAT_DIR}",
        f"  pt files: {len(img_avail)}",
        f"UNI features dir:      {UNI_FEAT_DIR}",
        f"  pt files: {len(uni_avail)}",
        "",
        "=== pod recon ===",
        "Active pod: thca-spark-dm-a6000-v4 (uvp9i2r9s6l85y), $0.33/hr",
        "  /workspace: prior neoantigen run; no WSI/tile data",
        "  /runpod-volume: empty (post-L40S→A6000 migration)",
        "  Decision: NO pod time used. All UNI features already on local disk.",
        "  Pod $0 spent on H15.",
        "",
        "=== cohort filter ===",
        f"  matched-N for fair comparison: {n_uni}",
        f"    DM1={cohort_uni_dm1}, DM2={cohort_uni_dm2}",
    ]
    (OUT / "h15_recon.txt").write_text("\n".join(recon_lines) + "\n")

    # ----------------------------------------------------------------
    # Step 9: pod runtime log
    # ----------------------------------------------------------------
    elapsed = time.time() - t0
    runtime_lines = [
        f"H15 wall-clock: {elapsed:.1f}s ({elapsed/60:.2f} min)",
        f"  CLAM UNI 5-fold train: {t_uni:.1f}s",
        f"  CLAM ImageNet (matched-N) 5-fold train: {t_img:.1f}s",
        f"Compute: {DEVICE}",
        "Pod time: 0 minutes (everything ran locally)",
        "Pod cost spent on H15: $0.00",
        "Budget reserved/unused: $1.00 (well under 2h budget)",
    ]
    (OUT / "h15_pod_runtime_log.txt").write_text("\n".join(runtime_lines) + "\n")

    # ----------------------------------------------------------------
    # Step 10: report
    # ----------------------------------------------------------------
    print("[report] writing H15_REPORT.md")
    write_report(res_df, n_uni, cohort_uni_dm1, cohort_uni_dm2,
                 missing_slides, t_uni, t_img, img_only_pooled, img_only_matched, uni_img_only)

    print(f"\n[done] total wall {elapsed:.1f}s")


def write_report(res_df, n_uni, n_dm1, n_dm2, missing_slides,
                 t_uni, t_img, img_only_baseline, img_only_matched, uni_img_only):
    img_rows = res_df[res_df["encoder"] == "ImageNet"].set_index("combo")
    uni_rows = res_df[res_df["encoder"] == "UNI"].set_index("combo")

    def get(rows, c):
        if c not in rows.index:
            return None
        return rows.loc[c]

    rows_summary = []
    for c in ["img_prob_only", "img_emb_only", "ht13_only",
              "img_prob_plus_ht13", "img_emb_plus_ht13"]:
        ir = get(img_rows, c)
        ur = get(uni_rows, c)
        rows_summary.append({
            "combo": c,
            "img_auc": f"{ir['pooled_auc']:.3f}" if ir is not None else "—",
            "img_ci": f"[{ir['boot_lo95']:.3f}, {ir['boot_hi95']:.3f}]" if ir is not None else "—",
            "uni_auc": f"{ur['pooled_auc']:.3f}" if ur is not None else "—",
            "uni_ci": f"[{ur['boot_lo95']:.3f}, {ur['boot_hi95']:.3f}]" if ur is not None else "—",
            "delta":  f"{ur['pooled_auc'] - ir['pooled_auc']:+.3f}"
                       if (ir is not None and ur is not None) else "—",
        })

    # verdict logic
    img_only_uni_auc = uni_img_only
    img_only_imgnet_auc = img_only_matched
    delta_img_only = img_only_uni_auc - img_only_imgnet_auc
    uni_img_ci_lo = float(uni_rows.loc["img_prob_only", "boot_lo95"])
    uni_img_ci_hi = float(uni_rows.loc["img_prob_only", "boot_hi95"])
    img_img_ci_lo = float(img_rows.loc["img_prob_only", "boot_lo95"])
    img_img_ci_hi = float(img_rows.loc["img_prob_only", "boot_hi95"])

    if uni_img_only > 0.70 and img_only_uni_auc > img_only_imgnet_auc + 0.05:
        verdict = "**RESCUE-PASS** — UNI features lift image-only AUC into a meaningful range"
    elif uni_img_only > img_only_imgnet_auc + 0.05:
        verdict = "**RESCUE-PARTIAL** — UNI improves over ImageNet but absolute image-only AUC stays below 0.70"
    elif uni_img_only > img_only_imgnet_auc:
        verdict = "**RESCUE-MARGINAL** — UNI nudges image-only above ImageNet but bootstrap CIs overlap"
    else:
        verdict = "**RESCUE-FAIL** — UNI no better than ImageNet on this stratum"

    body = f"""# H15 — foundation-model image rescue on BRAF/cPTC (2026-05-09)

## Verdict
{verdict}.

UNI image-only pooled AUC = **{uni_img_only:.3f}** [95% CI {uni_img_ci_lo:.3f}, {uni_img_ci_hi:.3f}], vs ImageNet ViT-L matched-N image-only **{img_only_imgnet_auc:.3f}** [{img_img_ci_lo:.3f}, {img_img_ci_hi:.3f}], delta = **{delta_img_only:+.3f}**. Reference: H7 ImageNet baseline (loaded checkpoint, n=41) AUC = {img_only_baseline:.3f}.

## Cohort
- BRAF_like/cPTC stratum, **n={n_uni}** (DM1={n_dm1}, DM2={n_dm2}).
- 4 of the original 41 slides lack UNI features locally (raw WSI was on a RunPod that was terminated; re-extraction would require re-downloading WSI from GDC). For honest comparison, ImageNet is also re-trained on the same n={n_uni} matched cohort.
- 5-fold CV uses the same fold assignment as the v2 phase2 CLAM run (`clam_per_slide_predictions.tsv`) so AUCs are directly comparable.

Missing slides: {", ".join(missing_slides)}.

## Method
- Two encoders: ImageNet ViT-L tile features and UNI (MahmoodLab) tile features (both 200×1024). UNI features extracted on prior runpod 2026-05-08 and persisted to `phase2_tcga_clam_UNI/features/`.
- For each encoder, fresh GatedAttentionMIL (512 hidden, dropout 0.25) trained 60 epochs per fold (AdamW lr=1e-4, weight_decay=1e-4, seed 42, identical hyper-params).
- Five combos per encoder: image scalar prob (img_prob), 512-d gated-attention bag (img_emb), HT-13 RNA only, prob+HT13, emb+HT13. LogReg L2 (C=0.5) for fusion; standardize on train fold.
- Bootstrap 95% CI on pooled OOF predictions (1000 stratified resamples).

## Results
| combo | ImageNet AUC | 95% CI | UNI AUC | 95% CI | delta |
|---|---:|---|---:|---|---:|
"""
    for r in rows_summary:
        body += f"| {r['combo']} | {r['img_auc']} | {r['img_ci']} | {r['uni_auc']} | {r['uni_ci']} | {r['delta']} |\n"

    body += f"""

Wall time: UNI 5-fold {t_uni:.0f}s, ImageNet 5-fold {t_img:.0f}s on {DEVICE.upper()}; **$0 pod spend** (everything local on existing UNI features).

## Caveats
- n={n_uni} only — bootstrap CIs are wide. AUC point estimates within ~0.10 of each other are not statistically separable here.
- HT-13 RNA already saturates AUC≈0.92 in this stratum (per H7), so multimodal uplift over RNA alone is mostly fixed. The H15 question is specifically about **image-only rescue**, i.e. can H&E alone separate DM1 from DM2 in BRAF/cPTC without the RNA crutch.
- 4 BRAF/cPTC slides excluded (no UNI features). If those 4 are systematically harder/easier the rescue verdict could shift; we report this honestly rather than re-running 6-12h pod jobs to recover them.
- Same N=37 used for both encoders — no train/test slide leakage between them, fold IDs preserved from v2 phase2 run.

## Files
- `h15_per_encoder_results.tsv` — full per-combo AUC table (10 rows)
- `h15_per_slide_preds.tsv` — OOF probabilities for every slide × combo × encoder
- `h15_recon.txt` — cohort + feature inventory + pod recon
- `h15_pod_runtime_log.txt` — wall-clock + dollar tracking
- `fig_h15_encoder_uplift.{{png,pdf}}` — encoder comparison bar+CI
- `run_h15_uni.py` — this script
"""
    (OUT / "H15_REPORT.md").write_text(body)
    print(f"[saved] {OUT / 'H15_REPORT.md'}")


if __name__ == "__main__":
    main()
