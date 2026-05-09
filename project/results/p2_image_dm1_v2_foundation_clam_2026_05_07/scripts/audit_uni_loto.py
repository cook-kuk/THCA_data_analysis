#!/usr/bin/env python3
"""UNI-final leave-one-TSS-out audit.

This directly stress-tests the 54-slide UNI-final feature set under
center holdout. It is the missing counterpart to the ViT-L/phase2 LOTO audit.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).parent))
from audit_clam_train_lib import load_features, train_one_fold  # noqa: E402


ROOT = Path(os.environ.get(
    "THCA_P2_ROOT",
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_image_dm1_v2_foundation_clam_2026_05_07",
))
META = Path(os.environ.get(
    "THCA_META_PATH",
    "/home/seungho/personal/THCA_data_analysis/project/metadata/"
    "sample_master_v3.tsv",
))
FEAT_DIR = ROOT / "phase2_tcga_clam_UNI" / "features"
MANIFEST = ROOT / "phase2_tcga_clam" / "slide_manifest.tsv"
PRED = ROOT / "phase2_tcga_clam_UNI" / "clam_per_slide_predictions.tsv"
AUD = ROOT / "analysis_supp" / "audit_uni_loto"
DEVICE = os.environ.get("THCA_DEVICE", "cpu")
EPOCHS = int(os.environ.get("THCA_EPOCHS", "30"))


def load_dataset():
    man = pd.read_csv(MANIFEST, sep="\t")
    man = man[man["dm"].isin(["DM1", "DM2"])].copy()
    man["has_feat"] = man["file_id"].apply(
        lambda f: (FEAT_DIR / f"{f}.pt").exists())
    man = man[man["has_feat"]].reset_index(drop=True)
    man["label"] = man["dm"].map({"DM2": 0, "DM1": 1}).astype(int)
    man["tss"] = man["submitter_id"].str.extract(r"TCGA-([^-]+)-")

    pred = pd.read_csv(PRED, sep="\t")
    pred_slides = set(pred["slide"])
    man = man[man["file_id"].isin(pred_slides)].reset_index(drop=True)

    meta = pd.read_csv(META, sep="\t", low_memory=False)
    meta = meta[meta["dataset"] == "TCGA-THCA"].copy()
    meta["case_short"] = meta["sample_id"].str.extract(r"^(TCGA-[^-]+-[^-]+)")
    meta_case = (meta.sort_values("sample_id")
                     .drop_duplicates("case_short", keep="first")
                     [["case_short", "histology_subtype", "molecular_subtype",
                       "sex", "age", "ajcc_stage_group"]])
    man = man.merge(meta_case, left_on="submitter_id",
                    right_on="case_short", how="left")
    bags = load_features(FEAT_DIR, man["file_id"].tolist())
    labels = man["label"].to_numpy()
    return man, bags, labels, pred


def safe_auc(y, p):
    if len(set(y)) < 2:
        return np.nan
    return float(roc_auc_score(y, p))


def leave_one_tss_out(man, bags, labels):
    tss_list = sorted(man["tss"].dropna().unique())
    oof = np.full(len(labels), np.nan)
    rows = []
    for tss in tss_list:
        va_idx = np.where(man["tss"].to_numpy() == tss)[0]
        tr_idx = np.where(man["tss"].to_numpy() != tss)[0]
        tr_bags = [bags[i] for i in tr_idx]
        va_bags = [bags[i] for i in va_idx]
        tr_lab = [int(labels[i]) for i in tr_idx]
        va_lab = [int(labels[i]) for i in va_idx]
        t0 = time.time()
        _, probs = train_one_fold(
            tr_bags, tr_lab, va_bags, va_lab,
            epochs=EPOCHS, device=DEVICE, seed=42,
        )
        for j, idx in enumerate(va_idx):
            oof[idx] = probs[j]
        ras_mask = man.iloc[va_idx]["molecular_subtype"].eq("RAS_like")
        rows.append({
            "tss_held_out": tss,
            "n_test": int(len(va_idx)),
            "n_pos": int(sum(va_lab)),
            "n_ras_test": int(ras_mask.sum()),
            "n_ras_pos_test": int(labels[va_idx][ras_mask.to_numpy()].sum()),
            "tss_holdout_auc": safe_auc(va_lab, probs),
            "seconds": round(time.time() - t0, 1),
        })
        print(
            f"TSS={tss}: n={len(va_idx)} pos={sum(va_lab)} "
            f"RAS={int(ras_mask.sum())} auc={rows[-1]['tss_holdout_auc']}"
        )
    return pd.DataFrame(rows), oof


def make_figure(summary: dict, per_tss: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    labels = ["Original UNI OOF", "UNI LOTO pooled", "UNI LOTO RAS-like"]
    vals = [
        summary["original_uni_overall_auc"],
        summary["pooled_overall_auc"],
        summary["pooled_ras_like_auc"],
    ]
    colors = ["#3b82f6", "#f59e0b", "#ef4444"]
    axes[0].bar(labels, vals, color=colors)
    axes[0].axhline(0.5, color="#333333", lw=1, ls="--")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_ylabel("AUC")
    axes[0].set_title("UNI OOF vs center holdout")
    axes[0].tick_params(axis="x", rotation=25)
    for i, v in enumerate(vals):
        axes[0].text(i, v + 0.025, f"{v:.3f}", ha="center", fontsize=10)

    plot_df = per_tss.copy()
    plot_df["auc_plot"] = plot_df["tss_holdout_auc"].astype(float)
    axes[1].bar(plot_df["tss_held_out"], plot_df["auc_plot"],
                color="#64748b")
    axes[1].axhline(0.5, color="#333333", lw=1, ls="--")
    axes[1].set_ylim(0, 1.05)
    axes[1].set_ylabel("AUC")
    axes[1].set_title("Per-TSS held-out AUC")
    for i, row in enumerate(plot_df.itertuples()):
        if not np.isnan(row.auc_plot):
            axes[1].text(i, row.auc_plot + 0.025, f"{row.auc_plot:.2f}",
                         ha="center", fontsize=8)

    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("UNI-final LOTO audit: center-holdout stress test",
                 fontsize=13, fontweight="bold")
    fig.savefig(AUD / "fig_H_uni_loto.png", dpi=220)
    fig.savefig(AUD / "fig_H_uni_loto.pdf")
    plt.close(fig)


def main() -> None:
    AUD.mkdir(parents=True, exist_ok=True)
    man, bags, labels, pred = load_dataset()
    print(f"[load] UNI-final dataset: {len(man)} slides, "
          f"{int(labels.sum())} DM1 / {int(len(labels)-labels.sum())} DM2")

    per_tss, oof = leave_one_tss_out(man, bags, labels)
    out_pred = man[["file_id", "submitter_id", "label", "tss",
                    "histology_subtype", "molecular_subtype", "sex"]].copy()
    out_pred["prob_DM1_uni_loto"] = oof
    out_pred.to_csv(AUD / "uni_loto_oof_predictions.tsv",
                    sep="\t", index=False)
    per_tss.to_csv(AUD / "e8_uni_leave_one_tss_out.tsv",
                   sep="\t", index=False)

    original_auc = float(roc_auc_score(pred["label"], pred["prob_DM1"]))
    ras_mask = man["molecular_subtype"].eq("RAS_like").to_numpy()
    braf_mask = man["molecular_subtype"].eq("BRAF_like").to_numpy()
    female_mask = man["sex"].eq("Female").to_numpy()
    male_mask = man["sex"].eq("Male").to_numpy()
    summary = {
        "n_slides": int(len(man)),
        "n_pos": int(labels.sum()),
        "n_neg": int(len(labels) - labels.sum()),
        "original_uni_overall_auc": original_auc,
        "pooled_overall_auc": safe_auc(labels, oof),
        "pooled_ras_like_auc": safe_auc(labels[ras_mask], oof[ras_mask]),
        "pooled_braf_like_auc": safe_auc(labels[braf_mask], oof[braf_mask]),
        "pooled_female_auc": safe_auc(labels[female_mask], oof[female_mask]),
        "pooled_male_auc": safe_auc(labels[male_mask], oof[male_mask]),
        "per_tss": per_tss.to_dict("records"),
        "epochs": EPOCHS,
        "device": DEVICE,
        "feature_dir": str(FEAT_DIR),
    }
    (AUD / "UNI_LOTO_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, default=str)
    )
    make_figure(summary, per_tss)

    md = [
        "# UNI-final leave-one-TSS-out audit",
        "",
        "| Metric | Value |",
        "|---|---:|",
        f"| n slides | {summary['n_slides']} |",
        f"| Original UNI OOF AUC | {summary['original_uni_overall_auc']:.3f} |",
        f"| UNI LOTO pooled overall AUC | {summary['pooled_overall_auc']:.3f} |",
        f"| UNI LOTO pooled RAS-like AUC | {summary['pooled_ras_like_auc']:.3f} |",
        f"| UNI LOTO pooled BRAF-like AUC | {summary['pooled_braf_like_auc']:.3f} |",
        f"| UNI LOTO Female AUC | {summary['pooled_female_auc']:.3f} |",
        f"| UNI LOTO Male AUC | {summary['pooled_male_auc']:.3f} |",
        "",
        "## Per-TSS held-out folds",
        "",
        per_tss.to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "This is the center-holdout counterpart to the UNI-final OOF result. "
        "It should be cited before making any UNI image-DM1 claim from TCGA.",
        "",
    ]
    (AUD / "UNI_LOTO_REPORT.md").write_text("\n".join(md))
    print(f"[done] {AUD / 'UNI_LOTO_SUMMARY.json'}")
    print(f"[done] {AUD / 'UNI_LOTO_REPORT.md'}")
    print(f"[done] {AUD / 'fig_H_uni_loto.png'}")


if __name__ == "__main__":
    main()
