"""
Phase 1b — UNI embedding ↔ DM1 score per-slide correlation + spatial overlay figure.

Steps:
  1. Load UNI embeddings (per-spot, 1024-d) + spot_scores (DM1_like, RAI_8, etc).
  2. Per slide: PCA on UNI embeddings → top 3 PCs.
  3. Spearman ρ (PC ↔ DM1_like_score) per slide.
  4. Spatial overlay: tile coordinates × DM1 score / PC1 / PC2 maps.

Output:
  uni_dm1_correlation_per_slide.tsv
  uni_dm1_overlay_per_slide_*.png  (spatial heatmap per slide)
  PHASE1_REPORT.md (auto-generated)
"""
from __future__ import annotations
import argparse, os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import spearmanr

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--embed_npz", type=Path,
                        default=Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/uni_embeddings_size224.npz"))
    parser.add_argument("--meta_tsv", type=Path,
                        default=Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521/uni_embed_metadata_size224.tsv"))
    parser.add_argument("--scores_tsv", type=Path,
                        default=Path("project/results/01_spatial_score/all_spots_scored.tsv.gz"))
    parser.add_argument("--out_dir", type=Path,
                        default=Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase1_gse250521"))
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # --- load ---
    embeds = np.load(args.embed_npz)["embeddings"]
    meta = pd.read_csv(args.meta_tsv, sep="\t")
    scores = pd.read_csv(args.scores_tsv, sep="\t")
    print(f"[load] embeds {embeds.shape}, meta {meta.shape}, scores {scores.shape}")

    # --- merge (slide+spot_id key) ---
    meta["key"] = meta["slide"].astype(str) + "_" + meta["spot_id"].astype(str)
    scores["key"] = scores["sample_id"].astype(str) + "_" + scores["spot_id"].astype(str)
    merged = meta.merge(scores[["key", "DM1_like_score", "RAI_8_score", "TDS_like_score",
                                  "pxl_row_in_fullres", "pxl_col_in_fullres",
                                  "stage", "stage_metadata_raw", "in_tissue"]],
                         on="key", how="left")
    merged["embed_idx"] = range(len(merged))
    in_tissue = merged["in_tissue"] == 1
    print(f"[merge] in_tissue n={in_tissue.sum()}")

    # --- per slide: PCA + corr + figure ---
    rows = []
    slides = merged["slide"].dropna().unique()
    for slide in slides:
        m = merged[(merged["slide"] == slide) & in_tissue & merged["DM1_like_score"].notna()]
        if len(m) < 50:
            print(f"  [skip] {slide} n={len(m)} too few")
            continue
        e = embeds[m["embed_idx"].values]
        pca = PCA(n_components=3)
        pcs = pca.fit_transform(e)
        for i in range(3):
            rho, p = spearmanr(pcs[:, i], m["DM1_like_score"])
            rows.append({"slide": slide, "n_spots": len(m),
                         "pc": i+1, "explained_var": float(pca.explained_variance_ratio_[i]),
                         "rho_DM1": float(rho), "p_DM1": float(p),
                         "stage": m["stage"].iloc[0]})

        # spatial overlay
        fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
        x = m["pxl_col_in_fullres"].values
        y = -m["pxl_row_in_fullres"].values
        for ax, vals, title in zip(axes,
                                    [m["DM1_like_score"].values, pcs[:, 0], pcs[:, 1]],
                                    [f"DM1_like_score (n={len(m)})",
                                     f"UNI PC1 (var {pca.explained_variance_ratio_[0]:.2%})",
                                     f"UNI PC2 (var {pca.explained_variance_ratio_[1]:.2%})"]):
            sc = ax.scatter(x, y, c=vals, cmap="RdBu_r", s=15, edgecolor="none")
            plt.colorbar(sc, ax=ax, fraction=0.04)
            ax.set_aspect("equal"); ax.set_title(title); ax.axis("off")
        rho1 = next((r["rho_DM1"] for r in rows if r["slide"] == slide and r["pc"] == 1), float("nan"))
        rho2 = next((r["rho_DM1"] for r in rows if r["slide"] == slide and r["pc"] == 2), float("nan"))
        plt.suptitle(f"{slide} (stage={m['stage'].iloc[0]}) — Spearman PC1↔DM1 ρ={rho1:.2f} | PC2↔DM1 ρ={rho2:.2f}",
                     fontsize=12, fontweight="bold")
        plt.tight_layout()
        plt.savefig(args.out_dir / f"uni_dm1_overlay_{slide}.png", dpi=140, bbox_inches="tight")
        plt.close()

    df = pd.DataFrame(rows)
    df.to_csv(args.out_dir / "uni_dm1_correlation_per_slide.tsv", sep="\t", index=False)
    print(df.to_string(index=False))

    # --- kill switch eval ---
    pc1 = df[df["pc"] == 1]
    pass_count = (pc1["rho_DM1"].abs() > 0.3).sum()
    abs_max = pc1["rho_DM1"].abs().max()
    n_slides = len(pc1)
    pass_rate = pass_count / max(n_slides, 1)

    verdict = ("PASS" if pass_rate >= 0.5 else
               ("MARGINAL" if pass_rate >= 0.25 else "FAIL"))

    rep = []
    rep.append(f"# Phase 1 report — UNI ↔ DM1 spatial correlation\n")
    rep.append(f"- Slides: {n_slides}")
    rep.append(f"- |ρ|>0.3 in PC1: {pass_count}/{n_slides} ({pass_rate:.0%})")
    rep.append(f"- max |ρ|: {abs_max:.3f}")
    rep.append(f"- **VERDICT: {verdict}**")
    rep.append("")
    rep.append("- Kill-switch criterion: PASS ≥ 50% slides with |ρ| > 0.3 in PC1")
    rep.append("")
    rep.append("## Per-slide ρ (PC1)")
    rep.append(pc1[["slide","stage","n_spots","rho_DM1","p_DM1"]].to_string(index=False))
    (args.out_dir / "PHASE1_REPORT.md").write_text("\n".join(rep))
    print("\n" + "\n".join(rep))


if __name__ == "__main__":
    main()
