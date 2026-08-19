#!/usr/bin/env python3
"""R17 master synthesis figure — 8 layers, single multi-panel visualization."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
R17 = ROOT / "project/results/r17_tcga_panel_d4p2_reconciliation"
OUT = R17 / "master_synthesis"
OUT.mkdir(parents=True, exist_ok=True)

ZONE_COLORS = {"BRAF-like": "#1f4f88", "RAS-like": "#d18b1f",
               "dark-matter": "#7d3c98", "WT-like": "#888888"}


def main() -> None:
    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(4, 4, figure=fig, hspace=0.55, wspace=0.45)

    # === Panel A: L1 bulk TCGA scatter ===
    ax = fig.add_subplot(gs[0, 0:2])
    L1 = pd.read_csv(R17 / "r17_per_sample_merged.tsv", sep="\t")
    driver_colors = {"BRAF V600E": "#1f4f88", "RAS": "#d18b1f",
                     "BRAF·RAS-neg": "#2c8a3a", "Other/unspecified": "#888"}
    for d, sub in L1.groupby("driver"):
        ax.scatter(sub["RAI_8"], sub["sig_score"], c=driver_colors.get(d, "#888"),
                   s=14, alpha=0.65, edgecolors="white", linewidths=0.3,
                   label=f"{d} (n={len(sub)})")
    ax.axhline(0, color="#444", lw=0.5, ls="--"); ax.axvline(0, color="#444", lw=0.5, ls="--")
    ax.set_xlabel("Panel RAI_8 z"); ax.set_ylabel("d4p2 sig_score")
    ax.set_title("A · L1 TCGA bulk RNA  n=500  concordance 14.8%", fontsize=10.5)
    ax.legend(fontsize=7, loc="lower right")

    # === Panel B: L2 Lu 2023 sc zone fractions ===
    ax = fig.add_subplot(gs[0, 2])
    L2 = pd.read_csv(R17 / "sc/lu2023_sc_zone_fraction_by_histology.tsv", sep="\t").set_index("histology")
    L2 = L2[["BRAF-like zone (silenced, no HT)", "RAS-like zone (preserved + HT)",
             "dark-matter zone (silenced + HT)", "wild-type-like (preserved, no HT)"]]
    L2.columns = ["BRAF-like", "RAS-like", "dark-matter", "WT-like"]
    L2.plot(kind="bar", stacked=True, ax=ax,
            color=[ZONE_COLORS[c] for c in L2.columns], legend=False)
    ax.set_ylabel("% malignant cells"); ax.set_xlabel("")
    ax.set_xticklabels(L2.index, rotation=0, fontsize=9)
    ax.set_title("B · L2 Lu 2023 sc  n=14,624", fontsize=10.5)

    # === Panel C: L3 Korean GSE286332 ===
    ax = fig.add_subplot(gs[0, 3])
    L3 = pd.read_csv(R17 / "cross_ethnic/gse286332_zone_fraction_by_group.tsv", sep="\t").set_index("group")
    L3_cols = []
    for c in ["BRAF-like", "RAS-like (preserved+HT)", "dark-matter (silenced+HT)", "WT-like (preserved,no HT)"]:
        for col in L3.columns:
            if col.startswith(c.split(" ")[0]):
                L3_cols.append(col); break
    L3 = L3[[c for c in L3.columns if c in L3_cols]]
    L3.columns = [c.split(" ")[0] for c in L3.columns]
    L3.plot(kind="bar", stacked=True, ax=ax,
            color=[ZONE_COLORS.get(c, "#888") if c in ZONE_COLORS else "#888"
                   for c in [{"BRAF-like":"BRAF-like","RAS-like":"RAS-like","dark-matter":"dark-matter","WT-like":"WT-like"}.get(x, x) for x in L3.columns]], legend=False)
    ax.set_ylabel("% samples"); ax.set_xlabel("")
    ax.set_xticklabels(L3.index, rotation=0, fontsize=9)
    ax.set_title("C · L3 Korean GSE286332  n=18", fontsize=10.5)

    # === Panel D: L4 Mun protein zone fractions ===
    ax = fig.add_subplot(gs[1, 0:2])
    L4 = pd.read_csv(R17 / "proteome_mun2025/mun2025_zone_fraction_by_group.tsv", sep="\t").set_index("group")
    L4 = L4.reindex(["PTC", "PDTC", "ATC"])
    L4 = L4[["BRAF-like", "RAS-like", "dark-matter", "WT-like"]]
    L4.plot(kind="bar", stacked=True, ax=ax,
            color=[ZONE_COLORS[c] for c in L4.columns], legend=False)
    ax.set_ylabel("% samples"); ax.set_xlabel("")
    ax.set_xticklabels(L4.index, rotation=0, fontsize=10)
    ax.set_title("D · L4 Mun 2025 proteome  n=336  ·  ATC vs PTC dark-matter Fisher OR=8.54  p=2.7e-15", fontsize=10.5)

    # === Panel E: L5 TERT enrichment by zone ===
    ax = fig.add_subplot(gs[1, 2:4])
    L5 = pd.read_csv(R17 / "tert_interaction/r17_tert_fraction_by_zone.tsv", sep="\t")
    L5 = L5.sort_values("zone", key=lambda s: s.map({"BRAF-like":0, "dark-matter":1, "WT-like":2, "RAS-like":3}))
    x = np.arange(len(L5))
    ax.bar(x, L5["pct_TERT_pos"], color=[ZONE_COLORS[z] for z in L5["zone"]])
    for xi, v, n, z in zip(x, L5["pct_TERT_pos"], L5["n"], L5["zone"]):
        ax.text(xi, v + 0.8, f"{v:.1f}%\n(n={int(n)})", ha="center", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(L5["zone"], rotation=10, fontsize=10)
    ax.set_ylabel("% TERT promoter mutation +")
    ax.set_title("E · L5 TCGA TERT+ by zone  ·  dark-matter OR=2.34  p=0.016", fontsize=10.5)
    ax.set_ylim(0, max(L5["pct_TERT_pos"].max() + 5, 15))

    # === Panel F: L6 Cox forest (PFI only) ===
    ax = fig.add_subplot(gs[2, 0:2])
    L6 = pd.read_csv(R17 / "zone_survival/r17_zone_cox_results.tsv", sep="\t")
    L6 = L6[(L6["model"] == "univariate") & (L6["outcome"] == "PFI") & (L6["term"].str.startswith("zone_"))]
    term_to_label = {"zone_BRAF": "BRAF-like", "zone_dark": "dark-matter", "zone_RAS": "RAS-like"}
    term_to_color = {"zone_BRAF": ZONE_COLORS["BRAF-like"], "zone_dark": ZONE_COLORS["dark-matter"], "zone_RAS": ZONE_COLORS["RAS-like"]}
    ax.axvline(1, color="#444", lw=0.6, ls="--")
    for i, (_, row) in enumerate(L6.iterrows()):
        hr, lo, hi, p = row["HR"], row["HR_lower_95"], row["HR_upper_95"], row["p"]
        ax.errorbar(hr, i, xerr=[[hr - lo], [hi - hr]],
                    fmt='o', color=term_to_color[row["term"]], markersize=10, capsize=4)
        ax.text(min(hi * 1.05, 10), i, f"  HR {hr:.2f} [{lo:.2f}-{hi:.2f}]  p={p:.3f}",
                va="center", fontsize=8.5)
    ax.set_yticks(range(len(L6))); ax.set_yticklabels([term_to_label[t] for t in L6["term"]])
    ax.set_xscale("log"); ax.set_xlim(0.1, 12)
    ax.set_xlabel("PFI HR (ref = WT-like)")
    ax.set_title("F · L6 TCGA Cox PFI  n=563, 63 events  ·  BRAF/dark-matter HR=1.89 p≈0.05", fontsize=10.5)

    # === Panel G: L7 RNA z heatmap ===
    ax = fig.add_subplot(gs[2, 2:4])
    L7 = pd.read_csv(R17 / "gene_profile/r17_zone_gene_z_means.tsv", sep="\t", index_col=0)
    L7 = L7[["BRAF-like", "RAS-like", "dark-matter", "WT-like"]]
    im = ax.imshow(L7.values, cmap="RdBu_r", vmin=-1.2, vmax=1.2, aspect="auto")
    ax.set_xticks(range(L7.shape[1])); ax.set_xticklabels(L7.columns, rotation=10, fontsize=10)
    ax.set_yticks(range(L7.shape[0])); ax.set_yticklabels(L7.index, fontsize=9)
    for i in range(L7.shape[0]):
        for j in range(L7.shape[1]):
            v = L7.values[i, j]
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=8,
                    color="white" if abs(v) > 0.8 else "black")
    fig.colorbar(im, ax=ax, label="RNA z", fraction=0.046, pad=0.04)
    ax.set_title("G · L7 RNA per-zone 8-gene z  n=522", fontsize=10.5)

    # === Panel H: L8 HM450 β heatmap ===
    ax = fig.add_subplot(gs[3, 0:2])
    L8 = pd.read_csv(R17 / "hm450_beta/r17_zone_hm450_beta_means.tsv", sep="\t", index_col=0)
    L8 = L8[["BRAF-like", "RAS-like", "dark-matter", "WT-like"]]
    im2 = ax.imshow(L8.values, cmap="Reds", vmin=0, vmax=0.9, aspect="auto")
    ax.set_xticks(range(L8.shape[1])); ax.set_xticklabels(L8.columns, rotation=10, fontsize=10)
    ax.set_yticks(range(L8.shape[0])); ax.set_yticklabels(L8.index, fontsize=9)
    for i in range(L8.shape[0]):
        for j in range(L8.shape[1]):
            v = L8.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8,
                    color="white" if v > 0.5 else "black")
    fig.colorbar(im2, ax=ax, label="HM450 β", fraction=0.046, pad=0.04)
    ax.set_title("H · L8 HM450 β per zone  n=518  ·  dark-matter mean_8g β 0.41 vs WT 0.28", fontsize=10.5)

    # === Panel I: synthesis text panel ===
    ax = fig.add_subplot(gs[3, 2:4])
    ax.axis("off")
    txt = (
        "R17 synthesis · 8 layers · 4 pillars (RNA · sc · protein · methylation)\n\n"
        "• Two-axis decomposition: HT-overlap (d4p2 sig.DM) vs RAI-lineage (8-gene panel)\n"
        "  separates BRAF-like (silenced RAI, no HT), RAS-like (preserved + HT),\n"
        "  dark-matter (silenced + HT), WT-like (neither) zones.\n\n"
        "• 14.8% RNA concordance is driver-pattern biology, not a labeling bug.\n\n"
        "• Dark-matter zone is uniquely:\n"
        "    – the cellular substrate at sc level (38.3% of ATC malignant cells)\n"
        "    – the dominant ATC zone at protein level (58.4%, OR=8.54 p=2.7e-15)\n"
        "    – enriched for TERT+ (OR=2.34 p=0.016) → cellular context for HR=7.57\n"
        "    – the most hypermethylated (mean_8g β 0.41 vs WT 0.28)\n"
        "    – the deepest RNA silencing (mean z = -0.54)\n"
        "    – borderline-prognostic (PFI HR=1.89 p=0.050)\n\n"
        "• BRAF-like vs dark-matter use different gene-level silencing programs\n"
        "  (TPO/DIO1 vs TG/PAX8/TSHR) — methylation-mediated for half the panel.\n\n"
        "Voice-protected manuscript sections untouched. Compiled 2026-05-21."
    )
    ax.text(0.0, 1.0, txt, fontsize=9.5, va="top", ha="left", family="monospace",
            transform=ax.transAxes)

    fig.suptitle(
        "R17 two-axis reconciliation · 8-layer 4-pillar synthesis (TCGA + Lu sc + Korean + Mun protein + TERT + Cox + gene + HM450 β)",
        fontsize=13, y=0.995,
    )
    fig.savefig(OUT / "fig_r17_master_synthesis.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT / "fig_r17_master_synthesis.pdf", bbox_inches="tight")
    plt.close(fig)

    print("Wrote:")
    for f in sorted(OUT.iterdir()):
        print(f" - {f.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
