#!/usr/bin/env python3
"""R17 × TERT promoter status interaction in TCGA THCA.

Joins R17 two-axis zone (from r17_per_sample_merged.tsv) with TERT promoter status
(from dm_master). Tests whether TERT+ enriches in any zone; runs Cox PFI within
zones if N permits.
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
R17 = ROOT / "project/results/r17_tcga_panel_d4p2_reconciliation/r17_per_sample_merged.tsv"
DMM = ROOT / "project/results/dark_matter_phase2/web/data/tcga_dm_master_with_pfi.tsv"
OUT = ROOT / "project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction"
OUT.mkdir(parents=True, exist_ok=True)


def zone_of(row):
    panel_silenced = row["DM_call"] == "DM1"
    ht_high = row["sig_score"] > 0   # d4p2 HT-axis high
    if panel_silenced and ht_high:
        return "dark-matter"
    if panel_silenced and not ht_high:
        return "BRAF-like"
    if (not panel_silenced) and ht_high:
        return "RAS-like"
    return "WT-like"


def main() -> None:
    r17 = pd.read_csv(R17, sep="\t")
    r17["zone"] = r17.apply(zone_of, axis=1)
    dmm = pd.read_csv(DMM, sep="\t")
    merged = r17.merge(dmm[["tcga_short", "tert_pos", "PFI", "PFI.time"]], on="tcga_short", how="left")
    merged["tert_known"] = merged["tert_pos"].notna()
    merged.to_csv(OUT / "r17_tert_per_sample.tsv", sep="\t", index=False)

    # per-zone TERT+ fraction (samples with known TERT only)
    sub = merged[merged["tert_known"]].copy()
    sub["tert_bool"] = sub["tert_pos"].astype(bool)

    rows = []
    for zone, g in sub.groupby("zone"):
        n = len(g)
        n_tert = int(g["tert_bool"].sum())
        rows.append({"zone": zone, "n": n, "n_TERT_pos": n_tert,
                     "pct_TERT_pos": round(100 * n_tert / n, 1) if n else 0.0})
    by_zone = pd.DataFrame(rows).sort_values("n", ascending=False)
    overall_pct = round(100 * sub["tert_bool"].mean(), 1)
    by_zone.to_csv(OUT / "r17_tert_fraction_by_zone.tsv", sep="\t", index=False)

    # fisher: each zone vs rest
    fisher_rows = []
    for zone in by_zone["zone"]:
        in_zone = sub["zone"] == zone
        a = int(((sub["tert_bool"]) & in_zone).sum())
        b = int(((~sub["tert_bool"]) & in_zone).sum())
        c = int(((sub["tert_bool"]) & ~in_zone).sum())
        d = int(((~sub["tert_bool"]) & ~in_zone).sum())
        odds, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        fisher_rows.append({"zone": zone, "TERT+_in_zone": a, "TERT-_in_zone": b,
                            "TERT+_rest": c, "TERT-_rest": d,
                            "OR": round(float(odds), 3), "p": float(f"{p:.3g}")})
    fisher_df = pd.DataFrame(fisher_rows)
    fisher_df.to_csv(OUT / "r17_tert_zone_fisher.tsv", sep="\t", index=False)

    # figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    zone_colors = {"BRAF-like": "#1f4f88", "RAS-like": "#d18b1f",
                   "dark-matter": "#7d3c98", "WT-like": "#888888"}

    ax = axes[0]
    by_zone_sorted = by_zone.sort_values("pct_TERT_pos", ascending=False)
    x = np.arange(len(by_zone_sorted))
    ax.bar(x, by_zone_sorted["pct_TERT_pos"], color=[zone_colors[z] for z in by_zone_sorted["zone"]])
    ax.axhline(overall_pct, color="#444", lw=1, ls="--", label=f"overall {overall_pct}%")
    for xi, z, n, p in zip(x, by_zone_sorted["zone"], by_zone_sorted["n"], by_zone_sorted["pct_TERT_pos"]):
        ax.text(xi, p + 1.0, f"{p:.1f}%\n(n={int(n)})", ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(by_zone_sorted["zone"], rotation=12, fontsize=10)
    ax.set_ylabel("% TERT promoter mutation +")
    ax.set_title(f"TERT+ fraction by R17 zone  (n={len(sub)} with known TERT)")
    ax.legend(fontsize=9, loc="upper right")
    ax.set_ylim(0, max(by_zone_sorted["pct_TERT_pos"].max() + 8, 30))

    # PFI by zone × TERT
    ax = axes[1]
    pfi_sub = sub.dropna(subset=["PFI", "PFI.time"]).copy()
    pfi_sub["zone_tert"] = pfi_sub["zone"] + " · TERT=" + pfi_sub["tert_bool"].map({True: "+", False: "−"})
    group_rates = pfi_sub.groupby(["zone", "tert_bool"], observed=True).agg(
        n=("PFI", "size"),
        events=("PFI", "sum"),
    )
    group_rates["event_rate_pct"] = (group_rates["events"] / group_rates["n"] * 100).round(1)
    group_rates.to_csv(OUT / "r17_pfi_by_zone_tert.tsv", sep="\t")

    # bar grouped
    zones_o = ["BRAF-like", "dark-matter", "RAS-like", "WT-like"]
    widths = 0.4
    for i, tert_state in enumerate([False, True]):
        ys = []
        ns = []
        for z in zones_o:
            try:
                rate = float(group_rates.loc[(z, tert_state), "event_rate_pct"])
                n = int(group_rates.loc[(z, tert_state), "n"])
            except KeyError:
                rate, n = 0.0, 0
            ys.append(rate); ns.append(n)
        offset = (-widths/2 if not tert_state else widths/2)
        bars = ax.bar(np.arange(len(zones_o)) + offset, ys, width=widths,
                      color="#999" if not tert_state else "#c0392b",
                      label=f"TERT={'−' if not tert_state else '+'}")
        for xi, y, n in zip(np.arange(len(zones_o)) + offset, ys, ns):
            if n > 0:
                ax.text(xi, y + 1, f"n={n}", ha="center", fontsize=8)
    ax.set_xticks(np.arange(len(zones_o)))
    ax.set_xticklabels(zones_o, rotation=12, fontsize=10)
    ax.set_ylabel("PFI event rate (%)")
    ax.set_title("PFI event rate by zone × TERT status")
    ax.legend(fontsize=9, loc="upper right")

    fig.suptitle("R17 × TERT  —  TERT+ enrichment and PFI events across two-axis zones (TCGA THCA)",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT / "fig_r17_tert_zone_interaction.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT / "fig_r17_tert_zone_interaction.pdf", bbox_inches="tight")
    plt.close(fig)

    print("Per-zone TERT+:")
    print(by_zone)
    print("\nFisher:")
    print(fisher_df)
    print("\nGroup rates (zone × TERT):")
    print(group_rates)

    # markdown
    lines = [
        "# R17 × TERT — TCGA THCA two-axis zone × TERT promoter status interaction\n",
        f"Samples with known TERT status: **n = {len(sub)}** of {len(merged)} R17-classified TCGA samples.",
        f"Overall TERT+ rate: **{overall_pct}%**.\n",
        "## TERT+ fraction by R17 zone\n",
        by_zone.to_markdown(index=False),
        "\n## Fisher exact (zone vs rest)\n",
        fisher_df.to_markdown(index=False),
        "\n## PFI event rate by zone × TERT status\n",
        group_rates.reset_index().to_markdown(index=False),
        "\n## Interpretation",
        "* TERT+ is an established late event in dedifferentiation. The R17 zones predict the histologic",
        "  / molecular state along two orthogonal axes (RAI silencing × HT overlap), so the question is",
        "  whether TERT+ enriches in any single zone.",
        "* The dark-matter zone (silenced + HT) is the predicted high-aggression cell-state intersection.",
        "  TERT+ enrichment here would tighten Paper 1's HR=7.57 TERT survival claim into a specific",
        "  cellular context.",
        "* PFI event rate by zone × TERT row shows the joint effect; n is small in some cells so treat",
        "  as exploratory.",
        "\n## Caveats",
        "* TCGA THCA is overwhelmingly low-risk PTC, so TERT+ rate is low overall (~10-13%); large-effect",
        "  zone enrichment may be hard to detect with this prior.",
        "* PFI event count is small per zone × TERT cell; survival inference within zones is underpowered.",
    ]
    (OUT / "R17_tert_zone_interaction.md").write_text("\n".join(lines))

    print("\nWrote:")
    for f in sorted(OUT.iterdir()):
        print(f" - {f.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
