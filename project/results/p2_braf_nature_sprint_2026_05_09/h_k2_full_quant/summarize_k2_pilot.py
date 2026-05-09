#!/usr/bin/env python3
"""Summarize the completed n=9 K2 kallisto pilot.

This is a reviewer-reserve audit, not a full PRJEB11591 result: only the local
FASTQ files were quantified.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats


OUT = Path(__file__).resolve().parent
SCORES = OUT / "k2_full_panel_scores.tsv"
AUCS = OUT / "k2_panel_aucs.tsv"
SUMMARY = OUT / "k2_full_summary.json"


def fmt(x: float | int | None, nd: int = 3) -> str:
    if x is None or pd.isna(x):
        return "NA"
    return f"{float(x):.{nd}f}"


def main() -> None:
    scores = pd.read_csv(SCORES, sep="\t")
    aucs = pd.read_csv(AUCS, sep="\t")
    summary = json.loads(SUMMARY.read_text())

    panel_cols = ["HT_13", "FA_12", "MAPK_9"]
    scores["pair"] = scores["sample_alias"].str.replace("-N$", "", regex=True)

    paired_records = []
    for pair, group in scores.groupby("pair"):
        tumor = group[group["tissue"] == "tumor"]
        normal = group[group["tissue"] == "normal_matched"]
        if len(tumor) != 1 or len(normal) != 1:
            continue
        rec = {"pair": pair}
        for panel in panel_cols:
            rec[f"{panel}_tumor"] = float(tumor[panel].iloc[0])
            rec[f"{panel}_normal"] = float(normal[panel].iloc[0])
            rec[f"{panel}_tumor_minus_normal"] = float(
                tumor[panel].iloc[0] - normal[panel].iloc[0]
            )
        paired_records.append(rec)
    paired = pd.DataFrame(paired_records)
    paired_path = OUT / "k2_paired_tumor_normal_deltas.tsv"
    paired.to_csv(paired_path, sep="\t", index=False)

    # Small n: report descriptive paired deltas only.
    paired_summary = []
    for panel in panel_cols:
        col = f"{panel}_tumor_minus_normal"
        vals = paired[col].dropna() if col in paired else pd.Series(dtype=float)
        if len(vals) >= 2:
            stat, p = stats.wilcoxon(vals, alternative="two-sided")
        else:
            stat, p = None, None
        paired_summary.append(
            {
                "panel": panel,
                "n_pairs": int(len(vals)),
                "median_tumor_minus_normal": float(vals.median()) if len(vals) else None,
                "mean_tumor_minus_normal": float(vals.mean()) if len(vals) else None,
                "wilcoxon_p": p,
            }
        )
    paired_summary_df = pd.DataFrame(paired_summary)
    paired_summary_path = OUT / "k2_paired_tumor_normal_summary.tsv"
    paired_summary_df.to_csv(paired_summary_path, sep="\t", index=False)

    corr = scores[panel_cols + ["p_DM2"]].corr(method="spearman")
    corr_path = OUT / "k2_panel_spearman.tsv"
    corr.to_csv(corr_path, sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), constrained_layout=True)
    tissue_color = {"normal_matched": "#8fb3ff", "tumor": "#f2b36d"}
    for ax, panel in zip(axes, panel_cols, strict=True):
        for pair, group in scores.groupby("pair"):
            group = group.sort_values("tissue")
            x = [0 if t == "normal_matched" else 1 for t in group["tissue"]]
            ax.plot(x, group[panel], color="#3b465c", alpha=0.5, linewidth=1)
            ax.scatter(
                x,
                group[panel],
                c=[tissue_color[t] for t in group["tissue"]],
                s=34,
                edgecolor="#111827",
                linewidth=0.4,
                zorder=3,
            )
        ax.axhline(0, color="#7b8497", linewidth=0.8, linestyle="--")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["normal", "tumor"], fontsize=8)
        ax.set_title(panel.replace("_", "-"), fontsize=11)
        ax.set_ylabel("within-K2 mean z" if panel == "HT_13" else "")
        ax.grid(axis="y", color="#e7ebf3", linewidth=0.6)
    fig.suptitle(
        "K2 PRJEB11591 local full-transcriptome pilot: n=9, all DM2 calls",
        fontsize=12,
    )
    fig_path = OUT / "fig_k2_pilot_panel_scores.png"
    fig_pdf = OUT / "fig_k2_pilot_panel_scores.pdf"
    fig.savefig(fig_path, dpi=220)
    fig.savefig(fig_pdf)
    plt.close(fig)

    tumor_rows = aucs[aucs["contrast"] == "tumor_vs_normal"].copy()
    tumor_rows = tumor_rows[["panel", "auc", "mw_p"]].to_dict("records")

    md = f"""# K2 full-transcriptome pilot verdict

Build: 2026-05-09 KST.

## Verdict

**Not a breakthrough result. Use as reviewer-reserve / negative-control only.**

The local kallisto run completed successfully for all available FASTQs, but the
available local data are only **n=9**: **5 tumors and 4 matched normals**. All
9 samples inherit the prior 8-gene mini-index call **DM2**, so the completed
full-transcriptome pilot cannot test DM1-vs-DM2 biology in K2.

## What succeeded

- Quantification completed: {summary["n_samples_quanted"]}/9 local FASTQ runs.
- Gene-level matrix: {summary["n_genes_aggregated"]:,} gene symbols.
- Panel coverage: HT-13 {summary["panel_coverage"]["HT_13"]["present"]}/13,
  FA-12 {summary["panel_coverage"]["FA_12"]["present"]}/12, MAPK-9
  {summary["panel_coverage"]["MAPK_9"]["present"]}/9.
- kallisto index: `{summary["kallisto_index"]}`.

## What it says

Tumor-vs-normal signal is weak and not deployable at n=9:

| Panel | AUC tumor vs normal | MW p |
|---|---:|---:|
"""
    for rec in tumor_rows:
        md += f"| {rec['panel']} | {fmt(rec['auc'])} | {fmt(rec['mw_p'])} |\n"

    md += f"""
Paired tumor-minus-normal medians across the 4 complete pairs:

| Panel | n pairs | median tumor-normal | Wilcoxon p |
|---|---:|---:|---:|
"""
    for rec in paired_summary:
        md += (
            f"| {rec['panel']} | {rec['n_pairs']} | "
            f"{fmt(rec['median_tumor_minus_normal'])} | {fmt(rec['wilcoxon_p'])} |\n"
        )

    ht_mapk_rho = corr.loc["HT_13", "MAPK_9"]
    fa_mapk_rho = corr.loc["FA_12", "MAPK_9"]
    md += f"""
Spearman correlations across the 9 samples are descriptive only:
HT-13 vs MAPK-9 rho={fmt(ht_mapk_rho)}; FA-12 vs MAPK-9 rho={fmt(fa_mapk_rho)}.
The signs are compatible with the broader two-axis story, but n=9 and all-DM2
labels make this non-deployable as evidence.

## Disposition

- Do **not** call this "K2 full cohort validation".
- Do **not** cite it as Korean DM1-vs-DM2 replication.
- Safe use: one reviewer-reserve sentence saying that the partial local K2
  re-quant was technically successful but underpowered/all-DM2 and therefore
  did not alter the Paper 1/2 decision.

## Files

- `{SCORES.name}`
- `{AUCS.name}`
- `{paired_path.name}`
- `{paired_summary_path.name}`
- `{corr_path.name}`
- `{fig_path.name}`
"""
    report_path = OUT / "K2_PILOT_VERDICT.md"
    report_path.write_text(md)

    print(f"WROTE {report_path}")
    print(f"WROTE {fig_path}")
    print(f"WROTE {paired_summary_path}")
    print(f"WROTE {corr_path}")


if __name__ == "__main__":
    main()
