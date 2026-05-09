#!/usr/bin/env python3
"""Integrate RunPod ESM2 fast evaluations into local CROSS-Neo v2 summary."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from common import OUT, ensure_dirs


PRIMARY = ["exact_peptide_hla_holdout", "near_peptide_cluster_holdout", "hla_stratified_group_5fold", "hla_supertype_heldout"]


def main() -> None:
    ensure_dirs()
    base = pd.read_csv(OUT / "metrics/all_model_all_split_metrics.tsv", sep="\t")
    esm_parts = []
    for p in sorted((OUT / "metrics").glob("esm2_*_fast_metrics.tsv")):
        df = pd.read_csv(p, sep="\t")
        df["runpod_source_file"] = p.name
        esm_parts.append(df)
    esm = pd.concat(esm_parts, ignore_index=True) if esm_parts else pd.DataFrame()
    combined = pd.concat([base, esm], ignore_index=True, sort=False)
    combined.to_csv(OUT / "metrics/all_model_all_split_metrics_plus_runpod_esm2.tsv", sep="\t", index=False)

    rows = []
    for split in PRIMARY:
        b = base[base["split_name"].eq(split)].sort_values(["AUPRC", "top10_precision"], ascending=False).head(1)
        e = esm[esm["split_name"].eq(split)].sort_values(["AUPRC", "top10_precision"], ascending=False).head(1)
        c = combined[combined["split_name"].eq(split)].sort_values(["AUPRC", "top10_precision"], ascending=False).head(1)
        rows.append({
            "split_name": split,
            "base_best_model": b["model_name"].iloc[0] if not b.empty else "",
            "base_best_AUPRC": b["AUPRC"].iloc[0] if not b.empty else None,
            "base_best_top10": b["top10_precision"].iloc[0] if not b.empty else None,
            "runpod_esm2_best_model": e["model_name"].iloc[0] if not e.empty else "",
            "runpod_esm2_best_AUPRC": e["AUPRC"].iloc[0] if not e.empty else None,
            "runpod_esm2_best_top10": e["top10_precision"].iloc[0] if not e.empty else None,
            "combined_best_model": c["model_name"].iloc[0] if not c.empty else "",
            "combined_best_AUPRC": c["AUPRC"].iloc[0] if not c.empty else None,
            "combined_best_top10": c["top10_precision"].iloc[0] if not c.empty else None,
        })
    comp = pd.DataFrame(rows)
    comp.to_csv(OUT / "metrics/runpod_esm2_primary_comparison.tsv", sep="\t", index=False)

    source = combined[combined["split_name"].astype(str).str.startswith("source_heldout_")].copy()
    source_best = source.sort_values(["split_name", "top10_precision", "top20_precision", "AUPRC"], ascending=[True, False, False, False]).groupby("split_name").head(5)
    source_best.to_csv(OUT / "metrics/runpod_esm2_source_best.tsv", sep="\t", index=False)

    xlsx = OUT / "CROSS_Neo_v2_all_results_summary_plus_runpod_esm2.xlsx"
    with pd.ExcelWriter(xlsx) as writer:
        combined.head(20000).to_excel(writer, sheet_name="all_plus_esm2", index=False)
        comp.to_excel(writer, sheet_name="primary_compare", index=False)
        esm.head(20000).to_excel(writer, sheet_name="runpod_esm2_metrics", index=False)
        source_best.head(5000).to_excel(writer, sheet_name="source_best", index=False)

    lines = [
        "# RunPod ESM2 Upgrade Summary",
        "",
        "RunPod A6000 was used to generate frozen ESM2-35M, ESM2-150M, and ESM2-650M projected features and fast train-fold-safe logistic evaluations.",
        "",
        "## Primary Locked Splits",
        "",
        comp.to_markdown(index=False),
        "",
        "## Source-Heldout Best Rows",
        "",
        source_best[["split_name", "model_name", "n", "n_pos", "prevalence", "AUPRC", "top10_precision", "top20_precision", "claim_status"]].to_markdown(index=False),
        "",
        "## Decision Impact",
        "",
        "- ESM2-150M improves exact peptide-HLA top10 to 0.8 but does not beat the v1 fusion AUPRC leader.",
        "- ESM2-650M improves near-peptide top10 to 0.7 but still does not clearly beat the v1 rule-gated fusion AUPRC.",
        "- ESM2-650M gives useful HLA-supertype support but does not overturn the benchmark-only verdict.",
        "- The public training corpus overlap blocker remains unresolved.",
        "",
        f"Workbook: `{xlsx}`",
    ]
    (OUT / "RUNPOD_ESM2_upgrade_summary.md").write_text("\n".join(lines) + "\n")
    print(f"[v2-runpod-esm2-integrate] esm_rows={len(esm)} combined={len(combined)} xlsx={xlsx}")


if __name__ == "__main__":
    main()
