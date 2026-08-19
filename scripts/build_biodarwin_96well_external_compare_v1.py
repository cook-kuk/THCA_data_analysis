#!/usr/bin/env python3
"""Combine the BioDarwin 96-well interpreter with external paper baselines.

This report is intentionally split into two comparable layers:
1. local 96-well wetlab plate algorithms
2. external benchmark winners from prior paper baselines

External methods are benchmark winners on their own datasets, not row-wise
predictions on the 96-well plate.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
LOCAL = ROOT / "project/results/biodarwin_96well_interpreter_2026_05_11"
EXT = ROOT / "project/results/biodarwin_external_mega_comparison_2026_05_11"
OUT = ROOT / "project/results/biodarwin_96well_external_compare_2026_05_11"
OUT.mkdir(parents=True, exist_ok=True)
HUB = ROOT / "project/papers_hub_2026_05_04"
LIVE_HUB = Path("/var/www/papers/papers_hub_2026_05_04")


def fmt_float(x: object, nd: int = 4) -> str:
    if pd.isna(x):
        return ""
    return f"{float(x):.{nd}f}"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def build_local_table(df: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "algorithm",
        "AUPRC",
        "AUROC",
        "top5_hits",
        "top5_precision",
        "top10_hits",
        "top10_precision",
        "top24_hits",
        "top24_precision",
    ]
    out = df[keep].copy()
    out.insert(1, "source", "96-well local plate")
    out.insert(2, "family", "BioDarwin / local comparators")
    out.insert(3, "dataset", "TG4050 96-well locked-control")
    out.insert(4, "claim_status", "wetlab routing comparator")
    return out


def build_external_table(df: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "dataset",
        "algorithm",
        "AUPRC",
        "AUROC",
        "top5_hits",
        "top5_precision",
        "top10_hits",
        "top10_precision",
        "top24_hits",
        "top24_precision",
        "training_overlap_fraction",
        "overlap_warning",
        "benchmark_family",
        "claim_status",
    ]
    out = df[cols].copy()
    out.insert(1, "source", "external paper baseline")
    return out


def to_html(df: pd.DataFrame) -> str:
    show = df.copy()
    for col in show.columns:
        if pd.api.types.is_float_dtype(show[col]):
            show[col] = show[col].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
    return show.to_html(index=False, escape=True, classes="data")


def write_html(summary: dict[str, object], local_table: pd.DataFrame, ext_table: pd.DataFrame, combined: pd.DataFrame) -> None:
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>BioDarwin 96-well + external baseline comparison v1</title>
<style>
body {{ margin:0; background:#0f1418; color:#e9edf0; font:14px/1.55 system-ui,-apple-system,Segoe UI,sans-serif; }}
header {{ padding:34px 5vw 20px; background:#0c1014; border-bottom:1px solid #27313a; }}
h1 {{ margin:.2rem 0 .45rem; font-size:34px; line-height:1.05; }}
.lead {{ color:#cfd8dd; max-width:1200px; }}
main {{ padding:24px 5vw 56px; }}
.note {{ background:#171e24; border-left:3px solid #c59b3b; padding:10px 14px; margin:0 0 18px; }}
table.data {{ width:100%; border-collapse:collapse; font-size:13px; margin:10px 0 28px; }}
table.data th, table.data td {{ border-bottom:1px solid #2a3640; padding:7px 8px; text-align:left; vertical-align:top; }}
table.data th {{ color:#f1d58b; background:#141b21; position:sticky; top:0; }}
</style>
</head>
<body>
<header>
  <div style="color:#c59b3b;text-transform:uppercase;letter-spacing:.12em;font-weight:700;font-size:12px;">BioDarwin 96-well + external baselines</div>
  <h1>Local wetlab plate and published baselines in one view</h1>
  <p class="lead">The 96-well plate uses row-level local scoring. The external paper baselines are benchmark winners on their own datasets, included here for comparison and reviewer framing.</p>
</header>
<main>
  <div class="note">External methods are not scored row-by-row on the 96-well plate here. They are pulled from already computed benchmark summaries and should be read as dataset winners, not direct wetlab plate predictions.</div>
  <h2>Summary</h2>
  <pre>{json.dumps(summary, indent=2)}</pre>
  <h2>96-well local algorithms</h2>
  {to_html(local_table)}
  <h2>External paper baselines</h2>
  {to_html(ext_table)}
  <h2>Combined winner map</h2>
  {to_html(combined)}
</main>
</body>
</html>
"""
    path = HUB / "biodarwin_96well_external_compare_v1.html"
    path.write_text(html, encoding="utf-8")
    try:
        shutil.copy2(path, LIVE_HUB / path.name)
    except Exception as exc:
        print(f"deploy warning: {exc}")


def main() -> None:
    local = read_tsv(LOCAL / "biodarwin_96well_interpreter_v1_metrics.tsv")
    ext = read_tsv(EXT / "biodarwin_reviewer_safe_best_per_dataset.tsv")

    local = local.copy()
    local["claim_status"] = "wetlab routing comparator"
    local["overlap_warning"] = False
    local["top5_precision"] = pd.to_numeric(local["top5_hits"], errors="coerce") / 5.0
    local["top10_precision"] = pd.to_numeric(local["top10_hits"], errors="coerce") / 10.0
    local["top24_precision"] = pd.to_numeric(local["top24_hits"], errors="coerce") / 24.0
    ext = ext.copy()
    for col in ["top5_hits", "top10_hits", "top24_hits"]:
        if col not in ext.columns:
            ext[col] = np.nan

    local_table = build_local_table(local)
    ext_table = build_external_table(ext)

    combined = pd.concat(
        [
            local.assign(layer="local 96-well"),
            ext.rename(columns={"dataset": "dataset", "algorithm": "algorithm"}).assign(layer="external benchmark"),
        ],
        ignore_index=True,
        sort=False,
    )
    combined = combined[
        [
            "layer",
            "dataset",
            "algorithm",
            "AUPRC",
            "AUROC",
            "top5_hits",
            "top10_hits",
            "top24_hits",
            "top5_precision",
            "top10_precision",
            "top24_precision",
            "claim_status",
            "overlap_warning",
        ]
    ].copy()
    combined["AUPRC"] = combined["AUPRC"].map(fmt_float)
    combined["AUROC"] = combined["AUROC"].map(fmt_float)
    for c in ["top5_hits", "top10_hits", "top24_hits"]:
        if c in combined.columns:
            combined[c] = combined[c].fillna("").astype(str)

    summary = {
        "local_wells": int(96),
        "local_locked_metric_rows": int(25),
        "external_datasets": int(ext["dataset"].nunique()),
        "external_rows": int(len(ext)),
        "local_winner": local.sort_values("AUPRC", ascending=False).iloc[0].to_dict(),
        "external_winners": ext.sort_values(["AUPRC", "AUROC"], ascending=False).head(5).to_dict(orient="records"),
        "outputs": {
            "html": str(HUB / "biodarwin_96well_external_compare_v1.html"),
            "local_metrics": str(LOCAL / "biodarwin_96well_interpreter_v1_metrics.tsv"),
            "external_review_safe": str(EXT / "biodarwin_reviewer_safe_best_per_dataset.tsv"),
        },
    }
    (OUT / "biodarwin_96well_external_compare_v1_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    local_table.to_csv(OUT / "biodarwin_96well_local_algorithms.tsv", sep="\t", index=False)
    ext_table.to_csv(OUT / "biodarwin_external_baselines.tsv", sep="\t", index=False)
    combined.to_csv(OUT / "biodarwin_96well_external_combined.tsv", sep="\t", index=False)
    (OUT / "BIODARWIN_96WELL_EXTERNAL_COMPARE_V1.md").write_text(
        "\n".join(
            [
                "# BioDarwin 96-well + external baseline comparison v1",
                "",
                pd.DataFrame([summary]).to_markdown(index=False),
                "",
                "## Local 96-well algorithms",
                local_table.to_markdown(index=False),
                "",
                "## External paper baselines",
                ext_table.to_markdown(index=False),
                "",
                "## Combined winner map",
                combined.to_markdown(index=False),
            ]
        ),
        encoding="utf-8",
    )
    write_html(summary, local_table, ext_table, combined)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
