#!/usr/bin/env python3
"""Build a separate benchmark pack for MHC-II scout rows.

The purpose is to keep the Class-I BAR-Neo / CLEAN-NeoBench story clean while
still preserving the Class-II scout output as a separate benchmark contract.
No pooled Class-I/Class-II claim is made here.
"""

from __future__ import annotations

import argparse
import html
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import ensure_dir, update_manifest, write_tsv


HTML_NAME = "mhc_ii_separate_benchmark_results_2026_05_11.html"
OUTPUT_FILES = [
    "mhc_ii_separate_benchmark.tsv",
    "mhc_ii_separate_benchmark_summary.tsv",
    "MHC_II_SEPARATE_BENCHMARK_KR.md",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    return p.parse_args()


def load_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def safe_md(df: pd.DataFrame, max_rows: int = 20) -> str:
    return "_No rows available._" if df.empty else df.head(max_rows).to_markdown(index=False)


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    fig_dir = output_root / "figures"
    ensure_dir(fig_dir)

    all_rows = load_tsv(output_root / "ga_rl_barneo_unique_candidate_summary.tsv")
    if all_rows.empty:
        raise FileNotFoundError(f"Missing GA/RL unique table: {output_root / 'ga_rl_barneo_unique_candidate_summary.tsv'}")

    ii = all_rows[all_rows.get("mhc_class", pd.Series(dtype=str)).astype(str).eq("II")].copy()
    if ii.empty:
        ii = all_rows[all_rows.get("ga_rl_track", pd.Series(dtype=str)).astype(str).str.contains("classII", case=False, na=False)].copy()
    ii = ii.sort_values(["ga_rl_score", "ga_rl_barneo_priority_score"], ascending=[False, False]).reset_index(drop=True)
    ii["separate_benchmark_rank"] = range(1, len(ii) + 1)

    summary_rows = []
    for cols in [["ga_rl_track"], ["source_name"], ["hla_allele_4digit"], ["ga_rl_barneo_decision"]]:
        if not all(c in ii.columns for c in cols):
            continue
        grouped = ii.groupby(cols, dropna=False)
        for keys, g in grouped:
            if not isinstance(keys, tuple):
                keys = (keys,)
            row = {c: v for c, v in zip(cols, keys)}
            row.update(
                {
                    "group_by": "+".join(cols),
                    "n_candidates": int(len(g)),
                    "mean_ga_rl_score": float(pd.to_numeric(g.get("ga_rl_score", pd.Series(dtype=float)), errors="coerce").mean()),
                    "mean_priority_score": float(pd.to_numeric(g.get("ga_rl_barneo_priority_score", pd.Series(dtype=float)), errors="coerce").mean()),
                }
            )
            summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)

    out_path = output_root / "mhc_ii_separate_benchmark.tsv"
    sum_path = output_root / "mhc_ii_separate_benchmark_summary.tsv"
    write_tsv(ii, out_path)
    write_tsv(summary, sum_path)

    fig_path = fig_dir / "fig_mhc_ii_separate_benchmark_scores.png"
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.hist(pd.to_numeric(ii.get("ga_rl_score", pd.Series(dtype=float)), errors="coerce").dropna(), bins=20, color="#60a5fa", edgecolor="black", linewidth=0.4)
    ax.set_title("MHC-II separate benchmark score distribution")
    ax.set_xlabel("GA/RL score")
    ax.set_ylabel("Candidate count")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    md_path = output_root / "MHC_II_SEPARATE_BENCHMARK_KR.md"
    md = f"""# MHC-II Separate Benchmark KR

## 한 줄 결론

Class-II scout rows are kept as a **separate benchmark contract**. They are not pooled with the Class-I BAR-Neo / CLEAN-NeoBench claim.

## Claim boundary

Allowed:

- MHC-II separate benchmark
- Class-II scout analysis
- separate contract reporting

Forbidden:

- pooled Class-I/Class-II predictor claim
- clinical vaccine selection
- external validation proven
- new SOTA claim

## Summary

{safe_md(summary, 20)}

## Top rows

{safe_md(ii[[c for c in ['separate_benchmark_rank','candidate_id','ga_rl_track','source_name','hla_allele_4digit','ga_rl_score','ga_rl_barneo_priority_score','ga_rl_barneo_decision','ga_rl_barneo_reason','locked_public','do_not_train','manual_QA_required'] if c in ii.columns]], 20)}

## Interpretation

- Rows: {len(ii)}
- Sources: {ii['source_name'].nunique() if 'source_name' in ii.columns else 0}
- HLA alleles: {ii['hla_allele_4digit'].nunique() if 'hla_allele_4digit' in ii.columns else 0}
- Decisions: {ii['ga_rl_barneo_decision'].nunique() if 'ga_rl_barneo_decision' in ii.columns else 0}

The class-II scout set is therefore a separate benchmark resource, not a Class-I comparator.

## Output files

- `{out_path}`
- `{sum_path}`
- `{fig_path}`
"""
    md_path.write_text(md)

    html_path = output_root / HTML_NAME
    html_page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>MHC-II Separate Benchmark</title>
<style>
body{{margin:0;background:#0b1020;color:#e5e7eb;font-family:Inter,system-ui,sans-serif;line-height:1.5}}
main{{max-width:1200px;margin:0 auto;padding:28px 24px 48px}}
.card{{background:#111827;border:1px solid #243042;border-radius:8px;padding:16px 18px;margin:0 0 18px}}
h1,h2{{margin:0 0 12px}}
.muted{{color:#9ca3af}}
img{{max-width:100%;height:auto;border-radius:8px;border:1px solid #243042}}
table{{border-collapse:collapse;width:100%;margin-top:10px}}
th,td{{border:1px solid #243042;padding:6px 8px;font-size:13px;vertical-align:top}}
th{{background:#0f172a;position:sticky;top:0}}
</style></head><body><main>
<div class="card"><h1>MHC-II Separate Benchmark</h1><p class="muted">Class-II scout rows are not pooled with the Class-I BAR-Neo benchmark.</p></div>
<div class="card"><img src="assets/ga_rl_barneo/{fig_path.name}" alt="MHC-II separate benchmark scores"></div>
<div class="card"><h2>Summary</h2>{summary.to_html(index=False, escape=False) if not summary.empty else '<p>No rows.</p>'}</div>
<div class="card"><h2>Rows</h2>{ii[[c for c in ['separate_benchmark_rank','candidate_id','ga_rl_track','source_name','hla_allele_4digit','ga_rl_score','ga_rl_barneo_priority_score','ga_rl_barneo_decision','ga_rl_barneo_reason'] if c in ii.columns]].to_html(index=False, escape=False) if not ii.empty else '<p>No rows.</p>'}</div>
<div class="card"><h2>Source</h2><p class="muted">{html.escape(str(out_path))}</p></div>
</main></body></html>"""
    html_path.write_text(html_page)

    update_manifest(
        output_root,
        "mhc_ii_separate_benchmark_pack",
        {
            "n_candidates": int(len(ii)),
            "n_sources": int(ii["source_name"].nunique()) if "source_name" in ii.columns else 0,
            "n_hla_alleles": int(ii["hla_allele_4digit"].nunique()) if "hla_allele_4digit" in ii.columns else 0,
            "output_files": OUTPUT_FILES + [HTML_NAME, f"figures/{fig_path.name}"],
            "warnings": [],
        },
    )

    print(f"[mhc-ii-separate-benchmark] rows={len(ii)}")


if __name__ == "__main__":
    main()
