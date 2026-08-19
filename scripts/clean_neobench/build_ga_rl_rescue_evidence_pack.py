#!/usr/bin/env python3
"""Build a rescue evidence pack for recoverable GA/RL watchlist candidates.

This pack focuses on the small set of rows that are not clean T1 survivors yet
but are not hard overlap blocks either. It turns them into a concrete evidence
intake queue by listing the exact metadata missing for translational uplift.
"""

from __future__ import annotations

import argparse
import html
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd

from common import ensure_dir, update_manifest, write_tsv


HTML_NAME = "ga_rl_rescue_evidence_results_2026_05_11.html"
OUTPUT_FILES = [
    "ga_rl_rescue_evidence.tsv",
    "ga_rl_rescue_evidence_summary.tsv",
    "GA_RL_RESCUE_EVIDENCE_KR.md",
]

REQUIRED_METADATA = [
    "wt_peptide",
    "gene",
    "mutation_id",
    "source_protein_window",
    "expression_tpm",
    "vaf",
    "clonality",
    "patient_id",
    "disease_context",
    "hla_loh",
    "b2m_status",
    "antigen_processing_status",
    "tumor_stage",
    "treatment_context",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def is_missing(value: Any) -> bool:
    if value is None:
        return True
    s = str(value).strip()
    return s == "" or s.lower() in {"na", "nan", "none", "null"}


def safe_md(df: pd.DataFrame, max_rows: int = 20) -> str:
    return "_No rows available._" if df.empty else df.head(max_rows).to_markdown(index=False)


def build_evidence(master: pd.DataFrame, rescue: pd.DataFrame) -> pd.DataFrame:
    if rescue.empty:
        return pd.DataFrame()
    enrich_cols = ["candidate_id", "source_name", "hla_allele_4digit", "peptide", "label", "ga_rl_barneo_priority_score", "top_nn_similarity"]
    enrich = rescue[[c for c in enrich_cols if c in rescue.columns]].copy()
    rows = []
    for _, row in enrich.iterrows():
        cid = str(row.get("candidate_id", ""))
        m = master[master["candidate_id"].astype(str).eq(cid)]
        if m.empty:
            continue
        mr = m.iloc[0]
        source_prev = master[master["source_name"].astype(str).eq(str(mr.get("source_name", "")))]["label"].mean()
        hla_prev = master[master["hla_allele_4digit"].astype(str).eq(str(mr.get("hla_allele_4digit", "")))]["label"].mean()
        missing = [f for f in REQUIRED_METADATA if is_missing(mr.get(f, None))]
        present = [f for f in REQUIRED_METADATA if f not in missing]
        evidence_grade = "A" if len(missing) <= 3 else ("B" if len(missing) <= 6 else "C")
        rows.append(
            {
                "candidate_id": cid,
                "peptide": mr.get("peptide", ""),
                "hla_allele_4digit": mr.get("hla_allele_4digit", ""),
                "source_name": mr.get("source_name", ""),
                "label": mr.get("label", pd.NA),
                "ga_rl_barneo_priority_score": row.get("ga_rl_barneo_priority_score", pd.NA),
                "top_nn_similarity": row.get("top_nn_similarity", pd.NA),
                "source_positive_prevalence": source_prev,
                "hla_positive_prevalence": hla_prev,
                "present_metadata_count": len(present),
                "missing_metadata_count": len(missing),
                "present_metadata": ";".join(present) if present else "none",
                "missing_metadata": ";".join(missing) if missing else "none",
                "evidence_grade": evidence_grade,
                "evidence_action": "intake missing metadata before T1 promotion",
                "evidence_boundary": "research triage only; not clinical selection",
                "rescue_note": "High-GA watchlist candidate without exact/near/public overlap.",
            }
        )
    return pd.DataFrame(rows).sort_values(["evidence_grade", "ga_rl_barneo_priority_score"], ascending=[True, False]).reset_index(drop=True)


def render_missingness_plot(df: pd.DataFrame, path: Path) -> None:
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(8.5, 3.8))
    ax.bar(df["candidate_id"], df["missing_metadata_count"], color="#f59e0b")
    ax.set_ylabel("Missing required fields")
    ax.set_title("GA/RL rescue evidence gaps")
    ax.set_ylim(0, max(1, int(df["missing_metadata_count"].max()) + 1))
    for i, r in enumerate(df.itertuples()):
        ax.text(i, float(r.missing_metadata_count) + 0.1, r.evidence_grade, ha="center", va="bottom", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    fig_dir = output_root / "figures"
    ensure_dir(fig_dir)

    master = read_tsv(output_root / "clean_neobench_master.tsv")
    rescue = read_tsv(output_root / "ga_rl_rescue_quadrant.tsv")
    if master.empty or rescue.empty:
        raise FileNotFoundError("Missing master or rescue quadrant table")

    rescue_watch = rescue[rescue["rescue_category"].astype(str).eq("recoverable_watchlist")].copy()
    evidence = build_evidence(master, rescue_watch)
    summary = (
        evidence.groupby("evidence_grade", dropna=False)
        .agg(
            n_candidates=("candidate_id", "count"),
            mean_missing=("missing_metadata_count", "mean"),
            mean_priority=("ga_rl_barneo_priority_score", "mean"),
            mean_source_prev=("source_positive_prevalence", "mean"),
            mean_hla_prev=("hla_positive_prevalence", "mean"),
        )
        .reset_index()
    )

    out_path = output_root / "ga_rl_rescue_evidence.tsv"
    sum_path = output_root / "ga_rl_rescue_evidence_summary.tsv"
    write_tsv(evidence, out_path)
    write_tsv(summary, sum_path)

    fig_path = fig_dir / "fig_ga_rl_rescue_evidence_missingness.png"
    render_missingness_plot(evidence, fig_path)

    md_path = output_root / "GA_RL_RESCUE_EVIDENCE_KR.md"
    md = f"""# GA/RL Rescue Evidence KR

## 한 줄 결론

회복 큐 3개는 exact/near/public overlap 없이 남았지만, 아직 patient/translational claim에는 필요한 metadata가 비어 있다. 이 패키지는 그 missingness를 바로 보여준다.

## Claim boundary

Allowed:

- rescue evidence intake
- metadata completion planning
- review-ready uplift queue

Forbidden:

- clinical vaccine selection
- new SOTA predictor
- public clean comparator claim without audit

## Summary

{safe_md(summary, 10)}

## Evidence queue

{safe_md(evidence, 20)}

## Missingness note

The rescue watchlist remains promising because it has no exact/near/public overlap, but the missing patient-level and biology fields keep it in research triage, not T1 translational claim.

## Output files

- `{out_path}`
- `{sum_path}`
- `{fig_path}`
"""
    md_path.write_text(md)

    html_path = output_root / HTML_NAME
    html_page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>GA/RL Rescue Evidence</title>
<style>
body{{margin:0;background:#0b1020;color:#e5e7eb;font-family:Inter,system-ui,sans-serif;line-height:1.5}}
main{{max-width:1100px;margin:0 auto;padding:28px 24px 48px}}
.card{{background:#111827;border:1px solid #243042;border-radius:8px;padding:16px 18px;margin:0 0 18px}}
h1,h2{{margin:0 0 12px}}
.muted{{color:#9ca3af}}
img{{max-width:100%;height:auto;border-radius:8px;border:1px solid #243042}}
table{{border-collapse:collapse;width:100%;margin-top:10px}}
th,td{{border:1px solid #243042;padding:6px 8px;font-size:13px;vertical-align:top}}
th{{background:#0f172a;position:sticky;top:0}}
</style></head><body><main>
<div class="card"><h1>GA/RL Rescue Evidence</h1><p class="muted">Missing metadata for the three recoverable watchlist candidates.</p></div>
<div class="card"><img src="assets/ga_rl_barneo/{fig_path.name}" alt="rescue evidence"></div>
<div class="card"><h2>Summary</h2>{summary.to_html(index=False, escape=False) if not summary.empty else '<p>No rows.</p>'}</div>
<div class="card"><h2>Evidence queue</h2>{evidence.to_html(index=False, escape=False) if not evidence.empty else '<p>No rows.</p>'}</div>
<div class="card"><h2>Source</h2><p class="muted">{html.escape(str(out_path))}</p></div>
</main></body></html>"""
    html_path.write_text(html_page)

    update_manifest(
        output_root,
        "ga_rl_rescue_evidence_pack",
        {
            "n_candidates": int(len(evidence)),
            "n_summary_rows": int(len(summary)),
            "output_files": OUTPUT_FILES + [HTML_NAME, f"figures/{fig_path.name}"],
            "warnings": [],
        },
    )

    print(f"[ga-rl-rescue-evidence] rows={len(evidence)}")


if __name__ == "__main__":
    main()
