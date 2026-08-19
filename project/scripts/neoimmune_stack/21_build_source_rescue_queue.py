#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from neoimmune_common import ensure_run_dir, safe_read_table, write_md, write_tsv


TARGET_SOURCES = [
    "TESLA",
    "TESLA_mmc4",
    "TESLA_mmc7_validation",
    "dbPepNeo2_MHCI",
    "dbPepNeo2_MHCII",
    "CEDAR",
    "IMPROVE_Neoepitopes_CEDAR_benchmark_data",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)

    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    high = safe_read_table(outdir / "predictions" / "high_disagreement_candidates.tsv")
    rescue = safe_read_table(outdir / "predictions" / "local_rescued_public_missed_positives.tsv")
    traps = safe_read_table(outdir / "predictions" / "binding_only_false_positive_traps.tsv")

    # Normalize source names to the canonical table.
    canon = canon[[
        "candidate_id",
        "dataset_source",
        "source_study",
        "patient_id",
        "tumor_type",
        "gene",
        "peptide_mut",
        "peptide_wt",
        "hla_allele",
        "label_immunogenicity",
    ]].copy()

    # Priority queue starts from high-disagreement rows and extends with rescued positives.
    d = pd.concat([high, rescue], ignore_index=True, sort=False)
    if "candidate_id" in d.columns:
        d = d.drop_duplicates(subset=["candidate_id"], keep="first")
    d = d.merge(canon, on="candidate_id", how="left", suffixes=("", "_canon"))
    d["dataset_source"] = d["dataset_source"].fillna(d["dataset_source_canon"])
    d["label_immunogenicity"] = pd.to_numeric(d.get("label_immunogenicity"), errors="coerce")
    d["clean_science_score"] = pd.to_numeric(d.get("clean_science_score"), errors="coerce")
    d["production_stack_score"] = pd.to_numeric(d.get("production_stack_score"), errors="coerce")
    d["public_best"] = pd.to_numeric(d.get("public_best"), errors="coerce")
    d["model_disagreement_score"] = pd.to_numeric(d.get("model_disagreement_score"), errors="coerce")
    d["local_rescue_delta_vs_public"] = pd.to_numeric(d.get("local_rescue_delta_vs_public"), errors="coerce")

    target_mask = d["dataset_source"].fillna("").astype(str).isin(TARGET_SOURCES)
    q = d[target_mask].copy()
    if q.empty:
        raise RuntimeError("No rescue candidates matched the target sources.")

    q["rescue_priority"] = (
        0.45 * q["clean_science_score"].fillna(0)
        + 0.20 * q["production_stack_score"].fillna(0)
        + 0.15 * q["model_disagreement_score"].fillna(0)
        + 0.20 * q["local_rescue_delta_vs_public"].fillna(0).clip(lower=0)
    )
    q = q.sort_values(["label_immunogenicity", "rescue_priority", "clean_science_score"], ascending=[False, False, False])

    summary = (
        q.groupby("dataset_source", dropna=False)
        .agg(
            rows=("candidate_id", "size"),
            positives=("label_immunogenicity", "sum"),
            top_priority=("rescue_priority", "max"),
            median_priority=("rescue_priority", "median"),
            mean_clean=("clean_science_score", "mean"),
            mean_public_best=("public_best", "mean"),
            mean_delta=("local_rescue_delta_vs_public", "mean"),
        )
        .reset_index()
        .sort_values("top_priority", ascending=False)
    )

    cols = [
        "candidate_id",
        "dataset_source",
        "source_study",
        "patient_id",
        "tumor_type",
        "gene",
        "peptide_mut",
        "peptide_wt",
        "hla_allele",
        "label_immunogenicity",
        "public_best",
        "clean_science_score",
        "production_stack_score",
        "local_rescue_delta_vs_public",
        "model_disagreement_score",
        "rescue_priority",
    ]
    for c in cols:
        if c not in q.columns:
            q[c] = np.nan
    write_tsv(q[cols], outdir / "predictions" / "source_rescue_queue.tsv")
    write_tsv(summary, outdir / "metrics" / "source_rescue_queue_summary.tsv")

    # Simple figure: top 20 rescue priority items.
    try:
        import matplotlib.pyplot as plt

        top = q.head(20).copy()
        fig, ax = plt.subplots(figsize=(12, 7), facecolor="#080b12")
        ax.set_facecolor("#111827")
        ax.barh(top["candidate_id"], top["rescue_priority"], color=["#26d9ff" if x > 0 else "#fb7185" for x in top["label_immunogenicity"].fillna(0)])
        ax.invert_yaxis()
        ax.set_xlabel("rescue priority")
        ax.set_title("Source rescue queue: high-disagreement candidates in TESLA/dbPepNeo/CEDAR/IMPROVE", color="white", fontsize=14)
        ax.tick_params(colors="white", labelsize=8)
        ax.xaxis.label.set_color("white")
        for spine in ax.spines.values():
            spine.set_color("#263348")
        ax.grid(axis="x", color="#263348", alpha=0.5)
        fig.tight_layout()
        fig.savefig(outdir / "figures" / "source_rescue_queue.png", dpi=220)
        plt.close(fig)
    except Exception:
        pass

    md = [
        "# Source Rescue Queue",
        "",
        "This queue is the practical handoff for source families where public reconstruction is still thin.",
        "",
        "## Reading rule",
        "- Positive labels with high rescue priority are the candidates most worth discussing with wet-lab collaborators.",
        "- Negative labels can still be useful as failure cases, but they are not rescue targets.",
        "",
        summary.to_markdown(index=False) if not summary.empty else "_No summary available._",
        "",
        "## Priority rule",
        "- High `clean_science_score` plus positive `local_rescue_delta_vs_public` means the candidate is rescued by local branch evidence beyond public predictors.",
        "- High `model_disagreement_score` means the candidate is a good audit case even when wet-lab validation is not immediate.",
        "",
        "## Claim boundary",
        "- This is a retrospective prioritization queue only.",
        "- No candidate is a validated vaccine target without orthogonal assay support.",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "SOURCE_RESCUE_QUEUE.md")
    print(outdir / "reports" / "SOURCE_RESCUE_QUEUE.md")


if __name__ == "__main__":
    main()
