#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from neoimmune_common import ensure_run_dir, safe_read_table, write_md, write_tsv


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)

    rescue = safe_read_table(outdir / "predictions" / "source_rescue_queue.tsv")
    traps = safe_read_table(outdir / "predictions" / "binding_only_false_positive_traps.tsv")
    gaps = safe_read_table(outdir / "metrics" / "source_transfer_gap_summary.tsv")

    rescue = rescue.copy()
    rescue["dataset_source"] = rescue["dataset_source"].fillna("unknown")
    rescue["label_immunogenicity"] = pd.to_numeric(rescue["label_immunogenicity"], errors="coerce").fillna(0).astype(int)
    rescue["rescue_priority"] = pd.to_numeric(rescue["rescue_priority"], errors="coerce")
    rescue = rescue.sort_values(["label_immunogenicity", "rescue_priority", "clean_science_score"], ascending=[False, False, False])

    positive_handoff = (
        rescue[rescue["label_immunogenicity"].eq(1)]
        .groupby("dataset_source", group_keys=False)
        .head(3)
        .copy()
    )
    negative_controls = (
        rescue[rescue["label_immunogenicity"].eq(0)]
        .groupby("dataset_source", group_keys=False)
        .head(2)
        .copy()
    )

    def pick_cols(df: pd.DataFrame) -> pd.DataFrame:
        cols = [
            "candidate_id",
            "dataset_source",
            "source_study",
            "patient_id",
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
            if c not in df.columns:
                df[c] = np.nan
        return df[cols]

    positive_out = pick_cols(positive_handoff)
    negative_out = pick_cols(negative_controls)

    write_tsv(positive_out, outdir / "predictions" / "collaborator_positive_handoff.tsv")
    write_tsv(negative_out, outdir / "predictions" / "collaborator_control_handoff.tsv")

    summary = (
        positive_out.groupby("dataset_source")
        .agg(
            n=("candidate_id", "size"),
            top_priority=("rescue_priority", "max"),
            mean_priority=("rescue_priority", "mean"),
            mean_clean=("clean_science_score", "mean"),
            mean_delta=("local_rescue_delta_vs_public", "mean"),
        )
        .reset_index()
        .sort_values("top_priority", ascending=False)
    )
    write_tsv(summary, outdir / "metrics" / "collaborator_handoff_summary.tsv")

    try:
        import matplotlib.pyplot as plt

        figdir = outdir / "figures"
        figdir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(12, 6), facecolor="#080b12")
        ax.set_facecolor("#111827")
        d = summary.sort_values("top_priority", ascending=True)
        ax.barh(d["dataset_source"], d["top_priority"], color="#26d9ff")
        ax.set_xlabel("top rescue priority")
        ax.set_title("Collaborator handoff: positive source-rescue shortlist", color="white", fontsize=14)
        ax.tick_params(colors="white")
        ax.xaxis.label.set_color("white")
        for spine in ax.spines.values():
            spine.set_color("#263348")
        ax.grid(axis="x", color="#263348", alpha=0.5)
        fig.tight_layout()
        fig.savefig(figdir / "collaborator_positive_handoff.png", dpi=220)
        plt.close(fig)
    except Exception:
        pass

    md = [
        "# Collaborator Handoff",
        "",
        "This is the minimal wet-lab-facing shortlist.",
        "",
        "## Positive handoff rule",
        "- Top 3 positive candidates per source family.",
        "- Prefer large clean score, positive local rescue delta, and high model disagreement.",
        "",
        positive_out.to_markdown(index=False) if not positive_out.empty else "_No positive candidates available._",
        "",
        "## Negative control handoff",
        "- Use the top 2 negative rows per source family as failure controls or calibration controls.",
        "",
        negative_out.to_markdown(index=False) if not negative_out.empty else "_No negative control candidates available._",
        "",
        "## Interpretation",
        "- CEDAR and dbPepNeo2_MHCI are the most immediately useful collaborator handoff families.",
        "- TESLA rows are present but the current rescue signal is weaker; they should be handled as secondary validation, not lead claims.",
        "",
        "## Claim boundary",
        "- Retrospective handoff only.",
        "- No candidate is validated without orthogonal assay support.",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "COLLABORATOR_HANDOFF.md")
    print(outdir / "reports" / "COLLABORATOR_HANDOFF.md")


if __name__ == "__main__":
    main()
