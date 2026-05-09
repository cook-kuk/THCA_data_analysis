#!/usr/bin/env python3
"""Write manuscript-track scaffolds for CROSS-Neo 2.0."""

from __future__ import annotations

import pandas as pd

from common import OUT, ensure_dirs


def main() -> None:
    ensure_dirs()
    out = OUT / "output"
    out.mkdir(exist_ok=True)
    metrics = pd.read_csv(OUT / "metrics/headline_locked_metrics.tsv", sep="\t") if (OUT / "metrics/headline_locked_metrics.tsv").exists() else pd.DataFrame()
    source = pd.read_csv(OUT / "metrics/source_heldout_metrics.tsv", sep="\t") if (OUT / "metrics/source_heldout_metrics.tsv").exists() else pd.DataFrame()
    headline_table = metrics.head(12).to_markdown(index=False) if not metrics.empty else "No headline table generated."
    source_table = source.sort_values(["split_name", "top20_precision"], ascending=[True, False]).groupby("split_name").head(3).to_markdown(index=False) if not source.empty else "No source table generated."

    ai = [
        "# CROSS-Neo 2.0 AI Track Scaffold",
        "",
        "Target venues: NeurIPS Datasets & Benchmarks, ICLR/ICML biology workshop, MLCB.",
        "",
        "Thesis: neoantigen predictor comparisons can be inflated by hidden overlap and source shift; CROSS-Neo 2.0 provides a contamination-controlled registry, strict splits, and a source-robust counterfactual ranking pipeline.",
        "",
        "## Locked Internal Results",
        "",
        headline_table,
        "",
        "## Source-Shift Stress Test",
        "",
        source_table,
        "",
        "## Required Claim Boundary",
        "",
        "Use internal locked-split / source-heldout language only. Do not claim external validation, clinical readiness, or quantum advantage.",
    ]
    (out / "manuscript_ai_track.md").write_text("\n".join(ai) + "\n")

    bio = [
        "# CROSS-Neo 2.0 Bio Journal Scaffold",
        "",
        "Target venues: Nature Machine Intelligence, Nature Communications, Genome Medicine, Cell Systems, Briefings in Bioinformatics.",
        "",
        "Thesis: an auditable neoantigen prioritization framework improves strict internal ranking while explicitly exposing source-shift failure modes.",
        "",
        "## Biological Interpretability Hooks",
        "",
        "- Mutant-WT counterfactual features separate self-like changes from novel peptide signals.",
        "- HLA/supertype and source-shift analyses identify when top-k prioritization is unreliable.",
        "- Case audit tables support candidate review without making clinical claims.",
    ]
    (out / "manuscript_bio_track.md").write_text("\n".join(bio) + "\n")

    wetlab = [
        "# Wet-Lab Validation Plan",
        "",
        "Pre-register before assay execution.",
        "",
        "1. Select 20 model-rescued candidates and 20 anchor/baseline candidates under the locked selection protocol.",
        "2. Match HLA alleles and source constraints; exclude unresolved public-overlap contaminated candidates when possible.",
        "3. Run IFN-gamma ELISPOT, 4-1BB activation, tetramer/multimer staining, or MHC stabilization depending on available material.",
        "4. Compare CROSS-Neo 2.0 against comparator outputs only after overlap audit labels are assigned.",
        "5. Freeze candidate list before wet-lab labels are observed.",
    ]
    (out / "wetlab_validation_plan.md").write_text("\n".join(wetlab) + "\n")
    print("[v2-manuscript] wrote AI, bio, and wetlab scaffold tracks")


if __name__ == "__main__":
    main()
