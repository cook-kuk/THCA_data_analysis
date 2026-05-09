#!/usr/bin/env python3
"""Annotate peptide, MHC, beta2m, and TCR chains for MD complexes."""

from __future__ import annotations

import pandas as pd

from common_md import OUT, annotate_chains, chain_ids, discover_runs, ensure_dirs, write_markdown_table


def main() -> None:
    ensure_dirs()
    rows = []
    seq_rows = []
    for run in discover_runs():
        if not run.topology:
            continue
        ann = annotate_chains(run.topology, run.peptide)
        if ann.empty:
            continue
        peptide_ids = chain_ids(ann, "peptide")
        mhc_ids = chain_ids(ann, "mhc_heavy_chain")
        b2m_ids = chain_ids(ann, "beta2m")
        tcr_ids = chain_ids(ann, "tcr_candidate")
        for _, r in ann.iterrows():
            seq_rows.append(
                {
                    "run_id": run.run_id,
                    "candidate": run.candidate,
                    "topology": str(run.topology),
                    **r.to_dict(),
                }
            )
        rows.append(
            {
                "run_id": run.run_id,
                "candidate": run.candidate,
                "expected_peptide": run.peptide,
                "peptide_chains": "|".join(peptide_ids),
                "mhc_heavy_chain": "|".join(mhc_ids),
                "beta2m_chain": "|".join(b2m_ids),
                "tcr_candidate_chains": "|".join(tcr_ids),
                "tcr_available": bool(tcr_ids),
                "peptide_sequence_from_structure": "|".join(ann.loc[ann["role"].eq("peptide"), "sequence"].astype(str).tolist()),
                "peptide_match": bool(ann["peptide_sequence_match"].any()),
                "mutation_position": "unknown",
                "manual_review_flag": "" if peptide_ids and mhc_ids else "ambiguous_peptide_or_mhc_chain",
                "topology": str(run.topology),
            }
        )
    chain_df = pd.DataFrame(rows)
    seq_df = pd.DataFrame(seq_rows)
    chain_df.to_csv(OUT / "complex_chain_annotation.tsv", sep="\t", index=False)
    seq_df.to_csv(OUT / "complex_sequence_report.tsv", sep="\t", index=False)

    lines = [
        "# Chain Annotation Report",
        "",
        "Chains are assigned by sequence and length heuristics: exact expected short peptide sequence, longest HLA-like chain, beta2m-like 80-120 residue chain, and 150-260 residue TCR candidate chains.",
        "",
        write_markdown_table(
            chain_df,
            [
                "run_id",
                "candidate",
                "expected_peptide",
                "peptide_sequence_from_structure",
                "peptide_match",
                "peptide_chains",
                "mhc_heavy_chain",
                "beta2m_chain",
                "tcr_candidate_chains",
                "manual_review_flag",
            ],
            max_rows=100,
        ),
        "",
        "Mutation residue positions are marked unknown unless a validated mutant/WT mapping is supplied. Do not claim mutant-site contact specificity without that mapping.",
    ]
    (OUT / "chain_annotation_report.md").write_text("\n".join(lines) + "\n")
    print(f"[md-annotate] runs={len(chain_df)} out={OUT}")
    if not chain_df.empty:
        print(chain_df[["run_id", "peptide_chains", "mhc_heavy_chain", "beta2m_chain", "tcr_candidate_chains", "peptide_match"]].to_string(index=False))


if __name__ == "__main__":
    main()
