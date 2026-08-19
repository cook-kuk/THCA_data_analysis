#!/usr/bin/env python3
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from neoimmune_common import ensure_run_dir, patient_topn_metrics, safe_read_table, write_md, write_tsv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    prod = safe_read_table(outdir / "production_track" / "production_stack_predictions.tsv")
    clean = safe_read_table(outdir / "clean_track" / "clean_science_predictions.tsv")
    d = canon.merge(prod[["candidate_id", "production_stack_score", "model_disagreement_score", "abstention_coverage"]], on="candidate_id", how="left")
    d = d.merge(clean[["candidate_id", "clean_science_score"]], on="candidate_id", how="left")
    d["patient_id"] = d["patient_id"].fillna("unknown_patient").astype(str)
    d["production_stack_score"] = pd.to_numeric(d["production_stack_score"], errors="coerce").fillna(0.5)
    d["clean_science_score"] = pd.to_numeric(d["clean_science_score"], errors="coerce").fillna(0.5)
    d["label_immunogenicity"] = pd.to_numeric(d["label_immunogenicity"], errors="coerce")
    d["moderna_like_top34_rank"] = d.groupby("patient_id")["production_stack_score"].rank(method="first", ascending=False)
    top34 = d[d["moderna_like_top34_rank"] <= 34].copy().sort_values(["patient_id", "moderna_like_top34_rank"])
    cols = [
        "patient_id",
        "moderna_like_top34_rank",
        "candidate_id",
        "dataset_source",
        "tumor_type",
        "gene",
        "mutation_id",
        "peptide_mut",
        "peptide_wt",
        "hla_allele",
        "label_immunogenicity",
        "production_stack_score",
        "clean_science_score",
        "model_disagreement_score",
        "expression_tpm",
        "vaf",
        "hla_loh_status",
        "b2m_status",
        "source_study",
    ]
    for c in cols:
        if c not in top34.columns:
            top34[c] = np.nan
    write_tsv(top34[cols], outdir / "predictions" / "patient_top34_candidates_moderna_like.tsv")

    rows = []
    eval_df = d[d["label_immunogenicity"].notna()].copy()
    for score_col in ["production_stack_score", "clean_science_score"]:
        m = patient_topn_metrics(eval_df, score_col, "label_immunogenicity", 34)
        m.update(
            {
                "track": "production_stack_track" if score_col == "production_stack_score" else "clean_science_track",
                "score_col": score_col,
                "N": 34,
                "rationale": "V940-public-analog cap; Moderna reports up to 34 neoantigens encoded per patient.",
            }
        )
        rows.append(m)
    met = pd.DataFrame(rows)
    write_tsv(met, outdir / "metrics" / "patient_top34_moderna_like_metrics.tsv")
    md = [
        "# Moderna-like top34 report",
        "",
        "V940/mRNA-4157 is publicly described as encoding up to 34 patient-specific neoantigens. This report adds a top34 cap to NeoImmune-Stack outputs.",
        "",
        "## Boundary",
        "- This is a ranking cap inspired by public V940 design logic.",
        "- It is not an mRNA construct, manufacturing protocol, or proprietary Moderna algorithm.",
        "- It is not a clinical efficacy claim.",
        "",
        "## Outputs",
        "- `predictions/patient_top34_candidates_moderna_like.tsv`",
        "- `metrics/patient_top34_moderna_like_metrics.tsv`",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "moderna_like_top34_report.md")
    print(outdir / "predictions" / "patient_top34_candidates_moderna_like.tsv")


if __name__ == "__main__":
    main()
