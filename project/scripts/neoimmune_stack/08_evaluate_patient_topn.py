#!/usr/bin/env python3
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from neoimmune_common import ensure_run_dir, patient_topn_metrics, safe_read_table, write_md, write_tsv


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    prod = safe_read_table(outdir / "production_track" / "production_stack_predictions.tsv")
    clean = safe_read_table(outdir / "clean_track" / "clean_science_predictions.tsv")
    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")

    d = canon.merge(prod[["candidate_id", "production_stack_score", "model_disagreement_score", "abstention_coverage"]], on="candidate_id", how="left")
    d = d.merge(clean[["candidate_id", "clean_science_score"]], on="candidate_id", how="left")
    d["patient_id"] = d["patient_id"].fillna("unknown_patient").astype(str)
    d["production_stack_score"] = pd.to_numeric(d["production_stack_score"], errors="coerce").fillna(0.5)
    d["clean_science_score"] = pd.to_numeric(d["clean_science_score"], errors="coerce").fillna(0.5)
    d["label_immunogenicity"] = pd.to_numeric(d["label_immunogenicity"], errors="coerce")

    top = d.sort_values(["patient_id", "production_stack_score"], ascending=[True, False]).groupby("patient_id", dropna=False).head(20)
    top_cols = [
        "patient_id",
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
        "abstention_coverage",
        "expression_tpm",
        "vaf",
        "hla_loh_status",
        "b2m_status",
        "source_study",
    ]
    for c in top_cols:
        if c not in top.columns:
            top[c] = np.nan
    write_tsv(top[top_cols], outdir / "predictions" / "patient_top20_candidates.tsv")

    rows = []
    eval_df = d[d["label_immunogenicity"].notna()].copy()
    for score_col in ["production_stack_score", "clean_science_score"]:
        for n in [10, 20, 50]:
            m = patient_topn_metrics(eval_df, score_col, "label_immunogenicity", n)
            m.update({"track": "production_stack_track" if score_col == "production_stack_score" else "clean_science_track", "score_col": score_col, "N": n})
            rows.append(m)
    met = pd.DataFrame(rows)
    write_tsv(met, outdir / "metrics" / "patient_topN_metrics.tsv")
    real_patient = d["patient_id"].ne("unknown_patient")
    patient_warning = ""
    if d.loc[real_patient, "patient_id"].nunique() < 3:
        patient_warning = (
            "WARNING: public/integrated tables contain fewer than 3 real patient IDs. "
            "Patient-level top-N files are structurally generated but are not a valid patient-level endpoint yet."
        )
    md = [
        "# Patient top-N report",
        "",
        f"- Patients represented: {d['patient_id'].nunique():,}",
        f"- Real non-placeholder patients represented: {d.loc[real_patient, 'patient_id'].nunique():,}",
        f"- Top-20 candidate rows written: {len(top):,}",
        f"- Interpretation warning: {patient_warning or 'patient-level endpoint has enough patient IDs for basic reporting'}",
        "",
        "## Primary endpoint",
        "Use patient-level Recall@20 and patient-level hit rate as the primary ranking endpoints when patient labels are available.",
        "",
        "## Boundary",
        "Top-20 candidates are prioritization hypotheses for discussion, not validated vaccine candidates.",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "patient_topN_report.md")


if __name__ == "__main__":
    main()
