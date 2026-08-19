#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import defaultdict

import numpy as np
import pandas as pd

from neoimmune_common import ensure_run_dir, safe_read_table, write_md, write_tsv


def overlap_count(train: pd.Series, test: pd.Series) -> int:
    a = set(train.dropna().astype(str))
    b = set(test.dropna().astype(str))
    a.discard("")
    b.discard("")
    return len(a & b)


def leakage_flags(df: pd.DataFrame) -> pd.DataFrame:
    out = df[["candidate_id", "peptide_mut", "peptide_wt", "hla_allele", "dataset_source", "patient_id", "source_study", "train_test_group"]].copy()
    out["peptide_hla_pair"] = out["peptide_mut"].fillna("").astype(str) + "|" + out["hla_allele"].fillna("").astype(str)
    out["mut_wt_pair"] = out["peptide_mut"].fillna("").astype(str) + "|" + out["peptide_wt"].fillna("").astype(str)
    for c in ["exact_peptide_train_overlap", "exact_peptide_hla_train_overlap", "near_peptide_train_overlap", "source_protein_window_train_overlap", "study_train_overlap", "patient_train_overlap", "public_tool_training_overlap_any"]:
        if c in df.columns:
            out[c] = df[c]
    risk = pd.Series(False, index=out.index)
    for c in out.columns:
        if "overlap" in c:
            risk = risk | out[c].astype(str).str.lower().isin(["true", "1", "yes", "high", "medium"])
    if "leakage_risk_level" in df.columns:
        out["leakage_risk_level"] = df["leakage_risk_level"]
        risk = risk | df["leakage_risk_level"].astype(str).str.lower().str.contains("high")
    if "test_set_safety" in df.columns:
        out["test_set_safety"] = df["test_set_safety"]
        risk = risk | df["test_set_safety"].astype(str).str.contains("TRAINING_OVERLAP|PARTIAL_OVERLAP", case=False, na=False)
    out["any_existing_overlap_flag"] = risk
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    df = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    flags = leakage_flags(df)
    write_tsv(flags, outdir / "metrics" / "leakage_candidate_flags.tsv")

    rows = []
    rows.append({"leakage_control": "existing_overlap_flags", "n_flagged": int(flags["any_existing_overlap_flag"].sum()), "denominator": len(flags), "fraction": float(flags["any_existing_overlap_flag"].mean())})
    for split_col in ["train_test_group", "split", "split_source_heldout", "split_study_heldout", "split_patient_heldout", "split_hla_heldout", "split_exact_phla"]:
        if split_col not in df.columns:
            continue
        s = df[split_col].fillna("").astype(str)
        train = df[s.str.contains("train", case=False, na=False)]
        test = df[s.str.contains("test|heldout|valid|validation", case=False, na=False)]
        if train.empty or test.empty:
            continue
        controls = {
            "exact_peptide_leakage": ("peptide_mut",),
            "peptide_hla_pair_leakage": ("peptide_mut", "hla_allele"),
            "mutant_wt_pair_leakage": ("peptide_mut", "peptide_wt"),
            "study_leakage": ("source_study",),
            "patient_leakage": ("patient_id",),
            "hla_allele_leakage": ("hla_allele",),
        }
        for name, cols in controls.items():
            tr_key = train[list(cols)].fillna("").astype(str).agg("|".join, axis=1)
            te_key = test[list(cols)].fillna("").astype(str).agg("|".join, axis=1)
            rows.append({"leakage_control": f"{name}::{split_col}", "n_flagged": overlap_count(tr_key, te_key), "denominator": int(test.shape[0]), "fraction": np.nan})
    if "source_window_30aa" in df.columns:
        rows.append({"leakage_control": "source_protein_window_available", "n_flagged": int(df["source_window_30aa"].notna().sum()), "denominator": len(df), "fraction": float(df["source_window_30aa"].notna().mean())})

    summary = pd.DataFrame(rows)
    write_tsv(summary, outdir / "metrics" / "leakage_summary.tsv")

    md = [
        "# Leakage audit",
        "",
        "## Controls implemented",
        "- exact peptide leakage",
        "- peptide-HLA pair leakage",
        "- mutant-WT pair leakage",
        "- source protein/window leakage where source windows are available",
        "- study leakage",
        "- patient leakage",
        "- HLA allele leakage",
        "- external tool training contamination risk flags when provided by existing artifacts",
        "",
        "## Summary",
    ]
    for _, r in summary.iterrows():
        frac = "" if pd.isna(r.get("fraction")) else f" ({float(r['fraction']):.3f})"
        md.append(f"- {r['leakage_control']}: {int(r['n_flagged']):,}/{int(r['denominator']):,}{frac}")
    md += [
        "",
        "## Reviewer boundary",
        "Random split performance is smoke-test evidence only. Clean-science claims should be read from strict, source-heldout, patient-heldout, or no-overlap views.",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "leakage_audit.md")


if __name__ == "__main__":
    main()

