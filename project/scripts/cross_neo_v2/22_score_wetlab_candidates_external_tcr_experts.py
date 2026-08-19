#!/usr/bin/env python3
"""Attach pMTnet/TEPCAM diagnostic scores to top wetlab candidate PMHC rows."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from common import OUT, SEED


TCR_OUT = OUT / "tcr_extension"
MODEL_OUT = TCR_OUT / "model_discovery" / "wetlab_external_tcr_expert_scores"
PMTNET = Path("/data/thca/tcr_model_repos/pMTnet")
TEPCAM = Path("/data/thca/tcr_model_repos/TEPCAM")
TOP_N = int(os.environ.get("CROSS_NEO_TCR_WETLAB_TOP_N", "25"))
MAX_TCR_PER_PMHC = int(os.environ.get("CROSS_NEO_TCR_MAX_TCR_PER_PMHC", "20"))


def run(cmd: list[str], cwd: Path, stdout_path: Path, stderr_path: Path) -> int:
    proc = subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True, check=False)
    stdout_path.write_text(proc.stdout)
    stderr_path.write_text(proc.stderr)
    return int(proc.returncode)


def clean_hla_for_pmtnet(hla: object) -> str:
    return str(hla or "").replace("HLA-", "").strip()


def select_pairs() -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = pd.read_csv(TCR_OUT / "tcr_wetlab_candidate_prioritization_unique_pmhc.tsv", sep="\t")
    candidates = candidates[
        candidates["peptide"].fillna("").astype(str).str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{8,15}")
        & candidates["hla_4digit"].fillna("").astype(str).str.contains(r"HLA-[ABC]\*", regex=True)
    ].copy()
    candidates = candidates.sort_values("wetlab_priority_score_evidence_adjusted", ascending=False).head(TOP_N)
    candidates["pmhc_key"] = candidates["peptide"].astype(str) + "|" + candidates["hla_4digit"].astype(str)

    tcr = pd.read_parquet(
        TCR_OUT / "tcr_registry.parquet",
        columns=[
            "row_id",
            "source_dataset",
            "peptide",
            "hla_4digit",
            "cdr3_alpha",
            "cdr3_beta",
            "paired_tcr_available",
            "binding_label_binary",
            "is_cancer_context",
            "is_pathogen_context",
            "antigen_source",
        ],
    )
    tcr = tcr[
        tcr["cdr3_beta"].fillna("").astype(str).str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")
        & tcr["peptide"].fillna("").astype(str).str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{8,15}")
        & tcr["hla_4digit"].fillna("").astype(str).str.contains(r"HLA-[ABC]\*", regex=True)
    ].copy()
    tcr["pmhc_key"] = tcr["peptide"].astype(str) + "|" + tcr["hla_4digit"].astype(str)
    tcr = tcr[tcr["pmhc_key"].isin(set(candidates["pmhc_key"]))].copy()
    for col in ["paired_tcr_available", "is_cancer_context", "is_pathogen_context"]:
        tcr[col] = tcr[col].fillna(False).astype(bool)
    tcr["binding_label_binary"] = pd.to_numeric(tcr["binding_label_binary"], errors="coerce")
    tcr["priority_sort"] = (
        tcr["binding_label_binary"].fillna(0).astype(float) * 4
        + tcr["is_cancer_context"].astype(float) * 2
        + tcr["paired_tcr_available"].astype(float)
    )
    selected = []
    for pmhc_key, sub in tcr.sort_values("priority_sort", ascending=False).groupby("pmhc_key"):
        sub = sub.drop_duplicates(["cdr3_beta", "peptide", "hla_4digit"])
        if len(sub) > MAX_TCR_PER_PMHC:
            # Keep the highest-priority rows, then diversify sources.
            head = sub.head(MAX_TCR_PER_PMHC * 2)
            pieces = []
            per_source = max(1, int(np.ceil(MAX_TCR_PER_PMHC / max(1, head["source_dataset"].nunique()))))
            for _, g in head.groupby("source_dataset"):
                pieces.append(g.head(per_source))
            sub = pd.concat(pieces, ignore_index=True).head(MAX_TCR_PER_PMHC)
        selected.append(sub)
    pairs = pd.concat(selected, ignore_index=True) if selected else pd.DataFrame()
    pairs = pairs.merge(candidates[["row_id", "pmhc_key", "wetlab_priority_score_evidence_adjusted", "wetlab_priority_tier"]].rename(columns={"row_id": "candidate_row_id"}), on="pmhc_key", how="left")
    pairs["eval_id"] = [f"wetlab_ext_{i:05d}" for i in range(len(pairs))]
    pairs["hla_for_pmtnet"] = pairs["hla_4digit"].map(clean_hla_for_pmtnet)
    return candidates, pairs


def run_pmtnet(pairs: pd.DataFrame) -> pd.DataFrame:
    uniq = pairs[["cdr3_beta", "peptide", "hla_for_pmtnet"]].drop_duplicates()
    inp = MODEL_OUT / "pmtnet_input_unique.csv"
    uniq.rename(columns={"cdr3_beta": "CDR3", "peptide": "Antigen", "hla_for_pmtnet": "HLA"}).to_csv(inp, index=False)
    shutil.rmtree(MODEL_OUT / "pmtnet_run", ignore_errors=True)
    (MODEL_OUT / "pmtnet_run").mkdir(parents=True, exist_ok=True)
    rc = run(
        [
            "python",
            str(PMTNET / "pMTnet.py"),
            "-input",
            str(inp),
            "-library",
            str(PMTNET / "library"),
            "-output",
            str(MODEL_OUT / "pmtnet_run"),
            "-output_log",
            str(MODEL_OUT / "pmtnet_output.log"),
        ],
        PMTNET,
        MODEL_OUT / "pmtnet_stdout.txt",
        MODEL_OUT / "pmtnet_stderr.txt",
    )
    if rc != 0:
        pairs["pmtnet_score"] = np.nan
        return pairs
    pred = pd.read_csv(MODEL_OUT / "pmtnet_run/prediction.csv")
    pred = pred.rename(columns={"CDR3": "cdr3_beta", "Antigen": "peptide", "HLA": "hla_for_pmtnet"})
    pred["pmtnet_rank"] = pd.to_numeric(pred["Rank"], errors="coerce")
    pred["pmtnet_score"] = 1.0 - pred["pmtnet_rank"]
    return pairs.merge(pred[["cdr3_beta", "peptide", "hla_for_pmtnet", "pmtnet_rank", "pmtnet_score"]], on=["cdr3_beta", "peptide", "hla_for_pmtnet"], how="left")


def run_tepcam(pairs: pd.DataFrame) -> pd.DataFrame:
    uniq = pairs[["cdr3_beta", "peptide"]].drop_duplicates().reset_index(drop=True)
    inp = MODEL_OUT / "tepcam_input_unique.csv"
    out = MODEL_OUT / "tepcam_output_unique.csv"
    metrics = MODEL_OUT / "tepcam_metrics_raw.txt"
    for p in [out, metrics]:
        if p.exists():
            p.unlink()
    pd.DataFrame({"TCR": uniq["cdr3_beta"], "epitope": uniq["peptide"], "Label": 0}).to_csv(inp, index=False)
    rc = run(
        [
            "python",
            "scripts/test.py",
            "--file_path",
            str(inp),
            "--model_path",
            "./ckpts/tepcam_test.pt",
            "--output_file",
            str(out),
            "--metric_file",
            str(metrics),
            "--batch_size",
            "256",
        ],
        TEPCAM,
        MODEL_OUT / "tepcam_stdout.txt",
        MODEL_OUT / "tepcam_stderr.txt",
    )
    if rc != 0:
        pairs["tepcam_score"] = np.nan
        return pairs
    pred = pd.read_csv(out)
    if len(pred) != len(uniq):
        pairs["tepcam_score"] = np.nan
        return pairs
    pred = pd.DataFrame({"cdr3_beta": uniq["cdr3_beta"], "peptide": uniq["peptide"], "tepcam_score": pd.to_numeric(pred["prediction"], errors="coerce")})
    return pairs.merge(pred, on=["cdr3_beta", "peptide"], how="left")


def aggregate(candidates: pd.DataFrame, pairs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pairs["external_tcr_expert_mean"] = pairs[["pmtnet_score", "tepcam_score"]].mean(axis=1, skipna=True)
    agg = (
        pairs.groupby("pmhc_key")
        .agg(
            external_tcr_pairs_scored=("external_tcr_expert_mean", lambda s: int(s.notna().sum())),
            external_tcr_pairs_total=("row_id", "count"),
            external_pmtnet_mean=("pmtnet_score", "mean"),
            external_pmtnet_max=("pmtnet_score", "max"),
            external_tepcam_mean=("tepcam_score", "mean"),
            external_tepcam_max=("tepcam_score", "max"),
            external_tcr_expert_mean=("external_tcr_expert_mean", "mean"),
            external_tcr_expert_max=("external_tcr_expert_mean", "max"),
            external_positive_label_count=("binding_label_binary", lambda s: int(pd.to_numeric(s, errors="coerce").eq(1).sum())),
            external_negative_label_count=("binding_label_binary", lambda s: int(pd.to_numeric(s, errors="coerce").eq(0).sum())),
            external_paired_count=("paired_tcr_available", "sum"),
            external_cancer_context_count=("is_cancer_context", "sum"),
            external_pathogen_context_count=("is_pathogen_context", "sum"),
            external_example_tcrb=("cdr3_beta", lambda s: "|".join(s.astype(str).head(5))),
            external_sources=("source_dataset", lambda s: "|".join(sorted(set(s.astype(str))))),
        )
        .reset_index()
    )
    out = candidates.merge(agg, on="pmhc_key", how="left")
    out["wetlab_priority_score_external_adjusted"] = (
        pd.to_numeric(out["wetlab_priority_score_evidence_adjusted"], errors="coerce").fillna(0)
        + 0.35 * pd.to_numeric(out["external_tcr_expert_mean"], errors="coerce").fillna(0)
        + 0.15 * pd.to_numeric(out["external_tcr_expert_max"], errors="coerce").fillna(0)
    )
    out = out.sort_values("wetlab_priority_score_external_adjusted", ascending=False)
    return out, pairs


def write_report(scored: pd.DataFrame, pairs: pd.DataFrame) -> None:
    top = scored.head(15)
    def fmt(x: object) -> str:
        return "NA" if pd.isna(x) else f"{float(x):.4f}"
    def fmt_int(x: object) -> str:
        return "0" if pd.isna(x) else str(int(x))
    lines = [
        "# Wetlab Candidate External TCR Expert Scores",
        "",
        f"Top {TOP_N} PMHC candidates were rescored using exact peptide-HLA TCR rows, up to {MAX_TCR_PER_PMHC} TCRs per PMHC.",
        "",
        "Claim boundary: these are diagnostic expert scores over available public TCR evidence. They are not transferred neoantigen labels.",
        "",
        "| Rank | Peptide | HLA | Old priority | External mean | External max | External-adjusted priority | TCR pairs | Sources |",
        "|---:|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for i, r in enumerate(top.itertuples(index=False), start=1):
        lines.append(
            f"| {i} | {r.peptide} | {r.hla_4digit} | {r.wetlab_priority_score_evidence_adjusted:.4f} | "
            f"{fmt(getattr(r, 'external_tcr_expert_mean', np.nan))} | {fmt(getattr(r, 'external_tcr_expert_max', np.nan))} | "
            f"{fmt(r.wetlab_priority_score_external_adjusted)} | {fmt_int(getattr(r, 'external_tcr_pairs_total', 0))} | "
            f"{'' if pd.isna(getattr(r, 'external_sources', np.nan)) else getattr(r, 'external_sources', '')} |"
        )
    lines.extend(
        [
            "",
            "Decision: prioritize candidates that remain high after external expert adjustment, especially when pMTnet and TEPCAM agree and exact TCR evidence is cancer-context or paired.",
        ]
    )
    (MODEL_OUT / "05_wetlab_external_tcr_expert_scores.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    MODEL_OUT.mkdir(parents=True, exist_ok=True)
    candidates, pairs = select_pairs()
    pairs.to_csv(MODEL_OUT / "wetlab_external_tcr_expert_input_pairs.tsv", sep="\t", index=False, na_rep="NA")
    if pairs.empty:
        raise SystemExit("No wetlab candidate TCR pairs found for external scoring.")
    pairs = run_pmtnet(pairs)
    pairs = run_tepcam(pairs)
    scored, pair_scores = aggregate(candidates, pairs)
    pair_scores.to_csv(MODEL_OUT / "wetlab_external_tcr_expert_pair_scores.tsv", sep="\t", index=False, na_rep="NA")
    scored.to_csv(MODEL_OUT / "wetlab_candidates_external_tcr_expert_ranked.tsv", sep="\t", index=False, na_rep="NA")
    write_report(scored, pair_scores)
    print(scored[["peptide", "hla_4digit", "external_tcr_pairs_total", "external_tcr_expert_mean", "external_tcr_expert_max", "wetlab_priority_score_external_adjusted"]].head(15).to_string(index=False))
    print(f"[wetlab-external-tcr-experts] candidates={len(scored)} pairs={len(pair_scores)} out={MODEL_OUT}")


if __name__ == "__main__":
    main()
