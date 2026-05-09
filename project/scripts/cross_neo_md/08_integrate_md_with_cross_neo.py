#!/usr/bin/env python3
"""Join CROSS-Neo prediction outputs with explicit-solvent MD audit evidence.

The join is intentionally conservative: MD is treated as structural audit
evidence, not a replacement label and not proof of immunogenicity.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
CROSS_OUT = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
INTEGRATED = MD_OUT / "integrated"


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def norm_peptide(x: object) -> str:
    if pd.isna(x):
        return ""
    return "".join(ch for ch in str(x).upper().strip() if ch.isalpha())


def norm_hla(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).upper().strip().replace("HLA-", "")
    s = s.replace("HLA_", "").replace("_", "*")
    if "*" not in s and len(s) >= 5:
        s = s[0] + "*" + s[1:]
    if ":" not in s and "*" in s:
        gene, allele = s.split("*", 1)
        digits = "".join(ch for ch in allele if ch.isdigit())
        if len(digits) >= 4:
            s = f"{gene}*{digits[:2]}:{digits[2:4]}"
    return "HLA-" + s if s else ""


def split_candidate(candidate: object) -> tuple[str, str]:
    text = "" if pd.isna(candidate) else str(candidate)
    if "/" not in text:
        return norm_peptide(text), ""
    peptide, hla = text.split("/", 1)
    return norm_peptide(peptide), norm_hla(hla)


def add_keys(df: pd.DataFrame, peptide_col: str = "peptide", hla_col: str = "hla_4digit") -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out["peptide_norm"] = out.get(peptide_col, "").map(norm_peptide)
    out["hla_norm"] = out.get(hla_col, "").map(norm_hla)
    out["pmhc_key"] = out["peptide_norm"] + "|" + out["hla_norm"]
    return out


def candidate_md_evidence() -> pd.DataFrame:
    md = read_tsv(MD_OUT / "md_evidence_scores.tsv")
    if md.empty:
        return pd.DataFrame()

    keys = md["candidate"].map(split_candidate)
    md["peptide_norm"] = [p for p, _ in keys]
    md["hla_norm"] = [h for _, h in keys]
    md["pmhc_key"] = md["peptide_norm"] + "|" + md["hla_norm"]
    md["is_completed_primary_10ns"] = (
        md["condition"].astype(str).str.contains("explicit_cuda_10ns", na=False)
        & (pd.to_numeric(md.get("runtime_fraction", 0), errors="coerce") >= 0.99)
        & ~md["MD_evidence_label"].astype(str).str.contains("INSUFFICIENT", na=False)
    )
    md["is_usable_md_evidence"] = ~md["MD_evidence_label"].astype(str).str.contains("INSUFFICIENT|FAIL", na=False)

    live = read_tsv(MD_OUT / "md_status_summary.tsv")
    if not live.empty:
        lk = live["candidate"].map(split_candidate)
        live["peptide_norm"] = [p for p, _ in lk]
        live["hla_norm"] = [h for _, h in lk]
        live["pmhc_key"] = live["peptide_norm"] + "|" + live["hla_norm"]
        live_agg = (
            live.groupby("pmhc_key", as_index=False)
            .agg(
                live_time_ps=("time_ps", "max"),
                live_target_ns=("target_ns", "max"),
                live_temperature_k=("temperature_k", "last"),
                live_speed_ns_per_day=("speed_ns_per_day", "max"),
            )
        )
    else:
        live_agg = pd.DataFrame(columns=["pmhc_key"])

    def pick_primary(g: pd.DataFrame) -> pd.Series:
        usable = g[g["is_completed_primary_10ns"]]
        if usable.empty:
            usable = g[g["is_usable_md_evidence"]]
        if usable.empty:
            usable = g
        chosen = usable.sort_values(
            ["is_completed_primary_10ns", "runtime_fraction", "MD_evidence_score"],
            ascending=[False, False, False],
        ).iloc[0]
        return pd.Series(
            {
                "candidate": chosen["candidate"],
                "md_representative_run_id": chosen["run_id"],
                "md_representative_condition": chosen["condition"],
                "MD_evidence_score": chosen["MD_evidence_score"],
                "MD_evidence_label": chosen["MD_evidence_label"],
                "pMHC_stability_score": chosen["pMHC_stability_score"],
                "TCR_recognition_score": chosen["TCR_recognition_score"],
                "counterfactual_specificity_score": chosen["counterfactual_specificity_score"],
                "replicate_confidence_score": chosen["replicate_confidence_score"],
                "simulation_qc_score": chosen["simulation_qc_score"],
                "runtime_fraction": chosen["runtime_fraction"],
                "peptide_rmsd_final_nm": chosen["peptide_rmsd_final_nm"],
                "pmhc_native_q_tail": chosen["pmhc_native_q_tail"],
                "tcr_peptide_contacts_tail": chosen["tcr_peptide_contacts_tail"],
                "any_completed_primary_10ns": bool(g["is_completed_primary_10ns"].any()),
                "n_md_runs_seen": int(len(g)),
                "best_md_evidence_score_seen": pd.to_numeric(g["MD_evidence_score"], errors="coerce").max(),
            }
        )

    agg = md.groupby("pmhc_key").apply(pick_primary, include_groups=False).reset_index()
    if not live_agg.empty:
        agg = agg.merge(live_agg, on="pmhc_key", how="left")
    agg["md_interpretation"] = np.select(
        [
            agg["any_completed_primary_10ns"].fillna(False) & (agg["MD_evidence_label"].eq("MD_MODERATE") | agg["MD_evidence_label"].eq("MD_STRONG") | agg["MD_evidence_label"].eq("MD_VERY_STRONG")),
            agg["MD_evidence_label"].astype(str).str.contains("INSUFFICIENT", na=False),
            ~agg["any_completed_primary_10ns"].fillna(False),
        ],
        [
            "usable_structural_support_not_immunogenicity_proof",
            "insufficient_runtime_or_missing_trajectory",
            "partial_or_replicate_only_evidence",
        ],
        default="structural_audit_evidence_only",
    )
    return agg


def load_predictions_with_registry() -> pd.DataFrame:
    registry = read_tsv(CROSS_OUT / "canonical_registry.tsv")
    if registry.empty:
        return pd.DataFrame()
    keep_cols = [
        c
        for c in [
            "row_id",
            "source_dataset",
            "study_id",
            "assay_type",
            "cancer_type",
            "peptide",
            "wildtype_peptide",
            "hla_4digit",
            "hla_supertype",
            "peptide_length",
            "source_protein",
            "label_binary",
            "public_overlap_flags",
        ]
        if c in registry.columns
    ]
    registry = add_keys(registry[keep_cols])

    pred_paths = [
        CROSS_OUT / "predictions/all_predictions.tsv",
        CROSS_OUT / "predictions/selective_ensemble_v2_1_predictions.tsv",
        CROSS_OUT / "predictions/fast_esm2_qk_gate_predictions.tsv",
    ]
    frames = []
    for path in pred_paths:
        pred = read_tsv(path)
        if pred.empty or "row_id" not in pred.columns:
            continue
        pred = pred.copy()
        pred["prediction_file"] = path.name
        frames.append(pred.merge(registry, on="row_id", how="left", suffixes=("", "_registry")))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def summarize_joined(joined: pd.DataFrame, wetlab: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if joined.empty:
        cols = ["row_id", "peptide", "hla_4digit", "score", "MD_evidence_score"]
        return pd.DataFrame(columns=cols), pd.DataFrame(columns=cols), pd.DataFrame(columns=cols)

    joined["score_numeric"] = pd.to_numeric(joined.get("score"), errors="coerce")
    joined["MD_evidence_score_numeric"] = pd.to_numeric(joined.get("MD_evidence_score"), errors="coerce")
    joined["rank_numeric"] = pd.to_numeric(joined.get("rank"), errors="coerce")
    joined["rank_pct_numeric"] = pd.to_numeric(joined.get("rank_pct"), errors="coerce")

    supported = joined[
        (joined["score_numeric"] >= 0.5)
        & (joined["MD_evidence_score_numeric"] >= 0.5)
        & joined["MD_evidence_label"].astype(str).isin(["MD_MODERATE", "MD_STRONG", "MD_VERY_STRONG"])
    ].copy()
    supported = supported.sort_values(["MD_evidence_score_numeric", "score_numeric"], ascending=False)

    high_md_low = joined[
        (joined["score_numeric"] >= 0.75)
        & (
            (joined["MD_evidence_score_numeric"] < 0.45)
            | joined["MD_evidence_label"].astype(str).str.contains("INSUFFICIENT|FAIL", na=False)
        )
    ].copy()
    high_md_low = high_md_low.sort_values(["score_numeric"], ascending=False)

    low_md_high = joined[
        (joined["score_numeric"] < 0.35)
        & (joined["MD_evidence_score_numeric"] >= 0.55)
        & joined["MD_evidence_label"].astype(str).isin(["MD_MODERATE", "MD_STRONG", "MD_VERY_STRONG"])
    ].copy()
    low_md_high = low_md_high.sort_values(["MD_evidence_score_numeric", "score_numeric"], ascending=[False, True])

    if not wetlab.empty:
        wetlab_cols = [
            c
            for c in [
                "row_id",
                "pmhc_key",
                "wetlab_priority_score",
                "wetlab_priority_score_evidence_adjusted",
                "wetlab_priority_tier",
                "pmhc_score_mean",
                "tcr_augmented_score_mean",
                "score_delta_mean",
                "tcr_link_claim_status",
                "tcr_link_warning",
            ]
            if c in wetlab.columns
        ]
        for frame_name, frame in [("supported", supported), ("high", high_md_low), ("low", low_md_high)]:
            if not frame.empty:
                merged = frame.merge(wetlab[wetlab_cols], on=["row_id", "pmhc_key"], how="left", suffixes=("", "_wetlab"))
                if frame_name == "supported":
                    supported = merged
                elif frame_name == "high":
                    high_md_low = merged
                else:
                    low_md_high = merged

    return supported, high_md_low, low_md_high


def write_report(joined: pd.DataFrame, md: pd.DataFrame, supported: pd.DataFrame, high_low: pd.DataFrame, low_high: pd.DataFrame) -> None:
    lines = [
        "# CROSS-Neo + MD Integration Audit",
        "",
        f"Generated from `{MD_OUT}` and `{CROSS_OUT}`.",
        "",
        "## Boundary",
        "",
        "MD evidence is used as structural audit evidence only. It is not an immunogenicity label, clinical validation, or SOTA proof.",
        "",
        "## Candidate-Level MD Evidence",
        "",
    ]
    if md.empty:
        lines.append("No MD evidence table was available.")
    else:
        cols = [
            "candidate",
            "md_representative_run_id",
            "MD_evidence_label",
            "MD_evidence_score",
            "pMHC_stability_score",
            "TCR_recognition_score",
            "runtime_fraction",
            "live_time_ps",
            "md_interpretation",
        ]
        lines.append(md[[c for c in cols if c in md.columns]].to_markdown(index=False))
    lines += [
        "",
        "## Joined Prediction Rows",
        "",
        f"- prediction rows with MD candidate key: {len(joined)}",
        f"- MD-supported model-high rows: {len(supported)}",
        f"- model-high / MD-low-or-insufficient rows: {len(high_low)}",
        f"- model-low / MD-high rows: {len(low_high)}",
        "",
        "## Interpretation",
        "",
        "- `GADGVGKSAL / HLA-C*08:02` has completed 10 ns explicit-solvent evidence and a 1 ns replicate screen, so it can be discussed as structural plausibility/wetlab prioritization evidence if the contact and QC plots remain stable.",
        "- `HMTEVVRHC / HLA-A*02:01` is still running in the primary 10 ns job; state traces are useful for live QC, but trajectory-derived contact claims must wait for DCD sync and rerun of the analysis.",
        "- No WT or scrambled-decoy trajectory is present yet, so mutant-specific TCR recognition and cross-reactivity claims remain forbidden.",
        "",
        "## Outputs",
        "",
        "- `integrated/cross_neo_md_joined.tsv`",
        "- `integrated/md_supported_top_candidates.tsv`",
        "- `integrated/model_high_md_low_cases.tsv`",
        "- `integrated/model_low_md_high_cases.tsv`",
    ]
    (INTEGRATED / "md_case_audit_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    INTEGRATED.mkdir(parents=True, exist_ok=True)
    md = candidate_md_evidence()
    pred = load_predictions_with_registry()
    if md.empty or pred.empty:
        joined = pd.DataFrame()
    else:
        joined = pred.merge(md, on="pmhc_key", how="inner", suffixes=("", "_md"))

    wetlab = read_tsv(CROSS_OUT / "tcr_extension/tcr_wetlab_candidate_prioritization_unique_pmhc.tsv")
    if not wetlab.empty:
        wetlab = add_keys(wetlab)
        wetlab = wetlab.drop_duplicates(["row_id", "pmhc_key"])
        if not joined.empty:
            cols = [
                c
                for c in [
                    "row_id",
                    "pmhc_key",
                    "pmhc_score_mean",
                    "tcr_augmented_score_mean",
                    "score_delta_mean",
                    "wetlab_priority_score",
                    "wetlab_priority_score_evidence_adjusted",
                    "wetlab_priority_tier",
                    "tcr_link_claim_status",
                    "tcr_link_warning",
                ]
                if c in wetlab.columns
            ]
            joined = joined.merge(wetlab[cols], on=["row_id", "pmhc_key"], how="left", suffixes=("", "_wetlab"))

    supported, high_low, low_high = summarize_joined(joined, wetlab)

    joined.to_csv(INTEGRATED / "cross_neo_md_joined.tsv", sep="\t", index=False)
    supported.to_csv(INTEGRATED / "md_supported_top_candidates.tsv", sep="\t", index=False)
    high_low.to_csv(INTEGRATED / "model_high_md_low_cases.tsv", sep="\t", index=False)
    low_high.to_csv(INTEGRATED / "model_low_md_high_cases.tsv", sep="\t", index=False)

    summary = {
        "joined_prediction_rows": int(len(joined)),
        "md_supported_top_candidate_rows": int(len(supported)),
        "model_high_md_low_rows": int(len(high_low)),
        "model_low_md_high_rows": int(len(low_high)),
        "candidate_md_rows": int(len(md)),
    }
    (INTEGRATED / "cross_neo_md_integration_summary.json").write_text(json.dumps(summary, indent=2))
    write_report(joined, md, supported, high_low, low_high)
    print("[md-integrate]", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
