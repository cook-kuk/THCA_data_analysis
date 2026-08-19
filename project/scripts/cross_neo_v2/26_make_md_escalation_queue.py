#!/usr/bin/env python3
"""Create an MD escalation queue for prediction-fragile neoantigen candidates."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from common import OUT


TCR_OUT = OUT / "tcr_extension"
MD_OUT = TCR_OUT / "md_escalation"


def num(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce")


def minmax(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    lo, hi = s.min(skipna=True), s.max(skipna=True)
    if pd.isna(lo) or pd.isna(hi) or hi <= lo:
        return pd.Series(0.0, index=s.index)
    return (s - lo) / (hi - lo)


def md_protocol(row: pd.Series) -> tuple[str, str, str]:
    exact = "exact_paired" in str(row.get("tcr_link_claim_status", ""))
    ext_scored = pd.notna(row.get("external_tcr_expert_mean"))
    has_wt = pd.notna(row.get("wildtype_peptide")) and str(row.get("wildtype_peptide")).strip() not in {"", "nan", "None"}
    pep_len = len(str(row.get("peptide", "")))
    long_pep = pep_len > 11
    has_full_chain = bool(row.get("md_full_chain_ready", False))

    blockers: list[str] = []
    if exact and not has_full_chain:
        blockers.append("recover_full_TCR_and_MHC_sequences")
    if not exact:
        blockers.append("no_exact_paired_TCR_for_TCR-pMHC_MD")
    if not has_wt:
        blockers.append("no_WT_peptide_for_mutant-WT_delta")

    if exact and ext_scored:
        if has_full_chain:
            protocol = "TCR-pMHC_MD_3x100ns_then_optional_enhanced_MD"
        else:
            protocol = "chain_recovery_then_TCR-pMHC_MD_3x100ns"
        tier = "P0_MD_TCR_pMHC"
        rationale = "external expert support plus exact paired TCR evidence; best use of expensive MD after chain recovery"
    elif long_pep:
        protocol = "pMHC_bulge_stability_MD_3x50ns"
        tier = "P1_MD_pMHC_bulge"
        rationale = "long class-I peptide likely has bulge/conformation uncertainty"
    elif row.get("external_tcr_pairs_total", 0) == 0 or pd.isna(row.get("external_tcr_pairs_total")):
        protocol = "pMHC_stability_MD_3x50ns_plus_TCR_search"
        tier = "P1_MD_pMHC_uncertain_TCR"
        rationale = "high-priority candidate lacks modelable exact TCR rows; use pMHC stability and seek TCR evidence"
    else:
        protocol = "short_pMHC_MD_or_hold"
        tier = "P2_MD_diagnostic"
        rationale = "candidate has diagnostic value but lower immediate MD leverage"

    if has_wt and exact:
        protocol += "_plus_mutant_WT_counterfactual"
        rationale += "; WT peptide enables cross-reactivity risk test"
    return tier, protocol, rationale + "; blockers=" + ("|".join(blockers) if blockers else "none")


def main() -> None:
    MD_OUT.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(
        TCR_OUT / "model_discovery/wetlab_external_tcr_expert_scores/wetlab_candidates_external_tcr_expert_ranked.tsv",
        sep="\t",
    )
    manifest = pd.read_csv(TCR_OUT / "structure_jobs/tcr_structure_job_manifest.tsv", sep="\t")
    manifest_counts = (
        manifest.groupby("neo_row_id")
        .agg(
            md_structure_jobs=("job_id", "count"),
            md_p0_jobs=("priority", lambda s: int((s == "P0").sum())),
            md_blocked_jobs=("claim_status", lambda s: int(s.astype(str).str.contains("blocked", na=False).sum())),
            md_full_chain_ready=("tool", lambda s: bool((~s.astype(str).str.contains("blocked", na=False)).any())),
            md_example_tcr_row=("tcr_row_id", lambda s: "|".join(s.astype(str).head(3))),
        )
        .reset_index()
        .rename(columns={"neo_row_id": "row_id"})
    )
    df = candidates.merge(manifest_counts, on="row_id", how="left")
    for c in ["md_structure_jobs", "md_p0_jobs", "md_blocked_jobs"]:
        df[c] = num(df, c).fillna(0).astype(int)
    df["md_full_chain_ready"] = df["md_full_chain_ready"].apply(lambda x: bool(x) if not pd.isna(x) else False).astype(bool)

    df["value_score"] = minmax(num(df, "wetlab_priority_score_external_adjusted").fillna(num(df, "wetlab_priority_score_evidence_adjusted")))
    df["score_delta_abs"] = num(df, "score_delta_mean").abs().fillna(0)
    df["external_disagreement"] = (num(df, "external_pmtnet_mean") - num(df, "external_tepcam_mean")).abs()
    df["peptide_length"] = df["peptide"].astype(str).str.len()

    uncertainty = pd.Series(0.0, index=df.index)
    uncertainty += np.where(df["external_tcr_expert_mean"].isna(), 0.22, 0.05)
    uncertainty += np.clip(df["score_delta_abs"].fillna(0), 0, 1) * 0.22
    uncertainty += np.clip(df["external_disagreement"].fillna(0), 0, 1) * 0.18
    uncertainty += np.where(~df["tcr_link_claim_status"].astype(str).str.contains("exact_paired", na=False), 0.15, 0.03)
    uncertainty += np.where(df["md_full_chain_ready"], 0.00, 0.10)
    uncertainty += np.where(df["wildtype_peptide"].notna(), 0.00, 0.05)
    uncertainty += np.where(df["peptide_length"].gt(11) | df["peptide_length"].lt(8), 0.08, 0.00)
    uncertainty += np.where((num(df, "label_binary").eq(0)) & num(df, "tcr_augmented_score_mean").gt(0.90), 0.08, 0.00)
    df["prediction_fragility_score"] = uncertainty.clip(0, 1)

    recognition_question = pd.Series(0.0, index=df.index)
    recognition_question = np.maximum(recognition_question, np.where(df["external_tcr_expert_mean"].notna(), 0.85, 0.0))
    recognition_question = np.maximum(recognition_question, np.where(num(df, "paired_tcr_evidence_count").gt(0), 0.75, 0.0))
    recognition_question = np.maximum(recognition_question, np.where(num(df, "score_delta_mean").abs().gt(0.40), 0.65, 0.0))
    recognition_question = np.maximum(recognition_question, np.where(df["peptide_length"].gt(11), 0.55, 0.0))
    df["recognition_question_score"] = recognition_question
    df["md_escalation_score"] = (
        0.52 * df["value_score"].fillna(0)
        + 0.28 * df["prediction_fragility_score"].fillna(0)
        + 0.20 * df["recognition_question_score"].fillna(0)
    )

    protocols = df.apply(md_protocol, axis=1, result_type="expand")
    protocols.columns = ["md_tier", "recommended_md_protocol", "md_rationale"]
    df = pd.concat([df, protocols], axis=1)
    df = df.sort_values(["md_tier", "md_escalation_score"], ascending=[True, False])

    cols = [
        "row_id",
        "peptide",
        "hla_4digit",
        "label_binary",
        "pmhc_score_mean",
        "tcr_augmented_score_mean",
        "score_delta_mean",
        "wetlab_priority_score_external_adjusted",
        "external_tcr_expert_mean",
        "external_pmtnet_mean",
        "external_tepcam_mean",
        "tcr_evidence_count",
        "paired_tcr_evidence_count",
        "external_tcr_pairs_total",
        "structure_evidence_count",
        "md_structure_jobs",
        "md_blocked_jobs",
        "md_full_chain_ready",
        "value_score",
        "prediction_fragility_score",
        "recognition_question_score",
        "md_escalation_score",
        "md_tier",
        "recommended_md_protocol",
        "md_rationale",
    ]
    df[cols].to_csv(MD_OUT / "md_escalation_queue.tsv", sep="\t", index=False, na_rep="NA")
    df[cols].head(20).to_csv(MD_OUT / "md_escalation_queue_top20.tsv", sep="\t", index=False, na_rep="NA")

    tiers = pd.DataFrame(
        [
            {
                "tier": "P0_MD_TCR_pMHC",
                "when_to_use": "exact paired TCR evidence plus external pMTnet/TEPCAM support",
                "protocol": "recover full TCR/MHC chains; model TCR-pMHC; run 3 replicas x 100 ns; escalate to enhanced MD if stable",
                "claim": "diagnostic interface support only",
            },
            {
                "tier": "P1_MD_pMHC_bulge",
                "when_to_use": "long class-I peptide or conformation uncertainty",
                "protocol": "pMHC model; 3 replicas x 50 ns; monitor peptide RMSD, anchor stability, groove contacts",
                "claim": "presentation-stability diagnostic",
            },
            {
                "tier": "P1_MD_pMHC_uncertain_TCR",
                "when_to_use": "high-priority PMHC but no modelable exact TCR rows",
                "protocol": "pMHC stability MD and TCR evidence search; do not infer recognition",
                "claim": "prioritization support only",
            },
            {
                "tier": "P2_MD_diagnostic",
                "when_to_use": "lower immediate leverage or incomplete evidence",
                "protocol": "short pMHC MD or hold until additional TCR/WT evidence exists",
                "claim": "future work",
            },
        ]
    )
    tiers.to_csv(MD_OUT / "md_simulation_tiers.tsv", sep="\t", index=False)

    top = df.head(10)
    lines = [
        "# CROSS-Neo-TCR MD Escalation Queue",
        "",
        "Purpose: identify prediction-fragile, high-value neoantigen candidates that justify expensive MD or enhanced MD.",
        "",
        "Important boundary: MD is an escalation diagnostic, not proof of immunogenicity. It should be run only after structure inputs are credible.",
        "",
        "## Runtime Reality",
        "",
        "- Local GROMACS/OpenMM/MDTraj/MDAnalysis were not detected in the current Python environment.",
        "- Most TCR-pMHC structure jobs are currently blocked by missing full TCR and MHC sequences.",
        "- Therefore the immediate action is a prioritized MD queue plus chain-recovery requirements, not blind production MD.",
        "",
        "## Top MD Escalation Candidates",
        "",
        "| Rank | Peptide | HLA | Tier | MD score | Fragility | External TCR mean | Protocol |",
        "|---:|---|---|---|---:|---:|---:|---|",
    ]
    for i, r in enumerate(top.itertuples(index=False), start=1):
        ext = "NA" if pd.isna(r.external_tcr_expert_mean) else f"{r.external_tcr_expert_mean:.3f}"
        lines.append(
            f"| {i} | {r.peptide} | {r.hla_4digit} | {r.md_tier} | {r.md_escalation_score:.3f} | "
            f"{r.prediction_fragility_score:.3f} | {ext} | {r.recommended_md_protocol} |"
        )
    lines.extend(
        [
            "",
            "## Recommended Use",
            "",
            "1. Start with P0 exact paired TCR candidates: HMTEVVRHC/HLA-A*02:01 and GADGVGKSAL/HLA-C*08:02.",
            "2. Recover full TCR alpha/beta and MHC sequences before TCR-pMHC MD.",
            "3. For long peptides without exact TCR evidence, run pMHC stability MD first and keep the recognition claim off.",
            "4. Only run mutant-WT counterfactual MD when WT peptide is available.",
            "",
            "## Outputs",
            "",
            "- `md_escalation_queue.tsv`",
            "- `md_escalation_queue_top20.tsv`",
            "- `md_simulation_tiers.tsv`",
        ]
    )
    (MD_OUT / "md_escalation_report.md").write_text("\n".join(lines) + "\n")
    print(f"[md-escalation] candidates={len(df)} out={MD_OUT}")
    print(df[["peptide", "hla_4digit", "md_tier", "md_escalation_score", "prediction_fragility_score", "recommended_md_protocol"]].head(12).to_string(index=False))


if __name__ == "__main__":
    main()
