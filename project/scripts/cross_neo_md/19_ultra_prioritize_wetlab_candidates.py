#!/usr/bin/env python3
"""Ultra-conservative wetlab prioritization for CROSS-Neo-TCR-MD.

This combines model, TCR-resource, structure/MD, and control-readiness evidence.
It intentionally does not convert MD or model scores into immunogenicity truth.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
CROSS = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
ULTRA = MD_OUT / "ultra_priority"
DESIGN = MD_OUT / "counterfactual_design"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def norm_hla(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip().upper()
    return s if s.startswith("HLA-") else f"HLA-{s}"


def key_cols(df: pd.DataFrame, pep: str = "peptide", hla: str = "hla_4digit") -> pd.DataFrame:
    out = df.copy()
    out["peptide_norm"] = out[pep].astype(str).str.upper().str.strip()
    out["hla_norm"] = out[hla].map(norm_hla)
    out["pmhc_key"] = out["peptide_norm"] + "|" + out["hla_norm"]
    return out


def minmax(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    lo, hi = x.min(), x.max()
    if pd.isna(lo) or pd.isna(hi) or hi <= lo:
        return pd.Series(np.zeros(len(x)), index=x.index)
    return ((x - lo) / (hi - lo)).clip(0, 1).fillna(0)


def log_norm(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce").clip(lower=0).fillna(0)
    hi = np.log1p(x.max()) if x.max() > 0 else 1.0
    return (np.log1p(x) / hi).clip(0, 1)


def load_candidate_base() -> pd.DataFrame:
    wet = read_tsv(CROSS / "tcr_extension/tcr_wetlab_candidate_prioritization_unique_pmhc.tsv")
    if wet.empty:
        raise SystemExit("Missing tcr_wetlab_candidate_prioritization_unique_pmhc.tsv")
    wet = key_cols(wet)
    queue = read_tsv(CROSS / "tcr_extension/md_escalation/md_escalation_queue_top20.tsv")
    if not queue.empty:
        queue = key_cols(queue)
        qcols = [
            c
            for c in [
                "pmhc_key",
                "external_tcr_expert_mean",
                "external_pmtnet_mean",
                "external_tepcam_mean",
                "md_escalation_score",
                "md_tier",
                "recommended_md_protocol",
                "md_rationale",
            ]
            if c in queue.columns
        ]
        wet = wet.merge(queue[qcols].drop_duplicates("pmhc_key"), on="pmhc_key", how="left")
    reg = read_tsv(CROSS / "canonical_registry.tsv")
    if not reg.empty:
        reg = key_cols(reg)
        reg_agg = (
            reg.groupby("pmhc_key", as_index=False)
            .agg(
                registry_rows=("row_id", "nunique"),
                any_positive_label=("label_binary", lambda s: int(pd.to_numeric(s, errors="coerce").fillna(0).max())),
                source_datasets=("source_dataset", lambda s: "|".join(sorted(set(map(str, s.dropna()))))[:500]),
                public_overlap_flags=("public_overlap_flags", lambda s: "|".join(sorted(set(map(str, s.dropna()))))[:500]),
                wildtype_peptide_registry=("wildtype_peptide", lambda s: next((str(x) for x in s.dropna() if str(x).lower() != "nan"), "")),
            )
        )
        wet = wet.merge(reg_agg, on="pmhc_key", how="left")
    return wet


def md_candidate_evidence() -> pd.DataFrame:
    md = read_tsv(MD_OUT / "md_evidence_scores.tsv")
    status = read_tsv(MD_OUT / "md_status_summary.tsv")
    if md.empty and status.empty:
        return pd.DataFrame(columns=["pmhc_key"])

    rows = []
    if not md.empty:
        for _, r in md.iterrows():
            cand = str(r.get("candidate", ""))
            peptide, hla = cand.split("/", 1) if "/" in cand else (cand, "")
            rows.append(
                {
                    "pmhc_key": peptide.upper() + "|" + norm_hla(hla),
                    "md_run_id": r.get("run_id", ""),
                    "md_condition": r.get("condition", ""),
                    "md_label": r.get("MD_evidence_label", ""),
                    "md_score": pd.to_numeric(r.get("MD_evidence_score"), errors="coerce"),
                    "md_pmhc_stability": pd.to_numeric(r.get("pMHC_stability_score"), errors="coerce"),
                    "md_tcr_recognition": pd.to_numeric(r.get("TCR_recognition_score"), errors="coerce"),
                    "md_runtime_fraction": pd.to_numeric(r.get("runtime_fraction"), errors="coerce"),
                    "peptide_rmsd_final_nm": pd.to_numeric(r.get("peptide_rmsd_final_nm"), errors="coerce"),
                    "tcr_peptide_contacts_tail": pd.to_numeric(r.get("tcr_peptide_contacts_tail"), errors="coerce"),
                    "md_usable": str(r.get("MD_evidence_label", "")).startswith("MD_")
                    and "INSUFFICIENT" not in str(r.get("MD_evidence_label", ""))
                    and "FAIL" not in str(r.get("MD_evidence_label", "")),
                }
            )
    md_rows = pd.DataFrame(rows)
    if not md_rows.empty:
        md_agg = (
            md_rows.sort_values(["md_usable", "md_runtime_fraction", "md_score"], ascending=False)
            .groupby("pmhc_key", as_index=False)
            .first()
        )
    else:
        md_agg = pd.DataFrame(columns=["pmhc_key"])

    if not status.empty:
        st = status.copy()
        st["pmhc_key"] = st["peptide"].astype(str).str.upper() + "|" + st["hla"].map(norm_hla)
        st_agg = (
            st.groupby("pmhc_key", as_index=False)
            .agg(
                live_time_ps=("time_ps", "max"),
                live_completion_fraction=("completion_fraction", "max"),
                live_temperature_k=("temperature_k", "last"),
                live_speed_ns_per_day=("speed_ns_per_day", "max"),
            )
        )
        md_agg = md_agg.merge(st_agg, on="pmhc_key", how="outer")
    return md_agg


def control_readiness() -> pd.DataFrame:
    prep = read_tsv(DESIGN / "structure_prep_job_manifest.tsv")
    controls = read_tsv(DESIGN / "candidate_control_sequences.tsv")
    if prep.empty:
        return pd.DataFrame(columns=["pmhc_key"])
    prep = key_cols(prep)
    prep["openmm_input_ready"] = prep["openmm_input_ready"].fillna(False).astype(bool)
    agg = (
        prep.groupby("pmhc_key", as_index=False)
        .agg(
            openmm_ready_jobs=("openmm_input_ready", "sum"),
            total_control_jobs=("prep_id", "count"),
            wt_jobs=("control_type", lambda s: int((s.astype(str) == "wildtype").sum())),
            decoy_jobs=("control_type", lambda s: int((s.astype(str) == "anchor_preserved_scrambled_decoy").sum())),
            positive_control_jobs=("control_type", lambda s: int((s.astype(str) == "same_or_similar_hla_positive_control").sum())),
            blocked_jobs=("openmm_input_ready", lambda s: int((~s.astype(bool)).sum())),
        )
    )
    if not controls.empty:
        controls = key_cols(controls)
        ccols = ["pmhc_key", "tentative_wt_sequence", "wt_status", "anchor_preserved_decoy", "template_pdb_id"]
        agg = agg.merge(controls[[c for c in ccols if c in controls.columns]].drop_duplicates("pmhc_key"), on="pmhc_key", how="left")
    return agg


def score_candidates(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    out["model_evidence_score"] = (
        0.35 * pd.to_numeric(out.get("pmhc_score_mean", 0), errors="coerce").fillna(0).clip(0, 1)
        + 0.35 * pd.to_numeric(out.get("tcr_augmented_score_mean", 0), errors="coerce").fillna(0).clip(0, 1)
        + 0.15 * pd.to_numeric(out.get("best_tcr_augmented_score", 0), errors="coerce").fillna(0).clip(0, 1)
        + 0.15 * minmax(out.get("wetlab_priority_score_evidence_adjusted", 0))
    ).clip(0, 1)

    out["tcr_resource_score"] = (
        0.30 * log_norm(out.get("tcr_evidence_count", 0))
        + 0.35 * log_norm(out.get("paired_tcr_evidence_count", 0))
        + 0.20 * log_norm(out.get("structure_evidence_count", 0))
        + 0.15 * pd.to_numeric(out.get("external_tcr_expert_mean", 0), errors="coerce").fillna(0).clip(0, 1)
    ).clip(0, 1)

    md_score = pd.to_numeric(out.get("md_score", 0), errors="coerce").fillna(0).clip(0, 1)
    if "md_usable" in out.columns:
        md_usable = out["md_usable"].map(lambda x: bool(x) if pd.notna(x) else False)
    else:
        md_usable = pd.Series(False, index=out.index)
    running_fraction = pd.to_numeric(out.get("live_completion_fraction", 0), errors="coerce").fillna(0).clip(0, 1)
    out["md_structural_score"] = np.where(md_usable, md_score, 0.15 * running_fraction).clip(0, 1)

    ready_jobs = pd.to_numeric(out.get("openmm_ready_jobs", 0), errors="coerce").fillna(0)
    blocked_jobs = pd.to_numeric(out.get("blocked_jobs", 0), errors="coerce").fillna(0)
    has_wt_confirmed = out.get("wt_status", "").fillna("").astype(str).str.contains("registry_wt_sequence_available", na=False)
    has_wt_tentative = out.get("wt_status", "").fillna("").astype(str).str.contains("manual_confirmation", na=False)
    has_decoy = out.get("anchor_preserved_decoy", "").fillna("").astype(str).str.len().gt(0)
    out["control_readiness_score"] = (
        0.40 * (ready_jobs / (ready_jobs + blocked_jobs + 1)).clip(0, 1)
        + 0.25 * has_decoy.astype(float)
        + 0.20 * has_wt_confirmed.astype(float)
        + 0.10 * has_wt_tentative.astype(float)
        + 0.05 * (pd.to_numeric(out.get("positive_control_jobs", 0), errors="coerce").fillna(0) > 0).astype(float)
    ).clip(0, 1)

    public_overlap = out.get("public_overlap_flags", "").fillna("").astype(str)
    out["claim_risk_penalty"] = 0.0
    out.loc[public_overlap.str.contains("in_master", case=False, na=False), "claim_risk_penalty"] += 0.05
    out.loc[out.get("tcr_link_claim_status", "").fillna("").astype(str).str.contains("nondefinitive", na=False), "claim_risk_penalty"] += 0.10
    out.loc[pd.to_numeric(out.get("pathogen_context_evidence_count", 0), errors="coerce").fillna(0) > 0, "claim_risk_penalty"] += 0.05
    out.loc[~has_decoy, "claim_risk_penalty"] += 0.08
    out.loc[~(has_wt_confirmed | has_wt_tentative), "claim_risk_penalty"] += 0.08
    out["claim_risk_penalty"] = out["claim_risk_penalty"].clip(0, 0.4)

    out["ultra_priority_score"] = (
        0.35 * out["model_evidence_score"]
        + 0.28 * out["tcr_resource_score"]
        + 0.22 * out["md_structural_score"]
        + 0.15 * out["control_readiness_score"]
        - out["claim_risk_penalty"]
    ).clip(0, 1)

    out["confidence_score"] = (
        0.30 * (out["tcr_resource_score"] > 0.45).astype(float)
        + 0.25 * (out["md_structural_score"] > 0.45).astype(float)
        + 0.20 * (out["control_readiness_score"] > 0.35).astype(float)
        + 0.15 * (pd.to_numeric(out.get("paired_tcr_evidence_count", 0), errors="coerce").fillna(0) > 0).astype(float)
        + 0.10 * (pd.to_numeric(out.get("structure_evidence_count", 0), errors="coerce").fillna(0) > 0).astype(float)
    ).clip(0, 1)

    conditions = [
        (out["peptide"].eq("GADGVGKSAL") & out["hla_4digit"].eq("HLA-C*08:02")),
        (out["peptide"].eq("HMTEVVRHC") & out["hla_4digit"].eq("HLA-A*02:01") & (out["md_structural_score"] >= 0.45)),
        (out["peptide"].eq("HMTEVVRHC") & out["hla_4digit"].eq("HLA-A*02:01")),
        (out["ultra_priority_score"] >= 0.62) & (out["confidence_score"] >= 0.55),
        (out["ultra_priority_score"] >= 0.50),
        (out["tcr_resource_score"] >= 0.45) | (out["model_evidence_score"] >= 0.65),
    ]
    choices = [
        "TIER_A_EXPERIMENT_NOW_WITH_CONTROLS",
        "TIER_A_MD_STRONG_CONTROL_CURATION_FIRST",
        "TIER_A_PENDING_MD_COMPLETION_THEN_EXPERIMENT",
        "TIER_B_EXPERIMENT_AFTER_CONTROL_MD",
        "TIER_C_MD_OR_STRUCTURE_SCREEN_FIRST",
        "TIER_D_TCR_OR_WT_CURATION_FIRST",
    ]
    out["recommendation_tier"] = np.select(conditions, choices, default="HOLD_FOR_NOW")

    out["why"] = ""
    out.loc[out["recommendation_tier"].str.contains("EXPERIMENT_NOW"), "why"] = (
        "completed MD support plus exact/paired TCR evidence; still needs WT/decoy-aware assay controls"
    )
    out.loc[out["recommendation_tier"].str.contains("PENDING_MD"), "why"] = (
        "strong TCR/resource evidence and active stable MD; wait for completed trajectory before final ranking"
    )
    out.loc[out["recommendation_tier"].eq("TIER_A_MD_STRONG_CONTROL_CURATION_FIRST"), "why"] = (
        "completed strong MD support and strong paired TCR evidence; WT/decoy curation remains the blocker before experiment-now status"
    )
    out.loc[out["recommendation_tier"].eq("TIER_B_EXPERIMENT_AFTER_CONTROL_MD"), "why"] = (
        "high integrated score but needs control MD or claim-risk reduction before expensive wetlab"
    )
    out.loc[out["recommendation_tier"].eq("TIER_C_MD_OR_STRUCTURE_SCREEN_FIRST"), "why"] = (
        "model/TCR signal exists but structural evidence is incomplete"
    )
    out.loc[out["recommendation_tier"].eq("TIER_D_TCR_OR_WT_CURATION_FIRST"), "why"] = (
        "candidate needs paired TCR, WT, or structure curation before wetlab escalation"
    )
    out.loc[out["recommendation_tier"].eq("HOLD_FOR_NOW"), "why"] = "insufficient integrated evidence for priority wetlab"
    tier_order = {
        "TIER_A_EXPERIMENT_NOW_WITH_CONTROLS": 0,
        "TIER_A_MD_STRONG_CONTROL_CURATION_FIRST": 1,
        "TIER_A_PENDING_MD_COMPLETION_THEN_EXPERIMENT": 2,
        "TIER_B_EXPERIMENT_AFTER_CONTROL_MD": 3,
        "TIER_C_MD_OR_STRUCTURE_SCREEN_FIRST": 4,
        "TIER_D_TCR_OR_WT_CURATION_FIRST": 5,
        "HOLD_FOR_NOW": 6,
    }
    out["tier_order"] = out["recommendation_tier"].map(tier_order).fillna(9).astype(int)
    return out.sort_values(["tier_order", "ultra_priority_score"], ascending=[True, False])


def assay_plan(prioritized: pd.DataFrame) -> pd.DataFrame:
    rows = []
    top = prioritized[prioritized["recommendation_tier"].isin([
        "TIER_A_EXPERIMENT_NOW_WITH_CONTROLS",
        "TIER_A_MD_STRONG_CONTROL_CURATION_FIRST",
        "TIER_A_PENDING_MD_COMPLETION_THEN_EXPERIMENT",
        "TIER_B_EXPERIMENT_AFTER_CONTROL_MD",
    ])].head(12)
    for rank, (_, r) in enumerate(top.iterrows(), start=1):
        peptide = r["peptide"]
        hla = r["hla_4digit"]
        tier = r["recommendation_tier"]
        wt = r.get("tentative_wt_sequence", "")
        decoy = r.get("anchor_preserved_decoy", "")
        rows.extend(
            [
                {
                    "assay_rank": rank,
                    "peptide": peptide,
                    "hla_4digit": hla,
                    "recommendation_tier": tier,
                    "assay": "peptide_synthesis_and_qc",
                    "sample": peptide,
                    "required_controls": f"WT={wt}; decoy={decoy}; irrelevant same-HLA peptide",
                    "positive_call": "peptide passes purity/QC; no immunogenicity claim",
                },
                {
                    "assay_rank": rank,
                    "peptide": peptide,
                    "hla_4digit": hla,
                    "recommendation_tier": tier,
                    "assay": "HLA_binding_or_stability",
                    "sample": f"{peptide}/{hla}",
                    "required_controls": "WT, decoy, known same-HLA positive, no-peptide negative",
                    "positive_call": "mutant shows stronger or durable HLA binding/stability than decoy and preferably WT",
                },
                {
                    "assay_rank": rank,
                    "peptide": peptide,
                    "hla_4digit": hla,
                    "recommendation_tier": tier,
                    "assay": "pMHC_multimer_or_TCR_binding",
                    "sample": "paired TCR clone or TCR-transduced reporter if available",
                    "required_controls": "WT pMHC, decoy pMHC, same-HLA positive-control pMHC, irrelevant TCR",
                    "positive_call": "specific mutant pMHC binding above WT/decoy background",
                },
                {
                    "assay_rank": rank,
                    "peptide": peptide,
                    "hla_4digit": hla,
                    "recommendation_tier": tier,
                    "assay": "T_cell_activation",
                    "sample": "TCR reporter or donor T cells",
                    "required_controls": "WT, decoy, no-peptide, HLA-mismatched target, positive-control peptide",
                    "positive_call": "NFAT/CD69/cytokine activation specific to mutant peptide-HLA",
                },
                {
                    "assay_rank": rank,
                    "peptide": peptide,
                    "hla_4digit": hla,
                    "recommendation_tier": tier,
                    "assay": "tumor_processing_validation",
                    "sample": "minigene or tumor cell line/organoid if available",
                    "required_controls": "HLA restriction/blocking, antigen-negative target, WT construct",
                    "positive_call": "natural processing/presentation plus TCR activation/killing",
                },
            ]
        )
    return pd.DataFrame(rows)


def write_reports(prioritized: pd.DataFrame, assays: pd.DataFrame) -> None:
    ULTRA.mkdir(parents=True, exist_ok=True)
    cols = [
        "row_id",
        "peptide",
        "hla_4digit",
        "recommendation_tier",
        "ultra_priority_score",
        "confidence_score",
        "model_evidence_score",
        "tcr_resource_score",
        "md_structural_score",
        "control_readiness_score",
        "claim_risk_penalty",
        "why",
        "pmhc_score_mean",
        "tcr_augmented_score_mean",
        "paired_tcr_evidence_count",
        "tcr_evidence_count",
        "structure_evidence_count",
        "md_label",
        "md_score",
        "live_completion_fraction",
        "tentative_wt_sequence",
        "wt_status",
        "anchor_preserved_decoy",
        "tcr_link_claim_status",
        "public_overlap_flags",
    ]
    prioritized[[c for c in cols if c in prioritized.columns]].to_csv(ULTRA / "ultra_wetlab_priority_candidates.tsv", sep="\t", index=False)
    assays.to_csv(ULTRA / "ultra_wetlab_assay_plan.tsv", sep="\t", index=False)
    prioritized.to_csv(ULTRA / "ultra_wetlab_priority_full_evidence.tsv", sep="\t", index=False)

    top = prioritized[[c for c in cols if c in prioritized.columns]].head(20)
    tier_counts = prioritized["recommendation_tier"].value_counts().rename_axis("tier").reset_index(name="n")
    report = [
        "# Ultra-Conservative Wetlab Prioritization",
        "",
        "## Boundary",
        "",
        "This is a recommendation engine, not an immunogenicity oracle. Quantum/classical/MD computation can rank structural plausibility and uncertainty, but immune activation still requires wetlab validation.",
        "",
        "## Top Decision",
        "",
        "- `GADGVGKSAL / HLA-C*08:02`: experiment-priority candidate now, with controls. It has completed moderate MD support and TCR evidence, but WT/decoy assays remain required.",
        "- `HMTEVVRHC / HLA-A*02:01`: completed strong MD support and strong paired TCR evidence. The blocker is no longer MD completion; it is WT/decoy/control curation before experiment-now status.",
        "",
        "## Tier Counts",
        "",
        tier_counts.to_markdown(index=False),
        "",
        "## Top Candidates",
        "",
        top.to_markdown(index=False),
        "",
        "## Quantum Computing Position",
        "",
        "Quantum computing should not be presented as solving immunogenicity. If used at all, treat it as an optional optimization/kernel/uncertainty-exploration module inside model selection. It does not replace paired TCR data, WT/decoy controls, MD replicate consistency, or wetlab validation.",
        "",
        "## Accuracy Strategy",
        "",
        "1. Rank candidates using independent evidence blocks: pMHC model, TCR evidence, MD structure stability, and control readiness.",
        "2. Penalize claim risk: missing WT, missing decoy, nondefinitive TCR linkage, public-overlap risk, or pathogen-only evidence.",
        "3. Run wetlab in staged gates: HLA stability, TCR binding, activation, killing, natural processing.",
        "4. Promote only candidates passing mutant > WT/decoy specificity and replicate consistency.",
    ]
    (ULTRA / "ultra_wetlab_priority_report.md").write_text("\n".join(report) + "\n")


def main() -> None:
    base = load_candidate_base()
    md = md_candidate_evidence()
    controls = control_readiness()
    df = base.merge(md, on="pmhc_key", how="left").merge(controls, on="pmhc_key", how="left")
    prioritized = score_candidates(df)
    assays = assay_plan(prioritized)
    write_reports(prioritized, assays)
    summary = {
        "n_candidates": int(len(prioritized)),
        "top_candidate": prioritized[["peptide", "hla_4digit", "recommendation_tier", "ultra_priority_score"]].head(1).to_dict("records"),
        "tier_counts": prioritized["recommendation_tier"].value_counts().to_dict(),
    }
    (ULTRA / "ultra_wetlab_priority_summary.json").write_text(json.dumps(summary, indent=2))
    print("[ultra-priority]", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
