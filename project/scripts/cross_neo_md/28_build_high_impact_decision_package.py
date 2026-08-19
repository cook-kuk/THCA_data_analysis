#!/usr/bin/env python3
"""High-impact decision package for CROSS-Neo DL/TCR/structure/MD.

This package goes beyond "possible" and creates decision-ready artifacts:
- no-false-positive top candidate list from optimized thresholds
- source-stratified performance audit
- wetlab validation plan
- claim ladder and reviewer risk register
- high-impact summary figures
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[3]
SCRIPT_DIR = REPO / "project/scripts/cross_neo_md"
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
OUT = MD_OUT / "high_impact_decision_package"
FIG = MD_OUT / "figures"


def load_mod(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def numeric(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def load_frame_and_presets() -> tuple[pd.DataFrame, list[dict]]:
    opt = load_mod(SCRIPT_DIR / "25_optimize_dl_funnel_thresholds.py", "threshold_opt")
    df = opt.prepare_frame()
    preset_path = MD_OUT / "dl_threshold_optimization/threshold_recommended_presets.json"
    presets = json.loads(preset_path.read_text()) if preset_path.exists() else []
    return df, presets


def pred_mask(df: pd.DataFrame, cfg: dict, call_rule: str) -> pd.Series:
    pep_len = numeric(df, "peptide_length")
    hla = df["hla_4digit"].fillna("").astype(str).str.upper()
    is_class_i = hla.str.match(r"^HLA-[ABC]\*")
    supp = pd.Series(
        np.where(
            is_class_i,
            (pep_len >= cfg["class_i_min_len"]) & (pep_len <= cfg["class_i_max_len"]),
            (pep_len >= cfg["other_min_len"]) & (pep_len <= cfg["other_max_len"]),
        ),
        index=df.index,
    )
    dropout_unstable = (numeric(df, "dropout_sensitivity") >= cfg["dropout_sens"]) | (numeric(df, "perturb_width_90") >= cfg["perturb_width"])
    main_dl = numeric(df, "main_dl_score")
    stage1 = supp & ((main_dl >= cfg["main_dl"]) | (numeric(df, "ensemble_q95") >= cfg["ensemble_upper"]) | (numeric(df, "pmhc_score_mean") >= cfg["pmhc"]))
    stage2 = stage1 & (
        ((numeric(df, "bayes_mean") >= cfg["bayes_mean"]) & (numeric(df, "perturb_prob_gt_050") >= cfg["perturb_prob"]) & (~dropout_unstable))
        | (numeric(df, "bayes_q95") >= cfg["bayes_upper"])
        | (main_dl >= cfg["robust_main"])
    )
    stage3 = stage2 & (~dropout_unstable) & (
        (numeric(df, "bayes_mean") >= cfg["robust_bayes"])
        | (numeric(df, "perturb_prob_gt_050") >= cfg["robust_perturb"])
        | (main_dl >= cfg["robust_main"])
    )
    source_tcr = (
        (numeric(df, "tcr_evidence_count") >= cfg["min_tcr_evidence"])
        & (numeric(df, "cancer_context_evidence_count") >= cfg["min_cancer_context"])
        & (numeric(df, "pathogen_context_evidence_count") <= cfg["max_pathogen_context"])
    )
    paired_tcr = source_tcr & (numeric(df, "paired_tcr_evidence_count") >= cfg["min_paired_tcr"])
    tcr_support = (numeric(df, "tcr_augmented_score_mean") >= cfg["tcr_branch"]) | (numeric(df, "best_tcr_augmented_score") >= cfg["tcr_rescue"])
    tcr_rescue = (
        supp
        & paired_tcr
        & (numeric(df, "tcr_augmented_score_mean") >= cfg["tcr_rescue"])
        & (numeric(df, "best_tcr_augmented_score") >= cfg["tcr_rescue"])
        & (numeric(df, "bayes_q95") >= cfg["tcr_rescue_bayes_upper"])
    )
    baker_pass = numeric(df, "baker_structural_score") >= cfg["min_baker_structural"]
    md_score = pd.concat([numeric(df, "md_structural_score"), numeric(df, "md_score")], axis=1).max(axis=1)
    md_pass = (md_score >= cfg["min_md_structural"]) | (numeric(df, "live_completion_fraction") >= cfg["min_live_completion"])
    primary_paired = stage3 & paired_tcr & tcr_support
    tcr_context = stage3 & source_tcr & tcr_support
    discordant = stage3 & source_tcr & (~tcr_support)
    cull = (~supp) | (~stage1) | (stage1 & (~stage2) & (~tcr_rescue)) | (dropout_unstable & (~paired_tcr))
    wt_ready = df.get("wt_or_decoy_ready", pd.Series(False, index=df.index)).fillna(False).astype(bool)

    if call_rule == "md_escalation":
        return primary_paired | tcr_rescue
    if call_rule == "wetlab_shortlist":
        return primary_paired & wt_ready
    if call_rule == "structure_md_supported":
        return stage3 & source_tcr & tcr_support & baker_pass & md_pass
    if call_rule == "tcr_supported":
        return tcr_context
    if call_rule == "robust_dl":
        return stage3
    if call_rule == "non_culled":
        return (~cull) & (~discordant)
    return primary_paired | tcr_rescue


def preset_candidate_lists(df: pd.DataFrame, presets: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    perf_rows = []
    label = numeric(df, "label_binary").astype(int)
    for p in presets:
        mask = pred_mask(df, p["thresholds"], p["call_rule"])
        sub = df[mask].copy()
        sub["preset_name"] = p["name"]
        sub["call_rule"] = p["call_rule"]
        sub["actual_label"] = sub["label_binary"].map({1: "positive", 0: "negative"}).fillna("unknown")
        sub["prediction_outcome"] = np.where(sub["label_binary"].astype(int) == 1, "TP", "FP")
        sub["priority_order"] = (
            numeric(sub, "tcr_augmented_score_mean")
            + numeric(sub, "main_dl_score")
            + numeric(sub, "bayes_mean")
            + numeric(sub, "baker_structural_score")
            + numeric(sub, "md_structural_score")
        ).rank(ascending=False, method="first")
        rows.append(sub)
        for source, idx in df.groupby("source_dataset").groups.items():
            sm = mask.loc[idx]
            lab = label.loc[idx]
            tp = int((sm & (lab == 1)).sum())
            fp = int((sm & (lab == 0)).sum())
            tn = int(((~sm) & (lab == 0)).sum())
            fn = int(((~sm) & (lab == 1)).sum())
            perf_rows.append(
                {
                    "preset_name": p["name"],
                    "call_rule": p["call_rule"],
                    "source_dataset": source,
                    "n": int(len(idx)),
                    "positive": int((lab == 1).sum()),
                    "negative": int((lab == 0).sum()),
                    "called_positive": int(sm.sum()),
                    "TP": tp,
                    "TN": tn,
                    "FP": fp,
                    "FN": fn,
                    "precision": tp / (tp + fp) if tp + fp else 0.0,
                    "recall": tp / (tp + fn) if tp + fn else 0.0,
                    "FPR": fp / (fp + tn) if fp + tn else 0.0,
                }
            )
    all_candidates = pd.concat(rows, ignore_index=True, sort=False) if rows else pd.DataFrame()
    return all_candidates, pd.DataFrame(perf_rows)


def build_wetlab_plan(top13: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for i, r in top13.sort_values("priority_order").reset_index(drop=True).iterrows():
        peptide = r.get("peptide", "")
        hla = r.get("hla_4digit", "")
        has_wt = str(r.get("wt_status", "")).lower() not in {"", "nan", "missing_wt_sequence_blocked"}
        has_decoy = str(r.get("anchor_preserved_decoy", "")).lower() not in {"", "nan"}
        tcr_count = int(float(r.get("paired_tcr_evidence_count", 0) or 0))
        tier = "TIER_1_IMMEDIATE" if i < 4 else ("TIER_2_BACKUP" if i < 10 else "TIER_3_RESERVE")
        rows.append(
            {
                "plate_order": i + 1,
                "tier": tier,
                "row_id": r.get("row_id", ""),
                "peptide": peptide,
                "hla_4digit": hla,
                "source_dataset": r.get("source_dataset", ""),
                "label_binary": r.get("label_binary", ""),
                "recommended_assays": "HLA binding/stability; pMHC multimer if TCR available; IFN-gamma ELISpot/ICS; killing assay only after positive screen",
                "required_controls": "WT peptide" + (" available/tentative" if has_wt else " must be curated") + "; anchor-preserved decoy" + (" available" if has_decoy else " must be generated") + "; irrelevant peptide; known HLA-matched positive control",
                "sample_requirement": "HLA-matched donor/PBMC or patient material; paired TCR clone if available",
                "go_no_go": "Promote only if peptide-HLA binding plus T-cell activation exceeds WT/decoy controls",
                "why_this_candidate": f"optimized no-FP preset; paired TCR evidence={tcr_count}; mainDL={float(r.get('main_dl_score', 0) or 0):.3f}; TCRbranch={float(r.get('tcr_augmented_score_mean', 0) or 0):.3f}; MD={r.get('md_label', 'NA')}",
            }
        )
    return pd.DataFrame(rows)


def refresh_current_decision_annotations(top13: pd.DataFrame) -> pd.DataFrame:
    """Replace stale pending-MD text after new trajectory evidence is merged."""
    if top13.empty:
        return top13
    top13 = top13.copy()
    md_support = top13["md_label"].astype(str).isin(["MD_MODERATE", "MD_STRONG", "MD_VERY_STRONG"])
    paired = pd.to_numeric(top13.get("paired_tcr_evidence_count", 0), errors="coerce").fillna(0) > 0
    top_rank = pd.to_numeric(top13.get("priority_order", 999), errors="coerce").fillna(999) <= 4
    wt_missing = top13.get("wt_status", pd.Series("", index=top13.index)).astype(str).str.contains("missing|blocked|nan", case=False, regex=True)
    supported = md_support & paired & top_rank
    top13.loc[supported & wt_missing, "recommendation_tier"] = "TIER_A_EXPERIMENT_NOW_WITH_WT_CURATION"
    top13.loc[supported & (~wt_missing), "recommendation_tier"] = "TIER_A_EXPERIMENT_NOW_WITH_CONTROLS"
    top13.loc[supported, "why"] = (
        "completed MD structural support plus paired TCR evidence; run WT/decoy-controlled wetlab assays before immunogenicity claim"
    )
    top13.loc[supported & wt_missing, "next_action"] = "curate_wt_decoy_then_specificity_md_and_wetlab"
    top13.loc[supported & (~wt_missing), "next_action"] = "run_wt_decoy_specificity_md_and_wetlab"
    return top13


def build_claim_ladder() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "claim_level": "L1_ALLOWED_NOW",
                "claim": "DL-first/TCR-aware triage enriches current labeled candidates and narrows expensive testing.",
                "supporting_artifacts": "stage_label_counts; threshold presets; slider dashboard",
                "required_before_stronger_claim": "external validation on heldout source/study",
                "risk": "optimized on current labels",
            },
            {
                "claim_level": "L2_ALLOWED_FOR_CASE_STUDY",
                "claim": "GADGVGKSAL and HMTEVVRHC are structurally auditable high-priority candidates.",
                "supporting_artifacts": "candidate board; Baker fallback; MD status/evidence",
                "required_before_stronger_claim": "WT/decoy controls and replicate MD/wetlab",
                "risk": "structure/MD do not prove recognition",
            },
            {
                "claim_level": "L3_ALLOWED_AFTER_WETLAB",
                "claim": "A candidate is immunogenic in the tested assay context.",
                "supporting_artifacts": "ELISpot/ICS/killing assay with controls",
                "required_before_stronger_claim": "orthogonal assay and source-independent replication",
                "risk": "assay-specific positivity may not generalize",
            },
            {
                "claim_level": "L4_FORBIDDEN_NOW",
                "claim": "Universal TCR-aware neoantigen prediction or clinical utility.",
                "supporting_artifacts": "none sufficient yet",
                "required_before_stronger_claim": "large paired TCR datasets, prospective validation, clinical endpoints",
                "risk": "overclaim",
            },
        ]
    )


def build_external_validation_plan() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "validation_axis": "source-heldout",
                "dataset_needed": "independent cancer neoantigen source not used in threshold optimization",
                "metric": "precision@k, AUPRC, FP rate",
                "success_criterion": "no-FP/high-precision preset remains enriched over prevalence",
                "owner": "model validation",
            },
            {
                "validation_axis": "TCR-heldout",
                "dataset_needed": "paired alpha/beta TCR-pMHC rows with non-overlapping TCRs",
                "metric": "TCR branch concordance, rescue/harm rate",
                "success_criterion": "TCR branch improves precision without hidden TCR leakage",
                "owner": "TCR expert track",
            },
            {
                "validation_axis": "WT/decoy specificity",
                "dataset_needed": "WT peptide and anchor-preserved decoy for top candidates",
                "metric": "mutant > WT/decoy activation and interface stability",
                "success_criterion": "mutant-specific response under matched controls",
                "owner": "wetlab/structure",
            },
            {
                "validation_axis": "MD/structure reproducibility",
                "dataset_needed": "3x10 ns per mutant/WT/decoy, selected 50-100 ns follow-up",
                "metric": "contact persistence, RMSD, interface delta",
                "success_criterion": "consistent structural plausibility and no WT cross-reactivity warning",
                "owner": "MD audit",
            },
        ]
    )


def build_reviewer_risk_register() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "risk": "threshold overfit",
                "why_reviewer_will_ask": "optimizer used current labels",
                "defense": "present as decision-support preset, not final clinical threshold",
                "next_action": "run source-heldout/external validation",
            },
            {
                "risk": "label leakage",
                "why_reviewer_will_ask": "TCR and peptide-HLA rows can overlap public resources",
                "defense": "report public overlap and strict split constraints",
                "next_action": "deduplicate peptide-HLA/TCR/CDR3 clusters before benchmark claim",
            },
            {
                "risk": "MD overclaim",
                "why_reviewer_will_ask": "short simulations do not prove binding or immunogenicity",
                "defense": "MD is only a structural audit layer",
                "next_action": "add WT/decoy/replicate controls",
            },
            {
                "risk": "TCR missingness",
                "why_reviewer_will_ask": "most neoantigen data lack paired TCRs",
                "defense": "TCR branch is optional expert and diagnostic branch",
                "next_action": "separate main pMHC model from TCR-aware manuscript track",
            },
        ]
    )


def build_evidence_to_claim_matrix(top13: pd.DataFrame, perf: pd.DataFrame, wetlab: pd.DataFrame) -> pd.DataFrame:
    top_n = int(len(top13))
    tp = int((top13.get("prediction_outcome", pd.Series(dtype=str)).astype(str).eq("TP")).sum())
    fp = int((top13.get("prediction_outcome", pd.Series(dtype=str)).astype(str).eq("FP")).sum())
    paired = int((pd.to_numeric(top13.get("paired_tcr_evidence_count", pd.Series(dtype=float)), errors="coerce").fillna(0) > 0).sum())
    md_supported = int(top13.get("md_label", pd.Series(dtype=str)).astype(str).str.contains("MODERATE|STRONG|VERY_STRONG", regex=True, na=False).sum())
    hm = top13[top13.get("peptide", pd.Series(dtype=str)).astype(str).eq("HMTEVVRHC")]
    gad = top13[top13.get("peptide", pd.Series(dtype=str)).astype(str).eq("GADGVGKSAL")]
    source_rows = perf[perf.get("preset_name", pd.Series(dtype=str)).astype(str).eq("NO_FALSE_POSITIVE_MAX_TP")]
    fp_by_source = int(pd.to_numeric(source_rows.get("FP", pd.Series(dtype=float)), errors="coerce").fillna(0).sum()) if not source_rows.empty else fp
    rows = [
        {
            "evidence_layer": "strict current-label operating point",
            "actual_data": f"Top-{top_n}; TP={tp}; FP={fp}",
            "supports_claim": "decision-support triage can create a compact wetlab board",
            "allowed_strength_now": "high for internal/current-label prioritization",
            "blocked_overclaim": "not external SOTA, not clinical utility, not universal precision",
            "figure_or_table": "fig_md27; no_false_positive_top13_candidates.tsv",
        },
        {
            "evidence_layer": "source-stratified current-label audit",
            "actual_data": f"{len(source_rows)} source rows; aggregate FP={fp_by_source}",
            "supports_claim": "the strict preset is not explained by one visible false-positive-heavy source in the current table",
            "allowed_strength_now": "reviewer-defense audit",
            "blocked_overclaim": "not a source-heldout external validation",
            "figure_or_table": "fig_md28; preset_source_stratified_performance.tsv",
        },
        {
            "evidence_layer": "paired TCR-resource evidence",
            "actual_data": f"{paired}/{top_n} strict Top candidates have paired TCR evidence",
            "supports_claim": "recognition-aware ranking is inspectable rather than peptide-only",
            "allowed_strength_now": "mechanistic plausibility / prioritization",
            "blocked_overclaim": "not direct T-cell activation proof",
            "figure_or_table": "fig_md15; fig_md26; no_false_positive_top13_candidates.tsv",
        },
        {
            "evidence_layer": "OpenMM structural audit",
            "actual_data": f"{md_supported}/{top_n} strict Top candidates have MD_MODERATE-or-better labels",
            "supports_claim": "completed MD supports structural plausibility for the top case studies",
            "allowed_strength_now": "case-study structural support",
            "blocked_overclaim": "not immunogenicity proof and not a substitute for WT/decoy controls",
            "figure_or_table": "fig_md12; fig_md13; md_evidence_scores.tsv",
        },
        {
            "evidence_layer": "GADGVGKSAL case",
            "actual_data": "absent" if gad.empty else f"mainDL={float(gad.iloc[0].get('main_dl_score', 0)):.3f}; TCR={float(gad.iloc[0].get('tcr_augmented_score_mean', 0)):.3f}; MD={gad.iloc[0].get('md_label', 'NA')}",
            "supports_claim": "balanced DL+TCR+MD candidate ready for controlled validation",
            "allowed_strength_now": "lead wetlab candidate",
            "blocked_overclaim": "needs WT/decoy assay readout before immunogenicity claim",
            "figure_or_table": "fig_md26; fig_md27; wetlab_validation_plate_plan.tsv",
        },
        {
            "evidence_layer": "HMTEVVRHC case",
            "actual_data": "absent" if hm.empty else f"mainDL={float(hm.iloc[0].get('main_dl_score', 0)):.3f}; TCR={float(hm.iloc[0].get('tcr_augmented_score_mean', 0)):.3f}; MD={hm.iloc[0].get('md_label', 'NA')}",
            "supports_claim": "TCR+MD rescue case despite weaker main DL",
            "allowed_strength_now": "high-value structural case study",
            "blocked_overclaim": "WT sequence curation remains a blocker before specificity claims",
            "figure_or_table": "fig_md26; fig_md27; md_evidence_score_card_HMTEVVRHC_HLA-A0201.md",
        },
        {
            "evidence_layer": "wetlab control plan",
            "actual_data": f"{len(wetlab)} plate-plan rows with WT/decoy/go-no-go columns",
            "supports_claim": "model output is immediately translatable into controlled assays",
            "allowed_strength_now": "experimental design readiness",
            "blocked_overclaim": "the assays have not been run yet",
            "figure_or_table": "fig_md30; wetlab_validation_plate_plan.tsv",
        },
    ]
    return pd.DataFrame(rows)


def build_manuscript_figure_plan() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "slot": "Main Figure A",
                "asset": "fig_md24_story_overview_flow",
                "purpose": "end-to-end candidate funnel overview",
                "data_used": "stage_label_composition_for_slides.tsv; threshold_preset_summary_for_slides.tsv",
                "claim_boundary": "prioritization workflow, not immunogenicity proof",
            },
            {
                "slot": "Main Figure B",
                "asset": "fig_md27_no_fp_top13_candidate_evidence",
                "purpose": "strict operating point candidate evidence board",
                "data_used": "no_false_positive_top13_candidates.tsv",
                "claim_boundary": "current-label result; external validation required",
            },
            {
                "slot": "Main Figure C",
                "asset": "fig_md26_representative_candidate_board",
                "purpose": "GADGVGKSAL vs HMTEVVRHC case comparison",
                "data_used": "representative_candidate_evidence_table.tsv; md_evidence_scores.tsv",
                "claim_boundary": "case-study structural audit",
            },
            {
                "slot": "Main Figure D",
                "asset": "fig_md46_evidence_to_claim_matrix",
                "purpose": "what evidence supports which claim and what remains blocked",
                "data_used": "evidence_to_claim_matrix.tsv",
                "claim_boundary": "prevents overclaiming",
            },
            {
                "slot": "Supplementary Figure",
                "asset": "fig_md12_evidence_score_components; fig_md13_visual_summary_dashboard",
                "purpose": "MD score components and run-level QC dashboard",
                "data_used": "md_evidence_scores.tsv; md_status_summary.tsv; qc/*.tsv; contacts/*.tsv",
                "claim_boundary": "structural plausibility only",
            },
            {
                "slot": "Supplementary Table",
                "asset": "data_usage_manifest.tsv; high_impact_figure_table_caption_guide.tsv",
                "purpose": "data provenance and figure/table explanation",
                "data_used": "decision_storyboard/*.tsv; high_impact_decision_package/*.tsv",
                "claim_boundary": "documentation/reproducibility layer",
            },
        ]
    )


def plot_evidence_to_claim_matrix(matrix: pd.DataFrame) -> None:
    if matrix.empty:
        return
    strength_map = {
        "high for internal/current-label prioritization": 3,
        "reviewer-defense audit": 2,
        "mechanistic plausibility / prioritization": 2,
        "case-study structural support": 2,
        "lead wetlab candidate": 3,
        "high-value structural case study": 3,
        "experimental design readiness": 2,
    }
    scores = matrix["allowed_strength_now"].map(strength_map).fillna(1).to_numpy()[:, None]
    fig, ax = plt.subplots(figsize=(12.2, 5.8))
    ax.imshow(scores, aspect="auto", cmap="YlGnBu", vmin=0, vmax=3)
    ax.set_xticks([0])
    ax.set_xticklabels(["claim support\nstrength"])
    ax.set_yticks(range(len(matrix)))
    ax.set_yticklabels(matrix["evidence_layer"], fontsize=9)
    for i, row in matrix.reset_index(drop=True).iterrows():
        ax.text(0, i, str(scores[i, 0]), ha="center", va="center", color="#07111f", weight="bold")
        ax.text(0.64, i, row["actual_data"], ha="left", va="center", fontsize=8, transform=ax.get_yaxis_transform())
        ax.text(1.03, i, row["blocked_overclaim"], ha="left", va="center", fontsize=8, color="#9b1c31", transform=ax.get_yaxis_transform())
    ax.set_title("Evidence-to-claim matrix: impact with explicit claim boundaries")
    ax.text(0.64, 1.04, "actual data", transform=ax.transAxes, fontsize=9, weight="bold")
    ax.text(1.03, 1.04, "blocked overclaim", transform=ax.transAxes, fontsize=9, weight="bold", color="#9b1c31")
    fig.savefig(FIG / "fig_md46_evidence_to_claim_matrix.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md46_evidence_to_claim_matrix.pdf", bbox_inches="tight")
    plt.close(fig)


def write_kakao_summary(top13: pd.DataFrame, matrix: pd.DataFrame) -> None:
    gad = top13[top13.get("peptide", pd.Series(dtype=str)).astype(str).eq("GADGVGKSAL")]
    hm = top13[top13.get("peptide", pd.Series(dtype=str)).astype(str).eq("HMTEVVRHC")]
    gad_md = gad.iloc[0].get("md_label", "NA") if not gad.empty else "NA"
    hm_md = hm.iloc[0].get("md_label", "NA") if not hm.empty else "NA"
    msg = (
        "CROSS-Neo cancer vaccine 쪽은 cheap DL로 649개 후보를 먼저 줄이고, Bayesian/dropout uncertainty, paired TCR evidence, "
        "Baker/static structure, OpenMM MD까지 단계적으로 붙인 decision-support 파이프라인입니다. 현재 label 기준 strict preset은 "
        f"Top-{len(top13)}에서 TP {int((top13.get('prediction_outcome', pd.Series(dtype=str)).astype(str).eq('TP')).sum())} / "
        f"FP {int((top13.get('prediction_outcome', pd.Series(dtype=str)).astype(str).eq('FP')).sum())}이고, "
        f"GADGVGKSAL/HLA-C*08:02는 DL+TCR+MD가 균형 잡힌 lead({gad_md}), HMTEVVRHC/HLA-A*02:01은 main DL은 약하지만 "
        f"TCR+10ns MD가 강한 rescue case({hm_md})입니다. 결론은 '면역원성 증명'이 아니라 WT/decoy-controlled wetlab으로 바로 넘길 "
        "하이임팩트 후보선별/구조감사 패키지입니다."
    )
    (OUT / "KAKAO_ONE_SHOT_HIGH_IMPACT_KR.md").write_text(msg + "\n")


def plot_no_fp_candidates(top13: pd.DataFrame) -> None:
    if top13.empty:
        return
    show = top13.sort_values("priority_order").head(13).copy()
    labels = show["peptide"].astype(str) + "\n" + show["hla_4digit"].astype(str)
    metrics = [
        ("main_dl_score", "main DL"),
        ("bayes_mean", "Bayes"),
        ("tcr_augmented_score_mean", "TCR branch"),
        ("perturb_prob_gt_050", "perturb"),
        ("baker_structural_score", "Baker"),
        ("md_structural_score", "MD"),
    ]
    mat = []
    for col, _ in metrics:
        mat.append(pd.to_numeric(show.get(col, 0), errors="coerce").fillna(0).clip(0, 1).to_numpy())
    mat = np.vstack(mat).T
    fig, ax = plt.subplots(figsize=(10.8, 6.5))
    im = ax.imshow(mat, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_yticks(range(len(show)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([m[1] for m in metrics], rotation=25, ha="right")
    ax.set_title("No-false-positive optimized candidate panel")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    fig.savefig(FIG / "fig_md27_no_fp_top13_candidate_evidence.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md27_no_fp_top13_candidate_evidence.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_source_stratified(perf: pd.DataFrame) -> None:
    if perf.empty:
        return
    sub = perf[perf["preset_name"].isin(["NO_FALSE_POSITIVE_MAX_TP", "BALANCED_F1", "RECALL_PRESERVING"])].copy()
    if sub.empty:
        sub = perf.copy()
    sources = sorted(sub["source_dataset"].astype(str).unique())
    presets = list(sub["preset_name"].drop_duplicates())
    x = np.arange(len(sources))
    width = 0.8 / max(len(presets), 1)
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    for i, preset in enumerate(presets):
        vals = []
        for source in sources:
            row = sub[(sub["preset_name"] == preset) & (sub["source_dataset"].astype(str) == source)]
            vals.append(float(row["precision"].iloc[0]) if not row.empty else 0.0)
        ax.bar(x + (i - len(presets) / 2) * width + width / 2, vals, width=width, label=preset.replace("_", "\n"))
    ax.set_xticks(x)
    ax.set_xticklabels(sources, rotation=20, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Precision")
    ax.set_title("Source-stratified precision by threshold preset")
    ax.legend(frameon=False, fontsize=7)
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(FIG / "fig_md28_source_stratified_preset_performance.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md28_source_stratified_preset_performance.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_claim_ladder(claims: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    ax.axis("off")
    colors = ["#2fbf71", "#4d96ff", "#f2c46d", "#ff7b72"]
    for i, r in claims.iterrows():
        y = 0.82 - i * 0.22
        ax.add_patch(plt.Rectangle((0.05, y - 0.075), 0.90, 0.15, color=colors[i], alpha=0.9, ec="#111"))
        ax.text(0.07, y + 0.025, r["claim_level"], weight="bold", fontsize=11, color="#07111f")
        ax.text(0.07, y - 0.025, r["claim"], fontsize=10, color="#07111f")
        ax.text(0.62, y - 0.025, f"Risk: {r['risk']}", fontsize=9, color="#07111f")
    ax.set_title("Claim ladder: what can be said now versus later", fontsize=16, weight="bold")
    fig.savefig(FIG / "fig_md29_claim_ladder.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md29_claim_ladder.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_wetlab_workflow(plan: pd.DataFrame) -> None:
    steps = ["Peptide\nsynthesis", "HLA binding\n/stability", "pMHC multimer\nif TCR", "ELISpot/ICS", "WT/decoy\nspecificity", "kill assay\nfollow-up"]
    fig, ax = plt.subplots(figsize=(12, 4.6))
    ax.axis("off")
    x = np.linspace(0.08, 0.92, len(steps))
    for i, (xi, step) in enumerate(zip(x, steps)):
        ax.add_patch(plt.Circle((xi, 0.58), 0.065, color="#183a5a", ec="#edf5ff", lw=1.2))
        ax.text(xi, 0.58, str(i + 1), ha="center", va="center", color="#f2c46d", weight="bold", fontsize=16)
        ax.text(xi, 0.37, step, ha="center", va="center", fontsize=10)
        if i < len(steps) - 1:
            ax.annotate("", xy=(x[i + 1] - 0.075, 0.58), xytext=(xi + 0.075, 0.58), arrowprops=dict(arrowstyle="->", lw=2))
    ax.text(0.5, 0.88, f"Wetlab validation workflow for {len(plan)} optimized candidates", ha="center", fontsize=16, weight="bold")
    ax.text(0.5, 0.12, "Advance only if mutant response exceeds WT/decoy and assay controls", ha="center", fontsize=12, color="#9b1c31")
    fig.savefig(FIG / "fig_md30_wetlab_validation_workflow.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md30_wetlab_validation_workflow.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(
    top13: pd.DataFrame,
    perf: pd.DataFrame,
    wetlab: pd.DataFrame,
    claims: pd.DataFrame,
    validation: pd.DataFrame,
    risks: pd.DataFrame,
    evidence_matrix: pd.DataFrame,
    figure_plan: pd.DataFrame,
) -> None:
    lines = [
        "# High-Impact CROSS-Neo Decision Package",
        "",
        "## Decision Headline",
        "",
        "The strongest current operating point finds 13 labeled positives with 0 labeled false positives in the current candidate table. This is a decision-support result that must be externally validated before strong claims.",
        "",
        "## Top No-FP Candidates",
        "",
        top13[["row_id", "peptide", "hla_4digit", "source_dataset", "label_binary", "main_dl_score", "bayes_mean", "tcr_augmented_score_mean", "paired_tcr_evidence_count", "baker_structural_score", "md_label"]].to_markdown(index=False)
        if not top13.empty
        else "No candidates.",
        "",
        "## Evidence-To-Claim Matrix",
        "",
        evidence_matrix.to_markdown(index=False),
        "",
        "## Manuscript Figure/Table Plan",
        "",
        figure_plan.to_markdown(index=False),
        "",
        "## Source-Stratified Performance",
        "",
        perf[["preset_name", "source_dataset", "n", "positive", "negative", "called_positive", "TP", "FP", "precision", "recall", "FPR"]].to_markdown(index=False),
        "",
        "## Wetlab Plan",
        "",
        wetlab[["plate_order", "tier", "peptide", "hla_4digit", "recommended_assays", "required_controls", "go_no_go"]].to_markdown(index=False),
        "",
        "## Claim Ladder",
        "",
        claims.to_markdown(index=False),
        "",
        "## External Validation Plan",
        "",
        validation.to_markdown(index=False),
        "",
        "## Reviewer Risk Register",
        "",
        risks.to_markdown(index=False),
    ]
    (OUT / "HIGH_IMPACT_DECISION_REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    df, presets = load_frame_and_presets()
    all_candidates, source_perf = preset_candidate_lists(df, presets)
    top13 = all_candidates[all_candidates["preset_name"].eq("NO_FALSE_POSITIVE_MAX_TP")].copy()
    top13 = top13.sort_values("priority_order")
    top13 = refresh_current_decision_annotations(top13)
    wetlab = build_wetlab_plan(top13)
    claims = build_claim_ladder()
    validation = build_external_validation_plan()
    risks = build_reviewer_risk_register()
    evidence_matrix = build_evidence_to_claim_matrix(top13, source_perf, wetlab)
    figure_plan = build_manuscript_figure_plan()

    all_candidates.to_csv(OUT / "optimized_preset_candidate_lists.tsv", sep="\t", index=False)
    top13.to_csv(OUT / "no_false_positive_top13_candidates.tsv", sep="\t", index=False)
    source_perf.to_csv(OUT / "preset_source_stratified_performance.tsv", sep="\t", index=False)
    wetlab.to_csv(OUT / "wetlab_validation_plate_plan.tsv", sep="\t", index=False)
    claims.to_csv(OUT / "claim_ladder.tsv", sep="\t", index=False)
    validation.to_csv(OUT / "external_validation_plan.tsv", sep="\t", index=False)
    risks.to_csv(OUT / "reviewer_risk_register.tsv", sep="\t", index=False)
    evidence_matrix.to_csv(OUT / "evidence_to_claim_matrix.tsv", sep="\t", index=False)
    figure_plan.to_csv(OUT / "manuscript_figure_table_plan.tsv", sep="\t", index=False)

    plot_no_fp_candidates(top13)
    plot_source_stratified(source_perf)
    plot_claim_ladder(claims)
    plot_wetlab_workflow(wetlab)
    plot_evidence_to_claim_matrix(evidence_matrix)
    write_kakao_summary(top13, evidence_matrix)
    write_report(top13, source_perf, wetlab, claims, validation, risks, evidence_matrix, figure_plan)

    summary = {
        "n_top_no_fp": int(len(top13)),
        "top_no_fp_positive": int(pd.to_numeric(top13.get("label_binary", pd.Series(dtype=float)), errors="coerce").fillna(0).sum()) if not top13.empty else 0,
        "n_wetlab_plan": int(len(wetlab)),
        "figures": [
            "fig_md27_no_fp_top13_candidate_evidence",
            "fig_md28_source_stratified_preset_performance",
            "fig_md29_claim_ladder",
            "fig_md30_wetlab_validation_workflow",
            "fig_md46_evidence_to_claim_matrix",
        ],
    }
    (OUT / "high_impact_summary.json").write_text(json.dumps(summary, indent=2))
    print("[high-impact]", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
