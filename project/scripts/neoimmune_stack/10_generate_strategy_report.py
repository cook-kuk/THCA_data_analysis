#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from neoimmune_common import HUB_ROOT, REPO_ROOT, ensure_run_dir, json_load, safe_read_table, write_md


def top_table(df: pd.DataFrame, cols: list[str], n: int = 8) -> str:
    if df.empty:
        return "_No rows available._"
    use = [c for c in cols if c in df.columns]
    return df[use].head(n).to_markdown(index=False)


def answer_quantum(clean_lb: pd.DataFrame) -> str:
    q = clean_lb[clean_lb["model"].astype(str).str.contains("quantum|GP_quantum|W7A_QK", case=False, na=False)] if not clean_lb.empty else pd.DataFrame()
    if q.empty:
        return "Not proven in this integrated run. Quantum artifacts are registered, but no source-heldout survivorship claim should be made until their split-specific scores are present."
    best = q.sort_values("AUPRC", ascending=False).head(1).iloc[0]
    return f"Best integrated quantum-labeled row is `{best['model']}` with AUPRC={best.get('AUPRC')}; this is still artifact-level evidence unless the row's split is strict/source-heldout."


def answer_wave8(clean_lb: pd.DataFrame) -> str:
    w = clean_lb[clean_lb["model"].astype(str).str.contains("Wave8|TCR|SelfSim", case=False, na=False)] if not clean_lb.empty else pd.DataFrame()
    if w.empty:
        return "Not measurable yet in this integrated run; Wave8/TCR-self-similarity remains a required branch, but the collector did not find a directly named score table."
    best = w.sort_values("AUPRC", ascending=False).head(1).iloc[0]
    return f"Wave8/TCR-self-similarity evidence exists as `{best['model']}` with AUPRC={best.get('AUPRC')}; interpret only within its stated split."


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    audit = outdir / "audit" / "repo_audit.md"
    registry = safe_read_table(outdir / "registry" / "model_registry.tsv")
    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    clean_lb = safe_read_table(outdir / "metrics" / "leaderboard_clean_track.tsv") if (outdir / "metrics" / "leaderboard_clean_track.tsv").exists() else pd.DataFrame()
    prod_lb = safe_read_table(outdir / "metrics" / "leaderboard_production_track.tsv") if (outdir / "metrics" / "leaderboard_production_track.tsv").exists() else pd.DataFrame()
    ga_metrics = safe_read_table(outdir / "metrics" / "ga_rl_algorithm_metrics.tsv") if (outdir / "metrics" / "ga_rl_algorithm_metrics.tsv").exists() else pd.DataFrame()
    ga_top34 = safe_read_table(outdir / "predictions" / "patient_top34_candidates_kg_ga.tsv") if (outdir / "predictions" / "patient_top34_candidates_kg_ga.tsv").exists() else pd.DataFrame()
    imneo_metrics = safe_read_table(outdir / "metrics" / "imneo_public_reconstruction_leaderboard.tsv") if (outdir / "metrics" / "imneo_public_reconstruction_leaderboard.tsv").exists() else pd.DataFrame()
    imneo_sources = safe_read_table(outdir / "reports" / "imneo_public_source_inventory.tsv") if (outdir / "reports" / "imneo_public_source_inventory.tsv").exists() else pd.DataFrame()
    source_gap = safe_read_table(outdir / "metrics" / "source_transfer_gap_summary.tsv") if (outdir / "metrics" / "source_transfer_gap_summary.tsv").exists() else pd.DataFrame()
    rescue_summary = safe_read_table(outdir / "metrics" / "source_rescue_queue_summary.tsv") if (outdir / "metrics" / "source_rescue_queue_summary.tsv").exists() else pd.DataFrame()
    collab_pos = safe_read_table(outdir / "predictions" / "collaborator_positive_handoff.tsv") if (outdir / "predictions" / "collaborator_positive_handoff.tsv").exists() else pd.DataFrame()
    collab_ctrl = safe_read_table(outdir / "predictions" / "collaborator_control_handoff.tsv") if (outdir / "predictions" / "collaborator_control_handoff.tsv").exists() else pd.DataFrame()
    leak = safe_read_table(outdir / "metrics" / "leakage_summary.tsv") if (outdir / "metrics" / "leakage_summary.tsv").exists() else pd.DataFrame()
    top20 = safe_read_table(outdir / "predictions" / "patient_top20_candidates.tsv") if (outdir / "predictions" / "patient_top20_candidates.tsv").exists() else pd.DataFrame()
    hub_leak = json_load(HUB_ROOT / "logs/leakage_audit.json") or {}
    qx = json_load(HUB_ROOT / "logs/qx_analyses.json") or {}
    esm = json_load(HUB_ROOT / "logs/esm2_full_analyses.json") or {}

    usable_public = registry[(registry["source"].astype(str).str.contains("external", na=False)) & (registry["production_track_allowed"].astype(str).isin(["True", "true", "1"]))]
    local = registry[registry["source"].astype(str).eq("local")]
    best_clean = clean_lb.head(1).to_dict("records")
    best_prod = prod_lb.head(1).to_dict("records")
    real_patient_rows = 0
    real_patient_unique = 0
    if "patient_id" in canon.columns:
        real_mask = canon["patient_id"].fillna("unknown_patient").astype(str).ne("unknown_patient")
        real_patient_rows = int(real_mask.sum())
        real_patient_unique = int(canon.loc[real_mask, "patient_id"].astype(str).nunique())

    md = [
        "# NEOIMMUNE-STACK STRATEGY REPORT",
        "",
        "**System names:** NeoImmune-Stack for business framing; CLEAN-Neo++ for paper framing.",
        "",
        "**Core decision:** do not build a new foundation model yet. Integrate local algorithms, public frozen predictors, leakage-aware benchmarking, and patient-level top-N candidate ranking.",
        "",
        "## Brutal executive decision",
        "- **대박 가능성:** yes, as an integration/product strategy layer. It already unifies local algorithms, public comparators, leakage audit, and candidate handoff.",
        f"- **논문 claim:** upgraded. The integrated public tables now contain {real_patient_unique:,} real non-placeholder patient IDs ({real_patient_rows:,} rows), so a retrospective public patient-level ranking scaffold is claimable.",
        "- **Still not claimable:** prospective hospital utility, clinical vaccine efficacy, presentation without MS, or immunogenicity without T-cell assay.",
        "- **Most defensible paper angle now:** leakage-aware patient-level immunogenicity ranking using public cohorts, plus a clear bridge to hospital-grade validation.",
        "- **Must-have next data:** stronger per-patient expression/VAF/clonality, HLA typing QC, APM/HLA-LOH, and orthogonal presentation/immunogenicity assay labels.",
        "",
        "## Scientific framing",
        "`CLEAN-Neo++: leakage-aware integration of peptide-HLA presentation, TCR-visible immunogenicity, structure proxies, and quantum-kernel features for patient-level cancer vaccine candidate ranking.`",
        "",
        "## Business framing",
        "`NeoImmune-Stack: an AI operating layer for patient-specific cancer vaccine candidate prioritization.`",
        "",
        "## Latest LLM strategy",
        "- Use the latest available language model only as an **analysis infrastructure layer**: evidence summarization, candidate rationale drafting, contradiction detection, wet-lab question generation, and reviewer-response audit.",
        "- Do not use LLM text output as a clean-science training feature.",
        "- Do not let the LLM override leakage flags, labels, or wet-lab validation requirements.",
        "- Add model name/config at runtime through `NEOIMMUNE_LLM_MODEL`; keep disabled by default until API/provider policy is selected.",
        "",
        "## Latest GA/RL algorithm integration",
        "The just-discovered `KG_GA_evolved_controller` has been integrated as a **production/experiment-priority branch**. It is deliberately not admitted to the clean science track because the controller includes public-predictor and product-value components.",
        "",
        top_table(ga_metrics, ["algorithm", "n", "positives", "AUPRC", "AUROC", "Precision@34", "Recall@34", "patient_hit_rate@34", "patient_recall@34", "patients_evaluated", "clean_track_allowed"], 8),
        "",
        f"- GA/RL patient top-34 rows: {len(ga_top34):,}" if not ga_top34.empty else "- GA/RL patient top-34 rows: not generated.",
        f"- GA/RL patient coverage: {ga_top34['patient_id'].nunique():,} patients" if not ga_top34.empty and "patient_id" in ga_top34 else "- GA/RL patient coverage: not available.",
        f"- Mapping modes: `{ga_top34['kg_ga_mapping_mode'].value_counts().to_dict()}`" if not ga_top34.empty and "kg_ga_mapping_mode" in ga_top34 else "- Mapping modes: not available.",
        "- Claim-safe wording: GA/RL found a retrospective candidate-prioritization controller ready for prospective assay validation, not a clinically validated vaccine selection model.",
        "",
        "## Competitor-facing imNEO public-source reconstruction",
        "This is the narrow comparison requested against the public-source space visible from CG Invites/imNEO patent/presentation context. It is not the proprietary imNEO peptide list and not a reverse-engineered company model.",
        "",
        top_table(imneo_metrics[imneo_metrics["subset"].astype(str).eq("imNEO_core_head_to_head_common")] if not imneo_metrics.empty and "subset" in imneo_metrics else imneo_metrics, ["model", "track", "clean_track_allowed", "n", "positives", "real_patients", "patient_level_claimable", "AUPRC", "AUROC", "Precision@34", "Recall@34"], 14),
        "",
        "Coverage of disclosed public sources in the local canonical table:",
        top_table(imneo_sources, ["source_name", "status", "present_rows", "positive_labels", "nmi_locked_rows", "public_comparator_rows", "claim_boundary"], 10),
        "",
        "Brutal boundary: NMI beats frozen public comparators on the current NEPdb common head-to-head subset, but patient-level top-N is not claimable there because real patient IDs are absent/collapsed. TESLA and dbPepNeo need NMI locked-branch coverage before broad competitor claims.",
        "",
        "## Source transfer gap",
        top_table(source_gap, ["source_family", "rows", "positives", "nmi_best_model", "nmi_best_auprc", "public_best_model", "public_best_auprc", "delta_nmi_minus_public", "patient_claimable"], 10),
        "",
        "Readout: NEPdb is already positive for NMI vs public comparators; TESLA/dbPepNeo remain coverage-limited and must not be oversold.",
        "",
        "## Source rescue queue",
        top_table(rescue_summary, ["dataset_source", "rows", "positives", "top_priority", "median_priority", "mean_clean", "mean_public_best", "mean_delta"], 10),
        "",
        "Interpretation: use the rescue queue for wet-lab discussion and failure-case triage. CEDAR and dbPepNeo rows dominate the current priority list because they have the strongest clean score and disagreement signal, while TESLA_mmc7_validation remains weaker and should stay in the control bucket unless new evidence appears.",
        "",
        "## Collaborator handoff",
        "Positive rows:",
        top_table(collab_pos, ["candidate_id", "dataset_source", "peptide_mut", "hla_allele", "label_immunogenicity", "public_best", "clean_science_score", "production_stack_score", "local_rescue_delta_vs_public", "model_disagreement_score", "rescue_priority"], 18),
        "",
        "Negative controls:",
        top_table(collab_ctrl, ["candidate_id", "dataset_source", "peptide_mut", "hla_allele", "label_immunogenicity", "public_best", "clean_science_score", "production_stack_score", "local_rescue_delta_vs_public", "model_disagreement_score", "rescue_priority"], 12),
        "",
        "## 1. What existing public models are usable now?",
        top_table(usable_public, ["model_name", "type", "runnable", "license_note", "role"]),
        "",
        "Existing artifacts were found for several public comparators, especially MHCflurry, BigMHC-IM, PRIME, and some NetMHCpan-formatted outputs. NetMHCpan/PRIME family tools remain license-sensitive and must be treated as frozen comparators/features, not clean training features.",
        "",
        "## 2. Which local algorithms already exist in this repo?",
        top_table(local, ["model_name", "type", "role", "runnable", "clean_track_allowed", "production_track_allowed"], 20),
        "",
        "## 3. Which local algorithm contributes most under strict no-leakage evaluation?",
        top_table(clean_lb, ["model", "split", "n", "positives", "AUPRC", "AUROC", "patient_hit_rate@20", "patient_recall@20"], 10),
        "",
        "Brutal boundary: the integrated run can rank available local artifacts, but a publishable clean-algorithm claim requires the row to be explicitly strict/no-overlap/source-heldout and free of public predictor scores.",
        "",
        "## 4. Does quantum_kernel_no_anchor_gamma1.0 survive source-heldout validation?",
        answer_quantum(clean_lb),
        "",
        "Historical JSON summaries show strong random/CV signal can degrade under leave-source-out conditions; therefore source-heldout survival is a must-pass gate, not an assumption.",
        "",
        "## 5. Does Wave8_TCR_SelfSim_full improve immunogenicity ranking beyond HLA presentation?",
        answer_wave8(clean_lb),
        "",
        "This is the most important scientific branch because vaccine ranking bottleneck is immunogenicity, not binding alone. But it must be shown against presentation-only comparators under patient/source leakage controls.",
        "",
        "## 6. Does production stacking beat BigMHC-IM / PRIME / MHCflurry / NetMHCpan?",
        top_table(prod_lb, ["model", "split", "n", "positives", "AUPRC", "AUROC", "Precision@20", "Recall@20", "patient_hit_rate@20"], 12),
        "",
        "Boundary: if production stack wins, the claim is practical integration, not a new clean predictor. Public predictor scores are frozen features and may carry training-corpus contamination risk.",
        "",
        "## 7. Which patients have top-20 rescued candidates?",
        top_table(top20, ["patient_id", "candidate_id", "dataset_source", "gene", "peptide_mut", "hla_allele", "production_stack_score", "clean_science_score", "model_disagreement_score"], 20),
        "",
        "## 8. Which candidates should be discussed with wet-lab collaborators?",
        "Prioritize top-20 candidates with high production score, adequate expression/VAF, low leakage risk, model agreement, and a positive local/TCR/structure branch. Exclude or down-rank HLA-LOH/APM-loss and high-leakage candidates unless the goal is methods debugging.",
        "",
        "## 9. What is the next minimal experiment to validate?",
        "- Assemble a small patient-level candidate set with matched tumor expression, VAF/clonality, HLA typing, and APM/HLA-LOH metadata.",
        "- Run frozen public predictors and local CLEAN-Neo++ branches.",
        "- Select top 20 plus disagreement controls per patient.",
        "- Validate presentation with MS where feasible and immunogenicity with T-cell assay where labels are claimed.",
        "",
        "## 10. What should be the paper/사업 framing?",
        "- Paper: leakage-aware patient-level immunogenicity ranking; endpoint is patient-level Recall@20/hit rate, secondary AUPRC under strict no-overlap/source-heldout splits.",
        "- Business: operating layer that standardizes public tools, local algorithms, leakage audit, candidate explanation, and wet-lab handoff.",
        "",
        "## Leakage summary",
        top_table(leak, ["leakage_control", "n_flagged", "denominator", "fraction"], 20),
        "",
        "## Data scale",
        f"- Canonical candidate rows: {len(canon):,}",
        f"- Sources: {canon['dataset_source'].nunique() if 'dataset_source' in canon else 0:,}",
        f"- Patients: {canon['patient_id'].nunique() if 'patient_id' in canon else 0:,}",
        f"- Real non-placeholder patient rows: {real_patient_rows:,}",
        f"- Real non-placeholder unique patients: {real_patient_unique:,}",
        "",
        "## Existing public leakage context from hub",
        f"- Hub total rows: {hub_leak.get('total_rows', 'NA')}",
        f"- Hub leakage-free rows: {hub_leak.get('leakage_free_n', 'NA')}",
        f"- Hub test-set safety counts: `{hub_leak.get('test_set_safety_counts', {})}`",
        "",
        "## Existing representation context",
        f"- QX 5-fold CV summary: `{qx.get('R_5fold_cv', {})}`",
        f"- ESM2 full summary: `{esm.get('QQ_esm2_full', {})}`",
        "",
        "## Command examples",
        "```bash",
        "OUT=project/results/neoimmune_stack_2026_05_10",
        "source /home/seungho/personal/THCA_data_analysis/.venv/bin/activate",
        "python project/scripts/neoimmune_stack/00_repo_audit.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/01_build_model_registry.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/02_build_canonical_candidate_table.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/03_run_external_adapters.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/04_collect_local_model_outputs.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/05_leakage_audit.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/06_train_clean_ranker.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/07_train_production_stack.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/08_evaluate_patient_topn.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/09_failure_case_audit.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/17_apply_ga_rl_algorithm.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/18_train_nmi_clean_method.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/19_build_imneo_public_reconstruction.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/20_source_transfer_gap_audit.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/21_build_source_rescue_queue.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/22_build_collaborator_handoff.py --outdir $OUT",
        "python project/scripts/neoimmune_stack/10_generate_strategy_report.py --outdir $OUT",
        "```",
        "",
        "## Claim boundaries",
        "- No clinical vaccine efficacy claim.",
        "- No antigen-presentation claim without MS immunopeptidomics or equivalent evidence.",
        "- No immunogenicity claim without T-cell assay labels.",
        "- HLA binding is not immunogenicity.",
        "- Public predictor features are excluded from clean-science training.",
        "- Leakage problems are surfaced, not hidden.",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "NEOIMMUNE_STACK_STRATEGY_REPORT.md")


if __name__ == "__main__":
    main()
