#!/usr/bin/env python3
"""Build a compact HTML visual dossier for CROSS-Neo MD audit figures."""

from __future__ import annotations

import html
import shutil
from pathlib import Path

import pandas as pd

from common_md import OUT, REPO


HUB = REPO / "project/papers_hub_2026_05_04"
WEB = Path("/var/www/papers/papers_hub_2026_05_04")
ASSET = HUB / "assets/cross_neo_md_audit"
WEB_ASSET = WEB / "assets/cross_neo_md_audit"
PAGE = HUB / "cross_neo_md_audit_dossier.html"
WEB_PAGE = WEB / "cross_neo_md_audit_dossier.html"


def esc(x) -> str:
    return html.escape("" if pd.isna(x) else str(x))


def fmt(x, digits=3) -> str:
    try:
        if pd.isna(x):
            return "NA"
        return f"{float(x):.{digits}f}"
    except Exception:
        return esc(x)


def table(df: pd.DataFrame, cols: list[str], max_rows: int = 20) -> str:
    if df.empty:
        return "<p class='muted'>No rows.</p>"
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return "<p class='muted'>Requested columns are absent.</p>"
    rows = ["<table><thead><tr>"]
    rows += [f"<th>{esc(c.replace('_', ' '))}</th>" for c in cols]
    rows.append("</tr></thead><tbody>")
    for _, r in df.head(max_rows).iterrows():
        rows.append("<tr>")
        for c in cols:
            v = r.get(c, "")
            rows.append(f"<td>{fmt(v) if isinstance(v, (int, float)) else esc(v)}</td>")
        rows.append("</tr>")
    rows.append("</tbody></table>")
    return "\n".join(rows)


def main() -> None:
    ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    for p in sorted((OUT / "figures").glob("*.png")):
        shutil.copy2(p, ASSET / p.name)
        shutil.copy2(p, WEB_ASSET / p.name)
    for p in sorted((OUT / "figures").glob("fig_md*.pdf")):
        shutil.copy2(p, ASSET / p.name)
        shutil.copy2(p, WEB_ASSET / p.name)
    for p in [
        OUT / "CROSS_Neo_MD_audit_decision_report.md",
        OUT / "MD_one_page_paper_summary.md",
        OUT / "MD_one_page_korean_summary_for_Seungho.md",
        OUT / "integrated/md_case_audit_report.md",
        OUT / "integrated/md_supported_top_candidates.tsv",
        OUT / "integrated/model_high_md_low_cases.tsv",
        OUT / "integrated/model_low_md_high_cases.tsv",
        OUT / "next_simulation_batch_recommendation.tsv",
        OUT / "counterfactual_design/counterfactual_md_batch_plan.md",
        OUT / "counterfactual_design/counterfactual_md_batch_manifest.tsv",
        OUT / "counterfactual_design/candidate_control_sequences.tsv",
        OUT / "counterfactual_design/candidate_decoy_controls.fasta",
        OUT / "counterfactual_design/positive_control_selection_report.md",
        OUT / "counterfactual_design/positive_control_candidates.tsv",
        OUT / "counterfactual_design/positive_control_selected_manifest_rows.tsv",
        OUT / "counterfactual_design/positive_control_pdb_extraction_report.md",
        OUT / "counterfactual_design/positive_control_pdb_extraction.tsv",
        OUT / "counterfactual_design/structure_prep_job_report.md",
        OUT / "counterfactual_design/structure_prep_job_manifest.tsv",
        OUT / "counterfactual_design/openmm_ready_structure_jobs.tsv",
        OUT / "counterfactual_design/openmm_ready_command_templates.tsv",
        OUT / "MD_control_readiness_report.md",
        OUT / "ultra_priority/ultra_wetlab_priority_report.md",
        OUT / "ultra_priority/ultra_wetlab_priority_candidates.tsv",
        OUT / "ultra_priority/ultra_wetlab_assay_plan.tsv",
        OUT / "ultra_priority/ultra_wetlab_priority_summary.json",
        OUT / "small_dataset_uncertainty/small_dataset_uncertainty_funnel_report.md",
        OUT / "small_dataset_uncertainty/small_dataset_uncertainty_funnel.tsv",
        OUT / "small_dataset_uncertainty/small_dataset_wetlab_survivors.tsv",
        OUT / "small_dataset_uncertainty/small_dataset_culled_candidates.tsv",
        OUT / "small_dataset_uncertainty/culling_decision_counts.tsv",
        OUT / "small_dataset_uncertainty/bayesian_posterior_scores.tsv",
        OUT / "small_dataset_uncertainty/dropout_perturbation_summary.tsv",
        OUT / "dl_first_funnel/dl_first_funnel_report.md",
        OUT / "dl_first_funnel/dl_first_candidate_funnel.tsv",
        OUT / "dl_first_funnel/dl_first_stage_counts.tsv",
        OUT / "dl_first_funnel/dl_first_decision_counts.tsv",
        OUT / "dl_first_funnel/dl_first_survivors_for_tcr_structure_review.tsv",
        OUT / "dl_first_funnel/dl_first_md_escalation_candidates.tsv",
        OUT / "dl_first_funnel/dl_first_wetlab_shortlist.tsv",
        OUT / "dl_first_funnel/dl_first_culled_candidates.tsv",
        OUT / "watchers/hmtevvrhc_autowatch_latest.json",
    ]:
        if p.exists():
            shutil.copy2(p, ASSET / p.name)
            shutil.copy2(p, WEB_ASSET / p.name)

    qc = pd.read_csv(OUT / "qc/md_qc_summary.tsv", sep="\t")
    scores = pd.read_csv(OUT / "md_evidence_scores.tsv", sep="\t")
    anchors = pd.read_csv(OUT / "contacts/pmhc_anchor_contacts.tsv", sep="\t")
    tcr_tail = pd.read_csv(OUT / "contacts/tcr_peptide_contact_timeseries.tsv", sep="\t").groupby("run_id", as_index=False).tail(1)
    integrated = pd.read_csv(OUT / "integrated/md_supported_top_candidates.tsv", sep="\t") if (OUT / "integrated/md_supported_top_candidates.tsv").exists() else pd.DataFrame()
    high_low = pd.read_csv(OUT / "integrated/model_high_md_low_cases.tsv", sep="\t") if (OUT / "integrated/model_high_md_low_cases.tsv").exists() else pd.DataFrame()
    next_batch = pd.read_csv(OUT / "next_simulation_batch_recommendation.tsv", sep="\t") if (OUT / "next_simulation_batch_recommendation.tsv").exists() else pd.DataFrame()
    cf_batch = pd.read_csv(OUT / "counterfactual_design/counterfactual_md_batch_manifest.tsv", sep="\t") if (OUT / "counterfactual_design/counterfactual_md_batch_manifest.tsv").exists() else pd.DataFrame()
    pos_controls = pd.read_csv(OUT / "counterfactual_design/positive_control_pdb_extraction.tsv", sep="\t") if (OUT / "counterfactual_design/positive_control_pdb_extraction.tsv").exists() else pd.DataFrame()
    prep_jobs = pd.read_csv(OUT / "counterfactual_design/openmm_ready_structure_jobs.tsv", sep="\t") if (OUT / "counterfactual_design/openmm_ready_structure_jobs.tsv").exists() else pd.DataFrame()
    ultra = pd.read_csv(OUT / "ultra_priority/ultra_wetlab_priority_candidates.tsv", sep="\t") if (OUT / "ultra_priority/ultra_wetlab_priority_candidates.tsv").exists() else pd.DataFrame()
    funnel = pd.read_csv(OUT / "small_dataset_uncertainty/small_dataset_uncertainty_funnel.tsv", sep="\t") if (OUT / "small_dataset_uncertainty/small_dataset_uncertainty_funnel.tsv").exists() else pd.DataFrame()
    cull_counts = pd.read_csv(OUT / "small_dataset_uncertainty/culling_decision_counts.tsv", sep="\t") if (OUT / "small_dataset_uncertainty/culling_decision_counts.tsv").exists() else pd.DataFrame()
    dl_first = pd.read_csv(OUT / "dl_first_funnel/dl_first_candidate_funnel.tsv", sep="\t") if (OUT / "dl_first_funnel/dl_first_candidate_funnel.tsv").exists() else pd.DataFrame()
    dl_first_stages = pd.read_csv(OUT / "dl_first_funnel/dl_first_stage_counts.tsv", sep="\t") if (OUT / "dl_first_funnel/dl_first_stage_counts.tsv").exists() else pd.DataFrame()
    dl_first_counts = pd.read_csv(OUT / "dl_first_funnel/dl_first_decision_counts.tsv", sep="\t") if (OUT / "dl_first_funnel/dl_first_decision_counts.tsv").exists() else pd.DataFrame()
    figs = [p.name for p in sorted((OUT / "figures").glob("fig_md*.png"))]
    fig_cards = "\n".join(
        f"<article><img src='assets/cross_neo_md_audit/{esc(name)}'><h3>{esc(name.replace('.png',''))}</h3><a href='assets/cross_neo_md_audit/{esc(name)}'>PNG</a></article>"
        for name in figs
    )

    html_text = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CROSS-Neo MD Audit</title>
<style>
body{{margin:0;background:#07111f;color:#edf5ff;font-family:Inter,Arial,sans-serif;line-height:1.55}}
a{{color:#39d4b5}} .hero{{padding:48px 32px;background:#0f2035;border-bottom:1px solid #26364d}}
.wrap{{max-width:1280px;margin:0 auto;padding:28px 30px 70px}} h1{{font-size:48px;margin:0 0 10px;font-family:Georgia,serif}}
.lead{{max-width:980px;color:#cbd7e7;font-size:17px}} .grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:18px 0}}
.stat{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:14px}} .stat b{{display:block;color:#f2c46d;font-size:26px}}
section{{border-bottom:1px solid #26364d;padding:24px 0}} h2{{font-family:Georgia,serif;font-size:30px;color:#fff2d0}}
table{{border-collapse:collapse;width:100%;font-size:12px}} th,td{{border:1px solid #293b55;padding:7px;vertical-align:top}} th{{color:#f2c46d;background:#13243a}}
.figs{{display:grid;grid-template-columns:repeat(2,1fr);gap:16px}} article{{background:#101d30;border:1px solid #293b55;border-radius:8px;padding:12px}}
img{{width:100%;background:white;border-radius:6px}} .warn{{border-left:4px solid #ff7b72;background:#111d2e;padding:12px 14px;margin:14px 0}}
.links{{display:flex;gap:10px;flex-wrap:wrap}} .pill{{border:1px solid #395170;border-radius:999px;padding:7px 10px;background:#101d30}}
.muted{{color:#9fb0c7}} @media(max-width:900px){{.grid,.figs{{grid-template-columns:1fr}}}}
</style></head><body>
<header class="hero"><div class="wrap">
<h1>CROSS-Neo MD Audit</h1>
<p class="lead">Structure-aware MD visualization layer for TCR-pMHC neoantigen candidates. This page summarizes synced OpenMM trajectories, QC, contact persistence, and claim boundaries.</p>
<div class="grid">
<div class="stat"><b>{len(qc)}</b><span>runs indexed</span></div>
<div class="stat"><b>{int(qc['trajectory_available'].sum())}</b><span>trajectories parsed</span></div>
<div class="stat"><b>{fmt(qc['total_time_ps'].max()/1000, 2)}</b><span>max ns synced</span></div>
<div class="stat"><b>{len(figs)}</b><span>main figure panels</span></div>
</div>
<p><a href="index.html">Hub index</a></p>
</div></header>
<main class="wrap">
<section><h2>Current Verdict</h2>
<div class="warn">MD is a structural audit layer, not immunogenicity proof. Completed 10 ns evidence is currently available for GADGVGKSAL / HLA-C*08:02; HMTEVVRHC 10 ns was still partial at the latest sync.</div>
<div class="links">
<a class="pill" href="assets/cross_neo_md_audit/CROSS_Neo_MD_audit_decision_report.md">Decision report</a>
<a class="pill" href="assets/cross_neo_md_audit/MD_one_page_paper_summary.md">Paper summary</a>
<a class="pill" href="assets/cross_neo_md_audit/MD_one_page_korean_summary_for_Seungho.md">Korean summary</a>
<a class="pill" href="assets/cross_neo_md_audit/md_case_audit_report.md">Model + MD audit</a>
<a class="pill" href="assets/cross_neo_md_audit/counterfactual_md_batch_plan.md">Counterfactual batch</a>
<a class="pill" href="assets/cross_neo_md_audit/positive_control_selection_report.md">Positive controls</a>
<a class="pill" href="assets/cross_neo_md_audit/structure_prep_job_report.md">Structure prep jobs</a>
<a class="pill" href="assets/cross_neo_md_audit/MD_control_readiness_report.md">Control readiness</a>
<a class="pill" href="assets/cross_neo_md_audit/ultra_wetlab_priority_report.md">Ultra priority</a>
<a class="pill" href="assets/cross_neo_md_audit/small_dataset_uncertainty_funnel_report.md">Bayesian/dropout funnel</a>
<a class="pill" href="assets/cross_neo_md_audit/dl_first_funnel_report.md">DL-first funnel</a>
<a class="pill" href="assets/cross_neo_md_audit/hmtevvrhc_autowatch_latest.json">6VRN watcher JSON</a>
</div>
</section>
<section><h2>QC Summary</h2>{table(qc, ['run_id','candidate','condition','n_frames','total_time_ps','peptide_rmsd_final_nm','peptide_com_drift_final_nm','qc_status'], 20)}</section>
<section><h2>MD Evidence Scores</h2>{table(scores, ['run_id','candidate','MD_evidence_score','MD_evidence_label','pMHC_stability_score','TCR_recognition_score','runtime_fraction'], 20)}</section>
<section><h2>Anchor Contacts</h2>{table(anchors, ['run_id','candidate','peptide_position','peptide_resname','is_anchor','max_contact_occupancy','sum_contact_occupancy'], 30)}</section>
<section><h2>TCR Contact Tail</h2>{table(tcr_tail, ['run_id','candidate','time_ps','tcr_peptide_heavy_atom_contacts'], 20)}</section>
<section><h2>CROSS-Neo + MD Agreement</h2>
<p class="muted">Rows where model score and completed MD evidence both support candidate prioritization.</p>
{table(integrated, ['row_id','source_dataset','peptide','hla_4digit','label','score','model_name','split_name','MD_evidence_label','MD_evidence_score','wetlab_priority_score_evidence_adjusted'], 12)}
<p class="muted">High model score but MD is low or still insufficient.</p>
{table(high_low, ['row_id','source_dataset','peptide','hla_4digit','label','score','model_name','MD_evidence_label','md_interpretation'], 12)}
</section>
<section><h2>Next Simulation Batch</h2>{table(next_batch, ['priority_rank','row_id','peptide','hla_4digit','md_tier','recommended_next_batch','reason'], 12)}</section>
<section><h2>Counterfactual Controls</h2>
<p class="muted">WT rows are blocked or manual-confirmation required unless a registry WT peptide is available. Decoys preserve P2 and C-terminal anchors.</p>
{table(cf_batch, ['batch_id','priority_rank','peptide','hla_4digit','control_type','complex_kind','sequence','sequence_status','readiness','blocked_by'], 18)}
</section>
<section><h2>Positive Control PDBs</h2>
<p class="muted">These are analysis sanity controls only, not cancer specificity evidence.</p>
{table(pos_controls, ['target_peptide','control_peptide','control_hla_4digit','pdb_id','status','selected_chains','contacts_tcr_peptide_initial','extracted_pdb'], 12)}
</section>
<section><h2>OpenMM-Ready Structure Jobs</h2>
<p class="muted">Command templates are generated, but no new MD runs are launched from this dossier.</p>
{table(prep_jobs, ['prep_id','batch_id','peptide','hla_4digit','control_type','complex_kind','sequence','structure_route','input_template_pdb'], 16)}
</section>
<section><h2>Ultra Wetlab Priority</h2>
<p class="muted">This is a conservative recommendation score, not an immunogenicity probability.</p>
{table(ultra, ['row_id','peptide','hla_4digit','recommendation_tier','ultra_priority_score','confidence_score','model_evidence_score','tcr_resource_score','md_structural_score','control_readiness_score','claim_risk_penalty','why'], 15)}
</section>
<section><h2>Small-Dataset Uncertainty Funnel</h2>
<p class="muted">Bayesian posterior and dropout/weight perturbation are used to cull weak or unstable candidates, not to certify positives.</p>
{table(cull_counts, ['decision','n'], 12)}
{table(funnel, ['row_id','peptide','hla_4digit','funnel_decision','funnel_priority_score','bayes_mean','bayes_q05','bayes_q95','perturb_median','perturb_width_90','dropout_sensitivity','funnel_reason'], 15)}
</section>
<section><h2>DL-First Candidate Funnel</h2>
<p class="muted">Cheap model ensemble and uncertainty gates run before TCR rescue, structure, and MD. TCR-discordant high-model rows are audit/control candidates, not positives.</p>
{table(dl_first_stages, ['stage','n'], 12)}
{table(dl_first_counts, ['decision','n'], 12)}
{table(dl_first, ['row_id','peptide','hla_4digit','dl_first_decision','next_action','dl_first_priority_score','main_dl_score','bayes_mean','bayes_q05','bayes_q95','perturb_prob_gt_050','tcr_augmented_score_mean','paired_tcr_evidence_count','md_label'], 18)}
</section>
<section><h2>Figure Gallery</h2><div class="figs">{fig_cards}</div></section>
<section><h2>Source Paths</h2><p class="muted">{esc(OUT)}</p></section>
</main></body></html>"""
    PAGE.write_text(html_text)
    WEB.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(PAGE, WEB_PAGE)
    print(f"[md-dossier] wrote {PAGE}")
    print(f"[md-dossier] deployed {WEB_PAGE}")


if __name__ == "__main__":
    main()
