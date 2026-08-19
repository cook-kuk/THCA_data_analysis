#!/usr/bin/env python3
"""Build a presentation/manuscript storyboard package for CROSS-Neo DL/TCR/MD.

Outputs include:
- end-to-end overview figure
- data usage matrix
- representative candidate board
- figure/table explanation master
- data lineage manifest
- Korean overview report
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
CROSS = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
MD_OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
OUT = MD_OUT / "decision_storyboard"
FIG = MD_OUT / "figures"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fmt_pct(x: float) -> str:
    return f"{100*x:.1f}%"


def nrows(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        return len(pd.read_csv(path, sep="\t"))
    except Exception:
        return 0


def build_data_usage_manifest() -> pd.DataFrame:
    rows = [
        {
            "artifact": "TCR wetlab candidate table",
            "path": str(CROSS / "tcr_extension/tcr_wetlab_candidate_prioritization_unique_pmhc.tsv"),
            "rows": nrows(CROSS / "tcr_extension/tcr_wetlab_candidate_prioritization_unique_pmhc.tsv"),
            "key_columns": "row_id, peptide, hla_4digit, label_binary, pmhc_score_mean, tcr_augmented_score_mean",
            "used_for": "candidate universe, labels, pMHC/TCR branch scores",
            "pipeline_stage": "input_registry",
            "claim_boundary": "labels are dataset labels, not direct clinical immunogenicity proof",
        },
        {
            "artifact": "CROSS-Neo prediction ensemble",
            "path": str(CROSS / "predictions/all_predictions.tsv"),
            "rows": nrows(CROSS / "predictions/all_predictions.tsv"),
            "key_columns": "row_id, score, model_name, model_family, split_name, rank_pct",
            "used_for": "main DL ensemble aggregation and rank evidence",
            "pipeline_stage": "DL_first_screen",
            "claim_boundary": "model score is prioritization evidence, not a calibrated probability",
        },
        {
            "artifact": "Bayesian/dropout uncertainty tables",
            "path": str(MD_OUT / "small_dataset_uncertainty"),
            "rows": nrows(MD_OUT / "small_dataset_uncertainty/small_dataset_uncertainty_funnel.tsv"),
            "key_columns": "bayes_mean, bayes_q05, bayes_q95, perturb_prob_gt_050, dropout_sensitivity",
            "used_for": "remove unstable or low-posterior candidates before expensive layers",
            "pipeline_stage": "uncertainty_gate",
            "claim_boundary": "uncertainty summaries are internal triage metrics",
        },
        {
            "artifact": "DL-first funnel with labels",
            "path": str(MD_OUT / "dl_first_funnel_ui/dl_first_stage_label_counts.tsv"),
            "rows": nrows(MD_OUT / "dl_first_funnel_ui/dl_first_stage_label_counts.tsv"),
            "key_columns": "stage, n, positive, negative, positive_fraction",
            "used_for": "show how positive/negative composition changes across gates",
            "pipeline_stage": "dashboard_and_reporting",
            "claim_boundary": "stage enrichment is retrospective on available labels",
        },
        {
            "artifact": "Threshold optimizer presets",
            "path": str(MD_OUT / "dl_threshold_optimization/threshold_recommended_presets.tsv"),
            "rows": nrows(MD_OUT / "dl_threshold_optimization/threshold_recommended_presets.tsv"),
            "key_columns": "preset_name, call_rule, TP, TN, FP, FN, precision, recall",
            "used_for": "label-aware operating point selection for slider presets",
            "pipeline_stage": "threshold_optimization",
            "claim_boundary": "optimized on current labels; external validation required",
        },
        {
            "artifact": "TCR/structure/MD escalation candidates",
            "path": str(MD_OUT / "dl_first_funnel/dl_first_md_escalation_candidates.tsv"),
            "rows": nrows(MD_OUT / "dl_first_funnel/dl_first_md_escalation_candidates.tsv"),
            "key_columns": "row_id, peptide, hla_4digit, label_binary, dl_first_decision, md_label",
            "used_for": "identify candidates to send into structure/MD or wetlab controls",
            "pipeline_stage": "escalation",
            "claim_boundary": "escalation means higher priority, not confirmed positive",
        },
        {
            "artifact": "Baker/Rosetta fallback interface scores",
            "path": str(MD_OUT / "baker_rosetta_filter/baker_rosetta_fallback_interface_scores.tsv"),
            "rows": nrows(MD_OUT / "baker_rosetta_filter/baker_rosetta_fallback_interface_scores.tsv"),
            "key_columns": "fallback_structural_score, pmhc_residue_pair_contacts_4A, tcr_peptide_residue_pair_contacts_4A",
            "used_for": "cheap structure sanity check before Rosetta/MD",
            "pipeline_stage": "cheap_structural_filter",
            "claim_boundary": "static contacts are not Rosetta energies or immunogenicity proof",
        },
        {
            "artifact": "OpenMM MD status and evidence",
            "path": str(MD_OUT / "md_status_summary.tsv") + "; " + str(MD_OUT / "md_evidence_scores.tsv"),
            "rows": nrows(MD_OUT / "md_evidence_scores.tsv"),
            "key_columns": "candidate, time_ps, completion_fraction, MD_evidence_score, MD_evidence_label, peptide_rmsd_final_nm, tcr_peptide_contacts_tail",
            "used_for": "late structural audit of peptide-MHC/TCR-pMHC stability",
            "pipeline_stage": "MD_audit",
            "claim_boundary": "MD supports structural plausibility, not immune recognition proof",
        },
    ]
    return pd.DataFrame(rows)


def build_figure_table_manifest() -> pd.DataFrame:
    rows = [
        {
            "id": "Fig. A / fig_md24",
            "title": "End-to-end DL-first recognition-aware funnel",
            "source_files": "dl_first_stage_label_counts.tsv; threshold_recommended_presets.tsv",
            "what_it_shows": "649 candidates are narrowed by cheap DL, uncertainty, TCR, structure, and MD gates.",
            "main_message": "The pipeline turns a broad labeled candidate set into a tiny, auditable wetlab/MD shortlist.",
            "use_in_slide": "Opening overview figure",
            "claim_boundary": "Enrichment and prioritization, not proof of immunogenicity.",
        },
        {
            "id": "Fig. B / fig_md25",
            "title": "Data usage and evidence lineage matrix",
            "source_files": "data_usage_manifest.tsv",
            "what_it_shows": "Which data artifacts feed each model, uncertainty, TCR, structure, MD, optimizer, and UI layer.",
            "main_message": "Every claim has a traceable source table and a stated boundary.",
            "use_in_slide": "Methods/reproducibility slide",
            "claim_boundary": "Lineage map does not validate model performance by itself.",
        },
        {
            "id": "Fig. C / fig_md26",
            "title": "Representative candidate evidence board",
            "source_files": "dl_first_md_escalation_candidates.tsv; baker_rosetta_fallback_interface_scores.tsv; md_evidence_scores.tsv; md_status_summary.tsv",
            "what_it_shows": "GADGVGKSAL and HMTEVVRHC across main DL, Bayesian uncertainty, TCR branch, Baker fallback, and completed MD.",
            "main_message": "HMTEVVRHC now has the strongest completed MD signal, while GADGVGKSAL remains the more balanced DL+TCR+MD candidate with controls queued.",
            "use_in_slide": "Candidate selection/case study slide",
            "claim_boundary": "Candidate evidence does not establish clinical utility.",
        },
        {
            "id": "Table 1",
            "title": "Stage-by-stage label composition",
            "source_files": "dl_first_stage_label_counts.tsv",
            "what_it_shows": "Counts of positive/negative rows after each gate.",
            "main_message": "Positive fraction rises from 24.8% in the full set to 100% in MD escalation/wetlab shortlist under current settings.",
            "use_in_slide": "Funnel performance table",
            "claim_boundary": "Current-label enrichment only; external holdout needed for final claims.",
        },
        {
            "id": "Table 2",
            "title": "Threshold operating presets",
            "source_files": "threshold_recommended_presets.tsv",
            "what_it_shows": "TP/TN/FP/FN, precision, recall, and F1 for optimized slider presets.",
            "main_message": "A no-false-positive preset found 13 labeled positives with 0 false positives in the current table.",
            "use_in_slide": "Decision threshold slide",
            "claim_boundary": "Preset is optimized on current labels and must be validated externally.",
        },
        {
            "id": "Table 3",
            "title": "Data usage manifest",
            "source_files": "data_usage_manifest.tsv",
            "what_it_shows": "Where each dataset/result is used and what claims it can support.",
            "main_message": "The pipeline is auditable from input label table to final UI.",
            "use_in_slide": "Supplement/methods table",
            "claim_boundary": "Documentation table, not performance evidence.",
        },
    ]
    return pd.DataFrame(rows)


def build_representative_candidate_table() -> pd.DataFrame:
    esc = read_tsv(MD_OUT / "dl_first_funnel/dl_first_md_escalation_candidates.tsv")
    baker = read_tsv(MD_OUT / "baker_rosetta_filter/baker_rosetta_fallback_interface_scores.tsv")
    md = read_tsv(MD_OUT / "md_status_summary.tsv")
    md_scores = read_tsv(MD_OUT / "md_evidence_scores.tsv")
    if esc.empty:
        return pd.DataFrame()
    baker_agg = (
        baker.groupby("row_id", as_index=False)
        .agg(
            baker_structural_score=("fallback_structural_score", "max"),
            pmhc_contacts=("pmhc_residue_pair_contacts_4A", "max"),
            tcr_peptide_contacts=("tcr_peptide_residue_pair_contacts_4A", "max"),
        )
        if not baker.empty
        else pd.DataFrame(columns=["row_id"])
    )
    md["pmhc_key"] = md["peptide"].astype(str).str.upper() + "|" + md["hla"].astype(str).str.upper()
    esc["pmhc_key"] = esc["peptide"].astype(str).str.upper() + "|" + esc["hla_4digit"].astype(str).str.upper()
    md_agg = (
        md.groupby("pmhc_key", as_index=False)
        .agg(time_ps=("time_ps", "max"), completion_fraction=("completion_fraction", "max"), temperature_k=("temperature_k", "last"))
        if not md.empty
        else pd.DataFrame(columns=["pmhc_key"])
    )
    if not md_scores.empty:
        md_scores = md_scores.copy()
        md_scores["pmhc_key"] = md_scores["candidate"].astype(str).str.replace("/", "|", regex=False).str.upper()
        md_scores["_runtime"] = pd.to_numeric(md_scores.get("runtime_fraction", 0), errors="coerce").fillna(0)
        md_scores["_score"] = pd.to_numeric(md_scores.get("MD_evidence_score", 0), errors="coerce").fillna(0)
        md_scores = md_scores.sort_values(["pmhc_key", "_runtime", "_score"], ascending=[True, False, False])
        md_scores = md_scores.groupby("pmhc_key", as_index=False).head(1)
        md_scores = md_scores.rename(
            columns={
                "run_id": "current_md_run_id",
                "MD_evidence_score": "current_md_evidence_score",
                "MD_evidence_label": "current_md_evidence_label",
                "runtime_fraction": "current_md_runtime_fraction",
                "peptide_rmsd_final_nm": "current_md_peptide_rmsd_final_nm",
                "tcr_peptide_contacts_tail": "current_md_tcr_peptide_contacts_tail",
            }
        )
        md_scores = md_scores[
            [
                "pmhc_key",
                "current_md_run_id",
                "current_md_evidence_score",
                "current_md_evidence_label",
                "current_md_runtime_fraction",
                "current_md_peptide_rmsd_final_nm",
                "current_md_tcr_peptide_contacts_tail",
            ]
        ]
    else:
        md_scores = pd.DataFrame(columns=["pmhc_key"])
    out = esc.merge(baker_agg, on="row_id", how="left").merge(md_agg, on="pmhc_key", how="left").merge(md_scores, on="pmhc_key", how="left")
    cols = [
        "row_id",
        "peptide",
        "hla_4digit",
        "label_binary",
        "dl_first_decision",
        "main_dl_score",
        "bayes_mean",
        "bayes_q05",
        "bayes_q95",
        "tcr_augmented_score_mean",
        "paired_tcr_evidence_count",
        "baker_structural_score",
        "pmhc_contacts",
        "tcr_peptide_contacts",
        "md_label",
        "current_md_evidence_label",
        "current_md_evidence_score",
        "current_md_run_id",
        "current_md_peptide_rmsd_final_nm",
        "current_md_tcr_peptide_contacts_tail",
        "time_ps",
        "completion_fraction",
        "dl_first_reason",
    ]
    return out[[c for c in cols if c in out.columns]]


def md_summary_for(candidates: pd.DataFrame, peptide: str, hla: str) -> tuple[str, str]:
    if candidates.empty:
        return "NA", "NA"
    mask = (candidates["peptide"].astype(str).str.upper() == peptide.upper()) & (
        candidates["hla_4digit"].astype(str).str.upper() == hla.upper()
    )
    if not mask.any():
        return "NA", "NA"
    row = candidates[mask].iloc[0]
    label = row.get("current_md_evidence_label") or row.get("md_label") or "NA"
    try:
        score = f"{float(row.get('current_md_evidence_score')):.3f}"
    except Exception:
        score = "NA"
    return str(label), score


def make_overview_figure(stage_counts: pd.DataFrame, presets: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.axis("off")
    stages = [
        ("Candidate\nregistry", "649\n161+ / 488-"),
        ("DL length\nsupport", "491\n161+ / 330-"),
        ("Main DL\nbroad pass", "280\n113+ / 167-"),
        ("Bayesian +\ndropout pass", "97\n56+ / 41-"),
        ("TCR branch\nsupported", "18\n11+ / 7-"),
        ("Baker / MD\nsupported", "1-2\n1-2+ / 0-"),
        ("Wetlab / MD\nshortlist", "1-2\npriority cases"),
    ]
    x = np.linspace(0.06, 0.94, len(stages))
    y = 0.58
    colors = ["#13243a", "#183a5a", "#245276", "#2a6f87", "#2f8f7c", "#5c8f3a", "#9f7c2f"]
    for i, ((title, value), xi) in enumerate(zip(stages, x)):
        ax.add_patch(plt.Rectangle((xi - 0.062, y - 0.11), 0.124, 0.22, color=colors[i], ec="#edf5ff", lw=1.2))
        ax.text(xi, y + 0.035, title, ha="center", va="center", color="white", fontsize=11, weight="bold")
        ax.text(xi, y - 0.055, value, ha="center", va="center", color="#f2c46d", fontsize=12)
        if i < len(stages) - 1:
            ax.annotate("", xy=(x[i + 1] - 0.075, y), xytext=(xi + 0.075, y), arrowprops=dict(arrowstyle="->", lw=2, color="#edf5ff"))
    nofp = presets[presets.get("preset_name", "").eq("NO_FALSE_POSITIVE_MAX_TP")] if not presets.empty else pd.DataFrame()
    if not nofp.empty:
        r = nofp.iloc[0]
        txt = f"Optimized preset: {int(r.TP)} TP / {int(r.FP)} FP\\nprecision {r.precision:.3f}, recall {r.recall:.3f}"
    else:
        txt = "Optimized presets unavailable"
    ax.text(0.5, 0.20, txt, ha="center", va="center", fontsize=18, color="#07111f", bbox=dict(boxstyle="round,pad=0.45", fc="#f2c46d", ec="#f2c46d"))
    ax.text(0.5, 0.91, "CROSS-Neo DL-first, TCR-aware, structure-audited candidate funnel", ha="center", fontsize=20, weight="bold", color="#edf5ff")
    fig.patch.set_facecolor("#07111f")
    fig.savefig(FIG / "fig_md24_story_overview_flow.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md24_story_overview_flow.pdf", bbox_inches="tight")
    plt.close(fig)


def make_data_usage_matrix(manifest: pd.DataFrame) -> None:
    stages = ["DL", "Uncertainty", "TCR", "Threshold", "Baker", "MD", "UI/Wetlab"]
    usage = np.zeros((len(manifest), len(stages)))
    for i, r in manifest.iterrows():
        used = f"{r['used_for']} {r['pipeline_stage']}".lower()
        for j, s in enumerate(stages):
            key = s.lower().split("/")[0]
            usage[i, j] = 1 if key in used else 0
        if "dashboard" in used or "ui" in used or "wetlab" in used:
            usage[i, stages.index("UI/Wetlab")] = 1
        if "bayesian" in used or "dropout" in used:
            usage[i, stages.index("Uncertainty")] = 1
        if "structure" in used or "rosetta" in used:
            usage[i, stages.index("Baker")] = 1
    fig, ax = plt.subplots(figsize=(11.8, 5.8))
    ax.imshow(usage, cmap="YlGnBu", vmin=0, vmax=1)
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(stages, rotation=25, ha="right")
    ax.set_yticks(range(len(manifest)))
    ax.set_yticklabels(manifest["artifact"].str.replace(" table", "").str.replace(" scores", ""), fontsize=8)
    for i in range(usage.shape[0]):
        for j in range(usage.shape[1]):
            ax.text(j, i, "●" if usage[i, j] else "", ha="center", va="center", color="#07111f")
    ax.set_title("Data usage matrix: where each artifact enters the pipeline")
    fig.savefig(FIG / "fig_md25_data_usage_matrix.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md25_data_usage_matrix.pdf", bbox_inches="tight")
    plt.close(fig)


def make_candidate_board(candidates: pd.DataFrame) -> None:
    if candidates.empty:
        return
    show = candidates.copy()
    labels = show["peptide"].astype(str) + "\n" + show["hla_4digit"].astype(str)
    metrics = [
        ("main_dl_score", "main DL"),
        ("bayes_mean", "Bayes"),
        ("tcr_augmented_score_mean", "TCR branch"),
        ("baker_structural_score", "Baker/static"),
        ("current_md_evidence_score", "MD evidence"),
    ]
    fig, ax = plt.subplots(figsize=(11, 5.8))
    y = np.arange(len(show))
    width = 0.15
    for k, (col, name) in enumerate(metrics):
        vals = pd.to_numeric(show.get(col, 0), errors="coerce").fillna(0).clip(0, 1)
        ax.barh(y + (k - 2) * width, vals, height=width, label=name)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Normalized evidence score")
    ax.set_title("Representative candidate evidence board")
    ax.legend(frameon=False, ncol=3)
    ax.grid(axis="x", alpha=0.25)
    fig.savefig(FIG / "fig_md26_representative_candidate_board.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIG / "fig_md26_representative_candidate_board.pdf", bbox_inches="tight")
    plt.close(fig)


def write_reports(manifest: pd.DataFrame, figtab: pd.DataFrame, candidates: pd.DataFrame, stage: pd.DataFrame, presets: pd.DataFrame) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(OUT / "data_usage_manifest.tsv", sep="\t", index=False)
    figtab.to_csv(OUT / "figure_table_explanation_master.tsv", sep="\t", index=False)
    candidates.to_csv(OUT / "representative_candidate_evidence_table.tsv", sep="\t", index=False)
    stage.to_csv(OUT / "stage_label_composition_for_slides.tsv", sep="\t", index=False)
    presets.to_csv(OUT / "threshold_preset_summary_for_slides.tsv", sep="\t", index=False)

    nofp = presets[presets.get("preset_name", "").eq("NO_FALSE_POSITIVE_MAX_TP")] if not presets.empty else pd.DataFrame()
    nofp_txt = ""
    if not nofp.empty:
        r = nofp.iloc[0]
        nofp_txt = f"- No-FP optimizer preset: {int(r.TP)} TP, {int(r.FP)} FP, precision {r.precision:.3f}, recall {r.recall:.3f}."
    hm_label, hm_score = md_summary_for(candidates, "HMTEVVRHC", "HLA-A*02:01")
    gad_label, gad_score = md_summary_for(candidates, "GADGVGKSAL", "HLA-C*08:02")
    lines = [
        "# CROSS-Neo Decision Storyboard Package",
        "",
        "## 한 줄 결론",
        "",
        "이 패키지는 CROSS-Neo 후보를 cheap DL, uncertainty, TCR branch, Baker/static structure, MD audit 순서로 줄이는 전체 흐름을 발표/논문용으로 설명한다.",
        "",
        "## 대표 숫자",
        "",
        "- 시작 후보: 649개 (positive 161, negative 488).",
        f"- HMTEVVRHC: 10 ns explicit-solvent MD 완료, current MD label `{hm_label}` (score {hm_score}).",
        f"- GADGVGKSAL: 10 ns explicit-solvent MD 완료, current MD label `{gad_label}` (score {gad_score}), 추가 10 ns replicate/control MD도 완료.",
        "- 현재 default MD escalation: 2개 (TP 2, FP 0).",
        "- 현재 default wetlab shortlist: 1개 (TP 1, FP 0).",
        nofp_txt,
        "",
        "## 추천 Figure 흐름",
        "",
        "1. `fig_md24_story_overview_flow`: 전체 후보 축소 흐름과 핵심 숫자.",
        "2. `fig_md18_dl_first_candidate_funnel`: DL-first gate별 후보 감소.",
        "3. `fig_md22_threshold_precision_recall_frontier`: threshold preset 탐색 결과.",
        "4. `fig_md26_representative_candidate_board`: GADGVGKSAL/HMTEVVRHC 대표 후보 근거와 최신 MD evidence.",
        "5. `fig_md25_data_usage_matrix`: 어떤 데이터가 어느 단계에 쓰였는지.",
        "",
        "## 데이터 사용처",
        "",
        manifest[["artifact", "rows", "used_for", "pipeline_stage", "claim_boundary"]].to_markdown(index=False),
        "",
        "## Figure/Table 설명",
        "",
        figtab[["id", "title", "main_message", "claim_boundary"]].to_markdown(index=False),
        "",
        "## Claim Boundary",
        "",
        "- 이 결과는 후보 우선순위/실험 설계 근거다.",
        "- MD와 구조 점수는 면역원성 증명이 아니다.",
        "- threshold preset은 현재 label table에서 최적화된 값이므로 외부 validation 전까지 clinical threshold로 주장하면 안 된다.",
    ]
    (OUT / "CROSS_Neo_decision_storyboard_overview_KR.md").write_text("\n".join([x for x in lines if x is not None]) + "\n")

    (OUT / "storyboard_summary.json").write_text(
        json.dumps(
            {
                "n_data_artifacts": int(len(manifest)),
                "n_figure_table_items": int(len(figtab)),
                "n_representative_candidates": int(len(candidates)),
                "figures": [
                    "fig_md24_story_overview_flow",
                    "fig_md25_data_usage_matrix",
                    "fig_md26_representative_candidate_board",
                ],
            },
            indent=2,
        )
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    manifest = build_data_usage_manifest()
    figtab = build_figure_table_manifest()
    candidates = build_representative_candidate_table()
    stage = read_tsv(MD_OUT / "dl_first_funnel_ui/dl_first_stage_label_counts.tsv")
    presets = read_tsv(MD_OUT / "dl_threshold_optimization/threshold_recommended_presets.tsv")
    make_overview_figure(stage, presets)
    make_data_usage_matrix(manifest)
    make_candidate_board(candidates)
    write_reports(manifest, figtab, candidates, stage, presets)
    print("[storyboard]", (OUT / "storyboard_summary.json").read_text())


if __name__ == "__main__":
    main()
