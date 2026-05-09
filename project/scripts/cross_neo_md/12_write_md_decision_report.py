#!/usr/bin/env python3
"""Write final decision reports for the CROSS-Neo MD audit layer."""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
CROSS = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"


def read_tsv(path: Path) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame()


def fmt(x, nd=3) -> str:
    if pd.isna(x):
        return "NA"
    try:
        return f"{float(x):.{nd}f}"
    except Exception:
        return str(x)


def md_table(df: pd.DataFrame, cols: list[str], max_rows: int = 12) -> str:
    if df.empty:
        return "No rows available."
    keep = [c for c in cols if c in df.columns]
    if not keep:
        return "No requested columns available."
    return df[keep].head(max_rows).to_markdown(index=False)


def candidate_metrics() -> dict[str, dict[str, object]]:
    status = read_tsv(OUT / "md_status_summary.tsv")
    scores = read_tsv(OUT / "md_evidence_scores.tsv")
    qc = read_tsv(OUT / "qc/md_qc_summary.tsv")
    pmhc = read_tsv(OUT / "contacts/pmhc_anchor_contacts.tsv")
    tcr = read_tsv(OUT / "contacts/tcr_peptide_contact_timeseries.tsv")

    out: dict[str, dict[str, object]] = {}
    for candidate in sorted(set(status.get("candidate", pd.Series(dtype=str)).dropna()) | set(scores.get("candidate", pd.Series(dtype=str)).dropna())):
        s = status[status.get("candidate", "") == candidate] if not status.empty else pd.DataFrame()
        sc = scores[scores.get("candidate", "") == candidate] if not scores.empty else pd.DataFrame()
        q = qc[qc.get("candidate", "") == candidate] if not qc.empty else pd.DataFrame()
        p = pmhc[pmhc.get("candidate", "") == candidate] if not pmhc.empty else pd.DataFrame()
        tc = tcr[tcr.get("candidate", "") == candidate] if not tcr.empty else pd.DataFrame()

        primary = sc[sc.get("condition", "").astype(str).str.contains("explicit_cuda_10ns", na=False)] if not sc.empty else pd.DataFrame()
        if primary.empty and not sc.empty:
            primary = sc.sort_values("runtime_fraction", ascending=False).head(1)
        row = primary.iloc[0].to_dict() if not primary.empty else {}
        if "time_ps" in s.columns and not s.empty:
            s_for_live = s.copy()
            s_for_live["time_ps_numeric"] = pd.to_numeric(s_for_live["time_ps"], errors="coerce")
            live_row = s_for_live.sort_values("time_ps_numeric", ascending=False).iloc[0].to_dict()
        else:
            live_row = {}
        out[candidate] = {
            "max_time_ps": live_row.get("time_ps", pd.NA),
            "temperature_k": live_row.get("temperature_k", pd.NA),
            "speed_ns_day": live_row.get("speed_ns_per_day", pd.NA),
            "status_rows": len(s),
            "md_label": row.get("MD_evidence_label", "NA"),
            "md_score": row.get("MD_evidence_score", pd.NA),
            "pmhc_score": row.get("pMHC_stability_score", pd.NA),
            "tcr_score": row.get("TCR_recognition_score", pd.NA),
            "runtime_fraction": row.get("runtime_fraction", pd.NA),
            "peptide_rmsd_final_nm": q["peptide_rmsd_final_nm"].dropna().iloc[0] if "peptide_rmsd_final_nm" in q.columns and q["peptide_rmsd_final_nm"].notna().any() else row.get("peptide_rmsd_final_nm", pd.NA),
            "peptide_drift_final_nm": q["peptide_com_drift_final_nm"].dropna().iloc[0] if "peptide_com_drift_final_nm" in q.columns and q["peptide_com_drift_final_nm"].notna().any() else pd.NA,
            "anchor_rows": len(p),
            "anchor_max_contact": p["max_contact_occupancy"].max() if "max_contact_occupancy" in p.columns and not p.empty else pd.NA,
            "tcr_contacts_tail": row.get("tcr_peptide_contacts_tail", pd.NA),
            "tcr_timeseries_rows": len(tc),
        }
    return out


def write_next_batch() -> pd.DataFrame:
    queue = read_tsv(CROSS / "tcr_extension/md_escalation/md_escalation_queue_top20.tsv")
    if queue.empty:
        rows = []
    else:
        use = queue.head(12).copy()
        rows = []
        for _, r in use.iterrows():
            tier = r.get("md_tier", "P1")
            peptide = r.get("peptide", "")
            hla = r.get("hla_4digit", "")
            if peptide == "GADGVGKSAL" and hla == "HLA-C*08:02":
                action = "P0: add WT, scrambled decoy, same-HLA positive control; run 3x10ns before any 50-100ns promotion"
                reason = "10ns mutant TCR-pMHC completed with moderate structural support, but no counterfactual control exists"
            elif peptide == "HMTEVVRHC" and hla == "HLA-A*02:01":
                action = "P0: let current 10ns finish, sync DCD, rerun full analysis, then launch WT/decoy 3x10ns"
                reason = "highest TCR evidence and live state is stable, but trajectory-derived contacts are pending"
            elif "P1_MD_pMHC_bulge" in str(tier):
                action = "P1: pMHC-only 3x10ns bulge/stability screen; no TCR-recognition claim"
                reason = "long class-I peptide or conformation uncertainty without paired TCR"
            else:
                action = "P1: pMHC 3x10ns plus TCR search; escalate only after TCR/template evidence"
                reason = "candidate is model-fragile or recognition-uncertain"
            rows.append(
                {
                    "priority_rank": len(rows) + 1,
                    "row_id": r.get("row_id", ""),
                    "peptide": peptide,
                    "hla_4digit": hla,
                    "label_binary": r.get("label_binary", ""),
                    "md_tier": tier,
                    "recommended_next_batch": action,
                    "reason": reason,
                    "md_escalation_score": r.get("md_escalation_score", ""),
                    "paired_tcr_evidence_count": r.get("paired_tcr_evidence_count", ""),
                    "tcr_evidence_count": r.get("tcr_evidence_count", ""),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "next_simulation_batch_recommendation.tsv", sep="\t", index=False)
    return df


def command_block() -> str:
    return """```bash
# Update live status and regenerate all currently possible MD outputs.
python project/scripts/cross_neo_md/00_parse_openmm_status.py
python project/scripts/cross_neo_md/01_annotate_complex.py
python project/scripts/cross_neo_md/02_md_qc.py
python project/scripts/cross_neo_md/03_analyze_pmhc_contacts.py
python project/scripts/cross_neo_md/04_analyze_tcr_contacts.py
python project/scripts/cross_neo_md/05_counterfactual_md_analysis.py
python project/scripts/cross_neo_md/06_replicate_consistency.py
python project/scripts/cross_neo_md/07_md_evidence_score.py
python project/scripts/cross_neo_md/08_integrate_md_with_cross_neo.py
python project/scripts/cross_neo_md/09_generate_md_figures.py
python project/scripts/cross_neo_md/11_generate_md_extra_visuals.py
python project/scripts/cross_neo_md/10_build_md_visual_dossier.py
python project/scripts/cross_neo_md/12_write_md_decision_report.py

# While 6VRN is running, sync state only.
OUT=project/results/cross_neo_md_audit_2026_05_10/remote_sync/pod1_thca_neo_bayesian_aux
mkdir -p "$OUT"
ssh -i ~/.runpod/ssh/RunPod-Key-Go -p 20878 root@135.84.176.142 \\
  'cd /workspace/openmm_pilot_10ns_package && tar --exclude=trajectory.dcd --exclude=final.chk -czf - prod_10ns_6VRN_1fs300K' \\
  | tar -xzf - -C "$OUT"

# After 6VRN finishes, sync full trajectory and rerun the analysis stack above.
rsync -avP -e 'ssh -i ~/.runpod/ssh/RunPod-Key-Go -p 20878' \\
  root@135.84.176.142:/workspace/openmm_pilot_10ns_package/prod_10ns_6VRN_1fs300K/ \\
  project/results/cross_neo_md_audit_2026_05_10/remote_sync/pod1_thca_neo_bayesian_aux/prod_10ns_6VRN_1fs300K/
```"""


def write_reports(next_batch: pd.DataFrame) -> None:
    metrics = candidate_metrics()
    status = read_tsv(OUT / "md_status_summary.tsv")
    qc = read_tsv(OUT / "qc/md_qc_summary.tsv")
    scores = read_tsv(OUT / "md_evidence_scores.tsv")
    anchor = read_tsv(OUT / "contacts/pmhc_anchor_contacts.tsv")
    tcr = read_tsv(OUT / "contacts/tcr_peptide_contact_timeseries.tsv")
    integ = read_tsv(OUT / "integrated/cross_neo_md_joined.tsv")
    supported = read_tsv(OUT / "integrated/md_supported_top_candidates.tsv")
    high_low = read_tsv(OUT / "integrated/model_high_md_low_cases.tsv")
    low_high = read_tsv(OUT / "integrated/model_low_md_high_cases.tsv")

    gad = metrics.get("GADGVGKSAL/HLA-C*08:02", {})
    hm = metrics.get("HMTEVVRHC/HLA-A*02:01", {})

    report = [
        "# CROSS-Neo MD Audit Decision Report",
        "",
        "## 1. Executive Verdict",
        "",
        "- `GADGVGKSAL / HLA-C*08:02`: **PROMOTE_TO_WETLAB_PRIORITIZATION_TOOL**, with MD evidence label `MD_MODERATE`. The completed 10 ns run supports structural plausibility, not immunogenicity proof.",
        "- `HMTEVVRHC / HLA-A*02:01`: **HOLD_AS_RUNNING_P0_MD_AUDIT** until the 10 ns DCD is synced and contact/QC analysis is rerun. Live state is stable so far, but trajectory-derived claims are not ready.",
        "- Overall: keep MD as a **structure-aware diagnostic audit layer**. Do not move it into the main CROSS-Neo claim without WT/decoy controls, replicates, and external validation.",
        "",
        "## 2. Which Simulations Completed",
        "",
        md_table(status, ["run_id", "candidate", "condition", "time_ps", "target_ns", "temperature_k", "speed_ns_per_day"]),
        "",
        "Completed primary evidence: `GADGVGKSAL / HLA-C*08:02` reached 10,000 ps. `HMTEVVRHC / HLA-A*02:01` was still partial at the latest sync, at "
        + f"{fmt(hm.get('max_time_ps'), 1)} ps.",
        "",
        "## 3. Which Simulations Are Still Partial",
        "",
        f"- HMTEVVRHC primary 10 ns OpenMM CUDA run: latest synced time {fmt(hm.get('max_time_ps'), 1)} ps, temperature {fmt(hm.get('temperature_k'), 2)} K, speed {fmt(hm.get('speed_ns_day'), 1)} ns/day.",
        "- HMTEVVRHC alternate replicate watcher should start after the primary exits; do not score replicate consistency until those DCDs are present.",
        "",
        "## 4. QC Status",
        "",
        md_table(qc, ["run_id", "candidate", "condition", "n_frames", "total_time_ps", "peptide_rmsd_final_nm", "peptide_com_drift_final_nm", "qc_status"]),
        "",
        f"GADGVGKSAL primary peptide RMSD final: {fmt(gad.get('peptide_rmsd_final_nm'), 3)} nm. Peptide COM drift final: {fmt(gad.get('peptide_drift_final_nm'), 3)} nm. This is compatible with a stable pilot, not a definitive binding claim.",
        "",
        "## 5. pMHC Stability Evidence",
        "",
        md_table(anchor, ["run_id", "candidate", "peptide_position", "peptide_resname", "is_anchor", "max_contact_occupancy", "sum_contact_occupancy"], 16),
        "",
        "GADGVGKSAL anchor contacts remained high in the completed primary and 1 ns replicate screens. HMTEVVRHC anchor evidence currently comes only from smoke/partial status until full trajectory analysis is available.",
        "",
        "## 6. TCR-Peptide Recognition Evidence",
        "",
        md_table(scores, ["run_id", "candidate", "condition", "TCR_recognition_score", "tcr_peptide_contacts_tail", "MD_evidence_label"]),
        "",
        "The GADGVGKSAL primary run has persistent TCR-peptide contact signal in the audit tables, but this remains structural compatibility evidence. It does not prove T-cell activation.",
        "",
        "## 7. Mutant-WT Or Decoy Counterfactual Evidence",
        "",
        "No WT, scrambled decoy, or same-HLA positive-control trajectory has been discovered yet. Mutant-specific recognition, specificity, and WT cross-reactivity claims are therefore blocked.",
        "",
        "## 8. Replicate Consistency",
        "",
        "GADGVGKSAL has one completed 1 ns alternate screen and one 10 ns primary run. The replicate is useful as an early instability filter only. A serious decision needs 3x10 ns for mutant, WT, and scrambled controls before any 50-100 ns escalation.",
        "",
        "## 9. Integration With CROSS-Neo Predictions",
        "",
        f"- Joined prediction rows with MD candidate keys: {len(integ)}",
        f"- MD-supported model-high rows: {len(supported)}",
        f"- model-high but MD-low-or-insufficient rows: {len(high_low)}",
        f"- model-low but MD-high rows: {len(low_high)}",
        "",
        "Top MD-supported rows:",
        "",
        md_table(supported, ["row_id", "source_dataset", "peptide", "hla_4digit", "label", "score", "model_name", "split_name", "MD_evidence_label", "MD_evidence_score", "wetlab_priority_score_evidence_adjusted"], 8),
        "",
        "## 10. Allowed Claims",
        "",
        "- Explicit-solvent MD was used as a structural audit layer for selected CROSS-Neo/TCR-aware candidates.",
        "- GADGVGKSAL / HLA-C*08:02 showed moderate structural audit support in a completed 10 ns pilot, including peptide-MHC stability and persistent TCR-peptide contact signal.",
        "- MD can help prioritize wetlab candidates and explain model/model-expert disagreements.",
        "",
        "## 11. Forbidden Claims",
        "",
        "- Do not claim MD proves immunogenicity, clinical efficacy, or SOTA.",
        "- Do not claim one short trajectory proves binding.",
        "- Do not claim mutant-specific recognition without WT/decoy comparison.",
        "- Do not treat AlphaFold/TCRdock/template structures as ground truth.",
        "- Do not use peptide-only TCR evidence as definitive neoantigen ground truth.",
        "",
        "## 12. Wetlab Candidate Recommendation",
        "",
        "- Immediate wetlab-priority candidate: `GADGVGKSAL / HLA-C*08:02`, as structural audit support now agrees with several CROSS-Neo/TCR-aware model rows. Use as prioritization evidence only.",
        "- Keep `HMTEVVRHC / HLA-A*02:01` as the highest-priority pending MD candidate because it has strong TCR evidence and an active 10 ns run, but wait for completed trajectory contacts before moving it above GADGVGKSAL on MD evidence.",
        "- Do not advance candidates lacking paired TCR/template support to TCR-pMHC MD claims; use pMHC-only stability screens for those.",
        "",
        "## 13. Next Simulation Batch Recommendation",
        "",
        md_table(next_batch, ["priority_rank", "row_id", "peptide", "hla_4digit", "md_tier", "recommended_next_batch", "reason"], 12),
        "",
        "Minimum serious batch per candidate: mutant pMHC, WT pMHC, mutant TCR-pMHC if paired TCR/template exists, WT TCR-pMHC, scrambled peptide negative control, same/similar-HLA positive control, 3 replicates each, 10 ns each. Promote only the strongest and most consistent candidates to 50-100 ns.",
        "",
        "## Exact Rerun Commands",
        "",
        command_block(),
    ]
    main_report = OUT / "CROSS_Neo_MD_audit_decision_report.md"
    main_report.write_text("\n".join(report) + "\n")
    shutil.copyfile(main_report, REPO / "CROSS_Neo_MD_audit_decision_report.md")

    paper_summary = [
        "# One-Page Paper Summary: CROSS-Neo MD Audit Layer",
        "",
        "We added an explicit-solvent OpenMM MD audit layer for selected CROSS-Neo/TCR-aware neoantigen candidates. The layer is designed to test structural plausibility and diagnostic consistency, not to replace the main pMHC ranking model and not to prove immunogenicity.",
        "",
        f"`GADGVGKSAL / HLA-C*08:02` completed a 10 ns CUDA explicit-solvent run and one 1 ns alternate screen. The primary run reached an MD evidence label of `MD_MODERATE`, with final peptide RMSD {fmt(gad.get('peptide_rmsd_final_nm'), 3)} nm and persistent TCR-peptide contact signal. This supports wetlab prioritization and case interpretation, while remaining below the threshold for any immunogenicity or clinical claim.",
        "",
        f"`HMTEVVRHC / HLA-A*02:01` remains a high-priority active run. The latest synced state reached {fmt(hm.get('max_time_ps'), 1)} ps at {fmt(hm.get('temperature_k'), 2)} K, but full DCD-based RMSD/contact analysis is pending. Claims for this candidate should wait until the 10 ns trajectory and follow-up replicate screens are parsed.",
        "",
        f"Integration with CROSS-Neo predictions found {len(supported)} MD-supported model-high rows and {len(low_high)} model-low/MD-high disagreement rows. These cases are useful for model audit and wetlab triage. No WT or decoy trajectories are available yet, so mutant-specific recognition and WT cross-reactivity claims remain blocked.",
    ]
    (OUT / "MD_one_page_paper_summary.md").write_text("\n".join(paper_summary) + "\n")

    korean = [
        "# Seungho용 1페이지 요약: CROSS-Neo MD 결과",
        "",
        "현재 결론은 명확합니다. MD는 면역원성 증명이 아니라 구조 audit 레이어입니다. 그래서 논문 메인 claim으로 쓰기보다, wetlab 후보 우선순위와 false positive/false negative 해석에 쓰는 것이 안전합니다.",
        "",
        f"`GADGVGKSAL / HLA-C*08:02`는 10 ns explicit-solvent CUDA run이 완료됐고, MD evidence는 `MD_MODERATE`입니다. peptide RMSD final은 {fmt(gad.get('peptide_rmsd_final_nm'), 3)} nm이고, peptide-MHC anchor와 TCR-peptide contact가 유지되는 쪽으로 보입니다. 이건 wetlab 우선순위 근거로 좋습니다.",
        "",
        f"`HMTEVVRHC / HLA-A*02:01`는 아직 running입니다. 최신 동기화 기준 {fmt(hm.get('max_time_ps'), 1)} ps, 온도 {fmt(hm.get('temperature_k'), 2)} K라서 run 자체는 정상입니다. 하지만 DCD 전체 분석 전에는 contact/RMSD claim을 하지 않는 게 맞습니다.",
        "",
        "다음 액션은 3개입니다. 1) HMTEVVRHC 10 ns 끝나면 DCD sync 후 전체 분석 재실행. 2) GADGVGKSAL과 HMTEVVRHC에 WT/decoy/same-HLA positive control을 붙여 3x10 ns batch 실행. 3) 그 다음에만 50-100 ns 장기 MD로 승격합니다.",
        "",
        "절대 하면 안 되는 말: MD가 면역원성을 증명했다, 임상 효능을 보였다, SOTA를 검증했다, 짧은 trajectory 하나로 binding을 증명했다. 지금 가능한 말은 `구조적으로 그럴듯해서 wetlab 후보로 올릴 가치가 있다`입니다.",
    ]
    (OUT / "MD_one_page_korean_summary_for_Seungho.md").write_text("\n".join(korean) + "\n")

    (OUT / "rerun_md_analysis_commands.md").write_text("# Rerun Commands\n\n" + command_block() + "\n")
    print(f"[md-decision] wrote {main_report}")
    print(f"[md-decision] wrote {OUT / 'MD_one_page_paper_summary.md'}")
    print(f"[md-decision] wrote {OUT / 'MD_one_page_korean_summary_for_Seungho.md'}")


def main() -> None:
    next_batch = write_next_batch()
    write_reports(next_batch)


if __name__ == "__main__":
    main()
