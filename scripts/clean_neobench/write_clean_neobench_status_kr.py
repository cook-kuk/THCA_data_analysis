#!/usr/bin/env python3
"""Write a compact Korean current-status report for CLEAN-NeoBench/BAR-Neo."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from common import dataframe_to_markdown, update_manifest


OUTPUT = "CLEAN_NEOBENCH_CURRENT_STATUS_KR.md"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def fmt(value: Any, digits: int = 3) -> str:
    try:
        if pd.isna(value):
            return "NA"
        if isinstance(value, bool):
            return str(value)
        x = float(value)
        if x.is_integer():
            return f"{int(x):,}"
        if abs(x) >= 1000:
            return f"{x:,.0f}"
        return f"{x:.{digits}f}"
    except Exception:
        return str(value)


def metric(summary: dict[str, Any], key: str) -> str:
    return fmt(summary.get(key, "NA"))


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root

    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    summary = manifest.get("summary", {}) if isinstance(manifest.get("summary", {}), dict) else {}

    winloss = read_tsv(output_root / "clean_neobench_winloss_method_summary.tsv")
    leaderboard = read_tsv(output_root / "clean_neobench_leaderboard.tsv")
    contextual = read_tsv(output_root / "barneo_contextual_bma_context_summary.tsv")
    challenge = read_tsv(output_root / "clean_neobench_challenge_axis_summary.tsv")
    patient_req = read_tsv(output_root / "patient_gated_clean_neo_metadata_requirements.tsv")
    public_audit = read_tsv(output_root / "clean_neobench_public_tool_overlap_audit.tsv")

    clean_winloss = winloss[
        winloss.get("method_role", pd.Series(dtype=str)).isin(["anchor", "internal_candidate"])
    ].copy() if not winloss.empty else winloss
    qk_winloss = winloss[
        winloss.get("method_role", pd.Series(dtype=str)).eq("bounded_fallback")
    ].copy() if not winloss.empty else winloss
    public_winloss = winloss[
        winloss.get("method_role", pd.Series(dtype=str)).eq("caveated_public_comparator")
    ].copy() if not winloss.empty else winloss

    blocked = patient_req[
        patient_req.get("metadata_status", pd.Series(dtype=str)).isin(["absent_from_master", "present_but_empty"])
    ].copy() if not patient_req.empty else patient_req
    public_unresolved = public_audit[
        public_audit.get("training_overlap_audit_status", pd.Series(dtype=str)).astype(str).str.contains("unresolved|missing", case=False, na=False)
    ].copy() if not public_audit.empty else public_audit

    text = f"""# CLEAN-NeoBench / BAR-Neo Current Status KR

## 한 줄 결론

현재 위치는 `new SOTA predictor`가 아니라 **leakage-aware AI neoantigen benchmark + benchmark-adaptive reliability/abstention controller**다. 실용적으로는 agent가 여러 expert/ensemble을 불러오고, BAR-Neo-BMA가 benchmark 성능과 실패 패턴으로 가중치를 조절한 뒤, claim-safe layer가 위험한 후보를 abstain/manual review로 보낸다.

## 현재 run 숫자

| 항목 | 값 |
|---|---:|
| candidates | {metric(summary, "n_candidates")} |
| methods | {metric(summary, "n_methods")} |
| BMA weighted methods | {metric(summary, "n_bma_methods_weighted")} |
| split metric rows | {metric(summary, "n_metric_rows")} |
| BAR-Neo-BMA candidates | {metric(summary, "n_bma_candidates")} |
| BAR-Neo-BMA abstain | {metric(summary, "n_bma_abstain")} |
| BAR-Neo-BMA non-abstain | {metric(summary, "n_bma_nonabstain")} |
| contextual clean claims allowed | {metric(summary, "n_contextual_bma_clean_claim_allowed")} |
| BAR-Neo-X priority review candidates | {metric(summary, "n_barneo_x_priority_review_candidates")} |
| BAR-Neo-X top claim-safe score | {metric(summary, "barneo_x_top_claim_safe_score")} |
| challenge rows | {metric(summary, "n_challenge_pack_rows")} |
| visual dashboard figures | {metric(summary, "n_visual_dashboard_figures")} |
| public clean comparators allowed | {metric(summary, "n_public_clean_comparators_allowed")} |

## 우리 알고리즘 위치

- `Structure_LR`: honest local anchor.
- `W7B_stacked`, `W7A_full`, source-balanced/PU RF: 현재 clean internal 후보군.
- `BAR-Neo`: 후보 intrinsic + method rank + reliability + uncertainty + leakage/source/HLA context로 candidate reliability score와 abstention reason을 낸다.
- `BAR-Neo-BMA`: 실용 agent ensemble controller. expert utility, top-k, calibration, source collapse, role prior, public caveat penalty, QK bounded fallback penalty를 가중치에 반영한다.
- `BAR-Neo-X`: leakage/claim-safe reranking과 설명 layer. 지금 clean claim은 막고 priority review만 만든다.
- QK 계열: 성능 좋은 slice가 있어도 **bounded fallback/fusion component**다. headline 또는 quantum advantage claim이 아니다.

## Structure_LR 대비 승패 핵심표

{dataframe_to_markdown(winloss[["method_name", "method_role", "method_disposition", "n_matched_anchor_splits", "n_wins_vs_anchor", "n_losses_vs_anchor", "n_ties_vs_anchor", "win_rate_vs_anchor", "median_delta_AUPRC_vs_anchor", "source_heldout_median_delta_AUPRC", "hla_heldout_median_delta_AUPRC"]].head(20) if not winloss.empty else winloss, max_rows=20)}

## Clean internal 후보만 보면

{dataframe_to_markdown(clean_winloss[["method_name", "n_wins_vs_anchor", "n_losses_vs_anchor", "win_rate_vs_anchor", "median_delta_AUPRC_vs_anchor", "source_heldout_median_delta_AUPRC", "hla_heldout_median_delta_AUPRC"]].head(15) if not clean_winloss.empty else clean_winloss, max_rows=15)}

## QK / fallback 계열 판단

{dataframe_to_markdown(qk_winloss[["method_name", "n_wins_vs_anchor", "n_losses_vs_anchor", "win_rate_vs_anchor", "median_delta_AUPRC_vs_anchor", "source_heldout_median_delta_AUPRC", "hla_heldout_median_delta_AUPRC"]].head(12) if not qk_winloss.empty else qk_winloss, max_rows=12)}

판단: `W7A_QK_only`는 여러 matched split에서 Structure_LR을 이기지만, 역할은 fallback/fusion이다. 보고서에서는 이 결과를 ranking headline이나 quantum claim으로 쓰지 않는다.

## Public pretrained tool 판단

{dataframe_to_markdown(public_winloss[["method_name", "n_wins_vs_anchor", "n_losses_vs_anchor", "win_rate_vs_anchor", "median_delta_AUPRC_vs_anchor", "method_disposition"]].head(12) if not public_winloss.empty else public_winloss, max_rows=12)}

{dataframe_to_markdown(public_unresolved[["method_name", "training_overlap_audit_status", "clean_comparator_allowed_after_audit", "reviewer_disposition"]].head(12) if not public_unresolved.empty else public_unresolved, max_rows=12)}

판단: public tool이 일부 split에서 좋아도 row-level training-corpus overlap audit 전까지는 caveated comparator다.

## 어디서 틀리는지

- source shift: 일부 RF/PU/source-balanced 계열은 overall AUPRC가 높지만 source-heldout median delta가 크게 음수다.
- low prevalence: TESLA-like 저 prevalence slice에서 false-positive pressure가 커진다.
- rare / underrepresented HLA: HLA-heldout 또는 Korean-HLA focus에서 method별 collapse가 갈린다.
- public/internal disagreement: public pretrained predictor가 높고 clean internal support가 약하면 manual review 또는 abstain.
- patient gate: PAAD/THCA patient-level fields가 부족해서 현재 patient-gated demo는 clinical-use가 아니라 research triage only다.

## Challenge / distribution audit

{dataframe_to_markdown(challenge[["challenge_axis", "n_unique_candidates", "positive_prevalence", "recommended_split_contract", "success_metric"]].head(20) if not challenge.empty else challenge, max_rows=20)}

## Patient-gated blocker

{dataframe_to_markdown(blocked[["field", "gate_group", "required_for", "metadata_status", "reviewer_safe_default_if_missing"]].head(25) if not blocked.empty else blocked, max_rows=25)}

## 바로 볼 파일

- HTML dossier: `project/papers_hub_2026_05_04/clean_neobench_barneo_dossier_2026_05_09.html`
- Visual dashboard: `project/papers_hub_2026_05_04/clean_neobench_visual_dashboard_2026_05_10.html`
- Win/loss report: `project/results/clean_neobench_barneo_2026_05_09/CLEAN_NEOBENCH_WINLOSS_REPORT_KR.md`
- Win/loss table: `project/results/clean_neobench_barneo_2026_05_09/clean_neobench_winloss_method_summary.tsv`
- Split winner board: `project/results/clean_neobench_barneo_2026_05_09/clean_neobench_split_winner_board.tsv`
- Distribution error audit: `project/results/clean_neobench_barneo_2026_05_09/CLEAN_NEOBENCH_DISTRIBUTION_ERROR_AUDIT_KR.md`
- Challenge pack: `project/results/clean_neobench_barneo_2026_05_09/CLEAN_NEOBENCH_CHALLENGE_PACK_KR.md`

## 다음 우선순위

1. Public tool training-corpus row-level overlap audit.
2. Source-heldout collapse가 큰 RF/PU/source-balanced 모델의 source reweighting 또는 abstention 강화.
3. Rare/Korean-HLA focus split에서 W7B/W7A/source-balanced RF를 다시 stress-test.
4. PAAD/THCA patient-gated demo에 실제 disease timing, expression/clonality, immune context, safety metadata 연결.
5. MHC-II는 별도 benchmark로 분리.

## Claim boundary

Allowed: leakage-aware AI neoantigen predictor benchmarking framework, benchmark-adaptive reliability ranking, calibrated candidate prioritization with abstention, reviewer-safe comparison, research triage framework.

Forbidden: clinical vaccine selection, new SOTA predictor, external validation proven, quantum advantage, public tools as clean baselines without overlap audit, Class I and Class II unified predictor.
"""

    path = output_root / OUTPUT
    path.write_text(text.strip() + "\n")
    update_manifest(
        output_root,
        "clean_neobench_current_status_kr",
        {
            "outputs": [OUTPUT],
            "warnings": ["Current status report is a reviewer-safe project snapshot, not a SOTA or clinical claim."],
        },
    )
    print(f"[clean-neobench-status-kr] wrote {path}")


if __name__ == "__main__":
    main()
