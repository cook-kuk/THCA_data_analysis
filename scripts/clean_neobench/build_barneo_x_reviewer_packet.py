#!/usr/bin/env python3
"""Build a portable BAR-Neo-X reviewer packet."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import zipfile
from pathlib import Path
from typing import Any


PAGE_NAME = "cancer_vaccine_barneo_x_interpretability_2026_05_09.html"
PACKET_NAME = "barneo_x_reviewer_packet_2026_05_09"
KAKAO_TEMPLATE = (
    "Cancer vaccine neoantigen benchmark에서 BAR-Neo 기반으로 BAR-Neo-X를 만들었고, "
    "후보별로 BAR-Neo evidence/BMA consensus/internal predictor/confidence와 "
    "leakage-overlap-uncertainty penalty를 분해해서 설명 가능하게 만들었습니다. "
    "raw BAR-Neo는 apparent AUPRC {raw_auprc}지만 top10이 전부 high-leakage라 SOTA 주장엔 위험하고, "
    "BAR-Neo-X claim-safe score는 AUPRC {claim_auprc}으로 보수화되는 대신 top10 high-leakage를 "
    "{raw_top10_leak}->{claim_top10_leak}로 제거하면서 top10 precision {claim_top10_precision}를 유지합니다. "
    "ablation에서도 positive-only는 top10 leakage {positive_top10_leak}라 위험하고, leakage gate 제거 시 top20 leakage가 "
    "{no_leak_gate_top20_leak}까지 올라갑니다. 안 될 때도 이유가 분해됩니다: {hard_blocked}개는 leakage/identity overlap으로 "
    "clean claim 불가, {metadata_rescuable}개는 patient/disease metadata 보강 시 구제 가능, "
    "{priority_if_metadata}개는 metadata 완성 시 priority threshold 통과 가능, 현재 priority-now는 {priority_now}개입니다. "
    "즉 무조건 SOTA/clinical/quantum advantage가 아니라, "
    "리뷰어 방어 가능한 해석형-누수차단 neoantigen triage layer입니다."
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench BAR-Neo output directory")
    parser.add_argument("--hub-root", default="project/papers_hub_2026_05_04", help="HTML hub directory")
    parser.add_argument("--page-name", default=PAGE_NAME)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def metric(rows: list[dict[str, str]], score: str, field: str) -> str:
    for row in rows:
        if row.get("score") == score:
            return fmt(row.get(field, "NA"))
    return "NA"


def fmt(value: object) -> str:
    try:
        return f"{float(value):.3f}"
    except Exception:
        return str(value)


def copy_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return True


def write_file_manifest(packet_dir: Path) -> None:
    rows = []
    for path in sorted(p for p in packet_dir.rglob("*") if p.is_file()):
        rows.append(
            {
                "relative_path": str(path.relative_to(packet_dir)),
                "bytes": str(path.stat().st_size),
            }
        )
    out = packet_dir / "FILE_MANIFEST.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "bytes"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def update_manifest(output_root: Path, packet_zip: Path, packet_dir: Path) -> None:
    path = output_root / "run_manifest.json"
    manifest = read_json(path)
    manifest.setdefault("output_files", [])
    for name in [packet_zip.name, packet_dir.name]:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("stages", {})
    manifest["stages"]["barneo_x_reviewer_packet"] = {
        "outputs": [packet_dir.name, packet_zip.name],
        "warnings": [
            "Reviewer packet is a portable evidence bundle, not a clinical decision package.",
            "BAR-Neo-X remains a claim-safe triage layer, not a public SOTA or quantum-advantage claim.",
        ],
    }
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def zip_dir(packet_dir: Path, zip_path: Path) -> None:
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(p for p in packet_dir.rglob("*") if p.is_file()):
            zf.write(path, path.relative_to(packet_dir.parent))


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    hub_root = Path(args.hub_root)
    packet_dir = output_root / PACKET_NAME
    packet_dir.mkdir(parents=True, exist_ok=True)
    stale_run_manifest = packet_dir / "run_manifest.json"
    if stale_run_manifest.exists():
        stale_run_manifest.unlink()

    metrics = read_tsv(output_root / "barneo_x_metric_audit.tsv")
    ablation = read_tsv(output_root / "barneo_x_ablation_audit.tsv")
    summary = read_json(output_root / "BAR_NEO_X_INTERPRETABILITY_BOOST_SUMMARY.json")
    why_summary = read_json(output_root / "BAR_NEO_X_WHY_NOT_SUMMARY.json")

    claim_card = {
        "project": "CLEAN-NeoBench BAR-Neo-X",
        "date": "2026-05-09",
        "allowed_claim": "Explainable, leakage-aware neoantigen research triage for reviewer-facing candidate review.",
        "forbidden_claims": [
            "public SOTA",
            "clinical vaccine selection",
            "quantum advantage",
            "clean public baseline without row-level overlap audit",
        ],
        "candidate_count": summary.get("n_candidates"),
        "component_attribution_rows": summary.get("n_component_attribution_rows"),
        "raw_barneo": {
            "apparent_auprc": metric(metrics, "barneo_score", "apparent_auprc"),
            "top10_precision": metric(metrics, "barneo_score", "top10_precision"),
            "top10_high_leakage_fraction": metric(metrics, "barneo_score", "top10_high_leakage_fraction"),
        },
        "barneo_x_claim_safe": {
            "apparent_auprc": metric(metrics, "barneo_x_claim_safe_score", "apparent_auprc"),
            "top10_precision": metric(metrics, "barneo_x_claim_safe_score", "top10_precision"),
            "top10_high_leakage_fraction": metric(metrics, "barneo_x_claim_safe_score", "top10_high_leakage_fraction"),
        },
        "ablation_boundary": {
            "positive_only_top10_high_leakage_fraction": metric(ablation, "barneo_x_positive_only_score", "top10_high_leakage_fraction"),
            "no_leakage_gate_top20_high_leakage_fraction": metric(ablation, "barneo_x_no_leakage_gate_score", "top20_high_leakage_fraction"),
        },
        "why_not_audit": {
            "hard_claim_blocked_by_leakage_or_identity_overlap": why_summary.get("n_hard_claim_blocked"),
            "metadata_rescuable": why_summary.get("n_metadata_rescuable"),
            "priority_review_now": why_summary.get("n_priority_now"),
            "would_be_priority_if_metadata_complete": why_summary.get("n_would_be_priority_if_metadata_complete"),
        },
    }

    kakao = KAKAO_TEMPLATE.format(
        raw_auprc=claim_card["raw_barneo"]["apparent_auprc"],
        claim_auprc=claim_card["barneo_x_claim_safe"]["apparent_auprc"],
        raw_top10_leak=claim_card["raw_barneo"]["top10_high_leakage_fraction"],
        claim_top10_leak=claim_card["barneo_x_claim_safe"]["top10_high_leakage_fraction"],
        claim_top10_precision=claim_card["barneo_x_claim_safe"]["top10_precision"],
        positive_top10_leak=claim_card["ablation_boundary"]["positive_only_top10_high_leakage_fraction"],
        no_leak_gate_top20_leak=claim_card["ablation_boundary"]["no_leakage_gate_top20_high_leakage_fraction"],
        hard_blocked=claim_card["why_not_audit"]["hard_claim_blocked_by_leakage_or_identity_overlap"],
        metadata_rescuable=claim_card["why_not_audit"]["metadata_rescuable"],
        priority_if_metadata=claim_card["why_not_audit"]["would_be_priority_if_metadata_complete"],
        priority_now=claim_card["why_not_audit"]["priority_review_now"],
    )

    (packet_dir / "BAR_NEO_X_CLAIM_CARD.json").write_text(json.dumps(claim_card, indent=2, sort_keys=True) + "\n")
    (packet_dir / "BAR_NEO_X_KAKAO_ONE_SHOT.txt").write_text(kakao + "\n")
    (packet_dir / "README.md").write_text(
        "\n".join(
            [
                "# BAR-Neo-X Reviewer Packet",
                "",
                "Portable evidence bundle for BAR-Neo-X claim-safe interpretability.",
                "",
                "Allowed claim: explainable, leakage-aware neoantigen research triage.",
                "",
                "Forbidden claims: public SOTA, clinical vaccine selection, quantum advantage, or clean public baseline without row-level overlap audit.",
                "",
                "Key boundary: BAR-Neo-X deliberately trades apparent AUPRC for top-k leakage control.",
                "",
                f"Raw BAR-Neo AUPRC: {claim_card['raw_barneo']['apparent_auprc']}; raw top10 high-leakage: {claim_card['raw_barneo']['top10_high_leakage_fraction']}.",
                f"BAR-Neo-X claim-safe AUPRC: {claim_card['barneo_x_claim_safe']['apparent_auprc']}; claim-safe top10 high-leakage: {claim_card['barneo_x_claim_safe']['top10_high_leakage_fraction']}.",
                f"Why-not audit: {claim_card['why_not_audit']['hard_claim_blocked_by_leakage_or_identity_overlap']} hard claim-blocked; {claim_card['why_not_audit']['metadata_rescuable']} metadata-rescuable; {claim_card['why_not_audit']['would_be_priority_if_metadata_complete']} would be priority if metadata were complete; {claim_card['why_not_audit']['priority_review_now']} priority now.",
                "",
            ]
        )
    )

    zip_path = output_root / f"{PACKET_NAME}.zip"
    copies = [
        output_root / "BAR_NEO_X_INTERPRETABILITY_BOOST_REPORT.md",
        output_root / "BAR_NEO_X_INTERPRETABILITY_BOOST_SUMMARY.json",
        output_root / "BAR_NEO_X_WHY_NOT_REPORT.md",
        output_root / "BAR_NEO_X_WHY_NOT_SUMMARY.json",
        output_root / "barneo_x_metric_audit.tsv",
        output_root / "barneo_x_topk_safety_audit.tsv",
        output_root / "barneo_x_ablation_audit.tsv",
        output_root / "barneo_x_candidate_scores.tsv",
        output_root / "barneo_x_candidate_explanations.tsv",
        output_root / "barneo_x_component_attributions.tsv",
        output_root / "barneo_x_why_not_audit.tsv",
        output_root / "barneo_x_why_not_reason_summary.tsv",
        output_root / "barneo_x_primary_blocker_summary.tsv",
        output_root / "barneo_x_rescue_lane_summary.tsv",
        hub_root / args.page_name,
        hub_root / "assets" / "barneo_x" / "barneo_x_tradeoff.png",
    ]
    for src in copies:
        copy_if_exists(src, packet_dir / src.name)

    write_file_manifest(packet_dir)
    zip_dir(packet_dir, zip_path)
    hub_zip = hub_root / "assets" / "barneo_x" / zip_path.name
    copy_if_exists(zip_path, hub_zip)

    print(json.dumps({"packet_dir": str(packet_dir), "zip": str(zip_path), "hub_zip": str(hub_zip)}, indent=2))


if __name__ == "__main__":
    main()
