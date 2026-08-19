#!/usr/bin/env python3
"""Build a compact reviewer packet for the CROSS-Neo-TCR extension."""

from __future__ import annotations

import shutil
import zipfile
from pathlib import Path

from common import OUT, REPO


TCR_OUT = OUT / "tcr_extension"
PACKET = TCR_OUT / "reviewer_packet"
ZIP_PATH = TCR_OUT / "cross_neo_tcr_reviewer_packet_2026_05_09.zip"
HUB_ASSET = REPO / "project/papers_hub_2026_05_04/assets/cross_neo_tcr_extension"
WEB_ASSET = Path("/var/www/papers/papers_hub_2026_05_04/assets/cross_neo_tcr_extension")


FILES: list[tuple[Path, str]] = [
    (REPO / "CROSS_Neo_TCR_extension_decision_report.md", "reports/CROSS_Neo_TCR_extension_decision_report.md"),
    (REPO / "output/manuscript_tcr_aware_track.md", "reports/manuscript_tcr_aware_track.md"),
    (TCR_OUT / "00_tcr_lit_tool_audit.md", "reports/00_tcr_lit_tool_audit.md"),
    (TCR_OUT / "tcr_registry_report.md", "reports/tcr_registry_report.md"),
    (TCR_OUT / "tcr_neo_linkage_report.md", "reports/tcr_neo_linkage_report.md"),
    (TCR_OUT / "tcr_feature_report.md", "reports/tcr_feature_report.md"),
    (TCR_OUT / "tcr_structure_pipeline_report.md", "reports/tcr_structure_pipeline_report.md"),
    (TCR_OUT / "model_discovery/01_tcr_model_discovery.md", "reports/01_tcr_model_discovery.md"),
    (TCR_OUT / "model_discovery/02_tcr_model_runtime_audit.md", "reports/02_tcr_model_runtime_audit.md"),
    (TCR_OUT / "model_discovery/external_tcr_expert_benchmark/04_external_tcr_expert_benchmark.md", "reports/04_external_tcr_expert_benchmark.md"),
    (TCR_OUT / "model_discovery/wetlab_external_tcr_expert_scores/05_wetlab_external_tcr_expert_scores.md", "reports/05_wetlab_external_tcr_expert_scores.md"),
    (TCR_OUT / "md_escalation/md_escalation_report.md", "reports/md_escalation_report.md"),
    (TCR_OUT / "md_escalation/md_jobs/md_job_stub_report.md", "reports/md_job_stub_report.md"),
    (TCR_OUT / "md_escalation/p0_structures/p0_structure_pilot_report.md", "reports/p0_structure_pilot_report.md"),
    (TCR_OUT / "md_escalation/p0_structures/openmm_minimized/openmm_p0_minimization_report.md", "reports/openmm_p0_minimization_report.md"),
    (TCR_OUT / "md_escalation/p0_structures/openmm_pilot_10ns_package/p0_openmm_smoke_report.md", "reports/p0_openmm_smoke_report.md"),
    (TCR_OUT / "md_escalation/p0_structures/openmm_pilot_10ns_package/README.md", "reports/openmm_p0_pilot_package_README.md"),
    (TCR_OUT / "tcr_extension_figure_report.md", "reports/tcr_extension_figure_report.md"),
    (TCR_OUT / "TCR_EXTENSION_REVIEWER_QA.md", "reports/TCR_EXTENSION_REVIEWER_QA.md"),
    (TCR_OUT / "tcr_registry_source_counts.tsv", "tables/tcr_registry_source_counts.tsv"),
    (TCR_OUT / "tcr_registry_missingness.tsv", "tables/tcr_registry_missingness.tsv"),
    (TCR_OUT / "tcr_registry_label_balance.tsv", "tables/tcr_registry_label_balance.tsv"),
    (TCR_OUT / "tcr_neo_linkage_summary.tsv", "tables/tcr_neo_linkage_summary.tsv"),
    (TCR_OUT / "model_discovery/tcr_model_discovery_matrix.tsv", "tables/tcr_model_discovery_matrix.tsv"),
    (TCR_OUT / "model_discovery/external_tcr_model_pilot_compare/external_tcr_model_pilot_comparison.tsv", "tables/external_tcr_model_pilot_comparison.tsv"),
    (TCR_OUT / "model_discovery/external_tcr_expert_benchmark/external_tcr_expert_benchmark_metrics.tsv", "tables/external_tcr_expert_benchmark_metrics.tsv"),
    (TCR_OUT / "model_discovery/external_tcr_expert_benchmark/external_tcr_expert_source_heldout_metrics.tsv", "tables/external_tcr_expert_source_heldout_metrics.tsv"),
    (TCR_OUT / "model_discovery/wetlab_external_tcr_expert_scores/wetlab_candidates_external_tcr_expert_ranked.tsv", "tables/wetlab_candidates_external_tcr_expert_ranked.tsv"),
    (TCR_OUT / "model_discovery/wetlab_external_tcr_expert_scores/wetlab_external_tcr_expert_pair_scores.tsv", "tables/wetlab_external_tcr_expert_pair_scores.tsv"),
    (TCR_OUT / "md_escalation/md_escalation_queue.tsv", "tables/md_escalation_queue.tsv"),
    (TCR_OUT / "md_escalation/md_escalation_queue_top20.tsv", "tables/md_escalation_queue_top20.tsv"),
    (TCR_OUT / "md_escalation/md_simulation_tiers.tsv", "tables/md_simulation_tiers.tsv"),
    (TCR_OUT / "md_escalation/md_jobs/md_job_manifest.tsv", "tables/md_job_manifest.tsv"),
    (TCR_OUT / "md_escalation/md_jobs/md_engine_availability.tsv", "tables/md_engine_availability.tsv"),
    (TCR_OUT / "md_escalation/md_jobs/submit_md_jobs.template.sh", "tables/submit_md_jobs.template.sh"),
    (TCR_OUT / "md_escalation/p0_structures/p0_exact_pdb_registry_rows.tsv", "tables/p0_exact_pdb_registry_rows.tsv"),
    (TCR_OUT / "md_escalation/p0_structures/p0_md_pilot_ready_complexes.tsv", "tables/p0_md_pilot_ready_complexes.tsv"),
    (TCR_OUT / "md_escalation/p0_structures/p0_pdb_chain_qc.tsv", "tables/p0_pdb_chain_qc.tsv"),
    (TCR_OUT / "md_escalation/p0_structures/openmm_minimized/openmm_p0_minimization_qc.tsv", "tables/openmm_p0_minimization_qc.tsv"),
    (TCR_OUT / "md_escalation/p0_structures/openmm_pilot_10ns_package/openmm_p0_pilot_manifest.tsv", "tables/openmm_p0_pilot_manifest.tsv"),
    (TCR_OUT / "md_escalation/p0_structures/openmm_pilot_10ns_package/p0_openmm_smoke_qc.tsv", "tables/p0_openmm_smoke_qc.tsv"),
    (TCR_OUT / "md_escalation/p0_structures/openmm_pilot_10ns_package/run_openmm_pilot.py", "scripts/run_openmm_pilot.py"),
    (TCR_OUT / "md_escalation/p0_structures/openmm_pilot_10ns_package/analyze_openmm_pilot.py", "scripts/analyze_openmm_pilot.py"),
    (TCR_OUT / "md_escalation/p0_structures/openmm_pilot_10ns_package/submit_p0_10ns_pilots.sh", "scripts/submit_p0_10ns_pilots.sh"),
    (TCR_OUT / "md_escalation/p0_structures/openmm_pilot_10ns_package/runpod_setup_openmm.sh", "scripts/runpod_setup_openmm.sh"),
    (TCR_OUT / "md_escalation/p0_structures/openmm_p0_10ns_pilot_package_2026_05_09.tar.gz", "packages/openmm_p0_10ns_pilot_package_2026_05_09.tar.gz"),
]


def copy_file(src: Path, rel: str) -> tuple[str, int]:
    dst = PACKET / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return rel, dst.stat().st_size


def main() -> None:
    shutil.rmtree(PACKET, ignore_errors=True)
    PACKET.mkdir(parents=True, exist_ok=True)
    copied: list[tuple[str, int]] = []
    missing: list[str] = []

    for src, rel in FILES:
        if src.exists():
            copied.append(copy_file(src, rel))
        else:
            missing.append(str(src))

    fig_dir = TCR_OUT / "figures"
    for src in sorted(fig_dir.glob("fig_tcr*.*")):
        copied.append(copy_file(src, f"figures/{src.name}"))

    readme = [
        "# CROSS-Neo-TCR Reviewer Packet",
        "",
        "Compact packet for reviewing the TCR-aware diagnostic extension.",
        "",
        "Decision: PROMOTE_TO_DIAGNOSTIC_CASE_STUDY and PROMOTE_TO_WETLAB_PRIORITIZATION_TOOL. Hold main-method promotion until strict external paired-TCR benchmarks and structure/wetlab validation are complete.",
        "",
        "Key headline:",
        "",
        "- pMTnet + TEPCAM mean ensemble: AUPRC 0.7841, AUROC 0.7652 on six fast challenge panels.",
        "- Top external-supported wetlab candidate: HMTEVVRHC / HLA-A*02:01.",
        "- Second external-supported candidate: GADGVGKSAL / HLA-C*08:02.",
        "- MD escalation queue: P0 TCR-pMHC MD for HMTEVVRHC/HLA-A*02:01 and GADGVGKSAL/HLA-C*08:02; existing PDB pilots are available for both.",
        "- OpenMM repair/minimization sanity check succeeded for representative P0 complexes 6VRN and 6UON.",
        "- Local OpenMM smoke dynamics succeeded for both P0 representatives; production still requires CUDA explicit-solvent pilot MD.",
        "",
        "Contents:",
        "",
        "- `reports/`: decision report, literature/tool audit, registry/linkage reports, model audit, benchmark reports, manuscript-track scaffold.",
        "- `tables/`: source counts, linkage summary, external model metrics, source-heldout metrics, wetlab ranked candidates, MD escalation queue.",
        "- `figures/`: seven TCR extension figures as PNG and PDF.",
        "",
        "Claim boundary:",
        "",
        "- This packet supports an optional TCR-aware diagnostic/wetlab-prioritization branch.",
        "- It does not support universal TCR-aware neoantigen prediction or clinical utility.",
        "- Synthetic decoy benchmarks are model-selection diagnostics, not final external SOTA evidence.",
        "",
    ]
    if missing:
        readme.extend(["Missing expected files:", ""])
        readme.extend([f"- `{m}`" for m in missing])
        readme.append("")
    readme.extend(["File manifest:", "", "| path | bytes |", "|---|---:|"])
    for rel, size in sorted(copied):
        readme.append(f"| `{rel}` | {size} |")
    (PACKET / "README.md").write_text("\n".join(readme) + "\n")
    copied.append(("README.md", (PACKET / "README.md").stat().st_size))

    if ZIP_PATH.exists():
        ZIP_PATH.unlink()
    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(PACKET.rglob("*")):
            if p.is_file():
                zf.write(p, p.relative_to(PACKET))

    HUB_ASSET.mkdir(parents=True, exist_ok=True)
    WEB_ASSET.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ZIP_PATH, HUB_ASSET / ZIP_PATH.name)
    shutil.copy2(ZIP_PATH, WEB_ASSET / ZIP_PATH.name)

    print(f"[tcr-reviewer-packet] files={len(copied)} zip={ZIP_PATH} size={ZIP_PATH.stat().st_size}")
    print(f"[tcr-reviewer-packet] deployed={HUB_ASSET / ZIP_PATH.name}")


if __name__ == "__main__":
    main()
