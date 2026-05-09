#!/usr/bin/env python3
"""Input/status scaffold for extending the TCGA pathology-AI audit.

This scaffold is deliberately conservative: it does not fabricate LUAD/BRCA
labels or assume that WSI features exist. It validates the expected files,
writes a machine-readable status report, and creates a small markdown summary
that can be filled once the per-slide UNI features and manifests are present.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path

import pandas as pd


ROOT = Path("/home/seungho/personal/THCA_data_analysis")
BASE = ROOT / "project/results/p2_image_dm1_v2_foundation_clam_2026_05_07"
OUT = BASE / "analysis_supp/audit_multicohort_generalization"


@dataclass
class CohortConfig:
    cohort: str
    proposed_task: str
    feature_dir: str
    manifest_tsv: str
    prediction_tsv: str
    metadata_tsv: str
    label_column: str = "label"
    positive_label: str = "1"
    tss_column: str = "tss"
    clinical_columns: tuple[str, ...] = ("histology", "sex", "tss")
    subgroup_columns: tuple[str, ...] = ("subtype", "histology")


CONFIGS = {
    "LUAD": CohortConfig(
        cohort="TCGA-LUAD",
        proposed_task=(
            "candidate binary pathology task to lock after metadata assembly "
            "(e.g. EGFR-mutant vs KRAS/STK11/other or TCGA subtype label)"
        ),
        feature_dir=str(BASE / "phase2_tcga_luad_UNI/features"),
        manifest_tsv=str(BASE / "phase2_tcga_luad_UNI/slide_manifest.tsv"),
        prediction_tsv=str(BASE / "phase2_tcga_luad_UNI/clam_per_slide_predictions.tsv"),
        metadata_tsv=str(ROOT / "project/metadata/tcga_luad_audit_metadata.tsv"),
    ),
    "BRCA": CohortConfig(
        cohort="TCGA-BRCA",
        proposed_task=(
            "candidate binary pathology task to lock after metadata assembly "
            "(e.g. basal-like vs luminal or ER-negative vs ER-positive)"
        ),
        feature_dir=str(BASE / "phase2_tcga_brca_UNI/features"),
        manifest_tsv=str(BASE / "phase2_tcga_brca_UNI/slide_manifest.tsv"),
        prediction_tsv=str(BASE / "phase2_tcga_brca_UNI/clam_per_slide_predictions.tsv"),
        metadata_tsv=str(ROOT / "project/metadata/tcga_brca_audit_metadata.tsv"),
    ),
}


REQUIRED_MANIFEST_COLUMNS = {
    "file_id",
    "submitter_id",
    "label",
    "tss",
}
REQUIRED_PRED_COLUMNS = {
    "slide",
    "fold",
    "label",
    "prob_pos",
    "n_tiles",
}


def inspect_tsv(path: Path, required: set[str]) -> dict:
    if not path.exists():
        return {
            "exists": False,
            "n_rows": 0,
            "columns": [],
            "missing_columns": sorted(required),
        }
    df = pd.read_csv(path, sep="\t", nrows=5)
    cols = list(df.columns)
    return {
        "exists": True,
        "n_rows": int(sum(1 for _ in path.open()) - 1),
        "columns": cols,
        "missing_columns": sorted(required.difference(cols)),
    }


def inspect_config(config: CohortConfig) -> dict:
    feat = Path(config.feature_dir)
    pt_files = sorted(feat.glob("*.pt")) if feat.exists() else []
    status = {
        "config": asdict(config),
        "feature_dir": {
            "exists": feat.exists(),
            "n_pt_files": len(pt_files),
            "example_pt_files": [p.name for p in pt_files[:5]],
        },
        "manifest": inspect_tsv(Path(config.manifest_tsv), REQUIRED_MANIFEST_COLUMNS),
        "predictions": inspect_tsv(Path(config.prediction_tsv), REQUIRED_PRED_COLUMNS),
        "metadata": inspect_tsv(Path(config.metadata_tsv), set()),
    }
    ready = (
        status["feature_dir"]["n_pt_files"] > 0
        and status["manifest"]["exists"]
        and not status["manifest"]["missing_columns"]
        and status["metadata"]["exists"]
    )
    status["ready_for_12_test_audit"] = bool(ready)
    blockers = []
    if status["feature_dir"]["n_pt_files"] == 0:
        blockers.append("missing per-slide UNI .pt feature files")
    if not status["manifest"]["exists"]:
        blockers.append("missing slide_manifest.tsv")
    elif status["manifest"]["missing_columns"]:
        blockers.append(
            "slide_manifest.tsv missing columns: "
            + ", ".join(status["manifest"]["missing_columns"])
        )
    if not status["metadata"]["exists"]:
        blockers.append("missing cohort audit metadata TSV")
    if not status["predictions"]["exists"]:
        blockers.append(
            "missing OOF predictions; needed after first CLAM run, not for input staging"
        )
    status["blockers"] = blockers
    return status


def write_markdown(status: dict) -> Path:
    cfg = status["config"]
    cohort_slug = cfg["cohort"].replace("-", "_")
    out = OUT / f"{cohort_slug}_input_status.md"
    lines = [
        f"# {cfg['cohort']} pathology-AI audit input status",
        "",
        f"**Proposed task:** {cfg['proposed_task']}",
        "",
        f"**Ready for 12-test audit:** `{status['ready_for_12_test_audit']}`",
        "",
        "## Required inputs",
        "",
        "| Input | Expected path | Status |",
        "|---|---|---|",
        [
            "Feature dir",
            cfg["feature_dir"],
            f"{status['feature_dir']['n_pt_files']} .pt files",
        ],
        [
            "Manifest",
            cfg["manifest_tsv"],
            "exists" if status["manifest"]["exists"] else "missing",
        ],
        [
            "Metadata",
            cfg["metadata_tsv"],
            "exists" if status["metadata"]["exists"] else "missing",
        ],
        [
            "OOF predictions",
            cfg["prediction_tsv"],
            "exists" if status["predictions"]["exists"] else "missing",
        ],
    ]
    table_rows = lines[8:]
    lines = lines[:8]
    for row in table_rows:
        lines.append("| " + " | ".join(str(x) for x in row) + " |")
    lines.extend([
        "",
        "## Blockers",
        "",
    ])
    if status["blockers"]:
        lines.extend(f"- {b}" for b in status["blockers"])
    else:
        lines.append("- none")
    lines.extend([
        "",
        "## Expected manifest schema",
        "",
        "`file_id`, `submitter_id`, `label`, `tss`; optional but recommended: "
        "`histology`, `sex`, `subtype`, `age`, `stage`, `site`.",
        "",
        "## Next command after feature extraction",
        "",
        "```bash",
        f"python3 {BASE}/scripts/audit_tcga_{cfg['cohort'].split('-')[-1].lower()}.py",
        "```",
        "",
    ])
    out.write_text("\n".join(lines))
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cohort", choices=sorted(CONFIGS))
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    status = inspect_config(CONFIGS[args.cohort])
    cohort_slug = status["config"]["cohort"].replace("-", "_")
    json_out = OUT / f"{cohort_slug}_input_status.json"
    json_out.write_text(json.dumps(status, indent=2))
    md_out = write_markdown(status)
    print(f"WROTE {json_out}")
    print(f"WROTE {md_out}")
    if not status["ready_for_12_test_audit"]:
        print("NOT READY:")
        for blocker in status["blockers"]:
            print(f"  - {blocker}")


if __name__ == "__main__":
    main()
