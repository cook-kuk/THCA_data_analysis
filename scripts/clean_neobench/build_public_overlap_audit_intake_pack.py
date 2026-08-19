#!/usr/bin/env python3
"""Build an intake pack for public-tool row-level training overlap audit."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


OUTPUTS = [
    "clean_neobench_public_overlap_audit_intake.tsv",
    "CLEAN_NEOBENCH_PUBLIC_OVERLAP_AUDIT_INTAKE.md",
    "CLEAN_NEOBENCH_PUBLIC_OVERLAP_AUDIT_INTAKE_KR.md",
    "public_tool_overlap_audit_intake/README.md",
    "public_tool_overlap_audit_intake/public_training_corpus_template.tsv",
]

TEMPLATE_COLUMNS = [
    "public_tool",
    "row_id",
    "peptide",
    "hla",
    "hla_allele_4digit",
    "source_dataset",
    "source_name",
    "train_split",
    "label",
    "assay_type",
    "mhc_class",
    "publication_or_url",
    "license_or_terms",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    parser.add_argument("--data-drop-root", default="project/data/public_training_corpora", help="Where row-level public tool corpora should be dropped")
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def build_intake(requirements: pd.DataFrame, audit: pd.DataFrame, data_drop_root: Path) -> pd.DataFrame:
    if requirements.empty:
        return pd.DataFrame(
            columns=[
                "public_tool",
                "audit_status",
                "expected_drop_file",
                "template_file",
                "minimum_required_fields",
                "clean_pass_rule",
                "current_disposition",
                "source_urls",
            ]
        )
    rows = []
    audit_idx = audit.drop_duplicates("method_name").set_index("method_name") if not audit.empty and "method_name" in audit.columns else pd.DataFrame()
    for _, r in requirements.iterrows():
        tool = str(r.get("public_tool", "unknown_public_tool"))
        audit_row = audit_idx.loc[tool] if not audit_idx.empty and tool in audit_idx.index else pd.Series(dtype=object)
        expected = data_drop_root / f"{slug(tool)}_training.tsv"
        template = Path("public_tool_overlap_audit_intake") / f"{slug(tool)}_training_template.tsv"
        rows.append(
            {
                "public_tool": tool,
                "audit_status": audit_row.get("training_overlap_audit_status", "unresolved_no_row_level_training_corpus"),
                "expected_drop_file": str(expected),
                "template_file": str(template),
                "filename_tokens": r.get("filename_tokens", ""),
                "minimum_required_fields": r.get("minimum_columns", ""),
                "required_key": r.get("required_key", ""),
                "clean_pass_rule": r.get("clean_pass_rule", ""),
                "current_clean_allowed": bool(audit_row.get("clean_comparator_allowed_after_audit", False)),
                "current_disposition": audit_row.get("reviewer_disposition", "caveated_public_comparator_only"),
                "training_overlap_risk": audit_row.get("training_overlap_risk", ""),
                "source_urls": audit_row.get("source_urls", ""),
            }
        )
    return pd.DataFrame(rows)


def write_templates(output_root: Path, intake_dir: Path, intake: pd.DataFrame, data_drop_root: Path) -> list[str]:
    ensure_dir(intake_dir)
    ensure_dir(data_drop_root)
    template = pd.DataFrame(columns=TEMPLATE_COLUMNS)
    write_tsv(template, intake_dir / "public_training_corpus_template.tsv")
    created = ["public_tool_overlap_audit_intake/public_training_corpus_template.tsv"]
    for tool in intake.get("public_tool", pd.Series(dtype=str)).astype(str):
        path = intake_dir / f"{slug(tool)}_training_template.tsv"
        tool_template = pd.DataFrame([{c: "" for c in TEMPLATE_COLUMNS}])
        tool_template.loc[0, "public_tool"] = tool
        write_tsv(tool_template, path)
        created.append(str(path.relative_to(output_root)))
    data_readme = data_drop_root / "README_PUBLIC_TRAINING_CORPORA.md"
    data_readme.write_text(
        """# Public Training Corpora Drop Directory

Place row-level public-tool training corpora here as TSV files.

Required minimum columns:
- `peptide`
- `hla` or `hla_allele_4digit` when available
- `public_tool`
- `source_dataset` or `source_name` when available
- `train_split`, `label`, `assay_type`, and `publication_or_url` when available

After adding files, rerun:

```bash
python scripts/clean_neobench/run_clean_neobench_pipeline.py \\
  --repo-root . \\
  --output-root project/results/clean_neobench_barneo_2026_05_09
```

Public pretrained tools remain caveated until row-level exact peptide, exact peptide-HLA, and near-peptide overlap audit passes.
"""
    )
    return created


def write_reports(output_root: Path, intake: pd.DataFrame, data_drop_root: Path, created_templates: list[str]) -> None:
    text = f"""# CLEAN-NeoBench Public Overlap Audit Intake

## Purpose

This intake pack defines exactly what row-level public-tool training corpora must be supplied before any public pretrained tool can be treated as a clean comparator. Until these files are provided and audited, public tools remain caveated comparators.

## Drop Location

`{data_drop_root}`

## Intake Table

{dataframe_to_markdown(intake, max_rows=40)}

## Required Template Columns

`{'`, `'.join(TEMPLATE_COLUMNS)}`

## Created Templates

{dataframe_to_markdown(pd.DataFrame({"template": created_templates}), max_rows=40)}

## Pass Rule

For a public tool to become a clean comparator, row-level audit must show no exact peptide-HLA training overlap and no unresolved exact peptide overlap. Near-peptide overlap must be summarized separately. Documentation-only provenance is not enough.

## Claim Boundary

Allowed: public tool caveated comparison, row-level overlap audit plan, clean-comparator eligibility after audit.

Forbidden: public pretrained tools as clean baselines without row-level corpus audit.
"""
    (output_root / "CLEAN_NEOBENCH_PUBLIC_OVERLAP_AUDIT_INTAKE.md").write_text(text.strip() + "\n")

    kr = f"""# CLEAN-NeoBench Public Overlap Audit Intake KR

## 한 줄 결론

public tool을 clean comparator로 올리려면 row-level training corpus 파일이 필요하다. 현재는 파일이 0개라 public clean comparator allowed는 0이 맞다.

## 어디에 넣나

`{data_drop_root}`

## 무엇을 넣나

{dataframe_to_markdown(intake[["public_tool", "expected_drop_file", "template_file", "required_key", "clean_pass_rule", "current_disposition"]], max_rows=40)}

## 최소 컬럼

`peptide`, `hla` 또는 `hla_allele_4digit`, `public_tool`이 핵심이다. 가능하면 `source_dataset`, `train_split`, `label`, `assay_type`, `publication_or_url`도 넣는다.

## 다시 실행

```bash
python scripts/clean_neobench/run_clean_neobench_pipeline.py \\
  --repo-root . \\
  --output-root project/results/clean_neobench_barneo_2026_05_09
```

## 정책

- exact peptide-HLA overlap 있으면 clean comparator 불가.
- exact peptide overlap이 unresolved이면 clean comparator 불가.
- near-peptide overlap은 별도 caveat로 보고한다.
- public tool이 점수상 이겨도 audit 전에는 caveated comparator다.
"""
    (output_root / "CLEAN_NEOBENCH_PUBLIC_OVERLAP_AUDIT_INTAKE_KR.md").write_text(kr.strip() + "\n")

    intake_readme = output_root / "public_tool_overlap_audit_intake" / "README.md"
    intake_readme.write_text(
        f"""# Public Tool Overlap Audit Intake

Use these templates to normalize row-level training corpora for public pretrained tools.

Drop completed training corpus TSVs into:

`{data_drop_root}`

Then rerun the CLEAN-NeoBench pipeline. Public tools remain caveated until the row-level audit passes.
"""
    )


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    data_drop_root = Path(args.data_drop_root)
    if not data_drop_root.is_absolute():
        data_drop_root = repo_root / data_drop_root
    ensure_dir(output_root)

    requirements = read_tsv(output_root / "clean_neobench_public_training_corpus_requirements.tsv")
    audit = read_tsv(output_root / "clean_neobench_public_tool_overlap_audit.tsv")
    intake = build_intake(requirements, audit, data_drop_root)
    intake_dir = output_root / "public_tool_overlap_audit_intake"
    created_templates = write_templates(output_root, intake_dir, intake, data_drop_root)

    write_tsv(intake, output_root / "clean_neobench_public_overlap_audit_intake.tsv")
    write_reports(output_root, intake, data_drop_root, created_templates)

    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("output_files", [])
    for name in OUTPUTS:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_public_overlap_intake_tools": int(len(intake)),
            "public_training_corpora_drop_root": str(data_drop_root),
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "public_overlap_audit_intake_pack",
        {
            "outputs": OUTPUTS + created_templates,
            "n_tools": int(len(intake)),
            "data_drop_root": str(data_drop_root),
            "warnings": ["Public pretrained tools remain caveated until row-level training corpus overlap audit passes."],
        },
    )
    print(f"[clean-neobench-public-intake] tools={len(intake)} drop_root={data_drop_root}")


if __name__ == "__main__":
    main()
