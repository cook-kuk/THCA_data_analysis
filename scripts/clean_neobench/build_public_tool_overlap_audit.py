#!/usr/bin/env python3
"""Build public-tool training overlap audit scaffold for CLEAN-NeoBench.

This script is deliberately conservative. Documentation-level provenance is
useful, but it is not a row-level training-corpus overlap audit. Public
pretrained tools remain caveated unless an explicit row-level audit artifact is
present.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from common import PUBLIC_METHOD_TOKENS, ensure_dir, method_family_and_role, normalize_empty, update_manifest, write_tsv


SEARCH_DIRS = [
    "project/results/p_neo_bayesian_2026_05_09",
    "project/results/p_neo_bayesian_2026_05_09/curation_2026_05_09",
    "project/results/p_neo_bayesian_2026_05_09/clean_neo_v0_2026_05_09",
    "project/results/cross_neo_v0",
    "project/results/cross_neo_v1",
    "project/results/cross_neo_v1_lockdown",
]

PUBLIC_TOOL_NOTES = {
    "mhcflurry": "public pretrained MHC-I presentation/binding tool; row-level training overlap unresolved",
    "netmhc": "public pretrained NetMHC-family tool; row-level training overlap unresolved",
    "bigmhc": "public deep presentation/immunogenicity model; row-level training overlap unresolved",
    "prime": "public immunogenicity predictor; row-level training overlap unresolved",
    "mixmhc": "public MHC ligand predictor; row-level training overlap unresolved",
    "netmhcstab": "public peptide-MHC stability predictor; row-level training overlap unresolved",
    "immunostruct": "public/pretrained structure-aware comparator; row-level training overlap unresolved",
    "neomhci": "public/pretrained neoantigen comparator; row-level training overlap unresolved",
    "unipmt": "public/pretrained pMHC/TCR-style comparator; row-level training overlap unresolved",
    "deephla": "public deep HLA/presentation comparator; row-level training overlap unresolved",
    "deepimmuno": "public immunogenicity model; row-level training overlap unresolved",
    "deephimmuno": "public immunogenicity model; row-level training overlap unresolved",
    "transphla": "public pHLA binding transformer; row-level training overlap unresolved",
    "mhcnuggets": "public MHC binding predictor; row-level training overlap unresolved",
    "tscape": "public multidomain immunogenicity/TCR model; row-level training overlap unresolved",
    "titanian": "public multidomain immunogenicity/TCR model; row-level training overlap unresolved",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output directory")
    return parser.parse_args()


def boolish(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def load_optional(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def discover_audit_files(repo_root: Path) -> list[Path]:
    files: list[Path] = []
    for rel in SEARCH_DIRS:
        base = repo_root / rel
        if not base.exists():
            continue
        for p in base.rglob("*training*audit*.tsv"):
            if p.is_file():
                files.append(p)
        for p in base.rglob("*overlap*audit*.tsv"):
            if p.is_file():
                files.append(p)
    return sorted(set(files))


def load_provenance(audit_files: list[Path]) -> pd.DataFrame:
    frames = []
    for path in audit_files:
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            continue
        if df.empty:
            continue
        method_col = next((c for c in ["method_name", "method", "tool", "predictor"] if c in df.columns), None)
        if method_col is None:
            continue
        keep = df.copy()
        keep["method_name_audit"] = keep[method_col].astype(str)
        keep["audit_source_file"] = str(path)
        frames.append(keep)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def provenance_for_method(method: str, provenance: pd.DataFrame) -> dict[str, str]:
    if provenance.empty:
        return {}
    low = method.lower()
    names = provenance["method_name_audit"].astype(str).str.lower()
    mask = names.eq(low) | names.apply(lambda x: x in low or low in x)
    if not mask.any():
        tokens = [tok for tok in PUBLIC_METHOD_TOKENS if tok in low]
        if tokens:
            token = tokens[0]
            mask = names.str.contains(token, regex=False)
    if not mask.any():
        return {}
    row = provenance.loc[mask].iloc[0]
    out = {}
    for src, dst in [
        ("training_data_summary", "training_data_summary"),
        ("training_overlap_risk", "training_overlap_risk"),
        ("paper_disposition", "paper_disposition"),
        ("provenance_confidence", "provenance_confidence"),
        ("training_data_note", "training_data_note"),
        ("source_urls", "source_urls"),
        ("audit_source_file", "audit_source_file"),
    ]:
        if src in row.index:
            out[dst] = normalize_empty(row[src])
    return out


def public_note(method_name: str) -> str:
    low = method_name.lower()
    for token, note in PUBLIC_TOOL_NOTES.items():
        if token in low:
            return note
    return "public/pretrained comparator detected; row-level training overlap unresolved"


def build_audit(output_root: Path, provenance: pd.DataFrame, audit_files: list[Path]) -> pd.DataFrame:
    method_scores = load_optional(output_root / "clean_neobench_method_scores.tsv")
    leaderboard = load_optional(output_root / "clean_neobench_leaderboard.tsv")
    if method_scores.empty and leaderboard.empty:
        return pd.DataFrame()
    if not method_scores.empty:
        method_table = method_scores.drop_duplicates("method_name").copy()
    else:
        method_table = leaderboard.drop_duplicates("method_name").copy()

    rows = []
    row_level_available = False
    row_level_files = []
    for path in audit_files:
        try:
            cols = pd.read_csv(path, sep="\t", nrows=2).columns
        except Exception:
            continue
        lower_cols = {c.lower() for c in cols}
        if {"candidate_id", "method_name"}.issubset(lower_cols) or {"peptide", "hla", "method_name"}.issubset(lower_cols):
            row_level_available = True
            row_level_files.append(str(path))

    for _, r in method_table.iterrows():
        method_name = str(r.get("method_name", ""))
        family, role, uses_public_detected, _, _, caveat_default = method_family_and_role(method_name)
        uses_public = boolish(r.get("uses_public_pretraining", uses_public_detected)) or uses_public_detected
        role = str(r.get("method_role", role))
        family = str(r.get("method_family", family))
        prov = provenance_for_method(method_name, provenance)

        if uses_public:
            audit_status = "unresolved_no_row_level_training_corpus"
            clean_allowed = False
            disposition = "caveated_public_comparator_only"
            caveat = public_note(method_name)
            overlap_detail = "No row-level public-tool training corpus matched to candidate_id/peptide-HLA was found."
        else:
            audit_status = "internal_or_local_known_split"
            clean_allowed = boolish(r.get("clean_comparator_allowed", True))
            disposition = "eligible_internal_comparator_if_split_contract_passes"
            caveat = normalize_empty(r.get("caveat", caveat_default))
            overlap_detail = "Not a public pretrained tool under current token rules."

        if uses_public and row_level_available:
            # The artifact exists, but this scaffold does not promote until an
            # explicit overlap row count for this method is present.
            audit_status = "row_level_audit_artifact_present_but_method_unresolved"
            overlap_detail = "A row-level audit artifact exists, but no method-specific clean pass was established by this scaffold."

        rows.append(
            {
                "method_name": method_name,
                "method_family": family,
                "method_role": "caveated_public_comparator" if uses_public else role,
                "uses_public_pretraining": uses_public,
                "training_overlap_audit_status": audit_status,
                "training_overlap_audited_row_level": False if uses_public else True,
                "row_level_audit_files_detected": "; ".join(row_level_files),
                "public_tool_training_overlap_any": True if uses_public else False,
                "candidate_overlap_rows": "NA",
                "candidate_overlap_fraction": "NA",
                "clean_comparator_allowed_after_audit": clean_allowed,
                "reviewer_disposition": disposition,
                "caveat": caveat,
                "overlap_detail": overlap_detail,
                "training_data_summary": prov.get("training_data_summary", ""),
                "training_overlap_risk": prov.get("training_overlap_risk", ""),
                "provenance_confidence": prov.get("provenance_confidence", ""),
                "audit_source_file": prov.get("audit_source_file", ""),
                "source_urls": prov.get("source_urls", ""),
            }
        )
    audit = pd.DataFrame(rows).sort_values(["uses_public_pretraining", "method_name"], ascending=[False, True])
    return audit.replace({"": "NA"})


def update_method_score_flags(output_root: Path, audit: pd.DataFrame) -> None:
    if audit.empty:
        return
    audit_idx = audit.set_index("method_name")
    for name in ["clean_neobench_method_scores.tsv", "clean_neobench_leaderboard.tsv"]:
        path = output_root / name
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        if "method_name" not in df.columns:
            continue
        public_mask = df["method_name"].isin(audit_idx.index[audit_idx["uses_public_pretraining"].astype(bool)])
        if public_mask.any():
            df.loc[public_mask, "uses_public_pretraining"] = True
            df.loc[public_mask, "training_overlap_audited"] = False
            df.loc[public_mask, "clean_comparator_allowed"] = False
            if "method_role" in df.columns:
                df.loc[public_mask, "method_role"] = "caveated_public_comparator"
            if "caveat" in df.columns:
                df.loc[public_mask, "caveat"] = "public pretrained comparator; row-level training overlap unresolved"
        write_tsv(df, path)


def write_report(output_root: Path, audit: pd.DataFrame, audit_files: list[Path]) -> None:
    public = audit[audit["uses_public_pretraining"].astype(bool)] if not audit.empty else pd.DataFrame()
    internal = audit[~audit["uses_public_pretraining"].astype(bool)] if not audit.empty else pd.DataFrame()
    public_table = public[
        [
            "method_name",
            "method_role",
            "training_overlap_audit_status",
            "clean_comparator_allowed_after_audit",
            "reviewer_disposition",
            "caveat",
        ]
    ].to_markdown(index=False) if not public.empty else "No public pretrained methods detected."
    files_md = "\n".join(f"- `{p}`" for p in audit_files) if audit_files else "- No training/overlap audit TSV artifacts discovered."
    text = f"""# CLEAN-NeoBench Public Tool Overlap Audit

## Summary

This audit scaffold separates documentation-level provenance from row-level training-corpus overlap auditing.

- Public pretrained methods detected: {len(public)}
- Internal/local methods detected: {len(internal)}
- Public methods promoted to clean comparators: 0
- Row-level clean-pass status: unresolved unless a future method-specific candidate/peptide-HLA audit table proves otherwise.

## Public Method Disposition

{public_table}

## Audit Artifacts Discovered

{files_md}

## Rule

Public pretrained outputs may be reported as caveated comparators. They must not be used as clean training features or described as clean external baselines until row-level training-corpus overlap is audited.

## Forbidden Claims

- Public tools are clean baselines without overlap audit.
- External validation is proven by public pretrained tool agreement.
- A public pretrained score can rescue a high-leakage candidate into a clean benchmark claim.
"""
    (output_root / "CLEAN_NEOBENCH_PUBLIC_OVERLAP_AUDIT.md").write_text(text.rstrip() + "\n")


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    ensure_dir(output_root)

    audit_files = discover_audit_files(repo_root)
    provenance = load_provenance(audit_files)
    audit = build_audit(output_root, provenance, audit_files)
    write_tsv(audit, output_root / "clean_neobench_public_tool_overlap_audit.tsv")
    update_method_score_flags(output_root, audit)
    write_report(output_root, audit, audit_files)

    public_count = int(audit["uses_public_pretraining"].astype(bool).sum()) if not audit.empty else 0
    update_manifest(
        output_root,
        "public_tool_overlap_audit",
        {
            "n_methods_audited": int(len(audit)),
            "n_public_methods": public_count,
            "n_public_clean_comparators_allowed": 0,
            "audit_files_discovered": [str(p) for p in audit_files],
            "warnings": [
                "Public pretrained methods remain caveated until row-level training-corpus overlap audit exists."
            ]
            if public_count
            else [],
        },
    )

    manifest_path = output_root / "run_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    else:
        manifest = {}
    manifest.setdefault("output_files", [])
    for name in ["clean_neobench_public_tool_overlap_audit.tsv", "CLEAN_NEOBENCH_PUBLIC_OVERLAP_AUDIT.md"]:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"]["n_public_methods_overlap_unresolved"] = public_count
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"[clean-neobench-public-audit] methods={len(audit)} public_unresolved={public_count}")


if __name__ == "__main__":
    main()
