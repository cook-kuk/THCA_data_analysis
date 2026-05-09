#!/usr/bin/env python3
"""Row-level public training-corpus overlap audit for CLEAN-NeoBench.

The script is designed to be useful even before public training corpora are
available locally. If corpora are present, it compares benchmark candidates to
public training rows by exact peptide-HLA, exact peptide, and near-peptide
similarity. If corpora are absent, it writes an explicit requirements table and
keeps public tools unresolved.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

from common import ensure_dir, near_similarity, normalize_empty, normalize_hla, update_manifest, write_tsv


PUBLIC_TOOL_SPECS = [
    {
        "public_tool": "MHCflurry",
        "tokens": ["mhcflurry"],
        "required_key": "peptide + HLA allele preferred; peptide-only accepted as weaker audit",
        "clean_pass_rule": "no exact peptide-HLA and no exact peptide overlap; near-peptide overlap summarized separately",
    },
    {
        "public_tool": "NetMHCpan_4.1",
        "tokens": ["netmhcpan", "netmhc"],
        "required_key": "peptide + HLA allele",
        "clean_pass_rule": "no exact peptide-HLA overlap against BA/EL training rows",
    },
    {
        "public_tool": "BigMHC_IM",
        "tokens": ["bigmhc"],
        "required_key": "peptide + HLA allele + split/source if available",
        "clean_pass_rule": "no exact peptide-HLA or exact peptide overlap against im_train/el_train rows",
    },
    {
        "public_tool": "PRIME",
        "tokens": ["prime"],
        "required_key": "peptide + HLA allele preferred",
        "clean_pass_rule": "no exact peptide-HLA and no exact peptide overlap against immunogenicity training rows",
    },
    {
        "public_tool": "DeepImmuno",
        "tokens": ["deepimmuno", "deephimmuno"],
        "required_key": "peptide + HLA allele",
        "clean_pass_rule": "no exact peptide-HLA overlap against IEDB-derived training rows",
    },
    {
        "public_tool": "MHCnuggets_2",
        "tokens": ["mhcnuggets"],
        "required_key": "peptide + HLA allele",
        "clean_pass_rule": "no exact peptide-HLA overlap against IEDB-derived binding rows",
    },
    {
        "public_tool": "NetMHCstabpan",
        "tokens": ["netmhcstab"],
        "required_key": "peptide + HLA allele",
        "clean_pass_rule": "no exact peptide-HLA overlap against stability training rows",
    },
    {
        "public_tool": "TransPHLA",
        "tokens": ["transphla"],
        "required_key": "peptide + HLA allele",
        "clean_pass_rule": "no exact peptide-HLA overlap against pHLA training rows",
    },
    {
        "public_tool": "TSCAPE_TITANiAN",
        "tokens": ["tscape", "titanian"],
        "required_key": "peptide + HLA allele + TCR fields when available",
        "clean_pass_rule": "no exact peptide-HLA overlap; TCR/pMHC train rows audited separately when available",
    },
]

DEFAULT_CORPUS_DIRS = [
    "project/data/public_training_corpora",
    "project/data/neoantigen_public_training",
    "project/data/public_tools",
    "project/results/clean_neobench_barneo_2026_05_09/public_training_corpora",
    "project/results/p_neo_bayesian_2026_05_09/public_training_corpora",
    "project/results/cross_neo_v1/public_training_corpora",
]
PEPTIDE_COLS = [
    "peptide",
    "epitope",
    "sequence",
    "peptide_sequence",
    "ligand",
    "mt_peptide",
    "mut_peptide",
    "mutant_peptide",
    "antigen_peptide",
]
HLA_COLS = ["hla", "allele", "mhc", "hla_allele", "mhc_allele", "restriction", "hla_type"]
LABEL_COLS = ["label", "immunogenic", "binder", "target", "y", "response", "measurement_value"]
OVERLAP_COLUMNS = [
    "candidate_id",
    "candidate_peptide",
    "candidate_hla",
    "label",
    "source_name",
    "public_peptide",
    "public_hla",
    "public_tool",
    "public_corpus_file",
    "public_row_index",
    "public_training_label",
    "peptide_column",
    "hla_column",
    "label_column",
    "overlap_type",
    "near_similarity",
]
SUMMARY_COLUMNS = [
    "public_tool",
    "public_corpus_file",
    "n_candidate_overlaps",
    "candidate_overlap_fraction",
    "n_exact_peptide_hla",
    "n_exact_peptide",
    "n_near_peptide",
    "audit_status",
    "clean_comparator_allowed_after_row_audit",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="CLEAN-NeoBench output directory")
    parser.add_argument("--public-corpus-root", action="append", default=[], help="Optional directory with public training corpus TSV/CSV files")
    parser.add_argument("--near-threshold", type=float, default=0.80, help="Near-peptide similarity threshold")
    parser.add_argument("--max-file-mb", type=float, default=250.0, help="Skip individual corpus files larger than this")
    return parser.parse_args()


def infer_tool(path: Path) -> str:
    low = path.name.lower()
    for spec in PUBLIC_TOOL_SPECS:
        if any(tok in low for tok in spec["tokens"]):
            return spec["public_tool"]
    return "unknown_public_training_corpus"


def discover_corpus_files(repo_root: Path, extra_roots: list[str], max_file_mb: float) -> tuple[list[Path], list[str]]:
    warnings = []
    roots = [repo_root / rel for rel in DEFAULT_CORPUS_DIRS] + [Path(p).expanduser() for p in extra_roots]
    files: list[Path] = []
    for root in roots:
        if not root.exists():
            continue
        for suffix in ("*.tsv", "*.csv", "*.txt"):
            for path in root.rglob(suffix):
                if not path.is_file():
                    continue
                low = path.name.lower()
                if not any(tok in low for spec in PUBLIC_TOOL_SPECS for tok in spec["tokens"]) and "training" not in low and "train" not in low:
                    continue
                size_mb = path.stat().st_size / (1024 * 1024)
                if size_mb > max_file_mb:
                    warnings.append(f"Skipped oversized corpus file: {path} ({size_mb:.1f} MB)")
                    continue
                files.append(path)
    return sorted(set(files)), warnings


def choose_col(columns: list[str], candidates: list[str]) -> str:
    lower = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    for col in columns:
        low = col.lower()
        if any(cand.lower() in low for cand in candidates):
            return col
    return ""


def normalize_peptide(value: Any) -> str:
    s = normalize_empty(value).upper()
    s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s)
    return s


def read_corpus(path: Path) -> tuple[pd.DataFrame, list[str]]:
    warnings: list[str] = []
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    try:
        df = pd.read_csv(path, sep=sep, low_memory=False)
    except Exception as exc:
        return pd.DataFrame(), [f"Failed to read {path}: {exc}"]
    if df.empty:
        return df, [f"Corpus file has no rows: {path}"]
    pep_col = choose_col(list(df.columns), PEPTIDE_COLS)
    hla_col = choose_col(list(df.columns), HLA_COLS)
    label_col = choose_col(list(df.columns), LABEL_COLS)
    if not pep_col:
        return pd.DataFrame(), [f"Skipped {path}: no peptide-like column found"]
    out = pd.DataFrame(
        {
            "public_row_index": df.index.astype(int),
            "public_peptide": df[pep_col].map(normalize_peptide),
            "public_hla": df[hla_col].map(normalize_hla) if hla_col else "",
            "public_training_label": df[label_col] if label_col else "",
        }
    )
    out["public_tool"] = infer_tool(path)
    out["public_corpus_file"] = str(path)
    out["peptide_column"] = pep_col
    out["hla_column"] = hla_col
    out["label_column"] = label_col
    out = out[out["public_peptide"].astype(str).ne("")].copy()
    return out, warnings


def master_keys(master: pd.DataFrame) -> pd.DataFrame:
    out = master[["candidate_id", "peptide", "hla", "hla_allele_4digit", "label", "source_name"]].copy()
    out["candidate_peptide"] = out["peptide"].map(normalize_peptide)
    out["candidate_hla"] = out["hla_allele_4digit"].map(normalize_hla)
    out.loc[out["candidate_hla"].eq(""), "candidate_hla"] = out.loc[out["candidate_hla"].eq(""), "hla"].map(normalize_hla)
    return out[out["candidate_peptide"].ne("")].copy()


def exact_overlaps(master: pd.DataFrame, corpus: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if corpus.empty:
        return pd.DataFrame()
    m = master.copy()
    c = corpus.copy()
    if c["public_hla"].astype(str).ne("").any():
        phla = m.merge(
            c[c["public_hla"].astype(str).ne("")],
            left_on=["candidate_peptide", "candidate_hla"],
            right_on=["public_peptide", "public_hla"],
            how="inner",
        )
        if not phla.empty:
            phla["overlap_type"] = "exact_peptide_hla"
            phla["near_similarity"] = 1.0
            rows.append(phla)
    pep = m.merge(c, left_on="candidate_peptide", right_on="public_peptide", how="inner")
    if not pep.empty:
        pep["overlap_type"] = "exact_peptide"
        pep["near_similarity"] = 1.0
        rows.append(pep)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True, sort=False)


def near_overlaps(master: pd.DataFrame, corpus: pd.DataFrame, threshold: float) -> pd.DataFrame:
    if corpus.empty:
        return pd.DataFrame()
    master_u = master[["candidate_id", "candidate_peptide", "candidate_hla", "label", "source_name"]].drop_duplicates()
    corpus_u = corpus[["public_peptide", "public_hla", "public_tool", "public_corpus_file"]].drop_duplicates()
    rows = []
    corpus_by_len: dict[int, pd.DataFrame] = {
        length: g for length, g in corpus_u.groupby(corpus_u["public_peptide"].str.len(), dropna=False)
    }
    for _, mr in master_u.iterrows():
        pep = mr["candidate_peptide"]
        if not pep:
            continue
        candidates = []
        for length in range(max(1, len(pep) - 1), len(pep) + 2):
            if length in corpus_by_len:
                candidates.append(corpus_by_len[length])
        if not candidates:
            continue
        sub = pd.concat(candidates, ignore_index=True)
        for _, cr in sub.iterrows():
            pub_pep = cr["public_peptide"]
            if pub_pep == pep:
                continue
            sim = near_similarity(pep, pub_pep)
            if sim >= threshold:
                rows.append(
                    {
                        "candidate_id": mr["candidate_id"],
                        "candidate_peptide": pep,
                        "candidate_hla": mr["candidate_hla"],
                        "label": mr["label"],
                        "source_name": mr["source_name"],
                        "public_peptide": pub_pep,
                        "public_hla": cr["public_hla"],
                        "public_tool": cr["public_tool"],
                        "public_corpus_file": cr["public_corpus_file"],
                        "public_row_index": "",
                        "public_training_label": "",
                        "peptide_column": "",
                        "hla_column": "",
                        "label_column": "",
                        "overlap_type": "near_peptide",
                        "near_similarity": sim,
                    }
                )
    return pd.DataFrame(rows)


def summarize(overlaps: pd.DataFrame, master_n: int, corpus_files: list[Path]) -> pd.DataFrame:
    if overlaps.empty:
        return pd.DataFrame(
            [
                {
                    "public_tool": infer_tool(path),
                    "public_corpus_file": str(path),
                    "n_candidate_overlaps": 0,
                    "candidate_overlap_fraction": 0.0,
                    "n_exact_peptide_hla": 0,
                    "n_exact_peptide": 0,
                    "n_near_peptide": 0,
                    "audit_status": "row_level_audit_completed_no_overlap",
                    "clean_comparator_allowed_after_row_audit": True,
                }
                for path in corpus_files
            ]
        )
    rows = []
    for (tool, path), g in overlaps.groupby(["public_tool", "public_corpus_file"], dropna=False):
        rows.append(
            {
                "public_tool": tool,
                "public_corpus_file": path,
                "n_candidate_overlaps": int(g["candidate_id"].nunique()),
                "candidate_overlap_fraction": float(g["candidate_id"].nunique() / max(1, master_n)),
                "n_exact_peptide_hla": int(g[g["overlap_type"].eq("exact_peptide_hla")]["candidate_id"].nunique()),
                "n_exact_peptide": int(g[g["overlap_type"].eq("exact_peptide")]["candidate_id"].nunique()),
                "n_near_peptide": int(g[g["overlap_type"].eq("near_peptide")]["candidate_id"].nunique()),
                "audit_status": "row_level_audit_completed_overlap_found",
                "clean_comparator_allowed_after_row_audit": False,
            }
        )
    seen_files = set(overlaps["public_corpus_file"].astype(str))
    for path in corpus_files:
        if str(path) not in seen_files:
            rows.append(
                {
                    "public_tool": infer_tool(path),
                    "public_corpus_file": str(path),
                    "n_candidate_overlaps": 0,
                    "candidate_overlap_fraction": 0.0,
                    "n_exact_peptide_hla": 0,
                    "n_exact_peptide": 0,
                    "n_near_peptide": 0,
                    "audit_status": "row_level_audit_completed_no_overlap",
                    "clean_comparator_allowed_after_row_audit": True,
                }
            )
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS).sort_values(["n_candidate_overlaps", "public_tool"], ascending=[False, True])


def requirements_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "public_tool": spec["public_tool"],
                "filename_tokens": ", ".join(spec["tokens"]),
                "required_key": spec["required_key"],
                "minimum_columns": "peptide plus HLA when available; source/split/label optional but recommended",
                "clean_pass_rule": spec["clean_pass_rule"],
                "drop_location": "project/data/public_training_corpora/<tool>_training.tsv",
            }
            for spec in PUBLIC_TOOL_SPECS
        ]
    )


def update_public_audit(output_root: Path, row_summary: pd.DataFrame) -> None:
    audit_path = output_root / "clean_neobench_public_tool_overlap_audit.tsv"
    if not audit_path.exists() or row_summary.empty:
        return
    audit = pd.read_csv(audit_path, sep="\t")
    if "method_name" not in audit.columns:
        return
    by_tool = row_summary.groupby("public_tool").agg(
        candidate_overlap_rows=("n_candidate_overlaps", "sum"),
        clean_allowed=("clean_comparator_allowed_after_row_audit", "all"),
        files=("public_corpus_file", lambda x: "; ".join(sorted(set(map(str, x))))),
    )
    for idx, row in audit.iterrows():
        method = str(row["method_name"])
        matches = [tool for tool in by_tool.index if tool.lower() in method.lower() or method.lower() in tool.lower()]
        if not matches:
            continue
        tool = matches[0]
        n_overlap = int(by_tool.loc[tool, "candidate_overlap_rows"])
        clean_allowed = bool(by_tool.loc[tool, "clean_allowed"]) and n_overlap == 0
        audit.loc[idx, "training_overlap_audit_status"] = (
            "row_level_audit_completed_no_overlap" if clean_allowed else "row_level_audit_completed_overlap_found"
        )
        audit.loc[idx, "training_overlap_audited_row_level"] = True
        audit.loc[idx, "candidate_overlap_rows"] = n_overlap
        audit.loc[idx, "candidate_overlap_fraction"] = n_overlap / max(1, len(pd.read_csv(output_root / "clean_neobench_master.tsv", sep="\t", usecols=["candidate_id"])))
        audit.loc[idx, "clean_comparator_allowed_after_audit"] = clean_allowed
        audit.loc[idx, "row_level_audit_files_detected"] = by_tool.loc[tool, "files"]
    write_tsv(audit, audit_path)


def write_report(output_root: Path, corpus_files: list[Path], overlaps: pd.DataFrame, summary: pd.DataFrame, warnings: list[str]) -> None:
    req = requirements_table()
    if corpus_files:
        files_md = "\n".join(f"- `{p}`" for p in corpus_files)
    else:
        files_md = "- No public training corpus files found in configured search directories."
    overlap_md = overlaps.head(50).to_markdown(index=False) if not overlaps.empty else "No overlap rows found or no corpora available."
    summary_md = summary.to_markdown(index=False) if not summary.empty else "No row-level audit summary available."
    req_md = req.to_markdown(index=False)
    warnings_md = "\n".join(f"- {w}" for w in warnings) if warnings else "- None."
    text = f"""# CLEAN-NeoBench Public Training Row-Level Audit

## Purpose

This report audits whether benchmark candidates overlap public pretrained tool training corpora at row level. Documentation-level provenance is not enough. A method becomes a cleaner comparator only when the relevant training corpus is present and exact peptide-HLA/peptide overlap is audited.

## Corpus Files Discovered

{files_md}

## Row-Level Summary

{summary_md}

## Example Overlap Rows

{overlap_md}

## Required Corpus Inputs

{req_md}

## Warnings

{warnings_md}

## Claim Rule

If no corpus file is available for a public pretrained tool, its training overlap status remains unresolved. If overlap is found, the method remains caveated. If a method-specific corpus is present and no overlap is found, it can be marked row-audited for this benchmark, subject to the corpus being complete and version-matched.
"""
    (output_root / "CLEAN_NEOBENCH_PUBLIC_TRAINING_ROW_AUDIT.md").write_text(text + "\n")


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    ensure_dir(output_root)

    master_path = output_root / "clean_neobench_master.tsv"
    master = pd.read_csv(master_path, sep="\t") if master_path.exists() else pd.DataFrame()
    master_norm = master_keys(master) if not master.empty else pd.DataFrame()
    corpus_files, warnings = discover_corpus_files(repo_root, args.public_corpus_root, args.max_file_mb)

    corpus_frames = []
    for path in corpus_files:
        df, w = read_corpus(path)
        warnings.extend(w)
        if not df.empty:
            corpus_frames.append(df)
    corpus = pd.concat(corpus_frames, ignore_index=True, sort=False) if corpus_frames else pd.DataFrame()
    exact = exact_overlaps(master_norm, corpus) if not master_norm.empty else pd.DataFrame()
    near = near_overlaps(master_norm, corpus, args.near_threshold) if not master_norm.empty else pd.DataFrame()
    overlaps = pd.concat([exact, near], ignore_index=True, sort=False) if not exact.empty or not near.empty else pd.DataFrame()
    summary = summarize(overlaps, len(master_norm), corpus_files) if corpus_files else pd.DataFrame()
    req = requirements_table()

    if overlaps.empty:
        overlaps = pd.DataFrame(columns=OVERLAP_COLUMNS)
    else:
        overlaps = overlaps[[c for c in OVERLAP_COLUMNS if c in overlaps.columns] + [c for c in overlaps.columns if c not in OVERLAP_COLUMNS]]
    if summary.empty:
        summary = pd.DataFrame(columns=SUMMARY_COLUMNS)
    write_tsv(overlaps, output_root / "clean_neobench_public_training_row_overlap.tsv")
    write_tsv(summary, output_root / "clean_neobench_public_training_row_overlap_summary.tsv")
    write_tsv(req, output_root / "clean_neobench_public_training_corpus_requirements.tsv")
    update_public_audit(output_root, summary)
    write_report(output_root, corpus_files, overlaps, summary, warnings)

    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("output_files", [])
    for name in [
        "clean_neobench_public_training_row_overlap.tsv",
        "clean_neobench_public_training_row_overlap_summary.tsv",
        "clean_neobench_public_training_corpus_requirements.tsv",
        "CLEAN_NEOBENCH_PUBLIC_TRAINING_ROW_AUDIT.md",
    ]:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_public_training_corpus_files": int(len(corpus_files)),
            "n_public_training_overlap_rows": int(len(overlaps)),
            "n_public_training_overlap_candidates": int(overlaps["candidate_id"].nunique()) if not overlaps.empty else 0,
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    stage_warnings = warnings[:]
    if not corpus_files:
        stage_warnings.append("No public training corpus files found; public pretrained tools remain unresolved.")
    update_manifest(
        output_root,
        "public_training_row_overlap_audit",
        {
            "n_corpus_files": int(len(corpus_files)),
            "n_overlap_rows": int(len(overlaps)),
            "n_overlap_candidates": int(overlaps["candidate_id"].nunique()) if not overlaps.empty else 0,
            "warnings": stage_warnings,
        },
    )
    print(f"[public-training-row-audit] corpus_files={len(corpus_files)} overlap_rows={len(overlaps)}")


if __name__ == "__main__":
    main()
