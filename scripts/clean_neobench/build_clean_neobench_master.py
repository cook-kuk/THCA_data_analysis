#!/usr/bin/env python3
"""Build the unified CLEAN-NeoBench candidate master table."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    KOREAN_HLA_ALLELES,
    MASTER_COLUMNS,
    ensure_dir,
    hla_allele_4digit,
    hla_gene,
    hla_supertype,
    normalize_empty,
    normalize_hla,
    peptide_cluster_fallback,
    read_table,
    stable_hash,
    update_manifest,
    write_tsv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="Output directory")
    return parser.parse_args()


def discover_inputs(repo: Path) -> list[Path]:
    roots = [
        repo / "project/results/p_neo_bayesian_2026_05_09",
        repo / "project/results/p_neo_bayesian_2026_05_09/curation_2026_05_09",
        repo / "project/results/p_neo_bayesian_2026_05_09/clean_neo_v0_2026_05_09",
        repo / "project/results/cross_neo_v0",
        repo / "project/results/cross_neo_v1",
        repo / "project/results/cross_neo_v1_lockdown",
    ]
    files: list[Path] = []
    for root in roots:
        if not root.exists():
            continue
        for suffix in ("*.tsv", "*.csv", "*.json", "*.md"):
            files.extend(root.glob(suffix))
    return sorted(set(files))


def first_present(row: pd.Series, names: list[str], default: str = "") -> str:
    for name in names:
        if name in row.index:
            val = normalize_empty(row.get(name))
            if val:
                return val
    return default


def normalize_label(value: object) -> object:
    s = normalize_empty(value)
    if s == "":
        return np.nan
    low = s.lower()
    if low in {"1", "true", "yes", "positive", "pos", "responder", "immunogenic"}:
        return 1
    if low in {"0", "false", "no", "negative", "neg", "nonresponder", "non-immunogenic"}:
        return 0
    try:
        f = float(s)
        if f in {0.0, 1.0}:
            return int(f)
    except Exception:
        pass
    return np.nan


def add_split_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["split_exact_phla"] = out.apply(
        lambda r: f"phla_fold_{stable_hash(str(r['peptide']) + '|' + str(r['hla_allele_4digit']))}",
        axis=1,
    )
    if "split_near_peptide_cluster" not in out.columns or out["split_near_peptide_cluster"].isna().all():
        out["split_near_peptide_cluster"] = out["peptide"].map(peptide_cluster_fallback)
    out["split_source_heldout"] = out["source_name"].replace("", "UNKNOWN")
    out["split_hla_heldout"] = out["hla_allele_4digit"].replace("", "UNKNOWN")
    out["split_supertype_heldout"] = out["hla_supertype"].replace("", "UNKNOWN")
    out["split_patient_heldout"] = out["patient_id"].replace("", "UNAVAILABLE")
    out["split_study_heldout"] = out["study_id"].replace("", "UNKNOWN")
    source_prev = out.groupby("source_name")["label"].mean(numeric_only=True).to_dict()
    source_npos = out.groupby("source_name")["label"].sum(numeric_only=True).to_dict()
    out["split_low_prevalence"] = out["source_name"].map(
        lambda s: "low_prevalence" if source_prev.get(s, 1.0) < 0.10 or source_npos.get(s, 99) < 10 else "not_low_prevalence"
    )
    out["split_korean_hla_focus"] = out["hla_allele_4digit"].isin(KOREAN_HLA_ALLELES).map({True: "korean_hla_focus", False: "not_korean_hla_focus"})
    return out


def normalize_from_cross_master(path: Path) -> pd.DataFrame:
    raw = read_table(path)
    rows = []
    for _, r in raw.iterrows():
        pep = first_present(r, ["peptide_mut", "mut_peptide", "peptide", "mt_peptide", "mutant_peptide"])
        hla = normalize_hla(first_present(r, ["hla", "allele", "mhc", "HLA_norm"]))
        sample_id = first_present(r, ["sample_id", "candidate_id"])
        source_name = first_present(r, ["study", "source", "source_name", "dataset"], "unknown")
        window = first_present(r, ["source_window_15aa", "source_window_30aa", "source_protein_window"])
        label = normalize_label(r.get("label", np.nan))
        rows.append(
            {
                "candidate_id": sample_id or f"CNB_{len(rows):06d}",
                "sample_id": sample_id,
                "patient_id": first_present(r, ["patient_id", "patient", "donor"]),
                "study_id": source_name,
                "source_name": source_name,
                "source_dataset": source_name,
                "cancer_type": "",
                "disease_context": "",
                "peptide": pep,
                "mut_peptide": pep,
                "wt_peptide": first_present(r, ["peptide_wt", "wt_peptide", "wildtype_peptide"]),
                "peptide_length": len(pep),
                "hla": hla,
                "hla_gene": hla_gene(hla),
                "hla_allele_4digit": hla_allele_4digit(hla),
                "hla_supertype": first_present(r, ["hla_supertype"]) or hla_supertype(hla),
                "gene": first_present(r, ["gene", "source_protein"]),
                "protein_id": first_present(r, ["protein_id", "source_protein"]),
                "mutation_id": first_present(r, ["mutation_id"]),
                "source_protein_window": window,
                "label": label,
                "label_type": "binary_immunogenicity" if pd.notna(label) else "",
                "mhc_class": "I",
                "exact_peptide_train_overlap": False,
                "exact_peptide_hla_train_overlap": False,
                "near_peptide_train_overlap": False,
                "source_protein_window_train_overlap": False,
                "study_train_overlap": False,
                "patient_train_overlap": False,
                "public_tool_training_overlap_any": bool(r.get("in_master", False)),
                "public_tool_training_overlap_detail": first_present(r, ["public_overlap_flags"]),
                "leakage_risk_level": "unknown",
                "split_near_peptide_cluster": first_present(r, ["near_peptide_cluster"]) or peptide_cluster_fallback(pep),
                "_original_split": first_present(r, ["split"]),
                "_strict_set_flag": bool(r.get("strict_set_flag", False)),
                "_source_file": str(path),
            }
        )
    return pd.DataFrame(rows)


def normalize_from_strict_seed(path: Path) -> pd.DataFrame:
    raw = read_table(path)
    rows = []
    for i, r in raw.iterrows():
        pep = first_present(r, ["peptide", "mut_peptide", "mt_peptide", "mutant_peptide"])
        hla = normalize_hla(first_present(r, ["hla", "allele", "mhc"]))
        label = normalize_label(r.get("label", np.nan))
        rows.append(
            {
                "candidate_id": f"CNB_STRICT_{i:05d}",
                "sample_id": f"CNB_STRICT_{i:05d}",
                "patient_id": "",
                "study_id": first_present(r, ["split"], "ITSNdb_strict"),
                "source_name": first_present(r, ["split"], "ITSNdb_strict"),
                "source_dataset": "clean_neo_v0_classI_strict_seed",
                "cancer_type": "",
                "disease_context": "",
                "peptide": pep,
                "mut_peptide": pep,
                "wt_peptide": "",
                "peptide_length": len(pep),
                "hla": hla,
                "hla_gene": hla_gene(hla),
                "hla_allele_4digit": hla_allele_4digit(hla),
                "hla_supertype": hla_supertype(hla),
                "gene": "",
                "protein_id": "",
                "mutation_id": "",
                "source_protein_window": "",
                "label": label,
                "label_type": "binary_immunogenicity" if pd.notna(label) else "",
                "mhc_class": "I",
                "public_tool_training_overlap_any": bool(r.get("in_master", False)),
                "public_tool_training_overlap_detail": ";".join(
                    [
                        f"mhcflurry={normalize_empty(r.get('mhcflurry_training_overlap_status'))}",
                        f"bigmhc={normalize_empty(r.get('bigmhc_training_overlap_status'))}",
                        f"wave8={normalize_empty(r.get('wave8_reference_status'))}",
                    ]
                ),
                "leakage_risk_level": "unknown",
                "split_near_peptide_cluster": peptide_cluster_fallback(pep),
                "_original_split": first_present(r, ["split"]),
                "_strict_set_flag": True,
                "_source_file": str(path),
            }
        )
    return pd.DataFrame(rows)


def synthetic_master() -> pd.DataFrame:
    peptides = [
        "SLYNTVATL",
        "GILGFVFTL",
        "NLVPMVATV",
        "LLFGYPVYV",
        "KLVVVGAVGV",
        "VVVGADGVGK",
        "GLATEKSRW",
        "AAGIGILTV",
        "ELAGIGILTV",
        "YLQPRTFLL",
        "FIAGLIAIV",
        "KIADYNYKL",
        "TPRVTGGGAM",
        "SPRWYFYYL",
        "RPHERNGFTV",
        "FPRPWLHGL",
        "ALWGPDPAAA",
        "IMDQVPFSV",
        "VLHDDLLEA",
        "YMDGTMSQV",
    ]
    hlas = ["HLA-A*02:01", "HLA-A*11:01", "HLA-A*24:02", "HLA-B*15:01"]
    rows = []
    for i, pep in enumerate(peptides):
        hla = hlas[i % len(hlas)]
        label = 1 if i in {0, 2, 4, 7, 12} else 0
        source = "synthetic_source_A" if i < 10 else "synthetic_source_B"
        rows.append(
            {
                "candidate_id": f"CNB_SYN_{i:05d}",
                "sample_id": f"CNB_SYN_{i:05d}",
                "patient_id": f"SYN_PAT_{i // 4:02d}",
                "study_id": source,
                "source_name": source,
                "source_dataset": "synthetic_fallback",
                "peptide": pep,
                "mut_peptide": pep,
                "wt_peptide": "",
                "peptide_length": len(pep),
                "hla": hla,
                "hla_gene": hla_gene(hla),
                "hla_allele_4digit": hla_allele_4digit(hla),
                "hla_supertype": hla_supertype(hla),
                "label": label,
                "label_type": "binary_immunogenicity",
                "mhc_class": "I",
                "source_protein_window": "",
                "public_tool_training_overlap_any": False,
                "public_tool_training_overlap_detail": "synthetic",
                "leakage_risk_level": "unknown",
                "split_near_peptide_cluster": peptide_cluster_fallback(pep),
                "_original_split": "train" if i < 12 else "test",
                "_strict_set_flag": True,
                "_source_file": "synthetic_fallback",
            }
        )
    return pd.DataFrame(rows)


def finalize_master(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in MASTER_COLUMNS:
        if col not in df.columns:
            df[col] = np.nan
    for col in ["candidate_id", "sample_id", "patient_id", "study_id", "source_name", "source_dataset", "peptide", "hla"]:
        df[col] = df[col].map(normalize_empty)
    df = df[df["peptide"].astype(str).str.len() > 0].copy()
    df = df[df["hla"].astype(str).str.len() > 0].copy()
    df["candidate_id"] = df["candidate_id"].where(df["candidate_id"].astype(str).str.len() > 0, [f"CNB_{i:06d}" for i in range(len(df))])
    df["candidate_id"] = df["candidate_id"].astype(str)
    duplicated = df["candidate_id"].duplicated(keep=False)
    if duplicated.any():
        counts: dict[str, int] = {}
        new_ids = []
        for cid in df["candidate_id"]:
            counts[cid] = counts.get(cid, 0) + 1
            new_ids.append(cid if counts[cid] == 1 else f"{cid}_dup{counts[cid]}")
        df["candidate_id"] = new_ids
    df["sample_id"] = df["sample_id"].where(df["sample_id"].astype(str).str.len() > 0, df["candidate_id"])
    df["label"] = df["label"].map(normalize_label)
    df["peptide_length"] = pd.to_numeric(df["peptide_length"], errors="coerce").fillna(df["peptide"].astype(str).str.len()).astype(int)
    df["hla"] = df["hla"].map(normalize_hla)
    df["hla_gene"] = df["hla"].map(hla_gene)
    df["hla_allele_4digit"] = df["hla"].map(hla_allele_4digit)
    df["hla_supertype"] = df["hla"].map(hla_supertype)
    df["mhc_class"] = df["mhc_class"].fillna("I").replace("", "I")
    df = add_split_columns(df)
    df = df[MASTER_COLUMNS + [c for c in df.columns if c.startswith("_")]]
    return df.reset_index(drop=True)


def main() -> None:
    args = parse_args()
    repo = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    ensure_dir(output_root)
    warnings: list[str] = []
    input_files = discover_inputs(repo)
    pieces = []

    cross_master = repo / "project/results/cross_neo_v0/master_table.tsv"
    strict_seed = repo / "project/results/p_neo_bayesian_2026_05_09/clean_neo_v0_2026_05_09/clean_neo_v0_classI_strict_seed.tsv"
    if cross_master.exists():
        pieces.append(normalize_from_cross_master(cross_master))
    else:
        warnings.append(f"missing preferred master: {cross_master}")
    if not pieces and strict_seed.exists():
        pieces.append(normalize_from_strict_seed(strict_seed))
    elif strict_seed.exists():
        # Keep strict-only rows only if they are not already represented by peptide-HLA.
        strict = normalize_from_strict_seed(strict_seed)
        existing_keys = set(zip(pieces[0]["peptide"], pieces[0]["hla"]))
        strict = strict[~strict.apply(lambda r: (r["peptide"], r["hla"]) in existing_keys, axis=1)]
        if len(strict):
            pieces.append(strict)

    if not pieces:
        warnings.append("no real candidate table found; using synthetic fallback")
        pieces.append(synthetic_master())

    master = finalize_master(pd.concat(pieces, ignore_index=True, sort=False))
    write_tsv(master.drop(columns=[c for c in master.columns if c.startswith("_")], errors="ignore"), output_root / "clean_neobench_master.tsv")
    update_manifest(
        output_root,
        "build_clean_neobench_master",
        {
            "input_files_discovered": [str(p) for p in input_files],
            "input_files_used": sorted(master["_source_file"].dropna().astype(str).unique().tolist()) if "_source_file" in master.columns else [],
            "n_candidates": int(len(master)),
            "n_labeled": int(master["label"].notna().sum()),
            "n_positive": int(pd.to_numeric(master["label"], errors="coerce").fillna(0).sum()),
            "source_counts": master["source_name"].value_counts().to_dict(),
            "warnings": warnings,
        },
    )
    print(f"[clean-neobench-master] rows={len(master)} output={output_root / 'clean_neobench_master.tsv'}")


if __name__ == "__main__":
    main()
