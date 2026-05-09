#!/usr/bin/env python3
"""CROSS-Neo v1 public/local overlap audit with explicit missing-data manifest."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from cross_neo_v1_common import REPO, V1, ensure_v1_dirs, load_master, seq_similarity


LOCAL_SOURCE_MAP = {
    "CEDAR": "master_table train source",
    "NEPdb": "master_table train source",
    "TESLA": "master_table train source TESLA_mmc4/TESLA_mmc7_validation",
}

MISSING_PUBLIC = [
    ("IEDB", "https://www.iedb.org/database_export_v3.php", "peptide, HLA allele, assay label", "project/data/external/iedb/"),
    ("BigMHC release", "BigMHC GitHub/release training data", "peptide, HLA, immunogenicity/presentation label", "project/data/external/bigmhc/"),
    ("PRIME training data", "PRIME/PRIME2 public training data", "peptide, HLA, immunogenicity label", "project/data/external/prime/"),
    ("MHCflurry training data", "MHCflurry downloads / local package data", "peptide, allele, affinity/presentation source", "project/data/external/mhcflurry_training/"),
    ("NetMHCpan public references", "NetMHCpan BA/EL public references", "peptide, allele, BA/EL label", "project/data/external/netmhcpan_training/"),
    ("MixMHCpred ligand corpus", "MixMHCpred public ligand training corpus", "peptide, allele/source sample", "project/data/external/mixmhcpred_training/"),
]


def best_near(pep: str, corpus: pd.DataFrame) -> tuple[float, str]:
    best = 0.0
    best_hla = ""
    for _, r in corpus.iterrows():
        s = seq_similarity(pep, str(r["peptide_mut"]))
        if s > best:
            best = s
            best_hla = str(r["hla"])
    return best, best_hla


def existing_hints() -> dict[str, str]:
    hints = {}
    patterns = {
        "IEDB": ["project/papers_hub_2026_05_04/neoantigen_hub_data/*iedb*"],
        "BigMHC": ["project/results/p_neo_bayesian_2026_05_09/wave3_bigmhc/*", "project/papers_hub_2026_05_04/neoantigen_hub_data/*bigmhc*"],
        "PRIME": ["project/results/p_neo_bayesian_2026_05_09/wave3_prime/*", "project/results/p_neo_bayesian_2026_05_09/wave11/_prime_work"],
        "MHCflurry": ["project/results/p_neo_bayesian_2026_05_09/wave5b/mhcflurry_raw.parquet", "project/results/p_neo_bayesian_2026_05_09/wave3_mhcflurry/*"],
        "NetMHCpan": ["project/results/p_neo_bayesian_2026_05_09/wave3_netmhcpan/*"],
        "TESLA_raw": ["project/results/proteogenomic_v1/_raw/wells2020_TESLA/*"],
    }
    for name, pats in patterns.items():
        hits = []
        for pat in pats:
            hits.extend(str(p) for p in REPO.glob(pat))
        hints[name] = ";".join(sorted(hits)[:20])
    return hints


def main() -> None:
    ensure_v1_dirs()
    master = load_master()
    strict = master[master["strict_set_flag"].astype(bool)].copy()
    local = master[master["split"] == "train"].copy()
    rows = []
    for corpus_name, note in LOCAL_SOURCE_MAP.items():
        if corpus_name == "TESLA":
            corpus = local[local["study"].str.startswith("TESLA", na=False)].copy()
        else:
            corpus = local[local["study"] == corpus_name].copy()
        pair_set = set(zip(corpus["peptide_mut"].astype(str), corpus["hla"].astype(str)))
        pep_set = set(corpus["peptide_mut"].astype(str))
        for _, r in strict.iterrows():
            pep, hla = str(r["peptide_mut"]), str(r["hla"])
            near, near_hla = best_near(pep, corpus) if len(corpus) else (0.0, "")
            rows.append(
                {
                    "audit_corpus": corpus_name,
                    "audit_status": "local_row_audited" if len(corpus) else "missing",
                    "corpus_note": note,
                    "sample_id": r["sample_id"],
                    "peptide_mut": pep,
                    "hla": hla,
                    "label": int(r["label"]),
                    "exact_peptide_hla": int((pep, hla) in pair_set),
                    "exact_peptide": int(pep in pep_set),
                    "near_peptide_similarity": near,
                    "near_peptide_hla": near_hla,
                    "near_hit_ge_0p75": int(near >= 0.75),
                    "source_window_auditable": 0,
                    "mutation_wt_pair_auditable": 0,
                    "tcr_motif_auditable": 0,
                }
            )
    audit = pd.DataFrame(rows)
    audit.to_csv(V1 / "public_overlap_audit_v1.tsv", sep="\t", index=False)
    hints = existing_hints()
    miss = []
    for dataset, url, cols, target in MISSING_PUBLIC:
        key = dataset.split()[0].replace("release", "")
        present_hint = hints.get(key, "")
        miss.append(
            {
                "dataset": dataset,
                "url_or_search_term": url,
                "expected_columns": cols,
                "local_target_path": target,
                "local_hint_files": present_hint,
                "status": "needs_training_corpus_row_audit" if not present_hint else "local_prediction_or_status_files_exist_but_training_corpus_not_confirmed",
            }
        )
    manifest = pd.DataFrame(miss)
    manifest.to_csv(V1 / "public_overlap_download_manifest_needed.tsv", sep="\t", index=False)
    summary = audit.groupby("audit_corpus")[["exact_peptide_hla", "exact_peptide", "near_hit_ge_0p75"]].sum().reset_index()
    lines = [
        "# CROSS-Neo v1 Public / Local Overlap Audit",
        "",
        "Local CEDAR/NEPdb/TESLA overlap is row-audited from the current master table. Public pretrained predictor training corpora remain unresolved unless their training rows are downloaded.",
        "",
        "## Local Audit Summary",
        summary.to_markdown(index=False),
        "",
        "## Missing / Unresolved Public Training Corpora",
        manifest.to_markdown(index=False),
        "",
        "No public predictor scores are used as model features in v1.",
    ]
    (V1 / "public_overlap_audit_v1.md").write_text("\n".join(lines) + "\n")
    print(f"[v1-overlap] rows={len(audit)} unresolved={len(manifest)}")


if __name__ == "__main__":
    main()
