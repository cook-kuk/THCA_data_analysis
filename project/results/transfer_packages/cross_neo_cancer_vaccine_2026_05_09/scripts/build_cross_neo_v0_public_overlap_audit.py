#!/usr/bin/env python3
"""Public/local overlap audit scaffold for CROSS-Neo v0.

This does not pretend to have every public tool's hidden training corpus. It
records what is locally auditable now, and labels the rest unresolved.
"""

from __future__ import annotations

import pandas as pd

from cross_neo_v0_common import OUT, ensure_dirs, load_master, seq_similarity


PUBLIC_SOURCES = [
    ("MHCflurry", "unresolved_public_pretrained", "IEDB/MS-ligand/affinity training corpus not row-audited locally"),
    ("NetMHCpan", "unresolved_public_pretrained", "BA/EL public training corpus not row-audited locally"),
    ("BigMHC", "unresolved_public_pretrained", "presentation/immunogenicity release exists but not downloaded into this audit"),
    ("PRIME", "unresolved_public_pretrained", "public immunogenicity training set not row-audited locally"),
    ("MixMHCpred", "unresolved_public_pretrained", "public ligand training corpus not row-audited locally"),
    ("CEDAR", "local_train_pool_present", "exact/near overlap against local train pool auditable"),
    ("TESLA", "local_train_pool_present", "exact/near overlap against local train pool auditable"),
    ("NEPdb", "local_train_pool_present", "exact/near overlap against local train pool auditable"),
]


def max_sim_with_source(pep: str, train: pd.DataFrame) -> tuple[float, str]:
    best = 0.0
    best_src = ""
    for _, r in train.iterrows():
        s = seq_similarity(pep, str(r["peptide_mut"]))
        if s > best:
            best = s
            best_src = str(r["study"])
    return best, best_src


def main() -> None:
    ensure_dirs()
    master = load_master()
    strict = master[master["strict_set_flag"].astype(bool)].copy()
    train = master[master["split"] == "train"].copy()
    pair_set = set(zip(train["peptide_mut"], train["hla"]))
    pep_set = set(train["peptide_mut"])
    rows = []
    for _, r in strict.iterrows():
        pep = str(r["peptide_mut"])
        hla = str(r["hla"])
        near, near_src = max_sim_with_source(pep, train)
        rows.append(
            {
                "sample_id": r["sample_id"],
                "peptide_mut": pep,
                "hla": hla,
                "label": int(r["label"]),
                "local_train_exact_peptide_hla": int((pep, hla) in pair_set),
                "local_train_exact_peptide": int(pep in pep_set),
                "local_train_near_similarity_max": near,
                "local_train_near_similarity_source": near_src,
                "local_train_near_hit_ge_0p75": int(near >= 0.75),
                "public_predictor_overlap_status": "unresolved_for_MHCflurry_NetMHCpan_BigMHC_PRIME_MixMHCpred",
            }
        )
    audit = pd.DataFrame(rows)
    audit.to_csv(OUT / "public_overlap_audit.tsv", sep="\t", index=False)
    src = pd.DataFrame(PUBLIC_SOURCES, columns=["source", "status", "note"])
    src.to_csv(OUT / "public_overlap_source_status.tsv", sep="\t", index=False)
    lines = [
        "# Public / Local Overlap Audit",
        "",
        "This audit is intentionally conservative.",
        "",
        f"Strict rows audited: {len(audit)}",
        f"Local exact peptide-HLA hits: {int(audit['local_train_exact_peptide_hla'].sum())}",
        f"Local exact peptide hits: {int(audit['local_train_exact_peptide'].sum())}",
        f"Local near hits similarity >=0.75: {int(audit['local_train_near_hit_ge_0p75'].sum())}",
        "",
        "Public predictor training overlap remains unresolved unless each public training corpus is downloaded and row-audited.",
        "",
        src.to_markdown(index=False),
    ]
    (OUT / "public_overlap_audit.md").write_text("\n".join(lines) + "\n")
    print(f"[public-overlap] strict={len(audit)} near>=0.75={int(audit['local_train_near_hit_ge_0p75'].sum())}")


if __name__ == "__main__":
    main()
