#!/usr/bin/env python3
"""Build compact CLEAN-Neo v0 seed tables from locked wave artifacts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
STRICT = ROOT / "wave8b_strict_structure_2026_05_09"
CURATION = ROOT / "curation_2026_05_09"
OUT = ROOT / "clean_neo_v0_2026_05_09"


KEEP = [
    "peptide",
    "hla",
    "label",
    "split",
    "in_master",
    "tcr_motif_score",
    "self_exact_match",
    "Structure_LR",
    "MHCflurry",
    "BigMHC_IM",
    "PRIME",
    "NetMHCpan_4.1",
    "Wave8_TCR_SelfSim_full",
    "Wave8_TCR_SelfSim_no_exact",
    "Wave8_TCR_motif_only",
    "ESM2_Bayesian",
    "ESMFold_3D_CV_score",
    "mean_pLDDT_peptide",
    "min_pLDDT_peptide",
    "interface_contacts_8A",
    "interface_contacts_10A",
    "mean_min_pep_to_hla_CA_dist",
    "radius_of_gyration_peptide",
    "peptide_helicity_proxy",
]


def main() -> None:
    OUT.mkdir(exist_ok=True)
    strict = pd.read_csv(STRICT / "strict_esmfold_merged_predictions.tsv", sep="\t")
    strict = strict[[c for c in KEEP if c in strict.columns]].copy()
    strict["benchmark_tier"] = "classI_strict_no_overlap_no_tcr_self_exact"
    strict["mhcflurry_training_overlap_status"] = "unresolved_public_pretrained"
    strict["bigmhc_training_overlap_status"] = "unresolved_public_pretrained"
    strict["wave8_reference_status"] = "exact_tcr_self_removed"
    strict["structure_status"] = "esmfold_pseudo_complex_proxy"
    strict.to_csv(OUT / "clean_neo_v0_classI_strict_seed.tsv", sep="\t", index=False)

    all_methods = pd.read_csv(CURATION / "class1_class2_all_method_comparison.tsv", sep="\t")
    all_methods.to_csv(OUT / "locked_classI_classII_method_comparison.tsv", sep="\t", index=False)

    method = pd.read_csv(STRICT / "strict_esmfold_method_comparison.tsv", sep="\t")
    method.to_csv(OUT / "locked_classI_strict_method_comparison.tsv", sep="\t", index=False)

    sources = pd.DataFrame(
        [
            ("IEDB", "peptide immunogenicity / MHC ligand source", "needed", "time-split and peptide-HLA overlap"),
            ("CEDAR", "cancer epitope resource", "needed", "public-tool training contamination risk"),
            ("TESLA", "blinded neoantigen consortium benchmark", "needed", "small but high-value positives"),
            ("NEPdb", "neoepitope database", "needed", "public immunogenicity overlap"),
            ("MHCflurry", "public pretrained class-I presentation", "partial_local_provenance", "treat as upper-bound comparator"),
            ("NetMHCpan/NetMHCIIpan", "DTU public predictors", "needed", "versioned training overlap hard to fully exclude"),
            ("BigMHC", "presentation + immunogenicity transfer", "needed", "public release overlap audit"),
            ("PRIME/MixMHCpred", "presentation/TCR-recognition baseline", "needed", "public training overlap audit"),
            ("IMMREP23/25", "TCR-pMHC challenge", "needed", "unseen pMHC split design"),
        ],
        columns=["source", "role", "local_status", "audit_use"],
    )
    sources.to_csv(OUT / "public_corpus_overlap_todo.tsv", sep="\t", index=False)

    readme = f"""# CLEAN-Neo v0 Seed Pack

Generated from locked wave artifacts.

## Files

- `clean_neo_v0_classI_strict_seed.tsv`: strict class-I rows for immediate v0 modeling; n={len(strict)}, positives={int(strict['label'].sum())}.
- `locked_classI_strict_method_comparison.tsv`: same-row strict comparator table.
- `locked_classI_classII_method_comparison.tsv`: full current method table; class II is intentionally marked not comparable.
- `public_corpus_overlap_todo.tsv`: contamination-audit source checklist.

## Immediate Rule

Use `Structure_LR` as the clean local anchor, keep `MHCflurry`/`BigMHC` caveated, and do not treat `ESMFold_3D_CV_score` as external evidence.
"""
    (OUT / "README.md").write_text(readme)


if __name__ == "__main__":
    main()
