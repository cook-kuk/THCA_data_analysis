#!/usr/bin/env python3
"""Build pMHC geometry features from existing ESMFold pseudo-complex outputs."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cross_neo_v0_common import INPUT, OUT, ensure_dirs, load_master


GEOM_BASE = [
    "mean_pLDDT_peptide",
    "min_pLDDT_peptide",
    "mean_pLDDT_HLA",
    "anchor_pLDDT",
    "interface_contacts_8A",
    "interface_contacts_10A",
    "n_buried_residues_8A",
    "mean_min_pep_to_hla_CA_dist",
    "max_min_pep_to_hla_CA_dist",
    "radius_of_gyration_peptide",
    "end_to_end_CA_dist",
    "peptide_helicity_proxy",
]


def main() -> None:
    ensure_dirs()
    master = load_master()
    strict_geom = pd.read_csv(INPUT / "wave8b_strict_structure_2026_05_09/strict_esmfold_features.tsv", sep="\t")
    geom = master[["sample_id", "peptide_mut", "hla", "label", "strict_set_flag"]].copy()
    geom = geom.merge(
        strict_geom[["peptide", "hla", *GEOM_BASE]],
        left_on=["peptide_mut", "hla"],
        right_on=["peptide", "hla"],
        how="left",
    )
    for c in GEOM_BASE:
        geom[c] = geom[c].astype(float)
    pep_len = geom["peptide_mut"].astype(str).str.len().clip(lower=1).astype(float)
    geom["peptide_bulge_proxy"] = geom["radius_of_gyration_peptide"] / pep_len
    geom["anchor_residue_burial_proxy"] = geom["n_buried_residues_8A"] / pep_len
    geom["tcr_facing_exposure_proxy"] = 1.0 / (1.0 + geom["interface_contacts_8A"].fillna(0))
    geom["contact_density"] = geom["interface_contacts_10A"] / pep_len
    geom["hla_groove_contact_distribution_proxy"] = geom["max_min_pep_to_hla_CA_dist"] - geom["mean_min_pep_to_hla_CA_dist"]
    geom["mutant_wt_surface_delta"] = np.nan
    geom["ensemble_variance"] = np.nan
    geom["structure_missing"] = geom["mean_pLDDT_peptide"].isna().astype(int)
    geom["structure_low_confidence"] = ((geom["mean_pLDDT_peptide"].fillna(0) < 50) | (geom["structure_missing"] == 1)).astype(int)
    out_cols = [
        "sample_id",
        "peptide_mut",
        "hla",
        "label",
        "strict_set_flag",
        *GEOM_BASE,
        "peptide_bulge_proxy",
        "anchor_residue_burial_proxy",
        "tcr_facing_exposure_proxy",
        "contact_density",
        "hla_groove_contact_distribution_proxy",
        "mutant_wt_surface_delta",
        "ensemble_variance",
        "structure_missing",
        "structure_low_confidence",
    ]
    geom[out_cols].to_csv(OUT / "structure_geometry_features.tsv", sep="\t", index=False)
    lines = [
        "# Structure Geometry Report",
        "",
        f"Rows with ESMFold pseudo-complex geometry: {int((geom['structure_missing'] == 0).sum())} / {len(geom)}.",
        "Source: existing strict ESMFold peptide + GGGGS + HLA-pseudo outputs.",
        "Boltz/AF3 structures were not found locally and were skipped.",
        "SaProt/ProSST structure-token embeddings were not generated in v0.",
        "These are proxy features, not physical pMHC structure claims.",
    ]
    (OUT / "structure_geometry_report.md").write_text("\n".join(lines) + "\n")
    print(f"[structure] wrote {OUT / 'structure_geometry_features.tsv'}")


if __name__ == "__main__":
    main()
