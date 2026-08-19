#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from neoimmune_common import HUB_ROOT, ensure_run_dir, safe_read_table, write_md, write_tsv


PUBLIC_SOURCES = [
    {
        "dataset": "IMPROVE neoepitope immunogenicity",
        "public_status": "open paper + GitHub; local ingest present",
        "patient_level": "yes in original publication; local ingest currently lacks patient_id",
        "what_it_has": "T-cell recognition labels, peptide/HLA, broad cancer patient cohort",
        "best_use": "immunogenicity benchmark and model training/evaluation with leakage controls",
        "local_path": "/data/neoantigen_vaccine_hub/data_processed/ingest_improve.tsv",
        "source_url": "https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1360281/full",
        "access_note": "Paper states data/scripts are on GitHub; use original tables to recover patient_id if not in local ingest.",
    },
    {
        "dataset": "TESLA public data",
        "public_status": "Synapse public access page; some raw files may require Synapse/AWS permissions",
        "patient_level": "yes",
        "what_it_has": "somatic variants, HLA calls, patient IDs, validated neoantigen candidate lists",
        "best_use": "patient-level benchmark, top-N ranking, variant-to-peptide audit",
        "local_path": "/data/neoantigen_vaccine_hub/data_processed/tesla_with_vdjdb_distance.tsv",
        "source_url": "https://www.synapse.org/TESLA_public",
        "access_note": "Synapse page lists Somatic variants, HLA Calls, Patient ID, and validated candidates.",
    },
    {
        "dataset": "NEPdb",
        "public_status": "open curated database; local ingest present",
        "patient_level": "yes in local ingest",
        "what_it_has": "neoantigen/neopeptide, HLA, cancer type, T-cell/clinical response labels, patient_id in local ingest",
        "best_use": "patient-level immunogenicity benchmark with caution about curation heterogeneity",
        "local_path": "/data/neoantigen_vaccine_hub/data_processed/ingest_nepdb.tsv",
        "source_url": "http://nep.whu.edu.cn",
        "access_note": "Local ingest has patient_id; curated positives/negatives need leakage/source audit.",
    },
    {
        "dataset": "CEDAR / IEDB cancer T-cell assays",
        "public_status": "open database-derived assay labels; local CEDAR partial present",
        "patient_level": "usually weak or absent",
        "what_it_has": "peptide, HLA, assay labels, source antigen/study metadata",
        "best_use": "assay-level immunogenicity labels and external stress test",
        "local_path": "/data/neoantigen_vaccine_hub/data_raw/cedar/cedar.tsv",
        "source_url": "https://www.iedb.org/",
        "access_note": "Good for labels; weaker for patient-level top-N unless original study IDs can be reconstructed.",
    },
    {
        "dataset": "dbPepNeo",
        "public_status": "open curated database",
        "patient_level": "limited",
        "what_it_has": "experimentally verified and MS-screened human tumor neoantigen peptides",
        "best_use": "positive evidence, MS/presentation support, external validation",
        "local_path": "/data/neoantigen_vaccine_hub/data_processed/ingest_dbpepneo.tsv",
        "source_url": "https://academic.oup.com/database/article/doi/10.1093/database/baaa004/5747759",
        "access_note": "Curated peptide evidence; not a full hospital-style patient table.",
    },
    {
        "dataset": "caAtlas cancer antigen atlas",
        "public_status": "open portal/public atlas",
        "patient_level": "sample-level immunopeptidomics; not full mutation/immunogenicity labels",
        "what_it_has": "MS-detected MHC-bound peptides across published cancer immunopeptidomics datasets",
        "best_use": "presentation evidence prior and normal/tumor peptide support",
        "local_path": "",
        "source_url": "http://www.zhang-lab.org/caatlas/",
        "access_note": "Use as presentation atlas, not T-cell immunogenicity labels.",
    },
    {
        "dataset": "Ligand.MHC Atlas",
        "public_status": "open recent large-scale immunopeptidomics atlas",
        "patient_level": "many patient tissue samples; mutation labels not guaranteed",
        "what_it_has": "HLA-I/HLA-II ligands across thousands of immunopeptidomic samples",
        "best_use": "presentation prior, HLA allele coverage, tumor/normal ligand context",
        "local_path": "",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12239622/",
        "access_note": "Use for MS ligand evidence; does not replace patient mutation/immunogenicity table.",
    },
    {
        "dataset": "CPTAC / GDC proteogenomics",
        "public_status": "public and/or controlled-access components",
        "patient_level": "yes",
        "what_it_has": "patient-level tumor genomics, proteomics, clinical context across cancer types",
        "best_use": "derive patient candidate universe, expression/protein context, tumor type validation",
        "local_path": "",
        "source_url": "https://gdc.cancer.gov/about-gdc/contributed-genomic-data-cancer-research/clinical-proteomic-tumor-analysis-consortium-cptac",
        "access_note": "Not a T-cell assay benchmark by itself; access and harmonization required.",
    },
]


def local_counts(path: str) -> dict[str, object]:
    if not path:
        return {"local_exists": False, "rows": 0, "patient_id_rows": 0, "unique_patients": 0, "hla_rows": 0}
    p = Path(path)
    if not p.exists():
        return {"local_exists": False, "rows": 0, "patient_id_rows": 0, "unique_patients": 0, "hla_rows": 0}
    try:
        sep = "\t" if p.suffix == ".tsv" else ","
        df = pd.read_csv(p, sep=sep, low_memory=False)
    except Exception:
        return {"local_exists": True, "rows": -1, "patient_id_rows": 0, "unique_patients": 0, "hla_rows": 0}
    pcols = [c for c in df.columns if "patient" in c.lower() or "subject" in c.lower() or "donor" in c.lower()]
    hcols = [c for c in df.columns if "hla" in c.lower()]
    patient_id_rows = 0
    unique_patients = 0
    if pcols:
        s = df[pcols[0]].dropna().astype(str)
        s = s[s.ne("") & s.ne("unknown_patient")]
        patient_id_rows = int(s.shape[0])
        unique_patients = int(s.nunique())
    hla_rows = int(df[hcols[0]].notna().sum()) if hcols else 0
    return {
        "local_exists": True,
        "rows": int(df.shape[0]),
        "patient_id_rows": patient_id_rows,
        "unique_patients": unique_patients,
        "hla_rows": hla_rows,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    rows = []
    for src in PUBLIC_SOURCES:
        r = dict(src)
        r.update(local_counts(src["local_path"]))
        rows.append(r)
    inv = pd.DataFrame(rows)
    write_tsv(inv, outdir / "reports" / "public_patient_level_data_inventory.tsv")

    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    real_mask = canon["patient_id"].fillna("unknown_patient").astype(str).ne("unknown_patient") if "patient_id" in canon else pd.Series(False, index=canon.index)
    by_source = (
        canon.assign(real_patient=real_mask)
        .groupby("dataset_source", dropna=False)
        .agg(rows=("candidate_id", "size"), real_patient_rows=("real_patient", "sum"), unique_patients=("patient_id", lambda x: x.fillna("unknown_patient").astype(str).loc[lambda s: s.ne("unknown_patient")].nunique()))
        .reset_index()
        .sort_values(["unique_patients", "real_patient_rows", "rows"], ascending=False)
    )
    write_tsv(by_source, outdir / "reports" / "canonical_patient_coverage_by_source.tsv")

    md = [
        "# Public Patient-Level Data Hunt",
        "",
        "## Bottom Line",
        "Yes. Open/public patient-level neoantigen data exist. The strongest immediate local sources are NEPdb and the integrated benchmark-ready neoantigen table; TESLA and IMPROVE are the next recovery targets for richer original patient-level fields.",
        "",
        "## What Changed",
        f"- Canonical table now has {int(real_mask.sum()):,} real patient-ID rows.",
        f"- Unique real non-placeholder patients: {canon.loc[real_mask, 'patient_id'].astype(str).nunique():,}.",
        "- Patient-level top-N metrics are now meaningful as a retrospective public benchmark scaffold.",
        "- This still does not prove prospective hospital utility or clinical vaccine efficacy.",
        "",
        "## Dataset Inventory",
        inv.to_markdown(index=False),
        "",
        "## Canonical Coverage By Source",
        by_source.to_markdown(index=False),
        "",
        "## Immediate Data Actions",
        "1. Recover original IMPROVE GitHub tables to restore patient IDs absent from current `ingest_improve.tsv`.",
        "2. Use TESLA Synapse patient/HLA/variant/validated-candidate tables as the main patient top-N benchmark.",
        "3. Keep NEPdb as the quickest open patient-ID benchmark already local.",
        "4. Use caAtlas/Ligand.MHC/CPTAC as presentation/proteogenomics support, not immunogenicity labels.",
        "5. Update clean/prod leaderboards with per-source and leave-patient-out splits after TESLA/IMPROVE recovery.",
        "",
        "## Claim Boundary",
        "- Public patient-level retrospective ranking scaffold: allowed.",
        "- Prospective hospital-grade performance: not yet.",
        "- Presentation confirmed: only with MS evidence.",
        "- Immunogenicity confirmed: only with T-cell assay labels.",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "public_patient_level_data_hunt.md")
    print(outdir / "reports" / "public_patient_level_data_hunt.md")


if __name__ == "__main__":
    main()
