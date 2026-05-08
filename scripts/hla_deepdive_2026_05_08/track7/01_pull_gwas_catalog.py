#!/usr/bin/env python3
"""Track 7 — pull GWAS Catalog associations for thyroid autoimmunity traits.

Strategy:
1. For each MONDO/EFO short form with disease/trait downloads support, hit the
   `/gwas/api/search/downloads` endpoint and save TSV.
2. For traits where shortForm-driven download fails, paginate the REST
   `/efoTraits/{id}/associations` endpoint with the `associationByEfoTrait`
   projection and collect SNP-level rows.
3. Merge into one harmonised TSV.

Boundary: thyroid autoimmunity / thyroid function only — not cancer.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track7_gwas_catalog_mhc")
RAW = OUT / "raw"
TAB = OUT / "tables"
RAW.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

# Trait registry: shortForm, friendly label, expected mode (download or rest)
TRAITS = [
    ("MONDO_0005364", "Graves disease",                "download"),
    ("MONDO_0007699", "Hashimoto thyroiditis",         "download"),
    ("MONDO_0005623", "Autoimmune thyroid disease",    "download"),
    ("MONDO_0005420", "Hypothyroidism",                "download"),
    ("MONDO_0004425", "Hyperthyroidism",               "download"),
    ("EFO_0004296",   "Thyroid function (TSH/T4)",     "rest"),
    ("EFO_0010050",   "Thyroglobulin measurement",     "rest"),
    ("EFO_0021632",   "L-Thyroxine measurement",       "rest"),
]

DOWNLOAD_URL = "https://www.ebi.ac.uk/gwas/api/search/downloads"
REST_BASE = "https://www.ebi.ac.uk/gwas/rest/api"

DOWNLOAD_COLS = [
    "DATE ADDED TO CATALOG", "PUBMEDID", "FIRST AUTHOR", "DATE", "JOURNAL",
    "LINK", "STUDY", "DISEASE/TRAIT", "INITIAL SAMPLE SIZE",
    "REPLICATION SAMPLE SIZE", "REGION", "CHR_ID", "CHR_POS",
    "REPORTED GENE(S)", "MAPPED_GENE", "STRONGEST SNP-RISK ALLELE",
    "SNPS", "MERGED", "SNP_ID_CURRENT", "CONTEXT", "INTERGENIC",
    "RISK ALLELE FREQUENCY", "P-VALUE", "PVALUE_MLOG",
    "P-VALUE (TEXT)", "OR or BETA", "95% CI (TEXT)", "PLATFORM [SNPS PASSING QC]",
    "CNV", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "STUDY ACCESSION",
    "GENOTYPING TECHNOLOGY",
]


def pull_download(short: str, label: str, attempts: int = 5) -> pd.DataFrame | None:
    params = {
        "q": f'shortForm:"{short}"',
        "efo": "true",
        "facet": "association",
        "pvalfilter": "",
        "orfilter": "",
        "betafilter": "",
        "datefilter": "",
        "genomicfilter": "",
        "genotypingfilter[]": "",
        "traitfilter[]": "",
        "dateaddedfilter": "",
    }
    for k in range(attempts):
        try:
            r = requests.get(DOWNLOAD_URL, params=params, timeout=120)
        except Exception as e:
            print(f"[{short}] attempt {k} exception {e}", flush=True)
            time.sleep(2 * (k + 1))
            continue
        if r.status_code == 200 and r.text.startswith("DATE"):
            break
        print(f"[{short}] attempt {k} status={r.status_code}", flush=True)
        time.sleep(2 * (k + 1))
    else:
        print(f"[{short}/{label}] download exhausted retries", flush=True)
        return None
    p = RAW / f"{short}_associations.tsv"
    p.write_text(r.text)
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    df["__source_id"] = short
    df["__source_label"] = label
    print(f"[{short}/{label}] download mode {len(df)} rows", flush=True)
    return df


def pull_rest(short: str, label: str, page_size: int = 100, max_pages: int = 50) -> pd.DataFrame | None:
    rows = []
    page = 0
    while page < max_pages:
        url = f"{REST_BASE}/efoTraits/{short}/associations"
        r = requests.get(
            url,
            params={"projection": "associationByEfoTrait", "size": page_size, "page": page},
            timeout=60,
        )
        if r.status_code != 200:
            print(f"[{short}/{label}] rest page {page} failed status={r.status_code}", flush=True)
            break
        d = r.json()
        embed = d.get("_embedded", {}).get("associations", [])
        if not embed:
            break
        for a in embed:
            study = a.get("study", {}) or {}
            ancestries = study.get("ancestries", []) or []
            anc_str = "; ".join(
                f"{x.get('type','?')}:{x.get('numberOfIndividuals','?')}:{x.get('ancestralGroups',[{}])[0].get('ancestralGroup','?') if x.get('ancestralGroups') else '?'}"
                for x in ancestries
            )
            loci = a.get("loci", []) or []
            snp_repr = []
            for L in loci:
                for s in L.get("strongestRiskAlleles", []) or []:
                    snp_repr.append(s.get("riskAlleleName", ""))
            rows.append(
                {
                    "shortForm": short,
                    "label": label,
                    "study_accession": study.get("accessionId"),
                    "pubmed": study.get("publicationInfo", {}).get("pubmedId") if isinstance(study.get("publicationInfo"), dict) else None,
                    "title": study.get("publicationInfo", {}).get("title") if isinstance(study.get("publicationInfo"), dict) else None,
                    "first_author": (study.get("publicationInfo", {}) or {}).get("author", {}).get("fullname"),
                    "initialSampleSize": study.get("initialSampleSize"),
                    "replicationSampleSize": study.get("replicationSampleSize"),
                    "ancestries": anc_str,
                    "snp": "; ".join(snp_repr),
                    "pvalue": a.get("pvalue"),
                    "orPerCopyNum": a.get("orPerCopyNum"),
                    "betaNum": a.get("betaNum"),
                    "betaUnit": a.get("betaUnit"),
                    "betaDirection": a.get("betaDirection"),
                    "standardError": a.get("standardError"),
                    "riskFrequency": a.get("riskFrequency"),
                    "range": a.get("range"),
                }
            )
        if len(embed) < page_size:
            break
        page += 1
        time.sleep(0.2)
    if not rows:
        return None
    df = pd.DataFrame(rows)
    p = RAW / f"{short}_rest_associations.tsv"
    df.to_csv(p, sep="\t", index=False)
    print(f"[{short}/{label}] rest mode {len(df)} rows", flush=True)
    return df


def harmonise_download(df: pd.DataFrame) -> pd.DataFrame:
    """Convert the GWAS Catalog download TSV into a unified schema."""
    keep = pd.DataFrame(
        {
            "trait":            df["DISEASE/TRAIT"],
            "trait_short":      df["__source_id"],
            "trait_label":      df["__source_label"],
            "mapped_trait":     df.get("MAPPED_TRAIT"),
            "mapped_trait_uri": df.get("MAPPED_TRAIT_URI"),
            "study":            df["STUDY"],
            "study_accession":  df.get("STUDY ACCESSION"),
            "pubmed":           df["PUBMEDID"],
            "first_author":     df["FIRST AUTHOR"],
            "year":             df["DATE"],
            "journal":          df["JOURNAL"],
            "ancestry_initial": df["INITIAL SAMPLE SIZE"],
            "ancestry_replication": df["REPLICATION SAMPLE SIZE"],
            "snp":              df["SNPS"],
            "strongest_risk_allele": df["STRONGEST SNP-RISK ALLELE"],
            "chr":              df["CHR_ID"],
            "pos":              df["CHR_POS"],
            "region":           df["REGION"],
            "mapped_gene":      df["MAPPED_GENE"],
            "context":          df.get("CONTEXT"),
            "risk_allele_freq": df["RISK ALLELE FREQUENCY"],
            "pvalue":           df["P-VALUE"],
            "pvalue_mlog":      df["PVALUE_MLOG"],
            "or_or_beta":       df["OR or BETA"],
            "ci_text":          df["95% CI (TEXT)"],
            "platform":         df.get("PLATFORM [SNPS PASSING QC]"),
            "source":           "gwas_catalog_download",
        }
    )
    return keep


def harmonise_rest(df: pd.DataFrame) -> pd.DataFrame:
    keep = pd.DataFrame(
        {
            "trait":            df["label"],
            "trait_short":      df["shortForm"],
            "trait_label":      df["label"],
            "mapped_trait":     df["label"],
            "mapped_trait_uri": "http://www.ebi.ac.uk/efo/" + df["shortForm"],
            "study":            df["title"],
            "study_accession":  df["study_accession"],
            "pubmed":           df["pubmed"],
            "first_author":     df["first_author"],
            "year":             None,
            "journal":          None,
            "ancestry_initial": df["initialSampleSize"],
            "ancestry_replication": df["replicationSampleSize"],
            "snp":              df["snp"],
            "strongest_risk_allele": df["snp"],
            "chr":              None,
            "pos":              None,
            "region":           None,
            "mapped_gene":      None,
            "context":          None,
            "risk_allele_freq": df["riskFrequency"],
            "pvalue":           df["pvalue"],
            "pvalue_mlog":      None,
            "or_or_beta":       df["orPerCopyNum"].fillna(df["betaNum"]),
            "ci_text":          df["range"],
            "platform":         None,
            "source":           "gwas_catalog_rest",
        }
    )
    return keep


def main() -> None:
    parts: list[pd.DataFrame] = []
    for short, label, mode in TRAITS:
        if mode == "download":
            df = pull_download(short, label)
            if df is not None:
                parts.append(harmonise_download(df))
        else:
            df = pull_rest(short, label)
            if df is not None:
                parts.append(harmonise_rest(df))
        time.sleep(0.4)

    big = pd.concat(parts, ignore_index=True)
    # Coerce numeric pvalue
    big["pvalue_num"] = pd.to_numeric(big["pvalue"], errors="coerce")
    big["chr_num"] = pd.to_numeric(big["chr"], errors="coerce")
    big["pos_num"] = pd.to_numeric(big["pos"], errors="coerce")

    out = TAB / "T01_gwas_catalog_thyroid_autoimmunity_associations.tsv"
    big.to_csv(out, sep="\t", index=False)
    print(f"WROTE {out} rows={len(big)}")

    # Trait summary
    summary = (
        big.groupby(["trait_short", "trait_label"])
        .agg(
            n_rows=("trait", "size"),
            n_snp=("snp", lambda s: s.dropna().nunique()),
            n_studies=("study", lambda s: s.dropna().nunique()),
            min_p=("pvalue_num", "min"),
            n_genome_wide=("pvalue_num", lambda s: int((s < 5e-8).sum())),
        )
        .reset_index()
    )
    out = TAB / "T01b_per_trait_summary.tsv"
    summary.to_csv(out, sep="\t", index=False)
    print(f"WROTE {out}")


if __name__ == "__main__":
    main()
