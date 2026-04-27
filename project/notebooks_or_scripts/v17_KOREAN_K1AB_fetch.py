#!/usr/bin/env python3
"""
v17 KOREAN K1-A + K1-B — PRJEB11591 ENA metadata + RNASeq-er count matrix.

K1-A: ENA portal API for sample metadata
K1-B: ENA RNASeq-er API for pre-quantified gene-level counts (skips alignment)

Fallback: if RNASeq-er has no data for PRJEB11591, log and skip alignment
(STAR alignment of 120 PE samples = 5-10 hours, deferred to manual run).

Outputs:
  - results/v17_korean/K1A_prjeb11591_runs.tsv
  - results/v17_korean/K1A_summary.json
  - results/v17_korean/K1B_count_matrix.tsv (if available)
  - results/v17_korean/K1B_summary.json
"""
from __future__ import annotations
import io
import json
import sys
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.parse import urlencode

import pandas as pd

OUT = Path("/opt/thyroid-dash/project/results/v17_korean")
OUT.mkdir(parents=True, exist_ok=True)

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
RNASEQER_API = "https://www.ebi.ac.uk/fg/rnaseq/api/json/getRunsByStudy"
RNASEQER_QUANT = "https://www.ebi.ac.uk/fg/rnaseq/api/json/getStudyDescription/PRJEB11591"


def fetch(url, timeout=120):
    req = Request(url, headers={"Accept": "application/json,text/plain", "User-Agent": "v17-korean/1.0"})
    with urlopen(req, timeout=timeout) as r:
        return r.read()


def k1a_metadata():
    print("\n[K1-A] ENA metadata for PRJEB11591", flush=True)
    params = {
        "accession": "PRJEB11591",
        "result": "read_run",
        "fields": "run_accession,sample_alias,sample_title,read_count,base_count,fastq_ftp,fastq_md5,library_strategy,library_layout,instrument_model,study_accession,experiment_accession,scientific_name,tax_id",
        "format": "tsv",
    }
    url = f"{ENA_PORTAL}?{urlencode(params)}"
    print(f"  GET {url}", flush=True)
    try:
        text = fetch(url).decode("utf-8")
        df = pd.read_csv(io.StringIO(text), sep="\t")
        print(f"  ✓ {len(df)} runs returned", flush=True)
    except Exception as e:
        print(f"  ✗ ENA fetch failed: {e}", flush=True)
        return {"error": str(e)}

    df.to_csv(OUT / "K1A_prjeb11591_runs.tsv", sep="\t", index=False)

    # heuristic histology mapping from sample_title (may be opaque alias)
    df["histology_guess"] = df.get("sample_title", "").astype(str).str.lower().apply(
        lambda s: ("ATC" if "atc" in s or "anaplastic" in s
                   else "PDTC" if "pdtc" in s or "poorly" in s
                   else "FVPTC" if "fvptc" in s or "follicular variant" in s
                   else "FTC" if "ftc" in s and "fvptc" not in s
                   else "FA" if " fa " in s or "adenoma" in s
                   else "Normal" if "normal" in s or "ctrl" in s
                   else "PTC" if "ptc" in s or "papillary" in s or "thyroid" in s
                   else "?")
    )
    histo = df["histology_guess"].value_counts().to_dict()

    summary = {
        "study": "PRJEB11591",
        "n_runs": int(len(df)),
        "n_unique_samples": int(df["sample_alias"].nunique()) if "sample_alias" in df.columns else int(len(df)),
        "library_strategy": df["library_strategy"].value_counts().to_dict() if "library_strategy" in df.columns else {},
        "library_layout": df["library_layout"].value_counts().to_dict() if "library_layout" in df.columns else {},
        "instrument": df["instrument_model"].value_counts().to_dict() if "instrument_model" in df.columns else {},
        "histology_guess_counts": histo,
        "median_read_count": int(df["read_count"].median()) if "read_count" in df.columns else None,
        "fastq_ftp_present": int(df["fastq_ftp"].astype(str).str.contains("ebi.ac.uk").sum()) if "fastq_ftp" in df.columns else 0,
    }
    with open(OUT / "K1A_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"  histology guess: {histo}", flush=True)
    return summary


def k1b_rnaseqer():
    print("\n[K1-B] RNASeq-er pre-quantified counts for PRJEB11591", flush=True)
    url = f"{RNASEQER_API}/PRJEB11591"
    print(f"  GET {url}", flush=True)
    try:
        data = fetch(url, timeout=60)
        runs = json.loads(data.decode("utf-8"))
        if isinstance(runs, dict) and runs.get("ERROR"):
            return {"error": "RNASeq-er study not indexed", "msg": runs.get("ERROR")}
        if not runs:
            return {"error": "RNASeq-er returned empty"}
        print(f"  ✓ {len(runs)} runs found in RNASeq-er", flush=True)
    except Exception as e:
        return {"error": f"RNASeq-er API failed: {e}"}

    # Try to fetch quantification per run
    # RNASeq-er provides ftp_files including FPKMs / counts
    pd.DataFrame(runs).to_csv(OUT / "K1B_rnaseqer_runs.tsv", sep="\t", index=False)
    sample = runs[0] if runs else {}
    summary = {
        "n_runs_indexed": len(runs),
        "first_run_keys": list(sample.keys())[:30],
        "first_run_sample": {k: v for k, v in list(sample.items())[:10]} if sample else {},
    }
    # Look for FTP URL fields
    for key in sample:
        if "ftp" in key.lower() or "url" in key.lower() or "file" in key.lower():
            summary[f"_field_{key}"] = sample[key]

    with open(OUT / "K1B_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    return summary


def main():
    a = k1a_metadata()
    b = k1b_rnaseqer()
    print("\n=== K1 SUMMARY ===")
    print("K1-A:", json.dumps({k: v for k, v in a.items() if k != "histology_guess_counts"}, indent=2, default=str))
    print("K1-B:", json.dumps(b, indent=2, default=str)[:1500])
    return 0


if __name__ == "__main__":
    sys.exit(main())
