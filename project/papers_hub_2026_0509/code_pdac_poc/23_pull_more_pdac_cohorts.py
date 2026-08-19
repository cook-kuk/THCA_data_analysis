"""
Aggressively pull more PUBLIC PDAC cohorts with bulk-RNA + OS for the
vaccine+ICI synergy meta-analysis.

Targets (immediate, public, no DACO):
  - GSE71729 (Moffitt 2015 Nat Genet, n=357 microarray + LCM, OS in clinical)
  - GSE57495 (Chen 2015 PLoS One, n=63 microarray + OS)
  - GSE21501 (Stratford 2010 PLoS Med, n=102 microarray + OS)
  - GSE62452 (Yang 2016 Cancer Res, n=130 microarray + OS)
  - GSE15471 (Pei 2009 Cancer Cell, n=39, normal-paired)

Strategy:
  Use NCBI GEO Series Matrix (.txt.gz) — single file per series,
  contains expression + sample-level clinical fields.
  Parse with GEOparse if installed, else regex on the raw matrix.

Outputs:
  raw/extra_geo/<gse>.tsv               canonicalized clinical
  raw/extra_geo/<gse>_matrix.tsv        expression matrix (subset to Moffitt + paradox genes)
  results/synergy/extra_cohorts.json    summary (n, OS coverage, allele coverage)
"""

import gzip
import json
import re
import sys
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import numpy as np
import requests

OUT = Path("/data/pdac_poc/raw/extra_geo")
OUT.mkdir(parents=True, exist_ok=True)
RES = Path("/data/pdac_poc/results/synergy")
RES.mkdir(parents=True, exist_ok=True)

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

# Genes we need (Moffitt + paradox modules)
MOFFITT_BASAL = ["VGLL1","UCA1","S100A2","LY6D","SPRR3","SPRR1B","LEMD1","LYPD3",
                 "KRT15","DHRS9","AREG","CST6","SERPINB3","KRT6A","SERPINB4",
                 "FAM83A","SCEL","FGFBP1","KRT7","KRT17","GPR87","TNS4","SLC2A1"]
MOFFITT_CLASS = ["BTNL8","FAM3D","AGR3","CTSE","LYZ","TFF2","TFF1","ANXA10",
                 "LGALS4","ECT2","CLRN3","MYO1A","CLDN18","LRRC31","TFF3","CDX2",
                 "SERPINA10","VSIG2","TSPAN8","ST6GALNAC1","AGR2","TOX3"]
TME_GENES = ["CXCL13","IFNG","HLA-DRA","HLA-DRB1","CD8A","CD68","CD163","CSF1R",
             "CD274","PDCD1","TIGIT","LAG3","HAVCR2","TOX",
             "ACTA2","FAP","POSTN","COL1A1","IL6","CXCL12","PDPN","S100A4"]
ALL_GENES = sorted(set(MOFFITT_BASAL + MOFFITT_CLASS + TME_GENES))


COHORTS = [
    ("GSE71729", "Moffitt 2015 Nat Genet"),
    ("GSE57495", "Chen 2015 PLoS One"),
    ("GSE21501", "Stratford 2010 PLoS Med"),
    ("GSE62452", "Yang 2016 Cancer Res"),
]


def gse_matrix_url(gse):
    """e.g. GSE71729 -> .../GSE71nnn/GSE71729/matrix/GSE71729_series_matrix.txt.gz"""
    prefix = re.sub(r'\d{1,3}$', 'nnn', gse)
    return f"{GEO_FTP}/{prefix}/{gse}/matrix/{gse}_series_matrix.txt.gz"


def download(url, dest):
    if Path(dest).exists() and Path(dest).stat().st_size > 1000:
        return True
    try:
        r = requests.get(url, timeout=180, stream=True)
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(8192):
                f.write(chunk)
        return Path(dest).stat().st_size > 1000
    except Exception as e:
        print(f"   download failed: {e}")
        return False


def parse_series_matrix(path):
    """Parse a GEO series_matrix.txt.gz file.
    Returns (clinical_df, expression_df).
    """
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    # Phenotype lines start with !Sample_
    pheno = {}
    sample_ids = None
    matrix_start = None
    matrix_end = None
    for i, line in enumerate(lines):
        if line.startswith("!Sample_"):
            key = line.split("\t",1)[0].lstrip("!").lower()
            vals = line.strip().split("\t")[1:]
            vals = [v.strip().strip('"') for v in vals]
            pheno.setdefault(key, []).append(vals)
        if line.startswith("!series_matrix_table_begin"):
            matrix_start = i + 1
        if line.startswith("!series_matrix_table_end"):
            matrix_end = i

    # Sample IDs
    if "sample_geo_accession" in pheno:
        sample_ids = pheno["sample_geo_accession"][0]

    # Build clinical DF: take all phenos that vary across samples
    clin_cols = {}
    for k, lol in pheno.items():
        # multiple lines for the same key (e.g. characteristics_ch1)
        for li, vals in enumerate(lol):
            col = f"{k}_{li}" if len(lol) > 1 else k
            if len(vals) == len(sample_ids or []):
                clin_cols[col] = vals
    clin = pd.DataFrame(clin_cols)
    if sample_ids:
        clin.insert(0, "sample_id", sample_ids)

    # Parse matrix portion
    if matrix_start and matrix_end:
        mat_lines = lines[matrix_start:matrix_end]
        if mat_lines:
            from io import StringIO
            mat_text = "".join(mat_lines)
            mat = pd.read_csv(StringIO(mat_text), sep="\t", index_col=0)
        else:
            mat = pd.DataFrame()
    else:
        mat = pd.DataFrame()
    return clin, mat


def extract_os_from_chars(clin):
    """Look for OS columns in characteristics_ch1_* fields."""
    out = {}
    for col in clin.columns:
        if col.startswith("characteristics") or "title" in col.lower():
            sample_vals = clin[col].astype(str).tolist()[:5]
            joined = " | ".join(sample_vals).lower()
            # OS time
            if "survival" in joined or ("os " in joined and ("month" in joined or "day" in joined)):
                out["os_col"] = col
            if "death" in joined or "vital" in joined or "status" in joined:
                out["event_col"] = col
            if "kras" in joined:
                out["kras_col"] = col
    # Return all phenotype columns for inspection
    return out


def fetch_one(args):
    gse, label = args
    print(f"[23] {gse} ({label}) …")
    url = gse_matrix_url(gse)
    dest = OUT / f"{gse}_series_matrix.txt.gz"
    if not download(url, dest):
        return {"gse": gse, "ok": False, "reason": "download failed", "url": url}
    print(f"     downloaded {dest.name} ({dest.stat().st_size // 1024} KB)")
    try:
        clin, mat = parse_series_matrix(dest)
    except Exception as e:
        return {"gse": gse, "ok": False, "reason": f"parse error: {e}"}
    # Save raw clin
    clin.to_csv(OUT / f"{gse}_clinical.tsv", sep="\t", index=False)
    # Subset matrix to genes we care about (try probe-id mapping by gene symbol)
    sub = pd.DataFrame()
    if not mat.empty:
        # The matrix index is probe IDs typically; we don't have probe→gene mapping
        # without GPL annotation download. Save the raw matrix dimensions.
        sub_genes_mask = mat.index.astype(str).str.upper().isin([g.upper() for g in ALL_GENES])
        sub = mat[sub_genes_mask]
    if not sub.empty:
        sub.to_csv(OUT / f"{gse}_expression_subset.tsv", sep="\t")

    # Inspect phenotype columns to find OS-like fields
    fields_found = extract_os_from_chars(clin)

    return {
        "gse": gse, "label": label, "ok": True,
        "n_samples": int(len(clin)),
        "matrix_genes_total": int(mat.shape[0]) if not mat.empty else 0,
        "matrix_samples": int(mat.shape[1]) if not mat.empty else 0,
        "panel_genes_recovered": int(len(sub)),
        "clinical_columns": list(clin.columns),
        "fields_found": fields_found,
    }


def main():
    t0 = time.time()
    results = []
    with ThreadPoolExecutor(max_workers=4) as ex:
        for r in ex.map(fetch_one, COHORTS):
            results.append(r)
            if r.get("ok"):
                print(f"     ✓ n={r['n_samples']} · matrix {r['matrix_genes_total']} probes × {r['matrix_samples']} samples · panel-genes recovered {r['panel_genes_recovered']}")
                if r["fields_found"]:
                    print(f"     phenotype fields: {r['fields_found']}")
            else:
                print(f"     ✗ {r.get('reason')}")

    summary = {"generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
               "elapsed_s": round(time.time()-t0, 1),
               "results": results,
               "n_cohorts_pulled": sum(1 for r in results if r.get("ok")),
               "total_samples_indexed": sum(r.get("n_samples",0) for r in results if r.get("ok"))}
    json.dump(summary, open(RES/"extra_cohorts.json","w"), indent=2, default=str)

    print("\n" + "="*70)
    print(f"[23] DONE · {summary['elapsed_s']}s")
    print(f"   pulled {summary['n_cohorts_pulled']}/{len(COHORTS)} cohorts")
    print(f"   total samples indexed: {summary['total_samples_indexed']}")
    print(f"   files in {OUT}/")
    for r in results:
        if r.get("ok"):
            print(f"     • {r['gse']:10s} ({r['label']:30s}) n={r['n_samples']}  panel-genes-found={r['panel_genes_recovered']}")


if __name__ == "__main__":
    main()
