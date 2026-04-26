#!/usr/bin/env python3
"""v5.1 Phase 1 — REAL TCGA + GEO downloads.

Sequential per cancer. Concurrent file fetch (max_workers=4) within cancer.
No synthetic fallback. Logs errors + emits v5p1_cohort_availability.tsv.
"""
from __future__ import annotations

import io
import os
import re
import sys
import json
import gzip
import time
import hashlib
import urllib.request
import urllib.error
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v5p1_common import (
    PROJECT, DATA_RAW_V5, RESULTS_V5, LOGS, GDC_CACHE, COHORTS, CANCERS, log_line
)

LOG = LOGS / "v5p1_download.log"
ERR_LOG = RESULTS_V5 / "v5p1_download_errors.log"
GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data/"

# Budgets to keep runtime bounded
MAX_RNASEQ_FILES = int(os.environ.get("V5P1_MAX_RNASEQ", "260"))  # per cancer
# MAFs: one file per (patient, aliquot) — TCGA ~1 maf/patient for primary
# cohorts. Need to grab all to get full patient coverage.
MAX_MAF_FILES = int(os.environ.get("V5P1_MAX_MAF", "600"))
DOWNLOAD_TIMEOUT = 180
CONCURRENT_DL = 6


def _ua_request(url: str, data: Optional[bytes] = None, headers: Optional[Dict] = None, timeout: int = 60):
    req = urllib.request.Request(url, data=data, headers=headers or {})
    req.add_header("User-Agent", "v5p1-dial-downloader/1.0")
    return urllib.request.urlopen(req, timeout=timeout)


# ============================================================
# TCGA: GDC API query
# ============================================================
def gdc_query_files(project_id: str, data_type: str, workflow_type: str, size: int = 5000) -> List[Dict]:
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project_id]}},
            {"op": "in", "content": {"field": "data_type", "value": [data_type]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": [workflow_type]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
        ],
    }
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type",
        "format": "JSON",
        "size": str(size),
    }
    body = urllib.parse.urlencode(params).encode("utf-8")
    try:
        with _ua_request(GDC_FILES, data=body, timeout=120) as r:
            resp = json.loads(r.read().decode("utf-8"))
        hits = resp.get("data", {}).get("hits", [])
        return hits
    except Exception as e:
        log_line(LOG, f"[gdc-query-err] {project_id} {data_type}: {e}")
        with open(ERR_LOG, "a") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\tGDC_QUERY\t{project_id}\t{data_type}\t{e}\n")
        return []


def _verify_cached(dest: Path) -> bool:
    """Return True iff cached file passes integrity. For .gz: full gunzip scan.
    Otherwise: size > 0. On failure, unlink so caller re-downloads."""
    if not (dest.exists() and dest.stat().st_size > 0):
        return False
    if dest.name.endswith(".gz"):
        try:
            with gzip.open(dest, "rb") as gf:
                while gf.read(1024 * 1024):
                    pass
            return True
        except Exception as e:
            with open(ERR_LOG, "a") as fh:
                fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\tINTEGRITY\t{dest.name}\tpartial gz, removing: {e}\n")
            try:
                dest.unlink()
            except Exception:
                pass
            return False
    return True


def _download_file(file_id: str, dest: Path, retries: int = 2) -> bool:
    if _verify_cached(dest):
        return True
    last_err = None
    for _ in range(retries + 1):
        try:
            url = GDC_DATA + file_id
            with _ua_request(url, timeout=DOWNLOAD_TIMEOUT) as r:
                data = r.read()
            dest.parent.mkdir(parents=True, exist_ok=True)
            tmp = dest.with_suffix(dest.suffix + ".part")
            with open(tmp, "wb") as fh:
                fh.write(data)
            tmp.replace(dest)
            return True
        except Exception as e:
            last_err = e
            time.sleep(1.0)
    with open(ERR_LOG, "a") as fh:
        fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\tGDC_DOWNLOAD\t{file_id}\t{last_err}\n")
    return False


def _submitter_to_sample_id(case: Dict) -> str:
    """Extract 16-char TCGA sample barcode."""
    samples = case.get("samples", [])
    if samples:
        sid = samples[0].get("submitter_id", "")
        if sid:
            return sid[:16]
    # Fallback: case submitter_id + sample_type code
    case_sid = case.get("submitter_id", "")
    if case_sid and samples:
        st = samples[0].get("sample_type", "")
        code = {"Primary Tumor": "01A", "Metastatic": "06A",
                "Recurrent Tumor": "02A", "Solid Tissue Normal": "11A"}.get(st, "01A")
        return f"{case_sid}-{code}"[:16]
    return case_sid[:16] if case_sid else ""


def download_tcga_rnaseq(cancer: str, max_files: int = MAX_RNASEQ_FILES) -> Dict:
    """Download TCGA STAR counts and build gene x sample matrix."""
    project_id = COHORTS[cancer]["tcga"]
    out_dir = DATA_RAW_V5 / cancer / "tcga_rnaseq"
    out_dir.mkdir(parents=True, exist_ok=True)

    hits = gdc_query_files(project_id, "Gene Expression Quantification", "STAR - Counts")
    log_line(LOG, f"[{cancer}] RNA-seq hits: {len(hits)}")

    # Subsample: keep primary tumor samples first
    def is_primary(h):
        cases = h.get("cases", [])
        if not cases:
            return False
        for s in cases[0].get("samples", []):
            if "Primary" in s.get("sample_type", ""):
                return True
        return False

    primary = [h for h in hits if is_primary(h)]
    other = [h for h in hits if not is_primary(h)]
    chosen = (primary + other)[:max_files]

    file_to_sample = {}
    for h in chosen:
        fid = h["file_id"]
        cases = h.get("cases", [])
        if not cases:
            continue
        case = cases[0]
        sid = _submitter_to_sample_id(case)
        if not sid:
            continue
        file_to_sample[fid] = sid

    log_line(LOG, f"[{cancer}] downloading {len(file_to_sample)} STAR-count files with {CONCURRENT_DL} workers")

    dest_paths = {fid: out_dir / f"{fid}.tsv" for fid in file_to_sample}
    ok = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=CONCURRENT_DL) as ex:
        futs = {ex.submit(_download_file, fid, dest_paths[fid]): fid for fid in file_to_sample}
        for fut in as_completed(futs):
            if fut.result():
                ok += 1
            if ok > 0 and ok % 30 == 0:
                log_line(LOG, f"[{cancer}] RNA-seq downloaded {ok}/{len(file_to_sample)} (elapsed {time.time()-t0:.0f}s)")

    log_line(LOG, f"[{cancer}] RNA-seq DONE {ok}/{len(file_to_sample)} files in {time.time()-t0:.0f}s")

    # Parse files into gene x sample matrix using stranded_second
    counts_matrix = {}
    gene_names_ref = None
    parsed_samples = []
    for fid, sid in file_to_sample.items():
        f = dest_paths[fid]
        if not f.exists() or f.stat().st_size == 0:
            continue
        try:
            df = pd.read_csv(f, sep="\t", comment="#", skiprows=1)
            df = df[~df["gene_id"].str.startswith("N_", na=False)]
            if gene_names_ref is None:
                gene_names_ref = df["gene_name"].values
                gene_ids_ref = df["gene_id"].values
            col_name = sid
            suffix = 1
            while col_name in counts_matrix:
                col_name = f"{sid}_d{suffix}"
                suffix += 1
            counts_matrix[col_name] = df["stranded_second"].values
            parsed_samples.append(col_name)
        except Exception as e:
            with open(ERR_LOG, "a") as fh:
                fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\tGDC_PARSE\t{fid}\t{e}\n")

    if not counts_matrix:
        return {"status": "failed", "n_samples": 0}

    mat = pd.DataFrame(counts_matrix, index=gene_names_ref)
    # collapse duplicate gene_name by sum
    mat = mat.groupby(level=0).sum()
    # atomic save: write to .part then rename, so SIGTERM can't leave partial gz
    out_path = DATA_RAW_V5 / cancer / "tcga_counts.tsv.gz"
    tmp_out = out_path.with_suffix(out_path.suffix + ".part")
    mat.to_csv(tmp_out, sep="\t", compression="gzip")
    tmp_out.replace(out_path)
    pd.DataFrame({"sample_id": parsed_samples, "patient_id": [s[:12] for s in parsed_samples]}).to_csv(
        DATA_RAW_V5 / cancer / "tcga_samples.tsv", sep="\t", index=False)

    log_line(LOG, f"[{cancer}] RNA-seq matrix: {mat.shape[0]} genes x {mat.shape[1]} samples saved")
    return {"status": "ok", "n_samples": mat.shape[1], "n_genes": mat.shape[0]}


def download_tcga_maf(cancer: str, max_files: int = MAX_MAF_FILES) -> Dict:
    project_id = COHORTS[cancer]["tcga"]
    out_dir = DATA_RAW_V5 / cancer / "tcga_maf"
    out_dir.mkdir(parents=True, exist_ok=True)

    hits = gdc_query_files(
        project_id, "Masked Somatic Mutation",
        "Aliquot Ensemble Somatic Variant Merging and Masking",
    )
    log_line(LOG, f"[{cancer}] MAF hits: {len(hits)}")
    chosen = hits[:max_files]

    dest_paths = {}
    for h in chosen:
        fid = h["file_id"]
        dest_paths[fid] = out_dir / f"{fid}.maf.gz"

    ok = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=CONCURRENT_DL) as ex:
        futs = {ex.submit(_download_file, fid, dest_paths[fid]): fid for fid in dest_paths}
        for fut in as_completed(futs):
            if fut.result():
                ok += 1

    log_line(LOG, f"[{cancer}] MAF DONE {ok}/{len(dest_paths)} in {time.time()-t0:.0f}s")

    # Parse all MAFs, concatenate
    maf_rows = []
    wanted_cols = ["Hugo_Symbol", "Tumor_Sample_Barcode", "HGVSp_Short", "Variant_Classification"]
    for fid, dest in dest_paths.items():
        if not dest.exists() or dest.stat().st_size == 0:
            continue
        try:
            df = pd.read_csv(dest, sep="\t", comment="#", compression="gzip", low_memory=False)
            keep = [c for c in wanted_cols if c in df.columns]
            if len(keep) < 3:
                continue
            maf_rows.append(df[keep])
        except Exception as e:
            with open(ERR_LOG, "a") as fh:
                fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\tMAF_PARSE\t{fid}\t{e}\n")

    if not maf_rows:
        return {"status": "failed"}
    maf = pd.concat(maf_rows, ignore_index=True)
    maf["patient_id"] = maf["Tumor_Sample_Barcode"].astype(str).str[:12]
    maf_out = DATA_RAW_V5 / cancer / "tcga_maf_combined.tsv.gz"
    tmp_maf = maf_out.with_suffix(maf_out.suffix + ".part")
    maf.to_csv(tmp_maf, sep="\t", index=False, compression="gzip")
    tmp_maf.replace(maf_out)
    log_line(LOG, f"[{cancer}] MAF rows: {len(maf)} patients: {maf['patient_id'].nunique()}")
    return {"status": "ok", "n_rows": len(maf), "n_patients": int(maf["patient_id"].nunique())}


# ============================================================
# Driver extraction per cancer
# ============================================================
def extract_driver_labels(cancer: str) -> Optional[pd.DataFrame]:
    """Return DataFrame with columns [patient_id, label] for class_A vs class_B."""
    maf_path = DATA_RAW_V5 / cancer / "tcga_maf_combined.tsv.gz"
    if not maf_path.exists():
        return None
    maf = pd.read_csv(maf_path, sep="\t", compression="gzip", low_memory=False)

    def has_mut(df_gene, hgvsp_regexes=None, classes=None):
        if df_gene.empty:
            return set()
        if hgvsp_regexes:
            hits = df_gene["HGVSp_Short"].astype(str).str.contains(
                "|".join(hgvsp_regexes), regex=True, na=False)
            return set(df_gene.loc[hits, "patient_id"].tolist())
        if classes:
            hits = df_gene["Variant_Classification"].isin(classes)
            return set(df_gene.loc[hits, "patient_id"].tolist())
        return set(df_gene["patient_id"].tolist())

    all_patients = set(maf["patient_id"].unique())

    # Note: current GDC MAF uses newer BRAF transcript (NM_001374258) with
    # offset numbering: canonical V600E = V640E, V600K = V640K. Match both.
    if cancer == "SKCM":
        braf = has_mut(maf[maf["Hugo_Symbol"] == "BRAF"],
                       ["V600E", "V600K", "V640E", "V640K", "V640M"])
        nras = has_mut(maf[maf["Hugo_Symbol"] == "NRAS"], ["Q61R", "Q61K", "Q61L"])
        class_A, class_B = braf, nras
        name_A, name_B = "BRAF", "NRAS"
    elif cancer == "LGG":
        # IDH1 canonical R132 = R172 in newer transcript; check both
        idh1 = has_mut(maf[maf["Hugo_Symbol"] == "IDH1"], ["R132", "R172"])
        idh2 = has_mut(maf[maf["Hugo_Symbol"] == "IDH2"], ["R172", "R212"])
        idhmut = idh1 | idh2
        idhwt = all_patients - idhmut
        class_A, class_B = idhmut, idhwt
        name_A, name_B = "IDHmut", "IDHwt"
    elif cancer == "LUAD":
        kras = has_mut(maf[maf["Hugo_Symbol"] == "KRAS"], ["G12", "G13"])
        egfr_l858 = has_mut(maf[maf["Hugo_Symbol"] == "EGFR"], ["L858R", "L898R"])
        egfr_exon19 = has_mut(
            maf[(maf["Hugo_Symbol"] == "EGFR") & (maf["Variant_Classification"] == "In_Frame_Del")]
        )
        egfr = egfr_l858 | egfr_exon19
        class_A, class_B = kras, egfr
        name_A, name_B = "KRAS", "EGFR"
    elif cancer == "COAD":
        braf = has_mut(maf[maf["Hugo_Symbol"] == "BRAF"],
                       ["V600", "V640"])
        kras = has_mut(maf[maf["Hugo_Symbol"] == "KRAS"], ["G12", "G13"])
        class_A, class_B = braf, kras
        name_A, name_B = "BRAF", "KRAS"
    else:
        return None

    # drop ambiguous (both) or neither
    only_A = class_A - class_B
    only_B = class_B - class_A
    rows = []
    for p in only_A:
        rows.append((p, name_A))
    for p in only_B:
        rows.append((p, name_B))
    df = pd.DataFrame(rows, columns=["patient_id", "label"])
    return df


# ============================================================
# GEO download via series matrix
# ============================================================
def _gse_ftp_url(gse: str) -> str:
    stub = gse[:-3] + "nnn" if len(gse) > 3 else gse
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/{gse}_series_matrix.txt.gz"


LABEL_PATTERNS = [
    re.compile(r'braf[- _]?(mutation|status|mut)', re.I),
    re.compile(r'nras', re.I),
    re.compile(r'\bidh[12]?\b', re.I),
    re.compile(r'\bkras\b', re.I),
    re.compile(r'\begfr\b', re.I),
    re.compile(r'subtype', re.I),
    re.compile(r'cms[1-4]', re.I),
    re.compile(r'\bmsi\b', re.I),
    re.compile(r'\bmss\b', re.I),
    re.compile(r'mutation status', re.I),
    re.compile(r'genotype', re.I),
    re.compile(r'wild[- ]?type', re.I),
    re.compile(r'\bv600\b', re.I),
    re.compile(r'driver', re.I),
]


def parse_series_matrix_header(raw_text: str) -> Dict:
    """Extract !Sample_* header attributes from series matrix."""
    attrs = {}
    sample_ids = []
    for line in raw_text.split("\n"):
        if not line.startswith("!"):
            continue
        if line.startswith("!series_matrix_table_begin"):
            break
        parts = line.split("\t")
        key = parts[0].lstrip("!").strip()
        vals = [p.strip().strip('"') for p in parts[1:]]
        if key == "Sample_geo_accession":
            sample_ids = vals
        if key.startswith("Sample_"):
            attrs.setdefault(key, []).append(vals)
    return {"attrs": attrs, "sample_ids": sample_ids}


def extract_expression_matrix(raw_text: str) -> Optional[pd.DataFrame]:
    lines = raw_text.split("\n")
    start = None
    end = None
    for i, line in enumerate(lines):
        if line.startswith("!series_matrix_table_begin"):
            start = i + 1
        elif line.startswith("!series_matrix_table_end"):
            end = i
            break
    if start is None or end is None:
        return None
    table = "\n".join(lines[start:end])
    try:
        df = pd.read_csv(io.StringIO(table), sep="\t", index_col=0)
        return df
    except Exception:
        return None


def assign_label_from_attrs(attrs: Dict, sample_ids: List[str], cancer: str) -> Optional[pd.Series]:
    """Return a Series index=sample_id, value=class_A/class_B/None.

    Builds per-sample strings by joining all Sample_characteristics_ch1 rows
    element-wise, then applies cancer-specific regex.
    """
    if not sample_ids:
        return None
    n = len(sample_ids)

    rows = []
    for key, kr in attrs.items():
        kl = key.lower()
        if ("characteristic" in kl) or ("source_name" in kl) or ("title" in kl) or ("description" in kl):
            for r in kr:
                if len(r) == n:
                    rows.append(r)

    if not rows:
        return None

    combined = [""] * n
    for r in rows:
        for i, v in enumerate(r):
            combined[i] = (combined[i] + " || " + str(v)).strip(" |")

    class_a = COHORTS[cancer]["class_a"]
    class_b = COHORTS[cancer]["class_b"]

    labels = []
    for text in combined:
        t = text.lower()
        a = False
        b = False
        if cancer == "THCA":
            # GSE27155 / GSE33630 / GSE29265 "braf t1799a mutation: Y/N"
            # and "kras, nras, or hras mutation: Y/N"
            if re.search(r'braf[^:]{0,40}mutation[: ]*y\b', t):
                a = True
            if re.search(r'(kras|nras|hras)[^:]{0,40}mutation[: ]*y\b', t):
                b = True
            # also explicit "braf v600e" as positive
            if re.search(r'braf[- ]?v?600e', t):
                a = True
        elif cancer == "SKCM":
            # GSE22153 "braf/nras: BRAF-V600E" or "NRAS-Q61K"
            if re.search(r'braf[- ]?(v600|mut|pos|mutant)', t):
                a = True
            if re.search(r'nras[- ]?(mut|q61|pos|mutant)', t):
                b = True
        elif cancer == "LGG":
            # Look for IDH status. Values typically "IDH1 mutation: mut" etc.
            if re.search(r'idh[12]?[^a-z0-9]{0,8}(mut|r132|r172|pos|positive|mutant|present)', t):
                a = True
            if re.search(r'idh[12]?[^a-z0-9]{0,8}(wt|wild[- ]?type|neg|negative|not[- ]?detect|absent)', t):
                b = True
        elif cancer == "LUAD":
            # GSE72094: "kras_status: Mut/WT" and "egfr_status: Mut/WT"
            m_kras = re.search(r'kras[_ ]status[: ]*(mut|wt|wild[- ]?type)', t)
            m_egfr = re.search(r'egfr[_ ]status[: ]*(mut|wt|wild[- ]?type)', t)
            # GSE31210: "gene alteration status: KRAS mutation +" or "EGFR mutation +"
            m_alt = re.search(r'(kras|egfr|alk)[ -]*(mutation|fusion)\s*\+', t)
            if m_kras:
                if m_kras.group(1) == "mut":
                    a = True
                # don't mark WT as B unless no EGFR
            if m_egfr:
                if m_egfr.group(1) == "mut":
                    b = True
            if m_alt:
                g = m_alt.group(1)
                if g == "kras":
                    a = True
                if g == "egfr":
                    b = True
            # also direct amino-acid substring
            if re.search(r'kras[- ]?(g12|g13|mutant)', t):
                a = True
            if re.search(r'egfr[- ]?(l858|exon[ -]?19|mutant)', t):
                b = True
        elif cancer == "COAD":
            # GSE39582: "braf.mutation: M/WT" and "kras.mutation: M/WT"
            # class_A=BRAF, class_B=KRAS
            m_braf = re.search(r'braf[.\s_]*mutation[: ]*([a-z/]+)', t)
            m_kras = re.search(r'kras[.\s_]*mutation[: ]*([a-z/]+)', t)
            if m_braf and m_braf.group(1).startswith("m"):
                a = True
            if m_kras and m_kras.group(1).startswith("m"):
                b = True
            # also support cms/msi tags
            if re.search(r'cms1|msi[- ]?h|msi-high', t):
                a = True
            if re.search(r'cms2|mss\b', t):
                b = True
        if a and not b:
            labels.append(class_a)
        elif b and not a:
            labels.append(class_b)
        else:
            labels.append(None)
    s = pd.Series(labels, index=sample_ids, name="label")
    return s


def download_geo_series(cancer: str, gse: str) -> Dict:
    out_dir = DATA_RAW_V5 / cancer / "geo"
    out_dir.mkdir(parents=True, exist_ok=True)
    mat_path = out_dir / f"{gse}_series_matrix.txt.gz"
    # download
    if not mat_path.exists() or mat_path.stat().st_size == 0:
        url = _gse_ftp_url(gse)
        try:
            with _ua_request(url, timeout=DOWNLOAD_TIMEOUT) as r:
                data = r.read()
            with open(mat_path, "wb") as fh:
                fh.write(data)
        except Exception as e:
            with open(ERR_LOG, "a") as fh:
                fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\tGEO_DL\t{gse}\t{e}\n")
            return {"status": "excluded_download_fail", "error": str(e)}

    try:
        with gzip.open(mat_path, "rt") as fh:
            raw = fh.read()
    except Exception as e:
        return {"status": "excluded_download_fail", "error": str(e)}

    parsed = parse_series_matrix_header(raw)
    labels = assign_label_from_attrs(parsed["attrs"], parsed["sample_ids"], cancer)
    if labels is None or labels.notna().sum() < 10:
        return {"status": "excluded_no_labels", "n_samples": len(parsed["sample_ids"]),
                "n_with_label": int(labels.notna().sum()) if labels is not None else 0}

    expr = extract_expression_matrix(raw)
    if expr is None or expr.shape[1] < 20:
        return {"status": "excluded_no_labels", "note": "no expression matrix"}

    # Map probes -> gene symbols using platform
    platform = None
    for row in parsed["attrs"].get("Sample_platform_id", []):
        if row:
            platform = row[0]
            break
    gene_expr = _probes_to_symbols(expr, platform)
    if gene_expr is None or gene_expr.shape[0] < 2000:
        # try keep original symbols if expr index looks like gene symbols
        if expr.index.astype(str).str.match(r'^[A-Z][A-Z0-9]{1,10}$').mean() > 0.5:
            gene_expr = expr
        else:
            return {"status": "excluded_no_labels", "note": f"probe->symbol mapping failed ({platform})"}

    # save
    expr_path = out_dir / f"{gse}_expr.tsv"
    pheno_path = out_dir / f"{gse}_pheno.tsv"
    gene_expr.to_csv(expr_path, sep="\t")
    pheno = pd.DataFrame({"sample_id": parsed["sample_ids"], "label": labels.values})
    pheno.to_csv(pheno_path, sep="\t", index=False)

    n_a = int((pheno["label"] == COHORTS[cancer]["class_a"]).sum())
    n_b = int((pheno["label"] == COHORTS[cancer]["class_b"]).sum())
    return {
        "status": "ok",
        "n_samples": len(parsed["sample_ids"]),
        "n_class_A": n_a,
        "n_class_B": n_b,
        "n_genes": gene_expr.shape[0],
        "platform": platform,
    }


def _probes_to_symbols(expr: pd.DataFrame, platform: Optional[str]) -> Optional[pd.DataFrame]:
    """Resolve probe IDs to gene symbols.

    Uses GEOparse-cached annotations if available, else local platform cache,
    else NCBI GEO platform annotation FTP. Aggregates probes -> symbol by mean.
    """
    if platform is None:
        return None
    plat_dir = PROJECT / "data_raw" / "geo" / "platforms"
    plat_dir.mkdir(parents=True, exist_ok=True)
    annot_path = plat_dir / f"{platform}.annot.csv"
    mapping = None
    if annot_path.exists():
        try:
            mapping = pd.read_csv(annot_path)
        except Exception:
            mapping = None
    if mapping is None:
        # Try to download .annot.gz from NCBI
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{platform[:-3]}nnn/{platform}/annot/{platform}.annot.gz"
        try:
            with _ua_request(url, timeout=DOWNLOAD_TIMEOUT) as r:
                data = r.read()
            with gzip.open(io.BytesIO(data), "rt") as fh:
                # .annot format has header lines starting with "!" then "ID\tGene symbol\t..."
                lines = []
                start = False
                for line in fh:
                    if not start:
                        if line.startswith("!platform_table_begin"):
                            start = True
                        continue
                    if line.startswith("!platform_table_end"):
                        break
                    lines.append(line)
            if lines:
                df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", low_memory=False)
                # normalize column names
                sym_col = None
                for c in df.columns:
                    cl = c.lower()
                    if cl in ("gene symbol", "gene_symbol", "symbol", "gene"):
                        sym_col = c
                        break
                if sym_col is None:
                    for c in df.columns:
                        if "symbol" in c.lower():
                            sym_col = c
                            break
                if sym_col is None:
                    return None
                mapping = pd.DataFrame({"probe": df["ID"].astype(str), "symbol": df[sym_col].astype(str)})
                mapping.to_csv(annot_path, index=False)
        except Exception as e:
            # try platform data page (GPL*.txt) which some platforms only have
            alt = plat_dir / f"{platform}.txt"
            if alt.exists():
                try:
                    with open(alt, "r", errors="ignore") as fh:
                        txt = fh.read()
                    # look for header row with ID and Gene symbol
                    lines = txt.split("\n")
                    hdr_idx = None
                    for i, line in enumerate(lines):
                        if line.startswith("ID\t") or line.startswith('"ID"\t'):
                            hdr_idx = i
                            break
                    if hdr_idx is not None:
                        tdf = pd.read_csv(io.StringIO("\n".join(lines[hdr_idx:])),
                                          sep="\t", low_memory=False)
                        sym_col = None
                        for c in tdf.columns:
                            if "symbol" in c.lower() or c.lower() == "gene_symbol":
                                sym_col = c
                                break
                        if sym_col is None:
                            for c in tdf.columns:
                                if "gene" in c.lower():
                                    sym_col = c
                                    break
                        if sym_col is not None:
                            mapping = pd.DataFrame({"probe": tdf["ID"].astype(str),
                                                    "symbol": tdf[sym_col].astype(str)})
                            mapping.to_csv(annot_path, index=False)
                except Exception as e2:
                    with open(ERR_LOG, "a") as fh:
                        fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\tGEO_PLATFORM\t{platform}\t{e2}\n")
                    return None
            else:
                with open(ERR_LOG, "a") as fh:
                    fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\tGEO_PLATFORM\t{platform}\t{e}\n")
                return None
    if mapping is None:
        return None
    mapping = mapping.dropna(subset=["symbol"])
    mapping = mapping[mapping["symbol"].astype(str).str.strip() != ""]
    mapping = mapping[~mapping["symbol"].astype(str).str.contains("///")]  # ambiguous
    mapping["probe"] = mapping["probe"].astype(str)

    expr = expr.copy()
    expr.index = expr.index.astype(str)
    m = mapping.drop_duplicates("probe").set_index("probe")["symbol"]
    expr = expr[expr.index.isin(m.index)]
    expr["__sym"] = expr.index.map(m)
    expr = expr.dropna(subset=["__sym"])
    # collapse to symbol by mean
    expr = expr.groupby("__sym").mean()
    expr = expr[~expr.index.astype(str).str.contains("^nan$", case=False, regex=True, na=False)]
    return expr


# ============================================================
# THCA shortcut: reuse existing real processed data
# ============================================================
def thca_shortcut() -> Dict:
    from v5p1_common import PROJECT
    expr_path = PROJECT / "data_processed" / "bulk_rnaseq" / "TCGA-THCA_rnaseq_expression_log2.tsv"
    master_path = PROJECT / "metadata" / "sample_master_v3_merged.tsv"
    if not expr_path.exists():
        return {"status": "excluded_download_fail", "error": f"{expr_path} missing"}
    m = pd.read_csv(master_path, sep="\t", dtype=str)
    m = m.set_index("sample_id", drop=False)
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    samples = expr.columns.tolist()
    labels = []
    for sid in samples:
        if sid not in m.index:
            labels.append(None); continue
        anchor = str(m.loc[sid].get("driver_anchor", "")).upper()
        if "BRAF" in anchor:
            labels.append("BRAF")
        elif "RAS" in anchor:
            labels.append("RAS")
        else:
            labels.append(None)
    s = pd.Series(labels, index=samples, name="label")
    keep = s.notna()
    n_a = int((s == "BRAF").sum())
    n_b = int((s == "RAS").sum())
    # save in v5p1 format
    out_dir = DATA_RAW_V5 / "THCA"
    out_dir.mkdir(parents=True, exist_ok=True)
    # THCA is log2 TPM already — use directly
    expr.loc[:, keep.index[keep].tolist()].to_csv(out_dir / "tcga_log2tpm.tsv.gz",
                                                  sep="\t", compression="gzip")
    pd.DataFrame({"sample_id": keep.index[keep].tolist(),
                  "patient_id": [x[:12] for x in keep.index[keep].tolist()],
                  "label": s.loc[keep.index[keep].tolist()].values
                  }).to_csv(out_dir / "tcga_samples_labels.tsv", sep="\t", index=False)
    return {"status": "ok", "n_samples": int(keep.sum()), "n_class_A": n_a, "n_class_B": n_b,
            "n_genes": int(expr.shape[0])}


# ============================================================
# Main
# ============================================================
def _tcga_usable_intersection(cancer: str, name_A: str, name_B: str,
                              expr_name: str, lab_name: str) -> Dict:
    """Recompute total_downloaded / label_universe / usable_intersection / class_A/B.

    Joins the per-cancer expression matrix (column = TCGA aliquot id) with the
    driver/subtype label table on patient barcode (TCGA-XX-XXXX). Used by the
    cohort_availability writer so the table is arithmetically consistent
    (audit fix F2): class_A + class_B == usable_intersection by construction.
    """
    expr_p = DATA_RAW_V5 / cancer / expr_name
    lab_p = DATA_RAW_V5 / cancer / lab_name
    if not expr_p.exists() or not lab_p.exists():
        return dict(total_downloaded=0, label_universe=0,
                    usable_intersection=0, class_A=0, class_B=0)
    open_fn = gzip.open if str(expr_p).endswith(".gz") else open
    with open_fn(expr_p, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    samples = header[1:]
    expr_patients = set(s[:12] for s in samples)
    lab = pd.read_csv(lab_p, sep="\t")
    lab_A_set = set(lab.loc[lab["label"] == name_A, "patient_id"].astype(str))
    lab_B_set = set(lab.loc[lab["label"] == name_B, "patient_id"].astype(str))
    label_universe = lab_A_set | lab_B_set
    usable = expr_patients & label_universe
    return dict(
        total_downloaded=len(samples),
        label_universe=len(label_universe),
        usable_intersection=len(usable),
        class_A=len(usable & lab_A_set),
        class_B=len(usable & lab_B_set),
    )


def _row_tcga(cancer: str, cohort_id: str, name_A: str, name_B: str,
              expr_name: str, lab_name: str, n_genes: int,
              status: str, err: str) -> Dict:
    f = _tcga_usable_intersection(cancer, name_A, name_B, expr_name, lab_name)
    # New schema columns + back-compat aliases (n_samples/n_class_A/n_class_B
    # remain as aliases of total_downloaded/class_A/class_B so existing call
    # sites keep working).
    return dict(
        cancer=cancer, cohort_type="tcga", cohort_id=cohort_id,
        total_downloaded=f["total_downloaded"],
        label_universe=f["label_universe"],
        usable_intersection=f["usable_intersection"],
        class_A=f["class_A"], class_B=f["class_B"],
        n_shared_genes=n_genes,
        status=status, download_error=err, semi_synthetic=False,
        n_samples=f["total_downloaded"],
        n_class_A=f["class_A"], n_class_B=f["class_B"],
    )


def _row_geo(cancer: str, gse: str, g: Dict) -> Dict:
    """GEO row. Labels are derived from the series matrix itself, so
    label_universe == usable_intersection == n_class_A + n_class_B."""
    n = int(g.get("n_samples", 0) or 0)
    a = int(g.get("n_class_A", 0) or 0)
    b = int(g.get("n_class_B", 0) or 0)
    status = "included" if g.get("status") == "ok" else g["status"]
    err = g.get("error", g.get("note", "")) or ""
    return dict(
        cancer=cancer, cohort_type="geo", cohort_id=gse,
        total_downloaded=n,
        label_universe=a + b,
        usable_intersection=a + b,
        class_A=a, class_B=b,
        n_shared_genes=g.get("n_genes", 0),
        status=status, download_error=err, semi_synthetic=False,
        n_samples=n, n_class_A=a, n_class_B=b,
    )


COHORT_AVAIL_COLS = [
    "cancer", "cohort_type", "cohort_id",
    "total_downloaded", "label_universe", "usable_intersection",
    "class_A", "class_B",
    "n_shared_genes", "status", "download_error", "semi_synthetic",
    # legacy aliases retained for downstream consumers
    "n_samples", "n_class_A", "n_class_B",
]


def main():
    log_line(LOG, "=== v5p1 Phase 1 START ===")
    rows = []

    # --- THCA shortcut ---
    r = thca_shortcut()
    thca_status = "included" if r["status"] == "ok" else r["status"]
    rows.append(_row_tcga(
        cancer="THCA", cohort_id="TCGA-THCA",
        name_A="BRAF", name_B="RAS",
        expr_name="tcga_log2tpm.tsv.gz", lab_name="tcga_samples_labels.tsv",
        n_genes=r.get("n_genes", 0),
        status=thca_status, err="",
    ))
    log_line(LOG, f"[THCA] {r}")

    # THCA GEO cohorts (GSE27155, GSE33630, GSE29265) — real labels
    for gse in COHORTS["THCA"]["geo"]:
        # check cached local first (project/data_raw/geo), fall back to fetch
        local = PROJECT / "data_raw" / "geo" / gse / f"{gse}_series_matrix.txt.gz"
        target = DATA_RAW_V5 / "THCA" / "geo" / f"{gse}_series_matrix.txt.gz"
        target.parent.mkdir(parents=True, exist_ok=True)
        if local.exists() and not target.exists():
            import shutil
            shutil.copy(local, target)
        g = download_geo_series("THCA", gse)
        rows.append(_row_geo("THCA", gse, g))
        log_line(LOG, f"[THCA GEO {gse}] {g}")

    # --- 4 other cancers ---
    for cancer in ["SKCM", "LGG", "LUAD", "COAD"]:
        log_line(LOG, f"=== {cancer} START ===")
        # RNA-seq
        rna = download_tcga_rnaseq(cancer)
        # MAF
        maf = download_tcga_maf(cancer)
        # Driver labels
        labels_df = extract_driver_labels(cancer) if maf.get("status") == "ok" else None

        tcga_n_a = int((labels_df["label"] == COHORTS[cancer]["class_a"]).sum()) if labels_df is not None else 0
        tcga_n_b = int((labels_df["label"] == COHORTS[cancer]["class_b"]).sum()) if labels_df is not None else 0

        tcga_status = "included"
        tcga_err = ""
        if rna.get("status") != "ok":
            tcga_status = "excluded_download_fail"
            tcga_err = "RNA-seq download failed"
        elif maf.get("status") != "ok":
            tcga_status = "excluded_download_fail"
            tcga_err = "MAF download failed"
        elif labels_df is None or tcga_n_a < 30 or tcga_n_b < 30:
            tcga_status = "excluded_small_n"
            tcga_err = f"n_A={tcga_n_a} n_B={tcga_n_b} (<30)"

        # Save labels TSV for this cancer
        if labels_df is not None:
            labels_df.to_csv(DATA_RAW_V5 / cancer / f"tcga_driver_labels.tsv",
                             sep="\t", index=False)

        # NOTE: tcga_n_a / tcga_n_b above are MAF-derived patient label counts
        # (label_universe-side numbers). The cohort_availability row uses
        # _row_tcga which joins the expression matrix with the label set so
        # class_A + class_B = usable_intersection holds (audit fix F2).
        rows.append(_row_tcga(
            cancer=cancer, cohort_id=COHORTS[cancer]["tcga"],
            name_A=COHORTS[cancer]["class_a"], name_B=COHORTS[cancer]["class_b"],
            expr_name="tcga_counts.tsv.gz", lab_name="tcga_driver_labels.tsv",
            n_genes=rna.get("n_genes", 0),
            status=tcga_status, err=tcga_err,
        ))
        log_line(LOG, f"[{cancer} TCGA] status={tcga_status} maf_n_A={tcga_n_a} maf_n_B={tcga_n_b}")

        # GEO cohorts
        for gse in COHORTS[cancer]["geo"]:
            g = download_geo_series(cancer, gse)
            rows.append(_row_geo(cancer, gse, g))
            log_line(LOG, f"[{cancer} GEO {gse}] {g}")

    df = pd.DataFrame(rows)
    # Reorder columns to the canonical schema (new columns first, legacy aliases last).
    df = df[[c for c in COHORT_AVAIL_COLS if c in df.columns]]
    out_path = RESULTS_V5 / "v5p1_cohort_availability.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    log_line(LOG, f"wrote {out_path}")

    # Also emit driver-label summary for TCGA
    driver_rows = []
    for cancer in ["SKCM", "LGG", "LUAD", "COAD"]:
        lp = DATA_RAW_V5 / cancer / "tcga_driver_labels.tsv"
        if lp.exists():
            lab = pd.read_csv(lp, sep="\t")
            for lname, cnt in lab["label"].value_counts().items():
                driver_rows.append(dict(cancer=cancer, label=lname, n_patients=int(cnt)))
    pd.DataFrame(driver_rows).to_csv(RESULTS_V5 / "v5p1_tcga_mutation_groups.tsv",
                                     sep="\t", index=False)
    log_line(LOG, "=== v5p1 Phase 1 DONE ===")


if __name__ == "__main__":
    main()
