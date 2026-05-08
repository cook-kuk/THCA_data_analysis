#!/usr/bin/env python3
"""Phase C v2 — harmonize 4 ICI cohorts to gene-symbol expression matrices + metadata.

Inputs (each cohort under phase_C_ICI/raw/{cohort}/):
  IMvigor210:        extracted/{counts,sample_phenotype,feature_annotation}.tsv
  GSE176307:         _log_trans_normalized_RNAseq.csv.gz  +  series_matrix.txt.gz  +  sample key
  riaz_GSE91061:     _hg19KnownGene.fpkm.csv.gz           +  series_matrix.txt.gz
  gide_PRJEB23709:   GSE115821_MGH_counts.csv.gz          +  GSE115821_family.soft.gz
                     (substituted MGH cohort — see logs/cohort_substitution_note.md)

Outputs (under phase_C_ICI/processed/{cohort}/):
  expr_matrix.tsv    log2-normalized expression, gene symbols × sample_id
  metadata.tsv       per-sample harmonized metadata (response/survival/etc.)

Metadata schema (subset present per cohort):
  sample_id patient_id cohort cancer_type therapy timepoint
  response_raw response_binary_CRPR_vs_SD_PD disease_control_CRPRSD_vs_PD
  os_time os_event pfs_time pfs_event source_file notes
"""
from __future__ import annotations
import gzip, re, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/data/thca/repo_results/paper11_pancancer/phase_C_ICI")
RAW = ROOT / "raw"
PROC = ROOT / "processed"
LOGS = ROOT / "logs"
LOGS.mkdir(exist_ok=True, parents=True)

LOG = []
def log(msg):
    print(msg)
    LOG.append(msg)


# ------------------------------------------------------------------ utils
def log2_cpm(counts_df: pd.DataFrame, pseudo: float = 1.0) -> pd.DataFrame:
    """counts_df: genes × samples (numeric). Returns log2(CPM+pseudo)."""
    lib = counts_df.sum(axis=0)
    cpm = counts_df.div(lib, axis=1) * 1e6
    return np.log2(cpm + pseudo)


def collapse_dups(df: pd.DataFrame, by: str = "index", how: str = "max_var") -> pd.DataFrame:
    """Collapse duplicate gene rows by max-variance (keeps the most informative copy)."""
    if df.index.duplicated().any():
        var = df.var(axis=1)
        # for each duplicate group, keep highest-var row
        keep = var.groupby(df.index).idxmax()
        df = df.loc[keep.values]
    return df


# ------------------------------------------------------------------ Entrez→Symbol lookup (NCBI Homo_sapiens.gene_info subset)
# We piggy-back on the IMvigor feature_annotation.tsv (entrez_id → symbol) since both
# IMvigor and Riaz use the same Entrez Gene IDs. For genes Riaz has but IMvigor doesn't,
# we fall back to mygene.info via local cache or skip with a warning.
def load_entrez_symbol_map() -> dict[str, str]:
    fa = RAW / "IMvigor210/extracted/feature_annotation.tsv"
    if not fa.exists():
        log(f"!!! missing {fa} — entrez map will be empty")
        return {}
    df = pd.read_csv(fa, sep="\t", dtype=str)
    sym_col = "symbol" if "symbol" in df.columns else "Symbol"
    m = dict(zip(df["entrez_id"].astype(str), df[sym_col].astype(str)))
    # drop entries where symbol is NA / empty
    m = {k: v for k, v in m.items() if v and v != "nan"}
    log(f"[entrez_map] loaded {len(m)} entrez→symbol pairs")
    return m


# ------------------------------------------------------------------ IMvigor210
def harmonize_imvigor() -> None:
    cohort = "IMvigor210"
    out = PROC / cohort
    out.mkdir(exist_ok=True, parents=True)
    extracted = RAW / cohort / "extracted"

    counts = pd.read_csv(extracted / "counts.tsv", sep="\t", index_col="entrez_id")
    feat = pd.read_csv(extracted / "feature_annotation.tsv", sep="\t", dtype=str).set_index("entrez_id")
    pheno = pd.read_csv(extracted / "sample_phenotype.tsv", sep="\t")

    log(f"[IMvigor] counts {counts.shape}, pheno {pheno.shape}")

    # --- expression: CPM + log2, then collapse to symbols (max-var)
    sym_col = "symbol" if "symbol" in feat.columns else "Symbol"
    counts.index = counts.index.astype(str)
    feat.index = feat.index.astype(str)
    sym_map = feat[sym_col].dropna().astype(str)
    counts_sym = counts.loc[counts.index.isin(sym_map.index)].copy()
    counts_sym.index = sym_map.loc[counts_sym.index].values
    # drop NA / empty symbols
    counts_sym = counts_sym.loc[(counts_sym.index.notna()) & (counts_sym.index != "") & (counts_sym.index != "nan")]
    counts_sym = collapse_dups(counts_sym)
    expr = log2_cpm(counts_sym)
    expr.index.name = "gene_symbol"
    expr.to_csv(out / "expr_matrix.tsv", sep="\t")
    log(f"[IMvigor] expr_matrix.tsv {expr.shape} (genes × samples)")

    # --- metadata
    md = pd.DataFrame()
    md["sample_id"] = pheno["sample_id"]
    md["patient_id"] = pheno.get("ANONPT_ID", pheno["sample_id"])
    md["cohort"] = cohort
    md["cancer_type"] = "urothelial_carcinoma"
    md["therapy"] = "atezolizumab_anti_PD_L1"
    md["timepoint"] = "pre"   # IMvigor210 phase II uses pre-treatment biopsy
    md["response_raw"] = pheno["Best Confirmed Overall Response"]
    bin_resp = pheno["binaryResponse"].astype(str).str.strip()
    md["response_binary_CRPR_vs_SD_PD"] = bin_resp.map(
        {"CR/PR": 1, "SD/PD": 0}).where(bin_resp.isin(["CR/PR", "SD/PD"]), np.nan)
    md["disease_control_CRPRSD_vs_PD"] = pheno["Best Confirmed Overall Response"].map(
        {"CR": 1, "PR": 1, "SD": 1, "PD": 0})
    md["os_time"] = pd.to_numeric(pheno["os"], errors="coerce")  # months
    md["os_event"] = pd.to_numeric(pheno["censOS"], errors="coerce")  # 0=alive, 1=dead
    md["pfs_time"] = np.nan
    md["pfs_event"] = np.nan
    md["source_file"] = "cds.RData (IMvigor210CoreBiologies_1.0.0)"
    md["notes"] = (
        "FMOne TMB: " + pheno.get("FMOne mutation burden per MB", pd.Series(["NA"]*len(pheno))).astype(str)
        + "; IC: " + pheno.get("IC Level", pd.Series([""]*len(pheno))).astype(str)
        + "; TCGA_subtype: " + pheno.get("TCGA Subtype", pd.Series([""]*len(pheno))).astype(str)
    )
    md.to_csv(out / "metadata.tsv", sep="\t", index=False)
    log(f"[IMvigor] metadata.tsv n={len(md)}, n_with_resp={md['response_binary_CRPR_vs_SD_PD'].notna().sum()}")


# ------------------------------------------------------------------ GSE176307 (urothelial real-world)
def harmonize_gse176307() -> None:
    cohort = "GSE176307"
    out = PROC / cohort
    out.mkdir(exist_ok=True, parents=True)
    raw = RAW / cohort

    # log_trans is already log2-normalized; rows = BACI samples, cols = symbols
    log_df = pd.read_csv(raw / "GSE176307_BACI_log_trans_normalized_RNAseq.csv.gz",
                         index_col=0, compression="gzip", low_memory=False)
    log(f"[GSE176307] log_trans shape: {log_df.shape}  (samples × genes)")
    expr = log_df.T  # gene × sample
    expr = collapse_dups(expr)
    expr.index.name = "gene_symbol"
    expr.to_csv(out / "expr_matrix.tsv", sep="\t")
    log(f"[GSE176307] expr_matrix.tsv {expr.shape}")

    # series_matrix → metadata
    sm_path = raw / "GSE176307_series_matrix.txt.gz"
    sm = {}
    with gzip.open(sm_path, "rt", errors="ignore") as f:
        for line in f:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0].lstrip("!")
            vals = [v.strip().strip('"') for v in parts[1:]]
            sm.setdefault(key, []).append(vals)

    titles = sm["Sample_title"][0]
    geos   = sm["Sample_geo_accession"][0]
    chars = sm["Sample_characteristics_ch1"]  # list-of-lists

    # parse characteristics into named cols
    n = len(titles)
    char_cols = {}
    for row in chars:
        keys = [c.split(":")[0].strip() if ":" in c else "unknown" for c in row]
        # use most frequent key as the column name
        if not keys: continue
        key = max(set(keys), key=keys.count)
        vals = [c.split(":", 1)[1].strip() if ":" in c else c for c in row]
        if key in char_cols:
            key = key + "_alt"
        char_cols[key] = vals

    md = pd.DataFrame({"sample_id": titles, "geo_accession": geos})
    for k, v in char_cols.items():
        md[k.replace(".", "_").replace(" ", "_").lower()] = v

    # The BACI patient ID in titles is "Patient sample BACI107" → BACI107
    md["patient_id"] = md["sample_id"].str.extract(r"(BACI\d+)", expand=False)
    md["sample_id"] = md["patient_id"]  # use BACI as sample id (matches expr_matrix row labels)
    md["cohort"] = cohort
    md["cancer_type"] = "urothelial_carcinoma"
    md["therapy"] = md.get("io_therapy", "ICB")
    md["timepoint"] = "pre"

    # response
    resp = md.get("io_response", pd.Series([np.nan]*n))
    md["response_raw"] = resp
    bin_map = {"CR": 1, "PR": 1, "SD": 0, "PD": 0}
    md["response_binary_CRPR_vs_SD_PD"] = resp.map(bin_map)
    md["disease_control_CRPRSD_vs_PD"] = resp.map({"CR": 1, "PR": 1, "SD": 1, "PD": 0})

    md["os_time"] = pd.to_numeric(md.get("overall_survival", pd.Series([np.nan]*n)), errors="coerce")
    alive = md.get("alive", pd.Series([""]*n)).astype(str).str.lower()
    md["os_event"] = (alive == "no").astype(int)
    md["pfs_time"] = pd.to_numeric(md.get("pfs", pd.Series([np.nan]*n)), errors="coerce")
    progressed = md.get("progressed", pd.Series([""]*n)).astype(str).str.lower()
    md["pfs_event"] = (progressed == "yes").astype(int)

    md["source_file"] = "GSE176307_BACI_log_trans_normalized_RNAseq.csv.gz"
    md["notes"] = "real-world urothelial ICB cohort"

    keep_cols = ["sample_id", "patient_id", "cohort", "cancer_type", "therapy", "timepoint",
                 "response_raw", "response_binary_CRPR_vs_SD_PD", "disease_control_CRPRSD_vs_PD",
                 "os_time", "os_event", "pfs_time", "pfs_event", "source_file", "notes",
                 "geo_accession", "io_therapy"]
    md = md[[c for c in keep_cols if c in md.columns]]

    # filter expr to samples that exist in metadata
    common = [s for s in expr.columns if s in md["sample_id"].values]
    expr_filtered = expr[common]
    expr_filtered.to_csv(out / "expr_matrix.tsv", sep="\t")  # overwrite with filtered
    md = md[md["sample_id"].isin(common)]
    md.to_csv(out / "metadata.tsv", sep="\t", index=False)
    log(f"[GSE176307] metadata n={len(md)}, n_with_resp={md['response_binary_CRPR_vs_SD_PD'].notna().sum()}, expr {expr_filtered.shape}")


# ------------------------------------------------------------------ Riaz GSE91061
def harmonize_riaz(entrez_map: dict[str, str]) -> None:
    cohort = "riaz_GSE91061"
    out = PROC / cohort
    out.mkdir(exist_ok=True, parents=True)
    raw = RAW / cohort

    fpkm = pd.read_csv(raw / "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz",
                       index_col=0, low_memory=False)
    log(f"[Riaz] fpkm raw: {fpkm.shape}")
    fpkm.index = fpkm.index.astype(str)
    # Map entrez → symbol
    sym = pd.Series(fpkm.index).map(entrez_map)
    keep = sym.notna() & (sym != "") & (sym != "nan")
    fpkm_sym = fpkm.loc[keep.values].copy()
    fpkm_sym.index = sym[keep].values
    fpkm_sym = collapse_dups(fpkm_sym)
    log(f"[Riaz] fpkm w/ symbols: {fpkm_sym.shape}")

    # log2(FPKM + 1)
    expr = np.log2(fpkm_sym + 1)
    expr.index.name = "gene_symbol"

    # series_matrix metadata
    sm_path = raw / "GSE91061_series_matrix.txt.gz"
    sm = {}
    with gzip.open(sm_path, "rt", errors="ignore") as f:
        for line in f:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0].lstrip("!")
            vals = [v.strip().strip('"') for v in parts[1:]]
            sm.setdefault(key, []).append(vals)
    titles = sm["Sample_title"][0]
    geos = sm["Sample_geo_accession"][0]
    chars = sm["Sample_characteristics_ch1"]

    n = len(titles)
    char_cols = {}
    for row in chars:
        keys = [c.split(":")[0].strip() if ":" in c else "unknown" for c in row]
        if not keys: continue
        key = max(set(keys), key=keys.count)
        vals = [c.split(":", 1)[1].strip() if ":" in c else c for c in row]
        if key in char_cols:
            key = key + "_alt"
        char_cols[key] = vals

    md = pd.DataFrame({"sample_id": titles, "geo_accession": geos})
    for k, v in char_cols.items():
        md[k.replace(".", "_").replace(" ", "_").lower().replace("(", "").replace(")", "")] = v

    md["patient_id"] = md["sample_id"].str.extract(r"^(Pt\d+)", expand=False)
    md["cohort"] = cohort
    md["cancer_type"] = "melanoma"
    md["therapy"] = "nivolumab_anti_PD_1"
    visit_col = next((c for c in md.columns if c.startswith("visit")), None)
    md["timepoint"] = md[visit_col].str.lower() if visit_col else "unknown"

    resp = md.get("response", pd.Series([np.nan]*n)).astype(str)
    md["response_raw"] = resp
    md["response_binary_CRPR_vs_SD_PD"] = resp.map({"PRCR": 1, "PD": 0, "SD": 0}).where(resp.isin(["PRCR","PD","SD"]), np.nan)
    md["disease_control_CRPRSD_vs_PD"] = resp.map({"PRCR": 1, "SD": 1, "PD": 0}).where(resp.isin(["PRCR","PD","SD"]), np.nan)

    md["os_time"] = np.nan  # not in series matrix; would need patient-level supplemental
    md["os_event"] = np.nan
    md["pfs_time"] = np.nan
    md["pfs_event"] = np.nan
    md["source_file"] = "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"
    md["notes"] = "Riaz 2017; pre/on treatment paired RNA-seq"

    keep_cols = ["sample_id", "patient_id", "cohort", "cancer_type", "therapy", "timepoint",
                 "response_raw", "response_binary_CRPR_vs_SD_PD", "disease_control_CRPRSD_vs_PD",
                 "os_time", "os_event", "pfs_time", "pfs_event", "source_file", "notes",
                 "geo_accession"]
    md = md[[c for c in keep_cols if c in md.columns]]

    # restrict expr to samples in metadata (titles are sample_ids)
    common = [s for s in expr.columns if s in md["sample_id"].values]
    expr = expr[common]
    expr.to_csv(out / "expr_matrix.tsv", sep="\t")
    md = md[md["sample_id"].isin(common)]
    md.to_csv(out / "metadata.tsv", sep="\t", index=False)
    log(f"[Riaz] expr {expr.shape}, metadata n={len(md)} (Pre={sum(md.timepoint=='pre')}, On={sum(md.timepoint=='on')}), "
        f"n_with_resp={md['response_binary_CRPR_vs_SD_PD'].notna().sum()}")


# ------------------------------------------------------------------ MGH GSE115821 (substituted for Gide)
def harmonize_mgh_gse115821() -> None:
    cohort = "MGH_GSE115821"
    out = PROC / cohort
    out.mkdir(exist_ok=True, parents=True)
    raw = RAW / "gide_PRJEB23709"

    counts = pd.read_csv(raw / "GSE115821_MGH_counts.csv.gz", low_memory=False)
    # Geneid is the first col, then Chr/Start/End/Strand/Length, then samples
    sample_cols = [c for c in counts.columns if c.endswith(".bam")]
    counts_only = counts[["Geneid"] + sample_cols].copy()
    counts_only["Geneid"] = counts_only["Geneid"].astype(str).str.strip().str.replace("﻿", "", regex=False)
    counts_only = counts_only.set_index("Geneid")
    counts_only = counts_only.apply(pd.to_numeric, errors="coerce")
    counts_only = counts_only.dropna(how="all")
    counts_only = collapse_dups(counts_only)
    log(f"[MGH] counts: {counts_only.shape} (genes × bams)")

    expr = log2_cpm(counts_only)
    expr.index.name = "gene_symbol"

    # SOFT metadata
    soft = raw / "GSE115821_family.soft.gz"
    rows = []
    cur = {}
    with gzip.open(soft, "rt", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                if cur:
                    rows.append(cur)
                cur = {"gsm": line.split("=", 1)[1].strip()}
            elif line.startswith("!Sample_title"):
                cur["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1"):
                v = line.split("=", 1)[1].strip()
                if ":" in v:
                    k, val = [x.strip() for x in v.split(":", 1)]
                    cur[k.lower().replace(" ", "_")] = val
    if cur:
        rows.append(cur)
    md_soft = pd.DataFrame(rows)
    log(f"[MGH] SOFT parsed: {len(md_soft)} samples, cols={list(md_soft.columns)[:10]}")

    # Map bam columns to titles. Sample titles like "115-031814" should match a bam stem
    title_to_bam = {}
    for bam in sample_cols:
        stem = bam.replace(".bam", "")
        # find best match in titles
        for t in md_soft["title"]:
            if stem == t or stem.replace("_", "-") == t or stem == t.replace("-", "_"):
                title_to_bam[t] = bam
                break
            # try fuzzy: starts-with
            if stem.startswith(str(t)) or str(t).startswith(stem):
                title_to_bam.setdefault(t, bam)
    log(f"[MGH] mapped {len(title_to_bam)} of {len(md_soft)} titles to bams")

    md_soft["bam_col"] = md_soft["title"].map(title_to_bam)
    matched = md_soft.dropna(subset=["bam_col"]).reset_index(drop=True)

    # rename expr columns from bam → title
    bam_to_title = {v: k for k, v in title_to_bam.items()}
    expr = expr.rename(columns=bam_to_title)
    expr = expr[[c for c in expr.columns if c in matched["title"].values]]

    # build harmonized metadata
    md = pd.DataFrame()
    md["sample_id"] = matched["title"]
    md["patient_id"] = matched.get("patient_id", matched["title"])
    md["cohort"] = cohort
    md["cancer_type"] = "melanoma"
    therapy_col = matched.get("antibody", pd.Series([""]*len(matched))).astype(str).str.lower()
    md["therapy"] = therapy_col
    ts = matched.get("treatment_state", pd.Series(["?"]*len(matched))).astype(str).str.lower()
    md["timepoint"] = ts.where(ts.str.contains("pre"), "on").where(~ts.str.contains("pre"), "pre")
    resp = matched.get("response", pd.Series([np.nan]*len(matched))).astype(str)
    md["response_raw"] = resp
    md["response_binary_CRPR_vs_SD_PD"] = resp.map({"R": 1, "NR": 0}).where(resp.isin(["R","NR"]), np.nan)
    md["disease_control_CRPRSD_vs_PD"] = md["response_binary_CRPR_vs_SD_PD"]  # binary R/NR collapses to single binary
    md["os_time"] = np.nan
    md["os_event"] = np.nan
    md["pfs_time"] = np.nan
    md["pfs_event"] = np.nan
    md["source_file"] = "GSE115821_MGH_counts.csv.gz"
    md["notes"] = "MGH (Auslander/Liu et al.) melanoma ICB; substituted for Gide PRJEB23709 (FASTQ-only)"

    expr.to_csv(out / "expr_matrix.tsv", sep="\t")
    md.to_csv(out / "metadata.tsv", sep="\t", index=False)
    log(f"[MGH] expr {expr.shape}, metadata n={len(md)}, n_with_resp={md['response_binary_CRPR_vs_SD_PD'].notna().sum()}")


# ------------------------------------------------------------------ main
def main() -> None:
    entrez_map = load_entrez_symbol_map()
    harmonize_imvigor()
    harmonize_gse176307()
    harmonize_riaz(entrez_map)
    harmonize_mgh_gse115821()

    (LOGS / "01_harmonize_log.txt").write_text("\n".join(LOG))
    log(f"\n[done] log → {LOGS/'01_harmonize_log.txt'}")


if __name__ == "__main__":
    main()
