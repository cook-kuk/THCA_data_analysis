#!/usr/bin/env python3
"""
H4 — HT-13 panel cross-cohort meta-direction validation.

HT-13 panel = HLA-DRA, HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1,
              CD79A, CD79B, MS4A1, AICDA, CXCL13, CCR6, IFNG.

For each cohort, compute mean(z(panel)) per sample, then Cohen's d (DM1/HT-overlap/aggressive vs comparator).

Cohorts:
  1. TCGA-THCA full primary tumor   (DM1 vs DM2 + by molecular_subtype)
  2. GSE286332 Korean PTC vs PTC+HT
  3. GSE76039 Landa ATC vs PDTC
  4. Mun 2025 proteomics (ATC vs other)
  5. (optional) GSE151179 post-RAI vs pre-RAI
"""
import gzip
import io
import json
import math
import os
import re
import sys
import warnings
from collections import OrderedDict, defaultdict

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings('ignore')

OUT = "/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/h4_meta_4cohort"
os.makedirs(OUT, exist_ok=True)

PANEL = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
         "CD79A", "CD79B", "MS4A1", "AICDA", "CXCL13", "CCR6", "IFNG"]

# canonical underscore form (Mun protein matrix uses '.' or '_' — handle robustly)
PANEL_ALT = {g: [g, g.replace("-", "."), g.replace("-", "_")] for g in PANEL}

# ----------------- helpers -----------------

def cohens_d(a, b):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan, len(a), len(b)
    s2a, s2b = a.var(ddof=1), b.var(ddof=1)
    sp = math.sqrt(((len(a)-1)*s2a + (len(b)-1)*s2b) / (len(a)+len(b)-2))
    if sp == 0:
        return np.nan, np.nan, len(a), len(b)
    d = (a.mean() - b.mean()) / sp
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return d, p, len(a), len(b)


def zscore_rows(df):
    """Row-wise z-score of a gene x sample matrix."""
    mu = df.mean(axis=1)
    sd = df.std(axis=1, ddof=0).replace(0, np.nan)
    return df.subtract(mu, axis=0).divide(sd, axis=0)


def panel_mean_score(expr_mat, panel_genes):
    """expr_mat: gene x sample (already row-zscored or comparable). Returns sample-level mean of available panel members."""
    avail = [g for g in panel_genes if g in expr_mat.index]
    sub = expr_mat.loc[avail]
    return sub.mean(axis=0), avail


def per_gene_cohens_d(expr_mat, group_a_samples, group_b_samples, genes):
    rows = []
    for g in genes:
        if g not in expr_mat.index:
            rows.append({"gene": g, "d": np.nan, "p": np.nan, "n_a": 0, "n_b": 0})
            continue
        a = expr_mat.loc[g, expr_mat.columns.intersection(group_a_samples)].values
        b = expr_mat.loc[g, expr_mat.columns.intersection(group_b_samples)].values
        d, p, na, nb = cohens_d(a, b)
        rows.append({"gene": g, "d": d, "p": p, "n_a": na, "n_b": nb})
    return rows


# ----------------- cohort 1: TCGA-THCA -----------------

def run_tcga():
    print("[TCGA-THCA] loading...")
    expr_path = "/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv"
    df = pd.read_csv(expr_path, sep="\t", index_col=0)
    # filter to panel genes only (faster)
    panel_present = [g for g in PANEL if g in df.index]
    expr = df.loc[panel_present].copy()
    print(f"  panel members present: {len(panel_present)}/{len(PANEL)}: {panel_present}")

    master = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv",
                         sep="\t")
    master = master[master["sample_id"].isin(expr.columns)].copy()

    score, _ = panel_mean_score(expr, PANEL)
    score = score.to_frame("ht13_score")
    score = score.merge(master[["sample_id", "dm", "molecular_subtype", "histology_subtype"]],
                        left_index=True, right_on="sample_id")

    out_rows = []
    pg_rows = []

    # overall DM1 vs DM2
    a = score.loc[score["dm"] == "DM1", "ht13_score"].values
    b = score.loc[score["dm"] == "DM2", "ht13_score"].values
    d, p, na, nb = cohens_d(a, b)
    out_rows.append({
        "cohort": "TCGA-THCA",
        "stratum": "all_primary_tumor",
        "contrast": "DM1_vs_DM2",
        "n_high": na, "n_low": nb,
        "mean_d_HT13": d, "p": p,
        "sign": "+" if d > 0.1 else ("-" if d < -0.1 else "0"),
    })
    a_ids = score.loc[score["dm"] == "DM1", "sample_id"].tolist()
    b_ids = score.loc[score["dm"] == "DM2", "sample_id"].tolist()
    for r in per_gene_cohens_d(expr, a_ids, b_ids, PANEL):
        r["cohort"] = "TCGA-THCA"; r["stratum"] = "all"; r["contrast"] = "DM1_vs_DM2"
        pg_rows.append(r)

    # by molecular_subtype
    for subt in sorted(score["molecular_subtype"].dropna().unique()):
        sub = score[score["molecular_subtype"] == subt]
        a = sub.loc[sub["dm"] == "DM1", "ht13_score"].values
        b = sub.loc[sub["dm"] == "DM2", "ht13_score"].values
        d, p, na, nb = cohens_d(a, b)
        if na < 3 or nb < 3:
            continue
        out_rows.append({
            "cohort": "TCGA-THCA",
            "stratum": f"subtype={subt}",
            "contrast": "DM1_vs_DM2",
            "n_high": na, "n_low": nb,
            "mean_d_HT13": d, "p": p,
            "sign": "+" if d > 0.1 else ("-" if d < -0.1 else "0"),
        })
    return out_rows, pg_rows


# ----------------- cohort 2: GSE286332 -----------------

def run_gse286332():
    print("[GSE286332] loading...")
    raw = "/data/thca/v17_korean/GSE286332/GSE286332_all_sample_rawdata.txt.gz"
    df = pd.read_csv(raw, sep="\t", compression="gzip", low_memory=False)
    # Use TPM columns
    tpm_cols = [c for c in df.columns if c.endswith("_TPM")]
    samples = [c.replace("_TPM", "") for c in tpm_cols]
    expr = df.set_index("Gene_Symbol")[tpm_cols].copy()
    expr.columns = samples
    # drop duplicate symbols (keep first)
    expr = expr[~expr.index.duplicated(keep="first")]
    # log2(tpm+1)
    expr = np.log2(expr.astype(float) + 1)
    # row z
    expr_z = zscore_rows(expr)

    panel_present = [g for g in PANEL if g in expr_z.index]
    print(f"  panel members present: {len(panel_present)}/{len(PANEL)}: {panel_present}")
    score, _ = panel_mean_score(expr_z, PANEL)

    # NG = PTC, TH = PTC+HT (per P3_summary.json: n_PTC=9, n_PTC_HT=9)
    a_ids = [s for s in samples if s.startswith("TH_")]   # PTC+HT
    b_ids = [s for s in samples if s.startswith("NG_")]   # PTC alone
    a = score.loc[score.index.intersection(a_ids)].values
    b = score.loc[score.index.intersection(b_ids)].values
    d, p, na, nb = cohens_d(a, b)

    out = [{
        "cohort": "GSE286332",
        "stratum": "Korean_PTC_HT_overlap",
        "contrast": "PTC+HT_vs_PTC",
        "n_high": na, "n_low": nb,
        "mean_d_HT13": d, "p": p,
        "sign": "+" if d > 0.1 else ("-" if d < -0.1 else "0"),
    }]
    pg = []
    for r in per_gene_cohens_d(expr_z, a_ids, b_ids, PANEL):
        r["cohort"] = "GSE286332"; r["stratum"] = "all"; r["contrast"] = "PTC+HT_vs_PTC"
        pg.append(r)
    return out, pg


# ----------------- cohort 3: GSE76039 -----------------

def run_gse76039():
    print("[GSE76039] loading...")
    expr = pd.read_csv("/data/thca/data_processed/microarray/GSE76039_microarray_expression_zscore.tsv",
                       sep="\t", index_col=0)
    # already row-zscored? check
    print(f"  expr shape: {expr.shape}, mean of row 0: {expr.iloc[0].mean():.3f}")

    meta = pd.read_csv("/data/thca/repo_results/p_deconv_2026_05_08/v11_GSE76039_sample_scores.tsv",
                       sep="\t")
    meta = meta.set_index("sample")
    panel_present = [g for g in PANEL if g in expr.index]
    print(f"  panel members present: {len(panel_present)}/{len(PANEL)}: {panel_present}")

    score, _ = panel_mean_score(expr, PANEL)

    atc_ids = meta.index[meta["hist_short"] == "ATC"].tolist()
    pdtc_ids = meta.index[meta["hist_short"] == "PDTC"].tolist()
    a = score.loc[score.index.intersection(atc_ids)].values
    b = score.loc[score.index.intersection(pdtc_ids)].values
    d, p, na, nb = cohens_d(a, b)

    out = [{
        "cohort": "GSE76039",
        "stratum": "Landa_dediff",
        "contrast": "ATC_vs_PDTC",
        "n_high": na, "n_low": nb,
        "mean_d_HT13": d, "p": p,
        "sign": "+" if d > 0.1 else ("-" if d < -0.1 else "0"),
    }]
    pg = []
    for r in per_gene_cohens_d(expr, atc_ids, pdtc_ids, PANEL):
        r["cohort"] = "GSE76039"; r["stratum"] = "all"; r["contrast"] = "ATC_vs_PDTC"
        pg.append(r)
    return out, pg


# ----------------- cohort 4: Mun 2025 proteomics -----------------

def run_mun2025():
    print("[Mun2025] loading proteomics S1C...")
    import openpyxl
    wb = openpyxl.load_workbook(
        "/data/thca/repo_results/proteogenomic_v1/_raw/mun2025_MOESM3.xlsx",
        read_only=True, data_only=True)
    ws = wb["Table S1C"]
    rows = list(ws.iter_rows(values_only=True))
    wb.close()

    header = list(rows[0])  # ['Sample', sample IDs...]
    group_row = list(rows[1])  # ['Group', 'ATC', 'ATC', ...]
    # build df
    sample_ids = header[1:]
    groups = group_row[1:]
    print(f"  total samples: {len(sample_ids)}; group counts:", pd.Series(groups).value_counts().to_dict())

    # build gene matrix
    data_rows = rows[2:]
    gene_idx = []
    mat = np.zeros((len(data_rows), len(sample_ids)), dtype=float)
    for i, r in enumerate(data_rows):
        gene_idx.append(r[0])
        # values
        for j, v in enumerate(r[1:1+len(sample_ids)]):
            try:
                mat[i, j] = float(v) if v is not None and v != "" else np.nan
            except (TypeError, ValueError):
                mat[i, j] = np.nan
    expr = pd.DataFrame(mat, index=gene_idx, columns=sample_ids)
    expr = expr[~expr.index.duplicated(keep="first")]

    # gene names in S1C use '.' instead of '-' for HLA genes (e.g. 'HLA.DRA')
    rename = {}
    for g in expr.index:
        if isinstance(g, str) and g.startswith("HLA."):
            rename[g] = g.replace("HLA.", "HLA-", 1)
    expr = expr.rename(index=rename)

    # log2 transform if values look like ratios (Mun S1C is intensity ratio per memory; values look ~0.0-3.0)
    # Just z-score directly across samples for each protein
    expr_z = zscore_rows(np.log2(expr.replace(0, np.nan).fillna(method="ffill") + 1e-3))
    # safer: skip log if too many zeros — use raw z
    expr_z = zscore_rows(expr.fillna(expr.median(axis=1)))

    panel_present = [g for g in PANEL if g in expr_z.index]
    panel_missing = [g for g in PANEL if g not in expr_z.index]
    print(f"  panel proteins detected: {len(panel_present)}/{len(PANEL)}: {panel_present}")
    print(f"  missing: {panel_missing}")

    score, _ = panel_mean_score(expr_z, PANEL)
    grp = pd.Series(groups, index=sample_ids)
    out = []
    pg = []

    # contrast 1: ATC vs PTC
    for cmp in [("ATC", "PTC"), ("ATC", "Normal"), ("PDTC", "PTC")]:
        a_ids = grp.index[grp == cmp[0]].tolist()
        b_ids = grp.index[grp == cmp[1]].tolist()
        if len(a_ids) < 3 or len(b_ids) < 3:
            continue
        a = score.loc[score.index.intersection(a_ids)].values
        b = score.loc[score.index.intersection(b_ids)].values
        d, p, na, nb = cohens_d(a, b)
        out.append({
            "cohort": "Mun2025_protein",
            "stratum": "ATC_proteomics",
            "contrast": f"{cmp[0]}_vs_{cmp[1]}",
            "n_high": na, "n_low": nb,
            "mean_d_HT13": d, "p": p,
            "sign": "+" if d > 0.1 else ("-" if d < -0.1 else "0"),
        })
        for r in per_gene_cohens_d(expr_z, a_ids, b_ids, PANEL):
            r["cohort"] = "Mun2025_protein"; r["stratum"] = "all"; r["contrast"] = f"{cmp[0]}_vs_{cmp[1]}"
            pg.append(r)
    return out, pg


# ----------------- cohort 5 (optional): GSE151179 RAI -----------------

def run_gse151179():
    print("[GSE151179] parsing series matrix...")
    sm_path = "/data/thca/repo_results/paper3_ici_public_data/GSE151179/GSE151179_series_matrix.txt.gz"
    # read header and expression block
    samples = None
    with gzip.open(sm_path, "rt") as fh:
        lines = fh.readlines()
    # find table begin
    for i, line in enumerate(lines):
        if line.startswith("!Sample_geo_accession"):
            samples = [s.strip().strip('"') for s in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!series_matrix_table_begin"):
            data_start = i + 1
        if line.startswith("!series_matrix_table_end"):
            data_end = i
            break
    expr_lines = lines[data_start:data_end]
    print(f"  rows: {len(expr_lines)-1}, samples: {len(samples)}")
    # first row is header
    expr_io = io.StringIO("".join(expr_lines))
    df = pd.read_csv(expr_io, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.strip('"')
    # df is probe x sample. Need to map probes to gene symbols.
    # Use existing aggressive_sprint probe map to find probes for our 13 genes via SOFT family file.
    # quicker: grep family.soft.gz for our panel symbols.
    soft = "/data/thca/repo_results/aggressive_sprint_2026_05_06/geo/GSE151179_family.soft.gz"
    if not os.path.exists(soft):
        # fallback: just probe file from probe_map.tsv has only 8-gene mappings
        return [], []
    # parse soft for platform table — Clariom D uses long descriptive 10th column
    # we match the panel symbol when present as `(SYMBOL),` (RefSeq description) or `Acc:HGNC:...` near description
    print("  parsing SOFT for probe -> gene symbol map...")
    target_genes = set(PANEL)
    probe2gene = {}
    # regex: gene symbol immediately after "(" and before "), mRNA" — RefSeq format
    sym_re = re.compile(r"\(([A-Z][A-Z0-9-]*)\), mRNA")
    in_table = False
    with gzip.open(soft, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                _ = next(fh)  # header
                continue
            if in_table:
                if line.startswith("!platform_table_end"):
                    in_table = False
                    break
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                probe = f[0]
                desc = f[9]  # SPOT_ID full description
                if not desc:
                    continue
                # gather all matched gene symbols from the long string
                syms = set(sym_re.findall(desc))
                hits = syms & target_genes
                for g in hits:
                    probe2gene.setdefault(g, []).append(probe)
    print(f"  probes mapped per panel gene: " +
          ", ".join(f"{g}={len(v)}" for g, v in probe2gene.items()))

    # build gene-level expression by mean of probes
    gene_expr = {}
    for g, probes in probe2gene.items():
        rows = df.index.intersection(probes)
        if len(rows) == 0:
            continue
        gene_expr[g] = df.loc[rows].mean(axis=0)
    if not gene_expr:
        return [], []
    expr = pd.DataFrame(gene_expr).T
    expr.columns = [c.strip('"') for c in expr.columns]
    # row z
    expr_z = zscore_rows(expr)
    panel_present = list(expr_z.index)
    print(f"  panel members usable: {len(panel_present)}/{len(PANEL)}: {panel_present}")

    score, _ = panel_mean_score(expr_z, PANEL)

    # use rai_scores meta to split
    meta = pd.read_csv("/data/thca/repo_results/aggressive_sprint_2026_05_06/gse151179_rai_scores.tsv",
                       sep="\t", index_col=0)
    # 'collection_before_after_rai'
    col = "collection_before_after_rai"
    after_ids = meta.index[meta[col].astype(str).str.lower().str.contains("after")].tolist()
    before_ids = meta.index[meta[col].astype(str).str.lower().str.contains("before")].tolist()
    a = score.loc[score.index.intersection(after_ids)].values
    b = score.loc[score.index.intersection(before_ids)].values
    d, p, na, nb = cohens_d(a, b)
    out = [{
        "cohort": "GSE151179",
        "stratum": "RAI_axis",
        "contrast": "post-RAI_vs_pre-RAI",
        "n_high": na, "n_low": nb,
        "mean_d_HT13": d, "p": p,
        "sign": "+" if d > 0.1 else ("-" if d < -0.1 else "0"),
    }]
    pg = []
    for r in per_gene_cohens_d(expr_z, after_ids, before_ids, PANEL):
        r["cohort"] = "GSE151179"; r["stratum"] = "all"; r["contrast"] = "post-RAI_vs_pre-RAI"
        pg.append(r)
    return out, pg


# ----------------- main -----------------

def main():
    all_meta = []
    all_pg = []

    for fn in [run_tcga, run_gse286332, run_gse76039, run_mun2025, run_gse151179]:
        try:
            m, pg = fn()
            all_meta.extend(m)
            all_pg.extend(pg)
        except Exception as e:
            print(f"  !! {fn.__name__} failed: {type(e).__name__}: {e}")
            import traceback; traceback.print_exc()

    meta_df = pd.DataFrame(all_meta)
    pg_df = pd.DataFrame(all_pg)

    meta_df.to_csv(os.path.join(OUT, "h4_meta_table.tsv"), sep="\t", index=False)
    pg_df.to_csv(os.path.join(OUT, "h4_per_gene_d.tsv"), sep="\t", index=False)

    print("\n========== META TABLE ==========")
    print(meta_df.to_string(index=False))

    # sign tally on the primary contrast per cohort
    primary = meta_df[meta_df["stratum"].isin(
        ["all_primary_tumor", "Korean_PTC_HT_overlap", "Landa_dediff", "ATC_proteomics", "RAI_axis"])]
    pos = (primary["mean_d_HT13"] > 0.1).sum()
    neg = (primary["mean_d_HT13"] < -0.1).sum()
    null = (primary["mean_d_HT13"].between(-0.1, 0.1)).sum()
    n = len(primary)
    print(f"\nPRIMARY-CONTRAST SIGN TALLY (HT panel up in DM1 / HT-overlap / aggressive):")
    print(f"  positive (+): {pos}/{n}")
    print(f"  negative (-): {neg}/{n}")
    print(f"  null    (0): {null}/{n}")

    summary = {
        "panel": PANEL,
        "n_cohorts_primary": n,
        "n_positive": int(pos),
        "n_negative": int(neg),
        "n_null": int(null),
    }
    with open(os.path.join(OUT, "h4_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
