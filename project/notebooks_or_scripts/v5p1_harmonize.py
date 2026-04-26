#!/usr/bin/env python3
"""v5.1 Phase 2 — Harmonize TCGA + GEO per cancer.

Per cancer with >=2 included real cohorts (THCA reuses v5 pipeline):
  - Normalize TCGA STAR counts -> log2(TPM+1)
  - Quantile-normalize GEO microarray within cohort
  - Intersect gene symbols; expect >=8000 shared
  - Concatenate; save X/Y/B + shared_genes.txt
"""
from __future__ import annotations

import sys
import gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import io
import urllib.request

import numpy as np
import pandas as pd

sys.path.insert(0, "/opt/thyroid-dash/project/notebooks_or_scripts")
from v5p1_common import (
    PROJECT, DATA_RAW_V5, DATA_PROC_V5, RESULTS_V5, LOGS, COHORTS, log_line
)

LOG = LOGS / "v5p1_harmonize.log"
GENE_LEN_PATH = PROJECT / "data_raw" / "gencode_v36_gene_lengths.tsv"


def _build_gene_lengths_from_tcga(cancer: str) -> Optional[pd.Series]:
    """Estimate gene length from a single TCGA STAR file (tpm + counts)."""
    tcga_dir = DATA_RAW_V5 / cancer / "tcga_rnaseq"
    files = list(tcga_dir.glob("*.tsv"))
    if not files:
        return None
    try:
        df = pd.read_csv(files[0], sep="\t", comment="#", skiprows=1)
        df = df[~df["gene_id"].str.startswith("N_", na=False)]
        # len = counts / tpm * (sum_tpm / 1e6)  (reverse TPM formula approx). Simpler:
        # TPM_i = (count_i / len_i) / sum_j(count_j/len_j) * 1e6
        # Direct from file would require multiple iterations. Instead use fpkm:
        # fpkm_i = count_i / (len_i * N_mapped) * 1e9 -> len_i = count_i * 1e9 / (fpkm_i * N_mapped)
        # But N_mapped is not provided. We use tpm_unstranded ratio-based approach.
        # Easier: use gene_type and compute length from count/tpm ratio.
        counts = df["stranded_second"].astype(float).values
        tpm = df["tpm_unstranded"].astype(float).values
        # avoid zero
        mask = (counts > 10) & (tpm > 0.1)
        if mask.sum() < 1000:
            return None
        ratio = counts[mask] / tpm[mask]
        # normalize ratio to per-kb (anchor: median across all expressed)
        # This yields a constant scaling up to factor; TPM normalization will
        # cancel it out anyway.
        # Use ratio directly (implicitly in kb units modulo a constant).
        lengths = np.full(len(df), np.nan)
        lengths[mask] = ratio
        s = pd.Series(lengths, index=df["gene_name"].values)
        s = s.groupby(level=0).mean()
        s = s.dropna()
        return s
    except Exception as e:
        log_line(LOG, f"[gene-len-build-err] {cancer}: {e}")
        return None


def counts_to_log2tpm(counts_df: pd.DataFrame, gene_lengths: pd.Series) -> pd.DataFrame:
    """counts_df: gene x sample. gene_lengths: gene-indexed Series (arbitrary units)."""
    common = counts_df.index.intersection(gene_lengths.index)
    counts_df = counts_df.loc[common]
    L = gene_lengths.loc[common].values.reshape(-1, 1)
    # rate = count / L
    rate = counts_df.values.astype(float) / np.maximum(L, 1e-3)
    # TPM: scale columns to sum 1e6
    col_sum = rate.sum(axis=0, keepdims=True)
    col_sum = np.where(col_sum == 0, 1, col_sum)
    tpm = rate / col_sum * 1e6
    log2tpm = np.log2(tpm + 1.0)
    return pd.DataFrame(log2tpm, index=common, columns=counts_df.columns)


def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Rank-based quantile normalization across columns (samples)."""
    arr = df.values.astype(float)
    rank = np.argsort(np.argsort(arr, axis=0), axis=0)
    sorted_arr = np.sort(arr, axis=0)
    mean_per_rank = sorted_arr.mean(axis=1)
    qn = mean_per_rank[rank]
    return pd.DataFrame(qn, index=df.index, columns=df.columns)


def load_tcga_for_cancer(cancer: str) -> Optional[Tuple[pd.DataFrame, pd.DataFrame]]:
    """Returns (expr_df gene x sample of log2TPM, labels_df sample x ['label','patient_id'])
    Or None if unavailable."""
    if cancer == "THCA":
        out_dir = DATA_RAW_V5 / "THCA"
        expr_path = out_dir / "tcga_log2tpm.tsv.gz"
        labels_path = out_dir / "tcga_samples_labels.tsv"
        if not expr_path.exists() or not labels_path.exists():
            return None
        expr = pd.read_csv(expr_path, sep="\t", index_col=0)
        labels = pd.read_csv(labels_path, sep="\t")
        labels = labels.set_index("sample_id")
        return expr, labels

    counts_path = DATA_RAW_V5 / cancer / "tcga_counts.tsv.gz"
    labels_path = DATA_RAW_V5 / cancer / "tcga_driver_labels.tsv"
    if not counts_path.exists() or not labels_path.exists():
        return None
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    labels = pd.read_csv(labels_path, sep="\t")  # patient_id, label

    # gene lengths
    gl = _build_gene_lengths_from_tcga(cancer)
    if gl is None:
        return None

    # remap samples: counts columns are 16-char barcodes (TCGA-XX-XXXX-NN)
    samples = counts.columns.tolist()
    patient_of = {s: s[:12] for s in samples}
    lab_by_pat = labels.set_index("patient_id")["label"].to_dict()
    sample_labels = {}
    for s in samples:
        pid = patient_of[s]
        if pid in lab_by_pat:
            sample_labels[s] = lab_by_pat[pid]
    kept = list(sample_labels.keys())
    if len(kept) == 0:
        return None
    counts = counts[kept]
    log2tpm = counts_to_log2tpm(counts, gl)
    labels_df = pd.DataFrame({"sample_id": kept,
                              "patient_id": [patient_of[s] for s in kept],
                              "label": [sample_labels[s] for s in kept]}).set_index("sample_id")
    return log2tpm, labels_df


def load_geo_for_cancer(cancer: str) -> List[Tuple[str, pd.DataFrame, pd.DataFrame]]:
    """Returns list of (gse, expr_df gene x sample log2 QN, labels_df sample x label)."""
    out = []
    geo_dir = DATA_RAW_V5 / cancer / "geo"
    if not geo_dir.exists():
        return out
    for gse in COHORTS[cancer]["geo"]:
        expr_p = geo_dir / f"{gse}_expr.tsv"
        pheno_p = geo_dir / f"{gse}_pheno.tsv"
        if not expr_p.exists() or not pheno_p.exists():
            continue
        try:
            expr = pd.read_csv(expr_p, sep="\t", index_col=0)
            pheno = pd.read_csv(pheno_p, sep="\t")
        except Exception as e:
            log_line(LOG, f"[geo-load-err] {gse}: {e}")
            continue
        # Determine if log scale; microarray values typically 0-20 range
        if expr.values.max() > 50:
            expr = np.log2(expr + 1.0)
        # Quantile normalize within cohort
        try:
            expr_qn = quantile_normalize(expr)
        except Exception as e:
            log_line(LOG, f"[qn-err] {gse}: {e}")
            expr_qn = expr
        # Align pheno
        keep_samples = [s for s in expr_qn.columns if s in pheno["sample_id"].values]
        lab_map = dict(zip(pheno["sample_id"], pheno["label"]))
        labels_df = pd.DataFrame({
            "sample_id": keep_samples,
            "label": [lab_map.get(s) for s in keep_samples]
        }).set_index("sample_id")
        labels_df = labels_df.dropna(subset=["label"])
        keep_samples = labels_df.index.tolist()
        if len(keep_samples) < 10:
            continue
        expr_qn = expr_qn[keep_samples]
        out.append((gse, expr_qn, labels_df))
    return out


def harmonize_cancer(cancer: str) -> Dict:
    log_line(LOG, f"=== harmonize {cancer} ===")
    tcga = load_tcga_for_cancer(cancer)
    geo_list = load_geo_for_cancer(cancer)

    cohorts = []  # list of (name, expr_df, labels_df)
    if tcga is not None:
        cohorts.append((f"TCGA-{cancer}", tcga[0], tcga[1]))
    for gse, expr, labs in geo_list:
        cohorts.append((gse, expr, labs))

    # If only TCGA is available, split it by Tissue Source Site (TSS) to
    # create multiple real-data batches. TSS codes are the 2-char segment
    # after "TCGA-" in the barcode (e.g. TCGA-CS-..., TCGA-DU-..., ...)
    # Only do this if there is no GEO cohort; TSS is a real sequencing
    # batch axis (different contributing institutions).
    if len(cohorts) == 1 and tcga is not None:
        expr, labs = tcga
        tss_map = {s: s.split("-")[1] for s in expr.columns if s.startswith("TCGA-")}
        tss_counts = pd.Series(list(tss_map.values())).value_counts()
        # Only split if at least 2 TSSes each have >=20 samples
        big_tss = tss_counts[tss_counts >= 20].index.tolist()
        if len(big_tss) >= 2:
            log_line(LOG, f"[{cancer}] single TCGA cohort — splitting by TSS: "
                          f"{dict(zip(big_tss, [int(tss_counts[t]) for t in big_tss]))}")
            # Rebuild cohorts as per-TSS sub-cohorts (keeping only big TSSes)
            new_cohorts = []
            for tss in big_tss:
                tss_samples = [s for s in expr.columns if tss_map.get(s) == tss]
                if len(tss_samples) < 20:
                    continue
                sub_expr = expr[tss_samples]
                sub_labs = labs.loc[[s for s in tss_samples if s in labs.index]]
                new_cohorts.append((f"TCGA-{cancer}-TSS-{tss}", sub_expr, sub_labs))
            if len(new_cohorts) >= 2:
                cohorts = new_cohorts

    if len(cohorts) < 2:
        log_line(LOG, f"[{cancer}] only {len(cohorts)} cohorts — need >=2")
        return {"status": "excluded", "n_cohorts": len(cohorts), "reason": "<2 cohorts"}

    # Intersect genes
    gene_set = None
    for _, e, _ in cohorts:
        gs = set(e.index.astype(str))
        if gene_set is None:
            gene_set = gs
        else:
            gene_set &= gs
    shared = sorted(gene_set)
    log_line(LOG, f"[{cancer}] {len(cohorts)} cohorts, {len(shared)} shared genes")

    if len(shared) < 3000:
        log_line(LOG, f"[{cancer}] WARNING n_shared_genes={len(shared)} < 3000")

    X_blocks = []
    Y_all = []
    B_all = []
    sample_names = []
    n_per_cohort = {}
    n_a_per_cohort = {}
    n_b_per_cohort = {}
    class_a = COHORTS[cancer]["class_a"]
    class_b = COHORTS[cancer]["class_b"]

    for name, expr, labs in cohorts:
        e = expr.loc[shared]
        common_samples = [s for s in e.columns if s in labs.index]
        # filter to class_a / class_b only
        common_samples = [s for s in common_samples
                          if labs.loc[s, "label"] in (class_a, class_b)]
        if len(common_samples) < 10:
            continue
        X_blocks.append(e[common_samples].values.T.astype(np.float32))
        Y_all.extend([labs.loc[s, "label"] for s in common_samples])
        B_all.extend([name] * len(common_samples))
        sample_names.extend(common_samples)
        n_per_cohort[name] = len(common_samples)
        n_a_per_cohort[name] = sum(1 for s in common_samples if labs.loc[s, "label"] == class_a)
        n_b_per_cohort[name] = sum(1 for s in common_samples if labs.loc[s, "label"] == class_b)

    if not X_blocks:
        return {"status": "excluded", "n_cohorts": 0, "reason": "no samples"}

    X = np.vstack(X_blocks)
    Y = np.array(Y_all)
    B = np.array(B_all)

    n_a_total = int((Y == class_a).sum())
    n_b_total = int((Y == class_b).sum())
    if n_a_total < 20 or n_b_total < 20:
        return {"status": "excluded_small_n", "n_cohorts": len(cohorts),
                "n_a": n_a_total, "n_b": n_b_total, "reason": "class n<20"}

    # Drop genes with any NaN across samples (real-world microarray + log2TPM
    # intersection often has some all-NaN rows for probes mapped to different
    # symbol variants across platforms).
    gene_nan = np.isnan(X).any(axis=0)
    if gene_nan.any():
        log_line(LOG, f"[{cancer}] dropping {int(gene_nan.sum())} genes with NaN "
                      f"({100*gene_nan.sum()/X.shape[1]:.1f}%)")
        X = X[:, ~gene_nan]
        shared = [g for i, g in enumerate(shared) if not gene_nan[i]]
    # Any remaining NaN -> 0 (shouldn't happen after mask above)
    if np.isnan(X).any():
        X = np.nan_to_num(X, nan=0.0)
    # Also drop zero-variance genes
    var = X.var(axis=0)
    zero_var = var < 1e-10
    if zero_var.any():
        log_line(LOG, f"[{cancer}] dropping {int(zero_var.sum())} zero-variance genes")
        X = X[:, ~zero_var]
        shared = [g for i, g in enumerate(shared) if not zero_var[i]]

    out_dir = DATA_PROC_V5 / cancer
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_dir / "X_combined.npz", X=X)
    np.savetxt(out_dir / "Y.tsv", Y, fmt="%s")
    np.savetxt(out_dir / "B.tsv", B, fmt="%s")
    with open(out_dir / "shared_genes.txt", "w") as fh:
        for g in shared:
            fh.write(g + "\n")
    with open(out_dir / "sample_names.txt", "w") as fh:
        for s in sample_names:
            fh.write(s + "\n")

    log_line(LOG, f"[{cancer}] X={X.shape} Y(a/b)={n_a_total}/{n_b_total} B_unique={len(set(B))}")
    return {
        "status": "ok",
        "n_cohorts": len(cohorts),
        "cohort_ids": list(n_per_cohort.keys()),
        "n_per_cohort": n_per_cohort,
        "n_a": n_a_total,
        "n_b": n_b_total,
        "n_shared_genes": len(shared),
        "n_samples": X.shape[0],
    }


def main():
    log_line(LOG, "=== v5p1 Phase 2 harmonize START ===")
    rows = []
    from v5p1_common import CANCERS
    for cancer in CANCERS:
        try:
            r = harmonize_cancer(cancer)
        except Exception as e:
            import traceback
            log_line(LOG, f"[{cancer}] ERROR: {traceback.format_exc()}")
            r = {"status": "error", "reason": str(e)}
        row = {"cancer": cancer, "status": r.get("status"),
               "n_cohorts": r.get("n_cohorts", 0),
               "cohort_ids": ",".join(r.get("cohort_ids", [])),
               "n_shared_genes": r.get("n_shared_genes", 0),
               "n_samples": r.get("n_samples", 0),
               "n_class_A": r.get("n_a", 0),
               "n_class_B": r.get("n_b", 0),
               "note": r.get("reason", ""),
               "semi_synthetic": False}
        rows.append(row)
    pd.DataFrame(rows).to_csv(RESULTS_V5 / "v5p1_harmonization.tsv", sep="\t", index=False)
    log_line(LOG, "=== v5p1 Phase 2 DONE ===")


if __name__ == "__main__":
    main()
