#!/usr/bin/env python3
"""
Track 2 / Step 03 — Concordance of K2 allele frequencies with:
  (A) Kim 2014 Korean reference panel (n=413 unrelated, 6 loci)
  (B) AFND South Korea pool (5 alleles only — supplementary)
  (C) AFND East Asian aggregate (Japan + China + Taiwan via 04_hla_long.csv)
  (D) Lee 2024 GSE213647 normal-only RNA-seq HLA atlas (matched RNA-seq method)

All comparisons are descriptive: report Pearson r, Spearman rho, total
variation distance, KL divergence (with smoothing), and a per-allele table.
NO disease association tests.
"""
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

T_DIR = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "hla_deepdive_2026_05_08/track2_k2_baseline/tables"
)
KIM_FRQ = Path(
    "/home/seungho/personal/THCA_data_analysis/project/external_refs/"
    "kim_korean_hla_ref/KOR_REF/Kim_KOR_HLA.FRQ.frq"
)
AFND_JSON = Path(
    "/home/seungho/personal/THCA_data_analysis/project/manuscript_p2_brief/"
    "lit_enrich_2026_05_02/data/afnd_alleles.json"
)
AFND_LONG = Path(
    "/home/seungho/personal/THCA_data_analysis/project/manuscript_p2_brief/"
    "lit_enrich_2026_05_02/deep_analysis/04_hla_long.csv"
)
LEE_TSV = Path(
    "/data/thca/_repo_offload/arcasHLA_GSE213647/GSE213647_arcasHLA_genotypes.tsv"
)
LEE_META = Path(
    "/data/thca/v17_korean/GSE213647/sample_sheet_clinical.tsv"
)


def parse_kim_panel():
    """Read Kim 2014 .frq file; map HLA_<LOCUS>_<DIGITS> -> standard 'L*XX:YY' form.

    Keep only 4-digit (2-field) alleles (length-4 digits after locus). Skip
    serological 2-digit summaries (length-2 digits).
    """
    rows = []
    with open(KIM_FRQ) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 6 or not parts[1].startswith("HLA_"):
                continue
            tag = parts[1]  # e.g., HLA_A_0201
            try:
                maf = float(parts[4])
                n_chrobs = int(parts[5])
            except ValueError:
                continue
            m = re.match(r"HLA_([A-Z0-9]+)_(\d+)$", tag)
            if not m:
                continue
            locus, digits = m.group(1), m.group(2)
            if len(digits) != 4:
                continue
            allele = f"{locus}*{digits[:2]}:{digits[2:]}"
            rows.append({
                "locus": locus,
                "allele": allele,
                "kim2014_freq": maf,
                "kim2014_n_chrom": n_chrobs,
            })
    return pd.DataFrame(rows)


def tvd(p, q):
    """Total variation distance for two prob vectors (same length)."""
    return 0.5 * np.sum(np.abs(np.array(p) - np.array(q)))


def kl_div(p, q, eps=1e-6):
    p = np.clip(np.array(p), eps, None)
    q = np.clip(np.array(q), eps, None)
    p = p / p.sum()
    q = q / q.sum()
    return float(np.sum(p * np.log(p / q)))


def pop_aggregate_from_afnd_long(target_pops):
    """Build per-allele weighted frequency for selected populations from AFND long
    CSV. Weights = study sample size (n).
    """
    df = pd.read_csv(AFND_LONG)
    df = df[df["country"].isin(target_pops)].copy()
    if df.empty:
        return pd.DataFrame(columns=["allele", "weighted_freq", "total_n"])
    df["wf"] = df["frequency"] * df["n"]
    out = df.groupby("allele").apply(
        lambda g: pd.Series({
            "weighted_freq": g["wf"].sum() / g["n"].sum(),
            "total_n": g["n"].sum(),
            "n_studies": len(g),
        }), include_groups=False
    ).reset_index()
    return out


def get_lee_normal_genotypes():
    """Load Lee 2024 GSE213647 arcasHLA genotypes restricted to normal-tissue samples.

    Map SRR -> GSM via PRJNA882018 ENA run table, then GSM -> tissue_type via
    sample_sheet_clinical.tsv ('Normal' rows).
    """
    geno = pd.read_csv(LEE_TSV, sep="\t")
    run_table = T_DIR / "T_aux_GSE213647_run_table.tsv"
    if run_table.exists() and LEE_META.exists():
        rmap = pd.read_csv(run_table, sep="\t")
        meta = pd.read_csv(LEE_META, sep="\t")
        normal_gsm = set(meta.loc[meta["tissue_type"].astype(str).str.lower()
                                  == "normal", "gsm"].astype(str))
        normal_srr = set(rmap.loc[rmap["sample_alias"].astype(str).isin(normal_gsm),
                                  "run_accession"].astype(str))
        n_pre = len(geno)
        geno = geno[geno["run"].astype(str).isin(normal_srr)].copy()
        print(f"[INFO] Lee normal mask via tissue_type=='Normal': "
              f"{len(geno)}/{n_pre} runs (n_normal_GSM={len(normal_gsm)})")
        geno.attrs["lee_subset"] = "normal_only"
    else:
        geno.attrs["lee_subset"] = "all_samples_(no_meta_or_run_table)"
        print("[WARN] missing meta or run table; using ALL Lee samples")
    return geno


def lee_allele_freqs(geno):
    """Compute allele freqs per locus from Lee genotype TSV (already 4-digit form).
    """
    locus_cols = {
        "A": ("A_a1_4d", "A_a2_4d"),
        "B": ("B_a1_4d", "B_a2_4d"),
        "C": ("C_a1_4d", "C_a2_4d"),
        "DRB1": ("DRB1_a1_4d", "DRB1_a2_4d"),
        "DPA1": ("DPA1_a1_4d", "DPA1_a2_4d"),
        "DQA1": ("DQA1_a1_4d", "DQA1_a2_4d"),
        "DQB1": ("DQB1_a1_4d", "DQB1_a2_4d"),
        "DPB1": ("DPB1_a1_4d", "DPB1_a2_4d"),
    }
    rows = []
    for locus, (c1, c2) in locus_cols.items():
        if c1 not in geno.columns:
            continue
        a = geno[[c1, c2]].dropna(how="all")
        chrom = pd.concat([a[c1], a[c2]]).dropna()
        n_chrom = len(chrom)
        if n_chrom == 0:
            continue
        vc = chrom.value_counts()
        for allele, k in vc.items():
            rows.append({
                "locus": locus,
                "allele": allele,
                "lee_freq": k / n_chrom,
                "lee_n_chrom": n_chrom,
                "lee_count": int(k),
            })
    return pd.DataFrame(rows)


def main():
    af = pd.read_csv(T_DIR / "T03_K2_allele_freq_per_locus.tsv", sep="\t")
    af_keep = af[["locus", "allele", "allele_freq", "af_ci_lo", "af_ci_hi",
                  "allele_count", "n_chrom"]].rename(
        columns={"allele_freq": "k2_freq", "n_chrom": "k2_n_chrom",
                 "allele_count": "k2_count"}
    )

    # ---------- (A) Kim 2014 Korean reference ----------
    kim = parse_kim_panel()
    kim_loci = set(kim["locus"])
    print(f"[INFO] Kim 2014 panel loci: {sorted(kim_loci)}; alleles n={len(kim)}")

    merged_kim = af_keep.merge(kim, on=["locus", "allele"], how="outer")
    merged_kim["k2_freq"] = merged_kim["k2_freq"].fillna(0)
    merged_kim["kim2014_freq"] = merged_kim["kim2014_freq"].fillna(0)
    merged_kim.to_csv(T_DIR / "T05_K2_vs_Kim2014.tsv", sep="\t", index=False)

    # Per-locus concordance metrics on alleles >=1% in either
    rows = []
    for locus in sorted(merged_kim["locus"].dropna().unique()):
        sub = merged_kim[merged_kim["locus"] == locus]
        sub = sub[(sub["k2_freq"] >= 0.01) | (sub["kim2014_freq"] >= 0.01)]
        if len(sub) < 3:
            continue
        x = sub["k2_freq"].values
        y = sub["kim2014_freq"].values
        r, p_pearson = stats.pearsonr(x, y)
        rho, p_spear = stats.spearmanr(x, y)
        rows.append({
            "locus": locus,
            "n_alleles_compared": int(len(sub)),
            "pearson_r": float(r),
            "pearson_p": float(p_pearson),
            "spearman_rho": float(rho),
            "spearman_p": float(p_spear),
            "tvd": float(tvd(x / x.sum() if x.sum() else x,
                             y / y.sum() if y.sum() else y)),
            "kl_k2_to_kim": float(kl_div(x, y)),
            "kl_kim_to_k2": float(kl_div(y, x)),
        })
    pd.DataFrame(rows).to_csv(
        T_DIR / "T06_K2_vs_Kim2014_concordance_metrics.tsv",
        sep="\t", index=False)
    print("Kim 2014 concordance:")
    print(pd.DataFrame(rows).to_string(index=False))

    # ---------- (B) AFND South Korea pool (5 alleles only) ----------
    afnd_kor = pop_aggregate_from_afnd_long(["Korea"])
    afnd_kor = afnd_kor.rename(columns={"weighted_freq": "afnd_kor_freq",
                                        "total_n": "afnd_kor_n",
                                        "n_studies": "afnd_kor_studies"})
    afnd_kor["locus"] = afnd_kor["allele"].str.split("*").str[0]
    afnd_kor_merge = af_keep.merge(afnd_kor, on=["locus", "allele"], how="inner")
    afnd_kor_merge.to_csv(T_DIR / "T07_K2_vs_AFND_Korea_5allele.tsv",
                          sep="\t", index=False)

    # ---------- (C) AFND East Asian aggregate (Japan + China + Taiwan) ----------
    afnd_ea = pop_aggregate_from_afnd_long(["Japan", "China", "Taiwan"])
    afnd_ea = afnd_ea.rename(columns={"weighted_freq": "afnd_ea_freq",
                                      "total_n": "afnd_ea_n",
                                      "n_studies": "afnd_ea_studies"})
    afnd_ea["locus"] = afnd_ea["allele"].str.split("*").str[0]
    afnd_ea_merge = af_keep.merge(afnd_ea, on=["locus", "allele"], how="inner")
    afnd_ea_merge.to_csv(T_DIR / "T08_K2_vs_AFND_EastAsian_5allele.tsv",
                         sep="\t", index=False)

    # ---------- (D) Lee 2024 GSE213647 normal-only ----------
    lee_geno = get_lee_normal_genotypes()
    lee_subset_label = lee_geno.attrs.get("lee_subset", "unknown")
    lee_af = lee_allele_freqs(lee_geno)
    print(f"[INFO] Lee subset = {lee_subset_label}; n_samples={lee_geno['run'].nunique()}; "
          f"n_alleles={len(lee_af)}")
    merged_lee = af_keep.merge(lee_af, on=["locus", "allele"], how="outer")
    merged_lee["k2_freq"] = merged_lee["k2_freq"].fillna(0)
    merged_lee["lee_freq"] = merged_lee["lee_freq"].fillna(0)
    merged_lee["lee_subset_used"] = lee_subset_label
    merged_lee.to_csv(T_DIR / "T09_K2_vs_Lee2024.tsv", sep="\t", index=False)

    rows_l = []
    for locus in sorted(merged_lee["locus"].dropna().unique()):
        sub = merged_lee[merged_lee["locus"] == locus]
        sub = sub[(sub["k2_freq"] >= 0.01) | (sub["lee_freq"] >= 0.01)]
        if len(sub) < 3:
            continue
        x = sub["k2_freq"].values
        y = sub["lee_freq"].values
        r, p_pearson = stats.pearsonr(x, y)
        rho, p_spear = stats.spearmanr(x, y)
        rows_l.append({
            "locus": locus,
            "n_alleles_compared": int(len(sub)),
            "pearson_r": float(r),
            "pearson_p": float(p_pearson),
            "spearman_rho": float(rho),
            "spearman_p": float(p_spear),
            "tvd": float(tvd(x / x.sum() if x.sum() else x,
                             y / y.sum() if y.sum() else y)),
            "kl_k2_to_lee": float(kl_div(x, y)),
            "kl_lee_to_k2": float(kl_div(y, x)),
        })
    pd.DataFrame(rows_l).to_csv(
        T_DIR / "T10_K2_vs_Lee2024_concordance_metrics.tsv",
        sep="\t", index=False)
    print("Lee 2024 concordance:")
    print(pd.DataFrame(rows_l).to_string(index=False))


if __name__ == "__main__":
    main()
