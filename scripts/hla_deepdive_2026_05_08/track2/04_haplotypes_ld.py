#!/usr/bin/env python3
"""
Track 2 / Step 04 — Inferred haplotypes (counting approximation, no phasing)
and pairwise linkage disequilibrium (D', r^2) at the genotype level.

CAVEAT (stated in narrative): arcasHLA does NOT report phased haplotypes.
We use the "haplotype proxy" of (allele_a from one chromosome, allele_b from
the matching chromosome) where the assignment is arbitrary across loci. So
"haplotypes" here are effectively *allele co-occurrence pairs* averaged over
both chromosomal assignments. This is descriptive only.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

T_DIR = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "hla_deepdive_2026_05_08/track2_k2_baseline/tables"
)
LOCI = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


def carrier_set_per_sample(geno: pd.DataFrame, locus: str):
    """For each sample with both alleles called at locus, return unordered allele set."""
    out = {}
    for _, r in geno.iterrows():
        a1, a2 = r[f"{locus}_a1"], r[f"{locus}_a2"]
        if pd.notna(a1) and pd.notna(a2):
            out[r["run_accession"]] = (a1, a2)
    return out


def ld_pair(geno: pd.DataFrame, l1: str, l2: str, min_freq: float = 0.05):
    """Approximate D' and r^2 for two HLA loci by treating each locus's
    most-common allele as the "reference" (binary). For multi-allelic LD we
    expand to all (allele_i, allele_j) combinations with allele_i, allele_j
    both at >= min_freq, computing per-pair D, D', r^2 with EM-free counting:

        2N samples contribute four chromosomes pairs: (a1@l1, a1@l2),
        (a1@l1, a2@l2), (a2@l1, a1@l2), (a2@l1, a2@l2).

    Without phasing we assume random pairing — D estimate = chromosome-pair
    co-occurrence frequency minus product of marginals over those 4*N pseudo-
    haplotypes. This is the "genotype-pair counting approximation" from the
    spec; report it as such.
    """
    pairs1 = carrier_set_per_sample(geno, l1)
    pairs2 = carrier_set_per_sample(geno, l2)
    common = sorted(set(pairs1) & set(pairs2))
    n = len(common)
    if n == 0:
        return pd.DataFrame()

    # Per-locus marginal allele freqs (chrom-level)
    chr1 = []
    chr2 = []
    for s in common:
        chr1.extend(pairs1[s])
        chr2.extend(pairs2[s])
    n_chr = len(chr1)
    af1 = pd.Series(chr1).value_counts() / n_chr
    af2 = pd.Series(chr2).value_counts() / n_chr

    # Pseudo-haplotype frequencies via 4-way assignment (random pairing assumption)
    hap = Counter()
    for s in common:
        a1a, a1b = pairs1[s]
        a2a, a2b = pairs2[s]
        # All 4 pairings, each with weight 0.25 of the 2 chromosomes
        for x in (a1a, a1b):
            for y in (a2a, a2b):
                hap[(x, y)] += 1
    total = sum(hap.values())  # = 4*n  (== 2*n_chr)
    rows = []
    for a1, p1 in af1.items():
        if p1 < min_freq:
            continue
        for a2, p2 in af2.items():
            if p2 < min_freq:
                continue
            f12 = hap.get((a1, a2), 0) / total
            d = f12 - p1 * p2
            if d >= 0:
                d_max = min(p1 * (1 - p2), (1 - p1) * p2)
            else:
                d_max = min(p1 * p2, (1 - p1) * (1 - p2))
            d_prime = d / d_max if d_max > 0 else np.nan
            r2 = (d ** 2) / (p1 * (1 - p1) * p2 * (1 - p2)) if min(
                p1 * (1 - p1) * p2 * (1 - p2), 1) > 0 else np.nan
            rows.append({
                "locus_pair": f"{l1}-{l2}",
                "allele_l1": a1,
                "allele_l2": a2,
                "freq_l1": p1,
                "freq_l2": p2,
                "f12_obs": f12,
                "D": d,
                "D_prime": d_prime,
                "r2": r2,
                "n_samples": n,
            })
    return pd.DataFrame(rows)


def main():
    geno = pd.read_csv(T_DIR / "T01_K2_genotypes_2field.tsv", sep="\t")

    # ---------- Inferred haplotypes (top 10 per class) ----------
    # CLASS-I trio A-B-C, CLASS-II trio DRB1-DQB1-DPB1
    class1_loci = ["A", "B", "C"]
    class2_loci = ["DRB1", "DQB1", "DPB1"]

    def hap_top(locs, k=20):
        # Keep samples with all 3 loci typed
        mask = geno[[f"{l}_a1" for l in locs] + [f"{l}_a2" for l in locs]].notna().all(axis=1)
        sub = geno[mask]
        haps = Counter()
        for _, r in sub.iterrows():
            # 8 pseudo-haplotype combinations per sample
            for x in (r[f"{locs[0]}_a1"], r[f"{locs[0]}_a2"]):
                for y in (r[f"{locs[1]}_a1"], r[f"{locs[1]}_a2"]):
                    for z in (r[f"{locs[2]}_a1"], r[f"{locs[2]}_a2"]):
                        haps[(x, y, z)] += 1
        total = sum(haps.values())
        rows = []
        for hap, cnt in haps.most_common(k):
            rows.append({
                "haplotype": "~".join(hap),
                "count": cnt,
                "freq_pseudo": cnt / total,
                "n_samples_typed": int(mask.sum()),
            })
        return pd.DataFrame(rows), int(mask.sum())

    h1, n1 = hap_top(class1_loci, k=20)
    h2, n2 = hap_top(class2_loci, k=20)
    h1.to_csv(T_DIR / "T11_K2_classI_ABC_haplotypes_top20.tsv", sep="\t", index=False)
    h2.to_csv(T_DIR / "T12_K2_classII_DRB1_DQB1_DPB1_haplotypes_top20.tsv", sep="\t", index=False)
    print(f"Class-I top10 haplotypes (n_typed={n1}):")
    print(h1.head(10).to_string(index=False))
    print(f"Class-II top10 haplotypes (n_typed={n2}):")
    print(h2.head(10).to_string(index=False))

    # ---------- LD: class-I (A-B, B-C, A-C); class-II (DRB1-DQB1, DRB1-DPB1, DQB1-DPB1) ----------
    pairs_to_test = [("A", "B"), ("B", "C"), ("A", "C"),
                     ("DRB1", "DQB1"), ("DRB1", "DPB1"), ("DQB1", "DPB1")]
    all_ld = []
    summary = []
    for l1, l2 in pairs_to_test:
        df = ld_pair(geno, l1, l2, min_freq=0.05)
        all_ld.append(df)
        if not df.empty:
            summary.append({
                "locus_pair": f"{l1}-{l2}",
                "n_allele_pairs": len(df),
                "median_D_prime": float(df["D_prime"].median()),
                "max_D_prime": float(df["D_prime"].max()),
                "median_r2": float(df["r2"].median()),
                "max_r2": float(df["r2"].max()),
                "top_pair_by_r2": df.sort_values("r2", ascending=False).iloc[0][
                    ["allele_l1", "allele_l2", "r2", "D_prime", "f12_obs"]].to_dict(),
            })
    full_ld = pd.concat(all_ld, ignore_index=True)
    full_ld.to_csv(T_DIR / "T13_K2_LD_per_allele_pair.tsv", sep="\t", index=False)
    pd.DataFrame(summary).to_csv(T_DIR / "T14_K2_LD_pair_summary.tsv",
                                 sep="\t", index=False)
    print("LD summary:")
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
