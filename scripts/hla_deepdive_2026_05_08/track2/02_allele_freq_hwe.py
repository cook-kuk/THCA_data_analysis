#!/usr/bin/env python3
"""
Track 2 / Step 02 — Allele frequencies (95% Clopper-Pearson CI) + Hardy-Weinberg
exact test per locus on the K2 PRJEB11591 cohort. Descriptive only.
"""
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

T_DIR = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "hla_deepdive_2026_05_08/track2_k2_baseline/tables"
)

LOCI = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


def clopper_pearson(k: int, n: int, alpha: float = 0.05):
    """Two-sided 95% Clopper-Pearson interval for a binomial proportion."""
    if n == 0:
        return (np.nan, np.nan)
    if k == 0:
        lo = 0.0
    else:
        lo = stats.beta.ppf(alpha / 2, k, n - k + 1)
    if k == n:
        hi = 1.0
    else:
        hi = stats.beta.ppf(1 - alpha / 2, k + 1, n - k)
    return float(lo), float(hi)


def hwe_collapse_chisq(genotypes_pairs):
    """HWE chi-square test for a multi-allelic locus.

    Strategy: collapse alleles with expected genotype count < 5 into a single
    'OTHER' allele, then do a chi-square comparing observed vs HW expected
    genotype counts on the collapsed table. dof = (g-1) - (k-1) where g =
    distinct collapsed genotypes used and k = collapsed alleles.

    Also reports the inbreeding coefficient F = 1 - H_obs / H_exp (Wright's F)
    which is robust to allele rarity, and a permutation p-value for F by
    randomly pairing chromosomes.
    """
    from collections import Counter
    pairs = [(a, b) for a, b in genotypes_pairs if a and b]
    n_geno = len(pairs)
    if n_geno == 0:
        return {"chi2": np.nan, "dof": np.nan, "p_chisq": np.nan,
                "F_obs": np.nan, "p_perm_F": np.nan,
                "n_genotypes": 0, "n_alleles": 0,
                "H_obs": np.nan, "H_exp": np.nan}

    # Allele counts
    ac = Counter()
    for a, b in pairs:
        ac[a] += 1
        ac[b] += 1
    n_chrom = sum(ac.values())
    n_alleles_total = len(ac)

    # Observed heterozygosity
    h_obs = sum(1 for a, b in pairs if a != b) / n_geno
    p_init = {a: ac[a] / n_chrom for a in ac}
    h_exp = 1 - sum(v ** 2 for v in p_init.values())
    F_obs = 1 - h_obs / h_exp if h_exp > 0 else np.nan

    # Identify rare alleles for collapse: q < 0.05 (so 2*p*q*N < 5 for most pairs)
    rare_thresh = 0.05
    rare = {a for a, q in p_init.items() if q < rare_thresh}

    def collapse(a):
        return "OTHER" if a in rare else a

    pairs_c = [(collapse(a), collapse(b)) for a, b in pairs]
    ac_c = Counter()
    for a, b in pairs_c:
        ac_c[a] += 1
        ac_c[b] += 1
    n_chr = sum(ac_c.values())
    p = {a: ac_c[a] / n_chr for a in ac_c}
    alleles = sorted(ac_c)
    k = len(alleles)

    gc = Counter()
    for a, b in pairs_c:
        gc[tuple(sorted([a, b]))] += 1

    obs, exp = [], []
    for i, ai in enumerate(alleles):
        for aj in alleles[i:]:
            if ai == aj:
                e = n_geno * p[ai] ** 2
            else:
                e = n_geno * 2 * p[ai] * p[aj]
            o = gc.get(tuple(sorted([ai, aj])), 0)
            obs.append(o)
            exp.append(e)
    obs = np.array(obs, float)
    exp = np.array(exp, float)
    # drop zero-expected
    keep = exp > 0
    obs = obs[keep]
    exp = exp[keep]
    chi2 = ((obs - exp) ** 2 / exp).sum() if exp.sum() else np.nan
    n_geno_classes = len(obs)
    # standard HW test dof = (n_genotype_classes) - (k_alleles) = k(k-1)/2 in textbook
    dof = max(1, n_geno_classes - k)
    p_chisq = 1 - stats.chi2.cdf(chi2, dof) if not np.isnan(chi2) else np.nan

    # Permutation p-value for F (random chromosome pairing under HWE).
    rng = np.random.default_rng(13)
    n_perm = 2000
    chrom = []
    for a, b in pairs:
        chrom.extend([a, b])
    chrom = np.array(chrom)
    F_null = np.empty(n_perm)
    for i in range(n_perm):
        idx = rng.permutation(len(chrom))
        c = chrom[idx]
        h_perm = np.mean(c[0::2] != c[1::2])
        F_null[i] = 1 - h_perm / h_exp if h_exp > 0 else np.nan
    p_perm = float(np.mean(np.abs(F_null) >= np.abs(F_obs)))

    return {"chi2": float(chi2), "dof": int(dof), "p_chisq": float(p_chisq),
            "F_obs": float(F_obs), "p_perm_F": p_perm,
            "n_genotypes": n_geno, "n_alleles": int(n_alleles_total),
            "n_alleles_collapsed": int(k),
            "H_obs": float(h_obs), "H_exp": float(h_exp)}


def main() -> None:
    geno = pd.read_csv(T_DIR / "T01_K2_genotypes_2field.tsv", sep="\t")
    n_total = len(geno)
    print(f"[INFO] N = {n_total} K2 samples")

    freq_rows = []
    hwe_rows = []
    pairs_per_locus = {}

    for locus in LOCI:
        a1 = geno[f"{locus}_a1"].dropna().tolist()
        a2 = geno[f"{locus}_a2"].dropna().tolist()
        # Carrier-level chr counts (each sample contributes up to 2 chromosomes)
        all_alleles = a1 + a2
        n_chrom = len(all_alleles)
        # Sample with both alleles present (typed sample)
        both_mask = geno[f"{locus}_a1"].notna() & geno[f"{locus}_a2"].notna()
        n_typed = int(both_mask.sum())
        n_chrom_typed = 2 * n_typed

        from collections import Counter
        ac = Counter(all_alleles)

        # Carrier counts (number of samples carrying ≥1 copy)
        carrier = Counter()
        for _, r in geno[both_mask].iterrows():
            s = {r[f"{locus}_a1"], r[f"{locus}_a2"]}
            for a in s:
                if a:
                    carrier[a] += 1

        for allele, k in ac.items():
            lo, hi = clopper_pearson(k, n_chrom_typed)
            car = carrier.get(allele, 0)
            car_lo, car_hi = clopper_pearson(car, n_typed)
            freq_rows.append({
                "locus": locus,
                "allele": allele,
                "allele_count": k,
                "n_chrom": n_chrom_typed,
                "allele_freq": k / n_chrom_typed if n_chrom_typed else np.nan,
                "af_ci_lo": lo,
                "af_ci_hi": hi,
                "carriers": car,
                "n_samples_typed": n_typed,
                "carrier_freq": car / n_typed if n_typed else np.nan,
                "carrier_ci_lo": car_lo,
                "carrier_ci_hi": car_hi,
            })

        pairs = [(r[f"{locus}_a1"], r[f"{locus}_a2"]) for _, r in geno[both_mask].iterrows()]
        pairs_per_locus[locus] = pairs
        hwe = hwe_collapse_chisq(pairs)
        hwe["locus"] = locus
        hwe_rows.append(hwe)

    af = pd.DataFrame(freq_rows).sort_values(
        ["locus", "allele_freq"], ascending=[True, False]
    ).reset_index(drop=True)
    af.to_csv(T_DIR / "T03_K2_allele_freq_per_locus.tsv", sep="\t", index=False)
    print(f"[OK] wrote T03_K2_allele_freq_per_locus.tsv ({len(af)} alleles)")

    hwe = pd.DataFrame(hwe_rows)[
        ["locus", "n_genotypes", "n_alleles", "n_alleles_collapsed",
         "H_obs", "H_exp", "F_obs",
         "chi2", "dof", "p_chisq", "p_perm_F"]
    ]
    hwe["hwe_pass_chisq_p_gt_0p05"] = hwe["p_chisq"].apply(
        lambda x: bool(x > 0.05) if pd.notna(x) else False
    )
    hwe["hwe_pass_perm_p_gt_0p05"] = hwe["p_perm_F"].apply(
        lambda x: bool(x > 0.05) if pd.notna(x) else False
    )
    hwe.to_csv(T_DIR / "T04_K2_HWE_per_locus.tsv", sep="\t", index=False)
    print("HWE summary:")
    print(hwe.to_string(index=False))


if __name__ == "__main__":
    main()
