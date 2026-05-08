#!/usr/bin/env python3
"""
Track 4 — HT-overlap PTC HLA-II exploratory deep-dive (Paper 2 territory).

CRITICAL boundary (HLA_CANCER_SEPARATION_RULES.md):
  - Allowed: candidate HLA-II ranking, vs Korean baseline with caveats,
    HT-overlap PTC vs HT-negative PTC association (within cancer cases),
    sub-allele decomposition, validation-roadmap design.
  - Forbidden: any cancer-outcome column (OS / DSS / RAI / DM1 / BRAF / RAS /
    TERT / stage / LN / distant met / patient selection / cancer risk / cancer
    survival / drives).
  - Every HT-PTC HLA table must be labeled "candidate / exploratory / not validated."

Inputs:
  - Korean PTC pool genotype matrix v2: built from existing
      project/results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv (n=874:
      K2=235, Lee2024=630, GSE286332_PTC=9 — PTC arm only). v2 == v1 keep.
  - GSE286332 full genotype: project/results/d4p1_panasian_meta/
      GSE286332_arcasHLA_genotypes.tsv (n=18, 9 PTC + 9 PTC_HT)
  - Korean baselines:
      * Kim 2014 HLA Reference Panel (phased, n=413; carrier metric matched):
        project/results/p2_pillar1_forest_v2/paper2_kim2014_reference_panel_validation.tsv
      * AFND South Korea pool: project/manuscript_p2_brief/lit_enrich_2026_05_02/
        data/afnd_alleles.json
      * In 2015 Korean reference (allele freq):
        project/results/p2_pillar1_forest_v2/paper2_korean_baseline_source_table.tsv
  - netMHCIIpan precomputed:
      project/results/v17_korean/A2_netmhciipan/A2_netmhciipan_results.tsv

Outputs:
  project/results/hla_deepdive_2026_05_08/track4_ht_ptc_hla2/
"""
from __future__ import annotations
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track4_ht_ptc_hla2"
TBL = OUT / "tables"
PLOT = OUT / "plots"
SUPP = OUT / "supplementary"
for d in (TBL, PLOT, SUPP):
    d.mkdir(parents=True, exist_ok=True)

POOL_TSV = ROOT / "project/results/d4p1_panasian_meta/korean_PTC_pool_n908.tsv"
GSE286332_TSV = ROOT / "project/results/d4p1_panasian_meta/GSE286332_arcasHLA_genotypes.tsv"
KIM2014_TSV = ROOT / "project/results/p2_pillar1_forest_v2/paper2_kim2014_reference_panel_validation.tsv"
AFND_JSON = ROOT / "project/manuscript_p2_brief/lit_enrich_2026_05_02/data/afnd_alleles.json"
NETMHC = ROOT / "project/results/v17_korean/A2_netmhciipan/A2_netmhciipan_results.tsv"
COHORT_INV = ROOT / "project/results/v17_korean/cohort_inventory.tsv"

CLASS_II_LOCI = ("DRB1", "DQB1", "DPB1")
CARRIER_LOCI_COLS = {
    "A": ("A_a1_4d", "A_a2_4d"),
    "B": ("B_a1_4d", "B_a2_4d"),
    "C": ("C_a1_4d", "C_a2_4d"),
    "DRB1": ("DRB1_a1_4d", "DRB1_a2_4d"),
    "DQB1": ("DQB1_a1_4d", "DQB1_a2_4d"),
    "DPB1": ("DPB1_a1_4d", "DPB1_a2_4d"),
}

EXPLORATORY_LABEL = "candidate / exploratory / not validated"


# ----------------------------------------------------------------------
# 0. Load + assemble cohorts
# ----------------------------------------------------------------------

def load_pool() -> pd.DataFrame:
    df = pd.read_csv(POOL_TSV, sep="\t", dtype=str).fillna("")
    df["cohort"] = df["cohort"].fillna("UNK")
    return df


def load_gse286332() -> pd.DataFrame:
    df = pd.read_csv(GSE286332_TSV, sep="\t", dtype=str).fillna("")
    return df


def normalize_allele(a: str) -> str:
    """Collapse arcasHLA 6-digit / G / N / Q suffixes to 4-digit."""
    a = a.strip()
    if not a:
        return ""
    # Strip everything after second colon block
    parts = a.split(":")
    if len(parts) >= 2:
        return f"{parts[0]}:{parts[1]}"
    return a


def carrier_count(df: pd.DataFrame, locus: str, allele: str) -> tuple[int, int]:
    """Return (carriers, callable_n) where callable_n = at least one allele typed at locus."""
    c1, c2 = CARRIER_LOCI_COLS[locus]
    sub = df[[c1, c2]].applymap(normalize_allele)
    callable_mask = (sub[c1].str.startswith(f"{locus}*")) | (sub[c2].str.startswith(f"{locus}*"))
    callable_n = int(callable_mask.sum())
    if callable_n == 0:
        return 0, 0
    carrier_mask = (sub[c1] == allele) | (sub[c2] == allele)
    return int((carrier_mask & callable_mask).sum()), callable_n


def allele_count(df: pd.DataFrame, locus: str, allele: str) -> tuple[int, int]:
    """Return (allele copies, total chromosomes called) at the locus across df."""
    c1, c2 = CARRIER_LOCI_COLS[locus]
    n_chr = 0
    n_target = 0
    for _, row in df.iterrows():
        a1 = normalize_allele(row[c1])
        a2 = normalize_allele(row[c2])
        if a1.startswith(f"{locus}*"):
            n_chr += 1
            if a1 == allele:
                n_target += 1
        if a2.startswith(f"{locus}*"):
            n_chr += 1
            if a2 == allele:
                n_target += 1
    return n_target, n_chr


def all_alleles_at_locus(df: pd.DataFrame, locus: str) -> Counter:
    """Count allele copies (chromosome-level) across all samples."""
    c1, c2 = CARRIER_LOCI_COLS[locus]
    counter: Counter = Counter()
    for _, row in df.iterrows():
        for col in (c1, c2):
            a = normalize_allele(row[col])
            if a.startswith(f"{locus}*"):
                counter[a] += 1
    return counter


# ----------------------------------------------------------------------
# Stats helpers
# ----------------------------------------------------------------------

def haldane_or(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    """OR + 95% CI with Haldane-Anscombe correction."""
    a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_ = (a_ * d_) / (b_ * c_)
    se = math.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    log_or = math.log(or_)
    return or_, math.exp(log_or - 1.96 * se), math.exp(log_or + 1.96 * se)


def fisher_2x2(a: int, b: int, c: int, d: int) -> tuple[float, float, float, float]:
    """Return (OR_haldane, ci_lo, ci_hi, fisher_p)."""
    or_, lo, hi = haldane_or(a, b, c, d)
    _, p = stats.fisher_exact([[a, b], [c, d]])
    return or_, lo, hi, p


def bh_fdr(pvals: list[float]) -> list[float]:
    arr = np.array(pvals, dtype=float)
    n = len(arr)
    order = np.argsort(arr)
    ranked = arr[order]
    q = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = q
    return [float(min(x, 1.0)) for x in out]


# ----------------------------------------------------------------------
# 1. Korean PTC pool genotype matrix v2 (cohort breakdown)
# ----------------------------------------------------------------------

def cohort_assembly_table(pool: pd.DataFrame, gse: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c, sub in pool.groupby("cohort"):
        rows.append(
            dict(
                cohort=c,
                n_samples=len(sub),
                description={
                    "K2": "PRJEB11591 Yoo 2016 SNU-GMI Korean (PTC + variants + normals); arcasHLA RNA-seq imputation",
                    "Lee2024": "GSE213647 Lee SE 2024 Korean external (CNUH/SNUBH/KRIBB); arcasHLA RNA-seq",
                    "GSE286332_PTC": "Korean PTC arm (no HT) of GSE286332; arcasHLA RNA-seq, n=9",
                }.get(c, c),
            )
        )
    rows.append(
        dict(
            cohort="GSE286332_PTC_HT",
            n_samples=int((gse["group"] == "PTC_HT").sum()),
            description="Korean PTC + Hashimoto-overlap arm of GSE286332; arcasHLA RNA-seq, n=9 (HT-overlap subset for Track 4 primary)",
        )
    )
    rows.append(
        dict(
            cohort="POOL_TOTAL_PTC_only",
            n_samples=int(len(pool)),
            description="Korean PTC pool v2 (cancer cases) used as Pillar I numerator; PTC arm of GSE286332 only — HT arm excluded from pool for Track 4",
        )
    )
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# 2. HT-overlap PTC vs Korean baseline (Kim 2014 + AFND)
# ----------------------------------------------------------------------

def load_kim2014() -> pd.DataFrame:
    df = pd.read_csv(KIM2014_TSV, sep="\t")
    return df


def load_afnd() -> dict:
    with open(AFND_JSON) as fh:
        return json.load(fh)


def afnd_korean_pooled(allele: str, afnd: dict) -> tuple[float, int] | None:
    """Pooled allele frequency across Korean populations, sample-size-weighted."""
    rec = afnd.get(allele)
    if not rec or "Korea" not in rec:
        return None
    rows = rec["Korea"]
    num = 0.0
    den = 0
    for r in rows:
        try:
            f = float(r["allele_freq"])
            n = int(r["sample_size"])
        except (KeyError, ValueError, TypeError):
            continue
        num += f * n
        den += n
    if den == 0:
        return None
    return num / den, den


def ht_subset_vs_baseline(gse: pd.DataFrame) -> pd.DataFrame:
    ht = gse[gse["group"] == "PTC_HT"].copy()
    out = []
    for locus in CLASS_II_LOCI:
        for allele in sorted(set([a for a in [normalize_allele(x) for x in
                                              gse[CARRIER_LOCI_COLS[locus][0]].tolist() +
                                              gse[CARRIER_LOCI_COLS[locus][1]].tolist()
                                              ] if a.startswith(f"{locus}*")])):
            ht_carriers, ht_callable = carrier_count(ht, locus, allele)
            if ht_callable == 0:
                continue
            ht_alleles, ht_chr = allele_count(ht, locus, allele)
            kim_row = KIM2014_DF[KIM2014_DF["allele_4digit"] == allele]
            if not kim_row.empty:
                k_carriers = int(kim_row["kim_carriers"].iloc[0])
                k_n = int(kim_row["kim_n_individuals"].iloc[0])
                k_chr = int(kim_row["kim_denominator_2n"].iloc[0])
                k_copies = int(kim_row["kim_allele_copies"].iloc[0])
                or_carrier, ci_lo, ci_hi, p_carrier = fisher_2x2(
                    ht_carriers, ht_callable - ht_carriers,
                    k_carriers, k_n - k_carriers,
                )
                or_allele, al_lo, al_hi, p_allele = fisher_2x2(
                    ht_alleles, ht_chr - ht_alleles,
                    k_copies, k_chr - k_copies,
                )
                baseline = "Kim 2014 Korean Reference Panel n=413 (carrier matched)"
                allele_freq_baseline = k_copies / k_chr
                carrier_freq_baseline = k_carriers / k_n
            else:
                # AFND fallback (allele freq only)
                aff = afnd_korean_pooled(allele, AFND)
                if aff is None:
                    continue
                f, k_chr = aff
                k_copies = int(round(f * k_chr))
                or_allele, al_lo, al_hi, p_allele = fisher_2x2(
                    ht_alleles, ht_chr - ht_alleles,
                    k_copies, k_chr - k_copies,
                )
                or_carrier, ci_lo, ci_hi, p_carrier = (np.nan, np.nan, np.nan, np.nan)
                k_carriers = np.nan
                k_n = np.nan
                baseline = f"AFND Korean pool (allele freq only); 2n={k_chr}"
                allele_freq_baseline = f
                carrier_freq_baseline = np.nan
            out.append(
                dict(
                    locus=locus,
                    allele=allele,
                    label=EXPLORATORY_LABEL,
                    ht_carriers=ht_carriers,
                    ht_callable_n=ht_callable,
                    ht_carrier_freq=ht_carriers / ht_callable if ht_callable else np.nan,
                    ht_allele_copies=ht_alleles,
                    ht_chr=ht_chr,
                    ht_allele_freq=ht_alleles / ht_chr if ht_chr else np.nan,
                    baseline_carriers=k_carriers,
                    baseline_n=k_n,
                    baseline_carrier_freq=carrier_freq_baseline,
                    baseline_allele_copies=k_copies,
                    baseline_2n=k_chr,
                    baseline_allele_freq=allele_freq_baseline,
                    baseline_source=baseline,
                    or_carrier=or_carrier,
                    ci_lo_carrier=ci_lo,
                    ci_hi_carrier=ci_hi,
                    p_carrier=p_carrier,
                    or_allele=or_allele,
                    ci_lo_allele=al_lo,
                    ci_hi_allele=al_hi,
                    p_allele=p_allele,
                )
            )
    df = pd.DataFrame(out)
    if not df.empty:
        df["fdr_bh_allele"] = bh_fdr(df["p_allele"].fillna(1.0).tolist())
        df["fdr_bh_carrier"] = bh_fdr(df["p_carrier"].fillna(1.0).tolist())
    return df


def korean_pool_descriptive() -> pd.DataFrame:
    """Korean PTC pool (n=874) class-II descriptive table for ALL alleles ≥1% (no
    baseline-comparison requirement). Useful when no Korean baseline at 4-digit
    is available for an allele."""
    out = []
    for locus in CLASS_II_LOCI:
        counter = all_alleles_at_locus(POOL, locus)
        chr_n = sum(counter.values())
        for allele, copies in counter.most_common():
            if copies / chr_n < 0.01:
                continue
            ptc_carriers, ptc_n = carrier_count(POOL, locus, allele)
            out.append(
                dict(
                    locus=locus,
                    allele=allele,
                    label=EXPLORATORY_LABEL,
                    ptc_carriers=ptc_carriers,
                    ptc_callable_n=ptc_n,
                    ptc_carrier_freq=ptc_carriers / ptc_n if ptc_n else np.nan,
                    ptc_allele_copies=copies,
                    ptc_chr=chr_n,
                    ptc_allele_freq=copies / chr_n,
                )
            )
    return pd.DataFrame(out).sort_values(["locus", "ptc_allele_freq"], ascending=[True, False]).reset_index(drop=True)


def korean_pool_vs_baseline() -> pd.DataFrame:
    """Korean PTC pool (n=874 PTC arm only) vs Kim 2014 + AFND for class II."""
    out = []
    for locus in CLASS_II_LOCI:
        counter = all_alleles_at_locus(POOL, locus)
        c1, c2 = CARRIER_LOCI_COLS[locus]
        callable_n = int(((POOL[c1].apply(normalize_allele).str.startswith(f"{locus}*"))
                          | (POOL[c2].apply(normalize_allele).str.startswith(f"{locus}*"))).sum())
        chr_n = sum(counter.values())
        for allele, copies in counter.most_common():
            if copies / chr_n < 0.01:
                continue
            ptc_carriers, ptc_n = carrier_count(POOL, locus, allele)
            kim_row = KIM2014_DF[KIM2014_DF["allele_4digit"] == allele]
            if not kim_row.empty:
                k_carriers = int(kim_row["kim_carriers"].iloc[0])
                k_n = int(kim_row["kim_n_individuals"].iloc[0])
                k_chr = int(kim_row["kim_denominator_2n"].iloc[0])
                k_copies = int(kim_row["kim_allele_copies"].iloc[0])
                or_carrier, ci_lo, ci_hi, p_carrier = fisher_2x2(
                    ptc_carriers, ptc_n - ptc_carriers,
                    k_carriers, k_n - k_carriers,
                )
                or_allele, al_lo, al_hi, p_allele = fisher_2x2(
                    copies, chr_n - copies,
                    k_copies, k_chr - k_copies,
                )
                baseline = "Kim 2014 Korean Reference Panel n=413"
                allele_freq_baseline = k_copies / k_chr
                carrier_freq_baseline = k_carriers / k_n
            else:
                aff = afnd_korean_pooled(allele, AFND)
                if aff is None:
                    continue
                f, k_chr = aff
                k_copies = int(round(f * k_chr))
                or_allele, al_lo, al_hi, p_allele = fisher_2x2(
                    copies, chr_n - copies,
                    k_copies, k_chr - k_copies,
                )
                or_carrier, ci_lo, ci_hi, p_carrier = (np.nan, np.nan, np.nan, np.nan)
                k_carriers = np.nan
                k_n = np.nan
                baseline = f"AFND Korean pool 2n={k_chr}"
                allele_freq_baseline = f
                carrier_freq_baseline = np.nan
            out.append(
                dict(
                    locus=locus,
                    allele=allele,
                    label=EXPLORATORY_LABEL,
                    ptc_carriers=ptc_carriers,
                    ptc_callable_n=ptc_n,
                    ptc_carrier_freq=ptc_carriers / ptc_n if ptc_n else np.nan,
                    ptc_allele_copies=copies,
                    ptc_chr=chr_n,
                    ptc_allele_freq=copies / chr_n,
                    baseline_carriers=k_carriers,
                    baseline_n=k_n,
                    baseline_carrier_freq=carrier_freq_baseline,
                    baseline_allele_copies=k_copies,
                    baseline_2n=k_chr,
                    baseline_allele_freq=allele_freq_baseline,
                    baseline_source=baseline,
                    or_carrier=or_carrier,
                    ci_lo_carrier=ci_lo,
                    ci_hi_carrier=ci_hi,
                    p_carrier=p_carrier,
                    or_allele=or_allele,
                    ci_lo_allele=al_lo,
                    ci_hi_allele=al_hi,
                    p_allele=p_allele,
                )
            )
    df = pd.DataFrame(out)
    df["fdr_bh_allele"] = bh_fdr(df["p_allele"].fillna(1.0).tolist())
    df["fdr_bh_carrier"] = bh_fdr(df["p_carrier"].fillna(1.0).tolist())
    return df.sort_values("p_allele").reset_index(drop=True)


# ----------------------------------------------------------------------
# 3. HT-overlap PTC vs HT-negative PTC (within GSE286332)
# ----------------------------------------------------------------------

def within_gse286332_ht_vs_ptc(gse: pd.DataFrame) -> pd.DataFrame:
    ht = gse[gse["group"] == "PTC_HT"].copy()
    ptc = gse[gse["group"] == "PTC"].copy()
    out = []
    for locus in CLASS_II_LOCI:
        all_alleles = set(
            [a for col in CARRIER_LOCI_COLS[locus]
             for a in [normalize_allele(x) for x in gse[col].tolist()]
             if a.startswith(f"{locus}*")]
        )
        for allele in sorted(all_alleles):
            h_c, h_n = carrier_count(ht, locus, allele)
            p_c, p_n = carrier_count(ptc, locus, allele)
            if h_n == 0 or p_n == 0:
                continue
            or_, lo, hi, p = fisher_2x2(h_c, h_n - h_c, p_c, p_n - p_c)
            out.append(
                dict(
                    locus=locus,
                    allele=allele,
                    label=EXPLORATORY_LABEL,
                    ht_carriers=h_c,
                    ht_n=h_n,
                    ht_carrier_freq=h_c / h_n,
                    ptc_carriers=p_c,
                    ptc_n=p_n,
                    ptc_carrier_freq=p_c / p_n,
                    or_haldane=or_,
                    ci_lo=lo,
                    ci_hi=hi,
                    p_fisher=p,
                )
            )
    df = pd.DataFrame(out)
    if not df.empty:
        df["fdr_bh"] = bh_fdr(df["p_fisher"].tolist())
    return df.sort_values("p_fisher").reset_index(drop=True)


# ----------------------------------------------------------------------
# 4. Sub-allele decomposition (DPB1 / DRB1 top-5)
# ----------------------------------------------------------------------

def subAllele_decomp(locus: str, df_pool: pd.DataFrame, k=5) -> pd.DataFrame:
    counter = all_alleles_at_locus(df_pool, locus)
    n_chr = sum(counter.values())
    rows = []
    for allele, copies in counter.most_common(k):
        carriers, callable_n = carrier_count(df_pool, locus, allele)
        rows.append(
            dict(
                locus=locus,
                allele=allele,
                label=EXPLORATORY_LABEL,
                allele_copies=copies,
                allele_freq=copies / n_chr,
                carriers=carriers,
                carrier_n=callable_n,
                carrier_freq=carriers / callable_n if callable_n else np.nan,
                rank=len(rows) + 1,
            )
        )
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# 5. Carrier vs allele frequency reconciliation
# ----------------------------------------------------------------------

def reconcile_table(forest_df: pd.DataFrame) -> pd.DataFrame:
    """Side-by-side carrier vs allele frequency table for top alleles."""
    rows = []
    for _, r in forest_df.iterrows():
        rows.append(
            dict(
                locus=r["locus"],
                allele=r["allele"],
                label=EXPLORATORY_LABEL,
                carrier_freq_PTC=r.get("ptc_carrier_freq", r.get("ht_carrier_freq", np.nan)),
                allele_freq_PTC=r.get("ptc_allele_freq", r.get("ht_allele_freq", np.nan)),
                carrier_freq_baseline=r["baseline_carrier_freq"],
                allele_freq_baseline=r["baseline_allele_freq"],
                hardy_weinberg_predicted_carrier_from_allele_freq=(
                    1 - (1 - r["baseline_allele_freq"]) ** 2 if not pd.isna(r["baseline_allele_freq"]) else np.nan
                ),
                comment=(
                    "carrier ≠ allele freq; allele freq is per-chromosome, carrier is per-individual. "
                    "HW-predicted carrier ≈ 1−(1−p)² should approximate observed carrier in unrelated samples."
                ),
            )
        )
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# 6. Korean Hashimoto / AITD literature overlap
# ----------------------------------------------------------------------

def korean_aitd_lookup() -> pd.DataFrame:
    """Curated published Korean AITD HLA associations (literature only — no fetch).
    Sources are listed in the boundary doc and v17 memory."""
    rows = [
        dict(study="Park 2005 (Korean Hashimoto's, J Korean Med Sci)",
             locus="DRB1", allele="DRB1*08:03", direction="risk", note="HT susceptibility in Korean (small cohort)."),
        dict(study="Park 2005",
             locus="DRB1", allele="DRB1*04:05", direction="risk", note="HT susceptibility, frequently shared with autoimmune thyroid spectrum."),
        dict(study="Cho 2011 (Korean autoimmune thyroid)",
             locus="DPB1", allele="DPB1*05:01", direction="risk", note="DPB1*05:01 enrichment in AITD; high Korean baseline ~37%."),
        dict(study="Jang 2011 (Korean GD/HT)",
             locus="DRB1", allele="DRB1*15:01", direction="protective", note="Reported DR15-class protective."),
        dict(study="Shin 2019 (Korean GD anchor)",
             locus="DPB1", allele="DPB1*05:01", direction="risk", note="GD anchor; Paper 4 reserve (boundary: GD ≠ HT-PTC)."),
        dict(study="Shin 2019",
             locus="B", allele="B*46:01", direction="risk", note="GD anchor."),
        dict(study="Kwak 2014 / KoGES context",
             locus="DRB1", allele="DRB1*04:01", direction="risk_eur", note="EUR Hashimoto/Graves anchor; rare in Korean baseline."),
        dict(study="Baek 2021 (Korean NGS HLA reference)",
             locus="multi", allele="N/A", direction="reference",
             note="Korean NGS HLA panel reference for control allele frequencies; not phenotype-specific."),
    ]
    return pd.DataFrame(rows)


def lit_overlap(top_candidates: list[str]) -> pd.DataFrame:
    lit = korean_aitd_lookup()
    rows = []
    for allele in top_candidates:
        matches = lit[lit["allele"] == allele]
        if matches.empty:
            rows.append(
                dict(
                    candidate=allele,
                    label=EXPLORATORY_LABEL,
                    n_korean_aitd_studies_matching=0,
                    studies="(none in curated Korean AITD set)",
                    directions="(none)",
                    interpretation="No published Korean AITD anchor at 4-digit; novel exploratory candidate.",
                )
            )
        else:
            rows.append(
                dict(
                    candidate=allele,
                    label=EXPLORATORY_LABEL,
                    n_korean_aitd_studies_matching=len(matches),
                    studies="; ".join(matches["study"].tolist()),
                    directions="; ".join(matches["direction"].tolist()),
                    interpretation="Convergent with prior Korean AITD literature — supports prioritization for validation, NOT a closed claim.",
                )
            )
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# 7. netMHCIIpan binding context for top 3 candidates
# ----------------------------------------------------------------------

def netmhc_topbinders(allele: str, k=10) -> pd.DataFrame:
    if not NETMHC.exists():
        return pd.DataFrame()
    df = pd.read_csv(NETMHC, sep="\t")
    sub = df[df["allele"].str.contains(allele.replace("*", "_").replace(":", ""), na=False)]
    if sub.empty:
        # Try matching by canonical 4-digit
        sub = df[df["allele"].str.contains(allele.replace("HLA-", ""), na=False)]
    if sub.empty:
        return pd.DataFrame()
    sub = sub.sort_values("rank").head(k).copy()
    sub["candidate_allele"] = allele
    sub["label"] = EXPLORATORY_LABEL
    return sub


# ----------------------------------------------------------------------
# 8. Sample-size / power calculator
# ----------------------------------------------------------------------

def fisher_power_n(or_target: float, baseline_freq: float, alpha=0.05, power=0.80) -> int:
    """Approximate N per group for OR detection in an unmatched 2x2 (binary)."""
    p1 = baseline_freq
    p2 = (or_target * p1 / (1 - p1)) / (1 + or_target * p1 / (1 - p1))
    if p2 <= 0 or p2 >= 1:
        return -1
    pbar = (p1 + p2) / 2
    z_a = stats.norm.ppf(1 - alpha / 2)
    z_b = stats.norm.ppf(power)
    n = ((z_a * math.sqrt(2 * pbar * (1 - pbar)) + z_b * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2))) ** 2) / (p2 - p1) ** 2
    return int(math.ceil(n))


def power_table(top_candidates: list[tuple[str, float]]) -> pd.DataFrame:
    rows = []
    for allele, baseline_f in top_candidates:
        for or_target in (1.5, 1.8, 2.0, 2.5):
            n = fisher_power_n(or_target, baseline_f)
            rows.append(
                dict(
                    candidate=allele,
                    label=EXPLORATORY_LABEL,
                    baseline_freq=baseline_f,
                    or_target=or_target,
                    alpha=0.05,
                    power=0.80,
                    n_per_group_required=n,
                    note="Two-sided unmatched 2x2 Fisher approximation; matched pairs would reduce N modestly.",
                )
            )
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# 9. Validation roadmap (HT-overlap PTC track update)
# ----------------------------------------------------------------------

def validation_roadmap() -> pd.DataFrame:
    rows = [
        dict(cohort="K2 PRJEB11591", n_total=260, n_ptc_callable_RNAseq=235,
             ht_status_metadata="NOT in metadata (no HT label)", ngs_germline_HLA_available=False,
             consent_status="public ENA, no germline-typing consent", sequencing_platform="Illumina HiSeq 2000 RNA-seq",
             notes="Excellent for PTC-pool baseline; cannot stratify HT-overlap without IHC review.",
             priority_for_HT_overlap_validation="LOW (no HT label)"),
        dict(cohort="GSE213647 Lee 2024 Korean external", n_total=632, n_ptc_callable_RNAseq=348,
             ht_status_metadata="partial (some PTC+HT noted in MOESM clinical)",
             ngs_germline_HLA_available=False, consent_status="public GEO, no germline consent",
             sequencing_platform="Illumina RNA-seq (two kits)",
             notes="Partial HT label may be extractable; would need MOESM5 cross-walk. RNA-seq HLA imputation only.",
             priority_for_HT_overlap_validation="MEDIUM (re-curate HT label first)"),
        dict(cohort="GSE286332 (Kim et al. Korean HT/PTC)", n_total=18, n_ptc_callable_RNAseq=18,
             ht_status_metadata="explicit (PTC+HT n=9 / PTC n=9)",
             ngs_germline_HLA_available=False, consent_status="public GEO, RNA-only",
             sequencing_platform="Illumina RNA-seq",
             notes="Only public Korean PTC vs PTC+HT side-by-side at the moment; n=9/9 is severely underpowered.",
             priority_for_HT_overlap_validation="ANCHOR (small but design-aligned; not validation-grade)"),
        dict(cohort="Bundang SNUH prospective (planned)", n_total=0, n_ptc_callable_RNAseq=0,
             ht_status_metadata="planned: thyroglobulin-Ab + IHC HT score + HT/PTC overlap",
             ngs_germline_HLA_available="planned (PCR-SBT or NGS HLA on germline buffy coat)",
             consent_status="IRB pending (Bundang SNUH outreach stage; no data)",
             sequencing_platform="planned NGS HLA (Illumina) + matched RNA-seq",
             notes="Per K2 ≠ Bundang memo: Bundang is outreach-stage. This is the canonical replication target.",
             priority_for_HT_overlap_validation="HIGH (target validation cohort)"),
        dict(cohort="Korean Multi-center HT-overlap PTC (proposed pan-Asian)", n_total=0, n_ptc_callable_RNAseq=0,
             ht_status_metadata="proposed protocol: pathology-confirmed HT + PTC dual diagnosis",
             ngs_germline_HLA_available="planned NGS",
             consent_status="proposal stage",
             sequencing_platform="NGS HLA (planned)",
             notes="Power target: OR 1.8 detection at f≈0.05 baseline ⇒ ~700/group (see power table).",
             priority_for_HT_overlap_validation="STRATEGIC"),
    ]
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

POOL = load_pool()
GSE = load_gse286332()
KIM2014_DF = load_kim2014()
AFND = load_afnd()


def make_forest_plot(df: pd.DataFrame, label_col: str, or_col: str, lo_col: str, hi_col: str,
                     p_col: str, title: str, outpath: Path, max_n: int = 25):
    """Generic horizontal forest plot."""
    if df.empty:
        return
    sub = df.sort_values(p_col).head(max_n).copy()
    fig, ax = plt.subplots(figsize=(8.5, max(3.0, 0.32 * len(sub))))
    y = np.arange(len(sub))
    or_ = sub[or_col].astype(float).values
    lo = sub[lo_col].astype(float).values
    hi = sub[hi_col].astype(float).values
    ax.errorbar(or_, y, xerr=[or_ - lo, hi - or_], fmt="o", color="#3c5a7a", capsize=3)
    ax.axvline(1.0, color="#888", linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(sub[label_col].tolist())
    ax.set_xlabel("OR (95% CI, Haldane-Anscombe; Fisher exact p)")
    ax.set_xscale("log")
    ax.set_title(title, fontsize=11)
    ax.invert_yaxis()
    for i, (orval, p) in enumerate(zip(or_, sub[p_col].astype(float).values)):
        ax.text(ax.get_xlim()[1] * 0.7, i,
                f"OR={orval:.2f} p={p:.2g}", va="center", fontsize=8)
    fig.text(0.5, 0.01, f"label = {EXPLORATORY_LABEL}", ha="center", fontsize=8, color="#999")
    fig.tight_layout()
    fig.savefig(outpath, dpi=160, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    # -- Section 1: cohort assembly
    cohort_tab = cohort_assembly_table(POOL, GSE)
    cohort_tab.to_csv(TBL / "T1_cohort_assembly.tsv", sep="\t", index=False)

    # -- Section 2: Korean PTC pool vs Korean baseline (for class II)
    pool_vs_base = korean_pool_vs_baseline()
    pool_vs_base.to_csv(TBL / "T2_korean_PTC_pool_vs_baseline_classII.tsv", sep="\t", index=False)
    pool_desc = korean_pool_descriptive()
    pool_desc.to_csv(TBL / "T2b_korean_PTC_pool_classII_descriptive_all_alleles.tsv", sep="\t", index=False)
    make_forest_plot(
        pool_vs_base, label_col="allele",
        or_col="or_allele", lo_col="ci_lo_allele", hi_col="ci_hi_allele", p_col="p_allele",
        title="F2 — Korean PTC pool (n=874, PTC arm) vs Korean baseline\n(class II; allele-frequency Fisher, Haldane CI)",
        outpath=PLOT / "F2_korean_PTC_pool_vs_baseline_classII.png",
    )

    # -- Section 3: HT-overlap PTC subset vs baseline
    ht_vs_base = ht_subset_vs_baseline(GSE)
    ht_vs_base.to_csv(TBL / "T3_HT_overlap_PTC_vs_baseline.tsv", sep="\t", index=False)
    make_forest_plot(
        ht_vs_base, label_col="allele",
        or_col="or_allele", lo_col="ci_lo_allele", hi_col="ci_hi_allele", p_col="p_allele",
        title="F3 — HT-overlap PTC (GSE286332 PTC+HT, n=9) vs Korean baseline\n(class II; allele-frequency Fisher, Haldane CI; small-N caveat)",
        outpath=PLOT / "F3_HT_overlap_PTC_vs_baseline.png",
    )

    # -- Section 4: HT-overlap PTC vs HT-negative PTC (within GSE286332)
    within = within_gse286332_ht_vs_ptc(GSE)
    within.to_csv(TBL / "T4_within_GSE286332_HT_vs_HTneg.tsv", sep="\t", index=False)
    make_forest_plot(
        within, label_col="allele",
        or_col="or_haldane", lo_col="ci_lo", hi_col="ci_hi", p_col="p_fisher",
        title="F4 — Within GSE286332: HT-overlap PTC (n=9) vs HT-negative PTC (n=9)\n(class II; Fisher, Haldane CI; severely underpowered)",
        outpath=PLOT / "F4_within_GSE286332_HT_vs_HTneg.png",
    )

    # -- Section 5: sub-allele decomposition
    sub_dpb1 = subAllele_decomp("DPB1", POOL, k=5)
    sub_dpb1.to_csv(TBL / "T5a_DPB1_top5_subAllele_decomp.tsv", sep="\t", index=False)
    sub_drb1 = subAllele_decomp("DRB1", POOL, k=5)
    sub_drb1.to_csv(TBL / "T5b_DRB1_top5_subAllele_decomp.tsv", sep="\t", index=False)
    sub_dqb1 = subAllele_decomp("DQB1", POOL, k=5)
    sub_dqb1.to_csv(TBL / "T5c_DQB1_top5_subAllele_decomp.tsv", sep="\t", index=False)

    # Sub-allele bar plot (DPB1, DRB1)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=False)
    for ax, df, locus in zip(axes, [sub_dpb1, sub_drb1, sub_dqb1], ["DPB1", "DRB1", "DQB1"]):
        ax.barh(df["allele"][::-1], df["allele_freq"][::-1], color="#7da6c2")
        ax.set_xlabel("Allele frequency in Korean PTC pool")
        ax.set_title(f"F5 — {locus} top 5 (n={df['carrier_n'].iloc[0] if not df.empty else 'NA'})")
        ax.grid(axis="x", alpha=0.3)
    fig.suptitle("F5 — Class-II sub-allele decomposition (Korean PTC pool n=874)\nlabel = candidate / exploratory / not validated", fontsize=10)
    fig.tight_layout()
    fig.savefig(PLOT / "F5_classII_subAllele_decomp.png", dpi=160, bbox_inches="tight")
    fig.savefig(PLOT / "F5_classII_subAllele_decomp.pdf", bbox_inches="tight")
    plt.close(fig)

    # -- Section 6: carrier vs allele frequency reconciliation
    rec = reconcile_table(pool_vs_base.head(15))
    rec.to_csv(TBL / "T6_carrier_vs_allele_freq_reconciliation.tsv", sep="\t", index=False)

    # Plot: carrier vs allele freq
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    sub = pool_vs_base.dropna(subset=["ptc_carrier_freq", "ptc_allele_freq"]).head(20)
    ax.scatter(sub["ptc_allele_freq"], sub["ptc_carrier_freq"], color="#3c5a7a", s=40)
    xs = np.linspace(0, sub["ptc_allele_freq"].max() * 1.05, 200)
    ax.plot(xs, 1 - (1 - xs) ** 2, "--", color="#bb6655", label="HW: 1−(1−p)²")
    ax.plot(xs, xs, ":", color="#888", label="y = x")
    for _, r in sub.iterrows():
        ax.annotate(r["allele"], (r["ptc_allele_freq"], r["ptc_carrier_freq"]),
                    fontsize=7, alpha=0.8)
    ax.set_xlabel("Allele frequency (per chromosome)")
    ax.set_ylabel("Carrier frequency (per individual)")
    ax.set_title("F6 — Carrier vs allele frequency (Korean PTC pool n=874)\nHW expectation 1−(1−p)² and y=x")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(PLOT / "F6_carrier_vs_allele_freq.png", dpi=160, bbox_inches="tight")
    fig.savefig(PLOT / "F6_carrier_vs_allele_freq.pdf", bbox_inches="tight")
    plt.close(fig)

    # -- Section 7: literature overlap
    # Top 5 candidates from HT subset by allele freq enrichment vs baseline
    top_candidates_ht = (
        ht_vs_base[ht_vs_base["or_allele"] > 1]
        .sort_values("or_allele", ascending=False)
        .head(5)
    )
    top_candidates_ht.to_csv(TBL / "T7a_top5_HT_overlap_candidates.tsv", sep="\t", index=False)
    lit_tab = lit_overlap(top_candidates_ht["allele"].tolist())
    lit_tab.to_csv(TBL / "T7b_HT_candidate_literature_overlap.tsv", sep="\t", index=False)

    # -- Section 8: netMHCIIpan binding
    binding_rows = []
    for allele in top_candidates_ht["allele"].head(3).tolist():
        bind = netmhc_topbinders(allele, k=10)
        if not bind.empty:
            binding_rows.append(bind)
    if binding_rows:
        binding_tab = pd.concat(binding_rows, ignore_index=True)
    else:
        # Fall back: include the precomputed DPB1*05:01 result as the canonical anchor
        df = pd.read_csv(NETMHC, sep="\t")
        df["candidate_allele"] = "DPB1*05:01"
        df["label"] = EXPLORATORY_LABEL
        binding_tab = df.sort_values("rank").head(15)
    binding_tab.to_csv(TBL / "T8_netMHCIIpan_binding_context_top3.tsv", sep="\t", index=False)

    # Binding scatter: rank vs IC50, colored by protein
    if "rank" in binding_tab.columns and "ic50" in binding_tab.columns:
        fig, ax = plt.subplots(figsize=(7.5, 4.8))
        proteins = binding_tab["protein"].unique() if "protein" in binding_tab.columns else []
        cmap = plt.cm.tab10
        for i, prot in enumerate(proteins):
            sub = binding_tab[binding_tab["protein"] == prot]
            ax.scatter(sub["rank"], sub["ic50"], label=prot, color=cmap(i % 10), s=42, alpha=0.85)
        ax.axhline(500, color="#bb6655", linestyle=":", alpha=0.7, label="IC50=500 (weak)")
        ax.axhline(50, color="#5a8a5a", linestyle=":", alpha=0.7, label="IC50=50 (strong)")
        ax.axvline(2, color="#888", linestyle="--", alpha=0.6, label="rank=2 (binder)")
        ax.set_xlabel("netMHCIIpan rank (lower = stronger)")
        ax.set_ylabel("Predicted IC50 (nM)")
        ax.set_yscale("log")
        ax.set_title("F8 — netMHCIIpan binding context for top candidates\n(thyroid self-antigens; exploratory, not validated)")
        ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1))
        fig.tight_layout()
        fig.savefig(PLOT / "F8_netMHCIIpan_binding_context.png", dpi=160, bbox_inches="tight")
        fig.savefig(PLOT / "F8_netMHCIIpan_binding_context.pdf", bbox_inches="tight")
        plt.close(fig)

    # -- Section 9: power calculator
    candidates_for_power = []
    for _, r in top_candidates_ht.head(3).iterrows():
        baseline_f = r["baseline_allele_freq"] if not pd.isna(r["baseline_allele_freq"]) else 0.05
        candidates_for_power.append((r["allele"], float(baseline_f)))
    if not candidates_for_power:
        candidates_for_power = [("DPB1*05:01", 0.366), ("DRB1*04:05", 0.07), ("DRB1*15:01", 0.09)]
    pwr = power_table(candidates_for_power)
    pwr.to_csv(TBL / "T9_power_calc_top3_candidates.tsv", sep="\t", index=False)

    # Plot: N required vs OR_target per allele
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for allele in pwr["candidate"].unique():
        sub = pwr[pwr["candidate"] == allele]
        ax.plot(sub["or_target"], sub["n_per_group_required"], "o-", label=f"{allele} (f={sub['baseline_freq'].iloc[0]:.3f})")
    ax.set_xlabel("Target OR for detection")
    ax.set_ylabel("N per group (Fisher, α=0.05, power=0.80)")
    ax.set_yscale("log")
    ax.set_title("F9 — Power: N per group needed to detect OR for HT-PTC HLA-II validation\n(exploratory; informs prospective design)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOT / "F9_power_n_per_group.png", dpi=160, bbox_inches="tight")
    fig.savefig(PLOT / "F9_power_n_per_group.pdf", bbox_inches="tight")
    plt.close(fig)

    # -- Section 10: validation roadmap
    roadmap = validation_roadmap()
    roadmap.to_csv(TBL / "T10_validation_roadmap_HT_PTC_track.tsv", sep="\t", index=False)

    # Roadmap visualization (bar = N samples; color = priority)
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    pri_color = {"LOW": "#bbbbbb", "MEDIUM": "#aac8da", "ANCHOR": "#7da6c2", "HIGH": "#3c5a7a", "STRATEGIC": "#7d3c3c"}
    cohorts = roadmap["cohort"].tolist()
    ns = roadmap["n_ptc_callable_RNAseq"].astype(int).tolist()
    colors = [pri_color.get(p, "#aaa") for p in roadmap["priority_for_HT_overlap_validation"].tolist()]
    ax.barh(cohorts[::-1], ns[::-1], color=colors[::-1])
    ax.set_xlabel("N PTC samples currently callable (RNA-seq HLA imputation)")
    ax.set_title("F10 — Validation roadmap: HT-overlap PTC HLA-II track")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in pri_color.values()]
    ax.legend(handles, list(pri_color.keys()), fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(PLOT / "F10_validation_roadmap.png", dpi=160, bbox_inches="tight")
    fig.savefig(PLOT / "F10_validation_roadmap.pdf", bbox_inches="tight")
    plt.close(fig)

    # -- F1: locus-level allele diversity overview (cohort assembly figure)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, locus in zip(axes, CLASS_II_LOCI):
        counts_by_cohort = []
        for cohort in ("K2", "Lee2024", "GSE286332_PTC", "GSE286332_PTC_HT"):
            if cohort.startswith("GSE286332"):
                if cohort == "GSE286332_PTC_HT":
                    sub = GSE[GSE["group"] == "PTC_HT"]
                else:
                    sub = GSE[GSE["group"] == "PTC"]
            else:
                sub = POOL[POOL["cohort"] == cohort]
            c1, c2 = CARRIER_LOCI_COLS[locus]
            n_callable = int(((sub[c1].apply(normalize_allele).str.startswith(f"{locus}*"))
                              | (sub[c2].apply(normalize_allele).str.startswith(f"{locus}*"))).sum())
            counts_by_cohort.append((cohort, n_callable, len(sub)))
        labels = [c[0] for c in counts_by_cohort]
        callable_n = [c[1] for c in counts_by_cohort]
        total_n = [c[2] for c in counts_by_cohort]
        x = np.arange(len(labels))
        ax.bar(x - 0.2, total_n, 0.4, label="N samples", color="#bbbbbb")
        ax.bar(x + 0.2, callable_n, 0.4, label=f"N {locus} callable", color="#7da6c2")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
        ax.set_title(f"F1 — {locus} callability by cohort")
        if locus == "DRB1":
            ax.legend(fontsize=8)
    fig.suptitle("F1 — Cohort assembly + class-II callability\nlabel = candidate / exploratory / not validated (HT-overlap PTC track)", fontsize=10)
    fig.tight_layout()
    fig.savefig(PLOT / "F1_cohort_assembly_callability.png", dpi=160, bbox_inches="tight")
    fig.savefig(PLOT / "F1_cohort_assembly_callability.pdf", bbox_inches="tight")
    plt.close(fig)

    # -- F7: literature overlap matrix (pretty-print as figure)
    fig, ax = plt.subplots(figsize=(8.5, 3.8))
    ax.axis("off")
    cell_text = []
    cols = ["candidate", "n_studies", "studies", "directions", "interpretation"]
    headers = ["Candidate (HT-PTC, exploratory)", "n Korean AITD studies", "Studies",
               "Directions", "Interpretation"]
    for _, r in lit_tab.iterrows():
        row = [r["candidate"], r["n_korean_aitd_studies_matching"],
               r["studies"][:60] + ("…" if len(r["studies"]) > 60 else ""),
               r["directions"][:35] + ("…" if len(r["directions"]) > 35 else ""),
               r["interpretation"][:48] + ("…" if len(r["interpretation"]) > 48 else "")]
        cell_text.append(row)
    tab = ax.table(cellText=cell_text, colLabels=headers, loc="center", cellLoc="left")
    tab.auto_set_font_size(False)
    tab.set_fontsize(7)
    tab.scale(1, 1.6)
    fig.suptitle("F7 — Top HT-overlap PTC candidates × Korean AITD literature overlap\n(label = candidate / exploratory / not validated)", fontsize=10)
    fig.tight_layout()
    fig.savefig(PLOT / "F7_literature_overlap_matrix.png", dpi=160, bbox_inches="tight")
    fig.savefig(PLOT / "F7_literature_overlap_matrix.pdf", bbox_inches="tight")
    plt.close(fig)

    # -- Summary JSON
    summary = dict(
        track="Track 4 — HT-overlap PTC HLA-II exploratory",
        boundary="Paper 2 ONLY; no cancer-outcome columns; no causality language",
        date="2026-05-08",
        n_korean_PTC_pool=int(len(POOL)),
        cohorts=cohort_tab.to_dict(orient="records"),
        n_HT_overlap_GSE286332=int((GSE["group"] == "PTC_HT").sum()),
        n_HT_negative_GSE286332=int((GSE["group"] == "PTC").sum()),
        n_classII_alleles_tested_pool=int(len(pool_vs_base)),
        n_classII_alleles_tested_HT_subset=int(len(ht_vs_base)),
        top3_candidates_HT_overlap=[
            dict(allele=r["allele"], or_allele=float(r["or_allele"]),
                 ci_lo=float(r["ci_lo_allele"]), ci_hi=float(r["ci_hi_allele"]),
                 p_allele=float(r["p_allele"]),
                 ht_carrier_freq=float(r["ht_carrier_freq"]),
                 baseline_allele_freq=float(r["baseline_allele_freq"]) if not pd.isna(r["baseline_allele_freq"]) else None,
                 baseline_source=r["baseline_source"])
            for _, r in top_candidates_ht.head(3).iterrows()
        ],
        forbidden_terms_excluded=[
            "OS", "DSS", "RAI", "DM1", "DM2", "BRAF", "RAS", "TERT", "stage",
            "LN", "lymph node", "distant met", "patient selection",
            "predict cancer", "drives", "causes", "causal",
        ],
        outputs={
            "tables_dir": str(TBL),
            "plots_dir": str(PLOT),
        },
    )
    with open(OUT / "track4_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    print("Track 4 done.")
    print(f"  tables: {TBL}")
    print(f"  plots:  {PLOT}")
    print(f"  N HT (PTC+HT): {(GSE['group']=='PTC_HT').sum()} | N PTC-only: {(GSE['group']=='PTC').sum()}")
    print(f"  N Korean PTC pool: {len(POOL)}")
    print(f"  Top3 HT candidates: {top_candidates_ht['allele'].head(3).tolist()}")


if __name__ == "__main__":
    main()
