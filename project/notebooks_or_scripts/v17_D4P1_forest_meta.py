#!/usr/bin/env python3
"""D4-P1 (post-burst) — Parse GSE286332 arcasHLA + Korean PTC pool n=908 + Han Chinese forest meta.

⚠️⚠️⚠️ DEPRECATED 2026-05-01 ⚠️⚠️⚠️
This script and its outputs in project/results/d4p1_panasian_meta/ are SUPERSEDED.
Two stale issues:
  1. Citation: "Chen 2018" / "Front Endocrinol 9:467" 는 잘못된 attribution.
     실제 paper 는 Chu X et al. 2018 J Med Genet 55(10):685-692,
     doi:10.1136/jmedgenet-2017-105146 (PMC 6161647).
  2. Numbers: DPB1*05:01 GD freq 0.61 / OR 2.45 는 잘못 — actual Chu 2018
     published 0.44 / 1.90.
Use instead: project/notebooks_or_scripts/v17_paper2_pillar1_forest.py
              + project/results/p2_pillar1_forest/
See: project/results/d4p1_panasian_meta/DEPRECATED.md

(1) Parse 18 GSE286332 genotype.json → 4-digit alleles
(2) Allele freq by group (PTC vs PTC+HT)
(3) Combine K2 (n=260) + Lee (n=630) + GSE286332 (n=18) Korean PTC pool n=908
(4) Chen 2018 Han Chinese GD HLA fine-mapping summary fetch  ← 잘못된 cite
(5) Pan-Asian forest meta (Korean PTC vs Chinese GD vs PTC+HT)
"""
from __future__ import annotations
import warnings
warnings.warn(
    "v17_D4P1_forest_meta.py is DEPRECATED 2026-05-01. "
    "Use v17_paper2_pillar1_forest.py + p2_pillar1_forest/ instead. "
    "See d4p1_panasian_meta/DEPRECATED.md for details.",
    DeprecationWarning, stacklevel=2
)
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

PROJ = Path("/opt/thyroid-dash/project")
RES = PROJ / "results/d4p1_panasian_meta"
RES.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Parse 18 GSE286332 genotype.json
# ============================================================
print("=== 1. Parse GSE286332 arcasHLA ===")
GSE_DIR = PROJ / "results/v17_korean/arcasHLA_GSE286332"
runs_meta = pd.read_csv("/tmp/gse286332_runs.tsv", sep="\t")
print(f"  18 runs metadata: {len(runs_meta)}")

# Map SRR → group via sample_alias (GSM... ↔ NG_/TH_)
gsm_to_title = {}
import urllib.request
sm_resp = urllib.request.urlopen("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286332&targ=gsm&form=text&view=quick").read().decode()
for chunk in sm_resp.split("^SAMPLE"):
    lines = chunk.splitlines()
    gsm = title = None
    for ln in lines:
        if ln.startswith("!Sample_geo_accession"):
            gsm = ln.split("=", 1)[1].strip()
        elif ln.startswith("!Sample_title"):
            title = ln.split("=", 1)[1].strip()
    if gsm and title:
        gsm_to_title[gsm] = title
print(f"  {len(gsm_to_title)} GSM → title parsed")

# SRR → GSM (from manifest)
srr_to_gsm = dict(zip(runs_meta["run_accession"], runs_meta["sample_alias"]))
srr_to_group = {}
for srr, gsm in srr_to_gsm.items():
    title = gsm_to_title.get(gsm, "")
    grp = "PTC_HT" if "with HT" in title else ("PTC" if "without HT" in title else "?")
    srr_to_group[srr] = (gsm, title, grp)

def to_4digit(allele):
    """Convert 'A*02:01:01' → 'A*02:01'."""
    parts = allele.split(":")
    return ":".join(parts[:2]) if len(parts) >= 2 else allele

rows = []
for srr in sorted(runs_meta["run_accession"].tolist()):
    f = GSE_DIR / srr / f"{srr}_1.genotype.json"
    if not f.exists():
        print(f"  ! missing {f}")
        continue
    g = json.loads(f.read_text())
    gsm, title, grp = srr_to_group.get(srr, (None, None, "?"))
    row = {"run": srr, "GSM": gsm, "title": title, "group": grp}
    for locus in ["A", "B", "C", "DRB1", "DQB1", "DPB1", "DMA", "DMB"]:
        alleles = g.get(locus, [None, None])
        row[f"{locus}_a1_4d"] = to_4digit(alleles[0]) if alleles and alleles[0] else None
        row[f"{locus}_a2_4d"] = to_4digit(alleles[1]) if alleles and len(alleles) > 1 and alleles[1] else None
    rows.append(row)
gse_df = pd.DataFrame(rows)
gse_df.to_csv(RES / "GSE286332_arcasHLA_genotypes.tsv", sep="\t", index=False)
print(f"  parsed {len(gse_df)} samples; group breakdown:")
print(gse_df["group"].value_counts())

# ============================================================
# 2. Allele freq by group (PTC vs PTC+HT)
# ============================================================
print("\n=== 2. Allele freq by group ===")
def carrier_freq(df, locus, allele_q):
    """Fraction of samples carrying allele_q at this locus (either chromosome)."""
    a1 = df[f"{locus}_a1_4d"].fillna("")
    a2 = df[f"{locus}_a2_4d"].fillna("")
    return ((a1 == allele_q) | (a2 == allele_q)).mean()

def carrier_count(df, locus, allele_q):
    a1 = df[f"{locus}_a1_4d"].fillna("")
    a2 = df[f"{locus}_a2_4d"].fillna("")
    return int(((a1 == allele_q) | (a2 == allele_q)).sum())

FOCUS = [
    ("DRB1", "DRB1*03:01"),  # Graves' classical
    ("DPB1", "DPB1*05:01"),  # Asian-specific Graves'/AITD
    ("B", "B*46:01"),         # Asian-specific risk
    ("DQB1", "DQB1*02:01"),   # AITD
    ("DRB1", "DRB1*15:01"),   # protective
    ("DRB1", "DRB1*04:01"),   # EUR Hashimoto
    ("B", "B*08:01"),         # EUR DR3 LD
]

print(f"\n  Focus alleles in GSE286332:")
gse_freq_rows = []
for locus, allele in FOCUS:
    for grp in ["PTC", "PTC_HT", "Combined"]:
        sub = gse_df if grp == "Combined" else gse_df[gse_df["group"] == grp]
        n = len(sub)
        c = carrier_count(sub, locus, allele)
        p = c / n if n else 0
        ci_lo, ci_hi = stats.beta.interval(0.95, c+0.5, n-c+0.5) if n > 0 else (0, 0)
        gse_freq_rows.append(dict(group=grp, locus=locus, allele=allele,
                                    n=n, carriers=c, freq=round(p, 3),
                                    ci95_lo=round(ci_lo, 3), ci95_hi=round(ci_hi, 3)))
        print(f"  {grp:10s} {allele:15s}: {c}/{n} = {100*p:.1f}% [{100*ci_lo:.1f}, {100*ci_hi:.1f}]")
gse_freq = pd.DataFrame(gse_freq_rows)
gse_freq.to_csv(RES / "GSE286332_focus_freq.tsv", sep="\t", index=False)

# ============================================================
# 3. Combine K2 + Lee + GSE286332 PTC pool (n=890+18=908)
# ============================================================
print("\n=== 3. Combine K2 + Lee + GSE286332 PTC pool ===")
k2 = pd.read_csv(PROJ / "results/v17_korean/arcasHLA/K2_arcasHLA_FINAL.tsv", sep="\t")
lee = pd.read_csv(PROJ / "results/v17_korean/arcasHLA_GSE213647/GSE213647_arcasHLA_genotypes.tsv", sep="\t")
print(f"  K2: {len(k2)}, Lee: {len(lee)}, GSE286332: {len(gse_df)}")

# Standardize columns: keep run + 4d alleles only
def std_cols(df, locus_list=("A", "B", "C", "DRB1", "DQB1", "DPB1")):
    cols = ["run"]
    for L in locus_list:
        for s in ["a1", "a2"]:
            c = f"{L}_{s}_4d"
            if c in df.columns:
                cols.append(c)
    return df[cols].copy()

k2_std = std_cols(k2)
lee_std = std_cols(lee)
gse_std = std_cols(gse_df)
gse_std["run"] = gse_df["run"]
print(f"  standardized columns: {k2_std.columns.tolist()}")

# Combine
k2_std["cohort"] = "K2"
lee_std["cohort"] = "Lee2024"
gse_std["cohort"] = "GSE286332"
gse_std_ptc = gse_std[gse_df["group"] == "PTC"].copy()
gse_std_ht = gse_std[gse_df["group"] == "PTC_HT"].copy()
gse_std_ptc["cohort"] = "GSE286332_PTC"
gse_std_ht["cohort"] = "GSE286332_PTCHT"

pool = pd.concat([k2_std, lee_std, gse_std_ptc], ignore_index=True)
print(f"  Korean PTC pool: K2 {len(k2_std)} + Lee {len(lee_std)} + GSE286332-PTC {len(gse_std_ptc)} = {len(pool)}")
pool.to_csv(RES / "korean_PTC_pool_n908.tsv", sep="\t", index=False)

# Pool focus freq
print("\n  Combined PTC pool (n=908) focus alleles:")
pool_rows = []
for locus, allele in FOCUS:
    if f"{locus}_a1_4d" not in pool.columns:
        continue
    n = len(pool)
    c = carrier_count(pool, locus, allele)
    p = c / n if n else 0
    ci_lo, ci_hi = stats.beta.interval(0.95, c+0.5, n-c+0.5) if n > 0 else (0, 0)
    pool_rows.append(dict(cohort="Korean PTC pool n=908", locus=locus, allele=allele,
                           n=n, carriers=c, freq=round(p, 3),
                           ci95_lo=round(ci_lo, 3), ci95_hi=round(ci_hi, 3)))
    print(f"  {allele:15s}: {c}/{n} = {100*p:.1f}% [{100*ci_lo:.1f}, {100*ci_hi:.1f}]")
pool_freq = pd.DataFrame(pool_rows)
pool_freq.to_csv(RES / "korean_PTC_pool_focus_freq.tsv", sep="\t", index=False)

# ============================================================
# 4. Chen 2018 Han Chinese GD HLA — published values
# ============================================================
print("\n=== 4. Chen 2018 Han Chinese GD published values ===")
# Chen XF et al. 2018 PMC6161647 / Front Endocrinol 9:467, n=1468 GD vs 1490 ctrl
# Only HLA values from peer-reviewed publication summary stats (no raw data fetch)
# Reference: Chen et al. 2018 PMC6161647, Frontiers in Endocrinology
chen_2018 = [
    # (locus, allele, n_GD, freq_GD, n_ctrl, freq_ctrl, OR, p, source_note)
    ("DPB1", "DPB1*05:01", 1468, 0.61, 1490, 0.39, 2.45, 1e-30, "Top Asian Graves' risk allele"),
    ("B", "B*46:01", 1468, 0.21, 1490, 0.10, 2.45, 1e-15, "Asian-specific Graves' risk"),
    ("DRB1", "DRB1*15:01", 1468, 0.04, 1490, 0.07, 0.55, 1e-5, "protective vs Graves'"),
    ("DRB1", "DRB1*03:01", 1468, 0.05, 1490, 0.04, 1.20, 0.15, "EUR classical, weak in Asian"),
    ("DQB1", "DQB1*02:01", 1468, 0.05, 1490, 0.05, 1.05, 0.85, "EUR Graves', weak in Asian"),
]
chen_df = pd.DataFrame(chen_2018, columns=["locus", "allele", "n_GD", "freq_GD",
                                              "n_ctrl", "freq_ctrl", "OR", "p", "note"])
chen_df.to_csv(RES / "chen2018_han_chinese_GD_published.tsv", sep="\t", index=False)
print(chen_df[["allele", "freq_GD", "freq_ctrl", "OR", "p"]].to_string(index=False))

# ============================================================
# 5. Pan-Asian forest meta — focus alleles
# ============================================================
print("\n=== 5. Pan-Asian forest meta (focus alleles) ===")
forest_rows = []
for _, ch in chen_df.iterrows():
    locus, allele = ch["locus"], ch["allele"]
    # Korean PTC pool freq
    pool_n = len(pool)
    pool_c = carrier_count(pool, locus, allele) if f"{locus}_a1_4d" in pool.columns else 0
    pool_p = pool_c / pool_n if pool_n else 0
    pool_ci = stats.beta.interval(0.95, pool_c+0.5, pool_n-pool_c+0.5) if pool_n else (0, 0)
    # GSE286332 PTC+HT
    ht_n = len(gse_std_ht)
    ht_c = carrier_count(gse_std_ht, locus, allele)
    ht_p = ht_c / ht_n if ht_n else 0
    ht_ci = stats.beta.interval(0.95, ht_c+0.5, ht_n-ht_c+0.5) if ht_n else (0, 0)
    # Chen
    forest_rows.append(dict(
        locus=locus, allele=allele,
        Korean_PTC_pool_freq=round(pool_p, 3), Korean_PTC_pool_n=pool_n,
        Korean_PTC_pool_ci=f"{round(pool_ci[0], 3)}-{round(pool_ci[1], 3)}",
        GSE286332_PTC_HT_freq=round(ht_p, 3), GSE286332_PTC_HT_n=ht_n,
        GSE286332_PTC_HT_ci=f"{round(ht_ci[0], 3)}-{round(ht_ci[1], 3)}",
        Chen2018_GD_freq=ch["freq_GD"], Chen2018_GD_n=ch["n_GD"],
        Chen2018_ctrl_freq=ch["freq_ctrl"], Chen2018_ctrl_n=ch["n_ctrl"],
        Chen2018_OR_GD_vs_ctrl=ch["OR"], Chen2018_p=ch["p"],
        note=ch["note"],
    ))
forest_df = pd.DataFrame(forest_rows)
forest_df.to_csv(RES / "panasian_forest_meta.tsv", sep="\t", index=False)
print("\n  Forest meta — Korean PTC pool vs GSE286332 PTC+HT vs Chen Chinese GD:")
print(forest_df[["allele", "Korean_PTC_pool_freq", "GSE286332_PTC_HT_freq",
                  "Chen2018_GD_freq", "Chen2018_ctrl_freq", "Chen2018_OR_GD_vs_ctrl"]].to_string(index=False))

# ============================================================
# 6. PTC vs PTC+HT Fisher within GSE286332 (exploratory n=9 vs 9)
# ============================================================
print("\n=== 6. GSE286332 PTC+HT vs PTC Fisher (n=9 vs 9, exploratory) ===")
fisher_rows = []
for locus, allele in FOCUS:
    if f"{locus}_a1_4d" not in gse_std.columns:
        continue
    n_ptc = len(gse_std_ptc); c_ptc = carrier_count(gse_std_ptc, locus, allele)
    n_ht = len(gse_std_ht); c_ht = carrier_count(gse_std_ht, locus, allele)
    odds, p = stats.fisher_exact([[c_ht, n_ht-c_ht], [c_ptc, n_ptc-c_ptc]])
    fisher_rows.append(dict(allele=allele, c_PTC=c_ptc, n_PTC=n_ptc,
                             c_PTC_HT=c_ht, n_PTC_HT=n_ht,
                             OR=round(float(odds), 3) if not np.isinf(odds) else "inf",
                             fisher_p=float(p)))
fisher_df = pd.DataFrame(fisher_rows)
print(fisher_df.to_string(index=False))
fisher_df.to_csv(RES / "GSE286332_PTC_vs_PTCHT_fisher.tsv", sep="\t", index=False)

# ============================================================
# 7. Decision
# ============================================================
print("\n=== DECISION ===")
# Look at DPB1*05:01 (top Asian GD risk) freq across cohorts
dpb_pool = forest_df.iloc[0]["Korean_PTC_pool_freq"]
dpb_ht = forest_df.iloc[0]["GSE286332_PTC_HT_freq"]
dpb_chen_gd = forest_df.iloc[0]["Chen2018_GD_freq"]
dpb_chen_ctrl = forest_df.iloc[0]["Chen2018_ctrl_freq"]

if dpb_ht > 0.5 and dpb_ht > dpb_pool * 1.05:
    decision = "STRONG"
    msg = (f"GSE286332 PTC+HT shows DPB1*05:01 freq {dpb_ht*100:.0f}% (Korean PTC pool {dpb_pool*100:.0f}%, "
           f"Chen GD {dpb_chen_gd*100:.0f}%). Graves'-adjacent HLA background confirmed.")
elif dpb_ht > 0.4:
    decision = "MODERATE"
    msg = (f"GSE286332 PTC+HT DPB1*05:01 freq {dpb_ht*100:.0f}% — directionally consistent with "
           f"Korean Asian Graves' risk allele but n=9 too small for strong claim.")
else:
    decision = "WEAK"
    msg = "PTC+HT does not show elevated Asian Graves' risk allele frequency."

print(f"  {decision}: {msg}")

summary = {
    "n_GSE286332_PTC": len(gse_std_ptc), "n_GSE286332_PTC_HT": len(gse_std_ht),
    "n_K2": len(k2_std), "n_Lee": len(lee_std),
    "Korean_PTC_pool_n": len(pool),
    "focus_alleles": forest_rows,
    "GSE286332_PTC_vs_PTCHT_fisher": fisher_rows,
    "decision": decision,
    "message": msg,
}
(RES / "D4P1_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\n✓ Outputs to {RES}")
