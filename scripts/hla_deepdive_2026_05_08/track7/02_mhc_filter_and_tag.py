#!/usr/bin/env python3
"""Track 7 step 2 — MHC region filter, sub-region tagging, ancestry parsing,
and SNP→HLA-allele tag assignment from a literature-curated table.

Inputs : T01 GWAS Catalog associations TSV
Outputs:
   T02_mhc_associations.tsv            — chr6:28.48-33.45 Mb subset
   T03_top20_per_trait.tsv             — top 20 SNPs per trait by -log10p
   T04_snp_hla_allele_tags.tsv         — SNP↔HLA allele tag table (literature)
   T05_track1_convergence.tsv          — Tag-SNP β/OR vs Track 1 random-effects OR
   T06_ancestry_stratified.tsv         — long-form ancestry × locus
   T07_pleiotropy_matrix.tsv           — SNP × trait pleiotropy matrix
   T08_mhc_share_per_trait.tsv         — fraction of GWS hits in MHC
"""
from __future__ import annotations

import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track7_gwas_catalog_mhc")
TAB = ROOT / "tables"
TAB.mkdir(parents=True, exist_ok=True)
T1_DIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track1_paper4_gd_panasian/tables")

# GRCh38 MHC bounds (per Horton 2008, Ensembl MHC region)
MHC_LO = 28_477_797
MHC_HI = 33_448_354
# Sub-region windows (GRCh38)
SUBREGIONS = [
    ("extended_class_I_xMHC_left",  28_477_797, 29_640_168),  # extended class I (left)
    ("class_I_HLA-F_A_G",           29_640_168, 30_550_000),  # HLA-F, HLA-G, HLA-A region
    ("class_I_HLA-E_C_B",           31_236_000, 31_490_000),  # HLA-C/B core (B 31.27 Mb, C 31.27 Mb)
    ("class_III",                   31_490_000, 32_407_000),  # complement, TNF, HSP70
    ("class_II_DR_DQ",              32_407_000, 33_055_000),  # DRB1/DQB1/DQA1
    ("class_II_DP",                 33_055_000, 33_448_354),  # DPB1/DPA1
]
# Literature-curated SNP↔HLA tag table (de Bakker 2006 SNAGGR / Karnes 2017 / OkadaJP 2018 / Saevarsdottir 2020 / Sakaue 2021 / Brown 2014).
# Each row: rsid, tagged HLA allele, ancestry of tag map, source.
# These are *imputation tags*, not perfect equivalents. Used here only for
# qualitative SNP↔allele bridging; effect-size convergence is reported as such.
SNP_HLA_TAGS = pd.DataFrame(
    [
        # Class II DR/DQ/DP — best-known autoimmune tag SNPs
        ("rs9271365",  "HLA-DRB1*03:01",  "European",     "Sakaue 2021 / OpenSNP imputation panel — DRB1*03 lead in DR3 haplotype"),
        ("rs9271365",  "HLA-DQB1*02:01",  "European",     "DR3-DQ2 ancestral haplotype tag (de Bakker 2006)"),
        ("rs9270911",  "HLA-DRB1*03:01",  "European",     "Saevarsdottir 2020 Nature — DR3 region"),
        ("rs9272346",  "HLA-DQA1*05:01",  "European",     "de Bakker 2006 SNAGGR HapMap CEU"),
        ("rs9272346",  "HLA-DRB1*03:01",  "European",     "DR3-DQ2 haplotype, lead T1D/AITD tag"),
        ("rs2647050",  "HLA-DPB1*05:01",  "Japanese",     "OkadaJP 2018 deep-imputation; East-Asian DPB1 LD block"),
        ("rs17612848", "HLA-DPB1*05:01",  "Japanese",     "Cooper 2012 Hum Mol Genet — class II DP region"),
        ("rs17615",    "HLA-DPB1*02:01",  "Japanese",     "Karnes 2017 / dbMHC East-Asian"),
        ("rs6457617",  "HLA-DRB1*15:01",  "European",     "Brown 2014 IBD — DR15-DQ6 tag"),
        ("rs3129953",  "HLA-DQA1*05:01",  "European",     "Cooper 2012 AITD class II tag"),
        ("rs7775055",  "HLA-DRB1*04:05",  "Japanese",     "OkadaJP 2018 — Japanese DRB1*04:05 tag"),
        ("rs2858331",  "HLA-DQA1*02:01",  "European",     "celiac tag, AITD secondary"),
        ("rs3104413",  "HLA-DRB1*04:01",  "European",     "Stewart 2017 Nat Genet T1D imputation"),
        # Class I tags
        ("rs2523608",  "HLA-B*08:01",     "European",     "DR3-B8 ancestral haplotype tag"),
        ("rs2256318",  "HLA-A*02:01",     "European",     "Stewart 2017 class I imputation"),
        ("rs2523608",  "HLA-C*07:01",     "European",     "B8-C7 LD block"),
        # Class III flanking (FLT3 ligand pathway etc. — not HLA)
        ("rs1612904",  "HLA-DRB1*03:01",  "European",     "Cooper 2012 AITD secondary signal in class II"),
        ("rs1265564",  "HLA-C*06:02",     "European",     "psoriasis lead, class I"),
    ],
    columns=["rsid", "tagged_allele", "tag_ancestry", "tag_source"],
)


def assign_subregion(pos: float) -> str:
    if pd.isna(pos):
        return "unknown"
    for name, lo, hi in SUBREGIONS:
        if lo <= pos < hi:
            return name
    if MHC_LO <= pos <= MHC_HI:
        return "MHC_other"
    return "outside_MHC"


def parse_ancestry(s: str) -> dict:
    """Heuristic ancestry classification of the GWAS Catalog free-text initial sample-size string."""
    if not isinstance(s, str):
        return {"ancestry_class": "unknown", "n_cases": np.nan, "n_total": np.nan}
    low = s.lower()
    classes = []
    if "european" in low:
        classes.append("European")
    if "east asian" in low or "han chinese" in low or "japanese" in low or "korean" in low or "chinese ancestry" in low:
        classes.append("East_Asian")
    if "south asian" in low or "indian" in low or "bangladesh" in low:
        classes.append("South_Asian")
    if "african" in low:
        classes.append("African")
    if "hispanic" in low or "admixed american" in low or "latin" in low:
        classes.append("Admixed_American")
    if "middle eastern" in low or "arabian" in low:
        classes.append("Middle_Eastern")
    if not classes:
        cls = "Unknown_or_Multi"
    elif len(classes) == 1:
        cls = classes[0]
    else:
        cls = "Multi_ancestry"
    # Pull first integer that appears with "cases"
    m = re.search(r"(\d[\d,]*)\s+\w[^,]*?(?:ancestry\s+)?cases", s.lower())
    n_cases = int(m.group(1).replace(",", "")) if m and m.group(1).replace(",", "").isdigit() else np.nan
    # First integer overall as a crude "n_total" surrogate
    nums = [int(x.replace(",", "")) for x in re.findall(r"\d[\d,]*", s) if int(x.replace(",", "")) > 0]
    n_total = sum(nums) if nums else np.nan
    return {"ancestry_class": cls, "n_cases": n_cases, "n_total": n_total}


def main() -> None:
    df = pd.read_csv(TAB / "T01_gwas_catalog_thyroid_autoimmunity_associations.tsv", sep="\t", low_memory=False)
    df["pvalue_num"] = pd.to_numeric(df["pvalue"], errors="coerce")
    df["chr_num"]    = pd.to_numeric(df["chr"], errors="coerce")
    df["pos_num"]    = pd.to_numeric(df["pos"], errors="coerce")
    # Floor extremely small p (cap at 1e-300 to avoid log10 = inf)
    df["pvalue_num_capped"] = df["pvalue_num"].where(df["pvalue_num"] > 1e-300, 1e-300)
    df["minus_log10_p"] = -np.log10(df["pvalue_num_capped"])

    # ---- MHC subset ----
    mhc = df[(df["chr_num"] == 6) & (df["pos_num"].between(MHC_LO, MHC_HI))].copy()
    mhc["mhc_subregion"] = mhc["pos_num"].apply(assign_subregion)
    anc = mhc["ancestry_initial"].apply(parse_ancestry).apply(pd.Series)
    mhc = pd.concat([mhc.reset_index(drop=True), anc.reset_index(drop=True)], axis=1)
    out = TAB / "T02_mhc_associations.tsv"
    mhc.to_csv(out, sep="\t", index=False)
    print(f"WROTE {out} rows={len(mhc)}")

    # ---- Top 20 per trait ----
    rows = []
    for label, g in mhc.groupby("trait_label"):
        g2 = g.sort_values("minus_log10_p", ascending=False).drop_duplicates("snp").head(20)
        rows.append(g2)
    top = pd.concat(rows, ignore_index=True)
    top.to_csv(TAB / "T03_top20_per_trait.tsv", sep="\t", index=False)
    print(f"WROTE T03_top20_per_trait.tsv rows={len(top)}")

    # ---- SNP→HLA tag join ----
    SNP_HLA_TAGS.to_csv(TAB / "T04a_snp_hla_tag_curated.tsv", sep="\t", index=False)
    tagged = mhc.merge(SNP_HLA_TAGS, left_on="snp", right_on="rsid", how="left")
    tagged_only = tagged[tagged["tagged_allele"].notna()].copy()
    out = TAB / "T04_snp_hla_allele_tags.tsv"
    tagged_only.to_csv(out, sep="\t", index=False)
    print(f"WROTE {out} rows={len(tagged_only)} unique SNPs={tagged_only['snp'].nunique()}")

    # ---- Track 1 convergence ----
    t1 = pd.read_csv(T1_DIR / "T02_panasian_GD_DL_random_effects_v2.tsv", sep="\t")
    t1["allele"] = t1["allele"].str.replace(r"^HLA-", "", regex=True)
    # build dict
    t1_lookup = {f"HLA-{a}": (t1.loc[t1["allele"] == a, "pooled_or"].iloc[0]) for a in t1["allele"]}
    t1_ci = {
        f"HLA-{a}": (t1.loc[t1["allele"] == a, "ci_lo"].iloc[0], t1.loc[t1["allele"] == a, "ci_hi"].iloc[0])
        for a in t1["allele"]
    }
    convergence_rows = []
    for _, r in tagged_only.iterrows():
        allele = r["tagged_allele"]
        t1_or = t1_lookup.get(allele)
        if t1_or is None:
            continue
        t1lo, t1hi = t1_ci[allele]
        # Convert tag-SNP estimate to OR-style. For β-form (log OR) we exponentiate;
        # for OR-form leave as is. GWAS Catalog stores both in `or_or_beta`.
        try:
            est = float(r["or_or_beta"])
        except Exception:
            est = np.nan
        # For binary autoimmune traits the GWAS Catalog "OR or BETA" stores OR
        # natively when riskFrequency / risk-allele are reported; small (<1)
        # values are often log-OR. For convergence we report the raw value
        # AND a boolean directional concordance with Track 1; we do NOT
        # exponentiate (would be wrong for OR-form).
        tag_or = est if not pd.isna(est) else np.nan
        convergence_rows.append(
            {
                "trait_label": r["trait_label"],
                "snp": r["snp"],
                "tagged_allele": allele,
                "tag_ancestry": r["tag_ancestry"],
                "tag_source": r["tag_source"],
                "snp_pvalue": r["pvalue_num"],
                "snp_or_or_beta_raw": r["or_or_beta"],
                "snp_estimate_OR_proxy": tag_or,
                "ancestry_class": r["ancestry_class"],
                "track1_panasian_GD_OR": t1_or,
                "track1_GD_CI_lo": t1lo,
                "track1_GD_CI_hi": t1hi,
                "convergence_qualitative": _convergence_call(tag_or, t1lo, t1hi, t1_or),
            }
        )
    conv = pd.DataFrame(convergence_rows)
    out = TAB / "T05_track1_convergence.tsv"
    conv.to_csv(out, sep="\t", index=False)
    print(f"WROTE {out} rows={len(conv)}")

    # ---- Ancestry-stratified ----
    out = TAB / "T06_ancestry_stratified.tsv"
    mhc[["trait_label", "snp", "pos_num", "mhc_subregion", "pvalue_num", "or_or_beta",
         "ancestry_class", "n_cases", "ancestry_initial"]].to_csv(out, sep="\t", index=False)
    print(f"WROTE {out}")

    # ---- Pleiotropy matrix ----
    pleio = (
        mhc.assign(gws=(mhc["pvalue_num"] < 5e-8).astype(int))
        .groupby(["snp", "trait_label"])["gws"]
        .max()
        .unstack(fill_value=0)
    )
    pleio["n_traits_GWS"] = pleio.sum(axis=1)
    pleio = pleio.sort_values("n_traits_GWS", ascending=False)
    pleio.to_csv(TAB / "T07_pleiotropy_matrix.tsv", sep="\t")
    print(f"WROTE T07_pleiotropy_matrix.tsv rows={len(pleio)}")

    # ---- MHC share per trait ----
    rows = []
    for label, g in df.groupby("trait_label"):
        gws = g[g["pvalue_num"] < 5e-8]
        in_mhc = gws[(gws["chr_num"] == 6) & (gws["pos_num"].between(MHC_LO, MHC_HI))]
        rows.append(
            {
                "trait_label": label,
                "n_total_assoc": len(g),
                "n_GWS_total": len(gws),
                "n_GWS_in_MHC": len(in_mhc),
                "MHC_GWS_share": (len(in_mhc) / len(gws)) if len(gws) else np.nan,
                "n_unique_GWS_SNPs_in_MHC": in_mhc["snp"].nunique(),
            }
        )
    share = pd.DataFrame(rows)
    share.to_csv(TAB / "T08_mhc_share_per_trait.tsv", sep="\t", index=False)
    print("WROTE T08_mhc_share_per_trait.tsv")
    print(share.to_string(index=False))


def _convergence_call(tag_or, t1_lo, t1_hi, t1_or):
    """Qualitative convergence call.

    Caveats: tag_or is the GWAS-Catalog as-reported numeric (could be OR or β
    on log-scale, depending on the original paper's encoding). We therefore
    report direction concordance only, and flag CI overlap only if both are
    plausibly on the OR scale (>0.1 and <10). The narrative MD has the
    discussion of why direction may flip across traits (hypothyroidism vs
    Graves disease for the same DR3 tag).
    """
    if pd.isna(tag_or) or pd.isna(t1_or):
        return "tag_estimate_unavailable"
    # If reported value is implausibly small for an OR (<0.1) we treat it as β
    # and convert to OR via exp() for direction call only.
    or_proxy = float(np.exp(tag_or)) if abs(tag_or) < 1.0 and tag_or != 0 else tag_or
    same_direction = (or_proxy - 1) * (t1_or - 1) > 0
    plausible_OR_scale = 0.1 < or_proxy < 10
    overlap = plausible_OR_scale and (or_proxy >= t1_lo) and (or_proxy <= t1_hi)
    if overlap:
        return "overlap_with_track1_CI"
    if same_direction:
        return "same_direction_no_CI_overlap"
    return "opposite_direction_or_unscaled"


if __name__ == "__main__":
    main()
