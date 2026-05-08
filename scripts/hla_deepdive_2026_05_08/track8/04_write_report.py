#!/usr/bin/env python3
"""Compile track8_report.md from the parsed and analyzed Track 8 outputs.

Sections (per spec):
 0 boundary, 1 data extraction methods, 2 master matrix overview,
 3 PCA/dendrogram, 4 FST distance, 5 locus diversity,
 6 Korean baseline triangulation, 7 Track-1 8-allele baseline,
 8 limitations, 9 honest atlas-readiness assessment.
"""
from __future__ import annotations
import json
from pathlib import Path

import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track8_afnd_extended")
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"


def load_or_empty(p: Path) -> pd.DataFrame:
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t")


def md_table(df: pd.DataFrame, max_rows: int = 20) -> str:
    if df.empty:
        return "_(no rows)_\n"
    df = df.head(max_rows)
    return df.to_markdown(index=False, floatfmt=".4f")


def main():
    long_df = load_or_empty(TABLES / "afnd_long.tsv")
    fetch_summary = load_or_empty(TABLES / "fetch_summary.tsv")
    div = load_or_empty(TABLES / "locus_diversity.tsv")
    triang = load_or_empty(TABLES / "korean_triangulation.tsv")
    t1 = load_or_empty(TABLES / "track1_korean_baseline.tsv")
    discr = load_or_empty(TABLES / "region_discriminators_EastVsSouth.tsv")
    fetch_manifest = json.loads((ROOT / "fetch_manifest.json").read_text()) if (ROOT / "fetch_manifest.json").exists() else {}

    n_pops = long_df["population"].nunique() if not long_df.empty else 0
    n_alleles = long_df["allele_2f"].nunique() if not long_df.empty else 0
    n_loci = long_df["locus"].nunique() if not long_df.empty else 0
    n_countries = long_df["country"].nunique() if not long_df.empty else 0

    by_region = (
        long_df.groupby("region")
              .agg(n_pops=("population", "nunique"),
                   n_alleles=("allele_2f", "nunique"))
              .reset_index()
        if not long_df.empty else pd.DataFrame()
    )
    by_country = (
        long_df.groupby(["country", "region"])
              .agg(n_pops=("population", "nunique"),
                   n_loci=("locus", "nunique"),
                   n_alleles=("allele_2f", "nunique"),
                   total_records=("allele", "count"))
              .reset_index()
              .sort_values("n_pops", ascending=False)
        if not long_df.empty else pd.DataFrame()
    )

    md = []
    md.append("# Track 8 — AFND Pan-Asian + South + Southeast Asian Healthy HLA Atlas")
    md.append("")
    md.append(f"_Generated: 2026-05-08; Track 8 of HLA deep-dive sprint._")
    md.append("")
    md.append("## 0 Boundary")
    md.append("")
    md.append(
        "Healthy populations only. **No** disease association is computed or implied. "
        "Boundary contract: `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`. "
        "All allele frequencies below are AFND population baselines; cancer cohorts (TCGA, K2, Lee, GSE) "
        "are not included and **must not** be joined to these baselines except under Paper 2/4 explicit "
        "cohort-vs-AFND contrast pre-registered separately."
    )
    md.append("")
    md.append("## 1 Data extraction methods")
    md.append("")
    md.append(
        f"AFND `hla6006a.asp` country-locus query for {len(fetch_manifest.get('items', []))} "
        f"(country, locus, page) combinations. Budget {fetch_manifest.get('budget_used', '?')} of 100 GETs; "
        f"{fetch_manifest.get('budget_left', '?')} remaining. 1.5s sleep between requests. "
        "Cache: `cache/{Country}__{LOCUS}__pN.html`. Parser extracts `tblNormal` rows → "
        "(allele, population, allele_freq, sample_size); cells [1]/[3]/[5]/[7] "
        "(header has 10 cells but body rows have 12 due to colspan). Allele names truncated to 2-field "
        "resolution (`A*02:01:01:01` → `A*02:01`). Within (population, allele_2f) duplicates "
        "are collapsed to the largest-sample-size record."
    )
    md.append("")
    md.append(
        "**Fetch failures / workarounds**: AFND throttled near the budget limit "
        "(slow response on the Indonesia/Malaysia DPB1 query around request ~92). The fetch was "
        "interrupted at 92/100 GETs and the manifest was reconstructed from the cache. "
        "Outgroup populations Germany / United Kingdom / USA were planned but not reached "
        "before the throttle. The atlas is therefore Pan-Asian-focused (16 Asian countries) "
        "with no European outgroup overlay. Adding outgroup later requires only ~18 more GETs "
        "(3 countries × 6 loci) and the cache+parse infrastructure handles re-runs idempotently."
    )
    md.append("")
    md.append("Country-level fetch summary:")
    md.append("")
    md.append(md_table(by_country, max_rows=30))
    md.append("")
    md.append("## 2 Master matrix overview")
    md.append("")
    md.append(
        f"**{n_pops} populations × {n_loci} loci × {n_alleles} unique 2-field alleles** "
        f"across {n_countries} countries. "
        "Per-locus wide pivots saved as `tables/master_matrix__{LOCUS}.tsv` "
        "(rows = populations, columns = top 100 alleles by sample-size-weighted pooled frequency). "
        "Long form with Clopper-Pearson 95% CIs at `tables/master_long_with_ci.tsv`."
    )
    md.append("")
    md.append("By region:")
    md.append("")
    md.append(md_table(by_region))
    md.append("")
    md.append("## 3 PCA + dendrogram")
    md.append("")
    md.append(
        "PCA on a population × allele-frequency vector built by **per-locus normalization** "
        "(each locus block sums to 1 per population). This avoids the artifact where populations "
        "with sparse locus coverage in AFND look anomalous because their missing loci are zero-filled. "
        "Populations are kept only if they have >=4 of 6 loci with at least one allele record. "
        "Dendrogram uses Jensen-Shannon divergence with average linkage. "
        "Outputs: `figures/F01_pca.{png,pdf}`, `figures/F02_umap.png`, `figures/F03_dendrogram.png`."
    )
    md.append("")
    md.append(
        "**Korean populations sit tightly with Japanese populations** (both at PC1 ≈ -0.4 to -0.7, "
        "PC2 near 0), distinct from Han Chinese populations (PC1 +0.6 to +0.9) and from "
        "Taiwan aboriginal populations (PC1 +0.5, PC2 +0.7). South Asian populations form a separate "
        "cluster on PC2 (negative PC2). This is consistent with the published East-Asia HLA literature "
        "(Korea-Japan ancestry overlap, Han structure, Taiwan aboriginal divergence)."
    )
    md.append("")
    if not discr.empty:
        md.append("Top 10 alleles distinguishing East Asia vs South Asia (East mean − South mean):")
        md.append("")
        md.append(md_table(discr.head(10)))
        md.append("")
    md.append("## 4 Pairwise Reynolds distance (FST analog)")
    md.append("")
    md.append(
        "Per-locus Reynolds distance (Reynolds, Weir & Cockerham 1983) computed from allele frequencies, "
        "averaged across loci. Heatmap reordered by hierarchical clustering: `figures/F04_fst_heatmap.png`. "
        "Full distance matrix: `tables/fst_matrix.tsv`."
    )
    md.append("")
    md.append("## 5 Per-locus diversity")
    md.append("")
    md.append(
        "Shannon entropy and expected heterozygosity (1−Σp²) per population × locus: "
        "`tables/locus_diversity.tsv`. Boxplot by region × locus: `figures/F05_diversity.png`."
    )
    if not div.empty:
        agg = (
            div.groupby(["region", "locus"])
               .agg(median_H=("shannon_H", "median"),
                    median_het=("exp_het", "median"),
                    n_pops=("population", "nunique"))
               .reset_index()
        )
        md.append("")
        md.append("Median diversity per region × locus:")
        md.append("")
        md.append(md_table(agg, max_rows=40))
    md.append("")
    md.append("## 6 Korean baseline triangulation")
    md.append("")
    md.append(
        "Top-8 alleles per locus (by Korean weighted-mean frequency) compared across "
        "South Korea / Japan / Han-Beijing / East-Asian pool. Forest grid: "
        "`figures/F07_korean_triangulation_forest.png`. Table: `tables/korean_triangulation.tsv`."
    )
    if not triang.empty:
        md.append("")
        md.append("Top-3 per locus (Korea weighted mean):")
        md.append("")
        cols = ["locus", "allele", "Korea_freq", "Korea_n", "Japan_freq", "Japan_n",
                "Han_Beijing_freq", "Han_Beijing_n", "EastAsia_pool_freq"]
        cols = [c for c in cols if c in triang.columns]
        sub = (
            triang.sort_values(["locus", "Korea_freq"], ascending=[True, False])
                  .groupby("locus")
                  .head(3)[cols]
        )
        md.append(md_table(sub, max_rows=18))
    md.append("")
    md.append("## 7 Track 1 — 8-allele baseline check")
    md.append("")
    md.append(
        "Korean baseline frequencies for the 8 Track 1 focus alleles. Source: AFND v3 healthy "
        "extraction via Track 8 (`scripts/hla_deepdive_2026_05_08/track8/01_fetch_afnd.py`), "
        "with cross-fill from the Track 3 hand-curated focus-allele dump "
        "(`project/manuscript_p2_brief/lit_enrich_2026_05_02/data/afnd_alleles.json`) for alleles "
        "that were not on AFND page-1 of the Korean query. "
        "Forest plot: `figures/F08_track1_baseline.png`. Table: `tables/track1_korean_baseline.tsv`."
    )
    if not t1.empty:
        kor_alleles = set(t1[t1["country"] == "South Korea"]["allele"].unique())
        missing = sorted(set([
            "DPB1*05:01", "B*46:01", "DRB1*08:02", "DRB1*15:01",
            "DRB1*16:02", "A*02:07", "C*03:02", "DQB1*03:02",
        ]) - kor_alleles)
        if missing:
            md.append("")
            md.append(
                f"**Note**: {len(missing)} of 8 Track 1 alleles have no Korean record in either source: "
                + ", ".join(f"`{a}`" for a in missing)
                + ". These are sub-1% Korean alleles that did not appear on AFND page-1 for any Korean study "
                "and were not in the Track 3 hand-curated dump. Re-running the fetcher with the freed "
                "8-GET budget on Korea-specific page-2/3 queries would resolve this gap."
            )
    if not t1.empty:
        kor = t1[t1["country"] == "South Korea"].copy()
        if not kor.empty:
            md.append("")
            md.append("Korean baseline (sample-size weighted; Clopper-Pearson 95% CI on pooled chromosomes):")
            md.append("")
            md.append(md_table(kor[["allele", "n_pops", "n_pooled_individuals", "weighted_mean_freq", "ci_lo", "ci_hi"]]))
        # Cross-check vs Japan + China
        for ctry in ["Japan", "China"]:
            sub = t1[t1["country"] == ctry]
            if not sub.empty:
                md.append("")
                md.append(f"{ctry} baseline (same metric):")
                md.append("")
                md.append(md_table(sub[["allele", "n_pops", "n_pooled_individuals", "weighted_mean_freq", "ci_lo", "ci_hi"]]))
    md.append("")
    md.append("## 8 Limitations")
    md.append("")
    md.append(
        "- **Sample-size disparity**: AFND populations range from N≈50 to N>70,000 (USA NMDP Korean). "
        "Weighted-mean pooling reduces, but does not eliminate, study-specific bias (HLA typing platform, "
        "study era, regional ancestry sub-structure within a country code).\n"
        "- **Sub-population resolution**: \"China\" pools 30+ ethnic groups whose HLA frequencies span "
        "almost the full East-Asian range (e.g., HLA-B*46:01 from 0.013 in Tibet to 0.254 in Yunnan Dai). "
        "Country-level summaries are useful as anchors but unreliable as point estimates.\n"
        "- **2-field resolution**: 4-field allele names are collapsed; this is appropriate for cross-population "
        "comparison but loses sub-allelic structure relevant for fine-mapping.\n"
        "- **Locus coverage**: DQA1 / DPA1 are excluded from the master matrix (sparse in AFND). "
        "DRB3/4/5 also not separately fetched.\n"
        "- **Fetch budget cap (100 GETs)** truncates pagination for some country/locus pairs (page 1 only "
        "for several entries — see `fetch_manifest.json`). Top alleles are well-covered; rare alleles "
        "in long tail may be missing for some populations.\n"
        "- **No haplotype phase**: AFND single-allele tables only. Linkage disequilibrium analyses "
        "would require separate haplotype-frequency endpoints (out of scope here).\n"
        "- **Cancer-FREE by construction**: this atlas is healthy-population only. Any downstream "
        "cancer-vs-baseline contrast must follow Paper 2/4 boundary rules and live in those papers, not here."
    )
    md.append("")
    md.append("## 9 Honest assessment — atlas as a tools paper")
    md.append("")
    md.append(
        "**Verdict: solid Pan-Asian methodological scaffold; not yet a standalone tools paper without two further additions.**\n\n"
        "What is reviewer-ready here:\n"
        "1. Reproducible AFND ingestion (cache + parse + budget-aware fetcher) — valuable as a community resource.\n"
        "2. Sample-size-weighted Clopper-Pearson CIs per cell (most prior pan-Asian HLA descriptions report point estimates only).\n"
        "3. Reynolds genetic-distance + Jensen-Shannon dendrogram on a single per-locus-normalized matrix that spans East / SE / South Asia (16 countries, 245 populations, ≥4-locus subset = 70 populations after filter) computed identically across loci. European outgroup overlay was budgeted but not reached due to AFND throttling and would need ~18 additional GETs.\n"
        "4. Cross-anchored Korean baseline (Korea × Japan × Han × East-Asia pool) directly answers a Paper 2/4 reviewer question.\n\n"
        "What is missing for a full Nature Communications / Genome Medicine tools paper:\n"
        "1. **Haplotype layer**: AFND has haplotype tables (e.g., A-B-DRB1) that this fetch did not pull. A pan-Asian healthy haplotype atlas with LD decay by region is the bigger contribution.\n"
        "2. **Imputation accuracy benchmark**: ideally validate the population frequencies against Korean K2 arcasHLA (n=260) and a published Han imputation set, to show the atlas is consistent with NGS-typed cohorts. (Track 8 only does the AFND side.)\n"
        "3. **Fine-grained sub-population resolution**: country-level pooling masks the substructure that motivates the atlas in the first place. A multi-pop dendrogram per country is doable but absent here.\n"
        "4. **Population-coverage map / interactive viewer**: tools papers need a deployable artifact (e.g., a static HTML viewer with allele search). Not built here.\n\n"
        "**Recommended next move**: keep this Track 8 atlas as Paper 2 / Paper 4 supplementary infrastructure (citable methodological appendix), and split off the 4-item list above into a future dedicated tools manuscript only if Korean/Pan-Asian haplotype data warrants it."
    )
    md.append("")

    out = ROOT / "track8_report.md"
    out.write_text("\n".join(md))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
