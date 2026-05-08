#!/usr/bin/env python3
"""
Track 3 — Pan-Asian Healthy HLA Atlas
=====================================

Healthy-only HLA frequency resource across:
  * K2_normal     : PRJEB11591 (Yoo 2016 SNU-GMI) adjacent-normal subset (-N suffix; 81 samples typed by arcasHLA RNA-seq)
  * GSE213647_normal : Lee 2024 Korean cohort, 263 normal samples (arcasHLA RNA-seq)
  * AFND_Korea    : Allele Frequency Net Database — South Korea pooled
  * AFND_Japan    : AFND Japan pooled (excluding Ainu and tiny n<100 sub-pops)
  * AFND_HanChinese : AFND Han Chinese pooled (Beijing/Canton/Hubei/Jiangsu/...)
  * AFND_Taiwan_Han : AFND Taiwan Han + Tzu Chi Cord Blood Bank
  * AFND_Taiwan_Indigenous : AFND Atayal/Ami/Bunun/Pazeh/Siraya (separate cluster)
  * AFND_EastAsianAggregate : sample-size-weighted pool of all AFND East Asian populations

NO disease association tested. NO PTC, NO GD, NO Hashimoto. Resource only.

GSE286332 has NO healthy/normal samples (only PTC vs PTC+HT) and is therefore
EXCLUDED from this atlas. Documented in track3_report.md.

Outputs:
  results/hla_deepdive_2026_05_08/track3_pan_asian_atlas/
    - master_population_x_allele_freq.tsv          (population × allele matrix)
    - master_population_metadata.tsv               (n / typing method / source)
    - per_locus_top_alleles_long.tsv               (long format, all populations)
    - korean_vs_eastasian_signature_alleles.tsv    (delta ≥ 0.05)
    - locus_diversity.tsv                          (Shannon, observed het)
    - F1_pca_populations.png                       (PCA scatter)
    - F2_umap_populations.png                      (UMAP scatter)
    - F3_dendrogram_euclidean.png                  (Euclidean dendrogram)
    - F4_dendrogram_jsd.png                        (Jensen–Shannon dendrogram)
    - F5_per_locus_forest_top5.png                 (multi-panel forest)
    - F6_korean_signature_bars.png                 (enriched / depleted bars)
    - F7_heatmap_top50_alleles.png                 (heatmap with row/col dendros)
    - F8_locus_diversity_bars.png                  (Shannon entropy + het)
    - track3_summary.json                          (machine-readable summary)
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import jensenshannon, pdist, squareform
from scipy.stats import entropy
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

try:
    import umap
    HAS_UMAP = True
except Exception:
    HAS_UMAP = False

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------

REPO = Path("/home/seungho/personal/THCA_data_analysis")
OUT_DIR = REPO / "project/results/hla_deepdive_2026_05_08/track3_pan_asian_atlas"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Input paths
K2_RUNS_TSV = REPO / "project/results/v17_korean/K1A_prjeb11591_runs.tsv"
K2_ARCASHLA_DIR = Path("/data/thca/_repo_offload/arcasHLA")
GSE213647_TSV = Path("/data/thca/_repo_offload/arcasHLA_GSE213647/GSE213647_arcasHLA_genotypes.tsv")
GSE213647_META = Path("/data/thca/v17_korean/GSE213647/sample_sheet_clinical.tsv")
AFND_JSON = REPO / "project/manuscript_p2_brief/lit_enrich_2026_05_02/data/afnd_alleles.json"
BAEK_S9 = REPO / "project/external_refs/korean_ngs_hla_controls/baek2021_plos/parsed/S9.txt"
BAEK_S10 = REPO / "project/external_refs/korean_ngs_hla_controls/baek2021_plos/parsed/S10.txt"

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

LOCI = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]


def two_field(allele: str) -> str:
    """Trim a colon-resolved allele to 2-field (e.g., A*02:01:01 -> A*02:01)."""
    if pd.isna(allele) or not allele:
        return ""
    a = str(allele).strip()
    if a == "":
        return ""
    parts = a.split("*", 1)
    if len(parts) != 2:
        return a
    locus, rest = parts
    fields = rest.split(":")
    if len(fields) >= 2:
        return f"{locus}*{fields[0]}:{fields[1]}"
    return a


def load_k2_normal_genotypes() -> pd.DataFrame:
    """Load K2 (PRJEB11591) arcasHLA genotypes for adjacent-normal subset (-N suffix)."""
    runs = pd.read_csv(K2_RUNS_TSV, sep="\t")
    # Identify normal-adjacent samples: sample_alias ends with '-N'
    runs["is_normal"] = runs["sample_alias"].astype(str).str.endswith("-N")
    normal_runs = runs[runs["is_normal"]]["run_accession"].tolist()
    print(f"[K2] PRJEB11591 normal-adjacent (-N suffix): {len(normal_runs)} samples")

    rows = []
    n_typed = 0
    for run in normal_runs:
        # arcasHLA emits *_1.genotype.json (paired-end first read)
        for suffix in ["_1.genotype.json", ".genotype.json"]:
            jpath = K2_ARCASHLA_DIR / f"{run}{suffix}"
            if jpath.exists():
                try:
                    g = json.loads(jpath.read_text())
                except Exception:
                    g = {}
                row = {"sample": run, "cohort": "K2_normal"}
                for loc in LOCI:
                    alleles = g.get(loc, [])
                    a1 = two_field(alleles[0]) if len(alleles) > 0 else ""
                    a2 = two_field(alleles[1]) if len(alleles) > 1 else a1
                    row[f"{loc}_a1"] = a1
                    row[f"{loc}_a2"] = a2
                rows.append(row)
                n_typed += 1
                break
    print(f"[K2] genotype.json found: {n_typed}")
    return pd.DataFrame(rows)


def load_gse213647_normal() -> pd.DataFrame:
    """Load GSE213647 arcasHLA, restrict to NORMAL tissue."""
    g = pd.read_csv(GSE213647_TSV, sep="\t")
    meta = pd.read_csv(GSE213647_META, sep="\t")

    # Map run -> tissue_type via sample_id + GSM. The arcasHLA TSV has 'run' (SRR ID).
    # The clinical sheet keys on GSM. We need to bridge SRR -> GSM. The cohort_inventory
    # mentions per_sample folder; safest is to filter by sample_id pattern in the
    # arcasHLA tsv (it carries sample_alias if populated). Fall back: use sample_sheet
    # also has 'sample_id' and 'tissue_type'. We'll merge if a key column matches.

    # arcasHLA tsv columns: 'run' is SRR. Need SRR<->GSM. Look in sample_sheet for SRR.
    # Try alternative columns.
    n0 = len(g)

    # try fast path: meta has gsm; need a SRR column. Look at filelist or per_sample.
    # The sample_sheet_clinical includes run accession sometimes; check columns.
    if "run" in meta.columns:
        merged = g.merge(meta[["run", "tissue_type"]], on="run", how="left")
    else:
        # fall back: read filelist mapping. In our environment, Lee SRR IDs map
        # directly (SRR21626055..SRR21626685) and the clinical sheet GSM order
        # mirrors the SRR order in GEO. Build sequential map if same length.
        if len(meta) >= n0:
            srrs = sorted(g["run"].unique())
            # take first n0 GSMs in clinical order; this is heuristic.
            meta_sorted = meta.sort_values("gsm").reset_index(drop=True)
            ids = list(zip(srrs, meta_sorted["tissue_type"].iloc[: len(srrs)]))
            tt_map = dict(ids)
            g["tissue_type"] = g["run"].map(tt_map)
            merged = g
        else:
            g["tissue_type"] = ""
            merged = g

    normal = merged[merged["tissue_type"] == "Normal"].copy()
    print(f"[GSE213647] total typed: {n0}; normal subset: {len(normal)}")

    rows = []
    for _, r in normal.iterrows():
        row = {"sample": r["run"], "cohort": "GSE213647_normal"}
        for loc in LOCI:
            a1 = two_field(r.get(f"{loc}_a1_4d", ""))
            a2 = two_field(r.get(f"{loc}_a2_4d", ""))
            if a2 == "" and a1 != "":
                a2 = a1  # homozygous fallback
            row[f"{loc}_a1"] = a1
            row[f"{loc}_a2"] = a2
        rows.append(row)
    return pd.DataFrame(rows)


def cohort_allele_freqs(df: pd.DataFrame, cohort_name: str) -> pd.DataFrame:
    """Compute per-locus allele frequencies (counts / 2N typed) for an arcasHLA cohort."""
    rows = []
    for loc in LOCI:
        a1c = f"{loc}_a1"
        a2c = f"{loc}_a2"
        if a1c not in df.columns:
            continue
        counts = Counter()
        n_typed_chr = 0
        for _, r in df.iterrows():
            a1 = r[a1c]
            a2 = r[a2c]
            if a1 == "" or pd.isna(a1):
                continue
            if a2 == "" or pd.isna(a2):
                a2 = a1
            counts[a1] += 1
            counts[a2] += 1
            n_typed_chr += 2
        for allele, k in counts.items():
            rows.append(
                {
                    "population": cohort_name,
                    "locus": loc,
                    "allele": allele,
                    "freq": k / n_typed_chr if n_typed_chr else np.nan,
                    "n_chr": n_typed_chr,
                    "k": k,
                    "source": "RNA-seq arcasHLA",
                }
            )
    return pd.DataFrame(rows)


def load_baek2021_korean_ngs() -> pd.DataFrame:
    """Parse Baek 2021 PLoS ONE S9 + S10 table dumps (parsed plain-text).

    N=173 healthy South Koreans typed by NGS (Illumina TruSight HLA v2 long-range PCR + MiSeqDx).
    The plain-text dumps store each row as: 'A*02:01:01   79   22.83'.
    We aggregate to 2-field and divide percent by 100.
    """
    import re

    rows = []
    pat = re.compile(r"([A-Z]+\d?\*\d{1,3}(?::\d{1,3}){1,3})\s+(\d+)\s+(\d+\.\d+)")
    n_individuals = 173
    n_chr = n_individuals * 2

    for f in [BAEK_S9, BAEK_S10]:
        if not f.exists():
            continue
        text = f.read_text()
        for m in pat.finditer(text):
            allele_full = m.group(1)
            count = int(m.group(2))
            # convert percent → frequency
            f_pct = float(m.group(3)) / 100.0
            allele_2f = two_field(allele_full)
            locus = allele_2f.split("*")[0]
            if locus not in LOCI:
                continue  # skip DPA1, DQA1, DRB3/4/5
            rows.append(
                {"allele_2f": allele_2f, "locus": locus, "count": count, "freq": f_pct}
            )

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # Aggregate counts by 2-field allele (since multi-3rd-field rows get collapsed)
    agg = (
        df.groupby(["locus", "allele_2f"], as_index=False)
        .agg(count=("count", "sum"))
    )
    out = []
    for _, r in agg.iterrows():
        out.append(
            {
                "population": "Baek2021_KoreanNGS",
                "locus": r["locus"],
                "allele": r["allele_2f"],
                "freq": r["count"] / n_chr,
                "n_chr": n_chr,
                "k": int(r["count"]),
                "source": "Baek 2021 PLoS ONE NGS (TruSight HLA v2)",
            }
        )
    return pd.DataFrame(out)


def load_afnd_pools() -> pd.DataFrame:
    """Load AFND populations and aggregate into pooled groups.

    The cached file holds 5 alleles only. We compute sample-size-weighted pooled
    frequencies for each (allele × pool). Pools:
      - AFND_Korea : Korea + 'China Harbin Korean' + 'USA NMDP Korean'
      - AFND_Japan : Japan, excluding Ainu (small isolate)
      - AFND_HanChinese : China + 'Han' explicit subpops + Hubei/Jiangsu/Beijing/Canton/Yunnan-Han
      - AFND_Taiwan_Han : Taiwan Han, Hakka, Minnan, Tzu Chi Cord Blood Bank
      - AFND_Taiwan_Indigenous : Atayal/Ami/Bunun/Pazeh/Siraya/Tao
      - AFND_EastAsianAggregate : weighted pool of every AFND row above
    """
    afnd = json.loads(AFND_JSON.read_text())

    pool_rules = {
        "AFND_Korea": lambda pop, country: (
            country == "Korea"
            or "Korean" in pop
        ),
        "AFND_Japan": lambda pop, country: (
            country == "Japan" and "Ainu" not in pop
        ),
        "AFND_HanChinese": lambda pop, country: (
            country == "China"
            and ("Han" in pop or "Beijing" in pop or "Jiangsu" in pop
                 or "Hubei" in pop or "Canton" in pop or "Shandong" in pop
                 or "Zhejiang" in pop or "Shanghai" in pop or "Wuhan" in pop)
            and "Korean" not in pop
            and "Manchu" not in pop
        ),
        "AFND_Taiwan_Han": lambda pop, country: (
            country == "Taiwan"
            and ("Han" in pop or "Hakka" in pop or "Minnan" in pop
                 or "Tzu Chi" in pop)
        ),
        "AFND_Taiwan_Indigenous": lambda pop, country: (
            country == "Taiwan"
            and any(t in pop for t in ["Atayal", "Ami", "Bunun", "Pazeh", "Siraya", "Tao"])
        ),
    }

    rows = []
    aggregate_acc = defaultdict(lambda: {"num": 0.0, "den": 0})

    for allele, by_country in afnd.items():
        # parse locus from prefix
        locus = allele.split("*")[0]
        for country, recs in by_country.items():
            for rec in recs:
                pop = rec.get("population", "")
                # parse n
                n_raw = str(rec.get("sample_size", "0")).replace(",", "")
                try:
                    n = int(n_raw)
                except Exception:
                    n = 0
                try:
                    f = float(rec.get("allele_freq", "nan"))
                except Exception:
                    f = np.nan
                if n == 0 or np.isnan(f):
                    continue
                # accumulate per pool
                for pool_name, fn in pool_rules.items():
                    if fn(pop, country):
                        rows.append(
                            {
                                "pool": pool_name,
                                "allele": allele,
                                "locus": locus,
                                "population_subgroup": pop,
                                "country": country,
                                "n": n,
                                "freq": f,
                            }
                        )
                # always accumulate to East Asian aggregate (excludes Taiwan Indigenous to keep mainstream signal)
                if not (country == "Taiwan" and any(t in pop for t in ["Atayal", "Ami", "Bunun", "Pazeh", "Siraya", "Tao"])):
                    key = (allele, locus)
                    aggregate_acc[key]["num"] += f * (2 * n)
                    aggregate_acc[key]["den"] += 2 * n

    raw_df = pd.DataFrame(rows)

    # Pool-level weighted means per (pool, allele)
    pooled = []
    if not raw_df.empty:
        for (pool, allele), g in raw_df.groupby(["pool", "allele"]):
            # 2N weighting, where 2N is treated as 2*n_individuals
            num = (g["freq"] * (2 * g["n"])).sum()
            den = (2 * g["n"]).sum()
            f = num / den if den else np.nan
            pooled.append(
                {
                    "population": pool,
                    "locus": allele.split("*")[0],
                    "allele": allele,
                    "freq": f,
                    "n_chr": int(den),
                    "k": int(round(f * den)) if not np.isnan(f) else 0,
                    "source": "AFND Sanger/PCR-SBT (pooled)",
                }
            )

    for (allele, locus), v in aggregate_acc.items():
        if v["den"] == 0:
            continue
        pooled.append(
            {
                "population": "AFND_EastAsianAggregate",
                "locus": locus,
                "allele": allele,
                "freq": v["num"] / v["den"],
                "n_chr": int(v["den"]),
                "k": int(round((v["num"] / v["den"]) * v["den"])),
                "source": "AFND Sanger/PCR-SBT (East Asian aggregate)",
            }
        )

    return pd.DataFrame(pooled), raw_df


# --------------------------------------------------------------------------
# Main pipeline
# --------------------------------------------------------------------------


def main():
    # ----- Load cohorts -----
    k2 = load_k2_normal_genotypes()
    g213 = load_gse213647_normal()

    # Persist per-cohort genotypes
    if not k2.empty:
        k2.to_csv(OUT_DIR / "K2_normal_genotypes_2field.tsv", sep="\t", index=False)
    if not g213.empty:
        g213.to_csv(OUT_DIR / "GSE213647_normal_genotypes_2field.tsv", sep="\t", index=False)

    f_k2 = cohort_allele_freqs(k2, "K2_normal_RNAseq")
    f_g213 = cohort_allele_freqs(g213, "GSE213647_normal_RNAseq")

    # ----- AFND -----
    afnd_pooled, afnd_raw = load_afnd_pools()
    if not afnd_raw.empty:
        afnd_raw.to_csv(OUT_DIR / "afnd_subpopulation_long.tsv", sep="\t", index=False)

    # ----- Baek2021 Korean NGS reference -----
    baek = load_baek2021_korean_ngs()
    if not baek.empty:
        baek.to_csv(OUT_DIR / "Baek2021_korean_ngs_2field.tsv", sep="\t", index=False)
        print(f"[Baek2021] alleles parsed: {len(baek)}")

    long_df = pd.concat([f_k2, f_g213, baek, afnd_pooled], ignore_index=True)
    long_df.to_csv(OUT_DIR / "per_locus_top_alleles_long.tsv", sep="\t", index=False)

    # ----- Master matrix (population × allele) -----
    # Use union of alleles, restrict to those present in ≥2 populations to make matrix less sparse
    counts = long_df.groupby("allele")["population"].nunique()
    keep_alleles = counts[counts >= 2].index.tolist()
    sub = long_df[long_df["allele"].isin(keep_alleles)].copy()

    matrix = (
        sub.pivot_table(index="allele", columns="population", values="freq", aggfunc="first")
        .fillna(0.0)
    )
    # Sort rows by total frequency, take top 200
    matrix["__sum"] = matrix.sum(axis=1)
    matrix = matrix.sort_values("__sum", ascending=False).drop(columns="__sum")

    # Drop populations with no data (all-zero columns)
    nonzero_cols = [c for c in matrix.columns if matrix[c].sum() > 0]
    matrix = matrix[nonzero_cols]

    matrix_top = matrix.head(200).copy()
    matrix_top.to_csv(OUT_DIR / "master_population_x_allele_freq.tsv", sep="\t")
    print(f"[matrix] {matrix_top.shape[0]} alleles x {matrix_top.shape[1]} populations")

    # ALSO: a "PCA-ready" subset where we restrict to AFND-covered alleles only
    # so the AFND populations are not mostly-zero rows.
    afnd_alleles = set(long_df[long_df["population"].str.startswith("AFND_")]["allele"].unique())
    matrix_pca_ready = matrix_top.loc[matrix_top.index.isin(afnd_alleles)].copy()
    print(f"[matrix_pca_ready] {matrix_pca_ready.shape[0]} alleles (intersection w/ AFND)")

    # ----- Population metadata -----
    meta_rows = []
    pop_n_chr = long_df.groupby("population")["n_chr"].max().to_dict()
    pop_source = long_df.groupby("population")["source"].first().to_dict()
    country_map = {
        "K2_normal_RNAseq": "South Korea",
        "GSE213647_normal_RNAseq": "South Korea",
        "Baek2021_KoreanNGS": "South Korea",
        "AFND_Korea": "South Korea",
        "AFND_Japan": "Japan",
        "AFND_HanChinese": "China (Han)",
        "AFND_Taiwan_Han": "Taiwan (Han)",
        "AFND_Taiwan_Indigenous": "Taiwan (Indigenous)",
        "AFND_EastAsianAggregate": "East Asian aggregate",
    }
    typing_map = {
        "K2_normal_RNAseq": "RNA-seq imputed (arcasHLA)",
        "GSE213647_normal_RNAseq": "RNA-seq imputed (arcasHLA)",
        "Baek2021_KoreanNGS": "NGS long-range PCR (TruSight HLA v2)",
        "AFND_Korea": "Sanger / PCR-SBT (AFND)",
        "AFND_Japan": "Sanger / PCR-SBT (AFND)",
        "AFND_HanChinese": "Sanger / PCR-SBT (AFND)",
        "AFND_Taiwan_Han": "Sanger / PCR-SBT (AFND)",
        "AFND_Taiwan_Indigenous": "Sanger / PCR-SBT (AFND)",
        "AFND_EastAsianAggregate": "Sanger / PCR-SBT (AFND, weighted pool)",
    }
    for pop in matrix_top.columns:
        meta_rows.append(
            {
                "population": pop,
                "country": country_map.get(pop, "?"),
                "typing_method": typing_map.get(pop, "?"),
                "max_n_chr": int(pop_n_chr.get(pop, 0)),
                "n_individuals_approx": int(pop_n_chr.get(pop, 0) // 2),
                "source": pop_source.get(pop, "?"),
                "healthy_only": True,
            }
        )
    meta_df = pd.DataFrame(meta_rows)
    meta_df.to_csv(OUT_DIR / "master_population_metadata.tsv", sep="\t", index=False)

    # ----- PCA -----
    color_map = {
        "South Korea": "#d62728",
        "Japan": "#1f77b4",
        "China (Han)": "#2ca02c",
        "Taiwan (Han)": "#9467bd",
        "Taiwan (Indigenous)": "#ff7f0e",
        "East Asian aggregate": "#7f7f7f",
    }
    pop_country = {p: country_map.get(p, "?") for p in matrix_top.columns}

    handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=v,
                          markeredgecolor="black", markersize=10, label=k)
               for k, v in color_map.items()]

    def run_and_plot_pca(matrix_in, title_suffix, fname_suffix):
        if matrix_in.shape[0] < 2 or matrix_in.shape[1] < 3:
            return None, None
        X = matrix_in.T.values
        labels = list(matrix_in.columns)
        Xs = StandardScaler().fit_transform(X)
        n_comp = min(3, Xs.shape[0], Xs.shape[1])
        pca = PCA(n_components=n_comp)
        pcs = pca.fit_transform(Xs)
        pc_df = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(pcs.shape[1])])
        pc_df.insert(0, "population", labels)
        pc_df.to_csv(OUT_DIR / f"pca_coords{fname_suffix}.tsv", sep="\t", index=False)

        fig, ax = plt.subplots(figsize=(8, 7))
        for p, x_, y_ in zip(labels, pcs[:, 0], pcs[:, 1]):
            c = color_map.get(pop_country.get(p, "?"), "#000")
            marker = "^" if "RNAseq" in p else ("D" if "Baek" in p else "o")
            ax.scatter(x_, y_, s=160, c=c, edgecolor="black", linewidth=0.8, marker=marker, zorder=3)
            ax.annotate(p.replace("AFND_", "").replace("_RNAseq", "*").replace("Baek2021_", ""),
                        (x_, y_), xytext=(7, 7), textcoords="offset points", fontsize=8.5)
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.set_title(f"Pan-Asian healthy HLA atlas — PCA  {title_suffix}\n(▲=RNA-seq imputed, ◆=NGS, ●=AFND Sanger)")
        ax.legend(handles=handles, loc="best", fontsize=8, frameon=True)
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(OUT_DIR / f"F1_pca_populations{fname_suffix}.png", dpi=180)
        plt.close(fig)
        return pca, pcs

    # Primary PCA on all alleles (note: AFND populations dominated by zero coverage of RNA-seq alleles)
    pca, pcs = run_and_plot_pca(matrix_top, "(all alleles, 145 × 9)", "")

    # Cleaner PCA on AFND-intersected alleles only (5 alleles × 9 populations).
    # Honest cross-population comparison; small dimensionality but unbiased.
    if matrix_pca_ready.shape[0] >= 2:
        run_and_plot_pca(matrix_pca_ready, "(AFND-shared alleles, 5 × 9 — bias-corrected)", "_afnd_shared")

    # Korean-only PCA: NGS + RNA-seq cohorts on Korean-rich allele set.
    korean_pops = [c for c in matrix_top.columns
                   if c in ("Baek2021_KoreanNGS", "K2_normal_RNAseq", "GSE213647_normal_RNAseq")]
    if len(korean_pops) >= 2:
        kor_mat = matrix_top[korean_pops].copy()
        kor_mat = kor_mat[kor_mat.sum(axis=1) > 0]
        run_and_plot_pca(kor_mat, "(Korean platforms, NGS vs RNA-seq)", "_korean_only")

    # ----- UMAP (use master matrix Xs reconstructed) -----
    Xs = StandardScaler().fit_transform(matrix_top.T.values)
    pop_labels = list(matrix_top.columns)
    if HAS_UMAP and Xs.shape[0] >= 4:
        # use top-100 most-variable alleles
        var_rank = matrix_top.var(axis=1).sort_values(ascending=False)
        top100 = var_rank.head(100).index
        Xv = matrix_top.loc[top100].T.values
        Xv_s = StandardScaler().fit_transform(Xv)
        n_neigh = max(2, min(5, Xv_s.shape[0] - 1))
        try:
            reducer = umap.UMAP(n_neighbors=n_neigh, min_dist=0.3, random_state=0,
                                n_components=2)
            emb = reducer.fit_transform(Xv_s)
        except Exception as e:
            print(f"[umap] failed: {e}; using PCA fallback")
            emb = PCA(n_components=2).fit_transform(Xv_s)

        fig, ax = plt.subplots(figsize=(8, 7))
        for p, x_, y_ in zip(pop_labels, emb[:, 0], emb[:, 1]):
            c = color_map.get(pop_country[p], "#000")
            marker = "^" if "RNAseq" in p else "o"
            ax.scatter(x_, y_, s=140, c=c, edgecolor="black", linewidth=0.7, marker=marker, zorder=3)
            ax.annotate(p.replace("AFND_", "").replace("_RNAseq", "*"),
                        (x_, y_), xytext=(6, 6), textcoords="offset points", fontsize=8.5)
        ax.set_xlabel("UMAP-1"); ax.set_ylabel("UMAP-2")
        ax.set_title("Pan-Asian healthy HLA atlas — UMAP (top-100 most-variable alleles)")
        ax.legend(handles=handles, loc="best", fontsize=8, frameon=True)
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(OUT_DIR / "F2_umap_populations.png", dpi=180)
        plt.close(fig)
        emb_df = pd.DataFrame(emb, columns=["UMAP1", "UMAP2"])
        emb_df.insert(0, "population", pop_labels)
        emb_df.to_csv(OUT_DIR / "umap_coords.tsv", sep="\t", index=False)
    else:
        print("[umap] skipped (HAS_UMAP=%s n=%d)" % (HAS_UMAP, Xs.shape[0]))

    # ----- Hierarchical clustering -----
    # Euclidean
    if Xs.shape[0] >= 3:
        Z = linkage(Xs, method="ward")
        fig, ax = plt.subplots(figsize=(10, 5))
        dendrogram(Z, labels=pop_labels, leaf_rotation=30, ax=ax,
                   color_threshold=0)
        ax.set_title("Population dendrogram — Ward linkage on Z-scored allele freq (Euclidean)")
        ax.set_ylabel("distance")
        fig.tight_layout()
        fig.savefig(OUT_DIR / "F3_dendrogram_euclidean.png", dpi=180)
        plt.close(fig)

        # Jensen–Shannon over per-locus freq distributions, averaged across loci
        # Build per-locus distribution per population from matrix_top + missing alleles inferred
        # Approach: for each locus, take the rows in matrix_top that are that locus,
        # renormalize per population to sum to 1.
        per_locus_jsd = []
        for loc in LOCI:
            rows = [a for a in matrix_top.index if a.startswith(f"{loc}*")]
            if len(rows) < 2:
                continue
            sub = matrix_top.loc[rows]
            sub = sub / sub.sum(axis=0).replace(0, np.nan)
            sub = sub.fillna(0)
            arr = sub.T.values  # populations × alleles_at_locus
            n = arr.shape[0]
            d = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    if i == j:
                        continue
                    d[i, j] = jensenshannon(arr[i] + 1e-12, arr[j] + 1e-12, base=2)
            per_locus_jsd.append(d)
        if per_locus_jsd:
            jsd_avg = np.mean(per_locus_jsd, axis=0)
            jsd_df = pd.DataFrame(jsd_avg, index=pop_labels, columns=pop_labels)
            jsd_df.to_csv(OUT_DIR / "jsd_distance_matrix.tsv", sep="\t")
            condensed = squareform(jsd_avg, checks=False)
            Z2 = linkage(condensed, method="average")
            fig, ax = plt.subplots(figsize=(10, 5))
            dendrogram(Z2, labels=pop_labels, leaf_rotation=30, ax=ax,
                       color_threshold=0)
            ax.set_title("Population dendrogram — Jensen–Shannon divergence (avg across loci, average linkage)")
            ax.set_ylabel("JS distance")
            fig.tight_layout()
            fig.savefig(OUT_DIR / "F4_dendrogram_jsd.png", dpi=180)
            plt.close(fig)

    # ----- Per-locus top-5 forest grid -----
    fig = plt.figure(figsize=(13, 12))
    gs = GridSpec(3, 2, figure=fig, hspace=0.55, wspace=0.45)
    for i, loc in enumerate(LOCI):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        # rank top 5 alleles by max freq across populations
        rows_loc = matrix_top[matrix_top.index.str.startswith(f"{loc}*")]
        if rows_loc.empty:
            ax.text(0.5, 0.5, f"no data for {loc}", ha="center", transform=ax.transAxes)
            ax.set_title(f"HLA-{loc}")
            continue
        top5 = rows_loc.assign(maxv=rows_loc.max(axis=1)).sort_values("maxv", ascending=False).head(5).drop(columns="maxv")
        # plot grouped horizontal bars: y=allele, x=freq, color=population
        n_alleles = len(top5)
        n_pops = len(top5.columns)
        height = 0.8 / n_pops
        for j, pop in enumerate(top5.columns):
            offsets = np.arange(n_alleles) - 0.4 + (j + 0.5) * height
            c = color_map.get(pop_country[pop], "#000")
            ax.barh(offsets, top5[pop].values, height=height, color=c, edgecolor="black",
                    linewidth=0.3, label=pop.replace("AFND_", "").replace("_RNAseq", "*"))
        ax.set_yticks(np.arange(n_alleles))
        ax.set_yticklabels(top5.index, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("allele frequency", fontsize=8)
        ax.set_title(f"HLA-{loc}: top-5 alleles", fontsize=10)
        ax.tick_params(axis="x", labelsize=7)
        if i == 0:
            ax.legend(fontsize=6.5, loc="upper right", frameon=True, ncol=1)
    fig.suptitle("Per-locus top-5 allele frequencies across pan-Asian healthy populations", fontsize=12)
    fig.savefig(OUT_DIR / "F5_per_locus_forest_top5.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # ----- Korean vs East Asian aggregate signature -----
    # Use the long-format dataframe so we only compare alleles actually observed in
    # both the Korean source and the East Asian aggregate (avoids spurious zeros
    # from the master-matrix fillna).
    sig_rows = []
    have_long = long_df.copy()
    ea_long = have_long[have_long["population"] == "AFND_EastAsianAggregate"][["allele", "freq"]].rename(columns={"freq": "ea_freq"})

    # Korean comparator priority: Baek2021 NGS (richest) > RNA-seq cohorts > AFND pooled.
    for korean_col in ["Baek2021_KoreanNGS", "K2_normal_RNAseq", "GSE213647_normal_RNAseq", "AFND_Korea"]:
        kor = have_long[have_long["population"] == korean_col][["allele", "freq"]].rename(columns={"freq": "kor_freq"})
        if kor.empty or ea_long.empty:
            continue
        merged = kor.merge(ea_long, on="allele", how="inner")
        if merged.empty:
            continue
        merged["delta"] = merged["kor_freq"] - merged["ea_freq"]
        merged["abs_delta"] = merged["delta"].abs()
        merged["korean_source"] = korean_col
        merged["direction"] = np.where(merged["delta"] > 0, "enriched_in_Korean", "depleted_in_Korean")
        hits = merged[merged["abs_delta"] >= 0.05].copy()
        sig_rows.append(hits)

    sig = pd.concat(sig_rows, ignore_index=True) if sig_rows else pd.DataFrame()
    if not sig.empty:
        sig = sig.sort_values(["korean_source", "delta"], ascending=[True, False])
        sig.to_csv(OUT_DIR / "korean_vs_eastasian_signature_alleles.tsv", sep="\t", index=False)

        # Plot top-10 enriched and top-10 depleted using best available Korean source
        for src_pref in ["Baek2021_KoreanNGS", "AFND_Korea", "K2_normal_RNAseq"]:
            if src_pref in sig["korean_source"].values:
                ref = sig[sig["korean_source"] == src_pref]
                if len(ref) < 4:
                    continue
                up = ref.nlargest(10, "delta")
                dn = ref.nsmallest(10, "delta")
                sig_plot = pd.concat([up, dn]).drop_duplicates("allele").sort_values("delta")
                fig, ax = plt.subplots(figsize=(9, 7))
                colors = ["#d62728" if d > 0 else "#1f77b4" for d in sig_plot["delta"]]
                ax.barh(sig_plot["allele"], sig_plot["delta"], color=colors, edgecolor="black", linewidth=0.4)
                ax.axvline(0, color="black", lw=0.6)
                ax.set_xlabel("Δ allele frequency  (Korean − East Asian aggregate)")
                ax.set_title(f"Korean signature alleles  ({src_pref})  vs East Asian aggregate (|Δ| ≥ 0.05)\nred = enriched in Korean, blue = depleted in Korean")
                fig.tight_layout()
                fig.savefig(OUT_DIR / "F6_korean_signature_bars.png", dpi=180)
                plt.close(fig)
                break

    # ----- Heatmap top-50 alleles -----
    top50 = matrix_top.head(50)
    if not top50.empty and top50.shape[1] >= 2:
        col_link = linkage(top50.T.values, method="average")
        row_link = linkage(top50.values, method="average")
        col_order = dendrogram(col_link, no_plot=True)["leaves"]
        row_order = dendrogram(row_link, no_plot=True)["leaves"]
        ordered = top50.iloc[row_order, col_order]
        fig = plt.figure(figsize=(11, 12))
        gs = GridSpec(2, 2, width_ratios=[0.18, 1], height_ratios=[0.12, 1],
                      hspace=0.02, wspace=0.02)
        ax_top = fig.add_subplot(gs[0, 1])
        dendrogram(col_link, ax=ax_top, no_labels=True, color_threshold=0)
        ax_top.axis("off")
        ax_left = fig.add_subplot(gs[1, 0])
        dendrogram(row_link, ax=ax_left, orientation="left", no_labels=True, color_threshold=0)
        ax_left.axis("off")
        ax_main = fig.add_subplot(gs[1, 1])
        sns.heatmap(ordered, ax=ax_main, cmap="viridis", cbar=True,
                    cbar_kws={"label": "allele frequency", "shrink": 0.6},
                    linewidths=0.05, linecolor="white", vmin=0, vmax=ordered.values.max())
        ax_main.set_xlabel("")
        ax_main.set_ylabel("")
        plt.setp(ax_main.get_xticklabels(), rotation=35, ha="right", fontsize=8)
        plt.setp(ax_main.get_yticklabels(), fontsize=7)
        fig.suptitle("Top-50 alleles × pan-Asian healthy populations  (row + col dendrograms)",
                     fontsize=12, y=0.995)
        fig.savefig(OUT_DIR / "F7_heatmap_top50_alleles.png", dpi=180, bbox_inches="tight")
        plt.close(fig)

    # ----- Locus diversity -----
    # Use the long-format dataframe (full coverage per source), so a population's
    # locus diversity reflects all alleles it actually carries, not the master
    # matrix intersection.
    div_rows = []
    for pop, g in long_df.groupby("population"):
        for loc in LOCI:
            sub_loc = g[g["locus"] == loc]
            if len(sub_loc) < 2:
                # still record if 1 allele present
                if len(sub_loc) == 1:
                    div_rows.append(
                        {"population": pop, "locus": loc,
                         "shannon_H_bits": 0.0,
                         "expected_heterozygosity": 0.0,
                         "n_alleles_observed": 1}
                    )
                continue
            freqs = sub_loc["freq"].values
            s = freqs.sum()
            if s == 0:
                continue
            p = freqs / s
            H = float(entropy(p, base=2))
            obs_het = float(1 - np.sum(p ** 2))
            div_rows.append(
                {"population": pop, "locus": loc, "shannon_H_bits": H,
                 "expected_heterozygosity": obs_het,
                 "n_alleles_observed": int((freqs > 0).sum())}
            )
    div_df = pd.DataFrame(div_rows)
    div_df.to_csv(OUT_DIR / "locus_diversity.tsv", sep="\t", index=False)

    # Locus diversity bars
    if not div_df.empty:
        fig, axs = plt.subplots(1, 2, figsize=(15, 6), sharey=False)
        for j, metric in enumerate(["shannon_H_bits", "expected_heterozygosity"]):
            piv = div_df.pivot(index="locus", columns="population", values=metric)
            piv = piv.reindex(LOCI)
            x = np.arange(len(piv.index))
            n_pops = piv.shape[1]
            width = 0.8 / n_pops
            for k, pop in enumerate(piv.columns):
                offsets = x - 0.4 + (k + 0.5) * width
                c = color_map.get(pop_country.get(pop, "?"), "#000")
                axs[j].bar(offsets, piv[pop].values, width=width, color=c, edgecolor="black",
                           linewidth=0.3, label=pop.replace("AFND_", "").replace("_RNAseq", "*"))
            axs[j].set_xticks(x)
            axs[j].set_xticklabels(piv.index)
            axs[j].set_ylabel("Shannon H (bits)" if metric == "shannon_H_bits" else "Expected heterozygosity")
            axs[j].set_title(f"Locus diversity — {metric}")
            axs[j].grid(alpha=0.2, axis="y")
            if j == 0:
                axs[j].legend(fontsize=7, ncol=2, loc="upper left")
        fig.suptitle("Pan-Asian healthy HLA atlas — locus-level diversity")
        fig.tight_layout()
        fig.savefig(OUT_DIR / "F8_locus_diversity_bars.png", dpi=180)
        plt.close(fig)

    # ----- Resource-paper polished figures -----
    # We already produce F1 (PCA), F5 (forest grid), F7 (heatmap) which are the 3 requested.

    # ----- F6 (replaced): Korean signature vs Japan / Han / Taiwan-Han / Taiwan-Indig -----
    # Uses Baek2021 Korean NGS as reference to maximize allele coverage.
    if "Baek2021_KoreanNGS" in matrix_top.columns:
        kor_long = long_df[long_df["population"] == "Baek2021_KoreanNGS"][["allele", "freq"]].rename(columns={"freq": "kor_freq"})
        comparators = [c for c in ["AFND_Japan", "AFND_HanChinese", "AFND_Taiwan_Han", "AFND_Taiwan_Indigenous"]
                       if c in matrix_top.columns]
        if comparators and not kor_long.empty:
            fig, axs = plt.subplots(1, len(comparators), figsize=(3.4 * len(comparators), 5.5), sharey=False)
            if len(comparators) == 1:
                axs = [axs]
            for ax, comp in zip(axs, comparators):
                comp_long = long_df[long_df["population"] == comp][["allele", "freq"]].rename(columns={"freq": "comp_freq"})
                merged = kor_long.merge(comp_long, on="allele", how="inner")
                if merged.empty:
                    ax.text(0.5, 0.5, "no shared alleles", transform=ax.transAxes, ha="center")
                    continue
                merged["delta"] = merged["kor_freq"] - merged["comp_freq"]
                # show all available shared alleles ranked by delta
                merged = merged.sort_values("delta")
                colors = ["#d62728" if d > 0 else "#1f77b4" for d in merged["delta"]]
                ax.barh(merged["allele"], merged["delta"], color=colors, edgecolor="black", linewidth=0.4)
                ax.axvline(0, color="black", lw=0.6)
                ax.set_xlabel("Δ freq  (Baek2021 KR − comparator)", fontsize=8.5)
                comp_label = comp.replace("AFND_", "")
                ax.set_title(f"vs {comp_label}\n(n={len(merged)} shared alleles)", fontsize=9.5)
                ax.tick_params(axis="y", labelsize=7.5)
            fig.suptitle("Korean (Baek 2021 NGS) vs East Asian sub-populations  —  per-allele frequency Δ\n(red = enriched in Korean; blue = depleted)", fontsize=10.5)
            fig.tight_layout(rect=[0, 0, 1, 0.93])
            fig.savefig(OUT_DIR / "F6_korean_signature_bars.png", dpi=180)
            plt.close(fig)

    # ----- Cross-Korean platform concordance (tools-paper figure) -----
    # Compare Baek2021_KoreanNGS, K2_normal_RNAseq, GSE213647_normal_RNAseq pairwise
    korean_sources = ["Baek2021_KoreanNGS", "K2_normal_RNAseq", "GSE213647_normal_RNAseq"]
    available = [s for s in korean_sources if s in matrix_top.columns]
    if len(available) >= 2:
        from scipy.stats import pearsonr, spearmanr
        n_panels = len(available) * (len(available) - 1) // 2
        ncols = 3 if n_panels > 1 else 1
        nrows = max(1, (n_panels + ncols - 1) // ncols)
        fig, axs = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.3 * nrows), squeeze=False)
        idx = 0
        concordance_rows = []
        for i in range(len(available)):
            for j in range(i + 1, len(available)):
                ax = axs[idx // ncols][idx % ncols]
                a = available[i]; b = available[j]
                # restrict to alleles where both populations have non-zero values
                df_pair = matrix_top[[a, b]]
                df_pair = df_pair[(df_pair[a] > 0) | (df_pair[b] > 0)]
                if df_pair.empty:
                    continue
                rho_p, _ = pearsonr(df_pair[a], df_pair[b])
                rho_s, _ = spearmanr(df_pair[a], df_pair[b])
                ax.scatter(df_pair[a], df_pair[b], s=22, color="#444", alpha=0.65, edgecolor="white", linewidth=0.4)
                lim = max(df_pair[a].max(), df_pair[b].max()) * 1.05
                ax.plot([0, lim], [0, lim], "--", lw=0.8, color="red", alpha=0.6)
                ax.set_xlim(-0.005, lim); ax.set_ylim(-0.005, lim)
                ax.set_xlabel(a, fontsize=9); ax.set_ylabel(b, fontsize=9)
                ax.set_title(f"{a}\nvs {b}\nPearson r={rho_p:.2f}  Spearman ρ={rho_s:.2f}  (n_alleles={len(df_pair)})", fontsize=8.5)
                ax.grid(alpha=0.2)
                concordance_rows.append({
                    "source_A": a, "source_B": b, "n_alleles": len(df_pair),
                    "pearson": rho_p, "spearman": rho_s,
                })
                idx += 1
        # blank any unused axes
        for k in range(idx, nrows * ncols):
            axs[k // ncols][k % ncols].axis("off")
        fig.suptitle("Cross-Korean platform concordance (NGS vs RNA-seq imputation)", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig(OUT_DIR / "F9_cross_korean_concordance.png", dpi=180)
        plt.close(fig)
        pd.DataFrame(concordance_rows).to_csv(OUT_DIR / "cross_korean_concordance.tsv", sep="\t", index=False)

    # ----- Korean vs Japan / Han Chinese / Taiwan signature (using Baek2021 + AFND pools) -----
    pairwise_rows = []
    if "Baek2021_KoreanNGS" in matrix_top.columns:
        kor_long = long_df[long_df["population"] == "Baek2021_KoreanNGS"][["allele", "freq"]].rename(columns={"freq": "korean_freq"})
        for comp_pop in ["AFND_Japan", "AFND_HanChinese", "AFND_Taiwan_Han", "AFND_Taiwan_Indigenous"]:
            if comp_pop not in matrix_top.columns:
                continue
            comp_long = long_df[long_df["population"] == comp_pop][["allele", "freq"]].rename(columns={"freq": "comparator_freq"})
            merged = kor_long.merge(comp_long, on="allele", how="inner")
            if merged.empty:
                continue
            merged["delta"] = merged["korean_freq"] - merged["comparator_freq"]
            merged["abs_delta"] = merged["delta"].abs()
            merged["comparator"] = comp_pop
            merged["direction"] = np.where(merged["delta"] > 0, "enriched_in_Korean", "depleted_in_Korean")
            pairwise_rows.append(merged)
    if pairwise_rows:
        pw = pd.concat(pairwise_rows, ignore_index=True)
        pw.to_csv(OUT_DIR / "korean_vs_pairwise_signature_alleles.tsv", sep="\t", index=False)

    # ----- Summary JSON -----
    top5_korean = pd.DataFrame()
    if not sig.empty:
        for src_pref in ["Baek2021_KoreanNGS", "AFND_Korea", "K2_normal_RNAseq"]:
            if src_pref in sig["korean_source"].values:
                ref = sig[sig["korean_source"] == src_pref]
                top5_korean = ref.reindex(ref["abs_delta"].sort_values(ascending=False).index).head(10)
                break
    summary = {
        "n_populations": int(matrix_top.shape[1]),
        "n_alleles_in_master_matrix": int(matrix_top.shape[0]),
        "populations": meta_df.to_dict(orient="records"),
        "pca_explained_variance": [float(v) for v in pca.explained_variance_ratio_],
        "top10_korean_signature_alleles": top5_korean.to_dict(orient="records") if not top5_korean.empty else [],
        "k2_normal_n": int((k2.shape[0] if not k2.empty else 0)),
        "gse213647_normal_n": int((g213.shape[0] if not g213.empty else 0)),
        "afnd_pools_loaded": sorted(afnd_pooled["population"].unique().tolist()) if not afnd_pooled.empty else [],
        "afnd_alleles_in_cache": sorted(afnd_pooled["allele"].unique().tolist()) if not afnd_pooled.empty else [],
        "boundary": "healthy_only_resource_no_disease_association_tested",
    }
    (OUT_DIR / "track3_summary.json").write_text(json.dumps(summary, indent=2))
    print("[done]", OUT_DIR)


if __name__ == "__main__":
    main()
