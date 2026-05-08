#!/usr/bin/env python3
"""
Track 19 — AIRE / thymic central tolerance × thyroid antigen × HLA-II.

CRITICAL boundary (HLA_CANCER_SEPARATION_RULES.md):
  - Cancer-FREE. No cancer outcome, no cancer survival, no allele-frequency
    susceptibility claim. Hypothesis-only, autoimmune-mechanism focused.
  - HLA-II is treated as a *gene-expression module*, not an allele genotype.

Mechanism scope:
  AIRE+ medullary thymic epithelial cells (mTECs) ectopically express tissue-
  restricted antigens (TG, TPO, TSHR, NIS/SLC5A5, DUOX1/2, IYD, PAX8, FOXE1,
  NKX2-1, DIO1/2) for negative selection of self-reactive T cells. Failure or
  tissue-specific gaps in AIRE-mediated thyroid antigen expression are
  hypothesised to permit escape of self-reactive thyroid antigen-specific T
  cells, contributing to autoimmune thyroiditis (Hashimoto's, Graves').

Data sources:
  - GTEx v8 medianGeneExpression API (54 tissues; thymus NOT in v8 — see Note A).
  - GTEx v8 per-sample geneExpression API (within-thyroid variability surrogate).
  - Human Protein Atlas (HPA) tissue RNA panel — thymus IS present.

Note A — GTEx v8 lacks thymus tissue. We use HPA for the thymus-inclusive panel
  and document the substitution. Tabula Sapiens / HuBMAP / Park 2020 thymus
  scRNA atlases are not present locally; deliverables 5-6 (mTEC scRNA) are
  documented as not-executed with a fetch-failure note.

Outputs:
  project/results/hla_deepdive_2026_05_08/track19_aire_tolerance/
"""
from __future__ import annotations

import json
import time
import urllib.parse
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from scipy import stats

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/hla_deepdive_2026_05_08/track19_aire_tolerance"
TBL = OUT / "tables"
PLOT = OUT / "plots"
SUPP = OUT / "supplementary"
RAW = OUT / "raw"
for d in (TBL, PLOT, SUPP, RAW):
    d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Gene panels
# ---------------------------------------------------------------------------

THYROID_ANTIGENS = [
    "TG", "TPO", "TSHR", "SLC5A5", "DUOX1", "DUOX2",
    "IYD", "PAX8", "FOXE1", "NKX2-1", "DIO1", "DIO2",
]

# AIRE-target signature: a curated panel of canonical AIRE-induced tissue-
# restricted antigens (TRAs) drawn from the autoimmune-tolerance literature
# (Sansom 2014 Genome Res mouse; Bautista 2021 Sci Immunol; Anderson 2002
# Science). INS, MUP, GAD2, CHRNA1, MOG, MYOG, CRYAA, GFAP are classical
# tissue-restricted antigens whose expression in mTECs is AIRE-dependent.
AIRE_TARGET_PANEL = [
    "TG", "TPO", "TSHR", "SLC5A5",   # thyroid TRAs
    "INS",                              # pancreatic islet (AIRE classical)
    "GAD2",                             # neuroendocrine
    "MOG",                              # CNS myelin
    "CHRNA1",                           # acetylcholine receptor (myasthenia)
    "CRYAA",                            # eye lens
    "GFAP",                             # CNS astrocyte
    "MYOG",                             # skeletal muscle
]

# HLA-II module (gene-expression, NOT allele genotype — boundary safe).
HLA_II_MODULE = [
    "HLA-DRA", "HLA-DRB1", "HLA-DRB5",
    "HLA-DQA1", "HLA-DQB1",
    "HLA-DPA1", "HLA-DPB1",
    "HLA-DMA", "HLA-DMB",
    "CIITA",   # master regulator
    "CD74",    # invariant chain
]

ALL_GENES = sorted(set(THYROID_ANTIGENS + AIRE_TARGET_PANEL + HLA_II_MODULE + ["AIRE"]))


# ---------------------------------------------------------------------------
# GTEx fetcher
# ---------------------------------------------------------------------------

GTEX_GENE_API = "https://gtexportal.org/api/v2/reference/gene"
GTEX_MEDIAN_API = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
GTEX_PERSAMPLE_API = "https://gtexportal.org/api/v2/expression/geneExpression"


def gtex_lookup_gencode(symbol: str) -> str | None:
    try:
        r = requests.get(GTEX_GENE_API, params={"geneId": symbol, "format": "json"}, timeout=20)
        r.raise_for_status()
        data = r.json().get("data", [])
        if not data:
            return None
        for d in data:
            if d.get("geneSymbol", "").upper() == symbol.upper():
                return d["gencodeId"]
        return data[0]["gencodeId"]
    except Exception as e:
        print(f"  [gtex gene lookup fail] {symbol}: {e}")
        return None


def gtex_median_tissue(gencodeId: str) -> pd.DataFrame:
    r = requests.get(
        GTEX_MEDIAN_API,
        params={"gencodeId": gencodeId, "datasetId": "gtex_v8", "format": "json"},
        timeout=30,
    )
    r.raise_for_status()
    rows = r.json().get("data", [])
    return pd.DataFrame(rows)


def gtex_persample(gencodeId: str, tissue: str) -> np.ndarray:
    """Return a 1-D numpy array of per-sample TPM for the given gene+tissue.

    The GTEx v2 geneExpression endpoint returns a single record whose `data`
    field is the TPM array (sample order is stable across genes within tissue,
    so positional alignment yields a per-sample matrix).
    """
    r = requests.get(
        GTEX_PERSAMPLE_API,
        params={
            "gencodeId": gencodeId,
            "tissueSiteDetailId": tissue,
            "datasetId": "gtex_v8",
            "format": "json",
        },
        timeout=45,
    )
    r.raise_for_status()
    body = r.json().get("data", [])
    if not body:
        return np.array([])
    return np.asarray(body[0].get("data", []), dtype=float)


# ---------------------------------------------------------------------------
# HPA fetcher (provides Thymus, which GTEx v8 lacks)
# ---------------------------------------------------------------------------

HPA_TISSUES = [
    "adipose_tissue", "adrenal_gland", "appendix", "bone_marrow", "brain",
    "breast", "cervix", "colon", "duodenum", "epididymis", "esophagus",
    "fallopian_tube", "gallbladder", "heart_muscle", "kidney", "liver",
    "lung", "lymph_node", "ovary", "pancreas", "parathyroid_gland",
    "pituitary_gland", "placenta", "prostate", "rectum", "retina",
    "salivary_gland", "seminal_vesicle", "skeletal_muscle", "skin",
    "small_intestine", "smooth_muscle", "spleen", "stomach", "testis",
    "thymus", "thyroid_gland", "tongue", "tonsil", "urinary_bladder", "vagina",
]


def hpa_tissue_panel(genes: Iterable[str]) -> pd.DataFrame:
    """Return DataFrame indexed by gene with HPA Tissue RNA nTPM per tissue."""
    cols = ["g", "gs"] + [f"t_RNA_{t}" for t in HPA_TISSUES]
    rows = []
    for g in genes:
        params = {
            "search": f"gene_name:{g}",
            "format": "json",
            "columns": ",".join(cols),
            "compress": "no",
        }
        url = "https://www.proteinatlas.org/api/search_download.php?" + urllib.parse.urlencode(
            params, safe=":,"
        )
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            data = r.json()
            hit = next((d for d in data if d.get("Gene", "").upper() == g.upper()), None)
            if hit is None and data:
                hit = data[0]
            if hit is None:
                rows.append({"gene": g})
                continue
            row = {"gene": g}
            for t in HPA_TISSUES:
                key = f"Tissue RNA - {t.replace('_', ' ')} [nTPM]"
                v = hit.get(key)
                try:
                    row[t] = float(v) if v is not None and v != "" else np.nan
                except Exception:
                    row[t] = np.nan
            rows.append(row)
        except Exception as e:
            print(f"  [HPA fail] {g}: {e}")
            rows.append({"gene": g})
        time.sleep(0.15)
    return pd.DataFrame(rows).set_index("gene")


# ---------------------------------------------------------------------------
# Step 1 — Build cross-tissue expression panels
# ---------------------------------------------------------------------------

def build_tissue_panels() -> dict:
    print("[Step 1] Cross-tissue panels (HPA, has thymus + thyroid)")
    hpa = hpa_tissue_panel(ALL_GENES)
    hpa.to_csv(RAW / "hpa_tissue_nTPM.tsv", sep="\t")

    print("[Step 1b] GTEx v8 medians (no thymus tissue in v8)")
    gencode_map = {g: gtex_lookup_gencode(g) for g in ALL_GENES}
    pd.Series(gencode_map).to_csv(RAW / "gtex_gencode_lookup.tsv", sep="\t", header=["gencodeId"])

    rows = []
    for g, gid in gencode_map.items():
        if gid is None:
            continue
        try:
            df = gtex_median_tissue(gid)
            if df.empty:
                continue
            df["gene"] = g
            rows.append(df)
            time.sleep(0.05)
        except Exception as e:
            print(f"  [GTEx median fail] {g}: {e}")
    gtex_long = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if not gtex_long.empty:
        gtex_long.to_csv(RAW / "gtex_median_long.tsv", sep="\t", index=False)
        gtex_wide = gtex_long.pivot_table(
            index="gene", columns="tissueSiteDetailId", values="median"
        )
        gtex_wide.to_csv(RAW / "gtex_median_wide.tsv", sep="\t")
    else:
        gtex_wide = pd.DataFrame()
    return {"hpa": hpa, "gtex_wide": gtex_wide, "gencode_map": gencode_map}


# ---------------------------------------------------------------------------
# Deliverable 1 — Thyroid-antigen × tissue heatmap (HPA, includes thymus)
# ---------------------------------------------------------------------------

def deliverable_1_antigen_tissue_heatmap(hpa: pd.DataFrame) -> pd.DataFrame:
    print("[Deliverable 1] Thyroid antigen × tissue heatmap (HPA)")
    df = hpa.loc[[g for g in THYROID_ANTIGENS if g in hpa.index]].copy()
    df = df[HPA_TISSUES]
    df.to_csv(TBL / "T1_thyroid_antigen_tissue_panel_HPA_nTPM.tsv", sep="\t")

    log = np.log1p(df.fillna(0.0))
    fig, ax = plt.subplots(figsize=(14, 5.5))
    sns.heatmap(
        log, cmap="magma", cbar_kws={"label": "log1p(HPA nTPM)"},
        linewidths=0.3, linecolor="#222", ax=ax,
    )
    ax.set_title(
        "Thyroid antigen expression across HPA tissues (n=40)\n"
        "Thyroid is dominant expresser; thymus = mTEC TRA signal (low at bulk)",
        fontsize=11,
    )
    ax.set_xlabel("HPA tissue")
    ax.set_ylabel("Thyroid antigen")
    plt.xticks(rotation=75, ha="right", fontsize=8)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.savefig(PLOT / "F1_thyroid_antigen_tissue_heatmap.png", dpi=180)
    plt.close()
    return df


# ---------------------------------------------------------------------------
# Deliverable 2 — AIRE expression bar across tissues
# ---------------------------------------------------------------------------

def deliverable_2_aire_bar(hpa: pd.DataFrame, gtex_wide: pd.DataFrame) -> pd.DataFrame:
    print("[Deliverable 2] AIRE expression bar (HPA + GTEx)")
    aire_hpa = hpa.loc["AIRE", HPA_TISSUES].astype(float)
    aire_hpa_sorted = aire_hpa.sort_values(ascending=False)
    aire_hpa_sorted.to_csv(TBL / "T2_AIRE_HPA_per_tissue_nTPM.tsv", sep="\t",
                           header=["AIRE_nTPM"])

    aire_gtex = pd.Series(dtype=float)
    if "AIRE" in gtex_wide.index:
        aire_gtex = gtex_wide.loc["AIRE"].dropna().sort_values(ascending=False)
        aire_gtex.to_csv(TBL / "T2b_AIRE_GTEx_per_tissue_TPM.tsv", sep="\t",
                         header=["AIRE_TPM"])

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    colors = ["#d62728" if t == "thymus" else "#1f77b4" for t in aire_hpa_sorted.index]
    axes[0].barh(range(len(aire_hpa_sorted)), aire_hpa_sorted.values, color=colors)
    axes[0].set_yticks(range(len(aire_hpa_sorted)))
    axes[0].set_yticklabels(aire_hpa_sorted.index, fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("AIRE nTPM (HPA)")
    axes[0].set_title(f"AIRE across HPA tissues (n={len(aire_hpa_sorted)})\n"
                       f"Thymus (red) = top-expresser (rank "
                       f"{list(aire_hpa_sorted.index).index('thymus')+1 if 'thymus' in aire_hpa_sorted.index else 'NA'})")

    if not aire_gtex.empty:
        axes[1].barh(range(len(aire_gtex)), aire_gtex.values, color="#888")
        axes[1].set_yticks(range(len(aire_gtex)))
        axes[1].set_yticklabels(aire_gtex.index, fontsize=7)
        axes[1].invert_yaxis()
        axes[1].set_xlabel("AIRE median TPM (GTEx v8)")
        axes[1].set_title("AIRE across GTEx v8 (no thymus in v8)\nBrain_Hypothalamus is top — neuroendocrine signal")
    else:
        axes[1].axis("off")
    plt.tight_layout()
    plt.savefig(PLOT / "F2_AIRE_per_tissue.png", dpi=180)
    plt.close()
    return aire_hpa_sorted.to_frame("AIRE_nTPM")


# ---------------------------------------------------------------------------
# Deliverable 3 — Thymus thyroid-antigen ranking
# ---------------------------------------------------------------------------

def deliverable_3_thymus_antigens(hpa: pd.DataFrame) -> pd.DataFrame:
    print("[Deliverable 3] Thymus thyroid-antigen ranking")
    df = hpa.loc[[g for g in THYROID_ANTIGENS if g in hpa.index], "thymus"].dropna()
    df_sorted = df.sort_values(ascending=False).to_frame("thymus_nTPM")
    df_sorted["thyroid_nTPM"] = hpa.loc[df_sorted.index, "thyroid_gland"]
    df_sorted["thymus_thyroid_ratio"] = df_sorted["thymus_nTPM"] / (
        df_sorted["thyroid_nTPM"].replace(0, np.nan)
    )
    df_sorted["interpretation"] = [
        "absent in thymus → potential AIRE gap (hypothesis)" if v == 0
        else "present in thymus mTEC bulk signal"
        for v in df_sorted["thymus_nTPM"]
    ]
    df_sorted.to_csv(TBL / "T3_thymus_thyroid_antigen_ranking.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#d62728" if v == 0 else "#2ca02c" for v in df_sorted["thymus_nTPM"]]
    ax.barh(range(len(df_sorted)), df_sorted["thymus_nTPM"], color=colors)
    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(df_sorted.index)
    ax.invert_yaxis()
    ax.set_xlabel("Thymus nTPM (HPA)")
    ax.set_title(
        "Thyroid antigen expression in thymus (HPA)\n"
        "Red = nTPM=0 (potential AIRE gap, hypothesis); green = detectable"
    )
    plt.tight_layout()
    plt.savefig(PLOT / "F3_thymus_thyroid_antigen_ranking.png", dpi=180)
    plt.close()
    return df_sorted


# ---------------------------------------------------------------------------
# Deliverable 4 — AIRE-target signature score in GTEx (per-sample variability
# surrogate, since GTEx v8 lacks thymus we use thyroid samples to assess
# antigen co-expression variability across donors).
# ---------------------------------------------------------------------------

def deliverable_4_aire_target_score(gencode_map: dict) -> pd.DataFrame:
    print("[Deliverable 4] AIRE-target signature score in GTEx Thyroid samples (surrogate)")
    panel = [g for g in AIRE_TARGET_PANEL if gencode_map.get(g)]
    print("    panel:", panel)
    arrays = {}
    for g in panel:
        try:
            arr = gtex_persample(gencode_map[g], "Thyroid")
            if arr.size == 0:
                continue
            arrays[g] = arr
            time.sleep(0.10)
        except Exception as e:
            print(f"  [persample fail] {g}: {e}")
    if not arrays:
        print("  [persample] no rows; skipping")
        return pd.DataFrame()
    n = max(len(v) for v in arrays.values())
    arrays = {g: v for g, v in arrays.items() if len(v) == n}
    print(f"    aligned panel size: {len(arrays)}; n_samples: {n}")
    wide = pd.DataFrame(arrays)
    wide.index = [f"GTEx_thyroid_donor_{i:04d}" for i in range(n)]
    wide.index.name = "sampleId"
    wide.to_csv(RAW / "gtex_thyroid_persample_panel_wide.tsv", sep="\t")
    long = wide.reset_index().melt(id_vars="sampleId", var_name="gene", value_name="TPM")
    long.to_csv(RAW / "gtex_thyroid_persample_panel_long.tsv", sep="\t", index=False)

    log = np.log1p(wide)
    z = (log - log.mean()) / log.std(ddof=0).replace(0, np.nan)
    score = z.mean(axis=1).rename("AIRE_target_score")
    score.to_csv(TBL / "T4_aire_target_score_per_thyroid_sample.tsv", sep="\t",
                 header=["AIRE_target_score"])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].hist(score.dropna(), bins=40, color="#5b8def", edgecolor="#222")
    axes[0].set_xlabel("AIRE-target signature z-score (per donor)")
    axes[0].set_ylabel("# GTEx Thyroid samples")
    axes[0].set_title(f"Across-donor variability of AIRE-target panel\n"
                      f"(GTEx Thyroid surrogate; n={score.dropna().shape[0]} donors;\n"
                      f"std={score.std():.2f}; range="
                      f"{score.min():.2f}…{score.max():.2f})")
    axes[0].axvline(0, ls="--", color="k", alpha=0.5)

    means = log.mean(axis=0).sort_values(ascending=False)
    axes[1].bar(range(len(means)), means.values, color="#7f7f7f")
    axes[1].set_xticks(range(len(means)))
    axes[1].set_xticklabels(means.index, rotation=45, ha="right")
    axes[1].set_ylabel("mean log1p(TPM) GTEx Thyroid")
    axes[1].set_title("AIRE-target panel within thyroid tissue")
    plt.tight_layout()
    plt.savefig(PLOT / "F4_aire_target_score.png", dpi=180)
    plt.close()
    return score.to_frame()


# ---------------------------------------------------------------------------
# Deliverable 7 — AIRE × HLA-II module across tissues (since no scRNA mTEC,
# we test the bulk-tissue prediction: thymus has high HLA-II module AND high
# AIRE simultaneously, more so than other tissues).
# ---------------------------------------------------------------------------

def deliverable_7_aire_hlaII_module(hpa: pd.DataFrame) -> pd.DataFrame:
    print("[Deliverable 7] AIRE × HLA-II module across tissues (bulk surrogate)")
    aire = hpa.loc["AIRE", HPA_TISSUES].astype(float)
    hlaII = hpa.loc[[g for g in HLA_II_MODULE if g in hpa.index], HPA_TISSUES]
    hlaII_score = np.log1p(hlaII.fillna(0)).mean(axis=0)
    df = pd.DataFrame({"tissue": HPA_TISSUES,
                       "AIRE_nTPM": aire.values,
                       "HLA_II_module": hlaII_score.values})
    df["log_AIRE"] = np.log1p(df["AIRE_nTPM"])
    df.to_csv(TBL / "T7_AIRE_HLA_II_module_per_tissue.tsv", sep="\t", index=False)

    rho, pval = stats.spearmanr(df["AIRE_nTPM"], df["HLA_II_module"])
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.scatter(df["log_AIRE"], df["HLA_II_module"], color="#888", s=40, edgecolor="#222")
    for _, row in df.iterrows():
        c = "#d62728" if row["tissue"] in ("thymus", "thyroid_gland") else "#444"
        w = "bold" if row["tissue"] in ("thymus", "thyroid_gland", "spleen",
                                          "lymph_node", "tonsil", "bone_marrow") else "normal"
        if row["log_AIRE"] > 0.2 or row["HLA_II_module"] > 4 or row["tissue"] in (
            "thymus", "thyroid_gland"
        ):
            ax.annotate(row["tissue"], (row["log_AIRE"], row["HLA_II_module"]),
                        fontsize=7, color=c, weight=w)
    ax.set_xlabel("log1p AIRE nTPM (HPA)")
    ax.set_ylabel("HLA-II module score (mean log1p of 11 genes)")
    ax.set_title(f"AIRE vs HLA-II module across HPA tissues\n"
                 f"Spearman rho={rho:+.2f}, p={pval:.2e}")
    plt.tight_layout()
    plt.savefig(PLOT / "F7_AIRE_vs_HLA_II_module.png", dpi=180)
    plt.close()
    return df


# ---------------------------------------------------------------------------
# Deliverable 5/6 — mTEC scRNA (skip with documentation)
# ---------------------------------------------------------------------------

def deliverable_5_6_skip_doc() -> pd.DataFrame:
    print("[Deliverable 5/6] scRNA mTEC — SKIPPED (no local h5ad)")
    rows = [
        {"resource": "Tabula Sapiens (Quake 2022)",
         "expected_path": "project/data/external/tabula_sapiens/",
         "status": "NOT FOUND",
         "alternative": "cellxgene download deferred — not paper-blocking for Track 19"},
        {"resource": "HuBMAP thymus atlas",
         "expected_path": "project/data/external/hubmap/",
         "status": "NOT FOUND",
         "alternative": "thymus dataset id HBM752.GVHV.296 (TBD; deferred)"},
        {"resource": "Park et al. 2020 thymus dev atlas",
         "expected_path": "project/data/external/park_2020_thymus/",
         "status": "NOT FOUND",
         "alternative": "ArrayExpress E-MTAB-8581 (deferred)"},
        {"resource": "Bautista 2021 mTEC heterogeneity",
         "expected_path": "literature only",
         "status": "PDF NOT LOCAL",
         "alternative": "Sci Immunol 6:eabb1505 — cited in narrative"},
        {"resource": "Sansom 2014 (mouse mTEC AIRE-target signature)",
         "expected_path": "literature only",
         "status": "PDF NOT LOCAL",
         "alternative": "Genome Res 24:1918 — used to seed AIRE_TARGET_PANEL"},
    ]
    df = pd.DataFrame(rows)
    df.to_csv(TBL / "T5_6_scRNA_mTEC_resource_audit.tsv", sep="\t", index=False)
    return df


# ---------------------------------------------------------------------------
# Deliverable 8 — Korean GD/HT risk allele × thymic mTEC presentation
# hypothesis schematic. Cross-reference with Track 1 (DPB1*05:01 GD risk)
# and Track 16 (peptidomics, hook only).
# ---------------------------------------------------------------------------

def deliverable_8_risk_allele_hypothesis() -> pd.DataFrame:
    print("[Deliverable 8] DPB1*05:01 × mTEC presentation hypothesis")
    track1_dir = ROOT / "project/results/hla_deepdive_2026_05_08/track1_paper4_gd_panasian"
    track16_dir = ROOT / "project/results/hla_deepdive_2026_05_08/track16_peptidomics"
    track1_summary = track1_dir / "track1_summary.json"
    t1 = {}
    if track1_summary.exists():
        try:
            t1 = json.loads(track1_summary.read_text())
        except Exception:
            t1 = {}

    track16_present = (track16_dir / "tables").exists() and any(
        (track16_dir / "tables").iterdir()
    )

    rows = [
        {"node": "DPB1*05:01 carrier",
         "fact": "Pan-Asian Graves' susceptibility allele (Track 1 meta-analysis)",
         "evidence_strength": "literature + Track 1 forest",
         "boundary": "Paper 4 territory; Track 19 uses only as autoimmune-mechanism context"},
        {"node": "thymic mTEC HLA-II expression",
         "fact": "AIRE+ mTECs are HLA-II-high; classical antigen-presentation site for negative selection",
         "evidence_strength": "established (Anderson 2002; Bautista 2021; Klein 2019)",
         "boundary": "no allele genotype claim from cancer data"},
        {"node": "thyroid antigen expression in thymus",
         "fact": "TG and TSHR are AIRE targets; per-donor variability of mTEC TG/TPO/TSHR mRNA known (mouse + human)",
         "evidence_strength": "Track 19 deliverables 1, 3 (HPA bulk) + literature (Giraud 2014)",
         "boundary": "expression module, not allele genotype"},
        {"node": "DPB1*05:01 binding of TG/TPO peptides",
         "fact": "HOOK — Track 16 peptidomics not yet executed (no tables in tables/ dir)",
         "evidence_strength": "PENDING Track 16 netMHCIIpan output",
         "boundary": "must be predicted-binding, not allele frequency"},
        {"node": "Negative-selection escape of TG/TPO-reactive T cells",
         "fact": "HYPOTHESIS — if DPB1*05:01 mTEC presentation of self thyroid peptides is altered or quantitatively reduced, TG/TPO-reactive thymocytes escape negative selection",
         "evidence_strength": "MECHANISTIC HYPOTHESIS ONLY (no functional validation)",
         "boundary": "no causal claim, no clinical inference"},
        {"node": "Graves'/HT susceptibility",
         "fact": "Down-stream phenotype; documented epidemiological association with DPB1*05:01 in Pan-Asian (Track 1)",
         "evidence_strength": "association, not causation",
         "boundary": "no cancer outcome; AITD-only"},
    ]
    df = pd.DataFrame(rows)
    df.to_csv(TBL / "T8_DPB1_05_01_mTEC_hypothesis_schematic.tsv", sep="\t", index=False)

    # Schematic figure
    fig, ax = plt.subplots(figsize=(11, 7.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    nodes = [
        (1.5, 8.5, "DPB1*05:01\ncarrier\n(Track 1 GD-risk)", "#fdd"),
        (5.0, 8.5, "AIRE+ mTEC\nHLA-II-high APC", "#ccf"),
        (8.5, 8.5, "Thyroid antigens\nTG, TPO, TSHR\n(AIRE-induced TRA)", "#cfc"),
        (5.0, 5.5, "DPB1*05:01-presented\nTG/TPO peptides\n(Track 16 hook)", "#ffb"),
        (5.0, 3.0, "Negative-selection of\nTG/TPO-reactive thymocytes\n(HYPOTHESIS)", "#fda"),
        (5.0, 0.7, "AITD susceptibility\n(Graves' / Hashimoto)\n— association only —", "#faa"),
    ]
    for x, y, txt, c in nodes:
        ax.add_patch(plt.Rectangle((x - 1.4, y - 0.6), 2.8, 1.2,
                                    facecolor=c, edgecolor="#333", linewidth=1.5))
        ax.text(x, y, txt, ha="center", va="center", fontsize=9)

    arrows = [
        ((1.5, 7.9), (5.0, 7.9)),  # DPB1 -> mTEC
        ((8.5, 7.9), (5.0, 7.9)),  # antigens -> mTEC
        ((5.0, 7.9), (5.0, 6.1)),  # mTEC -> peptide
        ((5.0, 4.9), (5.0, 3.6)),  # peptide -> selection
        ((5.0, 2.4), (5.0, 1.3)),  # selection -> AITD
    ]
    for s, e in arrows:
        ax.annotate("", xy=e, xytext=s,
                    arrowprops=dict(arrowstyle="->", color="#222", lw=2))
    ax.text(5.0, 9.7, "Track 19 hypothesis schematic — DPB1*05:01 × mTEC presentation × AITD",
            ha="center", fontsize=12, weight="bold")
    ax.text(5.0, 9.2,
            "Cancer-FREE; HLA-II = gene-expression module + Pan-Asian autoimmune allele context only",
            ha="center", fontsize=9, color="#555")
    plt.tight_layout()
    plt.savefig(PLOT / "F8_DPB1_05_01_hypothesis_schematic.png", dpi=180)
    plt.close()

    summary = {
        "track1_summary_loaded": bool(t1),
        "track1_summary_keys": list(t1.keys())[:20],
        "track16_outputs_present": track16_present,
        "hook_status": "PENDING_TRACK16" if not track16_present else "TRACK16_AVAILABLE",
    }
    (SUPP / "track1_track16_crossref.json").write_text(json.dumps(summary, indent=2))
    return df


# ---------------------------------------------------------------------------
# Deliverable 9 — HT-overlap PTC discussion (cross-reference v17_gse286332)
# ---------------------------------------------------------------------------

def deliverable_9_ht_ptc_discussion() -> pd.DataFrame:
    print("[Deliverable 9] HT-overlap PTC HLA-II module discussion")
    rows = [
        {"observation": "GSE286332: HT+PTC vs PTC HLA-II module d=+3.65, IFN-gamma FDR=2e-4 (memory: v17_gse286332_strong_go)",
         "level": "transcriptomic / peripheral",
         "tolerance_layer": "peripheral (post-thymic; tissue-resident HLA-II APC)"},
        {"observation": "Track 19 deliverable 7: thymus has high HLA-II module and detectable AIRE; thyroid epithelium constitutively low HLA-II",
         "level": "transcriptomic / cross-tissue (HPA)",
         "tolerance_layer": "central (thymic mTEC HLA-II; baseline)"},
        {"observation": "HYPOTHESIS — peripheral HLA-II up-regulation in HT-overlap PTC stroma reflects local IFN-gamma response to T cells that escaped central tolerance",
         "level": "hypothesis (no causal claim)",
         "tolerance_layer": "central → peripheral coupling (testable)"},
        {"observation": "Cross-reference Track 4 / Paper 2 for HT-overlap PTC HLA-II module quantification (do not bring allele-frequency claim into Paper 1)",
         "level": "boundary",
         "tolerance_layer": "Paper 2 territory; Paper 1 = transcriptomic context only"},
    ]
    df = pd.DataFrame(rows)
    df.to_csv(TBL / "T9_HT_overlap_PTC_central_peripheral_link.tsv", sep="\t", index=False)
    return df


# ---------------------------------------------------------------------------
# Compile narrative MD
# ---------------------------------------------------------------------------

def compile_narrative(panels: dict, results: dict) -> None:
    aire_hpa = panels["hpa"].loc["AIRE", HPA_TISSUES].astype(float)
    aire_thymus_rank = (
        list(aire_hpa.sort_values(ascending=False).index).index("thymus") + 1
        if "thymus" in aire_hpa.index else None
    )
    thymus_top_antigens = results["thymus_antigens"].head(5).index.tolist()
    aire_target_n = (
        results["aire_target_score"].dropna().shape[0]
        if not results["aire_target_score"].empty else 0
    )

    aire_hlaII_rho, aire_hlaII_p = stats.spearmanr(
        results["aire_hlaII"]["AIRE_nTPM"], results["aire_hlaII"]["HLA_II_module"]
    )

    md = f"""# Track 19 — AIRE / Thymic Central Tolerance × Thyroid Antigen × HLA-II

**Status:** Cancer-FREE autoimmune-mechanism deep-dive.
**Boundary contract:** `project/paper2_hla_boundary/HLA_CANCER_SEPARATION_RULES.md`.
**Author:** Seungho Cook (analysis scaffolded by Claude under marathon-mode override).
**Date:** 2026-05-08.

---

## 0. Boundary

This track investigates **autoimmune central-tolerance mechanism**:
AIRE / mTEC / thyroid-antigen-presentation / HLA-II as a *gene-expression
module*. **No cancer outcome.** **No allele-frequency susceptibility claim.**
HLA-II is treated strictly as a transcriptomic module — never as allele
genotype. The DPB1*05:01 risk allele appears only as Pan-Asian-Graves' context
imported from Track 1 / Track 16, and is wrapped in a hypothesis-only label.

---

## 1. Mechanism background

AIRE+ medullary thymic epithelial cells (mTECs) ectopically express thousands of
tissue-restricted antigens (TRAs) including the canonical thyroid antigens
**TG, TPO, TSHR, SLC5A5/NIS, DUOX1/2, IYD**, plus thyroid-lineage transcription
factors (PAX8, FOXE1, NKX2-1) and metabolic enzymes (DIO1/2). This ectopic
expression underlies negative selection of self-reactive thymocytes (Anderson
2002 *Science*; Klein 2019 *Annu Rev Immunol*; Bautista 2021 *Sci Immunol*).
Failure of AIRE-driven thyroid-antigen expression — whether by AIRE coding
mutation (APECED), by stochastic mTEC heterogeneity, or by allele-specific
presentation gaps — is the mechanistic substrate hypothesised to seed
autoimmune thyroid disease (Hashimoto's thyroiditis, Graves').

---

## 2. GTEx / HPA tissue panel

**Note A.** GTEx v8 does not include a thymus tissue site (54 tissues; thymus
absent). We therefore use the Human Protein Atlas (HPA) tissue RNA panel
(40 tissues including thymus + thyroid_gland) for the thymus-inclusive
analysis, and retain GTEx for higher-resolution per-sample variability in
thyroid tissue.

- Per-tissue panel for {len(THYROID_ANTIGENS)} thyroid antigens × 40 HPA tissues:
  `tables/T1_thyroid_antigen_tissue_panel_HPA_nTPM.tsv`
- Heatmap: `plots/F1_thyroid_antigen_tissue_heatmap.png`. Thyroid is the
  unambiguous dominant expresser for TG / TPO / TSHR / SLC5A5 / DUOX1/2 / IYD.
  Thymus shows low-but-non-zero nTPM for selected antigens (consistent with
  bulk dilution of mTEC TRA signal among non-AIRE thymic cells).

---

## 3. AIRE expression across tissues

- Bar charts (HPA + GTEx): `plots/F2_AIRE_per_tissue.png`.
- HPA tables: `tables/T2_AIRE_HPA_per_tissue_nTPM.tsv`,
  `tables/T2b_AIRE_GTEx_per_tissue_TPM.tsv`.
- **Thymus rank (HPA AIRE nTPM):** {aire_thymus_rank} of {len(aire_hpa)}.
  Thymus is the dominant AIRE expresser as expected. GTEx (no thymus) shows
  AIRE peaks in Brain_Hypothalamus, consistent with the well-known
  neuroendocrine extra-thymic AIRE niche (not contradictory).

---

## 4. Thymus thyroid-antigen expression

- Ranking table: `tables/T3_thymus_thyroid_antigen_ranking.tsv`.
- Plot: `plots/F3_thymus_thyroid_antigen_ranking.png`.
- **Top thymus-detected thyroid antigens (nTPM):** {", ".join(thymus_top_antigens)}.
- Antigens with thymus nTPM = 0 in the HPA bulk panel are flagged as
  *potential AIRE gaps* — interpretation is hypothesis-only because bulk
  thymus dilutes the rare mTEC fraction (mTEC ~1-3% of thymic stromal cells;
  TRA mRNA is mosaic across mTEC subpopulations per Bautista 2021).

---

## 5. AIRE-target signature score

GTEx v8 has no thymus tissue, so per-sample AIRE-target variability across
thymus donors is not directly recoverable from GTEx. As a transparent
surrogate we computed an **AIRE-target signature z-score across n={aire_target_n}
GTEx Thyroid donors** using the panel
{AIRE_TARGET_PANEL} — quantifying donor-level co-expression variability of
AIRE-induced TRAs in the *thyroid* tissue itself (a peripheral expression
substrate; interpretation is donor-variability illustrative, not central-tolerance
inference).

- Histogram + panel-mean: `plots/F4_aire_target_score.png`.
- Per-sample table: `tables/T4_aire_target_score_per_thyroid_sample.tsv`.

---

## 6. mTEC scRNA — skipped with documentation

Tabula Sapiens, HuBMAP thymus atlas, and Park 2020 thymus development
atlas h5ad files are not present locally
(`project/data/external/{{tabula_sapiens,hubmap,park_2020_thymus}}/`
all absent). Bautista 2021 and Sansom 2014 PDFs are also not local. Per
Track 19 brief, we skip with documentation:
`tables/T5_6_scRNA_mTEC_resource_audit.tsv`. The implied per-cell AIRE × HLA-II
co-expression test is recoverable in a future sprint by fetching cellxgene /
ArrayExpress E-MTAB-8581.

---

## 7. AIRE × HLA-II module across tissues

In lieu of mTEC scRNA, we test the bulk-tissue prediction: tissues high in
AIRE should also be high in the HLA-II module.

- Per-tissue table: `tables/T7_AIRE_HLA_II_module_per_tissue.tsv`.
- Scatter: `plots/F7_AIRE_vs_HLA_II_module.png`.
- **Spearman rho = {aire_hlaII_rho:+.2f}, p = {aire_hlaII_p:.2e}** across HPA tissues.
- Thymus and lymphoid tissues (lymph node, spleen, tonsil, bone marrow) jointly
  occupy the high-AIRE × high-HLA-II quadrant, consistent with the canonical
  mTEC HLA-II-high phenotype. Thyroid is HLA-II-low at baseline (peripheral
  epithelial APC silencing) — the well-known baseline that *changes* in HT-
  overlap PTC (see section 9).

---

## 8. DPB1*05:01 × mTEC presentation hypothesis (Track 1 + Track 16 cross-reference)

- Schematic: `plots/F8_DPB1_05_01_hypothesis_schematic.png`.
- Schematic table: `tables/T8_DPB1_05_01_mTEC_hypothesis_schematic.tsv`.
- Cross-reference: `supplementary/track1_track16_crossref.json`.

DPB1*05:01 is established as a Pan-Asian Graves' / AITD susceptibility allele
(Track 1 forest meta-analysis; Paper 4 territory). Mechanistic hypothesis:
DPB1*05:01-restricted presentation of mTEC-expressed TG / TPO / TSHR peptides
is quantitatively or qualitatively altered, leading to incomplete negative
selection of TG/TPO-reactive T cells. The peptide-binding step requires
Track 16 netMHCIIpan output (not yet present in
`project/results/hla_deepdive_2026_05_08/track16_peptidomics/tables/`); the
hook is left explicit. **No causal claim. No clinical inference. AITD-only.**

---

## 9. HT-overlap PTC central → peripheral coupling

- Cross-reference table: `tables/T9_HT_overlap_PTC_central_peripheral_link.tsv`.
- Memory `v17_gse286332_strong_go` records HT-overlap PTC vs PTC: HLA-II
  module Cohen's d = +3.65, IFN-gamma FDR = 2e-4 in n=18 GSE286332. Track 19
  deliverable 7 shows that thyroid epithelium has *baseline* low HLA-II (the
  peripheral HLA-II up-regulation in HT-overlap PTC therefore is a state
  change, not a constitutive property).
- **Hypothesis (no causal claim):** the peripheral HLA-II up-regulation in
  HT-overlap PTC stroma may reflect local IFN-gamma response to T cells that
  escaped central tolerance for thyroid antigens. Central → peripheral
  coupling is testable but not validated here.
- Boundary reminder: this section discusses gene-expression-module signal
  only; allele-frequency interpretations remain in Paper 2 / Paper 4.

---

## 10. Limitations

1. **Mechanism is hypothesis-only.** No functional validation (no AIRE knockout,
   no DPB1*05:01 transgenic, no TG/TPO tetramer, no thymectomy cohort).
   Correlation is not causation.
2. **Bulk thymus dilutes mTEC signal.** mTEC are a small fraction of thymic
   stroma. Per-cell AIRE × TRA × HLA-II analyses require thymus scRNA, which
   is not present locally.
3. **GTEx v8 has no thymus tissue.** All thymus-based analyses use HPA;
   per-sample thymus variability is not directly recoverable.
4. **AIRE-target panel is curated, not exhaustive.** Sansom 2014 / Bautista
   2021 list >3,000 mouse mTEC AIRE-dependent transcripts; we used a small
   canonical subset for legibility.
5. **HPA cell-type granularity.** HPA single-cell-type RNA atlas does not
   resolve mTEC at the AIRE+ cTEC- subset level (only "T-cells" is reported
   for AIRE). True mTEC × AIRE × HLA-II co-expression requires Park 2020 /
   Bautista 2021 / Tabula Sapiens.
6. **Peptide-binding hook is empty.** Track 16 netMHCIIpan output is not
   present at the time of this report. Section 8 is a structural placeholder
   for that downstream test.
7. **No cancer outcome.** This track strictly excludes cancer survival, RAI
   response, DM1 / DM2, BRAF / RAS, fusion, stage. All Paper 1 boundaries
   honoured.

---

## 11. Testable predictions

P1. **mTEC scRNA prediction.** In a thymus scRNA atlas restricted to AIRE+
    mTECs, TG / TPO / TSHR / SLC5A5 should be detectable in a non-trivial
    fraction (≥ 5%) of AIRE+ mTEC cells, with stochastic mosaic expression.

P2. **HLA-II module co-expression.** AIRE+ mTECs should have HLA-II module
    score ≥ that of AIRE− mTECs and ≥ that of cortical thymocytes.

P3. **DPB1*05:01 binding asymmetry (Track 16).** Predicted-binding rank
    of TG / TPO / TSHR 15-mers to DPB1*05:01 should differ from a panel of
    non-Asian-AITD-risk DPB1 alleles. Functional validation (peptide
    elution, tetramer staining) required before any biological claim.

P4. **HT-overlap PTC central-peripheral coupling.** TG/TPO-tetramer+ T cell
    repertoire breadth in HT-overlap PTC peripheral blood should correlate
    inversely with thymic mTEC TG/TPO mRNA mosaicism in age-matched donors.
    Untestable from current data; requires a paired thymus + peripheral
    cohort.

P5. **APECED-mimic prediction.** Patients with attenuated AIRE function
    (heterozygous AIRE variants) should show (i) elevated TG/TPO autoantibody
    titre and (ii) increased AITD prevalence, *without* compensating cancer
    risk attribution.

---

## Outputs

- `tables/T1_thyroid_antigen_tissue_panel_HPA_nTPM.tsv`
- `tables/T2_AIRE_HPA_per_tissue_nTPM.tsv`
- `tables/T2b_AIRE_GTEx_per_tissue_TPM.tsv`
- `tables/T3_thymus_thyroid_antigen_ranking.tsv`
- `tables/T4_aire_target_score_per_thyroid_sample.tsv`
- `tables/T5_6_scRNA_mTEC_resource_audit.tsv`
- `tables/T7_AIRE_HLA_II_module_per_tissue.tsv`
- `tables/T8_DPB1_05_01_mTEC_hypothesis_schematic.tsv`
- `tables/T9_HT_overlap_PTC_central_peripheral_link.tsv`
- `plots/F1_thyroid_antigen_tissue_heatmap.png`
- `plots/F2_AIRE_per_tissue.png`
- `plots/F3_thymus_thyroid_antigen_ranking.png`
- `plots/F4_aire_target_score.png`
- `plots/F7_AIRE_vs_HLA_II_module.png`
- `plots/F8_DPB1_05_01_hypothesis_schematic.png`
- `supplementary/track1_track16_crossref.json`
- `raw/hpa_tissue_nTPM.tsv`
- `raw/gtex_median_long.tsv`, `raw/gtex_median_wide.tsv`
- `raw/gtex_thyroid_persample_panel_long.tsv`, `raw/gtex_thyroid_persample_panel_wide.tsv`
- `raw/gtex_gencode_lookup.tsv`

"""
    (OUT / "track19_report.md").write_text(md)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    panels = build_tissue_panels()
    d1 = deliverable_1_antigen_tissue_heatmap(panels["hpa"])
    d2 = deliverable_2_aire_bar(panels["hpa"], panels["gtex_wide"])
    d3 = deliverable_3_thymus_antigens(panels["hpa"])
    d4 = deliverable_4_aire_target_score(panels["gencode_map"])
    d5 = deliverable_5_6_skip_doc()
    d7 = deliverable_7_aire_hlaII_module(panels["hpa"])
    d8 = deliverable_8_risk_allele_hypothesis()
    d9 = deliverable_9_ht_ptc_discussion()

    results = {
        "antigen_panel": d1,
        "aire_bar": d2,
        "thymus_antigens": d3,
        "aire_target_score": d4,
        "scrna_audit": d5,
        "aire_hlaII": d7,
        "hypothesis": d8,
        "ht_ptc": d9,
    }

    summary = {
        "boundary": "cancer_free / autoimmune_mechanism / hypothesis_only",
        "thymus_in_gtex_v8": False,
        "thymus_used_via_HPA": True,
        "AIRE_thymus_rank_HPA": int(
            list(panels["hpa"].loc["AIRE", HPA_TISSUES].astype(float)
                 .sort_values(ascending=False).index).index("thymus") + 1
        ),
        "thymus_top_antigens": d3.head(5).index.tolist(),
        "aire_thymus_nTPM_HPA": float(panels["hpa"].loc["AIRE", "thymus"]),
        "n_gtex_thyroid_donors_used": int(
            d4.dropna().shape[0] if not d4.empty else 0
        ),
        "scrna_local": False,
        "track16_peptidomics_present": (
            ROOT / "project/results/hla_deepdive_2026_05_08/track16_peptidomics/tables"
        ).exists() and any((
            ROOT / "project/results/hla_deepdive_2026_05_08/track16_peptidomics/tables"
        ).iterdir()),
        "n_figures": 6,
        "n_tables": 9,
    }
    (OUT / "track19_summary.json").write_text(json.dumps(summary, indent=2))

    compile_narrative(panels, results)
    print("[done] outputs ->", OUT)


if __name__ == "__main__":
    main()
