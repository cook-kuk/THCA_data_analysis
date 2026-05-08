#!/usr/bin/env python3
"""Paper 11 — Phase D: DepMap CRISPR essentiality × DM1 axis.

For each cell line in DepMap (~1200), score DM1-axis using available
gene expression (we use the existing CCLE DM1 score from
paper9_sl_first_pass/ccle_dm1_state.tsv when available, else fall back
to gene-level expression-based proxy not implemented here).

Tests: are top hub genes (LYN, NAMPT, TACSTD2, KCNN4, PAX8, JAK1, JAK2,
STAT3) essential in DM1-high vs DM1-low cell lines? If yes, those are
candidate therapeutic targets.

Outputs:
  - phase_D_depmap/dm1_high_low_dependency.tsv
  - phase_D_depmap/lineage_dependency_summary.tsv
  - phase_D_depmap/summary.json
"""
from __future__ import annotations
import argparse, json, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

DM_RAW = Path("/data/thca/repo_results/p3_p9_full_execution/paper9/raw")
CCLE_STATE = Path("/data/thca/repo_results/paper9_sl_first_pass/ccle_dm1_state.tsv")
OUT = Path("/data/thca/repo_results/paper11_pancancer/phase_D_depmap")

HUB_GENES = ["LYN", "NAMPT", "TACSTD2", "KCNN4", "PAX8",
             "JAK1", "JAK2", "STAT3", "FOXE1", "NKX2-1", "HHEX",
             "FYN", "SRC", "DNMT1", "EZH2", "MYC"]
PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]


def load_crispr_subset(path, genes):
    """CRISPRGeneEffect.csv: row=ModelID, col='SYMBOL (entrez)' format."""
    # Read header to map symbol→col
    header = pd.read_csv(path, nrows=0).columns.tolist()
    sym2col = {}
    for c in header:
        m = re.match(r"^([A-Z0-9\-]+)\s*\(\d+\)$", c)
        if m and m.group(1) in genes:
            sym2col[m.group(1)] = c
    # First column is unnamed (model ID)
    first_col = header[0]
    cols = [first_col] + [sym2col[g] for g in genes if g in sym2col]
    df = pd.read_csv(path, usecols=cols)
    df = df.rename(columns={first_col: "ModelID"})
    rename = {v: k for k, v in sym2col.items()}
    df = df.rename(columns=rename)
    return df


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    crispr_path = DM_RAW / "CRISPRGeneEffect.csv"
    model_path = DM_RAW / "Model.csv"
    if not crispr_path.exists():
        sys.exit(f"missing {crispr_path}")
    if not CCLE_STATE.exists():
        sys.exit(f"missing CCLE DM1 state {CCLE_STATE}; run paper9_sl_first_pass first")

    targets = HUB_GENES + PANEL_8
    print(f"[D] loading CRISPR effect for {len(targets)} target genes …")
    crispr = load_crispr_subset(crispr_path, set(targets))
    print(f"[D] CRISPR rows: {len(crispr)}, cols hit: "
          f"{[c for c in crispr.columns if c != 'ModelID']}")

    model = pd.read_csv(model_path, low_memory=False)[
        ["ModelID", "StrippedCellLineName", "OncotreeLineage",
         "OncotreePrimaryDisease", "OncotreeCode"]]

    # Use lineage-level DM1 proxy: each TCGA pancancer lineage has a median DM1
    # score (paper11/pancan_lineage_stats.tsv). Map DepMap OncotreeLineage to
    # the closest pancancer lineage and inherit its median DM1.
    pancan_stats = pd.read_csv(
        "/data/thca/repo_results/paper11_pancancer/pancan_lineage_stats.tsv", sep="\t")
    LINEAGE_MAP = {
        # DepMap OncotreeLineage : pancan lineage substring
        "Thyroid": "thyroid carcinoma",
        "Stomach/Esophagus": "stomach adenocarcinoma",
        "Pancreas": "pancreatic adenocarcinoma",
        "Bowel": "colon adenocarcinoma",
        "Skin": "skin cutaneous melanoma",
        "Breast": "breast invasive carcinoma",
        "Bladder/Urinary Tract": "bladder urothelial carcinoma",
        "Lung": "lung adenocarcinoma",
        "Soft Tissue": "sarcoma",
        "Testis": "testicular germ cell tumor",
        "Myeloid": "acute myeloid leukemia",
        "Lymphoid": "diffuse large B-cell lymphoma",
        "CNS/Brain": "brain lower grade glioma",
        "Pleura": "mesothelioma",
        "Esophagus/Stomach": "esophageal carcinoma",
        "Biliary Tract": "cholangiocarcinoma",
        "Prostate": "prostate adenocarcinoma",
        "Cervix": "cervical & endocervical cancer",
        "Uterus": "uterine corpus endometrioid carcinoma",
        "Adrenal Gland": "adrenocortical cancer",
        "Head and Neck": "head & neck squamous cell carcinoma",
        "Ovary/Fallopian Tube": "ovarian serous cystadenocarcinoma",
        "Kidney": "kidney clear cell carcinoma",
        "Eye": "uveal melanoma",
        "Liver": "liver hepatocellular carcinoma",
        "Thymus": "thymoma",
        "Peripheral Nervous System": "pheochromocytoma & paraganglioma",
    }
    pancan_med = dict(zip(pancan_stats["lineage"], pancan_stats["median_DM1"]))
    model["pancan_lineage"] = model["OncotreeLineage"].map(LINEAGE_MAP)
    model["DM1_like_score"] = model["pancan_lineage"].map(pancan_med)
    print(f"[D] mapped {model['DM1_like_score'].notna().sum()}/{len(model)} cell lines")
    df = model.dropna(subset=["DM1_like_score"]).merge(crispr, on="ModelID", how="inner")
    print(f"[D] cell lines with both DM1 score and CRISPR: {len(df)}")
    df.to_csv(OUT / "celllines_dm1_crispr.tsv", sep="\t", index=False)
    print(f"[D] saved per-cell-line table ({len(df)} rows)")

    # ---- DM1-high vs DM1-low dependency ----
    df["dm1_bin"] = pd.qcut(df["DM1_like_score"], 3,
                             labels=["low", "mid", "high"], duplicates="drop")
    eff_genes = [g for g in targets if g in df.columns]
    rows = []
    from scipy import stats as st
    for g in eff_genes:
        sub = df[[g, "dm1_bin"]].dropna()
        h = sub.loc[sub.dm1_bin == "high", g].values
        l = sub.loc[sub.dm1_bin == "low", g].values
        if len(h) < 5 or len(l) < 5: continue
        d = (h.mean() - l.mean()) / np.sqrt((h.std(ddof=1)**2 + l.std(ddof=1)**2) / 2 + 1e-12)
        try:
            t, p = st.ttest_ind(h, l, equal_var=False)
        except Exception:
            t, p = np.nan, np.nan
        rows.append({"gene": g, "n_high": len(h), "n_low": len(l),
                     "mean_eff_high": float(h.mean()), "mean_eff_low": float(l.mean()),
                     "delta": float(h.mean() - l.mean()),
                     "cohens_d": float(d), "p": float(p) if not np.isnan(p) else np.nan})
    eff_df = pd.DataFrame(rows)
    if not eff_df.empty:
        eff_df = eff_df.sort_values("cohens_d")
    eff_df.to_csv(OUT / "dm1_high_low_dependency.tsv", sep="\t", index=False)
    print(eff_df.to_string(index=False) if not eff_df.empty else "[D] no dependency rows")

    # ---- Lineage view: per OncotreeLineage, mean essentiality of top hubs ----
    lin_rows = []
    for g in eff_genes:
        for lineage, sub in df.groupby("OncotreeLineage"):
            v = sub[g].dropna()
            if len(v) < 5: continue
            lin_rows.append({"gene": g, "lineage": lineage, "n": len(v),
                             "median_effect": float(v.median()),
                             "pct_essential": float((v < -0.5).mean())})
    lin_df = pd.DataFrame(lin_rows)
    lin_df.to_csv(OUT / "lineage_dependency_summary.tsv", sep="\t", index=False)
    print(f"[D] lineage summary rows: {len(lin_df)}")

    summary = {
        "n_celllines_with_dm1_score": int(len(df)),
        "n_target_genes_with_crispr": len(eff_genes),
        "top_dm1high_dep": eff_df.head(5).to_dict(orient="records"),
        "top_dm1low_dep": eff_df.tail(5).to_dict(orient="records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"[D] summary → {OUT / 'summary.json'}")


if __name__ == "__main__":
    main()
