#!/usr/bin/env python3
"""Paper 11 — Phase E: lineage-portable DM1 reformulation.

The original DM1_like score uses thyroid-anchored signature (RAI_8 +
TF_collapse over PAX8/NKX2-1/FOXE1/HHEX), which causes saturation in
non-thyroid cancers (median ≈ +1.5). To make DM1 axis lineage-agnostic,
we substitute each cancer's canonical differentiation TF panel and
re-score using the Module 4 architecture identified in Paper 12 (P12).

Cancers covered (prognostic top-5 from Phase B Cox + LGG/LUAD/GBM):
  THCA: PAX8 / NKX2-1 / FOXE1 / TG / TPO / TSHR / SLC5A5 / DIO1
  LUAD: NKX2-1 / SFTPC / SFTPB / SFTPA1 / FOXA2 / ABCA3 / NAPSA
  LGG : OLIG1 / OLIG2 / SOX2 / SOX10 / NES / GFAP / S100B
  GBM : OLIG2 / SOX2 / NES / GFAP / S100B / VIM
  HNSC: TP63 / KRT5 / KRT14 / SOX2
  PRAD: AR / NKX3-1 / FOXA1 / KLK3 / TMPRSS2

Inflammation-axis (lineage-agnostic, shared):
  IFNG, STAT1, GZMB, PRF1, CXCL9, CXCL10, CD8A
  CD68, CD163, MRC1, CSF1R, CXCL8 (myeloid)

Metabolic / kinase / epigenetic (Module 4 core):
  MYC, NAMPT, JAK1, JAK2, STAT3, LYN, FYN, SRC, EZH2, DNMT1, HDAC1, HDAC2

Lineage-portable DM1 = - z(LineageTF panel mean)  +  z(Inflammation panel)
                        +  z(MYC/NAMPT/Module 4 core)

Output:
  - phase_E_lineage_specific/lineage_portable_dm1_per_sample.tsv
  - phase_E_lineage_specific/cohort_correlation_with_original.tsv
  - phase_E_lineage_specific/cox_per_cohort_portable.tsv
  - phase_E_lineage_specific/summary.json
"""
from __future__ import annotations
import argparse, json, gzip, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

P11 = Path("/data/thca/repo_results/paper11_pancancer")
PANCAN_EXPR = Path("/data/thca/repo_data/raw/TCGA_pancan/pancan_geneExp.gz")
PHENO = Path("/data/thca/repo_data/raw/TCGA_pancan/phenotype.tsv.gz")
SURV = Path("/data/thca/repo_data/raw/TCGA_pancan/survival.tsv")
SCORES = P11 / "pancan_dm1_scored.tsv"
OUT = P11 / "phase_E_lineage_specific"

LINEAGE_TF = {
    # solid epithelial — primary targets
    "thyroid carcinoma": ["PAX8", "NKX2-1", "FOXE1", "TG", "TPO", "TSHR", "SLC5A5", "DIO1"],
    "lung adenocarcinoma": ["NKX2-1", "SFTPC", "SFTPB", "SFTPA1", "FOXA2", "ABCA3", "NAPSA"],
    "lung squamous cell carcinoma": ["TP63", "KRT5", "KRT14", "SOX2", "PAX9", "KRT6A"],
    "brain lower grade glioma": ["OLIG1", "OLIG2", "SOX2", "SOX10", "NES", "GFAP", "S100B"],
    "glioblastoma multiforme": ["OLIG2", "SOX2", "NES", "GFAP", "S100B", "VIM"],
    "head & neck squamous cell carcinoma": ["TP63", "KRT5", "KRT14", "SOX2"],
    "prostate adenocarcinoma": ["AR", "NKX3-1", "FOXA1", "KLK3", "TMPRSS2"],
    "colon adenocarcinoma": ["CDX2", "MUC2", "VIL1", "KRT20"],
    "rectum adenocarcinoma": ["CDX2", "MUC2", "VIL1", "KRT20"],
    "stomach adenocarcinoma": ["CDX2", "MUC5AC", "MUC6", "PGC", "TFF1"],
    "esophageal carcinoma": ["TP63", "KRT5", "SOX2", "CDX2", "VIL1"],
    "breast invasive carcinoma": ["ESR1", "PGR", "FOXA1", "GATA3", "AR"],
    "pancreatic adenocarcinoma": ["PDX1", "GATA6", "HNF1A", "HNF4A"],
    "liver hepatocellular carcinoma": ["HNF4A", "HNF1A", "ALB", "AFP"],
    "cholangiocarcinoma": ["KRT7", "KRT19", "SOX17", "HNF1B"],
    "skin cutaneous melanoma": ["MITF", "TYR", "MLANA", "DCT", "PMEL"],
    "uveal melanoma": ["MITF", "TYR", "MLANA", "DCT", "PMEL"],
    "bladder urothelial carcinoma": ["GATA3", "FOXA1", "KRT20", "PPARG", "UPK1A", "UPK3A"],
    "cervical & endocervical cancer": ["TP63", "KRT5", "PAX8", "SOX2"],
    "uterine corpus endometrioid carcinoma": ["PGR", "FOXA1", "MUC1", "PAX2", "HAND2"],
    "uterine carcinosarcoma": ["PGR", "FOXA1", "MUC1", "PAX2"],
    "ovarian serous cystadenocarcinoma": ["PAX8", "MUC16", "WT1", "CA125"],
    "adrenocortical cancer": ["NR5A1", "CYP11A1", "CYP21A2", "STAR"],
    "testicular germ cell tumor": ["NANOG", "POU5F1", "SOX2", "LIN28A", "KIT"],
    "kidney clear cell carcinoma": ["PAX8", "PAX2", "CA9", "VHL", "HIF1A"],
    "kidney papillary cell carcinoma": ["PAX8", "PAX2", "MET", "VIM"],
    "kidney chromophobe": ["KIT", "PAX8"],
    "mesothelioma": ["WT1", "MSLN", "CALB2", "UPK3B"],
    "pheochromocytoma & paraganglioma": ["PNMT", "CHGA", "SYP", "DBH", "TH"],
    "thymoma": ["FOXN1", "EOMES", "CD3D"],
    "diffuse large B-cell lymphoma": ["PAX5", "CD19", "BCL6", "IRF4"],
    "acute myeloid leukemia": ["SPI1", "CEBPA", "RUNX1", "GATA2"],
    # heterogeneous — kept for completeness; expect noisier scoring
    "sarcoma": ["MYOD1", "MYOG", "SOX9"],
}

INFLAM_PANEL = ["IFNG", "STAT1", "GZMB", "PRF1", "CXCL9", "CXCL10", "CD8A",
                "CD68", "CD163", "MRC1", "CSF1R", "CXCL8"]
MOD4_CORE = ["MYC", "NAMPT", "JAK1", "JAK2", "STAT3", "LYN", "FYN", "SRC",
             "EZH2", "DNMT1", "HDAC1", "HDAC2"]


def stream_pancan_expr(path, target_genes):
    """Stream pancan_geneExp.gz tab-separated; first row = sample IDs."""
    print(f"[E] streaming {path} for {len(target_genes)} genes …")
    target = set(target_genes)
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        sample_ids = header[1:]
        out = {}
        n = 0
        for line in f:
            n += 1
            if n % 5000 == 0: print(f"  …{n}")
            tab = line.find("\t")
            if tab < 0: continue
            gene = line[:tab].strip()
            if gene in target:
                vals = line[tab + 1:].rstrip("\n").split("\t")
                out[gene] = pd.to_numeric(pd.Series(vals), errors="coerce").values
        print(f"[E] scanned {n} rows; matched {len(out)} genes")
    df = pd.DataFrame(out, index=sample_ids)
    df.index.name = "sample"
    return df


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if not PANCAN_EXPR.exists():
        sys.exit(f"missing {PANCAN_EXPR}")

    all_genes = set()
    for v in LINEAGE_TF.values(): all_genes.update(v)
    all_genes.update(INFLAM_PANEL)
    all_genes.update(MOD4_CORE)
    print(f"[E] target genes: {len(all_genes)}")

    expr = stream_pancan_expr(PANCAN_EXPR, all_genes)
    print(f"[E] expr shape: {expr.shape}")

    # Phenotype to map sample -> lineage
    pheno = pd.read_csv(PHENO, sep="\t", low_memory=False)
    print(f"[E] pheno: {len(pheno)} rows; cols: {list(pheno.columns)[:5]}")
    sid_col = pheno.columns[0]
    lin_col = [c for c in pheno.columns
               if "lineage" in c.lower() or "primary" in c.lower()
               or "cancer" in c.lower() or "disease" in c.lower()][0]
    pheno = pheno[[sid_col, lin_col]].rename(columns={sid_col: "sample", lin_col: "lineage"})
    df = expr.merge(pheno, on="sample", how="inner")
    print(f"[E] merged: {len(df)}")

    # z-score within lineage
    rows = []
    for lineage, sub in df.groupby("lineage"):
        if lineage not in LINEAGE_TF:
            continue
        tf = [g for g in LINEAGE_TF[lineage] if g in sub.columns]
        infl = [g for g in INFLAM_PANEL if g in sub.columns]
        mod4 = [g for g in MOD4_CORE if g in sub.columns]
        if len(tf) < 3 or len(infl) < 3 or len(mod4) < 3:
            print(f"[E] skip {lineage}: tf={len(tf)} infl={len(infl)} mod4={len(mod4)}")
            continue
        sub = sub.copy()
        # z-score within lineage
        for cols in (tf, infl, mod4):
            for g in cols:
                sd = sub[g].std()
                if sd > 0:
                    sub[g + "_z"] = (sub[g] - sub[g].mean()) / sd
                else:
                    sub[g + "_z"] = 0
        sub["lin_TF_score"] = sub[[g + "_z" for g in tf]].mean(axis=1)
        sub["inflam_score"] = sub[[g + "_z" for g in infl]].mean(axis=1)
        sub["mod4_score"]   = sub[[g + "_z" for g in mod4]].mean(axis=1)
        sub["DM1_portable"] = -sub["lin_TF_score"] + sub["inflam_score"] + sub["mod4_score"]
        rows.append(sub[["sample", "lineage", "lin_TF_score", "inflam_score",
                         "mod4_score", "DM1_portable"]])
    portable = pd.concat(rows, ignore_index=True)
    portable.to_csv(OUT / "lineage_portable_dm1_per_sample.tsv",
                    sep="\t", index=False)
    print(f"[E] portable scored: {len(portable)} samples × {portable['lineage'].nunique()} lineages")

    # Correlate with original DM1 (per cohort)
    orig = pd.read_csv(SCORES, sep="\t")
    merged = portable.merge(orig[["sample", "DM1_like"]], on="sample", how="left")
    corr_rows = []
    for lin, sub in merged.groupby("lineage"):
        valid = sub.dropna(subset=["DM1_like", "DM1_portable"])
        if len(valid) < 20: continue
        r, p = stats.pearsonr(valid["DM1_like"], valid["DM1_portable"])
        rs, _ = stats.spearmanr(valid["DM1_like"], valid["DM1_portable"])
        corr_rows.append({"lineage": lin, "n": len(valid),
                          "pearson_r": float(r), "spearman_r": float(rs),
                          "p": float(p)})
    corr_df = pd.DataFrame(corr_rows).sort_values("pearson_r", ascending=False)
    corr_df.to_csv(OUT / "cohort_correlation_with_original.tsv",
                   sep="\t", index=False)
    print(corr_df.to_string(index=False))

    # Cox per-cohort with portable score
    surv = pd.read_csv(SURV, sep="\t", low_memory=False)
    cox_rows = []
    try:
        from lifelines import CoxPHFitter
    except ImportError:
        print("[E] lifelines missing; skipping Cox", file=sys.stderr)
        cox_df = pd.DataFrame()
    else:
        for lin, sub in merged.groupby("lineage"):
            j = sub.merge(surv[["sample", "OS", "OS.time"]], on="sample", how="left")
            j = j.dropna(subset=["DM1_portable", "OS", "OS.time"])
            j = j[j["OS.time"] > 0]
            if len(j) < 30 or j["OS"].sum() < 10: continue
            try:
                cph = CoxPHFitter(penalizer=0.001)
                cph.fit(j[["OS", "OS.time", "DM1_portable"]],
                        duration_col="OS.time", event_col="OS")
                s = cph.summary.loc["DM1_portable"]
                cox_rows.append({"lineage": lin, "n": len(j),
                                 "events": int(j["OS"].sum()),
                                 "HR": float(s["exp(coef)"]),
                                 "p": float(s["p"]),
                                 "concordance": float(cph.concordance_index_)})
            except Exception as e:
                cox_rows.append({"lineage": lin, "error": str(e)[:80]})
        cox_df = pd.DataFrame(cox_rows).sort_values("p")
        cox_df.to_csv(OUT / "cox_per_cohort_portable.tsv", sep="\t", index=False)
        print(f"\n[E] Cox per cohort:\n{cox_df.to_string(index=False)}")

    summary = {
        "n_lineages_scored": int(portable["lineage"].nunique()),
        "n_samples": int(len(portable)),
        "median_correlation_with_original": float(corr_df["pearson_r"].median()) if not corr_df.empty else np.nan,
        "lineages_with_significant_cox_p05": int((cox_df["p"] < 0.05).sum()) if not cox_df.empty else 0,
        "top_cox_hits": cox_df.head(5).to_dict(orient="records") if not cox_df.empty else [],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
