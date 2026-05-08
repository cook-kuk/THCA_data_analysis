#!/usr/bin/env python3
"""
TCGA-THCA bulk RNA replicate of the CAF × RAI_lineage spatial avoidance signature.

Logic: in bulk RNA we cannot see spatial avoidance directly, but if CAF and RAI
modules are anatomically segregated within each tumor, their bulk-tumor
correlations should be NEGATIVE across patients (each tumor's mix of
CAF-rich vs RAI-rich tissue varies). If they're co-expressed within a cell, we
expect a positive bulk correlation. Negative bulk correlation = consistent with
spatial segregation.

Replicate also: CAF × RAI in DM1 vs DM2 dark-matter clusters.

Inputs:
  project/data/raw_data/.../tcga_thca_tpm.tsv (or similar) OR existing scored
  project/results/dark_matter_phase2/p2d_per_sample_classification.tsv

Output:
  spark_tcga_bulk_caf_rai_replicate.tsv
"""
from __future__ import annotations
import sys, glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent.parent.parent

CAF = ["FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "DCN", "POSTN"]
RAI = ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"]
CD8 = ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG", "NKG7"]
TLS = ["CXCL13", "CCL19", "CCR7", "LTB", "LTA"]
M12 = ["CD68", "CD163", "MRC1", "MARCO", "C1QA", "C1QB"]


PANCAN = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
PHENO = ROOT / "project/data/raw/TCGA_pancan/phenotype.tsv.gz"


def load_thca_expression(genes_needed):
    """Stream pancan_geneExp.gz, keep only THCA samples and target genes."""
    import gzip
    # Get THCA sample IDs from phenotype
    pheno = pd.read_csv(PHENO, sep="\t", compression="gzip")
    samp_col = next((c for c in ["sample", "sampleID", "submitter_id.samples"] if c in pheno.columns), None)
    proj_col = next((c for c in ["_PRIMARY_DISEASE", "primary_disease", "cancer type abbreviation",
                                  "_primary_disease"] if c in pheno.columns), None)
    print(f"phenotype cols: {pheno.columns[:8].tolist()} ; using sample={samp_col}, proj={proj_col}")
    if proj_col is None:
        # try grep "thyroid" anywhere
        text_cols = [c for c in pheno.columns if pheno[c].dtype == object]
        for c in text_cols:
            if pheno[c].astype(str).str.contains("thyroid", case=False, na=False).any():
                proj_col = c; print(f"  → using {c} as project col"); break
    thca_mask = pheno[proj_col].astype(str).str.contains("thyroid|THCA", case=False, na=False)
    thca_samples = set(pheno[samp_col][thca_mask].tolist())
    print(f"THCA samples in phenotype: {len(thca_samples)}")

    rows = []
    target_set = set(genes_needed)
    with gzip.open(PANCAN, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cols_to_keep = [i for i, c in enumerate(header) if c in thca_samples]
        kept_sample_names = [header[i] for i in cols_to_keep]
        print(f"matrix columns total {len(header)} ; THCA columns kept {len(cols_to_keep)}")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            g = parts[0]
            if g in target_set:
                vals = [float(parts[i]) if parts[i] not in {"", "NA"} else np.nan for i in cols_to_keep]
                rows.append([g] + vals)
            if len(rows) >= len(target_set):
                break
    df = pd.DataFrame(rows, columns=["gene"] + kept_sample_names).set_index("gene")
    return df


def main():
    all_genes = list(set(CAF + RAI + CD8 + TLS + M12))
    df = load_thca_expression(all_genes)
    print(f"matrix: {df.shape}")

    # Now rows = genes, columns = samples. Build module mean (z within gene).
    def module(genes):
        present = [g for g in genes if g in df.index]
        if len(present) < 2: return None
        sub = df.loc[present]
        z = sub.subtract(sub.mean(axis=1), axis=0).div(sub.std(axis=1) + 1e-6, axis=0)
        return z.mean(axis=0)

    caf = module(CAF); rai = module(RAI); cd8 = module(CD8); tls = module(TLS); m12 = module(M12)
    if caf is None or rai is None:
        print("CAF or RAI genes not found; skip"); sys.exit(0)

    rows = []
    for n1, m1 in [("CAF", caf), ("CD8", cd8), ("TLS", tls), ("M1_M2", m12)]:
        if m1 is None: continue
        for n2, m2 in [("RAI_lineage", rai)]:
            common = m1.index.intersection(m2.index)
            r = float(spearmanr(m1.loc[common], m2.loc[common]).statistic)
            rows.append({"a": n1, "b": n2, "spearman_r": r, "n": len(common)})

    out = pd.DataFrame(rows)
    out_path = ROOT / "project/results/03_pathology_poc/spark_tcga_bulk_caf_rai_replicate.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"\nwrote {out_path.relative_to(ROOT)}")
    print(out.to_string(index=False))

    # DM1 vs DM2 stratified
    dm = pd.read_csv(ROOT / "project/results/dark_matter_phase2/p2d_per_sample_classification.tsv", sep="\t")
    dm = dm.set_index("tcga_short")
    samp_dm = dm["v17_dark_cluster"].dropna()
    common_s = samp_dm.index.intersection(caf.index)
    print(f"\ncross-checked DM-labeled bulk samples: {len(common_s)}")
    if len(common_s) >= 30:
        for c1 in [("CAF", caf), ("CD8", cd8), ("M1_M2", m12), ("TLS", tls)]:
            if c1[1] is None: continue
            n1, m1 = c1
            for cluster in ["DM1", "DM2"]:
                ids = samp_dm[samp_dm == cluster].index.intersection(m1.index)
                ids = ids.intersection(rai.index)
                if len(ids) < 10: continue
                r = float(spearmanr(m1.loc[ids], rai.loc[ids]).statistic)
                print(f"  {n1} × RAI_lineage in {cluster} (n={len(ids)}): Spearman = {r:+.3f}")


if __name__ == "__main__":
    main()
