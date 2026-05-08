#!/usr/bin/env python3
"""
TCGA-THCA Cox survival of CCI signature modules (CAF, M1/M2, TLS, CD8, RAI_lineage).

For 5 module scores from bulk RNA, compute univariate Cox HR for PFI / OS.
Stratify by DM1/DM2 cluster (when label available).

Inputs:
  pancan_geneExp.gz + phenotype.tsv.gz + survival.tsv (TCGA pancan reference)
  project/results/dark_matter_phase2/p2d_per_sample_classification.tsv (DM labels)

Output:
  spark_tcga_cox_modules.tsv
"""
from __future__ import annotations
import gzip, sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent.parent
PANCAN = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
PHENO = ROOT / "project/data/raw/TCGA_pancan/phenotype.tsv.gz"
SURV = ROOT / "project/data/raw/TCGA_pancan/survival.tsv"
DM = ROOT / "project/results/dark_matter_phase2/p2d_per_sample_classification.tsv"
OUT = ROOT / "project/results/03_pathology_poc/spark_tcga_cox_modules.tsv"

MODULES = {
    "CAF":   ["FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "DCN", "POSTN"],
    "M1_M2": ["CD68", "CD163", "MRC1", "MARCO", "C1QA", "C1QB", "C1QC", "TYROBP"],
    "CD8":   ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG"],
    "TLS":   ["CXCL13", "CCL19", "CCR7", "LTB", "LTA"],
    "RAI":   ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"],
}


def load_thca_expr(target_genes):
    pheno = pd.read_csv(PHENO, sep="\t", compression="gzip")
    proj_col = next(c for c in pheno.columns if "primary" in c.lower() or "disease" in c.lower())
    samp_col = "sample"
    thca = set(pheno[samp_col][pheno[proj_col].astype(str).str.contains("thyroid|THCA", case=False, na=False)].tolist())
    rows = []
    target_set = set(target_genes)
    with gzip.open(PANCAN, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        keep_idx = [i for i, c in enumerate(header) if c in thca]
        keep_names = [header[i] for i in keep_idx]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in target_set:
                vals = [float(parts[i]) if parts[i] not in {"", "NA"} else np.nan for i in keep_idx]
                rows.append([parts[0]] + vals)
            if len(rows) >= len(target_set):
                break
    df = pd.DataFrame(rows, columns=["gene"] + keep_names).set_index("gene")
    return df


def cox_univariate(time, event, x):
    """Univariate Cox via lifelines if available; else logrank median split."""
    keep = ~(np.isnan(time) | np.isnan(event) | np.isnan(x))
    n = int(keep.sum())
    if n < 30:
        return {"n": n, "HR": float("nan"), "p": float("nan"), "method": "skip"}
    try:
        from lifelines import CoxPHFitter
        df = pd.DataFrame({"t": time[keep], "e": event[keep], "x": x[keep]})
        cph = CoxPHFitter().fit(df, "t", "e")
        s = cph.summary.iloc[0]
        return {"n": n, "HR": float(s["exp(coef)"]), "p": float(s["p"]), "method": "cox"}
    except Exception:
        # logrank median split fallback
        try:
            from lifelines.statistics import logrank_test
            med = np.median(x[keep])
            hi = x[keep] > med
            r = logrank_test(time[keep][hi], time[keep][~hi], event[keep][hi], event[keep][~hi])
            return {"n": n, "HR": float("nan"), "p": float(r.p_value), "method": "logrank-median"}
        except Exception as e:
            return {"n": n, "HR": float("nan"), "p": float("nan"), "method": str(e)[:30]}


def main():
    all_genes = list({g for v in MODULES.values() for g in v})
    expr = load_thca_expr(all_genes)
    print(f"loaded expr: {expr.shape}")

    # module score per sample (z within gene, then mean)
    mods = {}
    for name, genes in MODULES.items():
        present = [g for g in genes if g in expr.index]
        if len(present) < 2: continue
        sub = expr.loc[present]
        z = sub.subtract(sub.mean(axis=1), axis=0).div(sub.std(axis=1) + 1e-6, axis=0)
        mods[name] = z.mean(axis=0)

    surv = pd.read_csv(SURV, sep="\t")
    samp_col = "sample"
    print(f"survival cols: {surv.columns.tolist()[:8]}")
    surv_thca = surv[surv[samp_col].isin(expr.columns)].set_index(samp_col)
    print(f"THCA survival samples: {len(surv_thca)}")

    # PFI / OS
    rows = []
    for endpoint in ["OS", "PFI", "DSS", "DFI"]:
        time_col = f"{endpoint}.time"
        evt_col = endpoint
        if time_col not in surv_thca.columns or evt_col not in surv_thca.columns:
            continue
        time_v = surv_thca[time_col].astype(float)
        evt_v = surv_thca[evt_col].astype(float)
        for mname, mvec in mods.items():
            common = mvec.index.intersection(surv_thca.index)
            if len(common) < 30: continue
            time_a = time_v.loc[common].values
            evt_a = evt_v.loc[common].values
            x_a = mvec.loc[common].values
            res = cox_univariate(time_a, evt_a, x_a)
            res.update({"endpoint": endpoint, "module": mname, "n_used": len(common)})
            rows.append(res)
            print(f"  {endpoint} × {mname}: n={res['n']}, HR={res['HR']:.3g}, p={res['p']:.3g}, method={res['method']}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep="\t", index=False)
    print(f"\nwrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
