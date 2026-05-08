#!/usr/bin/env python3
"""
Driver class × CCI signature in TCGA-THCA bulk RNA.

Question: Is the CAF/M1-M2/TLS × RAI spatial avoidance signature universal across
driver mutations (BRAF, RAS, TERT, fusion, driver-neg) or specific to a subset?

Strata (from p2d_per_sample_classification.tsv):
  Class1_BRAF_V600E    n≈290 (TCGA majority)
  Class2_RAS           (HRAS / KRAS / NRAS)
  Class3_TERT_promoter
  Class4_RET_fusion
  Class5_other_fusion
  Class6_True_driver_neg

Output:
  spark_tcga_driver_class_CCI.tsv
  driver_class_cci_heatmap.png
"""
from __future__ import annotations
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent.parent.parent
PANCAN = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
PHENO = ROOT / "project/data/raw/TCGA_pancan/phenotype.tsv.gz"
DM = ROOT / "project/results/dark_matter_phase2/p2d_per_sample_classification.tsv"
RES = ROOT / "project/results/03_pathology_poc"

MODULES = {
    "CAF":   ["FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "DCN", "POSTN"],
    "M1_M2": ["CD68", "CD163", "MRC1", "MARCO", "C1QA", "C1QB", "C1QC", "TYROBP"],
    "CD8":   ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG"],
    "TLS":   ["CXCL13", "CCL19", "CCR7", "LTB", "LTA"],
    "RAI":   ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"],
}

DRIVER_LABEL = {
    "Class1_BRAF_V600E": "BRAF V600E",
    "Class1b_BRAF_other": "BRAF other",
    "Class2_RAS_HRAS": "RAS-H",
    "Class2_RAS_KRAS": "RAS-K",
    "Class2_RAS_NRAS": "RAS-N",
    "Class3_TERT_promoter": "TERT prom",
    "Class4_RET_fusion": "RET fusion",
    "Class5_other_fusion": "other fusion",
    "Class6_True_driver_neg": "driver-neg",
}


def load_thca_expr(target_genes):
    pheno = pd.read_csv(PHENO, sep="\t", compression="gzip")
    proj_col = next(c for c in pheno.columns if "primary" in c.lower() or "disease" in c.lower())
    thca = set(pheno["sample"][pheno[proj_col].astype(str).str.contains("thyroid|THCA", case=False, na=False)].tolist())
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


def main():
    all_genes = list({g for v in MODULES.values() for g in v})
    expr = load_thca_expr(all_genes)
    print(f"loaded expr: {expr.shape}")

    mods = {}
    for name, genes in MODULES.items():
        present = [g for g in genes if g in expr.index]
        if len(present) < 2: continue
        sub = expr.loc[present]
        z = sub.subtract(sub.mean(axis=1), axis=0).div(sub.std(axis=1) + 1e-6, axis=0)
        mods[name] = z.mean(axis=0)

    dm = pd.read_csv(DM, sep="\t")
    samp_to_driver = {}
    for _, row in dm.iterrows():
        case = row["tcga_short"]
        cls = row.get("driver_class")
        if pd.notna(cls):
            for sname in expr.columns:
                if sname.startswith(case):
                    samp_to_driver[sname] = cls
    print(f"driver-labeled bulk samples: {len(samp_to_driver)}")
    print("driver class distribution:")
    for cls, n in pd.Series(list(samp_to_driver.values())).value_counts().items():
        print(f"  {cls}: {n}")

    # Group small classes — use buckets to avoid n<20
    buckets = {
        "BRAF": ["Class1_BRAF_V600E", "Class1b_BRAF_other"],
        "RAS":  ["Class2_RAS_HRAS", "Class2_RAS_KRAS", "Class2_RAS_NRAS"],
        "TERT": ["Class3_TERT_promoter"],
        "FUSION": ["Class4_RET_fusion", "Class5_other_fusion"],
        "NEG":  ["Class6_True_driver_neg"],
    }

    rows = []
    for bucket, classes in buckets.items():
        ids = [s for s, c in samp_to_driver.items() if c in classes]
        if len(ids) < 15: continue
        for n1 in ["CAF", "M1_M2", "TLS", "CD8"]:
            common = mods[n1].index.intersection(mods["RAI"].index).intersection(ids)
            if len(common) < 15: continue
            r = float(spearmanr(mods[n1].loc[common], mods["RAI"].loc[common]).statistic)
            rows.append({"bucket": bucket, "module": n1, "n": len(common), "r_vs_RAI": r})

    # Also ALL for reference
    ids_all = list(samp_to_driver.keys())
    for n1 in ["CAF", "M1_M2", "TLS", "CD8"]:
        common = mods[n1].index.intersection(mods["RAI"].index).intersection(ids_all)
        r = float(spearmanr(mods[n1].loc[common], mods["RAI"].loc[common]).statistic)
        rows.append({"bucket": "ALL", "module": n1, "n": len(common), "r_vs_RAI": r})

    df = pd.DataFrame(rows)
    df.to_csv(RES / "spark_tcga_driver_class_CCI.tsv", sep="\t", index=False)
    print(f"\nwrote {(RES / 'spark_tcga_driver_class_CCI.tsv').relative_to(ROOT)}")
    print("\n=== Driver class × CCI heatmap ===")
    pivot = df.pivot(index="module", columns="bucket", values="r_vs_RAI")
    print(pivot.round(3))

    # Heatmap figure
    bucket_order = ["BRAF", "RAS", "TERT", "FUSION", "NEG", "ALL"]
    cols = [b for b in bucket_order if b in pivot.columns]
    pivot = pivot[cols]
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    vmax = float(np.nanmax(np.abs(pivot.values)))
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(pivot.columns))); ax.set_xticklabels(pivot.columns, fontsize=10)
    ax.set_yticks(range(len(pivot.index))); ax.set_yticklabels(pivot.index, fontsize=10)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.iloc[i, j]
            n_val = df[(df.module == pivot.index[i]) & (df.bucket == pivot.columns[j])]["n"]
            n = int(n_val.iloc[0]) if len(n_val) else 0
            txt = f"{v:+.3f}\n(n={n})"
            ax.text(j, i, txt, ha="center", va="center",
                    color="white" if abs(v) > vmax*0.5 else "#222", fontsize=8)
    fig.colorbar(im, ax=ax, label="Spearman r vs RAI lineage")
    ax.set_title("Driver class × CCI signature (TCGA-THCA bulk RNA)\n"
                 "Negative = spatial avoidance signature replicated in this stratum",
                 fontsize=11, fontweight="bold")
    out = RES / "summary_figs/F6_driver_class_cci.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
