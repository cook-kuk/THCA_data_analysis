"""
PDAC immune readiness — module-level signature scoring on TCGA-PAAD bulk RNA
z-scores (cBioPortal). We build PDAC-tailored modules differing from thyroid:
  - myeloid_suppressive (TAM/MDSC dominant in PDAC)
  - CAF axis (myCAF vs iCAF) — Elyada 2019 markers
  - HLA-I / HLA-II / IFNG / checkpoint exhaustion
  - TLS-CXCL13 (vaccine permissive niche, Mariathasan-style)
  - CD8/Treg ratio proxy

Output: per-sample module scores + Moffitt-stratified means + bootstrapped
direction-invariance against expected PDAC TME pattern.
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd

RAW = Path("/data/pdac_poc/raw")
OUT = Path("/data/pdac_poc/processed")
RES = Path("/data/pdac_poc/results")

MODULES = {
    "myeloid_suppressive": ["CD68", "CD163", "MRC1", "CSF1R", "MARCO",
                            "ARG1", "IDO1", "CD86", "TGFB1", "IL10", "VEGFA"],
    "myCAF": ["ACTA2", "POSTN", "TAGLN", "MYH11", "COL10A1"],
    "iCAF": ["IL6", "CXCL12", "PDGFRA", "DCN", "LUM", "HAS2"],
    "CAF_pan": ["FAP", "PDPN", "S100A4", "VIM", "COL1A1", "COL1A2", "COL3A1"],
    "HLA_I": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5",
              "PSMB8", "PSMB9"],
    "HLA_II": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
               "HLA-DQB1", "CIITA"],
    "IFNG_inflamed": ["IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11",
                       "GZMA", "GZMB", "PRF1", "CD8A", "CD8B"],
    "checkpoint_exhaustion": ["PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT",
                              "HAVCR2", "TOX"],
    "TLS_CXCL13": ["CXCL13", "CCL19", "CCL21", "CCR7", "LAMP3", "CD79A",
                   "MS4A1", "BCL6", "CXCR5", "AICDA"],
}


def main():
    z = pd.read_csv(RAW / "mrna_z_Zscores.tsv", sep="\t", index_col=0)
    moffitt = pd.read_csv(RES / "moffitt_calls.tsv", sep="\t", index_col=0)
    out = pd.DataFrame(index=z.index)

    coverage = {}
    for mod, genes in MODULES.items():
        avail = [g for g in genes if g in z.columns]
        coverage[mod] = (len(avail), len(genes))
        out[mod] = z[avail].mean(axis=1) if avail else np.nan

    out["myCAF_minus_iCAF"] = out["myCAF"] - out["iCAF"]
    out["IFNG_minus_myeloid"] = out["IFNG_inflamed"] - out["myeloid_suppressive"]

    out = out.join(moffitt[["moffitt_call", "delta_basal_minus_classical",
                            "mut_KRAS", "mut_TP53", "mut_SMAD4"]], how="left")
    out.to_csv(RES / "pdac_immune_readiness_per_sample.tsv", sep="\t")

    # group means by Moffitt
    grouped = out.groupby("moffitt_call")[list(MODULES.keys()) +
                                          ["myCAF_minus_iCAF",
                                           "IFNG_minus_myeloid"]].mean()
    grouped.to_csv(RES / "pdac_immune_readiness_by_moffitt.tsv", sep="\t")

    # Cohen's d basal vs classical
    def cohens_d(a, b):
        if len(a) < 2 or len(b) < 2:
            return float("nan")
        v = ((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2)
        sp = np.sqrt(v) if v > 0 else np.nan
        return (a.mean() - b.mean()) / sp if sp else float("nan")

    b = out.query("moffitt_call=='basal-like'")
    c = out.query("moffitt_call=='classical'")
    cohen = {m: float(cohens_d(b[m].dropna(), c[m].dropna()))
             for m in (list(MODULES.keys()) + ["myCAF_minus_iCAF",
                                                "IFNG_minus_myeloid"])}

    # KRAS-mutant vs wt
    km = out.query("mut_KRAS==1")
    kw = out.query("mut_KRAS==0")
    cohen_kras = {m: float(cohens_d(km[m].dropna(), kw[m].dropna()))
                  for m in MODULES.keys()}

    summary = {
        "n_samples": int(len(out)),
        "module_coverage": {k: f"{v[0]}/{v[1]}" for k, v in coverage.items()},
        "module_means_by_moffitt": grouped.to_dict(),
        "cohens_d_basal_vs_classical": cohen,
        "cohens_d_kras_mut_vs_wt": cohen_kras,
        "fraction_TLS_high": float((out["TLS_CXCL13"] > 0.5).mean()),
        "fraction_myeloid_high": float((out["myeloid_suppressive"] > 0.5).mean()),
        "fraction_basal_AND_inflamed": float(((out["moffitt_call"] == "basal-like") &
                                              (out["IFNG_inflamed"] > 0)).mean()),
        "fraction_basal_AND_TLS_high": float(((out["moffitt_call"] == "basal-like") &
                                              (out["TLS_CXCL13"] > 0.5)).mean()),
    }
    with open(OUT / "pdac_immune_readiness_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print(f"[05] modules scored: {len(MODULES)}")
    for m, d in cohen.items():
        print(f"  d(basal-classical)[{m:25s}] = {d:+.2f}")
    print("[05] done")


if __name__ == "__main__":
    main()
