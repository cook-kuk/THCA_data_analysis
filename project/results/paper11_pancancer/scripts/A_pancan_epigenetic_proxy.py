#!/usr/bin/env python3
"""Paper 11 — Phase A (revised): pan-cancer epigenetic dysregulation as
DM1 mechanism layer.

Direct pan-cancer HM450 methylation files are all 403/access-blocked
from Xena/GDC/PanCanAtlas hubs (programmatic fetch). As a proxy, we
score the **expression** of methylation/chromatin regulators (DNMTs,
HDACs, EZH2, KDM1A) and ask whether this epigenetic-dysregulation
index correlates with DM1 axis (original + portable) per cancer.

Mechanism interpretation: high DNMT/EZH2/HDAC expression → predicted
hyperactive de-novo methylation + chromatin compaction → silenced
lineage-TF promoters → DM1-high state. This was the THCA-specific
finding (TPO Δβ=+0.42 d=+2.30 p=1.9e-18 in TCGA HM450); we generalize
the prediction here without per-CpG access.

Outputs:
  - phase_A_epigenetic/epi_index_per_sample.tsv
  - phase_A_epigenetic/epi_index_dm1_corr_per_lineage.tsv
  - phase_A_epigenetic/summary.json
"""
from __future__ import annotations
import argparse, json, gzip, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

P11 = Path("/data/thca/repo_results/paper11_pancancer")
PANCAN_EXPR = Path("/data/thca/repo_data/raw/TCGA_pancan/pancan_geneExp.gz")
OUT = P11 / "phase_A_epigenetic"

EPI_PANEL = {
    "writers": ["DNMT1", "DNMT3A", "DNMT3B"],
    "PRC2": ["EZH2", "EZH1", "SUZ12", "EED"],
    "HDACs": ["HDAC1", "HDAC2", "HDAC3", "HDAC6"],
    "demethylase": ["KDM1A", "KDM6A", "TET1", "TET2"],
    "writers_marks": ["KMT2A", "KMT2D", "DOT1L"],
}
ALL_GENES = sorted({g for v in EPI_PANEL.values() for g in v})


def stream_pancan_expr(path, target_genes):
    print(f"[A] streaming {path} for {len(target_genes)} genes …")
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
        print(f"[A] scanned {n} rows; matched {len(out)} genes")
    df = pd.DataFrame(out, index=sample_ids)
    df.index.name = "sample"
    return df


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    expr = stream_pancan_expr(PANCAN_EXPR, ALL_GENES)
    print(f"[A] expr: {expr.shape}")

    # DM1 scores
    orig = pd.read_csv(P11 / "pancan_dm1_scored.tsv", sep="\t")
    portable = pd.read_csv(P11 / "phase_E_lineage_specific" /
                            "lineage_portable_dm1_per_sample.tsv", sep="\t")
    df = expr.merge(orig[["sample", "lineage", "DM1_like"]], on="sample", how="left")
    df = df.merge(portable[["sample", "DM1_portable"]], on="sample", how="left")
    print(f"[A] joined: {len(df)}")

    # z-score within lineage, sum into modules
    rows_long = []
    for lin, sub in df.groupby("lineage"):
        if len(sub) < 30: continue
        sub = sub.copy()
        for cls, genes in EPI_PANEL.items():
            present = [g for g in genes if g in sub.columns]
            if len(present) < 2: continue
            zsub = sub[present].apply(lambda c: (c - c.mean()) / (c.std() if c.std() > 0 else 1))
            sub[f"epi_{cls}"] = zsub.mean(axis=1)
        # composite epi index = mean of writers + PRC2 + HDACs (silencing machinery)
        cols = [c for c in ("epi_writers", "epi_PRC2", "epi_HDACs") if c in sub.columns]
        if cols:
            sub["epi_silencing_index"] = sub[cols].mean(axis=1)
        rows_long.append(sub)
    big = pd.concat(rows_long, ignore_index=True)
    keep = ["sample", "lineage", "DM1_like", "DM1_portable",
            "epi_writers", "epi_PRC2", "epi_HDACs", "epi_demethylase",
            "epi_writers_marks", "epi_silencing_index"]
    big = big[[c for c in keep if c in big.columns]]
    big.to_csv(OUT / "epi_index_per_sample.tsv", sep="\t", index=False)

    # Per-lineage Spearman corr: epi_silencing_index vs DM1_like, DM1_portable
    rows = []
    for lin, sub in big.groupby("lineage"):
        for axis in ("DM1_like", "DM1_portable"):
            valid = sub[["epi_silencing_index", axis]].dropna()
            if len(valid) < 20: continue
            rs, p = stats.spearmanr(valid["epi_silencing_index"], valid[axis])
            rows.append({"lineage": lin, "axis": axis, "n": len(valid),
                         "spearman_r": float(rs), "p": float(p)})
    res = pd.DataFrame(rows)
    if not res.empty:
        # FDR per axis
        out_rows = []
        for axis, g in res.groupby("axis"):
            g = g.copy().sort_values("p")
            valid = g["p"].notna()
            m = valid.sum()
            if m == 0:
                out_rows.append(g); continue
            ranks = np.arange(1, m + 1)
            pvals = g.loc[valid, "p"].values
            bh = pvals * m / ranks
            bh = np.minimum.accumulate(bh[::-1])[::-1]
            g.loc[valid, "fdr"] = bh
            out_rows.append(g)
        res = pd.concat(out_rows, ignore_index=True).sort_values(["axis", "fdr"])
    res.to_csv(OUT / "epi_index_dm1_corr_per_lineage.tsv", sep="\t", index=False)
    print(res.to_string(index=False))

    summary = {
        "n_lineages_tested": int(res["lineage"].nunique()) if not res.empty else 0,
        "n_DM1_like_pos_sig_fdr10": int(((res.axis == "DM1_like") &
                                         (res.spearman_r > 0) & (res.fdr < 0.1)).sum()),
        "n_DM1_portable_pos_sig_fdr10": int(((res.axis == "DM1_portable") &
                                              (res.spearman_r > 0) & (res.fdr < 0.1)).sum()),
        "median_r_DM1_like": float(res.loc[res.axis == "DM1_like", "spearman_r"].median()) if not res.empty else np.nan,
        "median_r_DM1_portable": float(res.loc[res.axis == "DM1_portable", "spearman_r"].median()) if not res.empty else np.nan,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))
    print(f"[A] note: full HM450 pancan fetch blocked from Xena/GDC; "
          f"epigenetic-regulator-expression index used as RNA-level proxy.")


if __name__ == "__main__":
    main()
