#!/usr/bin/env python3
"""Paper 11 — Phase G: pan-cancer MSigDB Hallmark × DM1 axis.

For each of 50 MSigDB Hallmark gene sets, score per sample
(z-mean of constituent genes within lineage), then correlate with
DM1_portable per lineage. Identify Hallmarks that consistently
co-vary with DM1 across cancers — these are the pan-cancer DM1
biology pillars.

Outputs:
  - phase_G_hallmark/hallmark_dm1_corr_long.tsv  (lineage × hallmark r,p)
  - phase_G_hallmark/hallmark_pancan_summary.tsv (per hallmark: median r,
                                                    n_sig_lineages_pos,
                                                    n_sig_lineages_neg)
  - phase_G_hallmark/heatmap.png
  - phase_G_hallmark/summary.json
"""
from __future__ import annotations
import argparse, gzip, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

P11 = Path("/data/thca/repo_results/paper11_pancancer")
PANCAN_EXPR = Path("/data/thca/repo_data/raw/TCGA_pancan/pancan_geneExp.gz")
HALLMARK_GMT = Path("/home/seungho/personal/THCA_data_analysis/project/metadata/genesets/MSigDB_Hallmark_2020.gmt")
OUT = P11 / "phase_G_hallmark"


def parse_gmt(path):
    """MSigDB Hallmark gmt: name<TAB>url<TAB>g1<TAB>g2<TAB>... per line."""
    sets = {}
    with open(path) as f:
        for L in f:
            parts = L.rstrip("\n").split("\t")
            if len(parts) < 3: continue
            name = parts[0]
            genes = [g for g in parts[2:] if g]
            sets[name] = genes
    return sets


def stream_pancan_expr(path, target_genes):
    print(f"[G] streaming {path} for {len(target_genes)} genes …")
    target = set(target_genes)
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        sample_ids = header[1:]
        out = {}
        n = 0
        for line in f:
            n += 1
            if n % 5000 == 0: print(f"  …{n}", file=sys.stderr)
            tab = line.find("\t")
            if tab < 0: continue
            gene = line[:tab].strip()
            if gene in target:
                vals = line[tab + 1:].rstrip("\n").split("\t")
                out[gene] = pd.to_numeric(pd.Series(vals), errors="coerce").values
        print(f"[G] scanned {n} rows; matched {len(out)} genes")
    df = pd.DataFrame(out, index=sample_ids)
    df.index.name = "sample"
    return df


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sets = parse_gmt(HALLMARK_GMT)
    print(f"[G] {len(sets)} Hallmark sets loaded")
    all_genes = sorted({g for genes in sets.values() for g in genes})
    expr = stream_pancan_expr(PANCAN_EXPR, all_genes)
    print(f"[G] expr: {expr.shape}")

    # join lineage + DM1
    portable = pd.read_csv(P11 / "phase_E_lineage_specific" /
                            "lineage_portable_dm1_per_sample.tsv", sep="\t")
    df = expr.merge(portable[["sample", "lineage", "DM1_portable"]],
                    on="sample", how="inner")
    print(f"[G] joined: {len(df)} samples × {df['lineage'].nunique()} lineages")

    # Per-lineage Hallmark module scoring + correlation
    rows = []
    for lin, sub in df.groupby("lineage"):
        if len(sub) < 30: continue
        sub = sub.copy()
        # z-score within lineage
        for col in [c for c in sub.columns if c not in ("sample", "lineage", "DM1_portable")]:
            std = sub[col].std()
            sub[col + "_z"] = (sub[col] - sub[col].mean()) / (std if std > 0 else 1)
        for hm, genes in sets.items():
            present = [g + "_z" for g in genes if g + "_z" in sub.columns]
            if len(present) < 5: continue
            mod_score = sub[present].mean(axis=1)
            valid = pd.DataFrame({"m": mod_score, "d": sub["DM1_portable"]}).dropna()
            if len(valid) < 20: continue
            r, p = stats.spearmanr(valid["m"], valid["d"])
            rows.append({"lineage": lin, "hallmark": hm,
                         "n": len(valid), "n_genes": len(present),
                         "spearman_r": float(r), "p": float(p)})

    res = pd.DataFrame(rows)
    if not res.empty:
        # FDR per lineage (50 hallmarks each)
        out_rows = []
        for lin, g in res.groupby("lineage"):
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
        res = pd.concat(out_rows, ignore_index=True).sort_values(
            ["lineage", "fdr"])
    res.to_csv(OUT / "hallmark_dm1_corr_long.tsv", sep="\t", index=False)
    print(f"[G] wrote {OUT / 'hallmark_dm1_corr_long.tsv'}: {len(res)} rows")

    # Per-hallmark summary
    pancan_sum = (res
                  .groupby("hallmark")
                  .agg(median_r=("spearman_r", "median"),
                       n_lineages=("lineage", "nunique"),
                       n_pos_fdr10=("spearman_r",
                                     lambda s: ((s > 0) & (res.loc[s.index, "fdr"] < 0.1)).sum()),
                       n_neg_fdr10=("spearman_r",
                                     lambda s: ((s < 0) & (res.loc[s.index, "fdr"] < 0.1)).sum()))
                  .sort_values("median_r", ascending=False)
                  .reset_index())
    pancan_sum.to_csv(OUT / "hallmark_pancan_summary.tsv",
                      sep="\t", index=False)
    print("[G] top hallmarks by median r:")
    print(pancan_sum.head(15).to_string(index=False))
    print("\n[G] bottom hallmarks (anti-DM1):")
    print(pancan_sum.tail(10).to_string(index=False))

    # Heatmap
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        wide = res.pivot_table(index="hallmark", columns="lineage",
                                values="spearman_r", aggfunc="first")
        # Order rows by median r
        wide = wide.loc[pancan_sum["hallmark"].values]
        fig, ax = plt.subplots(figsize=(max(8, 0.4 * wide.shape[1]),
                                          max(10, 0.22 * wide.shape[0])))
        im = ax.imshow(wide.values, aspect="auto", cmap="RdBu_r",
                       vmin=-0.5, vmax=0.5)
        ax.set_xticks(np.arange(wide.shape[1]))
        ax.set_xticklabels(wide.columns, rotation=90, fontsize=7)
        ax.set_yticks(np.arange(wide.shape[0]))
        ax.set_yticklabels([h.replace("HALLMARK_", "")[:30]
                             for h in wide.index], fontsize=7)
        plt.colorbar(im, ax=ax, label="Spearman r (DM1_portable)")
        ax.set_title("Hallmark × DM1_portable (per lineage)")
        plt.tight_layout()
        plt.savefig(OUT / "heatmap.png", dpi=140)
        plt.close()
        print(f"[G] heatmap → {OUT / 'heatmap.png'}")
    except Exception as e:
        print(f"[G] heatmap failed: {e}", file=sys.stderr)

    summary = {
        "n_lineages": int(res["lineage"].nunique()),
        "n_hallmarks": int(res["hallmark"].nunique()),
        "top10_pos_hallmarks": pancan_sum.head(10).to_dict(orient="records"),
        "bottom10_hallmarks": pancan_sum.tail(10).to_dict(orient="records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"[G] summary → {OUT / 'summary.json'}")


if __name__ == "__main__":
    main()
