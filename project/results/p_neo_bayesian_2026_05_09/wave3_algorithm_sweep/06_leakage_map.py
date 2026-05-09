#!/usr/bin/env python3
"""Leakage map: source × source peptide×HLA Jaccard overlap.

Uses bundle.tsv as the master peptide pool (8 sources, train + ext).
Per-source set = unique (peptide, HLA_norm) tuples.
Jaccard(A,B) = |A∩B| / |A∪B|.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent
df = pd.read_csv(ROOT.parent / "bundle.tsv", sep="\t")
print(f"bundle: {df.shape}, sources: {df['source'].nunique()}")
print(df["source"].value_counts())

# Build per-source set of (peptide, HLA_norm)
src_sets = {}
for src, grp in df.groupby("source"):
    src_sets[src] = set(zip(grp["peptide"], grp["HLA_norm"]))

# Also include peptide-only sets (for peptide-level leakage which is what most tools use)
src_pep = {src: set(grp["peptide"]) for src, grp in df.groupby("source")}

sources = sorted(src_sets.keys())
n = len(sources)
J_pmh = np.zeros((n, n))
J_pep = np.zeros((n, n))
counts = np.zeros((n, n), dtype=int)
for i, a in enumerate(sources):
    for j, b in enumerate(sources):
        A_pmh, B_pmh = src_sets[a], src_sets[b]
        A_p, B_p = src_pep[a], src_pep[b]
        if len(A_pmh | B_pmh) > 0:
            J_pmh[i, j] = len(A_pmh & B_pmh) / len(A_pmh | B_pmh)
        if len(A_p | B_p) > 0:
            J_pep[i, j] = len(A_p & B_p) / len(A_p | B_p)
        counts[i, j] = len(A_pmh & B_pmh)

# Save TSVs
pd.DataFrame(J_pmh, index=sources, columns=sources).to_csv(ROOT / "leakage_jaccard_peptideHLA.tsv", sep="\t")
pd.DataFrame(J_pep, index=sources, columns=sources).to_csv(ROOT / "leakage_jaccard_peptide.tsv", sep="\t")
pd.DataFrame(counts, index=sources, columns=sources).to_csv(ROOT / "leakage_overlap_counts.tsv", sep="\t")

# Top-3 overlapping pairs (excluding diagonal)
upper = []
for i in range(n):
    for j in range(i + 1, n):
        upper.append((sources[i], sources[j], J_pmh[i, j], J_pep[i, j], counts[i, j],
                       len(src_sets[sources[i]]), len(src_sets[sources[j]])))
upper_df = pd.DataFrame(upper, columns=["source_a", "source_b", "jaccard_pepHLA", "jaccard_pep_only", "n_overlap", "n_a", "n_b"])
upper_df = upper_df.sort_values("jaccard_pepHLA", ascending=False)
upper_df.to_csv(ROOT / "leakage_pairs_ranked.tsv", sep="\t", index=False)
print("\nTop 10 overlapping source pairs (peptide+HLA Jaccard):")
print(upper_df.head(10).to_string(index=False))

# Heatmap figure
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, M, title in [(axes[0], J_pmh, "Peptide+HLA Jaccard"),
                     (axes[1], J_pep, "Peptide-only Jaccard")]:
    im = ax.imshow(M, cmap="magma_r", vmin=0, vmax=max(0.05, M.max()))
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(sources, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(sources, fontsize=8)
    for i in range(n):
        for j in range(n):
            v = M[i, j]
            if v > 0.001:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", color="white" if v > 0.3 else "black", fontsize=7)
    ax.set_title(title)
    plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)

fig.suptitle(f"Wave 3 — Source × Source Leakage Map (n={n} sources)")
plt.tight_layout()
fig.savefig(ROOT / "fig_leakage_map.png", dpi=200)
fig.savefig(ROOT / "fig_leakage_map.pdf")
print(f"\nsaved fig_leakage_map.png/.pdf with {n} sources")
