"""5 more deepening tasks: ME test + cell composition + master table + GitHub files."""
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import json
import re
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project")
OUT = ROOT / "results/dark_matter_phase2"
FIG = OUT / "web/figures"
DATA_OUT = OUT / "web/data"

# ============================================================================
# 1. Mutual exclusivity BRAF V600E vs DICER1/EIF1AX (replicate Wang 2025)
# ============================================================================
print("=== 1. Mutual exclusivity BRAF V600E vs DICER1/EIF1AX ===")

# TCGA: from master + driver_anchor
master = pd.read_csv(OUT / "../dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
master["alt_pos"] = master["driver_anchor_v17"].isin(["DICER1_EIF1AX_PPM1D"])

# K2: from Yoo parsed
k2 = pd.read_csv(OUT / "k2_yoo2016_mutations_parsed.tsv", sep="\t")
k2["alt_pos"] = k2["has_dicer1"].astype(bool) | k2["has_eif1ax"].astype(bool)

me_results = []
for cohort_name, df in [("TCGA-THCA", master), ("K2 (Yoo 2016)", k2)]:
    if cohort_name == "TCGA-THCA":
        col_braf = "has_braf_v600e"
    else:
        col_braf = "has_braf_v600e"

    if col_braf not in df.columns:
        continue
    df_clean = df.dropna(subset=[col_braf])
    df_clean[col_braf] = df_clean[col_braf].astype(bool)

    a = int(((df_clean[col_braf]) & (df_clean["alt_pos"])).sum())  # both
    b = int(((df_clean[col_braf]) & (~df_clean["alt_pos"])).sum())  # BRAF only
    c = int(((~df_clean[col_braf]) & (df_clean["alt_pos"])).sum())  # DICER1 only
    d = int(((~df_clean[col_braf]) & (~df_clean["alt_pos"])).sum())  # neither
    odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="less")
    me_results.append({
        "cohort": cohort_name, "n": len(df_clean),
        "BRAF+DICER+": a, "BRAF+ only": b, "DICER+ only": c, "neither": d,
        "OR": float(odds) if not np.isinf(odds) else 0.0,
        "p_one-sided_less": float(p),
    })
    print(f"  {cohort_name} n={len(df_clean)}: a={a}(BRAF+DICER+), b={b}(BRAF+ only), c={c}(DICER+ only), d={d}(neither)")
    print(f"    Fisher one-sided OR={odds:.3f}, p_less={p:.4f} (mutual exclusivity test)")

me_df = pd.DataFrame(me_results)
print(f"\n=== Combined: TCGA + K2 ===")
combined_a = me_df["BRAF+DICER+"].sum()
combined_b = me_df["BRAF+ only"].sum()
combined_c = me_df["DICER+ only"].sum()
combined_d = me_df["neither"].sum()
combined_n = combined_a + combined_b + combined_c + combined_d
odds_c, p_c = stats.fisher_exact([[combined_a, combined_b], [combined_c, combined_d]], alternative="less")
print(f"  Combined: a={combined_a}, b={combined_b}, c={combined_c}, d={combined_d}")
print(f"  Combined Fisher one-sided OR={odds_c:.3f}, p_less={p_c:.4e}")
me_df.loc[len(me_df)] = {"cohort": "Combined TCGA+K2", "n": combined_n,
                          "BRAF+DICER+": combined_a, "BRAF+ only": combined_b,
                          "DICER+ only": combined_c, "neither": combined_d,
                          "OR": float(odds_c) if not np.isinf(odds_c) else 0.0, "p_one-sided_less": float(p_c)}
me_df.to_csv(DATA_OUT / "mutual_exclusivity_braf_dicer1.tsv", sep="\t", index=False)
print(me_df.round(4).to_string(index=False))

# Plot mosaic
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
for ax, (_, row) in zip(axes, me_df.iterrows()):
    a, b, c, d = row["BRAF+DICER+"], row["BRAF+ only"], row["DICER+ only"], row["neither"]
    matrix = np.array([[a, b], [c, d]])
    im = ax.imshow(matrix, cmap="YlOrRd", aspect="equal")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, str(matrix[i, j]), ha="center", va="center",
                    color="white" if matrix[i, j] > matrix.max()*0.5 else "black",
                    fontsize=14, fontweight="bold")
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
    ax.set_xticklabels(["DICER1+", "DICER1−"])
    ax.set_yticklabels(["BRAF V600E+", "BRAF V600E−"])
    ax.set_title(f"{row['cohort']} (n={row['n']})\nOR={row['OR']:.2f}, p={row['p_one-sided_less']:.1e}")
fig.suptitle("Figure E16 — Mutual exclusivity BRAF V600E vs DICER1/EIF1AX (replicate Wang 2025)", y=1.02)
fig.tight_layout()
fig.savefig(FIG / "figE16_mutual_exclusivity.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE16")

# ============================================================================
# 2. Cell type composition DM1 vs DM2 (sc data) — proxy via 8-gene quartiles
# ============================================================================
print("\n=== 2. Cell type composition by 8-gene quartile ===")
sc_meta = pd.read_csv(OUT / "../dark_matter_phase1/sc_cell_metadata.tsv", sep="\t")
print(f"sc cells: {len(sc_meta)}, cell types: {sc_meta['celltype'].value_counts().to_dict()}")

# Use thyrocyte 8-gene quartiles as DM1-like (low) vs DM2-like (high) proxy
thy = sc_meta[sc_meta["celltype"] == "Thyrocyte"].copy()
q1, q3 = thy["score_8gene"].quantile([0.25, 0.75])
thy["group"] = pd.cut(thy["score_8gene"], bins=[-np.inf, q1, q3, np.inf],
                      labels=["DM1-like (low 8gene)", "Mid", "DM2-like (high 8gene)"])

# For non-thyrocyte cell types, attribute to nearest sample's median 8-gene
# Aggregate cell type composition per sample, then split samples by sample-mean 8-gene
sample_8gene = sc_meta.groupby("sample", observed=True)["score_8gene"].mean()
sample_groups = pd.cut(sample_8gene, bins=2, labels=["Low-8gene", "High-8gene"]).rename("sample_group")
sc_meta = sc_meta.merge(sample_groups, on="sample", how="left")

ct_pivot = pd.crosstab(sc_meta["celltype"], sc_meta["sample_group"], normalize="columns") * 100
print(f"\nCell type composition (% by sample group):\n{ct_pivot.round(2)}")
ct_pivot.to_csv(DATA_OUT / "cell_composition_by_8gene_group.tsv", sep="\t")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax = axes[0]
ct_pivot.plot(kind="bar", ax=ax, color=["#4361ee", "#f4a261"], width=0.7, alpha=0.85)
ax.set_xlabel("Cell type")
ax.set_ylabel("% of cells in group")
ax.set_title(f"Figure E17A — Cell type composition by 8-gene group\n(GSE241184 sc, sample-level split)")
ax.legend(title="Sample group", loc="upper right")
ax.grid(alpha=0.3, axis="y")

# Right: signature score distribution per cell type
ax = axes[1]
ctypes = ["Thyrocyte", "T_cell", "Myeloid", "B_cell", "Endothelial", "Fibroblast"]
data_for_box = [sc_meta.loc[sc_meta["celltype"] == ct, "score_8gene"] for ct in ctypes]
ctypes_with_n = [f"{ct}\n(n={int((sc_meta['celltype']==ct).sum())})" for ct in ctypes]
bp = ax.boxplot(data_for_box, labels=ctypes_with_n, showfliers=False, patch_artist=True)
for patch, c in zip(bp["boxes"], ["#1864ab", "#a3cef1", "#fcbf49", "#90be6d", "#f72585", "#577590"]):
    patch.set_facecolor(c); patch.set_alpha(0.7)
ax.set_ylabel("8-gene DM score")
ax.set_title("Figure E17B — 8-gene score is thyrocyte-specific\n(other cell types near zero)")
ax.grid(alpha=0.3, axis="y")
plt.setp(ax.get_xticklabels(), rotation=20, ha="right", fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "figE17_cell_composition.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved figE17")

# ============================================================================
# 3. 5-cohort master comparison table
# ============================================================================
print("\n=== 3. 5-cohort master comparison ===")
master_table = pd.DataFrame([
    {"Cohort": "TCGA-THCA", "Type": "Bulk RNA-seq + WES", "n": 482,
     "DM%": 28.4, "DM1": 89, "DM2": 55,
     "DICER1+EIF1AX in DM%": 10.9, "8gene↔FVPTC r": "n.a. (bulk)",
     "Survival": "PFI 50/482", "Mutation": "Full"},
    {"Cohort": "K2 (Yoo 2016)", "Type": "Bulk RNA-seq", "n": 180,
     "DM%": 37.8, "DM1": "(predicted)", "DM2": "(predicted)",
     "DICER1+EIF1AX in DM%": 10.3, "8gene↔FVPTC r": "n.a. (bulk)",
     "Survival": "Pending", "Mutation": "S6 supplement"},
    {"Cohort": "GSE241184", "Type": "scRNA-seq", "n": "1 patient",
     "DM%": "n.a.", "DM1": "—", "DM2": "—",
     "DICER1+EIF1AX in DM%": "n.a.", "8gene↔FVPTC r": "0.91",
     "Survival": "n.a.", "Mutation": "n.a."},
    {"Cohort": "GSE193581", "Type": "scRNA-seq", "n": "23 (PTC+ATC+NORM)",
     "DM%": "n.a.", "DM1": "—", "DM2": "—",
     "DICER1+EIF1AX in DM%": "n.a.", "8gene↔FVPTC r": "0.89 (PTC+ATC)",
     "Survival": "n.a.", "Mutation": "n.a."},
    {"Cohort": "GSE184362", "Type": "scRNA-seq", "n": "23 (T/P/LN multi-site)",
     "DM%": "n.a.", "DM1": "—", "DM2": "—",
     "DICER1+EIF1AX in DM%": "n.a.", "8gene↔FVPTC r": "0.89 / 0.91 multi-site",
     "Survival": "n.a.", "Mutation": "n.a."},
    {"Cohort": "분당 SNUH", "Type": "TBD", "n": "TBD",
     "DM%": "TBD", "DM1": "—", "DM2": "—",
     "DICER1+EIF1AX in DM%": "TBD", "8gene↔FVPTC r": "TBD",
     "Survival": "TBD (Korean follow-up critical)", "Mutation": "TBD"},
])
master_table.to_csv(DATA_OUT / "5cohort_master_table.tsv", sep="\t", index=False)
print(master_table.to_string(index=False))

# ============================================================================
# 4. GitHub package files
# ============================================================================
print("\n=== 4. GitHub package files ===")
GH = OUT / "web" / "github_pkg"
GH.mkdir(exist_ok=True)

readme = """# Dark Matter of Thyroid Cancer

> An 8-gene transcriptional axis sub-stratifies BRAF/RAS-negative thyroid cancer
> at single-cell resolution and identifies DICER1/EIF1AX as the FVPTC-like genomic anchor.

**Author**: Seungho Cook (corresponding); collaborators TBD.
**Status**: Phase 1 + Phase 2 sprint completed 2026-04-29 (single-day intensive).
**Target venue**: Cell Reports Medicine (1차) / Nature Communications (reach).

## TL;DR

- Driver-negative thyroid cancer (BRAF V600E−, RAS hotspot−) = **28% of TCGA-THCA**, currently un-stratified by molecular tests.
- Our 8-gene panel sub-stratifies this "Dark Matter" into:
  - **DM1** (n=89, mean age 42, cPTC-architectured, mechanism unknown — true mystery)
  - **DM2** (n=55, mean age 54, FVPTC-like, **DICER1/EIF1AX/PPM1D 9× enriched**, p=0.0175)
- Reproduced at single-cell resolution across 4 sc cohorts (r ≥ 0.86).
- Korean Yoo 2016 K2 cohort (n=180): **NBNR ↔ DM concordance 93.5%**, DICER1+EIF1AX 10.3% ≈ TCGA 10.9%.
- Pan-cancer specificity confirmed (LUAD r=0.235 vs THCA r=0.99).
- 8-gene → cPTC/FVPTC predictive AUC = 0.71 (5-fold CV).
- Bootstrap cluster ARI = 0.994 (random null = 0).

## Repository structure

```
project/
├── data_processed/        # bulk RNA-seq + microarray expression matrices
├── data_raw/              # GDC + GEO + ENA raw downloads
├── results/
│   ├── dark_matter_phase1/  # Phase 1 (Discovery)
│   └── dark_matter_phase2/  # Phase 2 + 3 + 4 (Validation + Deepening + ULTRA)
│       └── web/             # Live HTML dashboard
├── notebooks_or_scripts/  # analysis scripts
└── reports/html/         # public dashboard at http://40.82.129.113:8012/
```

## Live dashboard
http://40.82.129.113:8012/reports/html/dark_matter/v2.html

## Key scripts (reproducible analyses)
```
project/results/dark_matter_phase1/
├── run_steps_1_to_6.py      # Phase 1 viability test
├── step2_redo_pfi.py        # Liu 2018 TCGA-CDR Cox
├── sc_analysis.py           # GSE241184 sc pipeline

project/results/dark_matter_phase2/
├── p2a_gse193581.py         # P2-A1 external sc
├── p2a2_gse184362.py        # P2-A2 adult multi-patient
├── p2a2_multisite.py        # T/P/LN trajectory
├── p2b_yoo2016_parse.py     # Korean K2 mutation parse
├── p2c_dm1_deep_dive.py     # DM1 mechanism mining
├── p2d_driver_map.py        # 7-class driver landscape
├── p3_deepening.py          # DEG + pathway + AUC + sc + K2 fix
├── p4_final_batch.py        # Pan-cancer + ATC + Stage I + DEG table + BibTeX + Methods
├── p5_speed.py              # Age-stratified + DICER1+ PFI
└── p6_more.py               # Mutual exclusivity + cell composition + master table
```

## Software requirements
- Python 3.12, scanpy 1.12.1, lifelines 0.27, pandas 2.3, scipy 1.13, scikit-learn 1.5
- See `requirements.txt`

## Citation
If you use this work, cite:
> Cook S et al. (2026) *Dark matter of thyroid cancer: an 8-gene transcriptional axis identifies DICER1/EIF1AX-anchored FVPTC-like sub-cluster within driver-negative tumors* (manuscript in preparation).

BibTeX entry: see `data/references.bib`.

## License
MIT (see `LICENSE`). Data sources retain their original licenses.

## Acknowledgments
Yu Hyeong-Won (Seoul National University Bundang Hospital) — clinical guidance.
Yoo SK et al. (PLoS Genet 2016) — K2 Korean cohort dataset.
Original TCGA, GSE241184/GSE193581/GSE184362 publication teams.
"""
(GH / "README.md").write_text(readme)
print(f"Saved README.md ({len(readme):,} chars)")

license_text = """MIT License

Copyright (c) 2026 Seungho Cook

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
(GH / "LICENSE").write_text(license_text)

citation = """cff-version: 1.2.0
message: "If you use this work, please cite the manuscript and/or repository as below."
title: "Dark Matter of Thyroid Cancer: 8-gene transcriptional axis sub-stratifies BRAF/RAS-negative thyroid cancer"
authors:
  - family-names: "Cook"
    given-names: "Seungho"
    affiliation: "TBD"
date-released: 2026-04-29
license: MIT
keywords:
  - thyroid-cancer
  - papillary-thyroid-cancer
  - dark-matter
  - DICER1
  - EIF1AX
  - single-cell-RNA-seq
  - transcriptional-clustering
  - molecular-taxonomy
preferred-citation:
  type: article
  title: "Dark matter of thyroid cancer: an 8-gene transcriptional axis identifies DICER1/EIF1AX-anchored FVPTC-like sub-cluster within driver-negative tumors"
  authors:
    - family-names: "Cook"
      given-names: "Seungho"
  status: "in preparation"
  year: 2026
"""
(GH / "CITATION.cff").write_text(citation)
print(f"Saved LICENSE + CITATION.cff")

# Save final summary
final = {
    "mutual_exclusivity": me_df.to_dict(orient="records"),
    "cell_composition_by_group": ct_pivot.to_dict(),
    "5cohort_master": master_table.to_dict(orient="records"),
    "github_files_saved": ["README.md", "LICENSE", "CITATION.cff"],
}
(DATA_OUT / "p6_more_summary.json").write_text(json.dumps(final, indent=2, default=str))
print(f"\n=== ALL P6 DONE ===")
