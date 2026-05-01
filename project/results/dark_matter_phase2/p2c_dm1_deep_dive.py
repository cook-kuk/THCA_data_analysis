"""P2-C DM1 deep dive — mechanism mining in true Dark Matter (driver-neg cPTC)."""
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

# Load Phase 1 master + clinical extended
df = pd.read_csv(ROOT / "results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
master_v17 = pd.read_csv(ROOT / "results/v17/tables/sample_master_v17_full.tsv", sep="\t", low_memory=False)
master_v17["tcga_short"] = master_v17["sample_id"].str[:12]
master_v17_tcga = master_v17[(master_v17["dataset"] == "TCGA-THCA") & (master_v17["normal_vs_tumor"] == "tumor")]
df = df.merge(master_v17_tcga[["tcga_short", "fusion_classes", "histology_subtype", "ajcc_stage_group",
                               "mutation_genes", "tcga12"]].drop_duplicates("tcga_short"),
              on="tcga_short", how="left", suffixes=("", "_m"))

# DM1 cohort = DM × cluster=DM1
dm1 = df[df["dm_status"] & (df["v17_dark_cluster"] == "DM1")].copy()
dm2 = df[df["dm_status"] & (df["v17_dark_cluster"] == "DM2")].copy()
print(f"DM1: {len(dm1)} | DM2: {len(dm2)}")

# === 1. mutation_genes content scan ===
# Parse mutation_genes column (format like "['BRAF']" or "['DICER1', 'EIF1AX']")
def parse_genes(s):
    if pd.isna(s) or s == "" or s == "[]":
        return []
    # Extract gene names from string-formatted list
    return re.findall(r"'([A-Z0-9\-]+)'", str(s))

dm1["genes"] = dm1["mutation_genes"].apply(parse_genes)
dm2["genes"] = dm2["mutation_genes"].apply(parse_genes)

# Count gene occurrences
def gene_freq(d):
    flat = [g for genes in d["genes"] for g in genes]
    return pd.Series(flat).value_counts() if flat else pd.Series([], dtype=int)

dm1_genes = gene_freq(dm1)
dm2_genes = gene_freq(dm2)
print(f"\nDM1 mutation_genes (top 20):\n{dm1_genes.head(20)}")
print(f"\nDM2 mutation_genes (top 20):\n{dm2_genes.head(20)}")

# === 2. Genes enriched in DM1 vs DM2 ===
all_genes = sorted(set(list(dm1_genes.index) + list(dm2_genes.index)))
enrich_rows = []
for g in all_genes:
    dm1_pos = int(dm1_genes.get(g, 0))
    dm2_pos = int(dm2_genes.get(g, 0))
    if dm1_pos + dm2_pos < 2:  # skip ultra-rare
        continue
    table = [[dm1_pos, len(dm1) - dm1_pos], [dm2_pos, len(dm2) - dm2_pos]]
    odds, p = stats.fisher_exact(table)
    enrich_rows.append({
        "gene": g, "DM1_n": dm1_pos, "DM1_pct": dm1_pos/len(dm1)*100,
        "DM2_n": dm2_pos, "DM2_pct": dm2_pos/len(dm2)*100,
        "OR_DM1_vs_DM2": float(odds), "p": float(p),
    })
enrich = pd.DataFrame(enrich_rows).sort_values("p")
enrich.to_csv(OUT / "p2c_dm1_vs_dm2_gene_enrichment.tsv", sep="\t", index=False)
print(f"\n=== Top 20 enriched genes (DM1 vs DM2) ===\n{enrich.head(20).to_string(index=False)}")

# === 3. Fusion landscape in DM1 ===
print(f"\n=== Fusion classes in DM1 ===")
print(dm1["fusion_classes"].value_counts(dropna=False))
print(f"\n=== Fusion classes in DM2 ===")
print(dm2["fusion_classes"].value_counts(dropna=False))

# === 4. Histology distribution within DM1 (just to confirm cPTC majority) ===
print(f"\nDM1 histology:\n{dm1['histology_subtype'].value_counts(dropna=False)}")
print(f"\nDM2 histology:\n{dm2['histology_subtype'].value_counts(dropna=False)}")

# === 5. Stage / age / clinical ===
print(f"\nDM1 stage:\n{dm1['ajcc_stage_group'].value_counts(dropna=False).head()}")
print(f"DM2 stage:\n{dm2['ajcc_stage_group'].value_counts(dropna=False).head()}")
print(f"\nAge: DM1 mean={pd.to_numeric(dm1['age'], errors='coerce').mean():.1f}, DM2 mean={pd.to_numeric(dm2['age'], errors='coerce').mean():.1f}")
print(f"PFI events: DM1={int(dm1['PFI'].sum())}/{int(dm1['PFI'].notna().sum())}, DM2={int(dm2['PFI'].sum())}/{int(dm2['PFI'].notna().sum())}")

# === 6. Save summary ===
summary = {
    "n_dm1": int(len(dm1)),
    "n_dm2": int(len(dm2)),
    "dm1_top_mutated_genes": dm1_genes.head(15).to_dict(),
    "dm2_top_mutated_genes": dm2_genes.head(15).to_dict(),
    "dm1_histology": dm1["histology_subtype"].value_counts(dropna=False).to_dict(),
    "dm2_histology": dm2["histology_subtype"].value_counts(dropna=False).to_dict(),
    "dm1_age_mean": float(pd.to_numeric(dm1["age"], errors="coerce").mean()),
    "dm2_age_mean": float(pd.to_numeric(dm2["age"], errors="coerce").mean()),
    "dm1_PFI_event_rate": float(dm1["PFI"].sum() / dm1["PFI"].notna().sum()) if dm1["PFI"].notna().sum() else None,
    "dm2_PFI_event_rate": float(dm2["PFI"].sum() / dm2["PFI"].notna().sum()) if dm2["PFI"].notna().sum() else None,
    "top_5_significant_enrichments": enrich.head(5).to_dict(orient="records"),
}
(OUT / "p2c_dm1_deep_dive_summary.json").write_text(json.dumps(summary, indent=2, default=str))
print(f"\nSaved p2c_dm1_deep_dive_summary.json")
