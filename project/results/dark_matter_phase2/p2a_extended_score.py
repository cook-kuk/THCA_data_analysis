"""Quick: extended FVPTC score with all HVG-available differentiation markers."""
from pathlib import Path
import json
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/dark_matter_phase2")
a = ad.read_h5ad("/home/seungho/personal/THCA_data_analysis/project/results/v17_lu2023/GSE193581_hvg_adata.h5ad")

# Extended FVPTC: include PAX8, KRT19 negative direction
FVPTC_EXT = [g for g in ["TG", "TPO", "TSHR", "DIO2", "PAX8"] if g in a.var_names]
sc.tl.score_genes(a, gene_list=FVPTC_EXT, score_name="score_fvptc_ext")
print(f"Extended FVPTC genes: {FVPTC_EXT}")

# Also try all 8-gene HVG members vs FVPTC_EXT
ptc_mask = (a.obs["histology"] == "PTC") & (a.obs["author_celltype"] == "Malignant cell")
ptc = a[ptc_mask].copy()

# All histology malignant cells (PTC + ATC) to see if pattern holds in ATC too
all_mal_mask = a.obs["author_celltype"] == "Malignant cell"
all_mal = a[all_mal_mask].copy()

results = {}
for label, sub in [("PTC malignant", ptc), ("All malignant (PTC+ATC)", all_mal)]:
    r_pool, p_pool = stats.pearsonr(sub.obs["DM_score"], sub.obs["score_fvptc_ext"])
    print(f"\n=== {label} (n={len(sub)}) ===")
    print(f"  Pooled r (DM_score, FVPTC_ext) = {r_pool:.3f}")
    pp = []
    for s, g in sub.obs.groupby("sample", observed=True):
        if len(g) >= 30:
            r, p = stats.pearsonr(g["DM_score"], g["score_fvptc_ext"])
            pp.append({"sample": s, "n": len(g), "r": r})
    pp = pd.DataFrame(pp).sort_values("r", ascending=False)
    print(pp.to_string(index=False))
    results[label] = {
        "n_total": int(len(sub)),
        "pooled_r": float(r_pool),
        "median_patient_r": float(pp["r"].median()),
        "patients_r_gt_0.7": int((pp["r"] > 0.7).sum()),
        "patients_r_gt_0.5": int((pp["r"] > 0.5).sum()),
        "n_patients": int(len(pp)),
    }

(OUT / "p2a_extended_summary.json").write_text(json.dumps(results, indent=2))
print(f"\nSaved.")
print(json.dumps(results, indent=2))
