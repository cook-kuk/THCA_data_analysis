"""v5 D NNLS-only writeout: compute DM d + 2-way concordance from existing NNLS fractions.

LinearSVR top-10K converged too slowly with 10K samples × 8 features in liblinear;
the 4/4 high-confidence axis direction match (NNLS Pu full vs Lu HVG nu-SVR v2)
is sufficient for the robustness claim, so finalize with NNLS only.
"""
from __future__ import annotations
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

celltypes_lu = ['B cell','Endothelial cell','Epithelial cell','Fibroblast',
                'Malignant cell','Myeloid cell','NK cell','T cell']

def cohen_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    n1, n2 = len(a), len(b)
    if n1<2 or n2<2 or (a.var()+b.var())==0: return 0
    pooled = np.sqrt(((n1-1)*a.var(ddof=1)+(n2-1)*b.var(ddof=1))/(n1+n2-2))
    return (a.mean()-b.mean())/pooled if pooled>0 else 0

# Load NNLS fractions
fr = pd.read_csv(OUT / "fractions_TCGA_full_Pu_NNLS.tsv", sep="\t", index_col=0)
fr['sample_id'] = fr.index
canonical = pd.read_csv("/data/thca/repo_results/v17p3/tables/A2_dm_score_full_cohort.tsv", sep="\t")
canonical = canonical[canonical['dataset']=='TCGA-THCA'].copy()
m = fr.merge(canonical[['sample_id','dm_like']], on='sample_id')
m = m[m['dm_like'].isin(['DM1_like','DM2_like'])]

print(f"NNLS Pu full DM merged: n_DM1={(m['dm_like']=='DM1_like').sum()}, n_DM2={(m['dm_like']=='DM2_like').sum()}")
ct_cols = [c for c in fr.columns if c not in ('sample_id','dm_like')]
rows = []
for ct in ct_cols:
    a = m.loc[m['dm_like']=='DM1_like', ct].values
    b = m.loc[m['dm_like']=='DM2_like', ct].values
    rows.append({'cell_type':ct,'cohen_d':cohen_d(a,b),'DM1_mean':a.mean(),'DM2_mean':b.mean()})
D_nnls = pd.DataFrame(rows).sort_values('cohen_d', key=abs, ascending=False)
D_nnls.to_csv(OUT / "v5D_dm1_dm2_full_pu_NNLS.tsv", sep="\t", index=False)
print("DM1 vs DM2 NNLS Pu full:")
print(D_nnls.round(4).to_string(index=False))

# 2-way concordance with v2 Lu HVG nu-SVR
v2 = pd.read_csv(OUT / "celltype_d_method_grid_v2.tsv", sep="\t")
v2_d = dict(zip(v2['cell_type'], v2['nu-SVR']))
nnls_d = dict(zip(D_nnls['cell_type'], D_nnls['cohen_d']))

rows = []
for ct in celltypes_lu:
    rows.append({
        'cell_type': ct,
        'd_LuHVG_nuSVR_v2': v2_d.get(ct, 0),
        'd_PuFull_NNLS': nnls_d.get(ct, 0),
    })
concord = pd.DataFrame(rows)
concord['informative'] = (
    (concord['d_LuHVG_nuSVR_v2'].abs() > 0.05)
    & (concord['d_PuFull_NNLS'].abs() > 0.05)
)
concord['sign_match'] = (
    np.sign(concord['d_LuHVG_nuSVR_v2'])
    == np.sign(concord['d_PuFull_NNLS'])
) | ~concord['informative']  # uninformative cells get sign-match True (no info to disagree)
concord.to_csv(OUT / "v5D_two_way_concordance.tsv", sep="\t", index=False)

print("\n=== 2-way concordance: Lu HVG nu-SVR (v2 primary) vs Pu full NNLS ===")
print(concord.round(3).to_string(index=False))
informative_match = (concord['informative'] & concord['sign_match']).sum()
informative_total = concord['informative'].sum()
print(f"\nInformative-cell sign-match: {informative_match}/{informative_total}")
print(f"Done — wrote v5D_dm1_dm2_full_pu_NNLS.tsv, v5D_two_way_concordance.tsv")
