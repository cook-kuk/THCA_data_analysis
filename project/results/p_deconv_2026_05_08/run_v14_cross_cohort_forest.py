"""v14: cross-cohort MAPK × {Panel-8, TDS-16, TDS_8only} forest.

Cohorts (5):
  TCGA-THCA           n=572  bulk RNA-seq z-score (already in v10/v13)
  Lee / GSE213647     n=632  bulk RNA-seq z-score (already in v10/v13; Ensembl→symbol)
  GSE126698 Korean    n=28   bulk RNA-seq z-score (already in v10)
  GSE286332 Korean+HT n=18   raw TPM → log2 → within-cohort z (NEW; PTC n=9 vs PTC+HT n=9)
  GSE76039 advanced   n=37   microarray z-score (PDTC n=17 + ATC n=20)

Outputs:
  v14_cross_cohort_forest.tsv      — ρ + 95% CI per cohort × score
  v14_cross_cohort_forest_long.tsv — long form for plotting
"""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

PANEL_8   = ['DIO1','FOXE1','NKX2-1','PAX8','SLC5A5','TG','TPO','TSHR']
TDS_16    = ['DIO1','DIO2','DUOX1','DUOX2','FOXE1','GLIS3','NKX2-1','PAX8',
             'SLC26A4','SLC5A5','SLC5A8','TG','THRA','THRB','TPO','TSHR']
TDS_8ONLY = sorted(set(TDS_16) - set(PANEL_8))
MAPK = ['DUSP4','DUSP5','DUSP6','SPRY2','SPRY4','ETV4','ETV5','PHLDA1','CCND1']

def fisher_ci(rho, n, alpha=0.05):
    """Fisher z-transform 95% CI for Spearman ρ."""
    if n < 4: return (np.nan, np.nan)
    z = np.arctanh(rho)
    se = 1.0 / np.sqrt(n - 3)
    zcrit = stats.norm.ppf(1 - alpha/2)
    lo = np.tanh(z - zcrit*se)
    hi = np.tanh(z + zcrit*se)
    return (lo, hi)

def score_corr_block(z, cohort_name, n_total):
    """z is genes × samples (rows=symbols). Compute MAPK × {Panel,TDS16,TDS8only} ρ."""
    avail_mapk = [g for g in MAPK if g in z.index]
    rows = []
    if len(avail_mapk) < 5:
        return rows
    mapk_score = z.loc[avail_mapk].mean(axis=0)
    for label, glist in [('Panel_8', PANEL_8), ('TDS_16', TDS_16), ('TDS_8only', TDS_8ONLY)]:
        avail = [g for g in glist if g in z.index]
        if len(avail) < max(3, len(glist)//2): continue
        s = z.loc[avail].mean(axis=0)
        rho, p = stats.spearmanr(mapk_score, s)
        lo, hi = fisher_ci(rho, len(s))
        rows.append({
            'cohort':cohort_name, 'n':len(s), 'n_total':n_total,
            'panel':label, 'genes_used':len(avail), 'genes_total':len(glist),
            'mapk_genes_used':len(avail_mapk),
            'spearman_rho':rho, 'spearman_p':p, 'ci_lo':lo, 'ci_hi':hi,
        })
    return rows

all_rows = []

# 1. TCGA-THCA z-score
print("=== TCGA-THCA ===")
tcga = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv",
                   sep="\t", index_col=0)
print(f"  shape: {tcga.shape}")
all_rows += score_corr_block(tcga, 'TCGA-THCA', tcga.shape[1])

# 2. Lee/GSE213647 z-score (Ensembl→symbol)
print("\n=== Lee/GSE213647 ===")
lee = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_zscore.tsv",
                  sep="\t", index_col=0)
gmap = pd.read_csv("/data/thca/repo_results/v17p3/tables/F1_gene_recovery_mapping.tsv", sep="\t")
ens_to_sym = dict(zip(gmap['ensembl'], gmap['symbol']))
target = sorted(set(PANEL_8 + TDS_16 + MAPK))
ens_target = [e for e, s in ens_to_sym.items() if s in target]
lee_sub = lee.loc[lee.index.intersection(ens_target)].copy()
lee_sub.index = [ens_to_sym[e] for e in lee_sub.index]
lee_sub = lee_sub.groupby(lee_sub.index).mean()
print(f"  recovered shape: {lee_sub.shape}")
all_rows += score_corr_block(lee_sub, 'Lee_GSE213647', lee.shape[1])

# 3. GSE126698 Korean PTC + HT z-score
print("\n=== GSE126698 ===")
gse126 = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/GSE126698_rnaseq_expression_zscore.tsv",
                     sep="\t", index_col=0)
print(f"  shape: {gse126.shape}")
all_rows += score_corr_block(gse126, 'GSE126698', gse126.shape[1])

# 4. GSE286332 Korean PTC vs PTC+HT (raw TPM → log2 → within-cohort z)
print("\n=== GSE286332 (NEW for v14) ===")
g286 = pd.read_csv("/data/thca/v17_korean/GSE286332/GSE286332_all_sample_rawdata.txt.gz",
                   sep="\t")
# TPM columns end with _TPM
tpm_cols = [c for c in g286.columns if c.endswith('_TPM')]
sym_col = 'Gene_Symbol'
g286 = g286[[sym_col] + tpm_cols].copy()
g286 = g286.dropna(subset=[sym_col])
g286 = g286.groupby(sym_col).sum()  # collapse duplicate symbols
g286.columns = [c.replace('_TPM','') for c in g286.columns]
mat = np.log2(g286.values + 1.0)
z = (mat - mat.mean(axis=1, keepdims=True)) / (mat.std(axis=1, ddof=1, keepdims=True) + 1e-9)
gse286_z = pd.DataFrame(z, index=g286.index, columns=g286.columns)
print(f"  z-score shape: {gse286_z.shape}; sample IDs: {list(gse286_z.columns)[:3]}...")
all_rows += score_corr_block(gse286_z, 'GSE286332', gse286_z.shape[1])

# 5. GSE76039 microarray z-score (advanced disease)
print("\n=== GSE76039 (PDTC + ATC) ===")
g76 = pd.read_csv("/data/thca/data_processed/microarray/GSE76039_microarray_expression_zscore.tsv",
                  sep="\t", index_col=0)
print(f"  shape: {g76.shape}")
all_rows += score_corr_block(g76, 'GSE76039_advanced', g76.shape[1])

# Compile + write
df = pd.DataFrame(all_rows)
print("\n=== Forest table ===")
print(df[['cohort','panel','n','spearman_rho','ci_lo','ci_hi','spearman_p']].round(4).to_string(index=False))

df.to_csv(OUT/"v14_cross_cohort_forest.tsv", sep="\t", index=False)
print(f"\nWrote {OUT}/v14_cross_cohort_forest.tsv  ({len(df)} rows)")

# Random-effects-ish summary on the Panel-8 row (DerSimonian-Laird; placeholder simple Fisher z mean)
panel_rows = df[df['panel']=='Panel_8'].copy()
panel_rows['z'] = np.arctanh(panel_rows['spearman_rho'])
panel_rows['se'] = 1.0 / np.sqrt(panel_rows['n'] - 3)
w = 1.0 / panel_rows['se']**2
mean_z = (w * panel_rows['z']).sum() / w.sum()
mean_rho = np.tanh(mean_z)
mean_se = 1.0 / np.sqrt(w.sum())
ci = np.tanh([mean_z - 1.96*mean_se, mean_z + 1.96*mean_se])
print(f"\nPooled (fixed-effect Fisher-z) Panel-8 ρ = {mean_rho:+.3f}  95% CI [{ci[0]:+.3f}, {ci[1]:+.3f}]  (5 cohorts, total n = {panel_rows['n'].sum()})")

print("\n=== DONE v14 ===")
