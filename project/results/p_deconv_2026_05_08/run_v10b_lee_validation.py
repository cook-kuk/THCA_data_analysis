"""v10b: Lee external validation of MAPK→panel correlation (Ensembl→symbol mapping)."""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

MAPK = ['DUSP4','DUSP5','DUSP6','SPRY2','SPRY4','ETV4','ETV5','PHLDA1','CCND1']
PANEL = ['DIO1','FOXE1','NKX2-1','PAX8','SLC5A5','TG','TPO','TSHR']

# Load Lee + map Ensembl -> symbol
lee = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_zscore.tsv", sep="\t", index_col=0)
gmap = pd.read_csv("/data/thca/repo_results/v17p3/tables/F1_gene_recovery_mapping.tsv", sep="\t")
ens_to_sym = dict(zip(gmap['ensembl'], gmap['symbol']))

target = list(set(MAPK + PANEL))
ens_target = [e for e,s in ens_to_sym.items() if s in target]
print(f"Lee Ensembl IDs in target panel: {len(ens_target)} (target={len(target)})")

lee_sub = lee.loc[lee.index.isin(ens_target)].copy()
lee_sub.index = [ens_to_sym[e] for e in lee_sub.index]
lee_sub = lee_sub.groupby(lee_sub.index).mean()  # collapse duplicates if any
print(f"  recovered symbols: {sorted(lee_sub.index)}")

available_mapk = [g for g in MAPK if g in lee_sub.index]
available_panel = [g for g in PANEL if g in lee_sub.index]
print(f"  MAPK: {len(available_mapk)}/{len(MAPK)}, Panel: {len(available_panel)}/{len(PANEL)}")

mapk_score = lee_sub.loc[available_mapk].mean(axis=0)
panel_score = lee_sub.loc[available_panel].mean(axis=0)
rho, p = stats.spearmanr(mapk_score, panel_score)
pr, pp = stats.pearsonr(mapk_score, panel_score)
print(f"\nLee n={len(mapk_score)}: MAPK × Panel Spearman ρ = {rho:.4f}, p = {p:.4g}")
print(f"  Pearson r = {pr:.4f}, p = {pp:.4g}")

# Stratify by tumor vs normal (Lee has Normal/PTC/PDTC/ATC labels per original metadata)
# Try to load metadata
try:
    lee_meta = pd.read_csv("/data/thca/repo_results/v17_korean/GSE213647_panel_score.tsv", sep="\t")
    print(f"\nLee metadata: {lee_meta.shape}")
    print(f"  cols: {list(lee_meta.columns)[:8]}")
    print(f"  hist: {lee_meta['histology'].value_counts().to_dict() if 'histology' in lee_meta.columns else 'no hist'}")
    sample_col = [c for c in lee_meta.columns if 'sample' in c.lower() or 'gsm' in c.lower()][0]
    print(f"  sample col: {sample_col}")
    samp_df = pd.DataFrame({'sample':mapk_score.index,'mapk':mapk_score.values,'panel':panel_score.values})
    samp_df = samp_df.merge(lee_meta, left_on='sample', right_on=sample_col, how='inner')
    print(f"  merged: {samp_df.shape}")
    if 'histology' in samp_df.columns:
        print("\nMAPK and Panel by histology (mean):")
        g = samp_df.groupby('histology').agg({'mapk':['mean','count'],'panel':['mean','count']}).round(3)
        print(g)
        g.to_csv(OUT/"v10b_lee_by_histology.tsv", sep="\t")
        # Per-histology Spearman
        for h, sub in samp_df.groupby('histology'):
            if len(sub)>=10:
                r, pv = stats.spearmanr(sub['mapk'], sub['panel'])
                print(f"  {h}: n={len(sub)}, ρ={r:.3f}, p={pv:.3g}")
except Exception as e:
    print(f"\nMetadata merge skipped: {e}")

pd.DataFrame([{'cohort':'Lee_GSE213647','n':len(mapk_score),
               'mapk_n':len(available_mapk),'panel_n':len(available_panel),
               'spearman_rho':rho,'spearman_p':p,'pearson_r':pr,'pearson_p':pp}]
).to_csv(OUT/"v10b_lee_correlation.tsv", sep="\t", index=False)

# Add GSE126698 Korean PTC+HT n=28 with proper symbol matching
print("\n=== GSE126698 detail ===")
gse = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/GSE126698_rnaseq_expression_zscore.tsv", sep="\t", index_col=0)
g_mapk = [g for g in MAPK if g in gse.index]
g_panel = [g for g in PANEL if g in gse.index]
print(f"  GSE126698 n={gse.shape[1]}, MAPK={len(g_mapk)}, Panel={len(g_panel)}")
m_s = gse.loc[g_mapk].mean(axis=0)
p_s = gse.loc[g_panel].mean(axis=0)
r, pv = stats.spearmanr(m_s, p_s)
print(f"  ρ={r:.3f}, p={pv:.3g} (n={len(m_s)})")

# Combine all 3 cohort summary
print("\n=== Combined cross-cohort MAPK × 8-gene panel Spearman ===")
final = pd.DataFrame([
    {'cohort':'TCGA-THCA','n':572,'spearman_rho':-0.291,'spearman_p':2.5e-12},
    {'cohort':'Lee_GSE213647','n':len(mapk_score),'spearman_rho':rho,'spearman_p':p},
    {'cohort':'GSE126698_Korean_PTC_HT','n':len(m_s),'spearman_rho':r,'spearman_p':pv},
])
print(final.round(4).to_string(index=False))
final.to_csv(OUT/"v10_final_cross_cohort.tsv", sep="\t", index=False)

print("\n=== DONE v10b ===")
