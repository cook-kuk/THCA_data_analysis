"""v9+v10 figure — MAPK pathway → 8-gene panel suppression mechanism."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, matplotlib as mpl
from scipy import stats
mpl.rcParams['pdf.fonttype']=42; mpl.rcParams['ps.fonttype']=42; mpl.rcParams['font.family']='DejaVu Sans'
OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

# Re-build TCGA scores for scatter
MAPK = ['DUSP4','DUSP5','DUSP6','SPRY2','SPRY4','ETV4','ETV5','PHLDA1','CCND1']
PANEL = ['DIO1','FOXE1','NKX2-1','PAX8','SLC5A5','TG','TPO','TSHR']

tcga = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv", sep="\t", index_col=0)
tcga_mapk = tcga.loc[MAPK].mean(axis=0)
tcga_panel = tcga.loc[PANEL].mean(axis=0)

lee_z = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/GSE213647_rnaseq_expression_zscore.tsv", sep="\t", index_col=0)
gmap = pd.read_csv("/data/thca/repo_results/v17p3/tables/F1_gene_recovery_mapping.tsv", sep="\t")
ens_to_sym = dict(zip(gmap['ensembl'], gmap['symbol']))
target = list(set(MAPK+PANEL))
ens_target = [e for e,s in ens_to_sym.items() if s in target]
lee_sub = lee_z.loc[lee_z.index.isin(ens_target)].copy()
lee_sub.index = [ens_to_sym[e] for e in lee_sub.index]
lee_sub = lee_sub.groupby(lee_sub.index).mean()
lee_mapk = lee_sub.loc[[g for g in MAPK if g in lee_sub.index]].mean(axis=0)
lee_panel = lee_sub.loc[[g for g in PANEL if g in lee_sub.index]].mean(axis=0)
lee_meta = pd.read_csv("/data/thca/repo_results/v17_korean/GSE213647_panel_score.tsv", sep="\t")
lee_df = pd.DataFrame({'sample':lee_mapk.index,'mapk':lee_mapk.values,'panel':lee_panel.values}).merge(
    lee_meta[['gsm','histology']], left_on='sample', right_on='gsm', how='inner')

# v9 driver methylation pivot
v9 = pd.read_csv(OUT/"v9_mapk_methylation_pivot.tsv", sep="\t", index_col=0)

fig = plt.figure(figsize=(15, 10))
gs = fig.add_gridspec(2, 3, hspace=0.4, wspace=0.35)

# Panel A: TCGA driver-class mean β bar
axA = fig.add_subplot(gs[0,0])
drivers = ['BRAF_V600E','RET_fusion','NTRK_fusion','RAS_mut','driver_neg']
ns = {'BRAF_V600E':290,'RET_fusion':33,'NTRK_fusion':10,'RAS_mut':53,'driver_neg':88}
betas = {'BRAF_V600E':0.372,'RET_fusion':0.385,'NTRK_fusion':None,'RAS_mut':0.272,'driver_neg':None}
# baseline = driver_neg overall mean (computed: ~0.30 typical)
clin = pd.read_csv("/data/thca/repo_results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv", sep="\t")
meth = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round5/r5_2_sample_methylation_8gene.tsv", sep="\t")
clin['patient'] = clin['tcga_short']; meth['patient'] = meth['sample_short']
mm = clin.merge(meth, on='patient', how='inner')
def patient_id(s): return "-".join(s.split("-")[:3])
cbio = pd.read_csv("/data/thca/repo_results/audit_2026_04_30/round3/cbio_sv_thca.tsv", sep="\t")
cbio['patient'] = cbio['sampleId'].apply(patient_id)
def fc(s):
    s=str(s) if pd.notna(s) else ''
    if 'RET' in s.upper(): return 'RET_fusion'
    if 'NTRK' in s.upper(): return 'NTRK_fusion'
    return None
cbio['fusion_class'] = cbio['eventInfo'].apply(fc)
fusion = cbio.dropna(subset=['fusion_class']).groupby('patient')['fusion_class'].first().reset_index()
mm = mm.merge(fusion, on='patient', how='left')
def driver(r):
    if r['has_braf_v600e']: return 'BRAF_V600E'
    if r['has_ras_mut']: return 'RAS_mut'
    if r['fusion_class']=='RET_fusion': return 'RET_fusion'
    if r['fusion_class']=='NTRK_fusion': return 'NTRK_fusion'
    return 'driver_neg'
mm['driver'] = mm.apply(driver, axis=1)
mean_betas = mm.groupby('driver')['mean_8g_beta'].agg(['mean','std','count']).reindex(drivers).fillna(0)
colors = {'BRAF_V600E':'#d8453a','RET_fusion':'#e08a3a','NTRK_fusion':'#e0bd3a','RAS_mut':'#3a78d8','driver_neg':'#888'}
xpos = np.arange(len(drivers))
axA.bar(xpos, mean_betas['mean'].values, yerr=mean_betas['std'].values/np.sqrt(mean_betas['count'].values),
        color=[colors[d] for d in drivers], edgecolor='black', linewidth=0.6, capsize=3)
axA.set_xticks(xpos); axA.set_xticklabels([f"{d}\n(n={int(mean_betas.loc[d,'count'])})" for d in drivers], fontsize=8)
axA.set_ylabel('mean 8-gene β (HM450)')
axA.set_title('A. Driver-class × 8-gene methylation\nMAPK-active (BRAF/RET/NTRK) hyper-methylate; RAS preserves', fontsize=9, loc='left')
axA.axhline(mean_betas.loc['driver_neg','mean'], ls='--', color='gray', alpha=0.5, lw=0.5)
axA.grid(axis='y', alpha=0.3)
# annotate d
axA.annotate(f"d(BRAF−RET)=−0.26", (0.5, mean_betas.loc['BRAF_V600E','mean']+0.04), fontsize=7, ha='center')
axA.annotate(f"d(BRAF−RAS)=+1.97", (1.5, 0.42), fontsize=7, ha='center', color='#a02020')

# Panel B: TCGA MAPK × Panel scatter
axB = fig.add_subplot(gs[0,1])
axB.scatter(tcga_mapk, tcga_panel, s=4, alpha=0.4, c='#3a78d8')
slope, intercept, _, _, _ = stats.linregress(tcga_mapk, tcga_panel)
xs = np.linspace(tcga_mapk.min(), tcga_mapk.max(), 100)
axB.plot(xs, slope*xs+intercept, 'r-', lw=1)
rho, p = stats.spearmanr(tcga_mapk, tcga_panel)
axB.set_xlabel('MAPK output mean z (DUSP/SPRY/ETV/PHLDA1/CCND1)')
axB.set_ylabel('8-gene panel mean z')
axB.set_title(f'B. TCGA-THCA n={len(tcga_mapk)}\nMAPK × Panel  ρ={rho:+.3f}  p={p:.2g}', fontsize=9, loc='left')
axB.grid(alpha=0.3)
axB.axhline(0, color='black', lw=0.4); axB.axvline(0, color='black', lw=0.4)

# Panel C: Lee MAPK × Panel by histology
axC = fig.add_subplot(gs[0,2])
hist_colors = {'Normal':'#3a78d8','PTC':'#d8453a','PDFP':'#a035a0','UTC/ATC':'#000000'}
for h, sub in lee_df.groupby('histology'):
    axC.scatter(sub['mapk'], sub['panel'], s=8, alpha=0.55, color=hist_colors.get(h,'gray'),
                label=f"{h} (n={len(sub)})", edgecolor='none')
slope, intercept, _, _, _ = stats.linregress(lee_df['mapk'], lee_df['panel'])
xs = np.linspace(lee_df['mapk'].min(), lee_df['mapk'].max(), 100)
axC.plot(xs, slope*xs+intercept, 'k--', lw=1, alpha=0.7)
rho, p = stats.spearmanr(lee_df['mapk'], lee_df['panel'])
axC.set_xlabel('MAPK output mean z')
axC.set_ylabel('8-gene panel mean z')
axC.set_title(f'C. Lee GSE213647 n={len(lee_df)}\nMAPK × Panel  ρ={rho:+.3f}  p={p:.2g}', fontsize=9, loc='left')
axC.legend(fontsize=7, loc='upper right')
axC.grid(alpha=0.3); axC.axhline(0, color='black', lw=0.4); axC.axvline(0, color='black', lw=0.4)

# Panel D: TCGA driver MAPK score bars
axD = fig.add_subplot(gs[1,0])
ds = pd.read_csv(OUT/"v10_driver_mapk_panel_means.tsv", sep="\t", header=[0,1])
print(ds.columns.tolist()[:8])
v10 = pd.read_csv(OUT/"v10_driver_mapk_panel_means.tsv", sep="\t", header=[0,1])
v10.columns = ['driver','mapk_mean','mapk_std','mapk_count','panel_mean','panel_std']
v10 = v10.set_index('driver').reindex(drivers).fillna(0)
xpos = np.arange(len(drivers))
w = 0.4
axD.bar(xpos-w/2, v10['mapk_mean'].values, yerr=v10['mapk_std'].values/np.sqrt(v10['mapk_count'].values),
        width=w, color='#d8453a', alpha=0.8, edgecolor='black', linewidth=0.5, capsize=3, label='MAPK output')
axD.bar(xpos+w/2, v10['panel_mean'].values, yerr=v10['panel_std'].values/np.sqrt(v10['mapk_count'].values),
        width=w, color='#3a78d8', alpha=0.8, edgecolor='black', linewidth=0.5, capsize=3, label='8-gene panel')
axD.set_xticks(xpos); axD.set_xticklabels([f"{d}\n(n={int(v10.loc[d,'mapk_count'])})" for d in drivers], fontsize=8)
axD.set_ylabel('mean z-score')
axD.set_title('D. Driver-class MAPK output vs 8-gene panel (TCGA)\nInverse pattern: MAPK↑ → Panel↓', fontsize=9, loc='left')
axD.legend(fontsize=8); axD.grid(axis='y', alpha=0.3)
axD.axhline(0, color='black', lw=0.4)

# Panel E: histology trajectory
axE = fig.add_subplot(gs[1,1])
hist_order = ['Normal','PTC','PDFP','UTC/ATC']
hd = lee_df.groupby('histology').agg({'mapk':['mean','std','count'],'panel':['mean','std']})
hd.columns = ['mapk_m','mapk_s','n','panel_m','panel_s']
hd = hd.reindex(hist_order)
xpos = np.arange(len(hist_order))
axE.errorbar(xpos, hd['mapk_m'], yerr=hd['mapk_s']/np.sqrt(hd['n']), fmt='-o', color='#d8453a', label='MAPK output', capsize=3)
axE.errorbar(xpos, hd['panel_m'], yerr=hd['panel_s']/np.sqrt(hd['n']), fmt='-o', color='#3a78d8', label='8-gene panel', capsize=3)
axE.set_xticks(xpos); axE.set_xticklabels([f"{h}\n(n={int(hd.loc[h,'n'])})" for h in hist_order], fontsize=8)
axE.set_ylabel('mean z-score')
axE.set_title('E. Lee histology trajectory\nMAPK rises in PTC; Panel monotonically falls Normal→ATC', fontsize=9, loc='left')
axE.axhline(0, color='black', lw=0.4)
axE.legend(fontsize=8); axE.grid(alpha=0.3)

# Panel F: text
axF = fig.add_subplot(gs[1,2])
axF.axis('off')
text = (
    "F. MAPK → 8-gene silencing — mechanism\n"
    "─────────────────────────────────\n"
    "Driver class × HM450 methylation (TCGA n=478):\n"
    "  BRAF V600E β=0.37, RET fusion β=0.39\n"
    "  RAS mut    β=0.27, driver-neg β baseline\n"
    "  d(BRAF−RET)=−0.26 (essentially equal)\n"
    "  d(BRAF−RAS)=+1.97, d(RET−RAS)=+2.17\n"
    "  → MAPK-output, NOT BRAF-specific\n"
    "\n"
    "RNA cross-cohort MAPK × Panel:\n"
    "  TCGA n=572:  ρ=−0.291  p=2.5e-12\n"
    "  Lee  n=632:  ρ=−0.395  p=4.8e-25\n"
    "  → bi-cohort replication\n"
    "\n"
    "Lee histology trajectory:\n"
    "  Normal (n=262): MAPK=−0.48, Panel=+0.45\n"
    "  PTC    (n=353): MAPK=+0.36, Panel=−0.28\n"
    "  ATC    (n=8):   Panel collapses to −1.57\n"
    "\n"
    "Mechanism: MAPK pathway output (DUSP/SPRY/\n"
    "ETV transcriptional program) → DNA-methylation\n"
    "reprogramming at thyrocyte-lineage TF promoters\n"
    "→ 8-gene silencing → DM1 phenotype.\n"
    "\n"
    "Why RAS retains panel: RAS in thyroid is a weak\n"
    "MEK/ERK activator (vs strong PI3K) — known\n"
    "biology; explains FVPTC histology overlap."
)
axF.text(0.0, 1.0, text, fontsize=7.5, va='top', family='DejaVu Sans Mono', transform=axF.transAxes)

fig.suptitle('Supp Fig SX (v9/v10) — MAPK-pathway-driven 8-gene silencing across drivers, methylation, RNA',
             fontsize=11, y=0.995)

for ext in ('png','pdf'):
    fig.savefig(OUT/f"Fig_SX_v9v10_mapk_mechanism.{ext}", dpi=180, bbox_inches='tight')
plt.close(fig)
print("Saved v9/v10 mechanism figure")
