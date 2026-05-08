"""v12: MAPK-axis vs HT-axis convergence to 8-gene silencing — two-path proof."""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
import matplotlib.pyplot as plt, matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42; mpl.rcParams['ps.fonttype']=42; mpl.rcParams['font.family']='DejaVu Sans'

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

MAPK = ['DUSP4','DUSP5','DUSP6','SPRY2','SPRY4','ETV4','ETV5','PHLDA1','CCND1']
PANEL = ['DIO1','FOXE1','NKX2-1','PAX8','SLC5A5','TG','TPO','TSHR']
# Hashimoto / B-cell-mediated immune signature
HT = ['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1',
      'CD79A','CD79B','MS4A1','AICDA','IGHM','IGKC','CXCL13','CCR6','IFNG']

def patient_id(s): return "-".join(s.split("-")[:3])

tcga = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv", sep="\t", index_col=0)
mapk_present = [g for g in MAPK if g in tcga.index]
panel_present = [g for g in PANEL if g in tcga.index]
ht_present = [g for g in HT if g in tcga.index]
print(f"MAPK={len(mapk_present)}/{len(MAPK)}, Panel={len(panel_present)}/{len(PANEL)}, HT={len(ht_present)}/{len(HT)}")
print(f"  HT hits: {ht_present}")

mapk_score = tcga.loc[mapk_present].mean(axis=0)
panel_score = tcga.loc[panel_present].mean(axis=0)
ht_score = tcga.loc[ht_present].mean(axis=0)

# Merge with sub-A/B
sub = pd.read_csv("/data/thca/repo_results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv", sep="\t")
sub.columns = ['sample_id','sub_cluster']
sub['patient'] = sub['sample_id'].apply(patient_id)

samp = pd.DataFrame({'sample':mapk_score.index,'mapk':mapk_score.values,
                     'panel':panel_score.values,'ht':ht_score.values})
samp['patient'] = samp['sample'].apply(patient_id)
m = samp.merge(sub, on='patient', how='inner')
print(f"\nsubA/B × tri-score: {len(m)}")

def cd(x,y):
    x,y=np.asarray(x),np.asarray(y)
    if len(x)<2 or len(y)<2: return 0
    p = np.sqrt(((len(x)-1)*x.var(ddof=1)+(len(y)-1)*y.var(ddof=1))/(len(x)+len(y)-2))
    return (x.mean()-y.mean())/p if p>0 else 0

a = m[m['sub_cluster']=='sub_A']; b = m[m['sub_cluster']=='sub_B']
print(f"\n  sub-A (n={len(a)}) vs sub-B (n={len(b)}):")
print(f"    MAPK  d(A-B) = {cd(a['mapk'].values, b['mapk'].values):+.3f}  (sub-A MAPK-active)")
print(f"    HT    d(A-B) = {cd(a['ht'].values,   b['ht'].values):+.3f}  (sub-B HT-active expected)")
print(f"    Panel d(A-B) = {cd(a['panel'].values, b['panel'].values):+.3f}  (similar — convergence)")
_, p_mapk = stats.mannwhitneyu(a['mapk'], b['mapk'])
_, p_ht = stats.mannwhitneyu(a['ht'], b['ht'])
_, p_pan = stats.mannwhitneyu(a['panel'], b['panel'])
print(f"    p MAPK={p_mapk:.3g}, HT={p_ht:.3g}, Panel={p_pan:.3g}")

# Whole-cohort: per-sample dual-axis classification
m['axis_class'] = pd.cut(m['mapk'], bins=[-np.inf, m['mapk'].median(), np.inf], labels=['mapk_low','mapk_high'])
m['ht_class']   = pd.cut(m['ht'],   bins=[-np.inf, m['ht'].median(),   np.inf], labels=['ht_low','ht_high'])
print("\n  Cross-tab axis × sub:")
ct = pd.crosstab(m['sub_cluster'], [m['axis_class'], m['ht_class']])
print(ct)

# Output
m.to_csv(OUT/"v12_subAB_mapk_ht_panel.tsv", sep="\t", index=False)
pd.DataFrame([
    {'metric':'MAPK','d_A_minus_B':cd(a['mapk'].values,b['mapk'].values),'p':p_mapk,'n_A':len(a),'n_B':len(b)},
    {'metric':'HT',  'd_A_minus_B':cd(a['ht'].values,b['ht'].values),'p':p_ht,  'n_A':len(a),'n_B':len(b)},
    {'metric':'Panel','d_A_minus_B':cd(a['panel'].values,b['panel'].values),'p':p_pan,'n_A':len(a),'n_B':len(b)},
]).to_csv(OUT/"v12_subAB_three_axes_d.tsv", sep="\t", index=False)

# === Figure: 2-axis scatter ===
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

ax = axes[0]
for s, color, label in [('sub_A','#d8453a','sub-A (driver+)'),('sub_B','#3a78d8','sub-B (HT-overlap)')]:
    sub_m = m[m['sub_cluster']==s]
    ax.scatter(sub_m['mapk'], sub_m['ht'], s=20, alpha=0.55, color=color, edgecolor='none',
               label=f"{label} n={len(sub_m)}")
ax.axhline(0, color='gray', lw=0.5); ax.axvline(0, color='gray', lw=0.5)
ax.set_xlabel('MAPK output z (DUSP/SPRY/ETV/PHLDA1/CCND1)')
ax.set_ylabel('HT/B-cell signature z (HLA-II/CD79/AICDA/CXCL13)')
ax.set_title('A. Two convergent axes to 8-gene silencing\nsub-A = MAPK-axis · sub-B = HT-axis', fontsize=10, loc='left')
ax.legend(fontsize=8, loc='upper right')
ax.grid(alpha=0.3)

ax = axes[1]
for s, color, label in [('sub_A','#d8453a','sub-A'),('sub_B','#3a78d8','sub-B')]:
    sub_m = m[m['sub_cluster']==s]
    ax.scatter(sub_m['mapk'], sub_m['panel'], s=20, alpha=0.55, color=color, edgecolor='none',
               label=f"{label} n={len(sub_m)}")
ax.axhline(0, color='gray', lw=0.5); ax.axvline(0, color='gray', lw=0.5)
ax.set_xlabel('MAPK output z')
ax.set_ylabel('8-gene panel z')
ax.set_title('B. Panel suppressed regardless of MAPK\n(both sub-A high-MAPK and sub-B low-MAPK)', fontsize=10, loc='left')
ax.legend(fontsize=8, loc='upper right')
ax.grid(alpha=0.3)

ax = axes[2]
for s, color, label in [('sub_A','#d8453a','sub-A'),('sub_B','#3a78d8','sub-B')]:
    sub_m = m[m['sub_cluster']==s]
    ax.scatter(sub_m['ht'], sub_m['panel'], s=20, alpha=0.55, color=color, edgecolor='none',
               label=f"{label} n={len(sub_m)}")
ax.axhline(0, color='gray', lw=0.5); ax.axvline(0, color='gray', lw=0.5)
ax.set_xlabel('HT/B-cell signature z')
ax.set_ylabel('8-gene panel z')
ax.set_title('C. Panel suppressed regardless of HT axis too\n→ panel = downstream sink of two upstream paths', fontsize=10, loc='left')
ax.legend(fontsize=8, loc='upper right')
ax.grid(alpha=0.3)

fig.suptitle('Supp Fig SX (v12) — DM1 sub-A/B convergent two-axis model: MAPK + HT both reach 8-gene silencing',
             fontsize=11, y=1.02)
fig.tight_layout()
for ext in ('png','pdf'):
    fig.savefig(OUT/f"Fig_SX_v12_two_axis.{ext}", dpi=180, bbox_inches='tight')
plt.close(fig)
print(f"\nFigure saved: Fig_SX_v12_two_axis.{{png,pdf}}")
print("\n=== DONE v12 ===")
