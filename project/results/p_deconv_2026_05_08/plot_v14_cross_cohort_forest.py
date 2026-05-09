"""v14 forest plot — MAPK × Panel-8 / TDS-16 / TDS_8only Spearman ρ across 5 cohorts."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, matplotlib as mpl
from scipy import stats
mpl.rcParams['pdf.fonttype']=42; mpl.rcParams['ps.fonttype']=42; mpl.rcParams['font.family']='DejaVu Sans'

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")
df = pd.read_csv(OUT/"v14_cross_cohort_forest.tsv", sep="\t")

C_PANEL = '#2c7fb8'
C_TDS16 = '#d8453a'
C_8O    = '#7a4ca8'

cohort_label = {
    'TCGA-THCA':         'TCGA-THCA  (US, primary)',
    'Lee_GSE213647':     'Lee 2024 / GSE213647  (Korean, primary)',
    'GSE126698':         'GSE126698  (Korean PTC + HT)',
    'GSE286332':         'GSE286332  (Korean PTC vs PTC+HT)',
    'GSE76039_advanced': 'GSE76039  (PDTC + ATC, advanced)',
}
cohort_order = ['TCGA-THCA','Lee_GSE213647','GSE126698','GSE286332','GSE76039_advanced']

def pooled(rho_list, n_list):
    """Fixed-effect Fisher-z pool."""
    rho_arr = np.asarray(rho_list); n_arr = np.asarray(n_list)
    z = np.arctanh(rho_arr)
    se2 = 1.0 / (n_arr - 3)
    w = 1.0 / se2
    mz = (w * z).sum() / w.sum()
    mse = np.sqrt(1.0 / w.sum())
    return np.tanh(mz), np.tanh(mz - 1.96*mse), np.tanh(mz + 1.96*mse)

def heterogeneity(rho_list, n_list):
    """Cochran Q + I²."""
    rho_arr = np.asarray(rho_list); n_arr = np.asarray(n_list)
    z = np.arctanh(rho_arr)
    se2 = 1.0 / (n_arr - 3)
    w = 1.0 / se2
    mz = (w * z).sum() / w.sum()
    Q = (w * (z - mz)**2).sum()
    df_q = len(z) - 1
    I2 = max(0, (Q - df_q) / Q) * 100 if Q > 0 else 0
    p_q = 1 - stats.chi2.cdf(Q, df_q) if df_q > 0 else 1
    return Q, df_q, I2, p_q

# Build forest figure
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)

for ax, panel, color, title in [
    (axes[0], 'Panel_8',   C_PANEL, "8-gene Panel (deployable)"),
    (axes[1], 'TDS_16',    C_TDS16, "TDS-16 (Yoo 2014 canonical)"),
    (axes[2], 'TDS_8only', C_8O,    "TDS−panel (8 disjoint genes)"),
]:
    sub = df[df['panel']==panel].set_index('cohort').reindex(cohort_order)
    y = np.arange(len(sub))
    rho = sub['spearman_rho'].values
    lo  = sub['ci_lo'].values
    hi  = sub['ci_hi'].values
    n   = sub['n'].values

    for yi, r, l, h, ni in zip(y, rho, lo, hi, n):
        ax.plot([l, h], [yi, yi], color=color, lw=1.6)
        # square size scales with sqrt(n)
        ms = 6 + 12*np.sqrt(min(ni,700)/700)
        ax.plot(r, yi, 's', color=color, ms=ms, mec='black', mew=0.5)
        ax.text(0.62, yi, f"ρ={r:+.3f}  n={int(ni)}", va='center', fontsize=8.5,
                family='monospace', color='#222' if abs(r)<0.5 else color)

    # pooled
    pr, plo, phi = pooled(rho, n)
    Q, dfq, I2, pq = heterogeneity(rho, n)
    py = len(sub)
    ax.plot([plo, phi], [py, py], color='black', lw=2)
    ax.plot([plo, plo], [py-0.18, py+0.18], color='black', lw=2)
    ax.plot([phi, phi], [py-0.18, py+0.18], color='black', lw=2)
    ax.plot(pr, py, 'D', color='gold', mec='black', mew=0.6, ms=11)
    ax.text(0.62, py, f"POOL  ρ={pr:+.3f}  [{plo:+.2f},{phi:+.2f}]  I²={I2:.0f}%",
            va='center', fontsize=8.5, family='monospace', weight='bold')

    ax.axvline(0, color='gray', lw=0.5, ls='--')
    ax.set_xlim(-0.65, 1.4)
    ax.set_yticks(list(y) + [py])
    ax.set_yticklabels([cohort_label[c] for c in cohort_order] + ['POOLED (5 cohorts)'],
                       fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("MAPK output × thyroid-score Spearman ρ", fontsize=9)
    ax.set_title(title, fontsize=10.5, weight='bold', loc='left')
    ax.grid(axis='x', alpha=0.3)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

fig.suptitle("Figure SX_v14. Cross-cohort forest — MAPK output × thyroid-score Spearman ρ.\n"
             "Two-axis convergence (v12) predicts decoupling in HT-route cohorts (GSE286332) "
             "and saturation in advanced disease (GSE76039).",
             fontsize=12, weight='bold', y=1.02)

fig.savefig(OUT/"Fig_SX_v14_cross_cohort_forest.png", dpi=200, bbox_inches='tight')
fig.savefig(OUT/"Fig_SX_v14_cross_cohort_forest.pdf", bbox_inches='tight')
print(f"Saved {OUT}/Fig_SX_v14_cross_cohort_forest.{{png,pdf}}")

# Print summary
for panel in ['Panel_8','TDS_16','TDS_8only']:
    sub = df[df['panel']==panel].set_index('cohort').reindex(cohort_order)
    pr, plo, phi = pooled(sub['spearman_rho'].values, sub['n'].values)
    Q, dfq, I2, pq = heterogeneity(sub['spearman_rho'].values, sub['n'].values)
    print(f"  {panel:10s} pooled ρ = {pr:+.3f}  95% CI [{plo:+.3f}, {phi:+.3f}]  "
          f"Q = {Q:.2f} (df={dfq}, p={pq:.3g}), I² = {I2:.1f}%")
