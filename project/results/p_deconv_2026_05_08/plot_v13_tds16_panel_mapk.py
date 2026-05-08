"""v13 figure — TDS-16 × Panel-8 × MAPK (Paper 1 Fig 8 / SX_v13). 6-panel composite.

Panels:
  A. Cross-cohort MAPK × thyroid-score Spearman ρ (TCGA + Lee × Panel/TDS16/TDS8only)
  B. Per-driver-class TCGA score means: MAPK / Panel / TDS-16 across drivers
  C. sub-A vs sub-B d(A−B) on MAPK / HT / Panel / TDS-16 / TDS_8only
  D. Per-gene MAPK × TDS-16 Spearman ρ heatmap (16 genes × {TCGA, Lee})
  E. MAPK-decile pseudotime: Panel / TDS-16 / TDS_8only score per decile (TCGA + Lee)
  F. ROC-AUC: MAPK-high vs MAPK-low classifier — Panel vs TDS-16 vs TDS_8only
"""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, matplotlib as mpl
from matplotlib.patches import Patch
mpl.rcParams['pdf.fonttype']=42; mpl.rcParams['ps.fonttype']=42; mpl.rcParams['font.family']='DejaVu Sans'

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

corr  = pd.read_csv(OUT/"v13_cross_cohort_corr.tsv",       sep="\t")
dmean = pd.read_csv(OUT/"v13_per_driver_class_means.tsv",  sep="\t")
abdf  = pd.read_csv(OUT/"v13_subAB_three_panels.tsv",      sep="\t")
pg    = pd.read_csv(OUT/"v13_per_gene_mapk_corr.tsv",      sep="\t")
dec   = pd.read_csv(OUT/"v13_decile_trajectory.tsv",       sep="\t")
auc   = pd.read_csv(OUT/"v13_auc_panel_vs_tds16.tsv",      sep="\t")

C_PANEL = '#2c7fb8'   # 8-gene
C_TDS16 = '#d8453a'   # TDS-16 canonical
C_8O    = '#7a4ca8'   # TDS-8only (TDS−panel)
C_MAPK  = '#e8a31a'   # MAPK
C_HT    = '#3a78d8'   # HT

fig = plt.figure(figsize=(14, 10.5))
gs  = fig.add_gridspec(3, 3, hspace=0.55, wspace=0.45)

# ============================================================
# A. Cross-cohort ρ bars
# ============================================================
ax = fig.add_subplot(gs[0, 0])
panels = ['Panel_8','TDS_16','TDS_8only']
labels = ['8-gene\n(deployable)','TDS-16\n(canonical Yoo)','TDS−panel\n(8 only)']
cohorts = ['TCGA-THCA','Lee_GSE213647']
x = np.arange(len(panels)); w = 0.38
for i, c in enumerate(cohorts):
    sub = corr[corr['cohort']==c].set_index('panel')
    vals = [sub.loc[p,'spearman_rho'] for p in panels]
    bars = ax.bar(x + (i-0.5)*w, vals, w,
                  color=[C_PANEL, C_TDS16, C_8O],
                  edgecolor='black', linewidth=0.6,
                  hatch=['','////','xxxx'][1] if c.startswith('Lee') else None,
                  alpha=0.95 if i==0 else 0.65,
                  label=c.replace('_',' '))
    for b, v in zip(bars, vals):
        ax.text(b.get_x()+b.get_width()/2, v - 0.02, f"{v:+.2f}",
                ha='center', va='top', fontsize=8)
ax.axhline(0, color='black', lw=0.5)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
ax.set_ylabel("MAPK × thyroid-score Spearman ρ", fontsize=9)
ax.set_ylim(-0.55, 0.05)
ax.set_title("A. MAPK × thyroid-score ρ — 8-gene ≈ TDS-16\n(canonical and compact panel agree, both cohorts)",
             fontsize=10, loc='left')
ax.legend(['TCGA n=572','Lee n=632'], fontsize=8, loc='lower right',
          handles=[Patch(facecolor='gray', alpha=0.95, label='TCGA n=572'),
                   Patch(facecolor='gray', alpha=0.65, label='Lee n=632')])
ax.grid(axis='y', alpha=0.3)

# ============================================================
# B. Per-driver-class score means
# ============================================================
ax = fig.add_subplot(gs[0, 1])
order = ['BRAF_V600E','RAS_mut','RET_fusion','NTRK_fusion','driver_neg']
dm = dmean.set_index('driver').reindex(order)
xpos = np.arange(len(order))
ax.plot(xpos, dm['mapk_mean'].values,  marker='o', color=C_MAPK,  lw=2,   label='MAPK output z')
ax.plot(xpos, dm['panel_mean'].values, marker='s', color=C_PANEL, lw=2,   label='Panel-8 z')
ax.plot(xpos, dm['tds16_mean'].values, marker='^', color=C_TDS16, lw=2,   label='TDS-16 z')
ax.plot(xpos, dm['tds8only_mean'].values, marker='v', color=C_8O, lw=1.5, ls='--', label='TDS−panel z')
ax.axhline(0, color='black', lw=0.5)
ax.set_xticks(xpos); ax.set_xticklabels([f"{c}\n(n={int(dm.loc[c,'n'])})" for c in order],
                                        rotation=20, ha='right', fontsize=8)
ax.set_ylabel("Mean score z", fontsize=9)
ax.set_title("B. Per-driver-class TCGA score means\n(TDS-16 tracks Panel-8 across all drivers)",
             fontsize=10, loc='left')
ax.legend(fontsize=7, loc='upper right')
ax.grid(alpha=0.3)

# ============================================================
# C. sub-A vs sub-B d(A-B) — 5 metrics
# ============================================================
ax = fig.add_subplot(gs[0, 2])
order_c = ['mapk','ht','panel','tds16','tds8only']
labels_c = ['MAPK','HT','Panel-8','TDS-16','TDS−panel']
colors_c = [C_MAPK, C_HT, C_PANEL, C_TDS16, C_8O]
ab = abdf.set_index('metric').reindex(order_c)
xp = np.arange(len(order_c))
bars = ax.bar(xp, ab['d_A_minus_B'].values, color=colors_c, edgecolor='black', linewidth=0.6)
for b, v, p in zip(bars, ab['d_A_minus_B'].values, ab['p'].values):
    star = '***' if p<1e-3 else ('**' if p<1e-2 else ('*' if p<0.05 else 'NS'))
    ax.text(b.get_x()+b.get_width()/2, v + (0.05 if v>=0 else -0.12),
            f"{v:+.2f}\n{star}", ha='center', va='bottom' if v>=0 else 'top',
            fontsize=7.5)
ax.axhline(0, color='black', lw=0.5)
ax.set_xticks(xp); ax.set_xticklabels(labels_c, rotation=20, ha='right', fontsize=8)
ax.set_ylabel(f"Cohen's d  (sub-A n={int(ab['n_A'].iloc[0])} − sub-B n={int(ab['n_B'].iloc[0])})", fontsize=9)
ax.set_ylim(-0.4, 2.2)
ax.set_title("C. Two-axis convergence holds for TDS-16\n(MAPK splits sub-A/B; thyroid panels both flat)",
             fontsize=10, loc='left')
ax.grid(axis='y', alpha=0.3)

# ============================================================
# D. Per-gene MAPK × gene Spearman ρ heatmap
# ============================================================
ax = fig.add_subplot(gs[1, :2])
genes = pg['gene'].tolist()
mat   = pg[['rho_TCGA','rho_Lee']].values.T   # 2 × 16
im = ax.imshow(mat, aspect='auto', cmap='RdBu_r', vmin=-0.7, vmax=0.7)
ax.set_yticks([0,1]); ax.set_yticklabels(['TCGA n=572','Lee n=632'], fontsize=9)
ax.set_xticks(range(len(genes)))
xlab = []
for g in genes:
    in8 = bool(pg.loc[pg['gene']==g,'in_panel_8'].iloc[0])
    xlab.append(("★ "+g) if in8 else g)
ax.set_xticklabels(xlab, rotation=45, ha='right', fontsize=8.5)
for i in range(2):
    for j, g in enumerate(genes):
        v = mat[i, j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:+.2f}", ha='center', va='center',
                    fontsize=7, color='white' if abs(v)>0.45 else 'black')
ax.set_title("D. Per-gene MAPK × thyroid-gene Spearman ρ (★ = 8-gene panel members)\n"
             "Negative = MAPK output silences gene. 8-gene members representative of full TDS-16 behavior.",
             fontsize=10, loc='left')
cb = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
cb.set_label("Spearman ρ", fontsize=9)
cb.ax.tick_params(labelsize=8)

# ============================================================
# E. MAPK-decile pseudotime
# ============================================================
ax = fig.add_subplot(gs[1, 2])
for cohort, ls in [('TCGA-THCA','-'), ('Lee_GSE213647','--')]:
    d = dec[dec['cohort']==cohort].sort_values('decile')
    ax.plot(d['decile'], d['panel_mean'], marker='s', ls=ls, color=C_PANEL, lw=1.6,
            label=f"Panel-8 ({cohort.split('_')[0]})")
    ax.plot(d['decile'], d['tds16_mean'], marker='^', ls=ls, color=C_TDS16, lw=1.6,
            label=f"TDS-16 ({cohort.split('_')[0]})")
    ax.plot(d['decile'], d['tds8only_mean'], marker='v', ls=ls, color=C_8O, lw=1.2,
            alpha=0.7, label=f"TDS−panel ({cohort.split('_')[0]})")
ax.axhline(0, color='black', lw=0.5)
ax.set_xlabel("MAPK-output decile (low → high)", fontsize=9)
ax.set_ylabel("Mean thyroid-score z per decile", fontsize=9)
ax.set_title("E. MAPK pseudotime: thyroid scores monotone-decline\n(Panel-8 and TDS-16 trace identical curves)",
             fontsize=10, loc='left')
ax.legend(fontsize=6.5, loc='lower left', ncol=2)
ax.grid(alpha=0.3)

# ============================================================
# F. AUC bars: Panel vs TDS-16 vs TDS_8only — MAPK-high classifier
# ============================================================
ax = fig.add_subplot(gs[2, 0])
xp = np.arange(2); w = 0.27
auc = auc.set_index('cohort').reindex(['TCGA-THCA','Lee_GSE213647'])
bars1 = ax.bar(xp - w, auc['AUC_Panel'].values,    w, color=C_PANEL, edgecolor='black', lw=0.6, label='Panel-8')
bars2 = ax.bar(xp,     auc['AUC_TDS16'].values,    w, color=C_TDS16, edgecolor='black', lw=0.6, label='TDS-16')
bars3 = ax.bar(xp + w, auc['AUC_TDS8only'].values, w, color=C_8O,    edgecolor='black', lw=0.6, label='TDS−panel')
for bars in (bars1, bars2, bars3):
    for b in bars:
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.005,
                f"{b.get_height():.3f}", ha='center', va='bottom', fontsize=7.5)
ax.axhline(0.5, color='gray', lw=0.5, ls='--')
ax.set_xticks(xp); ax.set_xticklabels(['TCGA n=572','Lee n=632'], fontsize=9)
ax.set_ylabel("ROC-AUC for MAPK-high vs MAPK-low\n(score sign-flipped, low score = high MAPK)", fontsize=8.5)
ax.set_ylim(0.5, 0.82)
delta_strs = [f"Δ TDS-16 − Panel = {auc.loc[c,'delta_TDS16_vs_Panel']:+.3f}" for c in auc.index]
ax.set_title(f"F. ROC-AUC: 8-gene Panel vs TDS-16 (NS)\n{delta_strs[0]} (TCGA), {delta_strs[1]} (Lee)",
             fontsize=10, loc='left')
ax.legend(fontsize=8, loc='lower right')
ax.grid(axis='y', alpha=0.3)

# ============================================================
# G. Methodology / interpretation panel
# ============================================================
ax = fig.add_subplot(gs[2, 1:])
ax.axis('off')
text = (
    "v13 — TDS-16 × 8-gene panel × MAPK (Paper 1 Fig 8 supporting analysis, 2026-05-09)\n\n"
    "Question. Is the MAPK→thyroid-silencing axis (v9–v12) specific to the deployable 8-gene\n"
    "compact readout, or does it operate on the full canonical Yoo 2014 TDS-16? If both panels\n"
    "behave identically, the 8-gene paper's compact-readout claim is bulletproof for the\n"
    "Reviewer Q 'why 8 genes not TDS-16 — cherry-picked from a 16-gene set?'.\n\n"
    "Methods. TCGA-THCA n=572 + Lee/GSE213647 n=632 z-scored bulk RNA-seq. Score = within-cohort\n"
    "z-mean of: MAPK output (DUSP4/5/6, SPRY2/4, ETV4/5, PHLDA1, CCND1; n=9), 8-gene panel\n"
    "(DIO1/FOXE1/NKX2-1/PAX8/SLC5A5/TG/TPO/TSHR), TDS-16 (Yoo 2014: + DIO2/DUOX1/DUOX2/GLIS3/\n"
    "SLC26A4/SLC5A8/THRA/THRB), TDS−panel (TDS-16 \\ Panel-8, n=8 disjoint genes), HT signature\n"
    "(HLA-II + B-cell + AICDA/CXCL13/IFNG, n=15). Driver class from cBioPortal SV/MAF anchors;\n"
    "DM1 sub-A/sub-B labels from d6p7. Mann-Whitney p where applicable.\n\n"
    "Findings.\n"
    "  • MAPK × Panel-8 ρ = −0.291 (TCGA) / −0.395 (Lee); MAPK × TDS-16 ρ = −0.306 / −0.435 —\n"
    "    indistinguishable across cohorts. TDS−panel disjoint set ρ = −0.310 / −0.449 (D).\n"
    "  • Per-driver d (BRAF V600E vs RAS): Panel-8 d = −1.615, TDS-16 d = −1.616 — identical (B).\n"
    "  • Sub-A vs sub-B convergence: MAPK d = +1.79 (p=4×10⁻¹⁷), TDS-16 d = +0.03 NS, Panel-8\n"
    "    d = +0.07 NS. Both compact and canonical panels are flat across the MAPK split — same\n"
    "    two-axis convergence pattern as v12, on the full TDS-16 (C).\n"
    "  • MAPK-high vs MAPK-low ROC-AUC: Panel-8 = 0.623 / 0.728, TDS-16 = 0.631 / 0.740, ΔAUC\n"
    "    +0.007 / +0.012 NS. Decile pseudotime (E) shows Panel-8 and TDS-16 trace identical\n"
    "    monotonic descent as MAPK-output rises.\n\n"
    "Verdict. The MAPK→silencing axis is a property of the canonical thyroid-differentiation\n"
    "program, not of the deployable 8-gene subset alone. Panel-8 is a compact lossless readout\n"
    "of the same biology — directly supporting the Paper 1 §2 / Reviewer-Q3 framing that 8-gene\n"
    "captures TDS-16 axis without cherry-pick (ΔAUC NS in both cohorts; per-gene heatmap shows\n"
    "8-gene members as median-rank within the 16, not selectively top-extreme)."
)
ax.text(0, 1, text, va='top', ha='left', fontsize=8.4, family='monospace',
        bbox=dict(facecolor='#f5f5f5', edgecolor='black', linewidth=0.4, pad=8))

fig.suptitle("Figure SX_v13. The MAPK→silencing axis operates on canonical TDS-16 = 8-gene compact panel.",
             fontsize=12, weight='bold', y=0.995)
fig.savefig(OUT/"Fig_SX_v13_TDS16_MAPK.png", dpi=200, bbox_inches='tight')
fig.savefig(OUT/"Fig_SX_v13_TDS16_MAPK.pdf", bbox_inches='tight')
print(f"Saved {OUT}/Fig_SX_v13_TDS16_MAPK.{{png,pdf}}")
