"""v5 E — Pseudotime trajectory along DM1/DM2 score with cell-type fractions.

Two cohorts:
  • TCGA-THCA: bulk × canonical 8-gene score (rai_score_recalc) — diffusion / PCA-1 trajectory
  • Lee/GSE213647: panel_z continuous score across PTC + ATC + PDTC

For each cohort: rank samples by canonical 8-gene score → bin into deciles →
plot mean cell-type fraction per decile. The decile trajectory IS the
pseudotime ordering (one-dim, score-driven).
"""
from __future__ import annotations
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'DejaVu Sans'

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

celltypes = ['B cell','Endothelial cell','Epithelial cell','Fibroblast',
             'Malignant cell','Myeloid cell','NK cell','T cell']
N_BINS = 10

def trajectory(frac_df, score_df, label):
    """frac_df indexed by sample_id with cell-type cols; score_df has [sample_id, score]"""
    m = frac_df.merge(score_df, on='sample_id', how='inner').dropna(subset=['score'])
    m = m.sort_values('score').reset_index(drop=True)
    m['bin'] = pd.qcut(m['score'], q=N_BINS, labels=False, duplicates='drop')
    bin_means = m.groupby('bin').agg({c: 'mean' for c in celltypes if c in m.columns})
    bin_means['mean_score'] = m.groupby('bin')['score'].mean()
    bin_means['n'] = m.groupby('bin').size()
    bin_means.to_csv(OUT / f"v5E_trajectory_{label}.tsv", sep="\t")
    print(f"{label}: n={len(m)}, bins={len(bin_means)}")
    print(bin_means.round(3).to_string())
    return bin_means, m

# === TCGA: bulk fractions × canonical RAI score ===
print("=" * 70); print("v5 E TCGA trajectory"); print("=" * 70)
tcga_frac = pd.read_csv(OUT / "fractions_nu_SVR.tsv", sep="\t", index_col=0)
tcga_frac['sample_id'] = tcga_frac.index
canonical = pd.read_csv("/data/thca/repo_results/v17p3/tables/A2_dm_score_full_cohort.tsv", sep="\t")
canonical = canonical[canonical['dataset']=='TCGA-THCA'].copy()
canonical = canonical.rename(columns={'rai_score_recalc':'score'})
tcga_traj, tcga_full = trajectory(tcga_frac, canonical[['sample_id','score']], "TCGA")

# === Lee/GSE213647: panel_z continuous score ===
print("\n" + "=" * 70); print("v5 E Lee/GSE213647 trajectory"); print("=" * 70)
try:
    lee_frac = pd.read_csv(OUT / "fractions_LEE_nu_SVR.tsv", sep="\t", index_col=0)
    lee_frac['sample_id'] = lee_frac.index
    # Find Lee panel_z source
    lee_score_candidates = [
        "/data/thca/repo_results/v17p3/tables/A2_dm_score_full_cohort.tsv",
    ]
    lee_score = None
    for p in lee_score_candidates:
        try:
            df = pd.read_csv(p, sep="\t")
            for ds in ['Lee_GSE213647','Lee','GSE213647']:
                sub = df[df.get('dataset','')==ds]
                if len(sub)>10:
                    lee_score = sub.copy()
                    if 'panel_z' in sub.columns: lee_score = sub.rename(columns={'panel_z':'score'})
                    elif 'rai_score_recalc' in sub.columns: lee_score = sub.rename(columns={'rai_score_recalc':'score'})
                    print(f"Lee score loaded from {p} for dataset={ds} (n={len(lee_score)})")
                    break
            if lee_score is not None: break
        except Exception as e:
            print(f"failed {p}: {e}")
    if lee_score is None or 'sample_id' not in lee_score.columns or 'score' not in lee_score.columns:
        # Use per-sample panel_z file
        alt = "/home/seungho/personal/THCA_data_analysis/project/results/v17_korean/GSE213647_panel_score.tsv"
        if Path(alt).exists():
            ds = pd.read_csv(alt, sep="\t")
            ds = ds.rename(columns={'panel_z':'score'})
            # Match sample_id to lee_frac index — frac uses GSM IDs typically
            print(f"Lee per-sample panel_z columns: {list(ds.columns)}; n={len(ds)}")
            print(f"frac index head: {list(lee_frac.index[:5])}; ds sample_id head: {list(ds['sample_id'].head())}; gsm head: {list(ds['gsm'].head())}")
            for sid_col in ['sample_id','gsm','sample_name']:
                if sid_col in ds.columns:
                    inter = set(lee_frac['sample_id']) & set(ds[sid_col].astype(str))
                    print(f"  intersect with frac sample_id via ds.{sid_col}: {len(inter)}")
                    if len(inter) > 50:
                        sub = ds[[sid_col, 'score']].copy()
                        sub.columns = ['sample_id', 'score']
                        lee_score = sub
                        break
    if lee_score is not None and 'sample_id' in lee_score.columns and 'score' in lee_score.columns:
        lee_traj, lee_full = trajectory(lee_frac, lee_score[['sample_id','score']], "Lee")
    else:
        print("Lee trajectory SKIPPED — no per-sample score source")
        lee_traj = None
except Exception as e:
    print(f"Lee trajectory error: {e}")
    lee_traj = None

# === Plot ===
n_panels = 2 if lee_traj is not None else 1
fig, axes = plt.subplots(n_panels, 1, figsize=(11, 4.0 * n_panels), squeeze=False)

palette = {
    'B cell':'#9467bd', 'Endothelial cell':'#8c564b', 'Epithelial cell':'#2ca02c',
    'Fibroblast':'#7f7f7f', 'Malignant cell':'#d62728', 'Myeloid cell':'#ff7f0e',
    'NK cell':'#bcbd22', 'T cell':'#1f77b4',
}

for idx, (traj, label) in enumerate([(tcga_traj, 'TCGA-THCA'), (lee_traj, 'Lee/GSE213647')]):
    if traj is None: continue
    ax = axes[idx, 0]
    bins = traj.index.values
    score_means = traj['mean_score'].values
    ct_in = [c for c in celltypes if c in traj.columns]
    for ct in ct_in:
        ax.plot(score_means, traj[ct].values, marker='o', linewidth=1.5,
                color=palette.get(ct, '#666'), label=ct, markersize=4.5)
    ax.set_xlabel(f'{label} mean canonical 8-gene score (per decile)')
    ax.set_ylabel('Mean cell-type fraction')
    ax.set_title(f"{label} — pseudotime trajectory along DM1↔DM2 axis (n_bins={len(bins)})",
                 fontsize=10, loc='left')
    ax.legend(fontsize=7, ncols=2, loc='best')
    ax.grid(alpha=0.3)
    # annotate Spearman ρ score-vs-fraction in legend (top 3 |ρ|)
    rs = []
    for ct in ct_in:
        if abs(traj[ct].std()) > 1e-9:
            r, _ = stats.spearmanr(traj['mean_score'], traj[ct])
        else: r = 0
        rs.append((ct, r))
    rs_sorted = sorted(rs, key=lambda x: abs(x[1]), reverse=True)
    text = "Top-3 score×fraction Spearman ρ (decile):\n"
    for ct, r in rs_sorted[:3]:
        text += f"  {ct}: ρ={r:+.2f}\n"
    ax.text(0.99, 0.05, text, transform=ax.transAxes, fontsize=7,
            ha='right', va='bottom',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.8))

fig.suptitle("Supp Fig SX (v5 E) — Cell-type composition trajectory along canonical 8-gene score",
             fontsize=11, y=0.995)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(OUT / f"Fig_SX_deconvolution_v5_e_trajectory.{ext}", dpi=180, bbox_inches='tight')
    print(f"saved Fig_SX_deconvolution_v5_e_trajectory.{ext}")
plt.close(fig)
print("DONE — v5 e trajectory")
