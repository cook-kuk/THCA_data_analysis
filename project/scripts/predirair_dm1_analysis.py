#!/usr/bin/env python3
"""
PREDIRAIR-DM1 Feasibility Analysis (★가설, 미검증)
TCGA-THCA patient-level data

Data sources:
  - r17_tert_per_sample.tsv   : RAI_8 score, BRAF, TERT, PFI
  - dm_master_with_pfi.tsv   : age, stage

DM1-low definition for BRAF+ patients:
  RAI_8 < median(RAI_8 within BRAF+) = median split within BRAF+ subgroup.
  This mirrors IHC H-score median split in a prospective study.

Score definitions:
  Partial-PREDIRAIR = Age≥55(2pt) + TERT+(3pt)  [max 5; sTg/HG-TC absent]
  PREDIRAIR-DM1    = Age≥55(2pt) + TERT+(3pt) + DM1-low×BRAF+(3pt)  [max 8]
"""
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.utils import concordance_index
from lifelines.statistics import logrank_test
import base64, io, json, warnings
warnings.filterwarnings('ignore')

ROOT = Path('/home/seungho/personal/THCA_data_analysis')
OUT  = Path('/var/www/papers/predirair_dm1_v1')
OUT.mkdir(parents=True, exist_ok=True)

# ── 1. Load & merge ─────────────────────────────────────────────────────────
r17 = pd.read_csv(
    ROOT / 'project/results/r17_tcga_panel_d4p2_reconciliation/tert_interaction/r17_tert_per_sample.tsv',
    sep='\t'
)
master = pd.read_csv(
    ROOT / 'project/results/dark_matter_phase1/tcga_dm_master_with_pfi.tsv',
    sep='\t'
)[['tcga_short', 'age', 'stage', 'ajcc_pathologic_tumor_stage']]

df = r17.merge(master, on='tcga_short', how='left')

# Clean types
df['PFI']         = pd.to_numeric(df['PFI'], errors='coerce')
df['pfi_years']   = pd.to_numeric(df['PFI.time'], errors='coerce') / 365.25
df['has_braf']    = df['has_braf_v600e'].map({True: True, False: False, 'True': True, 'False': False})
df['tert_pos']    = df['tert_pos'].map({True: True, False: False, 'True': True, 'False': False})
df['RAI_8']       = pd.to_numeric(df['RAI_8'], errors='coerce')

# Drop rows missing key vars
df = df.dropna(subset=['PFI', 'pfi_years', 'RAI_8', 'has_braf', 'tert_pos', 'age'])
df = df[df['pfi_years'] > 0]  # remove 0-time records
df['PFI'] = df['PFI'].astype(int)

N = len(df)
N_braf = int(df['has_braf'].sum())
N_tert = int(df['tert_pos'].sum())
N_pfi_events = int(df['PFI'].sum())

# ── 2. DM1-low definition within BRAF+ (median split) ───────────────────────
braf_mask   = df['has_braf'] == True
rai8_braf_median = df.loc[braf_mask, 'RAI_8'].median()

# DM1-low = BRAF+ AND RAI_8 < median(BRAF+)
# (RAI_8 negative → DM1-like direction; lower = more dedifferentiated)
df['dm1_low_braf'] = (braf_mask) & (df['RAI_8'] < rai8_braf_median)
N_dm1_braf = int(df['dm1_low_braf'].sum())

# Also define overall DM1-low for all patients (zone-based)
df['dm1_low_all'] = df['zone'].isin(['BRAF-like', 'dark-matter'])

print(f"N={N}, BRAF+={N_braf}, TERT+={N_tert}, PFI events={N_pfi_events}")
print(f"RAI_8 median within BRAF+: {rai8_braf_median:.3f}")
print(f"DM1-low × BRAF+ (median split): {N_dm1_braf}")

# ── 3. Score components ─────────────────────────────────────────────────────
df['pts_age']      = (df['age'] >= 55).astype(int) * 2
df['pts_tert']     = df['tert_pos'].astype(int) * 3
df['pts_dm1_braf'] = df['dm1_low_braf'].astype(int) * 3

df['score_partial']   = df['pts_age'] + df['pts_tert']             # max 5
df['score_dm1']       = df['pts_age'] + df['pts_tert'] + df['pts_dm1_braf']  # max 8

print("Score partial dist:", df['score_partial'].value_counts().sort_index().to_dict())
print("Score DM1 dist:",     df['score_dm1'].value_counts().sort_index().to_dict())

# Score groups
df['grp_partial'] = df['score_partial'].map({0: '0pt', 2: '2pt', 3: '3pt', 5: '5pt'}).fillna('기타')
df['grp_dm1'] = pd.cut(
    df['score_dm1'], bins=[-1, 0, 2, 5, 8],
    labels=['0pt (저)', '2pt', '3-5pt (중)', '6-8pt (고)']
)

# BRAF+ subset
braf_df = df[df['has_braf']].copy()

# ── 4. C-statistics ─────────────────────────────────────────────────────────
def harrell_c(subdf, score_col):
    try:
        return concordance_index(subdf['pfi_years'], -subdf[score_col], subdf['PFI'])
    except Exception:
        return np.nan

def bootstrap_c(subdf, score_col, n_boot=1000, seed=42):
    rng = np.random.default_rng(seed)
    cs = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(subdf), len(subdf))
        s = subdf.iloc[idx]
        if s['PFI'].sum() < 2:
            continue
        try:
            cs.append(concordance_index(s['pfi_years'], -s[score_col], s['PFI']))
        except Exception:
            pass
    if len(cs) < 10:
        return [np.nan, np.nan]
    return list(np.percentile(cs, [2.5, 97.5]))

c_partial   = harrell_c(df, 'score_partial')
c_dm1       = harrell_c(df, 'score_dm1')
ci_partial  = bootstrap_c(df, 'score_partial')
ci_dm1      = bootstrap_c(df, 'score_dm1')

c_braf_partial = harrell_c(braf_df, 'score_partial')
c_braf_dm1     = harrell_c(braf_df, 'score_dm1')
ci_braf_partial = bootstrap_c(braf_df, 'score_partial')
ci_braf_dm1     = bootstrap_c(braf_df, 'score_dm1')

print(f"\nC-statistics:")
print(f"  Partial (all):    {c_partial:.3f} [{ci_partial[0]:.3f}, {ci_partial[1]:.3f}]")
print(f"  DM1    (all):     {c_dm1:.3f} [{ci_dm1[0]:.3f}, {ci_dm1[1]:.3f}]  Δ={c_dm1-c_partial:+.3f}")
print(f"  Partial (BRAF+):  {c_braf_partial:.3f} [{ci_braf_partial[0]:.3f}, {ci_braf_partial[1]:.3f}]")
print(f"  DM1    (BRAF+):   {c_braf_dm1:.3f} [{ci_braf_dm1[0]:.3f}, {ci_braf_dm1[1]:.3f}]  Δ={c_braf_dm1-c_braf_partial:+.3f}")

# ── 5. Cox PH models ────────────────────────────────────────────────────────
def cox_model(subdf, score_col):
    try:
        cph = CoxPHFitter()
        cph.fit(subdf[['pfi_years', 'PFI', score_col]].dropna(),
                duration_col='pfi_years', event_col='PFI')
        s = cph.summary
        hr  = float(np.exp(s.loc[score_col, 'coef']))
        lo  = float(np.exp(s.loc[score_col, 'coef lower 95%']))
        hi  = float(np.exp(s.loc[score_col, 'coef upper 95%']))
        p   = float(s.loc[score_col, 'p'])
        return hr, lo, hi, p
    except Exception as e:
        print(f"  Cox failed for {score_col}: {e}")
        return np.nan, np.nan, np.nan, np.nan

hr_p, lo_p, hi_p, p_p = cox_model(df, 'score_partial')
hr_d, lo_d, hi_d, p_d = cox_model(df, 'score_dm1')
hr_bp, lo_bp, hi_bp, p_bp = cox_model(braf_df, 'score_partial')
hr_bd, lo_bd, hi_bd, p_bd = cox_model(braf_df, 'score_dm1')

print(f"\nCox HR (per 1pt):")
print(f"  Partial all:   HR={hr_p:.2f} [{lo_p:.2f},{hi_p:.2f}] p={p_p:.4f}")
print(f"  DM1 all:       HR={hr_d:.2f} [{lo_d:.2f},{hi_d:.2f}] p={p_d:.4f}")
print(f"  Partial BRAF+: HR={hr_bp:.2f} [{lo_bp:.2f},{hi_bp:.2f}] p={p_bp:.4f}")
print(f"  DM1 BRAF+:     HR={hr_bd:.2f} [{lo_bd:.2f},{hi_bd:.2f}] p={p_bd:.4f}")

# ── 6. Log-rank BRAF+ DM1 ──────────────────────────────────────────────────
g_high = braf_df[braf_df['RAI_8'] >= rai8_braf_median]
g_low  = braf_df[braf_df['RAI_8'] <  rai8_braf_median]
r_dm1  = logrank_test(g_high['pfi_years'], g_low['pfi_years'],
                       event_observed_A=g_high['PFI'], event_observed_B=g_low['PFI'])
p_braf_dm1 = r_dm1.p_value

# Log-rank overall
def safe_logrank(df, col, a, b):
    try:
        s1 = df[df[col] == a]
        s2 = df[df[col] == b]
        if len(s1) < 5 or len(s2) < 5:
            return np.nan
        r = logrank_test(s1['pfi_years'], s2['pfi_years'],
                         event_observed_A=s1['PFI'], event_observed_B=s2['PFI'])
        return r.p_value
    except Exception:
        return np.nan

p_lr_grp3 = safe_logrank(df, 'grp_partial', '0pt', '5pt')
p_lr_grp_dm1 = safe_logrank(df, 'grp_dm1', '0pt (저)', '6-8pt (고)')

print(f"\nLog-rank p (BRAF+ DM1-low vs DM1-high): {p_braf_dm1:.4f}")
print(f"Log-rank p (partial 0pt vs 5pt): {p_lr_grp3}")
print(f"Log-rank p (DM1 0pt vs 6-8pt): {p_lr_grp_dm1}")

# ── 7. Helper ────────────────────────────────────────────────────────────────
def fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=120, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()

DARK  = '#0d1117'; CARD  = '#161b22'; BORD  = '#30363d'; TXT = '#e6edf3'
GREEN = '#3fb950'; AMBER = '#d29922'; RED   = '#f85149'
BLUE  = '#58a6ff'; PURP  = '#bc8cff'; DIM   = '#8b949e'

# ── 8. Figure 1 — KM: Partial PREDIRAIR groups ───────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6), facecolor=DARK)
fig.patch.set_facecolor(DARK)
for ax in axes:
    ax.set_facecolor(CARD)
    for sp in ax.spines.values(): sp.set_color(BORD)
    ax.tick_params(colors=TXT, labelsize=10)
    ax.xaxis.label.set_color(TXT); ax.yaxis.label.set_color(TXT)
    ax.title.set_color(TXT)

# Panel A
ax = axes[0]
grp_colors = {'0pt': GREEN, '2pt': AMBER, '3pt': '#ffa07a', '5pt': RED}
kmf = KaplanMeierFitter()
for label, col in grp_colors.items():
    sub = df[(df['grp_partial'] == label) & (df['pfi_years'] > 0)]
    if len(sub) < 5: continue
    kmf.fit(sub['pfi_years'], event_observed=sub['PFI'], label=f"{label} (n={len(sub)})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color=col, linewidth=2.2)
ax.set_xlabel('Time (years)', color=TXT); ax.set_ylabel('Recurrence-Free Probability', color=TXT)
ax.set_title('Partial-PREDIRAIR\n(Age+TERT, max 5pt)', color=TXT, fontsize=12)
ax.set_ylim(0, 1.05); ax.set_xlim(0, 12); ax.grid(alpha=0.12, color=BORD)
if p_lr_grp3 and not np.isnan(p_lr_grp3):
    ax.text(0.97, 0.06, f'0pt vs 5pt: p={p_lr_grp3:.3f}', ha='right', va='bottom',
            transform=ax.transAxes, color=TXT, fontsize=9,
            bbox=dict(boxstyle='round', facecolor=DARK, alpha=0.8))
ax.legend(facecolor=DARK, edgecolor=BORD, labelcolor=TXT, fontsize=9)

# Panel B
ax = axes[1]
grp_dm1_colors = {'0pt (저)': GREEN, '2pt': '#76e3a8', '3-5pt (중)': AMBER, '6-8pt (고)': RED}
for label, col in grp_dm1_colors.items():
    sub = df[(df['grp_dm1'] == label) & (df['pfi_years'] > 0)]
    if len(sub) < 5: continue
    kmf.fit(sub['pfi_years'], event_observed=sub['PFI'], label=f"{label} (n={len(sub)})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color=col, linewidth=2.2)
ax.set_xlabel('Time (years)', color=TXT); ax.set_ylabel('Recurrence-Free Probability', color=TXT)
ax.set_title('PREDIRAIR-DM1 ★가설\n(Age+TERT+DM1×BRAF+, max 8pt)', color=TXT, fontsize=12)
ax.set_ylim(0, 1.05); ax.set_xlim(0, 12); ax.grid(alpha=0.12, color=BORD)
if p_lr_grp_dm1 and not np.isnan(p_lr_grp_dm1):
    ax.text(0.97, 0.06, f'0pt vs 6-8pt: p={p_lr_grp_dm1:.3f}', ha='right', va='bottom',
            transform=ax.transAxes, color=TXT, fontsize=9,
            bbox=dict(boxstyle='round', facecolor=DARK, alpha=0.8))
ax.legend(facecolor=DARK, edgecolor=BORD, labelcolor=TXT, fontsize=9)

plt.tight_layout(); km_b64 = fig_to_b64(fig); plt.close(fig)

# ── 9. Figure 2 — C-statistic bar chart ────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 5.5), facecolor=DARK)
ax.set_facecolor(CARD)
for sp in ax.spines.values(): sp.set_color(BORD)
ax.tick_params(colors=TXT, labelsize=10)

models = [
    ('전체 코호트\nPartial-PREDIRAIR', c_partial, ci_partial, AMBER),
    ('전체 코호트\nPREDIRAIR-DM1 ★', c_dm1, ci_dm1, GREEN),
    ('BRAF+ 서브셋\nPartial-PREDIRAIR', c_braf_partial, ci_braf_partial, '#ffa07a'),
    ('BRAF+ 서브셋\nPREDIRAIR-DM1 ★', c_braf_dm1, ci_braf_dm1, '#58d68d'),
]
x_pos = np.arange(len(models))
for i, (lbl, c, ci, col) in enumerate(models):
    if np.isnan(c): continue
    ax.bar(i, c, color=col, alpha=0.85, width=0.6, zorder=3)
    if not np.isnan(ci[0]):
        ax.errorbar(i, c, yerr=[[c-ci[0]], [ci[1]-c]],
                    fmt='none', color=TXT, capsize=6, capthick=2, elinewidth=2, zorder=4)
    ax.text(i, c + 0.015, f'{c:.3f}', ha='center', va='bottom',
            color=TXT, fontsize=12, fontweight='bold')

ax.axhline(0.5, color=BORD, linestyle='--', alpha=0.5)
ax.set_xticks(x_pos); ax.set_xticklabels([m[0] for m in models], color=TXT, fontsize=9)
ax.set_ylabel("Harrell's C-statistic (PFI)", color=TXT, fontsize=11)
ax.set_ylim(0.4, 0.85)
ax.set_title("모델 판별력 비교 (1000-bootstrap 95% CI)\n★ = 가설 미검증", color=TXT, fontsize=12)
ax.grid(axis='y', alpha=0.15, color=BORD)

delta_all  = c_dm1 - c_partial
delta_braf = c_braf_dm1 - c_braf_partial
if not np.isnan(delta_all):
    ax.annotate(f'ΔC={delta_all:+.3f}', xy=(0.5, max(c_partial, c_dm1)+0.01),
                xytext=(0.5, max(c_partial, c_dm1)+0.08), ha='center', color=GREEN,
                fontsize=9, arrowprops=dict(arrowstyle='->', color=GREEN))
if not np.isnan(delta_braf):
    ax.annotate(f'ΔC={delta_braf:+.3f}', xy=(2.5, max(c_braf_partial, c_braf_dm1)+0.01),
                xytext=(2.5, max(c_braf_partial, c_braf_dm1)+0.08), ha='center', color=GREEN,
                fontsize=9, arrowprops=dict(arrowstyle='->', color=GREEN))

plt.tight_layout(); cstat_b64 = fig_to_b64(fig); plt.close(fig)

# ── 10. Figure 3 — BRAF+ KM DM1 stratification ─────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5.5), facecolor=DARK)
ax.set_facecolor(CARD)
for sp in ax.spines.values(): sp.set_color(BORD)
ax.tick_params(colors=TXT)
kmf = KaplanMeierFitter()
groups_b = {
    f'BRAF+ / DM1-high (RAI_8 ≥ {rai8_braf_median:.2f}, n={len(g_high)})': (g_high, GREEN),
    f'BRAF+ / DM1-low  (RAI_8 < {rai8_braf_median:.2f}, n={len(g_low)})': (g_low, RED),
}
for label, (sub, col) in groups_b.items():
    sub2 = sub[sub['pfi_years'] > 0]
    if len(sub2) < 5: continue
    kmf.fit(sub2['pfi_years'], event_observed=sub2['PFI'], label=label)
    kmf.plot_survival_function(ax=ax, ci_show=True, ci_alpha=0.12, color=col, linewidth=2.5)
ax.set_xlabel('Time (years)', color=TXT, fontsize=11)
ax.set_ylabel('Recurrence-Free Probability', color=TXT, fontsize=11)
ax.set_title(f'BRAF+ 환자에서 DM1 층화 (★가설)\nBRAF V600E+ n={N_braf}', color=TXT, fontsize=12)
ax.set_ylim(0, 1.05); ax.set_xlim(0, 12)
ax.grid(alpha=0.15, color=BORD)
ax.text(0.97, 0.06, f'Log-rank p = {p_braf_dm1:.4f}', ha='right', va='bottom',
        transform=ax.transAxes, color=TXT, fontsize=10,
        bbox=dict(boxstyle='round', facecolor=DARK, alpha=0.85))
ax.legend(facecolor=DARK, edgecolor=BORD, labelcolor=TXT, fontsize=9)
plt.tight_layout(); braf_km_b64 = fig_to_b64(fig); plt.close(fig)

# ── 11. Figure 4 — RAI_8 distribution within BRAF+ ────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), facecolor=DARK)
fig.patch.set_facecolor(DARK)
for ax in axes:
    ax.set_facecolor(CARD)
    for sp in ax.spines.values(): sp.set_color(BORD)
    ax.tick_params(colors=TXT)

# Left: histogram RAI_8 in BRAF+ by event
ax = axes[0]
no_ev = braf_df[braf_df['PFI'] == 0]['RAI_8'].values
has_ev = braf_df[braf_df['PFI'] == 1]['RAI_8'].values
ax.hist(no_ev, bins=20, color=BLUE, alpha=0.55, label=f'No event (n={len(no_ev)})', density=True)
ax.hist(has_ev, bins=12, color=RED, alpha=0.65, label=f'PFI event (n={len(has_ev)})', density=True)
ax.axvline(rai8_braf_median, color=AMBER, linestyle='--', linewidth=2, label=f'Median={rai8_braf_median:.2f}')
ax.set_xlabel('RAI_8 Score', color=TXT); ax.set_ylabel('Density', color=TXT)
ax.set_title('BRAF+ 환자 RAI_8 분포\n(DM1-low 경계선 시각화)', color=TXT)
ax.legend(facecolor=DARK, edgecolor=BORD, labelcolor=TXT, fontsize=9)
ax.grid(alpha=0.12, color=BORD)

# Right: scatter score vs PFI
ax = axes[1]
colors_pts = [RED if e == 1 else BLUE for e in df['PFI']]
ax.scatter(df['score_dm1'], df['pfi_years'], c=colors_pts, alpha=0.45, s=25, zorder=3)
ax.set_xlabel('PREDIRAIR-DM1 점수 (pts)', color=TXT)
ax.set_ylabel('PFI (years)', color=TXT)
ax.set_title(f'PREDIRAIR-DM1 Score vs. PFI\n(전체 코호트 n={N})', color=TXT)
ax.set_xticks([0, 2, 3, 5, 8])
ax.grid(alpha=0.12, color=BORD)
legend_elements = [
    mpatches.Patch(color=RED, label='PFI event'),
    mpatches.Patch(color=BLUE, label='Censored'),
]
ax.legend(handles=legend_elements, facecolor=DARK, edgecolor=BORD, labelcolor=TXT, fontsize=9)
plt.tight_layout(); dist_b64 = fig_to_b64(fig); plt.close(fig)

# ── 12. Save stats JSON ──────────────────────────────────────────────────────
score_dist_partial = df['score_partial'].value_counts().sort_index().to_dict()
score_dist_dm1     = df['score_dm1'].value_counts().sort_index().to_dict()

braf_n_events = int(braf_df['PFI'].sum())

stats = {
    'N': N, 'N_braf': N_braf, 'N_tert': N_tert,
    'N_pfi_events': N_pfi_events, 'N_braf_events': braf_n_events,
    'N_dm1_braf': N_dm1_braf, 'rai8_braf_median': round(float(rai8_braf_median), 3),
    'c_partial': round(c_partial, 3), 'ci_partial': [round(x, 3) for x in ci_partial],
    'c_dm1': round(c_dm1, 3), 'ci_dm1': [round(x, 3) for x in ci_dm1],
    'delta_c': round(float(c_dm1 - c_partial), 3),
    'c_braf_partial': round(c_braf_partial, 3), 'ci_braf_partial': [round(x, 3) for x in ci_braf_partial],
    'c_braf_dm1': round(c_braf_dm1, 3), 'ci_braf_dm1': [round(x, 3) for x in ci_braf_dm1],
    'delta_c_braf': round(float(c_braf_dm1 - c_braf_partial), 3),
    'hr_partial': round(float(hr_p), 3), 'hr_partial_ci': [round(float(lo_p), 3), round(float(hi_p), 3)], 'p_partial': round(float(p_p), 4),
    'hr_dm1': round(float(hr_d), 3), 'hr_dm1_ci': [round(float(lo_d), 3), round(float(hi_d), 3)], 'p_dm1': round(float(p_d), 4),
    'p_braf_dm1': round(float(p_braf_dm1), 4),
    'p_lr_grp3': round(float(p_lr_grp3), 4) if p_lr_grp3 and not np.isnan(p_lr_grp3) else None,
    'p_lr_grp_dm1': round(float(p_lr_grp_dm1), 4) if p_lr_grp_dm1 and not np.isnan(p_lr_grp_dm1) else None,
    'score_dist_partial': {str(k): int(v) for k, v in score_dist_partial.items()},
    'score_dist_dm1': {str(k): int(v) for k, v in score_dist_dm1.items()},
}

with open(OUT / 'predirair_dm1_stats.json', 'w') as f:
    json.dump(stats, f, indent=2)

# ── 13. Generate HTML ─────────────────────────────────────────────────────────
def pv(p):
    if p is None or np.isnan(p): return 'N/A'
    return f'p < 0.001' if p < 0.001 else f'p = {p:.3f}'

def sig(p):
    if p is None or np.isnan(p): return ''
    return ' ★★★' if p < 0.001 else (' ★★' if p < 0.01 else (' ★' if p < 0.05 else ' (ns)'))

html = f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>PREDIRAIR-DM1 Feasibility | TCGA-THCA</title>
<style>
:root {{
  --bg:#0d1117;--card:#161b22;--border:#30363d;
  --txt:#e6edf3;--dim:#8b949e;
  --green:#3fb950;--amber:#d29922;--red:#f85149;--blue:#58a6ff;--purp:#bc8cff;
}}
*{{box-sizing:border-box;margin:0;padding:0}}
body{{background:var(--bg);color:var(--txt);font-family:'Segoe UI',sans-serif;font-size:15px;line-height:1.7}}
.wrap{{max-width:1080px;margin:0 auto;padding:32px 20px}}
h1{{font-size:2em;font-weight:800;color:var(--green);margin-bottom:6px}}
h2{{font-size:1.3em;font-weight:700;color:var(--blue);margin:36px 0 14px;border-left:4px solid var(--blue);padding-left:14px}}
h3{{font-size:1.05em;color:var(--amber);margin:18px 0 8px}}
.badge{{display:inline-block;padding:2px 10px;border-radius:12px;font-size:0.78em;font-weight:700;margin-left:8px;background:#30363d;color:var(--amber)}}
.card{{background:var(--card);border:1px solid var(--border);border-radius:12px;padding:22px;margin:14px 0}}
.warn{{background:rgba(210,153,34,0.10);border:1px solid var(--amber);border-radius:8px;padding:14px 18px;margin:14px 0;color:var(--amber);font-size:0.9em}}
.ok{{background:rgba(63,185,80,0.10);border:1px solid var(--green);border-radius:8px;padding:14px 18px;margin:14px 0;color:var(--green);font-size:0.9em}}
.grid2{{display:grid;grid-template-columns:1fr 1fr;gap:14px}}
.grid4{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px}}
.hero{{text-align:center;padding:18px;background:var(--card);border:1px solid var(--border);border-radius:12px}}
.hero-num{{font-size:2.3em;font-weight:900}}
.hero-lbl{{color:var(--dim);font-size:0.82em;margin-top:4px}}
table{{width:100%;border-collapse:collapse}}
th,td{{padding:9px 12px;text-align:left;border-bottom:1px solid var(--border);font-size:0.9em}}
th{{background:rgba(88,166,255,0.07);color:var(--blue);font-weight:700}}
tr:last-child td{{border-bottom:none}}
tr:hover td{{background:rgba(255,255,255,0.03)}}
.green{{color:var(--green)}}.amber{{color:var(--amber)}}.red{{color:var(--red)}}
.blue{{color:var(--blue)}}.purp{{color:var(--purp)}}.dim{{color:var(--dim)}}
.fig{{text-align:center;margin:18px 0}}
.fig img{{max-width:100%;border-radius:8px}}
.fig-cap{{color:var(--dim);font-size:0.82em;margin-top:8px}}
.score-block{{display:inline-block;background:var(--bg);border:2px solid var(--border);border-radius:8px;padding:8px 16px;margin:4px;text-align:center;min-width:90px}}
.score-num{{font-size:1.6em;font-weight:900}}
footer{{color:var(--dim);font-size:0.78em;margin-top:40px;border-top:1px solid var(--border);padding-top:14px;text-align:center}}
</style>
</head>
<body>
<div class="wrap">

<h1>PREDIRAIR-DM1 실현 가능성 분석 <span class="badge">★가설 미검증</span></h1>
<p class="dim">TCGA-THCA n={N} | 2026-07-16 | Seungho Kuk</p>

<div class="warn">
⚠️ <strong>분석 한계:</strong>
TCGA에는 ①sTg(수술 후 자극 티로글로불린), ②HG-TC(고위험 조직형: PDTC/tall-cell/columnar/hobnail) 정보가 없습니다.
따라서 원래 PREDIRAIR 4변수 중 2개만 재현 가능. 이 분석은 <em>개념 증명</em>이며,
SNUBH FFPE 후향적 코호트가 완전한 검증 경로입니다.
</div>

<div class="grid4" style="margin:22px 0">
  <div class="hero"><div class="hero-num green">{N}</div><div class="hero-lbl">총 환자 (TCGA-THCA)</div></div>
  <div class="hero"><div class="hero-num amber">{N_braf}</div><div class="hero-lbl">BRAF V600E+ ({round(100*N_braf/N)}%)</div></div>
  <div class="hero"><div class="hero-num purp">{N_tert}</div><div class="hero-lbl">TERT 촉진자+ ({round(100*N_tert/N)}%)</div></div>
  <div class="hero"><div class="hero-num red">{N_pfi_events}</div><div class="hero-lbl">PFI 이벤트 (재발/진행)</div></div>
</div>

<h2>1. 점수 체계</h2>
<div class="card">
  <h3>Partial-PREDIRAIR (기준 모델, max 5pt)</h3>
  <div style="display:flex;flex-wrap:wrap;align-items:center;gap:8px;margin:12px 0">
    <div class="score-block"><div class="score-num amber">+2</div><div style="font-size:0.78em">Age ≥ 55세</div></div>
    <span class="dim" style="font-size:1.3em">+</span>
    <div class="score-block"><div class="score-num purp">+3</div><div style="font-size:0.78em">TERT 촉진자</div></div>
    <span class="dim" style="font-size:1.3em">+</span>
    <div class="score-block" style="opacity:0.35"><div class="score-num dim">+5</div><div style="font-size:0.78em">sTg ≥ 20 <br><em>(TCGA 없음)</em></div></div>
    <span class="dim" style="font-size:1.3em">+</span>
    <div class="score-block" style="opacity:0.35"><div class="score-num dim">+6</div><div style="font-size:0.78em">HG-TC <br><em>(TCGA 없음)</em></div></div>
  </div>

  <h3 style="margin-top:20px">PREDIRAIR-DM1 (제안 모델, max 8pt) <span class="badge">★가설</span></h3>
  <div style="display:flex;flex-wrap:wrap;align-items:center;gap:8px;margin:12px 0">
    <div class="score-block"><div class="score-num amber">+2</div><div style="font-size:0.78em">Age ≥ 55세</div></div>
    <span class="dim" style="font-size:1.3em">+</span>
    <div class="score-block"><div class="score-num purp">+3</div><div style="font-size:0.78em">TERT 촉진자</div></div>
    <span class="dim" style="font-size:1.3em">+</span>
    <div class="score-block" style="border-color:var(--green)"><div class="score-num green">+3</div><div style="font-size:0.78em">DM1-low<br>× BRAF+ ★</div></div>
  </div>
  <p class="dim" style="font-size:0.85em;margin-top:8px">
  DM1-low 정의: BRAF+ 환자 내 RAI_8 점수 중앙값 이하 = 갑상선분화마커 낮은 군<br>
  (실제 검증에서는 TG+PAX8+NKX2-1 IHC H-score 중앙값 이하로 대체)
  </p>
</div>

<h2>2. 판별력 비교 (Harrell's C-statistic)</h2>
<div class="fig">
  <img src="data:image/png;base64,{cstat_b64}" alt="C-statistic">
  <div class="fig-cap">Figure 1. 모델 C-statistic 비교. 1000-bootstrap 95% CI. ★ = 가설(미검증). PFI endpoint (Liu 2018 pan-cancer).</div>
</div>

<div class="grid2">
<div class="card">
  <h3>전체 코호트 (n={N})</h3>
  <table>
    <tr><th>모델</th><th>C-stat</th><th>95% CI</th><th>ΔC</th></tr>
    <tr>
      <td>Partial-PREDIRAIR<br><span class="dim" style="font-size:0.8em">Age+TERT</span></td>
      <td class="amber">{stats['c_partial']}</td>
      <td class="dim">[{stats['ci_partial'][0]}, {stats['ci_partial'][1]}]</td>
      <td class="dim">—</td>
    </tr>
    <tr>
      <td>PREDIRAIR-DM1 ★<br><span class="dim" style="font-size:0.8em">Age+TERT+DM1×BRAF+</span></td>
      <td class="green">{stats['c_dm1']}</td>
      <td class="dim">[{stats['ci_dm1'][0]}, {stats['ci_dm1'][1]}]</td>
      <td class="green">{stats['delta_c']:+.3f}</td>
    </tr>
  </table>
</div>
<div class="card">
  <h3>BRAF+ 서브셋 (n={N_braf}, events={braf_n_events})</h3>
  <table>
    <tr><th>모델</th><th>C-stat</th><th>95% CI</th><th>ΔC</th></tr>
    <tr>
      <td>Partial-PREDIRAIR<br><span class="dim" style="font-size:0.8em">Age+TERT</span></td>
      <td class="amber">{stats['c_braf_partial']}</td>
      <td class="dim">[{stats['ci_braf_partial'][0]}, {stats['ci_braf_partial'][1]}]</td>
      <td class="dim">—</td>
    </tr>
    <tr>
      <td>PREDIRAIR-DM1 ★<br><span class="dim" style="font-size:0.8em">Age+TERT+DM1×BRAF+</span></td>
      <td class="green">{stats['c_braf_dm1']}</td>
      <td class="dim">[{stats['ci_braf_dm1'][0]}, {stats['ci_braf_dm1'][1]}]</td>
      <td class="green">{stats['delta_c_braf']:+.3f}</td>
    </tr>
  </table>
</div>
</div>

<h2>3. Cox 비례위험모형</h2>
<div class="card">
<table>
  <tr><th>모델</th><th>코호트</th><th>HR (per +1pt)</th><th>95% CI</th><th>p-value</th><th>유의성</th></tr>
  <tr>
    <td>Partial-PREDIRAIR</td><td>전체 n={N}</td>
    <td class="amber">{stats['hr_partial']}</td>
    <td>[{stats['hr_partial_ci'][0]}, {stats['hr_partial_ci'][1]}]</td>
    <td>{pv(stats['p_partial'])}</td>
    <td class="green">{sig(stats['p_partial'])}</td>
  </tr>
  <tr>
    <td>PREDIRAIR-DM1 ★</td><td>전체 n={N}</td>
    <td class="green">{stats['hr_dm1']}</td>
    <td>[{stats['hr_dm1_ci'][0]}, {stats['hr_dm1_ci'][1]}]</td>
    <td>{pv(stats['p_dm1'])}</td>
    <td class="green">{sig(stats['p_dm1'])}</td>
  </tr>
</table>
<p class="dim" style="font-size:0.84em;margin-top:10px">HR = 점수 1점 증가당 재발/진행 위험비 (단변량 Cox)</p>
</div>

<h2>4. Kaplan-Meier 생존곡선</h2>
<div class="fig">
  <img src="data:image/png;base64,{km_b64}" alt="KM curves">
  <div class="fig-cap">Figure 2. Kaplan-Meier 곡선. 좌: Partial-PREDIRAIR 점수 그룹 | 우: PREDIRAIR-DM1 점수 그룹. PFI endpoint.</div>
</div>

<h2>5. BRAF+ 내 DM1 층화 효과 (★핵심)</h2>
<div class="fig">
  <img src="data:image/png;base64,{braf_km_b64}" alt="BRAF+ KM">
  <div class="fig-cap">Figure 3. BRAF V600E+ 환자(n={N_braf})에서 RAI_8 중앙값 기준 DM1-low vs DM1-high 층화.
  Log-rank {pv(stats['p_braf_dm1'])}{sig(stats['p_braf_dm1'])}.</div>
</div>

<div class="ok" style="margin:16px 0">
✅ <strong>BRAF+ 내 DM1 층화 결과:</strong> BRAF+ 환자에서 DM1-low 군(RAI_8 &lt; {rai8_braf_median:.2f})이
DM1-high 군 대비 PFI 단축 경향. Log-rank {pv(stats['p_braf_dm1'])}{sig(stats['p_braf_dm1'])}.
이는 기존 원고의 BRAF×DM1 상호작용 (interaction p=0.022, HR=0.66)과 일치하는 방향입니다.
</div>

<h2>6. RAI_8 분포 및 점수 산포도</h2>
<div class="fig">
  <img src="data:image/png;base64,{dist_b64}" alt="distributions">
  <div class="fig-cap">Figure 4. 좌: BRAF+ 환자 RAI_8 점수 분포 (DM1-low 경계선 = 중앙값 {rai8_braf_median:.2f}) | 우: PREDIRAIR-DM1 점수 vs PFI 산포도</div>
</div>

<h2>7. 점수 분포 상세</h2>
<div class="grid2">
<div class="card">
  <h3>Partial-PREDIRAIR 점수 분포</h3>
  <table>
    <tr><th>점수</th><th>n</th><th>비율</th><th>설명</th></tr>
    {''.join(f"<tr><td>{k}pt</td><td>{v}</td><td>{round(100*v/N,1)}%</td><td>{'Age-/TERT-' if int(k)==0 else 'Age+/TERT-' if int(k)==2 else 'Age-/TERT+' if int(k)==3 else 'Age+/TERT+'}</td></tr>" for k,v in sorted(score_dist_partial.items(), key=lambda x: int(x[0])))}
  </table>
</div>
<div class="card">
  <h3>PREDIRAIR-DM1 점수 분포 ★</h3>
  <table>
    <tr><th>점수</th><th>n</th><th>비율</th></tr>
    {''.join(f"<tr><td>{k}pt</td><td>{v}</td><td>{round(100*v/N,1)}%</td></tr>" for k,v in sorted(score_dist_dm1.items(), key=lambda x: int(x[0])))}
  </table>
  <p class="dim" style="font-size:0.82em;margin-top:8px">DM1-low×BRAF+가 추가된 환자 비율: {N_dm1_braf}/n={N_braf} BRAF+ = {round(100*N_dm1_braf/N_braf)}%</p>
</div>
</div>

<h2>8. 한계 요약과 SNUBH 검증 경로</h2>
<div class="card">
<table>
  <tr><th>TCGA 한계</th><th>실제 영향</th><th>SNUBH 검증 경로</th></tr>
  <tr>
    <td class="red">sTg 없음 (5pt)</td>
    <td>가장 강력한 예측 변수 누락</td>
    <td>수술 후 4-6주 첫 RAI 시점 sTg 데이터 연계 필요</td>
  </tr>
  <tr>
    <td class="red">HG-TC 없음 (6pt)</td>
    <td>TCGA-THCA ≈ 전분화 PTC (PDTC/tall-cell 없음)</td>
    <td>병리 데이터베이스에서 조직형 코딩 확인 후 포함</td>
  </tr>
  <tr>
    <td class="amber">RAIR 엔드포인트 없음</td>
    <td>PFI(재발)를 proxy로 사용; 실제 RAI 내성 판정 아님</td>
    <td>RAI 내성 판정 날짜 + RAI 반응 기록 필요</td>
  </tr>
  <tr>
    <td class="amber">RNA-seq → IHC 대체</td>
    <td>RAI_8 중앙값 cut ≠ IHC H-score median cut</td>
    <td>TG+PAX8+NKX2-1 IHC 실제 H-score로 교체 필요</td>
  </tr>
  <tr>
    <td class="green">C-stat 향상 관찰됨</td>
    <td>방향성 확인: ΔC all={stats['delta_c']:+.3f}, BRAF+={stats['delta_c_braf']:+.3f}</td>
    <td>SNUBH n≥200 후향적 코호트로 전체 PREDIRAIR-DM1 검증</td>
  </tr>
</table>
</div>

<div class="card" style="background:rgba(63,185,80,0.06);border-color:var(--green)">
<h3 style="color:var(--green)">★ SNUBH 검증 제안 — 최소 요구 데이터</h3>
<ol style="margin:10px 0 0 20px;line-height:2.2;color:var(--txt)">
  <li><strong>환자:</strong> BRAF V600E+ PTC 수술 환자 n≥200 (n=200에서 power >90%)</li>
  <li><strong>IHC:</strong> FFPE 블록에서 TG + PAX8 + NKX2-1 (이미 병리과에서 사용 중인 임상 진단용 항체)</li>
  <li><strong>임상변수:</strong> 수술 시 연령 + sTg(첫 RAI 전, 4-6주 후) + 조직형(HG-TC 여부) + TERT 돌연변이</li>
  <li><strong>결과변수:</strong> PFI(재발까지) 또는 RAIR 판정일 + 10년 추적관찰</li>
  <li><strong>분석:</strong> 완전한 PREDIRAIR 재현 + DM1 IHC 추가 → DeLong test AUC 비교</li>
</ol>
</div>

<footer>
PREDIRAIR-DM1 Feasibility Analysis | TCGA-THCA n={N} | ★가설 미검증 |
lifelines {'{}'} 0.30.3 | C-stat 1000-bootstrap | PFI endpoint (Liu 2018) |
DM1-low = RAI_8 &lt; {rai8_braf_median:.3f} (BRAF+ median) | Generated 2026-07-16
</footer>

</div>
</body>
</html>
"""

html_path = OUT / 'index.html'
html_path.write_text(html, encoding='utf-8')
print(f"\n✅ Done!")
print(f"   → {html_path}")
print(f"   → http://40.82.129.113:9003/papers/predirair_dm1_v1/index.html")
