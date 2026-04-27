"""Fig4 — BRAF/RAS Orthogonality (4-panel: violin by driver, contingency, 17 outliers, correlation)."""
from v17p35_SYNTH1_common import RES, save_figure, apply_npj, DM1_COLOR, DM2_COLOR, DRIVER_PALETTE, log
import pandas as pd, numpy as np, plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats

def build():
    log('Fig4 start')
    df = pd.read_csv(RES / 'v17p3' / 'tables' / 'A2_dm_score_full_cohort.tsv', sep='\t', low_memory=False)

    # --- A: prob_dm2 by driver ---
    df['driver_simple'] = df['driver_anchor'].map(lambda x: x if x in ('BRAF', 'RAS') else ('Driver-neg' if pd.isna(x) or x == 'unknown' else 'Other'))

    # --- B: BRAF/RAS × DM contingency ---
    cont = pd.crosstab(df['driver_simple'], df['dm_like']).reindex(['BRAF', 'RAS', 'Other', 'Driver-neg']).fillna(0)
    cont = cont[['DM1_like', 'DM2_like']] if 'DM1_like' in cont.columns else cont

    # --- C: 17 outliers ---
    out_braf = df[(df['driver_anchor'] == 'BRAF') & (df['dm_like'] == 'DM2_like')].copy()
    out_ras = df[(df['driver_anchor'] == 'RAS') & (df['dm_like'] == 'DM1_like')].copy()
    outliers = pd.concat([out_braf, out_ras], ignore_index=True)

    # --- D: DM-driver correlation (Spearman) ---
    df['is_braf'] = (df['driver_anchor'] == 'BRAF').astype(int)
    df['is_ras'] = (df['driver_anchor'] == 'RAS').astype(int)
    rho_braf, p_braf = stats.spearmanr(df['prob_dm1'], df['is_braf'])
    rho_ras, p_ras = stats.spearmanr(df['prob_dm2'], df['is_ras'])

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            '<b>a</b>  P(DM2) by driver anchor',
            '<b>b</b>  Driver × DM contingency (n / row %)',
            f'<b>c</b>  17 outliers ({len(out_braf)} BRAF/DM2 + {len(out_ras)} RAS/DM1)',
            f'<b>d</b>  Spearman: BRAF↔P(DM1) ρ={rho_braf:.2f}; RAS↔P(DM2) ρ={rho_ras:.2f}',
        ),
        horizontal_spacing=0.13, vertical_spacing=0.16,
    )

    for grp, color in [('BRAF', '#D62728'), ('RAS', '#2CA02C'), ('Other', '#9467BD'), ('Driver-neg', '#7F7F7F')]:
        s = df[df['driver_simple'] == grp]
        if len(s) == 0: continue
        fig.add_trace(go.Violin(y=s['prob_dm2'], name=f'{grp} (n={len(s)})', marker_color=color,
                                box_visible=True, meanline_visible=True, opacity=0.8), row=1, col=1)

    fig.add_trace(go.Heatmap(z=cont.values, x=['DM1', 'DM2'], y=cont.index,
                             colorscale='Blues', text=cont.values, texttemplate='%{text}',
                             colorbar=dict(title='n', x=1.0, y=0.78, len=0.42, thickness=12)),
                  row=1, col=2)

    out_x = list(out_braf['prob_dm1']) + list(out_ras['prob_dm1'])
    out_y = list(out_braf['prob_dm2']) + list(out_ras['prob_dm2'])
    out_color = ['#D62728'] * len(out_braf) + ['#2CA02C'] * len(out_ras)
    out_text = ([f'BRAF/DM2 — {sid}' for sid in out_braf['sample_id']] +
                [f'RAS/DM1 — {sid}' for sid in out_ras['sample_id']])
    fig.add_trace(go.Scatter(x=out_x, y=out_y, mode='markers',
                             marker=dict(color=out_color, size=11, line=dict(width=1, color='black'), opacity=0.85),
                             text=out_text, hoverinfo='text', showlegend=False),
                  row=2, col=1)
    fig.add_shape(type='line', xref='x3', yref='y3', x0=0.5, x1=0.5, y0=0, y1=1,
                  line=dict(color='#888', dash='dash'))
    fig.add_shape(type='line', xref='x3', yref='y3', x0=0, x1=1, y0=0.5, y1=0.5,
                  line=dict(color='#888', dash='dash'))

    sub = df.dropna(subset=['prob_dm1', 'prob_dm2'])
    fig.add_trace(go.Scatter(x=sub.loc[sub['driver_simple'] == 'BRAF', 'prob_dm1'],
                             y=sub.loc[sub['driver_simple'] == 'BRAF', 'prob_dm2'],
                             mode='markers', name='BRAF', marker=dict(color='#D62728', size=4, opacity=0.5)),
                  row=2, col=2)
    fig.add_trace(go.Scatter(x=sub.loc[sub['driver_simple'] == 'RAS', 'prob_dm1'],
                             y=sub.loc[sub['driver_simple'] == 'RAS', 'prob_dm2'],
                             mode='markers', name='RAS', marker=dict(color='#2CA02C', size=4, opacity=0.5)),
                  row=2, col=2)
    fig.add_trace(go.Scatter(x=sub.loc[~sub['driver_simple'].isin(['BRAF', 'RAS']), 'prob_dm1'],
                             y=sub.loc[~sub['driver_simple'].isin(['BRAF', 'RAS']), 'prob_dm2'],
                             mode='markers', name='Other/neg', marker=dict(color='#7F7F7F', size=4, opacity=0.4)),
                  row=2, col=2)

    fig.update_yaxes(title_text='P(DM2)', range=[-0.05, 1.05], row=1, col=1)
    fig.update_xaxes(title_text='P(DM1)', row=2, col=1)
    fig.update_yaxes(title_text='P(DM2)', row=2, col=1)
    fig.update_xaxes(title_text='P(DM1)', row=2, col=2)
    fig.update_yaxes(title_text='P(DM2)', row=2, col=2)

    apply_npj(fig)
    fig.update_layout(showlegend=True, legend=dict(orientation='h', x=0.0, y=-0.05))
    sizes = save_figure(fig, 'Fig4', width=1300, height=950)
    log(f'Fig4 done {sizes}')

    # outliers table for supp
    outliers_path = RES / 'v17p35' / 'tables' / 'AMP_17_outliers.tsv'
    outliers[['sample_id', 'dataset', 'driver_anchor', 'dm_like', 'prob_dm1', 'prob_dm2',
              'tds_score', 'age', 'histology_subtype']].to_csv(outliers_path, sep='\t', index=False)
    log(f'Fig4 outliers table → {outliers_path}')
    return sizes

if __name__ == '__main__':
    build()
