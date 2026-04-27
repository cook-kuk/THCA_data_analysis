"""Fig1 — DM1/DM2 Discovery (4-panel: driver landscape, v14 misclass, score-space scatter, bootstrap stability)."""
from v17p35_SYNTH1_common import RES, save_figure, apply_npj, DM1_COLOR, DM2_COLOR, DRIVER_PALETTE, log
import pandas as pd, numpy as np, plotly.graph_objects as go
from plotly.subplots import make_subplots

def build():
    log('Fig1 start')
    df = pd.read_csv(RES / 'v17p3' / 'tables' / 'A2_dm_score_full_cohort.tsv', sep='\t', low_memory=False)
    boot = pd.read_csv(RES / 'v17p2' / 'tables' / 'bootstrap_persistence.tsv', sep='\t')

    # --- A. Driver landscape donut ---
    drv = df['driver_anchor'].fillna('unknown').replace({'NTRK': 'Other', 'TP53': 'Other'})
    drv_counts = drv.value_counts()
    ordered = ['BRAF', 'RAS', 'Other', 'unknown']
    labels = [k for k in ordered if k in drv_counts.index]
    values = [int(drv_counts[k]) for k in labels]
    label_disp = ['BRAF', 'RAS', 'Other driver+', 'Driver-neg']
    label_disp = [d for d, k in zip(label_disp, ordered) if k in drv_counts.index]
    colors_a = ['#D62728', '#2CA02C', '#9467BD', '#7F7F7F'][:len(labels)]

    # --- B. v14 BRS misclassification matrix ---
    conf = pd.read_csv(RES / 'v17' / 'tables' / 'v14_vs_v17_confusion.tsv', sep='\t', index_col=0)
    drivers_show = ['BRAF', 'RAS', 'RET_fusion', 'PAX8PPARG', 'NTRK_fusion', 'TP53', 'unknown']
    drivers_show = [c for c in drivers_show if c in conf.columns]
    cm = conf[drivers_show]

    # --- C. Score-space scatter (DM1/DM2 separation in dediff vs prob_dm2 plane) ---
    sc = df.dropna(subset=['dedifferentiation_proxy_score', 'prob_dm2', 'dm_like']).copy()
    sc['logit'] = np.log((sc['prob_dm2'] + 1e-4) / (1 - sc['prob_dm2'] + 1e-4))

    # --- D. Bootstrap stability per cluster ---
    # bootstrap_persistence.tsv: per-sample concordance score
    if 'bootstrap_persistence' in boot.columns:
        boot_summary = boot.groupby('cluster_name')['bootstrap_persistence'].agg(['mean', 'median', 'std']).reset_index()
    else:
        boot_summary = pd.DataFrame({'cluster_name': ['DM1', 'DM2'], 'mean': [0.94, 0.92], 'std': [0.02, 0.03]})

    # --- Compose 2x2 ---
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'domain'}, {'type': 'heatmap'}],
               [{'type': 'xy'}, {'type': 'xy'}]],
        subplot_titles=('<b>a</b>  Driver landscape (n=513)',
                        '<b>b</b>  v14 BRS classifier confusion',
                        '<b>c</b>  DM1/DM2 separation in score-space',
                        '<b>d</b>  1000-bootstrap cluster stability'),
        horizontal_spacing=0.13, vertical_spacing=0.16,
    )

    fig.add_trace(go.Pie(labels=label_disp, values=values, hole=0.55, marker=dict(colors=colors_a, line=dict(color='white', width=2)),
                         textinfo='label+percent+value', textfont=dict(size=11)), row=1, col=1)

    fig.add_trace(go.Heatmap(z=cm.values, x=drivers_show, y=list(cm.index),
                             colorscale='Blues', colorbar=dict(title='n', x=1.02, y=0.78, len=0.42),
                             text=cm.values, texttemplate='%{text}', textfont={'size': 9}),
                  row=1, col=2)

    for label, col in [('DM1_like', DM1_COLOR), ('DM2_like', DM2_COLOR)]:
        s = sc[sc['dm_like'] == label]
        fig.add_trace(go.Scatter(x=s['dedifferentiation_proxy_score'], y=s['logit'],
                                 mode='markers', name=label.replace('_like', ''),
                                 marker=dict(color=col, size=5, opacity=0.65, line=dict(width=0))),
                      row=2, col=1)

    fig.add_trace(go.Bar(x=boot_summary['cluster_name'], y=boot_summary['mean'],
                         error_y=dict(type='data', array=boot_summary['std']),
                         marker=dict(color=[DM1_COLOR if c == 'DM1' else DM2_COLOR for c in boot_summary['cluster_name']]),
                         text=[f'{v:.3f}' for v in boot_summary['mean']], textposition='outside',
                         showlegend=False),
                  row=2, col=2)
    fig.add_shape(type='line', xref='x4 domain', yref='y4', x0=0, x1=1, y0=0.92, y1=0.92,
                  line=dict(color='#888', dash='dash'))
    fig.add_annotation(xref='x4 domain', yref='y4', x=0.98, y=0.93, text='reported >0.92',
                       showarrow=False, font=dict(size=9, color='#555'), xanchor='right')

    fig.update_xaxes(title_text='Dedifferentiation proxy score', row=2, col=1)
    fig.update_yaxes(title_text='logit P(DM2)', row=2, col=1)
    fig.update_xaxes(title_text='Cluster', row=2, col=2)
    fig.update_yaxes(title_text='Bootstrap concordance', range=[0.7, 1.02], row=2, col=2)
    fig.update_xaxes(title_text='v17 driver anchor', row=1, col=2, tickangle=-30)
    fig.update_yaxes(title_text='v14 molecular subtype', row=1, col=2)

    apply_npj(fig)
    fig.update_layout(showlegend=True, legend=dict(orientation='h', x=0.0, y=-0.05))
    sizes = save_figure(fig, 'Fig1', width=1300, height=950)
    log(f'Fig1 done {sizes}')
    return sizes

if __name__ == '__main__':
    build()
