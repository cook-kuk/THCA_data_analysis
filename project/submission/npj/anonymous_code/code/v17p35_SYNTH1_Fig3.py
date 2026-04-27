"""Fig3 — Trajectory + Dedifferentiation (4-panel)."""
from v17p35_SYNTH1_common import RES, save_figure, apply_npj, DM1_COLOR, DM2_COLOR, HIST_PALETTE, log
import pandas as pd, numpy as np, plotly.graph_objects as go
from plotly.subplots import make_subplots

def build():
    log('Fig3 start')
    traj = pd.read_csv(RES / 'v17' / 'tables' / 'trajectory_pseudotime.tsv', sep='\t', low_memory=False)
    pdtc = pd.read_csv(RES / 'v17p3' / 'tables' / 'A6_pdtc_rai_validation.tsv', sep='\t', low_memory=False)
    rai = pd.read_csv(RES / 'v17p3' / 'tables' / 'A6_rai_uptake_score.tsv', sep='\t')
    amp4_summary = (RES / 'v17p35' / 'tables' / 'AMP4_summary.json')
    roc = pd.read_csv(RES / 'v17p35' / 'tables' / 'AMP4_roc_data.tsv', sep='\t')

    # --- A: trajectory pseudotime by histology (PTC vs PDTC vs ATC) ---
    traj['hist_simple'] = traj['histology_subtype'].map(
        lambda h: 'ATC' if str(h).strip().upper() in ('ATC', 'ANAPLASTIC THYROID TUMOR') else
                  ('PDTC' if 'PDTC' in str(h).upper() or 'POORLY' in str(h).upper() else
                   ('FVPTC' if 'FVPTC' in str(h).upper() or 'FOLLICULAR VARIANT' in str(h).upper() else
                    ('PTC' if 'PTC' in str(h).upper() or 'PAPILLARY' in str(h).upper() else 'other')))
    )

    # --- B: RAI by histology ---
    pdtc['hist_simple'] = pdtc['histology_subtype'].map(
        lambda h: 'ATC' if str(h).strip().upper() in ('ATC', 'ANAPLASTIC THYROID TUMOR') else
                  ('PDTC' if 'PDTC' in str(h).upper() or 'POORLY' in str(h).upper() else 'PTC')
    )
    rai_long = pdtc[['hist_simple', 'rai_score']].dropna()

    # --- C: DM1/DM2 along pseudotime (density-like ridge proxy = box per cluster) ---
    dm_traj = traj.dropna(subset=['v17_dark_cluster', 'pseudotime'])

    # --- D: AMP4 ROC ---
    # filter best model
    if 'model' in roc.columns:
        models = roc['model'].unique()
    else:
        models = []

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            '<b>a</b>  PTC → PDTC → ATC trajectory (pseudotime by histology)',
            '<b>b</b>  RAI score by histology (PDTC > ATC)',
            '<b>c</b>  DM1 along trajectory (ATC-proximal evidence)',
            '<b>d</b>  AMP4 8-gene model — ROC for high-risk vs low-risk',
        ),
        horizontal_spacing=0.13, vertical_spacing=0.16,
    )

    for h, color in [('PTC', '#1F77B4'), ('FVPTC', '#FF7F0E'), ('PDTC', '#9467BD'), ('ATC', '#D62728')]:
        s = traj[traj['hist_simple'] == h]
        if len(s) == 0: continue
        fig.add_trace(go.Box(y=s['pseudotime'], name=f'{h} (n={len(s)})', marker_color=color, boxmean=True),
                      row=1, col=1)

    for h, color in [('PTC', '#1F77B4'), ('PDTC', '#9467BD'), ('ATC', '#D62728')]:
        s = rai_long[rai_long['hist_simple'] == h]
        if len(s) == 0: continue
        fig.add_trace(go.Box(y=s['rai_score'], name=f'{h} RAI (n={len(s)})', marker_color=color, boxmean=True),
                      row=1, col=2)

    for cluster, color in [('DM1', DM1_COLOR), ('DM2', DM2_COLOR)]:
        s = dm_traj[dm_traj['v17_dark_cluster'] == cluster]
        if len(s) == 0: continue
        fig.add_trace(go.Histogram(x=s['pseudotime'], name=cluster, marker_color=color, opacity=0.6,
                                   nbinsx=30, histnorm='probability density'),
                      row=2, col=1)

    if len(models) > 0:
        for m in models[:3]:
            sub = roc[roc['model'] == m]
            fig.add_trace(go.Scatter(x=sub['fpr'], y=sub['tpr'], mode='lines', name=str(m), line=dict(width=2)),
                          row=2, col=2)
        fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode='lines', line=dict(color='#999', dash='dash'),
                                 showlegend=False), row=2, col=2)

    fig.update_yaxes(title_text='Pseudotime', row=1, col=1)
    fig.update_yaxes(title_text='RAI score', row=1, col=2)
    fig.update_xaxes(title_text='Pseudotime', row=2, col=1)
    fig.update_yaxes(title_text='Density', row=2, col=1)
    fig.update_xaxes(title_text='False positive rate', row=2, col=2)
    fig.update_yaxes(title_text='True positive rate', row=2, col=2)
    fig.update_layout(barmode='overlay')

    apply_npj(fig)
    fig.update_layout(showlegend=True, legend=dict(orientation='h', x=0.0, y=-0.05))
    sizes = save_figure(fig, 'Fig3', width=1300, height=950)
    log(f'Fig3 done {sizes}')
    return sizes

if __name__ == '__main__':
    build()
