"""Fig7 — Drug Actionability (4-panel: drug volcano, MOA enrichment, MEK/HMGCR class, BRAF cell-line concordance)."""
from v17p35_SYNTH1_common import RES, save_figure, apply_npj, DM1_COLOR, DM2_COLOR, log
import pandas as pd, numpy as np, plotly.graph_objects as go
from plotly.subplots import make_subplots

def build():
    log('Fig7 start')
    dm1 = pd.read_csv(RES / 'v17p35' / 'tables' / 'FIX1_top_drugs_dm1_selective_v2.tsv', sep='\t')
    dm2 = pd.read_csv(RES / 'v17p35' / 'tables' / 'FIX1_top_drugs_dm2_selective_v2.tsv', sep='\t')
    moa = pd.read_csv(RES / 'v17p35' / 'tables' / 'FIX1_moa_enrichment_v2.tsv', sep='\t')
    cell = pd.read_csv(RES / 'v17p35' / 'tables' / 'FIX1_celline_dm_scores_v2.tsv', sep='\t')

    # --- A: drug volcano ---
    all_drugs = pd.concat([dm1, dm2], ignore_index=True).drop_duplicates(subset=['compound', 'dose_uM'])
    all_drugs = all_drugs.dropna(subset=['delta_lfc', 'pvalue'])
    all_drugs['neglog_p'] = -np.log10(all_drugs['pvalue'].clip(lower=1e-10))

    # --- B: MOA bar (DM1 vs DM2 enriched) ---
    if 'dm1_count' in moa.columns:
        moa['delta'] = moa['dm1_count'] - moa['dm2_count']
        moa_top = moa.assign(absd=moa['delta'].abs()).sort_values('absd', ascending=True).tail(15)

    # --- C: MEK / HMGCR class selectivity ---
    targets_dm1 = ['MEK inhibitor', 'HMG-CoA reductase inhibitor', 'HMGCR inhibitor']
    moa_lower = moa['moa'].fillna('').astype(str).str.lower()
    is_mek = moa_lower.str.contains('mek')
    is_stat = moa_lower.str.contains('hmg|statin|reductase')
    sel_classes = moa[is_mek | is_stat].copy()
    if 'delta' not in sel_classes.columns:
        sel_classes['delta'] = sel_classes['dm1_count'] - sel_classes['dm2_count']

    # --- D: BRAF cell line concordance (DM1 prediction) ---
    cell['braf_label'] = cell['mutation_label'].astype(str).str.upper().str.contains('BRAF')

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            '<b>a</b>  DM1-selective drug volcano (Δ LFC vs −log10 p)',
            '<b>b</b>  MOA enrichment (Δ count: DM1 − DM2)',
            '<b>c</b>  Class-level selectivity (MEK + HMGCR)',
            '<b>d</b>  Cell-line DM1/DM2 axis (BRAF=red)',
        ),
        horizontal_spacing=0.13, vertical_spacing=0.16,
    )

    fig.add_trace(go.Scatter(x=all_drugs['delta_lfc'], y=all_drugs['neglog_p'],
                             mode='markers',
                             marker=dict(color=['#FF7F0E' if v < 0 else '#1F77B4' for v in all_drugs['delta_lfc']],
                                         size=7, opacity=0.7, line=dict(color='black', width=0.3)),
                             text=[f'{c}<br>MOA: {m}<br>p={p:.2g}' for c, m, p in
                                   zip(all_drugs['compound'], all_drugs['moa'], all_drugs['pvalue'])],
                             hoverinfo='text', showlegend=False),
                  row=1, col=1)
    fig.add_shape(type='line', xref='x1', yref='y1', x0=0, x1=0, y0=0,
                  y1=all_drugs['neglog_p'].max(), line=dict(color='#888', dash='dash'))

    if 'delta' in moa.columns:
        fig.add_trace(go.Bar(x=moa_top['delta'], y=moa_top['moa'], orientation='h',
                             marker_color=['#FF7F0E' if v > 0 else '#1F77B4' for v in moa_top['delta']],
                             showlegend=False), row=1, col=2)

    fig.add_trace(go.Bar(x=sel_classes['moa'], y=sel_classes['delta'],
                         marker_color=['#FF7F0E' if v > 0 else '#1F77B4' for v in sel_classes['delta']],
                         showlegend=False),
                  row=2, col=1)

    fig.add_trace(go.Scatter(x=cell.loc[~cell['braf_label'], 'dm1_score_z'],
                             y=cell.loc[~cell['braf_label'], 'dm2_score_z'],
                             mode='markers', name='non-BRAF',
                             marker=dict(color='#7F7F7F', size=8, opacity=0.6,
                                         line=dict(color='black', width=0.5))),
                  row=2, col=2)
    fig.add_trace(go.Scatter(x=cell.loc[cell['braf_label'], 'dm1_score_z'],
                             y=cell.loc[cell['braf_label'], 'dm2_score_z'],
                             mode='markers+text', name='BRAF',
                             text=cell.loc[cell['braf_label'], 'sample'],
                             textposition='top center', textfont=dict(size=9),
                             marker=dict(color='#D62728', size=11, line=dict(color='black', width=0.8))),
                  row=2, col=2)

    fig.update_xaxes(title_text='Δ LFC (DM1 − DM2)', row=1, col=1)
    fig.update_yaxes(title_text='−log10 p', row=1, col=1)
    fig.update_xaxes(title_text='Δ MOA count', row=1, col=2)
    fig.update_xaxes(title_text='MOA class', row=2, col=1, tickangle=-25)
    fig.update_yaxes(title_text='Δ count', row=2, col=1)
    fig.update_xaxes(title_text='DM1 score (z)', row=2, col=2)
    fig.update_yaxes(title_text='DM2 score (z)', row=2, col=2)

    apply_npj(fig)
    fig.update_layout(showlegend=True, legend=dict(orientation='h', x=0.0, y=-0.05))
    sizes = save_figure(fig, 'Fig7', width=1300, height=950)
    log(f'Fig7 done {sizes}')
    return sizes

if __name__ == '__main__':
    build()
