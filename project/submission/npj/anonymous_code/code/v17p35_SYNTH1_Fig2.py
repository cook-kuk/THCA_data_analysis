"""Fig2 — DM1/DM2 Biology + Clinical (4-panel: marker heatmap, age violin, TDS/RAI box, histology Fisher)."""
from v17p35_SYNTH1_common import RES, save_figure, apply_npj, DM1_COLOR, DM2_COLOR, log
import pandas as pd, numpy as np, plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats

def build():
    log('Fig2 start')
    cohort = pd.read_csv(RES / 'v17p3' / 'tables' / 'A2_dm_score_full_cohort.tsv', sep='\t', low_memory=False)
    markers = pd.read_csv(RES / 'v17' / 'tables' / 'dark_matter_cluster_markers.tsv', sep='\t')
    hist = pd.read_csv(RES / 'v17' / 'tables' / 'quad_group_histology_distribution.tsv', sep='\t')

    # --- A: Marker heatmap (top 8 per cluster) ---
    top_per = markers.sort_values('fdr').groupby('cluster').head(8)
    genes = top_per['gene'].drop_duplicates().tolist()
    mat = (markers[markers['gene'].isin(genes)]
           .pivot_table(index='gene', columns='cluster', values='log2fc', aggfunc='mean')
           .reindex(genes))
    mat.columns = [f'DM{int(c)+1}' for c in mat.columns]

    # --- B: Age violin DM1 vs DM2 ---
    age = cohort.dropna(subset=['age', 'dm_like'])
    a_dm1 = age.loc[age['dm_like'] == 'DM1_like', 'age']
    a_dm2 = age.loc[age['dm_like'] == 'DM2_like', 'age']
    age_t, age_p = stats.ttest_ind(a_dm1, a_dm2, equal_var=False, nan_policy='omit')
    age_diff = float(a_dm2.mean() - a_dm1.mean())

    # --- C: TDS + RAI box ---
    tds = cohort.dropna(subset=['tds_score', 'dm_like'])
    rai_col = 'rai_score_recalc'
    rai = cohort.dropna(subset=[rai_col, 'dm_like'])

    # --- D: Histology Fisher (PTC subtypes) ---
    hist2 = hist.set_index('quad_group')
    bg = hist2.loc[['A_braf_only', 'C_double_pos']].sum(axis=0) if 'C_double_pos' in hist2.index else hist2.iloc[0]
    rg = hist2.loc[['B_ras_only', 'D_neg']].sum(axis=0) if 'D_neg' in hist2.index else hist2.iloc[1]
    hist_df = pd.DataFrame({'DM1-proxy': bg, 'DM2-proxy': rg})

    # ---- Compose ----
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'heatmap'}, {'type': 'xy'}],
               [{'type': 'xy'}, {'type': 'xy'}]],
        subplot_titles=(
            '<b>a</b>  Cluster-defining genes (top 8 per cluster, log2FC)',
            f'<b>b</b>  Age distribution (Δ={age_diff:+.1f}y, p={age_p:.2e})',
            '<b>c</b>  Thyroid Differentiation + RAI uptake by DM',
            '<b>d</b>  PTC histology subtype enrichment',
        ),
        horizontal_spacing=0.13, vertical_spacing=0.16,
    )

    fig.add_trace(go.Heatmap(z=mat.values, x=mat.columns, y=mat.index,
                             colorscale='RdBu_r', zmid=0,
                             colorbar=dict(title='log2FC', x=0.46, y=0.78, len=0.42, thickness=12)),
                  row=1, col=1)

    fig.add_trace(go.Violin(y=a_dm1, name='DM1', side='negative', line_color=DM1_COLOR, fillcolor=DM1_COLOR,
                            opacity=0.7, box_visible=True, meanline_visible=True), row=1, col=2)
    fig.add_trace(go.Violin(y=a_dm2, name='DM2', side='positive', line_color=DM2_COLOR, fillcolor=DM2_COLOR,
                            opacity=0.7, box_visible=True, meanline_visible=True), row=1, col=2)

    fig.add_trace(go.Box(y=tds.loc[tds['dm_like'] == 'DM1_like', 'tds_score'], name='DM1·TDS',
                         marker_color=DM1_COLOR, boxmean=True), row=2, col=1)
    fig.add_trace(go.Box(y=tds.loc[tds['dm_like'] == 'DM2_like', 'tds_score'], name='DM2·TDS',
                         marker_color=DM2_COLOR, boxmean=True), row=2, col=1)
    fig.add_trace(go.Box(y=rai.loc[rai['dm_like'] == 'DM1_like', rai_col], name='DM1·RAI',
                         marker_color=DM1_COLOR, opacity=0.55, boxmean=True), row=2, col=1)
    fig.add_trace(go.Box(y=rai.loc[rai['dm_like'] == 'DM2_like', rai_col], name='DM2·RAI',
                         marker_color=DM2_COLOR, opacity=0.55, boxmean=True), row=2, col=1)

    for col_lbl, color in [('DM1-proxy', DM1_COLOR), ('DM2-proxy', DM2_COLOR)]:
        fig.add_trace(go.Bar(name=col_lbl, x=hist_df.index, y=hist_df[col_lbl],
                             marker_color=color, opacity=0.85), row=2, col=2)

    fig.update_yaxes(title_text='Age (years)', row=1, col=2)
    fig.update_yaxes(title_text='Score', row=2, col=1)
    fig.update_yaxes(title_text='n samples', row=2, col=2)
    fig.update_xaxes(title_text='Histology subtype', row=2, col=2)
    fig.update_layout(barmode='group')

    apply_npj(fig)
    fig.update_layout(showlegend=True, legend=dict(orientation='h', x=0.0, y=-0.05))
    sizes = save_figure(fig, 'Fig2', width=1300, height=950)
    log(f'Fig2 done {sizes}')
    return sizes

if __name__ == '__main__':
    build()
