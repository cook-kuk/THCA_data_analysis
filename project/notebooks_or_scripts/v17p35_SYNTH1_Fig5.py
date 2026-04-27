"""Fig5 — Hot/Cold Immune (4-panel: GSEA NES bar, scRNA immune stack, immune evasion heatmap, Hot/Cold composite)."""
from v17p35_SYNTH1_common import RES, save_figure, apply_npj, DM1_COLOR, DM2_COLOR, log
import pandas as pd, numpy as np, plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats

INFLAM_KEYS = ['Inflammatory_Response', 'Interferon_Gamma_Response', 'TNFA_Signaling_Via_NFKB',
               'Allograft_Rejection', 'IL6_JAK_STAT3_Signaling', 'Complement']

def cohens_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    pooled = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) / (len(a)+len(b)-2))
    return (np.mean(a) - np.mean(b)) / pooled if pooled > 0 else 0.0

def build():
    log('Fig5 start')
    gsea = pd.read_csv(RES / 'v17p3' / 'tables' / 'F2_gsea_hallmark_proper.tsv', sep='\t')
    sc = pd.read_csv(RES / 'v17p35' / 'tables' / 'FIX5_immune_cell_type_breakdown.tsv', sep='\t')
    evas = pd.read_csv(RES / 'v17p3' / 'tables' / 'A5_immune_evasion_genes.tsv', sep='\t')
    hm = pd.read_csv(RES / 'v17p2' / 'tables' / 'hallmark_score_per_sample.tsv', sep='\t')

    # --- A: GSEA top inflammatory ---
    # File has 'Name' column = workflow tag ('prerank') and 'pathway.1' = actual pathway name.
    # Pick whichever column holds the real pathway names.
    candidates = [c for c in ['pathway.1', 'pathway_name', 'name', 'gene_set']
                  if c in gsea.columns and gsea[c].astype(str).str.contains(
                      'Inflammatory|Interferon|TNFA|Allograft', case=False, regex=True).any()]
    if candidates:
        name_col = candidates[0]
    else:
        # fallback: try every string column to find one containing pathway names
        for c in gsea.columns:
            if gsea[c].astype(str).str.contains('Inflammatory|Interferon', case=False, regex=True).any():
                name_col = c
                break
        else:
            name_col = 'Name'
    nes_col = 'NES'
    fdr_col = 'FDR q-val' if 'FDR q-val' in gsea.columns else 'fdr'
    inflam = gsea[gsea[name_col].astype(str).str.contains(
        'Inflammatory|Interferon|TNFA|Allograft|IL6|Complement|IL2.STAT5|Coagulation|Apical', case=False, regex=True)].copy()
    # NES sign convention in this file is DM2-vs-DM1 (negative = DM1-up). Manuscript reports DM1-up (inflammatory hot)
    # so we flip sign to plot inflammatory pathways as positive (DM1-enriched).
    inflam['NES_dm1_up'] = -inflam[nes_col]
    inflam = inflam.sort_values('NES_dm1_up', key=abs, ascending=False).head(8)
    inflam_nes_col = 'NES_dm1_up'

    # --- B: scRNA immune cell type ---
    sc_pivot = sc.pivot_table(index='dominant_dm', columns='cell_type', values='fraction_within_dm',
                              aggfunc='sum', fill_value=0)

    # --- C: Immune evasion heatmap ---
    evas_genes = ['CD274', 'PDCD1', 'CTLA4', 'IDO1', 'HLA-A', 'HLA-B', 'HLA-C', 'B2M']
    evas_genes = [g for g in evas_genes if g in evas.columns]
    evas_means = evas.groupby('v17_dark_cluster')[evas_genes].mean()

    # --- D: Hot/Cold composite (cytolytic + IFN-γ + immune fraction proxy) ---
    cyto = evas[['sample_id', 'v17_dark_cluster', 'cytolytic_score']].copy()
    if 'Hallmark_Interferon_Gamma' in hm.columns:
        ifg = hm[['sample_id', 'cluster', 'Hallmark_Interferon_Gamma', 'Hallmark_Inflammatory_Response']].copy()
        merged = cyto.merge(ifg, on='sample_id', how='inner', suffixes=('', '_y'))
        merged['cluster_use'] = merged['v17_dark_cluster'].fillna(merged['cluster'])
        for col in ['cytolytic_score', 'Hallmark_Interferon_Gamma', 'Hallmark_Inflammatory_Response']:
            merged[col + '_z'] = (merged[col] - merged[col].mean()) / merged[col].std()
        merged['hot_cold_composite'] = merged[[c + '_z' for c in ['cytolytic_score', 'Hallmark_Interferon_Gamma', 'Hallmark_Inflammatory_Response']]].mean(axis=1)
        out_path = RES / 'v17p35' / 'tables' / 'AMP3_hot_cold_composite.tsv'
        merged[['sample_id', 'cluster_use', 'cytolytic_score_z', 'Hallmark_Interferon_Gamma_z',
                'Hallmark_Inflammatory_Response_z', 'hot_cold_composite']].to_csv(out_path, sep='\t', index=False)
        a = merged.loc[merged['cluster_use'] == 'DM1', 'hot_cold_composite'].dropna()
        b = merged.loc[merged['cluster_use'] == 'DM2', 'hot_cold_composite'].dropna()
        d_val = cohens_d(a, b)
        _, hc_p = stats.mannwhitneyu(a, b, alternative='two-sided') if len(a) and len(b) else (None, np.nan)
    else:
        merged = pd.DataFrame()
        d_val = float('nan'); hc_p = float('nan'); a = b = pd.Series(dtype=float)

    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'xy'}, {'type': 'xy'}],
               [{'type': 'heatmap'}, {'type': 'xy'}]],
        subplot_titles=(
            '<b>a</b>  GSEA: top inflammatory pathways (NES)',
            '<b>b</b>  scRNA immune cell breakdown by dominant DM',
            '<b>c</b>  Immune evasion gene expression (cluster mean)',
            f"<b>d</b>  Hot/Cold composite — Cohen's d={d_val:+.2f}, p={hc_p:.2e}",
        ),
        horizontal_spacing=0.13, vertical_spacing=0.16,
    )

    fig.add_trace(go.Bar(x=inflam[inflam_nes_col], y=inflam[name_col], orientation='h',
                         marker_color=['#D62728' if v > 0 else '#1F77B4' for v in inflam[inflam_nes_col]],
                         text=[f'FDR={f:.2g}' for f in inflam[fdr_col]], textposition='outside',
                         showlegend=False, name='GSEA NES'),
                  row=1, col=1)

    sc_long = sc_pivot.reset_index().melt(id_vars='dominant_dm', var_name='cell_type', value_name='frac')
    cell_palette = ['#3b82f6', '#f59e0b', '#10b981', '#a855f7', '#ec4899', '#8C564B', '#0ea5e9', '#7F7F7F', '#84cc16']
    for i, ct in enumerate(sc_pivot.columns):
        fig.add_trace(go.Bar(name=str(ct), x=sc_pivot.index, y=sc_pivot[ct],
                             marker_color=cell_palette[i % len(cell_palette)],
                             legendgroup='cell_types', legendgrouptitle_text='scRNA cell types (panel b)'),
                      row=1, col=2)

    fig.add_trace(go.Heatmap(z=evas_means.values, x=evas_genes, y=evas_means.index,
                             colorscale='Viridis',
                             colorbar=dict(title='Mean expr', x=0.46, y=0.22, len=0.42, thickness=12)),
                  row=2, col=1)

    if len(a) and len(b):
        fig.add_trace(go.Violin(y=a, name='DM1 cluster', marker_color=DM1_COLOR, box_visible=True,
                                meanline_visible=True, opacity=0.85,
                                legendgroup='dm_clusters', legendgrouptitle_text='DM-axis cluster (panel d)'),
                      row=2, col=2)
        fig.add_trace(go.Violin(y=b, name='DM2 cluster', marker_color=DM2_COLOR, box_visible=True,
                                meanline_visible=True, opacity=0.85,
                                legendgroup='dm_clusters'),
                      row=2, col=2)

    fig.update_xaxes(title_text='NES (DM1-up direction; positive = inflammatory hot in DM1)', row=1, col=1)
    fig.update_yaxes(title_text='Fraction (within dominant-DM)', row=1, col=2)
    fig.update_xaxes(title_text='Dominant DM', row=1, col=2)
    fig.update_yaxes(title_text='Hot/Cold composite (z)', row=2, col=2)
    fig.update_xaxes(title_text='DM-axis cluster', row=2, col=2)
    fig.update_xaxes(title_text='Immune-evasion / HLA gene', row=2, col=1)
    fig.update_yaxes(title_text='DM cluster', row=2, col=1)
    fig.update_layout(barmode='stack')

    apply_npj(fig)
    fig.update_layout(
        showlegend=True,
        legend=dict(orientation='v', x=1.04, y=1.0, xanchor='left', yanchor='top',
                    bgcolor='rgba(255,255,255,0.95)', bordercolor='#cbd5e1', borderwidth=1,
                    font=dict(size=10), groupclick='toggleitem'),
    )
    sizes = save_figure(fig, 'Fig5', width=1300, height=950)
    log(f'Fig5 done {sizes} d={d_val:+.3f} p={hc_p:.2e}')
    return sizes, d_val, hc_p

if __name__ == '__main__':
    build()
