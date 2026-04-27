"""Fig6 — already composed in submission/npj/figures/Fig6.html. Just export PNG/PDF.
Strategy: parse Fig6.html for embedded plotly JSON if present, otherwise re-render via screenshot.
Simplest robust path: load via plotly.io.read_json if Fig6 was saved as JSON-able. We try that first.
Fallback: just leave existing Fig6.html and create PNG/PDF from a tagged Plotly figure embedded in
results/v17p35/figs/AMP4_*.html (Fig6 components)."""
from v17p35_SYNTH1_common import FIG_DIR, RES, save_figure, apply_npj, DM1_COLOR, DM2_COLOR, log
import json, re, os, plotly.graph_objects as go, plotly.io as pio
import pandas as pd, numpy as np
from plotly.subplots import make_subplots

FIG6_HTML = FIG_DIR / 'Fig6.html'

def _parse_existing():
    """Try to find an embedded plotly figure in the existing Fig6.html."""
    if not FIG6_HTML.exists():
        return None
    txt = FIG6_HTML.read_text(errors='ignore')
    m = re.search(r'Plotly\.newPlot\(\s*"[^"]+"\s*,\s*(\[\{.*?\}\])\s*,\s*(\{.*?\})\s*,', txt, re.DOTALL)
    if not m:
        return None
    try:
        data = json.loads(m.group(1))
        layout = json.loads(m.group(2))
        fig = go.Figure(data=data, layout=layout)
        return fig
    except Exception as e:
        log(f'Fig6 parse failed: {e}')
        return None

def _rebuild():
    """Rebuild Fig6 as 4-panel AMP4 decision tool: feature importance, ROC, calibration, decision rule."""
    log('Fig6 rebuild from AMP4 sources')
    coef = pd.read_csv(RES / 'v17p35' / 'tables' / 'AMP4_8gene_model_coefficients.tsv', sep='\t')
    cv = pd.read_csv(RES / 'v17p35' / 'tables' / 'AMP4_cv_performance.tsv', sep='\t')
    roc = pd.read_csv(RES / 'v17p35' / 'tables' / 'AMP4_roc_data.tsv', sep='\t')
    summary_p = RES / 'v17p35' / 'tables' / 'AMP4_summary.json'
    summary = json.loads(summary_p.read_text()) if summary_p.exists() else {}

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            '<b>a</b>  AMP4 8-gene model — feature importance',
            '<b>b</b>  Cross-validated ROC',
            '<b>c</b>  CV fold AUC distribution',
            '<b>d</b>  Coefficient sign + magnitude',
        ),
        horizontal_spacing=0.13, vertical_spacing=0.16,
    )

    coef_sorted = coef.assign(absc=coef.get('logreg_coef', coef.iloc[:, 1]).abs()).sort_values('absc', ascending=True)
    fig.add_trace(go.Bar(x=coef_sorted['rf_importance'] if 'rf_importance' in coef_sorted else coef_sorted.iloc[:, 2],
                         y=coef_sorted['gene'], orientation='h',
                         marker_color='#9467BD', showlegend=False), row=1, col=1)

    if 'model' in roc.columns:
        for m in roc['model'].unique():
            sub = roc[roc['model'] == m]
            fig.add_trace(go.Scatter(x=sub['fpr'], y=sub['tpr'], mode='lines', name=m), row=1, col=2)
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode='lines', line=dict(color='#999', dash='dash'),
                             showlegend=False), row=1, col=2)

    auc_col = next((c for c in cv.columns if 'auc' in c.lower()), None)
    if auc_col:
        fig.add_trace(go.Box(y=cv[auc_col], name='CV AUC', marker_color='#2CA02C', boxmean=True,
                             showlegend=False), row=2, col=1)

    coef_col = 'logreg_coef' if 'logreg_coef' in coef.columns else coef.columns[1]
    fig.add_trace(go.Bar(x=coef['gene'], y=coef[coef_col],
                         marker_color=['#D62728' if v < 0 else '#1F77B4' for v in coef[coef_col]],
                         showlegend=False), row=2, col=2)

    fig.update_xaxes(title_text='RF importance', row=1, col=1)
    fig.update_xaxes(title_text='False positive rate', row=1, col=2)
    fig.update_yaxes(title_text='True positive rate', row=1, col=2)
    fig.update_yaxes(title_text='AUC (CV folds)', row=2, col=1)
    fig.update_xaxes(title_text='Gene', tickangle=-30, row=2, col=2)
    fig.update_yaxes(title_text='Logistic coefficient', row=2, col=2)

    apply_npj(fig)
    return fig

def build():
    log('Fig6 export start')
    fig = _parse_existing()
    if fig is None:
        fig = _rebuild()
    sizes = save_figure(fig, 'Fig6', width=1300, height=950)
    log(f'Fig6 done {sizes}')
    return sizes

if __name__ == '__main__':
    build()
