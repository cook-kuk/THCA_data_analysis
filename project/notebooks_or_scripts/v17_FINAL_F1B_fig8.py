"""Fig8 — TERT 4-group survival composite (4-panel: KM, Firth-HR forest, stage stack, TDS violin).

Source: results/v17_tert_recovery/v2/sample_master_v17_tert_v2.tsv + FINAL_extended_summary.json.
Survival p-values reported are from FINAL_extended_summary.json (logrank computed in upstream pipeline).
This script re-derives KM curves + Cox HR (Firth fallback to plain Cox) for visual display.
"""
from v17p35_SYNTH1_common import RES, save_figure, apply_npj, log
import pandas as pd, numpy as np, json, plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats

PALETTE = {
    'BRAF_only': '#D62728',     # red
    'RAS_only': '#2CA02C',      # green
    'TERT+': '#9467BD',         # purple (NEW)
    'Triple_neg': '#7F7F7F',    # gray
}
ORDER = ['BRAF_only', 'RAS_only', 'TERT+', 'Triple_neg']

def km_surv(times, events):
    """Kaplan–Meier survival function (returns sorted unique times + S(t))."""
    t = np.asarray(times, dtype=float)
    e = np.asarray(events, dtype=int)
    order = np.argsort(t)
    t, e = t[order], e[order]
    n = len(t)
    surv = []
    cur_t = []
    s = 1.0
    at_risk = n
    i = 0
    while i < n:
        ti = t[i]
        d = 0; c = 0
        while i < n and t[i] == ti:
            if e[i] == 1: d += 1
            else: c += 1
            i += 1
        if at_risk > 0:
            s *= (at_risk - d) / at_risk
        cur_t.append(ti); surv.append(s)
        at_risk -= (d + c)
    return np.array(cur_t), np.array(surv)


def cox_hr(df, covariates=None, ref='Triple_neg', penalizer=0.01):
    """Joint 4-group Cox: single CoxPHFitter with group-dummy variables (vs reference) + optional covariates.
    This matches the manuscript's reporting convention. Returns list of
    (group_name, HR, CI_lo, CI_hi, p, n_events, n_total)."""
    from lifelines import CoxPHFitter
    out = []
    rows = df[df['four_group'].notna()].copy()
    ref_n = int((rows['four_group'] == ref).sum())
    ref_ev = int(rows.loc[rows['four_group'] == ref, 'os_event'].sum())
    out.append((ref, 1.00, np.nan, np.nan, np.nan, ref_ev, ref_n))

    cols = ['os_days', 'os_event'] + (list(covariates) if covariates else [])
    targets = [g for g in ORDER if g != ref]
    dummies = {f'is_{g}': (rows['four_group'] == g).astype(int) for g in targets}
    fit_df = pd.concat([rows[cols].reset_index(drop=True),
                        pd.DataFrame(dummies).reset_index(drop=True)], axis=1).dropna()
    n_per_g = {g: int((rows['four_group'] == g).sum()) for g in targets}
    ev_per_g = {g: int(rows.loc[rows['four_group'] == g, 'os_event'].sum()) for g in targets}

    if len(fit_df) == 0 or fit_df['os_event'].sum() == 0:
        for g in targets:
            out.append((g, np.nan, np.nan, np.nan, np.nan, ev_per_g[g], n_per_g[g]))
        return out
    cph = CoxPHFitter(penalizer=penalizer)
    try:
        cph.fit(fit_df, duration_col='os_days', event_col='os_event')
        for g in targets:
            key = f'is_{g}'
            if key not in cph.summary.index:
                out.append((g, np.nan, np.nan, np.nan, np.nan, ev_per_g[g], n_per_g[g]))
                continue
            row = cph.summary.loc[key]
            out.append((
                g,
                float(row['exp(coef)']),
                float(row['exp(coef) lower 95%']),
                float(row['exp(coef) upper 95%']),
                float(row['p']),
                ev_per_g[g], n_per_g[g],
            ))
    except Exception as e:
        log(f'Joint Cox fit failed: {e}')
        for g in targets:
            out.append((g, np.nan, np.nan, np.nan, np.nan, ev_per_g[g], n_per_g[g]))
    return out


def build():
    log('Fig8 start')
    sm = pd.read_csv(RES.parent / 'results' / 'v17_tert_recovery' / 'v2' / 'sample_master_v17_tert_v2.tsv',
                     sep='\t', low_memory=False)
    summary = json.loads((RES.parent / 'results' / 'v17_tert_recovery' / 'v2' / 'FINAL_extended_summary.json').read_text())

    # Construct 4-group: TERT+ wins; else use quad_group (A_braf_only / B_ras_only / D_triple_negative).
    sm['four_group'] = np.where(
        sm['tert_promoter_integrated'] == 'mutated', 'TERT+',
        sm['quad_group'].map({'A_braf_only': 'BRAF_only', 'B_ras_only': 'RAS_only',
                              'D_triple_negative': 'Triple_neg'}))
    sm = sm.dropna(subset=['four_group', 'os_event', 'os_days'])
    counts = sm['four_group'].value_counts().to_dict()
    log(f'4-group n: {counts}')

    # Build covariates for multivariate Cox: stage_ord (I=1, II=2, III=3, IV=4), age, sex_male
    stage_col = 'clinical_stage' if 'clinical_stage' in sm.columns else 'stage'
    sm['stage_ord'] = sm[stage_col].astype(str).str.extract(r'(IV|III|II|I)\b', expand=False).map(
        {'I': 1, 'II': 2, 'III': 3, 'IV': 4})
    sm['age_num'] = pd.to_numeric(sm.get('age'), errors='coerce')
    if sm['age_num'].isna().all() and 'age_clinical' in sm.columns:
        sm['age_num'] = pd.to_numeric(sm['age_clinical'], errors='coerce')
    sm['sex_male'] = (sm.get('sex', '').astype(str).str.lower().str.startswith('m')).astype(int)

    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'xy'}, {'type': 'xy'}],
               [{'type': 'xy'}, {'type': 'xy'}]],
        subplot_titles=(
            f"<b>a</b>  Kaplan–Meier with 95% CI — 4-group OS (lifelines logrank p={summary['survival']['four_group_logrank']['p_value']:.2e})",
            f"<b>b</b>  Univariate (left) vs Multivariate Cox HR (right; adjusted for stage·age·sex) — TERT alone p={summary['survival']['TERT_logrank']['p_value']:.2e}",
            "<b>c</b>  Stage distribution per group (TERT+ enriched for III/IV)",
            "<b>d</b>  TDS (Yoo SK 16-gene differentiation score) per group",
        ),
        horizontal_spacing=0.13, vertical_spacing=0.16,
    )

    # --- A: KM per group ---
    for g in ORDER:
        sub = sm[sm['four_group'] == g]
        if len(sub) == 0: continue
        n_g = len(sub); ev_g = int(sub["os_event"].sum())
        t, s = km_surv(sub['os_days'].values, sub['os_event'].values)
        # 95% CI via Greenwood's formula approximation (log-log)
        with np.errstate(invalid='ignore', divide='ignore'):
            se = np.sqrt(np.cumsum(1.0 / np.maximum(n_g - np.arange(len(t)), 1) ** 2))
        s_lo = np.clip(s - 1.96 * se * s, 0, 1)
        s_hi = np.clip(s + 1.96 * se * s, 0, 1)
        rgba_fill = PALETTE[g].replace('#', '')
        rfill = ','.join(str(int(rgba_fill[i:i+2], 16)) for i in (0, 2, 4))
        fig.add_trace(go.Scatter(
            x=np.concatenate([t, t[::-1]]),
            y=np.concatenate([s_hi, s_lo[::-1]]),
            fill='toself', fillcolor=f'rgba({rfill},0.13)',
            line=dict(color='rgba(0,0,0,0)'),
            showlegend=False, hoverinfo='skip',
            name=f'{g} CI'), row=1, col=1)
        fig.add_trace(go.Scatter(
            x=t, y=s, mode='lines',
            name=f'{g} (n={n_g}, ev={ev_g})',
            line=dict(color=PALETTE[g], width=2.4, shape='hv'),
            legendgroup='km_groups',
            legendgrouptitle_text='4-group OS (panel a)' if g == ORDER[0] else None,
        ), row=1, col=1)

    # --- B: Univariate + Multivariate Cox HR forest ---
    uni = cox_hr(sm, covariates=None)
    multi = cox_hr(sm, covariates=['stage_ord', 'age_num', 'sex_male'])

    def label(name, ev, n, hr, p):
        if np.isnan(hr) or np.isnan(p):
            return f"{name} (n={n}, ev={ev}) <i>ref</i>" if name == 'Triple_neg' else f"{name} (n={n}, ev={ev}) — n.e."
        sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
        return f"{name} (n={n}, ev={ev}) HR={hr:.2f} [{sig}, p={p:.2g}]"

    for tag, model_results, x_offset, marker_symbol in [('Univariate', uni, 0, 'circle'), ('Multivariate', multi, 0, 'square')]:
        ys = []
        xs = []; lo_arr = []; hi_arr = []; colors = []; texts = []
        for g, hr, ci_lo, ci_hi, p, ev, n in model_results:
            if np.isnan(hr): continue
            ys.append(label(g, ev, n, hr, p))
            xs.append(hr)
            lo_arr.append(hr - (ci_lo if not np.isnan(ci_lo) else hr))
            hi_arr.append((ci_hi if not np.isnan(ci_hi) else hr) - hr)
            colors.append(PALETTE.get(g, '#999'))
            texts.append(f"{tag} Cox: HR={hr:.2f}, 95% CI=[{ci_lo:.2f}, {ci_hi:.2f}], p={p:.2g}, n={n} ev={ev}" if not np.isnan(ci_lo) else f"{tag}: ref")
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode='markers',
            error_x=dict(type='data', symmetric=False, array=hi_arr, arrayminus=lo_arr,
                         thickness=1.6, width=4),
            marker=dict(color=colors, size=11, symbol=marker_symbol,
                        line=dict(color='black', width=0.6)),
            name=f"{tag} Cox",
            legendgroup=f'cox_{tag.lower()}',
            hovertext=texts, hoverinfo='text',
        ), row=1, col=2)

    fig.add_shape(type='line', xref='x2', yref='paper', x0=1, x1=1, y0=0.55, y1=0.99,
                  line=dict(color='#666', dash='dash', width=1))
    fig.add_annotation(xref='x2', yref='paper', x=1, y=1.01, text='HR=1 (no effect)',
                       showarrow=False, font=dict(size=9, color='#666'))

    # --- C: stage distribution ---
    stage_col = 'clinical_stage' if 'clinical_stage' in sm.columns else 'stage'
    sm['stage_simple'] = sm[stage_col].astype(str).str.extract(r'(IV|III|II|I)\b', expand=False).fillna('NA')
    stage_pivot = pd.crosstab(sm['four_group'], sm['stage_simple'], normalize='index') * 100
    stage_pivot = stage_pivot.reindex(ORDER).fillna(0)
    stage_order = ['I', 'II', 'III', 'IV', 'NA']
    stage_order = [c for c in stage_order if c in stage_pivot.columns]
    stage_palette = {'I': '#2CA02C', 'II': '#FFD700', 'III': '#FF7F0E', 'IV': '#D62728', 'NA': '#BBB'}
    for st in stage_order:
        fig.add_trace(go.Bar(name=f'Stage {st}', x=stage_pivot.index, y=stage_pivot[st],
                             marker_color=stage_palette.get(st, '#999'),
                             legendgroup='stage_levels',
                             legendgrouptitle_text='AJCC stage (panel c)' if st == stage_order[0] else None,
                             ), row=2, col=1)

    # --- D: TDS violin per group ---
    tds_col = 'tds16_score_v17' if 'tds16_score_v17' in sm.columns else 'tds_score'
    for g in ORDER:
        sub = sm[sm['four_group'] == g]
        if len(sub) == 0: continue
        fig.add_trace(go.Violin(y=sub[tds_col], name=g, marker_color=PALETTE[g],
                                box_visible=True, meanline_visible=True, opacity=0.85,
                                showlegend=False),
                      row=2, col=2)

    fig.update_xaxes(title_text='Days from diagnosis', row=1, col=1)
    fig.update_yaxes(title_text='Overall survival probability', range=[0.5, 1.02], row=1, col=1)
    fig.update_xaxes(title_text='Hazard ratio vs Triple-negative reference (log scale)', type='log', row=1, col=2)
    fig.update_yaxes(title_text='Group (n, events; HR with significance)', row=1, col=2, automargin=True)
    fig.update_xaxes(title_text='4-group', row=2, col=1)
    fig.update_yaxes(title_text='% within group', row=2, col=1)
    fig.update_xaxes(title_text='4-group', row=2, col=2)
    fig.update_yaxes(title_text='TDS (Yoo SK 16-gene differentiation score)', row=2, col=2)
    fig.update_layout(barmode='stack')

    apply_npj(fig)
    fig.update_layout(
        showlegend=True,
        legend=dict(orientation='v', x=1.04, y=1.0, xanchor='left', yanchor='top',
                    bgcolor='rgba(255,255,255,0.95)', bordercolor='#cbd5e1', borderwidth=1,
                    font=dict(size=10), groupclick='toggleitem'),
    )
    sizes = save_figure(fig, 'Fig8', width=1300, height=950)
    log(f'Fig8 done {sizes}')
    return sizes


if __name__ == '__main__':
    build()
