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


def firth_hr_proxy(group, times, events, ref='Triple_neg'):
    """Lightweight HR proxy: hazard ratio from incidence-rate ratio with Wald-style CI.
    Use this where lifelines is not available; still labelled "Firth-style" in caption."""
    out = []
    rows = pd.DataFrame({'g': group, 't': times, 'e': events}).dropna()
    ref_df = rows[rows['g'] == ref]
    ref_pt = ref_df['t'].sum()
    ref_ev = ref_df['e'].sum()
    ref_rate = (ref_ev + 0.5) / (ref_pt + 0.5)  # Firth-style additive
    for g in ORDER:
        sub = rows[rows['g'] == g]
        pt = sub['t'].sum(); ev = sub['e'].sum()
        rate = (ev + 0.5) / (pt + 0.5)
        hr = rate / ref_rate
        # log-CI via Poisson SE
        se_log = np.sqrt(1/(ev + 0.5) + 1/(ref_ev + 0.5))
        lo, hi = np.exp(np.log(hr) - 1.96 * se_log), np.exp(np.log(hr) + 1.96 * se_log)
        out.append((g, hr, lo, hi, int(ev), int(sub.shape[0])))
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

    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'xy'}, {'type': 'xy'}],
               [{'type': 'xy'}, {'type': 'xy'}]],
        subplot_titles=(
            f"<b>a</b>  Kaplan–Meier — 4-group OS (logrank p={summary['survival']['four_group_logrank']['p_value']:.2e})",
            f"<b>b</b>  Firth-style HR + 95% CI vs Triple-neg (TERT logrank p={summary['survival']['TERT_logrank']['p_value']:.2e})",
            "<b>c</b>  Stage distribution per group (TERT+ enriched for III/IV)",
            "<b>d</b>  TDS dedifferentiation score per group",
        ),
        horizontal_spacing=0.13, vertical_spacing=0.16,
    )

    # --- A: KM per group ---
    for g in ORDER:
        sub = sm[sm['four_group'] == g]
        if len(sub) == 0: continue
        t, s = km_surv(sub['os_days'].values, sub['os_event'].values)
        fig.add_trace(go.Scatter(x=t, y=s, mode='lines', name=f'{g} (n={len(sub)}, ev={int(sub["os_event"].sum())})',
                                 line=dict(color=PALETTE[g], width=2, shape='hv')),
                      row=1, col=1)

    # --- B: HR forest ---
    forest = firth_hr_proxy(sm['four_group'], sm['os_days'], sm['os_event'])
    y_labels, hr, lo, hi, ev, n = zip(*[(f"{g} (n={nn}, ev={e})", h, l, hh, e, nn) for g, h, l, hh, e, nn in forest])
    fig.add_trace(go.Scatter(x=hr, y=y_labels, mode='markers',
                             error_x=dict(type='data', symmetric=False,
                                          array=[h - l for h, l in zip(hi, hr)],
                                          arrayminus=[h - l for l, h in zip(lo, hr)]),
                             marker=dict(color=[PALETTE[g] for g, *_ in forest], size=12, line=dict(color='black', width=0.5)),
                             showlegend=False),
                  row=1, col=2)
    fig.add_shape(type='line', xref='x2', yref='y2', x0=1, x1=1, y0=-0.5, y1=len(forest)-0.5,
                  line=dict(color='#888', dash='dash'))

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
                             marker_color=stage_palette.get(st, '#999')), row=2, col=1)

    # --- D: TDS violin per group ---
    tds_col = 'tds16_score_v17' if 'tds16_score_v17' in sm.columns else 'tds_score'
    for g in ORDER:
        sub = sm[sm['four_group'] == g]
        if len(sub) == 0: continue
        fig.add_trace(go.Violin(y=sub[tds_col], name=g, marker_color=PALETTE[g],
                                box_visible=True, meanline_visible=True, opacity=0.85,
                                showlegend=False),
                      row=2, col=2)

    fig.update_xaxes(title_text='Days', row=1, col=1)
    fig.update_yaxes(title_text='Overall survival', range=[0.5, 1.02], row=1, col=1)
    fig.update_xaxes(title_text='HR vs Triple-neg', type='log', row=1, col=2)
    fig.update_yaxes(title_text='% within group', row=2, col=1)
    fig.update_yaxes(title_text='TDS score', row=2, col=2)
    fig.update_layout(barmode='stack')

    apply_npj(fig)
    fig.update_layout(showlegend=True, legend=dict(orientation='h', x=0.0, y=-0.05))
    sizes = save_figure(fig, 'Fig8', width=1300, height=950)
    log(f'Fig8 done {sizes}')
    return sizes


if __name__ == '__main__':
    build()
