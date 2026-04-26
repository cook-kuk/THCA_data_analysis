// Interactive threshold explorer — binds to #threshold-explorer, loads model_curves_v2.json,
// renders a live-updating ROC + precision/recall panel and metrics as the user drags the
// threshold slider. Drop this script on any page that contains a #threshold-explorer host.
(function () {
  'use strict';
  const HOST_ID = 'threshold-explorer';
  const DATA_URL = '../assets/data/model_curves_v2.json';

  function fmt(n, digits) {
    if (n === null || n === undefined || Number.isNaN(n)) return '—';
    if (typeof n !== 'number') return String(n);
    return n.toFixed(digits ?? 3);
  }

  function makeSelect(id, label) {
    const wrap = document.createElement('label');
    wrap.className = 'te-label';
    wrap.setAttribute('for', id);
    wrap.innerHTML = `<span>${label}</span>`;
    const sel = document.createElement('select');
    sel.id = id;
    sel.className = 'te-select';
    wrap.appendChild(sel);
    return { wrap, sel };
  }

  function fillOptions(select, options, selected) {
    select.innerHTML = '';
    for (const opt of options) {
      const o = document.createElement('option');
      o.value = opt;
      o.textContent = opt;
      if (opt === selected) o.selected = true;
      select.appendChild(o);
    }
  }

  function build(host, data) {
    host.classList.add('te-root');
    host.innerHTML = `
      <div class="te-controls">
        <div class="te-control-row"></div>
        <div class="te-threshold-row">
          <label for="te-threshold" class="te-label te-slider-label"><span>Decision threshold</span></label>
          <input type="range" id="te-threshold" min="0" max="100" step="1" value="50" class="te-slider" aria-describedby="te-threshold-value">
          <span id="te-threshold-value" class="te-threshold-badge mono">0.50</span>
          <button type="button" id="te-use-youden" class="btn sm ghost">Youden 자동</button>
        </div>
      </div>
      <div class="te-grid">
        <div class="te-metrics" aria-live="polite">
          <div class="te-metric"><div class="te-mlabel">AUC (95% CI)</div><div class="te-mvalue mono" id="te-m-auc">—</div></div>
          <div class="te-metric"><div class="te-mlabel">PR-AUC</div><div class="te-mvalue mono" id="te-m-prauc">—</div></div>
          <div class="te-metric"><div class="te-mlabel">Brier</div><div class="te-mvalue mono" id="te-m-brier">—</div></div>
          <div class="te-metric"><div class="te-mlabel">Precision</div><div class="te-mvalue mono" id="te-m-prec">—</div></div>
          <div class="te-metric"><div class="te-mlabel">Recall (sens)</div><div class="te-mvalue mono" id="te-m-rec">—</div></div>
          <div class="te-metric"><div class="te-mlabel">F1</div><div class="te-mvalue mono" id="te-m-f1">—</div></div>
          <div class="te-metric"><div class="te-mlabel">Balanced acc</div><div class="te-mvalue mono" id="te-m-bacc">—</div></div>
          <div class="te-metric"><div class="te-mlabel">N (eval)</div><div class="te-mvalue mono" id="te-m-n">—</div></div>
        </div>
        <div id="te-chart" class="te-chart" role="img" aria-label="ROC with current threshold marker"></div>
      </div>
      <p class="subtle small">슬라이더를 움직이면 precision/recall/F1/balanced accuracy가 실시간으로 갱신됩니다. TCGA-THCA 내부는 CV 예측값, GSE27155/GSE126698은 TCGA 학습 모델의 외부 적용 점수입니다.</p>
    `;

    const row = host.querySelector('.te-control-row');
    const fs = makeSelect('te-fs', 'Feature set');
    const md = makeSelect('te-md', 'Model');
    const ds = makeSelect('te-ds', 'Eval cohort');
    row.appendChild(fs.wrap);
    row.appendChild(md.wrap);
    row.appendChild(ds.wrap);

    const state = { featureSet: 'TierA67_clean', model: null, dataset: 'TCGA-THCA_internal', tIdx: 50 };
    if (!data[state.featureSet]) state.featureSet = Object.keys(data)[0];

    function renderControls() {
      fillOptions(fs.sel, Object.keys(data), state.featureSet);
      const models = Object.keys(data[state.featureSet] || {});
      if (!models.includes(state.model)) state.model = models.includes('LogReg_l2') ? 'LogReg_l2' : models[0];
      fillOptions(md.sel, models, state.model);
      const evals = Object.keys(data[state.featureSet]?.[state.model] || {});
      if (!evals.includes(state.dataset)) state.dataset = evals[0];
      fillOptions(ds.sel, evals, state.dataset);
    }

    function currentEntry() {
      return data[state.featureSet]?.[state.model]?.[state.dataset] || null;
    }

    function refresh() {
      const e = currentEntry();
      const $ = (id) => host.querySelector('#' + id);
      if (!e) {
        ['te-m-auc','te-m-prauc','te-m-brier','te-m-prec','te-m-rec','te-m-f1','te-m-bacc','te-m-n'].forEach(id => $(id).textContent = '—');
        return;
      }
      const tc = e.threshold_curve || { thresholds: [], precision: [], recall: [], f1: [], balanced_accuracy: [] };
      const idx = Math.min(Math.max(0, state.tIdx), Math.max(0, tc.thresholds.length - 1));
      const t = tc.thresholds[idx];
      const slider = $('te-threshold');
      slider.max = Math.max(0, tc.thresholds.length - 1);
      slider.value = idx;
      $('te-threshold-value').textContent = t !== undefined ? t.toFixed(3) : '—';
      const ciLo = e.auc_ci_lo, ciHi = e.auc_ci_hi;
      const ciText = (ciLo != null && ciHi != null) ? ` (${fmt(ciLo)}–${fmt(ciHi)})` : '';
      $('te-m-auc').textContent = fmt(e.auc) + ciText;
      $('te-m-prauc').textContent = fmt(e.pr_auc);
      $('te-m-brier').textContent = e.brier != null ? fmt(e.brier, 4) : '—';
      $('te-m-prec').textContent = fmt(tc.precision?.[idx]);
      $('te-m-rec').textContent = fmt(tc.recall?.[idx]);
      $('te-m-f1').textContent = fmt(tc.f1?.[idx]);
      $('te-m-bacc').textContent = fmt(tc.balanced_accuracy?.[idx]);
      $('te-m-n').textContent = (e.n_samples != null) ? String(e.n_samples) : '—';
      renderChart(e, t);
    }

    function renderChart(entry, tCurrent) {
      const div = host.querySelector('#te-chart');
      if (!window.Plotly) { div.textContent = 'Plotly not loaded.'; return; }
      const fpr = entry.roc?.fpr || [];
      const tpr = entry.roc?.tpr || [];
      const thresholds = entry.threshold_curve?.thresholds || [];
      const scores = entry.scores || [];
      const y = entry.y_true || [];
      // Compute marker position on ROC at the current threshold
      let markerFpr = null, markerTpr = null;
      if (scores.length === y.length && scores.length) {
        let tp=0, fp=0, fn=0, tn=0;
        for (let i=0;i<scores.length;i++) {
          const p = scores[i] >= tCurrent ? 1 : 0;
          if (p === 1 && y[i] === 1) tp++;
          else if (p === 1 && y[i] === 0) fp++;
          else if (p === 0 && y[i] === 1) fn++;
          else tn++;
        }
        const fprV = (fp+tn) ? fp/(fp+tn) : 0;
        const tprV = (tp+fn) ? tp/(tp+fn) : 0;
        markerFpr = fprV; markerTpr = tprV;
      }
      const traces = [
        { x: fpr, y: tpr, mode: 'lines', name: 'ROC', line: { color: '#14b8a6', width: 2.5 }, hovertemplate: 'FPR=%{x:.3f}<br>TPR=%{y:.3f}<extra></extra>' },
        { x: [0,1], y: [0,1], mode: 'lines', name: 'Chance', line: { color: 'rgba(148,163,184,0.5)', dash: 'dash' }, hoverinfo: 'skip' },
      ];
      if (markerFpr !== null) {
        traces.push({ x: [markerFpr], y: [markerTpr], mode: 'markers', name: `threshold=${tCurrent.toFixed(3)}`, marker: { color: '#f59e0b', size: 12, line: { color: '#fff', width: 2 } }, hovertemplate: `FPR=${markerFpr.toFixed(3)}<br>TPR=${markerTpr.toFixed(3)}<extra></extra>` });
      }
      const layout = {
        margin: { l: 44, r: 12, t: 16, b: 36 },
        height: 280,
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        font: { family: 'Inter, ui-sans-serif, system-ui', color: '#94a3b8', size: 11 },
        xaxis: { title: 'False positive rate', range: [0,1], gridcolor: 'rgba(148,163,184,0.12)', zerolinecolor: 'rgba(148,163,184,0.2)' },
        yaxis: { title: 'True positive rate', range: [0,1], gridcolor: 'rgba(148,163,184,0.12)', zerolinecolor: 'rgba(148,163,184,0.2)' },
        legend: { orientation: 'h', x: 0.02, y: 1.12, bgcolor: 'rgba(0,0,0,0)' },
        showlegend: true,
      };
      window.Plotly.react(div, traces, layout, { displayModeBar: false, responsive: true });
    }

    fs.sel.addEventListener('change', () => { state.featureSet = fs.sel.value; renderControls(); refresh(); });
    md.sel.addEventListener('change', () => { state.model = md.sel.value; renderControls(); refresh(); });
    ds.sel.addEventListener('change', () => { state.dataset = ds.sel.value; refresh(); });
    const slider = host.querySelector('#te-threshold');
    slider.addEventListener('input', () => { state.tIdx = parseInt(slider.value, 10) || 0; refresh(); });
    host.querySelector('#te-use-youden').addEventListener('click', () => {
      const e = currentEntry();
      if (!e || !e.threshold_curve) return;
      const yt = e.youden_threshold;
      if (yt == null) return;
      const thrs = e.threshold_curve.thresholds;
      let best = 0, bestDiff = Infinity;
      for (let i = 0; i < thrs.length; i++) {
        const d = Math.abs(thrs[i] - yt);
        if (d < bestDiff) { bestDiff = d; best = i; }
      }
      state.tIdx = best;
      refresh();
    });

    renderControls();
    refresh();
  }

  function init() {
    const host = document.getElementById(HOST_ID);
    if (!host) return;
    host.innerHTML = '<div class="te-loading subtle small">Loading interactive threshold explorer…</div>';
    fetch(DATA_URL).then(r => r.json()).then(d => build(host, d)).catch(err => {
      host.innerHTML = `<div class="te-error subtle small">Failed to load ${DATA_URL}: ${err.message}</div>`;
    });
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
