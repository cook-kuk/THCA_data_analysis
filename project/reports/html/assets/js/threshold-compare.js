/*
 * threshold-compare.js — cohort compare mode for the Interactive threshold explorer.
 *
 * Augments 07_ml_baseline.html only. Adds a "비교 모드" toggle next to the
 * #threshold-explorer-section panel that, when enabled, mounts a second
 * independent threshold-explorer instance (feature_set, model, dataset) in a
 * side-by-side column. Reuses the same data source `model_curves_v2.json`.
 *
 * Does NOT import or touch threshold-explorer.js. It implements its own
 * compact viewer for side B so the two panels stay decoupled.
 */
(function () {
  'use strict';
  if (window.__thyroidThresholdCompareLoaded) return;
  window.__thyroidThresholdCompareLoaded = true;

  var SECTION_ID = 'threshold-explorer-section';
  var DATA_URL = '../assets/data/model_curves_v2.json';

  function fmt(n, digits) {
    if (n === null || n === undefined || Number.isNaN(n)) return '—';
    if (typeof n !== 'number') return String(n);
    return n.toFixed(digits == null ? 3 : digits);
  }

  function buildBuddy(host, data, initialKey) {
    host.innerHTML = '';
    host.classList.add('te-root', 'thyroid-compare-root');
    host.style.marginTop = '12px';
    host.innerHTML = [
      '<div class="te-controls">',
      '  <div class="te-control-row"></div>',
      '  <div class="te-threshold-row">',
      '    <label for="tec-threshold" class="te-label te-slider-label"><span>Decision threshold</span></label>',
      '    <input type="range" id="tec-threshold" min="0" max="100" step="1" value="50" class="te-slider">',
      '    <span id="tec-threshold-value" class="te-threshold-badge mono">0.50</span>',
      '  </div>',
      '</div>',
      '<div class="te-grid">',
      '  <div class="te-metrics" aria-live="polite">',
      '    <div class="te-metric"><div class="te-mlabel">AUC</div><div class="te-mvalue mono" id="tec-m-auc">—</div></div>',
      '    <div class="te-metric"><div class="te-mlabel">PR-AUC</div><div class="te-mvalue mono" id="tec-m-prauc">—</div></div>',
      '    <div class="te-metric"><div class="te-mlabel">Precision</div><div class="te-mvalue mono" id="tec-m-prec">—</div></div>',
      '    <div class="te-metric"><div class="te-mlabel">Recall</div><div class="te-mvalue mono" id="tec-m-rec">—</div></div>',
      '    <div class="te-metric"><div class="te-mlabel">F1</div><div class="te-mvalue mono" id="tec-m-f1">—</div></div>',
      '    <div class="te-metric"><div class="te-mlabel">Balanced acc</div><div class="te-mvalue mono" id="tec-m-bacc">—</div></div>',
      '  </div>',
      '  <div id="tec-chart" class="te-chart" role="img" aria-label="ROC panel B"></div>',
      '</div>',
      '<p class="subtle small">사이드-바이-사이드 비교용 독립 탐색기 (B panel). 좌측 패널과 다른 조합을 선택하세요.</p>'
    ].join('');

    function makeSelect(id, label) {
      var wrap = document.createElement('label');
      wrap.className = 'te-label';
      wrap.setAttribute('for', id);
      wrap.innerHTML = '<span>' + label + '</span>';
      var sel = document.createElement('select');
      sel.id = id;
      sel.className = 'te-select';
      wrap.appendChild(sel);
      return { wrap: wrap, sel: sel };
    }
    function fillOptions(sel, opts, selected) {
      sel.innerHTML = '';
      opts.forEach(function (o) {
        var opt = document.createElement('option');
        opt.value = o; opt.textContent = o;
        if (o === selected) opt.selected = true;
        sel.appendChild(opt);
      });
    }

    var row = host.querySelector('.te-control-row');
    var fs = makeSelect('tec-fs', 'Feature set');
    var md = makeSelect('tec-md', 'Model');
    var ds = makeSelect('tec-ds', 'Eval cohort');
    row.appendChild(fs.wrap); row.appendChild(md.wrap); row.appendChild(ds.wrap);

    var fsKeys = Object.keys(data);
    var state = {
      featureSet: initialKey && data[initialKey] ? initialKey : fsKeys[0],
      model: null,
      dataset: null,
      tIdx: 50
    };

    function renderControls() {
      fillOptions(fs.sel, fsKeys, state.featureSet);
      var models = Object.keys(data[state.featureSet] || {});
      // Prefer a model different from panel A to make the comparison meaningful.
      var aModel = null;
      try { aModel = (document.getElementById('te-md') || {}).value; } catch (e) {}
      if (!models.includes(state.model)) {
        if (aModel && models.indexOf(aModel) >= 0 && models.length > 1) {
          state.model = models.find(function (m) { return m !== aModel; }) || models[0];
        } else {
          state.model = models.includes('RandomForest') ? 'RandomForest' : models[0];
        }
      }
      fillOptions(md.sel, models, state.model);
      var evals = Object.keys(data[state.featureSet]?.[state.model] || {});
      if (!evals.includes(state.dataset)) state.dataset = evals[0];
      fillOptions(ds.sel, evals, state.dataset);
    }

    function currentEntry() {
      return (data[state.featureSet] || {})[state.model] && data[state.featureSet][state.model][state.dataset] || null;
    }

    function refresh() {
      var e = currentEntry();
      var $ = function (id) { return host.querySelector('#' + id); };
      var ids = ['tec-m-auc','tec-m-prauc','tec-m-prec','tec-m-rec','tec-m-f1','tec-m-bacc'];
      if (!e) { ids.forEach(function (id) { $(id).textContent = '—'; }); return; }
      var tc = e.threshold_curve || { thresholds: [], precision: [], recall: [], f1: [], balanced_accuracy: [] };
      var idx = Math.min(Math.max(0, state.tIdx), Math.max(0, tc.thresholds.length - 1));
      var t = tc.thresholds[idx];
      var slider = $('tec-threshold');
      slider.max = Math.max(0, tc.thresholds.length - 1);
      slider.value = idx;
      $('tec-threshold-value').textContent = (t !== undefined) ? t.toFixed(3) : '—';
      $('tec-m-auc').textContent = fmt(e.auc);
      $('tec-m-prauc').textContent = fmt(e.pr_auc);
      $('tec-m-prec').textContent = fmt(tc.precision && tc.precision[idx]);
      $('tec-m-rec').textContent = fmt(tc.recall && tc.recall[idx]);
      $('tec-m-f1').textContent = fmt(tc.f1 && tc.f1[idx]);
      $('tec-m-bacc').textContent = fmt(tc.balanced_accuracy && tc.balanced_accuracy[idx]);
      renderChart(e, t);
    }

    function renderChart(entry, tCurrent) {
      var div = host.querySelector('#tec-chart');
      if (!window.Plotly) { div.textContent = 'Plotly not loaded.'; return; }
      var fpr = (entry.roc && entry.roc.fpr) || [];
      var tpr = (entry.roc && entry.roc.tpr) || [];
      var scores = entry.scores || [];
      var y = entry.y_true || [];
      var markerFpr = null, markerTpr = null;
      if (scores.length && scores.length === y.length) {
        var tp=0, fp=0, fn=0, tn=0;
        for (var i=0;i<scores.length;i++) {
          var p = scores[i] >= tCurrent ? 1 : 0;
          if (p === 1 && y[i] === 1) tp++;
          else if (p === 1 && y[i] === 0) fp++;
          else if (p === 0 && y[i] === 1) fn++;
          else tn++;
        }
        markerFpr = (fp+tn) ? fp/(fp+tn) : 0;
        markerTpr = (tp+fn) ? tp/(tp+fn) : 0;
      }
      var traces = [
        { x: fpr, y: tpr, mode: 'lines', name: 'ROC (B)', line: { color: '#f59e0b', width: 2.5 } },
        { x: [0,1], y: [0,1], mode: 'lines', name: 'Chance', line: { color: 'rgba(148,163,184,0.5)', dash: 'dash' }, hoverinfo: 'skip' }
      ];
      if (markerFpr !== null) {
        traces.push({ x: [markerFpr], y: [markerTpr], mode: 'markers', name: 't=' + tCurrent.toFixed(3), marker: { color: '#ef4444', size: 12, line: { color: '#fff', width: 2 } } });
      }
      var layout = {
        margin: { l: 44, r: 12, t: 16, b: 36 },
        height: 280,
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        font: { family: 'Inter, ui-sans-serif, system-ui', color: '#94a3b8', size: 11 },
        xaxis: { title: 'FPR', range: [0,1], gridcolor: 'rgba(148,163,184,0.12)' },
        yaxis: { title: 'TPR', range: [0,1], gridcolor: 'rgba(148,163,184,0.12)' },
        legend: { orientation: 'h', x: 0.02, y: 1.12, bgcolor: 'rgba(0,0,0,0)' },
        showlegend: true
      };
      window.Plotly.react(div, traces, layout, { displayModeBar: false, responsive: true });
    }

    fs.sel.addEventListener('change', function () { state.featureSet = fs.sel.value; state.model = null; state.dataset = null; renderControls(); refresh(); });
    md.sel.addEventListener('change', function () { state.model = md.sel.value; state.dataset = null; renderControls(); refresh(); });
    ds.sel.addEventListener('change', function () { state.dataset = ds.sel.value; refresh(); });
    var slider = host.querySelector('#tec-threshold');
    slider.addEventListener('input', function () { state.tIdx = parseInt(slider.value, 10) || 0; refresh(); });

    renderControls();
    refresh();
  }

  var loadedData = null;
  function getData() {
    if (loadedData) return Promise.resolve(loadedData);
    return fetch(DATA_URL).then(function (r) { return r.json(); }).then(function (d) { loadedData = d; return d; });
  }

  function ensureToggle(section) {
    if (!section) return;
    if (section.querySelector('.thyroid-compare-toggle')) return;

    var controls = document.createElement('div');
    controls.className = 'thyroid-compare-controls';
    controls.setAttribute('data-print-hide', '1');
    controls.style.cssText = 'display:flex;align-items:center;gap:10px;margin:8px 0 12px 0';

    var label = document.createElement('label');
    label.className = 'te-label';
    label.style.cssText = 'display:inline-flex;align-items:center;gap:8px;cursor:pointer;font-size:.82rem';
    var cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.className = 'thyroid-compare-toggle';
    cb.id = 'thyroid-compare-toggle';
    var span = document.createElement('span');
    span.textContent = '비교 모드 (compare mode)';
    label.appendChild(cb); label.appendChild(span);
    controls.appendChild(label);

    var hint = document.createElement('span');
    hint.className = 'subtle small';
    hint.textContent = '두 번째 탐색기를 옆에 열어 (feature_set, model, dataset) 조합을 비교합니다.';
    controls.appendChild(hint);

    // Insert before the threshold-explorer host so the toggle lives above it.
    var host = section.querySelector('#threshold-explorer');
    if (host && host.parentNode) {
      host.parentNode.insertBefore(controls, host);
    } else {
      section.appendChild(controls);
    }

    cb.addEventListener('change', function () {
      if (cb.checked) enableCompare(section);
      else disableCompare(section);
    });
  }

  function enableCompare(section) {
    var host = section.querySelector('#threshold-explorer');
    if (!host) return;
    // Create a flex row wrapper on first activation
    var wrap = section.querySelector('.thyroid-compare-grid');
    var buddy = section.querySelector('#threshold-explorer-compare');
    if (!wrap) {
      wrap = document.createElement('div');
      wrap.className = 'thyroid-compare-grid';
      wrap.style.cssText = 'display:grid;grid-template-columns:1fr 1fr;gap:16px;align-items:start';
      host.parentNode.insertBefore(wrap, host);
      var colA = document.createElement('div');
      colA.className = 'thyroid-compare-col thyroid-compare-col-a';
      colA.appendChild(host);
      wrap.appendChild(colA);
    }
    if (!buddy) {
      var colB = document.createElement('div');
      colB.className = 'thyroid-compare-col thyroid-compare-col-b';
      buddy = document.createElement('div');
      buddy.id = 'threshold-explorer-compare';
      colB.appendChild(buddy);
      wrap.appendChild(colB);
      getData().then(function (d) { buildBuddy(buddy, d, null); }).catch(function (err) {
        buddy.innerHTML = '<div class="subtle small">Failed to load compare data: ' + err.message + '</div>';
      });
    } else {
      buddy.parentElement && (buddy.parentElement.style.display = '');
    }

    // responsive: stack on narrow viewports
    applyResponsive(wrap);
    window.addEventListener('resize', function () { applyResponsive(wrap); }, { passive: true });
  }

  function applyResponsive(wrap) {
    if (!wrap) return;
    if (window.innerWidth < 900) {
      wrap.style.gridTemplateColumns = '1fr';
    } else {
      wrap.style.gridTemplateColumns = '1fr 1fr';
    }
  }

  function disableCompare(section) {
    var wrap = section.querySelector('.thyroid-compare-grid');
    var host = section.querySelector('#threshold-explorer');
    var buddyCol = wrap && wrap.querySelector('.thyroid-compare-col-b');
    if (buddyCol) buddyCol.style.display = 'none';
    // move host back out of the grid so single-column look is preserved
    if (wrap && host && host.parentElement === wrap.querySelector('.thyroid-compare-col-a')) {
      wrap.parentNode.insertBefore(host, wrap);
      wrap.parentNode.removeChild(wrap);
    }
  }

  function init() {
    var section = document.getElementById(SECTION_ID);
    if (!section) return;
    ensureToggle(section);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init, { once: true });
  } else {
    init();
  }

  window.ThyroidDash = window.ThyroidDash || {};
  window.ThyroidDash.thresholdCompare = { enable: function () {
    var s = document.getElementById(SECTION_ID);
    var cb = s && s.querySelector('.thyroid-compare-toggle');
    if (cb && !cb.checked) { cb.checked = true; cb.dispatchEvent(new Event('change')); }
  }};
})();
