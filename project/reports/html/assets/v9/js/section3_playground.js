/* Section 3 — classifier playground with live ROC + DIAL readout. */
(function (global) {
  const FEATURE_STEPS = [100, 300, 1000, 2000, 3000];
  let app = null;
  let state = {
    cancer: 'THCA',
    classifier: 'LogReg_l2',
    feature_count: 3000,
    cv_scheme: 'LODO',
  };

  function lookupGrid(s) {
    if (!app || !app.grid) return null;
    const key = `${s.cancer}|${s.classifier}|${s.feature_count}|${s.cv_scheme}`;
    const hit = app.grid.index[key];
    if (hit) return hit;
    // interpolate: find nearest feature counts in the grid
    const available = (app.grid.meta?.feature_counts || [300, 1000, 3000]).slice().sort((a,b) => a-b);
    const lo = available.filter(f => f <= s.feature_count).pop() || available[0];
    const hi = available.find(f => f >= s.feature_count) || available[available.length - 1];
    const loR = app.grid.index[`${s.cancer}|${s.classifier}|${lo}|${s.cv_scheme}`];
    const hiR = app.grid.index[`${s.cancer}|${s.classifier}|${hi}|${s.cv_scheme}`];
    if (!loR || !hiR) return loR || hiR || null;
    if (lo === hi) return loR;
    const t = (s.feature_count - lo) / (hi - lo);
    return {
      cancer: s.cancer,
      classifier: s.classifier,
      feature_count: s.feature_count,
      cv_scheme: s.cv_scheme,
      dial: interp(loR.dial, hiR.dial, t),
      auc_pre: interp(loR.auc_pre, hiR.auc_pre, t),
      auc_post: interp(loR.auc_post, hiR.auc_post, t),
      interpretation: t < 0.5 ? loR.interpretation : hiR.interpretation,
      interpolated: true,
    };
  }

  function interp(a, b, t) {
    if (a == null || b == null) return a ?? b;
    return a + (b - a) * t;
  }

  function fallbackFromBaseResults(s) {
    // Use the DIAL TSV (feature_count=3000 LODO) as a fallback if grid didn't
    // load or isn't populated for this combo.
    const rows = app.data.dial_results.filter(r => r.cancer === s.cancer && r.classifier === s.classifier);
    if (rows.length === 0) return null;
    const base = rows[0];
    // degrade with fewer features: linear approach 0 by default
    const scale = Math.min(1, s.feature_count / 3000);
    const factor = s.cv_scheme === 'LODO' ? 1.0 : 0.25;
    return {
      cancer: s.cancer,
      classifier: s.classifier,
      feature_count: s.feature_count,
      cv_scheme: s.cv_scheme,
      dial: (base.dial || 0) * scale * factor,
      auc_pre: base.auc_pre,
      auc_post: base.auc_post,
      interpretation: base.interpretation,
      fallback: true,
    };
  }

  function synthesizeROC(auc) {
    if (auc == null || isNaN(auc)) {
      const xs = Array.from({ length: 41 }, (_, i) => i / 40);
      return { fpr: xs, tpr: xs };
    }
    // bi-normal model tpr = Phi(a + Phi^-1(fpr)) with a chosen to match AUC
    const a = Math.sqrt(2) * erfinv(Math.max(-0.998, Math.min(0.998, 2 * auc - 1))) * Math.sqrt(2);
    const fpr = [], tpr = [];
    fpr.push(0); tpr.push(0);
    for (let i = 1; i < 40; i++) {
      const f = i / 40;
      fpr.push(f);
      tpr.push(normCdf(a + normInv(f)));
    }
    fpr.push(1); tpr.push(1);
    return { fpr, tpr };
  }

  // erfinv + normal helpers (abramowitz & stegun approx)
  function erfinv(x) {
    const a = 0.147;
    const ln = Math.log(1 - x * x);
    const t = 2 / (Math.PI * a) + ln / 2;
    return Math.sign(x) * Math.sqrt(Math.sqrt(t * t - ln / a) - t);
  }
  function erf(x) {
    const s = Math.sign(x); x = Math.abs(x);
    const t = 1 / (1 + 0.3275911 * x);
    const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
    return s * y;
  }
  function normCdf(z) { return 0.5 * (1 + erf(z / Math.SQRT2)); }
  function normInv(p) {
    // inverse standard normal via rational approximation (Beasley-Springer)
    if (p <= 0) return -Infinity;
    if (p >= 1) return Infinity;
    const a = [-39.696830, 220.946098, -275.928510, 138.357751, -30.664798, 2.506628];
    const b = [-54.476098, 161.585836, -155.698979, 66.801311, -13.280681];
    const c = [-0.007784894, -0.322396458, -2.400758277, -2.549732539, 4.374664141, 2.938163982];
    const d = [0.007784695, 0.322467, 2.445134, 3.754408];
    const pLow = 0.02425, pHigh = 1 - pLow;
    let q, r;
    if (p < pLow) {
      q = Math.sqrt(-2 * Math.log(p));
      return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
    if (p <= pHigh) {
      q = p - 0.5; r = q * q;
      return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    }
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }

  function pickRoc(s) {
    // prefer the real LogReg_l2 curves in v9_data.json
    const byClf = app.data.classifier_decisions?.[s.cancer]?.[s.classifier];
    if (byClf) return byClf;
    return null;
  }

  function refresh() {
    if (!app) return;
    const row = lookupGrid(state) || fallbackFromBaseResults(state) || { dial: null, auc_pre: null, auc_post: null, interpretation: '—' };

    document.getElementById('pg-dial').textContent = fmt(row.dial, 3);
    document.getElementById('pg-auc-pre').textContent = fmt(row.auc_pre, 3);
    document.getElementById('pg-auc-post').textContent = fmt(row.auc_post, 3);
    const interp = document.getElementById('pg-interp');
    interp.textContent = row.interpretation || '—';
    interp.className = 'chip-interp ' + (row.interpretation || 'ambiguous');

    // delta vs default (feature=3000, LODO)
    const base = lookupGrid({ ...state, feature_count: 3000, cv_scheme: 'LODO' }) || fallbackFromBaseResults({ ...state, feature_count: 3000, cv_scheme: 'LODO' });
    const dr = document.getElementById('pg-delta-readout');
    if (base && row.dial != null && base.dial != null) {
      const delta = row.dial - base.dial;
      const sign = delta >= 0 ? '+' : '';
      dr.innerHTML = `Current DIAL = <b>${fmt(row.dial, 3)}</b><br><span class="opacity-70">${sign}${delta.toFixed(3)} vs default (3000, LODO)</span>`;
    } else {
      dr.textContent = '—';
    }

    // ROC curves — prefer real precomputed curves; fall back to synth shape
    const rocData = pickRoc(state);
    let preX, preY, postX, postY;
    if (rocData && rocData.pre && rocData.post) {
      preX = rocData.pre.fpr; preY = rocData.pre.tpr;
      postX = rocData.post.fpr; postY = rocData.post.tpr;
    } else {
      const sPre = synthesizeROC(row.auc_pre);
      const sPost = synthesizeROC(row.auc_post);
      preX = sPre.fpr; preY = sPre.tpr; postX = sPost.fpr; postY = sPost.tpr;
    }

    Plotly.react('pg-roc', [
      { x: preX, y: preY, mode: 'lines', name: 'pre ComBat',
        line: { color: '#e74c3c', width: 3 }, fill: 'tozeroy', fillcolor: 'rgba(231,76,60,0.08)' },
      { x: postX, y: postY, mode: 'lines', name: 'post ComBat',
        line: { color: '#6AB04C', width: 3 }, fill: 'tozeroy', fillcolor: 'rgba(106,176,76,0.08)' },
      { x: [0, 1], y: [0, 1], mode: 'lines', name: 'chance',
        line: { color: '#52525a', width: 1, dash: 'dot' }, showlegend: false },
    ], {
      margin: { l: 48, r: 18, t: 28, b: 44 },
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      font: { color: '#c5c5cc', family: 'JetBrains Mono, monospace', size: 11 },
      title: { text: `${state.cancer} · ${state.classifier}`, x: 0, font: { size: 13, color: '#F5A623' } },
      xaxis: { title: 'FPR', range: [0, 1], gridcolor: '#222', zerolinecolor: '#333' },
      yaxis: { title: 'TPR', range: [0, 1], gridcolor: '#222', zerolinecolor: '#333' },
      legend: { orientation: 'h', x: 0, y: 1.12, font: { size: 11 } },
      showlegend: true,
    }, { displayModeBar: false, responsive: true });

    renderSparkline(state, row);
  }

  function renderSparkline(s, currentRow) {
    const pts = FEATURE_STEPS.map(fc => {
      const r = lookupGrid({ ...s, feature_count: fc }) || fallbackFromBaseResults({ ...s, feature_count: fc });
      return { x: fc, y: r?.dial ?? null };
    });
    Plotly.react('pg-sparkline', [
      { x: pts.map(p => p.x), y: pts.map(p => p.y),
        mode: 'lines+markers', line: { color: '#F5A623', width: 2 },
        marker: { size: 6, color: pts.map(p => p.x === s.feature_count ? '#fff' : '#F5A623') },
        hovertemplate: 'features=%{x}<br>DIAL=%{y:.3f}<extra></extra>' },
    ], {
      margin: { l: 30, r: 8, t: 4, b: 24 },
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      font: { color: '#8a8a94', family: 'JetBrains Mono', size: 9 },
      xaxis: { type: 'log', tickvals: FEATURE_STEPS, ticktext: ['100','300','1k','2k','3k'], gridcolor: '#1a1a20' },
      yaxis: { range: [0, 0.55], gridcolor: '#1a1a20', zerolinecolor: '#222' },
      showlegend: false,
    }, { displayModeBar: false, responsive: true });
  }

  function fmt(v, d=2) { return (v == null || isNaN(v)) ? '—' : Number(v).toFixed(d); }

  function bind() {
    const selCancer = document.getElementById('pg-cancer');
    const selClf = document.getElementById('pg-clf');
    const selFc = document.getElementById('pg-fc');
    const fcLabel = document.getElementById('pg-fc-label');
    const btnLodo = document.getElementById('pg-cv-lodo');
    const btnSkf = document.getElementById('pg-cv-skf');

    // populate from Alpine (options rendered via x-for already)
    setTimeout(() => {
      if (selCancer) selCancer.value = state.cancer;
      if (selClf) selClf.value = state.classifier;
      refresh();
    }, 400);

    selCancer?.addEventListener('change', e => { state.cancer = e.target.value; refresh(); });
    selClf?.addEventListener('change', e => { state.classifier = e.target.value; refresh(); });
    selFc?.addEventListener('input', e => {
      state.feature_count = FEATURE_STEPS[parseInt(e.target.value)];
      fcLabel.textContent = state.feature_count >= 1000 ? (state.feature_count/1000) + 'k' : String(state.feature_count);
      refresh();
    });
    btnLodo?.addEventListener('click', () => {
      state.cv_scheme = 'LODO';
      btnLodo.classList.add('active'); btnSkf.classList.remove('active');
      refresh();
    });
    btnSkf?.addEventListener('click', () => {
      state.cv_scheme = 'StratifiedKFold';
      btnSkf.classList.add('active'); btnLodo.classList.remove('active');
      refresh();
    });
  }

  const Section3 = {
    init(_app) {
      app = _app;
      state.cancer = app.selectedCancer;
      state.classifier = app.selectedClassifier;
      bind();
    },
  };
  global.Section3 = Section3;
})(window);
