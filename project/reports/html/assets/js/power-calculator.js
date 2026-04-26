/* Power & sample-size calculator — pure vanilla JS, depends only on vendored Plotly.
 * Two-sample t-test (normal approximation):
 *   n_per_group = ((z_{α/2} + z_β)² * (σ1² + σ2²)) / d²
 * With σ=1 (standardized effect d = Cohen's d).
 *
 * Also back-computes *achievable power* given fixed sample sizes:
 *   z_β = sqrt(n * d² / 2) - z_{α/2}
 *   power = Φ(z_β)
 *
 * Fixed cohort sizes (per dashboard):
 *   TCGA: 333, GSE27155: 72, GSE126698: 12
 *
 * Renders live Plotly curve: power vs d for each cohort.
 */
(function(){
  'use strict';
  // --------- Numerical helpers --------------------------------------
  // erf via Abramowitz & Stegun 7.1.26 (max error ≈ 1.5e-7)
  function erf(x){
    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x);
    const a1= 0.254829592, a2=-0.284496736, a3= 1.421413741;
    const a4=-1.453152027, a5= 1.061405429, p = 0.3275911;
    const t = 1 / (1 + p * x);
    const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
    return sign * y;
  }
  function Phi(z){ return 0.5 * (1 + erf(z / Math.SQRT2)); }
  // Inverse normal via Beasley-Springer-Moro
  function invPhi(p){
    if (p <= 0 || p >= 1) return NaN;
    const a=[-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00];
    const b=[-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01];
    const c=[-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00];
    const d=[7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00];
    const plow = 0.02425, phigh = 1 - plow;
    let q, r;
    if (p < plow){ q = Math.sqrt(-2*Math.log(p));
      return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    } else if (p <= phigh){ q = p - 0.5; r = q*q;
      return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    } else { q = Math.sqrt(-2*Math.log(1-p));
      return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
  }

  function requiredN(d, alpha, power){
    const zA = invPhi(1 - alpha/2);
    const zB = invPhi(power);
    const num = Math.pow(zA + zB, 2) * 2; // σ1²+σ2² = 2 for σ=1
    return Math.ceil(num / Math.pow(d, 2));
  }
  function achievablePower(d, alpha, nPerGroup){
    if (nPerGroup < 2) return 0;
    const zA = invPhi(1 - alpha/2);
    const zB = Math.sqrt(nPerGroup * d * d / 2) - zA;
    return Math.max(0, Math.min(1, Phi(zB)));
  }

  // --------- Cohort sizes (per_group approximation = total/2) --------
  const COHORTS = [
    { name: 'TCGA',       total: 333, nPerGroup: 166 },
    { name: 'GSE27155',   total:  72, nPerGroup:  36 },
    { name: 'GSE126698',  total:  12, nPerGroup:   6 }
  ];

  // --------- State --------------------------------------------------
  const state = { d: 0.5, alpha: 0.05, power: 0.80 };

  // --------- Render -------------------------------------------------
  function fmt(x, n){ return (Number.isFinite(x) ? x.toFixed(n) : '—'); }

  function render(){
    const reqN = requiredN(state.d, state.alpha, state.power);

    // Update controls display
    document.getElementById('pc-d-val').textContent   = state.d.toFixed(2);
    document.getElementById('pc-a-val').textContent   = state.alpha.toFixed(3);
    document.getElementById('pc-p-val').textContent   = state.power.toFixed(2);

    document.getElementById('pc-reqN').textContent = Number.isFinite(reqN) ? reqN.toLocaleString() : '—';

    // Cohort power rows
    const tbody = document.getElementById('pc-cohort-body');
    tbody.innerHTML = '';
    COHORTS.forEach(c => {
      const pw = achievablePower(state.d, state.alpha, c.nPerGroup);
      const cls = pw >= 0.8 ? 'pc-ok' : (pw >= 0.5 ? 'pc-warn' : 'pc-bad');
      const mark = pw >= 0.8 ? '✓ powered' : (pw >= 0.5 ? '△ marginal' : '✗ under-powered');
      tbody.insertAdjacentHTML('beforeend',
        `<tr><td>${c.name}</td><td>${c.total}</td><td>${c.nPerGroup}</td><td class="${cls}">${(pw*100).toFixed(1)}%</td><td class="${cls}">${mark}</td></tr>`);
    });

    // Biomarker fraction message — of 2,773 novel biomarkers, how many reach power≥0.8 at this setup,
    // using d=state.d as a proxy (actual d distribution not known client-side; assume user tunes).
    const minD = (function(){
      // minimum d reliably detected at GSE27155 size, given current alpha & target power
      const zA = invPhi(1 - state.alpha/2);
      const zB = invPhi(0.8);
      const n = 36; // GSE27155 per_group
      return Math.sqrt(2 * Math.pow(zA + zB, 2) / n);
    })();
    const msg = `현재 α=${state.alpha.toFixed(3)} 에서 GSE27155 (n≈72) 로 80% power 를 확보하려면 |d| ≥ ${minD.toFixed(2)} 가 필요합니다. ` +
                `이보다 작은 효과크기를 가진 novel 바이오마커는 외부 검증에서 놓칠 위험이 큽니다.`;
    document.getElementById('pc-minD-msg').innerHTML = msg;

    // Plot: power-vs-d for each cohort, alpha fixed, plus target-power line
    const xs = [];
    for (let d = 0.05; d <= 3.0; d += 0.02) xs.push(+d.toFixed(2));
    const data = COHORTS.map(c => ({
      type: 'scatter', mode: 'lines', name: c.name + ' (n/group=' + c.nPerGroup + ')',
      x: xs, y: xs.map(d => achievablePower(d, state.alpha, c.nPerGroup)),
      line: { width: 2 }
    }));
    // Mark target d, target power
    data.push({
      type: 'scatter', mode: 'markers',
      x: [state.d], y: [state.power],
      marker: { size: 12, color: '#c084fc', symbol: 'star' },
      name: '현재 선택'
    });
    // horizontal target-power line
    const shapes = [{
      type: 'line', x0: 0.05, x1: 3.0, y0: state.power, y1: state.power,
      line: { color: 'rgba(192,132,252,.6)', dash: 'dash', width: 1 }
    }];
    const layout = {
      title: { text: 'Power vs Effect Size (|Cohen\'s d|)', font: { color: '#e2e8f0' } },
      paper_bgcolor: '#07132b', plot_bgcolor: '#0b1a3a',
      font: { color: '#cbd5e1' },
      xaxis: { title: '|Cohen\'s d|', range: [0, 3], gridcolor: 'rgba(148,163,184,.18)' },
      yaxis: { title: 'Statistical power', range: [0, 1.02], gridcolor: 'rgba(148,163,184,.18)' },
      legend: { orientation: 'h', y: 1.08 },
      shapes: shapes,
      margin: { l: 60, r: 20, t: 50, b: 50 }
    };
    if (window.Plotly) {
      Plotly.react('pc-plot', data, layout, { responsive: true, displaylogo: false });
    }
  }

  function bindSlider(id, key, parse){
    const el = document.getElementById(id);
    if (!el) return;
    const handler = () => { state[key] = parse(el.value); render(); };
    el.addEventListener('input', handler);
    el.addEventListener('change', handler);
  }

  document.addEventListener('DOMContentLoaded', function(){
    bindSlider('pc-d',     'd',     v => parseFloat(v));
    bindSlider('pc-alpha', 'alpha', v => parseFloat(v));
    bindSlider('pc-power', 'power', v => parseFloat(v));
    render();
  });
})();
