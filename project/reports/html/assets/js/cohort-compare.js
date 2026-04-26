/* Cohort compare — page 11. Toggle bar-stack grouping by dataset / histology / label. */
(function(){
  'use strict';
  async function init(){
    const target = document.getElementById('cc-toggle-plot');
    const toggles = document.querySelectorAll('[data-cc-group]');
    if(!target || !toggles.length) return;
    let samples = [];
    try { samples = await fetch('../assets/data/sample_master.json').then(r=>r.json()); }
    catch(e){ target.textContent='sample_master 로드 실패'; return; }
    const axisMap = {
      dataset: {x:'dataset', color:'molecular_subtype', label:'데이터셋'},
      histology: {x:'histology_subtype', color:'dataset', label:'조직학 분류'},
      label: {x:'label_confidence', color:'dataset', label:'라벨 신뢰도'}
    };
    function tabulate(rows, xk, ck){
      const xVals = new Set(), cVals = new Set(); const counts = {};
      rows.forEach(r => {
        const x = r[xk] || 'unknown', c = r[ck] || 'unknown';
        xVals.add(x); cVals.add(c);
        counts[c] = counts[c] || {}; counts[c][x] = (counts[c][x]||0)+1;
      });
      return { xs: Array.from(xVals).sort(), cs: Array.from(cVals).sort(), counts };
    }
    function render(mode){
      if(!window.Plotly){ target.textContent='Plotly 로드 실패'; return; }
      const cfg = axisMap[mode] || axisMap.dataset;
      const { xs, cs, counts } = tabulate(samples, cfg.x, cfg.color);
      const traces = cs.map(c => ({
        type:'bar', name:c, x:xs, y:xs.map(x => (counts[c]||{})[x] || 0)
      }));
      const layout = Object.assign({}, (window.THYROID_PLOTLY_TEMPLATE||{}).layout||{}, {
        barmode: 'stack', height: 440,
        title: { text: cfg.label+' 기준 분포', font:{size:14}},
        xaxis: { title: cfg.x }, yaxis: { title: 'Sample count' }
      });
      Plotly.newPlot(target, traces, layout, {displaylogo:false, responsive:true});
    }
    toggles.forEach(btn => {
      btn.addEventListener('click', ()=>{
        toggles.forEach(b => b.setAttribute('aria-pressed', 'false'));
        btn.setAttribute('aria-pressed','true');
        render(btn.dataset.ccGroup);
      });
    });
    const init = document.querySelector('[data-cc-group][aria-pressed="true"]') || toggles[0];
    init && init.click();
  }
  document.addEventListener('DOMContentLoaded', init);
})();
