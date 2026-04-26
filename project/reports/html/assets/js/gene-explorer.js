/* Gene explorer — page 10.
   Uses two new endpoints:
     assets/data/gene_universe.json   {symbols:[...], curated:[...]}
     assets/data/gene_expression.json {gene: {dataset: [values, ...]}}
   Autocomplete runs against the full symbol universe; violins come from
   the pre-computed expression cache (curated subset). If a searched gene
   is in the universe but not cached, show a clear "not cached — request via
   ..." message instead of the generic failure. */
(function(){
  'use strict';

  function installAutocomplete(input, symbols, curated, onPick){
    const cset = new Set(curated);
    const box = document.createElement('div');
    box.className = 'gx-suggest';
    box.style.cssText = 'position:absolute;top:100%;left:0;right:0;max-height:320px;overflow:auto;background:var(--card-solid,#0b1a3a);border:1px solid var(--line,rgba(148,163,184,.2));border-radius:12px;box-shadow:var(--shadow-lg,0 10px 30px rgba(0,0,0,.35));margin-top:6px;display:none;z-index:40';
    const wrap = document.createElement('div'); wrap.style.cssText='position:relative';
    input.parentNode.insertBefore(wrap, input); wrap.appendChild(input); wrap.appendChild(box);

    // Priority-sorted filter: prefix > substring, curated > non-curated
    function rank(g, q){
      if (!q) return cset.has(g) ? 1 : 3;
      if (g === q) return 0;
      if (g.startsWith(q)) return cset.has(g) ? 1 : 2;
      return cset.has(g) ? 3 : 4;
    }

    function render(q){
      q = (q||'').toUpperCase().trim();
      const hits = [];
      const Q = q.length;
      for (const g of symbols){
        if (Q && !g.includes(q)) continue;
        hits.push(g);
        if (hits.length >= 500) break;
      }
      hits.sort((a, b) => rank(a, q) - rank(b, q) || a.localeCompare(b));
      const shown = hits.slice(0, 20);
      box.innerHTML = shown.map(g => {
        const isCurated = cset.has(g);
        const badge = isCurated
          ? '<span style="float:right;font-size:.65rem;padding:2px 7px;border-radius:999px;background:rgba(20,184,166,.18);color:#5eead4">violin</span>'
          : '<span style="float:right;font-size:.65rem;padding:2px 7px;border-radius:999px;background:rgba(148,163,184,.12);color:#94a3b8">symbol</span>';
        return `<a href="#" data-gene="${g}" style="display:block;padding:8px 12px;color:var(--ink,#e2e8f0);font-size:.9rem;text-decoration:none;border-bottom:1px solid rgba(148,163,184,.08)"><strong>${g}</strong>${badge}</a>`;
      }).join('');
      box.style.display = shown.length ? 'block' : 'none';
      box.querySelectorAll('a').forEach(a => {
        a.addEventListener('click', e => {
          e.preventDefault();
          input.value = a.dataset.gene;
          box.style.display = 'none';
          onPick(a.dataset.gene);
        });
      });
    }
    input.addEventListener('input', e => render(e.target.value));
    input.addEventListener('focus', e => render(e.target.value));
    document.addEventListener('click', e => { if(!wrap.contains(e.target)) box.style.display='none'; });
  }

  async function loadUniverse(){
    try {
      const r = await fetch('../assets/data/gene_universe.json');
      if (!r.ok) throw new Error(r.status);
      return await r.json();
    } catch (e) {
      return { symbols: [], curated: [] };
    }
  }

  async function loadExpression(){
    try {
      const r = await fetch('../assets/data/gene_expression.json');
      if (!r.ok) throw new Error(r.status);
      return await r.json();
    } catch (e) {
      return {};
    }
  }

  function renderMissing(container, gene, reason){
    container.innerHTML = `<div class="subtle" style="padding:30px;text-align:center">
      <strong style="display:block;color:var(--ink-strong,#e2e8f0);margin-bottom:6px">${gene||''}</strong>
      ${reason}
    </div>`;
  }

  function plotGene(container, cache, gene){
    if(!window.Plotly){ container.textContent='Plotly 로드 실패'; return; }
    const G = (gene||'').toUpperCase().trim();
    const entry = cache[G];
    if (!entry){
      renderMissing(container, G,
        '이 유전자는 universe에는 있지만 violin 캐시(277 curated)에 포함되지 않습니다. 주요 panels (TDS16 · TierA67 · BRS71) + 상위 novel 후보 + driver genes 위주로 캐싱되어 있습니다. 다른 유전자는 필요 시 <code>build_gene_explorer_cache.py</code>에서 curated set을 확장하세요.');
      return;
    }
    const traces = Object.keys(entry).sort().map(ds => ({
      type: 'violin', name: ds, legendgroup: ds,
      y: entry[ds],
      x: entry[ds].map(() => ds),
      box: { visible: true }, meanline: { visible: true },
      points: 'all', jitter: 0.45, pointpos: 0,
      marker: { size: 4, opacity: 0.55 },
      line: { width: 1 },
      hovertemplate: `<b>${ds}</b><br>${G} log2 = %{y:.2f}<extra></extra>`,
    }));
    const layout = Object.assign({}, (window.THYROID_PLOTLY_TEMPLATE||{}).layout||{}, {
      title: { text: G + ' expression across datasets (log2)', font: { size: 14 } },
      yaxis: { title: 'log2 expression', zeroline: false },
      xaxis: { title: '', tickangle: -25 },
      showlegend: false, height: 460,
      margin: { l: 56, r: 14, t: 44, b: 64 },
    });
    Plotly.react(container, traces, layout, { displaylogo: false, responsive: true });
  }

  async function init(){
    const input = document.getElementById('gx-gene-search');
    const target = document.getElementById('gx-gene-plot');
    const status = document.getElementById('gx-gene-status');
    if(!input || !target) return;
    status && (status.textContent='유전자 목록을 불러오는 중…');
    const [uni, cache] = await Promise.all([loadUniverse(), loadExpression()]);
    const symbols = uni.symbols || [];
    const curated = uni.curated || [];
    if (!symbols.length){
      status && (status.textContent = '유전자 목록을 불러오지 못했습니다.');
      return;
    }
    status && (status.textContent =
      `${symbols.length.toLocaleString()}개 symbol 로드됨 (${curated.length}개 violin 가능). 이름을 입력하거나 Enter로 선택.`);
    installAutocomplete(input, symbols, curated, gene => plotGene(target, cache, gene));
    input.addEventListener('keydown', e => {
      if(e.key === 'Enter'){
        const v = (input.value || '').toUpperCase().trim();
        if (v) plotGene(target, cache, v);
      }
    });
    // Auto-plot a sensible default on page load
    const defaultGene = curated.find(g => g === 'TACSTD2') || curated[0];
    if (defaultGene){
      input.value = defaultGene;
      plotGene(target, cache, defaultGene);
    }
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
