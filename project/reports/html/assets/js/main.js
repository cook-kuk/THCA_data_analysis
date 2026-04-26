/* THCA Dashboard — main runtime (vanilla JS, no deps) */
(function(){
  'use strict';

  window.ThyroidDash = window.ThyroidDash || {};
  const TD = window.ThyroidDash;

  const $ = (s,r)=> (r||document).querySelector(s);
  const $$ = (s,r)=> Array.from((r||document).querySelectorAll(s));
  const byId = id => document.getElementById(id);

  const prefersReduced = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  const isPageContext = /\/pages\/[^/]+$/.test(window.location.pathname);

  function isLocalDocHref(href){
    return !!href && !/^(?:[a-z]+:|#|\/\/)/i.test(href);
  }

  function normalizeAppHref(href){
    if(!isLocalDocHref(href)) return href;
    const [path, hash=''] = href.split('#');
    if(!path || path.startsWith('../') || path.startsWith('./') || path.startsWith('/')) return href;
    if(!/\.html?$/.test(path) && !path.startsWith('pages/')) return href;

    let nextPath = path;
    if(isPageContext){
      if(path === 'index.html') nextPath = '../index.html';
      else if(path.startsWith('pages/')) nextPath = path.slice('pages/'.length);
    } else if(path !== 'index.html' && !path.startsWith('pages/')) {
      nextPath = 'pages/' + path;
    }

    return hash ? nextPath + '#' + hash : nextPath;
  }

  function initPathFixes(){
    $$('a[href]').forEach(el => {
      const href = el.getAttribute('href');
      const normalized = normalizeAppHref(href);
      if(normalized !== href) el.setAttribute('href', normalized);
    });
  }

  function ensureSiteMap(){
    if (!document.body) return;
    if (window.ThyraiSiteMap) {
      if (window.ThyraiSiteMap.init) window.ThyraiSiteMap.init();
      return;
    }
    const jsHref = isPageContext ? '../assets/js/site-map.js' : 'assets/js/site-map.js';
    if(!document.querySelector('script[src$="site-map.js"]')){
      const script = document.createElement('script');
      script.src = jsHref;
      script.defer = true;
      document.body.appendChild(script);
    }
  }

  TD.copyText = async function(txt){
    try { await navigator.clipboard.writeText(txt); return true; }
    catch(e){ const t=document.createElement('textarea'); t.value=txt; document.body.appendChild(t); t.select(); document.execCommand('copy'); t.remove(); return true; }
  };

  TD.toast = function(msg){
    let n = byId('td-toast');
    if(!n){ n=document.createElement('div'); n.id='td-toast';
      n.style.cssText='position:fixed;bottom:22px;left:50%;transform:translateX(-50%);padding:10px 16px;border-radius:12px;background:rgba(7,19,43,.92);border:1px solid rgba(45,212,191,.35);color:#e2e8f0;font-size:.85rem;z-index:6000;box-shadow:0 20px 48px rgba(3,8,20,.5);opacity:0;transition:opacity .2s ease';
      document.body.appendChild(n);
    }
    n.textContent = msg;
    requestAnimationFrame(()=>{ n.style.opacity='1'; });
    clearTimeout(TD.toast._t);
    TD.toast._t = setTimeout(()=>{ n.style.opacity='0'; }, 2000);
  };

  /* ───── theme ───── */
  window.toggleTheme = function(){
    const b = document.body;
    const next = b.classList.contains('light') ? 'dark' : 'light';
    b.classList.toggle('light', next==='light');
    localStorage.setItem('thyroid-theme', next);
  };
  function initTheme(){
    const saved = localStorage.getItem('thyroid-theme');
    if(saved==='light') document.body.classList.add('light');
  }

  /* ───── mobile drawer ───── */
  function initDrawer(){
    const btn = byId('td-menu-toggle');
    const drawer = byId('td-drawer');
    const backdrop = byId('td-drawer-backdrop');
    if(!btn || !drawer) return;
    const open = ()=>{ drawer.classList.add('open'); backdrop && backdrop.classList.add('open'); btn.setAttribute('aria-expanded','true'); drawer.setAttribute('aria-hidden','false'); };
    const close = ()=>{ drawer.classList.remove('open'); backdrop && backdrop.classList.remove('open'); btn.setAttribute('aria-expanded','false'); drawer.setAttribute('aria-hidden','true'); };
    btn.addEventListener('click',()=> drawer.classList.contains('open') ? close() : open());
    backdrop && backdrop.addEventListener('click', close);
    $$('#td-drawer a').forEach(a=> a.addEventListener('click', close));
    document.addEventListener('keydown', e => { if(e.key==='Escape') close(); });
  }

  /* ───── command palette (fuzzy) ───── */
  function fuzzyScore(q, text){
    q = q.toLowerCase(); text = text.toLowerCase();
    if(!q) return 0.001;
    if(text.indexOf(q)>=0) return 2 + (1/(1+text.indexOf(q)));
    let i=0, j=0, gaps=0, lastHit=-1;
    while(i<q.length && j<text.length){
      if(q[i]===text[j]){ if(lastHit>=0) gaps += (j-lastHit-1); lastHit=j; i++; }
      j++;
    }
    if(i<q.length) return 0;
    return 1/(1+gaps);
  }
  function highlight(text, q){
    if(!q) return text;
    const i = text.toLowerCase().indexOf(q.toLowerCase());
    if(i<0) return text;
    return text.slice(0,i)+'<mark>'+text.slice(i,i+q.length)+'</mark>'+text.slice(i+q.length);
  }
  function initPalette(){
    const input = byId('palette-input');
    const box = byId('palette-results');
    const modal = byId('command-palette');
    const payload = (TD.commandIndex || []).map(item => ({ ...item, href: normalizeAppHref(item.href) }));
    if(!input || !box || !modal) return;
    let cursor = 0;
    let rendered = [];
    function render(q){
      const scored = payload.map(x => ({...x, s: fuzzyScore(q, x.label + ' ' + (x.tags||''))}))
                            .filter(x => x.s > 0)
                            .sort((a,b) => b.s - a.s)
                            .slice(0, 40);
      rendered = scored;
      cursor = 0;
      box.innerHTML = scored.map((x,i)=>
        `<a href="${x.href}" data-idx="${i}" class="${i===0?'hl':''}">
           <strong>${highlight(x.label,q)}</strong>
           <div class="subtle">${x.tags||''}</div>
         </a>`
      ).join('') || '<div class="subtle" style="padding:12px">일치하는 항목 없음</div>';
    }
    function moveCursor(delta){
      if(!rendered.length) return;
      cursor = (cursor + delta + rendered.length) % rendered.length;
      $$('#palette-results a').forEach((a,i)=> a.classList.toggle('hl', i===cursor));
      const hl = $('#palette-results a.hl');
      if(hl) hl.scrollIntoView({block:'nearest'});
    }
    input.addEventListener('input', e => render(e.target.value));
    input.addEventListener('keydown', e => {
      if(e.key==='ArrowDown'){ e.preventDefault(); moveCursor(1); }
      if(e.key==='ArrowUp'){ e.preventDefault(); moveCursor(-1); }
      if(e.key==='Enter'){
        e.preventDefault();
        const a = $('#palette-results a.hl') || $('#palette-results a');
        if(a) location.href = a.getAttribute('href');
      }
    });
    document.addEventListener('keydown', e => {
      if((e.metaKey || e.ctrlKey) && e.key.toLowerCase()==='k'){
        e.preventDefault();
        modal.classList.add('open');
        input.focus(); input.select();
      }
      if(e.key==='Escape' && modal.classList.contains('open')){ modal.classList.remove('open'); }
      if(e.key.toLowerCase()==='p' && !['INPUT','TEXTAREA','SELECT'].includes(document.activeElement.tagName)){
        document.body.classList.toggle('presenter');
      }
    });
    modal.addEventListener('click', e => { if(e.target===modal) modal.classList.remove('open'); });
    render('');
  }

  /* ───── figure modal (enhanced) ───── */
  function initFigureModal(){
    const modal = byId('figure-modal');
    if(!modal) return;
    const frame = byId('figure-modal-frame');
    const label = byId('figure-modal-title');
    const openBtn = byId('figure-modal-open');
    const copyBtn = byId('figure-modal-copy');
    const pngBtn = byId('figure-modal-png');
    const tsvBtn = byId('figure-modal-tsv');
    let lastFocused = null;
    const focusable = ()=> $$('#figure-modal [data-figure-close],#figure-modal a,#figure-modal button').filter(el=>!el.hasAttribute('disabled'));

    function open(url, title, meta){
      lastFocused = document.activeElement;
      frame.src = url;
      if(openBtn) openBtn.href = url;
      if(label) label.textContent = title || '인터랙티브 Figure';
      const base = url.replace(/\.html$/, '');
      if(pngBtn) pngBtn.href = (meta && meta.png) || base.replace('/figs_interactive/','/figs_preview/') + '.png';
      if(tsvBtn) tsvBtn.href = (meta && meta.tsv) || base + '.tsv';
      if(copyBtn) copyBtn.onclick = ()=>{ TD.copyText(location.origin + location.pathname + '#fig=' + encodeURIComponent(title||'')); TD.toast('링크가 복사되었습니다'); };
      modal.classList.remove('hidden');
      modal.setAttribute('aria-hidden','false');
      document.body.style.overflow='hidden';
      setTimeout(()=>{ const fs = focusable(); fs[0] && fs[0].focus(); }, 30);
    }
    function close(){
      modal.classList.add('hidden');
      modal.setAttribute('aria-hidden','true');
      frame.src='about:blank';
      document.body.style.overflow='';
      lastFocused && lastFocused.focus && lastFocused.focus();
    }
    $$('[data-figure-open]').forEach(el => {
      el.addEventListener('click', e => {
        e.preventDefault();
        open(el.getAttribute('data-figure-open'),
             el.getAttribute('data-figure-title') || '인터랙티브 Figure',
             { png: el.getAttribute('data-figure-png'), tsv: el.getAttribute('data-figure-tsv') });
      });
    });
    $$('[data-figure-close]').forEach(el => el.addEventListener('click', close));
    document.addEventListener('keydown', e => {
      if(modal.classList.contains('hidden')) return;
      if(e.key==='Escape'){ close(); return; }
      if(e.key==='Tab'){
        const fs = focusable(); if(!fs.length) return;
        const first=fs[0], last=fs[fs.length-1];
        if(e.shiftKey && document.activeElement===first){ e.preventDefault(); last.focus(); }
        else if(!e.shiftKey && document.activeElement===last){ e.preventDefault(); first.focus(); }
      }
    });
    const fsBtn = byId('figure-modal-fullscreen');
    if(fsBtn){
      fsBtn.addEventListener('click', ()=> {
        const card = $('.figure-modal-card');
        if(!document.fullscreenElement){ card && card.requestFullscreen && card.requestFullscreen(); }
        else { document.exitFullscreen(); }
      });
    }
  }

  /* ───── share-view / copy-citation ───── */
  function initShareAndCopy(){
    $$('[data-copy-link]').forEach(el => el.addEventListener('click', async e => {
      e.preventDefault();
      await TD.copyText(location.href);
      TD.toast('URL이 클립보드에 복사되었습니다');
    }));
    $$('[data-copy-cite]').forEach(el => el.addEventListener('click', async e => {
      e.preventDefault();
      const cite = el.getAttribute('data-copy-cite') || `THCA Dashboard — ${document.title} — ${location.href}`;
      await TD.copyText(cite);
      TD.toast('인용 정보가 복사되었습니다');
    }));
    $$('[data-copy-code]').forEach(btn => btn.addEventListener('click', async () => {
      const pre = btn.closest('pre');
      if(!pre) return;
      const code = pre.innerText.replace(/Copy\s*$/,'').trim();
      await TD.copyText(code);
      btn.innerText = 'Copied'; setTimeout(()=>btn.innerText='Copy',1200);
    }));
  }

  /* ───── count-up (SVG progress + numeric) ───── */
  function initCountUp(){
    $$('[data-countup]').forEach(el => {
      const end = Number(el.getAttribute('data-countup'));
      if(Number.isNaN(end)) return;
      if(prefersReduced){ el.textContent = end.toLocaleString(); return; }
      const dur = 900, start = performance.now();
      function step(t){
        const p = Math.min(1, (t-start)/dur);
        const e = 1 - Math.pow(1-p, 3);
        el.textContent = Math.round(end * e).toLocaleString();
        if(p<1) requestAnimationFrame(step);
      }
      requestAnimationFrame(step);
    });
  }

  /* ───── AOS-like fade-in on scroll ───── */
  function initReveal(){
    if(prefersReduced) return;
    if(!('IntersectionObserver' in window)) return;
    const obs = new IntersectionObserver((entries)=>{
      entries.forEach(en => {
        if(en.isIntersecting){
          en.target.style.opacity='1';
          en.target.style.transform='translateY(0)';
          obs.unobserve(en.target);
        }
      });
    }, { threshold: 0.04, rootMargin: '0px 0px -60px 0px' });
    $$('[data-aos]').forEach(el=>{
      const r = el.getBoundingClientRect();
      // If the element is already in viewport on first paint, don't hide it — just leave as-is.
      if(r.top < (window.innerHeight || 800) && r.bottom > 0){
        el.style.opacity='1'; el.style.transform='none'; return;
      }
      el.style.opacity='0';
      el.style.transform='translateY(14px)';
      el.style.transition='opacity .5s cubic-bezier(.22,1,.36,1), transform .5s cubic-bezier(.22,1,.36,1)';
      obs.observe(el);
    });
  }

  /* ───── lazy-load JSON payloads ───── */
  TD.lazyLoadJSON = function(target, url, onLoad){
    const el = typeof target === 'string' ? byId(target) : target;
    if(!el){ return; }
    el.classList.add('skeleton');
    const shell = el.closest('.figure-preview-shell') || el;
    shell.classList.add('skeleton');
    const run = ()=>{
      fetch(url)
        .then(r => r.json())
        .then(data => {
          shell.classList.remove('skeleton');
          el.classList.remove('skeleton');
          try { onLoad && onLoad(data, el); } catch(err){ console.error(err); el.innerHTML='<p class="subtle">데이터 렌더에 실패했습니다.</p>'; }
        })
        .catch(err => {
          shell.classList.remove('skeleton');
          el.innerHTML = '<p class="subtle">데이터 로드 실패: '+(err && err.message || err)+'</p>';
        });
    };
    if(!('IntersectionObserver' in window)){ run(); return; }
    const obs = new IntersectionObserver((entries, o)=>{
      entries.forEach(en => {
        if(en.isIntersecting){ o.unobserve(en.target); run(); }
      });
    }, { rootMargin: '120px' });
    obs.observe(el);
  };

  /* ───── sortable + searchable tables (.dx-table) ───── */
  function initDxTables(){
    $$('.dx-table-wrap').forEach(wrap => {
      const table = $('.dx-table', wrap);
      const search = $('.dx-table-search', wrap);
      const status = $('.dx-table-status', wrap);
      if(!table) return;
      const headers = $$('thead th', table);
      const tbody = $('tbody', table);
      const rows = $$('tr', tbody);
      const cache = rows.map(r => ({ el: r, cells: $$('td', r).map(td => td.textContent) }));
      let sortCol = -1, sortDir = 1;
      function applySort(){
        if(sortCol<0) return;
        cache.sort((a,b)=>{
          const va = a.cells[sortCol] ?? '', vb = b.cells[sortCol] ?? '';
          const fa = parseFloat(va), fb = parseFloat(vb);
          if(!isNaN(fa) && !isNaN(fb) && String(fa)===va.trim() && String(fb)===vb.trim()) return (fa-fb)*sortDir;
          return va.localeCompare(vb, 'ko')*sortDir;
        });
        cache.forEach(c => tbody.appendChild(c.el));
      }
      headers.forEach((th, i) => {
        th.setAttribute('role','button'); th.tabIndex=0;
        const click = ()=>{
          if(sortCol===i) sortDir *= -1; else { sortCol=i; sortDir=1; }
          headers.forEach((h,j)=> h.setAttribute('aria-sort', j===sortCol ? (sortDir===1?'ascending':'descending') : 'none'));
          applySort();
        };
        th.addEventListener('click', click);
        th.addEventListener('keydown', e => { if(e.key==='Enter' || e.key===' '){ e.preventDefault(); click(); } });
      });
      function applyFilter(){
        const q = (search && search.value || '').trim().toLowerCase();
        let shown = 0;
        cache.forEach(({el, cells}) => {
          const hay = cells.join(' ').toLowerCase();
          const match = !q || hay.indexOf(q) >= 0;
          el.style.display = match ? '' : 'none';
          if(match){
            shown++;
            $$('td', el).forEach((td,i)=>{
              const raw = cells[i];
              if(q && raw.toLowerCase().indexOf(q)>=0){
                const re = new RegExp('('+q.replace(/[.*+?^${}()|[\]\\]/g,'\\$&')+')','gi');
                td.innerHTML = raw.replace(re,'<mark>$1</mark>');
              } else {
                td.textContent = raw;
              }
            });
          }
        });
        if(status){ status.textContent = shown + ' / ' + cache.length + ' 행'; }
      }
      if(search){ search.addEventListener('input', applyFilter); }
      applyFilter();
    });
  }

  /* ───── filter bar ───── */
  TD.initFilterBar = function(opts){
    const bar = byId(opts.barId);
    const target = byId(opts.tableId);
    if(!bar || !target) return;
    const rows = $$('tbody tr', target);
    const fields = opts.fields || [];
    const status = $('.filter-status', bar) || null;
    function apply(){
      let shown = 0;
      rows.forEach(tr => {
        let ok = true;
        fields.forEach(f => {
          const sel = bar.querySelector('[data-filter="'+f.id+'"]');
          if(!sel) return;
          const picked = Array.from(sel.selectedOptions || []).map(o=>o.value).filter(Boolean);
          if(picked.length && !picked.includes(tr.dataset[f.id])) ok = false;
        });
        tr.style.display = ok ? '' : 'none';
        if(ok) shown++;
      });
      if(status){ status.textContent = shown + ' / ' + rows.length + '건 표시'; }
    }
    fields.forEach(f => {
      const sel = bar.querySelector('[data-filter="'+f.id+'"]');
      if(sel){ sel.addEventListener('change', apply); sel.addEventListener('input', apply); }
    });
    apply();
  };

  /* ───── tabs ───── */
  function initTabs(){
    $$('[data-tabs]').forEach(group => {
      const buttons = $$('[role=tab]', group);
      const panels = $$('[role=tabpanel]', group);
      function activate(i){
        buttons.forEach((b,j)=> b.setAttribute('aria-selected', i===j));
        panels.forEach((p,j)=> p.classList.toggle('active', i===j));
      }
      buttons.forEach((b,i)=>{
        b.addEventListener('click', ()=>activate(i));
        b.addEventListener('keydown', e => {
          if(e.key==='ArrowRight'){ e.preventDefault(); activate((i+1)%buttons.length); buttons[(i+1)%buttons.length].focus(); }
          if(e.key==='ArrowLeft'){ e.preventDefault(); activate((i-1+buttons.length)%buttons.length); buttons[(i-1+buttons.length)%buttons.length].focus(); }
        });
      });
      activate(0);
    });
  }

  /* ───── prefetch on hover ───── */
  function initPrefetch(){
    const done = new Set();
    function queue(url){
      if(!url || done.has(url)) return;
      done.add(url);
      const l = document.createElement('link');
      l.rel = 'prefetch'; l.href = url; l.as = 'document';
      document.head.appendChild(l);
    }
    $$('.topnav a.nav-link, #td-drawer a').forEach(a => {
      let t;
      a.addEventListener('mouseenter', ()=>{ t = setTimeout(()=>queue(a.href), 80); });
      a.addEventListener('mouseleave', ()=>{ clearTimeout(t); });
      a.addEventListener('focus', ()=> queue(a.href));
    });
  }

  /* ───── markdown ───── */
  function initMarkdownRender(){
    $$('[data-md-source]').forEach(async el => {
      try {
        const src = el.getAttribute('data-md-source');
        const txt = await fetch(src).then(r=>r.text());
        if(window.marked){ el.innerHTML = window.marked.parse(txt); }
        else { el.innerHTML = '<pre>'+txt.replace(/[<&>]/g, ch=>({ '<':'&lt;','&':'&amp;','>':'&gt;' }[ch]))+'</pre>'; }
        $$('pre', el).forEach(pre => {
          if(!pre.querySelector('.copy-btn')){
            const b = document.createElement('button');
            b.className='copy-btn'; b.textContent='Copy'; b.setAttribute('data-copy-code','1');
            pre.appendChild(b);
          }
        });
        initShareAndCopy();
      } catch(e){
        el.innerHTML = '<p class="subtle">Markdown 렌더 실패.</p>';
      }
    });
  }

  /* ───── grid table fallback (used by existing page_js) ───── */
  window.renderGridTable = function(targetId, rows, columns, searchId){
    const target = byId(targetId); if(!target) return;
    target.innerHTML = '';
    const wrap = document.createElement('div'); wrap.className = 'dx-table-wrap';
    const status = document.createElement('div'); status.className = 'dx-table-status';
    const table = document.createElement('table'); table.className = 'dx-table';
    const thead = document.createElement('thead');
    const trh = document.createElement('tr');
    columns.forEach(c => { const th = document.createElement('th'); th.textContent = c.name || c.id; th.setAttribute('aria-sort','none'); trh.appendChild(th); });
    thead.appendChild(trh); table.appendChild(thead);
    const tbody = document.createElement('tbody');
    const limit = 300;
    rows.slice(0, limit).forEach(r => {
      const tr = document.createElement('tr');
      columns.forEach(c => { const td = document.createElement('td'); const v = r[c.id]; td.textContent = (v===undefined||v===null||Number.isNaN(v))?'':String(v); tr.appendChild(td); });
      tbody.appendChild(tr);
    });
    table.appendChild(tbody);
    wrap.appendChild(table); wrap.appendChild(status);
    target.appendChild(wrap);
    // Wire search input to wrap (if an external one exists)
    if(searchId){
      const ext = byId(searchId);
      if(ext){
        ext.classList.add('dx-table-search');
        wrap.insertBefore(ext.cloneNode(true), wrap.firstChild);
        // Sync: user types in original, proxy into internal
        const internal = wrap.querySelector('.dx-table-search');
        internal.id = searchId + '_proxy';
        internal.addEventListener('input', () => { ext.value = internal.value; ext.dispatchEvent(new Event('input')); });
        ext.addEventListener('input', () => { if(internal.value !== ext.value) internal.value = ext.value; });
      }
    }
    initDxTables();
    if(rows.length > limit && status){
      status.textContent = '처음 '+limit+'행을 표시 (전체 '+rows.length+')';
    }
  };

  /* ───── init ───── */
  document.addEventListener('DOMContentLoaded', () => {
    initTheme();
    initPathFixes();
    initDrawer();
    initPalette();
    initFigureModal();
    initShareAndCopy();
    initCountUp();
    initReveal();
    initDxTables();
    initTabs();
    initPrefetch();
    initMarkdownRender();
    ensureSiteMap();
  });
})();
