/*
 * tour.js — guided onboarding tour with SVG spotlight cutout.
 *
 * Auto-starts on first visit (localStorage flag `thyroid-tour-seen`). Can be
 * re-triggered via the `? 투어` button injected into `.nav-actions`.
 * Pure vanilla; CSP-safe (no eval, no inline handlers attached via string).
 */
(function () {
  'use strict';
  if (window.__thyroidTourLoaded) return;
  window.__thyroidTourLoaded = true;

  var STORAGE_KEY = 'thyroid-tour-seen';
  var SVG_NS = 'http://www.w3.org/2000/svg';

  // Page-scope can be: '*' (all pages) OR a substring of pathname.
  var STEPS = [
    {
      selector: '.topnav .nav-links',
      title: '상단 내비게이션 (Topnav)',
      body: '홈, 투자자, 연구자 뷰 간 빠른 전환이 가능합니다. 모든 섹션은 여기서 접근하세요.',
      scope: '*'
    },
    {
      selector: '.nav-actions',
      title: '명령 팔레트 (Command palette)',
      body: 'Cmd/Ctrl+K로 섹션과 유전자를 검색할 수 있습니다. "링크 복사" 버튼으로 현재 뷰를 공유하세요.',
      scope: '*'
    },
    {
      selector: '#threshold-explorer-section',
      title: '임계치 탐색기 (Interactive threshold explorer)',
      body: '슬라이더로 decision threshold를 조정하며 precision/recall/F1 변화를 실시간 확인합니다.',
      scope: '07_ml_baseline.html'
    },
    {
      selector: 'details.edu-toggle',
      title: '교육용 토글 (Educational toggles)',
      body: '페이지 곳곳의 💡 아이콘으로 시작하는 설명 박스를 열어 용어와 해석을 확인하세요.',
      scope: '*'
    },
    {
      selector: '.topnav .nav-link[href*="view_investor"], .topnav .nav-link[href*="view_researcher"]',
      title: '페르소나 랜딩 (Persona landings)',
      body: '투자자/연구자 뷰에서 각 페르소나에 맞춰 큐레이션된 엔트리 카드를 확인할 수 있습니다.',
      scope: '*'
    },
    {
      selector: '.topnav .nav-link[href*="99_glossary"]',
      title: '용어집 (Glossary)',
      body: '모든 기술 용어 정의는 용어집에서 검색할 수 있습니다. 약어 · 메트릭 · 데이터셋명 포함.',
      scope: '*'
    }
  ];

  function basename() {
    return (location.pathname.split('/').pop() || '').toLowerCase();
  }

  function filteredSteps() {
    var bn = basename();
    return STEPS.filter(function (s) {
      if (s.scope === '*') return true;
      return bn.indexOf(s.scope.toLowerCase()) >= 0;
    });
  }

  function prefersReducedMotion() {
    try { return window.matchMedia('(prefers-reduced-motion: reduce)').matches; }
    catch (e) { return false; }
  }

  var overlay = null;
  var tooltip = null;
  var currentIdx = 0;
  var active = false;
  var steps = [];

  function removeOverlay() {
    if (overlay && overlay.parentNode) overlay.parentNode.removeChild(overlay);
    if (tooltip && tooltip.parentNode) tooltip.parentNode.removeChild(tooltip);
    overlay = null;
    tooltip = null;
    active = false;
    document.removeEventListener('keydown', onKey);
  }

  function finish(flagSeen) {
    if (flagSeen) {
      try { localStorage.setItem(STORAGE_KEY, '1'); } catch (e) {}
    }
    removeOverlay();
  }

  function onKey(e) {
    if (!active) return;
    if (e.key === 'Escape') { e.preventDefault(); finish(true); }
    else if (e.key === 'ArrowRight') { e.preventDefault(); next(); }
    else if (e.key === 'ArrowLeft') { e.preventDefault(); prev(); }
  }

  function pickTarget(sel) {
    try {
      var candidates = document.querySelectorAll(sel);
      for (var i = 0; i < candidates.length; i++) {
        var el = candidates[i];
        var r = el.getBoundingClientRect();
        if (r.width > 0 && r.height > 0) return el;
      }
    } catch (e) {}
    return null;
  }

  function ensureOverlay() {
    if (overlay) return overlay;
    overlay = document.createElement('div');
    overlay.className = 'thyroid-tour-overlay';
    overlay.setAttribute('role', 'presentation');
    overlay.setAttribute('data-print-hide', '1');
    var svg = document.createElementNS(SVG_NS, 'svg');
    svg.setAttribute('width', '100%');
    svg.setAttribute('height', '100%');
    svg.setAttribute('preserveAspectRatio', 'none');
    var defs = document.createElementNS(SVG_NS, 'defs');
    var mask = document.createElementNS(SVG_NS, 'mask');
    mask.setAttribute('id', 'thyroidTourMask');
    var full = document.createElementNS(SVG_NS, 'rect');
    full.setAttribute('x', '0'); full.setAttribute('y', '0');
    full.setAttribute('width', '100%'); full.setAttribute('height', '100%');
    full.setAttribute('fill', 'white');
    mask.appendChild(full);
    var cut = document.createElementNS(SVG_NS, 'rect');
    cut.setAttribute('class', 'tour-cutout');
    cut.setAttribute('fill', 'black');
    cut.setAttribute('rx', '8'); cut.setAttribute('ry', '8');
    cut.setAttribute('x', '-999'); cut.setAttribute('y', '-999');
    cut.setAttribute('width', '1'); cut.setAttribute('height', '1');
    mask.appendChild(cut);
    defs.appendChild(mask);
    svg.appendChild(defs);
    var back = document.createElementNS(SVG_NS, 'rect');
    back.setAttribute('class', 'tour-backdrop');
    back.setAttribute('x', '0'); back.setAttribute('y', '0');
    back.setAttribute('width', '100%'); back.setAttribute('height', '100%');
    back.setAttribute('mask', 'url(#thyroidTourMask)');
    svg.appendChild(back);
    overlay.appendChild(svg);
    // Clicking backdrop dismisses
    overlay.addEventListener('click', function (ev) {
      if (ev.target === overlay || ev.target.tagName === 'svg' || ev.target === back) {
        finish(true);
      }
    });
    document.body.appendChild(overlay);
    overlay.__cut = cut;

    tooltip = document.createElement('div');
    tooltip.className = 'thyroid-tour-tooltip';
    tooltip.setAttribute('role', 'dialog');
    tooltip.setAttribute('aria-live', 'polite');
    tooltip.setAttribute('data-print-hide', '1');
    document.body.appendChild(tooltip);
    return overlay;
  }

  function positionTooltip(targetRect) {
    if (!tooltip) return;
    var pad = 12;
    var tipW = tooltip.offsetWidth || 300;
    var tipH = tooltip.offsetHeight || 160;
    var vw = window.innerWidth, vh = window.innerHeight;
    var left, top;
    // Prefer below the target, else above, else centered.
    if (!targetRect || (targetRect.width === 0 && targetRect.height === 0)) {
      left = Math.max(16, (vw - tipW) / 2);
      top = Math.max(16, (vh - tipH) / 2);
    } else if (targetRect.bottom + pad + tipH < vh) {
      top = targetRect.bottom + pad;
      left = Math.min(vw - tipW - 16, Math.max(16, targetRect.left));
    } else if (targetRect.top - pad - tipH > 0) {
      top = targetRect.top - pad - tipH;
      left = Math.min(vw - tipW - 16, Math.max(16, targetRect.left));
    } else {
      // Side-dock right if there's room, else left
      top = Math.max(16, Math.min(vh - tipH - 16, targetRect.top));
      if (targetRect.right + pad + tipW < vw) left = targetRect.right + pad;
      else left = Math.max(16, targetRect.left - tipW - pad);
    }
    tooltip.style.left = left + 'px';
    tooltip.style.top = top + 'px';
  }

  function render() {
    if (currentIdx < 0 || currentIdx >= steps.length) { finish(true); return; }
    ensureOverlay();
    var step = steps[currentIdx];
    var target = pickTarget(step.selector);
    var cut = overlay.__cut;
    if (target) {
      var r = target.getBoundingClientRect();
      var pad = 8;
      cut.setAttribute('x', Math.max(0, r.left - pad));
      cut.setAttribute('y', Math.max(0, r.top - pad));
      cut.setAttribute('width', Math.max(1, r.width + pad * 2));
      cut.setAttribute('height', Math.max(1, r.height + pad * 2));
      try { target.scrollIntoView({ behavior: prefersReducedMotion() ? 'auto' : 'smooth', block: 'center' }); } catch (e) {}
    } else {
      cut.setAttribute('x', '-999');
      cut.setAttribute('y', '-999');
      cut.setAttribute('width', '1');
      cut.setAttribute('height', '1');
    }
    // Tooltip content
    tooltip.innerHTML = '';
    var meta = document.createElement('div');
    meta.className = 'tour-step-count';
    meta.textContent = (currentIdx + 1) + ' / ' + steps.length;
    var h = document.createElement('h3');
    h.textContent = step.title;
    var p = document.createElement('p');
    p.textContent = step.body;
    var actions = document.createElement('div');
    actions.className = 'tour-actions';

    var skipBtn = document.createElement('button');
    skipBtn.type = 'button';
    skipBtn.textContent = '건너뛰기';
    skipBtn.addEventListener('click', function () { finish(true); });

    var prevBtn = document.createElement('button');
    prevBtn.type = 'button';
    prevBtn.textContent = '이전';
    prevBtn.disabled = currentIdx === 0;
    if (prevBtn.disabled) prevBtn.style.opacity = '0.4';
    prevBtn.addEventListener('click', prev);

    var nextBtn = document.createElement('button');
    nextBtn.type = 'button';
    nextBtn.className = 'primary';
    nextBtn.textContent = (currentIdx === steps.length - 1) ? '완료' : '다음';
    nextBtn.addEventListener('click', function () {
      if (currentIdx === steps.length - 1) finish(true);
      else next();
    });

    actions.appendChild(skipBtn);
    actions.appendChild(prevBtn);
    actions.appendChild(nextBtn);

    tooltip.appendChild(meta);
    tooltip.appendChild(h);
    tooltip.appendChild(p);
    tooltip.appendChild(actions);

    // Wait one frame for layout before positioning.
    requestAnimationFrame(function () {
      positionTooltip(target ? target.getBoundingClientRect() : null);
    });
  }

  function next() { currentIdx++; render(); }
  function prev() { currentIdx--; if (currentIdx < 0) currentIdx = 0; render(); }

  function start() {
    steps = filteredSteps();
    if (!steps.length) return;
    currentIdx = 0;
    active = true;
    document.addEventListener('keydown', onKey);
    render();
    window.addEventListener('resize', onResize, { passive: true });
    window.addEventListener('scroll', onResize, { passive: true });
  }

  function onResize() {
    if (!active || currentIdx < 0 || currentIdx >= steps.length) return;
    var step = steps[currentIdx];
    var target = pickTarget(step.selector);
    var cut = overlay && overlay.__cut;
    if (target && cut) {
      var r = target.getBoundingClientRect();
      var pad = 8;
      cut.setAttribute('x', Math.max(0, r.left - pad));
      cut.setAttribute('y', Math.max(0, r.top - pad));
      cut.setAttribute('width', Math.max(1, r.width + pad * 2));
      cut.setAttribute('height', Math.max(1, r.height + pad * 2));
      positionTooltip(r);
    }
  }

  function injectTriggerButton() {
    var actions = document.querySelector('.nav-actions');
    if (!actions) return;
    if (actions.querySelector('.thyroid-tour-btn')) return;
    var btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn sm ghost thyroid-tour-btn';
    btn.textContent = '? 투어';
    btn.setAttribute('aria-label', '가이드 투어 시작 (Start guided tour)');
    btn.addEventListener('click', function () { start(); });
    // Insert as the first child so it sits near theme toggle.
    actions.insertBefore(btn, actions.firstChild);
  }

  function init() {
    injectTriggerButton();
    var seen = false;
    try { seen = localStorage.getItem(STORAGE_KEY) === '1'; } catch (e) {}
    if (!seen) {
      // delay slightly so page has laid out
      setTimeout(start, 700);
    }
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init, { once: true });
  } else {
    init();
  }

  window.ThyroidDash = window.ThyroidDash || {};
  window.ThyroidDash.tour = { start: start, reset: function () { try { localStorage.removeItem(STORAGE_KEY); } catch (e) {} } };
  window.ThyroidDash.tourAvailable = true;
})();
