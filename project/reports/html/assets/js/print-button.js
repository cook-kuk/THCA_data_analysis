/*
 * print-button.js — floating "인쇄 / PDF 저장" button on report-style pages.
 *
 * The allowlist lives inline (no extra fetch). Button is only rendered on
 * pages whose basename matches. Idempotent: re-injection is a no-op.
 */
(function () {
  'use strict';
  if (window.__thyroidPrintBtnLoaded) return;
  window.__thyroidPrintBtnLoaded = true;

  var ALLOW = [
    '01_overview.html',
    '07_ml_baseline.html',
    '14_caveats.html',
    '15_drug_discovery.html',
    '16_quantum.html',
    '17_biomarker_insights.html',
    'view_investor.html',
    'view_researcher.html'
  ];

  function basename() {
    return (location.pathname.split('/').pop() || '').toLowerCase();
  }

  function init() {
    if (ALLOW.indexOf(basename()) < 0) return;
    if (document.querySelector('.thyroid-print-btn')) return;

    var btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn sm thyroid-print-btn';
    btn.textContent = '인쇄 / PDF 저장';
    btn.setAttribute('aria-label', '인쇄하거나 PDF로 저장 (Print / Save as PDF)');
    btn.setAttribute('data-print-hide', '1');
    btn.style.cssText = [
      'position:fixed',
      'right:18px',
      'bottom:18px',
      'z-index:9998',
      'padding:10px 14px',
      'border-radius:999px',
      'border:1px solid rgba(20,184,166,.6)',
      'background:rgba(13,148,136,0.92)',
      'color:#f8fafc',
      'font-weight:600',
      'font-size:0.82rem',
      'box-shadow:0 6px 24px rgba(15,23,42,.35)',
      'cursor:pointer',
      'backdrop-filter:blur(8px)'
    ].join(';');
    btn.addEventListener('click', function () {
      try { window.print(); } catch (e) { console.warn('print failed', e); }
    });
    btn.addEventListener('mouseenter', function () {
      btn.style.background = 'rgba(15,118,110,1)';
    });
    btn.addEventListener('mouseleave', function () {
      btn.style.background = 'rgba(13,148,136,0.92)';
    });
    document.body.appendChild(btn);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init, { once: true });
  } else {
    init();
  }

  window.ThyroidDash = window.ThyroidDash || {};
  window.ThyroidDash.printButton = { allowlist: ALLOW.slice() };
})();
