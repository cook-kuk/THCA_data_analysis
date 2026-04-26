/*
 * table-export.js — CSV/TSV export buttons for every table on every page.
 *
 * Scans on DOMContentLoaded + MutationObserver to cover async-rendered tables
 * (e.g. renderGridTable, Gene Explorer, view_researcher dynamic builds).
 *
 * Respects .dx-table filtering: rows with `display:none` are skipped.
 * Filename: `{page-slug}_{section-id}.csv`; UTF-8 with BOM for Excel-Korean.
 * No external dependencies, CSP-safe.
 */
(function () {
  'use strict';
  if (window.__thyroidTableExportLoaded) return;
  window.__thyroidTableExportLoaded = true;

  var BTN_MARK = 'data-csv-export-btn';
  var WIRED = 'data-csv-export-wired';
  var BOM = '﻿';

  function pageSlug() {
    var name = (location.pathname.split('/').pop() || 'page').replace(/\.html?$/i, '');
    if (!name || name === '') name = 'index';
    return name.replace(/[^A-Za-z0-9_-]/g, '_');
  }

  function nearestSectionId(table) {
    var el = table;
    while (el && el !== document.body) {
      if (el.id) return el.id;
      el = el.parentElement;
    }
    // fall back: walk up to <section>
    el = table.closest && table.closest('section,article,div[id]');
    return (el && el.id) ? el.id : 'table';
  }

  function cellText(td) {
    // Clone to avoid touching inline <mark> highlight state.
    var clone = td.cloneNode(true);
    // strip any export buttons accidentally included
    clone.querySelectorAll && clone.querySelectorAll('[' + BTN_MARK + ']').forEach(function (n) { n.remove(); });
    var txt = (clone.textContent || '').replace(/\s+/g, ' ').trim();
    return txt;
  }

  function csvEscape(v) {
    if (v == null) return '';
    var s = String(v);
    if (/[",\n\r]/.test(s)) return '"' + s.replace(/"/g, '""') + '"';
    return s;
  }

  function serialize(table) {
    var rows = [];
    var theadRows = table.tHead ? Array.prototype.slice.call(table.tHead.rows) : [];
    var bodyRows = [];
    var tbodies = table.tBodies ? Array.prototype.slice.call(table.tBodies) : [];
    tbodies.forEach(function (tb) {
      Array.prototype.slice.call(tb.rows).forEach(function (r) {
        var style = window.getComputedStyle(r);
        if (style.display === 'none') return;
        bodyRows.push(r);
      });
    });
    // If the table has no explicit thead, first row may live in tbody; keep as-is.
    theadRows.forEach(function (r) {
      var cells = Array.prototype.slice.call(r.cells).map(cellText);
      rows.push(cells.map(csvEscape).join(','));
    });
    bodyRows.forEach(function (r) {
      var cells = Array.prototype.slice.call(r.cells).map(cellText);
      rows.push(cells.map(csvEscape).join(','));
    });
    return rows.join('\r\n');
  }

  function triggerDownload(filename, csvText) {
    try {
      var blob = new Blob([BOM + csvText], { type: 'text/csv;charset=utf-8' });
      var url = URL.createObjectURL(blob);
      var a = document.createElement('a');
      a.href = url;
      a.download = filename;
      a.style.display = 'none';
      document.body.appendChild(a);
      a.click();
      setTimeout(function () {
        try { document.body.removeChild(a); } catch (e) {}
        URL.revokeObjectURL(url);
      }, 120);
    } catch (e) {
      console.warn('[table-export] download failed', e);
    }
  }

  function buildButton(table) {
    var btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn sm ghost';
    btn.setAttribute(BTN_MARK, '1');
    btn.setAttribute('aria-label', 'CSV 다운로드 (table CSV download)');
    btn.textContent = 'CSV 다운로드';
    btn.style.margin = '8px 0';
    btn.addEventListener('click', function (ev) {
      ev.preventDefault();
      var csv = serialize(table);
      var secId = nearestSectionId(table);
      var fname = pageSlug() + '_' + (secId || 'table') + '.csv';
      triggerDownload(fname, csv);
    });
    return btn;
  }

  function attach(table) {
    if (!table || table.getAttribute(WIRED) === '1') return;
    if (!table.tBodies || table.tBodies.length === 0) return;
    // Skip tables embedded inside already-exported content (paranoia).
    if (table.closest && table.closest('[' + BTN_MARK + ']')) return;
    var wrap = table.parentElement;
    if (!wrap) return;
    // If wrap already contains an export button for this table, skip.
    var existing = wrap.querySelector(':scope > [' + BTN_MARK + ']');
    if (existing) { table.setAttribute(WIRED, '1'); return; }
    var btn = buildButton(table);
    // Insert just before the table so it reads "download above the data".
    wrap.insertBefore(btn, table);
    table.setAttribute(WIRED, '1');
  }

  function scan(root) {
    var scope = root && root.querySelectorAll ? root : document;
    var tables = scope.querySelectorAll('table');
    tables.forEach(attach);
  }

  function init() {
    scan(document);
    // Watch for async-rendered tables (renderGridTable, gene-explorer, etc.).
    try {
      var mo = new MutationObserver(function (muts) {
        for (var i = 0; i < muts.length; i++) {
          var m = muts[i];
          for (var j = 0; j < m.addedNodes.length; j++) {
            var node = m.addedNodes[j];
            if (node.nodeType !== 1) continue;
            if (node.tagName === 'TABLE') attach(node);
            else if (node.querySelectorAll) scan(node);
          }
        }
      });
      mo.observe(document.body, { childList: true, subtree: true });
    } catch (e) { /* older browser */ }
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init, { once: true });
  } else {
    init();
  }

  // Expose a tiny hook for tests.
  window.ThyroidDash = window.ThyroidDash || {};
  window.ThyroidDash.tableExport = { scan: scan, attach: attach };
})();
