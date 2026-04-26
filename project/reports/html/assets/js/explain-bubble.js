/* Explain-bubble — global "자세히" button + modal system.
   Scans for elements with data-explain-id attributes, auto-injects a button
   that opens a modal with detailed explanation loaded from explain_content.json.
   Requires: assets/data/explain_content.json
   Usage: Add data-explain-id="my-key" to any element. Ensure matching entry exists in JSON.
   CSP-safe, no external deps. */
(function(){
  'use strict';
  var CONTENT_URL = '../assets/data/explain_content.json';
  // Try alternate paths for index.html (one level up):
  if (location.pathname.endsWith('/index.html') || location.pathname === '/reports/html/' || location.pathname === '/') {
    CONTENT_URL = 'assets/data/explain_content.json';
  }

  var content = {};
  var modalEl = null;

  function textify(node){
    if (!node) return '';
    if (node.cloneNode) {
      var clone = node.cloneNode(true);
      if (clone.querySelectorAll) {
        clone.querySelectorAll('.explain-btn').forEach(function(btn){ btn.remove(); });
      }
      return (clone.textContent || '').replace(/\s+/g, ' ').trim();
    }
    return (node.textContent || '').replace(/\s+/g, ' ').trim();
  }

  function slugify(text){
    return String(text || '')
      .toLowerCase()
      .replace(/[^a-z0-9가-힣]+/g, '-')
      .replace(/^-+|-+$/g, '')
      .slice(0, 80) || 'section';
  }

  function isHeading(el){
    return !!el && /^H[1-4]$/.test(el.tagName);
  }

  function getPrimaryHeading(host){
    if (!host || !host.querySelector) return null;
    if (isHeading(host)) return host;
    return host.querySelector(':scope > h1, :scope > h2, :scope > h3, :scope > h4, :scope > .section-title, :scope > figcaption, :scope > summary');
  }

  function hasNestedExplainTarget(host){
    if (!host || !host.querySelector) return false;
    return !!host.querySelector('[data-explain-id], [data-explain-auto]');
  }

  function collectSiblingText(host, selector, limit){
    var items = [];
    if (!host || !host.querySelectorAll) return items;
    host.querySelectorAll(selector).forEach(function(el){
      var txt = textify(el);
      if (txt) items.push(txt);
    });
    return items.slice(0, limit || 6);
  }

  function deriveEntryFromDom(id, host){
    var container = host;
    if (isHeading(host)) {
      container = host.closest('section, article, figure, .panel, .card, .glass, .dial-card, .v4a-section, .dd-panel, .cmap-card, .pim-card, .target-card, .biz-card, .rg-card') || host.parentElement || host;
    }
    var heading = isHeading(host) ? host : (getPrimaryHeading(host) || host);
    var title = textify(heading) || id;
    var paragraphs = collectSiblingText(container, 'p', 6);
    var bullets = collectSiblingText(container, 'li', 8);
    var tableHeads = collectSiblingText(container, 'th', 8);
    var controls = collectSiblingText(container, 'label, button, .chip, .section-label, .sec-anchor', 12)
      .filter(function(t){
        return t &&
          t !== '자세히' &&
          t !== 'Close' &&
          t.indexOf('자세히') === -1;
      });
    var tags = [];
    if (container.id) tags.push('section:' + container.id);
    if (container.className && typeof container.className === 'string') {
      container.className.split(/\s+/).filter(Boolean).slice(0, 4).forEach(function(c){ tags.push(c); });
    }

    var sections = [];
    if (paragraphs.length) {
      sections.push({
        heading: 'What This Section Shows',
        body: paragraphs.map(function(t){ return '- ' + t; }).join('\n')
      });
    }
    if (bullets.length) {
      sections.push({
        heading: 'Key Items On Screen',
        body: bullets.map(function(t){ return '- ' + t; }).join('\n')
      });
    }
    if (tableHeads.length) {
      sections.push({
        heading: 'How To Read The Table',
        body: '- The table columns visible in this block are: ' + tableHeads.join(' / ') + '.\n- Read left to right: identifier or cohort first, then metric or interpretation columns, then evidence or output.'
      });
    }
    if (controls.length) {
      sections.push({
        heading: 'Interactive Controls On Screen',
        body: controls.map(function(t){ return '- ' + t; }).join('\n')
      });
    }
    if (container.querySelector && (container.querySelector('canvas') || container.querySelector('svg') || container.querySelector('iframe'))) {
      sections.push({
        heading: 'How To Read This Visual',
        body: '- This block contains an interactive visual layer.\n- Start from the section title, then use the visible filters, legends, presets, or action buttons to change the view.\n- If the graphic appears blank, the underlying JS asset or figure payload may still be loading.'
      });
    }
    if (!sections.length) {
      sections.push({
        heading: 'Why This Matters',
        body: '- This block is part of the THYRAI evidence chain.\n- Use the title, nearby labels, and downloadable artefacts together when explaining the result to a reviewer or professor.\n- If this section looks empty, check whether the figure or table assets finished loading.'
      });
    }

    return {
      title: title,
      tagline: 'Auto-generated page explanation from the visible section content.',
      sections: sections,
      tags: tags.slice(0, 6)
    };
  }

  function getEntry(id, host){
    return content[id] || deriveEntryFromDom(id, host);
  }

  function ensureModal(){
    if (modalEl) return modalEl;
    modalEl = document.createElement('div');
    modalEl.className = 'explain-modal-backdrop';
    modalEl.setAttribute('role', 'dialog');
    modalEl.setAttribute('aria-modal', 'true');
    modalEl.setAttribute('aria-hidden', 'true');
    modalEl.innerHTML = '<div class="explain-modal" role="document">' +
      '<button class="explain-modal-close" aria-label="닫기" type="button">×</button>' +
      '<div class="explain-modal-body"></div>' +
      '</div>';
    document.body.appendChild(modalEl);
    // Close handlers
    modalEl.addEventListener('click', function(e){
      if (e.target === modalEl) closeModal();
    });
    modalEl.querySelector('.explain-modal-close').addEventListener('click', closeModal);
    document.addEventListener('keydown', function(e){
      if (e.key === 'Escape' && modalEl.classList.contains('open')) closeModal();
    });
    return modalEl;
  }

  function renderMarkdownLite(text){
    if (!text) return '';
    // Minimal markdown: **bold**, *italic*, `code`, [text](url), \n\n paragraphs
    return String(text)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
      .replace(/\*\*([^\*]+)\*\*/g, '<strong>$1</strong>')
      .replace(/\*([^\*]+)\*/g, '<em>$1</em>')
      .replace(/`([^`]+)`/g, '<code>$1</code>')
      .replace(/\[([^\]]+)\]\(([^)]+)\)/g, '<a href="$2" target="_blank" rel="noopener">$1</a>');
  }

  function openModal(id, host){
    var entry = getEntry(id, host);
    var m = ensureModal();
    var body = m.querySelector('.explain-modal-body');
    var html = '<h3>' + (entry.title || id) + '</h3>';
    if (entry.tagline) html += '<p class="explain-tagline">' + entry.tagline + '</p>';
    // Body: array of sections or plain body text
    if (Array.isArray(entry.sections)){
      entry.sections.forEach(function(sec){
        if (sec.heading) html += '<h4>' + sec.heading + '</h4>';
        if (sec.body) {
          // Paragraphs split on double newline
          sec.body.split(/\n\n+/).forEach(function(para){
            if (para.trim().startsWith('- ')){
              html += '<ul>';
              para.split('\n').forEach(function(line){
                var t = line.replace(/^-\s+/, '').trim();
                if (t) html += '<li>' + renderMarkdownLite(t) + '</li>';
              });
              html += '</ul>';
            } else {
              html += '<p>' + renderMarkdownLite(para) + '</p>';
            }
          });
        }
      });
    } else if (entry.body) {
      entry.body.split(/\n\n+/).forEach(function(para){
        if (para.trim().startsWith('- ')){
          html += '<ul>';
          para.split('\n').forEach(function(line){
            var t = line.replace(/^-\s+/, '').trim();
            if (t) html += '<li>' + renderMarkdownLite(t) + '</li>';
          });
          html += '</ul>';
        } else {
          html += '<p>' + renderMarkdownLite(para) + '</p>';
        }
      });
    }
    // Tags
    if (entry.tags && entry.tags.length){
      html += '<div class="explain-modal-tags">';
      entry.tags.forEach(function(t){ html += '<span class="explain-tag">' + t + '</span>'; });
      html += '</div>';
    }
    body.innerHTML = html;
    m.classList.add('open');
    m.setAttribute('aria-hidden', 'false');
    // Focus close button for keyboard users
    setTimeout(function(){ m.querySelector('.explain-modal-close').focus(); }, 50);
  }

  function closeModal(){
    if (!modalEl) return;
    modalEl.classList.remove('open');
    modalEl.setAttribute('aria-hidden', 'true');
  }

  function shouldSkipExplicitHost(host){
    if (!host || isHeading(host)) return false;
    if (host.getAttribute('data-explain-skip') === '1') return true;
    var primary = getPrimaryHeading(host);
    if (primary && primary.hasAttribute('data-explain-id')) return true;
    if (primary && primary.hasAttribute('data-explain-auto')) return true;
    return false;
  }

  function attachButton(host, id, mode){
    if (!host || !id) return;
    if (host.querySelector(':scope > .explain-btn')) return;
    var btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'explain-btn';
    if (mode === 'inline' || isHeading(host) || host.tagName === 'SUMMARY' || host.tagName === 'FIGCAPTION') {
      btn.classList.add('explain-btn--inline');
    }
    btn.setAttribute('aria-label', '자세히 보기: ' + id);
    btn.innerHTML = '<span class="icon">💡</span><span class="label-text">자세히</span>';
    btn.addEventListener('click', function(e){
      e.stopPropagation(); e.preventDefault();
      openModal(id, host);
    });
    if (btn.classList.contains('explain-btn--inline')) host.appendChild(btn);
    else host.appendChild(btn);
  }

  function injectButtons(){
    var hosts = document.querySelectorAll('[data-explain-id]');
    hosts.forEach(function(host){
      if (shouldSkipExplicitHost(host)) return;
      var id = host.getAttribute('data-explain-id');
      attachButton(host, id, isHeading(host) ? 'inline' : '');
    });

    var autoHosts = Array.prototype.filter.call(
      document.querySelectorAll('h1, h2, h3, h4, figcaption, details > summary'),
      function(host){
        return !host.closest('nav, footer, .topnav, .site-footer, .explain-modal, .explain-modal-backdrop, .toc, .sidebar');
      }
    );
    autoHosts.forEach(function(host){
      if (host.hasAttribute('data-explain-id')) return;
      if (host.hasAttribute('data-explain-auto')) return;
      if (host.closest('.explain-modal')) return;
      var parent = host.parentElement;
      if (parent && parent !== host && parent.hasAttribute('data-explain-id')) {
        // If parent already has explicit explanation, prefer the heading only.
        parent.setAttribute('data-explain-skip', '1');
      }
      host.setAttribute('data-explain-auto', '1');
      attachButton(host, 'auto-' + slugify(textify(host)), 'inline');
    });
  }

  function init(){
    fetch(CONTENT_URL, { cache: 'no-cache' })
      .then(function(r){ return r.ok ? r.json() : {}; })
      .then(function(data){ content = data || {}; injectButtons(); })
      .catch(function(){ content = {}; injectButtons(); /* still inject for "not yet" message */ });
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();

  // Expose for debug / other scripts
  window.ExplainBubble = { open: openModal, close: closeModal, reinjectButtons: injectButtons };
})();
