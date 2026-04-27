// v17 floating comment widget for 유형원 교수님 review
// Adds a floating "💬 의견 남기기" button to any page; saves to localStorage.
(function() {
  'use strict';

  // Section name = current page title (or override via window.SECTION_NAME)
  // For inline buttons we override per-click to use the nearest heading + page
  var pageTitle = document.title || (location.pathname.split('/').pop() || 'unknown');
  var sectionName = window.SECTION_NAME || pageTitle;
  var pageSlug = (location.pathname.split('/').pop() || 'unknown').replace(/\.html?$/, '');

  // Inject styles
  var style = document.createElement('style');
  style.textContent = `
    #yu-comment-fab {
      position: fixed; bottom: 24px; right: 24px; z-index: 99999;
      background: linear-gradient(135deg, #1e40af, #7c3aed);
      color: white; padding: 14px 20px; border-radius: 30px;
      font-size: 14px; font-weight: 700; cursor: pointer;
      box-shadow: 0 4px 20px rgba(0,0,0,0.25);
      font-family: -apple-system, "Segoe UI", "Apple SD Gothic Neo", "Noto Sans KR", sans-serif;
      border: none; transition: transform 0.15s;
    }
    #yu-comment-fab:hover { transform: scale(1.05); box-shadow: 0 6px 28px rgba(0,0,0,0.3); }
    #yu-comment-modal {
      position: fixed; bottom: 90px; right: 24px; z-index: 99999;
      background: white; border: 2px solid #1e40af; border-radius: 12px;
      padding: 20px; width: 380px; max-width: 92vw; display: none;
      box-shadow: 0 10px 40px rgba(0,0,0,0.3);
      font-family: -apple-system, "Segoe UI", "Apple SD Gothic Neo", "Noto Sans KR", sans-serif;
    }
    #yu-comment-modal.open { display: block; }
    #yu-comment-modal h3 {
      margin: 0 0 8px 0; font-size: 15px; color: #1e40af;
    }
    #yu-comment-modal .section-tag {
      display: inline-block; background: #ecfdf5; color: #065f46;
      padding: 3px 8px; border-radius: 4px; font-size: 11px;
      margin-bottom: 10px; font-weight: 600;
    }
    #yu-comment-modal textarea {
      width: 100%; min-height: 110px; padding: 10px;
      border: 1px solid #d1d5db; border-radius: 6px; resize: vertical;
      font-family: inherit; font-size: 13px; box-sizing: border-box;
    }
    #yu-comment-modal .btn-row {
      display: flex; gap: 8px; margin-top: 10px;
    }
    #yu-comment-modal button {
      flex: 1; padding: 8px 12px; border-radius: 6px; border: none;
      cursor: pointer; font-weight: 600; font-size: 13px;
    }
    #yu-comment-modal .btn-save {
      background: #10b981; color: white;
    }
    #yu-comment-modal .btn-cancel {
      background: #e5e7eb; color: #1f2937;
    }
    #yu-comment-modal .btn-view {
      background: #3730a3; color: white; flex: 0 0 auto; padding: 8px 14px;
    }
    #yu-comment-modal .recent {
      margin-top: 12px; padding-top: 10px; border-top: 1px solid #e5e7eb;
      font-size: 12px; color: #6b7280;
    }
    #yu-comment-modal .recent ul { margin: 4px 0; padding-left: 18px; }
    #yu-comment-modal .toast {
      position: fixed; top: 30px; right: 30px; background: #10b981;
      color: white; padding: 10px 16px; border-radius: 6px; z-index: 100000;
      font-size: 13px; font-weight: 600; display: none;
    }
    #yu-comment-modal .toast.show { display: block; animation: yu-slide 0.3s; }
    @keyframes yu-slide { from { opacity: 0; transform: translateY(-10px); } to { opacity: 1; transform: translateY(0); } }
  `;
  document.head.appendChild(style);

  // Helper: get/save comments
  function getComments() {
    try {
      return JSON.parse(localStorage.getItem('yu_comments_v17') || '[]');
    } catch (e) { return []; }
  }
  function saveComments(arr) {
    localStorage.setItem('yu_comments_v17', JSON.stringify(arr));
  }
  function addComment(section, text, author) {
    var entry = {
      section: section,
      text: text,
      author: author || '익명',
      timestamp: new Date().toISOString(),
      url: location.href,
    };
    var arr = getComments();
    arr.push(entry);
    saveComments(arr);
    // Also POST to server (best-effort; if server down, localStorage backup still saved)
    try {
      fetch('/api/comments', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify(entry),
      }).then(function(r) { return r.json(); }).then(function(data) {
        console.log('[comment] server save:', data);
        if (data.ok) {
          showToast('✓ 서버에도 저장됨 (총 ' + data.n_total + '개)');
        }
      }).catch(function(e) {
        console.warn('[comment] server save failed (localStorage backup OK):', e);
      });
    } catch (e) {
      console.warn('[comment] fetch unavailable:', e);
    }
    return arr;
  }
  function showToast(msg) {
    var t = document.getElementById('yu-toast');
    if (!t) {
      t = document.createElement('div');
      t.id = 'yu-toast';
      t.className = 'toast';
      document.body.appendChild(t);
      var s = document.createElement('style');
      s.textContent = '#yu-toast{position:fixed;top:30px;right:30px;background:#10b981;color:white;padding:10px 16px;border-radius:6px;z-index:100000;font-size:13px;font-weight:600;display:none;}#yu-toast.show{display:block;animation:yu-slide 0.3s;}@keyframes yu-slide{from{opacity:0;transform:translateY(-10px);}to{opacity:1;transform:translateY(0);}}';
      document.head.appendChild(s);
    }
    t.textContent = msg;
    t.classList.add('show');
    setTimeout(function() { t.classList.remove('show'); }, 2000);
  }

  // Build FAB button
  var fab = document.createElement('button');
  fab.id = 'yu-comment-fab';
  fab.innerHTML = '💬 의견 남기기';
  fab.title = '유형원 교수님 의견 남기기';
  document.body.appendChild(fab);

  // Build modal
  var modal = document.createElement('div');
  modal.id = 'yu-comment-modal';
  // Quick-pick author shortcuts
  var lastAuthor = localStorage.getItem('yu_comment_author') || '';
  modal.innerHTML = `
    <h3>💬 의견 / 코멘트</h3>
    <div class="section-tag">📍 ${escapeHtml(sectionName)}</div>
    <div style="margin:6px 0 12px 0;background:#f9fafb;padding:10px;border-radius:6px;border:1px solid #e5e7eb;">
      <label style="font-size:12px;color:#374151;display:block;margin-bottom:6px;font-weight:600;">👤 작성자 — 직접 입력 또는 빠른 선택</label>
      <input type="text" id="yu-author-input" value="${escapeHtml(lastAuthor)}" placeholder="✍️ 여기에 이름 직접 입력 (예: 김철수, OO과 OOO 교수)" style="width:100%;padding:10px;border:2px solid #1d4ed8;border-radius:5px;font-size:14px;box-sizing:border-box;font-weight:500;background:white;">
      <div style="font-size:11px;color:#6b7280;margin:6px 0 4px 0;">또는 빠른 선택:</div>
      <div style="display:flex;gap:5px;flex-wrap:wrap;">
        <button type="button" class="quick-author" data-name="유형원 교수님" style="background:#dbeafe;color:#1e40af;border:1px solid #93c5fd;padding:3px 9px;border-radius:14px;font-size:11px;cursor:pointer;">👨‍⚕️ 유형원 교수님</button>
        <button type="button" class="quick-author" data-name="국승호" style="background:#dcfce7;color:#065f46;border:1px solid #86efac;padding:3px 9px;border-radius:14px;font-size:11px;cursor:pointer;">👤 국승호</button>
        <button type="button" class="quick-author" data-name="외부 reviewer" style="background:#fef3c7;color:#92400e;border:1px solid #fcd34d;padding:3px 9px;border-radius:14px;font-size:11px;cursor:pointer;">📝 외부 reviewer</button>
        <button type="button" class="quick-author" data-name="익명" style="background:#f3f4f6;color:#374151;border:1px solid #d1d5db;padding:3px 9px;border-radius:14px;font-size:11px;cursor:pointer;">🕶 익명</button>
      </div>
      <div style="font-size:10px;color:#9ca3af;margin-top:6px;">💡 한 번 입력하면 자동 기억됩니다 (다음 의견에 자동 채워짐)</div>
    </div>
    <textarea id="yu-comment-input" placeholder="이 섹션에 대한 의견을 자유롭게 적어주세요..."></textarea>
    <div class="btn-row">
      <button class="btn-cancel" id="yu-cancel">취소</button>
      <button class="btn-view" id="yu-view">📋 모두 보기</button>
      <button class="btn-save" id="yu-save">저장</button>
    </div>
    <div class="recent" id="yu-recent"></div>
  `;
  document.body.appendChild(modal);

  // Quick author button handlers
  modal.querySelectorAll('.quick-author').forEach(function(btn) {
    btn.addEventListener('click', function() {
      var input = document.getElementById('yu-author-input');
      input.value = btn.getAttribute('data-name');
    });
  });

  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, function(c) {
      return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c];
    });
  }

  function updateRecent() {
    var arr = getComments();
    var thisSection = arr.filter(function(c) { return c.section === sectionName; });
    var recent = document.getElementById('yu-recent');
    if (thisSection.length === 0) {
      recent.innerHTML = '<i>이 섹션에 아직 의견 없음. 전체: ' + arr.length + '개 (모두 보기 클릭)</i>';
    } else {
      var lines = thisSection.slice(-3).reverse().map(function(c) {
        var t = new Date(c.timestamp).toLocaleString('ko-KR');
        var a = escapeHtml(c.author || '익명');
        return '<li><b>[' + a + ']</b> <span style="color:#6b7280;">' + escapeHtml(t) + '</span>: ' + escapeHtml(c.text.substring(0, 50)) + (c.text.length > 50 ? '...' : '') + '</li>';
      });
      recent.innerHTML = '이 섹션 최근 의견 (' + thisSection.length + '개 중 최대 3개):<ul>' + lines.join('') + '</ul>';
    }
  }

  // Open modal with a specific section context (used by inline buttons + FAB)
  function openModal(specificSection) {
    sectionName = specificSection || pageTitle;
    var tag = modal.querySelector('.section-tag');
    if (tag) tag.textContent = '📍 ' + sectionName;
    modal.classList.add('open');
    updateRecent();
    document.getElementById('yu-comment-input').focus();
  }

  // Wire events — main FAB opens modal with page-level section
  fab.addEventListener('click', function() {
    if (modal.classList.contains('open')) {
      modal.classList.remove('open');
    } else {
      openModal(pageTitle);
    }
  });

  // === Inline 💬 buttons on every H2 / H3 / details summary / figure / img ===
  function injectInlineButtons() {
    // Inject CSS for inline buttons
    var inlineStyle = document.createElement('style');
    inlineStyle.textContent = `
      .yu-inline-btn {
        display: inline-block; margin-left: 8px;
        background: rgba(124, 58, 237, 0.1); color: #7c3aed;
        border: 1px solid rgba(124, 58, 237, 0.3); border-radius: 12px;
        padding: 2px 9px; font-size: 11px; font-weight: 600; cursor: pointer;
        font-family: -apple-system, "Segoe UI", "Apple SD Gothic Neo", "Noto Sans KR", sans-serif;
        vertical-align: middle; transition: all 0.15s; line-height: 1.4;
      }
      .yu-inline-btn:hover { background: #7c3aed; color: white; transform: scale(1.05); }
      .yu-inline-btn.has-comments {
        background: #fef3c7; border-color: #f59e0b; color: #92400e;
      }
      .yu-inline-btn.has-comments:hover { background: #f59e0b; color: white; }
      .yu-inline-figbtn {
        position: absolute; top: 8px; right: 8px;
        background: rgba(124, 58, 237, 0.85); color: white;
        border: 0; border-radius: 14px; padding: 4px 10px;
        font-size: 11px; font-weight: 600; cursor: pointer;
        z-index: 10; backdrop-filter: blur(4px);
      }
      .yu-inline-figbtn:hover { background: #7c3aed; }
      .yu-fig-wrap { position: relative; display: inline-block; }
    `;
    document.head.appendChild(inlineStyle);

    // Get current comments to show count badges
    var commentsArr = getComments();
    function countForSection(s) {
      return commentsArr.filter(function(c) { return c.section === s; }).length;
    }

    function makeInlineBtn(sectionLabel, label) {
      var btn = document.createElement('button');
      btn.className = 'yu-inline-btn';
      btn.type = 'button';
      var n = countForSection(sectionLabel);
      btn.innerHTML = n > 0 ? ('💬 ' + n) : ('💬');
      btn.title = '이 섹션에 의견 남기기: ' + sectionLabel;
      if (n > 0) btn.classList.add('has-comments');
      btn.addEventListener('click', function(e) {
        e.preventDefault(); e.stopPropagation();
        openModal(sectionLabel);
      });
      return btn;
    }

    // Find headings (H2, H3) and inject inline buttons
    document.querySelectorAll('h2, h3').forEach(function(h) {
      // Skip if it's inside details summary or modal
      if (h.closest('#yu-comment-modal') || h.closest('.topnav')) return;
      var text = h.textContent.trim().substring(0, 80);
      if (!text) return;
      var sectionLabel = pageSlug + ' · ' + text;
      h.appendChild(makeInlineBtn(sectionLabel, text));
    });

    // Find details > summary and inject (more compact label)
    document.querySelectorAll('details > summary').forEach(function(s) {
      if (s.closest('#yu-comment-modal')) return;
      var text = s.textContent.trim().substring(0, 100);
      if (!text) return;
      var sectionLabel = pageSlug + ' · 펼침 · ' + text;
      // Append button — but be careful not to break click toggle
      var btn = makeInlineBtn(sectionLabel, text);
      s.appendChild(btn);
    });

    // Find <img> tags (figures) — wrap with relative container if needed
    document.querySelectorAll('img').forEach(function(img) {
      if (img.closest('.yu-fig-wrap') || img.closest('#yu-comment-modal') || img.closest('.sidebar')) return;
      var alt = (img.alt || img.title || 'figure').substring(0, 60);
      var src = (img.getAttribute('src') || '').split('/').pop().split('?')[0];
      var sectionLabel = pageSlug + ' · 🖼 ' + (src || alt);
      // Skip very small images (icons)
      if (img.width && img.width < 100) return;

      // Wrap if not already
      var parent = img.parentNode;
      if (parent && parent.classList && !parent.classList.contains('yu-fig-wrap')) {
        var wrap = document.createElement('span');
        wrap.className = 'yu-fig-wrap';
        wrap.style.display = 'inline-block';
        wrap.style.position = 'relative';
        parent.insertBefore(wrap, img);
        wrap.appendChild(img);
        var btn = document.createElement('button');
        btn.className = 'yu-inline-figbtn';
        btn.type = 'button';
        var n = countForSection(sectionLabel);
        btn.innerHTML = n > 0 ? ('💬 ' + n + '개 의견') : '💬 의견';
        btn.title = '이 figure 에 의견: ' + sectionLabel;
        btn.addEventListener('click', function(e) {
          e.preventDefault(); e.stopPropagation();
          openModal(sectionLabel);
        });
        wrap.appendChild(btn);
      }
    });
  }

  // Run after DOM ready (or immediately if already ready)
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', injectInlineButtons);
  } else {
    setTimeout(injectInlineButtons, 100); // small delay for late-rendered content
  }
  document.getElementById('yu-cancel').addEventListener('click', function() {
    modal.classList.remove('open');
    document.getElementById('yu-comment-input').value = '';
  });
  document.getElementById('yu-save').addEventListener('click', function() {
    var input = document.getElementById('yu-comment-input');
    var authorInput = document.getElementById('yu-author-input');
    var text = input.value.trim();
    var author = authorInput.value.trim() || '익명';
    if (!text) {
      showToast('빈 의견은 저장 안 됨');
      return;
    }
    addComment(sectionName, text, author);
    input.value = '';
    // Remember last-used author
    localStorage.setItem('yu_comment_author', author);
    showToast('✓ 저장됨 (' + author + ', ' + getComments().length + '개 누적)');
    updateRecent();
  });
  document.getElementById('yu-view').addEventListener('click', function() {
    location.href = 'comments.html';
  });

  // Close on outside click
  document.addEventListener('click', function(e) {
    if (!modal.contains(e.target) && e.target !== fab && modal.classList.contains('open')) {
      modal.classList.remove('open');
    }
  });
})();
