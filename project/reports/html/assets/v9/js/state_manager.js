/* State manager: URL state encoding + keyboard shortcuts. */
(function (global) {
  let app = null;

  function updateUrl(a = app) {
    if (!a) return;
    const q = new URLSearchParams({
      cancer: a.selectedCancer,
      classifier: a.selectedClassifier,
      flip: a.combatOn ? 'after' : 'before',
    });
    history.replaceState({}, '', location.pathname + '?' + q.toString());
  }

  function keydown(e) {
    // ignore if typing into input
    const tag = (e.target.tagName || '').toLowerCase();
    if (tag === 'input' || tag === 'textarea' || tag === 'select') return;
    const cancers = app.cancers;
    const key = e.key.toLowerCase();

    if (key === ' ') {
      e.preventDefault();
      app.toggleFlip();
      return;
    }
    if (/^[1-5]$/.test(e.key)) {
      const idx = parseInt(e.key, 10) - 1;
      if (cancers[idx]) app.selectCancer(cancers[idx]);
      return;
    }
    if (key === 'j' || key === 'arrowdown') {
      e.preventDefault();
      const i = app.sections.findIndex(s => s.id === app.activeSection);
      const next = app.sections[Math.min(app.sections.length - 1, i + 1)];
      if (next) app.goTo(next.id);
      return;
    }
    if (key === 'k' || key === 'arrowup') {
      e.preventDefault();
      const i = app.sections.findIndex(s => s.id === app.activeSection);
      const prev = app.sections[Math.max(0, i - 1)];
      if (prev) app.goTo(prev.id);
      return;
    }
    if (key === 'g') { app.goTo('section-hero'); return; }
    if (key === 's') { app.shareLink(); return; }
    if (e.key === '?') { app.showShortcuts = true; return; }
    if (e.key === 'Escape') {
      app.showShortcuts = false;
      app.aiOpen = false;
      return;
    }
  }

  const StateManager = {
    init(a) {
      app = a;
      window.addEventListener('keydown', keydown);
      updateUrl(a);
    },
    updateUrl,
  };
  global.StateManager = StateManager;
})(window);
