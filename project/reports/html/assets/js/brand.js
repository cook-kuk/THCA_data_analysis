/* THYRAI brand.js — hamburger toggle + pipeline hover tooltips.
   No external dependencies. CSP-safe. */
(function () {
  "use strict";

  function onReady(fn) {
    if (document.readyState !== "loading") fn();
    else document.addEventListener("DOMContentLoaded", fn);
  }

  function initHamburger() {
    var btn = document.querySelector(".topnav__hamburger");
    if (!btn) return;
    btn.addEventListener("click", function () {
      document.body.classList.toggle("nav-open");
      var expanded = document.body.classList.contains("nav-open");
      btn.setAttribute("aria-expanded", expanded ? "true" : "false");
    });
  }

  function initPipelineHover() {
    var rows = document.querySelectorAll(".pipeline-row[data-target]");
    rows.forEach(function (row) {
      // Tooltip content already rendered in .pipeline-tooltip
      // Hover styling is done in CSS; here we could add focus handling.
      row.addEventListener("focusin", function () { row.classList.add("is-hover"); });
      row.addEventListener("focusout", function () { row.classList.remove("is-hover"); });
    });
  }

  function initFigureLazyLoad() {
    // Defer setting iframe src until the element is near viewport for perf.
    var iframes = document.querySelectorAll("iframe[data-src]");
    if (!iframes.length) return;
    if (!("IntersectionObserver" in window)) {
      iframes.forEach(function (f) { f.src = f.dataset.src; });
      return;
    }
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting) {
          e.target.src = e.target.dataset.src;
          io.unobserve(e.target);
        }
      });
    }, { rootMargin: "200px" });
    iframes.forEach(function (f) { io.observe(f); });
  }

  function ensureExplainBubble() {
    if (!document.body || !document.body.classList.contains("thyrai")) return;
    if (window.ExplainBubble) {
      if (window.ExplainBubble.reinjectButtons) window.ExplainBubble.reinjectButtons();
      return;
    }
    var inPages = /\/pages\/[^/]+$/.test(window.location.pathname);
    var cssHref = inPages ? "../assets/css/explain-bubble.css" : "assets/css/explain-bubble.css";
    var jsHref = inPages ? "../assets/js/explain-bubble.js" : "assets/js/explain-bubble.js";

    if (!document.querySelector('link[href$="explain-bubble.css"]')) {
      var link = document.createElement("link");
      link.rel = "stylesheet";
      link.href = cssHref;
      document.head.appendChild(link);
    }
    if (!document.querySelector('script[src$="explain-bubble.js"]')) {
      var script = document.createElement("script");
      script.src = jsHref;
      script.defer = true;
      document.body.appendChild(script);
    }
  }

  function ensureSiteMap() {
    if (!document.body || !document.body.classList.contains("thyrai")) return;
    if (window.ThyraiSiteMap) {
      if (window.ThyraiSiteMap.init) window.ThyraiSiteMap.init();
      return;
    }
    var inPages = /\/pages\/[^/]+$/.test(window.location.pathname);
    var jsHref = inPages ? "../assets/js/site-map.js" : "assets/js/site-map.js";
    if (!document.querySelector('script[src$="site-map.js"]')) {
      var script = document.createElement("script");
      script.src = jsHref;
      script.defer = true;
      document.body.appendChild(script);
    }
  }

  onReady(function () {
    initHamburger();
    initPipelineHover();
    initFigureLazyLoad();
    ensureExplainBubble();
    ensureSiteMap();
  });
})();
