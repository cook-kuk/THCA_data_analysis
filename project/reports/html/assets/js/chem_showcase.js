/* chem_showcase.js — UI wiring for platform_chem_showcase.html
   Visual system inspired by modern drug-discovery platform conventions.
   Not affiliated with any commercial AI drug-discovery company.
*/
(function () {
  "use strict";

  var DATA_URL = "../assets/data/chem_targets.json";
  var BASE = "..";

  var state = {
    targets: [],
    byGene: {},
    currentGene: "CYP1B1",
  };

  // ---------- helpers ----------
  function $(sel, scope) { return (scope || document).querySelector(sel); }
  function $$(sel, scope) { return Array.from((scope || document).querySelectorAll(sel)); }

  function esc(s) {
    return String(s == null ? "" : s)
      .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;").replace(/'/g, "&#39;");
  }

  function setText(sel, text) { var el = $(sel); if (el) el.textContent = text; }

  function badgeClassForClassification(c) {
    if (!c) return "classification-badge";
    if (c === "validated_target") return "classification-badge classification-badge--validated";
    if (c === "emerging_target")  return "classification-badge classification-badge--emerging";
    return "classification-badge classification-badge--novel";
  }

  function cohensDFromLog2FC(log2fc) {
    // Rough proxy — readers see this as an approximate signal magnitude.
    if (log2fc == null) return null;
    return Math.abs(log2fc) * 0.6;   // heuristic mapping; caveat in tooltip
  }

  function replicationRate(bs) {
    var a = bs.replicated_27155 || 0;
    var b = bs.replicated_126698 || 0;
    return Math.round((a + b) / 2 * 100);
  }

  // ---------- hero: rdkit-js SVG render ----------
  function renderHero(smiles) {
    var slot = document.getElementById("hero-mol-svg");
    if (!slot) return;
    // Provide a default safe SVG immediately so layout doesn't jump.
    slot.innerHTML = _fallbackMoleculeSvg();
    if (!window.RDKit || !smiles) return;
    try {
      var mol = window.RDKit.get_mol(smiles);
      if (!mol) return;
      var svg = mol.get_svg_with_highlights(JSON.stringify({
        width: 520, height: 520, bondLineWidth: 2.2,
        backgroundColour: [0, 0, 0, 0],
        highlightColour: [0.96, 0.65, 0.14, 1.0],
      }));
      // Post-process: replace any pure black strokes with a lighter gray for
      // legibility on black background, and tint carbons amber.
      svg = svg
        .replace(/stroke='#000000'/g, "stroke='#E7E7E7'")
        .replace(/stroke="#000000"/g, 'stroke="#E7E7E7"')
        .replace(/fill='#000000'/g, "fill='#F5A623'")
        .replace(/fill="#000000"/g, 'fill="#F5A623"');
      slot.innerHTML = svg;
      mol.delete();
    } catch (e) {
      console.warn("[chem_showcase] RDKit render failed:", e);
    }
  }

  function _fallbackMoleculeSvg() {
    // Stylized hexagon-ring motif as a pre-RDKit placeholder.
    return [
      '<svg viewBox="0 0 200 200" xmlns="http://www.w3.org/2000/svg">',
      '<g fill="none" stroke="#F5A623" stroke-width="2" stroke-linecap="round">',
      '<polygon points="100,40 150,70 150,130 100,160 50,130 50,70" />',
      '<line x1="100" y1="40" x2="100" y2="10" />',
      '<line x1="150" y1="70" x2="178" y2="55" />',
      '<line x1="150" y1="130" x2="178" y2="145" />',
      '<line x1="50" y1="130" x2="22" y2="145" />',
      '<line x1="50" y1="70" x2="22" y2="55" />',
      '<circle cx="100" cy="40" r="4" fill="#F5A623" />',
      '<circle cx="150" cy="70" r="4" fill="#F5A623" />',
      '<circle cx="150" cy="130" r="4" fill="#F5A623" />',
      '<circle cx="100" cy="160" r="4" fill="#F5A623" />',
      '<circle cx="50" cy="130" r="4" fill="#F5A623" />',
      '<circle cx="50" cy="70" r="4" fill="#F5A623" />',
      '</g>',
      '</svg>'
    ].join("");
  }

  // ---------- viewer panel population ----------
  function renderTargetMeta(t) {
    if (!t) return;
    $("[data-meta-symbol]").textContent = t.gene;
    $("[data-meta-pref-name]").textContent = t.chembl_pref_name || "—";
    var badge = $("[data-meta-badge]");
    badge.className = badgeClassForClassification(t.classification);
    badge.textContent = (t.classification || "").replace("_", " ");

    var cohensD = cohensDFromLog2FC(t.biomarker_stats && t.biomarker_stats.log2FC);
    $("[data-meta-cohens]").textContent = cohensD == null ? "—" : cohensD.toFixed(2);
    var bs = t.biomarker_stats || {};
    $("[data-meta-replication]").textContent = replicationRate(bs) + "%";
    $("[data-meta-log2fc]").textContent = (bs.log2FC == null ? "—" : bs.log2FC.toFixed(2));
    $("[data-meta-fdr]").textContent = (bs.fdr == null ? "—" : bs.fdr.toExponential(1));
    $("[data-meta-novelty]").textContent = (t.novelty_score == null ? "—" : t.novelty_score.toFixed(1));

    var rcsbLink = $("[data-meta-rcsb]");
    if (t.structure && t.structure.source_kind === "pdb" && t.structure.identifier) {
      rcsbLink.href = "https://www.rcsb.org/structure/" + t.structure.identifier;
      rcsbLink.style.display = "";
      rcsbLink.innerHTML = "View in RCSB &nbsp;&rarr;";
    } else if (t.structure && t.structure.source_kind === "alphafold" && t.structure.identifier) {
      rcsbLink.href = "https://alphafold.ebi.ac.uk/entry/" + t.structure.identifier;
      rcsbLink.style.display = "";
      rcsbLink.innerHTML = "View in AlphaFold DB &nbsp;&rarr;";
    } else {
      rcsbLink.style.display = "none";
    }

    var rationaleEl = $("[data-meta-rationale]");
    if (rationaleEl) rationaleEl.textContent = t.rationale || "No rationale paragraph available — see the biomarker-to-drug report for full context.";
  }

  function renderCompoundList(t) {
    var ul = $("[data-compound-list]");
    ul.innerHTML = "";
    var topN = (t.compounds || []).slice(0, 3);
    if (topN.length === 0) {
      ul.innerHTML = '<li class="compound-list__item" style="opacity:0.6"><span class="compound-list__name">No ChEMBL compounds available</span></li>';
      return;
    }
    topN.forEach(function (c, idx) {
      var li = document.createElement("li");
      li.className = "compound-list__item";
      li.innerHTML =
        '<input type="radio" name="compound" value="' + esc(c.chembl_id) + '" ' + (idx === 0 ? "checked" : "") + '>' +
        '<span class="compound-list__name">' + esc(c.name || c.chembl_id) + '</span>' +
        '<span class="compound-list__pch">pCh ' + (c.pchembl == null ? "—" : c.pchembl) + '</span>';
      ul.appendChild(li);
    });
  }

  function renderBindingResidues(residues) {
    var bar = $("[data-residue-bar]");
    bar.innerHTML = "";
    if (!residues || !residues.length) {
      bar.innerHTML = '<span class="binding-residue-label" style="opacity:0.7">No bound ligand in reference structure</span>';
      return;
    }
    residues.forEach(function (r) {
      var span = document.createElement("span");
      span.className = "binding-residue-label";
      span.innerHTML = esc(r.resname) + esc(r.resi) + " <span class='dist'>" + r.distance_ang.toFixed(1) + "Å</span>";
      bar.appendChild(span);
    });
  }

  // ---------- charts (Plotly) ----------
  function drawPotencyChart() {
    if (!window.Plotly) return;
    var container = $("#potency-chart");
    if (!container) return;
    var focus = ["CYP1B1", "LDLR", "GABRB2", "TACSTD2"];
    var rows = focus.map(function (g) { return state.byGene[g]; }).filter(Boolean);
    var y = [], lo = [], hi = [], best = [], counts = [];
    rows.forEach(function (t) {
      var pchs = (t.compounds || []).slice(0, 3).map(function (c) { return c.pchembl; }).filter(function (v) { return v != null; });
      if (!pchs.length) {
        y.push(t.gene); lo.push(0); hi.push(0); best.push(null);
      } else {
        y.push(t.gene);
        lo.push(Math.min.apply(null, pchs));
        hi.push(Math.max.apply(null, pchs));
        best.push(Math.max.apply(null, pchs));
      }
      counts.push((t.compounds || []).length);
    });
    // Bars = hi - lo; use the base trick for Plotly horizontal.
    var traceRange = {
      type: "bar",
      orientation: "h",
      y: y,
      x: hi.map(function (v, i) { return v - lo[i]; }),
      base: lo,
      name: "pChEMBL range (top 3)",
      marker: { color: "rgba(245,166,35,0.6)", line: { color: "rgba(245,166,35,1)", width: 1 } },
      hovertemplate: "%{y}<br>pChEMBL range: %{base:.2f}–%{x:.2f}<extra></extra>",
    };
    var traceBest = {
      type: "scatter",
      mode: "markers+text",
      y: y,
      x: best,
      text: best.map(function (v, i) { return v == null ? "" : "best " + v.toFixed(2) + " · n=" + counts[i]; }),
      textposition: "middle right",
      textfont: { color: "#F5A623", family: "JetBrains Mono,monospace", size: 11 },
      name: "Top pChEMBL",
      marker: { size: 11, color: "#F5A623", symbol: "diamond", line: { color: "#000", width: 1 } },
      hovertemplate: "%{y}<br>Top pChEMBL: %{x:.2f}<extra></extra>",
    };
    var layout = {
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      margin: { l: 90, r: 60, t: 20, b: 40 },
      xaxis: { title: "pChEMBL (higher = more potent, approx. −log10[M])", color: "#CACACA", gridcolor: "rgba(255,255,255,0.06)", zerolinecolor: "rgba(255,255,255,0.12)", range: [4, 10] },
      yaxis: { color: "#CACACA", tickfont: { family: "JetBrains Mono,monospace" }, automargin: true },
      showlegend: false,
      font: { family: "Inter, Helvetica Neue, sans-serif" },
      bargap: 0.35,
    };
    Plotly.newPlot(container, [traceRange, traceBest], layout, { displayModeBar: false, responsive: true });
  }

  function drawSankey() {
    if (!window.Plotly) return;
    var el = $("#sankey-chart");
    if (!el) return;
    var labels = [
      "51,711 genes",                // 0
      "2,773 replicated biomarkers", // 1
      "Non-replicated (dropped)",    // 2
      "8 druggable targets",         // 3
      "No ChEMBL match",             // 4
      "15 compounds (pChEMBL ≥ 6)",  // 5
    ];
    var source = [0, 0, 1, 1, 3, 3];
    var target = [1, 2, 3, 4, 5, 4];
    var value  = [2773, 51711 - 2773, 8, 2773 - 8, 15, 8 - 3];
    var linkLabel = [
      "Welch t-test + BH-FDR passed (TCGA)",
      "Filtered by FDR / effect size threshold",
      "ChEMBL target match + thyroid literature anchor",
      "No chemical matter in ChEMBL",
      "pChEMBL ≥ 6.0 assayed compounds",
      "Greenfield targets (no compounds yet)",
    ];
    var trace = {
      type: "sankey",
      orientation: "h",
      node: {
        pad: 14,
        thickness: 16,
        line: { color: "rgba(255,255,255,0.12)", width: 1 },
        label: labels,
        color: ["#3A3A3A", "#F5A623", "#222", "#F5A623", "#222", "#F5A623"],
      },
      link: {
        source: source,
        target: target,
        value: value,
        label: linkLabel,
        color: [
          "rgba(245,166,35,0.4)", "rgba(255,255,255,0.06)",
          "rgba(245,166,35,0.4)", "rgba(255,255,255,0.06)",
          "rgba(245,166,35,0.4)", "rgba(255,255,255,0.06)",
        ],
        hovertemplate: "%{label}<br>count: %{value}<extra></extra>",
      },
    };
    var layout = {
      paper_bgcolor: "rgba(0,0,0,0)",
      font: { color: "#CACACA", family: "JetBrains Mono, monospace", size: 11 },
      margin: { l: 0, r: 0, t: 16, b: 16 },
      height: 420,
    };
    Plotly.newPlot(el, [trace], layout, { displayModeBar: false, responsive: true });
  }

  // ---------- init ----------
  function wireControls() {
    var sel = $("#target-select");
    sel.addEventListener("change", function () {
      state.currentGene = sel.value;
      var t = state.byGene[state.currentGene];
      renderTargetMeta(t);
      renderCompoundList(t);
      window.MolViewer.loadTarget(state.currentGene);
    });

    $$(".style-toggle button").forEach(function (btn) {
      btn.addEventListener("click", function () {
        $$(".style-toggle button").forEach(function (b) { b.setAttribute("aria-pressed", "false"); });
        btn.setAttribute("aria-pressed", "true");
        window.MolViewer.setStyle(btn.dataset.style);
      });
    });
    $$(".scheme-toggle button").forEach(function (btn) {
      btn.addEventListener("click", function () {
        $$(".scheme-toggle button").forEach(function (b) { b.setAttribute("aria-pressed", "false"); });
        btn.setAttribute("aria-pressed", "true");
        window.MolViewer.setColorScheme(btn.dataset.scheme);
      });
    });

    $("[data-action='spin']").addEventListener("click", function () {
      window.MolViewer.toggleSpin();
      this.setAttribute("aria-pressed", window.viewerState.spinning ? "true" : "false");
    });
    $("[data-action='reset']").addEventListener("click", function () { window.MolViewer.resetView(); });
    $("[data-action='png']").addEventListener("click", function () { window.MolViewer.downloadPNG(); });
    $("[data-action='fullscreen']").addEventListener("click", function () { window.MolViewer.fullscreen(); });

    // Mobile collapse for the selector drawer
    var drawer = $(".target-selector");
    var toggle = $(".target-selector__toggle");
    if (toggle && drawer) {
      toggle.addEventListener("click", function () {
        var collapsed = drawer.getAttribute("data-collapsed") === "true";
        drawer.setAttribute("data-collapsed", collapsed ? "false" : "true");
        toggle.setAttribute("aria-expanded", collapsed ? "true" : "false");
      });
    }

    // Mobile hamburger wiring (match brand.css behavior)
    var hamb = document.querySelector(".topnav__hamburger");
    if (hamb) {
      hamb.addEventListener("click", function () {
        document.body.classList.toggle("nav-open");
        var expanded = document.body.classList.contains("nav-open");
        hamb.setAttribute("aria-expanded", expanded ? "true" : "false");
      });
    }
  }

  function populateTargetSelect() {
    var sel = $("#target-select");
    sel.innerHTML = "";
    state.targets.forEach(function (t) {
      var opt = document.createElement("option");
      opt.value = t.gene;
      opt.textContent = t.gene + " — " + (t.chembl_pref_name || t.classification || "");
      if (t.gene === state.currentGene) opt.selected = true;
      sel.appendChild(opt);
    });
  }

  async function boot() {
    var resp = await fetch(DATA_URL, { cache: "no-cache" });
    if (!resp.ok) {
      console.error("[chem_showcase] payload fetch failed:", resp.status);
      return;
    }
    var data = await resp.json();
    state.targets = data.targets || [];
    state.byGene = {};
    state.targets.forEach(function (t) { state.byGene[t.gene] = t; });

    // Hero: top compound for CYP1B1 (CHEMBL3132932 — highest pChEMBL in the set)
    var heroTarget = state.byGene["CYP1B1"] || state.targets[0];
    var heroSmiles = heroTarget && heroTarget.compounds && heroTarget.compounds[0] && heroTarget.compounds[0].smiles;
    if (window.initRDKitModule) {
      window.initRDKitModule().then(function (RDKit) {
        window.RDKit = RDKit;
        renderHero(heroSmiles);
      }).catch(function () { renderHero(null); });
    } else {
      renderHero(null);
    }

    // Selector UI
    populateTargetSelect();
    renderTargetMeta(heroTarget);
    renderCompoundList(heroTarget);

    // Viewer
    window.MolViewer.init({
      container: document.getElementById("viewer-3d"),
      baseDataUrl: BASE,
      datasetByGene: state.byGene,
      onBindingResidues: renderBindingResidues,
    });
    window.MolViewer.loadTarget("CYP1B1");

    // Charts
    drawPotencyChart();
    drawSankey();

    // Controls
    wireControls();
  }

  document.addEventListener("DOMContentLoaded", boot);
})();
