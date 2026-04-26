/* mol_viewer.js — 3Dmol.js wrapper
   Exposes a tiny, state-aware API used by chem_showcase.js.
   Visual system inspired by modern drug-discovery platform conventions.
*/
(function () {
  "use strict";

  var HYDROPHOBIC_RES = new Set(["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "PRO", "TRP", "CYS"]);
  var POLAR_RES       = new Set(["SER", "THR", "ASN", "GLN", "TYR", "HIS"]);
  var CHARGED_POS     = new Set(["LYS", "ARG"]);
  var CHARGED_NEG     = new Set(["ASP", "GLU"]);

  // Shared state — read by the UI layer
  var state = {
    currentTarget: null,
    currentCompound: null,
    style: "cartoon",
    scheme: "spectrum",
    spinning: true,
  };

  // PDB string cache: target name -> raw PDB text
  var pdbCache = new Map();
  // SDF cache: path -> raw SDF text
  var sdfCache = new Map();

  var viewer = null;
  var container = null;
  var datasetByGene = {};
  var baseDataUrl = "";
  var onBindingResidues = null;   // hook back into the page

  function _fetchText(url) {
    if (pdbCache.has(url)) return Promise.resolve(pdbCache.get(url));
    return fetch(url, { cache: "force-cache" })
      .then(function (r) { if (!r.ok) throw new Error("HTTP " + r.status + " " + url); return r.text(); })
      .then(function (txt) { pdbCache.set(url, txt); return txt; });
  }

  function _fetchSdf(url) {
    if (sdfCache.has(url)) return Promise.resolve(sdfCache.get(url));
    return fetch(url, { cache: "force-cache" })
      .then(function (r) { if (!r.ok) throw new Error("HTTP " + r.status + " " + url); return r.text(); })
      .then(function (txt) { sdfCache.set(url, txt); return txt; });
  }

  function _applyProteinStyle() {
    if (!viewer) return;
    // Clear all protein styles first
    viewer.setStyle({ hetflag: false }, {});

    var colorfunc = _proteinColor();

    if (state.style === "cartoon" || state.style === "surface") {
      viewer.setStyle({ hetflag: false }, { cartoon: colorfunc });
    } else if (state.style === "stick") {
      viewer.setStyle({ hetflag: false }, { stick: { radius: 0.18, colorscheme: "chainHetatm" } });
    } else if (state.style === "sphere") {
      viewer.setStyle({ hetflag: false }, { sphere: { scale: 0.9, colorscheme: "chainHetatm" } });
    }

    // remove any lingering surface before re-adding
    viewer.removeAllSurfaces();
    if (state.style === "surface") {
      var surfOpts = { opacity: 0.55 };
      if (state.scheme === "hydrophobicity") {
        surfOpts.colorfunc = function (atom) { return _hydroColor(atom); };
      } else {
        surfOpts.color = "#F5A623";
        surfOpts.opacity = 0.35;
      }
      try { viewer.addSurface(window.$3Dmol.SurfaceType.VDW, surfOpts, { hetflag: false }); } catch (e) { /* noop */ }
    }
  }

  function _proteinColor() {
    switch (state.scheme) {
      case "secondary":       return { color: "spectrum" }; // fallback
      case "chain":           return { colorscheme: "chainHetatm" };
      case "hydrophobicity":  return { colorfunc: function (atom) { return _hydroColor(atom); } };
      case "spectrum":
      default:                return { color: "spectrum" };
    }
  }

  function _hydroColor(atom) {
    var rn = (atom && atom.resn) || "";
    if (HYDROPHOBIC_RES.has(rn)) return "#F5A623";
    if (CHARGED_POS.has(rn))     return "#4EA8FF";
    if (CHARGED_NEG.has(rn))     return "#FF5C7A";
    if (POLAR_RES.has(rn))       return "#C8C8C8";
    return "#6A6A6A";
  }

  function _applyLigandStyle(sdfModel) {
    if (!viewer || !sdfModel) return;
    // ligand as sticks with element colors + transparent spheres for glow
    viewer.setStyle({ model: sdfModel }, {
      stick: { radius: 0.15, colorscheme: "default" },
      sphere: { scale: 0.35, opacity: 0.18, color: "#F5A623" },
    });
  }

  function _applyCoLigandStyle() {
    if (!viewer) return;
    // co-crystal ligand lives as HETATM in the PDB; style separately
    viewer.setStyle({ hetflag: true, byres: false, resn: "HOH", invert: true },
      { stick: { radius: 0.18, colorscheme: "default" },
        sphere: { scale: 0.3, opacity: 0.22, color: "#F5A623" } });
    viewer.setStyle({ resn: "HOH" }, {}); // hide water
  }

  function _applyBindingSiteSticks(residues) {
    if (!viewer || !residues || !residues.length) return;
    residues.forEach(function (r) {
      var sel = { chain: r.chain, resi: r.resi };
      viewer.addStyle(sel, { stick: { radius: 0.12, color: "#F5A623" } });
      try {
        viewer.addLabel(r.resname + r.resi, {
          position: undefined,
          backgroundColor: "rgba(0,0,0,0.75)",
          backgroundOpacity: 0.75,
          fontColor: "#F5A623",
          fontSize: 10,
          borderThickness: 0,
          alignment: "centerCenter",
          inFront: true,
        }, sel);
      } catch (e) { /* labels optional */ }
    });
  }

  async function loadTarget(name) {
    var ds = datasetByGene[name];
    if (!ds) { console.warn("[mol_viewer] unknown target:", name); return; }
    state.currentTarget = name;
    state.currentCompound = (ds.compounds && ds.compounds[0] && ds.compounds[0].chembl_id) || null;

    var loadingEl = document.querySelector("[data-viewer-loading]");
    if (loadingEl) loadingEl.hidden = false;

    if (viewer) {
      try { viewer.clear(); } catch (e) { /* noop */ }
    }

    var bindingResidues = (ds.structure && ds.structure.binding_site_residues) || [];
    if (!ds.structure || !ds.structure.file) {
      if (loadingEl) loadingEl.hidden = false;
      if (loadingEl) loadingEl.textContent = "No structure available — compound-only view.";
      if (typeof onBindingResidues === "function") onBindingResidues([]);
      return;
    }

    try {
      var pdbText = await _fetchText(baseDataUrl + "/" + ds.structure.file.replace(/^assets\//, "assets/"));
      viewer.addModel(pdbText, "pdb");
      _applyProteinStyle();

      // Co-crystal ligand (already part of the PDB as HETATM)
      _applyCoLigandStyle();

      // Overlay the top compound's RDKit-generated 3D SDF when available
      if (ds.compounds && ds.compounds[0] && ds.compounds[0].sdf) {
        try {
          var sdfText = await _fetchSdf(baseDataUrl + "/" + ds.compounds[0].sdf);
          var preCount = viewer.getModel() ? 1 : 0;
          var sdfModel = viewer.addModel(sdfText, "sdf");
          // translate overlay to avoid clash — keep simple: leave at origin; users can orbit
          _applyLigandStyle(sdfModel);
        } catch (sdfErr) {
          console.warn("[mol_viewer] SDF overlay failed:", sdfErr);
        }
      }

      _applyBindingSiteSticks(bindingResidues);

      viewer.setBackgroundColor("#000000");
      viewer.zoomTo();
      if (state.spinning) viewer.spin("y", 0.6);
      viewer.render();
    } catch (err) {
      console.error("[mol_viewer] loadTarget error:", err);
      if (loadingEl) loadingEl.textContent = "Structure load failed.";
    } finally {
      if (loadingEl) loadingEl.hidden = true;
    }

    if (typeof onBindingResidues === "function") onBindingResidues(bindingResidues);
  }

  function setStyle(s) {
    state.style = s;
    _applyProteinStyle();
    _applyCoLigandStyle();
    var ds = datasetByGene[state.currentTarget];
    if (ds && ds.structure) _applyBindingSiteSticks(ds.structure.binding_site_residues || []);
    viewer && viewer.render();
  }

  function setColorScheme(s) {
    state.scheme = s;
    _applyProteinStyle();
    _applyCoLigandStyle();
    var ds = datasetByGene[state.currentTarget];
    if (ds && ds.structure) _applyBindingSiteSticks(ds.structure.binding_site_residues || []);
    viewer && viewer.render();
  }

  function toggleSpin() {
    state.spinning = !state.spinning;
    if (!viewer) return;
    if (state.spinning) viewer.spin("y", 0.6);
    else viewer.spin(false);
  }

  function resetView() {
    if (!viewer) return;
    viewer.zoomTo();
    viewer.render();
  }

  function downloadPNG() {
    if (!viewer) return;
    try {
      var dataUrl = viewer.pngURI();
      var a = document.createElement("a");
      a.href = dataUrl;
      a.download = "thyrai_chem_" + (state.currentTarget || "view") + ".png";
      document.body.appendChild(a); a.click(); a.remove();
    } catch (e) {
      console.warn("[mol_viewer] downloadPNG failed:", e);
    }
  }

  function fullscreen() {
    if (!container) return;
    if (document.fullscreenElement) {
      document.exitFullscreen && document.exitFullscreen();
    } else {
      container.requestFullscreen && container.requestFullscreen();
    }
  }

  function init(opts) {
    container = opts.container;
    baseDataUrl = opts.baseDataUrl || "..";
    datasetByGene = opts.datasetByGene || {};
    onBindingResidues = opts.onBindingResidues || null;
    if (!window.$3Dmol) {
      console.error("[mol_viewer] 3Dmol.js not loaded");
      return;
    }
    viewer = window.$3Dmol.createViewer(container, { backgroundColor: "#000000", antialias: true });
    window.addEventListener("resize", function () { if (viewer) viewer.resize(); });
    return viewer;
  }

  // Expose public API
  window.viewerState = state;
  window.MolViewer = {
    init: init,
    loadTarget: loadTarget,
    setStyle: setStyle,
    setColorScheme: setColorScheme,
    toggleSpin: toggleSpin,
    resetView: resetView,
    downloadPNG: downloadPNG,
    fullscreen: fullscreen,
    state: state,
  };
})();
