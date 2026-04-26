"""HTML + Markdown templates for v7_repurposing.py (keeps the main script lean)."""

DASHBOARD_HTML = """<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>v7 Drug Repurposing | THYROID DASH</title>
<link rel="icon" href="../assets/img/favicon.svg" type="image/svg+xml">
<link rel="stylesheet" href="../assets/css/main.css">
<link rel="stylesheet" href="../assets/css/v7_repurposing.css">
</head>
<body class="dash">
<main class="wrap">
<section class="v7-hero">
  <h1>v7 &middot; Drug Repurposing Atlas</h1>
  <p>End-to-end computational repurposing across 8 novel THCA druggable targets, triangulating Open Targets knownDrugs, DrugCentral interactions, and ChEMBL clinical-stage bioactivities. Each (target, drug) pair is scored on a weighted rubric covering development stage, potency, oral bioavailability, patent horizon, mechanism fit, and thyroid safety signals.</p>
</section>
<div class="v7-stat-row">
  <div class="v7-stat"><div class="n">__STAT_TOTAL__</div><div class="l">Total drug-target rows</div></div>
  <div class="v7-stat"><div class="n">__STAT_FDA__</div><div class="l">FDA-approved (phase 4)</div></div>
  <div class="v7-stat"><div class="n">__STAT_P23__</div><div class="l">Clinical phase 2-3</div></div>
  <div class="v7-stat"><div class="n">__STAT_GT05__</div><div class="l">Repurposing score &gt; 0.5</div></div>
</div>
<section class="v7-section">
  <h2>Ranked repurposing candidates</h2>
  <p class="lead">Sort any column; filter by typing a gene, drug, source, or indication.</p>
  <div class="v7-table-wrap">
    <div class="v7-table-tools"><input id="v7-filter" placeholder="filter rows..." autocomplete="off"><span id="v7-count" class="muted"></span></div>
    <table class="v7-tbl" id="v7-tbl">
      <thead><tr>
        <th data-k="target">target</th><th data-k="drug">drug</th><th data-k="score">score</th>
        <th data-k="phase">max phase</th><th data-k="indication">original indication</th>
        <th data-k="source">source</th><th data-k="rationale">rationale</th>
      </tr></thead><tbody></tbody>
    </table>
  </div>
</section>
<section class="v7-section">
  <h2>Top 10 deep-dive candidates</h2>
  <p class="lead">Clinical-trial footprint (ClinicalTrials.gov v2), PubMed activity since 2020 (thyroid OR cancer), and ChEMBL physicochemical/PK context.</p>
  <div class="v7-card-grid">__CARDS__</div>
</section>
<section class="v7-section v7-heat-wrap">
  <h2>Combination hypothesis heatmap</h2>
  <p class="lead">Upper triangle: synergy score across target pairs where both have a drug scoring &ge; 0.4. Lower triangle: number of shared Reactome pathways. Diagonal: pathway count per target.</p>
  <iframe src="__HEAT__" loading="lazy"></iframe>
</section>
<section class="v7-caveat">
  <h3>Caveats</h3>
  <p>All scores are <em>computational hypotheses</em> generated from public bioactivity/target/drug databases (Open Targets, DrugCentral, ChEMBL, ClinicalTrials.gov, PubMed, Reactome). The thyroid-safety term is a heuristic flag derived from text fields only and does NOT replace FDA label review, DDI checks, or QT/hepatotoxicity adjudication. ChEMBL clinical-stage filtering (max_phase &ge; 1) inherits whatever phase tags ChEMBL records, which lag ClinicalTrials.gov. Reactome pathway overlap is a coarse measure of shared biology, not of pharmacodynamic synergy. No combination listed here has been validated experimentally.</p>
</section>
</main>
<script>
const DATA = __TBL_JSON__;
const tbody = document.querySelector('#v7-tbl tbody');
const filterEl = document.getElementById('v7-filter');
const countEl = document.getElementById('v7-count');
let sortKey = 'score', sortDir = 'desc';
function fmtScore(v){ return Number(v).toFixed(3); }
function render(){
  const term = (filterEl.value || '').toLowerCase();
  const rows = DATA.filter(r => !term || JSON.stringify(r).toLowerCase().includes(term));
  rows.sort((a,b) => {
    const va = a[sortKey], vb = b[sortKey];
    const na = Number(va), nb = Number(vb);
    const bothNum = !isNaN(na) && !isNaN(nb);
    if (bothNum) return sortDir==='asc' ? na-nb : nb-na;
    return sortDir==='asc' ? String(va).localeCompare(String(vb)) : String(vb).localeCompare(String(va));
  });
  tbody.innerHTML = rows.map(r => `<tr><td>${r.target}</td><td><b>${r.drug}</b></td><td>${fmtScore(r.score)}</td><td>${r.phase}</td><td>${r.indication||''}</td><td>${r.source||''}</td><td>${r.rationale||''}</td></tr>`).join('');
  countEl.textContent = `${rows.length} / ${DATA.length} rows`;
  document.querySelectorAll('#v7-tbl th').forEach(th => {
    th.classList.remove('sorted-asc','sorted-desc');
    if (th.dataset.k===sortKey) th.classList.add(sortDir==='asc'?'sorted-asc':'sorted-desc');
  });
}
document.querySelectorAll('#v7-tbl th').forEach(th => {
  th.addEventListener('click', () => {
    const k = th.dataset.k;
    if (sortKey===k) sortDir = (sortDir==='asc'?'desc':'asc');
    else { sortKey = k; sortDir = (k==='drug'||k==='target'||k==='source'||k==='indication'||k==='rationale') ? 'asc' : 'desc'; }
    render();
  });
});
filterEl.addEventListener('input', render);
render();
</script>
</body></html>
"""

PAPER_MD = """# v7. Computational drug repurposing for novel THCA targets

## Methods

We assembled an end-to-end repurposing workflow around the eight top-novelty druggable targets identified in our THCA differential-expression and druggability screen (TACSTD2, TMPRSS4, PLEKHA6, CYP1B1, LDLR, GABRB2, B3GNT3, PTPRE). Drug-target evidence was triangulated from three public resources: (i) Open Targets Platform `knownDrugs` (GraphQL v4) after resolving each gene symbol to its Ensembl stable ID via the `search` query; (ii) DrugCentral `drug.target.interaction.tsv.gz` (cached locally, filtered to HGNC matches of the eight symbols); and (iii) ChEMBL `activity.json` with a `pchembl_value >= 6` threshold followed by a molecule-level fetch restricted to `max_phase >= 1` (i.e., at least phase-1 clinical development). Six commercially available thyroid-cancer tyrosine kinase inhibitors (lenvatinib, sorafenib, cabozantinib, vandetanib, selpercatinib, pralsetinib) were excluded so the pipeline surfaces genuinely repurposable assets.

Each (target, drug) record was scored with a weighted rubric: `score = 0.30 * (max_phase / 4) + 0.20 * potency + 0.15 * thyroid_safety + 0.15 * oral_bioavailability + 0.10 * off_patent + 0.10 * mechanism_fit`, with a 10x downweight for withdrawn drugs. Potency was mapped from pChEMBL bins (>=8 = 1.0, >=7 = 0.7, >=6 = 0.4, else 0.2; missing = 0.3). Thyroid safety was a text-derived heuristic flag on QT or hepatotoxicity keywords. Off-patent status used first-approval year cutoffs (<=2005, 2006-2015, >2015). Mechanism fit favored inhibitors/antagonists/blockers (all eight targets are overexpressed in THCA) and penalized agonists/activators. Where multiple sources reported the same drug-target pair we merged phase tags (maximum), pChEMBL values (maximum), and mechanism strings (union).

For the ten highest-scoring unique drugs we layered (a) ClinicalTrials.gov v2 study counts with a thyroid-indication flag, (b) ChEMBL-derived pharmacokinetic tags (oral/parenteral, first-approval year, withdrawal), and (c) a PubMed E-utilities count of papers since 2020 matching `<drug> AND (thyroid OR cancer)`. Combination hypotheses were enumerated across all C(8,2)=28 target pairs where both members had at least one candidate with score >=0.4; pathway overlap was drawn live from Reactome's UniProt-to-pathway mapping with a curated fallback for offline runs. A synergy score combined pathway overlap, FDA-approval status of both partners, and mechanism-class diversity.

## Top-10 repurposing candidates

__TOP_LINES__

## Combination hypotheses

__COMBO_LINES__

A complete 8x8 visualization of synergy scores (upper triangle) and Reactome pathway overlap counts (lower triangle) is provided as an interactive heatmap (`combination_heatmap.html`).

## Summary statistics

- Total drug-target interaction rows after deduplication: **__STAT_TOTAL__**
- FDA-approved (phase 4) hits: **__STAT_FDA__**
- Clinical phase 2-3 hits: **__STAT_P23__**
- Rows with repurposing score > 0.5: **__STAT_GT05__**

## Figures and tables referenced

- Figure v7-1: interactive repurposing table (`v7_repurposing_ranked.tsv`).
- Figure v7-2: 8x8 combination heatmap (`combination_heatmap.html`).
- Table v7-A: `v7_combination_hypotheses.tsv`.

## Limitations

The analysis is entirely computational and public-data-driven: (1) bioactivity measurements in ChEMBL span heterogeneous assays, so the potency term is a coarse binning rather than calibrated Ki/IC50; (2) our thyroid-safety term is a keyword heuristic and cannot replace FDA label review, DDI profiling, or in vivo QT/hepatotoxicity work; (3) Open Targets phase tags occasionally lag ClinicalTrials.gov, which is why STEP 5 re-queries CT.gov; (4) Reactome overlap is a necessary but not sufficient proxy for pharmacodynamic synergy, and all combination hypotheses require wet-lab validation before translational claims; (5) the eight targets are transcriptionally overexpressed in THCA but proteomic/surface-exposure confirmation is a separate workstream. Future work should layer PDX efficacy signals, patient-derived-organoid sensitivity data, and structure-based DTI refinement for the novel targets without co-crystal structures.
"""
