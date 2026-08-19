from __future__ import annotations

import csv
import html
import shutil
from pathlib import Path

ROOT = Path('/home/seungho/personal/THCA_data_analysis')
HUB = ROOT / 'project/papers_hub_2026_05_04'
LIVE_HUB = Path('/var/www/papers/papers_hub_2026_05_04')
RESULT = ROOT / 'project/results/cell_therapy_independent_series_2026_05_11'
ASSET_DIR = HUB / 'assets/cell_therapy_independent_series_2026_05_11'
LIVE_ASSET_DIR = LIVE_HUB / 'assets/cell_therapy_independent_series_2026_05_11'

PROMISING = ROOT / 'project/results/therapy_spectrum_promising_master_2026_05_11/therapy_spectrum_promising_master.tsv'
VALMAT = ROOT / 'project/results/therapy_spectrum_validation_matrix_2026_05_11/therapy_spectrum_validation_matrix.tsv'
CAR_T = ROOT / 'project/results/car_t_bundle_2026_05_11/car_t_bundle.tsv'
CAPTION = ROOT / 'project/results/cell_therapy_caption_atlas_2026_05_11/cell_therapy_caption_atlas.tsv'
HARD = ROOT / 'project/results/cell_therapy_hard_dossier_2026_05_11/cell_therapy_hard_dossier.tsv'
ULTRA = ROOT / 'project/results/cell_therapy_ultra_analysis_2026_05_11/cell_therapy_ultra_analysis.tsv'
MEGA = ROOT / 'project/results/cell_therapy_megapack_2026_05_11/cell_therapy_megapack.tsv'

LIVE_BASE = 'http://40.82.129.113/papers_hub_2026_05_04/'

FILES = [
    ('Hub index', HUB / 'index.html', LIVE_BASE + 'index.html'),
    ('CAR-T bundle', HUB / 'car_t_bundle_2026_05_11.html', LIVE_BASE + 'car_t_bundle_2026_05_11.html'),
    ('Caption atlas', HUB / 'cell_therapy_caption_atlas_2026_05_11.html', LIVE_BASE + 'cell_therapy_caption_atlas_2026_05_11.html'),
    ('Hard dossier', HUB / 'cell_therapy_hard_dossier_2026_05_11.html', LIVE_BASE + 'cell_therapy_hard_dossier_2026_05_11.html'),
    ('Ultra analysis', HUB / 'cell_therapy_ultra_analysis_2026_05_11.html', LIVE_BASE + 'cell_therapy_ultra_analysis_2026_05_11.html'),
    ('Mega pack', HUB / 'cell_therapy_megapack_2026_05_11.html', LIVE_BASE + 'cell_therapy_megapack_2026_05_11.html'),
]

MODS = {
    'CAR-T': {
        'rank': 1,
        'anchor': 'GSE290722',
        'title': 'CAR-T response durability',
        'verdict': 'strong',
        'score': '15',
        'tier': 'TRIPLE',
        'why': 'Deepest validation corridor, strongest translational claim, and best fit for response prediction plus persistence engineering.',
        'risk': 'Use as LBCL response framework, not universal CAR-T biology.',
        'task': 'response durability / native repertoire expansion / outcome prediction',
        'validation': ['GSE273170', 'GSE243973', 'GSE246342', 'GSE145007', 'GSE162975', 'GSE223655', 'GSE158676'],
        'figs': [
            ('Fig 1', 'Durable response versus early relapse', 'Contrast durable responders and early relapsers using longitudinal CAR-T sampling and tie the split to native cytotoxic repertoire expansion.'),
            ('Fig 2', 'Tonic signaling and binding-domain redesign', 'Show how construct-level architecture changes tonic signaling, antigen sensitivity, and fitness without changing the target antigen.'),
            ('Fig 3', 'Manufacturing QC and dose efficiency', 'Show how process changes affect yield, phenotype, and anti-lymphoma potency so the claim stays product-relevant.'),
            ('Fig 4', 'Persistence and translational rewiring', 'Connect translational control and lineage state to solid-tumor persistence and keep the interpretation at the state-transition level.'),
            ('Fig 5', 'Exhaustion rescue and memory skew', 'Use lineage-state results to explain why memory-like skew can improve potency.'),
            ('Fig 6', 'Armored payload localization', 'Show local payload delivery as a tumor-restricted activity strategy rather than blanket cytokine secretion.'),
            ('Fig 7', 'Validation corridor overview', 'Summarize the independent cohorts and orthogonal support that make this the only strong lane.'),
            ('Fig 8', 'Claim boundary schematic', 'Make explicit where the claim stops: construct-specific, disease-specific, and not universal CAR-T biology.'),
        ],
        'tabs': [
            ('Table 1', 'Anchor dataset map', 'List each CAR-T GEO record with score, best task, and a practical caveat.'),
            ('Table 2', 'Validation corridor', 'Pair the anchor with validation cohorts or orthogonal support and label the evidence depth.'),
            ('Table 3', 'Claim boundary table', 'Mark mouse-only, process-only, construct-specific, and antigen-low-dependent claims.'),
            ('Table 4', 'BD-fit table', 'Separate response prediction, potency/QC, and design guidance for partner-ready framing.'),
            ('Table 5', 'File map', 'List the companion pages that sit with the CAR-T story and what each one adds.'),
        ],
        'related': ['cell_therapy_ultra_analysis_2026_05_11.html', 'cell_therapy_hard_dossier_2026_05_11.html', 'car_t_bundle_2026_05_11.html'],
    },
    'CAR-NKT': {
        'rank': 2,
        'anchor': 'GSE270430',
        'title': 'CAR-NKT off-the-shelf platform',
        'verdict': 'moderate',
        'score': '14',
        'tier': 'DUAL',
        'why': 'Best platform pitch after CAR-T, centered on manufacturing simplification and myeloid-malignancy synergy.',
        'risk': 'Thin field; avoid disease-agnostic generalization.',
        'task': 'cell-therapy prioritization / off-the-shelf design',
        'validation': ['GSE223071', 'GSE235923'],
        'figs': [
            ('Fig 1', 'Allogeneic CAR-NKT architecture', 'Show how HSPC-engineered CAR-NKT is an off-the-shelf alternative to bespoke CAR-T manufacturing.'),
            ('Fig 2', 'Hypomethylating-agent synergy', 'Keep the focus on HMA co-treatment as the actual translation hook.'),
            ('Fig 3', 'Responder phenotype from the clinical series', 'Use the clinical readout to define the most relevant responder states.'),
            ('Fig 4', 'Cross-disease bridge', 'Treat AML/MDS and neuroblastoma as adjacent validation lanes, not one proof.'),
            ('Fig 5', 'Manufacturing and selection logic', 'Explain why the value is standardization and availability as much as antigen targeting.'),
            ('Fig 6', 'Safety and persistence balance', 'State the safety and persistence tradeoff so the platform pitch stays grounded.'),
        ],
        'tabs': [
            ('Table 1', 'Validation sources', 'List the primary CAR-NKT dataset and the adjacent clinical series in one place.'),
            ('Table 2', 'Platform caveats', 'Mark where the off-the-shelf story is direct evidence and where it is inference.'),
            ('Table 3', 'BD-fit table', 'Show the concrete partner value: manufacturing simplification and bridge data.'),
            ('Table 4', 'File map', 'List the companion pages that contextualize the CAR-NKT lane.'),
        ],
        'related': ['cell_therapy_ultra_analysis_2026_05_11.html', 'cell_therapy_hard_dossier_2026_05_11.html', 'cell_therapy_caption_atlas_2026_05_11.html'],
    },
    'TIL': {
        'rank': 3,
        'anchor': 'GSE303442',
        'title': 'TIL exhaustion rescue and potency',
        'verdict': 'moderate',
        'score': '14',
        'tier': 'DUAL',
        'why': 'Strongest non-CAR-T product-quality story; best framed as manufacturing rescue and potency control.',
        'risk': 'Orthogonal validation is supportive but not same-lane clinical proof.',
        'task': 'TIL expansion / exhaustion rescue / product potency',
        'validation': ['Organoids', 'Xenografts', 'TCR-seq diversity retention', 'Memory/costimulatory shift'],
        'figs': [
            ('Fig 1', 'RELB rescue of exhausted TILs', 'Show how RELB moves exhausted TILs toward a therapy-competent state.'),
            ('Fig 2', 'TCR diversity preservation', 'Use TCR-seq to show expansion without clonal collapse.'),
            ('Fig 3', 'Organoid killing and xenograft control', 'Move from expression to function with orthogonal functional readouts.'),
            ('Fig 4', 'Manufacturing QC logic', 'Translate the signal into release-assay logic and potency-positive states.'),
            ('Fig 5', 'Memory/costimulatory shift', 'Show the state distribution shift that makes the product more durable.'),
            ('Fig 6', 'Heterogeneity risk map', 'Flag where heterogeneity and exhaustion could break translation.'),
        ],
        'tabs': [
            ('Table 1', 'Rescue vs potency table', 'Pair proliferation, TCR diversity, killing, and persistence in one table.'),
            ('Table 2', 'Manufacturing risk table', 'Flag where heterogeneity or exhaustion could break the adoptive-cell-therapy path.'),
            ('Table 3', 'Orthogonal validation table', 'Keep organoid/xenograft support separate from the core GEO result.'),
            ('Table 4', 'File map', 'List the companion pages that contextualize the TIL lane.'),
        ],
        'related': ['cell_therapy_ultra_analysis_2026_05_11.html', 'cell_therapy_hard_dossier_2026_05_11.html', 'cell_therapy_caption_atlas_2026_05_11.html'],
    },
    'Blinatumomab': {
        'rank': 4,
        'anchor': 'GSE294241',
        'title': 'Blinatumomab scheduling and dysfunction',
        'verdict': 'moderate',
        'score': '13',
        'tier': 'SINGLE',
        'why': 'Most bounded lane; useful for schedule optimization and biomarker timing in B-ALL.',
        'risk': 'Disease- and age-specific; stay within B-ALL.',
        'task': 'response durability / T-cell dysfunction / dosing window',
        'validation': ['GSE196463', 'PMCID 7923532', 'B-ALL longitudinal response studies'],
        'figs': [
            ('Fig 1', 'T-cell dysfunction during continuous dosing', 'Show that persistent exposure erodes function and makes treatment duration biologically meaningful.'),
            ('Fig 2', 'Re-dosing and TCF7-high memory state', 'Show how re-dosing shifts the compartment toward a central-memory-like phenotype.'),
            ('Fig 3', 'Cytokine-risk versus response window', 'Place efficacy and cytokine risk on the same panel so the scheduling tradeoff is explicit.'),
            ('Fig 4', 'Clinical sample timing map', 'Lay out pre-, on-treatment, and re-dose timing for biomarker utility.'),
            ('Fig 5', 'Response-state biomarkers', 'Highlight the marker set most useful for who should be re-dosed and when.'),
            ('Fig 6', 'Safety window schematic', 'Make clear why cytokine-risk biology and efficacy should be interpreted separately.'),
        ],
        'tabs': [
            ('Table 1', 'Scheduling biomarkers', 'Convert response-state markers into a practical re-dose decision aid.'),
            ('Table 2', 'Safety boundary table', 'Separate cytokine-risk biology from efficacy claims.'),
            ('Table 3', 'B-ALL-specific caveat table', 'Mark the disease- and age-specific boundaries of the result.'),
            ('Table 4', 'File map', 'List the companion pages that contextualize the blinatumomab lane.'),
        ],
        'related': ['cell_therapy_ultra_analysis_2026_05_11.html', 'cell_therapy_hard_dossier_2026_05_11.html', 'cell_therapy_caption_atlas_2026_05_11.html'],
    },
}


def esc(x: object) -> str:
    return html.escape(str(x), quote=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open('r', encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def merge_maps(*paths: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    for path in paths:
        for row in read_tsv(path):
            acc = row.get('accession') or row.get('anchor')
            if not acc:
                continue
            prev = rows.get(acc, {})
            merged = dict(prev)
            for k, v in row.items():
                if v not in ('', None):
                    merged[k] = v
            rows[acc] = merged
    return rows


def caption_cards(items: list[tuple[str, str, str]], kind: str) -> str:
    return ''.join(
        f"<div class='cap'><div class='pill {kind.lower()}'>{esc(lbl)}</div><h3>{esc(title)}</h3><p>{esc(desc)}</p></div>"
        for lbl, title, desc in items
    )


def build_page(mod: str, rows: dict[str, dict[str, str]]) -> str:
    meta = MODS[mod]
    row = rows.get(meta['anchor'], {})
    file_map_rows = ''.join(
        f"<tr><td><b>{esc(name)}</b></td><td><code>{esc(str(path.relative_to(ROOT)))}</code></td><td><a href='{esc(url)}'>{esc(url)}</a></td></tr>"
        for name, path, url in FILES
    )
    related_rows = ''.join(
        f"<tr><td><code>{esc(rel)}</code></td><td><a href='{esc(LIVE_BASE + rel)}'>{esc(LIVE_BASE + rel)}</a></td></tr>"
        for rel in [f"{m}" for m in meta['related']]
    )
    validation = ''.join(f'<li>{esc(v)}</li>' for v in meta['validation'])
    accessions = ', '.join([meta['anchor']] + [r for r in row.get('validation_links', row.get('external_validation', '')).replace(';', ',').split(',') if r.strip()][:8])
    return f"""<!DOCTYPE html>
<html lang='en'>
<head>
<meta charset='UTF-8'/>
<meta name='viewport' content='width=device-width,initial-scale=1.0'/>
<title>{esc(mod)} independent dossier · cell therapy</title>
<style>
:root{{--bg:#050a11;--panel:#0f1726;--line:#22324b;--ink:#edf3fb;--muted:#9caec4;--gold:#ffd28a;--green:#35d39d;--blue:#7eb6ff;--red:#ff8a6b}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.6 Inter,"Noto Sans KR",sans-serif}} a{{color:inherit;text-decoration:none}} h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dc}} .hero{{padding:62px 34px 38px;background:radial-gradient(circle at top left,#182b4a 0%,#0a1220 58%,#050810 100%);border-bottom:1px solid var(--line)}} .hero-inner{{max-width:1500px;margin:0 auto}} .wrap{{max-width:1500px;margin:0 auto;padding:0 34px 80px}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}} h1{{font-size:62px;line-height:.95;margin:12px 0}} .lead{{max-width:1180px;color:#d4dfec;font:18px/1.68 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin-top:24px}} .stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.16);padding:13px 14px;border-radius:12px}} .stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}} .stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
section{{padding:30px 0;border-bottom:1px solid var(--line)}} h2{{font-size:38px;margin:0 0 6px}} h2 .num{{font:700 13px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;margin-right:10px}} .sub{{font:700 10px "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.12em;text-transform:uppercase;margin-bottom:16px}}
.grid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}} .grid2{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}} .cap{{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px 16px}} .cap h3{{margin:10px 0 6px}} .pill{{display:inline-block;padding:3px 9px;border-radius:999px;font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em}} .pill.figure{{background:rgba(53,211,157,.12);color:var(--green);border:1px solid rgba(53,211,157,.22)}} .pill.table{{background:rgba(126,182,255,.12);color:var(--blue);border:1px solid rgba(126,182,255,.22)}}
.box{{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 12px 12px 0;padding:14px 16px;margin:14px 0}} .box.good{{border-left-color:var(--green);background:#0a1a16}} .box.warn{{border-left-color:var(--red);background:#1d1310}} .small{{font-size:12px;color:var(--muted)}} .table{{width:100%;border-collapse:collapse;background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden}} .table th,.table td{{padding:10px 12px;border-bottom:1px solid rgba(255,255,255,.06);vertical-align:top;text-align:left}} .table th{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);text-transform:uppercase;letter-spacing:.08em}} code{{color:#dce7f4;background:#08101b;padding:2px 4px;border-radius:4px}}
@media(max-width:1100px){{.stats{{grid-template-columns:repeat(2,1fr)}} .grid,.grid2{{grid-template-columns:1fr}} h1{{font-size:40px}} .wrap{{padding:0 20px 60px}}}}
</style>
</head>
<body>
<header class='hero'>
  <div class='hero-inner'>
    <div class='kicker'>{esc(mod)} independent dossier · 2026-05-11</div>
    <h1>{esc(mod)}<br/><em>{esc(meta['title'])}</em></h1>
    <p class='lead'>This is the independent version: only one modality is shown, but the page still includes the strict claim boundary, validation ladder, expanded figure/table captions, and a file map so you can send it alone without losing context.</p>
    <div style='margin-top:18px;color:#9caec4;font:12px "JetBrains Mono",monospace;'>
      <a href='index.html' style='color:#35d39d'>Hub index</a> ·
      <a href='cell_therapy_megapack_2026_05_11.html' style='color:#ffd28a'>mega pack</a> ·
      <a href='cell_therapy_ultra_analysis_2026_05_11.html' style='color:#7eb6ff'>ultra analysis</a>
    </div>
    <div class='stats'>
      <div class='stat'><b>{meta['rank']}</b><span>Rank</span></div>
      <div class='stat'><b>{meta['score']}</b><span>Score</span></div>
      <div class='stat'><b>{len(meta['figs'])}</b><span>Figure captions</span></div>
      <div class='stat'><b>{len(meta['tabs'])}</b><span>Table captions</span></div>
      <div class='stat'><b>{meta['verdict']}</b><span>Verdict</span></div>
    </div>
  </div>
</header>
<div class='wrap'>
<section>
  <h2><span class='num'>01</span>Snapshot</h2>
  <div class='sub'>What this page says</div>
  <div class='box good'><b>Anchor accession:</b> {esc(meta['anchor'])}<br/><b>Best task:</b> {esc(meta['task'])}</div>
  <div class='grid2'>
    <div class='box good'><b>Why it matters.</b> {esc(meta['why'])}</div>
    <div class='box warn'><b>Claim boundary.</b> {esc(meta['risk'])}</div>
  </div>
  <div class='box'><b>Validation ladder.</b><ul>{validation}</ul></div>
</section>
<section>
  <h2><span class='num'>02</span>Analysis</h2>
  <div class='sub'>Independent interpretation</div>
  <div class='box'>
    <table class='table'>
      <tbody>
        <tr><th>Verdict</th><td>{esc(meta['verdict'])}</td></tr>
        <tr><th>Score</th><td>{esc(meta['score'])}</td></tr>
        <tr><th>Tier</th><td>{esc(meta['tier'])}</td></tr>
        <tr><th>Task</th><td>{esc(meta['task'])}</td></tr>
        <tr><th>Validation depth</th><td>{esc(', '.join(meta['validation']))}</td></tr>
      </tbody>
    </table>
  </div>
  <div class='box'><b>Summary.</b> {esc(meta['why'])}</div>
</section>
<section>
  <h2><span class='num'>03</span>Figure captions</h2>
  <div class='sub'>Expanded figure bank</div>
  <div class='grid'>{caption_cards(meta['figs'], 'figure')}</div>
</section>
<section>
  <h2><span class='num'>04</span>Table captions</h2>
  <div class='sub'>Expanded table bank</div>
  <div class='grid'>{caption_cards(meta['tabs'], 'table')}</div>
</section>
<section>
  <h2><span class='num'>05</span>File map</h2>
  <div class='sub'>Where this page sits relative to the rest of the bundle</div>
  <table class='table'>
    <thead><tr><th>Name</th><th>Local path</th><th>Live URL</th></tr></thead>
    <tbody>
      <tr><td><b>{esc(mod)} independent dossier</b></td><td><code>project/papers_hub_2026_05_04/{esc(mod.lower().replace('-','_'))}_independent_2026_05_11.html</code></td><td><a href='{esc(LIVE_BASE + mod.lower().replace('-', '_') + '_independent_2026_05_11.html')}'>{esc(LIVE_BASE + mod.lower().replace('-', '_') + '_independent_2026_05_11.html')}</a></td></tr>
      <tr><td><b>Cell therapy mega pack</b></td><td><code>project/papers_hub_2026_05_04/cell_therapy_megapack_2026_05_11.html</code></td><td><a href='{esc(LIVE_BASE + 'cell_therapy_megapack_2026_05_11.html')}'>{esc(LIVE_BASE + 'cell_therapy_megapack_2026_05_11.html')}</a></td></tr>
      <tr><td><b>Cell therapy ultra analysis</b></td><td><code>project/papers_hub_2026_05_04/cell_therapy_ultra_analysis_2026_05_11.html</code></td><td><a href='{esc(LIVE_BASE + 'cell_therapy_ultra_analysis_2026_05_11.html')}'>{esc(LIVE_BASE + 'cell_therapy_ultra_analysis_2026_05_11.html')}</a></td></tr>
      <tr><td><b>Cell therapy hard dossier</b></td><td><code>project/papers_hub_2026_05_04/cell_therapy_hard_dossier_2026_05_11.html</code></td><td><a href='{esc(LIVE_BASE + 'cell_therapy_hard_dossier_2026_05_11.html')}'>{esc(LIVE_BASE + 'cell_therapy_hard_dossier_2026_05_11.html')}</a></td></tr>
      <tr><td><b>Cell therapy caption atlas</b></td><td><code>project/papers_hub_2026_05_04/cell_therapy_caption_atlas_2026_05_11.html</code></td><td><a href='{esc(LIVE_BASE + 'cell_therapy_caption_atlas_2026_05_11.html')}'>{esc(LIVE_BASE + 'cell_therapy_caption_atlas_2026_05_11.html')}</a></td></tr>
      <tr><td><b>CAR-T bundle</b></td><td><code>project/papers_hub_2026_05_04/car_t_bundle_2026_05_11.html</code></td><td><a href='{esc(LIVE_BASE + 'car_t_bundle_2026_05_11.html')}'>{esc(LIVE_BASE + 'car_t_bundle_2026_05_11.html')}</a></td></tr>
    </tbody>
  </table>
</section>
<section>
  <h2><span class='num'>06</span>Related pages</h2>
  <div class='sub'>Companion pages for the same lane</div>
  <table class='table'>
    <thead><tr><th>File</th><th>URL</th></tr></thead>
    <tbody>{related_rows}</tbody>
  </table>
</section>
<section>
  <h2><span class='num'>07</span>Source files</h2>
  <div class='sub'>Where the evidence comes from</div>
  <table class='table'>
    <thead><tr><th>File</th><th>Contents</th></tr></thead>
    <tbody>
      <tr><td><code>{esc(str(PROMISING.relative_to(ROOT)))}</code></td><td>Score and prioritization master used to pick the anchor and task.</td></tr>
      <tr><td><code>{esc(str(VALMAT.relative_to(ROOT)))}</code></td><td>Validation matrix used to decide whether the lane is strong or moderate.</td></tr>
      <tr><td><code>{esc(str(CAPTION.relative_to(ROOT)))}</code></td><td>Caption atlas used as the base caption bank and modality grouping.</td></tr>
      <tr><td><code>{esc(str(HARD.relative_to(ROOT)))}</code></td><td>Hard dossier with claim boundary and evidence matrix.</td></tr>
      <tr><td><code>{esc(str(ULTRA.relative_to(ROOT)))}</code></td><td>Strict ranking layer across the four cell-therapy lanes.</td></tr>
      <tr><td><code>{esc(str(MEGA.relative_to(ROOT)))}</code></td><td>Expanded caption bank and file map for the full bundle.</td></tr>
    </tbody>
  </table>
</section>
</div>
</body>
</html>"""


def main() -> None:
    RESULT.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    LIVE_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    HUB.mkdir(parents=True, exist_ok=True)
    LIVE_HUB.mkdir(parents=True, exist_ok=True)
    rows = merge_maps(PROMISING, VALMAT, CAR_T, HARD, ULTRA)

    index_rows = []
    for mod in MODS:
        page_name = f"{mod.lower().replace('-', '_')}_independent_2026_05_11.html"
        html_text = build_page(mod, rows)
        html_path = HUB / page_name
        html_path.write_text(html_text, encoding='utf-8')
        page_tsv = RESULT / f"{mod.lower().replace('-', '_')}_independent.tsv"
        page_summary = RESULT / f"{mod.lower().replace('-', '_')}_independent_SUMMARY.md"
        meta = MODS[mod]
        with page_tsv.open('w', encoding='utf-8', newline='') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(['modality', 'anchor', 'verdict', 'score', 'tier', 'figures', 'tables', 'validation_depth'])
            w.writerow([mod, meta['anchor'], meta['verdict'], meta['score'], meta['tier'], len(meta['figs']), len(meta['tabs']), len(meta['validation'])])
        page_summary.write_text('\n'.join([
            f'# {mod} independent dossier',
            '',
            'Date: 2026-05-11',
            '',
            f'Independent page for {mod} with expanded figure and table captions plus a file map.',
            '',
            f'- Anchor: {meta["anchor"]}',
            f'- Verdict: {meta["verdict"]}',
            f'- Score: {meta["score"]}',
            f'- Tier: {meta["tier"]}',
            f'- Figure captions: {len(meta["figs"])}',
            f'- Table captions: {len(meta["tabs"])}',
            '',
            'Use:',
            f'- Send this when you only want the {mod} lane without the other modalities.',
        ]) + '\n', encoding='utf-8')
        for src in [html_path, page_tsv, page_summary]:
            shutil.copy2(src, LIVE_HUB / src.name)
            shutil.copy2(src, ASSET_DIR / src.name)
            shutil.copy2(src, LIVE_ASSET_DIR / src.name)
        index_rows.append((mod, page_name, LIVE_BASE + page_name))

    index = HUB / 'cell_therapy_independent_index_2026_05_11.html'
    index_html = ['<!DOCTYPE html>', '<html lang="en">', '<head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width,initial-scale=1.0"/><title>Cell Therapy Independent Index</title></head>', '<body style="margin:0;background:#050a11;color:#edf3fb;font:14px/1.6 Inter,sans-serif;">', '<div style="max-width:1200px;margin:0 auto;padding:40px 24px;">', '<h1 style="font-family:serif;">Cell Therapy Independent Index</h1>', '<p>Four stand-alone pages, each with its own analysis, figure bank, table bank, and file map.</p>', '<ul>']
    for mod, page_name, url in index_rows:
        index_html.append(f'<li><a href="{page_name}" style="color:#35d39d">{mod}</a> - <a href="{url}" style="color:#7eb6ff">live</a></li>')
    index_html.extend(['</ul>', '</div>', '</body>', '</html>'])
    index.write_text('\n'.join(index_html), encoding='utf-8')
    shutil.copy2(index, LIVE_HUB / index.name)
    shutil.copy2(index, ASSET_DIR / index.name)
    shutil.copy2(index, LIVE_ASSET_DIR / index.name)

    summary = RESULT / 'SUMMARY.md'
    summary.write_text('\n'.join([
        '# Cell Therapy Independent Series',
        '',
        'Date: 2026-05-11',
        '',
        'Standalone dossiers for CAR-T, CAR-NKT, TIL, and blinatumomab with expanded figure/table banks and file maps.',
        '',
        f'- Pages: {len(MODS)}',
        f'- Figure captions total: {sum(len(m["figs"]) for m in MODS.values())}',
        f'- Table captions total: {sum(len(m["tabs"]) for m in MODS.values())}',
        '',
        'Files:',
        '- car_t_independent_2026_05_11.html',
        '- car_nkt_independent_2026_05_11.html',
        '- til_independent_2026_05_11.html',
        '- blinatumomab_independent_2026_05_11.html',
    ]) + '\n', encoding='utf-8')

if __name__ == '__main__':
    main()
