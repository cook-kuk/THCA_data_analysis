from __future__ import annotations

import csv
import html
import shutil
from pathlib import Path

ROOT = Path('/home/seungho/personal/THCA_data_analysis')
RESULT = ROOT / 'project/results/cell_therapy_megapack_2026_05_11'
HUB = ROOT / 'project/papers_hub_2026_05_04'
LIVE_HUB = Path('/var/www/papers/papers_hub_2026_05_04')
ASSET_DIR = HUB / 'assets/cell_therapy_megapack_2026_05_11'
LIVE_ASSET_DIR = LIVE_HUB / 'assets/cell_therapy_megapack_2026_05_11'

HARD = ROOT / 'project/results/cell_therapy_hard_dossier_2026_05_11/cell_therapy_hard_dossier.tsv'
ULTRA = ROOT / 'project/results/cell_therapy_ultra_analysis_2026_05_11/cell_therapy_ultra_analysis.tsv'
CAPTION = ROOT / 'project/results/cell_therapy_caption_atlas_2026_05_11/cell_therapy_caption_atlas.tsv'
CAR_T = ROOT / 'project/results/car_t_bundle_2026_05_11/car_t_bundle.tsv'
PROMISING = ROOT / 'project/results/therapy_spectrum_promising_master_2026_05_11/therapy_spectrum_promising_master.tsv'
VALMAT = ROOT / 'project/results/therapy_spectrum_validation_matrix_2026_05_11/therapy_spectrum_validation_matrix.tsv'

LIVE_BASE = 'http://40.82.129.113/papers_hub_2026_05_04/'

FILE_MAP = [
    ('Hub index', 'project/papers_hub_2026_05_04/index.html', LIVE_BASE + 'index.html'),
    ('Cell therapy caption atlas', 'project/papers_hub_2026_05_04/cell_therapy_caption_atlas_2026_05_11.html', LIVE_BASE + 'cell_therapy_caption_atlas_2026_05_11.html'),
    ('Cell therapy hard dossier', 'project/papers_hub_2026_05_04/cell_therapy_hard_dossier_2026_05_11.html', LIVE_BASE + 'cell_therapy_hard_dossier_2026_05_11.html'),
    ('Cell therapy ultra analysis', 'project/papers_hub_2026_05_04/cell_therapy_ultra_analysis_2026_05_11.html', LIVE_BASE + 'cell_therapy_ultra_analysis_2026_05_11.html'),
    ('CAR-T bundle', 'project/papers_hub_2026_05_04/car_t_bundle_2026_05_11.html', LIVE_BASE + 'car_t_bundle_2026_05_11.html'),
    ('Therapy topic report', 'project/papers_hub_2026_05_04/therapy_spectrum_topic_report_2026_05_11.html', LIVE_BASE + 'therapy_spectrum_topic_report_2026_05_11.html'),
    ('Therapy validation matrix', 'project/papers_hub_2026_05_04/therapy_spectrum_validation_matrix_2026_05_11.html', LIVE_BASE + 'therapy_spectrum_validation_matrix_2026_05_11.html'),
    ('Therapy promising master', 'project/papers_hub_2026_05_04/therapy_spectrum_promising_master_2026_05_11.html', LIVE_BASE + 'therapy_spectrum_promising_master_2026_05_11.html'),
]

MODS = {
    'CAR-T': {
        'anchor': 'GSE290722',
        'rank': '1',
        'title': 'Response durability and engineering',
        'figs': [
            ('Fig 1', 'Durable response vs early relapse', 'Use longitudinal sampling to separate durable responders from early relapsers and anchor the rest of the panel around the native cytotoxic T-cell repertoire shift.'),
            ('Fig 2', 'Tonic signaling tradeoff map', 'Show how construct-level design affects tonic signaling, antigen sensitivity, and fitness, then point to the design knobs that can be tuned without changing the target.'),
            ('Fig 3', 'Manufacturing QC ladder', 'Make the process side explicit: expansion, phenotype preservation, and potency should be shown together so the claim does not overreach into clinical efficacy.'),
            ('Fig 4', 'Persistence and translational rewiring', 'Tie translational control to solid-tumor persistence and keep the interpretation at the level of state transition rather than universal biology.'),
            ('Fig 5', 'Exhaustion rescue and memory skew', 'Show how lineage-state and transcriptional brakes can move the product away from dysfunctional exhaustion and toward a more durable memory-like state.'),
            ('Fig 6', 'Armored payload localization', 'Use tumor-restricted payload delivery as the design principle and explicitly separate local activity from systemic cytokine exposure.'),
            ('Fig 7', 'Validation corridor overview', 'Summarize the independent cohorts and orthogonal support so the reader sees why this lane clears the strong threshold.'),
        ],
        'tabs': [
            ('Table 1', 'Anchor dataset map', 'Map each CAR-T accession to its best task, score, and practical caveat so the bundle does not collapse into “many CAR-T papers.”'),
            ('Table 2', 'Validation corridor table', 'Match each anchor to the external cohorts or orthogonal support and label whether the evidence is response-level, mechanistic, or process-level.'),
            ('Table 3', 'Claim boundary table', 'State what must not be claimed: universal CAR-T biology, mouse-only process, construct-only effects, or antigen-low dependence.'),
            ('Table 4', 'Translation-to-BD table', 'List which outputs are useful for response prediction, which for potency/QC, and which are only useful as design guidance.'),
        ],
        'analysis': [
            'Deepest validation corridor in the bundle',
            'Only lane that supports a strong translational claim',
            'Best for response prediction plus persistence engineering',
            'Useful as the central anchor for all sendable cell-therapy narratives',
        ],
    },
    'CAR-NKT': {
        'anchor': 'GSE270430',
        'rank': '2',
        'title': 'Off-the-shelf platform pitch',
        'figs': [
            ('Fig 1', 'Allogeneic CAR-NKT design', 'Show how HSPC-engineered CAR-NKT cells are positioned as an off-the-shelf alternative to bespoke CAR-T production.'),
            ('Fig 2', 'HMA synergy and anti-malignancy activity', 'Keep the focus on how hypomethylating-agent synergy converts a construct into a translational story.'),
            ('Fig 3', 'Responder phenotype from the clinical series', 'Use the clinical readout to show what responder-like behavior looks like and where the evidence is still thin.'),
            ('Fig 4', 'Cross-disease bridge', 'Make the AML/MDS to neuroblastoma bridge explicit so the platform pitch stays bounded.'),
            ('Fig 5', 'Manufacturing and selection logic', 'Explain why this lane is about standardization and availability as much as antigen recognition.'),
        ],
        'tabs': [
            ('Table 1', 'Validation source map', 'List the primary CAR-NKT dataset and adjacent clinical series in one place so the bridge and the gap are visible.'),
            ('Table 2', 'Platform caveat table', 'Mark where the allogeneic/off-the-shelf claim is directly supported and where it remains an inference.'),
            ('Table 3', 'BD-fit table', 'Specify what a partner would actually buy here: manufacturing simplification, rescue logic, and disease-specific bridge data.'),
        ],
        'analysis': [
            'Best platform pitch after CAR-T',
            'Off-the-shelf positioning is the core value',
            'Needs careful disease-specific framing',
            'Validation corridor is shorter than CAR-T but still coherent',
        ],
    },
    'TIL': {
        'anchor': 'GSE303442',
        'rank': '3',
        'title': 'Product potency and rescue',
        'figs': [
            ('Fig 1', 'RELB exhaustion rescue', 'Show how RELB shifts exhausted TILs toward a more therapy-competent state rather than only reporting a differential-expression result.'),
            ('Fig 2', 'TCR diversity preservation', 'Use TCR-seq to show that expansion does not have to mean clonal collapse, which is key for product-level trust.'),
            ('Fig 3', 'Organoid killing and xenograft control', 'Connect the molecular rescue to a functional readout so the reader sees the bridge from expression to function.'),
            ('Fig 4', 'Manufacturing QC logic', 'Make the release-assay implication explicit: which states should be considered potency-positive and which should be treated as liabilities.'),
            ('Fig 5', 'Memory/costimulatory shift', 'Translate the mechanistic signal into a practical manufacturing objective by showing how the cell state distribution moves.'),
        ],
        'tabs': [
            ('Table 1', 'Rescue vs potency table', 'Pair proliferation, TCR diversity, killing, and persistence into a single product-quality table.'),
            ('Table 2', 'Manufacturing risk table', 'Flag where heterogeneity, exhaustion, or expansion drift can break the adoptive-cell-therapy path.'),
            ('Table 3', 'Orthogonal validation table', 'Keep the organoid/xenograft support separate from the core GEO result so the evidence stack stays honest.'),
        ],
        'analysis': [
            'Strongest non-CAR-T product-quality story',
            'Best read as potency rescue rather than outcome proof',
            'Orthogonal support is useful but not fully same-lane validation',
            'Good partner pitch if manufacturing is the target',
        ],
    },
    'Blinatumomab': {
        'anchor': 'GSE294241',
        'rank': '4',
        'title': 'Dosing-window decision support',
        'figs': [
            ('Fig 1', 'T-cell dysfunction during continuous dosing', 'Show that duration of exposure is biologically meaningful because persistent treatment can erode T-cell function.'),
            ('Fig 2', 'Re-dosing and TCF7-high memory state', 'Use the re-dosing profile to show how the compartment shifts toward a more central-memory phenotype.'),
            ('Fig 3', 'Cytokine-risk versus response window', 'Put efficacy and cytokine risk on the same plot so the scheduling tradeoff becomes concrete.'),
            ('Fig 4', 'Clinical sample timing map', 'Lay out pre-, on-treatment, and re-dose timing so biomarker collection can be interpreted correctly.'),
            ('Fig 5', 'Response-state biomarker panel', 'Show the small set of markers that are most useful for deciding who should be re-dosed and when.'),
        ],
        'tabs': [
            ('Table 1', 'Scheduling biomarkers table', 'Convert response-state markers into a practical redose decision aid.'),
            ('Table 2', 'Safety boundary table', 'Keep cytokine-risk biology separate from efficacy claims so the page stays claim-safe.'),
            ('Table 3', 'B-ALL-specific caveat table', 'Mark that the result is disease- and age-specific and should not be generalized beyond the original setting.'),
        ],
        'analysis': [
            'Most bounded lane',
            'Best interpreted as schedule optimization',
            'Useful for biomarker timing and dysfunction tracking',
            'Do not generalize beyond B-ALL',
        ],
    },
}


def esc(x: object) -> str:
    return html.escape(str(x), quote=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open('r', encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def load_maps() -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    prom = {row['accession']: row for row in read_tsv(PROMISING)}
    val = {row['accession']: row for row in read_tsv(VALMAT)}
    return prom, val


def rows_for_anchor(prom: dict[str, dict[str, str]], val: dict[str, dict[str, str]], anchor: str) -> dict[str, str]:
    row = {}
    if anchor in prom:
        row.update(prom[anchor])
    if anchor in val:
        row.update({k: v for k, v in val[anchor].items() if v not in ('', None)})
    return row


def build_caption_block(items: list[tuple[str, str, str]], kind: str) -> str:
    return ''.join(
        f"<div class='cap'><div class='pill {kind.lower()}'>{esc(lbl)}</div><h3>{esc(title)}</h3><p>{esc(desc)}</p></div>"
        for lbl, title, desc in items
    )


def build_section(modality: str, prom: dict[str, dict[str, str]], val: dict[str, dict[str, str]]) -> str:
    meta = MODS[modality]
    row = rows_for_anchor(prom, val, meta['anchor'])
    return f"""
<section id='{modality.lower().replace('-','_')}'>
  <h2><span class='num'>{esc(meta['rank'])}</span>{esc(modality)}</h2>
  <div class='sub'>{esc(meta['anchor'])} | {esc(row.get('verdict','n/a'))} | score {esc(row.get('score_total','n/a'))} | tier {esc(row.get('tier', row.get('class','n/a')))}</div>
  <div class='grid2'>
    <div class='box good'><b>Why it matters.</b> {esc(row.get('why', row.get('title', meta['title'])))} </div>
    <div class='box warn'><b>Claim boundary.</b> {esc(row.get('caveat', meta['analysis'][-1]))}</div>
  </div>
  <div class='grid2'>
    <div class='box'><b>Best task.</b> {esc(row.get('best_task', meta['analysis'][0]))}</div>
    <div class='box'><b>Validation depth.</b> {esc(row.get('validation', row.get('external_validation', 'n/a')))}</div>
  </div>
  <div class='box'><b>Analysis notes.</b><ul>{''.join(f'<li>{esc(x)}</li>' for x in meta['analysis'])}</ul></div>
  <div class='subhead'>Figure captions</div>
  <div class='grid'>{build_caption_block(meta['figs'], 'figure')}</div>
  <div class='subhead'>Table captions</div>
  <div class='grid'>{build_caption_block(meta['tabs'], 'table')}</div>
</section>
"""


def build_html(prom: dict[str, dict[str, str]], val: dict[str, dict[str, str]]) -> str:
    counts = {
        'modalities': len(MODS),
        'figs': sum(len(v['figs']) for v in MODS.values()),
        'tabs': sum(len(v['tabs']) for v in MODS.values()),
    }
    cards = ''.join(
        f"<a class='topic' href='#{mod.lower().replace('-','_')}'><div class='k'>{esc(mod)}</div><h3>{esc(meta['anchor'])}</h3><p>{esc(meta['title'])}</p><p class='small'>figs {len(meta['figs'])} · tables {len(meta['tabs'])}</p></a>"
        for mod, meta in MODS.items()
    )
    files = ''.join(
        f"<tr><td><b>{esc(name)}</b></td><td><code>{esc(rel)}</code></td><td><a href='{esc(url)}'>{esc(url)}</a></td></tr>"
        for name, rel, url in FILE_MAP
    )
    sections = ''.join(build_section(mod, prom, val) for mod in MODS)
    return f"""<!DOCTYPE html>
<html lang='en'>
<head>
<meta charset='UTF-8'/>
<meta name='viewport' content='width=device-width,initial-scale=1.0'/>
<title>Cell Therapy Mega Pack · figures, tables, file map</title>
<style>
:root{{--bg:#050a11;--panel:#0f1726;--line:#22324b;--ink:#edf3fb;--muted:#9caec4;--gold:#ffd28a;--green:#35d39d;--blue:#7eb6ff;--red:#ff8a6b;--violet:#c79cff}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.6 Inter,"Noto Sans KR",sans-serif}} a{{color:inherit;text-decoration:none}} h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dc}} .hero{{padding:64px 34px 40px;background:radial-gradient(circle at top left,#182b4a 0%,#0a1220 58%,#050810 100%);border-bottom:1px solid var(--line)}} .hero-inner{{max-width:1540px;margin:0 auto}} .wrap{{max-width:1540px;margin:0 auto;padding:0 34px 80px}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}} h1{{font-size:64px;line-height:.95;margin:12px 0}} .lead{{max-width:1180px;color:#d4dfec;font:18px/1.68 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(6,1fr);gap:10px;margin-top:24px}} .stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.16);padding:13px 14px;border-radius:12px}} .stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}} .stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
section{{padding:30px 0;border-bottom:1px solid var(--line)}} h2{{font-size:38px;margin:0 0 6px}} h2 .num{{font:700 13px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.18em;margin-right:10px}} .sub{{font:700 10px "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.12em;text-transform:uppercase;margin-bottom:16px}}
.grid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}} .grid2{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}} .topic{{display:block;background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:16px 18px;min-height:122px}} .topic:hover{{border-color:var(--gold)}} .topic p{{margin:8px 0 0}} .k{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.2em;text-transform:uppercase}} .small{{font-size:12px;color:var(--muted)}}
.box{{background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 12px 12px 0;padding:14px 16px;margin:14px 0}} .box.good{{border-left-color:var(--green);background:#0a1a16}} .box.warn{{border-left-color:var(--red);background:#1d1310}} .subhead{{font:700 12px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.12em;margin:26px 0 12px;text-transform:uppercase}}
.cap{{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:14px 16px}} .cap h3{{margin:10px 0 6px}} .pill{{display:inline-block;padding:3px 9px;border-radius:999px;font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em}} .pill.figure{{background:rgba(53,211,157,.12);color:var(--green);border:1px solid rgba(53,211,157,.22)}} .pill.table{{background:rgba(126,182,255,.12);color:var(--blue);border:1px solid rgba(126,182,255,.22)}}
.table{{width:100%;border-collapse:collapse;background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden}} .table th,.table td{{padding:10px 12px;border-bottom:1px solid rgba(255,255,255,.06);vertical-align:top;text-align:left}} .table th{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);text-transform:uppercase;letter-spacing:.08em}} code{{color:#dfe8f4;background:#08101b;padding:2px 4px;border-radius:4px}}
@media(max-width:1100px){{.stats{{grid-template-columns:repeat(3,1fr)}} .grid,.grid2{{grid-template-columns:1fr}} h1{{font-size:42px}} .wrap{{padding:0 20px 60px}}}}
</style>
</head>
<body>
<header class='hero'>
  <div class='hero-inner'>
    <div class='kicker'>Cell therapy mega pack · 2026-05-11</div>
    <h1>Cell Therapy Mega Pack<br/><em>figure/table atlas + file map</em></h1>
    <p class='lead'>This is the expanded sendable pack. It puts the figure and table explanations in one place, then lists exactly where each companion file lives and what its live URL is. Use this when you need a full cell-therapy brief rather than a short bundle.</p>
    <div style='margin-top:18px;color:#9caec4;font:12px "JetBrains Mono",monospace;'>
      <a href='index.html' style='color:#35d39d'>Hub index</a> ·
      <a href='cell_therapy_caption_atlas_2026_05_11.html' style='color:#c79cff'>caption atlas</a> ·
      <a href='cell_therapy_hard_dossier_2026_05_11.html' style='color:#ffd28a'>hard dossier</a> ·
      <a href='cell_therapy_ultra_analysis_2026_05_11.html' style='color:#7eb6ff'>ultra analysis</a> ·
      <a href='car_t_bundle_2026_05_11.html' style='color:#35d39d'>CAR-T bundle</a>
    </div>
    <div class='stats'>
      <div class='stat'><b>{counts['modalities']}</b><span>Modalities</span></div>
      <div class='stat'><b>{counts['figs']}</b><span>Figure captions</span></div>
      <div class='stat'><b>{counts['tabs']}</b><span>Table captions</span></div>
      <div class='stat'><b>4</b><span>File map pages</span></div>
      <div class='stat'><b>GEO</b><span>Public anchors</span></div>
      <div class='stat'><b>Sendable</b><span>Expanded pack</span></div>
    </div>
  </div>
</header>
<div class='wrap'>
<section>
  <h2><span class='num'>01</span>File map</h2>
  <div class='sub'>Where everything lives</div>
  <table class='table'>
    <thead><tr><th>Name</th><th>Local path</th><th>Live URL</th></tr></thead>
    <tbody>{files}</tbody>
  </table>
</section>
<section>
  <h2><span class='num'>02</span>Ranking</h2>
  <div class='sub'>Strict lane order from the ultra analysis</div>
  <div class='grid'>{cards}</div>
  <div class='box good'><b>Bottom line:</b> CAR-T is the only `strong` lane. TIL is the strongest `moderate` lane. CAR-NKT stays a platform pitch, and blinatumomab remains a schedule-specific decision-support lane.</div>
</section>
<section>
  <h2><span class='num'>03</span>CAR-T</h2>
  <div class='sub'>GSE290722 plus adjacent CAR-T bundle accessions</div>
  <div class='box good'><b>Anchor accessions:</b> GSE290722, GSE313971, GSE308384, GSE236468, GSE304795, GSE312384, GSE292859, GSE284026</div>
  <div class='grid'>{build_caption_block(MODS['CAR-T']['figs'], 'figure')}</div>
  <div class='subhead'>Table captions</div>
  <div class='grid'>{build_caption_block(MODS['CAR-T']['tabs'], 'table')}</div>
</section>
<section>
  <h2><span class='num'>04</span>CAR-NKT</h2>
  <div class='sub'>Off-the-shelf platform pitch</div>
  <div class='box good'><b>Anchor accessions:</b> GSE270430, GSE223071</div>
  <div class='grid'>{build_caption_block(MODS['CAR-NKT']['figs'], 'figure')}</div>
  <div class='subhead'>Table captions</div>
  <div class='grid'>{build_caption_block(MODS['CAR-NKT']['tabs'], 'table')}</div>
</section>
<section>
  <h2><span class='num'>05</span>TIL</h2>
  <div class='sub'>Exhaustion rescue and product potency</div>
  <div class='box good'><b>Anchor accessions:</b> GSE303442</div>
  <div class='grid'>{build_caption_block(MODS['TIL']['figs'], 'figure')}</div>
  <div class='subhead'>Table captions</div>
  <div class='grid'>{build_caption_block(MODS['TIL']['tabs'], 'table')}</div>
</section>
<section>
  <h2><span class='num'>06</span>Blinatumomab</h2>
  <div class='sub'>Dose-window and T-cell dysfunction</div>
  <div class='box good'><b>Anchor accessions:</b> GSE294241, GSE292621</div>
  <div class='grid'>{build_caption_block(MODS['Blinatumomab']['figs'], 'figure')}</div>
  <div class='subhead'>Table captions</div>
  <div class='grid'>{build_caption_block(MODS['Blinatumomab']['tabs'], 'table')}</div>
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
    prom, val = load_maps()
    html_text = build_html(prom, val)
    html_path = HUB / 'cell_therapy_megapack_2026_05_11.html'
    html_path.write_text(html_text, encoding='utf-8')
    # companion result files for traceability
    summary = RESULT / 'SUMMARY.md'
    summary.write_text('\n'.join([
        '# Cell Therapy Mega Pack',
        '',
        'Date: 2026-05-11',
        '',
        'Expanded figure/table atlas with a file map for the cell-therapy bundle.',
        '',
        'Snapshot:',
        f'- Modalities: {len(MODS)}',
        f'- Figure captions: {sum(len(v["figs"]) for v in MODS.values())}',
        f'- Table captions: {sum(len(v["tabs"]) for v in MODS.values())}',
        f'- File map entries: {len(FILE_MAP)}',
        '',
        'Use:',
        '- Send this when you want the fullest cell-therapy caption bank and an exact map of where the companion files live.',
    ]) + '\n', encoding='utf-8')
    tsv = RESULT / 'cell_therapy_megapack.tsv'
    with tsv.open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['modality', 'rank', 'anchor', 'figures', 'tables', 'file_count'])
        for mod, meta in MODS.items():
            w.writerow([mod, meta['rank'], meta['anchor'], len(meta['figs']), len(meta['tabs']), len(FILE_MAP)])
    for src in [html_path, tsv, summary]:
        shutil.copy2(src, LIVE_HUB / src.name)
        shutil.copy2(src, ASSET_DIR / src.name)
        shutil.copy2(src, LIVE_ASSET_DIR / src.name)

if __name__ == '__main__':
    main()
