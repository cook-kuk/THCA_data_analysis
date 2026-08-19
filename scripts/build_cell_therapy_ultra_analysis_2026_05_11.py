from __future__ import annotations

import csv
import html
import shutil
from pathlib import Path

ROOT = Path('/home/seungho/personal/THCA_data_analysis')
PROMISING = ROOT / 'project/results/therapy_spectrum_promising_master_2026_05_11/therapy_spectrum_promising_master.tsv'
VALMAT = ROOT / 'project/results/therapy_spectrum_validation_matrix_2026_05_11/therapy_spectrum_validation_matrix.tsv'
HARD = ROOT / 'project/results/cell_therapy_hard_dossier_2026_05_11/cell_therapy_hard_dossier.tsv'
RESULT = ROOT / 'project/results/cell_therapy_ultra_analysis_2026_05_11'
HUB = ROOT / 'project/papers_hub_2026_05_04'
LIVE_HUB = Path('/var/www/papers/papers_hub_2026_05_04')
ASSET_DIR = HUB / 'assets/cell_therapy_ultra_analysis_2026_05_11'
LIVE_ASSET_DIR = LIVE_HUB / 'assets/cell_therapy_ultra_analysis_2026_05_11'

MODALITIES = ['CAR-T', 'CAR-NKT', 'TIL', 'Blinatumomab']

MODALITY_MAP = {
    'CAR-T': {
        'anchor': 'GSE290722',
        'name': 'CAR-T response durability',
        'claim': 'Response prediction and persistence engineering in LBCL CAR-T',
        'use': 'Sendable now; strongest lane in the pack.',
        'risk': 'Disease-specific; should not be framed as universal CAR-T biology.',
    },
    'CAR-NKT': {
        'anchor': 'GSE270430',
        'name': 'CAR-NKT / off-the-shelf cell therapy',
        'claim': 'Off-the-shelf design and myeloid-malignancy synergy',
        'use': 'Sendable as a platform pitch, but not as a broad efficacy claim.',
        'risk': 'Thin field and cross-disease extrapolation risk.',
    },
    'TIL': {
        'anchor': 'GSE303442',
        'name': 'TIL adoptive therapy / RELB',
        'claim': 'Exhaustion rescue and product potency',
        'use': 'Strong for manufacturing/potency story.',
        'risk': 'Mostly orthogonal validation; not direct clinical endpoint proof.',
    },
    'Blinatumomab': {
        'anchor': 'GSE294241',
        'name': 'blinatumomab / therapy response',
        'claim': 'Dosing-window and T-cell dysfunction',
        'use': 'Useful for schedule optimization in B-ALL.',
        'risk': 'Disease- and age-specific; should stay within B-ALL.',
    },
}


def esc(x: object) -> str:
    return html.escape(str(x), quote=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open('r', encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def load_rows() -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    for src in [PROMISING, VALMAT, HARD]:
        for row in read_tsv(src):
            acc = row['accession']
            prev = rows.get(acc, {})
            merged = dict(prev)
            merged.update({k: v for k, v in row.items() if v not in ('', None)})
            rows[acc] = merged
    return rows


def build_html(rows: dict[str, dict[str, str]]) -> str:
    modalities = []
    for mod in MODALITIES:
        meta = MODALITY_MAP[mod]
        row = rows.get(meta['anchor'], {})
        verdict = row.get('verdict', 'n/a')
        tier = row.get('tier', row.get('class', 'n/a'))
        score = row.get('score_total', 'n/a')
        paper = row.get('paperability', 'n/a')
        val = row.get('validation', 'n/a')
        bd = row.get('bd', 'n/a')
        val_mat = rows.get(meta['anchor'], {})
        ext = rows.get(meta['anchor'], {}).get('validation_links') or rows.get(meta['anchor'], {}).get('external_validation', 'n/a')
        modalities.append({
            'mod': mod,
            'anchor': meta['anchor'],
            'name': meta['name'],
            'claim': meta['claim'],
            'use': meta['use'],
            'risk': meta['risk'],
            'verdict': verdict,
            'tier': tier,
            'score': score,
            'paper': paper,
            'val': val,
            'bd': bd,
            'ext': ext,
        })

    ranked = sorted(modalities, key=lambda x: (x['verdict'] != 'strong', -int(x['score']) if str(x['score']).isdigit() else 0, -int(x['paper']) if str(x['paper']).isdigit() else 0))

    rows_html = ''
    for m in ranked:
        rows_html += f"""
        <tr>
          <td><b>{esc(m['mod'])}</b><br/><span class='muted'>{esc(m['anchor'])}</span></td>
          <td>{esc(m['name'])}</td>
          <td><span class='pill {('strong' if m['verdict']=='strong' else 'moderate')}'>${{verdict}}</span></td>
          <td>{esc(m['score'])}</td>
          <td>{esc(m['tier'])}</td>
          <td>{esc(m['paper'])}</td>
          <td>{esc(m['val'])}</td>
          <td>{esc(m['bd'])}</td>
          <td>{esc(m['ext'])}</td>
        </tr>
        """.replace('${verdict}', m['verdict'])

    body_cards = ''
    for m in ranked:
        body_cards += f"""
        <div class='card'>
          <div class='top'><span class='k'>{esc(m['mod'])}</span><span class='pill {('strong' if m['verdict']=='strong' else 'moderate')}'> {esc(m['verdict'])} </span></div>
          <h3>{esc(m['anchor'])}</h3>
          <p class='lead2'><b>Claim:</b> {esc(m['claim'])}</p>
          <p><b>Use:</b> {esc(m['use'])}</p>
          <p><b>Risk:</b> {esc(m['risk'])}</p>
          <table class='mini'>
            <tbody>
              <tr><th>Score</th><td>{esc(m['score'])}</td></tr>
              <tr><th>Tier</th><td>{esc(m['tier'])}</td></tr>
              <tr><th>Paperability</th><td>{esc(m['paper'])}</td></tr>
              <tr><th>Validation</th><td>{esc(m['val'])}</td></tr>
              <tr><th>BD</th><td>{esc(m['bd'])}</td></tr>
              <tr><th>External / orthogonal</th><td>{esc(m['ext'])}</td></tr>
            </tbody>
          </table>
        </div>
        """

    return f"""<!DOCTYPE html>
<html lang='en'>
<head>
<meta charset='UTF-8'/>
<meta name='viewport' content='width=device-width,initial-scale=1.0'/>
<title>Cell Therapy Ultra Analysis · hard comparison</title>
<style>
:root{{--bg:#050b12;--panel:#10192a;--panel2:#0a1322;--line:#22324b;--ink:#edf3fb;--muted:#9bacc2;--gold:#ffd28a;--green:#35d39d;--blue:#7eb6ff;--red:#ff8a6b;--violet:#c79cff}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 Inter,"Noto Sans KR",sans-serif}} a{{color:inherit;text-decoration:none}} h1,h2,h3{{font-family:"Cormorant Garamond","Newsreader",serif;letter-spacing:-.02em;color:#fff7dc}}
.hero{{padding:60px 34px 36px;background:radial-gradient(circle at top left,#17304f 0%,#0b1220 58%,#050810 100%);border-bottom:1px solid var(--line)}} .wrap{{max-width:1520px;margin:0 auto;padding:0 34px 80px}} .hero-inner{{max-width:1520px;margin:0 auto}}
.kicker{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.24em;text-transform:uppercase}} h1{{font-size:60px;line-height:.96;margin:10px 0 12px}} .lead{{max-width:1160px;color:#d2dcec;font:18px/1.65 "Newsreader",serif}}
.stats{{display:grid;grid-template-columns:repeat(6,1fr);gap:10px;margin-top:24px}} .stat{{background:rgba(255,255,255,.04);border:1px solid rgba(255,210,138,.16);padding:13px 14px;border-radius:12px}} .stat b{{display:block;color:var(--gold);font:700 28px "Cormorant Garamond",serif;line-height:1}} .stat span{{display:block;color:var(--muted);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em;margin-top:4px}}
section{{padding:30px 0;border-bottom:1px solid var(--line)}} h2{{font-size:36px;margin:0 0 6px}} .sub{{font:700 10px "JetBrains Mono",monospace;color:var(--muted);letter-spacing:.12em;text-transform:uppercase;margin-bottom:16px}} .grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}} .card{{background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:18px}} .top{{display:flex;justify-content:space-between;align-items:center;margin-bottom:6px}} .k{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);letter-spacing:.2em;text-transform:uppercase}} .pill{{display:inline-block;padding:4px 10px;border-radius:999px;font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em}} .pill.strong{{background:rgba(53,211,157,.12);color:var(--green);border:1px solid rgba(53,211,157,.22)}} .pill.moderate{{background:rgba(126,182,255,.12);color:var(--blue);border:1px solid rgba(126,182,255,.22)}}
.table{{width:100%;border-collapse:collapse;background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden}} .table th,.table td{{padding:10px 12px;border-bottom:1px solid rgba(255,255,255,.06);vertical-align:top;text-align:left}} .table th{{font:700 11px "JetBrains Mono",monospace;color:var(--gold);text-transform:uppercase;letter-spacing:.08em}} .muted{{color:var(--muted)}} .mini{{width:100%;border-collapse:collapse;margin-top:10px}} .mini th,.mini td{{padding:6px 8px;border-bottom:1px solid rgba(255,255,255,.06);text-align:left;vertical-align:top}} .mini th{{width:170px;color:var(--gold);font:700 10px "JetBrains Mono",monospace;text-transform:uppercase;letter-spacing:.08em}} .mini td{{color:#dce5f2}} .note{{background:var(--panel2);border:1px solid var(--line);border-left:4px solid var(--gold);border-radius:0 12px 12px 0;padding:14px 16px;margin:12px 0}}
@media(max-width:1100px){{.stats{{grid-template-columns:repeat(3,1fr)}} .grid{{grid-template-columns:1fr}} h1{{font-size:42px}} .wrap{{padding:0 20px 60px}}}}
</style>
</head>
<body>
<header class='hero'>
  <div class='hero-inner'>
    <div class='kicker'>Ultra analysis · 2026-05-11</div>
    <h1>Cell Therapy Ultra Analysis<br/><em>hard comparison across the sendable cell-therapy lanes</em></h1>
    <p class='lead'>This page is the stricter layer on top of the hard dossier: it compares the four lanes on score, tier, validation depth, paperability, BD-readiness, and claim safety. The goal is not to broaden the story; the goal is to identify the strongest lane and the exact boundary where each other lane remains useful.</p>
    <div style='margin-top:18px;color:#9aabc0;font:12px "JetBrains Mono",monospace;'>
      <a href='index.html' style='color:#35d39d'>Hub index</a> ·
      <a href='cell_therapy_hard_dossier_2026_05_11.html' style='color:#ffd28a'>hard dossier</a> ·
      <a href='cell_therapy_caption_atlas_2026_05_11.html' style='color:#c79cff'>caption atlas</a> ·
      <a href='car_t_bundle_2026_05_11.html' style='color:#35d39d'>CAR-T bundle</a>
    </div>
    <div class='stats'>
      <div class='stat'><b>4</b><span>Modalities</span></div>
      <div class='stat'><b>1</b><span>Strong lane</span></div>
      <div class='stat'><b>3</b><span>Moderate lanes</span></div>
      <div class='stat'><b>7</b><span>CAR-T validation links</span></div>
      <div class='stat'><b>4</b><span>Evidence matrices</span></div>
      <div class='stat'><b>Sendable</b><span>Execution pack</span></div>
    </div>
  </div>
</header>
<div class='wrap'>
<section>
  <h2><span class='num'>01</span>Hard ranking</h2>
  <div class='sub'>What survives a strict comparison</div>
  <div class='note'><b>Bottom line.</b> CAR-T remains the only lane that clears the full `strong` bar. CAR-NKT is the cleanest platform story after that. TIL is mechanistically strong but still mostly a product-quality lane. Blinatumomab is useful as a schedule/dysfunction decision-support lane, not as a broad platform story.</div>
  <table class='table'>
    <thead>
      <tr><th>Rank</th><th>Modality</th><th>Anchor</th><th>Verdict</th><th>Score</th><th>Tier</th><th>Paper</th><th>Validation</th><th>BD</th><th>External / orthogonal support</th></tr>
    </thead>
    <tbody>{rows_html}</tbody>
  </table>
</section>
<section>
  <h2><span class='num'>02</span>What each lane can honestly claim</h2>
  <div class='sub'>Claim-safe framing only</div>
  <div class='grid'>{body_cards}</div>
</section>
<section>
  <h2><span class='num'>03</span>Interpretation</h2>
  <div class='sub'>Why the ranking looks this way</div>
  <div class='note'>
    <b>CAR-T:</b> has the deepest validation corridor and the highest combined paper/validation/BD scores, so it is the only lane that can support both a mechanistic and translational story without overreach.<br/>
    <b>CAR-NKT:</b> has the strongest platform pitch after CAR-T because the off-the-shelf story is operationally clear, but the validation corridor is shorter.<br/>
    <b>TIL:</b> has the strongest productization story of the non-CAR-T lanes, but the evidence mostly supports potency rescue and manufacturing logic rather than direct clinical claim.<br/>
    <b>Blinatumomab:</b> is the most bounded lane; it is useful because it informs dosing-window and T-cell dysfunction, not because it supports a broad immunotherapy platform.
  </div>
  <div class='note'>
    <b>Claim-safety rule.</b> The more the evidence leans on longitudinal patient data and independent cohorts, the more you can say `response prediction`. The more it leans on orthogonal functional support or manufacturing readouts, the more you must stay in `potency / design / decision-support` territory.
  </div>
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

    rows = load_rows()
    html_text = build_html(rows)
    html_path = HUB / 'cell_therapy_ultra_analysis_2026_05_11.html'
    html_path.write_text(html_text, encoding='utf-8')

    tsv = RESULT / 'cell_therapy_ultra_analysis.tsv'
    ranked_modalities = sorted(
        MODALITIES,
        key=lambda mod: (
            MODALITY_MAP[mod]['anchor'] != 'GSE290722',
            -(int(rows.get(MODALITY_MAP[mod]['anchor'], {}).get('score_total', '0')) if rows.get(MODALITY_MAP[mod]['anchor'], {}).get('score_total', '0').isdigit() else 0),
            -(int(rows.get(MODALITY_MAP[mod]['anchor'], {}).get('paperability', '0')) if rows.get(MODALITY_MAP[mod]['anchor'], {}).get('paperability', '0').isdigit() else 0),
            -(int(rows.get(MODALITY_MAP[mod]['anchor'], {}).get('validation', '0')) if rows.get(MODALITY_MAP[mod]['anchor'], {}).get('validation', '0').isdigit() else 0),
            -(int(rows.get(MODALITY_MAP[mod]['anchor'], {}).get('bd', '0')) if rows.get(MODALITY_MAP[mod]['anchor'], {}).get('bd', '0').isdigit() else 0),
        ),
    )
    with tsv.open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['rank', 'modality', 'anchor', 'verdict', 'score_total', 'tier', 'paperability', 'validation', 'bd', 'external_support', 'claim', 'risk'])
        for i, m in enumerate(ranked_modalities, start=1):
            meta = MODALITY_MAP[m]
            row = rows.get(meta['anchor'], {})
            support = row.get('validation_links') or row.get('external_validation', '')
            w.writerow([i, m, meta['anchor'], row.get('verdict', ''), row.get('score_total', ''), row.get('tier', row.get('class', '')), row.get('paperability', ''), row.get('validation', ''), row.get('bd', ''), support, meta['claim'], meta['risk']])

    summary = RESULT / 'SUMMARY.md'
    summary.write_text('\n'.join([
        '# Cell Therapy Ultra Analysis',
        '',
        'Date: 2026-05-11',
        '',
        'Strict comparison layer on top of the hard dossier: ranks the four cell-therapy lanes by score, tier, validation depth, paperability, BD-readiness, and claim safety.',
        '',
        'Snapshot:',
        '- Modalities: 4',
        '- Strong lanes: 1',
        '- Moderate lanes: 3',
        '- Validation corridor: CAR-T has 7 external links; CAR-NKT 2; TIL 4 orthogonal supports; blinatumomab 3 bounded supports',
        '',
        'Use:',
        '- Send this when you need the strictest cell-therapy comparison available in the current bundle.',
    ]) + '\n', encoding='utf-8')

    for src in [html_path, tsv, summary]:
        shutil.copy2(src, LIVE_HUB / src.name)
        shutil.copy2(src, ASSET_DIR / src.name)
        shutil.copy2(src, LIVE_ASSET_DIR / src.name)

if __name__ == '__main__':
    main()
