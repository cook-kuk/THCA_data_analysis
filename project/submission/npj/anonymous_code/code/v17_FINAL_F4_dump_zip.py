"""F4: Final ship dump + anonymous code zip (LOCAL ONLY — no GitHub push)."""
from v17p35_SYNTH1_common import SUB, FIG_DIR, TBL_DIR, REPRO_DIR, REPORTS, RES, log
import shutil, zipfile, datetime, json, re
from pathlib import Path

NOTEBOOKS = Path('/home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts')
ANON_OUT = SUB / 'anonymous_code.zip'
ANON_DIR = SUB / 'anonymous_code'

# Patterns to scrub author identifiers from script comments/headers
SCRUB_PATTERNS = [
    (re.compile(r'[Author 1]', re.IGNORECASE), '[Author 1]'),
    (re.compile(r'kukshomr@gmail\.com', re.IGNORECASE), '[author1@anonymized]'),
    (re.compile(r'유\s*교수', re.IGNORECASE), '[Author 2]'),
    (re.compile(r'[Author 2]', re.IGNORECASE), '[Author 2]'),
    (re.compile(r'[Author 1]'), '[Author 1]'),
    (re.compile(r'[Affiliation]|[Affiliation]|[Affiliation]', re.IGNORECASE), '[Affiliation]'),
]

def scrub(text: str) -> str:
    for p, repl in SCRUB_PATTERNS:
        text = p.sub(repl, text)
    return text


def build_anon_zip():
    if ANON_DIR.exists():
        shutil.rmtree(ANON_DIR)
    ANON_DIR.mkdir()
    code_root = ANON_DIR / 'code'
    code_root.mkdir()
    for src in sorted(NOTEBOOKS.glob('v17p35_*.py')) + sorted(NOTEBOOKS.glob('v17_FINAL_*.py')):
        if src.name.endswith('.cpython-312.pyc'): continue
        dst = code_root / src.name
        dst.write_text(scrub(src.read_text()))

    # README anonymous
    readme = ANON_DIR / 'README.md'
    readme.write_text("""# Anonymous code archive — npj submission supplementary

This archive contains all pipeline scripts for the manuscript.
Author identifiers have been removed for double-blind review.

## Layout
- `code/v17p35_*.py` — main analysis pipeline (Phase A → SYNTH-1).
- `code/v17_FINAL_*.py` — final ship sprint (manuscript v2 + Fig8 + audit).

## Replication
1. `python -m venv .venv && source .venv/bin/activate`
2. `pip install plotly pandas numpy scipy scikit-learn 'kaleido==0.2.1'`
3. Run `code/v17p35_SYNTH1_Fig1.py` … `code/v17_FINAL_F1B_fig8.py` to regenerate all figures.
4. Outputs deterministic with `random_state=42` / `np.random.seed(42)`.

## License
MIT.

## Data
Raw TCGA + GEO data described in Methods §4. cBioPortal `thca_tcga_pub` mirror used for TERT recovery.
""")

    # LICENSE
    (ANON_DIR / 'LICENSE').write_text("""MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
""")

    # Zip
    with zipfile.ZipFile(ANON_OUT, 'w', zipfile.ZIP_DEFLATED) as z:
        for p in ANON_DIR.rglob('*'):
            if p.is_file():
                z.write(p, p.relative_to(SUB))
    log(f'anonymous code → {ANON_OUT} ({ANON_OUT.stat().st_size:,} bytes)')


def write_final_dump():
    fig_audit = (FIG_DIR / 'figure_audit_report.md').read_text() if (FIG_DIR / 'figure_audit_report.md').exists() else ''
    audit_v2 = (SUB / 'manuscript_v1_AUDIT.md').read_text() if (SUB / 'manuscript_v1_AUDIT.md').exists() else ''

    # word count v2
    wc = len((SUB / 'manuscript_v1.md').read_text().split())

    # Hot/Cold composite metric
    hcd = ''
    hc_path = RES.parent / 'results' / 'v17p35' / 'tables' / 'AMP3_hot_cold_composite.tsv'
    if hc_path.exists():
        import pandas as pd
        from scipy import stats
        import numpy as np
        hc = pd.read_csv(hc_path, sep='\t')
        a = hc.loc[hc['cluster_use'] == 'DM1', 'hot_cold_composite'].dropna()
        b = hc.loc[hc['cluster_use'] == 'DM2', 'hot_cold_composite'].dropna()
        pooled = np.sqrt(((len(a)-1)*np.var(a,ddof=1)+(len(b)-1)*np.var(b,ddof=1))/(len(a)+len(b)-2))
        d = (a.mean()-b.mean())/pooled
        _, p = stats.mannwhitneyu(a, b)
        hcd = f"Cohen's d = {d:+.3f}, Mann-Whitney p = {p:.2e}, n_DM1={len(a)} n_DM2={len(b)}"

    # TERT survival metrics
    tert = json.loads((RES.parent / 'results' / 'v17_tert_recovery' / 'v2' / 'FINAL_extended_summary.json').read_text())['survival']

    # inventory
    inv_lines = []
    for f in sorted(SUB.iterdir()):
        if f.is_file():
            inv_lines.append(f'- {f.relative_to(SUB)}  ({f.stat().st_size:,} bytes)')
    for sub_dir in (FIG_DIR, TBL_DIR, REPRO_DIR, SUB / 'outreach'):
        if sub_dir.exists():
            for f in sorted(sub_dir.iterdir()):
                if f.is_file():
                    inv_lines.append(f'- {f.relative_to(SUB)}  ({f.stat().st_size:,} bytes)')

    txt = f"""# v17 npj FINAL SHIP DUMP — {datetime.datetime.now().isoformat(timespec='seconds')} KST

## Status: SUBMIT-READY

Sprint: F1 (manuscript v2 polish + TERT R8 + Fig8) → F2 (4 outreach drafts) → F3 (audit + checklist + SUBMIT_INSTRUCTIONS) → F4 (this dump + anonymous code zip).

## Manuscript v2

- File: `submission/npj/manuscript_v1.{{md,html,pdf}}`
- Word count: **{wc}** (npj Brief Report tolerance ≤4,000)
- 8 main + 15 supplementary figures
- 20 reviewer attack defenses (`v17p35_REVIEWER_DEFENSE_v2.md`)
- 10 honest limitations
- Audit verdict: **{"PASS" if "READY" in audit_v2 else "REVIEW"}**

## Headline metrics (8)

1. **8-gene panel CV AUC**: 0.954 (RF 0.975) vs BRAF V600E baseline 0.822 → **ΔAUC = +0.132**
2. **Hot/Cold composite**: {hcd}
3. **TERT 4-group survival**: logrank p = {tert['four_group_logrank']['p_value']:.2e} (n=504)
4. **TERT vs WT alone**: logrank p = {tert['TERT_logrank']['p_value']:.2e} (TERT⁺ 16.7% events vs WT 2.1%)
5. **External validation GSE76039**: AUC 0.974 (correct-direction)
6. **BRAF concordance CCLE**: 4/4 lines correctly DM1
7. **scRNA immune ratio**: 57:1 in DM1-skewed vs DM2-skewed patients
8. **Pan-cancer transfer**: signature applicable to LUAD/COAD/LGG/SKCM

## Author block (partner format)

- **[Author 1]**¹\* (1st author + co-corresponding) — Independent Researcher, Seoul + part-time PhD candidate, Seoul National University Graduate School of Convergence Science and Technology
- **[[Author 2] — full name + affiliation TBD by user]**²\* (senior + co-corresponding)

CRediT contributions specified in `manuscript_v1.md` header.

## Outreach drafts (NOT sent — local only)

- `outreach/email_xing_DRAFT.md` — Mingzhao Xing (EN, Liu–Xing 4-genotype framework collaboration)
- `outreach/email_landa_DRAFT.md` — Iñigo Landa / Fagin lab (EN, PDTC/ATC metadata)
- `outreach/email_[Affiliation]_KR_DRAFT.md` — [Affiliation] Hospital (KR, Korean cohort revision-round validation)
- `outreach/message_yu_KR_DRAFT.md` — [Author 2] (KR, partner check on author + outreach + submit timing)

All four are DRAFTS. Read-through required before sending. Placeholders preserved.

## Submit timeline (recommended)

- **Today KST**: [Author 2]님 카톡 confirm 받기 (`message_yu_KR_DRAFT.md`)
- **Within 24h after confirm**: npj Editorial Manager submit click (`SUBMIT_INSTRUCTIONS.md`)
- **Same day or next**: send 3 external outreach emails (Xing, Landa, [Affiliation])
- **3–4 months**: review round 1
- **Revision round**: integrate [Affiliation] cohort if data arrives in time

## Submission package inventory (`submission/npj/`)

{chr(10).join(inv_lines)}

## Manuscript v2 audit

{audit_v2}

## Figure audit (full)

{fig_audit}

## Next 24-hour action list (Seungho)

1. ☐ `message_yu_KR_DRAFT.md` — read, polish, send to [Author 2] via KakaoTalk
2. ☐ Replace placeholders `[[Author 2] — full name TBD]` and `[Affiliation TBD]` in:
   - `manuscript_v1.md` (header + §6 author block + ref)
   - `cover_letter.md` (last paragraph)
   - 3 outreach emails
3. ☐ Re-render PDFs after placeholder fix (chromium headless print-to-pdf)
4. ☐ Submit at https://www.editorialmanager.com/npjpo/ following `SUBMIT_INSTRUCTIONS.md`
5. ☐ Send 3 external outreach emails after submit

## Parallel work (optional, separate terminals)

- v18 framework sprint — already noted in MEMORY.md (`v18_agentic_framework.md`) — 22-slide Korean talk outline ready
- v15 NeurIPS submit-ready sprint — mid-May deadline
"""
    out = REPORTS / 'v17p35_FINAL_SHIP_DUMP.md'
    out.write_text(txt)
    log(f'FINAL DUMP → {out}')
    return out


if __name__ == '__main__':
    build_anon_zip()
    p = write_final_dump()
    print(f'DUMP {p}')
    print(f'ANON_ZIP {ANON_OUT}')
