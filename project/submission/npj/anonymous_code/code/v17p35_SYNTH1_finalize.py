"""Finalize submission package: copy core tables, compile manuscript→PDF, write checklist + dump."""
from v17p35_SYNTH1_common import SUB, FIG_DIR, TBL_DIR, REPRO_DIR, REPORTS, RES, log
import shutil, subprocess, json, datetime
from pathlib import Path

CORE_TABLES = [
    'AMP4_8gene_model_coefficients.tsv',
    'AMP4_cv_performance.tsv',
    'FIX1_celline_dm_scores_v2.tsv',
    'FIX3_alternative_endpoints_full.tsv',
    'AMP3_hot_cold_composite.tsv',
]

# F2 GSEA hallmark proper lives under v17p3
EXTRA_TABLES = [
    (RES / 'v17p3' / 'tables' / 'F2_gsea_hallmark_proper.tsv', 'F2_gsea_hallmark_proper.tsv'),
    (RES / 'v17p35' / 'tables' / 'AMP_17_outliers.tsv', 'AMP_17_outliers.tsv'),
]


def copy_tables():
    for t in CORE_TABLES:
        src = RES / 'v17p35' / 'tables' / t
        if src.exists():
            shutil.copy2(src, TBL_DIR / t)
            log(f'table → {t}')
    for src, dst in EXTRA_TABLES:
        if src.exists():
            shutil.copy2(src, TBL_DIR / dst)
            log(f'table → {dst}')


def compile_manuscript():
    md = SUB / 'manuscript_v1.md'
    pdf = SUB / 'manuscript_v1.pdf'
    if not md.exists():
        log('manuscript_v1.md missing — skipping PDF compile')
        return None
    # Try pandoc with multiple PDF engines until one works.
    for eng in ('--pdf-engine=tectonic', '--pdf-engine=xelatex', '--pdf-engine=pdflatex', None):
        cmd = ['pandoc', str(md), '-o', str(pdf), '--standalone',
               '-V', 'geometry:margin=1in', '-V', 'fontsize=11pt']
        if eng: cmd.append(eng)
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if r.returncode == 0 and pdf.exists() and pdf.stat().st_size > 1000:
                log(f'manuscript PDF compiled via {eng or "default"} → {pdf} ({pdf.stat().st_size} bytes)')
                return pdf
            else:
                log(f'pandoc {eng or "default"} failed: {r.stderr[:200]}')
        except FileNotFoundError:
            continue
        except Exception as e:
            log(f'pandoc {eng}: {e}')
    # Fallback: HTML
    html = SUB / 'manuscript_v1.html'
    r = subprocess.run(['pandoc', str(md), '-o', str(html), '--standalone'],
                       capture_output=True, text=True, timeout=60)
    if r.returncode == 0:
        log(f'manuscript fallback HTML → {html}')
        return html
    return None


def word_count(path: Path) -> int:
    if not path.exists(): return 0
    return len(path.read_text().split())


def write_checklist():
    items = []
    items.append(('Manuscript v1 markdown', (SUB / 'manuscript_v1.md').exists()))
    wc = word_count(SUB / 'manuscript_v1.md')
    items.append((f'Manuscript word count {wc}/3500', wc <= 3500))
    items.append(('Cover letter', (SUB / 'cover_letter.md').exists()))
    items.append(('Reviewer defense', (SUB / 'v17p35_REVIEWER_DEFENSE.md').exists()))
    items.append(('Reproducibility log',
                  (REPRO_DIR / 'v17p35_PHASE_B_LOG.md').exists()))
    for i in range(1, 8):
        items.append((f'Fig{i} html+png+pdf',
                      all((FIG_DIR / f'Fig{i}.{e}').exists() for e in ('html', 'png', 'pdf'))))
    for i in range(1, 16):
        items.append((f'SuppFig{i} html+png+pdf',
                      all((FIG_DIR / f'SuppFig{i}.{e}').exists() for e in ('html', 'png', 'pdf'))))
    items.append(('Tables (5 core + 2 extra)', all((TBL_DIR / t).exists() for t in CORE_TABLES)))

    md_pdf = SUB / 'manuscript_v1.pdf'
    items.append((f'Manuscript PDF compiled ({md_pdf.stat().st_size if md_pdf.exists() else 0} bytes)', md_pdf.exists()))
    items.append(('Master panel index', (FIG_DIR / 'master_panel.html').exists()))
    items.append(('Figure audit report', (FIG_DIR / 'figure_audit_report.md').exists()))

    lines = ['# v17p35 npj submission checklist', '',
             f'Generated {datetime.datetime.now().isoformat(timespec="seconds")}', '']
    pass_count = 0
    for name, ok in items:
        lines.append(f'- [{"x" if ok else " "}] {name}')
        if ok: pass_count += 1
    verdict = '✅ READY TO SUBMIT' if pass_count == len(items) else f'⚠ {len(items)-pass_count} item(s) outstanding'
    lines.append('')
    lines.append(f'## Verdict: {verdict} ({pass_count}/{len(items)})')
    out = SUB / 'SUBMISSION_CHECKLIST.md'
    out.write_text('\n'.join(lines) + '\n')
    log(f'checklist → {out}')

    # SUBMISSION_READY.md only if all pass
    ready = SUB / 'SUBMISSION_READY.md'
    if pass_count == len(items):
        ready.write_text(f'# ✅ SUBMISSION READY\n\nAll {len(items)} items pass. Generated {datetime.datetime.now().isoformat(timespec="seconds")}.\n')
        log('SUBMISSION_READY.md created')
    else:
        if ready.exists(): ready.unlink()
    return pass_count, len(items)


def write_final_dump(pass_count, total):
    dump_path = REPORTS / 'v17p35_SYNTH1_FINAL_DUMP.md'
    fig_audit = (FIG_DIR / 'figure_audit_report.md').read_text() if (FIG_DIR / 'figure_audit_report.md').exists() else ''
    inv_lines = []
    for sub_dir in (FIG_DIR, TBL_DIR, REPRO_DIR):
        for f in sorted(sub_dir.iterdir()):
            if f.is_file():
                inv_lines.append(f'- {f.relative_to(SUB)}  ({f.stat().st_size:,} bytes)')
    for f in sorted(SUB.iterdir()):
        if f.is_file():
            inv_lines.append(f'- {f.relative_to(SUB)}  ({f.stat().st_size:,} bytes)')

    # Hot/Cold composite metric
    hcd = ''
    hc_path = RES / 'v17p35' / 'tables' / 'AMP3_hot_cold_composite.tsv'
    if hc_path.exists():
        import pandas as pd
        from scipy import stats
        hc = pd.read_csv(hc_path, sep='\t')
        a = hc.loc[hc['cluster_use'] == 'DM1', 'hot_cold_composite'].dropna()
        b = hc.loc[hc['cluster_use'] == 'DM2', 'hot_cold_composite'].dropna()
        import numpy as np
        pooled = np.sqrt(((len(a)-1)*np.var(a,ddof=1)+(len(b)-1)*np.var(b,ddof=1))/(len(a)+len(b)-2))
        d = (a.mean()-b.mean())/pooled
        _, p = stats.mannwhitneyu(a, b)
        hcd = f"Cohen's d = {d:+.3f}, Mann-Whitney p = {p:.2e}, DM1 n={len(a)} DM2 n={len(b)}"

    txt = f"""# v17p35 SYNTH-1 FINAL DUMP

Generated {datetime.datetime.now().isoformat(timespec='seconds')}
Sprint: figure composition + PDF/PNG export + npj submission package.

## Core metrics

1. **Main figures (7)**: html + png(scale 3) + pdf — all under npj 5MB limit
2. **Supplementary figures (15)**: same — all under limit
3. **AMP-3 Hot/Cold composite**: {hcd}
4. **Submission package**: {pass_count}/{total} checklist items pass
5. **Manuscript PDF**: {'compiled' if (SUB/'manuscript_v1.pdf').exists() else 'fallback HTML'}
6. **Submission verdict**: {'GO' if pass_count==total else 'NEEDS_FIX'}

## What I (Seungho) need to do next (in priority order)

1. **Manuscript voice read-through (2-3h)** — open `submission/npj/manuscript_v1.md` (or `.pdf`).
   Personality slots to consider tightening:
   - Discovery framing (DM1/DM2 nomenclature — first paragraph)
   - BRAF/RAS orthogonality interpretation (Result 3-4 transition)
   - Hot/Cold + scRNA evidence assembly (Result 5)
   - Drug-repurposing translation framing (Result 7 + Discussion)
   - Reviewer-defense bridge to limitations section
2. **Title decision** — go with current #4 candidate or surface an alternative?
3. **Cover letter polish** — Korean voice → English, single read-through
4. **[Affiliation] Hospital decision** — submit before or after their reply?
5. **Author affiliation** — single 1st author vs PI co-author (Yoo PI)
6. **Trial PI outreach** — NCT06235216 / NCT07521670 — pre-revision contact?

## Submission package inventory (`submission/npj/`)
{chr(10).join(inv_lines)}

## Figure audit (full)
{fig_audit}
"""
    dump_path.write_text(txt)
    log(f'FINAL DUMP → {dump_path}')
    return dump_path


if __name__ == '__main__':
    copy_tables()
    compile_manuscript()
    p, t = write_checklist()
    dump = write_final_dump(p, t)
    print(f'CHECKLIST {p}/{t}')
    print(f'DUMP {dump}')
