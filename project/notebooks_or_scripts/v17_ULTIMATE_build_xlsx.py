"""Consolidate all U1-U3 result TSVs + key JSONs into a single XLSX with pretty formatting."""
import json, os
from pathlib import Path
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

ROOT = Path('/home/seungho/personal/THCA_data_analysis')
RES = ROOT / 'project' / 'results' / 'v17_ultimate'
OUT = ROOT / 'project' / 'submission' / 'npj' / 'Supplementary_Tables_v6_ULTIMATE.xlsx'

# style helpers
HDR = Font(bold=True, color='FFFFFF', size=11)
HDR_FILL = PatternFill('solid', fgColor='1E3A8A')
ALT_FILL = PatternFill('solid', fgColor='F3F4F6')
BORDER = Border(left=Side(style='thin', color='D1D5DB'), right=Side(style='thin', color='D1D5DB'),
                top=Side(style='thin', color='D1D5DB'), bottom=Side(style='thin', color='D1D5DB'))


def style_sheet(ws, df, title=None, note=None):
    if title:
        ws['A1'] = title
        ws['A1'].font = Font(bold=True, size=14, color='1E3A8A')
        ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=max(2, len(df.columns)))
        start_row = 3
    else:
        start_row = 1
    if note:
        ws.cell(row=start_row - 1, column=1).value = note
        ws.cell(row=start_row - 1, column=1).font = Font(italic=True, size=10, color='6B7280')
        start_row += 1

    # write columns
    for j, col in enumerate(df.columns, 1):
        c = ws.cell(row=start_row, column=j, value=str(col))
        c.font = HDR
        c.fill = HDR_FILL
        c.alignment = Alignment(horizontal='left', vertical='center')
        c.border = BORDER

    # write data
    for i, row in enumerate(df.itertuples(index=False), 1):
        for j, val in enumerate(row, 1):
            c = ws.cell(row=start_row + i, column=j, value=val if not pd.isna(val) else '')
            c.alignment = Alignment(horizontal='left', vertical='top', wrap_text=False)
            c.border = BORDER
            if i % 2 == 0:
                c.fill = ALT_FILL
            # numeric formatting
            if isinstance(val, float):
                if abs(val) < 0.01 and val != 0:
                    c.number_format = '0.0000'
                else:
                    c.number_format = '0.000'

    # auto width (capped)
    for j, col in enumerate(df.columns, 1):
        max_len = max([len(str(col))] + [len(str(v)) for v in df[col].astype(str).head(50)])
        ws.column_dimensions[get_column_letter(j)].width = min(36, max(10, max_len + 2))

    # freeze header
    ws.freeze_panes = ws.cell(row=start_row + 1, column=1)


# read inputs
sheets = []  # list of (name, df, title, note)

# S1 — U1A 4-variant summary
df = pd.read_csv(RES / 'U1A_4variants_summary.tsv', sep='\t')
sheets.append(('S1_U1A_4variants', df,
    'Supplementary Table S1 · U1A 4 leak-free variants (de-circularization)',
    'Source: REALFIX R1-R2 / U1A consolidation; n=500 TCGA-THCA primary tumours; 5-fold CV; bootstrap 1000 95% CI.'))

# S2 — U1B per-cohort AUC (7-cohort meta)
df = pd.read_csv(RES / 'U1B_per_cohort_real_aucs.tsv', sep='\t')
sheets.append(('S2_U1B_7cohort_AUC', df,
    'Supplementary Table S2 · U1B 7-cohort external AUC (★ including RAI clinical)',
    '8 cohorts attempted (TCGA in-sample + 7 external); GSE192683 was non-thyroid (excluded); GSE151179 = ★ true RAI clinical ground truth.'))

# S3 — U1B RAI clinical detail
df = pd.read_csv(RES / 'U1B_RAI_clinical_validation.tsv', sep='\t')
sheets.append(('S3_RAI_clinical_GSE151179', df,
    'Supplementary Table S3 · ★ RAI clinical validation (GSE151179)',
    'Pu RAI series; refractory vs avid tumour annotation; n=39 specimens after primary/LNM split.'))

# S4 — U2A WHO 2022 mapping
df = pd.read_csv(RES / 'U2A_who2022_mapping.tsv', sep='\t')
sheets.append(('S4_U2A_WHO2022', df,
    'Supplementary Table S4 · U2A WHO 2022 morphologic mapping (n=476)',
    'TCGA primary_diagnosis → WHO 2022 binary axis; FVPTC mixed = RAS-like proxy at medium confidence (TCGA limitation).'))

# S5 — U2D K-metrics
df = pd.read_csv(RES / 'U2D_K_metrics.tsv', sep='\t')
sheets.append(('S5_U2D_K_selection', df,
    'Supplementary Table S5 · U2D K=2 vs K=3,4,5 internal-validity metrics',
    'Variant A leak-free (47 genes), KMeans random_state=42, n_init=20.'))

# S6 — U2D K=2 vs K=4 confusion (multi-block file → parse first block)
def _read_first_block(path, sep='\t'):
    rows = []
    in_header = False
    cols = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('#') or not line.strip():
                if cols and rows:
                    break
                continue
            parts = line.split(sep)
            if cols is None:
                cols = parts
            else:
                rows.append(parts)
    return pd.DataFrame(rows, columns=cols)

df = _read_first_block(RES / 'U2D_K2_vs_K4_confusion.tsv')
sheets.append(('S6_U2D_K2_K4_nest', df,
    'Supplementary Table S6 · U2D K=4 nests inside K=2',
    'Each K=4 partition → Pu 2022 canonical subtype name (Stromal / CNV-enriched / Immune-enriched / BRAF-enriched).'))

# S7 — U3B Korean TERT comparison
df = pd.read_csv(RES / 'U3B_korean_tert_comparison.tsv', sep='\t')
sheets.append(('S7_U3B_Korean_TERT', df,
    'Supplementary Table S7 · U3B 3-cohort TERT promoter prevalence',
    'Yang H et al. (ENM 2022, n=2,092 Korean) vs TCGA-THCA (n=513) vs MSK-IMPACT 2017 (advanced/refractory).'))

# S8 — U3C cBioPortal mutation counts
df = pd.read_csv(RES / 'U3C_per_study_mutation_counts.tsv', sep='\t')
sheets.append(('S8_U3C_cBioPortal_muts', df,
    'Supplementary Table S8 · U3C cBioPortal cross-cohort mutation prevalence',
    '6 thyroid studies, 2,194 total samples; TERT C228T/C250T separately; BRAF V600E; HRAS/KRAS/NRAS hotspots; TP53.'))

# S9 — U3C consolidated supplementary table
df = pd.read_csv(RES / 'U3C_supplementary_table_S20.tsv', sep='\t')
sheets.append(('S9_U3C_per_study_consolidated', df,
    'Supplementary Table S9 · U3C consolidated paper-ready cross-cohort table',
    'Rows = histology × gene × hotspot; columns = per-study prevalence + meta-pooled.'))

# S10 — U3A RAI landscape (cohort options)
df = pd.read_csv(RES / 'U3A_rai_landscape.tsv', sep='\t')
sheets.append(('S10_U3A_RAI_cohort_landscape', df,
    'Supplementary Table S10 · U3A RAI clinical ground-truth cohort landscape',
    'Mu Z 2024 controlled-access; GSE151179 = best open alternative (used in S2/S3).'))

# S11 — U2C Pu 2021 BRAF-like-B overlap
df = pd.read_csv(RES / 'U2C_braf_like_B_overlap.tsv', sep='\t')
sheets.append(('S11_U2C_Pu2021_BRAF_like_B', df,
    'Supplementary Table S11 · U2C Pu 2021 BRAF-like-B subtype overlap',
    'DM2 vs Pu 2021 dedifferentiation signature direction concordance; 88.9% (8/9 marker genes).'))

# S12 — U1D methylation cluster
df = pd.read_csv(RES / 'U1D_methylation_cluster.tsv', sep='\t')
sheets.append(('S12_U1D_methylation_cluster', df,
    'Supplementary Table S12 · U1D GSE97466 methylation cluster (n=141)',
    'Top-5000 β-values, K=2 KMeans random_state=42; 100% of aggressive histologies → met_DM2.'))

# === build workbook ===
wb = Workbook()
wb.remove(wb.active)

# README sheet first
ws = wb.create_sheet('README', 0)
ws['A1'] = 'Supplementary Tables — v17 ULTIMATE (npj Precision Oncology submission)'
ws['A1'].font = Font(bold=True, size=14, color='1E3A8A')
ws.merge_cells('A1:D1')

ws['A3'] = 'Author'; ws['A3'].font = Font(bold=True)
ws['B3'] = 'Seungho Cook (Independent Researcher, Seoul) · with [Yu Kyungho]'

ws['A4'] = 'Date'; ws['A4'].font = Font(bold=True)
ws['B4'] = '2026-04-27 KST'

ws['A5'] = 'Version'; ws['A5'].font = Font(bold=True)
ws['B5'] = 'v6 ULTIMATE (Scenario A · npj Precision Oncology)'

ws['A6'] = 'Manuscript'; ws['A6'].font = Font(bold=True)
ws['B6'] = 'manuscript_v6_ULTIMATE.{md,html,pdf,docx}'

ws['A8'] = 'Sheet'; ws['A8'].font = HDR; ws['A8'].fill = HDR_FILL
ws['B8'] = 'Title'; ws['B8'].font = HDR; ws['B8'].fill = HDR_FILL
ws['C8'] = 'Source'; ws['C8'].font = HDR; ws['C8'].fill = HDR_FILL
ws['D8'] = 'n rows'; ws['D8'].font = HDR; ws['D8'].fill = HDR_FILL

for i, (name, df, title, note) in enumerate(sheets, 1):
    r = 8 + i
    ws.cell(row=r, column=1, value=name).font = Font(name='JetBrains Mono', size=10, bold=True, color='1E3A8A')
    ws.cell(row=r, column=2, value=title).alignment = Alignment(wrap_text=True)
    ws.cell(row=r, column=3, value=note).alignment = Alignment(wrap_text=True)
    ws.cell(row=r, column=4, value=len(df))

ws.column_dimensions['A'].width = 28
ws.column_dimensions['B'].width = 60
ws.column_dimensions['C'].width = 70
ws.column_dimensions['D'].width = 8

# data sheets
for name, df, title, note in sheets:
    ws = wb.create_sheet(name)
    style_sheet(ws, df, title, note)

wb.save(OUT)
print(f'wrote {OUT} with {1 + len(sheets)} sheets')
print(f'size: {OUT.stat().st_size:,} bytes')
