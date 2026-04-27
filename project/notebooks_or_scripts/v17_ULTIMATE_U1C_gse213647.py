"""v17 ULTIMATE U1C — apply 8-gene RAI panel to GSE213647 (n=632, largest local cohort).

Steps:
  1. Load GSE213647 v3 (ComBat-corrected) log2 expression.
  2. Subset to GENE_8; report match count.
  3. Train 8-gene LogisticRegression on TCGA-THCA using R1A cluster_orig (DM1=1, DM2=0).
  4. Predict DM probabilities on GSE213647.
  5. Parse family.soft.gz for tumor/normal + histology subtype.
  6. Tumor-only DM1/DM2 distribution by histology.
  7. AUC: differentiated (PTC) vs dedifferentiated (ATC + PD) on tumors.
  8. Outputs:
       - U1C_gse213647_predictions.tsv
       - U1C_gse213647_summary.json
       - figures/U1C_gse213647_dm_by_histology.png
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

# Local helpers
SCR = Path('/home/seungho/personal/THCA_data_analysis/project/notebooks_or_scripts')
sys.path.insert(0, str(SCR))
from v17_ULTIMATE_common import GENE_8, load_tcga_expr, log, jdump, RES, DATA_PROC, DATA_RAW, ROOT  # noqa: E402

OUT_DIR = RES
FIG_DIR = RES / 'figures'
FIG_DIR.mkdir(parents=True, exist_ok=True)


def _build_ensg_to_symbol() -> dict:
    """Build Ensembl-stable-ID -> HGNC symbol mapping from a GDC gene-counts file.

    GDC GENCODE v36 covers ~60k genes; sufficient for the 8-gene RAI panel.
    """
    p = Path('/data/thca/data_raw/gdc/TCGA-THCA/counts/0e554902-e3ec-400b-82db-b2514ae50b99.tsv')
    m = pd.read_csv(p, sep='\t', skiprows=1)
    m = m[m['gene_id'].str.startswith('ENSG', na=False)].copy()
    m['ensg'] = m['gene_id'].str.split('.').str[0]
    # Drop duplicates (keep first)
    m = m.drop_duplicates(subset='ensg', keep='first')
    return dict(zip(m['ensg'], m['gene_name']))


def load_gse213647_expr() -> pd.DataFrame:
    """Prefer v3 ComBat-corrected log2; fall back to plain log2.

    The GSE213647 processed files have ENSG* IDs (despite the index name being
    'gene_symbol'). We translate to HGNC symbols using GDC GENCODE v36 mapping.
    """
    p_v3 = DATA_PROC / 'bulk_rnaseq_v3' / 'GSE213647_v3_log2.tsv'
    p_old = DATA_PROC / 'bulk_rnaseq' / 'GSE213647_rnaseq_expression_log2.tsv'
    if p_v3.exists():
        log(f'GSE213647: using v3 ComBat-corrected file -> {p_v3}')
        df = pd.read_csv(p_v3, sep='\t', index_col=0)
    else:
        log(f'GSE213647: v3 not found; falling back to {p_old}')
        df = pd.read_csv(p_old, sep='\t', index_col=0)
    # Translate Ensembl -> symbol if needed
    if df.index.astype(str).str.startswith('ENSG').mean() > 0.5:
        log('GSE213647: index is ENSG*; mapping to HGNC symbols via GDC v36')
        ensg2sym = _build_ensg_to_symbol()
        df.index = df.index.map(lambda e: ensg2sym.get(str(e).split('.')[0], None))
        df = df[~df.index.isna()]
        # Collapse duplicate symbols by mean (rare for protein-coding panel genes)
        if df.index.duplicated().any():
            df = df.groupby(level=0).mean()
        log(f'GSE213647: post-mapping shape {df.shape}')
    log(f'GSE213647 expr: {df.shape} (genes x samples)')
    return df


def parse_soft_metadata(soft_path: Path) -> pd.DataFrame:
    """Parse SOFT family file; return per-sample tumor/normal + histology subtype."""
    rows = []
    cur = None
    with gzip.open(soft_path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('^SAMPLE'):
                if cur is not None:
                    rows.append(cur)
                gsm = line.split('=', 1)[1].strip()
                cur = {'sample_id': gsm, 'title': '', 'source': '', 'cell_type': '', 'cell_subtype': '',
                       'genotype': '', 'treatment': ''}
            elif cur is not None:
                if line.startswith('!Sample_title'):
                    cur['title'] = line.split('=', 1)[1].strip()
                elif line.startswith('!Sample_source_name_ch1'):
                    cur['source'] = line.split('=', 1)[1].strip()
                elif line.startswith('!Sample_characteristics_ch1'):
                    val = line.split('=', 1)[1].strip()
                    m = re.match(r'^([^:]+):\s*(.*)$', val)
                    if m:
                        key = m.group(1).strip().lower().replace(' ', '_')
                        v = m.group(2).strip()
                        if key in cur:
                            cur[key] = v
        if cur is not None:
            rows.append(cur)
    meta = pd.DataFrame(rows)
    log(f'SOFT parsed: {meta.shape}; cell_type counts -> {meta["cell_type"].value_counts().to_dict()}')
    log(f'  cell_subtype counts -> {meta["cell_subtype"].value_counts().to_dict()}')
    return meta


def train_tcga_8gene_classifier():
    """Train LR on TCGA-THCA 8-gene panel using R1A cluster as label.

    Returns (scaler, model, gene_order_used).
    """
    tcga = load_tcga_expr()  # genes x samples
    labels_path = ROOT / 'project' / 'results' / 'v17_realfix' / 'R1A_cluster_labels.tsv'
    lbl = pd.read_csv(labels_path, sep='\t')
    # Map: DM1_A=1, DM2_A=0
    lbl['y'] = (lbl['cluster'] == 'DM1_A').astype(int)
    log(f'TCGA labels: {lbl.shape}; DM1={int(lbl.y.sum())}, DM2={int((1 - lbl.y).sum())}')

    common = [g for g in GENE_8 if g in tcga.index]
    log(f'TCGA: 8-gene match = {len(common)}/8 ({common})')
    X = tcga.loc[common].T  # samples x genes
    X = X.loc[X.index.intersection(lbl['sample_id'])]
    y = lbl.set_index('sample_id').loc[X.index, 'y'].values
    log(f'TCGA training matrix: {X.shape}, y_mean={y.mean():.3f}')

    scaler = StandardScaler().fit(X.values)
    Xs = scaler.transform(X.values)
    model = LogisticRegression(max_iter=5000, C=1.0, class_weight='balanced')
    model.fit(Xs, y)
    train_auc = roc_auc_score(y, model.predict_proba(Xs)[:, 1])
    log(f'TCGA train AUC (in-sample): {train_auc:.3f}')
    return scaler, model, common, float(train_auc)


def main():
    log('=== U1C: GSE213647 8-gene RAI panel ===')
    expr = load_gse213647_expr()  # genes x samples

    # Train on TCGA
    scaler, model, tcga_genes_used, train_auc = train_tcga_8gene_classifier()

    # Subset to GENE_8 in GSE213647 — only use genes present in BOTH sets
    common = [g for g in tcga_genes_used if g in expr.index]
    n_match = len(common)
    log(f'GSE213647: 8-gene match = {n_match}/8 ({common})')
    if n_match < 6:
        log('WARNING: <6 genes match in GSE213647; predictions may be unreliable.')

    # Re-train using only commonly available genes (so feature order is consistent)
    if set(common) != set(tcga_genes_used):
        log('Re-training on common-gene subset for cross-cohort consistency.')
        tcga = load_tcga_expr()
        labels_path = ROOT / 'project' / 'results' / 'v17_realfix' / 'R1A_cluster_labels.tsv'
        lbl = pd.read_csv(labels_path, sep='\t')
        lbl['y'] = (lbl['cluster'] == 'DM1_A').astype(int)
        X = tcga.loc[common].T
        X = X.loc[X.index.intersection(lbl['sample_id'])]
        y = lbl.set_index('sample_id').loc[X.index, 'y'].values
        scaler = StandardScaler().fit(X.values)
        model = LogisticRegression(max_iter=5000, C=1.0, class_weight='balanced')
        model.fit(scaler.transform(X.values), y)
        train_auc = roc_auc_score(y, model.predict_proba(scaler.transform(X.values))[:, 1])
        log(f'TCGA retrained AUC (n_genes={len(common)}): {train_auc:.3f}')

    # GSE213647 prediction
    Xg = expr.loc[common].T  # samples x genes
    Xg = Xg.dropna(axis=0, how='any')
    log(f'GSE213647 prediction matrix: {Xg.shape}')
    Xs = scaler.transform(Xg.values)
    proba = model.predict_proba(Xs)[:, 1]  # P(DM1)
    pred = (proba >= 0.5).astype(int)

    pred_df = pd.DataFrame({
        'sample_id': Xg.index,
        'n_genes_matched': n_match,
        'dm_score': proba,
        'dm_pred': pred,  # 1 = DM1, 0 = DM2
    })

    # Metadata parse
    meta = parse_soft_metadata(DATA_RAW / 'geo' / 'GSE213647_family.soft.gz')
    pred_df = pred_df.merge(meta[['sample_id', 'cell_type', 'cell_subtype', 'title']],
                            on='sample_id', how='left')
    pred_df['tumor_normal'] = pred_df['cell_type'].str.lower().map(
        {'tumor': 'tumor', 'normal': 'normal'}).fillna('unknown')
    pred_df = pred_df.rename(columns={'cell_subtype': 'histology'})

    # Tumor-only analysis
    tumors = pred_df[pred_df['tumor_normal'] == 'tumor'].copy()
    n_tumor = len(tumors)
    log(f'Tumor samples: {n_tumor}')
    hist_counts_all = pred_df['histology'].value_counts().to_dict()
    hist_counts_tumor = tumors['histology'].value_counts().to_dict()
    log(f'Tumor histology counts: {hist_counts_tumor}')

    # AUC: differentiated (PTC) vs dedifferentiated (ATC, PD) on TUMORS only
    diff_set = {'PTC'}
    dediff_set = {'ATC', 'PD'}
    auc_block = {}
    auc_value = None
    sub = tumors[tumors['histology'].isin(diff_set | dediff_set)].copy()
    if sub['histology'].nunique() >= 2 and sub['histology'].isin(dediff_set).sum() >= 3 \
            and sub['histology'].isin(diff_set).sum() >= 3:
        # DM2 (0) score = 1 - dm_score; here higher dm_score = DM1.
        # For RAI panel: DM2 = differentiated/RAI-avid, DM1 = dedifferentiated/RAI-refractory.
        # So expected: dediff (ATC/PD) -> high dm_score (DM1).
        y_true = sub['histology'].isin(dediff_set).astype(int).values  # 1 = dedifferentiated
        y_score = sub['dm_score'].values
        auc_value = float(roc_auc_score(y_true, y_score))
        # Symmetric (in case orientation is flipped)
        auc_flipped = float(roc_auc_score(y_true, -y_score))
        auc_block = {
            'comparison': 'PTC (differentiated) vs ATC+PD (dedifferentiated)',
            'n_diff': int((sub['histology'].isin(diff_set)).sum()),
            'n_dediff': int((sub['histology'].isin(dediff_set)).sum()),
            'auc_dm1_high_eq_dediff': auc_value,
            'auc_flipped': auc_flipped,
        }
        log(f'AUC PTC vs ATC+PD (DM1=dediff orientation): {auc_value:.3f}; flipped={auc_flipped:.3f}')
    else:
        log('Insufficient histology classes for AUC; descriptive only.')

    # Mean DM score per histology (tumors only)
    dm_by_hist = tumors.groupby('histology')['dm_score'].agg(['count', 'mean', 'std']).round(4)
    log(f'Tumor DM score by histology:\n{dm_by_hist}')

    # Save predictions tsv
    pred_path = OUT_DIR / 'U1C_gse213647_predictions.tsv'
    pred_df[['sample_id', 'n_genes_matched', 'dm_score', 'dm_pred',
             'histology', 'tumor_normal']].to_csv(pred_path, sep='\t', index=False)
    log(f'wrote {pred_path}')

    # JSON summary
    summary = {
        'cohort': 'GSE213647',
        'n_samples_total': int(len(pred_df)),
        'n_tumor': int(n_tumor),
        'n_normal': int((pred_df['tumor_normal'] == 'normal').sum()),
        'n_genes_matched': int(n_match),
        'genes_used': common,
        'genes_missing': [g for g in GENE_8 if g not in common],
        'tcga_train_auc': train_auc,
        'histology_counts_all': hist_counts_all,
        'histology_counts_tumor': hist_counts_tumor,
        'tumor_dm1_pred_count': int((tumors['dm_pred'] == 1).sum()),
        'tumor_dm2_pred_count': int((tumors['dm_pred'] == 0).sum()),
        'tumor_dm_score_mean_by_histology': {
            k: {'count': int(v['count']), 'mean': float(v['mean']),
                'std': float(v['std']) if not np.isnan(v['std']) else None}
            for k, v in dm_by_hist.to_dict('index').items()
        },
        'auc_block': auc_block,
        'note_braf': 'BRAF mutation status not provided in GSE213647 metadata; ΔAUC vs BRAF not computable.',
        'data_source_expr': 'GSE213647_v3_log2.tsv (ComBat v3)',
        'tcga_label_source': 'project/results/v17_realfix/R1A_cluster_labels.tsv (cluster column; DM1_A=1, DM2_A=0)',
    }
    if auc_value is None:
        summary['note_auc'] = 'no ground-truth label, descriptive only'
    jdump(summary, OUT_DIR / 'U1C_gse213647_summary.json')

    # Figure
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    # Left: violin/strip of dm_score by histology (tumors only)
    order = [h for h in ['PTC', 'PD', 'ATC', 'Follicular neoplasm', 'Nodular hyperplasia']
             if h in tumors['histology'].unique()]
    data = [tumors.loc[tumors['histology'] == h, 'dm_score'].values for h in order]
    parts = axes[0].violinplot(data, showmeans=True, showmedians=False)
    for pc in parts['bodies']:
        pc.set_facecolor('#4C72B0')
        pc.set_edgecolor('black')
        pc.set_alpha(0.7)
    axes[0].set_xticks(range(1, len(order) + 1))
    axes[0].set_xticklabels([f'{h}\nn={len(d)}' for h, d in zip(order, data)], fontsize=9)
    axes[0].set_ylabel('Predicted DM1 probability\n(higher = more dedifferentiated)')
    axes[0].set_title(f'GSE213647 tumors (n={n_tumor}) — 8-gene RAI panel DM score by histology')
    axes[0].axhline(0.5, ls='--', color='red', alpha=0.5)
    axes[0].grid(axis='y', alpha=0.3)

    # Right: stacked bar of DM1 vs DM2 calls per histology
    counts = tumors.groupby(['histology', 'dm_pred']).size().unstack(fill_value=0)
    counts = counts.reindex(order)
    if 0 not in counts.columns:
        counts[0] = 0
    if 1 not in counts.columns:
        counts[1] = 0
    counts = counts[[0, 1]].rename(columns={0: 'DM2 (RAI-avid)', 1: 'DM1 (RAI-refractory)'})
    counts.plot(kind='bar', stacked=True, ax=axes[1],
                color=['#55A868', '#C44E52'], edgecolor='black')
    axes[1].set_ylabel('Tumor samples')
    axes[1].set_xlabel('Histology subtype')
    axes[1].set_title('DM1/DM2 call counts by histology (tumors)')
    axes[1].set_xticklabels(counts.index, rotation=20, ha='right')
    axes[1].grid(axis='y', alpha=0.3)

    auc_txt = (f'\nAUC PTC vs ATC+PD: {auc_value:.3f}'
               if auc_value is not None else '\n(no AUC: insufficient histology contrast)')
    fig.suptitle(f'U1C: GSE213647 (n={len(pred_df)}, {n_match}/8 genes){auc_txt}',
                 fontsize=12, fontweight='bold')
    fig.tight_layout()
    fig_path = FIG_DIR / 'U1C_gse213647_dm_by_histology.png'
    fig.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    log(f'wrote {fig_path}')

    log('=== U1C done ===')
    return summary


if __name__ == '__main__':
    main()
