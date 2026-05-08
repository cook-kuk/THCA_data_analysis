"""Build predictions.tsv + auroc_summary.tsv for BigMHC IM on ITSNdb."""
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

BUNDLE = '/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/bundle.tsv'
PRED   = '/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/wave3_bigmhc/itsndb_bigmhc_im.csv'
OUTDIR = '/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09/wave3_bigmhc'

bundle = pd.read_csv(BUNDLE, sep='\t')
sub = bundle[bundle['split'].isin(['ext_itsndb_main','ext_itsndb_val'])].copy().reset_index(drop=True)

pred = pd.read_csv(PRED)
pred = pred.rename(columns={'mhc':'hla','pep':'peptide','BigMHC_IM':'score_bigmhc_im'})

# Merge — preserve order
m = sub.merge(pred[['hla','peptide','score_bigmhc_im']],
              left_on=['HLA_norm','peptide'], right_on=['hla','peptide'], how='left')
n_skipped = m['score_bigmhc_im'].isna().sum()
print(f'merged rows: {len(m)}, n_skipped (NaN scores): {n_skipped}')

out = m[['peptide','hla','label','in_master','score_bigmhc_im','source','split']].copy()
out.to_csv(f'{OUTDIR}/predictions.tsv', sep='\t', index=False)
print(f'wrote predictions.tsv ({len(out)} rows)')

# Bootstrap AUROC
def bootstrap_auc(y_true, y_score, n_boot=1000, seed=42):
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    if len(np.unique(y_true)) < 2:
        return np.nan, np.nan, np.nan
    auc = roc_auc_score(y_true, y_score)
    rng = np.random.default_rng(seed)
    n = len(y_true)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[idx], y_score[idx]))
    aucs = np.array(aucs)
    lo, hi = np.quantile(aucs, [0.025, 0.975])
    return auc, lo, hi

clean = out.dropna(subset=['score_bigmhc_im']).copy()

rows = []
def add(name, df):
    auc, lo, hi = bootstrap_auc(df['label'].values, df['score_bigmhc_im'].values)
    rows.append({
        'testset': name,
        'in_master': 'all' if name == 'ITSNdb_combined' else ('False' if 'no_overlap' in name else 'True'),
        'n': len(df),
        'n_pos': int(df['label'].sum()),
        'AUROC': round(auc, 4) if not np.isnan(auc) else np.nan,
        'AUROC_lo95': round(lo, 4) if not np.isnan(lo) else np.nan,
        'AUROC_hi95': round(hi, 4) if not np.isnan(hi) else np.nan,
    })

add('ITSNdb_combined', clean)
add('ITSNdb_no_overlap', clean[clean['in_master'] == False])
add('ITSNdb_in_master', clean[clean['in_master'] == True])

summary = pd.DataFrame(rows)
summary.to_csv(f'{OUTDIR}/auroc_summary.tsv', sep='\t', index=False)
print(summary.to_string(index=False))
