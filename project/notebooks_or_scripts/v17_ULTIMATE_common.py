"""Common helpers for v17 ULTIMATE FIX sprint."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path('/home/seungho/personal/THCA_data_analysis')
DATA_PROC = Path('/data/thca/data_processed')
DATA_RAW = Path('/data/thca/data_raw')
SCR = ROOT / 'project' / 'notebooks_or_scripts'
RES = ROOT / 'project' / 'results' / 'v17_ultimate'
SUB = ROOT / 'project' / 'submission' / 'npj'
META = ROOT / 'project' / 'metadata'
RES.mkdir(parents=True, exist_ok=True)
(RES / 'figures').mkdir(exist_ok=True)

GENE_8 = ['SLC5A5', 'TPO', 'TG', 'TSHR', 'PAX8', 'NKX2-1', 'FOXE1', 'DIO1']

# leak-free variant gene sets (zero overlap with GENE_8)
VARIANTS = {
    'A_TIERA67_minus8': None,  # populated below
    'B_MAPK_immune_EMT': [
        # MAPK
        'DUSP4', 'DUSP5', 'DUSP6', 'SPRY1', 'SPRY2', 'SPRY4', 'ETV4', 'ETV5', 'PHLDA1', 'FOSL1',
        # Immune
        'CD274', 'CD8A', 'FOXP3', 'IDO1', 'HLA-DRA', 'PRF1', 'GZMB', 'GZMA',
        # EMT
        'VIM', 'ZEB1', 'ZEB2', 'SNAI1', 'SNAI2', 'TWIST1', 'CDH1', 'CDH2',
        # Cell cycle
        'TP53', 'CDKN2A', 'CDKN2B', 'MKI67',
    ],
    'C_immune_5feature': ['CD8A', 'PRF1', 'GZMB', 'IFNG', 'CXCL10'],
    'D_BRS71_minus8': None,  # populated below
}


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


def load_tcga_expr() -> pd.DataFrame:
    """Return genes × samples dataframe for TCGA-THCA primary tumours."""
    p = DATA_PROC / 'bulk_rnaseq' / 'TCGA-THCA_rnaseq_expression_log2.tsv'
    df = pd.read_csv(p, sep='\t', index_col=0)
    log(f'TCGA-THCA loaded: {df.shape} (genes × samples)')
    return df


def load_dataset_master() -> pd.DataFrame:
    return pd.read_csv(META / 'dataset_master.tsv', sep='\t')


def load_brs71_genes() -> list[str]:
    p = META / 'brs71_genes.txt'
    return [l.strip() for l in p.read_text().splitlines()
            if l.strip() and not l.strip().startswith('#')]


def load_gene_panels() -> dict:
    return json.loads((META / 'gene_panels.json').read_text())


def load_tiera67_genes() -> list[str]:
    """Load TIERA67 used in original DM1/DM2 cluster definition."""
    for cand in (META / 'tierA67_genes.txt', META / 'v3_tierA67_clean_genes.txt'):
        if cand.exists():
            return [l.strip() for l in cand.read_text().splitlines()
                    if l.strip() and not l.strip().startswith('#')]
    return []


def load_orig_dm_labels() -> pd.DataFrame:
    """Load original v17 DM1/DM2 cluster labels for the 500-sample subset.

    Falls back to R1A which was generated from TIERA67-8gene (close proxy).
    Returns a dataframe with columns: sample_id, dm_orig (0=DM2, 1=DM1).
    """
    # The REALFIX R1 was the de-circularized version. The "original" label used in
    # R1*_summary.json is the v17 labels. We pull that from the consensus output.
    p = ROOT / 'project' / 'results' / 'v17_realfix' / 'R1A_cluster_labels.tsv'
    df = pd.read_csv(p, sep='\t')
    # cluster_orig is the v17 original; cluster_new is R1A's leak-free cluster.
    return df[['sample_id', 'cluster_orig']].rename(columns={'cluster_orig': 'dm_orig'})


def populate_variants():
    """Fill VARIANTS A and D with concrete gene lists."""
    if VARIANTS['A_TIERA67_minus8'] is None:
        tiera = load_tiera67_genes()
        if tiera:
            VARIANTS['A_TIERA67_minus8'] = [g for g in tiera if g not in GENE_8]
    if VARIANTS['D_BRS71_minus8'] is None:
        brs = load_brs71_genes()
        if brs:
            VARIANTS['D_BRS71_minus8'] = [g for g in brs if g not in GENE_8]


def jdump(obj, path):
    Path(path).write_text(json.dumps(obj, indent=2, default=str))
    log(f'  wrote {path}')


if __name__ == '__main__':
    populate_variants()
    log('=== variant gene counts ===')
    for k, v in VARIANTS.items():
        log(f'  {k}: {len(v) if v else "(unfilled)"} genes')
    log('=== dataset_master ===')
    dm = load_dataset_master()
    log(f'  {dm.shape}; cohorts: {dm.columns.tolist()[:6]}')
