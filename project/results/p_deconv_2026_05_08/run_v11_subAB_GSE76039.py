"""v11: sub-A vs sub-B MAPK + GSE76039 advanced disease MAPK gradient."""
from pathlib import Path
import warnings, gzip, re
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats

OUT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_deconv_2026_05_08")

MAPK = ['DUSP4','DUSP5','DUSP6','SPRY2','SPRY4','ETV4','ETV5','PHLDA1','CCND1']
PANEL = ['DIO1','FOXE1','NKX2-1','PAX8','SLC5A5','TG','TPO','TSHR']

def patient_id(s): return "-".join(s.split("-")[:3])

# === Part 1: sub-A vs sub-B MAPK score in TCGA ===
print("=== Part 1: sub-A vs sub-B MAPK ===")
tcga = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_zscore.tsv", sep="\t", index_col=0)
tcga_mapk = tcga.loc[MAPK].mean(axis=0)
tcga_panel = tcga.loc[PANEL].mean(axis=0)

sub = pd.read_csv("/data/thca/repo_results/d6p7_dm1_subcluster/dm1_subcluster_labels.tsv", sep="\t")
sub.columns = ['sample_id','sub_cluster']
sub['patient'] = sub['sample_id'].apply(patient_id)

sample_df = pd.DataFrame({'sample':tcga_mapk.index,'mapk':tcga_mapk.values,'panel':tcga_panel.values})
sample_df['patient'] = sample_df['sample'].apply(patient_id)
m = sample_df.merge(sub, on='patient', how='inner')
print(f"sub × TCGA RNA: {len(m)} rows, {m['sub_cluster'].value_counts().to_dict()}")

a = m.loc[m['sub_cluster']=='sub_A']
b = m.loc[m['sub_cluster']=='sub_B']
def cd(x,y):
    x,y = np.asarray(x), np.asarray(y)
    if len(x)<2 or len(y)<2: return 0
    p = np.sqrt(((len(x)-1)*x.var(ddof=1)+(len(y)-1)*y.var(ddof=1))/(len(x)+len(y)-2))
    return (x.mean()-y.mean())/p if p>0 else 0

print(f"\nsub-A (n={len(a)}) vs sub-B (n={len(b)}):")
print(f"  MAPK: A mean={a['mapk'].mean():.3f}, B mean={b['mapk'].mean():.3f}")
print(f"  d(A−B) MAPK = {cd(a['mapk'].values, b['mapk'].values):+.3f}")
print(f"  d(A−B) Panel = {cd(a['panel'].values, b['panel'].values):+.3f}")
_, p_mapk = stats.mannwhitneyu(a['mapk'], b['mapk'])
_, p_panel = stats.mannwhitneyu(a['panel'], b['panel'])
print(f"  p (MAPK) = {p_mapk:.3g}, p (Panel) = {p_panel:.3g}")

pd.DataFrame([
    {'metric':'MAPK','subA_mean':a['mapk'].mean(),'subB_mean':b['mapk'].mean(),'cohen_d':cd(a['mapk'].values,b['mapk'].values),'p':p_mapk},
    {'metric':'Panel','subA_mean':a['panel'].mean(),'subB_mean':b['panel'].mean(),'cohen_d':cd(a['panel'].values,b['panel'].values),'p':p_panel}
]).to_csv(OUT/"v11_subAB_mapk_panel.tsv", sep="\t", index=False)

# === Part 2: GSE76039 advanced disease (PDTC vs ATC) MAPK ===
print("\n=== Part 2: GSE76039 advanced disease ===")
import io
gz_path = "/data/thca/repo_results/paper3_ici_public_data/GSE76039/GSE76039_series_matrix.txt.gz"
with gzip.open(gz_path, 'rt') as f:
    lines = f.readlines()

# Find table block
start = None; end = None; source_line = None
for i, line in enumerate(lines):
    if line.startswith('!series_matrix_table_begin'): start = i+1
    elif line.startswith('!series_matrix_table_end'): end = i
    elif line.startswith('!Sample_source_name_ch1'): source_line = line
print(f"  Table: lines {start}–{end}")

table = pd.read_csv(io.StringIO(''.join(lines[start:end])), sep="\t", index_col=0)
table.index.name = 'probe'
print(f"  Expression: {table.shape}")

# Sample-level histology from source line
src_parts = source_line.strip().split('\t')[1:]
src_cleaned = [p.strip('"') for p in src_parts]
gsm_line = [l for l in lines if l.startswith('!Sample_geo_accession')][0]
gsms = [g.strip('"') for g in gsm_line.strip().split('\t')[1:]]
hist_map = dict(zip(gsms, src_cleaned))
print(f"  Histology: {pd.Series(list(hist_map.values())).value_counts().to_dict()}")

# Map probes (GPL570 affy hgu133plus2) → gene symbols
# Affy IDs need annotation — use conservative match: search NetAffx-style or use a probe→symbol file if present
probe_anno_path = Path("/data/thca/_probe_anno_GPL570.tsv")
if not probe_anno_path.exists():
    # Fallback: use TPM index gene names if any TCGA probes match (won't); skip and use direct text mining via annotation
    # Try Bioconductor-equivalent fallback: use hardcoded known probes
    PROBE_TO_GENE = {
        'DUSP4': ['204014_at','204015_s_at'],
        'DUSP5': ['209457_at','209458_x_at'],
        'DUSP6': ['208891_at','208892_s_at','208893_s_at'],
        'SPRY2': ['204011_at','204012_s_at'],
        'SPRY4': ['221489_s_at','221490_at'],
        'ETV4':  ['203829_s_at','211603_s_at'],
        'ETV5':  ['203349_s_at','203350_at'],
        'PHLDA1':['204602_at','217997_at','225842_at'],
        'CCND1': ['208711_s_at','208712_at','214019_at'],
        'DIO1':  ['206457_s_at','206458_x_at'],
        'FOXE1': ['206912_at'],
        'NKX2-1':['203929_s_at','205317_s_at'],
        'PAX8':  ['208795_s_at','208796_s_at'],
        'SLC5A5':['208486_at'],
        'TG':    ['205969_at','206924_s_at'],
        'TPO':   ['203689_s_at','205255_x_at'],
        'TSHR':  ['205455_at'],
    }
else:
    pa = pd.read_csv(probe_anno_path, sep="\t")
    PROBE_TO_GENE = pa.groupby('symbol')['probe'].apply(list).to_dict()

def gene_score(genes, table, probe_map):
    scores = []
    for g in genes:
        probes = probe_map.get(g, [])
        present = [p for p in probes if p in table.index]
        if present:
            scores.append(table.loc[present].mean(axis=0))
    if not scores: return None
    return pd.concat(scores, axis=1).mean(axis=1)

mapk_g = gene_score(MAPK, table, PROBE_TO_GENE)
panel_g = gene_score(PANEL, table, PROBE_TO_GENE)
if mapk_g is None or panel_g is None:
    print("  WARN: probe lookup failed, skipping GSE76039")
else:
    # z-score within cohort
    mapk_z = (mapk_g - mapk_g.mean())/mapk_g.std()
    panel_z = (panel_g - panel_g.mean())/panel_g.std()
    df = pd.DataFrame({'sample':mapk_z.index,'mapk':mapk_z.values,'panel':panel_z.values})
    df['hist'] = df['sample'].map(hist_map)
    df['hist_short'] = df['hist'].apply(lambda s: 'ATC' if 'Anaplastic' in str(s) else ('PDTC' if 'Poorly' in str(s) else 'OTHER'))
    print(f"\n  GSE76039 by histology:")
    g = df.groupby('hist_short').agg({'mapk':['mean','std','count'],'panel':['mean','std']})
    print(g.round(3))
    g.to_csv(OUT/"v11_GSE76039_by_histology.tsv", sep="\t")
    
    pdtc = df[df['hist_short']=='PDTC']
    atc = df[df['hist_short']=='ATC']
    print(f"\n  ATC (n={len(atc)}) vs PDTC (n={len(pdtc)}):")
    print(f"    d(ATC−PDTC) MAPK = {cd(atc['mapk'].values, pdtc['mapk'].values):+.3f}")
    print(f"    d(ATC−PDTC) Panel = {cd(atc['panel'].values, pdtc['panel'].values):+.3f}")
    rho, p = stats.spearmanr(df['mapk'], df['panel'])
    print(f"  Within GSE76039 MAPK × Panel ρ={rho:+.3f}, p={p:.3g}")
    df.to_csv(OUT/"v11_GSE76039_sample_scores.tsv", sep="\t", index=False)

print("\n=== DONE v11 ===")
