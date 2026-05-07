#!/usr/bin/env python3
"""Paper 3 ICI — Track B-lite analysis on already-downloaded public data.

Within marathon mode + Track A freeze constraints (G1/G2 still OPEN, but user
explicitly authorized "Track B" on 2026-05-06). Scope: limited to data already
on disk. NO new downloads of large scRNA RAW.tar. NO LOHHLA. NO NetMHCpan.
NO scVI/scANVI integration. NO publication-quality figures.

Outputs:
  scores_per_cohort/<cohort>_module_scores.tsv
  pooled_module_scores.tsv
  pca_pooled.tsv
  ecotype_sanity_hclust_assignments.tsv
  dial_lite_hugo_riaz.tsv
  module_subtype_summary_thyroid.tsv
"""
from __future__ import annotations
import gzip, io, os, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler

DATA_DIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_public_data")
REG_DIR  = Path("/home/seungho/personal/THCA_data_analysis/project/reports/paper3_ici")
OUT_DIR  = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_track_b_lite")
OUT_DIR.mkdir(parents=True, exist_ok=True)
SCORES_DIR = OUT_DIR / "scores_per_cohort"
SCORES_DIR.mkdir(exist_ok=True)

# ---------------- module gene list ----------------
mods_df = pd.read_csv(REG_DIR / "paper3_ici_module_gene_list.tsv", sep="\t")
MODULES = mods_df.groupby("module")["gene_symbol"].apply(list).to_dict()
ALL_MODULE_GENES = sorted({g for genes in MODULES.values() for g in genes})
print(f"[modules] {len(MODULES)} modules, {len(ALL_MODULE_GENES)} unique genes")

# ---------------- GPL570 probe -> symbol + entrez ----------------
def load_gpl570_subset(genes_of_interest):
    p = DATA_DIR / "_annotations" / "GPL570.annot.gz"
    rows = []
    in_table = False
    with gzip.open(p, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                header = next(fh).rstrip("\n").split("\t")
                col_id = header.index("ID")
                col_sym = header.index("Gene symbol")
                col_eid = header.index("Gene ID")
                continue
            if not in_table: continue
            if line.startswith("!platform_table_end"): break
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= max(col_id,col_sym,col_eid): continue
            sym_field = cols[col_sym]
            # symbol may be "MIR4640///DDR1" — split
            for sym in sym_field.split("///"):
                sym = sym.strip()
                if sym in genes_of_interest:
                    eid = cols[col_eid].split("///")[0].strip() if cols[col_eid] else ""
                    rows.append((cols[col_id], sym, eid))
    df = pd.DataFrame(rows, columns=["probe_id","gene_symbol","entrez_id"])
    return df

probe_map = load_gpl570_subset(set(ALL_MODULE_GENES))
print(f"[gpl570] probe rows for module genes: {len(probe_map)}; unique symbols: {probe_map.gene_symbol.nunique()}")
probe_map.to_csv(OUT_DIR / "gpl570_probe_symbol_map_module_subset.tsv", sep="\t", index=False)

# Build entrez -> symbol for module genes
entrez2sym = {}
for _, r in probe_map.iterrows():
    if r.entrez_id and r.entrez_id != "":
        entrez2sym[str(r.entrez_id)] = r.gene_symbol
print(f"[entrez map] {len(entrez2sym)} entrez->symbol covering modules")

# ---------------- generic module scoring helper ----------------
def score_modules(expr: pd.DataFrame) -> pd.DataFrame:
    """expr: rows=samples, cols=gene symbols. log2 already applied if needed.
       Returns rows=samples, cols=module names, z-scored within cohort."""
    rows = []
    for s in expr.index:
        rec = {"sample_id": s}
        for mod, genes in MODULES.items():
            present = [g for g in genes if g in expr.columns]
            if len(present) == 0:
                rec[mod] = np.nan
            else:
                vals = expr.loc[s, present].astype(float)
                rec[mod] = vals.mean()
        rows.append(rec)
    s = pd.DataFrame(rows).set_index("sample_id")
    # z-score within cohort
    z = (s - s.mean()) / s.std(ddof=0).replace(0, np.nan)
    return z, s  # zscored, raw

# ---------------- Cohort 1: GSE126698 (RNA-seq FPM) ----------------
def load_gse126698():
    p = DATA_DIR / "GSE126698" / "GSE126698_DE_Thyroid_totalRNA_all.csv.gz"
    df = pd.read_csv(p)
    fpm_cols = [c for c in df.columns if c.endswith(".FPM")]
    sub = df[["names"] + fpm_cols].copy()
    sub = sub[sub.names.isin(ALL_MODULE_GENES)]
    sub = sub.groupby("names").max(numeric_only=True)  # max if duplicate symbols
    expr = sub.T  # samples x genes
    expr.index = [c.replace(".FPM","") for c in expr.index]
    expr = np.log2(expr + 1.0)
    return expr

# ---------------- Cohort 2: Hugo GSE78220 (FPKM xlsx) ----------------
def load_hugo():
    import openpyxl
    p = DATA_DIR / "GSE78220" / "GSE78220_PatientFPKM.xlsx"
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    header = list(rows[0])
    body = rows[1:]
    df = pd.DataFrame(body, columns=header)
    df = df[df["Gene"].isin(ALL_MODULE_GENES)]
    df = df.groupby("Gene").max(numeric_only=True)
    expr = df.T
    expr.index = [str(c).replace(".baseline","") for c in expr.index]
    expr = expr.astype(float)
    expr = np.log2(expr + 1.0)
    return expr

# ---------------- Cohort 3: Riaz GSE91061 (FPKM, Entrez IDs) ----------------
def load_riaz_pre():
    p = DATA_DIR / "GSE91061" / "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"
    df = pd.read_csv(p, index_col=0)
    df.index = df.index.astype(str)
    # Map entrez -> symbol (subset)
    df = df.loc[df.index.intersection(entrez2sym.keys())]
    df["symbol"] = [entrez2sym[i] for i in df.index]
    df = df.groupby("symbol").max(numeric_only=True)
    expr_all = df.T.astype(float)
    expr_all = np.log2(expr_all + 1.0)
    # Filter to Pre-treatment only via series matrix
    # Match sample IDs by parsing PtX_Pre_*
    pre = [s for s in expr_all.index if "_Pre_" in s]
    return expr_all.loc[pre]

# ---------------- Cohort 4-8: GPL570 microarray ----------------
def parse_series_matrix_table(path: Path, probe_subset: set) -> pd.DataFrame:
    """Return rows=probe_id (subset only), cols=GSM."""
    in_table = False
    rows = []
    samples = []
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if not in_table:
                if line.startswith("!series_matrix_table_begin"):
                    in_table = True
                    header = next(fh).rstrip("\n").split("\t")
                    samples = [h.strip().strip('"') for h in header[1:]]
                continue
            if line.startswith("!series_matrix_table_end"): break
            cols = line.rstrip("\n").split("\t")
            pid = cols[0].strip().strip('"')
            if pid in probe_subset:
                vals = [float(x) if x not in ("","NA","NaN") else np.nan for x in cols[1:]]
                rows.append([pid]+vals)
    if not rows: return pd.DataFrame()
    df = pd.DataFrame(rows, columns=["probe_id"] + samples).set_index("probe_id")
    return df

def load_microarray(acc: str) -> pd.DataFrame:
    p = DATA_DIR / acc / f"{acc}_series_matrix.txt.gz"
    probe_subset = set(probe_map.probe_id)
    raw = parse_series_matrix_table(p, probe_subset)
    if raw.empty: return raw
    # Map probe -> symbol; aggregate by max across probes per symbol (robust default)
    raw = raw.copy()
    sym_map = dict(zip(probe_map.probe_id, probe_map.gene_symbol))
    raw["symbol"] = [sym_map.get(p,"") for p in raw.index]
    raw = raw[raw["symbol"]!=""]
    agg = raw.groupby("symbol").max(numeric_only=True)
    expr = agg.T  # samples x genes
    # microarray: assume already log-scale (typical RMA). Use as-is.
    return expr

# ---------------- Run all cohorts ----------------
def run_cohort(name, expr):
    if expr is None or expr.empty:
        print(f"  [{name}] empty — skip")
        return None, None
    z, raw = score_modules(expr)
    z = z.assign(cohort=name)
    raw = raw.assign(cohort=name)
    z.to_csv(SCORES_DIR / f"{name}_module_scores.tsv", sep="\t")
    print(f"  [{name}] n_samples={len(z)}  modules_present={z.drop(columns=['cohort']).notna().all().sum()}")
    return z, raw

print("\n[scoring cohorts]")
all_z, all_raw = [], []
for name, loader in [
    ("GSE126698", load_gse126698),
    ("Hugo_GSE78220", load_hugo),
    ("Riaz_GSE91061_pre", load_riaz_pre),
    ("GSE76039", lambda: load_microarray("GSE76039")),
    ("GSE65144", lambda: load_microarray("GSE65144")),
    ("GSE29265", lambda: load_microarray("GSE29265")),
    ("GSE33630", lambda: load_microarray("GSE33630")),
    ("GSE60542", lambda: load_microarray("GSE60542")),
]:
    try:
        expr = loader()
        z, raw = run_cohort(name, expr)
        if z is not None:
            all_z.append(z); all_raw.append(raw)
    except Exception as e:
        print(f"  [{name}] ERROR: {e}")

pooled = pd.concat(all_z, axis=0)
pooled.to_csv(OUT_DIR / "pooled_module_scores.tsv", sep="\t")
print(f"\n[pooled] n_samples={len(pooled)} cohorts={pooled.cohort.value_counts().to_dict()}")

# ---------------- PCA + hierarchical sanity ----------------
print("\n[PCA + hclust sanity]")
mod_cols = [c for c in pooled.columns if c != "cohort"]
X = pooled[mod_cols].dropna()
if len(X) >= 5:
    Xs = StandardScaler().fit_transform(X.values)
    pca = PCA(n_components=min(5, len(mod_cols)))
    P = pca.fit_transform(Xs)
    pca_df = pd.DataFrame(P, index=X.index, columns=[f"PC{i+1}" for i in range(P.shape[1])])
    pca_df["cohort"] = pooled.loc[X.index, "cohort"].values
    pca_df.to_csv(OUT_DIR / "pca_pooled.tsv", sep="\t")
    print(f"  PCA explained_var: {pca.explained_variance_ratio_.round(3).tolist()}")
    # hclust K=4..6
    for k in [3,4,5,6]:
        try:
            ag = AgglomerativeClustering(n_clusters=k, linkage="ward")
            lbl = ag.fit_predict(Xs)
            pca_df[f"hclust_K{k}"] = lbl
        except Exception as e:
            print(f"  hclust K={k} err: {e}")
    pca_df.to_csv(OUT_DIR / "ecotype_sanity_hclust_assignments.tsv", sep="\t")
else:
    print(f"  too few samples for PCA: {len(X)}")

# ---------------- DIAL-lite: Hugo + Riaz response ----------------
print("\n[DIAL-lite Hugo + Riaz]")
def parse_response_meta(acc, key="response"):
    p = DATA_DIR / acc / f"{acc}_series_matrix.txt.gz"
    samples=[]; ch_lines=[]
    with gzip.open(p, "rt", errors="replace") as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            tag = cols[0]
            if tag == '!Sample_title':
                samples = [c.strip().strip('"') for c in cols[1:]]
            elif tag == '!Sample_characteristics_ch1':
                ch_lines.append([c.strip().strip('"') for c in cols[1:]])
            if line.startswith("!series_matrix_table_begin"): break
    # find response row
    for row in ch_lines:
        if any(("response" in c.lower() and "anti-pd-1 response" in c.lower()) or
               c.lower().startswith("response:") for c in row if c):
            return dict(zip(samples, row))
    return {}

def hugo_response_class(s):
    s=(s or "").lower()
    if "complete response" in s or "partial response" in s: return 1
    if "progressive disease" in s: return 0
    if "stable disease" in s: return np.nan  # exclude
    return np.nan

def riaz_response_class(s):
    s=(s or "").split(":")[-1].strip().lower()
    if s in ("prcr","cr","pr"): return 1
    if s in ("pd",): return 0
    return np.nan

hugo_resp_map = parse_response_meta("GSE78220")
hugo_class = {k:hugo_response_class(v) for k,v in hugo_resp_map.items()}
riaz_resp_map = parse_response_meta("GSE91061")
riaz_class = {k:riaz_response_class(v) for k,v in riaz_resp_map.items()}
print(f"  Hugo responders/non-responders/excluded: {sum(v==1 for v in hugo_class.values())}/{sum(v==0 for v in hugo_class.values())}/{sum(np.isnan(v) for v in hugo_class.values())}")
print(f"  Riaz responders/non-responders/excluded: {sum(v==1 for v in riaz_class.values())}/{sum(v==0 for v in riaz_class.values())}/{sum(np.isnan(v) for v in riaz_class.values())}")

# match scores by GSM-aware mapping
# Hugo: scores indexed by Pt1, Pt2... ; map via series matrix sample title
hugo_scores = pooled[pooled.cohort=="Hugo_GSE78220"].drop(columns=["cohort"]).copy()
# hugo_resp_map keys are Pt1, Pt2... and hugo_scores index is also PtX after our cleanup
hugo_scores["response_pos"] = [hugo_class.get(s, np.nan) for s in hugo_scores.index]

# Riaz: index is PtX_Pre_AD..., keys in resp_map same
riaz_scores = pooled[pooled.cohort=="Riaz_GSE91061_pre"].drop(columns=["cohort"]).copy()
riaz_scores["response_pos"] = [riaz_class.get(s, np.nan) for s in riaz_scores.index]

dial_rows = []
for cname, scores in [("Hugo", hugo_scores), ("Riaz", riaz_scores)]:
    s = scores.dropna(subset=["response_pos"])
    if s.empty:
        continue
    pos = s[s.response_pos==1]
    neg = s[s.response_pos==0]
    for mod in MODULES:
        if mod not in s.columns: continue
        a = pos[mod].dropna().values
        b = neg[mod].dropna().values
        if len(a)<2 or len(b)<2: continue
        t,p = stats.ttest_ind(a, b, equal_var=False)
        d = (a.mean() - b.mean()) / np.sqrt((a.std(ddof=1)**2 + b.std(ddof=1)**2)/2 + 1e-12)
        dial_rows.append({
            "cohort": cname, "module": mod,
            "n_responders": len(a), "n_nonresponders": len(b),
            "mean_responders": float(a.mean()), "mean_nonresponders": float(b.mean()),
            "cohens_d": float(d), "t_stat": float(t), "p_value": float(p),
            "sign": "+" if d>0 else "-"
        })
ddf = pd.DataFrame(dial_rows)
ddf.to_csv(OUT_DIR / "dial_lite_hugo_riaz.tsv", sep="\t", index=False)
print(f"  rows: {len(ddf)}")

# Cross-cohort sign consistency
if not ddf.empty:
    pivot = ddf.pivot(index="module", columns="cohort", values="sign")
    pivot["consistent"] = pivot.apply(lambda r: "yes" if r.get("Hugo")==r.get("Riaz") else "no", axis=1)
    pivot.to_csv(OUT_DIR / "dial_lite_sign_consistency.tsv", sep="\t")
    print(pivot.to_string())

# ---------------- Thyroid subtype × module summary ----------------
print("\n[thyroid subtype summary]")
# GSE126698 sample naming: A=ATC, F=FTC, N=NT, P=PTC
def gse126698_subtype(s):
    if s.startswith("A"): return "ATC"
    if s.startswith("F"): return "FTC"
    if s.startswith("N"): return "NT"
    if s.startswith("P"): return "PTC"
    return "Unknown"
ts = []
gse126_z = pooled[pooled.cohort=="GSE126698"].drop(columns=["cohort"]).copy()
gse126_z["subtype"] = [gse126698_subtype(s) for s in gse126_z.index]
for mod in MODULES:
    if mod not in gse126_z.columns: continue
    for sub in ["ATC","PDTC","FTC","NT","PTC"]:
        sub_df = gse126_z[gse126_z.subtype==sub]
        if len(sub_df)>=2:
            ts.append({
                "cohort":"GSE126698","module":mod,"subtype":sub,
                "n":len(sub_df),"mean_z":float(sub_df[mod].mean()),"sd_z":float(sub_df[mod].std(ddof=0))
            })
# GSE76039: per series matrix characteristics — let's derive from sample title
# GSE76039 publication: 17 PDTC + 20 ATC. The series matrix should have subtype
# Defer rigorous subtype parsing — use overall placeholder
tdf = pd.DataFrame(ts)
tdf.to_csv(OUT_DIR / "module_subtype_summary_thyroid.tsv", sep="\t", index=False)
print(tdf.to_string())

print("\n[done]")
