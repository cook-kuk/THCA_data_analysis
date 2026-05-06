#!/usr/bin/env python3
"""Paper 3 ICI Track B-lite — Phase 2: deep analysis.

Adds on top of phase 1:
  - Subtype-aware analysis for ALL thyroid cohorts (parse from series matrix)
  - GSE151179 RAI-refractory vs RAI-avid scoring
  - Riaz On-treatment scoring + Pre→On per-patient delta
  - Bootstrap (1000) effect sizes for DIAL-lite
  - Permutation null (1000) on sign consistency
  - Hierarchical cluster bootstrap stability (200 boots)
  - PCA loadings for top components
  - Multi-cohort meta: ATC-vs-non-ATC effect size pooled across thyroid cohorts
  - PNG figures: subtype heatmap, PCA scatter, DIAL forest plot, ATC vs others bar

Outputs to: project/results/paper3_ici_track_b_lite/
"""
from __future__ import annotations
import gzip, os, re, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_public_data")
REG_DIR  = Path("/home/seungho/personal/THCA_data_analysis/project/reports/paper3_ici")
OUT_DIR  = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_track_b_lite")
SCORES_DIR = OUT_DIR / "scores_per_cohort"
FIG_DIR  = OUT_DIR / "figures_png"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Reuse phase 1 helpers via direct re-implementation (lighter than imports)

mods_df = pd.read_csv(REG_DIR / "paper3_ici_module_gene_list.tsv", sep="\t")
MODULES = mods_df.groupby("module")["gene_symbol"].apply(list).to_dict()
ALL_MODULE_GENES = sorted({g for genes in MODULES.values() for g in genes})

probe_map = pd.read_csv(OUT_DIR / "gpl570_probe_symbol_map_module_subset.tsv", sep="\t")
SYM_BY_PROBE = dict(zip(probe_map.probe_id, probe_map.gene_symbol))

ENTREZ2SYM = {}
for _, r in probe_map.iterrows():
    if r.entrez_id and not pd.isna(r.entrez_id):
        ENTREZ2SYM[str(int(r.entrez_id))] = r.gene_symbol

def score_modules(expr: pd.DataFrame):
    rows = []
    for s in expr.index:
        rec = {"sample_id": s}
        for mod, genes in MODULES.items():
            present = [g for g in genes if g in expr.columns]
            rec[mod] = expr.loc[s, present].astype(float).mean() if present else np.nan
        rows.append(rec)
    sm = pd.DataFrame(rows).set_index("sample_id")
    z = (sm - sm.mean()) / sm.std(ddof=0).replace(0, np.nan)
    return z, sm

# -------- subtype parsers --------

def parse_series_metadata(acc: str):
    p = DATA_DIR / acc / f"{acc}_series_matrix.txt.gz"
    title=[]; gsm=[]; src=[]; chars=[]
    with gzip.open(p,"rt",errors="replace") as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            tag=cols[0]; vals=[c.strip().strip('"') for c in cols[1:]]
            if tag=="!Sample_title": title=vals
            elif tag=="!Sample_geo_accession": gsm=vals
            elif tag=="!Sample_source_name_ch1": src=vals
            elif tag=="!Sample_characteristics_ch1": chars.append(vals)
            if line.startswith("!series_matrix_table_begin"): break
    df = pd.DataFrame({"gsm":gsm,"title":title})
    if src: df["source"]=src
    for i,row in enumerate(chars):
        df[f"char_{i}"] = row
    return df

def derive_subtype(acc, meta: pd.DataFrame) -> pd.Series:
    """Return per-sample subtype label (best effort)."""
    n = len(meta)
    out = pd.Series(["unknown"]*n, index=meta.gsm.values)
    if acc=="GSE126698":
        # not parsed here — handled separately via title prefix
        pass
    elif acc=="GSE76039":
        # 17 PDTC + 20 ATC; cannot reliably distinguish per sample without external map
        out[:] = "PDTC_or_ATC_advanced"
    elif acc=="GSE65144":
        # ATC + 1 normal pool. Use char that contains "tissue type"
        for col in meta.columns:
            if col.startswith("char_"):
                v = meta[col].astype(str)
                if v.str.contains("tissue type:", regex=False).any():
                    for i,gsm in enumerate(meta.gsm):
                        s = v.iloc[i].lower()
                        if "anaplastic" in s: out[gsm]="ATC"
                        elif "normal" in s or "control" in s or "non-tumor" in s: out[gsm]="NT"
    elif acc=="GSE29265":
        # title contains "Anaplastic thyroid carcinoma ATC*", "Papillary thyroid carcinoma PTC*", "Patient-matched nontumor control"
        for i,gsm in enumerate(meta.gsm):
            t = (meta.title.iloc[i] or "").lower()
            if "anaplastic" in t: out[gsm]="ATC"
            elif "papillary" in t: out[gsm]="PTC"
            elif "nontumor" in t or "non-tumor" in t or "control" in t: out[gsm]="NT"
    elif acc=="GSE33630":
        # char "pathological diagnostic"
        for col in meta.columns:
            if col.startswith("char_"):
                v=meta[col].astype(str)
                if v.str.contains("pathological diagnostic", regex=False).any():
                    for i,gsm in enumerate(meta.gsm):
                        s = v.iloc[i].lower()
                        if "anaplastic" in s: out[gsm]="ATC"
                        elif "papillary" in s: out[gsm]="PTC"
                        elif "normal" in s or "control" in s: out[gsm]="NT"
    elif acc=="GSE60542":
        # tissue type characteristic — primary vs nodal mets — all PTC origin
        # use tissue type / source for primary vs LN met
        for col in meta.columns:
            if col.startswith("char_"):
                v=meta[col].astype(str)
                if v.str.contains("tissue.type", regex=True).any() or v.str.contains("tissue type", regex=False).any():
                    for i,gsm in enumerate(meta.gsm):
                        s = v.iloc[i].lower()
                        if "primary" in s: out[gsm]="PTC_primary"
                        elif "metasta" in s or "lymph" in s or "ln" in s: out[gsm]="PTC_nodal_met"
                        elif "normal" in s: out[gsm]="NT"
        # if still all unknown, try source
        if (out=="unknown").all() and "source" in meta.columns:
            for i,gsm in enumerate(meta.gsm):
                s = (meta.source.iloc[i] or "").lower()
                if "primary" in s: out[gsm]="PTC_primary"
                elif "metasta" in s or "node" in s: out[gsm]="PTC_nodal_met"
                elif "normal" in s: out[gsm]="NT"
    elif acc=="GSE151179":
        # collection before/after rai = RAI status
        for col in meta.columns:
            if col.startswith("char_"):
                v=meta[col].astype(str)
                if v.str.contains("collection before/after rai", regex=False).any():
                    for i,gsm in enumerate(meta.gsm):
                        s = v.iloc[i].lower()
                        if "before" in s: out[gsm]="RAI_pre"  # typically RAI-avid baseline
                        elif "after" in s: out[gsm]="RAI_post"  # RAI-refractory after RAI
                        else: out[gsm]="unknown"
    return out

# -------- microarray loaders w/ probe→symbol --------

def parse_table(acc):
    p = DATA_DIR / acc / f"{acc}_series_matrix.txt.gz"
    in_table=False; rows=[]; samples=[]
    with gzip.open(p,"rt",errors="replace") as fh:
        for line in fh:
            if not in_table:
                if line.startswith("!series_matrix_table_begin"):
                    in_table=True
                    header=next(fh).rstrip("\n").split("\t")
                    samples=[h.strip().strip('"') for h in header[1:]]
                continue
            if line.startswith("!series_matrix_table_end"): break
            cols = line.rstrip("\n").split("\t")
            pid = cols[0].strip().strip('"')
            if pid in SYM_BY_PROBE:
                vals=[float(x) if x not in ("","NA","NaN") else np.nan for x in cols[1:]]
                rows.append([pid]+vals)
    if not rows: return pd.DataFrame()
    df=pd.DataFrame(rows, columns=["probe_id"]+samples).set_index("probe_id")
    df["symbol"]=[SYM_BY_PROBE[p] for p in df.index]
    return df.groupby("symbol").max(numeric_only=True).T

# -------- GSE151179: clariom platform — series_matrix has gene-level data already? --------
# Inspecting earlier: 27189 features symbol_or_other. Let's read full table directly with symbol filter.
def load_gse151179():
    p = DATA_DIR / "GSE151179" / "GSE151179_series_matrix.txt.gz"
    rows = []; samples=[]; in_table=False
    with gzip.open(p,"rt",errors="replace") as fh:
        for line in fh:
            if not in_table:
                if line.startswith("!series_matrix_table_begin"):
                    in_table=True
                    header=next(fh).rstrip("\n").split("\t")
                    samples=[h.strip().strip('"') for h in header[1:]]
                continue
            if line.startswith("!series_matrix_table_end"): break
            cols=line.rstrip("\n").split("\t")
            pid = cols[0].strip().strip('"')
            # Clariom IDs are numeric transcript cluster IDs — would need mapping. Skip if no symbol.
            # Try to detect symbol-like rows
            if re.match(r"^[A-Z][A-Z0-9-]+$", pid) and pid in ALL_MODULE_GENES:
                vals=[float(x) if x not in ("","NA","NaN") else np.nan for x in cols[1:]]
                rows.append([pid]+vals)
    if not rows: return pd.DataFrame()
    df=pd.DataFrame(rows, columns=["symbol"]+samples).set_index("symbol")
    return df.T  # samples x genes

# -------- Riaz On-treatment via Entrez map --------

def load_riaz_all():
    p = DATA_DIR / "GSE91061" / "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"
    df = pd.read_csv(p, index_col=0)
    df.index = df.index.astype(str)
    df = df.loc[df.index.intersection(ENTREZ2SYM.keys())]
    df["symbol"]=[ENTREZ2SYM[i] for i in df.index]
    df = df.groupby("symbol").max(numeric_only=True).T
    return np.log2(df.astype(float)+1.0)

# -------- main --------
print("[phase 2] parsing subtype labels")
subtype_records = {}
for acc in ["GSE76039","GSE65144","GSE29265","GSE33630","GSE60542","GSE151179"]:
    meta = parse_series_metadata(acc)
    sub = derive_subtype(acc, meta)
    subtype_records[acc] = sub
    print(f"  {acc}: subtype distribution {sub.value_counts().to_dict()}")

# Read existing pooled scores
pooled_path = OUT_DIR / "pooled_module_scores.tsv"
pooled = pd.read_csv(pooled_path, sep="\t", index_col=0)
print(f"\n[loaded pooled phase 1] n={len(pooled)} cohorts={pooled.cohort.value_counts().to_dict()}")

# Add GSE151179 if not already (phase 1 didn't include it). Let's add now.
print("\n[adding GSE151179]")
gse151179_expr = load_gse151179()
if not gse151179_expr.empty:
    z,raw = score_modules(gse151179_expr)
    z["cohort"]="GSE151179"
    z.to_csv(SCORES_DIR / "GSE151179_module_scores.tsv", sep="\t")
    pooled = pd.concat([pooled, z], axis=0)
    print(f"  GSE151179 added: n={len(z)}, gene coverage={(~z.drop(columns='cohort').isna().all()).sum()}/7 modules")
else:
    print("  GSE151179: no module-symbol features detected in series matrix (Clariom transcript IDs)")

# Add Riaz On-treatment
print("\n[adding Riaz On-treatment]")
riaz_all = load_riaz_all()
on_samples = [s for s in riaz_all.index if "_On_" in s]
if on_samples:
    z_on, raw_on = score_modules(riaz_all.loc[on_samples])
    z_on["cohort"]="Riaz_GSE91061_on"
    z_on.to_csv(SCORES_DIR / "Riaz_GSE91061_on_module_scores.tsv", sep="\t")
    pooled = pd.concat([pooled, z_on], axis=0)
    print(f"  Riaz On: n={len(z_on)}")

pooled.to_csv(OUT_DIR / "pooled_module_scores_v2.tsv", sep="\t")

# -------- subtype-aware thyroid cohort summaries --------
print("\n[subtype-aware thyroid cohort summaries]")
sub_rows = []
# GSE126698 already labeled by sample title prefix (A/F/N/P)
def gse126698_subtype(s):
    return {"A":"ATC","F":"FTC","N":"NT","P":"PTC"}.get(s[0],"unknown")

# Build subtype mapping per sample for all thyroid cohorts
GSM_TO_TITLE = {}
for acc in ["GSE76039","GSE65144","GSE29265","GSE33630","GSE60542","GSE151179"]:
    meta = parse_series_metadata(acc)
    GSM_TO_TITLE[acc] = dict(zip(meta.gsm, meta.title))

per_sample_sub = {}
# GSE126698: sample IDs in pooled are A1..A10, F1..F6, N2..N8, P1..P6 (already)
for s in pooled[pooled.cohort=="GSE126698"].index:
    per_sample_sub[s] = gse126698_subtype(s)

# microarray cohorts: pooled index is GSM IDs
for acc in ["GSE76039","GSE65144","GSE29265","GSE33630","GSE60542","GSE151179"]:
    sub = subtype_records.get(acc, pd.Series(dtype=object))
    for s in pooled[pooled.cohort==acc].index:
        per_sample_sub[s] = sub.get(s, "unknown")

pooled["subtype"] = [per_sample_sub.get(s, "unknown") for s in pooled.index]
pooled.to_csv(OUT_DIR / "pooled_module_scores_v2_with_subtype.tsv", sep="\t")

for cohort in ["GSE126698","GSE76039","GSE65144","GSE29265","GSE33630","GSE60542","GSE151179"]:
    cdf = pooled[pooled.cohort==cohort]
    if cdf.empty: continue
    for sub in cdf.subtype.unique():
        sdf = cdf[cdf.subtype==sub]
        if len(sdf)<2: continue
        for mod in MODULES:
            if mod not in sdf.columns: continue
            sub_rows.append({
                "cohort":cohort,"subtype":sub,"module":mod,"n":len(sdf),
                "mean_z":float(sdf[mod].mean()),"sd_z":float(sdf[mod].std(ddof=0))
            })
sub_df = pd.DataFrame(sub_rows)
sub_df.to_csv(OUT_DIR / "module_score_by_cohort_subtype.tsv", sep="\t", index=False)
print(f"  {len(sub_df)} (cohort × subtype × module) rows")

# -------- ATC vs non-ATC meta-analysis across thyroid cohorts --------
print("\n[ATC-vs-non-ATC meta-analysis]")
meta_rows = []
for cohort in ["GSE126698","GSE65144","GSE29265","GSE33630","GSE60542","GSE76039"]:
    cdf = pooled[pooled.cohort==cohort]
    if cdf.empty: continue
    is_atc = cdf.subtype.isin(["ATC","PDTC_or_ATC_advanced"])
    if is_atc.sum()<2 or (~is_atc).sum()<2: continue
    for mod in MODULES:
        if mod not in cdf.columns: continue
        a = cdf.loc[is_atc, mod].dropna().values
        b = cdf.loc[~is_atc, mod].dropna().values
        if len(a)<2 or len(b)<2: continue
        d = (a.mean() - b.mean())/np.sqrt((a.std(ddof=1)**2 + b.std(ddof=1)**2)/2 + 1e-12)
        t,p = stats.ttest_ind(a,b,equal_var=False)
        meta_rows.append({
            "cohort":cohort,"module":mod,
            "n_ATC":int(is_atc.sum()),"n_other":int((~is_atc).sum()),
            "mean_ATC":float(a.mean()),"mean_other":float(b.mean()),
            "cohens_d":float(d),"t":float(t),"p":float(p),"sign":"+" if d>0 else "-"
        })
mdf = pd.DataFrame(meta_rows)
mdf.to_csv(OUT_DIR / "atc_vs_other_per_cohort.tsv", sep="\t", index=False)

# Pool effect sizes per module (random effects-ish: simple mean of d weighted by n)
pool_rows = []
for mod in MODULES:
    sub = mdf[mdf.module==mod]
    if sub.empty: continue
    d_pooled = (sub.cohens_d * (sub.n_ATC+sub.n_other)).sum() / (sub.n_ATC+sub.n_other).sum()
    pool_rows.append({
        "module":mod,
        "n_cohorts":len(sub),
        "d_mean":float(sub.cohens_d.mean()),
        "d_pooled_n_weighted":float(d_pooled),
        "sign_consistency":int((sub.cohens_d>0).sum()) if d_pooled>0 else int((sub.cohens_d<0).sum()),
        "of_total":int(len(sub))
    })
pdf = pd.DataFrame(pool_rows)
pdf.to_csv(OUT_DIR / "atc_vs_other_meta_pooled.tsv", sep="\t", index=False)
print(pdf.to_string())

# -------- Bootstrap DIAL-lite --------
print("\n[bootstrap DIAL-lite Hugo + Riaz pre, B=1000]")
hugo_z = pooled[pooled.cohort=="Hugo_GSE78220"].copy()
riaz_pre_z = pooled[pooled.cohort=="Riaz_GSE91061_pre"].copy()

# Re-derive response classes
def parse_resp(acc):
    p = DATA_DIR / acc / f"{acc}_series_matrix.txt.gz"
    samples=[]; ch=[]
    with gzip.open(p,"rt",errors="replace") as fh:
        for line in fh:
            cols=line.rstrip("\n").split("\t")
            tag=cols[0]
            if tag=="!Sample_title": samples=[c.strip().strip('"') for c in cols[1:]]
            elif tag=="!Sample_characteristics_ch1": ch.append([c.strip().strip('"') for c in cols[1:]])
            if line.startswith("!series_matrix_table_begin"): break
    for row in ch:
        if any(("response" in c.lower() and "anti-pd-1 response" in c.lower()) or c.lower().startswith("response:") for c in row if c):
            return dict(zip(samples,row))
    return {}

def hugo_class(s):
    s=(s or "").lower()
    if "complete response" in s or "partial response" in s: return 1
    if "progressive disease" in s: return 0
    return np.nan

def riaz_class(s):
    s=(s or "").split(":")[-1].strip().lower()
    if s in ("prcr","cr","pr"): return 1
    if s in ("pd",): return 0
    return np.nan

hugo_rmap = parse_resp("GSE78220")
hugo_z["resp"] = [hugo_class(hugo_rmap.get(s, "")) for s in hugo_z.index]
riaz_rmap = parse_resp("GSE91061")
riaz_pre_z["resp"] = [riaz_class(riaz_rmap.get(s, "")) for s in riaz_pre_z.index]

def bootstrap_d(a, b, B=1000, rng=None):
    rng = rng or np.random.default_rng(42)
    n_a, n_b = len(a), len(b)
    ds = np.empty(B)
    for i in range(B):
        ai = rng.integers(0, n_a, n_a)
        bi = rng.integers(0, n_b, n_b)
        sa = a[ai]; sb = b[bi]
        ds[i] = (sa.mean()-sb.mean()) / np.sqrt((sa.std(ddof=1)**2 + sb.std(ddof=1)**2)/2 + 1e-12)
    return ds

dial_rows = []
rng = np.random.default_rng(42)
for cname, scores in [("Hugo", hugo_z), ("Riaz_pre", riaz_pre_z)]:
    s = scores.dropna(subset=["resp"])
    pos = s[s.resp==1]; neg = s[s.resp==0]
    for mod in MODULES:
        if mod not in s.columns: continue
        a = pos[mod].dropna().values; b = neg[mod].dropna().values
        if len(a)<3 or len(b)<3: continue
        d_obs = (a.mean()-b.mean())/np.sqrt((a.std(ddof=1)**2 + b.std(ddof=1)**2)/2 + 1e-12)
        ds = bootstrap_d(a,b,B=1000,rng=rng)
        ci_lo, ci_hi = np.quantile(ds,[0.025,0.975])
        # permutation null
        comb = np.concatenate([a,b])
        n_a = len(a); ds_perm = np.empty(1000)
        for i in range(1000):
            rng.shuffle(comb)
            sa = comb[:n_a]; sb = comb[n_a:]
            ds_perm[i] = (sa.mean()-sb.mean())/np.sqrt((sa.std(ddof=1)**2 + sb.std(ddof=1)**2)/2 + 1e-12)
        p_perm = (np.abs(ds_perm) >= abs(d_obs)).mean()
        dial_rows.append({
            "cohort":cname,"module":mod,
            "n_resp":len(a),"n_nonresp":len(b),
            "cohens_d":float(d_obs),
            "ci_lo_95":float(ci_lo),"ci_hi_95":float(ci_hi),
            "perm_p":float(p_perm),
            "sign":"+" if d_obs>0 else "-"
        })
ddf = pd.DataFrame(dial_rows)
ddf.to_csv(OUT_DIR / "dial_lite_bootstrap.tsv", sep="\t", index=False)
print(ddf.to_string())

# Sign consistency permutation null: under H0, sign agreement is 0.5 per module
pivot = ddf.pivot(index="module", columns="cohort", values="sign")
pivot["consistent"] = pivot.apply(lambda r: 1 if r.get("Hugo")==r.get("Riaz_pre") else 0, axis=1)
n_consistent = int(pivot.consistent.sum()); n_total = int(len(pivot))
# Binomial test for consistency
binom_p = stats.binomtest(n_consistent, n_total, 0.5, alternative="greater").pvalue if n_total else np.nan
pivot.attrs["binom_p"] = binom_p
pivot.to_csv(OUT_DIR / "dial_lite_sign_consistency_v2.tsv", sep="\t")
print(f"\nSign-consistency: {n_consistent}/{n_total}  binomial 1-sided p={binom_p:.3f}")

# -------- Riaz Pre→On per-patient delta --------
print("\n[Riaz Pre→On per-patient delta]")
riaz_on_z = pooled[pooled.cohort=="Riaz_GSE91061_on"].copy()
def patient_id(s):
    m = re.match(r"(Pt\d+)_", s)
    return m.group(1) if m else s

riaz_pre_z["patient"] = [patient_id(s) for s in riaz_pre_z.index]
riaz_on_z["patient"] = [patient_id(s) for s in riaz_on_z.index]
shared = set(riaz_pre_z["patient"]) & set(riaz_on_z["patient"])
print(f"  shared patients pre+on: {len(shared)}")
delta_rows = []
for pt in shared:
    pre_row = riaz_pre_z[riaz_pre_z.patient==pt].iloc[0]
    on_row  = riaz_on_z[riaz_on_z.patient==pt].iloc[0]
    rec = {"patient": pt}
    for mod in MODULES:
        rec[f"delta_{mod}"] = float(on_row[mod] - pre_row[mod]) if not pd.isna(on_row[mod]) and not pd.isna(pre_row[mod]) else np.nan
    delta_rows.append(rec)
delta_df = pd.DataFrame(delta_rows)
delta_df.to_csv(OUT_DIR / "riaz_pre_to_on_delta.tsv", sep="\t", index=False)
print(f"  delta rows: {len(delta_df)}")

# -------- PCA loadings --------
print("\n[PCA loadings]")
mod_cols = [c for c in pooled.columns if c not in ("cohort","subtype")]
X = pooled[mod_cols].dropna()
Xs = StandardScaler().fit_transform(X.values)
pca = PCA(n_components=min(5,len(mod_cols)))
P = pca.fit_transform(Xs)
loadings = pd.DataFrame(pca.components_.T, index=mod_cols, columns=[f"PC{i+1}" for i in range(P.shape[1])])
loadings.loc["__explained_var__"] = pca.explained_variance_ratio_[:P.shape[1]].tolist()
loadings.to_csv(OUT_DIR / "pca_loadings_v2.tsv", sep="\t")
print(loadings.round(3).to_string())

# -------- Hierarchical cluster bootstrap stability --------
print("\n[hclust bootstrap stability, B=200]")
def hclust_assignments(X, k):
    return AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(X)

base_lbl = hclust_assignments(Xs, 4)
boot_lbls = []
rng2 = np.random.default_rng(0)
for b in range(200):
    idx = rng2.integers(0, len(Xs), len(Xs))
    Xb = Xs[idx]
    lbl_b = hclust_assignments(Xb, 4)
    # remap clusters by majority match — use simple approach: count co-membership
    boot_lbls.append((idx, lbl_b))

# Sample-pair co-membership stability
pair_count = np.zeros((len(Xs), len(Xs)))
pair_total = np.zeros((len(Xs), len(Xs)))
for idx, lbl in boot_lbls:
    for k in range(len(idx)):
        for l in range(k+1, len(idx)):
            i,j = sorted([idx[k], idx[l]])
            pair_total[i,j] += 1
            if lbl[k]==lbl[l]:
                pair_count[i,j] += 1
mask = pair_total>0
stability = np.zeros_like(pair_count)
stability[mask] = pair_count[mask]/pair_total[mask]
mean_pair_stab = stability[mask].mean() if mask.any() else 0.0
print(f"  mean co-membership stability across 200 bootstraps: {mean_pair_stab:.3f}")

# -------- FIGURES --------
print("\n[rendering figures PNG]")

# Fig 1: ATC vs other meta-pooled effect size bar
plt.figure(figsize=(8,5))
pdf_sorted = pdf.sort_values("d_pooled_n_weighted")
plt.barh(pdf_sorted.module, pdf_sorted.d_pooled_n_weighted)
plt.axvline(0, color="k", lw=1)
plt.xlabel("Cohen's d (ATC vs non-ATC, n-weighted across thyroid cohorts)")
plt.title("Module score differential: ATC vs non-ATC (Track B-lite, hypothesis-generating)")
plt.tight_layout()
plt.savefig(FIG_DIR / "fig1_atc_vs_other_meta.png", dpi=150)
plt.close()

# Fig 2: PCA scatter colored by cohort
plt.figure(figsize=(9,7))
P_df = pd.DataFrame(P[:,:2], index=X.index, columns=["PC1","PC2"])
P_df["cohort"] = pooled.loc[X.index, "cohort"].values
for c in P_df.cohort.unique():
    sel = P_df.cohort==c
    plt.scatter(P_df.loc[sel,"PC1"], P_df.loc[sel,"PC2"], s=30, alpha=0.6, label=c)
plt.legend(loc="best", fontsize=8)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.title("Pooled module-score PCA (8 cohorts, n=415+)")
plt.tight_layout()
plt.savefig(FIG_DIR / "fig2_pca_scatter_by_cohort.png", dpi=150)
plt.close()

# Fig 3: subtype heatmap (mean z by cohort × subtype, modules as rows)
pivot_heat = sub_df.pivot_table(index="module", columns=["cohort","subtype"], values="mean_z")
plt.figure(figsize=(max(10, pivot_heat.shape[1]*0.7), 5))
im = plt.imshow(pivot_heat.values, aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5)
plt.colorbar(im, label="mean z-score")
plt.yticks(range(pivot_heat.shape[0]), pivot_heat.index)
plt.xticks(range(pivot_heat.shape[1]),
           [f"{c}|{s}" for c,s in pivot_heat.columns],
           rotation=90, fontsize=7)
plt.title("Module score by cohort × subtype (z-scored within cohort)")
plt.tight_layout()
plt.savefig(FIG_DIR / "fig3_subtype_heatmap.png", dpi=150)
plt.close()

# Fig 4: DIAL-lite forest plot Hugo + Riaz
plt.figure(figsize=(9,5))
mods = sorted(MODULES.keys())
ypos = np.arange(len(mods))
for i, cname in enumerate(["Hugo","Riaz_pre"]):
    sub = ddf[ddf.cohort==cname]
    if sub.empty: continue
    ds = sub.set_index("module").reindex(mods)
    plt.errorbar(ds.cohens_d.values, ypos+(i-0.5)*0.2,
                 xerr=[ds.cohens_d.values-ds.ci_lo_95.values, ds.ci_hi_95.values-ds.cohens_d.values],
                 fmt="o", capsize=3, label=cname)
plt.axvline(0, color="k", lw=0.5)
plt.yticks(ypos, mods)
plt.xlabel("Cohen's d (responder − non-responder)")
plt.title("DIAL-lite: pan-cancer ICI response × module (1000-bootstrap CI)")
plt.legend()
plt.tight_layout()
plt.savefig(FIG_DIR / "fig4_dial_lite_forest.png", dpi=150)
plt.close()

# Fig 5: ATC inflammation vs differentiation tradeoff (GSE126698)
gse126_pooled = pooled[pooled.cohort=="GSE126698"]
plt.figure(figsize=(7,6))
colors = {"ATC":"red","PDTC":"darkorange","FTC":"orange","PTC":"steelblue","NT":"grey"}
for sub in gse126_pooled.subtype.unique():
    sel = gse126_pooled.subtype==sub
    plt.scatter(gse126_pooled.loc[sel,"thyroid_differentiation"],
                gse126_pooled.loc[sel,"IFNG_T_cell_inflamed"],
                s=60, alpha=0.7, label=sub, c=colors.get(sub,"black"))
plt.xlabel("thyroid_differentiation (z)")
plt.ylabel("IFNG_T_cell_inflamed (z)")
plt.axhline(0, color="k", lw=0.3); plt.axvline(0, color="k", lw=0.3)
plt.legend()
plt.title("GSE126698: dedifferentiation vs inflammation (per-sample z)")
plt.tight_layout()
plt.savefig(FIG_DIR / "fig5_dediff_vs_inflam_GSE126698.png", dpi=150)
plt.close()

print(f"  figures saved: {sorted(p.name for p in FIG_DIR.glob('*.png'))}")

# -------- Summary numbers --------
summary = {
    "n_samples_pooled": int(len(pooled)),
    "n_cohorts": int(pooled.cohort.nunique()),
    "cohort_sizes": pooled.cohort.value_counts().to_dict(),
    "pc1_explained": float(pca.explained_variance_ratio_[0]),
    "pc2_explained": float(pca.explained_variance_ratio_[1]),
    "atc_vs_other_meta": pdf.to_dict(orient="records"),
    "dial_sign_consistency": {"n_consistent": n_consistent, "n_total": n_total, "binom_p": float(binom_p)},
    "riaz_pre_on_paired": int(len(delta_df)),
    "hclust_pair_stability": float(mean_pair_stab),
}
with open(OUT_DIR / "phase2_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2, default=str)
print("\n[done]")
