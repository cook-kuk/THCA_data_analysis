#!/usr/bin/env python3
"""Paper 3 ICI Track B-lite — Phase 4 continue.

(1) Build GPL23159 (Clariom S) probe→symbol map for our 60 module genes.
(2) Score GSE151179 modules (n=52), split by RAI status (before vs after RAI).
(3) Score GSE193581 ATC cell-line bulk (n=9 cell lines) — already have ENSG_SYMBOL.
(4) Update pooled scores with both new cohorts.
(5) Refresh figures with RAI overlay.
"""
from __future__ import annotations
import gzip, re, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_public_data")
REG_DIR  = Path("/home/seungho/personal/THCA_data_analysis/project/reports/paper3_ici")
OUT_DIR  = Path("/home/seungho/personal/THCA_data_analysis/project/results/paper3_ici_track_b_lite")
FIG_DIR  = OUT_DIR / "figures_png"

mods_df = pd.read_csv(REG_DIR / "paper3_ici_module_gene_list.tsv", sep="\t")
MODULES = mods_df.groupby("module")["gene_symbol"].apply(list).to_dict()
ALL_MODULE_GENES = sorted({g for genes in MODULES.values() for g in genes})

# (1) Parse GPL23159 platform table — extract probe → symbol for our 60 genes
print("[GPL23159 parse]")
gpl_path = DATA_DIR / "_annotations" / "GPL23159.full.txt"
mod_gene_re = re.compile(r"\b(" + "|".join(re.escape(g) for g in ALL_MODULE_GENES) + r")\b")

probe2sym = {}
in_table = False
header = None
with open(gpl_path, "r", errors="replace") as fh:
    for line in fh:
        if line.startswith("!platform_table_begin"):
            in_table = True; continue
        if line.startswith("!platform_table_end"):
            break
        if not in_table:
            continue
        if header is None:
            header = line.rstrip("\n").split("\t")
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9: continue
        probe_id = cols[0]
        annot_blob = "\t".join(cols[8:])  # SPOT_ID columns concatenated
        # Find module gene symbol mentions; pick the most frequent / first matching
        matches = mod_gene_re.findall(annot_blob)
        if not matches: continue
        # take majority
        from collections import Counter
        sym = Counter(matches).most_common(1)[0][0]
        probe2sym[probe_id] = sym
print(f"  GPL23159 probe→symbol map for {len(probe2sym)} probes covering "
      f"{len(set(probe2sym.values()))} module genes")

probe_map_clariom = pd.DataFrame(
    [(k,v) for k,v in probe2sym.items()], columns=["probe_id","gene_symbol"]
)
probe_map_clariom.to_csv(OUT_DIR / "gpl23159_clariom_probe_symbol_map.tsv", sep="\t", index=False)

# (2) Score GSE151179 — extract probe-level data, map → symbol → module score
print("\n[GSE151179 module scoring]")
def parse_table_w_probe_filter(acc, probe_set):
    p = DATA_DIR / acc / f"{acc}_series_matrix.txt.gz"
    in_table=False; rows=[]; samples=[]
    with gzip.open(p, "rt", errors="replace") as fh:
        for line in fh:
            if not in_table:
                if line.startswith("!series_matrix_table_begin"):
                    in_table=True
                    header=next(fh).rstrip("\n").split("\t")
                    samples=[h.strip().strip('"') for h in header[1:]]
                continue
            if line.startswith("!series_matrix_table_end"): break
            cols=line.rstrip("\n").split("\t")
            pid=cols[0].strip().strip('"')
            if pid in probe_set:
                vals=[float(x) if x not in ("","NA","NaN") else np.nan for x in cols[1:]]
                rows.append([pid]+vals)
    if not rows: return pd.DataFrame()
    df = pd.DataFrame(rows, columns=["probe_id"]+samples).set_index("probe_id")
    return df

gse151179_probe_df = parse_table_w_probe_filter("GSE151179", set(probe2sym.keys()))
print(f"  rows pulled: {len(gse151179_probe_df)} probes × {gse151179_probe_df.shape[1] if not gse151179_probe_df.empty else 0} samples")

if not gse151179_probe_df.empty:
    gse151179_probe_df["symbol"] = [probe2sym[p] for p in gse151179_probe_df.index]
    gse151179_expr = gse151179_probe_df.groupby("symbol").max(numeric_only=True).T  # samples × genes
    print(f"  expr shape: {gse151179_expr.shape} (samples × genes)")
    print(f"  module gene coverage: {gse151179_expr.shape[1]}/60")

    def score_modules(expr):
        rows=[]
        for s in expr.index:
            rec={"sample_id":s}
            for mod, genes in MODULES.items():
                present=[g for g in genes if g in expr.columns]
                rec[mod] = expr.loc[s, present].astype(float).mean() if present else np.nan
            rows.append(rec)
        sm=pd.DataFrame(rows).set_index("sample_id")
        z=(sm-sm.mean())/sm.std(ddof=0).replace(0,np.nan)
        return z

    z151179 = score_modules(gse151179_expr)

    # parse RAI status from series matrix
    p = DATA_DIR / "GSE151179" / "GSE151179_series_matrix.txt.gz"
    samples=[]; chars=[]
    with gzip.open(p,"rt",errors="replace") as fh:
        for line in fh:
            cols=line.rstrip("\n").split("\t")
            tag=cols[0]
            if tag=="!Sample_geo_accession": samples=[c.strip().strip('"') for c in cols[1:]]
            elif tag=="!Sample_characteristics_ch1": chars.append([c.strip().strip('"') for c in cols[1:]])
            if line.startswith("!series_matrix_table_begin"): break

    def find_field(needle):
        for row in chars:
            if any(needle in (c or "").lower() for c in row):
                return dict(zip(samples, row))
        return {}

    rai_map = find_field("collection before/after rai")
    tissue_map = find_field("tissue type")
    histo_map = find_field("histological variant")

    # Also check for "outcome" or "rai-refractory"
    refractory_map = {}
    for row in chars:
        if any("refractory" in (c or "").lower() or "rai" in (c or "").lower() for c in row):
            refractory_map = dict(zip(samples, row))
            break

    z151179["rai_collection"] = [rai_map.get(s,"") for s in z151179.index]
    z151179["tissue_type"] = [tissue_map.get(s,"") for s in z151179.index]
    z151179["histo_variant"] = [histo_map.get(s,"") for s in z151179.index]
    z151179["cohort"] = "GSE151179"
    z151179.to_csv(OUT_DIR / "scores_per_cohort/GSE151179_module_scores.tsv", sep="\t")

    # RAI before vs after
    print(f"\n  RAI distribution: {z151179.rai_collection.value_counts().to_dict()}")
    print(f"  tissue type: {z151179.tissue_type.value_counts().to_dict()}")

    # Per-module RAI before vs after Cohen's d
    rai_pre = z151179[z151179.rai_collection.str.lower().str.contains("before",na=False)]
    rai_post = z151179[z151179.rai_collection.str.lower().str.contains("after",na=False)]
    rai_rows=[]
    for mod in MODULES:
        if mod not in z151179.columns: continue
        a=rai_pre[mod].dropna().values; b=rai_post[mod].dropna().values
        if len(a)<2 or len(b)<2: continue
        d=(b.mean()-a.mean())/np.sqrt((a.std(ddof=1)**2+b.std(ddof=1)**2)/2+1e-12)
        t,p=stats.ttest_ind(b,a,equal_var=False)
        rai_rows.append({
            "module":mod,"n_before":len(a),"n_after":len(b),
            "mean_before":float(a.mean()),"mean_after":float(b.mean()),
            "cohens_d_after_minus_before":float(d),"t":float(t),"p":float(p),
            "sign":"+" if d>0 else "-"
        })
    rai_df = pd.DataFrame(rai_rows)
    rai_df.to_csv(OUT_DIR / "gse151179_rai_after_vs_before.tsv", sep="\t", index=False)
    print("\n  RAI after vs before per-module:")
    print(rai_df.round(3).to_string())

# (3) Score GSE193581 ATC cell line bulk — ENSG_SYMBOL combined IDs
print("\n[GSE193581 ATC cell-line bulk module scoring]")
p = DATA_DIR / "GSE193581" / "GSE193581_bulkcellline.gz"
df = pd.read_csv(p, sep="\t", index_col=0)
# Index format: ENSG00000000003_TSPAN6
sym = df.index.to_series().str.split("_").str[1]
df["symbol"] = sym.values
df_mod = df[df.symbol.isin(ALL_MODULE_GENES)].drop(columns=[c for c in df.columns if c=="symbol"])
df_mod["symbol"] = sym.loc[df_mod.index].values
df_mod = df_mod.groupby("symbol").max(numeric_only=True)
print(f"  ATC cell-line genes covered: {df_mod.shape[0]}/60")
expr_cl = df_mod.T.astype(float)  # samples (P1..P9) × genes
expr_cl = np.log2(expr_cl + 1.0)
def score_modules2(expr):
    rows=[]
    for s in expr.index:
        rec={"sample_id":s}
        for mod, genes in MODULES.items():
            present=[g for g in genes if g in expr.columns]
            rec[mod] = expr.loc[s, present].astype(float).mean() if present else np.nan
        rows.append(rec)
    sm=pd.DataFrame(rows).set_index("sample_id")
    z=(sm-sm.mean())/sm.std(ddof=0).replace(0,np.nan)
    return z
z_cl = score_modules2(expr_cl)
z_cl["cohort"] = "GSE193581_atc_cellline"
z_cl.to_csv(OUT_DIR / "scores_per_cohort/GSE193581_cellline_module_scores.tsv", sep="\t")
print(f"  cell-line samples: {len(z_cl)}")
print(z_cl.drop(columns=["cohort"]).round(2).to_string())

# (4) Update pooled scores
print("\n[updating pooled v3]")
pooled = pd.read_csv(OUT_DIR / "pooled_module_scores_v2_with_subtype.tsv", sep="\t", index_col=0)
# Add GSE151179 rows
if not gse151179_probe_df.empty:
    add = z151179.copy()
    add["subtype"] = add.rai_collection.apply(
        lambda v: "RAI_before" if "before" in v.lower() else ("RAI_after" if "after" in v.lower() else "unknown"))
    add = add[["HLA_class_I","HLA_class_II","IFNG_T_cell_inflamed","TLS_CXCL13_like",
               "checkpoint_exhaustion","myeloid_suppressive","thyroid_differentiation",
               "cohort","subtype"]]
    pooled = pd.concat([pooled, add], axis=0)

# Add GSE193581 cell line rows
add2 = z_cl.copy()
add2["subtype"] = "ATC_cell_line"
add2 = add2[pooled.columns.intersection(add2.columns).tolist()]
# ensure same column order
add2 = z_cl.assign(subtype="ATC_cell_line")[
    ["HLA_class_I","HLA_class_II","IFNG_T_cell_inflamed","TLS_CXCL13_like",
     "checkpoint_exhaustion","myeloid_suppressive","thyroid_differentiation","cohort","subtype"]]
pooled = pd.concat([pooled, add2], axis=0)
pooled.to_csv(OUT_DIR / "pooled_module_scores_v3.tsv", sep="\t")
print(f"  pooled v3 n_samples = {len(pooled)}; cohorts = {pooled.cohort.value_counts().to_dict()}")

# (5) RAI figure
print("\n[fig7 RAI before/after]")
plt.figure(figsize=(8,5))
mods = list(MODULES.keys())
ypos = np.arange(len(mods))
if not rai_df.empty:
    rai_d = rai_df.set_index("module").reindex(mods)
    plt.barh(mods, rai_d.cohens_d_after_minus_before,
             color=["red" if x>0 else "steelblue" for x in rai_d.cohens_d_after_minus_before])
    plt.axvline(0, color="k", lw=0.8)
    plt.xlabel("Cohen's d (post-RAI − pre-RAI)")
    plt.title(f"GSE151179 — RAI-refractory progression: post-RAI vs pre-RAI module shift\n"
              f"(n_before={rai_df.n_before.iloc[0]}, n_after={rai_df.n_after.iloc[0]})")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig7_rai_after_vs_before.png", dpi=150)
    plt.close()
    print("  fig7 saved")

# Cell-line vs tumor compare for GSE193581
print("\n[fig8 GSE193581 cell-line vs GSE126698 ATC tumor]")
gse126_atc = pooled[(pooled.cohort=="GSE126698") & (pooled.subtype=="ATC")]
cl = pooled[pooled.cohort=="GSE193581_atc_cellline"]
plt.figure(figsize=(8,5))
mod_means_atc = gse126_atc[mods].mean()
mod_means_cl  = cl[mods].mean()
x = np.arange(len(mods))
width = 0.35
plt.bar(x-width/2, mod_means_atc.values, width, label=f"GSE126698 ATC tumor (n={len(gse126_atc)})", color="#9b7cff")
plt.bar(x+width/2, mod_means_cl.values, width, label=f"GSE193581 ATC cell line (n={len(cl)})", color="#65a9ff")
plt.axhline(0, color="k", lw=0.5)
plt.xticks(x, mods, rotation=45, ha="right", fontsize=9)
plt.ylabel("mean within-cohort z-score")
plt.title("ATC tumor vs ATC cell line module score profile")
plt.legend()
plt.tight_layout()
plt.savefig(FIG_DIR / "fig8_atc_tumor_vs_cellline.png", dpi=150)
plt.close()
print("  fig8 saved")

# Updated montage with 8 panels
print("\n[updated 8-panel montage]")
import matplotlib.image as mpimg
fig_files = [
    ("fig1_atc_vs_other_meta.png", "Fig 1 — ATC vs non-ATC pooled meta"),
    ("fig2_pca_scatter_by_cohort.png", "Fig 2 — Pooled PCA"),
    ("fig3_subtype_heatmap.png", "Fig 3 — Subtype heatmap"),
    ("fig4_dial_lite_forest.png", "Fig 4 — DIAL-lite forest"),
    ("fig5_dediff_vs_inflam_GSE126698.png", "Fig 5 — Dedifferentiation × inflammation"),
    ("fig6_riaz_pre_to_on_delta.png", "Fig 6 — Riaz Pre→On"),
    ("fig7_rai_after_vs_before.png", "Fig 7 — RAI post vs pre"),
    ("fig8_atc_tumor_vs_cellline.png", "Fig 8 — ATC tumor vs cell line"),
]
fig, axes = plt.subplots(2,4, figsize=(22,10))
for ax, (fname, title) in zip(axes.flat, fig_files):
    p = FIG_DIR / fname
    if p.exists():
        img = mpimg.imread(p)
        ax.imshow(img)
    ax.set_title(title, fontsize=9)
    ax.axis("off")
plt.suptitle("Paper 3 ICI Track B-lite — 8-panel overview (hypothesis-generating)", fontsize=11, y=0.995)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig_montage_8panel.png", dpi=130, bbox_inches="tight")
plt.close()
print("  fig_montage_8panel.png saved")

# Summary
summary = {
    "phase4_cohorts_added": ["GSE151179 (RAI axis)", "GSE193581 ATC cell-line bulk"],
    "gse151179_module_gene_coverage": int(gse151179_expr.shape[1]) if not gse151179_probe_df.empty else 0,
    "gse151179_n": int(len(z151179)) if not gse151179_probe_df.empty else 0,
    "gse151179_rai_distribution": z151179.rai_collection.value_counts().to_dict() if not gse151179_probe_df.empty else {},
    "gse193581_cellline_n": int(len(z_cl)),
    "gse193581_cellline_module_coverage": int(df_mod.shape[0]),
    "pooled_v3_n_samples": int(len(pooled)),
    "pooled_v3_cohorts": pooled.cohort.value_counts().to_dict(),
}
with open(OUT_DIR / "phase4_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2, default=str)
print(f"\n[done]\nSummary: {json.dumps(summary, indent=2, default=str)}")
