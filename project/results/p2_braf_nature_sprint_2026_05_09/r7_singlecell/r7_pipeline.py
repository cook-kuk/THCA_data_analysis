#!/usr/bin/env python3
"""R7 — single-cell projection of HT-13 / FA-12 / MAPK-9 / PANEL-8 axes
on Lu 2023 GSE193581 thyroid scRNA atlas (n=67,678 cells, 23 samples).

Strategy:
  1. Load celltype_annotation.txt (sample × barcode → celltype).
  2. For each per-sample raw UMI file, pull rows for panel genes only,
     log-normalize per cell using total UMI library size (CP10k + log1p).
  3. Stack into one panel-gene × all-cells matrix; compute per-cell
     panel scores (mean log-normalized expression).
  4. Aggregate per cell type (median, IQR, mean).
  5. Decompose: per-cohort axis-score = sum_{ct} fraction_ct * mean_score_ct.
  6. LR-pair re-mining within HT-axis-high vs HT-axis-low spots — done at
     STAGE level using existing PantheonOS LR recurrence files.

Outputs in r7_singlecell/.
"""
from __future__ import annotations
import os, glob, gzip, json, time
import numpy as np
import pandas as pd

OUT = "/home/seungho/personal/THCA_data_analysis/project/results/p2_braf_nature_sprint_2026_05_09/r7_singlecell"
RAW_DIR = "/data/thca/v17_lu2023_GSE193581/raw"
ANN = "/data/thca/v17_lu2023_GSE193581/celltype_annotation.txt.gz"
PANEL_SCORES = "/data/thca/repo_results/v17_lu2023/GSE193581_cell_panel_score.tsv"
LR_RECUR = "/home/seungho/personal/THCA_data_analysis/project/results/pantheonos_demo/_AGGREGATE/lr_condition_recurrence.csv"
LR_DELTA = "/home/seungho/personal/THCA_data_analysis/project/results/pantheonos_demo/_AGGREGATE/cross_stage_delta.csv"

os.makedirs(OUT, exist_ok=True)

# Panels — broad, literature-grounded sets
HT13   = ["HLA-DRA","HLA-DRB1","HLA-DPB1","HLA-DPA1","HLA-DQA1","HLA-DQB1","CD74","CXCL13","CD79A","MS4A1","CD3D","CD3E","GZMK"]
FA12   = ["FABP4","SCD","FASN","ACACA","ELOVL6","PLIN2","APOE","APOC1","LIPE","ACSL1","ACADM","HADHA"]
MAPK9  = ["DUSP4","DUSP5","DUSP6","SPRY1","SPRY2","SPRY4","ETV4","ETV5","PHLDA1"]
PANEL8 = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
PANELS = {"HT13": HT13, "FA12": FA12, "MAPK9": MAPK9, "PANEL8": PANEL8}
ALL_GENES = sorted(set(HT13 + FA12 + MAPK9 + PANEL8))

t0 = time.time()
print(f"[r7] panels: HT13={len(HT13)} FA12={len(FA12)} MAPK9={len(MAPK9)} PANEL8={len(PANEL8)} unique={len(ALL_GENES)}")

# ------- 1. cell-type annotation
# File header is misaligned: 2-col header but 3-col data => cell_id became index,
# 'sample ID' actually holds sample name, 'celltype' holds celltype.
ann = pd.read_csv(ANN, sep="\t")
ann.index.name = "cell_id"
ann = ann.rename(columns={"sample ID": "sample"})
ann.index = ann.index.astype(str)
print(f"[r7] annotation rows={len(ann)}; unique celltypes={ann['celltype'].nunique()}")
print(ann["celltype"].value_counts().to_dict())

# ------- 2. per-sample raw UMI -> panel-gene expression (cache to parquet)
CACHE = os.path.join(OUT, "_panel_expr_cache.parquet")
if os.path.exists(CACHE):
    print(f"[r7] loading cached panel expr: {CACHE}")
    expr_df = pd.read_parquet(CACHE)
    genes_present_all = set([c for c in expr_df.columns if expr_df[c].notna().any()])
    print(f"[r7] cached expr matrix shape: {expr_df.shape}; genes_present={len(genes_present_all)}")
else:
    files = sorted(glob.glob(os.path.join(RAW_DIR, "GSM*_UMI.txt.gz")))
    print(f"[r7] raw UMI files: {len(files)}")
    cells_all = []
    genes_present_all = set()
    panel_mats = []
    for f in files:
        name = os.path.basename(f).replace(".txt.gz", "")
        sample_name = name.split("_", 2)[1]
        df = pd.read_csv(f, sep="\t", index_col=0)
        cells = df.columns.tolist()
        total_umi = df.sum(axis=0).values
        pres = [g for g in ALL_GENES if g in df.index]
        genes_present_all.update(pres)
        sub = df.loc[pres].values.astype(np.float32)
        safe_total = np.where(total_umi == 0, 1, total_umi).astype(np.float32)
        norm = (sub / safe_total[None, :]) * 1e4
        norm = np.log1p(norm)
        panel_mats.append((sample_name, cells, pres, norm))
        print(f"  [{sample_name}] cells={len(cells)} panel_found={len(pres)} elapsed={time.time()-t0:.1f}s")

    print(f"[r7] union genes present: {len(genes_present_all)}/{len(ALL_GENES)}")
    all_cell_ids, all_data = [], {g: [] for g in ALL_GENES}
    for sample_name, cells, pres, norm in panel_mats:
        n = len(cells)
        all_cell_ids.extend(cells)
        pres_set = set(pres)
        pres_idx = {g: i for i, g in enumerate(pres)}
        for g in ALL_GENES:
            if g in pres_set:
                all_data[g].append(norm[pres_idx[g]])
            else:
                all_data[g].append(np.full(n, np.nan, dtype=np.float32))
    expr_df = pd.DataFrame({g: np.concatenate(all_data[g]) for g in ALL_GENES},
                           index=all_cell_ids)
    expr_df.to_parquet(CACHE)
    print(f"[r7] expr matrix shape: {expr_df.shape}; cached -> {CACHE}")

# Merge with annotation
common = expr_df.index.intersection(ann.index)
print(f"[r7] cells with celltype: {len(common)} / {len(expr_df)}")
expr_df = expr_df.loc[common]
celltype = ann.loc[common, "celltype"]
sample = ann.loc[common, "sample"]

# Histology and DM info from existing panel_score (15331 malignant cells covered)
ps = pd.read_csv(PANEL_SCORES, sep="\t", index_col=0)
# For non-malignant cells we infer histology from sample prefix
def hist_from_sample(s):
    if s.startswith("PTC"): return "PTC"
    if s.startswith("ATC"): return "ATC"
    if s.startswith("NORM"): return "NORM"
    return "UNK"
hist = sample.apply(hist_from_sample)

# ------- 3. per-cell panel scores (mean of present, NaN-skip)
scores = {}
for name, genes in PANELS.items():
    cols = [g for g in genes if g in expr_df.columns]
    pres = [g for g in cols if expr_df[g].notna().any()]
    sub = expr_df[pres]
    mean_score = sub.mean(axis=1, skipna=True)  # per-cell avg log-norm
    scores[name] = mean_score
    print(f"[r7] panel={name} present_in_data={len(pres)}/{len(genes)}: {pres}")

score_df = pd.DataFrame(scores)
score_df["celltype"] = celltype.values
score_df["sample"] = sample.values
score_df["histology"] = hist.values
score_df.to_csv(os.path.join(OUT, "r7_per_cell_panel_scores.tsv"), sep="\t",
                float_format="%.4f")
print(f"[r7] per-cell panel scores written, {score_df.shape}")

# ------- 4. cell-type aggregate (median, IQR, mean, n)
def agg_panel(df, group_col):
    rows = []
    for grp, gd in df.groupby(group_col):
        for panel in PANELS.keys():
            v = gd[panel].dropna().values
            if len(v) < 5:
                continue
            rows.append({
                group_col: grp,
                "panel": panel,
                "n": len(v),
                "median": float(np.median(v)),
                "iqr_low": float(np.percentile(v, 25)),
                "iqr_hi":  float(np.percentile(v, 75)),
                "mean": float(np.mean(v)),
                "std": float(np.std(v)),
            })
    return pd.DataFrame(rows)

ct_agg = agg_panel(score_df, "celltype")
ct_agg.to_csv(os.path.join(OUT, "r7_celltype_panel_scores.tsv"), sep="\t",
              index=False, float_format="%.4f")
print(f"[r7] celltype aggregate -> r7_celltype_panel_scores.tsv ({len(ct_agg)} rows)")
print(ct_agg.pivot(index="celltype", columns="panel", values="median").round(3))

# Per histology × celltype
rows = []
for (hist_v, ct), gd in score_df.groupby(["histology", "celltype"]):
    for panel in PANELS.keys():
        v = gd[panel].dropna().values
        if len(v) < 5:
            continue
        rows.append({
            "histology": hist_v, "celltype": ct, "panel": panel,
            "n": len(v),
            "median": float(np.median(v)),
            "iqr_low": float(np.percentile(v, 25)),
            "iqr_hi":  float(np.percentile(v, 75)),
            "mean": float(np.mean(v)),
        })
hist_ct_agg = pd.DataFrame(rows)
hist_ct_agg.to_csv(os.path.join(OUT, "r7_hist_x_celltype_panel_scores.tsv"),
                   sep="\t", index=False, float_format="%.4f")
print(f"[r7] hist×celltype aggregate -> r7_hist_x_celltype_panel_scores.tsv "
      f"({len(hist_ct_agg)} rows)")

# ------- 5. axis decomposition
# bulk-equivalent score = sum_ct fraction_ct * mean_score_ct
# bulk-by-celltype-contribution = fraction_ct * mean_score_ct (gives % share)
def decompose(df, hist_v):
    sub = df[df["histology"] == hist_v]
    n_total = len(sub)
    rows = []
    for panel in PANELS.keys():
        # only cells with at least one panel gene present
        valid = sub[sub[panel].notna()]
        n_valid = len(valid)
        if n_valid == 0:
            continue
        bulk_mean = float(valid[panel].mean())
        for ct, gd in valid.groupby("celltype"):
            frac = len(gd) / n_valid
            ct_mean = float(gd[panel].mean())
            contribution = frac * ct_mean
            share = contribution / bulk_mean if bulk_mean != 0 else np.nan
            rows.append({
                "histology": hist_v, "panel": panel, "celltype": ct,
                "n_cells_ct": len(gd), "n_cells_hist": n_valid,
                "fraction": round(frac, 4),
                "ct_mean_score": round(ct_mean, 4),
                "bulk_mean_score": round(bulk_mean, 4),
                "contribution": round(contribution, 4),
                "share_of_bulk": round(share, 4),
            })
    return pd.DataFrame(rows)

decomp_all = pd.concat([decompose(score_df, h) for h in ["NORM", "PTC", "ATC"]],
                        ignore_index=True)
decomp_all.to_csv(os.path.join(OUT, "r7_axis_decomposition.tsv"),
                  sep="\t", index=False, float_format="%.4f")
print(f"[r7] axis decomposition -> r7_axis_decomposition.tsv ({len(decomp_all)} rows)")

# Quick view: which celltype dominates each panel in PTC?
for hist_v in ["NORM", "PTC", "ATC"]:
    print(f"\n  -- {hist_v} contribution share (top per panel) --")
    sub = decomp_all[decomp_all["histology"] == hist_v]
    for panel in PANELS.keys():
        ps = sub[sub["panel"] == panel].sort_values("share_of_bulk", ascending=False).head(3)
        if len(ps) == 0: continue
        print(f"    {panel}: " + " | ".join(
            f"{r.celltype}={r.share_of_bulk:.2f}" for r in ps.itertuples()))

# ------- 6. LR-pair re-mining: HT-axis-high stages vs low
# We re-use existing per-stage LR recurrence (PT/PTC/LPTC/ATC).
# HT-axis is highest in PTC + LPTC (B/T cell density), low in ATC (myeloid / loss)
# Compute delta_ht_high = mean(PTC,LPTC) - mean(PT,ATC)
lr = pd.read_csv(LR_RECUR)
ht_high_mean = lr[["PTC_recurrence", "LPTC_recurrence"]].mean(axis=1)
ht_low_mean  = lr[["PT_recurrence", "ATC_recurrence"]].mean(axis=1)
lr["delta_ht_high_minus_low"] = ht_high_mean - ht_low_mean
lr["ht_high_mean"] = ht_high_mean
lr["ht_low_mean"]  = ht_low_mean
top_ht = lr.sort_values("delta_ht_high_minus_low", ascending=False).head(20)
top_low = lr.sort_values("delta_ht_high_minus_low", ascending=True).head(20)
out_lr = pd.concat([
    top_ht.assign(direction="HT_high_enriched"),
    top_low.assign(direction="HT_low_enriched")
])
out_lr.to_csv(os.path.join(OUT, "r7_lr_pairs_HTaxis.tsv"), sep="\t",
              index=False, float_format="%.4f")
print(f"\n[r7] LR pairs -> r7_lr_pairs_HTaxis.tsv ({len(out_lr)} rows)")
print("Top 10 HT-axis-high LR pairs:")
print(top_ht[["lr_pair", "ht_high_mean", "ht_low_mean",
              "delta_ht_high_minus_low"]].head(10).to_string(index=False))

# ------- summary json
summary = {
    "n_cells_total": int(len(score_df)),
    "n_celltypes": int(score_df["celltype"].nunique()),
    "celltypes": sorted(score_df["celltype"].unique().tolist()),
    "histology_counts": score_df["histology"].value_counts().to_dict(),
    "panels_present_in_raw": {
        p: [g for g in genes if g in genes_present_all]
        for p, genes in PANELS.items()
    },
    "panels_total": {p: len(g) for p, g in PANELS.items()},
    "celltype_panel_medians": ct_agg.pivot(index="celltype", columns="panel",
                                            values="median").round(3).to_dict(),
    "elapsed_sec": round(time.time() - t0, 1),
}
with open(os.path.join(OUT, "r7_summary.json"), "w") as fh:
    json.dump(summary, fh, indent=2)
print(f"\n[r7] DONE in {time.time()-t0:.1f}s -> {OUT}")
