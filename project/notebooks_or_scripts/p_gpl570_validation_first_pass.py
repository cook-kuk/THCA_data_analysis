#!/usr/bin/env python3
"""
GPL570 thyroid dedifferentiation external validation pack — first pass.

Datasets (all GPL570 — Affymetrix HG-U133 Plus 2.0):
  1. GSE33630  Tomas et al  — 49 PTC + 11 ATC + 45 Normal (n=105)
  2. GSE29265  Tomas et al  — 20 PTC +  9 ATC + 20 Normal (n=49)
  3. GSE65144  von Roemeling — 12 ATC + 13 Normal (n=25)

Inputs (read-only):
  project/results/p_gpl570_validation/raw/<ACC>/<ACC>_series_matrix.txt.gz
  project/results/p_landa_2016/raw/GPL570_full.soft           (re-used annotation)

Outputs:
  project/results/p_gpl570_validation/sample_metadata.tsv
  project/results/p_gpl570_validation/<ACC>_expression_gene_log.tsv.gz
  project/results/p_gpl570_validation/<ACC>_scores.tsv
  project/results/p_gpl570_validation/probe_to_gene_panel.tsv
  project/results/p_gpl570_validation/gpl570_score_tests.tsv
  project/results/p_gpl570_validation/gpl570_meta_effect_summary.tsv
  project/results/p_gpl570_validation/gpl570_rai_lineage_boxplots.png
  project/results/p_gpl570_validation/gpl570_dm1_nonoverlap_scatter_grid.png
  project/results/p_gpl570_validation/gpl570_direction_consistency_forest.png

Marathon discipline: CPU-only. No download. No CEL re-normalisation. No TCGA classifier
transfer (within-cohort z-score per dataset only). No survival / mutation / pooled-batch
modelling (separate within-dataset stats only).
"""
import os, sys, gzip, json, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
OUT = ROOT / "project/results/p_gpl570_validation"
RAW = OUT / "raw"
SOFT = ROOT / "project/results/p_landa_2016/raw/GPL570_full.soft"  # reuse Landa 2016 download
RNG_SEED = 42

DATASETS = ["GSE33630", "GSE29265", "GSE65144"]

# ---------- gene panels (identical to GSE76039 first pass) ----------
RAI_8       = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
NONOVERLAP  = ["SLC26A4", "IYD", "DUOX1", "DUOX2", "TFF3", "HHEX", "GLIS3", "DIO2"]
TDS_LIKE    = sorted(set(RAI_8) | set(NONOVERLAP))
TF_COLLAPSE = ["FOXE1", "NKX2-1", "PAX8", "HHEX"]
STAT3_AP1_DNMT = ["STAT3", "FOSL1", "JUNB", "DNMT1", "DNMT3B"]
TACSTD2_GENE = "TACSTD2"
GENE_ALIASES = {"NKX2-1": ["NKX2-1", "TITF1"]}

PANEL_UNION = sorted(set(RAI_8) | set(NONOVERLAP) | set(TF_COLLAPSE) |
                     set(STAT3_AP1_DNMT) | {TACSTD2_GENE})


# ---------- 1. parse Series Matrix ----------
def parse_series_matrix(path):
    meta_lines, table_lines, in_table = [], [], False
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln == "!series_matrix_table_begin": in_table = True; continue
            if ln == "!series_matrix_table_end": in_table = False; continue
            (table_lines if in_table else meta_lines).append(ln)
    sample_id = title = src_name = platform = proc = None
    chars = []
    for ln in meta_lines:
        if ln.startswith("!Sample_geo_accession"):
            sample_id = [x.strip('"') for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_title"):
            title = [x.strip('"') for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_source_name_ch1"):
            src_name = [x.strip('"') for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1"):
            chars.append([x.strip('"') for x in ln.split("\t")[1:]])
        elif ln.startswith("!Sample_platform_id"):
            platform = [x.strip('"') for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_data_processing"):
            proc = [x.strip('"') for x in ln.split("\t")[1:]]
    n = len(sample_id)
    rows = []
    for i in range(n):
        row = {"sample_id": sample_id[i],
               "title": title[i] if title else "",
               "source_name": src_name[i] if src_name else "",
               "platform": platform[i] if platform else "",
               "data_processing": proc[i] if proc else ""}
        for j, c in enumerate(chars):
            v = c[i] if i < len(c) else ""
            if ":" in v:
                k, val = v.split(":", 1)
                row[f"char_{k.strip().replace(' ', '_').lower()}"] = val.strip()
            else:
                row[f"char_extra_{j}"] = v
        rows.append(row)
    meta_df = pd.DataFrame(rows)
    header = [c.strip('"') for c in table_lines[0].split("\t")]
    body = [ln.split("\t") for ln in table_lines[1:]]
    expr_df = pd.DataFrame(body, columns=header).rename(columns={header[0]: "probe_id"})
    expr_df["probe_id"] = expr_df["probe_id"].astype(str).str.strip('"')
    for c in expr_df.columns[1:]:
        expr_df[c] = pd.to_numeric(expr_df[c], errors="coerce")
    return meta_df, expr_df


# ---------- 2. dataset-specific histology assignment ----------
def assign_histology(meta_df, dataset):
    """
    Consistent histology label across datasets:
    Normal | PTC | ATC | UNKNOWN  (PDTC absent in these 3)
    """
    out = meta_df.copy()

    def label_GSE33630(row):
        diag = (row.get("char_pathological_diagnostic", "") or "").lower()
        if "anaplastic" in diag: return "ATC"
        if "papillary" in diag: return "PTC"
        if "non-tumor" in diag or "non tumor" in diag or "normal" in diag: return "Normal"
        return "UNKNOWN"

    def label_GSE29265(row):
        title = (row.get("title", "") or "").lower()
        if "anaplastic" in title: return "ATC"
        if "papillary" in title: return "PTC"
        if "non-tumor" in title or "nontumor" in title or "non tumor" in title or "normal" in title:
            return "Normal"
        return "UNKNOWN"

    def label_GSE65144(row):
        tt = (row.get("char_tissue_type", "") or "").lower()
        src = (row.get("source_name", "") or "").lower()
        if "anaplastic" in tt or "anaplastic" in src: return "ATC"
        if "normal" in tt or "normal" in src: return "Normal"
        return "UNKNOWN"

    fn = {"GSE33630": label_GSE33630, "GSE29265": label_GSE29265, "GSE65144": label_GSE65144}[dataset]
    out["histology"] = out.apply(fn, axis=1)
    out["disease_group"] = out["histology"]
    out["dataset"] = dataset
    return out


# ---------- 3. GPL570 annotation (parse once, reuse) ----------
def parse_gpl570_annotation(soft_path):
    in_tab = False; header = None; rows = []
    with open(soft_path, "rt", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln == "!platform_table_begin": in_tab = True; continue
            if ln == "!platform_table_end": break
            if in_tab:
                if header is None:
                    header = ln.split("\t")
                else:
                    parts = ln.split("\t")
                    if len(parts) < len(header):
                        parts = parts + [""] * (len(header) - len(parts))
                    rows.append(parts[:len(header)])
    df = pd.DataFrame(rows, columns=header)
    keep = ["ID", "Gene Symbol", "ENTREZ_GENE_ID"]
    return df[keep].rename(columns={"ID": "probe_id", "Gene Symbol": "gene_symbol"})


# ---------- 4. probe → gene mapping (best-mean probe collapse) ----------
def map_probes_to_genes(expr_df, ann_df, gene_list):
    sample_cols = [c for c in expr_df.columns if c != "probe_id"]
    ann_df = ann_df.copy()
    ann_df["gene_list"] = ann_df["gene_symbol"].fillna("").apply(
        lambda s: [x.strip() for x in re.split(r"\s*///\s*", s) if x.strip()])
    rows, matched = [], []
    for g in gene_list:
        aliases = GENE_ALIASES.get(g, [g])
        mask = ann_df["gene_list"].apply(lambda gl: any(a in gl for a in aliases))
        probes = ann_df.loc[mask, "probe_id"].tolist()
        e = expr_df[expr_df["probe_id"].isin(probes)]
        if e.empty:
            rows.append({"gene": g, "found": False, "n_probes": 0,
                         **{s: np.nan for s in sample_cols}})
            matched.append({"gene": g, "alias_used": ",".join(aliases),
                            "n_probes": 0, "best_probe": "", "probes_all": ""})
            continue
        means = e[sample_cols].mean(axis=1)
        best_idx = means.idxmax()
        best_row = e.loc[best_idx, sample_cols].to_dict()
        rows.append({"gene": g, "found": True, "n_probes": len(probes), **best_row})
        matched.append({"gene": g, "alias_used": ",".join(aliases),
                        "n_probes": len(probes),
                        "best_probe": e.loc[best_idx, "probe_id"],
                        "probes_all": ",".join(probes)})
    return pd.DataFrame(rows).set_index("gene"), pd.DataFrame(matched)


# ---------- 5. score helpers ----------
def zscore_within_cohort(values):
    v = np.asarray(values, dtype=float)
    m = np.nanmean(v); s = np.nanstd(v, ddof=1)
    if not np.isfinite(s) or s == 0: return np.full_like(v, np.nan)
    return (v - m) / s

def module_score(gx, sample_cols, genes):
    avail = [g for g in genes if g in gx.index and bool(gx.loc[g, "found"])]
    if not avail:
        return pd.Series([np.nan]*len(sample_cols), index=sample_cols), avail
    sub = gx.loc[avail, sample_cols].astype(float)
    z = sub.apply(zscore_within_cohort, axis=1, result_type="broadcast")
    return z.mean(axis=0), avail

def cohens_d(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2: return np.nan
    pooled = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) /
                     (len(a)+len(b)-2))
    if pooled == 0: return np.nan
    return (np.mean(a) - np.mean(b)) / pooled


# ---------- 6. process one dataset ----------
def process_dataset(dataset, ann_df):
    smatrix = RAW / dataset / f"{dataset}_series_matrix.txt.gz"
    print(f"\n========== {dataset} ==========")
    print(f"[load] {smatrix}")
    meta_df, expr_df = parse_series_matrix(smatrix)
    n = len(meta_df)
    sample_cols = [c for c in expr_df.columns if c != "probe_id"]
    assert n == len(sample_cols), (n, len(sample_cols))
    meta_df = assign_histology(meta_df, dataset)
    print(f"  n_samples={n}, histology={meta_df['histology'].value_counts().to_dict()}")

    gx, probe_map = map_probes_to_genes(expr_df, ann_df, PANEL_UNION)
    probe_map["dataset"] = dataset
    print("  per-gene found:", {g: int(bool(gx.loc[g, "found"])) for g in PANEL_UNION})

    gene_out = OUT / f"{dataset}_expression_gene_log.tsv.gz"
    gx[sample_cols].to_csv(gene_out, sep="\t", compression="gzip")

    # scores
    rai8, rai8_used = module_score(gx, sample_cols, RAI_8)
    nonol, nonol_used = module_score(gx, sample_cols, NONOVERLAP)
    tds, tds_used = module_score(gx, sample_cols, TDS_LIKE)
    tfc, tfc_used = module_score(gx, sample_cols, TF_COLLAPSE)
    sad, sad_used = module_score(gx, sample_cols, STAT3_AP1_DNMT)
    tac_z = (pd.Series(zscore_within_cohort(gx.loc[TACSTD2_GENE, sample_cols].values),
                       index=sample_cols)
             if bool(gx.loc[TACSTD2_GENE, "found"])
             else pd.Series(np.nan, index=sample_cols))
    scores = pd.DataFrame({
        "sample_id": sample_cols,
        "RAI_8_score": rai8.values,
        "DM1_like_score": -rai8.values,
        "THYROID_NONOVERLAP_score": nonol.values,
        "TDS_like_score": tds.values,
        "TF_collapse_score": tfc.values,
        "STAT3_AP1_DNMT_score": sad.values,
        "TACSTD2_z": tac_z.values,
    }).merge(meta_df[["sample_id", "title", "histology", "disease_group"]], on="sample_id")
    scores["dataset"] = dataset
    scores.to_csv(OUT / f"{dataset}_scores.tsv", sep="\t", index=False)

    coverage = {"dataset": dataset,
                "RAI_8_n_used": len(rai8_used),
                "NONOVERLAP_n_used": len(nonol_used),
                "TDS_like_n_used": len(tds_used),
                "TF_collapse_n_used": len(tfc_used),
                "STAT3_AP1_DNMT_n_used": len(sad_used),
                "TACSTD2_found": int(bool(gx.loc[TACSTD2_GENE, "found"]))}

    return meta_df, scores, probe_map, coverage


# ---------- 7. tests per dataset ----------
def run_tests(dataset, scores):
    rows = []
    def mw(label, x, y, lx, ly):
        x = np.asarray(x, float); y = np.asarray(y, float)
        x = x[~np.isnan(x)]; y = y[~np.isnan(y)]
        if len(x) < 2 or len(y) < 2:
            rows.append({"dataset": dataset, "test": label, "n_x": len(x), "n_y": len(y),
                         "label_x": lx, "label_y": ly,
                         "mean_x": np.nanmean(x) if len(x) else np.nan,
                         "mean_y": np.nanmean(y) if len(y) else np.nan,
                         "cohens_d": np.nan, "MW_U": np.nan, "MW_p_two_sided": np.nan,
                         "note": "insufficient n"})
            return
        u, p = mannwhitneyu(x, y, alternative="two-sided")
        rows.append({"dataset": dataset, "test": label, "n_x": len(x), "n_y": len(y),
                     "label_x": lx, "label_y": ly,
                     "mean_x": np.mean(x), "mean_y": np.mean(y),
                     "cohens_d": cohens_d(x, y), "MW_U": u, "MW_p_two_sided": p, "note": ""})

    pairs = [("ATC", "PTC"), ("ATC", "Normal"), ("PTC", "Normal"), ("ATC", "PDTC")]
    axes  = ["RAI_8_score", "DM1_like_score", "THYROID_NONOVERLAP_score", "TDS_like_score",
             "TF_collapse_score", "STAT3_AP1_DNMT_score", "TACSTD2_z"]
    for ga, gb in pairs:
        sa = scores[scores["histology"] == ga]
        sb = scores[scores["histology"] == gb]
        if len(sa) < 2 or len(sb) < 2: continue
        for ax in axes:
            mw(f"{ax} :: {ga} vs {gb}", sa[ax].values, sb[ax].values, ga, gb)

    # Spearmans (whole dataset, all available samples)
    def spr(label, a, b):
        if scores[a].isna().all() or scores[b].isna().all():
            rows.append({"dataset": dataset, "test": label, "rho": np.nan, "p": np.nan, "n": 0,
                         "note": "all NaN"})
            return
        rho, p = spearmanr(scores[a], scores[b], nan_policy="omit")
        rows.append({"dataset": dataset, "test": label,
                     "rho": rho, "p": p, "n": int(scores[a].notna().sum()), "note": ""})
    spr("Spearman :: DM1_like vs THYROID_NONOVERLAP", "DM1_like_score", "THYROID_NONOVERLAP_score")
    spr("Spearman :: TF_collapse vs DM1_like",         "TF_collapse_score", "DM1_like_score")
    spr("Spearman :: STAT3_AP1_DNMT vs DM1_like",      "STAT3_AP1_DNMT_score", "DM1_like_score")
    spr("Spearman :: TACSTD2 vs DM1_like",             "TACSTD2_z", "DM1_like_score")
    return rows


# ---------- 8. meta summary ----------
def build_meta_summary(all_tests, all_coverage, all_meta):
    """
    Per-dataset Cohen's d for the headline ATC contrasts and direction-consistency.
    """
    df = pd.DataFrame(all_tests)
    headline = df[df["test"].fillna("").str.contains(":: ATC vs ", regex=False)
                  & df["test"].fillna("").str.contains("RAI_8_score|DM1_like_score|THYROID_NONOVERLAP_score|TDS_like_score|TF_collapse_score|STAT3_AP1_DNMT_score|TACSTD2_z", regex=True)].copy()

    # consistency: lineage axes should be NEGATIVE in ATC contrasts (ATC < comparator);
    # STAT3_AP1_DNMT and TACSTD2_z should be POSITIVE (ATC > comparator).
    expected_sign = {
        "RAI_8_score": -1, "DM1_like_score": +1, "THYROID_NONOVERLAP_score": -1,
        "TDS_like_score": -1, "TF_collapse_score": -1,
        "STAT3_AP1_DNMT_score": +1, "TACSTD2_z": +1,
    }
    def axis_of(test):
        return test.split(" :: ")[0]
    headline["axis"] = headline["test"].map(axis_of)
    headline["expected_sign"] = headline["axis"].map(expected_sign)
    headline["observed_sign"] = np.sign(headline["cohens_d"])
    headline["direction_consistent"] = (
        (headline["observed_sign"] == headline["expected_sign"]) & headline["cohens_d"].notna()
    )
    summary = (headline.groupby(["dataset", "axis", "expected_sign"])
                       .agg(n_contrasts=("cohens_d", "count"),
                            median_d=("cohens_d", "median"),
                            min_d=("cohens_d", "min"),
                            max_d=("cohens_d", "max"),
                            n_consistent=("direction_consistent", "sum"))
                       .reset_index())

    # global summary by axis (across datasets)
    overall = (headline.groupby(["axis", "expected_sign"])
                       .agg(n_contrasts=("cohens_d", "count"),
                            median_d=("cohens_d", "median"),
                            n_datasets=("dataset", "nunique"),
                            n_consistent=("direction_consistent", "sum"))
                       .reset_index())
    overall.insert(0, "scope", "ALL_GPL570")
    return summary, overall, headline


# ---------- 9. figures ----------
def fig_boxplots(all_scores, out_path):
    panels = [("RAI_8_score","RAI_8 (8-gene)"),
              ("DM1_like_score","DM1_like (= -RAI_8)"),
              ("THYROID_NONOVERLAP_score","THYROID_NONOVERLAP"),
              ("TDS_like_score","TDS-like 16-gene"),
              ("TF_collapse_score","TF collapse (FOXE1/NKX2-1/PAX8/HHEX)"),
              ("STAT3_AP1_DNMT_score","STAT3 / AP-1 / DNMT axis"),
              ("TACSTD2_z","TACSTD2 (TROP2) z")]
    histo_order = ["Normal", "PTC", "ATC"]
    colors = {"Normal":"#2ca02c", "PTC":"#1f77b4", "ATC":"#d62728"}
    nrow, ncol = len(panels), 3
    fig, axes = plt.subplots(nrow, ncol,
                              figsize=(3.4*ncol, 2.4*nrow), sharey=False)
    for col_idx, ds in enumerate(DATASETS):
        sub_ds = all_scores[all_scores["dataset"] == ds]
        for row_idx, (col, title) in enumerate(panels):
            ax = axes[row_idx, col_idx]
            data = []
            labels = []
            for h in histo_order:
                vals = sub_ds.loc[sub_ds["histology"] == h, col].dropna().values
                if len(vals) >= 2:
                    data.append(vals); labels.append(f"{h}\n(n={len(vals)})")
            if not data:
                ax.text(0.5, 0.5, "—", ha="center", va="center"); ax.set_axis_off(); continue
            bp = ax.boxplot(data, labels=labels, showmeans=True, widths=0.6, patch_artist=True)
            for patch, lbl in zip(bp["boxes"], labels):
                key = lbl.split("\n")[0]
                patch.set_facecolor(colors.get(key, "lightgray"))
                patch.set_alpha(0.5)
            for i, vals in enumerate(data, start=1):
                jitter = np.random.RandomState(RNG_SEED+row_idx*7+col_idx*3+i).uniform(-0.08, 0.08, len(vals))
                ax.scatter(np.full_like(vals, i, float)+jitter, vals, s=10, alpha=0.7, color="black")
            if col_idx == 0:
                ax.set_ylabel(title, fontsize=9)
            if row_idx == 0:
                ax.set_title(ds, fontsize=10)
            ax.tick_params(axis="x", labelsize=8); ax.tick_params(axis="y", labelsize=8)
    fig.suptitle("GPL570 external validation pack — module scores by dataset (within-cohort z)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(out_path, dpi=160, bbox_inches="tight"); plt.close(fig)


def fig_scatter_grid(all_scores, all_tests, out_path):
    df_t = pd.DataFrame(all_tests)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4))
    colors = {"Normal":"#2ca02c", "PTC":"#1f77b4", "ATC":"#d62728"}
    for ax, ds in zip(axes, DATASETS):
        sub = all_scores[all_scores["dataset"] == ds]
        for h, ssub in sub.groupby("histology"):
            if h not in colors: continue
            ax.scatter(ssub["DM1_like_score"], ssub["THYROID_NONOVERLAP_score"],
                       c=colors[h], s=36, alpha=0.85,
                       label=f"{h} (n={len(ssub)})", edgecolor="white")
        # Spearman label
        rec = df_t[(df_t["dataset"]==ds) & (df_t["test"]=="Spearman :: DM1_like vs THYROID_NONOVERLAP")]
        if len(rec):
            rho = rec.iloc[0]["rho"]; p = rec.iloc[0]["p"]; n = rec.iloc[0]["n"]
            ax.set_title(f"{ds}\nSpearman ρ = {rho:.3f} (p = {p:.2e}, n={int(n)})", fontsize=10)
        else:
            ax.set_title(ds)
        ax.axhline(0, ls=":", c="gray", lw=0.7); ax.axvline(0, ls=":", c="gray", lw=0.7)
        ax.set_xlabel("DM1_like score"); ax.set_ylabel("THYROID_NONOVERLAP score")
        ax.legend(loc="best", fontsize=8)
    fig.suptitle("GPL570 external validation — DM1_like vs zero-overlap lineage score (within-cohort z)", fontsize=11)
    fig.tight_layout(rect=[0,0,1,0.93])
    fig.savefig(out_path, dpi=160, bbox_inches="tight"); plt.close(fig)


def fig_forest(headline, out_path):
    """
    Forest of Cohen's d for headline ATC contrasts, colored by axis, faceted by axis.
    Negative-expected lineage axes shown on left; positive-expected (STAT3/TACSTD2) on right.
    """
    axes_neg = ["RAI_8_score", "DM1_like_score", "THYROID_NONOVERLAP_score", "TDS_like_score", "TF_collapse_score"]
    axes_pos = ["STAT3_AP1_DNMT_score", "TACSTD2_z"]
    order_axes = axes_neg + axes_pos
    contrasts = [("ATC","PTC"), ("ATC","Normal")]
    rows = []
    for ax_name in order_axes:
        for ds in DATASETS:
            for ga, gb in contrasts:
                rec = headline[(headline["axis"]==ax_name) &
                               (headline["dataset"]==ds) &
                               (headline["label_x"]==ga) &
                               (headline["label_y"]==gb)]
                if len(rec):
                    rows.append({"axis": ax_name, "dataset": ds,
                                 "contrast": f"{ga} vs {gb}",
                                 "d": rec.iloc[0]["cohens_d"],
                                 "consistent": bool(rec.iloc[0]["direction_consistent"])})
    df = pd.DataFrame(rows)
    n = len(df)
    fig, ax = plt.subplots(figsize=(8.2, 0.30*n + 1.2))
    ypos = np.arange(n)
    color = ["#2ca02c" if c else "#d62728" for c in df["consistent"].values]
    ax.barh(ypos, df["d"].values, color=color, alpha=0.85, height=0.7)
    ax.axvline(0, color="black", lw=0.8)
    labels = [f"{r['axis']}  |  {r['dataset']}  |  {r['contrast']}" for _, r in df.iterrows()]
    ax.set_yticks(ypos); ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Cohen's d (positive = mean_x > mean_y)")
    consistent = int(df["consistent"].sum()); total = len(df)
    ax.set_title(f"GPL570 external validation — direction consistency forest "
                 f"(green = matches expected sign; {consistent}/{total} consistent)",
                 fontsize=10)
    fig.tight_layout(); fig.savefig(out_path, dpi=160, bbox_inches="tight"); plt.close(fig)


# ---------- main ----------
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"[load] GPL570 annotation (reused) {SOFT}")
    ann_df = parse_gpl570_annotation(SOFT)
    print(f"  n_probes_in_annot={len(ann_df)}")

    all_meta, all_scores, all_probe, all_coverage, all_tests = [], [], [], [], []
    for ds in DATASETS:
        meta, scores, probe_map, cov = process_dataset(ds, ann_df)
        all_meta.append(meta)
        all_scores.append(scores)
        all_probe.append(probe_map)
        all_coverage.append(cov)
        all_tests.extend(run_tests(ds, scores))

    meta_all = pd.concat(all_meta, axis=0, ignore_index=True, sort=False)
    scores_all = pd.concat(all_scores, axis=0, ignore_index=True, sort=False)
    probe_all = pd.concat(all_probe, axis=0, ignore_index=True, sort=False)

    # canonical metadata (subset of cols that are stable across datasets)
    keep = ["dataset", "sample_id", "title", "source_name", "platform",
            "histology", "disease_group"]
    raw_label_col_map = {"GSE33630": "char_pathological_diagnostic",
                         "GSE29265": "title",
                         "GSE65144": "char_tissue_type"}
    meta_all["histology_raw"] = meta_all.apply(
        lambda r: r.get(raw_label_col_map.get(r["dataset"], ""), ""), axis=1)
    keep.insert(3, "histology_raw")
    meta_all["source_file"] = meta_all["dataset"].apply(
        lambda d: f"project/results/p_gpl570_validation/raw/{d}/{d}_series_matrix.txt.gz")
    keep.append("source_file")
    meta_all = meta_all[keep + [c for c in meta_all.columns if c not in keep]]
    meta_all.to_csv(OUT / "sample_metadata.tsv", sep="\t", index=False)
    probe_all.to_csv(OUT / "probe_to_gene_panel.tsv", sep="\t", index=False)

    tests_df = pd.DataFrame(all_tests)
    tests_df.to_csv(OUT / "gpl570_score_tests.tsv", sep="\t", index=False)

    cov_df = pd.DataFrame(all_coverage)
    summary, overall, headline = build_meta_summary(all_tests, all_coverage, all_meta)
    meta_summary = pd.concat(
        [pd.DataFrame([{"scope": "_coverage_per_dataset"}]),
         cov_df,
         pd.DataFrame([{"scope": "_per_dataset_axis_summary"}]),
         summary,
         pd.DataFrame([{"scope": "_overall_per_axis_summary"}]),
         overall,
         pd.DataFrame([{"scope": "_headline_full_table"}]),
         headline],
        axis=0, ignore_index=True, sort=False)
    meta_summary.to_csv(OUT / "gpl570_meta_effect_summary.tsv", sep="\t", index=False)

    # figures
    fig_boxplots(scores_all, OUT / "gpl570_rai_lineage_boxplots.png")
    fig_scatter_grid(scores_all, all_tests, OUT / "gpl570_dm1_nonoverlap_scatter_grid.png")
    fig_forest(headline, OUT / "gpl570_direction_consistency_forest.png")

    # ---------- console summary ----------
    print("\n========== meta-summary (per dataset, ATC contrasts) ==========")
    for ds in DATASETS:
        sub = headline[headline["dataset"] == ds]
        print(f"\n--- {ds} ---")
        for _, r in sub.iterrows():
            mark = "OK" if r["direction_consistent"] else "WRONG"
            print(f"  [{mark:>5}] axis={r['axis']:<25} {r['label_x']} vs {r['label_y']:<8} "
                  f"d={r['cohens_d']:+.3f}  expected_sign={int(r['expected_sign']):+d}")
    print("\n========== overall per-axis (across all datasets) ==========")
    print(overall.to_string(index=False))
    print("\nDONE.")


if __name__ == "__main__":
    np.random.seed(RNG_SEED)
    main()
