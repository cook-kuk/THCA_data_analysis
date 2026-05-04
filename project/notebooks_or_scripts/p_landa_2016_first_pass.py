#!/usr/bin/env python3
"""
GSE76039 (Landa 2016 PDTC/ATC) first-pass DM1/RAI-lineage axis check.

Inputs (read-only):
  project/results/p_landa_2016/raw/GSE76039_series_matrix.txt.gz   (gcRMA log2)
  project/results/p_landa_2016/raw/GPL570_full.soft                (probe annotation)

Outputs:
  project/results/p_landa_2016/sample_metadata.tsv
  project/results/p_landa_2016/expression_matrix_log_gene.tsv.gz
  project/results/p_landa_2016/probe_to_gene_panel.tsv
  project/results/p_landa_2016/scores.tsv
  project/results/p_landa_2016/score_tests.tsv
  project/results/p_landa_2016/gse76039_dm1_lineage_boxplot.png
  project/results/p_landa_2016/gse76039_dm1_vs_nonoverlap_scatter.png
  project/results/p_landa_2016/gse76039_mechanism_heatmap.png

Marathon discipline: CPU-only, no download, no GPU/RunPod. No raw FASTQ alignment, no CEL re-normalisation. Uses only the gcRMA-normalised Series Matrix produced by the Landa lab.
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
RAW = ROOT / "project/results/p_landa_2016/raw"
OUT = ROOT / "project/results/p_landa_2016"
SOFT = RAW / "GPL570_full.soft"
SMATRIX = RAW / "GSE76039_series_matrix.txt.gz"
RNG_SEED = 42

# ---------- gene panels (per user prompt) ----------
RAI_8       = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
NONOVERLAP  = ["SLC26A4", "IYD", "DUOX1", "DUOX2", "TFF3", "HHEX", "GLIS3", "DIO2"]
TDS_LIKE    = sorted(set(RAI_8) | set(NONOVERLAP))   # broader thyroid-lineage
TF_COLLAPSE = ["FOXE1", "NKX2-1", "PAX8", "HHEX"]    # all expected DOWN in advanced disease (as in DM2)
STAT3_AP1_DNMT = ["STAT3", "FOSL1", "JUNB", "DNMT1", "DNMT3B"]
TACSTD2_GENE = "TACSTD2"

# NKX2-1 alias on older Affymetrix annotation:
GENE_ALIASES = {"NKX2-1": ["NKX2-1", "TITF1"]}


# ---------- 1. parse Series Matrix metadata + expression ----------
def parse_series_matrix(path):
    meta_lines = []
    table_lines = []
    in_table = False
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln == "!series_matrix_table_begin":
                in_table = True
                continue
            if ln == "!series_matrix_table_end":
                in_table = False
                continue
            if in_table:
                table_lines.append(ln)
            else:
                meta_lines.append(ln)
    # build sample metadata
    sample_id = None
    title = None
    src_name = None
    chars = []  # list of lists, one per characteristics line
    platform = None
    proc = None
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
        row = {
            "sample_id": sample_id[i],
            "title": title[i] if title else "",
            "source_name": src_name[i] if src_name else "",
            "platform": platform[i] if platform else "",
            "data_processing": proc[i] if proc else "",
        }
        # parse characteristics rows ("key: value")
        for j, c in enumerate(chars):
            v = c[i] if i < len(c) else ""
            if ":" in v:
                k, val = v.split(":", 1)
                row[f"char_{k.strip().replace(' ', '_').lower()}"] = val.strip()
            else:
                row[f"char_extra_{j}"] = v
        rows.append(row)
    meta_df = pd.DataFrame(rows)
    # convenience histology / disease_group column
    def histology_of(src):
        s = src.lower()
        if "anaplastic" in s:
            return "ATC"
        if "poorly" in s or "pdtc" in s:
            return "PDTC"
        if "papillary" in s or "ptc" in s:
            return "PTC"
        if "normal" in s:
            return "Normal"
        return "UNKNOWN"
    meta_df["histology"] = meta_df["source_name"].map(histology_of)
    meta_df["disease_group"] = meta_df["histology"]  # synonym for downstream

    # expression table
    header = table_lines[0].split("\t")
    header = [c.strip('"') for c in header]
    body_rows = [ln.split("\t") for ln in table_lines[1:]]
    expr_df = pd.DataFrame(body_rows, columns=header)
    expr_df = expr_df.rename(columns={header[0]: "probe_id"})
    # strip quotes from probe_id values (Series Matrix encloses both header and ID_REF column in quotes)
    expr_df["probe_id"] = expr_df["probe_id"].astype(str).str.strip('"')
    # numeric cast for sample columns
    for c in expr_df.columns[1:]:
        expr_df[c] = pd.to_numeric(expr_df[c], errors="coerce")
    return meta_df, expr_df


# ---------- 2. parse GPL570 annotation: probe -> gene symbol ----------
def parse_gpl570_annotation(soft_path):
    in_tab = False
    header = None
    rows = []
    with open(soft_path, "rt", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln == "!platform_table_begin":
                in_tab = True
                continue
            if ln == "!platform_table_end":
                in_tab = False
                break
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


# ---------- 3. map probes to gene panels; collapse multiple probes per gene by max mean ----------
def map_probes_to_genes(expr_df, ann_df, gene_list):
    """Return dataframe: rows=samples, cols=genes (or NaN where missing)."""
    # all probes that match any of the gene symbols (handle multiple genes per probe joined by ' /// ')
    needed = set()
    for g in gene_list:
        for alias in GENE_ALIASES.get(g, [g]):
            needed.add(alias)
    sample_cols = [c for c in expr_df.columns if c != "probe_id"]
    ann_df = ann_df.copy()
    ann_df["gene_list"] = ann_df["gene_symbol"].fillna("").apply(
        lambda s: [x.strip() for x in re.split(r"\s*///\s*", s) if x.strip()]
    )
    rows = []
    matched_table = []
    for g in gene_list:
        aliases = GENE_ALIASES.get(g, [g])
        mask = ann_df["gene_list"].apply(lambda gl: any(a in gl for a in aliases))
        probes_for_g = ann_df.loc[mask, "probe_id"].tolist()
        e = expr_df[expr_df["probe_id"].isin(probes_for_g)]
        if e.empty:
            rows.append({"gene": g, "found": False, "n_probes": 0,
                          **{s: np.nan for s in sample_cols}})
            matched_table.append({"gene": g, "alias_used": None, "n_probes": 0, "probes": ""})
            continue
        # collapse: keep probe with maximum mean expression across samples (standard for Affy)
        means = e[sample_cols].mean(axis=1)
        best_idx = means.idxmax()
        best_row = e.loc[best_idx, sample_cols].to_dict()
        best_probe = e.loc[best_idx, "probe_id"]
        rows.append({"gene": g, "found": True, "n_probes": len(probes_for_g), **best_row})
        matched_table.append({"gene": g,
                                "alias_used": ",".join(aliases),
                                "n_probes": len(probes_for_g),
                                "best_probe": best_probe,
                                "probes_all": ",".join(probes_for_g)})
    gx = pd.DataFrame(rows).set_index("gene")
    return gx, pd.DataFrame(matched_table)


# ---------- 4. score helpers ----------
def zscore_within_cohort(values):
    v = np.asarray(values, dtype=float)
    m = np.nanmean(v); s = np.nanstd(v, ddof=1)
    if not np.isfinite(s) or s == 0:
        return np.full_like(v, np.nan)
    return (v - m) / s

def module_score(gene_expr_df, sample_cols, genes):
    """gene_expr_df: rows=genes, cols=metadata+samples; returns per-sample mean z."""
    available = [g for g in genes if g in gene_expr_df.index and gene_expr_df.loc[g, "found"] == True]
    if not available:
        return pd.Series([np.nan]*len(sample_cols), index=sample_cols), available
    sub = gene_expr_df.loc[available, sample_cols].astype(float)
    z = sub.apply(zscore_within_cohort, axis=1, result_type="broadcast")
    score = z.mean(axis=0)
    return score, available

def cohens_d(a, b):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2: return np.nan
    pooled = np.sqrt(((len(a)-1)*np.var(a, ddof=1) + (len(b)-1)*np.var(b, ddof=1)) / (len(a)+len(b)-2))
    if pooled == 0: return np.nan
    return (np.mean(a) - np.mean(b)) / pooled


# ---------- 5. main ----------
def main():
    print(f"[load] series matrix {SMATRIX}")
    meta_df, expr_df = parse_series_matrix(SMATRIX)
    n_samples = len(meta_df)
    sample_cols = [c for c in expr_df.columns if c != "probe_id"]
    assert n_samples == len(sample_cols), (n_samples, len(sample_cols))
    print(f"  n_samples={n_samples}, n_probes={len(expr_df)}")
    print(f"  histology counts: {meta_df['histology'].value_counts().to_dict()}")

    # save metadata
    meta_out = OUT / "sample_metadata.tsv"
    meta_df.to_csv(meta_out, sep="\t", index=False)
    print(f"[save] {meta_out}")

    # parse annotation
    print(f"[load] GPL570 annotation {SOFT}")
    ann_df = parse_gpl570_annotation(SOFT)
    print(f"  n_probes_in_annot={len(ann_df)}")

    # build per-gene matrix for panel union
    panel_union = sorted(set(RAI_8) | set(NONOVERLAP) | set(TF_COLLAPSE) |
                          set(STAT3_AP1_DNMT) | {TACSTD2_GENE})
    print(f"[map] panel_union genes: {len(panel_union)} → {panel_union}")
    gene_expr_df, probe_map_df = map_probes_to_genes(expr_df, ann_df, panel_union)
    probe_map_out = OUT / "probe_to_gene_panel.tsv"
    probe_map_df.to_csv(probe_map_out, sep="\t", index=False)
    print(f"[save] {probe_map_out}")
    print("  per-gene found:")
    for g in panel_union:
        f = gene_expr_df.loc[g, "found"]
        n = gene_expr_df.loc[g, "n_probes"]
        print(f"    {g:>10s}  found={f}  n_probes={n}")

    # save the gene-level expression matrix (panel only)
    gene_out = OUT / "expression_matrix_log_gene.tsv.gz"
    gene_expr_df[sample_cols].to_csv(gene_out, sep="\t", compression="gzip")
    print(f"[save] {gene_out} (panel genes only, log2 gcRMA, GPL570 best-probe collapsed)")

    # ---------- compute scores ----------
    print("[score] computing module scores")
    rai8_score, rai8_used = module_score(gene_expr_df, sample_cols, RAI_8)
    nonol_score, nonol_used = module_score(gene_expr_df, sample_cols, NONOVERLAP)
    tds_score, tds_used = module_score(gene_expr_df, sample_cols, TDS_LIKE)
    tfcoll_score, tfcoll_used = module_score(gene_expr_df, sample_cols, TF_COLLAPSE)
    sad_score, sad_used = module_score(gene_expr_df, sample_cols, STAT3_AP1_DNMT)
    # TACSTD2 single-gene z
    if gene_expr_df.loc[TACSTD2_GENE, "found"]:
        tacstd2_z = pd.Series(zscore_within_cohort(gene_expr_df.loc[TACSTD2_GENE, sample_cols].values),
                                index=sample_cols)
    else:
        tacstd2_z = pd.Series(np.nan, index=sample_cols)

    scores_df = pd.DataFrame({
        "sample_id": sample_cols,
        "RAI_8_score": rai8_score.values,
        "DM1_like_score": -rai8_score.values,
        "THYROID_NONOVERLAP_score": nonol_score.values,
        "TDS_like_score": tds_score.values,
        "TF_collapse_score": tfcoll_score.values,         # FOXE1+NKX2-1+PAX8+HHEX (high → DIFFERENTIATED)
        "STAT3_AP1_DNMT_score": sad_score.values,         # STAT3+FOSL1+JUNB+DNMT1+DNMT3B
        "TACSTD2_z": tacstd2_z.values,
    })
    scores_df = scores_df.merge(meta_df[["sample_id", "histology", "disease_group"]], on="sample_id")
    scores_out = OUT / "scores.tsv"
    scores_df.to_csv(scores_out, sep="\t", index=False)
    print(f"[save] {scores_out}")

    # ---------- minimal tests ----------
    print("[tests] minimal")
    test_rows = []
    def mw_test(label, x, y, label_x, label_y):
        x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
        x = x[~np.isnan(x)]; y = y[~np.isnan(y)]
        if len(x) < 2 or len(y) < 2:
            test_rows.append({"test": label, "n_x": len(x), "n_y": len(y),
                              "mean_x": np.nanmean(x), "mean_y": np.nanmean(y),
                              "cohens_d": np.nan, "MW_p_two_sided": np.nan, "note": "insufficient n"})
            return
        u, p = mannwhitneyu(x, y, alternative="two-sided")
        d = cohens_d(x, y)
        test_rows.append({"test": label, "n_x": len(x), "n_y": len(y),
                          "label_x": label_x, "label_y": label_y,
                          "mean_x": np.nanmean(x), "mean_y": np.nanmean(y),
                          "cohens_d": d, "MW_U": u, "MW_p_two_sided": p, "note": ""})
    g_atc = scores_df[scores_df["histology"] == "ATC"]
    g_pdtc = scores_df[scores_df["histology"] == "PDTC"]
    # 1) ATC vs PDTC for each axis (the only within-cohort contrast available — no PTC controls)
    for axis in ["RAI_8_score", "DM1_like_score", "THYROID_NONOVERLAP_score", "TDS_like_score",
                 "TF_collapse_score", "STAT3_AP1_DNMT_score", "TACSTD2_z"]:
        mw_test(f"{axis} :: ATC vs PDTC",
                g_atc[axis].values, g_pdtc[axis].values,
                "ATC", "PDTC")
    tests_df = pd.DataFrame(test_rows)
    # 2) Spearman DM1_like vs THYROID_NONOVERLAP across all samples
    rho_all, p_all = spearmanr(scores_df["DM1_like_score"], scores_df["THYROID_NONOVERLAP_score"], nan_policy="omit")
    print(f"  Spearman(DM1_like, NONOVERLAP) all37: rho={rho_all:.3f}, p={p_all:.2e}")
    # 3) Spearman TF_collapse_score vs DM1_like
    rho_tf, p_tf = spearmanr(scores_df["TF_collapse_score"], scores_df["DM1_like_score"], nan_policy="omit")
    print(f"  Spearman(TF_collapse, DM1_like) all37: rho={rho_tf:.3f}, p={p_tf:.2e}")
    # 4) Spearman STAT3_AP1_DNMT vs DM1_like
    rho_s, p_s = spearmanr(scores_df["STAT3_AP1_DNMT_score"], scores_df["DM1_like_score"], nan_policy="omit")
    print(f"  Spearman(STAT3_AP1_DNMT, DM1_like) all37: rho={rho_s:.3f}, p={p_s:.2e}")
    # 5) Spearman TACSTD2_z vs DM1_like
    rho_t, p_t = spearmanr(scores_df["TACSTD2_z"], scores_df["DM1_like_score"], nan_policy="omit")
    print(f"  Spearman(TACSTD2, DM1_like) all37: rho={rho_t:.3f}, p={p_t:.2e}")
    extra = pd.DataFrame([
        {"test": "Spearman :: DM1_like vs THYROID_NONOVERLAP all37", "rho": rho_all, "p": p_all, "n": len(scores_df)},
        {"test": "Spearman :: TF_collapse vs DM1_like all37", "rho": rho_tf, "p": p_tf, "n": len(scores_df)},
        {"test": "Spearman :: STAT3_AP1_DNMT vs DM1_like all37", "rho": rho_s, "p": p_s, "n": len(scores_df)},
        {"test": "Spearman :: TACSTD2 vs DM1_like all37", "rho": rho_t, "p": p_t, "n": len(scores_df)},
    ])
    tests_out = OUT / "score_tests.tsv"
    pd.concat([tests_df, extra], axis=0, ignore_index=True).to_csv(tests_out, sep="\t", index=False)
    print(f"[save] {tests_out}")

    # ---------- figures ----------
    print("[fig] boxplot")
    plt.figure(figsize=(8, 4.5))
    panels = [
        ("RAI_8_score", "RAI_8 (8-gene)"),
        ("DM1_like_score", "DM1_like (= -RAI_8)"),
        ("THYROID_NONOVERLAP_score", "THYROID_NONOVERLAP (zero-overlap)"),
        ("TDS_like_score", "TDS-like (RAI_8 ∪ NONOVERLAP)"),
    ]
    fig, axes = plt.subplots(1, len(panels), figsize=(3.0*len(panels), 3.6), sharey=False)
    for ax, (col, title) in zip(axes, panels):
        data = [scores_df.loc[scores_df["histology"] == h, col].values for h in ["ATC", "PDTC"]]
        ax.boxplot(data, labels=["ATC", "PDTC"], showmeans=True, widths=0.5)
        ax.set_title(title, fontsize=9)
        ax.set_ylabel("within-cohort z" if "score" in col else "")
        # overlay points
        for i, d in enumerate(data, start=1):
            jitter = np.random.RandomState(RNG_SEED + i).uniform(-0.08, 0.08, size=len(d))
            ax.scatter(np.full_like(d, i, dtype=float) + jitter, d, s=14, alpha=0.7, color="black")
    fig.suptitle("GSE76039 (n=37; 20 ATC, 17 PDTC) — RAI / lineage axis distributions", fontsize=10)
    fig.tight_layout()
    fig_box = OUT / "gse76039_dm1_lineage_boxplot.png"
    fig.savefig(fig_box, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {fig_box}")

    # scatter DM1_like vs NONOVERLAP
    print("[fig] scatter")
    fig, ax = plt.subplots(figsize=(5, 5))
    colors = {"ATC": "#d62728", "PDTC": "#1f77b4"}
    for h, sub in scores_df.groupby("histology"):
        ax.scatter(sub["DM1_like_score"], sub["THYROID_NONOVERLAP_score"],
                    c=colors.get(h, "gray"), s=42, alpha=0.85, label=f"{h} (n={len(sub)})", edgecolor="white")
    ax.set_xlabel("DM1_like score (within-cohort z)")
    ax.set_ylabel("THYROID_NONOVERLAP score (within-cohort z)")
    ax.set_title(f"GSE76039 — DM1_like vs zero-overlap lineage score\nSpearman ρ = {rho_all:.3f} (p = {p_all:.2e}, n=37)", fontsize=10)
    ax.axhline(0, ls=":", c="gray", lw=0.7); ax.axvline(0, ls=":", c="gray", lw=0.7)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig_sc = OUT / "gse76039_dm1_vs_nonoverlap_scatter.png"
    fig.savefig(fig_sc, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {fig_sc}")

    # mechanism heatmap
    print("[fig] mechanism heatmap")
    mech_genes = ["FOXE1", "NKX2-1", "PAX8", "HHEX",
                   "DNMT1", "DNMT3B", "STAT3", "FOSL1", "JUNB",
                   "TACSTD2"]
    mech_avail = [g for g in mech_genes if gene_expr_df.loc[g, "found"]]
    if not mech_avail:
        print("  [skip] no mechanism genes available")
    else:
        # z within cohort per gene
        sub = gene_expr_df.loc[mech_avail, sample_cols].astype(float)
        z = sub.apply(zscore_within_cohort, axis=1, result_type="broadcast")
        # order samples by histology then DM1_like
        order_meta = scores_df.set_index("sample_id").loc[sample_cols][["histology", "DM1_like_score"]].copy()
        order_meta["hist_rank"] = order_meta["histology"].map({"ATC": 0, "PDTC": 1, "PTC": 2, "Normal": 3, "UNKNOWN": 9})
        sample_order = order_meta.sort_values(["hist_rank", "DM1_like_score"]).index.tolist()
        z = z[sample_order]
        fig, ax = plt.subplots(figsize=(11, 0.32*len(mech_avail) + 1.8))
        vmax = float(np.nanpercentile(np.abs(z.values), 95)); vmax = max(vmax, 1.5)
        im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_yticks(range(len(mech_avail)))
        ax.set_yticklabels(mech_avail)
        ax.set_xticks(range(len(sample_order)))
        # color sample labels by histology
        ax.set_xticklabels(sample_order, rotation=90, fontsize=6)
        for tl, sid in zip(ax.get_xticklabels(), sample_order):
            h = order_meta.loc[sid, "histology"]
            tl.set_color(colors.get(h, "black"))
        ax.set_title(f"GSE76039 mechanism panel — within-cohort z (samples ordered by histology then DM1_like; ATC red, PDTC blue)", fontsize=9)
        cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
        cbar.set_label("z (per gene, within cohort)")
        fig.tight_layout()
        fig_h = OUT / "gse76039_mechanism_heatmap.png"
        fig.savefig(fig_h, dpi=160, bbox_inches="tight")
        plt.close(fig)
        print(f"[save] {fig_h}")

    # ---------- summary printout ----------
    print("\n=== histology-level score means (within-cohort z) ===")
    for col in ["RAI_8_score", "DM1_like_score", "THYROID_NONOVERLAP_score", "TDS_like_score",
                 "TF_collapse_score", "STAT3_AP1_DNMT_score", "TACSTD2_z"]:
        m = scores_df.groupby("histology")[col].agg(["count", "mean", "std"]).round(3)
        print(f"\n--- {col} ---")
        print(m.to_string())

    print("\nDONE.")


if __name__ == "__main__":
    np.random.seed(RNG_SEED)
    main()
