"""
Cross-cohort gene-panel combination scan — TCGA-THCA × Korean cohorts.

Cohorts:
  - TCGA-THCA n=179 (DM1=152, DM2=27)              [DM1 vs DM2]
  - Korean K2 PRJEB11591 n=260 (DM1=14, DM2=246)   [DM1 vs DM2, RAI_8-subset panels only;
                                                    K2 kallisto index = RAI_8 + TTF1/TTF2 only]
  - Korean GSE286332 n=18 (PTC=9, PTC+HT=9)        [PTC vs PTC+HT, full-transcriptome FPKM,
                                                    secondary cross-cohort sanity]

For each panel we report:
  TCGA_AUC, TCGA_ARI, K2_AUC, K2_ARI (if all panel genes in K2 quant),
  GSE286332_AUC (PTC vs PTC+HT), AUC_min_DM (TCGA & K2),
  direction_consistent_DM_3cohort (sign of AUC-0.5 agrees across all evaluable cohorts).

Outputs:
  - panel_combos_cross_cohort.tsv
  - fig_panel_combos_cross_cohort.png
"""
from __future__ import annotations
import gzip, sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import roc_auc_score, adjusted_rand_score
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

OUT = Path(__file__).resolve().parent
ROOT = Path("/home/seungho/personal/THCA_data_analysis")

RAI_8        = ["SLC5A5","TPO","TG","TSHR","PAX8","NKX2-1","FOXE1","DIO1"]
NONOVERLAP_8 = ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"]
TF_4         = ["FOXE1","NKX2-1","PAX8","HHEX"]
RAI_5_eff    = ["SLC5A5","TPO","TG","TSHR","DIO1"]
TDS_16       = sorted(set(RAI_8 + NONOVERLAP_8))
TIERA_extra  = ["DIO2","SLC5A8","THRA","THRB"]
MECHANISM    = ["STAT3","FOSL1","JUNB","DNMT1","DNMT3B"]

PANELS = {
    "RAI_8 (canonical)":               RAI_8,
    "TDS-16 (full)":                   TDS_16,
    "NONOVERLAP_8":                    NONOVERLAP_8,
    "TF_4_backbone":                   TF_4,
    "RAI_5_effectors_only":            RAI_5_eff,
    "TF_3_within_RAI8":                ["FOXE1","NKX2-1","PAX8"],
    "RAI_8 + HHEX (9)":                RAI_8 + ["HHEX"],
    "RAI_8 + TF_4 (12 unique)":        sorted(set(RAI_8 + TF_4)),
    "RAI_8 + NONOVERLAP_8 (16)":       RAI_8 + NONOVERLAP_8,
    "TF_4 + RAI_5 (9)":                TF_4 + RAI_5_eff,
    "TF_3 + iodide (NIS+TPO+TG+SLC26A4+IYD)":
                                       ["FOXE1","NKX2-1","PAX8","SLC5A5","TPO","TG","SLC26A4","IYD"],
    "Iodide-handling 6":               ["SLC5A5","SLC26A4","TPO","DUOX1","DUOX2","IYD"],
    "Lineage TF + Effector minimal 6": ["FOXE1","NKX2-1","PAX8","TG","TPO","DIO1"],
    "Lineage TF + Effector minimal 4": ["FOXE1","NKX2-1","TG","TPO"],
    "TIERA TDS-extended 12":           sorted(set(RAI_8 + TIERA_extra)),
    "Mechanism arm 5 (STAT3+AP1+DNMT)":MECHANISM,
}
ALL_GENES = sorted({g for gl in PANELS.values() for g in gl})

# ---------- TCGA ----------
def load_tcga():
    print(f"[TCGA] streaming pancan for {len(ALL_GENES)} genes...", flush=True)
    rows = {}
    needed = set(ALL_GENES) | {"TITF1","NKX2_1"}
    with gzip.open(ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz","rt") as f:
        header = next(f).strip().split("\t")
        samples = header[1:]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            sym = parts[0].strip()
            if sym in needed:
                rows[sym] = [float(x) if x not in ("","NA","NaN") else np.nan for x in parts[1:]]
    expr = pd.DataFrame(rows, index=samples).T
    if "NKX2-1" not in expr.index:
        if "NKX2_1" in expr.index: expr = expr.rename(index={"NKX2_1":"NKX2-1"})
        elif "TITF1" in expr.index: expr = expr.rename(index={"TITF1":"NKX2-1"})
    expr = expr[~expr.index.duplicated(keep="first")]
    master = pd.read_csv(ROOT/"project/results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv", sep="\t")
    master["short"] = master["sample_id"].str[:15]
    expr.columns = [c[:15] for c in expr.columns]
    common = sorted(set(expr.columns) & set(master["short"]))
    expr_thca = expr[common]
    dm_lookup = master.set_index("short")["dm"]
    dm = dm_lookup.loc[common]
    keep = dm.isin(["DM1","DM2"])
    print(f"[TCGA] expr {expr_thca.shape}; DM1={int((dm=='DM1').sum())} DM2={int((dm=='DM2').sum())}", flush=True)
    return expr_thca, dm, keep

# ---------- K2 (8 RAI genes only — kallisto index restriction) ----------
def load_k2():
    expr = pd.read_csv(ROOT/"project/results/v17_korean/K2_8gene_tpm_matrix_v4.tsv", sep="\t", index_col=0)
    # transpose → genes × samples
    expr = expr.T
    preds = pd.read_csv(ROOT/"project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t")
    preds = preds.set_index("run")
    # dedup tech replicates by sample_alias if metadata exists
    meta_path = ROOT/"project/results/v17_korean/K2_subset_v2.tsv"
    alias = None
    if meta_path.exists():
        meta = pd.read_csv(meta_path, sep="\t")
        if "run_accession" in meta.columns and "sample_alias" in meta.columns:
            alias = meta.set_index("run_accession")["sample_alias"]
    common = sorted(set(expr.columns) & set(preds.index))
    expr_k2 = expr[common]
    dm = preds.loc[common, "DM_call"]
    # dedup
    if alias is not None:
        groups = alias.reindex(common).fillna(pd.Series(common, index=common))
        # collapse: pick first run per alias
        seen = {}
        keep_runs = []
        for r in common:
            a = groups.loc[r]
            if a not in seen:
                seen[a] = r
                keep_runs.append(r)
        expr_k2 = expr_k2[keep_runs]
        dm = dm.loc[keep_runs]
        print(f"[K2] deduplicated by sample_alias: {len(common)} runs -> {len(keep_runs)} unique samples", flush=True)
    keep = dm.isin(["DM1","DM2"])
    print(f"[K2] expr {expr_k2.shape}; DM1={int((dm=='DM1').sum())} DM2={int((dm=='DM2').sum())}", flush=True)
    return expr_k2, dm, keep

# ---------- GSE286332 (full transcriptome, PTC vs PTC+HT) ----------
def load_gse286332():
    p = Path("/data/thca/v17_korean/GSE286332/GSE286332_all_sample_rawdata.txt.gz")
    print(f"[GSE286332] streaming {p}", flush=True)
    df = pd.read_csv(p, sep="\t")
    fpkm_cols = [c for c in df.columns if c.endswith("_FPKM")]
    if not fpkm_cols:
        raise RuntimeError("FPKM columns missing in GSE286332")
    sym = df["Gene_Symbol"].astype(str)
    expr = df[fpkm_cols].copy()
    expr.index = sym
    expr.columns = [c.replace("_FPKM","") for c in fpkm_cols]
    expr = expr.groupby(level=0).sum()
    sub = expr.loc[[g for g in ALL_GENES if g in expr.index]]
    preds = pd.read_csv(ROOT/"project/results/p3_gse286332/dm12_predictions.tsv", sep="\t", index_col=0)
    common = sorted(set(sub.columns) & set(preds.index))
    sub = sub[common]
    grp = preds.loc[common, "group"]   # PTC vs PTC_HT
    keep = grp.isin(["PTC","PTC_HT"])
    print(f"[GSE286332] expr {sub.shape}; PTC={int((grp=='PTC').sum())} PTC_HT={int((grp=='PTC_HT').sum())}", flush=True)
    return sub, grp, keep

# ---------- scoring ----------
def score_panel(expr, label, keep, gene_list, pos_label, log_tpm=False, require_full=False):
    found = [g for g in gene_list if g in expr.index]
    n_panel = len(gene_list); n_found = len(found)
    if n_found < 2:
        return dict(AUC=np.nan, ARI=np.nan, n_found=n_found, n_panel=n_panel,
                    n_samples=int(keep.sum()), missing=",".join(sorted(set(gene_list)-set(found))))
    if require_full and n_found < n_panel:
        return dict(AUC=np.nan, ARI=np.nan, n_found=n_found, n_panel=n_panel,
                    n_samples=int(keep.sum()), missing=",".join(sorted(set(gene_list)-set(found))))
    idx_keep = keep[keep].index
    sub = expr.loc[found, idx_keep].copy()
    if log_tpm:
        sub = np.log1p(sub.clip(lower=0))
    sub = sub.replace([np.inf,-np.inf], np.nan).ffill(axis=1).bfill(axis=1).fillna(0)
    sub_z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    score = sub_z.mean(axis=0)
    ref = (label.loc[idx_keep] == pos_label).astype(int).values
    if ref.sum() == 0 or ref.sum() == len(ref) or score.values.std() == 0:
        return dict(AUC=np.nan, ARI=np.nan, n_found=n_found, n_panel=n_panel,
                    n_samples=int(keep.sum()), missing=",".join(sorted(set(gene_list)-set(found))))
    auc = roc_auc_score(ref, -score.values)
    if len(np.unique(sub_z.T.values, axis=0)) >= 2:
        km = KMeans(n_clusters=2, n_init=10, random_state=0).fit(sub_z.T.values)
        labels_km = km.labels_
        m0 = score.values[labels_km==0].mean(); m1 = score.values[labels_km==1].mean()
        if m0 < m1: cl = (labels_km==0).astype(int)
        else:       cl = (labels_km==1).astype(int)
        ari = adjusted_rand_score(ref, cl)
    else:
        ari = np.nan
    return dict(AUC=round(float(auc),3), ARI=round(float(ari),3) if np.isfinite(ari) else np.nan,
                n_found=n_found, n_panel=n_panel, n_samples=int(keep.sum()),
                missing=",".join(sorted(set(gene_list)-set(found))))

# ---------- main ----------
def main():
    expr_tcga, dm_tcga, keep_tcga = load_tcga()
    expr_k2,   dm_k2,   keep_k2   = load_k2()
    expr_gs,   grp_gs,  keep_gs   = load_gse286332()

    rows = []
    for name, gene_list in PANELS.items():
        t = score_panel(expr_tcga, dm_tcga, keep_tcga, gene_list, pos_label="DM1", log_tpm=False)
        k = score_panel(expr_k2,   dm_k2,   keep_k2,   gene_list, pos_label="DM1", log_tpm=True,
                        require_full=True)
        g = score_panel(expr_gs,   grp_gs,  keep_gs,   gene_list, pos_label="PTC_HT", log_tpm=True)
        # robust metrics — only over cohorts with finite AUC AND full panel coverage
        dm_aucs = [a for a in [t["AUC"], k["AUC"]] if np.isfinite(a)]
        auc_min_dm = round(float(np.min(dm_aucs)),3) if dm_aucs else np.nan
        auc_geomean_dm = round(float(np.sqrt(np.prod(dm_aucs))),3) if (len(dm_aucs)==2 and all(a>0 for a in dm_aucs)) else np.nan
        # direction-consistent: TCGA & K2 both same side of 0.5 (DM1 task), independent of GSE286332
        dir_ok = bool(np.isfinite(t["AUC"]) and np.isfinite(k["AUC"]) and ((t["AUC"]-0.5)*(k["AUC"]-0.5) > 0))
        rows.append({
            "panel": name, "n_panel": len(gene_list),
            "TCGA_AUC": t["AUC"], "TCGA_ARI": t["ARI"], "TCGA_n": t["n_samples"],
            "K2_AUC":   k["AUC"], "K2_ARI":   k["ARI"], "K2_n":   k["n_samples"], "K2_n_found": k["n_found"],
            "GSE286332_AUC": g["AUC"], "GSE286332_ARI": g["ARI"], "GSE286332_n": g["n_samples"], "GSE286332_n_found": g["n_found"],
            "AUC_min_DM_TCGA_K2": auc_min_dm,
            "AUC_geomean_DM_TCGA_K2": auc_geomean_dm,
            "direction_consistent_TCGA_K2": dir_ok,
            "K2_panel_coverage": "FULL" if k["n_found"]==len(gene_list) else f"PARTIAL ({k['n_found']}/{len(gene_list)})",
            "K2_missing_genes": k["missing"],
        })
    df = pd.DataFrame(rows).sort_values(["AUC_min_DM_TCGA_K2","TCGA_AUC"], ascending=[False, False],
                                        na_position="last")
    df.to_csv(OUT/"panel_combos_cross_cohort.tsv", sep="\t", index=False)
    print("\n[joint result]"); print(df.to_string(index=False))

    # ---------- plot ----------
    fig = plt.figure(figsize=(17, 11))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1], width_ratios=[1.05, 1])

    # (a) scatter TCGA vs K2 (only panels with FULL K2 coverage)
    ax0 = fig.add_subplot(gs[0, 0])
    full = df[df["K2_panel_coverage"]=="FULL"].dropna(subset=["TCGA_AUC","K2_AUC"])
    ax0.scatter(full["TCGA_AUC"], full["K2_AUC"], s=120, c="#244e73", edgecolor="black", zorder=3)
    for _, r in full.iterrows():
        ax0.annotate(r["panel"], (r["TCGA_AUC"], r["K2_AUC"]),
                     xytext=(6,4), textcoords="offset points", fontsize=8, color="#102033")
    ax0.axhline(0.5, color="grey", linestyle=":"); ax0.axvline(0.5, color="grey", linestyle=":")
    ax0.axhline(0.85, color="#426b50", linestyle="--", lw=0.8, label="AUC=0.85")
    ax0.axvline(0.85, color="#426b50", linestyle="--", lw=0.8)
    ax0.plot([0,1],[0,1], color="#b58534", linestyle="-.", lw=0.6, label="y=x")
    ax0.fill_between([0.85,1],0.85,1, color="#e9f3e8", alpha=0.5, zorder=0, label="cross-cohort sweet spot")
    ax0.set_xlim(0,1); ax0.set_ylim(0,1)
    ax0.set_xlabel("TCGA-THCA AUC (DM1 vs DM2)"); ax0.set_ylabel("Korean K2 AUC (DM1 vs DM2)")
    ax0.set_title("(a) TCGA × Korean K2 — only RAI_8-subset panels evaluable", fontsize=11)
    ax0.legend(fontsize=8, loc="lower right"); ax0.grid(alpha=0.25)

    # (b) horizontal bar — AUC_min ranking (cross-cohort robustness)
    ax1 = fig.add_subplot(gs[0, 1])
    full_sorted = full.sort_values("AUC_min_DM_TCGA_K2", ascending=True)
    y = np.arange(len(full_sorted))
    ax1.barh(y, full_sorted["AUC_min_DM_TCGA_K2"], color="#2c5e9c", edgecolor="black",
             label="AUC_min (worst-of TCGA, K2)")
    ax1.barh(y, full_sorted["AUC_geomean_DM_TCGA_K2"], color="#9bb6d4", edgecolor="none",
             alpha=0.45, label="AUC_geomean")
    ax1.axvline(0.85, color="#426b50", linestyle="--", lw=0.8)
    ax1.set_yticks(y); ax1.set_yticklabels(full_sorted["panel"], fontsize=8)
    ax1.set_xlim(0,1); ax1.set_xlabel("AUC")
    ax1.set_title("(b) Cross-cohort robustness — TCGA × K2 worst-of vs geomean", fontsize=11)
    ax1.legend(fontsize=8, loc="lower right"); ax1.grid(axis="x", alpha=0.25)

    # (c) 3-cohort side-by-side bars
    ax2 = fig.add_subplot(gs[1, :])
    bar_df = df.sort_values("TCGA_AUC", ascending=False)
    x = np.arange(len(bar_df)); w = 0.27
    ax2.bar(x-w, bar_df["TCGA_AUC"], width=w, color="#2c5e9c", edgecolor="black", label="TCGA (DM1 vs DM2)")
    # K2 — only plot where FULL coverage
    k2_vals = [v if cov=="FULL" else np.nan for v, cov in zip(bar_df["K2_AUC"], bar_df["K2_panel_coverage"])]
    ax2.bar(x, k2_vals, width=w, color="#c0392b", edgecolor="black", label="K2 (DM1 vs DM2; only RAI_8 subset)")
    ax2.bar(x+w, bar_df["GSE286332_AUC"], width=w, color="#b58534", edgecolor="black",
            label="GSE286332 (PTC vs PTC+HT; full transcriptome)")
    ax2.set_xticks(x); ax2.set_xticklabels(bar_df["panel"], rotation=35, ha="right", fontsize=8)
    ax2.axhline(0.5, color="grey", linestyle=":"); ax2.axhline(0.85, color="#426b50", linestyle="--", lw=0.8)
    ax2.set_ylim(0,1); ax2.set_ylabel("AUC")
    ax2.set_title("(c) Per-panel AUC across 3 cohorts (panels sorted by TCGA AUC)", fontsize=11)
    ax2.legend(fontsize=9, loc="upper right")
    ax2.grid(axis="y", alpha=0.25)

    plt.suptitle("Cross-cohort gene-panel scan — TCGA-THCA × Korean K2 (PRJEB11591) × Korean GSE286332\n"
                 "Reviewer-targeted: which 8-gene-style panels work on BOTH TCGA AND a Korean cohort?",
                 fontsize=12.5, y=1.00)
    plt.tight_layout()
    plt.savefig(OUT/"fig_panel_combos_cross_cohort.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[plot] {OUT}/fig_panel_combos_cross_cohort.png")

if __name__ == "__main__":
    main()
