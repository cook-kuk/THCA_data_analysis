"""
Cross-cohort panel scan v2 — adds K2 DM1 subset full-transcriptome (n=23).

Cohorts (v2):
  - TCGA-THCA n=179         (DM1 vs DM2, master r9_1)
  - K2 n=260, 8 RAI genes   (DM1 vs DM2; original K2 quant — RAI_8 panel only)
  - K2 n=23 full transcriptome subset
       = 14 DM1 (newly re-quanted from ENA FASTQ) + 9 DM2 (already done)
       (DM1 vs DM2; ALL 16-gene panels evaluable)
  - Korean GSE286332 n=18   (PTC vs PTC+HT; full transcriptome FPKM)

Outputs:
  - panel_combos_cross_cohort.tsv   (replaces v1)
  - fig_panel_combos_cross_cohort.png
"""
from __future__ import annotations
import gzip, sys, os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import roc_auc_score, adjusted_rand_score
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

OUT  = Path(__file__).resolve().parent
ROOT = Path("/home/seungho/personal/THCA_data_analysis")
K2_FULL_DIR = Path("/data/thca/v17_korean_k2_full_quant")

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
    print(f"[TCGA] streaming pancan", flush=True)
    rows = {}
    needed = set(ALL_GENES) | {"TITF1","NKX2_1"}
    with gzip.open(ROOT/"project/data/raw/TCGA_pancan/pancan_geneExp.gz","rt") as f:
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
    dm = master.set_index("short")["dm"].loc[common]
    keep = dm.isin(["DM1","DM2"])
    print(f"[TCGA] DM1={(dm=='DM1').sum()} DM2={(dm=='DM2').sum()}", flush=True)
    return expr_thca, dm, keep

# ---------- K2 8-gene 260-sample (existing TPM matrix) ----------
def load_k2_rai8():
    expr = pd.read_csv(ROOT/"project/results/v17_korean/K2_8gene_tpm_matrix_v4.tsv", sep="\t", index_col=0).T
    preds = pd.read_csv(ROOT/"project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t").set_index("run")
    common = sorted(set(expr.columns) & set(preds.index))
    expr_k2 = expr[common]
    dm = preds.loc[common, "DM_call"]
    keep = dm.isin(["DM1","DM2"])
    print(f"[K2_RAI8] DM1={(dm=='DM1').sum()} DM2={(dm=='DM2').sum()} (RAI_8 only)", flush=True)
    return expr_k2, dm, keep

# ---------- K2 full transcriptome (23-sample subset) ----------
def _agg_one_sample(args):
    run, path, gene_set = args
    out = {g: 0.0 for g in gene_set}
    try:
        with open(path) as f:
            next(f)
            for line in f:
                parts = line.split("\t")
                fields = parts[0].split("|")
                if len(fields) >= 6:
                    sym = fields[5]
                    if sym in gene_set:
                        try: out[sym] += float(parts[4])
                        except ValueError: pass
    except FileNotFoundError:
        return run, None
    return run, out

def load_k2_full():
    runs = sorted(d.name for d in K2_FULL_DIR.iterdir()
                  if (d/"abundance.tsv").exists() and (d/"abundance.tsv").stat().st_size > 0)
    print(f"[K2_full] aggregating {len(runs)} samples × {len(ALL_GENES)} genes (parallel)...", flush=True)
    cache = OUT/"k2_full_gene_tpm_matrix.tsv"
    if cache.exists():
        df_cache = pd.read_csv(cache, sep="\t", index_col=0)
        if set(runs).issubset(set(df_cache.columns)) and set(ALL_GENES).issubset(set(df_cache.index)):
            print(f"[K2_full] cached", flush=True)
            sub = df_cache.loc[ALL_GENES, runs]
            preds = pd.read_csv(ROOT/"project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t").set_index("run")
            dm = preds.loc[runs, "DM_call"]
            keep = dm.isin(["DM1","DM2"])
            print(f"[K2_full] DM1={(dm=='DM1').sum()} DM2={(dm=='DM2').sum()}", flush=True)
            return sub, dm, keep
    tasks = [(r, K2_FULL_DIR/r/"abundance.tsv", set(ALL_GENES)) for r in runs]
    res = {}
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
        for run, vals in ex.map(_agg_one_sample, tasks, chunksize=2):
            if vals: res[run] = vals
    expr = pd.DataFrame(res).reindex(ALL_GENES)
    expr.to_csv(cache, sep="\t")
    preds = pd.read_csv(ROOT/"project/results/v17_korean/K2_korean_predictions_v4.tsv", sep="\t").set_index("run")
    common = [r for r in runs if r in preds.index]
    expr = expr[common]
    dm = preds.loc[common, "DM_call"]
    keep = dm.isin(["DM1","DM2"])
    print(f"[K2_full] DM1={(dm=='DM1').sum()} DM2={(dm=='DM2').sum()}", flush=True)
    return expr, dm, keep

# ---------- GSE286332 ----------
def load_gse286332():
    p = Path("/data/thca/v17_korean/GSE286332/GSE286332_all_sample_rawdata.txt.gz")
    df = pd.read_csv(p, sep="\t")
    fpkm_cols = [c for c in df.columns if c.endswith("_FPKM")]
    sym = df["Gene_Symbol"].astype(str)
    expr = df[fpkm_cols].copy()
    expr.index = sym
    expr.columns = [c.replace("_FPKM","") for c in fpkm_cols]
    expr = expr.groupby(level=0).sum()
    sub = expr.loc[[g for g in ALL_GENES if g in expr.index]]
    preds = pd.read_csv(ROOT/"project/results/p3_gse286332/dm12_predictions.tsv", sep="\t", index_col=0)
    common = sorted(set(sub.columns) & set(preds.index))
    sub = sub[common]
    grp = preds.loc[common, "group"]
    keep = grp.isin(["PTC","PTC_HT"])
    print(f"[GSE286332] PTC={(grp=='PTC').sum()} PTC_HT={(grp=='PTC_HT').sum()}", flush=True)
    return sub, grp, keep

# ---------- scoring ----------
def score(expr, label, keep, gene_list, pos_label, log_tpm=False, require_full=False):
    found = [g for g in gene_list if g in expr.index]
    if len(found) < 2 or (require_full and len(found) < len(gene_list)):
        return dict(AUC=np.nan, ARI=np.nan, n_found=len(found))
    idx_keep = keep[keep].index
    sub = expr.loc[found, idx_keep].copy()
    if log_tpm: sub = np.log1p(sub.clip(lower=0))
    sub = sub.replace([np.inf,-np.inf], np.nan).ffill(axis=1).bfill(axis=1).fillna(0)
    sub_z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    score = sub_z.mean(axis=0)
    ref = (label.loc[idx_keep] == pos_label).astype(int).values
    if ref.sum() == 0 or ref.sum() == len(ref) or score.values.std() == 0:
        return dict(AUC=np.nan, ARI=np.nan, n_found=len(found))
    auc = roc_auc_score(ref, -score.values)
    try:
        if len(np.unique(sub_z.T.values, axis=0)) >= 2:
            km = KMeans(n_clusters=2, n_init=10, random_state=0).fit(sub_z.T.values)
            labels_km = km.labels_
            m0 = score.values[labels_km==0].mean(); m1 = score.values[labels_km==1].mean()
            cl = (labels_km==0).astype(int) if m0 < m1 else (labels_km==1).astype(int)
            ari = adjusted_rand_score(ref, cl)
        else:
            ari = np.nan
    except Exception:
        ari = np.nan
    return dict(AUC=round(float(auc),3), ARI=round(float(ari),3) if np.isfinite(ari) else np.nan,
                n_found=len(found))

# ---------- main ----------
def main():
    expr_t, dm_t, keep_t = load_tcga()
    expr_r, dm_r, keep_r = load_k2_rai8()
    expr_f, dm_f, keep_f = load_k2_full()
    expr_g, grp_g, keep_g = load_gse286332()

    rows = []
    for name, gene_list in PANELS.items():
        t = score(expr_t, dm_t, keep_t, gene_list, pos_label="DM1", log_tpm=False)
        r = score(expr_r, dm_r, keep_r, gene_list, pos_label="DM1", log_tpm=True, require_full=True)
        f = score(expr_f, dm_f, keep_f, gene_list, pos_label="DM1", log_tpm=True)
        g = score(expr_g, grp_g, keep_g, gene_list, pos_label="PTC_HT", log_tpm=True)
        # robust metrics: AUC_min across cohorts with finite, FULL-coverage AUC
        coverage_full_for_k2_rai8 = (r["n_found"] == len(gene_list))
        coverage_full_for_k2_full = (f["n_found"] == len(gene_list))
        dm_aucs = []
        if np.isfinite(t["AUC"]): dm_aucs.append(t["AUC"])
        if np.isfinite(f["AUC"]) and coverage_full_for_k2_full: dm_aucs.append(f["AUC"])
        if np.isfinite(r["AUC"]) and coverage_full_for_k2_rai8: dm_aucs.append(r["AUC"])
        auc_min = round(float(np.min(dm_aucs)),3) if dm_aucs else np.nan
        auc_geomean = round(float(np.prod(dm_aucs) ** (1.0/len(dm_aucs))),3) if (dm_aucs and all(a>0 for a in dm_aucs)) else np.nan
        # direction-consistency: all evaluable DM cohorts on same side of 0.5
        eval_auc = [t["AUC"]]
        if coverage_full_for_k2_full and np.isfinite(f["AUC"]): eval_auc.append(f["AUC"])
        if coverage_full_for_k2_rai8 and np.isfinite(r["AUC"]): eval_auc.append(r["AUC"])
        dir_ok = bool(len(eval_auc) >= 2 and (all(a > 0.5 for a in eval_auc) or all(a < 0.5 for a in eval_auc)))
        rows.append({
            "panel": name, "n_panel": len(gene_list),
            "TCGA_AUC": t["AUC"], "TCGA_ARI": t["ARI"],
            "K2_n260_RAI8_AUC": r["AUC"] if coverage_full_for_k2_rai8 else np.nan,
            "K2_n260_RAI8_coverage": "FULL" if coverage_full_for_k2_rai8 else f"PARTIAL({r['n_found']}/{len(gene_list)})",
            "K2_n23_full_AUC": f["AUC"] if coverage_full_for_k2_full else np.nan,
            "K2_n23_full_ARI": f["ARI"] if coverage_full_for_k2_full else np.nan,
            "K2_n23_full_coverage": "FULL" if coverage_full_for_k2_full else f"PARTIAL({f['n_found']}/{len(gene_list)})",
            "GSE286332_AUC": g["AUC"], "GSE286332_ARI": g["ARI"],
            "AUC_min_DM": auc_min,
            "AUC_geomean_DM": auc_geomean,
            "direction_consistent": dir_ok,
            "n_eval_DM_cohorts": len(eval_auc),
        })
    df = pd.DataFrame(rows).sort_values(["AUC_min_DM","TCGA_AUC"], ascending=[False, False], na_position="last")
    df.to_csv(OUT/"panel_combos_cross_cohort.tsv", sep="\t", index=False)
    print("\n[joint result v2]"); print(df.to_string(index=False))

    # ---------- plot ----------
    fig = plt.figure(figsize=(17, 11.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1], width_ratios=[1.05, 1])

    # (a) TCGA vs K2_full scatter (n=23 subset; now ALL panels evaluable)
    ax0 = fig.add_subplot(gs[0, 0])
    full = df.dropna(subset=["TCGA_AUC","K2_n23_full_AUC"])
    colors = ["#244e73" if r["direction_consistent"] else "#c0392b" for _, r in full.iterrows()]
    ax0.scatter(full["TCGA_AUC"], full["K2_n23_full_AUC"], s=130, c=colors, edgecolor="black", zorder=3)
    for _, r in full.iterrows():
        ax0.annotate(r["panel"], (r["TCGA_AUC"], r["K2_n23_full_AUC"]),
                     xytext=(6,4), textcoords="offset points", fontsize=7.5, color="#102033")
    ax0.axhline(0.5, color="grey", linestyle=":"); ax0.axvline(0.5, color="grey", linestyle=":")
    ax0.axhline(0.85, color="#426b50", linestyle="--", lw=0.8, label="AUC=0.85")
    ax0.axvline(0.85, color="#426b50", linestyle="--", lw=0.8)
    ax0.plot([0,1],[0,1], color="#b58534", linestyle="-.", lw=0.6, label="y=x")
    ax0.fill_between([0.85,1],0.85,1, color="#e9f3e8", alpha=0.5, zorder=0, label="cross-cohort sweet spot")
    ax0.set_xlim(0,1); ax0.set_ylim(0,1)
    ax0.set_xlabel("TCGA-THCA AUC (DM1 vs DM2, n=179)")
    ax0.set_ylabel("Korean K2 full-quant subset AUC (DM1 vs DM2, n=23)")
    ax0.set_title("(a) TCGA × K2 n=23 full-transcriptome — ALL 16 panels evaluable", fontsize=11)
    ax0.legend(fontsize=8, loc="lower right"); ax0.grid(alpha=0.25)

    # (b) AUC_min ranking with multi-cohort weighting
    ax1 = fig.add_subplot(gs[0, 1])
    ranked = df.dropna(subset=["AUC_min_DM"]).sort_values("AUC_min_DM", ascending=True)
    y = np.arange(len(ranked))
    ax1.barh(y, ranked["AUC_min_DM"], color="#2c5e9c", edgecolor="black",
             label="AUC_min across DM cohorts")
    ax1.barh(y, ranked["AUC_geomean_DM"], color="#9bb6d4", edgecolor="none",
             alpha=0.45, label="AUC_geomean")
    ax1.axvline(0.85, color="#426b50", linestyle="--", lw=0.8)
    ax1.set_yticks(y); ax1.set_yticklabels(ranked["panel"], fontsize=8)
    ax1.set_xlim(0,1); ax1.set_xlabel("AUC")
    ax1.set_title("(b) Cross-cohort robustness — AUC_min vs AUC_geomean", fontsize=11)
    ax1.legend(fontsize=8, loc="lower right"); ax1.grid(axis="x", alpha=0.25)

    # (c) 4-cohort side-by-side bars
    ax2 = fig.add_subplot(gs[1, :])
    bar_df = df.sort_values("TCGA_AUC", ascending=False)
    x = np.arange(len(bar_df)); w = 0.20
    ax2.bar(x-1.5*w, bar_df["TCGA_AUC"], width=w, color="#2c5e9c", edgecolor="black", label="TCGA n=179")
    k2r = [v if c=="FULL" else np.nan for v, c in zip(bar_df["K2_n260_RAI8_AUC"], bar_df["K2_n260_RAI8_coverage"])]
    ax2.bar(x-0.5*w, k2r, width=w, color="#c0392b", edgecolor="black", label="K2 n=260 (RAI_8 only)")
    k2f = [v if c=="FULL" else np.nan for v, c in zip(bar_df["K2_n23_full_AUC"], bar_df["K2_n23_full_coverage"])]
    ax2.bar(x+0.5*w, k2f, width=w, color="#f39c12", edgecolor="black", label="K2 n=23 full")
    ax2.bar(x+1.5*w, bar_df["GSE286332_AUC"], width=w, color="#b58534", edgecolor="black",
            label="GSE286332 n=18 (PTC vs PTC+HT)")
    ax2.set_xticks(x); ax2.set_xticklabels(bar_df["panel"], rotation=35, ha="right", fontsize=8)
    ax2.axhline(0.5, color="grey", linestyle=":"); ax2.axhline(0.85, color="#426b50", linestyle="--", lw=0.8)
    ax2.set_ylim(0,1); ax2.set_ylabel("AUC")
    ax2.set_title("(c) Per-panel AUC across 4 cohorts (sorted by TCGA AUC)", fontsize=11)
    ax2.legend(fontsize=9, loc="upper right", ncol=2)
    ax2.grid(axis="y", alpha=0.25)

    plt.suptitle("Cross-cohort gene-panel scan v2 — TCGA × K2 (full + RAI8-only) × GSE286332\n"
                 "K2 n=23 full-transcriptome subset (14 DM1 newly re-quanted + 9 DM2) enables real AUC for NONOVERLAP_8/TDS-16/TIERA",
                 fontsize=12.5, y=1.00)
    plt.tight_layout()
    plt.savefig(OUT/"fig_panel_combos_cross_cohort.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[plot] {OUT}/fig_panel_combos_cross_cohort.png")

if __name__ == "__main__":
    main()
