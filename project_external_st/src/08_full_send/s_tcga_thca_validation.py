#!/usr/bin/env python3
"""S — TCGA-THCA bulk RNA-seq validation of RAI_8 / DM1_like / THYROID_NONOVERLAP.

1. Load Pancan EBPP-adjusted geneExp (Entrez IDs × 11069 TCGA samples)
2. Filter to TCGA-THCA via clinical join with tcga_dm_master_with_pfi.tsv
3. Score 16 RAI/NONOVERLAP genes + Epithelial/Proliferation
4. Stage trend (Stage I→IVC ordinal)
5. Survival KM by DM1 quantile + Cox regression
6. Compare to GSE250521 stage trend
"""
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu

OUT = Path("project_external_st/results/extra")

# Entrez Gene IDs for our gene panels
RAI_8 = ["TPO","DIO1","TSHR","PAX8","TG","FOXE1","NKX2-1","SLC5A5"]
NONOVERLAP = ["SLC26A4","IYD","DUOX1","DUOX2","TFF3","HHEX","GLIS3","DIO2"]
EPITHELIAL = ["EPCAM","KRT8","KRT18","KRT19","TACSTD2"]
PROLIF = ["MKI67","TOP2A","PCNA","MCM2","MCM5","STMN1","CENPF"]
ALIASES = {"NKX2-1": ["NKX2-1","NKX2_1","TITF1"]}
ALL_GENES = set(RAI_8 + NONOVERLAP + EPITHELIAL + PROLIF) | {"NKX2_1","TITF1"}


def load_target_genes(gz_path: str, target_symbols: set[str]) -> pd.DataFrame:
    """Stream pancan TSV, keep only rows whose row label matches a target gene symbol."""
    rows = {}
    with gzip.open(gz_path, "rt") as f:
        header = next(f).strip().split("\t")
        samples = header[1:]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            sym = parts[0].strip()
            if sym in target_symbols:
                rows[sym] = [float(x) if x not in ("","NA","NaN") else np.nan for x in parts[1:]]
    df = pd.DataFrame(rows, index=samples).T  # genes × samples
    df.index.name = "gene_symbol"
    # consolidate aliases
    if "NKX2_1" in df.index and "NKX2-1" not in df.index:
        df = df.rename(index={"NKX2_1":"NKX2-1"})
    elif "TITF1" in df.index and "NKX2-1" not in df.index:
        df = df.rename(index={"TITF1":"NKX2-1"})
    return df


def main():
    print("Loading Pancan EBPP (streaming target genes only)...")
    expr = load_target_genes("project/data/raw/TCGA_pancan/pancan_geneExp.gz", ALL_GENES)
    print(f"  matrix: {expr.shape[0]} genes × {expr.shape[1]} TCGA pan-cancer samples")
    print(f"  found genes: {sorted(expr.index)}")

    # load clinical to identify THCA samples
    clin = pd.read_csv("project/results/dark_matter_phase2/web/data/tcga_dm_master_with_pfi.tsv", sep="\t")
    clin = clin.rename(columns={"tcga_short": "patient"})
    print(f"  THCA clinical: {len(clin)} patients")

    # match: pancan sample = TCGA-XX-XXXX-01 → patient = TCGA-XX-XXXX
    sample_to_patient = {s: "-".join(s.split("-")[:3]) for s in expr.columns}
    pancan_patients = pd.Series(sample_to_patient)
    keep_samples = pancan_patients[pancan_patients.isin(clin["patient"])].index.tolist()
    expr_thca = expr[keep_samples].copy()
    print(f"  THCA samples in pancan: {len(keep_samples)}")

    # within-cohort z-score per gene then mean → score
    z = expr_thca.sub(expr_thca.mean(axis=1), axis=0).div(expr_thca.std(axis=1), axis=0)

    def score(genes: list) -> pd.Series:
        avail = [g for g in genes if g in z.index]
        if not avail: return pd.Series(np.nan, index=z.columns)
        return z.loc[avail].mean(axis=0)

    scores = pd.DataFrame({
        "RAI_8":       score(RAI_8),
        "DM1_like":    -score(RAI_8),
        "NONOVERLAP":  score(NONOVERLAP),
        "Epithelial":  score(EPITHELIAL),
        "Proliferation": score(PROLIF),
    })
    scores.index.name = "sample"
    scores["patient"] = scores.index.map(sample_to_patient)
    merged = scores.merge(clin, on="patient", how="left")

    # ===== 1. Cross-validation: DM1 vs NONOVERLAP across n=507 =====
    rho, p = spearmanr(merged["DM1_like"], merged["NONOVERLAP"])
    from scipy.stats import pearsonr
    r, pp = pearsonr(merged["DM1_like"], merged["NONOVERLAP"])
    print(f"\n=== TCGA-THCA bulk DM1 vs NONOVERLAP (n={len(merged)}) ===")
    print(f"  Spearman ρ = {rho:.3f}, p = {p:.2e}")
    print(f"  Pearson r  = {r:.3f}, p = {pp:.2e}")

    # ===== 2. Stage trend =====
    stage_map = {"Stage I":1,"Stage II":2,"Stage III":3,
                 "Stage IVA":4,"Stage IVB":5,"Stage IVC":6}
    merged["stage_ord"] = merged["ajcc_pathologic_tumor_stage"].map(stage_map)
    stage_ok = merged.dropna(subset=["stage_ord", "DM1_like"])
    rho_s, p_s = spearmanr(stage_ok["stage_ord"], stage_ok["DM1_like"])
    print(f"\n=== TCGA-THCA stage trend (n={len(stage_ok)}) ===")
    print(f"  DM1_like vs stage_ord Spearman ρ = {rho_s:.3f}, p = {p_s:.2e}")
    rho_r, p_r = spearmanr(stage_ok["stage_ord"], stage_ok["RAI_8"])
    print(f"  RAI_8 vs stage_ord Spearman ρ    = {rho_r:.3f}, p = {p_r:.2e}")
    print("  per-stage mean DM1:")
    print(stage_ok.groupby("ajcc_pathologic_tumor_stage")["DM1_like"].agg(["mean","std","count"]).to_string())

    # ===== 3. DM1 vs DM cluster (Paper 1 framework) =====
    if "v17_dark_cluster" in merged.columns:
        dm_groups = merged.dropna(subset=["v17_dark_cluster", "DM1_like"])
        print(f"\n=== DM1 by Paper 1 dark-matter cluster ===")
        print(dm_groups.groupby("v17_dark_cluster")["DM1_like"].agg(["mean","std","count"]).to_string())

    # ===== 4. Survival (PFI = progression-free interval) =====
    surv_ok = merged.dropna(subset=["PFI","PFI.time","DM1_like"])
    print(f"\n=== TCGA-THCA PFI survival (n={len(surv_ok)}) ===")
    try:
        from lifelines import CoxPHFitter
        from lifelines.statistics import logrank_test
        # Cox: DM1_like continuous + age + stage_ord
        cph_df = surv_ok[["PFI.time","PFI","DM1_like"]].rename(columns={"PFI.time":"time","PFI":"event"})
        cph = CoxPHFitter().fit(cph_df, duration_col="time", event_col="event")
        hr = cph.hazard_ratios_["DM1_like"]
        p_cox = cph.summary.loc["DM1_like","p"]
        ci_lo, ci_hi = np.exp(cph.confidence_intervals_.loc["DM1_like"]).values
        print(f"  Cox HR for DM1_like (per 1 SD): HR = {hr:.2f}, 95% CI [{ci_lo:.2f}, {ci_hi:.2f}], p = {p_cox:.2e}")
        # KM by DM1 tertile
        surv_ok = surv_ok.copy()
        surv_ok["DM1_tertile"] = pd.qcut(surv_ok["DM1_like"], 3, labels=["DM1_low","DM1_mid","DM1_high"])
        lo = surv_ok[surv_ok["DM1_tertile"] == "DM1_low"]
        hi = surv_ok[surv_ok["DM1_tertile"] == "DM1_high"]
        lr = logrank_test(lo["PFI.time"], hi["PFI.time"], lo["PFI"], hi["PFI"])
        print(f"  Log-rank DM1_low vs DM1_high: p = {lr.p_value:.2e}")
    except ImportError:
        print("  lifelines not installed, skipping survival")

    # save
    merged.to_csv(OUT / "s_tcga_thca_scored.tsv", sep="\t", index=False)
    print(f"\n→ {OUT / 's_tcga_thca_scored.tsv'}")

    # KM plot
    try:
        from lifelines import KaplanMeierFitter
        fig, ax = plt.subplots(figsize=(8, 6))
        for label, sub, color in [("DM1_low (T1)", lo, "#3C6B4F"),
                                   ("DM1_mid (T2)", surv_ok[surv_ok["DM1_tertile"]=="DM1_mid"], "#B8893C"),
                                   ("DM1_high (T3)", hi, "#962E2E")]:
            kmf = KaplanMeierFitter()
            kmf.fit(sub["PFI.time"], sub["PFI"], label=f"{label} (n={len(sub)})")
            kmf.plot_survival_function(ax=ax, color=color, ci_show=True)
        ax.set_xlabel("PFI time (days)"); ax.set_ylabel("Progression-free survival")
        ax.set_title(f"TCGA-THCA PFI by DM1_like tertile (n={len(surv_ok)})\n"
                     f"Cox HR per 1 SD = {hr:.2f} [{ci_lo:.2f}, {ci_hi:.2f}], p = {p_cox:.2e}",
                     fontsize=11)
        ax.legend(loc="lower left", frameon=False)
        fig.tight_layout()
        fig.savefig(OUT / "s_tcga_km_dm1.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"→ {OUT / 's_tcga_km_dm1.png'}")
    except Exception as e:
        print(f"KM plot failed: {e}")


if __name__ == "__main__":
    main()
