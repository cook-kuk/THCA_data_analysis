#!/usr/bin/env python3
"""
TCGA-THCA multivariate Cox + DM1/DM2 stratified CCI.

Adds two paper-grade strengthening pieces:

(A) Multivariate Cox DSS — adjust M1/M2 + CD8 protective effect by age + stage + RAI score.
    Tests whether the immune protective signal survives canonical confounders.

(B) DM1/DM2 stratified CCI bulk replicate — does CAF/M1-M2/TLS × RAI negative
    correlation hold inside each dark-matter cluster (DM1 / DM2 / unlabeled)?
    If yes: avoidance is universal. If only DM1: HT-like-specific.

Inputs: pancan_geneExp.gz + survival.tsv + phenotype + p2d_per_sample_classification.tsv

Outputs:
  spark_tcga_multivariate_cox_DSS.tsv
  spark_tcga_DM_stratified_CCI.tsv
"""
from __future__ import annotations
import gzip, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent.parent.parent
PANCAN = ROOT / "project/data/raw/TCGA_pancan/pancan_geneExp.gz"
PHENO = ROOT / "project/data/raw/TCGA_pancan/phenotype.tsv.gz"
SURV = ROOT / "project/data/raw/TCGA_pancan/survival.tsv"
DM = ROOT / "project/results/dark_matter_phase2/p2d_per_sample_classification.tsv"
RES = ROOT / "project/results/03_pathology_poc"

MODULES = {
    "CAF":   ["FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "DCN", "POSTN"],
    "M1_M2": ["CD68", "CD163", "MRC1", "MARCO", "C1QA", "C1QB", "C1QC", "TYROBP"],
    "CD8":   ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG"],
    "TLS":   ["CXCL13", "CCL19", "CCR7", "LTB", "LTA"],
    "RAI":   ["TPO", "DIO1", "TSHR", "PAX8", "TG", "FOXE1", "NKX2-1", "SLC5A5"],
}


def load_thca_expr(target_genes):
    pheno = pd.read_csv(PHENO, sep="\t", compression="gzip")
    proj_col = next(c for c in pheno.columns if "primary" in c.lower() or "disease" in c.lower())
    thca = set(pheno["sample"][pheno[proj_col].astype(str).str.contains("thyroid|THCA", case=False, na=False)].tolist())
    rows = []
    target_set = set(target_genes)
    with gzip.open(PANCAN, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        keep_idx = [i for i, c in enumerate(header) if c in thca]
        keep_names = [header[i] for i in keep_idx]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in target_set:
                vals = [float(parts[i]) if parts[i] not in {"", "NA"} else np.nan for i in keep_idx]
                rows.append([parts[0]] + vals)
            if len(rows) >= len(target_set):
                break
    df = pd.DataFrame(rows, columns=["gene"] + keep_names).set_index("gene")
    return df


def stage_to_int(s):
    if not isinstance(s, str): return np.nan
    s = s.upper().replace("STAGE ", "")
    mapping = {"I":1, "IA":1, "IB":1, "II":2, "IIA":2, "IIB":2, "IIC":2,
               "III":3, "IIIA":3, "IIIB":3, "IIIC":3, "IV":4, "IVA":4, "IVB":4, "IVC":4}
    return mapping.get(s, np.nan)


def main():
    all_genes = list({g for v in MODULES.values() for g in v})
    expr = load_thca_expr(all_genes)
    print(f"loaded expr: {expr.shape}")

    # module scores
    mods = {}
    for name, genes in MODULES.items():
        present = [g for g in genes if g in expr.index]
        if len(present) < 2: continue
        sub = expr.loc[present]
        z = sub.subtract(sub.mean(axis=1), axis=0).div(sub.std(axis=1) + 1e-6, axis=0)
        mods[name] = z.mean(axis=0)

    # survival + phenotype
    surv = pd.read_csv(SURV, sep="\t").set_index("sample")
    pheno = pd.read_csv(PHENO, sep="\t", compression="gzip").set_index("sample")
    samples = list(set(surv.index) & set(pheno.index) & set(expr.columns))
    print(f"samples with surv+pheno+expr: {len(samples)}")

    # ============================================================
    # (A) Multivariate Cox DSS adjusting for age + stage
    # ============================================================
    from lifelines import CoxPHFitter
    df = pd.DataFrame(index=samples)
    df["time"] = surv.loc[samples, "DSS.time"].astype(float)
    df["event"] = surv.loc[samples, "DSS"].astype(float)
    df["age"] = surv.loc[samples, "age_at_initial_pathologic_diagnosis"].astype(float)
    df["stage"] = surv.loc[samples, "ajcc_pathologic_tumor_stage"].apply(stage_to_int)
    for name, mvec in mods.items():
        df[name] = mvec.reindex(samples).values
    print(f"\n=== (A) Multivariate Cox DSS — n={df.dropna().shape[0]} ===")

    rows = []
    # Univariate baseline (already done)
    for module in ["M1_M2", "CD8", "TLS", "CAF", "RAI"]:
        for adjust in [["age", "stage"], ["age", "stage", "RAI"]]:
            sub = df[["time", "event", module] + adjust].dropna()
            try:
                cph = CoxPHFitter(penalizer=0.01).fit(sub, "time", "event")
                if module not in cph.summary.index: continue
                s = cph.summary.loc[module]
                rows.append({
                    "module": module, "adjust": "+".join(adjust),
                    "n": int(len(sub)), "HR": float(s["exp(coef)"]),
                    "p": float(s["p"]),
                    "HR_lo": float(s["exp(coef) lower 95%"]),
                    "HR_hi": float(s["exp(coef) upper 95%"]),
                })
                marker = " ★★" if s["p"] < 0.01 else (" ★" if s["p"] < 0.05 else "")
                print(f"  {module} | adj={'+'.join(adjust):20s} n={len(sub)} HR={float(s['exp(coef)']):.3f} "
                      f"[{float(s['exp(coef) lower 95%']):.2f},{float(s['exp(coef) upper 95%']):.2f}] "
                      f"p={float(s['p']):.4f}{marker}")
            except Exception as e:
                print(f"  {module} | adj={'+'.join(adjust)} fail: {e}")

    pd.DataFrame(rows).to_csv(RES / "spark_tcga_multivariate_cox_DSS.tsv", sep="\t", index=False)
    print(f"wrote {(RES / 'spark_tcga_multivariate_cox_DSS.tsv').relative_to(ROOT)}")

    # ============================================================
    # (B) DM1/DM2 stratified CCI bulk
    # ============================================================
    print(f"\n=== (B) DM1/DM2 stratified CCI ===")
    dm = pd.read_csv(DM, sep="\t")
    # tcga_short typically TCGA-XX-XXXX; expr cols are TCGA-XX-XXXX-01 (sample with -01 suffix)
    # Build sample-to-DM lookup
    samp_to_dm = {}
    for _, row in dm.iterrows():
        case = row["tcga_short"]
        cluster = row.get("v17_dark_cluster")
        if pd.notna(cluster):
            for sname in expr.columns:
                if sname.startswith(case):
                    samp_to_dm[sname] = cluster
    print(f"DM-labeled bulk samples: {len(samp_to_dm)}")

    rows = []
    for cluster in ["DM1", "DM2", "ALL"]:
        if cluster == "ALL":
            ids = list(set(expr.columns) & set(samp_to_dm.keys()))
        else:
            ids = [s for s, c in samp_to_dm.items() if c == cluster]
        if len(ids) < 20: continue
        for n1 in ["CAF", "M1_M2", "TLS", "CD8"]:
            common = mods[n1].index.intersection(mods["RAI"].index).intersection(ids)
            if len(common) < 15: continue
            r = float(spearmanr(mods[n1].loc[common], mods["RAI"].loc[common]).statistic)
            rows.append({"cluster": cluster, "module": n1, "n": len(common), "r_vs_RAI": r})
            print(f"  {cluster:5s} (n={len(common):3d}): {n1} × RAI Spearman = {r:+.3f}")

    pd.DataFrame(rows).to_csv(RES / "spark_tcga_DM_stratified_CCI.tsv", sep="\t", index=False)
    print(f"wrote {(RES / 'spark_tcga_DM_stratified_CCI.tsv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
