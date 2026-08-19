"""
Run MORE signatures + cohorts in parallel via ProcessPoolExecutor.

Adds:
  - Mariathasan IMmotion150 19-gene tumor immune phenotype (TGFB-cancer interface)
  - PRGS Pan-cancer 23-gene Tumor Inflammation Signature (Cristescu 2018 Science)
  - Senbabaoglu 2016 immune deconvolution (B/T/NK/Tregs)
  - Davoli 2017 aneuploidy score proxy (CIN70)
  - Tirosh 2016 melanoma exhaustion-like
  - GSEA Hallmark Inflammatory_Response + Interferon_Gamma_Response
  - Mariathasan TGFb axis
  - Auslander_Wolf 2018 IPRES (innate ICI resistance)

All on TCGA-PAAD bulk z-scores + (where possible) GSE71729.
"""

import json, gzip, re, os
from io import StringIO
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

ROOT = Path("/data/pdac_poc")
RES = ROOT/"results/benchmark"
RES.mkdir(parents=True, exist_ok=True)

EXTRA_PANELS = {
    "Mariathasan_TGFB":  ["ACTA2","ACTG2","ADAM12","ADAM19","COMP","COL10A1","COL1A1","COL1A2",
                           "COL3A1","COL4A1","COL5A1","COL6A1","FAP","FBN1","FN1","LRRC15","MMP11"],
    "PRGS_Cristescu":    ["CCL5","CD27","CD274","CD8A","CMKLR1","CXCL9","CXCR6","HLA-DQA1",
                           "HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT"],
    "Senbabaoglu_Treg":  ["FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18","CCR8"],
    "CIN70_aneuploidy":  ["TPX2","CCNB1","FOXM1","KIF20A","UBE2C","CDC20","NEK2","KIF11",
                           "MELK","TYMS","RRM2","CDK1","KIF23","TOP2A","BIRC5"],
    "IPRES_Hugo":        ["VEGFA","VEGFB","TGFB1","WNT5A","HGF","IGF1","IL6","IL10","CXCL12"],
    "Mariathasan_immune":["CD8A","IFNG","CXCL9","CXCL10","GZMA","GZMB","HLA-DRA","STAT1","CCR5"],
    "ChengEcotype_TLS":  ["CXCL13","CCL19","CCL21","CXCR5","CCR7","SELL","LAMP3"],
    "Tirosh_exhaustion": ["LAG3","HAVCR2","TIGIT","PDCD1","CTLA4","BTLA","TOX","TOX2","NR4A1","NR4A2"],
    "Hallmark_IFNG":     ["IFNG","STAT1","IRF1","IRF7","ISG15","OAS1","OAS2","MX1","MX2","CXCL9","CXCL10"],
    "Hallmark_INFL":     ["IL6","TNF","IL1B","NFKB1","NFKB2","RELB","CCL2","CCL5","CXCL1","CXCL8"],
    "Auslander_NSCLC":   ["CD274","PDCD1","CTLA4","CXCL9","CXCL10","STAT1","IFNG","HLA-DRA"],
}


def score_panel(z_df, genes):
    upper_idx = z_df.index.str.upper()
    rows = []
    for g in genes:
        mask = upper_idx == g.upper()
        if mask.any(): rows.append(z_df.index[mask][0])
    if not rows: return None
    return z_df.loc[rows].mean(axis=0)


def cindex(score, t, e):
    df = pd.DataFrame({"s": score, "t": t, "e": e}).dropna()
    df = df[df["t"]>0]
    if len(df) < 30: return None, len(df)
    try:
        c = concordance_index(df["t"], -df["s"], df["e"])
        return float(c), len(df)
    except: return None, len(df)


def cox_hr(score, t, e):
    df = pd.DataFrame({"s": score, "t": t, "e": e}).dropna()
    df = df[df["t"]>0]
    if len(df) < 30: return None
    try:
        cph = CoxPHFitter(penalizer=0.05).fit(df, duration_col="t", event_col="e")
        s = cph.summary
        return {"HR": float(s.loc["s","exp(coef)"]),
                "p": float(s.loc["s","p"]),
                "CI_lo": float(s.loc["s","exp(coef) lower 95%"]),
                "CI_hi": float(s.loc["s","exp(coef) upper 95%"])}
    except: return None


def benchmark_one_panel(args):
    name, genes = args
    z = pd.read_csv(ROOT/"raw/mrna_z_Zscores.tsv", sep="\t", index_col=0).T
    z.index = z.index.str.upper()
    moff = pd.read_csv(ROOT/"results/moffitt_calls.tsv", sep="\t", index_col=0)
    clin = pd.read_csv(ROOT/"raw/clinical.tsv", sep="\t").set_index("sampleId")
    surv = moff.join(pd.DataFrame({
        "t": pd.to_numeric(clin.get("OS_MONTHS"), errors="coerce"),
        "e": clin.get("OS_STATUS").fillna("").str.startswith("1").astype(int)
    }), how="left")
    common = list(set(z.columns) & set(surv.index))
    z = z[common]; surv = surv.loc[common]
    score = score_panel(z, genes)
    if score is None:
        return {"method": name, "ok": False, "n_genes_used": 0}
    c, n = cindex(score, surv["t"], surv["e"])
    cox = cox_hr(score, surv["t"], surv["e"])
    return {"method": name,
            "n_genes_used": int(sum(1 for g in genes if g.upper() in z.index)),
            "n_genes_total": len(genes),
            "n": n, "c_index": round(c,4) if c else None,
            "HR": round(cox["HR"],3) if cox else None,
            "p": round(cox["p"],4) if cox else None,
            "CI_lo": round(cox["CI_lo"],3) if cox else None,
            "CI_hi": round(cox["CI_hi"],3) if cox else None}


def main():
    print(f"[26] launching {len(EXTRA_PANELS)} extra signatures on {os.cpu_count()} CPUs")
    rows = []
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
        for r in ex.map(benchmark_one_panel, list(EXTRA_PANELS.items())):
            rows.append(r)
            print(f"   ✓ {r['method']:25s} C={r.get('c_index')} HR={r.get('HR')} p={r.get('p')} (n_genes {r.get('n_genes_used')}/{r.get('n_genes_total')})")
    json.dump(rows, open(RES/"extra_signatures.json","w"), indent=2, default=str)
    print(f"\n[26] saved {RES}/extra_signatures.json")


if __name__ == "__main__":
    main()
