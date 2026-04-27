#!/usr/bin/env python3
"""
v17 KOREAN K2 v4 — Single-end mode (FAST, 153s/sample tested).

Process all available R1 FASTQs with kallisto in single-end mode.
Aggregate 8-gene TPM, apply TCGA-trained LogReg, compute Korean predictions.
"""
from __future__ import annotations
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/opt/thyroid-dash/project/results/v17_korean")
KAL = "/opt/thyroid-dash/project/.venv/lib/python3.12/site-packages/kb_python/bins/linux/kallisto/kallisto"
DATA_REF = Path("/data/thca/reference_kallisto")
DATA_FASTQ = Path("/data/thca/PRJEB11591_fastq")
QUANT = Path("/data/thca/PRJEB11591_quant_se")
QUANT.mkdir(parents=True, exist_ok=True)
WORK = Path("/data/thca/k2_work")
(WORK / "tmp").mkdir(parents=True, exist_ok=True)

GENE_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
GENE_8_ALIASES = {"SLC5A5": ["SLC5A5", "NIS"], "NKX2-1": ["NKX2-1", "NKX2_1", "TTF1", "TITF1"],
                  "FOXE1": ["FOXE1", "TTF2", "FKHL15"], "TPO": ["TPO"], "TG": ["TG"],
                  "TSHR": ["TSHR"], "PAX8": ["PAX8"], "DIO1": ["DIO1"]}


def get_t2g():
    fa = DATA_REF / "gencode.v44.8gene.fa"
    t2g = {}
    with open(fa) as f:
        for line in f:
            if line.startswith(">"):
                fields = line[1:].rstrip().split("|")
                if len(fields) >= 6:
                    tid = fields[0]; gname = fields[5]
                    canonical = None
                    for canon, aliases in GENE_8_ALIASES.items():
                        if gname.upper() in {a.upper() for a in aliases}:
                            canonical = canon; break
                    if canonical:
                        t2g[tid] = canonical
    return t2g


def quantify_se(idx, run, f1):
    out = QUANT / run
    out.mkdir(parents=True, exist_ok=True)
    if (out / "abundance.tsv").exists():
        return str(out / "abundance.tsv")
    t0 = time.time()
    try:
        r = subprocess.run([KAL, "quant", "-i", str(idx), "-o", str(out),
                            "--single", "-l", "200", "-s", "30", "-t", "4",
                            f1],
                           cwd=str(WORK), capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            print(f"  ! [{run}] quant rc={r.returncode}: {r.stderr[:200]}", flush=True)
            return None
        print(f"  ✓ [{run}] quantified in {time.time()-t0:.0f}s", flush=True)
        return str(out / "abundance.tsv")
    except Exception as e:
        print(f"  ! [{run}] exception: {e}", flush=True)
        return None


def main():
    idx = DATA_REF / "gencode.v44.8gene.kallisto.idx"
    if not idx.exists():
        print("[ERROR] index missing"); return 1
    t2g = get_t2g()
    print(f"  t2g: {len(t2g)} transcripts")

    # Find all R1 FASTQs that are >500MB (complete)
    r1_files = []
    for d in sorted(DATA_FASTQ.iterdir()):
        if not d.is_dir(): continue
        run = d.name
        f1 = next((f for f in d.iterdir() if f.name.endswith("_1.fastq.gz") and f.stat().st_size > 500_000_000), None)
        if f1:
            r1_files.append((run, str(f1)))
    print(f"  Found {len(r1_files)} R1 FASTQs (>500MB)")

    print("\n=== Sequential single-end kallisto quant ===", flush=True)
    summaries = []
    for i, (run, f1) in enumerate(r1_files):
        print(f"  [{i+1}/{len(r1_files)}] {run}: quant single-end...", flush=True)
        ab = quantify_se(idx, run, f1)
        if ab:
            summaries.append((run, ab))

    print(f"\n  ✓ {len(summaries)}/{len(r1_files)} quantified")

    # Build 8-gene TPM matrix
    print("\n=== 8-gene TPM matrix ===", flush=True)
    rows = []
    for run, ab in summaries:
        try:
            df_q = pd.read_csv(ab, sep="\t")
            sample_row = {"run": run}
            for g in GENE_8:
                txs = [t for t, gn in t2g.items() if gn == g]
                m = df_q[df_q["target_id"].isin(txs)]
                sample_row[g] = float(m["tpm"].sum()) if len(m) > 0 else 0.0
            rows.append(sample_row)
        except Exception as e:
            print(f"  ! {run}: {e}")
    if not rows:
        print("[ERROR] no quantifications"); return 1
    mat = pd.DataFrame(rows).set_index("run")[GENE_8]
    mat.to_csv(ROOT / "K2_8gene_tpm_matrix_v4.tsv", sep="\t")
    print(f"  matrix: {mat.shape}")
    print(mat.describe().round(2).to_string())

    # Apply TCGA LogReg
    print("\n=== Apply TCGA LogReg + Korean predictions ===", flush=True)
    tpm = pd.read_csv("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv",
                     sep="\t", index_col=0)
    lbl = pd.read_csv("/opt/thyroid-dash/project/results/v17_realfix/R1A_cluster_labels.tsv", sep="\t")
    lbl["y"] = lbl["cluster"].str.startswith("DM2").astype(int)
    g_have = [g for g in GENE_8 if g in tpm.index]
    X_tr = tpm.loc[g_have, [s for s in lbl["sample_id"] if s in tpm.columns]].T
    common = X_tr.index.intersection(lbl.set_index("sample_id").index)
    X_tr = X_tr.loc[common].values
    y_tr = lbl.set_index("sample_id").loc[common, "y"].values

    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    sc = StandardScaler().fit(X_tr)
    model = LogisticRegression(C=1.0, max_iter=2000, random_state=42).fit(sc.transform(X_tr), y_tr)
    print(f"  trained on TCGA n={len(common)}")

    Xk = np.log2(mat.values + 1.0)
    Xk_s = sc.transform(Xk)
    p = model.predict_proba(Xk_s)[:, 1]
    pred = pd.DataFrame({"run": mat.index, "p_DM2": p})
    pred["DM_call"] = ["DM2" if x > 0.5 else "DM1" for x in p]
    pred = pd.concat([pred.reset_index(drop=True), mat.reset_index(drop=True)], axis=1)
    pred = pred.loc[:, ~pred.columns.duplicated()]
    pred.to_csv(ROOT / "K2_korean_predictions_v4.tsv", sep="\t", index=False)
    print(f"  predictions: DM1={int((p<0.5).sum())}, DM2={int((p>=0.5).sum())}")

    # Distribution check
    summary = {
        "method": "Single-end kallisto pseudoalignment (-l 200 -s 30) on full R1 FASTQ; 8-gene-only index (93 transcripts)",
        "n_korean_samples": int(len(mat)),
        "panel_genes": GENE_8,
        "tpm_summary": mat.describe().round(2).to_dict(),
        "n_DM1": int((p < 0.5).sum()),
        "n_DM2": int((p >= 0.5).sum()),
        "mean_p_DM2": float(p.mean()),
        "median_p_DM2": float(np.median(p)),
        "p_DM2_distribution": {
            "p<0.1": int((p < 0.1).sum()),
            "0.1-0.3": int(((p >= 0.1) & (p < 0.3)).sum()),
            "0.3-0.5": int(((p >= 0.3) & (p < 0.5)).sum()),
            "0.5-0.7": int(((p >= 0.5) & (p < 0.7)).sum()),
            "0.7-0.9": int(((p >= 0.7) & (p < 0.9)).sum()),
            "p>=0.9": int((p >= 0.9).sum()),
        },
        "tcga_train_dm1_dm2": f"{int((y_tr==0).sum())}:{int((y_tr==1).sum())}",
        "korean_dm1_dm2_ratio": f"{int((p<0.5).sum())}:{int((p>=0.5).sum())}",
        "interpretation": "Cross-population prediction. DM2-skewed = preserved thyroid differentiation (NIS/TPO/TG/TSHR/DIO1 high); DM1-skewed = de-differentiated. Korean cohort distribution can be compared to TCGA (DM1 28%/DM2 72% in R1-A leak-free).",
        "expected": "If panel transfers correctly (Tier 1), Korean cohort should show similar DM1/DM2 ratio to TCGA (~30% DM1 / 70% DM2 for primary indolent PTC). Yoo 2016 BRAF rate 62-70% suggests slightly more BRAF-DM1 enrichment expected.",
    }
    with open(ROOT / "K2_korean_summary_v4.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print("\n=== K2 v4 SUMMARY ===")
    print(json.dumps(summary, indent=2, default=str)[:2500])
    return 0


if __name__ == "__main__":
    sys.exit(main())
