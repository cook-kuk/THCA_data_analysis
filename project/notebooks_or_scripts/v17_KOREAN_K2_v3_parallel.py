#!/usr/bin/env python3
"""
v17 KOREAN K2 v3 — Parallel download + sequential quant.
3 concurrent FASTQ downloads, kallisto quant after each.
Reuses 8-gene index from v2.
"""
from __future__ import annotations
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/opt/thyroid-dash/project/results/v17_korean")
KAL = "/opt/thyroid-dash/project/.venv/lib/python3.12/site-packages/kb_python/bins/linux/kallisto/kallisto"
DATA_REF = Path("/data/thca/reference_kallisto")
DATA_FASTQ = Path("/data/thca/PRJEB11591_fastq")
QUANT = Path("/data/thca/PRJEB11591_quant")
WORK = Path("/data/thca/k2_work")
(WORK / "tmp").mkdir(parents=True, exist_ok=True)
DATA_FASTQ.mkdir(parents=True, exist_ok=True)
QUANT.mkdir(parents=True, exist_ok=True)

GENE_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
GENE_8_ALIASES = {"SLC5A5": ["SLC5A5", "NIS"], "NKX2-1": ["NKX2-1", "NKX2_1", "TTF1", "TITF1"],
                  "FOXE1": ["FOXE1", "TTF2", "FKHL15"], "TPO": ["TPO"], "TG": ["TG"],
                  "TSHR": ["TSHR"], "PAX8": ["PAX8"], "DIO1": ["DIO1"]}
SUBSET_SIZE = 20
PARALLEL = 3


def get_t2g():
    fa = DATA_REF / "gencode.v44.8gene.fa"
    t2g = {}
    with open(fa) as f:
        for line in f:
            if line.startswith(">"):
                fields = line[1:].rstrip().split("|")
                if len(fields) >= 6:
                    tid = fields[0]
                    gname = fields[5]
                    canonical = None
                    for canon, aliases in GENE_8_ALIASES.items():
                        if gname.upper() in {a.upper() for a in aliases}:
                            canonical = canon; break
                    if canonical:
                        t2g[tid] = canonical
    return t2g


def download_one(run, ftp):
    urls = [u.strip() for u in str(ftp).split(";") if u.strip()]
    urls = ["http://" + u if not u.startswith("http") else u for u in urls]
    if len(urls) != 2: return None, None
    out_dir = DATA_FASTQ / run
    out_dir.mkdir(parents=True, exist_ok=True)
    f1 = out_dir / urls[0].split("/")[-1]
    f2 = out_dir / urls[1].split("/")[-1]
    # If both already exist with valid sizes, skip
    if f1.exists() and f2.exists() and f1.stat().st_size > 100_000_000 and f2.stat().st_size > 100_000_000:
        print(f"  [{run}] cached", flush=True)
        return str(f1), str(f2)
    t0 = time.time()
    try:
        for u, dest in [(urls[0], f1), (urls[1], f2)]:
            r = subprocess.run(["wget", "-q", "-c", "--timeout=600", "--tries=2",
                                "-O", str(dest), u],
                               capture_output=True, timeout=1200)
            if r.returncode != 0:
                return None, None
        sz = (f1.stat().st_size + f2.stat().st_size) // 1024 // 1024
        print(f"  ↓ [{run}] downloaded {sz} MB in {time.time()-t0:.0f}s", flush=True)
        return str(f1), str(f2)
    except Exception as e:
        print(f"  ✗ [{run}] download exception: {e}", flush=True)
        return None, None


def quantify_one(idx, run, f1, f2):
    out = QUANT / run
    out.mkdir(parents=True, exist_ok=True)
    if (out / "abundance.tsv").exists():
        return str(out / "abundance.tsv")
    t0 = time.time()
    try:
        r = subprocess.run([KAL, "quant", "-i", str(idx), "-o", str(out), "-t", "4",
                            f1, f2],
                           cwd=str(WORK), capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            print(f"  ! [{run}] quant rc={r.returncode}", flush=True)
            return None
        print(f"  ✓ [{run}] quantified in {time.time()-t0:.0f}s", flush=True)
        try:
            Path(f1).unlink(); Path(f2).unlink()
        except Exception: pass
        return str(out / "abundance.tsv")
    except Exception as e:
        print(f"  ! [{run}] quant exception: {e}", flush=True)
        return None


def main():
    idx = DATA_REF / "gencode.v44.8gene.kallisto.idx"
    if not idx.exists():
        print("[ERROR] index missing")
        return 1
    t2g = get_t2g()
    print(f"  t2g: {len(t2g)} transcripts")

    runs_tsv = ROOT / "K1A_prjeb11591_runs.tsv"
    df = pd.read_csv(runs_tsv, sep="\t").sort_values("run_accession").head(SUBSET_SIZE)
    print(f"  {len(df)} samples queued, parallel={PARALLEL}")

    print("\n=== Parallel pipeline ===", flush=True)
    summaries = []
    pending = [(row["run_accession"], row["fastq_ftp"]) for _, row in df.iterrows()]

    # Use one thread to drive downloads, main loop drives quant
    completed = 0
    futures = {}
    with ThreadPoolExecutor(max_workers=PARALLEL) as exe:
        for run, ftp in pending:
            futures[exe.submit(download_one, run, ftp)] = run
        for fut in as_completed(futures):
            run = futures[fut]
            f1, f2 = fut.result()
            if not f1:
                continue
            ab = quantify_one(idx, run, f1, f2)
            if ab:
                summaries.append((run, ab))
                completed += 1
                print(f"  [{completed}/{len(pending)}] {run} pipeline complete", flush=True)

    print(f"\n  ✓ {len(summaries)}/{len(pending)} samples processed")

    # 5. Extract 8-gene TPM matrix
    print("\n=== Extract 8-gene TPM ===", flush=True)
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
        print("[ERROR] no quantifications")
        return 1
    mat = pd.DataFrame(rows).set_index("run")[GENE_8]
    mat.to_csv(ROOT / "K2_8gene_tpm_matrix.tsv", sep="\t")
    print(f"  matrix shape {mat.shape}")
    print(mat.describe().to_string())

    # 6. Apply LogReg
    print("\n=== Apply TCGA-trained LogReg ===", flush=True)
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
    pred.to_csv(ROOT / "K2_korean_predictions.tsv", sep="\t", index=False)

    summary = {
        "n_korean_samples": int(len(mat)),
        "n_attempted": SUBSET_SIZE,
        "n_successful": int(len(summaries)),
        "n_DM1_predicted": int((p < 0.5).sum()),
        "n_DM2_predicted": int((p >= 0.5).sum()),
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
        "tcga_dm1_dm2_ratio_in_train": f"{int((y_tr==0).sum())}:{int((y_tr==1).sum())}",
        "korean_dm1_dm2_ratio": f"{int((p<0.5).sum())}:{int((p>=0.5).sum())}",
        "panel_genes": GENE_8,
        "method": "8-gene-only kallisto pseudoalignment (k=31, 93 transcripts), 3-way parallel download, full GENCODE v44 cdna had build segfault → small index workaround",
        "interpretation": "Relative TPM ranking valid; DM1/DM2 prediction distribution 비교 가능. AUC requires Yoo 2016 BRAF/RAS-like ground truth labels (revision-round outreach to Yoo SK).",
    }
    with open(ROOT / "K2_korean_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    print("\n=== K2 KOREAN SUMMARY ===")
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
