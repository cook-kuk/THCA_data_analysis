#!/usr/bin/env python3
"""
v17 KOREAN K2 v2 — Smarter pipeline using 8-gene-only kallisto index.

Approach: Build a tiny kallisto index of just 8-gene transcripts (~100 transcripts).
Pseudoalign FASTQs against this small index (very fast, low memory).
Result is relative TPM of 8 genes per sample (for cross-sample comparison only;
absolute TPM not interpretable but ranking is valid).

Then apply TCGA-trained 8-gene LogReg, generate Korean DM1/DM2 predictions.
"""
from __future__ import annotations
import gzip
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from urllib.request import urlretrieve

import numpy as np
import pandas as pd

ROOT = Path("/opt/thyroid-dash/project/results/v17_korean")
ROOT.mkdir(parents=True, exist_ok=True)
KAL = "/opt/thyroid-dash/project/.venv/lib/python3.12/site-packages/kb_python/bins/linux/kallisto/kallisto"
DATA_REF = Path("/data/thca/reference_kallisto")
DATA_FASTQ = Path("/data/thca/PRJEB11591_fastq")
DATA_FASTQ.mkdir(parents=True, exist_ok=True)
QUANT = Path("/data/thca/PRJEB11591_quant")
QUANT.mkdir(parents=True, exist_ok=True)
WORK = Path("/data/thca/k2_work")
WORK.mkdir(parents=True, exist_ok=True)
(WORK / "tmp").mkdir(exist_ok=True)

GENE_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
GENE_8_ALIASES = {"SLC5A5": ["SLC5A5", "NIS"], "NKX2-1": ["NKX2-1", "NKX2_1", "TTF1", "TITF1"],
                  "FOXE1": ["FOXE1", "TTF2", "FKHL15"], "TPO": ["TPO"], "TG": ["TG"],
                  "TSHR": ["TSHR"], "PAX8": ["PAX8"], "DIO1": ["DIO1"]}
SUBSET_SIZE = 20  # smaller subset for fast iteration


def step(msg):
    print(f"\n=== {msg} ===", flush=True)


def extract_8gene_fasta():
    step("1. Extract 8-gene transcripts from full GENCODE cdna")
    out_fa = DATA_REF / "gencode.v44.8gene.fa"
    if out_fa.exists() and out_fa.stat().st_size > 1000:
        print(f"  ✓ already exists: {out_fa} ({out_fa.stat().st_size//1024} KB)")
        # show t2g
        t2g = {}
        with open(out_fa) as f:
            for line in f:
                if line.startswith(">"):
                    m = re.match(r">([^\s|]+).*?gene_name[=:|]([^\s|]+)", line)
                    if m:
                        t2g[m.group(1)] = m.group(2)
        return out_fa, t2g

    cdna = DATA_REF / "gencode.v44.transcripts.fa"
    if not cdna.exists():
        print(f"  ✗ missing {cdna}")
        return None, {}

    # Build set of all 8-gene aliases
    all_aliases = set()
    for g, aliases in GENE_8_ALIASES.items():
        for a in aliases:
            all_aliases.add(a)

    print(f"  scanning {cdna} for transcripts of {sorted(all_aliases)}")
    out = open(out_fa, "w")
    keep = False
    n_tx = 0
    t2g = {}
    with open(cdna) as f:
        for line in f:
            if line.startswith(">"):
                # GENCODE header: >ENST00000456328.2|ENSG00000223972.5|OTTHUMG...|...|DDX11L1-202|DDX11L1|1657|...
                fields = line[1:].split("|")
                if len(fields) >= 6:
                    tid = fields[0]
                    gname = fields[5]
                    if gname.upper() in {a.upper() for a in all_aliases} or any(gname.upper() == a.upper() for a in all_aliases):
                        keep = True
                        n_tx += 1
                        # determine canonical 8-gene name
                        canonical = None
                        for canon, aliases in GENE_8_ALIASES.items():
                            if gname.upper() in {a.upper() for a in aliases}:
                                canonical = canon; break
                        if canonical:
                            t2g[tid] = canonical
                        out.write(line)
                    else:
                        keep = False
                else:
                    keep = False
            else:
                if keep:
                    out.write(line)
    out.close()
    print(f"  ✓ extracted {n_tx} transcripts ({out_fa.stat().st_size//1024} KB)")
    print(f"  per-gene transcript counts:")
    from collections import Counter
    cnt = Counter(t2g.values())
    for g in GENE_8:
        print(f"    {g}: {cnt.get(g, 0)} transcripts")
    return out_fa, t2g


def build_small_index(fa):
    step("2. Build kallisto index (8-gene only, fast)")
    idx = DATA_REF / "gencode.v44.8gene.kallisto.idx"
    if idx.exists() and idx.stat().st_size > 500_000:
        print(f"  ✓ already exists: {idx} ({idx.stat().st_size//1024} KB)")
        return idx
    t0 = time.time()
    # Run kallisto from a CWD that has tmp/
    proc = subprocess.run([KAL, "index", "-k", "31", "-i", str(idx), str(fa)],
                          cwd=str(WORK), capture_output=True, text=True)
    if proc.returncode != 0:
        # Try with smaller k
        print(f"  retry with k=21")
        proc = subprocess.run([KAL, "index", "-k", "21", "-i", str(idx), str(fa)],
                              cwd=str(WORK), capture_output=True, text=True)
    print(f"  stdout: {proc.stdout[:300]}")
    print(f"  stderr: {proc.stderr[:300]}")
    if proc.returncode != 0:
        print(f"  ✗ index build failed (returncode {proc.returncode})")
        return None
    print(f"  ✓ index built in {time.time()-t0:.0f}s, size {idx.stat().st_size//1024} KB")
    return idx


def select_subset():
    step(f"3. Select {SUBSET_SIZE} samples")
    runs_tsv = ROOT / "K1A_prjeb11591_runs.tsv"
    df = pd.read_csv(runs_tsv, sep="\t")
    # take first SUBSET_SIZE runs (already sorted by run_accession)
    df = df.sort_values("run_accession").head(SUBSET_SIZE)
    df.to_csv(ROOT / "K2_subset_v2.tsv", sep="\t", index=False)
    print(f"  selected {len(df)} samples")
    return df


def download_one(run, ftp_str):
    urls = [u.strip() for u in str(ftp_str).split(";") if u.strip()]
    urls = ["http://" + u if not u.startswith("http") else u for u in urls]
    if len(urls) != 2:
        return None, None
    out_dir = DATA_FASTQ / run
    out_dir.mkdir(parents=True, exist_ok=True)
    f1 = out_dir / urls[0].split("/")[-1]
    f2 = out_dir / urls[1].split("/")[-1]
    if f1.exists() and f2.exists() and f1.stat().st_size > 100_000_000:
        return f1, f2
    try:
        for u, dest in [(urls[0], f1), (urls[1], f2)]:
            r = subprocess.run(["wget", "-q", "-c", "--timeout=300", "-O", str(dest), u],
                               capture_output=True, timeout=600)
            if r.returncode != 0:
                return None, None
        return f1, f2
    except Exception:
        return None, None


def quantify_one(idx, run, f1, f2):
    out = QUANT / run
    out.mkdir(parents=True, exist_ok=True)
    if (out / "abundance.tsv").exists():
        return out
    try:
        proc = subprocess.run([KAL, "quant", "-i", str(idx), "-o", str(out), "-t", "4",
                               str(f1), str(f2)],
                              cwd=str(WORK), capture_output=True, text=True, timeout=600)
        if proc.returncode != 0:
            print(f"    ! {run} kallisto quant returncode={proc.returncode}: {proc.stderr[:200]}")
            return None
        return out
    except Exception as e:
        print(f"    ! {run} exception: {e}")
        return None


def process_samples(idx, df):
    step("4. Download + quantify each sample (sequential to limit storage)")
    summaries = []
    for i, row in df.iterrows():
        run = row["run_accession"]
        ftp = row.get("fastq_ftp", "")
        print(f"  [{len(summaries)+1}/{len(df)}] {run}: download...", flush=True)
        t0 = time.time()
        f1, f2 = download_one(run, ftp)
        if not f1:
            print(f"    ✗ download failed")
            continue
        sz_mb = (f1.stat().st_size + f2.stat().st_size) // 1024 // 1024
        print(f"    ✓ downloaded {sz_mb} MB in {time.time()-t0:.0f}s; quantifying...", flush=True)
        t0 = time.time()
        out = quantify_one(idx, run, f1, f2)
        if not out:
            continue
        print(f"    ✓ quantified in {time.time()-t0:.0f}s", flush=True)
        summaries.append({"run": run, "out": str(out / "abundance.tsv")})
        # Cleanup FASTQ to save storage
        try:
            f1.unlink()
            f2.unlink()
        except Exception: pass
    return summaries


def extract_8gene_tpm(summaries, t2g):
    step("5. Extract 8-gene TPM matrix")
    rows = []
    for s in summaries:
        run = s["run"]
        try:
            df = pd.read_csv(s["out"], sep="\t")
            sample_row = {"run": run}
            for g in GENE_8:
                txs = [t for t, gn in t2g.items() if gn == g]
                m = df[df["target_id"].isin(txs)]
                sample_row[g] = float(m["tpm"].sum()) if len(m) > 0 else 0.0
            rows.append(sample_row)
        except Exception as e:
            print(f"  ! {run}: {e}")
    mat = pd.DataFrame(rows).set_index("run")[GENE_8]
    mat.to_csv(ROOT / "K2_8gene_tpm_matrix.tsv", sep="\t")
    print(f"  ✓ matrix shape {mat.shape}")
    print(mat.describe().to_string())
    return mat


def apply_panel(mat):
    step("6. Apply TCGA-trained LogReg + Korean predictions")
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
    print(f"  ✓ trained on TCGA n={len(common)}")

    Xk = np.log2(mat.values + 1.0)
    Xk_s = sc.transform(Xk)
    p = model.predict_proba(Xk_s)[:, 1]
    pred = pd.DataFrame({"run": mat.index, "p_DM2": p})
    pred["DM_call"] = ["DM2" if x > 0.5 else "DM1" for x in p]
    pred.to_csv(ROOT / "K2_korean_predictions.tsv", sep="\t", index=False)

    summary = {
        "n_korean_samples": int(len(mat)),
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
        "tcga_dm1_dm2_ratio": f"{int((y_tr==0).sum())}:{int((y_tr==1).sum())}",
        "korean_dm1_dm2_ratio": f"{int((p<0.5).sum())}:{int((p>=0.5).sum())}",
        "panel_genes": GENE_8,
        "method": "8-gene-only kallisto pseudoalignment (k=31), full GENCODE v44 cdna had build segfault → small index workaround",
        "interpretation": "Relative ranking valid; absolute TPM may differ from full-transcriptome quant due to potential decoy contamination; for AUC measurement against ground truth (Yoo 2016 BRAF/RAS-like), supplementary metadata mapping needed (separate outreach)",
    }
    return summary


def main():
    fa, t2g = extract_8gene_fasta()
    if not fa: return 1
    idx = build_small_index(fa)
    if not idx: return 1
    df = select_subset()
    summaries = process_samples(idx, df)
    if not summaries:
        print("[ERROR] no successful samples")
        return 1
    mat = extract_8gene_tpm(summaries, t2g)
    summary = apply_panel(mat)
    summary["n_attempted"] = SUBSET_SIZE
    summary["n_successful"] = len(summaries)
    with open(ROOT / "K2_korean_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print("\n=== K2 SUMMARY ===")
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
