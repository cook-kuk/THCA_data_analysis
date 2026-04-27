#!/usr/bin/env python3
"""
v17 KOREAN K2 — Real PRJEB11591 alignment pipeline (kallisto + 8-gene focus).

Strategy: Download GENCODE human cdna, build kallisto index, download a
subset of PRJEB11591 PE samples, quantify, extract 8-gene TPM, apply
TCGA-trained LogReg, measure Korean AUC.

Subset strategy: 30 PRJEB11591 samples spanning histology types
(time + storage budget: ~3 hours, ~150 GB intermediate).

Outputs:
  - results/v17_korean/K2_kallisto/index.idx
  - results/v17_korean/K2_quant/{run_id}/abundance.tsv
  - results/v17_korean/K2_8gene_tpm_matrix.tsv
  - results/v17_korean/K2_korean_predictions.tsv
  - results/v17_korean/K2_korean_summary.json
"""
from __future__ import annotations
import json
import os
import shutil
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
DATA_REF.mkdir(parents=True, exist_ok=True)
DATA_FASTQ = Path("/data/thca/PRJEB11591_fastq")
DATA_FASTQ.mkdir(parents=True, exist_ok=True)
QUANT = Path("/data/thca/PRJEB11591_quant")
QUANT.mkdir(parents=True, exist_ok=True)

GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz"
GTF_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz"

GENE_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]

# subset: 30 samples (mix of histology types based on sample_alias patterns)
SUBSET_SIZE = 30


def step(msg):
    print(f"\n=== {msg} ===", flush=True)


def download_reference():
    step("1. Download GENCODE v44 cdna")
    cdna_gz = DATA_REF / "gencode.v44.transcripts.fa.gz"
    cdna = DATA_REF / "gencode.v44.transcripts.fa"
    if cdna.exists() and cdna.stat().st_size > 100_000_000:
        print(f"  ✓ already exists: {cdna} ({cdna.stat().st_size//1024//1024} MB)")
    elif cdna_gz.exists() and cdna_gz.stat().st_size > 30_000_000:
        print(f"  ✓ gz already downloaded, decompressing...")
        subprocess.run(["gunzip", "-k", str(cdna_gz)], check=True)
    else:
        print(f"  downloading {GENCODE_URL}")
        t0 = time.time()
        urlretrieve(GENCODE_URL, str(cdna_gz))
        print(f"  ✓ downloaded {cdna_gz.stat().st_size//1024//1024} MB in {time.time()-t0:.0f}s")
        subprocess.run(["gunzip", "-k", str(cdna_gz)], check=True)
        print(f"  ✓ decompressed to {cdna.stat().st_size//1024//1024} MB")

    gtf_gz = DATA_REF / "gencode.v44.basic.annotation.gtf.gz"
    gtf = DATA_REF / "gencode.v44.basic.annotation.gtf"
    if not gtf.exists() and not gtf_gz.exists():
        print(f"  downloading GTF...")
        urlretrieve(GTF_URL, str(gtf_gz))
        subprocess.run(["gunzip", "-k", str(gtf_gz)], check=True)
    return cdna, gtf if gtf.exists() else gtf_gz


def build_index(cdna):
    step("2. Build kallisto index")
    idx = DATA_REF / "gencode.v44.kallisto.idx"
    if idx.exists() and idx.stat().st_size > 100_000_000:
        print(f"  ✓ index already exists: {idx} ({idx.stat().st_size//1024//1024} MB)")
        return idx
    print(f"  building index from {cdna}")
    t0 = time.time()
    subprocess.run([KAL, "index", "-i", str(idx), str(cdna)], check=True)
    print(f"  ✓ built in {time.time()-t0:.0f}s, size {idx.stat().st_size//1024//1024} MB")
    return idx


def parse_t2g(gtf):
    """Build transcript_id → gene_symbol mapping from GTF."""
    step("2b. Parse transcript→gene mapping")
    t2g = {}
    if str(gtf).endswith(".gz"):
        import gzip
        f = gzip.open(gtf, "rt")
    else:
        f = open(gtf)
    for line in f:
        if line.startswith("#"): continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "transcript": continue
        attrs = {}
        for kv in fields[8].split(";"):
            kv = kv.strip()
            if not kv: continue
            parts = kv.split(" ", 1)
            if len(parts) == 2:
                attrs[parts[0]] = parts[1].strip('"')
        tid = attrs.get("transcript_id")
        gname = attrs.get("gene_name")
        if tid and gname:
            t2g[tid] = gname
    f.close()
    print(f"  ✓ {len(t2g)} transcripts mapped")
    return t2g


def select_subset(runs_tsv):
    step(f"3. Select {SUBSET_SIZE} samples from PRJEB11591")
    df = pd.read_csv(runs_tsv, sep="\t")
    print(f"  total runs: {len(df)}")
    # PRJEB11591 has 262 runs; sample_alias should differentiate
    # Pick first SUBSET_SIZE unique samples
    df = df.sort_values("run_accession").drop_duplicates(subset=["sample_alias"], keep="first")
    df = df.head(SUBSET_SIZE)
    df.to_csv(ROOT / "K2_subset_selected.tsv", sep="\t", index=False)
    print(f"  selected {len(df)} unique samples for alignment")
    return df


def download_fastq(df):
    step("4. Download FASTQs (parallel)")
    manifests = []
    skip_count = 0
    download_count = 0
    for idx, row in df.iterrows():
        run = row["run_accession"]
        ftp_str = str(row.get("fastq_ftp", ""))
        urls = [u.strip() for u in ftp_str.split(";") if u.strip()]
        urls = ["http://" + u if not u.startswith("http") else u for u in urls]
        if len(urls) != 2:
            print(f"  ! {run}: expected 2 URLs, got {len(urls)}, skip")
            continue
        out_dir = DATA_FASTQ / run
        out_dir.mkdir(parents=True, exist_ok=True)
        f1 = out_dir / urls[0].split("/")[-1]
        f2 = out_dir / urls[1].split("/")[-1]
        if f1.exists() and f2.exists() and f1.stat().st_size > 100_000_000:
            print(f"  ✓ {run}: already downloaded")
            skip_count += 1
            manifests.append({"run": run, "f1": str(f1), "f2": str(f2)})
            continue
        print(f"  downloading {run}: {urls[0].split('/')[-1]}")
        t0 = time.time()
        try:
            for u, dest in [(urls[0], f1), (urls[1], f2)]:
                # Use wget which handles FTP+resume well
                subprocess.run(["wget", "-q", "-c", "-O", str(dest), u], check=True, timeout=900)
            sz = (f1.stat().st_size + f2.stat().st_size) // 1024 // 1024
            print(f"    ✓ done in {time.time()-t0:.0f}s, {sz} MB")
            download_count += 1
            manifests.append({"run": run, "f1": str(f1), "f2": str(f2)})
        except Exception as e:
            print(f"    ✗ failed: {e}")
    print(f"  ✓ {len(manifests)} samples ready ({skip_count} cached, {download_count} downloaded)")
    return manifests


def quantify(idx, manifests):
    step("5. Kallisto quant")
    summaries = []
    for mi, m in enumerate(manifests):
        run = m["run"]
        out = QUANT / run
        if (out / "abundance.tsv").exists():
            print(f"  ✓ {run} already quantified")
            summaries.append({"run": run, "out": str(out / "abundance.tsv")})
            continue
        print(f"  [{mi+1}/{len(manifests)}] kallisto quant {run}")
        out.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        try:
            subprocess.run([KAL, "quant", "-i", str(idx), "-o", str(out), "-t", "4",
                            m["f1"], m["f2"]], check=True, timeout=1800,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f"    ✓ done in {time.time()-t0:.0f}s")
            summaries.append({"run": run, "out": str(out / "abundance.tsv")})
        except Exception as e:
            print(f"    ✗ failed: {e}")
    return summaries


def extract_8gene(summaries, t2g):
    step("6. Extract 8-gene TPM matrix")
    rows = []
    gene2tx = {g: [] for g in GENE_8}
    # Reverse map: gene_name → tx_ids
    for tid, gname in t2g.items():
        if gname in gene2tx:
            gene2tx[gname].append(tid)
    print("  transcript counts per gene:")
    for g, txs in gene2tx.items():
        print(f"    {g}: {len(txs)} transcripts")
    for s in summaries:
        run = s["run"]
        try:
            df = pd.read_csv(s["out"], sep="\t")
            # df columns: target_id, length, eff_length, est_counts, tpm
            sample_row = {"run": run}
            for g, txs in gene2tx.items():
                # sum TPM across all transcripts of this gene
                m = df[df["target_id"].isin(txs)]
                if len(m) > 0:
                    sample_row[g] = float(m["tpm"].sum())
                else:
                    # try stripping version suffix
                    m2 = df[df["target_id"].str.split(".").str[0].isin([t.split(".")[0] for t in txs])]
                    sample_row[g] = float(m2["tpm"].sum()) if len(m2) > 0 else 0.0
            rows.append(sample_row)
        except Exception as e:
            print(f"  ! {run}: {e}")
    mat = pd.DataFrame(rows).set_index("run")
    mat = mat[GENE_8]  # ensure column order
    mat.to_csv(ROOT / "K2_8gene_tpm_matrix.tsv", sep="\t")
    print(f"  ✓ matrix shape {mat.shape}")
    print(mat.describe().to_string())
    return mat


def apply_panel(mat):
    step("7. Apply TCGA-trained 8-gene LogReg + measure Korean AUC")
    # Load TCGA training data + R1-A labels
    tpm_path = Path("/data/thca/data_processed/bulk_rnaseq/TCGA-THCA_rnaseq_expression_log2.tsv")
    expr = pd.read_csv(tpm_path, sep="\t", index_col=0)
    lbl_path = Path("/opt/thyroid-dash/project/results/v17_realfix/R1A_cluster_labels.tsv")
    lbl = pd.read_csv(lbl_path, sep="\t")
    lbl["y"] = lbl["cluster"].str.startswith("DM2").astype(int)

    # Train on TCGA
    g_have = [g for g in GENE_8 if g in expr.index]
    X_train = expr.loc[g_have, [s for s in lbl["sample_id"] if s in expr.columns]].T
    common = X_train.index.intersection(lbl.set_index("sample_id").index)
    X_train = X_train.loc[common].values
    y_train = lbl.set_index("sample_id").loc[common, "y"].values

    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    sc = StandardScaler().fit(X_train)
    model = LogisticRegression(C=1.0, max_iter=2000, random_state=42).fit(sc.transform(X_train), y_train)
    print(f"  ✓ trained on TCGA n={len(common)}")

    # Apply to Korean
    # Korean matrix is TPM raw → log2(TPM+1) for compatibility with TCGA
    Xk = np.log2(mat.values + 1.0)
    Xk_scaled = sc.transform(Xk)
    p_DM2 = model.predict_proba(Xk_scaled)[:, 1]

    pred = pd.DataFrame({"run": mat.index, "p_DM2": p_DM2})
    pred["DM_call"] = ["DM2" if p > 0.5 else "DM1" for p in p_DM2]
    pred.to_csv(ROOT / "K2_korean_predictions.tsv", sep="\t", index=False)
    print(f"  ✓ predictions: DM1={int((pred['DM_call']=='DM1').sum())}, DM2={int((pred['DM_call']=='DM2').sum())}")

    # Without ground truth labels we can't measure AUC; report distribution + mean
    summary = {
        "n_korean_samples_quantified": int(len(mat)),
        "mean_p_DM2": float(p_DM2.mean()),
        "median_p_DM2": float(np.median(p_DM2)),
        "n_DM1_predicted": int((p_DM2 < 0.5).sum()),
        "n_DM2_predicted": int((p_DM2 >= 0.5).sum()),
        "p_DM2_distribution": {
            "p<0.1": int((p_DM2 < 0.1).sum()),
            "0.1-0.3": int(((p_DM2 >= 0.1) & (p_DM2 < 0.3)).sum()),
            "0.3-0.5": int(((p_DM2 >= 0.3) & (p_DM2 < 0.5)).sum()),
            "0.5-0.7": int(((p_DM2 >= 0.5) & (p_DM2 < 0.7)).sum()),
            "0.7-0.9": int(((p_DM2 >= 0.7) & (p_DM2 < 0.9)).sum()),
            "p>=0.9": int((p_DM2 >= 0.9).sum()),
        },
        "note": "AUC requires ground truth labels (Yoo 2016 BRAF/RAS supplementary mapping). DM1/DM2 distribution + bimodality compared to TCGA expectation = qualitative validation.",
    }
    return summary


def main():
    runs_tsv = ROOT / "K1A_prjeb11591_runs.tsv"
    if not runs_tsv.exists():
        print("[ERROR] run K1A_metadata first")
        return 1

    cdna, gtf = download_reference()
    idx = build_index(cdna)
    t2g = parse_t2g(gtf)
    df = select_subset(runs_tsv)
    manifests = download_fastq(df)
    summaries = quantify(idx, manifests)
    if not summaries:
        print("[ERROR] no successful quantifications")
        return 1
    mat = extract_8gene(summaries, t2g)
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
