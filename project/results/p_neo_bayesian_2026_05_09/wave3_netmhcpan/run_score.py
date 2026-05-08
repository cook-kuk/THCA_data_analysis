#!/usr/bin/env python
"""Score NetMHCpan-4.1 on leakage-stratified ITSNdb (wave3 single-algo agent).

Install path: docker image `kevinr9525/netmhcpan-4.1b.linux:latest`, which ships
the DTU NetMHCpan-4.1b binary at `/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan`.

HLA format conversion: bundle uses `HLA-A*02:01`; netMHCpan wants `HLA-A02:01`.

Per-allele batch run: -p (peptide file) + -a HLA-X02:01. Score = -%Rank_EL.
"""
from __future__ import annotations
import json
import re
import subprocess
import time
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
OUT = ROOT / "wave3_netmhcpan"
BUNDLE = ROOT / "bundle.tsv"
DOCKER_IMG = "kevinr9525/netmhcpan-4.1b.linux:latest"
NETMHC_BIN = "/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan"

OUT.mkdir(parents=True, exist_ok=True)
WORK = OUT / "_work"
WORK.mkdir(exist_ok=True)

t0 = time.time()
print("[load] bundle...", flush=True)
df = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master": "boolean"})

splits_keep = {"ext_itsndb_main", "ext_itsndb_val"}
df = df[df["split"].isin(splits_keep)].reset_index(drop=True).copy()
print(f"[load] kept {len(df)} ITSNdb rows", flush=True)
print(f"[load] split counts: {df['split'].value_counts().to_dict()}", flush=True)
print(f"[load] in_master counts: {df['in_master'].value_counts(dropna=False).to_dict()}", flush=True)

df["pep_len"] = df["peptide"].str.len()
mask_len = df["pep_len"].between(8, 11)
n_dropped_len = int((~mask_len).sum())
df = df[mask_len].copy()

AA = set("ACDEFGHIKLMNPQRSTVWY")
mask_aa = df["peptide"].apply(lambda s: set(s).issubset(AA))
n_dropped_aa = int((~mask_aa).sum())
df = df[mask_aa].copy()
print(f"[filter] dropped len={n_dropped_len}, non-std AA={n_dropped_aa}; {len(df)} rows remain", flush=True)


def hla_to_netmhc(s: str) -> str:
    # "HLA-A*02:01" -> "HLA-A02:01"
    return s.replace("*", "")


df["hla_netmhc"] = df["HLA_norm"].apply(hla_to_netmhc)


# Pre-fetch list of supported alleles
print("[netmhc] fetching supported allele list...", flush=True)
res = subprocess.run(
    [
        "docker", "run", "--rm",
        "-e", "NETMHCpan=/netMHCpan-4.1",
        "-e", "TMPDIR=/tmp",
        "--entrypoint", NETMHC_BIN,
        DOCKER_IMG,
        "-listMHC",
    ],
    capture_output=True, text=True, timeout=120,
)
supported = set()
for line in res.stdout.splitlines():
    line = line.strip()
    if line.startswith("HLA-"):
        supported.add(line)
print(f"[netmhc] {len(supported)} alleles supported", flush=True)

allele_skipped = {}


def is_supported(hla: str) -> bool:
    if hla in supported:
        return True
    allele_skipped[hla] = allele_skipped.get(hla, 0) + 1
    return False


df["hla_ok"] = df["hla_netmhc"].apply(is_supported)
n_skipped_total = int((~df["hla_ok"]).sum())
print(f"[netmhc] skipping {n_skipped_total} rows (unsupported HLA)", flush=True)
print(f"[netmhc] skipped allele table: {allele_skipped}", flush=True)

df_score = df[df["hla_ok"]].copy().reset_index(drop=True)
df_score["score_netmhcpan"] = np.nan
df_score["pct_rank_el"] = np.nan
df_score["pct_rank_ba"] = np.nan


# Header in netMHCpan output:
# Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity Score_EL %Rank_EL Score_BA %Rank_BA Aff(nM) [BindLevel]
def parse_netmhcpan_output(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-") or line.startswith("Pos"):
            continue
        if line.startswith("Protein") or line.startswith("Number of"):
            continue
        # Data lines start with a position integer
        parts = re.split(r"\s+", line)
        if not parts or not parts[0].isdigit():
            continue
        # need at least 16 fields for -BA output (BindLevel may be present as 'SB'/'WB' with optional <=)
        if len(parts) < 15:
            continue
        try:
            pos = int(parts[0])
            mhc = parts[1]
            peptide = parts[2]
            score_el = float(parts[11])
            rank_el = float(parts[12])
            score_ba = float(parts[13])
            rank_ba = float(parts[14])
            aff = float(parts[15])
        except (ValueError, IndexError):
            continue
        rows.append({
            "MHC": mhc,
            "peptide": peptide,
            "Score_EL": score_el,
            "Rank_EL": rank_el,
            "Score_BA": score_ba,
            "Rank_BA": rank_ba,
            "Aff_nM": aff,
        })
    return pd.DataFrame(rows)


t_pred = 0.0
alleles = sorted(df_score["hla_netmhc"].unique())
print(f"[predict] running per-allele on {len(alleles)} alleles ({len(df_score)} rows)...", flush=True)
for i, hla in enumerate(alleles):
    sub = df_score[df_score["hla_netmhc"] == hla]
    peps = sub["peptide"].tolist()
    pep_file = WORK / f"peps_{i:03d}.txt"
    pep_file.write_text("\n".join(peps) + "\n")

    t1 = time.time()
    proc = subprocess.run(
        [
            "docker", "run", "--rm",
            "-v", f"{WORK.resolve()}:/work",
            "-e", "NETMHCpan=/netMHCpan-4.1",
            "-e", "TMPDIR=/tmp",
            "--entrypoint", NETMHC_BIN,
            DOCKER_IMG,
            "-p", f"/work/{pep_file.name}",
            "-a", hla,
            "-BA",
        ],
        capture_output=True, text=True, timeout=600,
    )
    t_pred += time.time() - t1
    if proc.returncode != 0:
        print(f"[predict] {hla} FAILED rc={proc.returncode}; stderr head: {proc.stderr[:300]}", flush=True)
        continue

    parsed = parse_netmhcpan_output(proc.stdout)
    if parsed.empty:
        print(f"[predict] {hla} returned 0 parsed rows; raw head:\n{proc.stdout[:500]}", flush=True)
        continue
    # netMHCpan often outputs MHC like "HLA-A*02:01" (re-formatted with asterisk)
    parsed["pep_key"] = parsed["peptide"]
    # map back by peptide (peptides per allele are unique-or-deduped within sub; if not we average)
    score_map = parsed.groupby("pep_key").agg({
        "Rank_EL": "mean", "Rank_BA": "mean",
    })
    for idx in sub.index:
        pep = df_score.at[idx, "peptide"]
        if pep in score_map.index:
            df_score.at[idx, "pct_rank_el"] = score_map.at[pep, "Rank_EL"]
            df_score.at[idx, "pct_rank_ba"] = score_map.at[pep, "Rank_BA"]
    print(f"[predict] {hla} ({len(peps)} peps) ok in {time.time()-t1:.1f}s", flush=True)

# Score = -%Rank_EL (higher = better). Fall back to -%Rank_BA if EL missing.
df_score["score_netmhcpan"] = -df_score["pct_rank_el"]
mask_el_missing = df_score["pct_rank_el"].isna() & df_score["pct_rank_ba"].notna()
df_score.loc[mask_el_missing, "score_netmhcpan"] = -df_score.loc[mask_el_missing, "pct_rank_ba"]
n_no_score = int(df_score["score_netmhcpan"].isna().sum())
print(f"[score] {n_no_score} rows still missing a score", flush=True)

out_pred = df_score.rename(columns={"HLA_norm": "hla"})[
    ["peptide", "hla", "label", "in_master", "score_netmhcpan", "source", "split"]
].copy()
out_pred.to_csv(OUT / "predictions.tsv", sep="\t", index=False)
print(f"[write] {OUT/'predictions.tsv'} ({len(out_pred)} rows)", flush=True)


def auroc_with_ci(y, s, n_boot=1000, seed=0):
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    keep = ~np.isnan(s)
    y, s = y[keep], s[keep]
    if len(np.unique(y)) < 2 or len(y) < 2:
        return (float("nan"), float("nan"), float("nan"), int(len(y)))
    auc = roc_auc_score(y, s)
    rng = np.random.default_rng(seed)
    n = len(y)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yy = y[idx]
        ss = s[idx]
        if len(np.unique(yy)) < 2:
            continue
        aucs.append(roc_auc_score(yy, ss))
    if not aucs:
        return (float(auc), float("nan"), float("nan"), int(n))
    lo = float(np.percentile(aucs, 2.5))
    hi = float(np.percentile(aucs, 97.5))
    return (float(auc), lo, hi, int(n))


rows = []
its = out_pred.copy()

for testset, sub in [
    ("ITSNdb_combined", its),
    ("ITSNdb_no_overlap", its[its["in_master"] == False]),
    ("ITSNdb_in_master", its[its["in_master"] == True]),
]:
    y = sub["label"].astype(int)
    s = sub["score_netmhcpan"]
    auc, lo, hi, n_used = auroc_with_ci(y, s)
    rows.append(dict(
        testset=testset,
        in_master={"ITSNdb_combined": "any", "ITSNdb_no_overlap": "False", "ITSNdb_in_master": "True"}[testset],
        n=len(sub),
        n_pos=int(y.sum()),
        AUROC=auc,
        AUROC_lo95=lo,
        AUROC_hi95=hi,
    ))

auroc_df = pd.DataFrame(rows)
auroc_df.to_csv(OUT / "auroc_summary.tsv", sep="\t", index=False)
print(auroc_df.to_string(index=False), flush=True)

runtime_total = time.time() - t0
meta = dict(
    runtime_total_sec=runtime_total,
    runtime_predict_sec=t_pred,
    n_input_after_filter=len(df),
    n_scored=len(df_score),
    n_skipped_unsupported_hla=n_skipped_total,
    n_dropped_pep_len=n_dropped_len,
    n_dropped_non_std_aa=n_dropped_aa,
    n_no_score_after_predict=n_no_score,
    skipped_alleles=allele_skipped,
    install_path="docker image kevinr9525/netmhcpan-4.1b.linux:latest",
    netmhc_binary=NETMHC_BIN,
    n_alleles_run=len(alleles),
)
(OUT / "_run_meta.json").write_text(json.dumps(meta, indent=2))
print(f"[done] total {runtime_total:.1f}s", flush=True)
