#!/usr/bin/env python
"""Score NetMHCstabpan-1.0 (DTU; pMHC complex stability prediction) on leakage-stratified ITSNdb.

Container: kevinr9525/netmhcstabpan-1.0a.linux:latest
Score: -%Rank_Stab (lower rank = more stable; sign-flipped so higher = better)

Output: predictions_netmhcstabpan.tsv, auroc_netmhcstabpan.tsv
"""
from __future__ import annotations
import os, sys, time, json, subprocess, re, tempfile
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
OUT = ROOT / "wave9"
OUT.mkdir(parents=True, exist_ok=True)
BUNDLE = ROOT / "bundle.tsv"
DOCKER_IMG = "kevinr9525/netmhcstabpan-1.0a.linux:latest"
BIN = "/netMHCstabpan-1.0/Linux_x86_64/bin/netMHCstabpan"

t0 = time.time()
print("[stabpan] load bundle", flush=True)
df = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master": "boolean"})
df = df[df["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].reset_index(drop=True)
print(f"[stabpan] kept {len(df)} ITSNdb rows", flush=True)

AA = set("ACDEFGHIKLMNPQRSTVWY")
df = df[df["peptide"].apply(lambda s: set(s).issubset(AA))].reset_index(drop=True)
df = df[df["peptide"].str.len().between(8, 11)].reset_index(drop=True)
print(f"[stabpan] after AA/length filter: {len(df)}", flush=True)


def to_stabpan(hla: str) -> str:
    """HLA-A*02:01 -> HLA-A02:01"""
    return hla.replace("*", "")


df["hla_stab"] = df["HLA_norm"].map(to_stabpan)

# Get supported alleles
print("[stabpan] fetching supported allele list", flush=True)
r = subprocess.run(["docker", "run", "--rm", DOCKER_IMG, BIN, "-listMHC"],
                   capture_output=True, text=True)
supported = set(line.strip() for line in r.stdout.splitlines() if line.strip().startswith("HLA-"))
print(f"[stabpan] {len(supported)} alleles supported", flush=True)

df["hla_ok"] = df["hla_stab"].isin(supported)
n_skip = int((~df["hla_ok"]).sum())
skipped_alleles = df[~df["hla_ok"]]["HLA_norm"].value_counts().to_dict()
print(f"[stabpan] {n_skip} skipped (unsupported HLA): {skipped_alleles}", flush=True)
df_score = df[df["hla_ok"]].copy().reset_index(drop=True)
print(f"[stabpan] scoring {len(df_score)} rows over {df_score['hla_stab'].nunique()} alleles", flush=True)


def parse_stabpan_output(stdout: str) -> list[dict]:
    """Parse netMHCstabpan output. Format:
    pos  HLA  peptide  Identity  Pred  Thalf(h)  %Rank_Stab  BindLevel
    """
    rows = []
    for line in stdout.splitlines():
        if not re.match(r"^\s*\d+\s+HLA-", line):
            continue
        parts = line.split()
        if len(parts) < 7:
            continue
        try:
            hla = parts[1]  # HLA-A*02:01
            pep = parts[2]
            pred = float(parts[4])
            thalf = float(parts[5])
            rank_stab = float(parts[6])
        except (ValueError, IndexError):
            continue
        rows.append({"hla_full": hla, "peptide": pep,
                     "stab_pred": pred, "stab_thalf_h": thalf, "stab_rank": rank_stab})
    return rows


# Run per allele × length
all_scores = []
errors = []
unique_keys = df_score.groupby("hla_stab")
for hla_stab, sub in unique_keys:
    pep_lens = sorted(sub["peptide"].str.len().unique())
    peps = sub["peptide"].unique().tolist()
    with tempfile.TemporaryDirectory() as td:
        pep_path = Path(td) / "peps.txt"
        pep_path.write_text("\n".join(peps) + "\n")
        for L in pep_lens:
            sub_L_peps = [p for p in peps if len(p) == L]
            if not sub_L_peps:
                continue
            pep_L = Path(td) / f"peps_L{L}.txt"
            pep_L.write_text("\n".join(sub_L_peps) + "\n")
            cmd = ["docker", "run", "--rm", "-v", f"{td}:/work", DOCKER_IMG,
                   BIN, "-p", "-f", f"/work/{pep_L.name}", "-a", hla_stab, "-l", str(L)]
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if r.returncode != 0:
                errors.append({"hla": hla_stab, "L": L, "stderr": r.stderr[:500]})
                print(f"[stabpan] FAIL {hla_stab} L={L}: rc={r.returncode}", flush=True)
                continue
            parsed = parse_stabpan_output(r.stdout)
            all_scores.extend(parsed)
    print(f"[stabpan] {hla_stab}: {len([s for s in all_scores if s['hla_full'].replace('*','')==hla_stab])} scored", flush=True)

sc = pd.DataFrame(all_scores)
print(f"[stabpan] total parsed {len(sc)} rows", flush=True)
sc["hla_norm_match"] = sc["hla_full"]  # already HLA-A*02:01 form

# Merge: align by peptide + HLA_norm
df_score = df_score.merge(
    sc[["hla_norm_match", "peptide", "stab_pred", "stab_thalf_h", "stab_rank"]],
    left_on=["HLA_norm", "peptide"], right_on=["hla_norm_match", "peptide"],
    how="left"
).drop(columns=["hla_norm_match"])

# Score: lower rank_stab = more stable, flip
df_score["score_stabpan"] = -df_score["stab_rank"]
n_no_score = int(df_score["stab_rank"].isna().sum())
print(f"[stabpan] {n_no_score} rows lost score after merge", flush=True)

# Predictions output
out_pred = df_score.rename(columns={"HLA_norm": "hla"})[
    ["peptide", "hla", "label", "in_master", "split", "source",
     "stab_pred", "stab_thalf_h", "stab_rank", "score_stabpan"]
].copy()
out_pred.to_csv(OUT / "predictions_netmhcstabpan.tsv", sep="\t", index=False)
print(f"[stabpan] wrote {OUT/'predictions_netmhcstabpan.tsv'} ({len(out_pred)} rows)", flush=True)


# AUROC
from sklearn.metrics import roc_auc_score


def auroc_ci(y, s, n_boot=1000, seed=0):
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    if len(np.unique(y)) < 2:
        return (np.nan, np.nan, np.nan)
    auc = roc_auc_score(y, s)
    rng = np.random.default_rng(seed)
    n = len(y)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yy, ss = y[idx], s[idx]
        if len(np.unique(yy)) < 2:
            continue
        aucs.append(roc_auc_score(yy, ss))
    if not aucs:
        return (auc, np.nan, np.nan)
    return (float(auc), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5)))


rows = []
ok = out_pred[out_pred["score_stabpan"].notna()].copy()
for label, mask in [
    ("ITSNdb_combined", ok["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])),
    ("ITSNdb_no_overlap", (ok["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])) & (ok["in_master"] == False)),
    ("ITSNdb_in_master", (ok["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])) & (ok["in_master"] == True)),
]:
    sub = ok[mask]
    y = sub["label"].astype(int).to_numpy()
    s = sub["score_stabpan"].astype(float).to_numpy()
    auc, lo, hi = auroc_ci(y, s, n_boot=1000, seed=42)
    rows.append({"algorithm": "NetMHCstabpan", "testset": label,
                 "n": len(sub), "n_pos": int((sub["label"] == 1).sum()),
                 "n_neg": int((sub["label"] == 0).sum()),
                 "AUROC": auc, "AUROC_lo95": lo, "AUROC_hi95": hi,
                 "n_skipped_unsupported_hla": n_skip,
                 "score_col": "neg_rank_stab"})

ar = pd.DataFrame(rows)
ar.to_csv(OUT / "auroc_netmhcstabpan.tsv", sep="\t", index=False)
print(ar.to_string(index=False), flush=True)

(OUT / "_meta_netmhcstabpan.json").write_text(json.dumps({
    "n_input": int(len(df)),
    "n_scored": int(len(ok)),
    "n_skipped_unsupported_hla": n_skip,
    "skipped_alleles": skipped_alleles,
    "n_lost_after_merge": n_no_score,
    "n_alleles_run": int(df_score["hla_stab"].nunique()),
    "errors": errors,
    "runtime_sec": time.time() - t0,
}, indent=2))
print(f"[stabpan] total {time.time()-t0:.1f}s", flush=True)
