#!/usr/bin/env python
"""Score MHCnuggets-2 (Karchin Lab, Shao 2020) on leakage-stratified ITSNdb.

Output: predictions_mhcnuggets.tsv, auroc_mhcnuggets.tsv
"""
from __future__ import annotations
import os, sys, time, json, tempfile
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
OUT = ROOT / "wave9"
OUT.mkdir(parents=True, exist_ok=True)
BUNDLE = ROOT / "bundle.tsv"

t0 = time.time()
print("[mhcnuggets] load bundle", flush=True)
df = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master": "boolean"})
df = df[df["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].reset_index(drop=True)
print(f"[mhcnuggets] kept {len(df)} ITSNdb rows", flush=True)

AA = set("ACDEFGHIKLMNPQRSTVWY")
df = df[df["peptide"].apply(lambda s: set(s).issubset(AA))].reset_index(drop=True)
df = df[df["peptide"].str.len().between(8, 11)].reset_index(drop=True)
print(f"[mhcnuggets] after AA/length filter: {len(df)}", flush=True)

# Load supported alleles list (MHCnuggets format: HLA-A02:01, no asterisk)
sup_path = "/tmp/mhcnuggets/mhcnuggets/data/production/mhcI/alleles_with_trained_models.txt"
supported = set(open(sup_path).read().splitlines())
print(f"[mhcnuggets] {len(supported)} alleles supported", flush=True)


def to_mhcnuggets(hla: str) -> str:
    """HLA-A*02:01 -> HLA-A02:01"""
    return hla.replace("*", "")


df["mhc"] = df["HLA_norm"].map(to_mhcnuggets)
df["hla_ok"] = df["mhc"].isin(supported)
n_skip_hla = int((~df["hla_ok"]).sum())
skipped_alleles = df[~df["hla_ok"]]["HLA_norm"].value_counts().to_dict()
print(f"[mhcnuggets] {n_skip_hla} rows skipped (unsupported HLA): {skipped_alleles}", flush=True)

df_score = df[df["hla_ok"]].copy().reset_index(drop=True)
print(f"[mhcnuggets] scoring {len(df_score)} rows over {df_score['mhc'].nunique()} alleles",
      flush=True)

# Run mhcnuggets predict per allele (writes peptides.txt -> predictions.csv)
from mhcnuggets.src.predict import predict as mhc_predict

scores = []  # list of (peptide, mhc, ic50, score)
for mhc, sub in df_score.groupby("mhc"):
    peptides = sub["peptide"].tolist()
    with tempfile.TemporaryDirectory() as td:
        pep_path = Path(td) / "peptides.txt"
        out_path = Path(td) / "out.csv"
        pep_path.write_text("\n".join(peptides) + "\n")
        try:
            mhc_predict(class_="I", peptides_path=str(pep_path), mhc=mhc,
                        output=str(out_path), rank_output=True)
        except Exception as e:
            print(f"[mhcnuggets] FAILED {mhc}: {e}", flush=True)
            continue
        if not out_path.exists():
            print(f"[mhcnuggets] no output for {mhc}", flush=True)
            continue
        pr = pd.read_csv(out_path)
        # Output columns: peptide, ic50, human_proteome_rank (when rank_output)
        for _, r in pr.iterrows():
            scores.append({
                "peptide": r["peptide"],
                "mhc": mhc,
                "mhcnuggets_ic50": float(r["ic50"]),
                "mhcnuggets_rank": float(r.get("human_proteome_rank", np.nan)) if "human_proteome_rank" in pr.columns else np.nan,
            })
    print(f"[mhcnuggets] done {mhc}: {len(sub)} peptides", flush=True)

sc = pd.DataFrame(scores)
df_score = df_score.merge(sc, on=["peptide", "mhc"], how="left")

# Build predictions output. Score for AUROC: -log10(IC50) so higher = stronger binder
# Equivalent: -ic50 (negate so higher = better). MHCnuggets is binding affinity not immunogenicity
# but it's a recent (2024 maintenance) update of a classic so it's a reasonable comparator.
df_score["score_mhcnuggets"] = -df_score["mhcnuggets_ic50"]

out_pred = df_score.rename(columns={"HLA_norm": "hla"})[
    ["peptide", "hla", "label", "in_master", "split", "source",
     "mhcnuggets_ic50", "mhcnuggets_rank", "score_mhcnuggets"]
].copy()
out_pred.to_csv(OUT / "predictions_mhcnuggets.tsv", sep="\t", index=False)
print(f"[mhcnuggets] wrote {OUT/'predictions_mhcnuggets.tsv'} ({len(out_pred)} rows)", flush=True)


# --- AUROC w/ bootstrap ---
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
ok = out_pred[out_pred["score_mhcnuggets"].notna()].copy()
n_skipped_total = int(out_pred["score_mhcnuggets"].isna().sum() + n_skip_hla)

for label, mask in [
    ("ITSNdb_combined", ok["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])),
    ("ITSNdb_no_overlap", (ok["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])) & (ok["in_master"] == False)),
    ("ITSNdb_in_master", (ok["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])) & (ok["in_master"] == True)),
]:
    sub = ok[mask]
    y = sub["label"].astype(int).to_numpy()
    s = sub["score_mhcnuggets"].astype(float).to_numpy()
    auc, lo, hi = auroc_ci(y, s, n_boot=1000, seed=42)
    rows.append({"algorithm": "MHCnuggets", "testset": label,
                 "n": len(sub), "n_pos": int((sub["label"] == 1).sum()),
                 "n_neg": int((sub["label"] == 0).sum()),
                 "AUROC": auc, "AUROC_lo95": lo, "AUROC_hi95": hi,
                 "n_skipped_unsupported_hla": n_skip_hla,
                 "score_col": "neg_ic50"})

ar = pd.DataFrame(rows)
ar.to_csv(OUT / "auroc_mhcnuggets.tsv", sep="\t", index=False)
print(ar.to_string(index=False), flush=True)

(OUT / "_meta_mhcnuggets.json").write_text(json.dumps({
    "n_input": int(len(df)),
    "n_scored": int(len(out_pred[out_pred["score_mhcnuggets"].notna()])),
    "n_skipped_unsupported_hla": n_skip_hla,
    "skipped_alleles": skipped_alleles,
    "n_unique_alleles_scored": int(df_score["mhc"].nunique()),
    "runtime_sec": time.time() - t0,
}, indent=2))
print(f"[mhcnuggets] total {time.time()-t0:.1f}s", flush=True)
