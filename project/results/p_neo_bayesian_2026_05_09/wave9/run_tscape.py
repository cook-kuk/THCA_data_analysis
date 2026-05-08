#!/usr/bin/env python
"""Score T-SCAPE / TITANiAN (Kim et al 2025 bioRxiv) pmhc_im_neo on leakage-stratified ITSNdb.

Repo: https://github.com/seoklab/TITANiAN
Model: best_param/pmhc_im_neo/BigMHC_finalMedium_OAS_el-mlm_ADV1.0_bestvalloss.pt (cancer-immunotherapy variant)

Output: predictions_tscape.tsv, auroc_tscape.tsv
"""
from __future__ import annotations
import os, sys, time, json, shutil, subprocess
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
OUT = ROOT / "wave9"
OUT.mkdir(parents=True, exist_ok=True)
BUNDLE = ROOT / "bundle.tsv"
TSCAPE_DIR = Path("/tmp/TITANiAN")

# Stage HF-downloaded checkpoint into TITANiAN dir under best_param/pmhc_im_neo/
CKPT_NAME = "BigMHC_finalMedium_OAS_el-mlm_ADV1.0_bestvalloss.pt"
HF_CKPT = Path("/data/thca/_tmp/tscape/models--seoklab--T-SCAPE/snapshots/314dcb5d1ff0a081c6f50e4f4857918cf6d18586/best_param/pmhc_im_neo") / CKPT_NAME
target_dir = TSCAPE_DIR / "best_param" / "pmhc_im_neo"
target_dir.mkdir(parents=True, exist_ok=True)
target_ckpt = target_dir / CKPT_NAME
if not target_ckpt.exists():
    if target_ckpt.is_symlink():
        target_ckpt.unlink()
    target_ckpt.symlink_to(HF_CKPT)
print(f"[tscape] checkpoint at {target_ckpt} -> {target_ckpt.resolve()}", flush=True)

t0 = time.time()
print("[tscape] load bundle", flush=True)
df = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master": "boolean"})
df = df[df["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])].reset_index(drop=True)
print(f"[tscape] kept {len(df)} ITSNdb rows", flush=True)

AA = set("ACDEFGHIKLMNPQRSTVWY")
df = df[df["peptide"].apply(lambda s: set(s).issubset(AA))].reset_index(drop=True)
df = df[df["peptide"].str.len().between(8, 11)].reset_index(drop=True)
print(f"[tscape] after AA/length filter: {len(df)}", flush=True)

# Build T-SCAPE input csv: needs columns "Allele","peptide". Allele format: HLA-A*02:01 acceptable
# (mhc_pseudo_matching.modify_entry_2 strips HLA, * and :)
inp_csv = OUT / "_tscape_input.csv"
mod_csv = OUT / "_tscape_input_modified.csv"
out_csv = OUT / "_tscape_output.csv"
df_in = df[["HLA_norm", "peptide"]].rename(columns={"HLA_norm": "Allele"})
df_in.to_csv(inp_csv, index=False)
print(f"[tscape] wrote input csv: {inp_csv} ({len(df_in)} rows)", flush=True)

# Run pseudo matching to expand allele -> mhc pseudo-seq
cmd1 = [sys.executable, str(TSCAPE_DIR / "mhc_pseudo_matching.py"), "I", str(inp_csv), str(mod_csv)]
print(f"[tscape] running: {' '.join(cmd1)}", flush=True)
r1 = subprocess.run(cmd1, cwd=str(TSCAPE_DIR), capture_output=True, text=True)
print(r1.stdout, flush=True)
if r1.returncode != 0:
    print(f"[tscape] mhc_pseudo_matching FAILED: {r1.stderr}", flush=True)
    sys.exit(1)
print(f"[tscape] mhc_pseudo_matching stderr: {r1.stderr[:500]}", flush=True)

# Inspect modified csv counts
mod = pd.read_csv(mod_csv)
print(f"[tscape] modified csv shape: {mod.shape}; cols: {list(mod.columns)}", flush=True)

# Run inference
cmd2 = [sys.executable, str(TSCAPE_DIR / "inference_csv.py"),
        "--csv_path", str(mod_csv),
        "--inf_type", "pmhc_im_neo",
        "--output", str(out_csv)]
print(f"[tscape] running: {' '.join(cmd2)}", flush=True)
t1 = time.time()
r2 = subprocess.run(cmd2, cwd=str(TSCAPE_DIR), capture_output=True, text=True)
print(r2.stdout[-2000:], flush=True)
if r2.returncode != 0:
    print(f"[tscape] inference FAILED: {r2.stderr[-2000:]}", flush=True)
    sys.exit(1)
t_inf = time.time() - t1
print(f"[tscape] inference done in {t_inf:.1f}s", flush=True)

# Read T-SCAPE output and merge into bundle
ts_out = pd.read_csv(out_csv)
print(f"[tscape] output shape: {ts_out.shape}; cols: {list(ts_out.columns)}", flush=True)
# Output should have Allele, peptide, score
# Merge by peptide+Allele back into our df
# T-SCAPE strips Allele formatting; mod was already merged so we use the MODIFIED Allele from mod_csv to align
# Easier: ts_out should preserve the Allele column from mod_csv (which is the modified short form)
df["Allele_short"] = mod["Allele"][:len(df)].values  # mod row order matches df row order (mhc_pseudo_matching kept order)

# Robust merge: use index alignment if length matches, else by (Allele_short, peptide)
if len(ts_out) == len(df):
    df["score_tscape"] = ts_out["score"].values
else:
    # merge by Allele_short + peptide
    print(f"[tscape] WARNING length mismatch {len(ts_out)} vs {len(df)}, falling back to merge", flush=True)
    df = df.merge(ts_out[["Allele", "peptide", "score"]],
                  left_on=["Allele_short", "peptide"], right_on=["Allele", "peptide"],
                  how="left").rename(columns={"score": "score_tscape"})
n_skip = int(df["score_tscape"].isna().sum())
print(f"[tscape] {n_skip} rows missing score post-merge", flush=True)

# Save predictions
out_pred = df.rename(columns={"HLA_norm": "hla"})[
    ["peptide", "hla", "label", "in_master", "split", "source", "score_tscape"]
].copy()
out_pred.to_csv(OUT / "predictions_tscape.tsv", sep="\t", index=False)
print(f"[tscape] wrote {OUT/'predictions_tscape.tsv'} ({len(out_pred)} rows)", flush=True)

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
ok = out_pred[out_pred["score_tscape"].notna()].copy()
for label, mask in [
    ("ITSNdb_combined", ok["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])),
    ("ITSNdb_no_overlap", (ok["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])) & (ok["in_master"] == False)),
    ("ITSNdb_in_master", (ok["split"].isin(["ext_itsndb_main", "ext_itsndb_val"])) & (ok["in_master"] == True)),
]:
    sub = ok[mask]
    y = sub["label"].astype(int).to_numpy()
    s = sub["score_tscape"].astype(float).to_numpy()
    auc, lo, hi = auroc_ci(y, s, n_boot=1000, seed=42)
    rows.append({"algorithm": "T-SCAPE", "testset": label,
                 "n": len(sub), "n_pos": int((sub["label"] == 1).sum()),
                 "n_neg": int((sub["label"] == 0).sum()),
                 "AUROC": auc, "AUROC_lo95": lo, "AUROC_hi95": hi,
                 "n_skipped_unsupported_hla": n_skip,
                 "score_col": "tscape_pmhc_im_neo"})

ar = pd.DataFrame(rows)
ar.to_csv(OUT / "auroc_tscape.tsv", sep="\t", index=False)
print(ar.to_string(index=False), flush=True)

(OUT / "_meta_tscape.json").write_text(json.dumps({
    "n_input": int(len(df)),
    "n_scored": int(len(ok)),
    "n_skipped": int(n_skip),
    "model_variant": "pmhc_im_neo (cancer-immunotherapy)",
    "ckpt": CKPT_NAME,
    "runtime_inference_sec": t_inf,
    "runtime_total_sec": time.time() - t0,
}, indent=2))
print(f"[tscape] total {time.time()-t0:.1f}s", flush=True)
