#!/usr/bin/env python3
"""Wave 11 — score TransPHLA-AOMP on every test bundle."""
from __future__ import annotations
import os, sys, time, json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
WAVE11 = ROOT / "wave11"
BUN = WAVE11 / "test_bundles"
OUT = WAVE11 / "wave11_predictions"

TRANSPHLA_DIR = Path("/tmp/TransPHLA/TransPHLA-AOMP")
sys.path.insert(0, str(TRANSPHLA_DIR))
import scipy as _scipy
if not hasattr(_scipy, "interp"):
    _scipy.interp = np.interp

import torch
import torch.nn as _nn
_orig_cwd = os.getcwd()
os.chdir(TRANSPHLA_DIR)
from model import Transformer, read_predict_data
os.chdir(_orig_cwd)

def eval_step_safe(model, val_loader):
    model.eval()
    torch.manual_seed(19961231)
    y_prob_list = []
    with torch.no_grad():
        for pep_in, hla_in in val_loader:
            outputs, _, _, _ = model(pep_in, hla_in)
            y_prob = _nn.Softmax(dim=1)(outputs)[:, 1].cpu().detach().numpy()
            y_prob_list.extend(y_prob.tolist())
    return y_prob_list

print("[load] HLA db + model...", flush=True)
hla_db = pd.read_csv(TRANSPHLA_DIR / "common_hla_sequence.csv")
supported = set(hla_db["HLA"])

device = torch.device("cpu")
model = Transformer().to(device)
state = torch.load(str(TRANSPHLA_DIR / "pHLAIformer.pkl"), map_location="cpu")
model.load_state_dict(state, strict=True)
model.eval()
print(f"[load] model loaded; {len(supported)} HLA supported", flush=True)

AA = set("ACDEFGHIKLMNPQRSTVWY")
t0 = time.time()
for bp_path in sorted(BUN.glob("*.tsv")):
    bname = bp_path.stem
    df = pd.read_csv(bp_path, sep="\t")
    if len(df) == 0: continue
    df = df.copy()
    df["plen"] = df["peptide"].str.len()
    df = df[df["plen"].between(8, 15) & df["hla"].notna()
            & df["peptide"].apply(lambda s: set(s).issubset(AA))]
    df = df[df["hla"].isin(supported)]
    if len(df) == 0:
        print(f"[skip] {bname}: 0 eligible"); continue
    print(f"[transphla] {bname}: scoring n={len(df)}...", flush=True)
    t1 = time.time()
    input_df = pd.DataFrame({"HLA": df["hla"].values, "peptide": df["peptide"].values})
    os.chdir(TRANSPHLA_DIR)
    predict_data, pep_inputs, hla_inputs, predict_loader = read_predict_data(input_df, batch_size=1024)
    y_prob = eval_step_safe(model, predict_loader)
    os.chdir(_orig_cwd)
    predict_data = predict_data.copy()
    predict_data["score_transphla"] = y_prob
    score_map = predict_data.set_index(["HLA", "peptide"])["score_transphla"].to_dict()
    df["score_transphla"] = df.apply(
        lambda r: score_map.get((r["hla"], r["peptide"]), float("nan")), axis=1
    )
    out = df[["peptide", "hla", "label", "source", "length",
              "n_overlap_with_other_sources", "score_transphla"]].copy()
    out.to_csv(OUT / f"TransPHLA__{bname}.tsv", sep="\t", index=False)
    print(f"  done in {time.time()-t1:.1f}s; {len(out)} preds", flush=True)

print(f"[done] total {time.time()-t0:.1f}s")
