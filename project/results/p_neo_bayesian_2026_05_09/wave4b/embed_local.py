"""Compute ESM2-150M peptide embeddings locally on CPU.

Wave 4B prerequisite — recreates a minimal embeddings.pt covering ALL peptides
in bundle.tsv (train + ITSNdb + Venus already enumerated as windows).

Key differences vs the original Wave 1 embed_esm2.py:
  - CPU only (Wave 4A is on the pod)
  - HF cache redirected to /data/thca/_hf_cache (root disk near full)
  - mean-pool over attention mask (matches Wave 1 convention)
  - skip Venus window enumeration: Wave 4B targets ITSNdb only

Output: wave4b/embeddings_local.pt
"""
from __future__ import annotations
import os, sys, time
from pathlib import Path

os.environ.setdefault("HF_HOME", "/data/thca/_hf_cache")
os.environ.setdefault("TRANSFORMERS_CACHE", "/data/thca/_hf_cache/transformers")
os.environ.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")

import numpy as np
import pandas as pd
import torch

MODEL_NAME = "facebook/esm2_t30_150M_UR50D"
HIDDEN = 640
BATCH = 8  # CPU-friendly
ROOT = Path(__file__).resolve().parent
WAVE1 = ROOT.parent
BUNDLE = WAVE1 / "bundle.tsv"
OUT = ROOT / "embeddings_local.pt"


def load_model():
    from transformers import AutoTokenizer, AutoModel
    tok = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME)
    model.eval()
    for p in model.parameters():
        p.requires_grad_(False)
    return tok, model


@torch.no_grad()
def embed_batch(tok, model, seqs):
    enc = tok(list(seqs), padding=True, return_tensors="pt", add_special_tokens=True)
    out = model(**enc).last_hidden_state
    mask = enc["attention_mask"].unsqueeze(-1).float()
    pooled = (out * mask).sum(1) / mask.sum(1).clamp(min=1)
    return pooled.float().numpy()


def chunk(xs, n):
    for i in range(0, len(xs), n):
        yield xs[i:i + n]


def main():
    t0 = time.time()
    print(f"loading {MODEL_NAME} on CPU...")
    tok, model = load_model()
    print(f"  loaded in {time.time()-t0:.1f}s")

    bundle = pd.read_csv(BUNDLE, sep="\t")
    print(f"bundle rows: {len(bundle)}")

    # Skip Venus rows; we only need train + ITSNdb peptides
    keep_splits = {"train", "ext_itsndb_main", "ext_itsndb_val"}
    sub = bundle[bundle["split"].isin(keep_splits)].copy()
    pep_set = set()
    for p in sub["peptide"]:
        if isinstance(p, str) and 7 < len(p) <= 16 and all(c in "ACDEFGHIKLMNPQRSTVWY" for c in p):
            pep_set.add(p)
    pep_keys = sorted(pep_set)
    print(f"unique peptides to embed: {len(pep_keys)}")

    pep_emb = np.zeros((len(pep_keys), HIDDEN), dtype=np.float32)
    t1 = time.time()
    n_b = (len(pep_keys) + BATCH - 1) // BATCH
    for i, batch in enumerate(chunk(pep_keys, BATCH)):
        e = embed_batch(tok, model, batch)
        pep_emb[i * BATCH:i * BATCH + len(batch)] = e
        if (i + 1) % 25 == 0 or (i + 1) == n_b:
            elapsed = time.time() - t1
            done = i * BATCH + len(batch)
            rate = done / max(elapsed, 1e-9)
            print(f"  batch {i+1}/{n_b}: {done}/{len(pep_keys)} rate={rate:.1f}/s eta={(len(pep_keys)-done)/max(rate,1e-9):.0f}s")

    print(f"\npep embed: {time.time()-t1:.1f}s")

    payload = {
        "pep_emb": torch.from_numpy(pep_emb),
        "pep_keys": pep_keys,
        "model_name": MODEL_NAME,
        "hidden": HIDDEN,
    }
    torch.save(payload, OUT)
    print(f"saved -> {OUT}  ({OUT.stat().st_size/1e6:.1f} MB)")
    print(f"total: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
