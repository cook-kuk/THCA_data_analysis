#!/usr/bin/env python
"""Wave 5B — frozen ESM2-150M (640d) mean-pool peptide + HLA pseudo embeddings.

Only embeds the 8-11mer peptides + 34-aa HLA pseudo sequences that appear in
multitask_supervision.tsv. CPU-only. Saves to wave5b/embeddings.pt.
"""
from __future__ import annotations
import os, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch

# Route HF cache to /data (root disk is at 84%)
os.environ.setdefault("HF_HOME", "/data/thca/_repo_offload/hf_cache")
os.environ.setdefault("HF_HUB_CACHE", "/data/thca/_repo_offload/hf_cache/hub")

MODEL_NAME = "facebook/esm2_t30_150M_UR50D"
HIDDEN = 640
BATCH = 32
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
W5B = ROOT / "wave5b"
SUP = W5B / "multitask_supervision.tsv"
PSEU = ROOT / "hla_pseudo.tsv"
OUT = W5B / "embeddings.pt"


def chunk(xs, n):
    for i in range(0, len(xs), n):
        yield xs[i:i + n]


@torch.no_grad()
def embed_batch(tok, model, seqs):
    enc = tok(list(seqs), padding=True, return_tensors="pt", add_special_tokens=True)
    enc = {k: v.to(DEVICE) for k, v in enc.items()}
    out = model(**enc).last_hidden_state
    mask = enc["attention_mask"].unsqueeze(-1).float()
    pooled = (out * mask).sum(1) / mask.sum(1).clamp(min=1)
    return pooled.float().cpu().numpy()


def main():
    t0 = time.time()
    if OUT.exists():
        print(f"[skip] {OUT} exists ({OUT.stat().st_size/1e6:.1f} MB)")
        return

    print(f"[load] device={DEVICE} model={MODEL_NAME}")
    from transformers import AutoTokenizer, AutoModel
    tok = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME).to(DEVICE)
    model.eval()
    for p in model.parameters():
        p.requires_grad_(False)
    print(f"[load] done {time.time()-t0:.1f}s")

    sup = pd.read_csv(SUP, sep="\t")
    pep_keys = sorted(sup["peptide"].unique().tolist())
    hla_keys_in_sup = sorted(sup["hla"].unique().tolist())
    print(f"[data] peptides={len(pep_keys)}  hlas_in_sup={len(hla_keys_in_sup)}")

    pseu = pd.read_csv(PSEU, sep="\t")
    pseu = pseu[pseu["pseudo_seq"].astype(str).str.len() == 34].reset_index(drop=True)
    pseu_map = dict(zip(pseu["HLA_norm"], pseu["pseudo_seq"]))
    have_pseu = [h for h in hla_keys_in_sup if h in pseu_map]
    miss_pseu = [h for h in hla_keys_in_sup if h not in pseu_map]
    print(f"[data] hla with pseudo: {len(have_pseu)}; missing: {len(miss_pseu)}")

    # ---- peptide embed ----
    pep_emb = np.zeros((len(pep_keys), HIDDEN), dtype=np.float32)
    t1 = time.time()
    for i, batch in enumerate(chunk(pep_keys, BATCH)):
        e = embed_batch(tok, model, batch)
        pep_emb[i*BATCH:i*BATCH + len(batch)] = e
        if (i + 1) % 25 == 0:
            done = (i + 1) * BATCH
            rate = done / (time.time() - t1)
            print(f"  pep {min(done, len(pep_keys))}/{len(pep_keys)} rate={rate:.0f}/s")
    print(f"[pep] done {time.time()-t1:.1f}s")

    # ---- hla pseudo embed ----
    hla_emb = np.zeros((len(have_pseu), HIDDEN), dtype=np.float32)
    seqs = [pseu_map[h] for h in have_pseu]
    t2 = time.time()
    for i, batch in enumerate(chunk(seqs, BATCH)):
        e = embed_batch(tok, model, batch)
        hla_emb[i*BATCH:i*BATCH + len(batch)] = e
    print(f"[hla] done {time.time()-t2:.1f}s ({len(have_pseu)} alleles)")

    payload = {
        "pep_emb": torch.from_numpy(pep_emb),
        "pep_keys": pep_keys,
        "hla_emb": torch.from_numpy(hla_emb),
        "hla_keys": have_pseu,
        "model_name": MODEL_NAME,
        "hidden": HIDDEN,
    }
    torch.save(payload, OUT)
    print(f"[save] {OUT}  ({OUT.stat().st_size/1e6:.1f} MB)")
    print(f"[done] total {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
