"""Wave 5A — Step 2.5: Build ESM2-150M frozen embeddings cache for the
Wave 5A pool (distill_pool ∪ ITSNdb ∪ Wave1 train pool). Same architecture as
Wave 1's `embed_esm2.py` but extended pool.

Saves: wave5a_embeddings.pt
  pep_emb [N_pep, 640], pep_keys [N_pep]
  hla_emb [N_hla, 640], hla_keys [N_hla]
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch

torch.set_num_threads(8)

MODEL_NAME = "facebook/esm2_t30_150M_UR50D"
HIDDEN = 640
BATCH = 64

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
W5A = ROOT / "wave5a"
BUNDLE = ROOT / "bundle.tsv"
HLA_PSEUDO_TSV = ROOT / "hla_pseudo.tsv"
POOL_TSV = W5A / "distill_pool.tsv"
OUT = W5A / "wave5a_embeddings.pt"

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


def load_model():
    from transformers import AutoTokenizer, AutoModel
    tok = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME).to(DEVICE).eval()
    for p in model.parameters():
        p.requires_grad_(False)
    return tok, model


@torch.no_grad()
def embed_batch(tok, model, seqs):
    enc = tok(list(seqs), padding=True, return_tensors="pt", add_special_tokens=True)
    enc = {k: v.to(DEVICE) for k, v in enc.items()}
    out = model(**enc).last_hidden_state
    m = enc["attention_mask"].unsqueeze(-1).float()
    pooled = (out * m).sum(1) / m.sum(1).clamp(min=1)
    return pooled.float().cpu().numpy()


def chunk(xs, n):
    for i in range(0, len(xs), n):
        yield xs[i:i + n]


def main():
    t0 = time.time()
    print(f"Wave 5A — ESM2 embeddings on {DEVICE} (threads={torch.get_num_threads()})")
    tok, model = load_model()
    print(f"  loaded in {time.time()-t0:.1f}s")

    # Pool peptides
    pool = pd.read_csv(POOL_TSV, sep="\t")
    pool_peps = set(pool["peptide"].astype(str).tolist())

    # ITSNdb (eval set) peptides
    bundle = pd.read_csv(BUNDLE, sep="\t")
    eval_peps = set(bundle[bundle["split"].isin(
        ["ext_itsndb_main", "ext_itsndb_val"])]["peptide"].astype(str).tolist())

    # Wave 1 train peptides (already in pool but ensure)
    train_peps = set(bundle[bundle["split"] == "train"]["peptide"].astype(str).tolist())

    pep_set = pool_peps | eval_peps | train_peps
    pep_set = {p for p in pep_set
               if isinstance(p, str) and 7 < len(p) < 16 and
               all(c in "ACDEFGHIKLMNPQRSTVWY" for c in p)}
    pep_keys = sorted(pep_set)
    print(f"  unique peptides to embed: {len(pep_keys):,}")

    pep_emb = np.zeros((len(pep_keys), HIDDEN), dtype=np.float32)
    t1 = time.time()
    last_log = t1
    for i, batch in enumerate(chunk(pep_keys, BATCH)):
        e = embed_batch(tok, model, batch)
        pep_emb[i*BATCH:i*BATCH + len(batch)] = e
        if (i + 1) % 20 == 0 or time.time() - last_log > 30:
            elapsed = time.time() - t1
            done = (i + 1) * BATCH
            rate = done / max(elapsed, 1)
            eta = max(0, (len(pep_keys) - done) / max(rate, 1))
            print(f"  pep batch {i+1}: {done}/{len(pep_keys)}  rate={rate:.0f}/s  eta={eta:.0f}s",
                  flush=True)
            last_log = time.time()
    print(f"  peptide embedding done: {time.time()-t1:.1f}s")

    # HLA pseudo-seqs from existing table (Wave 1)
    pseu = pd.read_csv(HLA_PSEUDO_TSV, sep="\t")
    pseu = pseu[pseu["pseudo_seq"].astype(str).str.len() == 34].reset_index(drop=True)
    hla_keys = pseu["HLA_norm"].tolist()
    seqs = pseu["pseudo_seq"].tolist()
    hla_emb = np.zeros((len(hla_keys), HIDDEN), dtype=np.float32)
    t2 = time.time()
    for i, batch in enumerate(chunk(seqs, BATCH)):
        e = embed_batch(tok, model, batch)
        hla_emb[i*BATCH:i*BATCH + len(batch)] = e
    print(f"  HLA embedding done: {time.time()-t2:.1f}s  ({len(hla_keys)} alleles)")

    payload = {
        "pep_emb": torch.from_numpy(pep_emb),
        "pep_keys": pep_keys,
        "hla_emb": torch.from_numpy(hla_emb),
        "hla_keys": hla_keys,
        "model_name": MODEL_NAME,
        "hidden": HIDDEN,
    }
    torch.save(payload, OUT)
    print(f"\nsaved: {OUT}  ({OUT.stat().st_size/1e6:.1f} MB)  total {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
