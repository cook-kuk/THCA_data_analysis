"""Compute frozen ESM2-150M (facebook/esm2_t30_150M_UR50D, hidden=640) mean-pool
embeddings for all unique peptides + HLA pseudo-sequences in bundle.tsv.

Saves: embeddings.pt (dict)
  {"pep_emb": [N_pep, 640], "pep_keys": [N_pep],
   "hla_emb": [N_hla, 640], "hla_keys": [N_hla],
   "venus_windows": dict(protein_id -> [n_w, 640]) }

We enumerate 8-11mer windows for VenusVaccine proteins inline and embed each
window so the eval script can score (window × HLA) pairs.

Run on pod: A6000 ~3-5 min for ~30k sequences batch=64.
"""
from __future__ import annotations
import sys, time, json
from pathlib import Path
import numpy as np
import pandas as pd
import torch

MODEL_NAME = "facebook/esm2_t30_150M_UR50D"
HIDDEN = 640
BATCH = 64
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

ROOT = Path(__file__).parent
BUNDLE = ROOT / "bundle.tsv"
HLA_PSEUDO_TSV = ROOT / "hla_pseudo.tsv"
OUT = ROOT / "embeddings.pt"


def load_model():
    from transformers import AutoTokenizer, AutoModel
    tok = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME).to(DEVICE)
    model.eval()
    for p in model.parameters():
        p.requires_grad_(False)
    return tok, model


@torch.no_grad()
def embed_batch(tok, model, seqs):
    """Tokenize → forward → masked mean-pool over residues. Returns [B, 640]."""
    enc = tok(list(seqs), padding=True, return_tensors="pt", add_special_tokens=True)
    enc = {k: v.to(DEVICE) for k, v in enc.items()}
    out = model(**enc).last_hidden_state  # [B, L, H]
    mask = enc["attention_mask"].unsqueeze(-1).float()
    # Drop CLS/EOS tokens from pool: facebook/esm2 uses cls=0, eos=2 — easiest
    # is the attention mask which already excludes pads; we keep CLS/EOS in
    # for simplicity (this is what fair-esm + transformers people do, mean
    # over the full attended span). Effect on a 9-11mer is small.
    pooled = (out * mask).sum(1) / mask.sum(1).clamp(min=1)
    return pooled.float().cpu().numpy()


def chunk(xs, n):
    for i in range(0, len(xs), n):
        yield xs[i:i + n]


def enumerate_kmers(seq, ks=(8, 9, 10, 11)):
    s = str(seq).upper()
    out = []
    for k in ks:
        for i in range(0, len(s) - k + 1):
            sub = s[i:i + k]
            if all(c in "ACDEFGHIKLMNPQRSTVWY" for c in sub):
                out.append(sub)
    return out


def main():
    t0 = time.time()
    print(f"loading {MODEL_NAME} on {DEVICE}…")
    tok, model = load_model()
    print(f"loaded in {time.time()-t0:.1f}s")

    bundle = pd.read_csv(BUNDLE, sep="\t")
    print(f"bundle: {len(bundle)} rows, splits={bundle['split'].value_counts().to_dict()}")

    # ---- Peptides ------------------------------------------------------------
    # For non-Venus rows: the peptide is the actual peptide
    pep_set = set()
    nonvenus = bundle[~bundle["split"].isin(["ext_venus_test", "ext_venus_valid"])].copy()
    for p in nonvenus["peptide"]:
        if isinstance(p, str) and 7 < len(p) <= 16 and all(c in "ACDEFGHIKLMNPQRSTVWY" for c in p):
            pep_set.add(p)
    print(f"unique non-Venus peptides: {len(pep_set)}")

    # ---- Venus enumeration ---------------------------------------------------
    venus_rows = bundle[bundle["split"].isin(["ext_venus_test", "ext_venus_valid"])].copy()
    venus_windows = {}
    for _, r in venus_rows.iterrows():
        pid = r["protein_id"]
        kmers = enumerate_kmers(r["peptide"])
        venus_windows[pid] = kmers
        pep_set.update(kmers)
    n_venus_windows = sum(len(v) for v in venus_windows.values())
    print(f"venus proteins: {len(venus_windows)}, total windows: {n_venus_windows}")
    print(f"unique-peptide universe (incl venus windows): {len(pep_set)}")

    pep_keys = sorted(pep_set)
    pep_emb = np.zeros((len(pep_keys), HIDDEN), dtype=np.float32)

    t1 = time.time()
    for i, batch in enumerate(chunk(pep_keys, BATCH)):
        e = embed_batch(tok, model, batch)
        pep_emb[i*BATCH:i*BATCH + len(batch)] = e
        if (i + 1) % 50 == 0:
            elapsed = time.time() - t1
            done = (i + 1) * BATCH
            rate = done / elapsed
            print(f"  pep batch {i+1}: {done}/{len(pep_keys)}  rate={rate:.0f}/s  eta={(len(pep_keys)-done)/rate:.0f}s")
    print(f"peptide embedding done: {time.time()-t1:.1f}s")

    # ---- HLA pseudo-seqs -----------------------------------------------------
    pseu = pd.read_csv(HLA_PSEUDO_TSV, sep="\t")
    pseu = pseu[pseu["pseudo_seq"].astype(str).str.len() == 34].reset_index(drop=True)
    hla_keys = pseu["HLA_norm"].tolist()
    hla_emb = np.zeros((len(hla_keys), HIDDEN), dtype=np.float32)
    seqs = pseu["pseudo_seq"].tolist()
    t2 = time.time()
    for i, batch in enumerate(chunk(seqs, BATCH)):
        e = embed_batch(tok, model, batch)
        hla_emb[i*BATCH:i*BATCH + len(batch)] = e
    print(f"HLA embedding done: {time.time()-t2:.1f}s  ({len(hla_keys)} alleles)")

    # ---- Save -----------------------------------------------------------------
    payload = {
        "pep_emb": torch.from_numpy(pep_emb),
        "pep_keys": pep_keys,
        "hla_emb": torch.from_numpy(hla_emb),
        "hla_keys": hla_keys,
        "venus_windows": {k: v for k, v in venus_windows.items()},
        "model_name": MODEL_NAME,
        "hidden": HIDDEN,
    }
    torch.save(payload, OUT)
    print(f"\nsaved → {OUT}  ({OUT.stat().st_size/1e6:.1f} MB)")
    print(f"total wall: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
