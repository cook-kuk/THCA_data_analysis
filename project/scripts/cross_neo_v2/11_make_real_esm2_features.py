#!/usr/bin/env python3
"""Generate frozen ESM2-35M compact features for CROSS-Neo 2.0.

This is a CPU-safe first real-PLM pass. It uses local fair-esm checkpoints if
present and falls back cleanly if unavailable. No labels or split assignments
enter embedding generation; the projection is deterministic and label-free.
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import pandas as pd
import torch

from common import OUT, SEED, ensure_dirs, safe_parquet


AA = set("ACDEFGHIKLMNPQRSTVWY")
MODEL_CHOICES = {
    "35m": ("esm2_t12_35M_UR50D", 12, 480),
    "150m": ("esm2_t30_150M_UR50D", 30, 640),
    "650m": ("esm2_t33_650M_UR50D", 33, 1280),
}
REQUESTED_TAG = os.environ.get("CROSS_NEO_ESM2_TAG", os.environ.get("CROSS_NEO_ESM2_MODEL", "35m")).lower()
TAG = REQUESTED_TAG if REQUESTED_TAG in MODEL_CHOICES else "35m"
MODEL_NAME, LAYER, HIDDEN = MODEL_CHOICES[TAG]
PROJ_DIM = 192
BATCH = 24


def clean_seq(x: object) -> str:
    s = "".join([c for c in str(x or "").upper() if c in AA])
    return s if 0 < len(s) <= 80 else ""


def load_esm2():
    import esm

    model, alphabet = getattr(esm.pretrained, MODEL_NAME)()
    model.eval()
    for p in model.parameters():
        p.requires_grad_(False)
    return model, alphabet.get_batch_converter()


@torch.no_grad()
def embed_sequences(model, batch_converter, seqs: list[str]) -> dict[str, np.ndarray]:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    out: dict[str, np.ndarray] = {}
    valid = sorted({s for s in seqs if s})
    for i in range(0, len(valid), BATCH):
        batch = valid[i:i + BATCH]
        labels = [(f"s{j}", s) for j, s in enumerate(batch)]
        _, _, toks = batch_converter(labels)
        toks = toks.to(device)
        rep = model(toks, repr_layers=[LAYER], return_contacts=False)["representations"][LAYER].cpu()
        for j, s in enumerate(batch):
            out[s] = rep[j, 1:len(s) + 1].mean(0).numpy().astype(np.float32)
        if (i // BATCH) % 20 == 0:
            print(f"[esm2] embedded {min(i + BATCH, len(valid))}/{len(valid)}")
    return out


def cos(a: np.ndarray, b: np.ndarray) -> float:
    den = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / den) if den > 0 else 0.0


def main() -> None:
    ensure_dirs()
    t0 = time.time()
    reg = pd.read_csv(OUT / "canonical_registry.tsv", sep="\t")
    reg["mut_clean"] = reg["mutant_peptide"].map(clean_seq)
    reg["wt_clean"] = reg["wildtype_peptide"].map(clean_seq)
    reg["source_clean"] = reg["source_window"].map(clean_seq)
    seqs = reg["mut_clean"].tolist() + reg["wt_clean"].tolist() + reg["source_clean"].tolist()
    seqs = [s for s in seqs if s]

    status = {
        "model": MODEL_NAME,
        "layer": LAYER,
        "hidden": HIDDEN,
        "projection_dim": PROJ_DIM,
        "n_unique_sequences": len(set(seqs)),
        "device": "cuda" if torch.cuda.is_available() else "cpu",
        "requested_tag": REQUESTED_TAG,
        "resolved_tag": TAG,
        "status": "started",
    }
    try:
        torch.set_num_threads(max(1, min(8, torch.get_num_threads())))
        model, batch_converter = load_esm2()
        emb = embed_sequences(model, batch_converter, seqs)
        rng = np.random.default_rng(SEED)
        proj = rng.normal(0, 1 / np.sqrt(HIDDEN * 4), size=(HIDDEN * 4, PROJ_DIM)).astype(np.float32)
        rows = []
        zero = np.zeros(HIDDEN, dtype=np.float32)
        for _, r in reg.iterrows():
            mut = emb.get(r["mut_clean"], zero)
            wt = emb.get(r["wt_clean"], zero)
            src = emb.get(r["source_clean"], zero)
            delta = mut - wt
            full = np.concatenate([mut, wt, delta, src]).astype(np.float32)
            z = full @ proj
            row = {"row_id": r["row_id"]}
            for i, v in enumerate(z):
                row[f"esm2_{TAG}_proj_{i:03d}"] = float(v)
            row["esm2_mut_wt_cosine"] = cos(mut, wt)
            row["esm2_mut_source_cosine"] = cos(mut, src)
            row["esm2_wt_source_cosine"] = cos(wt, src)
            row["esm2_wt_missing"] = float(not bool(r["wt_clean"]))
            row["esm2_source_missing"] = float(not bool(r["source_clean"]))
            rows.append(row)
        feat = pd.DataFrame(rows).fillna(0)
        safe_parquet(feat, OUT / f"features/esm2_{TAG}_features.parquet")
        feat.to_csv(OUT / f"features/esm2_{TAG}_features.tsv", sep="\t", index=False)
        status |= {"status": "ok", "feature_shape": list(feat.shape), "wall_seconds": round(time.time() - t0, 2)}
    except Exception as exc:
        status |= {"status": "failed", "error": repr(exc), "wall_seconds": round(time.time() - t0, 2)}
        pd.DataFrame({"row_id": reg["row_id"], f"esm2_{TAG}_unavailable": 1.0}).to_csv(OUT / f"features/esm2_{TAG}_features.tsv", sep="\t", index=False)
    (OUT / f"features/esm2_{TAG}_report.json").write_text(json.dumps(status, indent=2) + "\n")
    (OUT / f"features/esm2_{TAG}_report.md").write_text(
        f"# ESM2-{TAG} Feature Report\n\n"
        f"- Status: {status['status']}\n"
        f"- Model: {MODEL_NAME}\n"
        f"- Device: {status['device']}\n"
        f"- Unique sequences: {status['n_unique_sequences']}\n"
        f"- Wall seconds: {status.get('wall_seconds')}\n"
        "- Claim boundary: frozen PLM features only; no fine-tuning and no public predictor scores.\n"
    )
    print(f"[v2-esm2-{TAG}] status={status['status']} wall={status.get('wall_seconds')}s")


if __name__ == "__main__":
    main()
