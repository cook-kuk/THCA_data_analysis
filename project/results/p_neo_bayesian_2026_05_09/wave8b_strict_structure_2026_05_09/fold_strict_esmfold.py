#!/usr/bin/env python3
"""Fold strict no-reference ITSNdb peptides with HF ESMFold.

Input rows must include:
  peptide, hla, label, pseudo_seq, joined_seq

The "complex" is a pragmatic single-chain proxy:
  peptide + GGGGS + 34-aa HLA contact pseudo-sequence

This is not a physical pMHC complex model. It is a structure-imputed proxy for
TCR-facing geometry when no TCR reference motif evidence is available.
"""
from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
from transformers.models.esm.openfold_utils.protein import Protein, to_pdb


LINKER = "GGGGS"


def to_pdb_strings(outputs: dict) -> list[str]:
    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
    out_np = {k: v.to("cpu").numpy() for k, v in outputs.items()}
    final_atom_positions = final_atom_positions.cpu().numpy()
    final_atom_mask = out_np["atom37_atom_exists"]
    pdbs = []
    for i in range(out_np["aatype"].shape[0]):
        pred = Protein(
            aatype=out_np["aatype"][i],
            atom_positions=final_atom_positions[i],
            atom_mask=final_atom_mask[i],
            residue_index=out_np["residue_index"][i] + 1,
            b_factors=out_np["plddt"][i],
            chain_index=out_np["chain_index"][i] if "chain_index" in out_np else None,
        )
        pdbs.append(to_pdb(pred))
    return pdbs


def normalize_plddt(raw: np.ndarray) -> np.ndarray:
    """HF has used both 0..1 and 0..100-like pLDDT conventions across APIs."""
    raw = np.asarray(raw, dtype=float)
    if np.nanmax(raw) <= 1.5:
        return raw * 100.0
    return raw


def extract_features(ca: np.ndarray, plddt: np.ndarray, pep_len: int, pseudo_len: int) -> dict:
    plddt = normalize_plddt(plddt)
    hla_start = pep_len + len(LINKER)
    hla_end = min(len(ca), hla_start + pseudo_len)
    pep_idx = np.arange(pep_len)
    hla_idx = np.arange(hla_start, hla_end)

    pep_ca = ca[pep_idx]
    hla_ca = ca[hla_idx] if len(hla_idx) else np.zeros((0, 3))
    pep_plddt = plddt[pep_idx]
    hla_plddt = plddt[hla_idx] if len(hla_idx) else np.array([np.nan])

    feats = {
        "mean_pLDDT_peptide": float(np.nanmean(pep_plddt)),
        "min_pLDDT_peptide": float(np.nanmin(pep_plddt)),
        "mean_pLDDT_HLA": float(np.nanmean(hla_plddt)),
        "anchor_pLDDT": float(np.nanmean([pep_plddt[1], pep_plddt[-1]])) if pep_len >= 2 else float(pep_plddt[0]),
    }

    if len(hla_ca):
        dist = np.linalg.norm(pep_ca[:, None, :] - hla_ca[None, :, :], axis=2)
        min_dist = dist.min(axis=1)
        feats["interface_contacts_8A"] = int((dist < 8.0).sum())
        feats["interface_contacts_10A"] = int((dist < 10.0).sum())
        feats["n_buried_residues_8A"] = int(((dist < 8.0).sum(axis=1) >= 3).sum())
        feats["mean_min_pep_to_hla_CA_dist"] = float(np.mean(min_dist))
        feats["max_min_pep_to_hla_CA_dist"] = float(np.max(min_dist))
    else:
        feats["interface_contacts_8A"] = 0
        feats["interface_contacts_10A"] = 0
        feats["n_buried_residues_8A"] = 0
        feats["mean_min_pep_to_hla_CA_dist"] = np.nan
        feats["max_min_pep_to_hla_CA_dist"] = np.nan

    center = pep_ca.mean(axis=0)
    centered = pep_ca - center
    feats["radius_of_gyration_peptide"] = float(np.sqrt(np.mean(np.sum(centered**2, axis=1))))
    feats["end_to_end_CA_dist"] = float(np.linalg.norm(pep_ca[0] - pep_ca[-1])) if pep_len >= 2 else 0.0
    if pep_len >= 4:
        i3 = np.linalg.norm(pep_ca[:-3] - pep_ca[3:], axis=1)
        feats["peptide_helicity_proxy"] = float(np.mean((i3 > 4.5) & (i3 < 6.5)))
    else:
        feats["peptide_helicity_proxy"] = 0.0
    return feats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--max-min", type=float, default=30.0)
    ap.add_argument("--save-pdb-n", type=int, default=12)
    args = ap.parse_args()

    inp = Path(args.input)
    out = Path(args.out)
    pdb_dir = Path(args.pdb_dir)
    pdb_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(inp, sep="\t")
    df["pep_len"] = df["peptide"].str.len()
    df["pseudo_len"] = df["pseudo_seq"].str.len()
    print(f"[input] rows={len(df)} pos={int(df['label'].sum())}", flush=True)

    print("[load] facebook/esmfold_v1", flush=True)
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model = model.cuda().eval()
    model.esm = model.esm.half()
    model.trunk.set_chunk_size(64)
    torch.set_grad_enabled(False)
    print(f"[load] done in {(time.time() - t0):.1f}s", flush=True)

    deadline = time.time() + args.max_min * 60
    rows = []
    fold_t0 = time.time()
    for i, row in df.iterrows():
        if time.time() > deadline:
            print(f"[stop] deadline at {i}/{len(df)}", flush=True)
            break
        seq = row["joined_seq"]
        try:
            tokens = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
            tokens = {k: v.cuda() for k, v in tokens.items()}
            with torch.no_grad():
                pred = model(tokens["input_ids"])
            pos14 = pred["positions"][-1, 0].detach().cpu().numpy()
            ca = pos14[:, 1, :]
            plddt_raw = pred["plddt"][0].detach().cpu().numpy()
            if plddt_raw.ndim > 1:
                plddt_raw = plddt_raw.mean(axis=-1)
            feats = extract_features(ca, plddt_raw, int(row["pep_len"]), int(row["pseudo_len"]))
            rec = {
                "peptide": row["peptide"],
                "hla": row["hla"],
                "label": int(row["label"]),
                "split": row["split"],
                "in_master": row["in_master"],
                "tcr_motif_score": row["tcr_motif_score"],
                "self_exact_match": row["self_exact_match"],
            }
            rec.update(feats)
            rows.append(rec)
            if len(rows) <= args.save_pdb_n:
                pdb = to_pdb_strings({k: v.detach() if torch.is_tensor(v) else v for k, v in pred.items()})[0]
                safe = f"{row['peptide']}_{row['hla'].replace('*','-').replace(':','-')}.pdb"
                (pdb_dir / safe).write_text(pdb)
            if len(rows) % 10 == 0:
                rate = (time.time() - fold_t0) / len(rows)
                print(f"[fold] {len(rows)}/{len(df)} avg={rate:.2f}s", flush=True)
                pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        except torch.cuda.OutOfMemoryError:
            print(f"[oom] {row['peptide']} {row['hla']}", flush=True)
            torch.cuda.empty_cache()
        except Exception as exc:
            print(f"[err] {row['peptide']} {row['hla']} {repr(exc)[:180]}", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(out, sep="\t", index=False)
    print(f"[done] folded={len(res)} wrote={out} total_min={(time.time() - t0)/60:.1f}", flush=True)


if __name__ == "__main__":
    main()
