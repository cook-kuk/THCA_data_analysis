#!/usr/bin/env python
"""
Wave 8B-aux: ESMFold pseudo-complex structure prediction for ITSNdb peptides.
Uses HuggingFace transformers EsmForProteinFolding (no openfold dep needed).
Single-chain: peptide + GGGGS + HLA_pseudo_seq (34aa NetMHCpan contacts).
"""
import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37


LINKER = "GGGGS"


def load_inputs(bundle_tsv, itsndb_tsv, hla_tsv, n_train_anchor=200, max_len=11, seed=0):
    bundle = pd.read_csv(bundle_tsv, sep="\t")
    its = pd.read_csv(itsndb_tsv, sep="\t")
    hla = pd.read_csv(hla_tsv, sep="\t").dropna(subset=["pseudo_seq"])
    hla_map = dict(zip(hla["HLA_norm"], hla["pseudo_seq"]))

    its = its[its["peptide"].str.len().le(max_len)].copy()
    its["set"] = "itsndb"
    its["priority"] = its["in_master"].map(lambda x: 0 if str(x) == "False" else 1)
    its = its.sort_values("priority").reset_index(drop=True)

    train = bundle[(bundle["split"] == "train") & bundle["peptide"].str.len().le(max_len)].copy()
    rng = np.random.default_rng(seed)
    if len(train) > n_train_anchor:
        idx = rng.choice(len(train), size=n_train_anchor, replace=False)
        train = train.iloc[idx].copy()
    train["set"] = "train"
    train["in_master"] = ""
    train["priority"] = 2

    cols = ["peptide", "HLA_norm", "label", "set", "in_master", "priority"]
    df = pd.concat([its[cols], train[cols]], ignore_index=True)
    df["pseudo_seq"] = df["HLA_norm"].map(hla_map)
    df = df.dropna(subset=["pseudo_seq"]).copy()
    df["joined_seq"] = df.apply(lambda r: r["peptide"] + LINKER + r["pseudo_seq"], axis=1)
    df["pep_len"] = df["peptide"].str.len()
    df["pseudo_len"] = df["pseudo_seq"].str.len()
    return df


def convert_outputs_to_pdb(outputs):
    """From HF docs: convert model output dict -> PDB string."""
    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
    outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
    final_atom_positions = final_atom_positions.cpu().numpy()
    final_atom_mask = outputs["atom37_atom_exists"]
    pdbs = []
    for i in range(outputs["aatype"].shape[0]):
        aa = outputs["aatype"][i]
        pred_pos = final_atom_positions[i]
        mask = final_atom_mask[i]
        resid = outputs["residue_index"][i] + 1
        pred = Protein(
            aatype=aa,
            atom_positions=pred_pos,
            atom_mask=mask,
            residue_index=resid,
            b_factors=outputs["plddt"][i],
            chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
        )
        pdbs.append(to_pdb(pred))
    return pdbs


def extract_features(coords_ca, plddt_per_res, pep_len, pseudo_len, linker_len=5):
    """coords_ca: (L,3) numpy CA coords. plddt_per_res: (L,) in [0,100]."""
    total = pep_len + linker_len + pseudo_len
    if len(coords_ca) < pep_len + 1:
        return None
    coords = coords_ca
    plddt = plddt_per_res / 100.0  # normalize to [0,1]

    pep_idx = np.arange(pep_len)
    hla_start = pep_len + linker_len
    hla_idx = np.arange(hla_start, min(len(coords), total))

    pep_coords = coords[pep_idx]
    hla_coords = coords[hla_idx] if len(hla_idx) > 0 else np.zeros((0, 3))
    pep_plddt = plddt[pep_idx]
    hla_plddt = plddt[hla_idx] if len(hla_idx) > 0 else np.array([np.nan])

    feats = {}
    feats["mean_pLDDT_peptide"] = float(np.mean(pep_plddt))
    feats["min_pLDDT_peptide"] = float(np.min(pep_plddt))
    feats["mean_pLDDT_HLA"] = float(np.mean(hla_plddt)) if len(hla_plddt) else np.nan
    if pep_len >= 2:
        feats["anchor_pLDDT"] = float(np.mean([pep_plddt[1], pep_plddt[-1]]))
    else:
        feats["anchor_pLDDT"] = float(pep_plddt[0])

    if len(hla_coords) > 0:
        d = np.linalg.norm(pep_coords[:, None, :] - hla_coords[None, :, :], axis=2)
        feats["interface_contacts_8A"] = int((d < 8.0).sum())
        feats["n_buried_residues"] = int(((d < 8.0).sum(axis=1) >= 3).sum())
    else:
        feats["interface_contacts_8A"] = 0
        feats["n_buried_residues"] = 0

    com = pep_coords.mean(axis=0)
    rg = np.sqrt(np.mean(np.sum((pep_coords - com) ** 2, axis=1)))
    feats["radius_of_gyration_peptide"] = float(rg)

    if pep_len >= 4:
        i3 = np.linalg.norm(pep_coords[:-3] - pep_coords[3:], axis=1)
        feats["peptide_helicity"] = float(np.mean((i3 > 4.5) & (i3 < 6.5)))
    else:
        feats["peptide_helicity"] = 0.0
    return feats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--itsndb", required=True)
    ap.add_argument("--hla", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-min", type=float, default=40.0)
    ap.add_argument("--n-anchor", type=int, default=200)
    args = ap.parse_args()

    print("[esm] loading HF EsmForProteinFolding (esmfold_v1)...", flush=True)
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model = model.cuda()
    model.esm = model.esm.half()
    model.trunk.set_chunk_size(64)
    model = model.eval()
    torch.set_grad_enabled(False)
    print(f"[esm] loaded in {time.time()-t0:.1f}s", flush=True)

    df = load_inputs(args.bundle, args.itsndb, args.hla, n_train_anchor=args.n_anchor)
    print(f"[input] n_total={len(df)}; n_itsndb={(df['set']=='itsndb').sum()}; n_train={(df['set']=='train').sum()}", flush=True)
    print(f"[input] no_overlap (in_master=False): {(df['in_master'].astype(str)=='False').sum()}", flush=True)

    out_rows = []
    pdb_dir = os.path.dirname(args.out) + "/pdb_files"
    os.makedirs(pdb_dir, exist_ok=True)

    deadline = time.time() + args.max_min * 60.0
    n_done = 0
    n_skip = 0
    t_fold0 = time.time()
    save_interval = 25  # snapshot tsv every 25 folds

    for i, row in df.iterrows():
        if time.time() > deadline:
            print(f"[stop] deadline hit at i={i}/{len(df)}", flush=True)
            break
        seq = row["joined_seq"]
        try:
            inputs = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
            inputs = {k: v.cuda() for k, v in inputs.items()}
            with torch.no_grad():
                out = model(inputs["input_ids"])
            # CA index = atom14 idx 1
            positions = out["positions"][-1, 0].cpu().numpy()  # (L, 14, 3)
            ca = positions[:, 1, :]  # CA
            plddt = out["plddt"][0].mean(dim=-1).cpu().numpy()  # (L,)

            feats = extract_features(ca, plddt, row["pep_len"], row["pseudo_len"])
            if feats is None:
                n_skip += 1
                continue
            rec = {"peptide": row["peptide"], "HLA_norm": row["HLA_norm"], "label": row["label"],
                   "set": row["set"], "in_master": row["in_master"]}
            rec.update(feats)
            out_rows.append(rec)

            if n_done < 5:
                # save first few PDBs for spot-check
                try:
                    pdb_str = convert_outputs_to_pdb({k: v.detach() if torch.is_tensor(v) else v for k, v in out.items()})[0]
                    fn = f"{pdb_dir}/{row['peptide']}_{row['HLA_norm'].replace(':','-').replace('*','-')}.pdb"
                    with open(fn, "w") as f:
                        f.write(pdb_str)
                except Exception:
                    pass
            n_done += 1
            if n_done % 10 == 0:
                rate = (time.time() - t_fold0) / n_done
                rem_budget = (deadline - time.time())
                print(f"[fold] {n_done} done | avg {rate:.2f}s/fold | budget {rem_budget/60:.1f}min", flush=True)
            if n_done % save_interval == 0:
                pd.DataFrame(out_rows).to_csv(args.out, sep="\t", index=False)
        except torch.cuda.OutOfMemoryError:
            print(f"[oom] {row['peptide']}/{row['HLA_norm']}", flush=True)
            n_skip += 1
            torch.cuda.empty_cache()
        except Exception as e:
            n_skip += 1
            if n_skip < 8:
                print(f"[err] {row['peptide']}/{row['HLA_norm']}: {repr(e)[:200]}", flush=True)

    feat_df = pd.DataFrame(out_rows)
    feat_df.to_csv(args.out, sep="\t", index=False)
    print(f"[done] folded={n_done} skipped={n_skip} | wrote {args.out}", flush=True)
    print(f"[time] total {(time.time()-t0)/60:.1f}min", flush=True)


if __name__ == "__main__":
    main()
