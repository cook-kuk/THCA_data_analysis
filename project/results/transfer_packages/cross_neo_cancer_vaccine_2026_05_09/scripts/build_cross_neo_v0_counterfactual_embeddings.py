#!/usr/bin/env python3
"""Frozen counterfactual peptide/HLA embeddings for CROSS-Neo v0.

ESM2 tensors exist in the repo, but this v0 stores a deterministic lightweight
fallback representation so the full pipeline stays CPU-only and reproducible.
No PLM is fine-tuned.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from cross_neo_v0_common import INPUT, OUT, aa_features, ensure_dirs, load_master


def main() -> None:
    ensure_dirs()
    master = load_master()
    pseu = pd.read_csv(INPUT / "hla_pseudo.tsv", sep="\t").fillna("")
    hla_to_pseudo = dict(zip(pseu["HLA_norm"], pseu["pseudo_seq"]))

    feature_names = None
    mat = []
    index_rows = []
    for i, r in master.iterrows():
        mut = str(r["peptide_mut"])
        wt = str(r.get("peptide_wt", "") or "")
        hla_seq = hla_to_pseudo.get(str(r["hla"]), "")
        mut_names, mut_vec = aa_features(mut, "mut")
        wt_names, wt_vec = aa_features(wt, "wt")
        hla_names, hla_vec = aa_features(hla_seq, "hla")
        delta = mut_vec - wt_vec
        abs_delta = np.abs(delta)
        names = mut_names + wt_names + [f"delta_{n.replace('mut_', '')}" for n in mut_names] + [
            f"abs_delta_{n.replace('mut_', '')}" for n in mut_names
        ] + hla_names + ["anchor_mutation_flag", "wt_missing_flag", "hla_pseudo_missing_flag"]
        vals = np.concatenate(
            [
                mut_vec,
                wt_vec,
                delta,
                abs_delta,
                hla_vec,
                np.asarray([0.0, float(wt == ""), float(hla_seq == "")], dtype=float),
            ]
        )
        if feature_names is None:
            feature_names = names
        mat.append(vals)
        index_rows.append(
            {
                "sample_id": r["sample_id"],
                "embedding_row": i,
                "peptide_mut": mut,
                "peptide_wt_available": int(wt != ""),
                "hla_pseudo_available": int(hla_seq != ""),
                "embedding_source": "deterministic_AAindex_kmer_fallback_no_finetune",
            }
        )
    arr = np.vstack(mat).astype(np.float32)
    np.save(OUT / "counterfactual_embeddings.npy", arr)
    pd.DataFrame(index_rows).to_csv(OUT / "counterfactual_feature_index.tsv", sep="\t", index=False)
    pd.DataFrame({"feature_index": range(len(feature_names or [])), "feature": feature_names or []}).to_csv(
        OUT / "counterfactual_feature_names.tsv", sep="\t", index=False
    )
    report = [
        "# Counterfactual Embedding Report",
        "",
        "No PLM was fine-tuned.",
        "WT peptide is missing in the current master table, so WT vectors are zero and `wt_missing_flag=1`.",
        "The representation is a deterministic AA-property + k-mer fallback with HLA pseudo-sequence encoding.",
        f"Matrix shape: {arr.shape}",
    ]
    (OUT / "counterfactual_embedding_report.md").write_text("\n".join(report) + "\n")
    print(f"[counterfactual] wrote {arr.shape}")


if __name__ == "__main__":
    main()
