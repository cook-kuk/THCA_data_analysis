#!/usr/bin/env python3
"""Generate CROSS-Neo 2.0 feature families."""

from __future__ import annotations

import json
import numpy as np
import pandas as pd

from common import OUT, V0, aa_physchem, ensure_dirs, safe_parquet, stable_hash, write_json


def hashed_plm(seq: object, prefix: str, n: int = 128) -> dict[str, float]:
    s = str(seq or "")
    arr = np.zeros(n, dtype=float)
    for k in (1, 2, 3, 4):
        if len(s) < k:
            continue
        for i in range(len(s) - k + 1):
            km = s[i:i+k]
            arr[stable_hash(f"{k}:{km}", n)] += 1
    arr = arr / max(1.0, np.linalg.norm(arr))
    return {f"{prefix}_emb_{i:03d}": float(v) for i, v in enumerate(arr)}


def main() -> None:
    ensure_dirs()
    reg = pd.read_csv(OUT / "canonical_registry.tsv", sep="\t")
    rows_cf, rows_plm, rows_proc = [], [], []
    for _, r in reg.iterrows():
        mut = r["mutant_peptide"]
        wt = "" if pd.isna(r["wildtype_peptide"]) else str(r["wildtype_peptide"])
        cf = {"row_id": r["row_id"]}
        mut_f = aa_physchem(mut, "mut")
        wt_f = aa_physchem(wt, "wt")
        cf |= mut_f | wt_f
        for k in mut_f:
            wk = k.replace("mut_", "wt_", 1)
            if wk in wt_f:
                cf[k.replace("mut_", "delta_", 1)] = mut_f[k] - wt_f[wk]
                cf[k.replace("mut_", "absdelta_", 1)] = abs(mut_f[k] - wt_f[wk])
        cf["hla_gene_hash"] = stable_hash(str(r["hla_gene"]), 17)
        cf["hla_supertype_hash"] = stable_hash(str(r["hla_supertype"]), 67)
        cf["peptide_length"] = float(r["peptide_length"])
        cf["wt_available"] = float(bool(wt))
        cf["anchor_position_change_unknown"] = 1.0
        rows_cf.append(cf)

        plm = {"row_id": r["row_id"]}
        plm |= hashed_plm(mut, "mut")
        plm |= hashed_plm(wt, "wt")
        for i in range(128):
            plm[f"delta_emb_{i:03d}"] = plm[f"mut_emb_{i:03d}"] - plm[f"wt_emb_{i:03d}"]
        plm |= hashed_plm(str(r["hla_4digit"]) + ":" + str(mut), "pmhc")
        sw = "" if pd.isna(r["source_window"]) else str(r["source_window"])
        plm |= hashed_plm(sw, "source_window")
        rows_plm.append(plm)

        rows_proc.append({
            "row_id": r["row_id"],
            "processing_has_source_window": float(pd.notna(r["source_window"])),
            "processing_has_source_protein": float(pd.notna(r["source_protein"])),
            "processing_has_expression": float(pd.notna(r.get("expression", pd.NA)) or pd.notna(r.get("TPM", pd.NA))),
            "processing_has_vaf": float(pd.notna(r.get("VAF", pd.NA))),
            "processing_len": float(r["peptide_length"]),
        })

    cf_df = pd.DataFrame(rows_cf)
    plm_df = pd.DataFrame(rows_plm)
    proc_df = pd.DataFrame(rows_proc)

    # Preserve v1/v0 counterfactual embedding features if available.
    if (V0 / "counterfactual_embeddings.npy").exists():
        arr = np.load(V0 / "counterfactual_embeddings.npy")
        idx = pd.read_csv(V0 / "counterfactual_feature_index.tsv", sep="\t")
        keep = min(160, arr.shape[1])
        emb = pd.DataFrame(arr[:, :keep], columns=[f"v1_cf_emb_{i:03d}" for i in range(keep)])
        emb.insert(0, "row_id", idx["sample_id"].astype(str).values)
        cf_df = cf_df.merge(emb, on="row_id", how="left")
    cf_df = cf_df.fillna(0)
    plm_df = plm_df.fillna(0)
    proc_df = proc_df.fillna(0)

    geom_path = V0 / "structure_geometry_features.tsv"
    if geom_path.exists():
        geom = pd.read_csv(geom_path, sep="\t").rename(columns={"sample_id": "row_id"})
        keep = ["row_id"] + [c for c in geom.columns if c not in {"row_id", "peptide_mut", "hla", "label", "strict_set_flag"}]
        struct = geom[keep].copy()
    else:
        struct = pd.DataFrame({"row_id": reg["row_id"], "structure_missing": 1})
    for c in struct.columns:
        if c != "row_id":
            struct[c] = pd.to_numeric(struct[c], errors="coerce").fillna(0)

    safe_parquet(cf_df, OUT / "features/counterfactual_features.parquet")
    safe_parquet(plm_df, OUT / "features/plm_embeddings.parquet")
    safe_parquet(struct, OUT / "features/structure_features.parquet")
    safe_parquet(proc_df, OUT / "features/processing_context_features.parquet")
    cf_df.to_csv(OUT / "features/counterfactual_features.tsv", sep="\t", index=False, na_rep="NA")
    plm_df.iloc[:, :80].to_csv(OUT / "features/plm_embeddings_preview.tsv", sep="\t", index=False, na_rep="NA")
    struct.to_csv(OUT / "features/structure_features.tsv", sep="\t", index=False, na_rep="NA")
    proc_df.to_csv(OUT / "features/processing_context_features.tsv", sep="\t", index=False, na_rep="NA")

    train_corpus = reg[["row_id", "peptide", "hla_4digit", "source_dataset", "strict_set_flag"]].copy()
    train_corpus["pmhc_text"] = train_corpus["hla_4digit"].astype(str) + " " + train_corpus["peptide"].astype(str)
    (OUT / "corpora/pMHC_mlm_train.txt").write_text("\n".join(train_corpus["pmhc_text"].tolist()) + "\n")
    train_corpus[["row_id", "peptide", "hla_4digit", "source_dataset"]].to_csv(OUT / "corpora/peptide_hla_pairs_train.tsv", sep="\t", index=False, na_rep="NA")
    train_corpus.groupby("source_dataset").size().reset_index(name="n").to_csv(OUT / "corpora/corpus_manifest.tsv", sep="\t", index=False, na_rep="NA")

    pd.DataFrame({
        "job_type": ["structure_prediction_manifest"],
        "n_rows": [len(reg)],
        "recommended_tools": ["Boltz/Chai/AF3/PANDORA if available"],
        "input_table": [str(OUT / "canonical_registry.tsv")],
    }).to_csv(OUT / "structure_jobs/boltz_or_chai_or_af3_jobs.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame({
        "job_type": ["template_structure_manifest"],
        "n_rows": [len(reg)],
        "recommended_tools": ["PANDORA/template pMHC if available"],
        "input_table": [str(OUT / "canonical_registry.tsv")],
    }).to_csv(OUT / "structure_jobs/pandora_or_template_jobs.tsv", sep="\t", index=False, na_rep="NA")

    write_json(OUT / "features/feature_manifest.json", {
        "counterfactual_features": cf_df.shape,
        "plm_embeddings": plm_df.shape,
        "structure_features": struct.shape,
        "processing_context_features": proc_df.shape,
        "plm_status": "hashed AA/kmer fallback; ESM2/ProtT5 job should replace for high-compute run",
    })
    (OUT / "structure_report.md").write_text("# Structure Report\n\nExisting v0 structure geometry was reused where available. Missing structures are explicit. Boltz/Chai/AF3 job manifests were created.\n")
    print(f"[v2-features] cf={cf_df.shape} plm={plm_df.shape} struct={struct.shape}")


if __name__ == "__main__":
    main()
