#!/usr/bin/env python3
"""Generate CROSS-Neo-TCR feature families from the local TCR registry.

This script intentionally uses lightweight, deterministic fallbacks. External
TCR language-model or structure-prediction scores can replace these files later,
but the missingness and uncertainty flags are already explicit.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from common import AA, OUT, aa_physchem, kmer_set, safe_parquet, seq_similarity, stable_hash, write_json


TCR_OUT = OUT / "tcr_extension"
FEATURE_OUT = TCR_OUT / "features"


def clean(value: object) -> str:
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value or "").strip().upper()
    return "" if text in {"", "NA", "NAN", "NONE", "<NA>"} else text


def module_available(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


def hashed_embedding(seq: object, prefix: str, dim: int = 64) -> dict[str, float]:
    s = clean(seq)
    arr = np.zeros(dim, dtype=np.float32)
    for k in (1, 2, 3, 4):
        if len(s) < k:
            continue
        for i in range(len(s) - k + 1):
            arr[stable_hash(f"{prefix}:{k}:{s[i:i+k]}", dim)] += 1.0
    norm = float(np.linalg.norm(arr))
    if norm > 0:
        arr /= norm
    return {f"{prefix}_fallback_emb_{i:03d}": float(v) for i, v in enumerate(arr)}


def gene_hash_features(value: object, prefix: str, bins: int = 24) -> dict[str, float]:
    out = {f"{prefix}_missing": float(clean(value) == "")}
    text = clean(value)
    for i in range(bins):
        out[f"{prefix}_hash_{i:02d}"] = 0.0
    if text:
        out[f"{prefix}_hash_{stable_hash(text, bins):02d}"] = 1.0
    return out


def shared_kmer_fraction(a: object, b: object, k: int) -> float:
    ka = kmer_set(clean(a), k)
    kb = kmer_set(clean(b), k)
    return float(len(ka & kb) / len(ka | kb)) if ka and kb else 0.0


def sequence_feature_row(row) -> dict[str, float | str]:
    alpha = clean(row.cdr3_alpha)
    beta = clean(row.cdr3_beta)
    peptide = clean(row.peptide)
    out: dict[str, float | str] = {"row_id": row.row_id}
    out.update(aa_physchem(alpha, "cdr3a"))
    out.update(aa_physchem(beta, "cdr3b"))
    out.update(aa_physchem(peptide, "pep"))
    out.update(gene_hash_features(row.tcr_alpha_v, "trav"))
    out.update(gene_hash_features(row.tcr_alpha_j, "traj"))
    out.update(gene_hash_features(row.tcr_beta_v, "trbv"))
    out.update(gene_hash_features(row.tcr_beta_j, "trbj"))
    out["paired_alpha_beta_available"] = float(bool(alpha and beta))
    out["alpha_only_available"] = float(bool(alpha and not beta))
    out["beta_only_available"] = float(bool(beta and not alpha))
    out["any_cdr3_available"] = float(bool(alpha or beta))
    out["no_tcr_sequence"] = float(not bool(alpha or beta))
    out["peptide_available"] = float(bool(peptide))
    out["peptide_hla_available"] = float(bool(peptide and clean(row.hla_4digit)))
    out["hla_supertype_hash"] = float(stable_hash(clean(row.hla_supertype), 67)) if clean(row.hla_supertype) else 0.0
    out["peptide_cdr3a_similarity"] = seq_similarity(peptide, alpha)
    out["peptide_cdr3b_similarity"] = seq_similarity(peptide, beta)
    for k in (2, 3, 4, 5):
        out[f"pep_cdr3a_k{k}_jaccard"] = shared_kmer_fraction(peptide, alpha, k)
        out[f"pep_cdr3b_k{k}_jaccard"] = shared_kmer_fraction(peptide, beta, k)
    out["cdr3a_cdr3b_similarity"] = seq_similarity(alpha, beta)
    out["cdr3_total_len"] = float(len(alpha) + len(beta))
    out["cdr3_len_abs_delta"] = float(abs(len(alpha) - len(beta))) if alpha and beta else 0.0
    return out


def embedding_row(row, dim: int) -> dict[str, float | str]:
    alpha = clean(row.cdr3_alpha)
    beta = clean(row.cdr3_beta)
    peptide = clean(row.peptide)
    hla = clean(row.hla_4digit)
    out: dict[str, float | str] = {"row_id": row.row_id}
    out.update(hashed_embedding(alpha, "cdr3a", dim))
    out.update(hashed_embedding(beta, "cdr3b", dim))
    out.update(hashed_embedding(peptide, "pep", dim))
    out.update(hashed_embedding(f"{hla}:{peptide}:{alpha}:{beta}", "tcr_pmhc_cross", dim))
    out["embedding_is_fallback_kmer"] = 1.0
    out["embedding_has_alpha"] = float(bool(alpha))
    out["embedding_has_beta"] = float(bool(beta))
    out["embedding_has_peptide_hla"] = float(bool(peptide and hla))
    return out


def structure_feature_row(row) -> dict[str, float | str]:
    pdb = clean(row.structure_pdb_id)
    paired = bool(clean(row.cdr3_alpha) and clean(row.cdr3_beta))
    peptide = clean(row.peptide)
    hla = clean(row.hla_4digit)
    return {
        "row_id": row.row_id,
        "structure_pdb_id_present": float(bool(pdb)),
        "no_structure": float(not bool(pdb)),
        "low_structure_confidence": float(not bool(pdb)),
        "structure_claimable_interface": 0.0,
        "structure_known_template_or_pdb": float(bool(pdb)),
        "structure_model_generated": 0.0,
        "tcr_pmhc_modelable_now": float(bool(paired and peptide and hla)),
        "paired_tcr_missing_for_structure": float(not paired),
        "peptide_hla_missing_for_structure": float(not bool(peptide and hla)),
        "interface_contact_count": np.nan,
        "peptide_tcr_contact_count": np.nan,
        "cdr3a_peptide_contact_count": np.nan,
        "cdr3b_peptide_contact_count": np.nan,
        "crossing_angle": np.nan,
        "mutation_residue_contact_flag": np.nan,
        "mutant_wt_interface_delta": np.nan,
        "wt_cross_reactivity_risk": np.nan,
    }


def read_registry() -> pd.DataFrame:
    path = TCR_OUT / "tcr_registry.parquet"
    if path.exists():
        return pd.read_parquet(path)
    return pd.read_csv(TCR_OUT / "tcr_registry.tsv", sep="\t", low_memory=False)


def numeric_frame(rows: list[dict[str, float | str]]) -> pd.DataFrame:
    df = pd.DataFrame(rows).fillna(0)
    for c in df.columns:
        if c != "row_id":
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(np.float32)
    return df


def write_chunked_parquet(
    source: pd.DataFrame,
    path: Path,
    row_builder,
    chunk_size: int,
    preview_path: Path,
) -> tuple[int, int]:
    if path.exists():
        path.unlink()
    writer: pq.ParquetWriter | None = None
    preview_written = False
    n_rows = 0
    n_cols = 0
    for start in range(0, len(source), chunk_size):
        chunk = source.iloc[start : start + chunk_size]
        frame = numeric_frame([row_builder(r) for r in chunk.itertuples(index=False)])
        table = pa.Table.from_pandas(frame, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(path, table.schema, compression="zstd")
            n_cols = len(frame.columns)
        writer.write_table(table)
        if not preview_written:
            frame.head(2000).to_csv(preview_path, sep="\t", index=False, na_rep="NA")
            preview_written = True
        n_rows += len(frame)
    if writer is not None:
        writer.close()
    else:
        numeric_frame([]).to_csv(preview_path, sep="\t", index=False, na_rep="NA")
    return n_rows, n_cols


def main() -> None:
    FEATURE_OUT.mkdir(parents=True, exist_ok=True)
    tcr = read_registry()
    for col in [
        "row_id",
        "source_dataset",
        "peptide",
        "hla_4digit",
        "hla_supertype",
        "tcr_alpha_v",
        "tcr_alpha_j",
        "cdr3_alpha",
        "tcr_beta_v",
        "tcr_beta_j",
        "cdr3_beta",
        "structure_pdb_id",
    ]:
        if col not in tcr:
            tcr[col] = ""

    emb_source = tcr[
        tcr["cdr3_alpha"].fillna("").astype(str).ne("")
        | tcr["cdr3_beta"].fillna("").astype(str).ne("")
    ].copy()
    chunk_size = 20000
    seq_shape = write_chunked_parquet(
        tcr,
        FEATURE_OUT / "tcr_sequence_features.parquet",
        sequence_feature_row,
        chunk_size,
        FEATURE_OUT / "tcr_sequence_features_preview.tsv",
    )
    emb_shape = write_chunked_parquet(
        emb_source,
        FEATURE_OUT / "tcr_plm_embeddings.parquet",
        lambda r: embedding_row(r, dim=64),
        chunk_size,
        FEATURE_OUT / "tcr_plm_embeddings_preview.tsv",
    )
    struct_shape = write_chunked_parquet(
        tcr,
        FEATURE_OUT / "tcr_structure_features.parquet",
        structure_feature_row,
        chunk_size,
        FEATURE_OUT / "tcr_structure_features_preview.tsv",
    )

    alpha = tcr["cdr3_alpha"].fillna("").astype(str).str.strip().ne("")
    beta = tcr["cdr3_beta"].fillna("").astype(str).str.strip().ne("")
    peptide = tcr["peptide"].fillna("").astype(str).str.strip().ne("")
    hla = tcr["hla_4digit"].fillna("").astype(str).str.strip().ne("")
    pdb = tcr["structure_pdb_id"].fillna("").astype(str).str.strip().ne("")

    availability = {
        "TCR_BERT_local_module": module_available("tcr_bert"),
        "TCRpeg_local_module": module_available("tcrpeg"),
        "transformers_local_module": module_available("transformers"),
        "torch_local_module": module_available("torch"),
        "esm_local_module": module_available("esm"),
        "prot_t5_status": "not_run; fallback hashed k-mer embeddings emitted",
        "embedding_dim_per_family": 64,
        "embedding_rows": int(emb_shape[0]),
        "sequence_rows": int(seq_shape[0]),
        "structure_rows": int(struct_shape[0]),
    }
    counts = {
        "registry_rows": int(len(tcr)),
        "paired_alpha_beta": int((alpha & beta).sum()),
        "beta_only": int((beta & ~alpha).sum()),
        "alpha_only": int((alpha & ~beta).sum()),
        "any_cdr3": int((alpha | beta).sum()),
        "peptide_hla_available": int((peptide & hla).sum()),
        "known_structure_or_pdb": int(pdb.sum()),
        "modelable_paired_tcr_pmhc": int((alpha & beta & peptide & hla).sum()),
    }
    write_json(
        FEATURE_OUT / "tcr_feature_manifest.json",
        {
            "outputs": {
                "tcr_sequence_features": str(FEATURE_OUT / "tcr_sequence_features.parquet"),
                "tcr_plm_embeddings": str(FEATURE_OUT / "tcr_plm_embeddings.parquet"),
                "tcr_structure_features": str(FEATURE_OUT / "tcr_structure_features.parquet"),
            },
            "counts": counts,
            "local_model_availability": availability,
            "claim_status": "diagnostic_feature_layer; fallback embeddings are not TCR-aware SOTA",
        },
    )

    lines = [
        "# CROSS-Neo-TCR Feature Report",
        "",
        "Generated deterministic TCR sequence, fallback k-mer embedding, and structure-missingness feature tables from the local TCR registry.",
        "",
        "## Counts",
        "",
        "| Field | Count |",
        "|---|---:|",
    ]
    lines.extend([f"| {k} | {v} |" for k, v in counts.items()])
    lines.extend(
        [
            "",
            "## Local PLM / Tool Status",
            "",
            "| Component | Status |",
            "|---|---|",
        ]
    )
    lines.extend([f"| {k} | {v} |" for k, v in availability.items()])
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "- `tcr_plm_embeddings.parquet` currently contains deterministic hashed k-mer fallback embeddings, not trained TCR-BERT/TCRpeg/ESM2/ProtT5 outputs.",
            "- Structure features are explicit missingness/QC placeholders unless a PDB/template exists; they should not be interpreted as interface evidence without parsed structures.",
            "- These features are suitable for diagnostic pilots and wetlab prioritization triage, not as a standalone claim of full TCR-aware neoantigen SOTA.",
        ]
    )
    (TCR_OUT / "tcr_feature_report.md").write_text("\n".join(lines) + "\n")
    print(f"[tcr-features] seq={seq_shape} emb={emb_shape} struct={struct_shape}")


if __name__ == "__main__":
    main()
