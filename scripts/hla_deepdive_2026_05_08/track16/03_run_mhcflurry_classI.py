#!/usr/bin/env python3
"""
Track 16 Step 3: MHCflurry-2.0 (presentation models) HLA-I binding for the 14
thyroid self-antigens × top-10 Korean HLA-A/B/C alleles (from K2 freq table).
Boundary: thyroid autoimmunity peptidomics — not cancer neoantigen prediction.
"""
from __future__ import annotations
import sys
import time
from pathlib import Path
import pandas as pd
from mhcflurry import Class1PresentationPredictor

OUTDIR = Path("/home/seungho/personal/THCA_data_analysis/project/results/hla_deepdive_2026_05_08/track16_peptidomics")
PEP_PATH = OUTDIR / "tables" / "T02_peptide_library_9mer.tsv"

# Top 10 per locus from K2 (T03_K2_allele_freq_per_locus.tsv) — all MHCflurry-supported
TOP_A = ["A*24:02", "A*33:03", "A*02:06", "A*02:01", "A*02:07", "A*31:01", "A*26:01", "A*11:01", "A*30:01", "A*32:01"]
TOP_B = ["B*51:01", "B*44:03", "B*15:01", "B*46:01", "B*40:02", "B*58:01", "B*35:01", "B*54:01", "B*40:01", "B*15:18"]
TOP_C = ["C*01:02", "C*14:02", "C*03:03", "C*08:01", "C*07:02", "C*03:04", "C*04:01", "C*06:02", "C*12:02", "C*15:02"]

# B*15:18 is sometimes labelled separately in MHCflurry; substitute B*55:02 if missing
TOP_ALLELES_RAW = TOP_A + TOP_B + TOP_C


def main() -> int:
    print("Loading MHCflurry presentation predictor ...")
    pred = Class1PresentationPredictor.load()
    supported = set(pred.supported_alleles)
    alleles = []
    for a in TOP_ALLELES_RAW:
        full = "HLA-" + a
        if full in supported:
            alleles.append(full)
        else:
            print(f"  WARNING: {full} not in MHCflurry — skipping", file=sys.stderr)
    print(f"Will predict on {len(alleles)} alleles.")

    pep_df = pd.read_csv(PEP_PATH, sep="\t")
    print(f"Loaded {len(pep_df):,} 9-mer peptides across {pep_df.gene.nunique()} antigens.")

    # MHCflurry handles peptide x allele "products" via predict()
    # Use predict_to_dataframe for efficiency
    t0 = time.time()
    chunks = []
    # Run per-allele to keep memory steady
    for a in alleles:
        ta = time.time()
        pred_df = pred.predict(
            peptides=pep_df["peptide"].tolist(),
            alleles=[a],
            include_affinity_percentile=True,
            verbose=0,
        )
        pred_df = pred_df.rename(columns={
            "best_allele": "allele",
            "affinity": "affinity_nM",
            "affinity_percentile": "affinity_percentile_rank",
        })
        # MHCflurry preserves input order in `peptide_num`; merge on peptide_num to be safe vs duplicate sequences
        pep_df_idx = pep_df.copy()
        pep_df_idx["peptide_num"] = range(len(pep_df_idx))
        merged = pep_df_idx.merge(
            pred_df[["peptide_num", "peptide", "allele", "affinity_nM", "affinity_percentile_rank",
                     "processing_score", "presentation_score", "presentation_percentile"]],
            on=["peptide_num", "peptide"],
            how="inner",
        ).drop(columns=["peptide_num"])
        chunks.append(merged)
        print(f"  {a}: n={len(merged):,} preds, dt={time.time()-ta:.1f}s")

    out = pd.concat(chunks, ignore_index=True)
    print(f"Total preds: {len(out):,} (elapsed {time.time()-t0:.1f}s)")

    # Tag binder bins on presentation_percentile (lower = stronger)
    out["binder_class"] = pd.cut(
        out["presentation_percentile"],
        bins=[-0.001, 0.5, 2.0, 100.0],
        labels=["SB", "WB", "NB"],
    )

    out_path = OUTDIR / "tables" / "T04_HLA_I_binding_long.tsv.gz"
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"wrote {out_path} ({out_path.stat().st_size/1e6:.1f} MB)")

    # Summary
    summ = (out.groupby(["allele", "binder_class"], observed=True)
              .size().unstack(fill_value=0).reset_index())
    summ_path = OUTDIR / "tables" / "T05_HLA_I_binders_per_allele.tsv"
    summ.to_csv(summ_path, sep="\t", index=False)
    print(f"wrote {summ_path}")

    by_ag = (out.groupby(["gene", "allele", "binder_class"], observed=True)
               .size().unstack(fill_value=0).reset_index())
    by_ag_path = OUTDIR / "tables" / "T06_HLA_I_binders_per_antigen_allele.tsv"
    by_ag.to_csv(by_ag_path, sep="\t", index=False)
    print(f"wrote {by_ag_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
