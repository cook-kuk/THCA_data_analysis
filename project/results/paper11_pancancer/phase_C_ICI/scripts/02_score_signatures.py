#!/usr/bin/env python3
"""Phase C v2 — Signature scoring per sample.

For each cohort:
  1. Load expr_matrix.tsv (gene_symbol × sample) and metadata.tsv
  2. z-score each gene WITHIN COHORT (so cross-cohort batch is not double-counted in score)
  3. For each module (from paper3_ici_module_gene_list.tsv), compute mean z across genes present
  4. Compute composite scores:
       DM1_inflam_composite = mean(IFNG_T_cell_inflamed, myeloid_suppressive,
                                   checkpoint_exhaustion, HLA_class_II)
       cytolytic_GZMA_PRF1  = mean( z(GZMA), z(PRF1) )    (Rooney et al. 2015)
       APM_score            = mean of HLA_class_I genes
  5. Lineage-portable DM1:
       LineageTF_z = mean z of cohort-appropriate lineage TF panel
                     (urothelial → GATA3, FOXA1, KRT20, PPARG, UPK1A, UPK3A, UPK3B)
                     (melanoma   → MITF, TYR, MLANA, DCT, PMEL, TYRP1)
       Inflam_z    = mean z of INFLAM_PANEL
       Mod4_z      = mean z of MOD4_CORE
       lineage_portable_DM1 = -LineageTF_z + Inflam_z + Mod4_z   (matches Phase E definition)

Output:
  results/tables/signature_scores_per_sample.tsv  (1 row per sample, columns = cohort + scores + metadata join)
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("/data/thca/repo_results/paper11_pancancer/phase_C_ICI")
PROC = ROOT / "processed"
OUT = ROOT / "results" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

PANEL = pd.read_csv("/home/seungho/personal/THCA_data_analysis/project/reports/paper3_ici/paper3_ici_module_gene_list.tsv", sep="\t")
MODULES = PANEL.groupby("module")["gene_symbol"].apply(list).to_dict()
print(f"[modules] {len(MODULES)} modules: {list(MODULES.keys())}")

COMPOSITE_INFLAM_MODULES = ["IFNG_T_cell_inflamed", "myeloid_suppressive",
                            "checkpoint_exhaustion", "HLA_class_II"]
CYTOLYTIC_GENES = ["GZMA", "PRF1"]

# Lineage-portable DM1 (matches phase_E semantics)
LINEAGE_TF = {
    "urothelial_carcinoma": ["GATA3", "FOXA1", "KRT20", "PPARG", "UPK1A", "UPK3A", "UPK3B"],
    "melanoma":             ["MITF", "TYR", "MLANA", "DCT", "PMEL", "TYRP1"],
}
INFLAM_PANEL = ["IFNG", "STAT1", "GZMB", "PRF1", "CXCL9", "CXCL10", "CD8A",
                "CD68", "CD163", "MRC1", "CSF1R", "CXCL8"]
MOD4_CORE    = ["MYC", "NAMPT", "JAK1", "JAK2", "STAT3", "LYN", "FYN", "SRC",
                "EZH2", "DNMT1", "HDAC1", "HDAC2"]


def zscore_genes_within_cohort(expr: pd.DataFrame) -> pd.DataFrame:
    """expr: gene × sample. Returns gene × sample z-scored across samples per gene."""
    mu = expr.mean(axis=1)
    sd = expr.std(axis=1).replace(0, np.nan)
    return expr.sub(mu, axis=0).div(sd, axis=0)


def score_module(z: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in z.index]
    if not present:
        return pd.Series(np.nan, index=z.columns)
    return z.loc[present].mean(axis=0)


def score_one_cohort(cohort_dir: Path) -> pd.DataFrame:
    expr = pd.read_csv(cohort_dir / "expr_matrix.tsv", sep="\t", index_col=0)
    md = pd.read_csv(cohort_dir / "metadata.tsv", sep="\t")
    cohort = cohort_dir.name

    # --- column dedup + restrict to metadata samples
    expr = expr.loc[:, ~expr.columns.duplicated()]
    expr = expr[[c for c in expr.columns if c in md["sample_id"].astype(str).values]]
    md = md[md["sample_id"].astype(str).isin(expr.columns)].copy()

    # --- coverage QC for module genes
    cov_rows = []
    for mod, genes in MODULES.items():
        present = [g for g in genes if g in expr.index]
        cov_rows.append({"module": mod, "n_total": len(genes), "n_present": len(present),
                         "fraction": round(len(present)/len(genes), 3),
                         "missing": ",".join(g for g in genes if g not in expr.index)})
    cov = pd.DataFrame(cov_rows).assign(cohort=cohort)
    print(f"[{cohort}] expr {expr.shape}, n_md={len(md)}; module coverage:")
    print(cov[["module","n_total","n_present","fraction"]].to_string(index=False))

    # --- z-score and module means
    z = zscore_genes_within_cohort(expr)
    out = pd.DataFrame(index=expr.columns)
    for mod in MODULES:
        out[mod] = score_module(z, MODULES[mod])
    out["DM1_inflam_composite"] = out[COMPOSITE_INFLAM_MODULES].mean(axis=1)
    out["cytolytic_GZMA_PRF1"] = score_module(z, CYTOLYTIC_GENES)

    # --- lineage-portable DM1
    cancer_type = md["cancer_type"].mode().iloc[0] if len(md) else "unknown"
    tf_panel = LINEAGE_TF.get(cancer_type, [])
    out["lineageTF_z"] = score_module(z, tf_panel)
    out["inflam_panel_z"] = score_module(z, INFLAM_PANEL)
    out["mod4_core_z"] = score_module(z, MOD4_CORE)
    out["lineage_portable_DM1"] = (-out["lineageTF_z"]
                                   + out["inflam_panel_z"]
                                   + out["mod4_core_z"])

    out = out.reset_index().rename(columns={"index": "sample_id"})
    out["cohort"] = cohort

    # join metadata
    out_full = out.merge(md, on=["sample_id", "cohort"], how="left", suffixes=("", "_md"))
    return out_full, cov


def main() -> None:
    all_scores = []
    all_cov = []
    for d in sorted(PROC.iterdir()):
        if not d.is_dir(): continue
        if not (d / "expr_matrix.tsv").exists(): continue
        sc, cov = score_one_cohort(d)
        all_scores.append(sc)
        all_cov.append(cov)
        print()
    sigs = pd.concat(all_scores, ignore_index=True, sort=False)
    sigs.to_csv(OUT / "signature_scores_per_sample.tsv", sep="\t", index=False)
    cov = pd.concat(all_cov, ignore_index=True)
    cov.to_csv(OUT / "module_gene_coverage_per_cohort.tsv", sep="\t", index=False)

    # quick verdict
    print("\n=== overall ===")
    print(f"signature_scores_per_sample.tsv {sigs.shape} → {OUT}")
    print(sigs.groupby("cohort")[["DM1_inflam_composite", "lineage_portable_DM1",
                                   "cytolytic_GZMA_PRF1", "response_binary_CRPR_vs_SD_PD"]]
          .agg(["count", "mean"]).round(3).to_string())


if __name__ == "__main__":
    main()
