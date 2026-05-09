#!/usr/bin/env python3
"""K2 full-transcriptome kallisto -> gene TPM matrix -> HT-13/FA-12/MAPK-9 panel scores.

Aggregates transcript-level abundance.tsv per sample to gene-level TPM by summing
within-gene transcript TPMs (standard tximport "no" length-correction; equivalent
for TPM since per-transcript TPM is already length-normalized). Within-K2 z-scores
per gene, then panel = mean z across panel members.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

QUANT_ROOT = Path("/data/thca/v17_korean_k2_full_quant")
TX2GENE = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_braf_nature_sprint_2026_05_09/h_k2_full_quant/tx2gene_v44.tsv"
)
OUT = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/"
    "p2_braf_nature_sprint_2026_05_09/h_k2_full_quant"
)
META = Path("/data/thca/repo_results/v17_korean/K1A_prjeb11591_runs.tsv")
TCGA_MASTER = Path(
    "/data/thca/repo_results/audit_2026_04_30/round9/r9_1_master_with_rai.tsv"
)
EXISTING_PRED = Path("/data/thca/repo_results/v17_korean/K2_korean_predictions_v4.tsv")

# ---- Panel definitions ------------------------------------------------------
HT13 = [
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
    "CD79A", "CD79B", "MS4A1", "AICDA", "CXCL13", "CCR6", "IFNG",
]
# FA-12: fatty-acid metabolism axis (down in DM1 per Round 5 GSEA hallmark).
# Take the canonical Hallmark FATTY_ACID_METABOLISM core 12.
FA12 = [
    "ACADL", "ACADM", "ACAA2", "HADHA", "HADHB", "ECHS1",
    "CPT1A", "CPT2", "ACOX1", "HMGCS2", "PPARA", "ACAT1",
]
# MAPK-9: canonical RAF-MEK-ERK output 9-gene (DUSPs + ETV4/5 + SPRY + immediate-early).
MAPK9 = [
    "DUSP4", "DUSP5", "DUSP6", "ETV4", "ETV5",
    "SPRY2", "SPRY4", "PHLDA1", "FOS",
]
PANELS = {"HT_13": HT13, "FA_12": FA12, "MAPK_9": MAPK9}


def read_one(sample_dir: Path) -> pd.Series | None:
    p = sample_dir / "abundance.tsv"
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t", usecols=["target_id", "tpm"])
    df = df.rename(columns={"tpm": sample_dir.name})
    df["target_id"] = df["target_id"].str.split("|").str[0]
    df = df.set_index("target_id")
    return df[sample_dir.name]


def main() -> None:
    samples = sorted([d for d in QUANT_ROOT.iterdir()
                      if d.is_dir() and (d / "abundance.tsv").exists()])
    print(f"[1/5] Loading {len(samples)} sample abundance.tsv ...")
    series = []
    for s in samples:
        ser = read_one(s)
        if ser is not None:
            series.append(ser)
    tx_tpm = pd.concat(series, axis=1)
    print(f"   tx-level TPM matrix: {tx_tpm.shape}")

    print("[2/5] Aggregating to gene-symbol TPM (sum per gene) ...")
    tx2g = pd.read_csv(TX2GENE, sep="\t",
                       names=["tx", "ensg", "symbol"])
    tx2g["tx"] = tx2g["tx"].astype(str)
    sym_map = tx2g.set_index("tx")["symbol"]
    tx_tpm.index = tx_tpm.index.astype(str)
    aligned = tx_tpm.join(sym_map.rename("symbol"), how="inner")
    gene_tpm = aligned.groupby("symbol").sum(numeric_only=True)
    print(f"   gene-level TPM matrix: {gene_tpm.shape}")
    gene_tpm_path = OUT / "k2_full_gene_tpm.tsv"
    gene_tpm.to_csv(gene_tpm_path, sep="\t")
    print(f"   wrote {gene_tpm_path}")

    print("[3/5] Within-K2 log2(TPM+1) z-score per gene ...")
    log_tpm = np.log2(gene_tpm + 1.0)
    z = log_tpm.sub(log_tpm.mean(axis=1), axis=0).div(
        log_tpm.std(axis=1).replace(0, np.nan), axis=0
    )

    print("[4/5] Panel scores (mean z across panel members)...")
    panel_records = []
    panel_coverage = {}
    for name, genes in PANELS.items():
        present = [g for g in genes if g in z.index]
        missing = [g for g in genes if g not in z.index]
        panel_coverage[name] = {"present": len(present), "missing": missing}
        if not present:
            continue
        panel_z = z.loc[present].mean(axis=0)
        for samp, v in panel_z.items():
            panel_records.append({"sample": samp, "panel": name, "z_score": v})
    panel_df = pd.DataFrame(panel_records)
    panel_wide = panel_df.pivot(index="sample", columns="panel", values="z_score")

    # Annotate with metadata
    meta = pd.read_csv(META, sep="\t")
    meta_subset = meta[["run_accession", "sample_alias", "sample_title"]].set_index(
        "run_accession"
    )
    panel_wide = panel_wide.join(meta_subset, how="left")
    # Tissue type from sample_alias suffix: "-N" = matched normal, no suffix = tumor
    panel_wide["tissue"] = np.where(
        panel_wide["sample_alias"].fillna("").str.endswith("-N"),
        "normal_matched", "tumor",
    )
    # Existing 8-gene DM call
    if EXISTING_PRED.exists():
        pred = pd.read_csv(EXISTING_PRED, sep="\t").set_index("run")
        panel_wide = panel_wide.join(pred[["DM_call", "p_DM2"]], how="left")
    panel_path = OUT / "k2_full_panel_scores.tsv"
    panel_wide.to_csv(panel_path, sep="\t")
    print(f"   wrote {panel_path}")

    print("[5/5] Bimodality + tumor-vs-normal contrast (no DM1/DM2 RNA labels for K2)...")
    auc_records = []
    from scipy import stats
    try:
        from sklearn.metrics import roc_auc_score
        sklearn_ok = True
    except Exception:
        sklearn_ok = False

    contrasts = []
    if (panel_wide["tissue"] == "normal_matched").any():
        contrasts.append(("tumor_vs_normal", panel_wide["tissue"] == "tumor"))
    if "DM_call" in panel_wide.columns and panel_wide["DM_call"].notna().any():
        # 8gene-derived DM_call exists from prior pipeline; treat DM1 as positive class
        dm_mask = panel_wide["DM_call"].isin(["DM1", "DM2"])
        if (panel_wide.loc[dm_mask, "DM_call"] == "DM1").any():
            contrasts.append(
                ("DM1_vs_DM2_8gene_label",
                 panel_wide.loc[dm_mask, "DM_call"] == "DM1")
            )

    for panel_name in PANELS.keys():
        if panel_name not in panel_wide.columns:
            continue
        x_all = panel_wide[panel_name].dropna()
        # GMM BIC bimodality
        try:
            from sklearn.mixture import GaussianMixture
            xx = x_all.values.reshape(-1, 1)
            g1 = GaussianMixture(n_components=1, random_state=0).fit(xx)
            g2 = GaussianMixture(n_components=2, random_state=0).fit(xx)
            bic_delta = g2.bic(xx) - g1.bic(xx)
        except Exception as e:
            bic_delta = np.nan
        # Hartigan dip via diptest if available
        try:
            from diptest import diptest
            _, dip_p = diptest(x_all.values)
        except Exception:
            dip_p = np.nan
        auc_records.append({
            "panel": panel_name,
            "contrast": "GMM_BIC_2_minus_1",
            "n": int(x_all.size),
            "value": bic_delta,
            "interpretation": "negative supports bimodal",
        })
        auc_records.append({
            "panel": panel_name,
            "contrast": "Hartigan_dip_p",
            "n": int(x_all.size),
            "value": dip_p,
            "interpretation": "p<0.05 supports non-unimodal",
        })

        for cname, mask in contrasts:
            sub = panel_wide.loc[mask.index] if hasattr(mask, "index") else panel_wide
            sub_y = mask.astype(int)
            sub_x = sub[panel_name]
            valid = sub_x.notna() & sub_y.reindex(sub_x.index).notna()
            if valid.sum() < 4 or sub_y.reindex(sub_x.index)[valid].nunique() < 2:
                continue
            if sklearn_ok:
                auc = roc_auc_score(sub_y.reindex(sub_x.index)[valid].values,
                                    sub_x[valid].values)
            else:
                auc = np.nan
            # Mann-Whitney U
            pos = sub_x[valid][sub_y.reindex(sub_x.index)[valid] == 1].values
            neg = sub_x[valid][sub_y.reindex(sub_x.index)[valid] == 0].values
            u, mwp = stats.mannwhitneyu(pos, neg, alternative="two-sided")
            auc_records.append({
                "panel": panel_name,
                "contrast": cname,
                "n": int(valid.sum()),
                "n_pos": int((sub_y.reindex(sub_x.index)[valid] == 1).sum()),
                "auc": auc,
                "mw_u": u,
                "mw_p": mwp,
            })

    auc_df = pd.DataFrame(auc_records)
    auc_path = OUT / "k2_panel_aucs.tsv"
    auc_df.to_csv(auc_path, sep="\t", index=False)
    print(f"   wrote {auc_path}")

    summary = {
        "n_samples_quanted": int(len(samples)),
        "n_genes_aggregated": int(gene_tpm.shape[0]),
        "panel_coverage": {
            k: {"present": v["present"], "missing": v["missing"]}
            for k, v in panel_coverage.items()
        },
        "n_tumor": int((panel_wide["tissue"] == "tumor").sum()),
        "n_normal": int((panel_wide["tissue"] == "normal_matched").sum()),
        "kallisto_index": "/data/thca/reference_kallisto/gencode.v44.kallisto.idx.NEW",
        "kallisto_version": "0.51.1",
    }
    with (OUT / "k2_full_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print(json.dumps(summary, indent=2, default=str))
    print("DONE")


if __name__ == "__main__":
    main()
