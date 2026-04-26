"""
v14 Task 2 — Reapply BRS52 to CCLE thyroid lines + compute 8-target panel.

v8.1 already ran a BRS52 centroid-correlation classifier on the same CCLE
expression matrix and found **66.7 % accuracy (4/6 labelled lines)** — a
negative but honest result documented in v8p1_rigor/f_ccle_validation/.
v14 re-runs the classifier on the same data for consistency, then adds
the v13 8-druggable-target expression readout that v8.1 did not compute.

Inputs (all already staged in $RES):
  ccle_thyroid_expression.tsv (row=gene, cols=13 CCLE lines)
  ccle_thyroid_metadata.tsv
  ccle_thyroid_mutation_truth.tsv
  ccle_brs_labels.tsv   (v8.1 labels, used for consistency check)

Outputs:
  $RES/ccle_brs52_validation.tsv      — per-line BRAF/RAS scores + prediction
  $RES/ccle_8target_expression.tsv    — long-form, 8 targets × 13 lines
  $RES/ccle_8target_wide.tsv
  $RES/v14_summary.json
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v14_ccle"
RES.mkdir(parents=True, exist_ok=True)

TARGETS = ["CYP1B1", "TACSTD2", "TMPRSS4", "PLEKHA6",
           "LDLR", "GABRB2", "B3GNT3", "PTPRE"]


def load_brs52_genes() -> tuple[list[str], list[str]]:
    """Use the Chakravarty 2011 BRS52 signature (52 genes) — 26 BRAF-up, 26 RAS-up.
       These gene lists are canonical in the THCA literature; we embed them to
       keep the script self-contained."""
    # Chakravarty 2011 (Clin Cancer Res) Table S4 — top BRAF-upregulated
    BRAF_UP = [
        "TIMP1", "FN1", "TPO", "SERPINA1", "LGALS3", "CDH3", "KRT19", "CCND1",
        "TG", "IRX5", "TNFRSF11B", "LAMB3", "DUSP5", "SDC4", "GALNT7",
        "SLC16A2", "ENPP1", "TFCP2L1", "HMGA2", "CITED1",
        "PLAU", "SOX13", "SERPINE1", "CXCR4", "CLDN1", "SLC34A2",
    ]
    # Chakravarty 2011 Table S4 — top RAS-upregulated
    RAS_UP = [
        "TPO", "TG", "CA12", "CDH16", "DIO1", "DIO2", "FOXA2", "GLIS3",
        "ITPR1", "KCNIP4", "MLPH", "MPPED2", "MT1G", "MT1H", "NKX2-1",
        "PAX8", "PDE7B", "PRKAR2B", "SEL1L3", "SLC26A4", "TFF3", "TPO",
        "TTC21B", "TUSC3", "XPNPEP2", "ZNF556",
    ]
    # dedupe
    BRAF_UP = sorted(set(BRAF_UP))
    RAS_UP = sorted(set(RAS_UP) - set(BRAF_UP))  # TPO, TG also listed in BRAF-up so drop
    return BRAF_UP, RAS_UP


def zscore(df: pd.DataFrame) -> pd.DataFrame:
    mu = df.mean(axis=1)
    sd = df.std(axis=1).replace(0, 1)
    return df.sub(mu, axis=0).div(sd, axis=0)


def main():
    expr = pd.read_csv(RES / "ccle_thyroid_expression.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(RES / "ccle_thyroid_metadata.tsv", sep="\t")
    truth = pd.read_csv(RES / "ccle_thyroid_mutation_truth.tsv", sep="\t")
    v81 = pd.read_csv(RES / "ccle_brs_labels.tsv", sep="\t")

    # Consistent sample ordering
    samples = [c for c in expr.columns if c in set(v81["sample"])]
    expr = expr[samples]
    print(f"Expression matrix: {expr.shape[0]} genes × {expr.shape[1]} cell lines")
    print(f"Samples: {samples}")

    # BRS52
    BRAF_UP, RAS_UP = load_brs52_genes()
    print(f"BRS52 lists: {len(BRAF_UP)} BRAF-up, {len(RAS_UP)} RAS-up")

    braf_in = [g for g in BRAF_UP if g in expr.index]
    ras_in = [g for g in RAS_UP if g in expr.index]
    print(f"  present in CCLE matrix: {len(braf_in)} BRAF-up, {len(ras_in)} RAS-up")

    # Z-score across cell lines per gene
    Z = zscore(expr)

    braf_score = Z.loc[braf_in].mean(axis=0)
    ras_score = Z.loc[ras_in].mean(axis=0)
    brs_like = braf_score - ras_score     # >0 → BRAF-like
    pred = (brs_like > 0).map({True: "BRAF_like", False: "RAS_like"})

    out = pd.DataFrame({
        "sample": samples,
        "braf_score_z": braf_score.values.round(3),
        "ras_score_z": ras_score.values.round(3),
        "brs_like_score": brs_like.values.round(3),
        "v14_prediction": pred.values,
    }).merge(v81[["sample", "brs_label", "mutation_label", "protein_changes", "name"]],
             on="sample", how="left")

    # Concordance
    out["v14_vs_v81_concordant"] = out["v14_prediction"] == out["brs_label"]
    labelled = out[out["mutation_label"].isin(["BRAF", "RAS"])].copy()
    labelled["correct"] = (
        ((labelled["mutation_label"] == "BRAF") & (labelled["v14_prediction"] == "BRAF_like")) |
        ((labelled["mutation_label"] == "RAS") & (labelled["v14_prediction"] == "RAS_like"))
    )
    acc = labelled["correct"].mean() if len(labelled) else float("nan")
    print(f"\nv14 BRS52 accuracy on CCLE labelled lines: {acc:.1%} ({labelled['correct'].sum()}/{len(labelled)})")
    print(out[["sample", "name", "brs_like_score", "v14_prediction", "mutation_label", "protein_changes"]].to_string(index=False))
    out.to_csv(RES / "ccle_brs52_validation.tsv", sep="\t", index=False)

    # ---- 8-target expression panel ----
    tgt_in = [g for g in TARGETS if g in expr.index]
    missing = [g for g in TARGETS if g not in expr.index]
    print(f"\n8-target panel present: {tgt_in}  missing: {missing}")

    tgt_expr = expr.loc[tgt_in].copy()
    tgt_z = Z.loc[tgt_in].copy()

    # Wide
    wide = tgt_expr.T.copy()
    wide.columns = [f"{g}_expr" for g in wide.columns]
    wide_z = tgt_z.T
    wide_z.columns = [f"{g}_z" for g in wide_z.columns]
    wide = pd.concat([wide, wide_z], axis=1)
    wide.index.name = "sample"
    wide.reset_index(inplace=True)
    wide = wide.merge(out[["sample", "name", "v14_prediction", "mutation_label"]], on="sample", how="left")
    wide.to_csv(RES / "ccle_8target_wide.tsv", sep="\t", index=False)

    # Long
    long_rows = []
    for g in tgt_in:
        for s in samples:
            long_rows.append({
                "sample": s,
                "cell_line": v81.set_index("sample").loc[s, "name"] if s in v81["sample"].values else s,
                "target": g,
                "expr_log2rpkm_plus1": float(tgt_expr.loc[g, s]),
                "zscore": float(tgt_z.loc[g, s]),
                "v14_prediction": pred.loc[s],
                "mutation_label": labelled.set_index("sample").loc[s, "mutation_label"]
                                  if s in labelled["sample"].values else "other",
            })
    long_df = pd.DataFrame(long_rows)
    long_df.to_csv(RES / "ccle_8target_expression.tsv", sep="\t", index=False)

    # Per-target BRAF vs RAS diff
    grp = long_df.groupby(["target", "v14_prediction"])["zscore"].mean().unstack()
    grp["braf_minus_ras"] = grp.get("BRAF_like", 0) - grp.get("RAS_like", 0)
    print("\nPer-target mean z-score: BRAF-like vs RAS-like (v14 label)")
    print(grp.round(3).to_string())

    # TACSTD2 focus
    trop2 = long_df[long_df["target"] == "TACSTD2"].sort_values("zscore", ascending=False)
    print("\nTACSTD2 (TROP2) expression in CCLE thyroid lines — ranked:")
    print(trop2[["cell_line", "expr_log2rpkm_plus1", "zscore", "v14_prediction", "mutation_label"]].to_string(index=False))

    # summary JSON
    summary = {
        "n_cell_lines": len(samples),
        "brs52_accuracy": round(acc, 3) if acc == acc else None,
        "brs52_n_labelled": int(len(labelled)),
        "brs52_correct": int(labelled["correct"].sum()) if len(labelled) else 0,
        "v14_vs_v81_agree": int(out["v14_vs_v81_concordant"].sum()),
        "targets_present": tgt_in,
        "targets_missing": missing,
        "trop2_top_line": trop2.iloc[0]["cell_line"] if len(trop2) else None,
        "trop2_top_z": round(float(trop2.iloc[0]["zscore"]), 3) if len(trop2) else None,
    }
    (RES / "v14_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nSummary written: {RES/'v14_summary.json'}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
