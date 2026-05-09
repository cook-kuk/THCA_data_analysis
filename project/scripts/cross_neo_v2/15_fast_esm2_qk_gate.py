#!/usr/bin/env python3
"""Fast vectorized ESM2+QK fold-safe gate."""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import OUT, ensure_dirs, metrics


PRIMARY = ["exact_peptide_hla_holdout", "near_peptide_cluster_holdout", "hla_stratified_group_5fold", "hla_supertype_heldout"]
BASES = [
    "prespecified_rf_qk_no_anchor_w0.5",
    "rule_gate_rf_qk_fallback_train_selected",
    "v2_cf_plm_lr",
    "v2_multimodal_lr",
]
ESMS = [
    "v2_cf_esm2_150m_lr_fast",
    "v2_cf_esm2_650m_lr_fast",
    "v2_multimodal_esm2_650m_lr_fast",
]
WEIGHTS = [0.25, 0.5, 0.75, 1.0]


def pct(train: np.ndarray, x: np.ndarray) -> np.ndarray:
    train = np.asarray(train, dtype=float)
    x = np.asarray(x, dtype=float)
    if len(train) < 10:
        return x
    return np.asarray([(train <= v).mean() for v in x], dtype=float)


def load_pred() -> pd.DataFrame:
    parts = [pd.read_csv(OUT / "predictions/all_predictions.tsv", sep="\t")]
    for p in sorted((OUT / "predictions").glob("esm2_*_fast_predictions.tsv")):
        parts.append(pd.read_csv(p, sep="\t"))
    df = pd.concat(parts, ignore_index=True, sort=False)
    keep = BASES + ESMS
    df = df[df["model_name"].isin(keep)].copy()
    df["row_id"] = df["row_id"].astype(str)
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df["label"] = pd.to_numeric(df["label"], errors="coerce").fillna(0).astype(int)
    return df.dropna(subset=["score"])


def main() -> None:
    ensure_dirs()
    pred = load_pred()
    fused_rows = []
    selected = []
    for split in PRIMARY:
        sdf = pd.read_csv(OUT / "splits" / f"{split}.tsv", sep="\t")
        p = pred[pred["split_name"].eq(split)].copy()
        wide = p.pivot_table(index="row_id", columns="model_name", values="score", aggfunc="mean")
        labels = p.groupby("row_id")["label"].first()
        for fold_id, fold in sdf.groupby("fold_id"):
            train_ids = fold.loc[fold["role"].eq("train"), "row_id"].astype(str)
            test_ids = fold.loc[fold["role"].eq("test"), "row_id"].astype(str)
            train_ids = [x for x in train_ids if x in wide.index]
            test_ids = [x for x in test_ids if x in wide.index]
            best = None
            for base in BASES:
                for esm in ESMS:
                    if base not in wide or esm not in wide:
                        continue
                    tr = wide.loc[train_ids, [base, esm]].dropna()
                    if len(tr) < 30:
                        continue
                    ytr = labels.loc[tr.index].astype(int).values
                    if len(np.unique(ytr)) < 2:
                        continue
                    btr = pct(tr[base].values, tr[base].values)
                    etr = pct(tr[esm].values, tr[esm].values)
                    for w in WEIGHTS:
                        sc = w * btr + (1 - w) * etr
                        m = metrics(ytr, sc, ks=(10,))
                        key = (m["AUPRC"], m["top10_precision"], w)
                        if best is None or key > best["key"]:
                            best = {"base": base, "esm": esm, "w": w, "key": key, "train_AUPRC": m["AUPRC"], "train_top10": m["top10_precision"]}
            if best is None:
                continue
            te = wide.loc[test_ids, [best["base"], best["esm"]]].dropna()
            if te.empty:
                continue
            train_ref = wide.loc[train_ids, [best["base"], best["esm"]]].dropna()
            b = pct(train_ref[best["base"]].values, te[best["base"]].values)
            e = pct(train_ref[best["esm"]].values, te[best["esm"]].values)
            score = best["w"] * b + (1 - best["w"]) * e
            for rid, sc in zip(te.index, score):
                fused_rows.append({
                    "split_name": split,
                    "fold_id": fold_id,
                    "row_id": rid,
                    "label": int(labels.loc[rid]),
                    "score": float(sc),
                    "model_name": "fast_nested_esm2_qk_gate",
                    "model_family": "esm2_qk_foldsafe_gate",
                    "claim_status": "reviewer_safe_internal_locked",
                    "selected_base": best["base"],
                    "selected_esm": best["esm"],
                    "base_weight": best["w"],
                    "train_selected_AUPRC": best["train_AUPRC"],
                    "train_selected_top10": best["train_top10"],
                })
            selected.append({"split_name": split, "fold_id": fold_id, **{k: best[k] for k in ["base", "esm", "w", "train_AUPRC", "train_top10"]}})
    fused = pd.DataFrame(fused_rows)
    fused["score"] = fused["score"].round(12)
    fused["rank"] = fused.groupby(["split_name", "fold_id", "model_name"])["score"].rank(method="first", ascending=False)
    fused["rank_pct"] = fused["rank"] / fused.groupby(["split_name", "fold_id", "model_name"])["row_id"].transform("size")
    fused.to_csv(OUT / "predictions/fast_esm2_qk_gate_predictions.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(selected).to_csv(OUT / "metrics/fast_esm2_qk_gate_selected.tsv", sep="\t", index=False, na_rep="NA")

    rows = []
    for split, g in fused.groupby("split_name"):
        rows.append({"split_name": split, "model_name": "fast_nested_esm2_qk_gate", "model_family": "esm2_qk_foldsafe_gate", "claim_status": "reviewer_safe_internal_locked", **metrics(g["label"].values, g["score"].values)})
    met = pd.DataFrame(rows)
    met.to_csv(OUT / "metrics/fast_esm2_qk_gate_metrics.tsv", sep="\t", index=False, na_rep="NA")
    base = pd.read_csv(OUT / "metrics/all_model_all_split_metrics_plus_runpod_esm2.tsv", sep="\t")
    combined = pd.concat([base, met], ignore_index=True, sort=False)
    combined.to_csv(OUT / "metrics/all_model_all_split_metrics_plus_runpod_esm2_fast_gate.tsv", sep="\t", index=False, na_rep="NA")
    comp = []
    for split in PRIMARY:
        b = base[base["split_name"].eq(split)].sort_values(["AUPRC", "top10_precision"], ascending=False).head(1)
        f = met[met["split_name"].eq(split)].head(1)
        c = combined[combined["split_name"].eq(split)].sort_values(["AUPRC", "top10_precision"], ascending=False).head(1)
        comp.append({
            "split_name": split,
            "pre_gate_best": b["model_name"].iloc[0],
            "pre_gate_AUPRC": b["AUPRC"].iloc[0],
            "pre_gate_top10": b["top10_precision"].iloc[0],
            "gate_AUPRC": f["AUPRC"].iloc[0] if not f.empty else np.nan,
            "gate_top10": f["top10_precision"].iloc[0] if not f.empty else np.nan,
            "combined_best": c["model_name"].iloc[0],
            "combined_AUPRC": c["AUPRC"].iloc[0],
            "combined_top10": c["top10_precision"].iloc[0],
        })
    comp = pd.DataFrame(comp)
    comp.to_csv(OUT / "metrics/fast_esm2_qk_gate_primary_comparison.tsv", sep="\t", index=False, na_rep="NA")
    xlsx = OUT / "CROSS_Neo_v2_all_results_summary_plus_runpod_esm2_fast_gate.xlsx"
    with pd.ExcelWriter(xlsx) as writer:
        combined.to_excel(writer, sheet_name="all_plus_fast_gate", index=False)
        comp.to_excel(writer, sheet_name="primary_compare", index=False)
        met.to_excel(writer, sheet_name="gate_metrics", index=False)
        pd.DataFrame(selected).to_excel(writer, sheet_name="selected", index=False)
    report = [
        "# Fast ESM2+QK Gate Report",
        "",
        comp.to_markdown(index=False),
        "",
        "Selection used same-split OOF predictions on outer-train rows only.",
        "",
        "## Selected Experts",
        "",
        pd.DataFrame(selected).to_markdown(index=False),
        "",
        f"Workbook: `{xlsx}`",
    ]
    (OUT / "FAST_ESM2_QK_gate_report.md").write_text("\n".join(report) + "\n")
    print(f"[fast-gate] rows={len(fused)}")
    print(comp.to_string(index=False))


if __name__ == "__main__":
    main()
