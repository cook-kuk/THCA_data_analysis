#!/usr/bin/env python3
"""Source-heldout collapse and rescue analysis for CROSS-Neo 2.0."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import OUT, ensure_dirs


SOURCES = ["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"]


def positive_ranks(pred: pd.DataFrame, split: str, model: str) -> pd.DataFrame:
    g = pred[(pred["split_name"].eq(split)) & (pred["model_name"].eq(model))].copy()
    if g.empty:
        return pd.DataFrame()
    g = g.sort_values(["fold_id", "score"], ascending=[True, False])
    g["rank"] = g.groupby("fold_id")["score"].rank(method="first", ascending=False)
    g["rank_pct"] = g["rank"] / g.groupby("fold_id")["row_id"].transform("size")
    return g[g["label"].astype(int).eq(1)].copy()


def hla_overlap(reg: pd.DataFrame, split_df: pd.DataFrame) -> float:
    train = set(split_df.loc[split_df["role"].eq("train"), "row_id"].astype(str))
    test = set(split_df.loc[split_df["role"].eq("test"), "row_id"].astype(str))
    tr = set(reg.loc[reg["row_id"].isin(train), "hla_4digit"].dropna().astype(str))
    te = reg.loc[reg["row_id"].isin(test), "hla_4digit"].dropna().astype(str)
    return float(te.isin(tr).mean()) if len(te) else np.nan


def source_shift_table(reg: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for src in SOURCES:
        split_file = OUT / "splits" / f"source_heldout_{src}.tsv"
        if not split_file.exists():
            continue
        sp = pd.read_csv(split_file, sep="\t")
        train = set(sp.loc[sp["role"].eq("train"), "row_id"].astype(str))
        test = set(sp.loc[sp["role"].eq("test"), "row_id"].astype(str))
        tr = reg[reg["row_id"].isin(train)]
        te = reg[reg["row_id"].isin(test)]
        train_proteins = set(tr["source_protein"].dropna().astype(str))
        test_protein_overlap = te["source_protein"].dropna().astype(str).isin(train_proteins).mean() if te["source_protein"].notna().any() else np.nan
        rows.append({
            "heldout_source": src,
            "n_test": len(te),
            "n_pos": int(te["label_binary"].sum()),
            "prevalence": float(te["label_binary"].mean()) if len(te) else np.nan,
            "test_mean_peptide_length": float(te["peptide_length"].mean()) if len(te) else np.nan,
            "train_mean_peptide_length": float(tr["peptide_length"].mean()) if len(tr) else np.nan,
            "length_shift": float(te["peptide_length"].mean() - tr["peptide_length"].mean()) if len(te) and len(tr) else np.nan,
            "wt_available_test": float(te["wildtype_peptide"].notna().mean()) if len(te) else np.nan,
            "wt_available_train": float(tr["wildtype_peptide"].notna().mean()) if len(tr) else np.nan,
            "source_window_available_test": float(te["source_window"].notna().mean()) if len(te) else np.nan,
            "source_window_available_train": float(tr["source_window"].notna().mean()) if len(tr) else np.nan,
            "hla_overlap_with_train": hla_overlap(reg, sp),
            "source_protein_overlap_with_train": float(test_protein_overlap) if pd.notna(test_protein_overlap) else np.nan,
        })
    return pd.DataFrame(rows)


def main() -> None:
    ensure_dirs()
    pred = pd.read_csv(OUT / "predictions/all_predictions.tsv", sep="\t")
    metrics = pd.read_csv(OUT / "metrics/source_heldout_metrics.tsv", sep="\t")
    reg = pd.read_csv(OUT / "canonical_registry.tsv", sep="\t")
    reg["label_binary"] = pd.to_numeric(reg["label_binary"], errors="coerce").fillna(0).astype(int)
    pred["label"] = pd.to_numeric(pred["label"], errors="coerce").fillna(0).astype(int)

    rescue_rows = []
    rank_rows = []
    for src in SOURCES:
        split = f"source_heldout_{src}"
        sub = metrics[metrics["split_name"].eq(split)].copy()
        if sub.empty:
            continue
        sub["source"] = src
        sub["nonzero_top10"] = sub["top10_precision"].fillna(0).gt(0)
        sub["nonzero_top20"] = sub["top20_precision"].fillna(0).gt(0)
        rescue_rows.extend(sub.to_dict("records"))
        best_models = sub.sort_values(["top10_precision", "top20_precision", "AUPRC"], ascending=False).head(6)["model_name"].tolist()
        for model in ["anchor_rf", "anchor_lr"] + best_models:
            pr = positive_ranks(pred, split, model)
            if not pr.empty:
                pr["source"] = src
                pr["model_name"] = model
                rank_rows.append(pr)
    rescue = pd.DataFrame(rescue_rows)
    ranks = pd.concat(rank_rows, ignore_index=True) if rank_rows else pd.DataFrame()
    shift = source_shift_table(reg)
    rescue.to_csv(OUT / "metrics/source_collapse_rescue_summary.tsv", sep="\t", index=False, na_rep="NA")
    ranks.to_csv(OUT / "metrics/source_positive_rank_positions_v2.tsv", sep="\t", index=False, na_rep="NA")
    shift.to_csv(OUT / "metrics/source_shift_summary.tsv", sep="\t", index=False, na_rep="NA")

    fig_dir = OUT / "figures"
    fig_dir.mkdir(exist_ok=True)
    if not shift.empty:
        mat_cols = ["prevalence", "length_shift", "wt_available_test", "source_window_available_test", "hla_overlap_with_train", "source_protein_overlap_with_train"]
        mat = shift.set_index("heldout_source")[mat_cols].astype(float)
        norm = (mat - mat.mean()) / mat.std(ddof=0).replace(0, 1)
        plt.figure(figsize=(9, 4.5))
        plt.imshow(norm.fillna(0), aspect="auto", cmap="coolwarm", vmin=-2, vmax=2)
        plt.xticks(range(len(mat_cols)), mat_cols, rotation=35, ha="right", fontsize=8)
        plt.yticks(range(len(norm.index)), norm.index)
        plt.colorbar(label="z-score shift")
        plt.title("Source-heldout internal feature shift (descriptive small-n)")
        plt.tight_layout()
        plt.savefig(fig_dir / "source_shift_heatmap.png", dpi=220)
        plt.savefig(fig_dir / "source_shift_heatmap.pdf")
        plt.close()

    if not ranks.empty:
        plt.figure(figsize=(9, 4.5))
        for i, (src, g) in enumerate(ranks.groupby("source")):
            vals = g["rank_pct"].astype(float).values
            x = np.full(len(vals), i) + np.linspace(-0.18, 0.18, max(1, len(vals)))
            plt.scatter(x, vals, s=18, alpha=0.75, label=src)
        plt.axhline(0.1, color="black", linestyle="--", lw=1, label="top10%")
        plt.xticks(range(len(SOURCES)), SOURCES, rotation=20, ha="right")
        plt.ylabel("positive rank percentile; lower is better")
        plt.title("Source-heldout internal positive rank positions")
        plt.tight_layout()
        plt.savefig(fig_dir / "positive_rank_shift_plots.png", dpi=220)
        plt.savefig(fig_dir / "positive_rank_shift_plots.pdf")
        plt.close()

    curve_rows = []
    source_pred = pred[pred["split_name"].astype(str).str.startswith("source_heldout_")].copy()
    for (split, model), g in source_pred.groupby(["split_name", "model_name"]):
        for cov in [1.0, 0.8, 0.6, 0.4, 0.25]:
            n_keep = max(1, int(np.ceil(len(g) * cov)))
            kept = g.sort_values("score", ascending=False).head(n_keep)
            top = kept.head(min(20, len(kept)))
            curve_rows.append({
                "split_name": split,
                "model_name": model,
                "coverage": cov,
                "top20_precision_after_score_abstention": float(top["label"].mean()) if len(top) else np.nan,
                "n_kept": len(kept),
            })
    curve = pd.DataFrame(curve_rows)
    curve.to_csv(OUT / "metrics/source_abstention_or_rescue_curve.tsv", sep="\t", index=False, na_rep="NA")
    if not curve.empty:
        show = curve[curve["model_name"].isin(["anchor_rf", "v2_source_balanced_cf_rf", "v2_rule_moe_cf_plm_structure", "v2_cf_plm_rf"])]
        plt.figure(figsize=(8.5, 4.5))
        for model, g in show.groupby("model_name"):
            plt.plot(g.groupby("coverage")["top20_precision_after_score_abstention"].mean().sort_index().index,
                     g.groupby("coverage")["top20_precision_after_score_abstention"].mean().sort_index().values,
                     marker="o", label=model)
        plt.xlabel("coverage retained")
        plt.ylabel("mean source-heldout top20 precision")
        plt.title("Source-heldout internal abstention/rescue curve")
        plt.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(fig_dir / "abstention_or_rescue_curve.png", dpi=220)
        plt.savefig(fig_dir / "abstention_or_rescue_curve.pdf")
        plt.close()

    lines = ["# Source Collapse Rescue Report", ""]
    if not rescue.empty:
        best = rescue.sort_values(["source", "top10_precision", "top20_precision", "AUPRC"], ascending=[True, False, False, False]).groupby("source").head(3)
        lines += ["## Best Source-Heldout Rows", "", best[["source", "model_name", "n", "n_pos", "prevalence", "AUPRC", "top10_precision", "top20_precision", "claim_status"]].to_markdown(index=False), ""]
        nep = rescue[rescue["source"].eq("NEPdb")]
        tesla = rescue[rescue["source"].astype(str).str.contains("TESLA")]
        lines.append(f"- NEPdb nonzero top10 models: {int(nep['nonzero_top10'].sum()) if not nep.empty else 0}")
        lines.append(f"- TESLA nonzero top20 models: {int(tesla['nonzero_top20'].sum()) if not tesla.empty else 0}")
    lines += [
        "",
        "## Interpretation",
        "",
        "- This remains an internal/source-heldout stress test, not external validation.",
        "- Imported v1 diagnostic rows can have different source-heldout cardinalities from v2 split files; compare source rows descriptively unless the n/test definition matches.",
        "- Nonzero top-k recovery is necessary but not sufficient for SOTA; public overlap remains a claim blocker unless resolved.",
        "- If TESLA positives remain in the low-score tail across model families, the result should be framed as source-shift benchmark evidence rather than hidden.",
    ]
    (OUT / "source_collapse_rescue_report.md").write_text("\n".join(lines) + "\n")
    print(f"[v2-source] rescue_rows={len(rescue)} rank_rows={len(ranks)}")


if __name__ == "__main__":
    main()
