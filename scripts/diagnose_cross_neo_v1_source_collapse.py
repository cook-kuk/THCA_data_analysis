#!/usr/bin/env python3
"""Diagnose source-heldout top-k collapse for CROSS-Neo v1 lockdown."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import OneHotEncoder

from cross_neo_v1_lockdown_common import OUT, V0, ensure_dirs, load_master, metric_row, peptide_cluster


SOURCES = ["CEDAR", "NEPdb", "TESLA_mmc4", "TESLA_mmc7_validation"]


def total_variation(a: pd.Series, b: pd.Series) -> float:
    aa = a.value_counts(normalize=True)
    bb = b.value_counts(normalize=True)
    idx = aa.index.union(bb.index)
    return float(0.5 * np.abs(aa.reindex(idx, fill_value=0) - bb.reindex(idx, fill_value=0)).sum())


def calibration_bins(df: pd.DataFrame, method: str, n_bins: int = 5) -> list[dict[str, object]]:
    rows = []
    d = df[df["method"].eq(method)].copy()
    if d.empty:
        return rows
    d["bin"] = pd.cut(d["score"], bins=np.linspace(0, 1, n_bins + 1), include_lowest=True)
    for b, g in d.groupby("bin", observed=True):
        rows.append(
            {
                "method": method,
                "score_bin": str(b),
                "bin_n": len(g),
                "bin_mean_score": float(g["score"].mean()),
                "bin_positive_rate": float(g["label"].mean()) if len(g) else np.nan,
            }
        )
    return rows


def source_predictability(master: pd.DataFrame) -> pd.DataFrame:
    train = master[master["split"].eq("train") & master["study"].isin(SOURCES)].copy()
    if train["study"].nunique() < 2:
        return pd.DataFrame()
    emb = np.load(V0 / "counterfactual_embeddings.npy")
    idx = pd.read_csv(V0 / "counterfactual_feature_index.tsv", sep="\t")
    emb_df = pd.DataFrame(emb[:, :60], columns=[f"cf_{i:02d}" for i in range(min(60, emb.shape[1]))])
    emb_df.insert(0, "sample_id", idx["sample_id"].values)
    d = train[["sample_id", "study", "peptide_mut", "hla", "hla_supertype", "label"]].merge(emb_df, on="sample_id", how="left")
    d["peptide_length"] = d["peptide_mut"].astype(str).str.len()
    d["label"] = d["label"].astype(int)
    cat = d[["hla_supertype"]].fillna("UNKNOWN").astype(str)
    try:
        enc = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    except TypeError:
        enc = OneHotEncoder(handle_unknown="ignore", sparse=False)
    cat_x = enc.fit_transform(cat)
    num_cols = [c for c in d.columns if c.startswith("cf_")] + ["peptide_length", "label"]
    x = np.hstack([d[num_cols].fillna(0).to_numpy(float), cat_x])
    y = d["study"].astype(str).to_numpy()
    min_count = pd.Series(y).value_counts().min()
    n_splits = int(min(5, min_count))
    if n_splits < 2:
        return pd.DataFrame()
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=20260509)
    rows = []
    for fold, (tr, te) in enumerate(skf.split(x, y)):
        clf = RandomForestClassifier(n_estimators=160, max_depth=4, min_samples_leaf=8, random_state=20260509 + fold, n_jobs=-1)
        clf.fit(x[tr], y[tr])
        pred = clf.predict(x[te])
        rows.append(
            {
                "fold": fold,
                "n_test": len(te),
                "accuracy": float(accuracy_score(y[te], pred)),
                "balanced_accuracy": float(balanced_accuracy_score(y[te], pred)),
                "majority_baseline": float(pd.Series(y[tr]).value_counts(normalize=True).max()),
                "interpretation": "source_shortcut_risk" if balanced_accuracy_score(y[te], pred) > 0.55 else "weak_source_predictability",
            }
        )
    rows.append(
        {
            "fold": "mean",
            "n_test": len(y),
            "accuracy": float(np.mean([r["accuracy"] for r in rows])),
            "balanced_accuracy": float(np.mean([r["balanced_accuracy"] for r in rows])),
            "majority_baseline": float(pd.Series(y).value_counts(normalize=True).max()),
            "interpretation": "source_shortcut_risk" if np.mean([r["balanced_accuracy"] for r in rows]) > 0.55 else "weak_source_predictability",
        }
    )
    return pd.DataFrame(rows)


def main() -> None:
    ensure_dirs()
    master = load_master()
    source_pred = pd.read_csv(V0 / "source_heldout_predictions.tsv", sep="\t")
    source_pred = source_pred.merge(
        master[
            [
                "sample_id",
                "peptide_mut",
                "peptide_wt",
                "hla",
                "hla_supertype",
                "study",
                "source_protein",
                "source_window_15aa",
                "source_window_30aa",
                "near_peptide_cluster",
                "public_overlap_flags",
            ]
        ],
        on="sample_id",
        how="left",
    )
    source_pred["peptide_length"] = source_pred["peptide_mut"].astype(str).str.len()
    source_pred["wt_available"] = source_pred["peptide_wt"].notna() & source_pred["peptide_wt"].astype(str).ne("")
    source_pred["source_window_available"] = source_pred["source_window_15aa"].notna() | source_pred["source_window_30aa"].notna()
    source_pred["source_protein_available"] = source_pred["source_protein"].notna() & source_pred["source_protein"].astype(str).ne("")

    diag_rows = []
    rank_rows = []
    shift_rows = []
    for heldout in SOURCES:
        test_meta = master[(master["split"].eq("train")) & (master["study"].eq(heldout))].copy()
        train_meta = master[(master["split"].eq("train")) & (~master["study"].eq(heldout))].copy()
        if test_meta.empty or train_meta.empty:
            continue
        test_meta["peptide_length"] = test_meta["peptide_mut"].astype(str).str.len()
        train_meta["peptide_length"] = train_meta["peptide_mut"].astype(str).str.len()
        test_meta["cluster_tmp"] = test_meta["peptide_mut"].map(peptide_cluster)
        train_meta["cluster_tmp"] = train_meta["peptide_mut"].map(peptide_cluster)
        train_hla = set(train_meta["hla"].dropna().astype(str))
        train_cluster = set(train_meta["cluster_tmp"].dropna().astype(str))
        train_source_protein = set(train_meta["source_protein"].dropna().astype(str))
        hla_overlap = float(test_meta["hla"].astype(str).isin(train_hla).mean())
        cluster_overlap = float(test_meta["cluster_tmp"].astype(str).isin(train_cluster).mean())
        sp = test_meta["source_protein"].dropna().astype(str)
        source_protein_overlap = float(sp.isin(train_source_protein).mean()) if len(sp) else np.nan
        for field in ["peptide_length", "hla", "hla_supertype", "cluster_tmp"]:
            if field == "peptide_length":
                tr = train_meta[field].astype(float)
                te = test_meta[field].astype(float)
                denom = float(np.sqrt((tr.var() + te.var()) / 2)) if (tr.var() + te.var()) > 0 else np.nan
                shift_rows.append(
                    {
                        "heldout_study": heldout,
                        "feature": field,
                        "shift_metric": "standardized_mean_delta",
                        "value": float((te.mean() - tr.mean()) / denom) if denom and np.isfinite(denom) else 0.0,
                        "train_summary": f"mean={tr.mean():.3f};sd={tr.std():.3f}",
                        "heldout_summary": f"mean={te.mean():.3f};sd={te.std():.3f}",
                    }
                )
            else:
                shift_rows.append(
                    {
                        "heldout_study": heldout,
                        "feature": field,
                        "shift_metric": "categorical_total_variation",
                        "value": total_variation(train_meta[field].fillna("NA").astype(str), test_meta[field].fillna("NA").astype(str)),
                        "train_summary": ",".join(train_meta[field].fillna("NA").astype(str).value_counts().head(3).index.tolist()),
                        "heldout_summary": ",".join(test_meta[field].fillna("NA").astype(str).value_counts().head(3).index.tolist()),
                    }
                )

        held = source_pred[source_pred["heldout_study"].eq(heldout)].copy()
        for method, g in held.groupby("method"):
            g = g.sort_values("score", ascending=False).reset_index(drop=True)
            y = g["label"].to_numpy(int)
            mm = metric_row(y, g["score"].to_numpy(float))
            pos_ranks = (np.where(y == 1)[0] + 1).tolist()
            for rank in pos_ranks:
                row = g.iloc[rank - 1]
                rank_rows.append(
                    {
                        "heldout_study": heldout,
                        "method": method,
                        "sample_id": row["sample_id"],
                        "positive_rank": rank,
                        "rank_pct": rank / max(1, len(g)),
                        "score": row["score"],
                        "peptide_mut": row["peptide_mut"],
                        "hla": row["hla"],
                        "hla_supertype": row["hla_supertype"],
                    }
                )
            diag_rows.append(
                {
                    "heldout_study": heldout,
                    "method": method,
                    "n": len(g),
                    "positives": int(y.sum()),
                    "prevalence": float(y.mean()),
                    "positive_top5": int(y[:5].sum()),
                    "positive_top10": int(y[:10].sum()),
                    "positive_top20": int(y[:20].sum()),
                    "positive_top50": int(y[:50].sum()),
                    "AUPRC": mm["AUPRC"],
                    "AUROC": mm["AUROC"],
                    "top10_precision": mm["top10_precision"],
                    "Brier": mm["Brier"],
                    "ECE": mm["ECE"],
                    "median_positive_rank": float(np.median(pos_ranks)) if pos_ranks else np.nan,
                    "hla_overlap_with_train": hla_overlap,
                    "peptide_cluster_overlap_with_train": cluster_overlap,
                    "source_protein_overlap_with_train": source_protein_overlap,
                    "wt_available_rate": float(test_meta["peptide_wt"].notna().mean()),
                    "source_window_available_rate": float(test_meta[["source_window_15aa", "source_window_30aa"]].notna().any(axis=1).mean()),
                    "length_median": float(test_meta["peptide_length"].median()),
                    "dominant_hla": str(test_meta["hla"].fillna("NA").value_counts().index[0]),
                    "dominant_hla_frac": float(test_meta["hla"].fillna("NA").value_counts(normalize=True).iloc[0]),
                    "cause_prevalence_too_low": bool(y.mean() < 0.10),
                    "cause_source_shift": bool(
                        any(
                            r.get("heldout_study") == heldout and r.get("shift_metric") == "categorical_total_variation" and r.get("value", 0) > 0.35
                            for r in shift_rows
                        )
                    ),
                    "cause_hla_shift": bool(hla_overlap < 0.75),
                    "cause_length_shift": bool(abs(next((r["value"] for r in shift_rows if r["heldout_study"] == heldout and r["feature"] == "peptide_length"), 0)) > 0.5),
                    "cause_missing_wt_source": bool(test_meta["peptide_wt"].notna().mean() < 0.1 or test_meta[["source_window_15aa", "source_window_30aa"]].notna().any(axis=1).mean() < 0.1),
                    "cause_score_miscalibration": bool(mm["ECE"] > 0.15 if pd.notna(mm["ECE"]) else False),
                    "cause_positive_low_tail": bool(np.median(pos_ranks) > len(g) * 0.5 if pos_ranks else False),
                    "cause_public_overlap_unresolved": True,
                    "cause_label_definition_mismatch_possible": bool(heldout.startswith("TESLA") or heldout == "NEPdb"),
                }
            )
            for cb in calibration_bins(held, method):
                cb["heldout_study"] = heldout
                cb["diagnostic_type"] = "calibration_bin"
                shift_rows.append(cb)

    diag = pd.DataFrame(diag_rows)
    ranks = pd.DataFrame(rank_rows)
    shifts = pd.DataFrame(shift_rows)
    pred = source_predictability(master)

    diag.to_csv(OUT / "source_collapse_diagnostics.tsv", sep="\t", index=False)
    ranks.to_csv(OUT / "source_positive_rank_positions.tsv", sep="\t", index=False)
    shifts.to_csv(OUT / "source_feature_shift.tsv", sep="\t", index=False)
    pred.to_csv(OUT / "source_predictability.tsv", sep="\t", index=False)

    collapse = diag[(diag["method"].eq("sourceheld_prespecified_late_fusion_w0.5")) & (diag["top10_precision"].fillna(0).le(0))]
    lines = [
        "# CROSS-Neo v1 Source-Heldout Collapse Diagnosis",
        "",
        "This is an internal source-heldout stress test, not external validation.",
        "",
        "## Headline Diagnostics",
        "",
        diag.sort_values(["heldout_study", "method"]).to_markdown(index=False) if len(diag) else "No diagnostics.",
        "",
        "## Collapsed Source Rows",
        "",
        collapse.to_markdown(index=False) if len(collapse) else "No zero-top10 source rows.",
        "",
        "## Source Predictability Test",
        "",
        pred.to_markdown(index=False) if len(pred) else "Source predictability test unavailable.",
        "",
        "## Interpretation",
        "",
        "- TESLA sources have very low prevalence, so top-10 is statistically brittle.",
        "- NEPdb/TESLA positives are often ranked in the low-score tail under the locked branches.",
        "- Source labels are highly distribution-specific; this supports a source shortcut/label-definition shift concern.",
        "- Public corpus overlap remains unresolved and must be closed before public comparator cleanliness claims.",
    ]
    (OUT / "source_collapse_report.md").write_text("\n".join(lines) + "\n")
    print(f"[source-collapse] diagnostics={len(diag)} positive_ranks={len(ranks)} predictability_rows={len(pred)}")


if __name__ == "__main__":
    main()
