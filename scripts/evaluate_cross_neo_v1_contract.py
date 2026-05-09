#!/usr/bin/env python3
"""CROSS-Neo v1 evaluation contract."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cross_neo_v1_common import V0, V1, bootstrap_metric_ci, ece_score, ensure_v1_dirs, metrics, read_v0_expert_predictions


def collect_predictions() -> pd.DataFrame:
    rows = []
    v0 = read_v0_expert_predictions()
    v0_keep = {
        "C_counterfactual_rf": "v0_C_counterfactual_RF",
        "A_structure_baseline_rf": "v0_Structure_baseline_RF",
        "qk_no_anchor_gamma1": "v0_QK_no_anchor_fallback",
        "qk_quantum_only_gamma1": "v0_QK_quantum_only_fallback",
        "late_prespecified_equal_weight_Cw0.5_qk_no_anchor_gamma1": "v0_fixed_late_fusion_C_QK_no_anchor",
        "late_prespecified_equal_weight_Cw0.5_qk_quantum_only_gamma1": "v0_fixed_late_fusion_C_QK_quantum",
    }
    for exp, model in v0_keep.items():
        sub = v0[v0["expert"] == exp].copy()
        sub["model"] = model
        sub["family"] = "v0_reference"
        rows.append(sub[["split_name", "fold_id", "sample_id", "label", "score", "model", "family"]])
    for path, family in [
        (V1 / "gated_moe_predictions.tsv", "v1_gated_moe"),
        (V1 / "pu_ranking_predictions.tsv", "v1_pu_ranking"),
        (V1 / "source_balanced_predictions.tsv", "v1_source_bias_corrected"),
        (V1 / "decoy_focal_predictions.tsv", "v1_decoy_focal"),
    ]:
        if path.exists():
            p = pd.read_csv(path, sep="\t")
            if "heldout_study" in p.columns:
                p["split_name"] = "source_heldout_" + p["heldout_study"].astype(str)
                p["fold_id"] = p["heldout_study"].astype(str)
            p["family"] = family
            rows.append(p[["split_name", "fold_id", "sample_id", "label", "score", "model", "family"]])
    return pd.concat(rows, ignore_index=True)


def metric_table(pred: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (family, model, split), sub in pred.groupby(["family", "model", "split_name"]):
        mm = metrics(sub["label"].to_numpy(), sub["score"].to_numpy())
        mm["ECE"] = ece_score(sub["label"].to_numpy(), sub["score"].to_numpy())
        rows.append({"family": family, "model": model, "split_name": split, **mm})
    return pd.DataFrame(rows)


def subset_metrics(pred: pd.DataFrame) -> pd.DataFrame:
    out = []
    retr_path = V0 / "retrieval_features_by_fold.tsv"
    if retr_path.exists():
        retr = pd.read_csv(retr_path, sep="\t")
        m = pred.merge(retr[["split_name", "fold_id", "sample_id", "retrieval_leakage_risk"]], on=["split_name", "fold_id", "sample_id"], how="inner")
        clean = m[m["retrieval_leakage_risk"] == "clean_no_reference"]
        for (family, model, split), sub in clean.groupby(["family", "model", "split_name"]):
            out.append({"subset": "retrieval_clean_only", "family": family, "model": model, "split_name": split, **metrics(sub["label"].to_numpy(), sub["score"].to_numpy())})
    audit_path = V1 / "public_overlap_audit_v1.tsv"
    if audit_path.exists():
        aud = pd.read_csv(audit_path, sep="\t")
        clean_ids = set(aud.groupby("sample_id")[["exact_peptide_hla", "exact_peptide", "near_hit_ge_0p75"]].sum().query("exact_peptide_hla == 0 and exact_peptide == 0 and near_hit_ge_0p75 == 0").index)
        clean = pred[pred["sample_id"].isin(clean_ids)]
        for (family, model, split), sub in clean.groupby(["family", "model", "split_name"]):
            out.append({"subset": "public_overlap_local_clean", "family": family, "model": model, "split_name": split, **metrics(sub["label"].to_numpy(), sub["score"].to_numpy())})
    return pd.DataFrame(out)


def bootstrap_table(pred: pd.DataFrame, comp: pd.DataFrame) -> pd.DataFrame:
    rows = []
    top = comp.sort_values(["AUPRC", "top10_precision"], ascending=False).head(60)
    for _, r in top.iterrows():
        sub = pred[(pred["family"] == r["family"]) & (pred["model"] == r["model"]) & (pred["split_name"] == r["split_name"])]
        lo, hi = bootstrap_metric_ci(sub["label"].to_numpy(), sub["score"].to_numpy(), "AUPRC", 250)
        tlo, thi = bootstrap_metric_ci(sub["label"].to_numpy(), sub["score"].to_numpy(), "top10_precision", 250)
        rows.append({"family": r["family"], "model": r["model"], "split_name": r["split_name"], "AUPRC_ci_low": lo, "AUPRC_ci_high": hi, "top10_ci_low": tlo, "top10_ci_high": thi})
    return pd.DataFrame(rows)


def permutation_check(pred: pd.DataFrame, comp: pd.DataFrame) -> pd.DataFrame:
    eligible = comp[(comp["family"].str.startswith("v1")) & (comp["split_name"] == "hla_stratified_group_5fold")].sort_values("AUPRC", ascending=False)
    if eligible.empty:
        return pd.DataFrame()
    r = eligible.iloc[0]
    sub = pred[(pred["family"] == r["family"]) & (pred["model"] == r["model"]) & (pred["split_name"] == r["split_name"])]
    y = sub["label"].to_numpy(int)
    s = sub["score"].to_numpy(float)
    obs = metrics(y, s)["AUPRC"]
    rng = np.random.default_rng(20260509)
    null = []
    for _ in range(300):
        yy = rng.permutation(y)
        if len(np.unique(yy)) > 1:
            null.append(metrics(yy, s)["AUPRC"])
    p = float((np.sum(np.asarray(null) >= obs) + 1) / (len(null) + 1)) if null else np.nan
    return pd.DataFrame([{"family": r["family"], "model": r["model"], "split_name": r["split_name"], "observed_AUPRC": obs, "permutation_p_ge_observed": p, "n_permutations": len(null)}])


def main() -> None:
    ensure_v1_dirs()
    pred = collect_predictions()
    pred.to_csv(V1 / "v1_all_model_predictions.tsv", sep="\t", index=False)
    comp = metric_table(pred).sort_values(["split_name", "AUPRC"], ascending=[True, False])
    subset = subset_metrics(pred)
    comp_full = pd.concat([comp.assign(subset="all_rows"), subset], ignore_index=True, sort=False)
    comp_full.to_csv(V1 / "v1_model_comparison.tsv", sep="\t", index=False)
    boot = bootstrap_table(pred, comp)
    boot.to_csv(V1 / "v1_model_comparison_bootstrap.tsv", sep="\t", index=False)
    perm = permutation_check(pred, comp)
    perm.to_csv(V1 / "v1_permutation_check.tsv", sep="\t", index=False)
    failures = []
    for split in ["source_heldout_NEPdb", "source_heldout_TESLA_mmc4", "source_heldout_TESLA_mmc7_validation"]:
        sub = comp[comp["split_name"] == split].sort_values("AUPRC", ascending=False)
        if len(sub):
            best = sub.iloc[0]
            failures.append({"failure_mode": split, "best_model": best["model"], "best_AUPRC": best["AUPRC"], "best_top10": best["top10_precision"], "prevalence": best["prevalence"], "interpretation": "source-heldout stress; not external validation"})
    if (V1 / "public_overlap_download_manifest_needed.tsv").exists():
        manifest = pd.read_csv(V1 / "public_overlap_download_manifest_needed.tsv", sep="\t")
        unresolved = int(manifest["status"].str.contains("needs|not_confirmed", case=False, regex=True).sum())
        failures.append({"failure_mode": "public_predictor_training_overlap", "best_model": "", "best_AUPRC": np.nan, "best_top10": np.nan, "prevalence": np.nan, "interpretation": f"{unresolved} public comparator training corpora still unresolved"})
    pd.DataFrame(failures).to_csv(V1 / "v1_failure_modes.tsv", sep="\t", index=False)
    lines = [
        "# CROSS-Neo v1 Evaluation Contract",
        "",
        "This compares v0 references and v1 candidates on locked/internal and source-heldout stress splits. It does not claim external validation.",
        "",
        "## Top HLA-Stratified Rows",
        comp_full[(comp_full["split_name"] == "hla_stratified_group_5fold") & (comp_full["subset"] == "all_rows")].sort_values("AUPRC", ascending=False).head(15).to_markdown(index=False),
        "",
        "## Top Retrieval-Clean Rows",
        subset[subset["subset"] == "retrieval_clean_only"].sort_values("AUPRC", ascending=False).head(12).to_markdown(index=False) if len(subset) else "Not available.",
        "",
        "## Source-Heldout Failures",
        pd.DataFrame(failures).to_markdown(index=False) if failures else "No source-heldout rows.",
        "",
        "## Permutation Check",
        perm.to_markdown(index=False) if len(perm) else "Not available.",
    ]
    (V1 / "v1_evaluation_contract.md").write_text("\n".join(lines) + "\n")
    print(f"[v1-eval] comparisons={len(comp_full)} bootstrap={len(boot)}")


if __name__ == "__main__":
    main()
