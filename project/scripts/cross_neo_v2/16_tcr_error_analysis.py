#!/usr/bin/env python3
"""TCR-aware case/error analysis for the CROSS-Neo-TCR pilot."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from common import OUT, safe_parquet


TCR_OUT = OUT / "tcr_extension"
PILOT_OUT = TCR_OUT / "pilot_subset_compare"


BASELINE_MODEL = "pmhc_counterfactual"
CONSERVATIVE_TCR_MODEL = "pmhc_counterfactual_plus_tcr_exsource_nolabel"
RAW_TCR_MODEL = "pmhc_counterfactual_plus_tcr"


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    view = df.copy()
    view = view.where(view.notna(), "NA").replace("", "NA")
    view.to_csv(path, sep="\t", index=False, na_rep="NA")


def read_linked() -> pd.DataFrame:
    path = TCR_OUT / "tcr_neo_linked_registry.parquet"
    if path.exists():
        return pd.read_parquet(path)
    return pd.read_csv(TCR_OUT / "tcr_neo_linked_registry.tsv", sep="\t", low_memory=False)


def read_registry_counts() -> dict[str, int]:
    path = TCR_OUT / "tcr_registry.parquet"
    if path.exists():
        tcr = pd.read_parquet(
            path,
            columns=[
                "paired_tcr_available",
                "peptide_hla_available",
                "binding_label_binary",
                "cdr3_alpha",
                "cdr3_beta",
            ],
        )
    else:
        tcr = pd.read_csv(TCR_OUT / "tcr_registry.tsv", sep="\t", low_memory=False)
    label = pd.to_numeric(tcr.get("binding_label_binary"), errors="coerce")
    peptide_hla = tcr.get("peptide_hla_available", False).fillna(False).astype(bool)
    paired = tcr.get("paired_tcr_available", False).fillna(False).astype(bool)
    alpha = tcr.get("cdr3_alpha", "").fillna("").astype(str).ne("")
    beta = tcr.get("cdr3_beta", "").fillna("").astype(str).ne("")
    any_tcr = alpha | beta
    return {
        "tcr_registry_rows": int(len(tcr)),
        "paired_alpha_beta_rows": int(paired.sum()),
        "any_tcr_peptide_hla_labeled_rows": int((any_tcr & peptide_hla & label.notna()).sum()),
        "paired_tcr_peptide_hla_labeled_rows": int((paired & peptide_hla & label.notna()).sum()),
    }


def model_comparison(pred: pd.DataFrame, comparison_model: str) -> pd.DataFrame:
    base = pred[pred["model_name"].eq(BASELINE_MODEL)].copy()
    comp = pred[pred["model_name"].eq(comparison_model)].copy()
    keep = ["split_name", "fold_id", "row_id", "label_binary", "tcr_link_category", "tcr_evidence_count", "paired_tcr_evidence_count"]
    base = base[keep + ["score"]].rename(columns={"score": "pmhc_score"})
    comp = comp[["split_name", "fold_id", "row_id", "score"]].rename(columns={"score": "tcr_augmented_score"})
    out = base.merge(comp, on=["split_name", "fold_id", "row_id"], how="inner")
    out["score_delta_tcr_minus_pmhc"] = out["tcr_augmented_score"] - out["pmhc_score"]
    out["comparison_model"] = comparison_model
    return out


def add_rank_columns(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (_, _), g in df.groupby(["split_name", "fold_id"], dropna=False):
        g = g.copy()
        n = max(1, len(g))
        g["pmhc_rank_pct"] = g["pmhc_score"].rank(ascending=False, method="average") / n
        g["tcr_augmented_rank_pct"] = g["tcr_augmented_score"].rank(ascending=False, method="average") / n
        g["rank_delta_pct_positive_is_rescue"] = g["pmhc_rank_pct"] - g["tcr_augmented_rank_pct"]
        rows.append(g)
    return pd.concat(rows, ignore_index=True) if rows else df


def annotate_cases(comp: pd.DataFrame, linked: pd.DataFrame) -> pd.DataFrame:
    annot_cols = [
        "row_id",
        "source_dataset",
        "study_id",
        "assay_type",
        "cancer_type",
        "peptide",
        "wildtype_peptide",
        "hla_4digit",
        "hla_supertype",
        "source_protein",
        "tcr_evidence_sources",
        "beta_only_tcr_evidence_count",
        "peptide_hla_tcr_label_count",
        "cancer_context_evidence_count",
        "pathogen_context_evidence_count",
        "tcr_binding_positive_count",
        "tcr_binding_negative_count",
        "structure_evidence_count",
        "example_tcr_registry_rows",
        "example_tcr_cdr3b",
        "example_tcr_cdr3a",
        "example_tcr_antigen_sources",
        "tcr_link_evidence_basis",
        "tcr_link_claim_status",
        "tcr_link_warning",
    ]
    annot_cols = [c for c in annot_cols if c in linked.columns]
    return comp.merge(linked[annot_cols], on="row_id", how="left")


def case_tables(comp: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    label = comp["label_binary"].astype(int)
    rescue = comp[
        label.eq(1)
        & (
            ((comp["pmhc_score"] < 0.35) & (comp["tcr_augmented_score"] >= 0.50))
            | (comp["score_delta_tcr_minus_pmhc"] >= 0.25)
            | (comp["rank_delta_pct_positive_is_rescue"] >= 0.30)
        )
    ].copy()
    rescue["case_reason"] = "positive_rescued_by_tcr_branch"

    harm = comp[
        label.eq(1)
        & (
            ((comp["pmhc_score"] >= 0.50) & (comp["tcr_augmented_score"] < 0.35))
            | (comp["score_delta_tcr_minus_pmhc"] <= -0.25)
            | (comp["rank_delta_pct_positive_is_rescue"] <= -0.30)
        )
    ].copy()
    harm["case_reason"] = "positive_harmed_by_tcr_branch"

    fp = comp[
        label.eq(0)
        & (
            (comp["pmhc_score"] >= 0.75)
            | (comp["tcr_augmented_score"] >= 0.75)
            | (comp["tcr_augmented_rank_pct"] <= 0.10)
        )
    ].copy()
    fp["case_reason"] = np.select(
        [
            (fp["pmhc_score"] >= 0.75) & (fp["tcr_augmented_score"] < 0.35),
            (fp["pmhc_score"] >= 0.75) & (fp["tcr_augmented_score"] >= 0.75),
            (fp["pmhc_score"] < 0.50) & (fp["tcr_augmented_score"] >= 0.75),
        ],
        ["main_high_tcr_low", "main_high_tcr_high", "tcr_high_only"],
        default="top_ranked_negative",
    )

    fn = comp[
        label.eq(1)
        & (
            ((comp["pmhc_score"] < 0.35) & (comp["tcr_augmented_score"] < 0.35))
            | ((comp["pmhc_rank_pct"] > 0.80) & (comp["tcr_augmented_rank_pct"] > 0.80))
        )
    ].copy()
    fn["case_reason"] = "positive_missed_by_both"

    rescue = rescue.sort_values(["score_delta_tcr_minus_pmhc", "rank_delta_pct_positive_is_rescue"], ascending=False)
    harm = harm.sort_values(["score_delta_tcr_minus_pmhc", "rank_delta_pct_positive_is_rescue"], ascending=True)
    fp = fp.sort_values(["tcr_augmented_score", "pmhc_score"], ascending=False)
    fn = fn.sort_values(["tcr_augmented_score", "pmhc_score"], ascending=True)
    return rescue, harm, fp, fn


def aggregate_candidate_table(comp: pd.DataFrame) -> pd.DataFrame:
    score = (
        comp.groupby("row_id", dropna=False)
        .agg(
            label_binary=("label_binary", "max"),
            n_eval=("row_id", "size"),
            pmhc_score_mean=("pmhc_score", "mean"),
            tcr_augmented_score_mean=("tcr_augmented_score", "mean"),
            score_delta_mean=("score_delta_tcr_minus_pmhc", "mean"),
            rank_delta_mean=("rank_delta_pct_positive_is_rescue", "mean"),
            best_tcr_augmented_score=("tcr_augmented_score", "max"),
            split_names=("split_name", lambda s: "|".join(sorted(set(map(str, s))))),
        )
        .reset_index()
    )
    score["wetlab_priority_score"] = (
        score["tcr_augmented_score_mean"].fillna(0)
        + 0.35 * score["score_delta_mean"].fillna(0)
        + 0.20 * score["rank_delta_mean"].fillna(0)
    )
    return score.sort_values("wetlab_priority_score", ascending=False)


def write_case_report(
    conservative: pd.DataFrame,
    raw: pd.DataFrame,
    rescue: pd.DataFrame,
    harm: pd.DataFrame,
    fp: pd.DataFrame,
    fn: pd.DataFrame,
    linked: pd.DataFrame,
    counts: dict[str, int],
) -> None:
    linked_any = int(linked["tcr_link_category"].ne("no_tcr_match").sum())
    linked_exact = int(linked["tcr_link_category"].eq("exact_tcr_pmhc_match").sum())
    linked_cancer = int(pd.to_numeric(linked.get("cancer_context_evidence_count", 0), errors="coerce").fillna(0).gt(0).sum())
    split_counts = (
        conservative.groupby("split_name")
        .agg(
            n=("row_id", "size"),
            positives=("label_binary", "sum"),
            mean_pmhc=("pmhc_score", "mean"),
            mean_tcr_aug=("tcr_augmented_score", "mean"),
            mean_delta=("score_delta_tcr_minus_pmhc", "mean"),
        )
        .reset_index()
    )
    lines = [
        "# CROSS-Neo-TCR Case/Error Analysis",
        "",
        "This is a diagnostic pilot over the current subset-comparison predictions. The conservative comparison uses source-excluded, no-public-label TCR evidence features.",
        "",
        "## Registry And Linkage Counts",
        "",
        f"- TCR registry rows: {counts['tcr_registry_rows']:,}",
        f"- Paired alpha/beta TCR rows: {counts['paired_alpha_beta_rows']:,}",
        f"- Any-TCR peptide-HLA labeled rows: {counts['any_tcr_peptide_hla_labeled_rows']:,}",
        f"- Paired-TCR peptide-HLA labeled rows: {counts['paired_tcr_peptide_hla_labeled_rows']:,}",
        f"- CROSS-Neo rows with any TCR-resource overlap: {linked_any:,} / {len(linked):,}",
        f"- CROSS-Neo rows with exact paired TCR-pMHC overlap: {linked_exact:,} / {len(linked):,}",
        f"- CROSS-Neo rows with cancer-context TCR evidence: {linked_cancer:,} / {len(linked):,}",
        "",
        "## Conservative Case Counts",
        "",
        f"- TCR rescue cases: {len(rescue):,}",
        f"- TCR harm cases: {len(harm):,}",
        f"- False-positive audit rows: {len(fp):,}",
        f"- False-negative audit rows: {len(fn):,}",
        "",
        "## Split-Level Mean Score Shift",
        "",
        "| Split | n | positives | mean pMHC | mean TCR-aug | mean delta |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in split_counts.itertuples(index=False):
        lines.append(f"| {r.split_name} | {int(r.n)} | {int(r.positives)} | {r.mean_pmhc:.3f} | {r.mean_tcr_aug:.3f} | {r.mean_delta:.3f} |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Raw TCR evidence remains useful as a diagnostic annotation, but can be inflated by source self-evidence.",
            "- The conservative source-excluded/no-label branch is the safer readout for claim discussions.",
            "- Handcrafted sequence/TCR-evidence features do not establish de novo TCR recognition. A true TCR expert still needs ERGO-II/NetTCR/pMTnet/PLM or structure-backed validation.",
        ]
    )
    (TCR_OUT / "tcr_case_report.md").write_text("\n".join(lines) + "\n")


def write_decision_report(linked: pd.DataFrame, counts: dict[str, int], aggregate: pd.DataFrame, deltas: pd.DataFrame) -> None:
    linked_any = int(linked["tcr_link_category"].ne("no_tcr_match").sum())
    linked_exact = int(linked["tcr_link_category"].eq("exact_tcr_pmhc_match").sum())
    linked_cancer = int(pd.to_numeric(linked.get("cancer_context_evidence_count", 0), errors="coerce").fillna(0).gt(0).sum())

    def delta_line(split: str, comparison: str = CONSERVATIVE_TCR_MODEL) -> str:
        sub = deltas[(deltas["split_name"].eq(split)) & (deltas["comparison_model"].eq(comparison))]
        if sub.empty:
            return "not available"
        r = sub.iloc[0]
        return f"AUPRC {r['baseline_auprc']:.3f} -> {r['comparison_auprc']:.3f} (delta {r['delta_auprc']:+.3f}), AUROC delta {r['delta_auroc']:+.3f}"

    decoy = aggregate[aggregate["pilot"].eq("paired_tcr_decoy_sequence_pilot")].copy()
    decoy_best = "not available"
    if not decoy.empty:
        best = decoy.sort_values("auprc_mean", ascending=False).iloc[0]
        decoy_best = f"{best['model_name']} AUPRC {best['auprc_mean']:.3f}, AUROC {best['auroc_mean']:.3f}"

    lines = [
        "# CROSS-Neo-TCR Extension Decision Report",
        "",
        "Decision: **PROMOTE_TO_DIAGNOSTIC_CASE_STUDY** and **PROMOTE_TO_WETLAB_PRIORITIZATION_TOOL** for the current sprint. Hold main-method promotion until paired TCR benchmarks, stronger TCR sequence models, and structure validation are complete.",
        "",
        "## Required Answers",
        "",
        f"1. Paired TCR alpha/beta rows: **{counts['paired_alpha_beta_rows']:,}** in the local TCR registry.",
        f"2. Peptide-HLA-TCR labeled rows: **{counts['any_tcr_peptide_hla_labeled_rows']:,}** with any TCR sequence; **{counts['paired_tcr_peptide_hla_labeled_rows']:,}** with paired alpha/beta.",
        f"3. Cancer neoantigen overlap: **{linked_any:,}/{len(linked):,}** CROSS-Neo rows have some TCR-resource overlap; **{linked_exact:,}** exact paired TCR-pMHC rows; **{linked_cancer:,}** rows have cancer-context evidence.",
        f"4. Top-k/ranking improvement: internal/exact/near pilots improve with TCR evidence. Conservative exact holdout: {delta_line('exact_peptide_hla_holdout')}; near holdout: {delta_line('near_peptide_cluster_holdout')}.",
        f"5. Structure over sequence-only: **not yet supported**. Structure files currently contain missingness/QC features unless PDB/template evidence exists.",
        f"6. Mutant-WT interface delta: **not yet supported**. Needs parsed mutant/WT TCR-pMHC structures.",
        f"7. Source-heldout rescue: partial. Conservative NEPdb source-heldout readout: {delta_line('source_heldout_NEPdb')}. Raw TCR evidence is stronger but treated as leakage-prone.",
        "8. False positives: case audit files identify main-high/TCR-low and TCR-high-only negatives, but explanation remains diagnostic until structures or external TCR assays support it.",
        "9. Main claim or supplement: **supplement/diagnostic branch now**, not the main CROSS-Neo ranking claim.",
        "10. Wetlab candidates: prioritize rows in `tcr_wetlab_candidate_prioritization.tsv` with high conservative TCR-augmented score, positive delta, exact paired or cancer-context TCR evidence, and low pathogen-only dependence.",
        "",
        "## Decoy Recognition Pilot",
        "",
        f"- Paired positive versus shuffled-TCR decoy pilot best result: {decoy_best}. This indicates the simple handcrafted TCR sequence features are not sufficient for de novo cognate-recognition prediction.",
        "",
        "## Claim Boundary",
        "",
        "- Do not replace the main pMHC neoantigen model with the TCR branch.",
        "- Do not transfer pathogen epitope labels to cancer neoantigens.",
        "- Do not claim clinical utility or universal TCR-aware prediction.",
        "- Safe current claim: optional TCR evidence/diagnostic layer plus wetlab prioritization scaffold.",
    ]
    (OUT / "CROSS_Neo_TCR_extension_decision_report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    pred = pd.read_parquet(PILOT_OUT / "cross_neo_tcr_evidence_pilot_predictions.parquet")
    linked = read_linked()
    counts = read_registry_counts()
    aggregate = pd.read_csv(PILOT_OUT / "quick_tcr_subset_compare_aggregate.tsv", sep="\t")
    deltas = pd.read_csv(PILOT_OUT / "quick_tcr_subset_compare_key_deltas.tsv", sep="\t")

    conservative = add_rank_columns(model_comparison(pred, CONSERVATIVE_TCR_MODEL))
    raw = add_rank_columns(model_comparison(pred, RAW_TCR_MODEL))
    conservative = annotate_cases(conservative, linked)
    raw = annotate_cases(raw, linked)

    rescue, harm, fp, fn = case_tables(conservative)
    candidates = aggregate_candidate_table(conservative)
    candidates = annotate_cases(candidates, linked)

    out_cols_first = [
        "split_name",
        "fold_id",
        "row_id",
        "case_reason",
        "label_binary",
        "pmhc_score",
        "tcr_augmented_score",
        "score_delta_tcr_minus_pmhc",
        "pmhc_rank_pct",
        "tcr_augmented_rank_pct",
        "rank_delta_pct_positive_is_rescue",
        "comparison_model",
        "source_dataset",
        "peptide",
        "wildtype_peptide",
        "hla_4digit",
        "source_protein",
        "tcr_link_category",
        "tcr_evidence_count",
        "paired_tcr_evidence_count",
        "cancer_context_evidence_count",
        "pathogen_context_evidence_count",
        "structure_evidence_count",
        "example_tcr_cdr3b",
        "example_tcr_cdr3a",
        "tcr_link_claim_status",
        "tcr_link_warning",
    ]

    def write_table(df: pd.DataFrame, name: str, limit: int | None = None) -> None:
        cols = [c for c in out_cols_first if c in df.columns] + [c for c in df.columns if c not in out_cols_first]
        view = df[cols].copy()
        if limit is not None:
            view = view.head(limit)
        write_tsv(view, TCR_OUT / name)

    write_table(rescue, "tcr_rescue_cases.tsv")
    write_table(harm, "tcr_harm_cases.tsv")
    write_table(fp, "tcr_false_positive_audit.tsv", limit=5000)
    write_table(fn, "tcr_false_negative_audit.tsv", limit=5000)
    write_tsv(candidates, TCR_OUT / "tcr_wetlab_candidate_prioritization.tsv")
    safe_parquet(conservative, TCR_OUT / "tcr_error_analysis_conservative_predictions.parquet")
    safe_parquet(raw, TCR_OUT / "tcr_error_analysis_raw_predictions.parquet")

    write_case_report(conservative, raw, rescue, harm, fp, fn, linked, counts)
    write_decision_report(linked, counts, aggregate, deltas)
    print(f"[tcr-errors] rescue={len(rescue)} harm={len(harm)} fp={len(fp)} fn={len(fn)} candidates={len(candidates)}")


if __name__ == "__main__":
    main()
