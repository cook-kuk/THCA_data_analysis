#!/usr/bin/env python3
from __future__ import annotations

import argparse
import numpy as np
import pandas as pd

from neoimmune_common import (
    REPO_ROOT,
    clean_hla,
    ensure_run_dir,
    metric_binary,
    patient_topn_metrics,
    precision_recall_at_k,
    safe_read_table,
    write_md,
    write_tsv,
)


GA_PATH = REPO_ROOT / "project/results/cross_neo_kg_ga_evolutionary_ensemble_2026_05_10/kg_ga_evolved_candidate_scores.tsv"
GA_BENCH = REPO_ROOT / "project/results/cross_neo_kg_ga_evolutionary_ensemble_2026_05_10/kg_ga_benchmark.tsv"
GA_WEIGHTS = REPO_ROOT / "project/results/cross_neo_kg_ga_evolutionary_ensemble_2026_05_10/kg_ga_best_feature_weights.tsv"
GA_GATES = REPO_ROOT / "project/results/cross_neo_kg_ga_evolutionary_ensemble_2026_05_10/kg_ga_best_gates_and_synergies.tsv"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    if not GA_PATH.exists():
        raise SystemExit(f"Missing GA evolved score file: {GA_PATH}")

    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    ga = safe_read_table(GA_PATH)
    c = canon.copy()
    c["peptide_key"] = c["peptide_mut"].fillna("").astype(str).str.strip()
    c["hla_key"] = c["hla_allele"].map(clean_hla)
    c["source_key"] = c["dataset_source"].fillna("").astype(str).str.strip()
    g = ga.copy()
    g["peptide_key"] = g["peptide"].fillna("").astype(str).str.strip()
    g["hla_key"] = g["hla_allele_4digit"].map(clean_hla)
    g["source_key"] = g["source_name"].fillna("").astype(str).str.strip()
    direct = c.merge(g, on="candidate_id", how="inner", suffixes=("", "_ga")) if "candidate_id" in g.columns else pd.DataFrame()
    direct["kg_ga_mapping_mode"] = "direct_candidate_id"
    mapped = c.merge(g, on=["peptide_key", "hla_key", "source_key"], how="inner", suffixes=("", "_ga"))
    if "candidate_id_ga" in mapped.columns:
        mapped["kg_ga_source_candidate_id"] = mapped["candidate_id_ga"]
    mapped["kg_ga_mapping_mode"] = "peptide_hla_source_patient_mapped"
    d = pd.concat([direct, mapped], ignore_index=True, sort=False)
    d = d.drop_duplicates(subset=["candidate_id", "kg_ga_evolved_score"], keep="first")
    d["label_immunogenicity"] = pd.to_numeric(d["label_immunogenicity"], errors="coerce")
    for c in [
        "kg_ga_evolved_score",
        "immunogenicity_discovery_score",
        "immunogenicity_claim_safe_score",
        "bigmhc_im_score",
        "bigmhc_el_score",
        "tcr_recognition_score_norm",
        "md_control_score_norm",
        "fitness_foreignness_proxy",
    ]:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")

    score_map = {
        "KG_GA_evolved_controller": "kg_ga_evolved_score",
        "KG_GA_fixed_integrated": "immunogenicity_discovery_score",
        "KG_GA_fixed_claimsafe": "immunogenicity_claim_safe_score",
        "BigMHC_IM_inside_GA": "bigmhc_im_score",
        "BigMHC_EL_inside_GA": "bigmhc_el_score",
        "GA_TCR_expert_component": "tcr_recognition_score_norm",
        "GA_MD_control_component": "md_control_score_norm",
        "GA_foreignness_component": "fitness_foreignness_proxy",
    }
    rows = []
    eval_df = d[d["label_immunogenicity"].notna()].copy()
    for name, col in score_map.items():
        if col not in eval_df.columns or eval_df[col].notna().sum() == 0:
            continue
        m = metric_binary(eval_df["label_immunogenicity"], eval_df[col])
        for k in [10, 20, 34, 50, 96]:
            m.update(precision_recall_at_k(eval_df, col, "label_immunogenicity", k))
        for n in [20, 34]:
            m.update(patient_topn_metrics(eval_df, col, "label_immunogenicity", n))
        m.update(
            {
                "algorithm": name,
                "score_column": col,
                "track": "production_stack_track",
                "clean_track_allowed": False,
                "reason_clean_disallowed": "GA controller includes public predictor/product-value components; use as production/experiment-priority branch only.",
            }
        )
        rows.append(m)
    metrics = pd.DataFrame(rows)
    order = [
        "track",
        "algorithm",
        "score_column",
        "clean_track_allowed",
        "n",
        "positives",
        "AUROC",
        "AUPRC",
        "Precision@10",
        "Precision@20",
        "Precision@34",
        "Precision@50",
        "Recall@20",
        "Recall@34",
        "Recall@96",
        "patient_hit_rate@20",
        "patient_recall@20",
        "patient_hit_rate@34",
        "patient_recall@34",
        "patients_evaluated",
        "reason_clean_disallowed",
    ]
    for c in order:
        if c not in metrics:
            metrics[c] = np.nan
    metrics = metrics[order + [c for c in metrics.columns if c not in order]].sort_values(["AUPRC", "patient_hit_rate@34"], ascending=False)
    write_tsv(metrics, outdir / "metrics" / "ga_rl_algorithm_metrics.tsv")

    d["kg_ga_rank_global"] = d["kg_ga_evolved_score"].rank(method="first", ascending=False)
    d["kg_ga_rank_within_patient"] = d.groupby("patient_id")["kg_ga_evolved_score"].rank(method="first", ascending=False)
    keep = [
        "candidate_id",
        "patient_id",
        "dataset_source",
        "source_name",
        "peptide_mut",
        "peptide",
        "peptide_wt",
        "hla_allele",
        "hla_allele_4digit",
        "label_immunogenicity",
        "kg_ga_evolved_score",
        "kg_ga_rank_global",
        "kg_ga_rank_within_patient",
        "immunogenicity_action",
        "immunogenicity_claim_blockers",
        "bigmhc_im_score",
        "bigmhc_el_score",
        "tcr_recognition_score_norm",
        "md_control_score_norm",
        "fitness_foreignness_proxy",
        "leakage_risk_level",
    ]
    for c in keep:
        if c not in d.columns:
            d[c] = np.nan
    keep.append("kg_ga_mapping_mode")
    write_tsv(d[keep].sort_values("kg_ga_evolved_score", ascending=False), outdir / "predictions" / "kg_ga_evolved_ranked_candidates.tsv")
    write_tsv(
        d[d["kg_ga_rank_within_patient"] <= 34][keep].sort_values(["patient_id", "kg_ga_rank_within_patient"]),
        outdir / "predictions" / "patient_top34_candidates_kg_ga.tsv",
    )
    write_tsv(
        d[d["kg_ga_rank_within_patient"] <= 20][keep].sort_values(["patient_id", "kg_ga_rank_within_patient"]),
        outdir / "predictions" / "patient_top20_candidates_kg_ga.tsv",
    )

    if GA_WEIGHTS.exists():
        weights = safe_read_table(GA_WEIGHTS)
        write_tsv(weights, outdir / "reports" / "kg_ga_feature_weights.tsv")
    else:
        weights = pd.DataFrame()
    if GA_GATES.exists():
        gates = safe_read_table(GA_GATES)
        write_tsv(gates, outdir / "reports" / "kg_ga_gates_and_synergies.tsv")
    else:
        gates = pd.DataFrame()
    if GA_BENCH.exists():
        bench = safe_read_table(GA_BENCH)
        write_tsv(bench, outdir / "metrics" / "kg_ga_original_benchmark.tsv")
    else:
        bench = pd.DataFrame()

    best = metrics.head(1).to_dict("records")
    best_line = ""
    if best:
        r = best[0]
        best_line = f"- Integrated GA/RL best row: {r['algorithm']} AUPRC={r['AUPRC']:.3f}, AUROC={r['AUROC']:.3f}, patient_hit_rate@34={r.get('patient_hit_rate@34', np.nan):.3f}."
    md = [
        "# GA/RL Algorithm Applied to NeoImmune-Stack",
        "",
        "## What was used",
        f"- GA score file: `{GA_PATH}`",
        "- Algorithm: `KG_GA_evolved_controller` from knowledge-graph genetic/evolutionary architecture search.",
        "- Related RL artifact: Darwin/RL blueprint remains a method-search blueprint, not a trained scoring model.",
        "",
        "## Boundary",
        "- This branch is **production/experiment-priority only**.",
        "- It is not allowed in the clean science track because the evolved controller includes BigMHC/public predictor and impact/product components.",
        "- Claim-safe wording: GA/RL found a retrospective candidate-prioritization controller ready for prospective assay validation.",
        "",
        "## Result",
        best_line,
        "",
        "## Feature weights",
        weights.head(12).to_markdown(index=False) if not weights.empty else "_No weights found._",
        "",
        "## Gates and synergies",
        gates.head(16).to_markdown(index=False) if not gates.empty else "_No gates found._",
        "",
        "## Outputs",
        "- `metrics/ga_rl_algorithm_metrics.tsv`",
        "- `predictions/kg_ga_evolved_ranked_candidates.tsv`",
        "- `predictions/patient_top20_candidates_kg_ga.tsv`",
        "- `predictions/patient_top34_candidates_kg_ga.tsv`",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "ga_rl_algorithm_applied_report.md")
    print(outdir / "reports" / "ga_rl_algorithm_applied_report.md")


if __name__ == "__main__":
    main()
