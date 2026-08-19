#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from neoimmune_common import (
    LOCAL_SCORE_COLUMNS,
    REPO_ROOT,
    build_candidate_lookup,
    clean_hla,
    ensure_run_dir,
    make_key,
    norm_series,
    rank_percentile,
    safe_read_table,
    write_md,
    write_tsv_gz,
)


def add_rows_from_score_table(rows: list[dict], path: Path, id_col: str, score_map: dict[str, str], split: str, scope: str) -> None:
    if not path.exists():
        return
    try:
        df = safe_read_table(path)
    except Exception:
        return
    for col, model in score_map.items():
        if col not in df.columns:
            continue
        norm = norm_series(df[col], True)
        for cid, raw, nrm in zip(df[id_col].astype(str), df[col], norm):
            rows.append(
                {
                    "candidate_id": cid,
                    "local_model_name": model,
                    "local_score": raw,
                    "normalized_score": nrm,
                    "fold_id": "",
                    "split_type": split,
                    "training_scope": scope,
                    "notes": f"collected from {path}",
                }
            )


def add_rows_from_kg_ga_mapped(rows: list[dict], path: Path, canon: pd.DataFrame) -> None:
    if not path.exists():
        return
    try:
        ga = safe_read_table(path)
    except Exception:
        return
    score_map = {
        "kg_ga_evolved_score": "KG_GA_evolved_controller",
        "immunogenicity_discovery_score": "KG_GA_fixed_integrated",
        "immunogenicity_claim_safe_score": "KG_GA_fixed_claimsafe",
        "tcr_recognition_score_norm": "KG_GA_TCR_expert_component",
        "md_control_score_norm": "KG_GA_MD_control_component",
        "fitness_foreignness_proxy": "KG_GA_foreignness_component",
        "impact_portfolio_score": "KG_GA_impact_portfolio_component",
    }
    c = canon[["candidate_id", "peptide_mut", "hla_allele", "dataset_source"]].copy()
    c["peptide_key"] = c["peptide_mut"].fillna("").astype(str).str.strip()
    c["hla_key"] = c["hla_allele"].map(clean_hla)
    c["source_key"] = c["dataset_source"].fillna("").astype(str).str.strip()
    g = ga.copy()
    g["peptide_key"] = g["peptide"].fillna("").astype(str).str.strip()
    g["hla_key"] = g["hla_allele_4digit"].map(clean_hla)
    g["source_key"] = g["source_name"].fillna("").astype(str).str.strip()
    mapped = c.merge(g, on=["peptide_key", "hla_key", "source_key"], how="inner", suffixes=("_canon", "_ga"))
    if "candidate_id_canon" in mapped.columns:
        mapped["candidate_id"] = mapped["candidate_id_canon"]
    direct = canon[["candidate_id"]].merge(g, on="candidate_id", how="inner") if "candidate_id" in g.columns else pd.DataFrame()
    for frame, mapping_note in [(direct, "direct_candidate_id"), (mapped, "peptide_hla_source_patient_mapped")]:
        if frame.empty:
            continue
        for col, model in score_map.items():
            if col not in frame.columns:
                continue
            norm = norm_series(frame[col], True)
            for cid, raw, nrm in zip(frame["candidate_id"].astype(str), frame[col], norm):
                rows.append(
                    {
                        "candidate_id": cid,
                        "local_model_name": model,
                        "local_score": raw,
                        "normalized_score": nrm,
                        "fold_id": "",
                        "split_type": "kg_ga_retrospective_architecture_search",
                        "training_scope": "production-only GA/evolutionary controller; contains public predictor and product-value components",
                        "notes": f"{mapping_note} from {path}",
                    }
                )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    lookup = build_candidate_lookup(canon)
    rows: list[dict] = []

    base = REPO_ROOT / "project/results/clean_neobench_barneo_2026_05_09"
    p = base / "clean_neobench_method_scores.tsv"
    if p.exists():
        df = safe_read_table(p)
        for _, r in df.iterrows():
            rows.append(
                {
                    "candidate_id": str(r.get("candidate_id", "")),
                    "local_model_name": str(r.get("method_name", "")),
                    "local_score": r.get("score_raw", np.nan),
                    "normalized_score": r.get("score_calibrated", np.nan),
                    "fold_id": "",
                    "split_type": str(r.get("score_context", "")),
                    "training_scope": str(r.get("method_role", "")),
                    "notes": str(r.get("caveat", "")),
                }
            )

    add_rows_from_score_table(rows, base / "barneo_candidate_scores.tsv", "candidate_id", {"barneo_score": "BAR-Neo", "patient_gated_score": "BAR-Neo_patient_gated", "confidence_score": "BAR-Neo_confidence"}, "production_existing", "patient-gated local score")
    add_rows_from_score_table(rows, base / "barneo_bma_candidate_scores.tsv", "candidate_id", {"bma_score": "BAR-Neo-BMA", "claim_safe_score": "BAR-Neo-BMA_claim_safe", "discovery_score": "BAR-Neo-BMA_discovery"}, "production_existing", "BMA local score")
    add_rows_from_score_table(rows, base / "barneo_x_candidate_scores.tsv", "candidate_id", {"barneo_x_discovery_score": "BAR-Neo-X", "barneo_x_claim_safe_score": "BAR-Neo-X_claim_safe"}, "production_existing", "reviewer-aware local score")
    add_rows_from_score_table(rows, REPO_ROOT / "project/results/cross_neo_kaggle_winner_playbook_2026_05_10/score_finetune_2026_05_10/fine_tuned_candidate_scores.tsv", "candidate_id", {"finetuned_clean_stack_score": "stress_guarded_clean_stack", "stress_guarded_discovery_score": "stress_guarded_discovery", "stress_guarded_claim_safe_score": "stress_guarded_claim_safe", "clean_contextual_bma_score": "clean_contextual_bma"}, "existing_stress_guarded", "existing local/product score")
    add_rows_from_kg_ga_mapped(rows, REPO_ROOT / "project/results/cross_neo_kg_ga_evolutionary_ensemble_2026_05_10/kg_ga_evolved_candidate_scores.tsv", canon)

    pred_dir = REPO_ROOT / "project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions"
    local_prefix = ("Structure_LR", "ESM2_Bayesian", "GP_quantum", "RF_biophys", "DeepImmuno")
    if pred_dir.exists():
        for p in sorted(pred_dir.glob("*.tsv")):
            model = p.name.split("__", 1)[0]
            if not model.startswith(local_prefix):
                continue
            try:
                df = safe_read_table(p)
            except Exception:
                continue
            score_cols = [c for c in df.columns if c.startswith("score_")]
            for sc in score_cols:
                norm = norm_series(df[sc], True)
                for _, r in df.iterrows():
                    cid = lookup.get(make_key(r.get("peptide"), r.get("hla"), r.get("source")), lookup.get(make_key(r.get("peptide"), r.get("hla"), None), ""))
                    if not cid:
                        continue
                    rows.append(
                        {
                            "candidate_id": cid,
                            "local_model_name": model,
                            "local_score": r.get(sc, np.nan),
                            "normalized_score": norm.loc[_],
                            "fold_id": "",
                            "split_type": str(r.get("source", "")),
                            "training_scope": "wave11_existing_prediction",
                            "notes": f"{sc} from {p}",
                        }
                    )

    out = pd.DataFrame(rows, columns=LOCAL_SCORE_COLUMNS).drop_duplicates()
    write_tsv_gz(out, outdir / "local_scores" / "all_local_scores.tsv.gz")
    md = [
        "# Local model output collection",
        "",
        f"- Rows collected: {len(out):,}",
        f"- Unique local model names: {out['local_model_name'].nunique() if not out.empty else 0:,}",
        "",
        "## Boundary",
        "These are existing local artifacts. A clean claim still requires split-aware leakage checks before a model is treated as a scientific contribution.",
    ]
    write_md("\n".join(md) + "\n", outdir / "local_scores" / "local_score_report.md")


if __name__ == "__main__":
    main()
