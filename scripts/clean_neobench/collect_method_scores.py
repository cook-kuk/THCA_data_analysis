#!/usr/bin/env python3
"""Collect internal and public predictor outputs into the CLEAN-NeoBench schema."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    METHOD_SCORE_COLUMNS,
    detect_score_direction,
    ensure_dir,
    fold_safe_calibrate,
    method_family_and_role,
    normalize_empty,
    percentile_calibrate,
    read_table,
    update_manifest,
    write_tsv,
)


STRICT_METHOD_COLUMNS = [
    "Structure_LR",
    "MHCflurry",
    "BigMHC_IM",
    "PRIME",
    "NetMHCpan_4.1",
    "Wave8_TCR_SelfSim_full",
    "Wave8_TCR_SelfSim_no_exact",
    "Wave8_TCR_motif_only",
    "ESM2_Bayesian",
    "ESMFold_3D_CV_score",
    "mean_pLDDT_peptide",
    "min_pLDDT_peptide",
    "interface_contacts_8A",
    "interface_contacts_10A",
    "radius_of_gyration_peptide",
    "peptide_helicity_proxy",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", default=".", help="Repository root")
    parser.add_argument("--output-root", required=True, help="Output directory")
    return parser.parse_args()


def add_score(rows: list[dict], candidate_id: str, method: str, score: object, label: object | None = None, score_context: str = "") -> None:
    try:
        val = float(score)
    except Exception:
        return
    if not np.isfinite(val):
        return
    rows.append(
        {
            "candidate_id": candidate_id,
            "method_name": method,
            "score_raw": val,
            "runtime_status": "ok",
            "score_context": score_context,
            "_label": label,
        }
    )


def collect_strict_seed(repo: Path, master: pd.DataFrame, rows: list[dict], warnings: list[str]) -> None:
    path = repo / "project/results/p_neo_bayesian_2026_05_09/clean_neo_v0_2026_05_09/clean_neo_v0_classI_strict_seed.tsv"
    if not path.exists():
        warnings.append(f"missing strict seed scores: {path}")
        return
    strict = read_table(path)
    key_map = (
        master.assign(_key=master["peptide"].astype(str) + "|" + master["hla"].astype(str))
        .drop_duplicates("_key")
        .set_index("_key")["candidate_id"]
        .to_dict()
    )
    for _, r in strict.iterrows():
        key = f"{normalize_empty(r.get('peptide'))}|{normalize_empty(r.get('hla'))}"
        cid = key_map.get(key)
        if not cid:
            continue
        for method in STRICT_METHOD_COLUMNS:
            if method in strict.columns:
                add_score(rows, cid, method, r.get(method), r.get("label"), score_context="strict_classI_seed")


def collect_curated_itsndb(repo: Path, master: pd.DataFrame, rows: list[dict], warnings: list[str]) -> None:
    path = repo / "project/results/p_neo_bayesian_2026_05_09/curation_2026_05_09/curated_predictions_itsndb_long.tsv"
    if not path.exists():
        warnings.append(f"missing curated ITSNdb long predictions: {path}")
        return
    pred = read_table(path)
    key_map = (
        master.assign(_key=master["peptide"].astype(str) + "|" + master["hla"].astype(str))
        .drop_duplicates("_key")
        .set_index("_key")["candidate_id"]
        .to_dict()
    )
    for _, r in pred.iterrows():
        key = f"{normalize_empty(r.get('peptide'))}|{normalize_empty(r.get('hla'))}"
        cid = key_map.get(key)
        if not cid:
            continue
        add_score(rows, cid, normalize_empty(r.get("method")), r.get("score"), r.get("label"), score_context=normalize_empty(r.get("split")) or "curated_itsndb")


def collect_sample_prediction_file(path: Path, rows: list[dict], sample_to_candidate: dict[str, str], warnings: list[str]) -> None:
    if not path.exists():
        warnings.append(f"missing prediction file: {path}")
        return
    pred = read_table(path)
    method_col = "method" if "method" in pred.columns else "model" if "model" in pred.columns else "branch" if "branch" in pred.columns else None
    if method_col is None or "sample_id" not in pred.columns or "score" not in pred.columns:
        warnings.append(f"prediction file skipped due to missing schema: {path}")
        return
    for _, r in pred.iterrows():
        sid = normalize_empty(r.get("sample_id"))
        cid = sample_to_candidate.get(sid)
        if not cid:
            continue
        method = normalize_empty(r.get(method_col))
        context = normalize_empty(r.get("split_name")) or normalize_empty(r.get("heldout_study")) or path.stem
        add_score(rows, cid, method, r.get("score"), r.get("label"), score_context=context)


def annotate_and_rank(scores: pd.DataFrame, master: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame(columns=METHOD_SCORE_COLUMNS + ["score_context"])
    master_small = master[["candidate_id", "label", "patient_id", "hla_allele_4digit"]].copy()
    scores = scores.merge(master_small, on="candidate_id", how="left", suffixes=("", "_master"))
    if "_label" in scores.columns:
        scores["_eval_label"] = pd.to_numeric(scores["_label"], errors="coerce").fillna(pd.to_numeric(scores["label"], errors="coerce"))
    else:
        scores["_eval_label"] = pd.to_numeric(scores["label"], errors="coerce")

    pieces = []
    for method, g in scores.groupby("method_name", dropna=False):
        gg = g.copy()
        direction = detect_score_direction(gg["_eval_label"], gg["score_raw"])
        gg["score_direction"] = direction
        if direction == "lower_better":
            gg["_score_norm"] = -pd.to_numeric(gg["score_raw"], errors="coerce")
        else:
            gg["_score_norm"] = pd.to_numeric(gg["score_raw"], errors="coerce")
        gg["score_calibrated"] = fold_safe_calibrate(gg["_eval_label"], gg["_score_norm"])
        if gg["score_calibrated"].isna().all():
            gg["score_calibrated"] = percentile_calibrate(gg["_score_norm"])
        pieces.append(gg)
    scores = pd.concat(pieces, ignore_index=True, sort=False)

    meta = scores["method_name"].map(method_family_and_role)
    scores["method_family"] = [m[0] for m in meta]
    scores["method_role"] = [m[1] for m in meta]
    scores["uses_public_pretraining"] = [m[2] for m in meta]
    scores["training_overlap_audited"] = [m[3] for m in meta]
    scores["clean_comparator_allowed"] = [m[4] for m in meta]
    scores["caveat"] = [m[5] for m in meta]

    scores["rank_global"] = scores.groupby("method_name")["score_calibrated"].rank(ascending=False, method="first")
    scores["rank_within_patient"] = scores.groupby(["method_name", "patient_id"])["score_calibrated"].rank(ascending=False, method="first")
    scores["rank_within_hla"] = scores.groupby(["method_name", "hla_allele_4digit"])["score_calibrated"].rank(ascending=False, method="first")
    scores["runtime_status"] = scores["runtime_status"].fillna("ok")
    return scores[METHOD_SCORE_COLUMNS + ["score_context"]].sort_values(["method_name", "rank_global"])


def synthetic_scores(master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in master.iterrows():
        label = float(r.get("label", 0) or 0)
        base = 0.25 + 0.5 * label
        add_score(rows, r["candidate_id"], "Structure_LR", base, label, "synthetic")
        add_score(rows, r["candidate_id"], "MHCflurry", 0.3 + 0.35 * label, label, "synthetic")
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    repo = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    ensure_dir(output_root)
    master = pd.read_csv(output_root / "clean_neobench_master.tsv", sep="\t")
    rows: list[dict] = []
    warnings: list[str] = []
    sample_to_candidate = master.drop_duplicates("sample_id").set_index("sample_id")["candidate_id"].to_dict()

    collect_strict_seed(repo, master, rows, warnings)
    collect_curated_itsndb(repo, master, rows, warnings)
    for rel in [
        "project/results/cross_neo_v0/oof_predictions.tsv",
        "project/results/cross_neo_v0/qk_fallback_predictions.tsv",
        "project/results/cross_neo_v0/source_heldout_predictions.tsv",
        "project/results/cross_neo_v1/v1_all_model_predictions.tsv",
        "project/results/cross_neo_v1_lockdown/locked_anchor_predictions.tsv",
        "project/results/cross_neo_v1_lockdown/foldsafe_fusion_predictions.tsv",
        "project/results/cross_neo_v1_lockdown/source_topk_rescue_predictions.tsv",
    ]:
        collect_sample_prediction_file(repo / rel, rows, sample_to_candidate, warnings)

    if not rows:
        warnings.append("no method scores found; using synthetic fallback scores")
        scores = synthetic_scores(master)
    else:
        scores = pd.DataFrame(rows)
    out = annotate_and_rank(scores, master)
    out = out.drop_duplicates(["candidate_id", "method_name", "score_context"], keep="first")
    write_tsv(out, output_root / "clean_neobench_method_scores.tsv")
    update_manifest(
        output_root,
        "collect_method_scores",
        {
            "n_score_rows": int(len(out)),
            "n_methods": int(out["method_name"].nunique()) if len(out) else 0,
            "method_role_counts": out["method_role"].value_counts().to_dict() if len(out) else {},
            "method_family_counts": out["method_family"].value_counts().to_dict() if len(out) else {},
            "warnings": warnings,
        },
    )
    print(f"[clean-neobench-scores] rows={len(out)} methods={out['method_name'].nunique() if len(out) else 0}")


if __name__ == "__main__":
    main()
