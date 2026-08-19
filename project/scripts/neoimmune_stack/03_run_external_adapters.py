#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from neoimmune_common import (
    EXTERNAL_MODELS,
    EXTERNAL_SCORE_COLUMNS,
    REPO_ROOT,
    build_candidate_lookup,
    ensure_run_dir,
    make_key,
    norm_series,
    rank_percentile,
    safe_read_table,
    write_md,
    write_tsv,
    write_tsv_gz,
)


def infer_model(path: Path) -> str:
    stem = path.name.split("__", 1)[0]
    mapping = {
        "MHCflurry": "MHCflurry_2.0_presentation",
        "BigMHC_IM": "BigMHC_IM",
        "BigMHC_EL": "BigMHC_EL",
        "PRIME": "PRIME",
        "NetMHCpan_4.1": "NetMHCpan_4.1_EL",
    }
    return mapping.get(stem, stem)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)
    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    lookup = build_candidate_lookup(canon)

    pred_dir = REPO_ROOT / "project/results/p_neo_bayesian_2026_05_09/wave11/wave11_predictions"
    rows = []
    files = sorted(pred_dir.glob("*.tsv")) if pred_dir.exists() else []
    for p in files:
        model = infer_model(p)
        if model not in EXTERNAL_MODELS and model not in {"DeepImmuno", "TransPHLA", "RF_biophys"}:
            continue
        try:
            df = safe_read_table(p)
        except Exception as e:
            rows.append(
                {
                    "candidate_id": "__adapter_error__",
                    "model_name": model,
                    "score_type": "external_existing_artifact",
                    "raw_score": np.nan,
                    "normalized_score": np.nan,
                    "rank_percentile": np.nan,
                    "higher_is_better": True,
                    "runnable_status": "error_reading_existing_artifact",
                    "error_message": repr(e),
                    "runtime_sec": 0.0,
                }
            )
            continue
        score_cols = [c for c in df.columns if c.startswith("score_")]
        for sc in score_cols:
            higher = not ("aff" in sc.lower() or "rank" in sc.lower() or "nm" in sc.lower())
            tmp = df[["peptide", "hla"] + (["source"] if "source" in df.columns else []) + [sc]].copy()
            tmp["candidate_id"] = [
                lookup.get(make_key(r.get("peptide"), r.get("hla"), r.get("source")), lookup.get(make_key(r.get("peptide"), r.get("hla"), None), ""))
                for _, r in tmp.iterrows()
            ]
            tmp = tmp[tmp["candidate_id"].ne("")]
            if tmp.empty:
                continue
            norm = norm_series(tmp[sc], higher)
            pct = rank_percentile(tmp[sc], higher)
            for cid, raw, nrm, rp in zip(tmp["candidate_id"], tmp[sc], norm, pct):
                rows.append(
                    {
                        "candidate_id": cid,
                        "model_name": model,
                        "score_type": sc,
                        "raw_score": raw,
                        "normalized_score": nrm,
                        "rank_percentile": rp,
                        "higher_is_better": higher,
                        "runnable_status": "existing_artifact_collected",
                        "error_message": "",
                        "runtime_sec": 0.0,
                    }
                )

    present = {r["model_name"] for r in rows}
    status_rows = []
    for m in EXTERNAL_MODELS:
        status = "existing_artifact_collected" if m in present else "pending_adapter_or_license_install"
        status_rows.append({"model_name": m, "status": status, "clean_track_allowed": False, "production_track_allowed": True})
        if m not in present:
            rows.append(
                {
                    "candidate_id": "__adapter_stub__",
                    "model_name": m,
                    "score_type": "pending",
                    "raw_score": np.nan,
                    "normalized_score": np.nan,
                    "rank_percentile": np.nan,
                    "higher_is_better": True,
                    "runnable_status": "pending_adapter_or_license_install",
                    "error_message": "No runnable local binary or compatible existing artifact found in this run.",
                    "runtime_sec": 0.0,
                }
            )
    out = pd.DataFrame(rows, columns=EXTERNAL_SCORE_COLUMNS).drop_duplicates()
    write_tsv_gz(out, outdir / "external_scores" / "all_external_scores.tsv.gz")
    write_tsv(pd.DataFrame(status_rows), outdir / "external_scores" / "external_adapter_status.tsv")
    md = [
        "# External adapter collection",
        "",
        f"- Existing score rows collected: {len(out[out['candidate_id'] != '__adapter_stub__']):,}",
        f"- Existing wave11 files scanned: {len(files):,}",
        "",
        "## Boundary",
        "External model scores are frozen comparators/features for the production stack only. They are not allowed as training features in the clean science track.",
    ]
    write_md("\n".join(md) + "\n", outdir / "external_scores" / "external_adapter_report.md")


if __name__ == "__main__":
    main()

