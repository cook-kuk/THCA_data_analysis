#!/usr/bin/env python3
"""cBioPortal sweep — Landa 2016 + GATCI 2024 + TCGA-PTC GDC.

Goal: extend Paper 1 trajectory PT(TCGA) -> PTC(TCGA) -> PDTC(Landa) -> ATC(Landa+GATCI)
with formal 8-gene scoring + DM1/DM2 calls + survival where available.

Output (results/cbioportal_sweep/):
  panel_expression_<study>.tsv   per-sample × 8-gene log2 expression
  dm_calls_<study>.tsv           per-sample DM1/DM2 + RAI_8 + DM1_like
  trajectory_summary.tsv         pooled stage-wise summary (PT/PTC/PDTC/ATC)
  pooled_meta.json               key numbers for manuscript drop-in
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import (  # noqa: E402
    PANEL_8,
    PANEL_8_ENTREZ,
    kmeans_dm,
    log,
    resolve_panel,
    score_panel,
    set_threads,
    write_status,
)

set_threads(2)

WORKER = "cbioportal_sweep"
OUT = Path(__file__).parent
API = "https://www.cbioportal.org/api"
SESSION = requests.Session()
SESSION.headers.update({"Accept": "application/json"})

STUDIES = [
    "thyroid_mskcc_2016",       # Landa 2016 PDTC + ATC
    "thyroid_gatci_2024",       # GATCI ATC, Cell Reports 2024
    "thpa_tcga_gdc",            # TCGA-PTC, GDC release
]


def _get(path: str, **params) -> Any:
    url = f"{API}{path}"
    r = SESSION.get(url, params=params, timeout=60)
    r.raise_for_status()
    return r.json()


def _post(path: str, payload: Any, **params) -> Any:
    url = f"{API}{path}"
    r = SESSION.post(url, params=params, json=payload, timeout=180)
    r.raise_for_status()
    return r.json()


def fetch_study(study_id: str) -> dict:
    log(WORKER, f"[{study_id}] fetching study metadata")
    study = _get(f"/studies/{study_id}")
    samples = _get(f"/studies/{study_id}/samples")
    profiles = _get(f"/studies/{study_id}/molecular-profiles")
    log(
        WORKER,
        f"[{study_id}] n_samples={len(samples)} n_profiles={len(profiles)}",
    )
    return {"study": study, "samples": samples, "profiles": profiles}


def pick_expression_profile(profiles: list[dict]) -> dict | None:
    """Prefer continuous mRNA RNA-seq, fallback to anything mRNA."""
    rna = [p for p in profiles if p.get("molecularAlterationType") == "MRNA_EXPRESSION"]
    if not rna:
        return None
    # priorities (substring match on molecularProfileId)
    priorities = (
        ("rsem", "zscore"),       # avoid zscore variants
        ("rna_seq", "zscore"),
        ("mrna", "zscore"),
        ("rsem",),
        ("rna_seq",),
        ("mrna",),
    )
    def has(p: dict, *needles: str) -> bool:
        s = p["molecularProfileId"].lower()
        return all(n in s for n in needles)
    def lacks(p: dict, *avoids: str) -> bool:
        s = p["molecularProfileId"].lower()
        return not any(a in s for a in avoids)
    # First pass: continuous (no zscore in id)
    for needles in (("rsem",), ("rna_seq_v2",), ("rna_seq",), ("mrna",)):
        for p in rna:
            if has(p, *needles) and lacks(p, "zscore"):
                return p
    # Fallback: anything mRNA
    return rna[0]


def fetch_expression(profile_id: str, sample_ids: list[str]) -> pd.DataFrame:
    """Fetch 8-gene × samples matrix from cBioPortal molecular-data endpoint."""
    payload = {
        "entrezGeneIds": list(PANEL_8_ENTREZ.values()),
        "sampleIds": sample_ids,
    }
    log(
        WORKER,
        f"[{profile_id}] fetching expr for {len(sample_ids)} samples × {len(PANEL_8)} genes",
    )
    rows = _post(
        f"/molecular-profiles/{profile_id}/molecular-data/fetch",
        payload,
        sampleListId=None,
    )
    if not rows:
        raise RuntimeError(f"empty response from {profile_id}")
    df = pd.DataFrame(rows)
    # Sample-by-gene matrix
    df["geneSymbol"] = df["entrezGeneId"].map({v: k for k, v in PANEL_8_ENTREZ.items()})
    pivot = df.pivot_table(
        index="sampleId", columns="geneSymbol", values="value", aggfunc="first"
    )
    return pivot


def score_study(study_id: str) -> dict | None:
    try:
        meta = fetch_study(study_id)
    except Exception as exc:
        log(WORKER, f"[{study_id}] FAIL fetch_study: {exc}")
        return None
    sample_ids = [s["sampleId"] for s in meta["samples"]]
    if not sample_ids:
        log(WORKER, f"[{study_id}] no samples — skip")
        return None
    profile = pick_expression_profile(meta["profiles"])
    if profile is None:
        log(WORKER, f"[{study_id}] no mRNA profile available — skip")
        return None
    pid = profile["molecularProfileId"]
    log(WORKER, f"[{study_id}] using profile={pid}")
    try:
        expr = fetch_expression(pid, sample_ids)
    except Exception as exc:
        log(WORKER, f"[{study_id}] FAIL fetch_expression: {exc}")
        return None
    # log2 if values look raw (rsem)
    if expr.values.max() > 50:
        expr = np.log2(expr.clip(lower=0) + 1)

    gene_map = resolve_panel(expr.columns)
    missing = [g for g in PANEL_8 if g not in gene_map]
    if missing:
        log(WORKER, f"[{study_id}] missing panel genes: {missing}")
    if len(gene_map) < 6:
        log(WORKER, f"[{study_id}] too few panel genes ({len(gene_map)}/8) — skip")
        return None

    expr = expr.dropna(thresh=int(0.6 * len(gene_map)))
    expr = expr.fillna(expr.median(numeric_only=True))
    scores = score_panel(expr, gene_map)
    dm = kmeans_dm(scores)

    out = pd.concat([expr, scores, dm.rename("DM_call")], axis=1)
    out.to_csv(OUT / f"panel_expression_{study_id}.tsv", sep="\t")
    dm_df = pd.concat([scores, dm.rename("DM_call")], axis=1)
    dm_df.to_csv(OUT / f"dm_calls_{study_id}.tsv", sep="\t")

    summary = {
        "study_id": study_id,
        "n_samples": int(out.shape[0]),
        "n_panel_genes_resolved": len(gene_map),
        "panel_gene_map": gene_map,
        "missing_panel_genes": missing,
        "dm1_count": int((dm == "DM1").sum()),
        "dm2_count": int((dm == "DM2").sum()),
        "rai_8_mean": float(scores["RAI_8"].mean()),
        "rai_8_std": float(scores["RAI_8"].std()),
        "dm1_like_mean": float(scores["DM1_like"].mean()),
    }
    log(WORKER, f"[{study_id}] OK n={summary['n_samples']} DM1={summary['dm1_count']} DM2={summary['dm2_count']}")
    return summary


def main() -> None:
    write_status(WORKER, "running", started_at=time.time())
    summaries = []
    for sid in STUDIES:
        s = score_study(sid)
        if s is not None:
            summaries.append(s)

    pooled = pd.DataFrame(summaries)
    pooled.to_csv(OUT / "trajectory_summary.tsv", sep="\t", index=False)
    (OUT / "pooled_meta.json").write_text(json.dumps(summaries, indent=2))

    write_status(
        WORKER,
        "done",
        n_studies=len(summaries),
        studies_succeeded=[s["study_id"] for s in summaries],
        finished_at=time.time(),
    )
    log(WORKER, f"DONE — {len(summaries)} studies scored")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(WORKER, f"FATAL: {exc}")
        write_status(WORKER, "error", error=str(exc))
        raise
