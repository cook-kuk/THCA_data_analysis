#!/usr/bin/env python3
"""CPTAC THCA proteogenomics — protein-level confirmation of DM1 hypermethylation.

Public sources:
  PDC (proteomic data commons): TCGA / CPTAC THCA dataset
  cBioPortal mirror: thca_tcga_pan_can_atlas_2018 has CPTAC RPPA + protein partial
  Easiest entry: cptac Python package (pip), or NIH PDC API.

We try cptac Python package first (clean), fall back to PDC GraphQL API.

Goal:
  - For each TCGA-THCA sample with protein abundance, score 8-gene panel from RNA
    (using existing cBioPortal TCGA pull from sibling worker; or self-pull RNA).
  - Check that DM1 vs DM2 also separate at protein level: mean panel protein abundance
    should be lower in DM1 (mirror of mRNA result).
  - Compute mRNA-protein correlation per gene.
"""
from __future__ import annotations

import json
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _common import (  # noqa: E402
    PANEL_8,
    kmeans_dm,
    log,
    resolve_panel,
    score_panel,
    set_threads,
    write_status,
)

set_threads(2)

WORKER = "cptac_thca"
OUT = Path(__file__).parent
RAW = OUT / "_raw"
RAW.mkdir(exist_ok=True)


def try_cptac_pkg() -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """Use the cptac Python package if installed/installable. Returns (rna, protein)."""
    try:
        import cptac  # noqa: F401
    except ImportError:
        log(WORKER, "cptac pkg not installed — attempting pip install (venv mode)")
        try:
            import subprocess
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "--quiet", "cptac"],
                timeout=300,
            )
            import cptac  # noqa: F401
        except Exception as exc:
            log(WORKER, f"cptac install failed: {exc}")
            return None

    import cptac
    cptac.download(dataset="THCA")  # pyright: ignore[reportAttributeAccessIssue]
    th = cptac.Thca()  # pyright: ignore[reportAttributeAccessIssue]
    rna = th.get_transcriptomics()
    prot = th.get_proteomics()
    return rna, prot


def fallback_pdc_graphql() -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """Try direct PDC GraphQL fetch.

    Note: the PDC GraphQL schema is non-trivial; this fallback is mainly a marker.
    If cptac pkg fails we surface a partial-success status and let the user retry.
    """
    log(WORKER, "PDC GraphQL fallback not implemented — leaving stub")
    return None


def main() -> None:
    write_status(WORKER, "running", started_at=time.time())

    pair = try_cptac_pkg()
    if pair is None:
        pair = fallback_pdc_graphql()
    if pair is None:
        msg = ("cptac pkg unavailable and PDC fallback not implemented; "
               "skipping CPTAC for this run. Re-run after `pip install cptac`.")
        log(WORKER, msg)
        write_status(WORKER, "skipped", reason=msg)
        return

    rna, prot = pair
    log(WORKER, f"rna shape={rna.shape} protein shape={prot.shape}")

    # CPTAC frames are sample × gene. Convert to log2 if raw (heuristic).
    if rna.values.max() > 50:
        rna = np.log2(rna.clip(lower=0) + 1)

    common_samples = rna.index.intersection(prot.index)
    log(WORKER, f"common samples (RNA ∩ protein): n={len(common_samples)}")
    if len(common_samples) < 10:
        write_status(WORKER, "thin", n_common=int(len(common_samples)))
        return

    rna_c = rna.loc[common_samples]
    prot_c = prot.loc[common_samples]

    gene_map = resolve_panel(rna_c.columns)
    if len(gene_map) < 6:
        log(WORKER, f"thin RNA panel coverage ({len(gene_map)}/8) — abort")
        write_status(WORKER, "thin", n_panel=len(gene_map))
        return

    rna_panel = rna_c.loc[:, list(gene_map.values())].dropna(thresh=int(0.7 * len(gene_map)))
    rna_panel = rna_panel.fillna(rna_panel.median(numeric_only=True))
    scores = score_panel(rna_panel, gene_map)
    dm = kmeans_dm(scores)

    # Protein-side panel (same gene_map)
    prot_panel = prot_c.loc[:, [g for g in gene_map.values() if g in prot_c.columns]]
    log(WORKER, f"protein panel coverage: {prot_panel.shape[1]}/{len(gene_map)} genes")
    if prot_panel.shape[1] < 4:
        log(WORKER, "too few protein-detected panel genes — partial output only")

    # Per-gene DM1 vs DM2 protein delta
    from scipy.stats import mannwhitneyu

    rows = []
    dm1_samples = [s for s in dm.index if dm.loc[s] == "DM1"]
    dm2_samples = [s for s in dm.index if dm.loc[s] == "DM2"]
    for g in prot_panel.columns:
        v1 = prot_panel.loc[[s for s in dm1_samples if s in prot_panel.index], g].dropna()
        v2 = prot_panel.loc[[s for s in dm2_samples if s in prot_panel.index], g].dropna()
        if len(v1) < 3 or len(v2) < 3:
            continue
        try:
            stat, p = mannwhitneyu(v1, v2, alternative="two-sided")
        except ValueError:
            continue
        sd = np.sqrt(((len(v1) - 1) * v1.var() + (len(v2) - 1) * v2.var()) / max(1, len(v1) + len(v2) - 2))
        d = float((v1.mean() - v2.mean()) / sd) if sd else float("nan")
        rows.append({
            "gene": g, "n_dm1": int(len(v1)), "n_dm2": int(len(v2)),
            "mean_dm1": float(v1.mean()), "mean_dm2": float(v2.mean()),
            "cohens_d": d, "mw_p": float(p),
        })
    diff = pd.DataFrame(rows).sort_values("mw_p")
    diff.to_csv(OUT / "protein_dm1_vs_dm2.tsv", sep="\t", index=False)

    # mRNA-protein per-gene correlation (cross-check)
    corr_rows = []
    for g in prot_panel.columns:
        rna_g = rna_panel.get(g)
        prot_g = prot_panel[g]
        if rna_g is None:
            continue
        joined = pd.concat([rna_g.rename("rna"), prot_g.rename("prot")], axis=1).dropna()
        if len(joined) < 10:
            continue
        spear = joined.corr(method="spearman").loc["rna", "prot"]
        pear = joined.corr(method="pearson").loc["rna", "prot"]
        corr_rows.append({"gene": g, "n": int(len(joined)), "spearman": float(spear), "pearson": float(pear)})
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(OUT / "rna_protein_correlation.tsv", sep="\t", index=False)

    summary = {
        "n_common_samples": int(len(common_samples)),
        "n_panel_resolved_rna": len(gene_map),
        "n_panel_resolved_protein": int(prot_panel.shape[1]),
        "n_dm1": int((dm == "DM1").sum()),
        "n_dm2": int((dm == "DM2").sum()),
        "panel_gene_map": gene_map,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    write_status(WORKER, "done", **{k: v for k, v in summary.items() if not isinstance(v, dict)}, finished_at=time.time())
    log(WORKER, "DONE")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(WORKER, f"FATAL: {exc}")
        write_status(WORKER, "error", error=str(exc))
        raise
