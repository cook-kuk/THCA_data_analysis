#!/usr/bin/env python3
"""CCLE thyroid cell line panel scoring via cBioPortal `ccle_broad_2025`.

This replaces the figshare-blocked DepMap path. Sources:
  - cBioPortal study: ccle_broad_2025 (Cancer Cell Line Encyclopedia, Nat Rev Cancer 2025).
  - Profile used: ccle_broad_2025_rna_seq_mrna (CONTINUOUS).
  - Sample filter: ONCOTREE in (THAP, THPA, THFO, THME) → 24 thyroid cell lines.

Outputs:
  thyroid_cell_panel.tsv     per-cell-line × {ONCOTREE, panel z, RAI_8, DM1_like, DM_call}
  oncotree_summary.tsv       grouped by ONCOTREE class
  thap_vs_thpa_cohens_d.json key effect-size headline
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

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

WORKER = "ccle_thyroid"
OUT = Path(__file__).parent
API = "https://www.cbioportal.org/api"
SESSION = requests.Session()
SESSION.headers.update({"Accept": "application/json"})

STUDY = "ccle_broad_2025"
PROFILE = "ccle_broad_2025_rna_seq_mrna"
THYROID_ONCOTREE = {"THAP", "THPA", "THFO", "THME", "THCO", "THYC", "THRA", "THRD"}


def main() -> None:
    write_status(WORKER, "running", started_at=time.time())

    log(WORKER, f"fetching ONCOTREE codes for {STUDY}")
    r = SESSION.get(f"{API}/studies/{STUDY}/clinical-data",
                    params={"clinicalDataType": "SAMPLE",
                            "attributeId": "ONCOTREE_CODE"},
                    timeout=60)
    r.raise_for_status()
    onc = pd.DataFrame(r.json())
    thy = onc[onc["value"].astype(str).str.upper().isin(THYROID_ONCOTREE)].copy()
    log(WORKER, f"thyroid cell lines: n={len(thy)}, ONCOTREE distribution = {dict(thy['value'].value_counts())}")
    sample_ids = thy["sampleId"].tolist()

    # Fetch panel expression for these 24 samples
    payload = {"entrezGeneIds": list(PANEL_8_ENTREZ.values()), "sampleIds": sample_ids}
    log(WORKER, f"fetching panel expression for {len(sample_ids)} cell lines × {len(PANEL_8)} genes")
    r = SESSION.post(f"{API}/molecular-profiles/{PROFILE}/molecular-data/fetch",
                     json=payload, timeout=180)
    r.raise_for_status()
    rows = r.json()
    df = pd.DataFrame(rows)
    df["geneSymbol"] = df["entrezGeneId"].map({v: k for k, v in PANEL_8_ENTREZ.items()})
    expr = df.pivot_table(index="sampleId", columns="geneSymbol", values="value", aggfunc="first")
    if expr.values.max() > 50:
        expr = np.log2(expr.clip(lower=0) + 1)

    gene_map = resolve_panel(expr.columns)
    missing = [g for g in PANEL_8 if g not in gene_map]
    log(WORKER, f"panel resolved {len(gene_map)}/{len(PANEL_8)} genes; missing={missing}")
    if len(gene_map) < 6:
        raise RuntimeError(f"thin panel coverage: {gene_map}")

    expr = expr.dropna(thresh=int(0.6 * len(gene_map)))
    expr = expr.fillna(expr.median(numeric_only=True))
    scores = score_panel(expr, gene_map)
    dm = kmeans_dm(scores)

    # Attach ONCOTREE class
    onc_map = thy.set_index("sampleId")["value"].to_dict()
    out = pd.concat([expr, scores, dm.rename("DM_call")], axis=1)
    out["ONCOTREE"] = out.index.map(onc_map)
    out.to_csv(OUT / "thyroid_cell_panel.tsv", sep="\t")
    log(WORKER, f"saved thyroid_cell_panel.tsv n={len(out)}")

    # Per-ONCOTREE summary
    summ = out.groupby("ONCOTREE").agg(
        n=("RAI_8", "count"),
        rai_8_mean=("RAI_8", "mean"),
        rai_8_std=("RAI_8", "std"),
        dm1_count=("DM_call", lambda x: (x == "DM1").sum()),
        dm2_count=("DM_call", lambda x: (x == "DM2").sum()),
        panel_log2_mean=("RAI_8", "mean"),  # placeholder — filled below
    )
    # Compute raw panel log2 means too
    panel_log2 = expr.loc[:, list(gene_map.values())].mean(axis=1)
    out["panel_log2"] = panel_log2
    summ["panel_log2_mean"] = out.groupby("ONCOTREE")["panel_log2"].mean()
    summ["panel_log2_std"] = out.groupby("ONCOTREE")["panel_log2"].std()
    summ.to_csv(OUT / "oncotree_summary.tsv", sep="\t")
    log(WORKER, f"summary by ONCOTREE:\n{summ.to_string()}")

    # ATC (THAP) vs PTC (THPA) headline contrast
    from scipy.stats import mannwhitneyu
    head = {}
    for code1, code2 in (("THAP", "THPA"), ("THAP", "THFO"), ("THPA", "THFO")):
        v1 = out.loc[out["ONCOTREE"] == code1, "panel_log2"].dropna()
        v2 = out.loc[out["ONCOTREE"] == code2, "panel_log2"].dropna()
        if len(v1) >= 2 and len(v2) >= 2:
            sd_pool = np.sqrt(((len(v1) - 1) * v1.var() + (len(v2) - 1) * v2.var())
                              / max(1, len(v1) + len(v2) - 2))
            d = float((v1.mean() - v2.mean()) / sd_pool) if sd_pool else float("nan")
            try:
                stat, p = mannwhitneyu(v1, v2, alternative="two-sided")
                p = float(p)
            except ValueError:
                p = float("nan")
            head[f"{code1}_vs_{code2}"] = {
                "n_first": int(len(v1)), "n_second": int(len(v2)),
                "mean_first_log2": float(v1.mean()),
                "mean_second_log2": float(v2.mean()),
                "delta_log2": float(v1.mean() - v2.mean()),
                "cohens_d": d,
                "mw_p": p,
            }
    (OUT / "thap_vs_thpa_cohens_d.json").write_text(json.dumps(head, indent=2))
    log(WORKER, f"headline contrasts: {json.dumps(head, indent=2)}")

    # Figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        order = ["THPA", "THFO", "THME", "THAP"]  # increasing aggressiveness
        present = [o for o in order if o in summ.index]
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        positions, labels = [], []
        for x, code in enumerate(present):
            v = out.loc[out["ONCOTREE"] == code, "panel_log2"].dropna()
            ax.boxplot(v, positions=[x], widths=0.7,
                       showmeans=True, meanline=True, patch_artist=True,
                       boxprops=dict(facecolor="#e6f0ff", edgecolor="#3b6ea8"))
            ax.scatter([x] * len(v), v, c="#b03a2e", s=22, zorder=3, alpha=0.7)
            positions.append(x)
            labels.append(f"{code}\n(n={len(v)})")
        ax.set_xticks(positions)
        ax.set_xticklabels(labels)
        ax.set_ylabel("8-gene panel mean log2 expression")
        ax.set_title("CCLE 2025 thyroid cell lines — 8-gene panel by ONCOTREE class\n(n=24 lines, cBioPortal `ccle_broad_2025`)")
        plt.tight_layout()
        (OUT.parent / "figures").mkdir(exist_ok=True)
        plt.savefig(OUT.parent / "figures" / "ccle_thyroid_panel.png", dpi=150)
        plt.close(fig)
        log(WORKER, "saved figures/ccle_thyroid_panel.png")
    except Exception as exc:
        log(WORKER, f"figure block failed: {exc}")

    write_status(WORKER, "done",
                 n_cell_lines=int(len(out)),
                 oncotree_breakdown={k: int(v) for k, v in thy["value"].value_counts().items()},
                 panel_genes_resolved=len(gene_map),
                 finished_at=time.time())
    log(WORKER, "DONE")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(WORKER, f"FATAL: {exc}")
        write_status(WORKER, "error", error=str(exc))
        raise
