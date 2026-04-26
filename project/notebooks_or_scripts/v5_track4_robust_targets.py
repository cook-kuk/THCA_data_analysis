#!/usr/bin/env python3
"""v5 Track 4 — Batch-robust novel target shortlist.

Combines Track 1 + Track 2 outputs with existing v3 replication /
drug-discovery / literature tables to produce:
  - Top 20 batch-robust shortlist (overall score)
  - Top 5 zero-literature novel (score high + 0 PubMed)
  - Top 5 fully-validated known (sanity check)

Score:
  0.30 * replication_rate
+ 0.25 * adversarial_shap_stability_rank (normalized)
+ 0.15 * correction_method_survival_count / n_methods
+ 0.15 * (ChEMBL pChEMBL>=6 present?)
+ 0.10 * (PDB or AlphaFold structure present?)
+ 0.05 * thyroid PubMed >=1 flag
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

sys.path.insert(0, str(Path(__file__).parent))
from v5_common import (  # noqa: E402
    RESULTS_ML,
    RESULTS_TABLES,
    ASSETS_DATA,
    FIGS_INTERACTIVE,
    safe_write_json,
)


def log(msg: str) -> None:
    print(f"[v5-track4] {msg}", flush=True)


def load_input(path: Path, **kw):
    if not path.exists():
        log(f"  missing {path.name}")
        return None
    try:
        return pd.read_csv(path, sep="\t", **kw)
    except Exception as e:
        log(f"  failed to read {path.name}: {e}")
        return None


def load_chem_targets():
    path = ASSETS_DATA / "chem_targets.json"
    if not path.exists():
        return {}
    try:
        obj = json.loads(path.read_text())
        if isinstance(obj, dict):
            if "targets" in obj:
                arr = obj["targets"]
            elif "rows" in obj:
                arr = obj["rows"]
            elif "items" in obj:
                arr = obj["items"]
            else:
                arr = next((v for v in obj.values() if isinstance(v, list)), [])
        else:
            arr = obj
        out = {}
        for item in arr:
            if not isinstance(item, dict):
                continue
            g = item.get("gene") or item.get("symbol") or item.get("target")
            if not g:
                continue
            out[g.upper()] = item
        return out
    except Exception as e:
        log(f"  chem_targets parse failed: {e}")
        return {}


def main() -> int:
    # ----- inputs -----
    val = load_input(RESULTS_TABLES / "biomarker_validated.tsv")
    if val is None:
        log("FATAL: biomarker_validated.tsv missing")
        return 1
    val["gene"] = val["gene"].astype(str).str.upper()
    log(f"validated biomarkers: {len(val)}")

    # replication_rate — assume 1.0 for now, refine via replicated flags if present
    replication_rate = pd.Series(1.0, index=val.index)
    for c in ["replicated_27155", "replicated_126698", "replicated_213647"]:
        if c in val.columns:
            replication_rate = replication_rate * 0
            break
    # build replication rate from per-cohort replication flags
    rep_cols = [c for c in val.columns if c.startswith("replicated_")]
    if rep_cols:
        rr = val[rep_cols].fillna(0).astype(float).mean(axis=1)
        val["replication_rate"] = rr
    else:
        val["replication_rate"] = 1.0

    # adversarial SHAP stability
    shap_stable = load_input(RESULTS_ML / "v5_track2_shap_stable_genes.tsv")
    if shap_stable is not None and "gene" in shap_stable.columns:
        shap_stable["gene"] = shap_stable["gene"].str.upper()
        max_n = shap_stable["n_lambdas_in_top10"].max() if len(shap_stable) else 1
        shap_map = dict(zip(shap_stable["gene"],
                            shap_stable["n_lambdas_in_top10"] / max(max_n, 1)))
    else:
        shap_map = {}
    val["adversarial_shap_stability_rank"] = val["gene"].map(
        lambda g: shap_map.get(g, 0.0))

    # Track 1 correction-survival: how many methods have non-nan LODO and kept the gene well
    # Approximate: a gene "survives" correction if it remains in top 500 variable genes post-correction
    t1_top = load_input(RESULTS_ML / "v5_track1_top_corrected_genes.tsv")
    correction_methods = load_input(RESULTS_ML / "v5_track1_correction_comparison.tsv")
    n_methods_ok = (correction_methods["status"] == "OK").sum() \
        if correction_methods is not None and "status" in correction_methods.columns else 5
    if t1_top is not None:
        survivors = set(t1_top["gene"].astype(str).str.upper().tolist())
    else:
        survivors = set()
    # We count survival as a binary for now (in top-stable genes). For more
    # methods, replicate the top-gene extraction across methods; here we use
    # a single representative (ComBat) + imply >=2 methods if present.
    val["correction_method_survival_count"] = val["gene"].apply(
        lambda g: 2 if g in survivors else 0)

    # Drug discovery: ChEMBL pChEMBL>=6
    dd = load_input(RESULTS_TABLES / "drug_discovery_targets.tsv")
    chembl_hits = set()
    if dd is not None:
        for _, r in dd.iterrows():
            try:
                pch = float(r.get("best_pchembl", 0) or 0)
            except Exception:
                pch = 0
            if pch >= 6:
                g = str(r["gene"]).upper()
                chembl_hits.add(g)
    val["has_chembl_pchembl_ge6"] = val["gene"].isin(chembl_hits).astype(int)

    # Structure (PDB/AlphaFold) — from chem_targets.json
    chem_map = load_chem_targets()
    log(f"chem_targets entries: {len(chem_map)}")
    def has_structure(g):
        item = chem_map.get(g)
        if not item:
            return 0
        for k in ["has_pdb", "has_alphafold", "has_structure", "pdb_id", "pdb",
                  "alphafold_id", "af_id"]:
            v = item.get(k) if isinstance(item, dict) else None
            if v:
                return 1
        return 0
    val["has_structure"] = val["gene"].apply(has_structure)

    # Thyroid literature count
    lit = load_input(RESULTS_TABLES / "drug_discovery_literature.tsv")
    lit_counts = {}
    if lit is not None:
        lit_counts = lit.groupby("target").size().to_dict()
    val["thyroid_pubmed_count"] = val["gene"].apply(
        lambda g: int(lit_counts.get(g, 0)))
    val["has_thyroid_pubmed"] = (val["thyroid_pubmed_count"] >= 1).astype(int)

    # Final score
    val["correction_survival_frac"] = val["correction_method_survival_count"] / max(n_methods_ok, 1)
    val["v5_robust_score"] = (
        0.30 * val["replication_rate"].fillna(0)
        + 0.25 * val["adversarial_shap_stability_rank"].fillna(0)
        + 0.15 * val["correction_survival_frac"].fillna(0)
        + 0.15 * val["has_chembl_pchembl_ge6"].fillna(0)
        + 0.10 * val["has_structure"].fillna(0)
        + 0.05 * val["has_thyroid_pubmed"].fillna(0)
    )

    score_cols = [
        "gene", "v5_robust_score", "replication_rate",
        "adversarial_shap_stability_rank", "correction_method_survival_count",
        "correction_survival_frac",
        "has_chembl_pchembl_ge6", "has_structure",
        "thyroid_pubmed_count", "has_thyroid_pubmed",
    ]
    ranked = val.sort_values("v5_robust_score", ascending=False)

    # Top 20 overall
    top20 = ranked.head(20)[score_cols].copy()
    # Top 5 zero-literature novel
    zero_lit = ranked[ranked["thyroid_pubmed_count"] == 0].head(5)[score_cols].copy()
    # Top 5 well-known
    known = ranked[ranked["thyroid_pubmed_count"] >= 2].head(5)[score_cols].copy()

    # Combine into a single shortlist with category
    top20["category"] = "batch-robust top-20"
    zero_lit["category"] = "zero-literature novel"
    known["category"] = "fully-validated known"
    shortlist = pd.concat([top20, zero_lit, known])
    shortlist.to_csv(RESULTS_TABLES / "v5_robust_target_shortlist.tsv",
                     sep="\t", index=False)
    log(f"shortlist: top20 n={len(top20)}, zero-lit n={len(zero_lit)}, "
        f"known n={len(known)}")

    # Payload for HTML
    payload = {
        "generated_at": pd.Timestamp.utcnow().isoformat(),
        "n_validated_total": int(len(val)),
        "top_20": top20.to_dict(orient="records"),
        "zero_literature_top_5": zero_lit.to_dict(orient="records"),
        "known_top_5": known.to_dict(orient="records"),
        "scoring_weights": {
            "replication_rate": 0.30,
            "adversarial_shap_stability_rank": 0.25,
            "correction_method_survival_count": 0.15,
            "chembl_pchembl_ge6": 0.15,
            "structure_pdb_or_alphafold": 0.10,
            "thyroid_pubmed_ge1": 0.05,
        },
    }
    (ASSETS_DATA / "v5_shortlist_payload.json").write_text(
        json.dumps(payload, indent=2, default=str))

    # Radar chart: top-8 targets × 6 score components
    top8 = ranked.head(8)
    components = [
        "replication_rate",
        "adversarial_shap_stability_rank",
        "correction_survival_frac",
        "has_chembl_pchembl_ge6",
        "has_structure",
        "has_thyroid_pubmed",
    ]
    fig = go.Figure()
    for _, r in top8.iterrows():
        vals = [float(r[c]) if not pd.isna(r[c]) else 0.0 for c in components]
        fig.add_trace(go.Scatterpolar(r=vals + [vals[0]],
                                       theta=components + [components[0]],
                                       fill="toself", name=str(r["gene"])))
    fig.update_layout(title="v5 Track 4 — top-8 batch-robust targets × 6 score components",
                      polar=dict(radialaxis=dict(range=[0, 1])))
    fig.write_html(FIGS_INTERACTIVE / "v5_shortlist_radar.html",
                   include_plotlyjs="cdn")

    # Literature count vs score scatter
    plot_df = ranked.head(200).copy()
    fig = px.scatter(plot_df, x="thyroid_pubmed_count", y="v5_robust_score",
                     hover_name="gene",
                     title="v5 Track 4 — thyroid PubMed count vs v5_robust_score (top-200)",
                     labels={"thyroid_pubmed_count": "thyroid PubMed count",
                             "v5_robust_score": "v5 robust score"})
    fig.add_hline(y=float(plot_df["v5_robust_score"].quantile(0.95)),
                  line_dash="dot", annotation_text="95th pct")
    fig.write_html(FIGS_INTERACTIVE / "v5_shortlist_literature_vs_score.html",
                   include_plotlyjs="cdn")

    # Summary JSON for synthesizer
    summary = {
        "n_validated_total": int(len(val)),
        "n_top_20": int(len(top20)),
        "n_zero_lit": int(len(zero_lit)),
        "n_known": int(len(known)),
        "top_3_robust_novel_genes": zero_lit["gene"].tolist()[:3],
        "top_3_robust_known_genes": known["gene"].tolist()[:3],
        "top_20_genes": top20["gene"].tolist(),
    }
    safe_write_json(RESULTS_ML / "v5_track4_summary.json", summary)
    log(f"done; top-20={summary['top_20_genes'][:5]}...")
    return 0


if __name__ == "__main__":
    sys.exit(main())
