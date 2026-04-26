"""
v13 Task 7 — Prioritize candidates with a composite score.

  score = 0.25 * binding_affinity_norm      (= max Tanimoto to reference binder)
        + 0.20 * admet_favorability
        + 0.15 * clinical_maturity           (from max_phase)
        + 0.15 * off_patent_status           (first_approval ≤ 2005)
        + 0.10 * thyroid_evidence            (from PubMed hits)
        + 0.10 * target_specificity          (from v13 modality fit)
        + 0.05 * novel_vs_repurpose          (0 for pure novel; 1 for repurpose)

Tiers:
  A (go now)       : ≥ 0.70 AND max_phase == 4 AND no red flag
  B (clinical dev) : 0.55–0.69 OR (phase 2–3 with no red flag)
  C (early)        : 0.40–0.54
  D (novel hyp.)   : < 0.40

Outputs:
  $RES/prioritization/ranked_candidates.tsv   (all 400)
  $RES/prioritization/tier_a_go_now.md        (executive summary)
  $RES/prioritization/per_target_top10.tsv
"""
from __future__ import annotations
from pathlib import Path
import json
import pandas as pd
import numpy as np
import requests, time

PROJECT = Path("/opt/thyroid-dash/project")
RES = PROJECT / "results" / "v13_drug_discovery"
ADMET = RES / "admet"
PRIO = RES / "prioritization"
PRIO.mkdir(parents=True, exist_ok=True)
MOD = RES / "modalities"

UA = {"User-Agent": "thca-v13 (kukshomr@gmail.com)"}


def thyroid_evidence_count(drug_name: str) -> int:
    """Small PubMed count proxy. Cached across calls."""
    cache = PRIO / "thyroid_pubmed_cache.json"
    if cache.exists():
        db = json.loads(cache.read_text())
    else:
        db = {}
    key = drug_name.lower().strip()
    if not key or key == "nan":
        return 0
    if key in db:
        return int(db[key])
    try:
        r = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params={"db": "pubmed", "term": f"{drug_name} AND thyroid", "retmax": 1, "retmode": "json"},
            headers=UA, timeout=20,
        )
        if r.status_code == 200:
            n = int(r.json().get("esearchresult", {}).get("count", 0))
        else:
            n = 0
    except Exception:
        n = 0
    db[key] = n
    cache.write_text(json.dumps(db))
    time.sleep(0.25)
    return n


def main():
    cand = pd.read_csv(ADMET / "top_candidates_admet.tsv", sep="\t")
    # Merge screen info (first_approval, max_phase etc. already present)
    print(f"Input: {len(cand)} candidates")

    # Load modality fit to upweight targets with matching modality
    mod_json = json.loads((MOD / "target_modality_matrix.json").read_text())

    def target_spec(row):
        # SM-focused: small-molecule fit is relevant
        g = row["target"]
        m_info = mod_json.get(g, {}).get("SM", {"score": 0})
        return m_info["score"] / 4  # 0..1

    cand["target_spec"] = cand.apply(target_spec, axis=1)

    # Thyroid evidence: PubMed count on drug + thyroid. Cache to avoid repeated calls.
    print("Fetching PubMed counts for unique drug names (cached)...")
    unique_names = cand["name"].dropna().unique()
    counts = {}
    for i, n in enumerate(unique_names):
        counts[n] = thyroid_evidence_count(str(n))
        if (i+1) % 50 == 0:
            print(f"  {i+1}/{len(unique_names)}")
    cand["thyroid_pubmed"] = cand["name"].map(counts).fillna(0).astype(int)

    # Normalizations
    cand["binding_norm"] = cand["max_sim_to_ref"].fillna(0).clip(0, 1)

    def maturity(p):
        try:
            p = float(p)
        except Exception:
            return 0.0
        return min(1.0, max(0.0, p / 4.0))
    cand["maturity_norm"] = cand["max_phase"].apply(maturity)

    # off_patent: first_approval present and ≤ 2005 → off patent
    def off_patent(fa):
        try:
            y = int(float(fa))
        except Exception:
            return 0.3  # unknown ≈ neutral
        if y == 0:
            return 0.3
        if y <= 2005:
            return 1.0
        if y <= 2015:
            return 0.5
        return 0.1
    cand["off_patent"] = cand["first_approval"].apply(off_patent)

    # thyroid evidence normalized (log-scale saturation)
    cand["thyroid_norm"] = cand["thyroid_pubmed"].apply(lambda n: min(1.0, np.log1p(n)/np.log1p(200)))

    # novel_vs_repurpose: all phase-4 FDA candidates are "repurpose" (=1.0)
    cand["repurpose_norm"] = cand["maturity_norm"].apply(lambda x: 1.0 if x >= 0.99 else 0.6)

    cand["composite_score"] = (
        0.25 * cand["binding_norm"] +
        0.20 * cand["admet_favorability"].fillna(0) +
        0.15 * cand["maturity_norm"] +
        0.15 * cand["off_patent"] +
        0.10 * cand["thyroid_norm"] +
        0.10 * cand["target_spec"] +
        0.05 * cand["repurpose_norm"]
    ).round(4)

    cand["red_flag"] = (
        (cand["herg_risk"] > 0.6) |
        (cand["pains_count"] > 0) |
        (cand["hepatotox_flags"].fillna("") != "")
    )

    def tier(row):
        s = row["composite_score"]
        if s >= 0.70 and row["max_phase"] == 4 and not row["red_flag"]:
            return "A"
        if s >= 0.55 or (row["max_phase"] in (2.0, 3.0) and not row["red_flag"]):
            return "B"
        if s >= 0.40:
            return "C"
        return "D"
    cand["tier"] = cand.apply(tier, axis=1)

    cand = cand.sort_values("composite_score", ascending=False)
    cand.to_csv(PRIO / "ranked_candidates.tsv", sep="\t", index=False)
    print(f"Wrote ranked_candidates.tsv ({len(cand)} rows)")

    # per-target top10
    per_tgt = cand.groupby("target").head(10)
    per_tgt.to_csv(PRIO / "per_target_top10.tsv", sep="\t", index=False)

    # Tier counts
    print("\nTier distribution:")
    print(cand["tier"].value_counts())

    # Tier A exec summary
    tier_a = cand[cand["tier"] == "A"].copy()
    lines = ["# v13 Tier A — 'Go Now' Candidates",
             "",
             f"Total Tier-A candidates: **{len(tier_a)}**",
             "",
             "Tier-A criteria: composite ≥ 0.70, FDA-approved (max_phase=4), no PAINS / hepatotox / hERG red flag.",
             "",
             "| Rank | Target | Drug | ChEMBL | Composite | Tanimoto | ADMET | Thyroid PubMed | First approval |",
             "|-----:|:-------|:-----|:-------|----------:|---------:|------:|---------------:|---------------:|"]
    for i, (_, r) in enumerate(tier_a.head(50).iterrows(), 1):
        lines.append(
            f"| {i} | {r['target']} | {r['name']} | {r['chembl_id']} | "
            f"{r['composite_score']:.3f} | {r['max_sim_to_ref']:.3f} | "
            f"{r['admet_favorability']:.3f} | {int(r['thyroid_pubmed'])} | "
            f"{r['first_approval'] if pd.notna(r['first_approval']) else ''} |"
        )
    (PRIO / "tier_a_go_now.md").write_text("\n".join(lines))
    print(f"Wrote tier_a_go_now.md  (top Tier A = {len(tier_a)})")

    # per-target top3 condensed
    print("\n=== Per-target top 3 ===")
    for g in sorted(cand["target"].unique()):
        top = cand[cand["target"] == g].head(3)
        print(f"\n[{g}]")
        for _, r in top.iterrows():
            print(f"  {r['name']:<30}  score={r['composite_score']:.3f}  phase={r['max_phase']}  sim={r['max_sim_to_ref']:.2f}  admet={r['admet_favorability']:.2f}  thyr_pubmed={int(r['thyroid_pubmed'])}")


if __name__ == "__main__":
    main()
