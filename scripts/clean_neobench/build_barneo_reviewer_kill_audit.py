#!/usr/bin/env python3
"""Build reviewer kill-audit for BAR-Neo high-impact leads."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from common import dataframe_to_markdown, ensure_dir, update_manifest, write_tsv


OUTPUTS = [
    "barneo_high_impact_reviewer_kill_audit.tsv",
    "barneo_high_impact_neighbor_audit.tsv",
    "barneo_high_impact_decoy_pressure.tsv",
    "BAR_NEO_REVIEWER_KILL_AUDIT.md",
    "BAR_NEO_REVIEWER_KILL_AUDIT_KR.md",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root", default=".", help="Repository root")
    p.add_argument("--output-root", required=True, help="CLEAN-NeoBench output root")
    p.add_argument("--near-threshold", type=float, default=0.80, help="Near-peptide similarity threshold")
    return p.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def num(s: pd.Series, default: float = 0.0) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").fillna(default)


def truth(x: object) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"true", "1", "yes"}


def levenshtein(a: str, b: str) -> int:
    if a == b:
        return 0
    if not a:
        return len(b)
    if not b:
        return len(a)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def edit_similarity(a: str, b: str) -> float:
    denom = max(len(a), len(b), 1)
    return 1.0 - levenshtein(a, b) / denom


def kmers(seq: str, k: int = 3) -> set[str]:
    seq = str(seq)
    if len(seq) < k:
        return {seq}
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def kmer_jaccard(a: str, b: str, k: int = 3) -> float:
    ka, kb = kmers(a, k), kmers(b, k)
    return len(ka & kb) / max(len(ka | kb), 1)


def peptide_similarity(a: str, b: str) -> float:
    return max(edit_similarity(a, b), kmer_jaccard(a, b))


def build_neighbor_audit(leads: pd.DataFrame, master: pd.DataFrame, threshold: float) -> pd.DataFrame:
    rows = []
    if leads.empty or master.empty:
        return pd.DataFrame()
    m = master.copy()
    m["label_num"] = pd.to_numeric(m.get("label"), errors="coerce")
    for _, lead in leads.iterrows():
        cid = lead["candidate_id"]
        pep = str(lead.get("peptide", ""))
        hla = str(lead.get("hla_allele_4digit", ""))
        source = str(lead.get("source_name", ""))
        others = m[m["candidate_id"].ne(cid)].copy()
        sims = []
        for _, row in others.iterrows():
            sim = peptide_similarity(pep, str(row.get("peptide", "")))
            if sim >= threshold or str(row.get("peptide", "")) == pep:
                sims.append(
                    {
                        "neighbor_id": row.get("candidate_id", ""),
                        "neighbor_peptide": row.get("peptide", ""),
                        "neighbor_hla": row.get("hla_allele_4digit", ""),
                        "neighbor_source": row.get("source_name", ""),
                        "neighbor_label": row.get("label", np.nan),
                        "similarity": sim,
                        "same_hla": str(row.get("hla_allele_4digit", "")) == hla,
                        "same_source": str(row.get("source_name", "")) == source,
                    }
                )
        ndf = pd.DataFrame(sims)
        if ndf.empty:
            rows.append(
                {
                    "candidate_id": cid,
                    "peptide": pep,
                    "hla_allele_4digit": hla,
                    "source_name": source,
                    "n_near_neighbors": 0,
                    "n_near_neighbors_other_source": 0,
                    "n_near_neighbors_same_hla": 0,
                    "nearest_neighbor_similarity": 0.0,
                    "nearest_neighbor_id": "",
                    "nearest_neighbor_peptide": "",
                    "near_neighbor_positive_prevalence": np.nan,
                    "neighbor_audit_note": "no near neighbors at threshold",
                }
            )
            continue
        ndf = ndf.sort_values("similarity", ascending=False)
        rows.append(
            {
                "candidate_id": cid,
                "peptide": pep,
                "hla_allele_4digit": hla,
                "source_name": source,
                "n_near_neighbors": int(len(ndf)),
                "n_near_neighbors_other_source": int((~ndf["same_source"]).sum()),
                "n_near_neighbors_same_hla": int(ndf["same_hla"].sum()),
                "nearest_neighbor_similarity": float(ndf["similarity"].iloc[0]),
                "nearest_neighbor_id": ndf["neighbor_id"].iloc[0],
                "nearest_neighbor_peptide": ndf["neighbor_peptide"].iloc[0],
                "near_neighbor_positive_prevalence": float(pd.to_numeric(ndf["neighbor_label"], errors="coerce").mean()),
                "neighbor_audit_note": "near-neighbor caveat" if len(ndf) else "no near neighbors at threshold",
            }
        )
    return pd.DataFrame(rows)


def build_decoy_pressure(leads: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if leads.empty or scores.empty:
        return pd.DataFrame()
    s = scores.copy()
    s["score"] = num(s["stress_guarded_final_review_score"])
    s["label_num"] = pd.to_numeric(s.get("label"), errors="coerce")
    for _, lead in leads.iterrows():
        cid = lead["candidate_id"]
        row = s[s["candidate_id"].eq(cid)]
        if row.empty:
            continue
        r = row.iloc[0]
        score = float(r["score"])
        source = str(r.get("source_name", ""))
        hla = str(r.get("hla_allele_4digit", ""))
        same_source = s[s["source_name"].astype(str).eq(source)]
        same_hla = s[s["hla_allele_4digit"].astype(str).eq(hla)]
        same_pair = s[s["source_name"].astype(str).eq(source) & s["hla_allele_4digit"].astype(str).eq(hla)]
        def pressure(df: pd.DataFrame, name: str) -> dict[str, object]:
            above = df[df["score"].ge(score)]
            neg_above = int(above["label_num"].eq(0).sum())
            pos_above = int(above["label_num"].eq(1).sum())
            return {
                f"{name}_n": int(len(df)),
                f"{name}_rank": int((df["score"].gt(score)).sum() + 1),
                f"{name}_negatives_at_or_above": neg_above,
                f"{name}_positives_at_or_above": pos_above,
                f"{name}_precision_at_or_above": pos_above / max(pos_above + neg_above, 1),
            }
        out = {
            "candidate_id": cid,
            "source_name": source,
            "hla_allele_4digit": hla,
            "score": score,
        }
        out.update(pressure(same_source, "source"))
        out.update(pressure(same_hla, "hla"))
        out.update(pressure(same_pair, "source_hla"))
        rows.append(out)
    return pd.DataFrame(rows)


def build_kill_audit(leads: pd.DataFrame, overlap: pd.DataFrame, neighbor: pd.DataFrame, decoy: pd.DataFrame, contrib: pd.DataFrame) -> pd.DataFrame:
    if leads.empty:
        return pd.DataFrame()
    out = leads.copy()
    if not overlap.empty:
        out = out.merge(overlap, on="candidate_id", how="left")
    if not neighbor.empty:
        out = out.merge(neighbor.drop(columns=["peptide", "hla_allele_4digit", "source_name"], errors="ignore"), on="candidate_id", how="left")
    if not decoy.empty:
        out = out.merge(decoy.drop(columns=["source_name", "hla_allele_4digit"], errors="ignore"), on="candidate_id", how="left")
    if not contrib.empty:
        c = contrib[contrib["candidate_id"].isin(out["candidate_id"])].copy()
        c["candidate_method_contribution"] = num(c["candidate_method_contribution"])
        agg = c.groupby("candidate_id").agg(
            top_method=("method_name", "first"),
            n_contrib_methods=("method_name", "nunique"),
            public_contrib_methods=("method_role", lambda x: int((x == "caveated_public_comparator").sum())),
            fallback_contrib_methods=("method_role", lambda x: int((x == "bounded_fallback").sum())),
        ).reset_index()
        out = out.merge(agg, on="candidate_id", how="left")

    dispositions = []
    reasons = []
    for _, r in out.iterrows():
        reason = []
        if str(r.get("benchmark_label_status", "")).startswith("not_positive"):
            reason.append("not benchmark-positive")
        if truth(r.get("exact_peptide_hla_train_overlap")) or truth(r.get("exact_peptide_train_overlap")):
            reason.append("exact train overlap")
        if str(r.get("leakage_risk_level", "")).lower() == "high":
            reason.append("high leakage risk")
        if truth(r.get("near_peptide_train_overlap")) or float(r.get("nearest_neighbor_similarity", 0) or 0) >= 0.80:
            reason.append("near-peptide/neighbor caveat")
        if float(r.get("public_weight_fraction", 0) or 0) > 0.05:
            reason.append("public support too high")
        if float(r.get("fallback_weight_fraction", 0) or 0) > 0.08:
            reason.append("fallback support too high")
        if int(r.get("source_negatives_at_or_above", 0) or 0) > 0 and str(r.get("manual_audit_priority", "")).startswith("A_"):
            reason.append("same-source negatives at or above score")
        if not reason:
            dispositions.append("passes_current_reviewer_kill_audit")
            reasons.append("passes current automated leakage/neighbor/support/decoy checks")
        elif "high leakage risk" in reason or "exact train overlap" in reason:
            dispositions.append("blocked_from_clean_claim")
            reasons.append("; ".join(reason))
        else:
            dispositions.append("manual_audit_required")
            reasons.append("; ".join(reason))
    out["reviewer_kill_disposition"] = dispositions
    out["reviewer_kill_reason"] = reasons
    keep = [
        "candidate_id",
        "manual_audit_priority",
        "stress_guarded_rank_global",
        "stress_guarded_final_review_score",
        "benchmark_label_status",
        "source_name",
        "hla_allele_4digit",
        "peptide",
        "label",
        "leakage_risk_level",
        "exact_peptide_train_overlap",
        "exact_peptide_hla_train_overlap",
        "near_peptide_train_overlap",
        "n_near_neighbors",
        "nearest_neighbor_similarity",
        "nearest_neighbor_id",
        "source_negatives_at_or_above",
        "source_hla_negatives_at_or_above",
        "public_weight_fraction",
        "fallback_weight_fraction",
        "top_method",
        "reviewer_kill_disposition",
        "reviewer_kill_reason",
    ]
    return out[[c for c in keep if c in out.columns]]


def write_reports(output_root: Path, kill: pd.DataFrame, neighbor: pd.DataFrame, decoy: pd.DataFrame) -> None:
    counts = kill["reviewer_kill_disposition"].value_counts().rename_axis("disposition").reset_index(name="n") if not kill.empty else pd.DataFrame()
    text = f"""# BAR-Neo Reviewer Kill Audit

## Purpose

This report stress-tests the high-impact BAR-Neo leads against common reviewer attacks: leakage, exact/near peptide overlap, neighbor dependence, public/fallback dependence, and same-source negative pressure.

## Disposition Counts

{dataframe_to_markdown(counts, max_rows=20)}

## Kill Audit

{dataframe_to_markdown(kill, max_rows=40)}

## Near-Neighbor Audit

{dataframe_to_markdown(neighbor, max_rows=40)}

## Decoy / Negative Pressure

{dataframe_to_markdown(decoy, max_rows=40)}

## Claim Boundary

Passing this audit supports manual review priority only. It is not clinical vaccine selection, SOTA validation, or external validation.
"""
    (output_root / "BAR_NEO_REVIEWER_KILL_AUDIT.md").write_text(text.strip() + "\n")
    kr = f"""# BAR-Neo Reviewer Kill Audit KR

## 한 줄 결론

하이임팩트 lead를 리뷰어 공격 기준으로 다시 걸렀다. 핵심은 `passes_current_reviewer_kill_audit`만 다음 clean manual-review lead로 남기고, near-neighbor/high-leakage/exact-overlap 후보는 claim-safe에서 낮추는 것이다.

## Disposition counts

{dataframe_to_markdown(counts, max_rows=20)}

## Kill audit table

{dataframe_to_markdown(kill, max_rows=40)}

## Near-neighbor audit

{dataframe_to_markdown(neighbor, max_rows=40)}

## Decoy / negative pressure

{dataframe_to_markdown(decoy, max_rows=40)}

## 해석

- high leakage 또는 exact train overlap은 clean claim 차단.
- near-peptide caveat는 manual audit 필요.
- public/QK/fallback weight가 낮으면 public pretrained 또는 fallback이 lead를 만든 것이 아니라는 방어 근거가 된다.
- 이 audit을 통과해도 clinical/SOTA claim은 금지이며, manual biological/source/overlap audit의 우선순위로만 쓴다.
"""
    (output_root / "BAR_NEO_REVIEWER_KILL_AUDIT_KR.md").write_text(kr.strip() + "\n")


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    ensure_dir(output_root)

    leads = read_tsv(output_root / "barneo_high_impact_lead_candidates.tsv")
    master = read_tsv(output_root / "clean_neobench_master.tsv")
    overlap = read_tsv(output_root / "clean_neobench_overlap_flags.tsv")
    scores = read_tsv(output_root / "barneo_stress_guarded_candidate_scores.tsv")
    contrib = read_tsv(output_root / "barneo_stress_guarded_expert_contributions.tsv")

    neighbor = build_neighbor_audit(leads, master, args.near_threshold)
    decoy = build_decoy_pressure(leads, scores)
    kill = build_kill_audit(leads, overlap, neighbor, decoy, contrib)

    write_tsv(kill, output_root / "barneo_high_impact_reviewer_kill_audit.tsv")
    write_tsv(neighbor, output_root / "barneo_high_impact_neighbor_audit.tsv")
    write_tsv(decoy, output_root / "barneo_high_impact_decoy_pressure.tsv")
    write_reports(output_root, kill, neighbor, decoy)

    manifest_path = output_root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("output_files", [])
    for name in OUTPUTS:
        if name not in manifest["output_files"]:
            manifest["output_files"].append(name)
    manifest.setdefault("summary", {})
    manifest["summary"].update(
        {
            "n_reviewer_kill_audit_rows": int(len(kill)),
            "n_reviewer_kill_pass": int(kill["reviewer_kill_disposition"].eq("passes_current_reviewer_kill_audit").sum()) if not kill.empty else 0,
            "n_reviewer_kill_blocked": int(kill["reviewer_kill_disposition"].eq("blocked_from_clean_claim").sum()) if not kill.empty else 0,
            "n_reviewer_kill_manual_audit": int(kill["reviewer_kill_disposition"].eq("manual_audit_required").sum()) if not kill.empty else 0,
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    update_manifest(
        output_root,
        "barneo_reviewer_kill_audit",
        {
            "outputs": OUTPUTS,
            "n_rows": int(len(kill)),
            "warnings": ["Reviewer kill audit is manual-review triage, not external validation or clinical selection."],
        },
    )
    print(
        "[barneo-kill-audit] "
        f"rows={len(kill)} pass={manifest['summary']['n_reviewer_kill_pass']} "
        f"blocked={manifest['summary']['n_reviewer_kill_blocked']} manual={manifest['summary']['n_reviewer_kill_manual_audit']}"
    )


if __name__ == "__main__":
    main()
