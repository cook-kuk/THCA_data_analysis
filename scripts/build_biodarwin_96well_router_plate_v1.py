#!/usr/bin/env python3
"""Build a 96-well BioDarwin router plate from the public industrial scout.

Plate composition:
- 67 high-confidence metric-ready scout rows
- 29 exploratory rows from the broader public industrial candidate pool

The goal is not to claim SOTA; it is to build a useful 96-well wetlab board
that mixes reliable positives/negatives with router-disagreement probes and
helper-lane candidates.
"""
from __future__ import annotations

import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/seungho/personal/THCA_data_analysis")
SCOUT = ROOT / "project/results/public_industrial_neoantigen_benchmark_scout_2026_05_10"
OUT = ROOT / "project/results/biodarwin_96well_router_plate_2026_05_11"
OUT.mkdir(parents=True, exist_ok=True)


def peptide_features(peptide: object) -> dict[str, float]:
    pep = str(peptide or "").upper()
    if not pep or pep == "NAN":
        return {
            "length_pref": 0.0,
            "diversity": 0.0,
            "hydrophobic": 0.0,
            "aromatic": 0.0,
            "charged_balance": 0.0,
            "peptide_quality_proxy": 0.0,
        }
    n = len(pep)
    length_pref = 1.0 if n in {9, 10, 11} else 0.75 if n in {8, 12, 13, 14, 15} else 0.0
    diversity = min(1.0, len(set(pep)) / 8.0)
    hydrophobic = sum(aa in "AILMVPFWY" for aa in pep) / n
    aromatic = sum(aa in "FWY" for aa in pep) / n
    charged = sum(aa in "KRDE" for aa in pep) / n
    charged_balance = 1.0 - min(1.0, abs(charged - 0.22) / 0.35)
    quality = 0.42 * length_pref + 0.22 * diversity + 0.16 * hydrophobic + 0.10 * aromatic + 0.10 * charged_balance
    return {
        "length_pref": float(length_pref),
        "diversity": float(diversity),
        "hydrophobic": float(hydrophobic),
        "aromatic": float(aromatic),
        "charged_balance": float(charged_balance),
        "peptide_quality_proxy": float(quality),
    }


def normalize_hla(hla: object) -> str:
    s = str(hla or "").upper().replace("HLA-", "").strip()
    return s


def assign_wells(df: pd.DataFrame) -> pd.DataFrame:
    wells = []
    rows = list("ABCDEFGH")
    cols = list(range(1, 13))
    for r in rows:
        for c in cols:
            wells.append(f"{r}{c:02d}")
    if len(df) > len(wells):
        raise ValueError("too many rows for 96-well plate")
    out = df.copy().reset_index(drop=True)
    out["plate_well"] = wells[: len(out)]
    return out


def pick_top(df: pd.DataFrame, n: int, score_col: str, allow_cols: dict[str, int] | None = None, ascending=False) -> pd.DataFrame:
    if df.empty or n <= 0:
        return df.head(0).copy()
    d = df.copy()
    if allow_cols:
        # Greedy cap by source to keep diversity while preserving ranking.
        d = d.sort_values(score_col, ascending=ascending).copy()
        chosen = []
        seen = defaultdict(int)
        for _, row in d.iterrows():
            key = row["source_id"]
            cap = allow_cols.get(key, np.inf)
            if seen[key] < cap:
                chosen.append(row)
                seen[key] += 1
            if len(chosen) >= n:
                break
        if len(chosen) < n:
            remaining = d.loc[~d.index.isin([r.name for r in chosen])]
            for _, row in remaining.iterrows():
                chosen.append(row)
                if len(chosen) >= n:
                    break
        return pd.DataFrame(chosen).head(n).copy()
    return d.sort_values(score_col, ascending=ascending).head(n).copy()


def main() -> None:
    candidate = pd.read_csv(SCOUT / "public_industrial_candidate_pairs.tsv", sep="\t")
    high_conf = pd.read_csv(SCOUT / "public_industrial_high_confidence_candidate_pairs.tsv", sep="\t")
    locked = pd.read_csv(SCOUT / "public_industrial_locked_testset_v0.tsv", sep="\t")

    def length_key(s: pd.Series) -> pd.Series:
        return pd.to_numeric(s, errors="coerce").fillna(-1).astype(int).astype(str)

    # Build a canonical key.
    candidate["key"] = (
        candidate["source_id"].astype(str)
        + "|"
        + candidate["peptide_or_sequence"].astype(str).str.upper()
        + "|"
        + candidate["hla_allele"].astype(str)
        + "|"
        + length_key(candidate["length"])
    )
    high_conf["key"] = (
        high_conf["source_id"].astype(str)
        + "|"
        + high_conf["peptide"].astype(str).str.upper()
        + "|"
        + high_conf["hla_allele"].astype(str)
        + "|"
        + length_key(high_conf["peptide_length"])
    )
    locked["key"] = (
        locked["source_id"].astype(str)
        + "|"
        + locked["peptide"].astype(str).str.upper()
        + "|"
        + locked["hla_allele"].astype(str)
        + "|"
        + length_key(locked["peptide_length"])
    )

    core = candidate.merge(high_conf[["key", "line_number", "response_label", "label_raw", "high_confidence_metric_ready"]], on="key", how="inner")
    core = core.merge(locked[["key", "label", "benchmark_use"]], on="key", how="left", suffixes=("", "_locked"))
    core["bucket"] = "core_high_confidence"
    core["selection_reason"] = np.where(core["key"].isin(set(locked["key"])), "locked_metric_row", "high_confidence_scout")
    core["selection_score"] = 1.0

    orphan_hc = high_conf.loc[~high_conf["key"].isin(set(core["key"]))].copy()
    orphan_hc["bucket"] = "core_high_confidence_orphan"
    orphan_hc["selection_reason"] = "high_confidence_orphan"
    orphan_hc["selection_score"] = 1.0

    # Exploratory pool: direct class-I / class-II ready rows not already in core.
    pool = candidate.loc[~candidate["key"].isin(set(core["key"]))].copy()
    pool = pool[pool["direct_classI_score_ready"] | pool["direct_classII_score_ready"] | pool["pan_vaccine_score_ready"]].copy()
    pf = pd.DataFrame([peptide_features(x) for x in pool["peptide_or_sequence"]], index=pool.index)
    pool = pd.concat([pool.reset_index(drop=True), pf.reset_index(drop=True)], axis=1)
    pool["source_count"] = pool.groupby("source_id")["source_id"].transform("count")
    pool["source_balance"] = 1.0 / np.sqrt(pool["source_count"].astype(float))
    pool["tier_bonus"] = np.where(pool["candidate_label_tier"].astype(str).eq("A_like"), 1.0, 0.0)
    pool["classI_bonus"] = np.where(pool["direct_classI_score_ready"].astype(bool), 1.0, 0.0)
    pool["classII_bonus"] = np.where(pool["direct_classII_score_ready"].astype(bool), 1.0, 0.0)
    pool["pan_bonus"] = np.where(pool["pan_vaccine_score_ready"].astype(bool), 1.0, 0.0)
    pool["explore_score"] = (
        0.30 * pool["peptide_quality_proxy"]
        + 0.20 * pool["length_pref"]
        + 0.15 * pool["diversity"]
        + 0.10 * pool["hydrophobic"]
        + 0.08 * pool["aromatic"]
        + 0.07 * pool["charged_balance"]
        + 0.10 * pool["source_balance"]
        + 0.05 * pool["tier_bonus"]
        + 0.04 * pool["classI_bonus"]
        + 0.03 * pool["classII_bonus"]
        + 0.03 * pool["pan_bonus"]
    )

    # Selection buckets for the remaining 29 wells.
    classI = pool[(pool["candidate_label_tier"].eq("A_like")) & (pool["direct_classI_score_ready"].astype(bool)) & (pool["mhc_class"].astype(str).eq("I"))].copy()
    classII = pool[(pool["direct_classII_score_ready"].astype(bool)) & (pool["mhc_class"].astype(str).eq("II"))].copy()
    controls = pool[(pool["candidate_label_tier"].eq("B_or_C_like")) & (pool["direct_classI_score_ready"].astype(bool)) & (pool["mhc_class"].astype(str).eq("I"))].copy()
    diverse = pool.copy()

    # Pull rows in a way that preserves source diversity.
    pick_map = defaultdict(int)

    def top_diverse(df: pd.DataFrame, n: int, cap_per_source: int, ascending=False) -> pd.DataFrame:
        if df.empty:
            return df.head(0).copy()
        d = df.sort_values("explore_score", ascending=ascending).copy()
        chosen = []
        seen = defaultdict(int)
        for _, row in d.iterrows():
            src = row["source_id"]
            if seen[src] >= cap_per_source:
                continue
            chosen.append(row)
            seen[src] += 1
            if len(chosen) >= n:
                break
        if len(chosen) < n:
            rest = d.loc[~d.index.isin([r.name for r in chosen])]
            for _, row in rest.iterrows():
                chosen.append(row)
                if len(chosen) >= n:
                    break
        return pd.DataFrame(chosen).head(n).copy()

    extra_a = top_diverse(classI, 12, cap_per_source=3)
    extra_ii = top_diverse(classII, 6, cap_per_source=2)
    extra_ctrl = controls.sort_values("explore_score", ascending=True).head(6).copy()
    extra_src = top_diverse(diverse, 5, cap_per_source=1)

    extras = pd.concat([extra_a, extra_ii, extra_ctrl, extra_src], ignore_index=True)
    extras = extras.loc[~extras["key"].isin(set(core["key"]))].copy()
    extras = extras.drop_duplicates(subset=["key"]).head(29).copy()

    if len(core) + len(extras) != 96:
        # Backfill from remaining pool by score if some buckets underfill.
        remaining = pool.loc[~pool["key"].isin(set(pd.concat([core, extras], ignore_index=True)["key"]))].copy()
        remaining = remaining.sort_values("explore_score", ascending=False)
        need = 96 - (len(core) + len(extras))
        if need > 0:
            extras = pd.concat([extras, remaining.head(need)], ignore_index=True)
        extras = extras.head(29).copy()

    core["plate_bucket"] = np.where(core["selection_reason"].eq("locked_metric_row"), "core_locked_metric", "core_high_confidence")
    extras["plate_bucket"] = np.where(extras["candidate_label_tier"].astype(str).eq("A_like"), "explore_A_like", "explore_controls")
    extras.loc[extras["direct_classII_score_ready"].astype(bool), "plate_bucket"] = "explore_classII"
    extras.loc[extras["candidate_label_tier"].eq("B_or_C_like"), "plate_bucket"] = "explore_controls"

    core_out = pd.DataFrame(
        {
            "source_id": core["source_id"],
            "line_number": core["line_number"],
            "peptide": core["peptide_or_sequence"],
            "hla_allele": core["hla_allele"],
            "mhc_class": core["mhc_class"],
            "candidate_label_tier": core["candidate_label_tier"],
            "response_label": core["response_label"],
            "label_raw": core["label_raw"],
            "plate_bucket": core["plate_bucket"],
            "selection_score": core["selection_score"],
            "selection_reason": core["selection_reason"],
            "analysis_role": np.where(core["selection_reason"].eq("locked_metric_row"), "locked_metric_control", "high_confidence_core"),
        }
    )
    orphan_out = pd.DataFrame(
        {
            "source_id": orphan_hc["source_id"],
            "line_number": orphan_hc["line_number"],
            "peptide": orphan_hc["peptide"],
            "hla_allele": orphan_hc["hla_allele"],
            "mhc_class": orphan_hc["mhc_class"],
            "candidate_label_tier": np.where(orphan_hc["response_label"].astype(float) > 0, "A_like", "B_or_C_like"),
            "response_label": orphan_hc["response_label"],
            "label_raw": orphan_hc["label_raw"],
            "plate_bucket": "core_high_confidence_orphan",
            "selection_score": orphan_hc["selection_score"],
            "selection_reason": orphan_hc["selection_reason"],
            "analysis_role": "high_confidence_core",
        }
    )
    extra_out = pd.DataFrame(
        {
            "source_id": extras["source_id"],
            "line_number": np.nan,
            "peptide": extras["peptide_or_sequence"],
            "hla_allele": extras["hla_allele"],
            "mhc_class": extras["mhc_class"],
            "candidate_label_tier": extras["candidate_label_tier"],
            "response_label": np.nan,
            "label_raw": extras["candidate_label_tier"] if "candidate_label_tier" in extras.columns else "",
            "plate_bucket": extras["plate_bucket"],
            "selection_score": extras["explore_score"],
            "selection_reason": np.where(extras["plate_bucket"].eq("explore_classII"), "helper_lane_probe",
                                         np.where(extras["plate_bucket"].eq("explore_controls"), "specificity_control", "disagreement_probe")),
            "analysis_role": np.where(extras["plate_bucket"].eq("explore_classII"), "helper_lane",
                                      np.where(extras["plate_bucket"].eq("explore_controls"), "negative_control", "exploratory_probe")),
        }
    )

    plate = pd.concat([core_out, orphan_out, extra_out], ignore_index=True)
    # Make a stable interleaved order to avoid clustering all cores first.
    order = []
    buckets = ["core_locked_metric", "core_high_confidence", "core_high_confidence_orphan", "explore_A_like", "explore_classII", "explore_controls"]
    bucket_frames = {b: plate[plate["plate_bucket"].eq(b)].sort_values(["selection_score", "source_id"], ascending=[False, True]).copy() for b in buckets}
    while sum(len(v) for v in bucket_frames.values()) > 0:
        for b in buckets:
            if len(bucket_frames[b]) == 0:
                continue
            row = bucket_frames[b].iloc[[0]]
            bucket_frames[b] = bucket_frames[b].iloc[1:].copy()
            order.append(row)
    plate = pd.concat(order, ignore_index=True)
    plate = assign_wells(plate)
    plate["rank"] = np.arange(1, len(plate) + 1)

    # Append summary fields useful for the wetlab interpreter.
    plate["hit_fail_interpretation"] = np.where(
        plate["analysis_role"].eq("locked_metric_control"),
        "known_locked_reference",
        np.where(
            plate["analysis_role"].eq("negative_control"),
            "expected_fail_reference",
            np.where(
                plate["analysis_role"].eq("helper_lane"),
                "helper_lane_candidate",
                "exploratory_candidate",
            ),
        ),
    )

    out_tsv = OUT / "biodarwin_96well_router_plate_v1.tsv"
    out_md = OUT / "BIODARWIN_96WELL_ROUTER_PLATE_V1.md"
    out_json = OUT / "biodarwin_96well_router_plate_v1_summary.json"
    plate.to_csv(out_tsv, sep="\t", index=False)

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "n_total": int(len(plate)),
        "n_core": int((plate["analysis_role"] == "locked_metric_control").sum() + (plate["analysis_role"] == "high_confidence_core").sum()),
        "n_explore": int((plate["analysis_role"] != "locked_metric_control").sum() - (plate["analysis_role"] == "high_confidence_core").sum()),
        "bucket_counts": plate["plate_bucket"].value_counts().to_dict(),
        "source_counts_top10": plate["source_id"].value_counts().head(10).to_dict(),
        "locked_metric_rows": int((plate["analysis_role"] == "locked_metric_control").sum()),
        "high_confidence_rows": int((plate["analysis_role"] == "high_confidence_core").sum()),
        "helper_lane_rows": int((plate["analysis_role"] == "helper_lane").sum()),
        "negative_control_rows": int((plate["analysis_role"] == "negative_control").sum()),
        "outputs": {
            "tsv": str(out_tsv),
            "md": str(out_md),
        },
    }
    out_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    report = []
    report.append("# BioDarwin 96-well router plate v1")
    report.append("")
    report.append(f"Generated: {summary['generated_at']}")
    report.append("")
    report.append("## Summary")
    report.append(f"- Total wells: {summary['n_total']}")
    report.append(f"- Locked metric rows: {summary['locked_metric_rows']}")
    report.append(f"- High-confidence core rows: {summary['high_confidence_rows']}")
    report.append(f"- Helper-lane rows: {summary['helper_lane_rows']}")
    report.append(f"- Negative-control rows: {summary['negative_control_rows']}")
    report.append("")
    report.append("## Bucket counts")
    report.append(pd.Series(summary["bucket_counts"]).to_frame("count").to_markdown())
    report.append("")
    report.append("## Top source counts")
    report.append(pd.Series(summary["source_counts_top10"]).to_frame("count").to_markdown())
    report.append("")
    report.append("## Claim boundary")
    report.append("- This is a wetlab routing manifest.")
    report.append("- Core rows are high-confidence public industrial rows.")
    report.append("- Extra rows are exploratory probes, helper-lane candidates, and specificity controls.")
    report.append("- No class-II locked metric claim is made; class-II rows are exploratory/helper only.")
    out_md.write_text("\n".join(report), encoding="utf-8")

    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
