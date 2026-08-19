#!/usr/bin/env python3
"""Interpret CROSS-Neo assay result sheets against preregistered endpoints.

The script intentionally works from candidate-level PASS/FAIL calls. Raw
instrument values and lab-specific calibration should be resolved upstream into
the result-entry fields before this endpoint-level interpreter is run.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import pandas as pd


PASS_VALUES = {"PASS", "POS", "POSITIVE", "YES", "Y", "TRUE", "1", "SUCCESS", "UNLOCK"}
FAIL_VALUES = {"FAIL", "NEG", "NEGATIVE", "NO", "N", "FALSE", "0", "MISS", "LOCK"}
PENDING_VALUES = {"", "NA", "N/A", "NAN", "NONE", "PENDING", "TBD", "UNKNOWN"}
EXCLUDE_VALUES = {"EXCLUDE", "EXCLUDED", "QC_FAIL", "INVALID"}


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def normalize_call(value: object) -> str:
    text = str(value).strip().upper()
    if text in PASS_VALUES:
        return "PASS"
    if text in FAIL_VALUES:
        return "FAIL"
    if text in EXCLUDE_VALUES:
        return "EXCLUDE"
    if text in PENDING_VALUES:
        return "PENDING"
    return "UNRECOGNIZED"


def first_present(row: pd.Series, names: list[str]) -> object:
    for name in names:
        if name in row.index:
            value = row.get(name, "")
            if str(value).strip() != "":
                return value
    return ""


def derive_candidate_call(row: pd.Series) -> tuple[str, str]:
    manual = normalize_call(first_present(row, ["candidate_success_call", "observed_candidate_success"]))
    if manual in {"PASS", "FAIL", "EXCLUDE"}:
        return manual, "manual_candidate_success_call"
    if manual == "UNRECOGNIZED":
        return "PENDING", "manual_candidate_success_call_unrecognized"

    arm = str(row.get("plate_v4_arm", ""))
    mut = normalize_call(first_present(row, ["mut_pmhc_call", "mutant_pmhc_call", "binding_call"]))
    specificity = normalize_call(first_present(row, ["specificity_pass_call", "observed_specificity_pass"]))
    tcr_or_prov = normalize_call(
        first_present(row, ["tcr_or_provenance_pass_call", "observed_tcr_or_provenance_pass", "tcr_readout_call", "provenance_call"])
    )
    qc = normalize_call(first_present(row, ["technical_qc_call", "qc_call"]))
    weak_negative = normalize_call(first_present(row, ["hard_negative_weak_call", "negative_remains_weak_call"]))

    if "D_specificity_moat" == arm:
        if weak_negative in {"PASS", "FAIL", "EXCLUDE"}:
            return weak_negative, "hard_negative_weak_call"
        if mut == "FAIL" and specificity in {"PASS", "PENDING"}:
            return "PASS", "derived_negative_remains_weak"
        return "PENDING", "needs_hard_negative_weak_call"

    if "E_positive_QC_control" == arm:
        if mut == "PASS" and qc in {"PASS", "PENDING"}:
            return "PASS", "derived_positive_control_binding"
        if mut == "FAIL":
            return "FAIL", "derived_positive_control_binding"
        return "PENDING", "needs_positive_control_binding"

    if "B_mechanism_TCR_MD" == arm:
        if mut == "PASS" and tcr_or_prov == "PASS":
            return "PASS", "derived_binding_plus_tcr"
        if mut == "FAIL" or tcr_or_prov == "FAIL":
            return "FAIL", "derived_binding_plus_tcr"
        return "PENDING", "needs_binding_and_tcr_readout"

    if "C_label_rescue" == arm:
        if mut == "PASS" and tcr_or_prov == "PASS":
            return "PASS", "derived_binding_plus_provenance"
        if mut == "FAIL" or tcr_or_prov == "FAIL":
            return "FAIL", "derived_binding_plus_provenance"
        return "PENDING", "needs_binding_and_provenance"

    if "F_model_boundary" == arm:
        if tcr_or_prov in {"PASS", "FAIL", "EXCLUDE"}:
            return tcr_or_prov, "manual_tcr_structure_resolution_call"
        return "PENDING", "needs_tcr_structure_resolution"

    if mut == "PASS" and specificity in {"PASS", "PENDING"}:
        return "PASS", "derived_mutant_binding_specificity"
    if mut == "FAIL" or specificity == "FAIL":
        return "FAIL", "derived_mutant_binding_specificity"
    return "PENDING", "needs_binding_specificity"


def endpoint_status(successes: int, completed: int, expected: int, threshold: int, direction: str) -> str:
    if completed == 0:
        return "PENDING_NO_DATA"
    if direction == "at_most":
        if completed < expected:
            if successes <= threshold and successes + (expected - completed) <= threshold:
                return "UNLOCKED_EVEN_IF_PENDING_FAILS"
            return "PENDING"
        return "UNLOCKED" if successes <= threshold else "NOT_UNLOCKED"

    if successes >= threshold:
        return "UNLOCKED"
    if completed < expected:
        remaining = expected - completed
        if successes + remaining < threshold:
            return "NOT_UNLOCKED_EVEN_IF_PENDING_PASS"
        return "PENDING"
    return "NOT_UNLOCKED"


def interpret(results: pd.DataFrame, endpoint: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    calls = results.copy()
    derived = calls.apply(derive_candidate_call, axis=1, result_type="expand")
    calls["normalized_candidate_success"] = derived[0]
    calls["candidate_call_source"] = derived[1]
    calls["endpoint_countable"] = calls["normalized_candidate_success"].isin(["PASS", "FAIL"])
    calls["endpoint_success"] = calls["normalized_candidate_success"].eq("PASS")

    rows = []
    for _, plan in endpoint.iterrows():
        arm = str(plan["plate_v4_arm"])
        arm_calls = calls[calls["plate_v4_arm"].astype(str).eq(arm)]
        expected = int(plan["n_candidates"])
        successes = int(arm_calls["endpoint_success"].sum())
        completed = int(arm_calls["endpoint_countable"].sum())
        excluded = int(arm_calls["normalized_candidate_success"].eq("EXCLUDE").sum())
        threshold = int(plan["threshold_successes"])
        confirmatory_threshold = int(plan.get("confirmatory_threshold_successes", threshold))
        direction = str(plan.get("threshold_direction", "at_least"))
        rows.append(
            {
                "plate_v4_arm": arm,
                "endpoint": plan.get("endpoint", ""),
                "expected_candidates": expected,
                "observed_candidates": int(len(arm_calls)),
                "completed_candidates": completed,
                "excluded_candidates": excluded,
                "success_count": successes,
                "threshold_successes": threshold,
                "confirmatory_threshold_successes": confirmatory_threshold,
                "exploratory_unlock_status": endpoint_status(successes, completed, expected, threshold, direction),
                "confirmatory_unlock_status": endpoint_status(successes, completed, expected, confirmatory_threshold, direction),
                "claim_layer_after_unlock": plan.get("claim_layer_after_unlock", ""),
                "main_text_claim_tier": plan.get("main_text_claim_tier", ""),
                "interpretation_if_unlocked": plan.get("interpretation_if_unlocked", ""),
                "interpretation_if_not_unlocked": plan.get("interpretation_if_not_unlocked", ""),
            }
        )
    endpoint_calls = pd.DataFrame(rows)
    return calls, endpoint_calls


def write_markdown(endpoint_calls: pd.DataFrame, out_path: Path) -> None:
    lines = [
        "# CROSS-Neo v6 assay endpoint calls",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "This file is an endpoint-level interpretation of candidate PASS/FAIL calls against the v5 preregistered thresholds.",
        "",
        "## Endpoint calls",
        "",
    ]
    lines.append(endpoint_calls.to_markdown(index=False))
    lines.extend(
        [
            "",
            "## Claim boundary",
            "",
            "Unlocked exploratory endpoints support the matching assay/diagnostic layer only. Confirmatory-ready claims require the stricter confirmatory threshold and still do not create a clinical vaccine-selection claim.",
            "",
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path, help="Candidate-level result-entry TSV.")
    parser.add_argument("--endpoint-plan", required=True, type=Path, help="v5 preregistered endpoint-plan TSV.")
    parser.add_argument("--out-dir", required=True, type=Path, help="Directory for interpreted calls.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    results = read_tsv(args.results)
    endpoint = read_tsv(args.endpoint_plan)

    required = {"plate_v4_arm", "candidate_id"}
    missing = sorted(required - set(results.columns))
    if missing:
        raise ValueError(f"Result sheet missing required columns: {missing}")

    candidate_calls, endpoint_calls = interpret(results, endpoint)
    candidate_path = args.out_dir / "candidate_calls_v6.tsv"
    endpoint_path = args.out_dir / "endpoint_calls_v6.tsv"
    report_path = args.out_dir / "claim_unlock_calls_v6.md"
    summary_path = args.out_dir / "interpretation_summary_v6.json"
    candidate_calls.to_csv(candidate_path, sep="\t", index=False)
    endpoint_calls.to_csv(endpoint_path, sep="\t", index=False)
    write_markdown(endpoint_calls, report_path)
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "results_input": str(args.results),
        "endpoint_plan_input": str(args.endpoint_plan),
        "candidate_count": int(len(candidate_calls)),
        "endpoint_count": int(len(endpoint_calls)),
        "exploratory_unlocked": int(endpoint_calls["exploratory_unlock_status"].eq("UNLOCKED").sum()),
        "confirmatory_unlocked": int(endpoint_calls["confirmatory_unlock_status"].eq("UNLOCKED").sum()),
        "pending_endpoints": int(endpoint_calls["exploratory_unlock_status"].str.startswith("PENDING").sum()),
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
