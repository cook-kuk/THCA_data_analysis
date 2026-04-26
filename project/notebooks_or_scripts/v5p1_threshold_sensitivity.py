"""
v5.1 audit F4: Threshold sensitivity analysis (+/- 0.05) for DIAL interpretation.

Re-classifies the 25 rows of v5p1_dial_all_cancers.tsv under a grid of
(dial_high, dial_low, auc_post_high, auc_post_low) thresholds and counts
label changes vs. the v5.1b baseline (0.30, 0.10, 0.70, 0.60).

For each perturbation we report:
  - total number of rows that change interpretation
  - number of THCA batch_entangled flips that downgrade
  - number of non-THCA rows that newly become batch_entangled
  - the list of (cancer, classifier) cells that changed

Output:
  results/v5/v5p1_threshold_sensitivity.tsv

Run:
  python -u notebooks_or_scripts/v5p1_threshold_sensitivity.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path("/opt/thyroid-dash/project")
DIAL_TSV = PROJECT / "results" / "v5" / "v5p1_dial_all_cancers.tsv"
OUT_TSV = PROJECT / "results" / "v5" / "v5p1_threshold_sensitivity.tsv"

# True-biology branch uses an extra dial_strict threshold (existing v5.1
# rule: auc_post > auc_post_high AND dial < 0.05). We hold this fixed at 0.05
# because the audit only varies the 4 named thresholds.
DIAL_STRICT = 0.05


def classify(dial_val: float, auc_post: float,
             dial_high: float, dial_low: float,
             auc_post_high: float, auc_post_low: float) -> str:
    """Apply the v5.1 interpretation rule with parameterised thresholds.

    Mirrors notebooks_or_scripts/v5p1_common.py classification logic.
    """
    if (dial_val is None or auc_post is None
            or (isinstance(dial_val, float) and math.isnan(dial_val))
            or (isinstance(auc_post, float) and math.isnan(auc_post))):
        return "ambiguous"
    if dial_val > dial_high:
        return "batch_entangled"
    if dial_val > dial_low:
        return "partial_batch"
    if auc_post > auc_post_high and dial_val < DIAL_STRICT:
        return "true_biology"
    if auc_post < auc_post_low:
        return "no_signal"
    return "ambiguous"


def relabel(df: pd.DataFrame, dh: float, dl: float, ah: float, al: float) -> pd.Series:
    return df.apply(
        lambda r: classify(float(r["dial"]), float(r["auc_post"]), dh, dl, ah, al),
        axis=1,
    )


def main() -> int:
    print(f"[F4] reading {DIAL_TSV}", flush=True)
    df = pd.read_csv(DIAL_TSV, sep="\t")
    print(f"[F4] loaded {len(df)} rows, columns: {list(df.columns)}", flush=True)

    baseline = (0.30, 0.10, 0.70, 0.60)
    perturbations = [
        ("baseline_v5.1b", 0.30, 0.10, 0.70, 0.60),
        ("dial_high-0.05", 0.25, 0.10, 0.70, 0.60),
        ("dial_high+0.05", 0.35, 0.10, 0.70, 0.60),
        ("dial_low-0.05",  0.30, 0.05, 0.70, 0.60),
        ("dial_low+0.05",  0.30, 0.15, 0.70, 0.60),
        ("auc_post_high-0.05", 0.30, 0.10, 0.65, 0.60),
        ("auc_post_high+0.05", 0.30, 0.10, 0.75, 0.60),
        ("auc_post_low-0.05",  0.30, 0.10, 0.70, 0.55),
        ("auc_post_low+0.05",  0.30, 0.10, 0.70, 0.65),
    ]

    base_labels = relabel(df, *baseline)
    # Sanity check: baseline should reproduce the labels in the table.
    file_labels = df["interpretation"].astype(str)
    # Note: file uses {batch_entangled, true_biology, no_signal, ambiguous}.
    # Our reclassifier may emit "partial_batch" (an additional bin in the
    # original code that does not appear in the file because no row has
    # 0.10 < dial <= 0.30). Confirm equality of labels that do appear.
    n_match = int((base_labels.values == file_labels.values).sum())
    print(f"[F4] baseline reclassification matches file in {n_match}/{len(df)} rows",
          flush=True)
    if n_match != len(df):
        diff_idx = np.where(base_labels.values != file_labels.values)[0]
        for i in diff_idx:
            print(f"  diff: cancer={df.iloc[i]['cancer']} clf={df.iloc[i]['classifier']} "
                  f"file={file_labels.iloc[i]} recomputed={base_labels.iloc[i]} "
                  f"dial={df.iloc[i]['dial']:.4f} auc_post={df.iloc[i]['auc_post']:.4f}",
                  flush=True)

    out_rows = []
    for name, dh, dl, ah, al in perturbations:
        new_labels = relabel(df, dh, dl, ah, al)
        diff_mask = (new_labels.values != base_labels.values)
        n_changed = int(diff_mask.sum())

        # THCA batch_entangled downgrade: rows where baseline=batch_entangled
        # and cancer=THCA but new label is anything else.
        thca_mask = (df["cancer"].values == "THCA")
        thca_be_base = thca_mask & (base_labels.values == "batch_entangled")
        thca_be_lost = int((thca_be_base & (new_labels.values != "batch_entangled")).sum())

        # Non-THCA newly batch_entangled.
        non_thca_mask = ~thca_mask
        non_thca_be_gained = int(
            (non_thca_mask
             & (base_labels.values != "batch_entangled")
             & (new_labels.values == "batch_entangled")).sum()
        )

        changed_pairs = []
        for i in np.where(diff_mask)[0]:
            changed_pairs.append(
                f"{df.iloc[i]['cancer']}/{df.iloc[i]['classifier']}:"
                f"{base_labels.iloc[i]}->{new_labels.iloc[i]}"
            )
        out_rows.append({
            "perturbation": name,
            "dial_high": dh,
            "dial_low": dl,
            "auc_post_high": ah,
            "auc_post_low": al,
            "total_label_changes": n_changed,
            "n_THCA_BE_lost": thca_be_lost,
            "n_nonTHCA_BE_gained": non_thca_be_gained,
            "list_of_changed_rows": "; ".join(changed_pairs) if changed_pairs else "",
        })
        print(f"[F4] {name:24s} (dh={dh}, dl={dl}, ah={ah}, al={al})"
              f"  changes={n_changed}  THCA_BE_lost={thca_be_lost}"
              f"  nonTHCA_BE_gained={non_thca_be_gained}", flush=True)
        if changed_pairs:
            for cp in changed_pairs:
                print(f"    - {cp}", flush=True)

    out_df = pd.DataFrame(out_rows)[
        ["perturbation", "dial_high", "dial_low", "auc_post_high", "auc_post_low",
         "total_label_changes", "n_THCA_BE_lost", "n_nonTHCA_BE_gained",
         "list_of_changed_rows"]
    ]
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"[F4] wrote {OUT_TSV}", flush=True)

    # Summary
    non_baseline = out_df[out_df["perturbation"] != "baseline_v5.1b"]
    max_changes = int(non_baseline["total_label_changes"].max())
    max_thca_lost = int(non_baseline["n_THCA_BE_lost"].max())
    any_nonthca_be = int(non_baseline["n_nonTHCA_BE_gained"].sum())
    print(f"\n[F4-SUMMARY] max single-perturbation label changes: {max_changes}",
          flush=True)
    print(f"[F4-SUMMARY] max THCA batch_entangled downgrades in any perturbation: "
          f"{max_thca_lost}", flush=True)
    print(f"[F4-SUMMARY] total non-THCA rows newly batch_entangled (across all perts): "
          f"{any_nonthca_be}", flush=True)
    print(f"[F4-SUMMARY] THCA flip count robust = "
          f"{'YES' if max_thca_lost == 0 else 'NO'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
