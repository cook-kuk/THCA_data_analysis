#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import norm

from pilot_utils import fdr_bh, pilot_root, setup_logging


EXPECTED = {
    ("ATC_vs_PTC", "RAI differentiation"): "down",
    ("ATC_vs_normal", "RAI differentiation"): "down",
    ("PTC_vs_normal", "RAI differentiation"): "down",
    ("PDTC_vs_PTC", "RAI differentiation"): "down",
    ("ATC_vs_PTC", "RAI-restorable proxy"): "down",
    ("ATC_vs_normal", "RAI-restorable proxy"): "down",
    ("ATC_vs_PTC", "Aggressive dedifferentiation"): "up",
    ("ATC_vs_normal", "Aggressive dedifferentiation"): "up",
    ("PDTC_vs_PTC", "Aggressive dedifferentiation"): "up",
    ("ATC_vs_PTC", "Drug-delivery failure proxy"): "up",
    ("ATC_vs_normal", "Drug-delivery failure proxy"): "up",
    ("ATC_vs_PTC", "Proliferation"): "up",
    ("ATC_vs_normal", "Proliferation"): "up",
    ("ATC_vs_PTC", "Hypoxia"): "up",
    ("ATC_vs_normal", "Hypoxia"): "up",
    ("ATC_vs_PTC", "Myeloid/CAF barrier"): "up",
    ("ATC_vs_normal", "Myeloid/CAF barrier"): "up",
}


def signed_z(effect: float, p: float, expected: str) -> float:
    if not np.isfinite(effect) or not np.isfinite(p) or p <= 0 or p > 1:
        return np.nan
    p = max(p, 1e-300)
    z = norm.isf(p / 2)
    if expected == "down":
        sign = 1 if effect < 0 else -1
    elif expected == "up":
        sign = 1 if effect > 0 else -1
    else:
        sign = 1 if effect > 0 else -1
    return sign * z


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=pilot_root())
    args = parser.parse_args()
    logger = setup_logging("21_external_meta_consistency")
    out = args.root / "results" / "extra_analyses"
    tables = out / "tables"
    figs = out / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)
    path = args.root / "results" / "public_validation_expansion" / "tables" / "exact_kthyro_external_bulk_contrasts.tsv"
    df = pd.read_csv(path, sep="\t")
    rows = []
    for (contrast, axis), expected in EXPECTED.items():
        sub = df[df["contrast"].eq(contrast) & df["axis_label"].eq(axis)].copy()
        if sub.empty:
            continue
        effects = pd.to_numeric(sub["effect_size_cohens_d"], errors="coerce")
        pvals = pd.to_numeric(sub["p_value"], errors="coerce")
        direction_match = (effects < 0) if expected == "down" else (effects > 0)
        zvals = [signed_z(e, p, expected) for e, p in zip(effects, pvals)]
        zvals = np.array([z for z in zvals if np.isfinite(z)])
        stouffer_z = zvals.sum() / np.sqrt(len(zvals)) if len(zvals) else np.nan
        stouffer_p = norm.sf(stouffer_z) if np.isfinite(stouffer_z) else np.nan
        rows.append(
            {
                "contrast": contrast,
                "axis_label": axis,
                "expected_direction": expected,
                "n_datasets": int(sub["dataset"].nunique()),
                "datasets": ",".join(sorted(sub["dataset"].unique())),
                "median_effect_cohens_d": float(effects.median()),
                "min_effect": float(effects.min()),
                "max_effect": float(effects.max()),
                "n_direction_match": int(direction_match.sum()),
                "direction_match_fraction": float(direction_match.mean()),
                "n_strong_or_moderate": int(sub["support_strength"].isin(["strong", "moderate"]).sum()),
                "stouffer_signed_z_expected_direction": stouffer_z,
                "stouffer_one_sided_p_expected_direction": stouffer_p,
                "claim_boundary": "Cross-cohort direction consistency supports reproducibility of expression axes, not therapy response.",
            }
        )
    outdf = pd.DataFrame(rows)
    outdf["stouffer_fdr_q"] = fdr_bh(outdf["stouffer_one_sided_p_expected_direction"])
    outdf["support_summary"] = np.select(
        [
            outdf["direction_match_fraction"].ge(0.75) & outdf["stouffer_fdr_q"].lt(0.05),
            outdf["direction_match_fraction"].ge(0.5) & outdf["stouffer_fdr_q"].lt(0.10),
            outdf["direction_match_fraction"].lt(0.5),
        ],
        ["consistent_strong", "consistent_moderate", "inconsistent"],
        default="weak_or_underpowered",
    )
    outdf.to_csv(tables / "external_exact_kthyro_meta_consistency.tsv", sep="\t", index=False)

    plot = outdf.sort_values(["support_summary", "direction_match_fraction", "median_effect_cohens_d"], ascending=[True, False, True])
    fig, ax = plt.subplots(figsize=(12.5, 7.5))
    sns.barplot(data=plot, x="direction_match_fraction", y="axis_label", hue="contrast", ax=ax)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of external cohorts matching expected direction")
    ax.set_ylabel("")
    ax.set_title("External cohort direction consistency for K-Thyro axes", weight="bold")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(figs / "external_exact_kthyro_meta_direction_consistency.png", dpi=300, bbox_inches="tight")
    fig.savefig(figs / "external_exact_kthyro_meta_direction_consistency.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("Wrote %d external meta-consistency rows.", len(outdf))


if __name__ == "__main__":
    main()
