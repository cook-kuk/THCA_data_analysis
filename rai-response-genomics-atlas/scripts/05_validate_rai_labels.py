#!/usr/bin/env python3
"""05 — statistical validation of panel score against parsed RAI labels.

Usage:  python3 scripts/05_validate_rai_labels.py <ACCESSION>

Inputs:
  data/processed/<ACCESSION>_panel_score.tsv  (from 04)
  data/interim/<ACCESSION>_metadata.tsv       (from 03)

Outputs:
  data/processed/<ACCESSION>_label_joined.tsv
  results/tables/<ACCESSION>_test_results.tsv
  results/reports/<ACCESSION>_validation_brief.md
"""
from __future__ import annotations
import re
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, kruskal
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
INTERIM = ROOT / "data" / "interim"
PROCESSED = ROOT / "data" / "processed"
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "results" / "reports"


def normalize_label(s: str) -> str:
    s = (s or "").lower()
    if "refractor" in s and "avid" in s:
        return "mixed"
    if "refractor" in s:
        return "refractory"
    if "avid" in s:
        return "avid"
    if "no uptake" in s or "non-avid" in s or "no-avid" in s:
        return "non_avid"
    if "remission" in s:
        return "remission"
    if "persisten" in s:
        return "persistent"
    return ""


def derive_label(meta: pd.DataFrame) -> pd.Series:
    """Prefer structured parsed columns (rai_response, disease_status, rai_uptake_met) over blob."""
    # Priority 1: rai_response column from 03's structured parser
    if "rai_response" in meta.columns:
        s = meta["rai_response"].fillna("").astype(str).str.strip().str.lower()
        # GSE299988 uses 'refractive', GSE151179 uses 'refractory'; normalize both.
        s = s.replace({"refractive": "refractory", "non-avid": "non_avid", "non avid": "non_avid"})
        s = s.where(s.isin(["avid", "refractory", "non_avid"]), other="")
        if s.replace("", pd.NA).dropna().shape[0] > 0:
            return s.rename("rai_label")
    # Fallback: blob scan
    char_cols = [c for c in meta.columns if "characteristics" in c.lower() or c.endswith("_joined") or c == "rai_label_blob"]
    if not char_cols:
        return pd.Series([""] * len(meta), index=meta.index, name="rai_label")
    blob = meta[char_cols].astype(str).agg(" | ".join, axis=1)
    return blob.map(normalize_label).rename("rai_label")


def derive_secondary_labels(meta: pd.DataFrame) -> pd.DataFrame:
    """Pull rai_uptake_met, disease_status, sample_type, lesion_class as analysis columns."""
    cols = ["rai_uptake_met", "disease_status", "sample_type", "collection_timing",
            "lesion_class", "lesion_driver", "patient_id", "histology_variant"]
    out = pd.DataFrame(index=meta.index)
    for c in cols:
        if c in meta.columns:
            out[c] = meta[c].fillna("").astype(str).str.strip().str.lower()
    return out


def two_group_test(a, b, name_test, name_ref):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) < 3 or len(b) < 3:
        return None
    stat, p = mannwhitneyu(a, b, alternative="two-sided")
    y = np.r_[np.zeros(len(b)), np.ones(len(a))]
    s = np.r_[-b, -a]   # convention: lower panel_z => "test" group
    try:
        auc = roc_auc_score(y, s)
    except Exception:
        auc = float("nan")
    # Cohen's d
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / max(len(a)+len(b)-2, 1))
    d = (a.mean() - b.mean()) / pooled if pooled > 0 else float("nan")
    return {
        "n_test": int(len(a)), "n_ref": int(len(b)),
        "median_test": float(np.median(a)), "median_ref": float(np.median(b)),
        "cohen_d_test_minus_ref": round(float(d), 3),
        "mwu_p": float(p),
        "AUC_test_low_panel": float(auc),
    }


def main() -> None:
    if len(sys.argv) < 2:
        print("usage: 05_validate_rai_labels.py <ACCESSION>", file=sys.stderr); sys.exit(1)
    acc = sys.argv[1].strip()
    panel = pd.read_csv(PROCESSED / f"{acc}_panel_score.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(INTERIM / f"{acc}_metadata.tsv", sep="\t", index_col=0)
    joined = panel.join(meta, how="left")
    joined["rai_label"] = derive_label(joined)
    secondary = derive_secondary_labels(joined)
    for c in secondary.columns:
        if c not in joined.columns:
            joined[c] = secondary[c]

    out_joined = PROCESSED / f"{acc}_label_joined.tsv"
    joined.to_csv(out_joined, sep="\t")
    print(f"# wrote {out_joined}")

    tests = []

    # 1) Sanity: tumor vs non-neoplastic / normal
    NORMAL_PATTERN = r"non-neoplastic|^normal$|adjacent"
    if "sample_type" in joined.columns:
        st = joined["sample_type"].fillna("").astype(str).str.lower()
        is_normal = st.str.contains(NORMAL_PATTERN, regex=True)
        is_tumor  = ~is_normal & st.ne("")
        if is_normal.sum() >= 3 and is_tumor.sum() >= 3:
            r = two_group_test(joined.loc[is_tumor, "panel_z"], joined.loc[is_normal, "panel_z"], "tumor", "non_neoplastic")
            if r:
                tests.append({"comparison": "tumor_vs_non_neoplastic", **r})

    # 2) Within-tumor primary: avid vs refractory
    if "sample_type" in joined.columns:
        tumor_only = joined[~joined["sample_type"].fillna("").str.lower().str.contains(NORMAL_PATTERN, regex=True)]
    else:
        tumor_only = joined
    g = tumor_only.groupby("rai_label")["panel_z"].apply(list).to_dict()
    g = {k: np.asarray(v, dtype=float) for k, v in g.items() if k and len(v) >= 2}
    print(f"# tumor-only label groups (n>=2): { {k: len(v) for k,v in g.items()} }")
    if "avid" in g and "refractory" in g:
        r = two_group_test(g["refractory"], g["avid"], "refractory", "avid")
        if r:
            tests.append({"comparison": "refractory_vs_avid_within_tumor", **r})

    # 3) Uptake yes vs no
    if "rai_uptake_met" in tumor_only.columns:
        uy = tumor_only.loc[tumor_only["rai_uptake_met"] == "yes", "panel_z"]
        un = tumor_only.loc[tumor_only["rai_uptake_met"] == "no", "panel_z"]
        r = two_group_test(un, uy, "no_uptake", "uptake")
        if r:
            tests.append({"comparison": "uptake_no_vs_yes_within_tumor", **r})

    # 4) Persistence vs remission
    if "disease_status" in tumor_only.columns:
        per = tumor_only.loc[tumor_only["disease_status"] == "persistence", "panel_z"]
        rem = tumor_only.loc[tumor_only["disease_status"] == "remission", "panel_z"]
        r = two_group_test(per, rem, "persistence", "remission")
        if r:
            tests.append({"comparison": "persistence_vs_remission_within_tumor", **r})

    # 5) Per-driver panel score within tumor
    if "lesion_class" in tumor_only.columns:
        drivers = ["brafv600e", "fusion", "ptert", "wt"]
        present = [d for d in drivers if (tumor_only["lesion_class"] == d).sum() >= 3]
        if len(present) >= 2:
            ks_arrays = [tumor_only.loc[tumor_only["lesion_class"] == d, "panel_z"].dropna().values for d in present]
            stat, p = kruskal(*ks_arrays)
            tests.append({
                "comparison": "kruskal_by_lesion_class_tumor_only",
                "n_groups": len(present),
                "labels": ",".join(f"{d}:{len(a)}" for d, a in zip(present, ks_arrays)),
                "stat": float(stat),
                "p": float(p),
            })

    # 6) Multi-class Kruskal on rai_label (tumor only)
    if len(g) >= 3:
        stat, p = kruskal(*g.values())
        tests.append({"comparison": "kruskal_all_rai_labels_tumor_only",
                       "n_groups": len(g),
                       "labels": ",".join(f"{k}:{len(v)}" for k, v in g.items()),
                       "stat": float(stat),
                       "p": float(p)})

    if not tests:
        print("WARN: not enough samples per label group for tests; check 03 parsing")
    TABLES.mkdir(parents=True, exist_ok=True)
    test_df = pd.DataFrame(tests)
    test_df.to_csv(TABLES / f"{acc}_test_results.tsv", sep="\t", index=False)
    print(f"# wrote {TABLES / (acc + '_test_results.tsv')}")

    # brief report
    REPORTS.mkdir(parents=True, exist_ok=True)
    lines = [
        f"# {acc} — RAI panel validation brief",
        "",
        f"- n samples merged: {joined.shape[0]}",
        f"- panel genes available per sample: {int(joined['n_genes_available'].iloc[0]) if 'n_genes_available' in joined.columns else 'n/a'}",
        f"- label group counts: { {k: len(v) for k,v in g.items()} }",
        "",
        "## Tests",
        "",
    ]
    if test_df.empty:
        lines.append("(no group had sufficient samples for tests)")
    else:
        lines.append(test_df.to_markdown(index=False))
    lines += [
        "",
        "## Wording rule",
        "",
        "Any claim derived from this report must declare the dataset tier "
        "(Tier 2 RAI avidity for GSE151179 / GSE299988; Tier 1 for Boucai 2023 if available). ",
        "Use 'associated with' or 'stratifies', never 'predicts RAI response'.",
    ]
    (REPORTS / f"{acc}_validation_brief.md").write_text("\n".join(lines))
    print(f"# wrote {REPORTS / (acc + '_validation_brief.md')}")


if __name__ == "__main__":
    main()
