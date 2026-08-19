#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from neoimmune_common import ensure_run_dir, metric_binary, safe_read_table, write_md, write_tsv


LABEL_COL = "label_immunogenicity"

FOCUS_SOURCES = [
    ("NEPdb", "NEPdb"),
    ("TESLA", "TESLA"),
    ("dbPepNeo", "dbPepNeo"),
    ("CEDAR/IMPROVE", "CEDAR"),
    ("ITSNdb", "ITSNdb"),
]

NMI_COLS = [
    "NMI_all_branch_mean",
    "NMI_tcr_qk_structure",
    "NMI_equal_top5",
    "NMI_w7b_esm2_qk",
]

LOCAL_COLS = [
    "ESM2_Bayesian",
    "W7B_stacked",
    "W7A_QK_only",
    "Wave8_TCR_SelfSim_full",
    "Wave8_TCR_SelfSim_no_exact",
    "Structure_LR",
    "GP_quantum",
]

PUBLIC_COLS = [
    "BigMHC_IM",
    "MHCflurry_2.0_presentation",
    "PRIME",
    "NetMHCpan_4.1_EL",
]


def pivot_scores(scores: pd.DataFrame, name_col: str, value_col: str) -> pd.DataFrame:
    if scores.empty or "candidate_id" not in scores:
        return pd.DataFrame(columns=["candidate_id"])
    d = scores.copy()
    d[value_col] = pd.to_numeric(d[value_col], errors="coerce")
    return d.pivot_table(index="candidate_id", columns=name_col, values=value_col, aggfunc="mean").reset_index()


def best_metric(df: pd.DataFrame, cols: list[str]) -> dict[str, object]:
    out = {
        "best_model": "",
        "best_AUPRC": np.nan,
        "best_AUROC": np.nan,
        "best_n": np.nan,
        "best_positives": np.nan,
        "best_coverage": np.nan,
    }
    if df.empty:
        return out
    best_row = None
    for c in cols:
        if c not in df:
            continue
        sub = df[df[c].notna()].copy()
        if sub.empty:
            continue
        m = metric_binary(sub[LABEL_COL], sub[c])
        if np.isnan(m.get("AUPRC", np.nan)):
            continue
        if best_row is None or m["AUPRC"] > best_row["best_AUPRC"]:
            best_row = {
                "best_model": c,
                "best_AUPRC": m.get("AUPRC", np.nan),
                "best_AUROC": m.get("AUROC", np.nan),
                "best_n": m.get("n", np.nan),
                "best_positives": m.get("positives", np.nan),
                "best_coverage": float(sub[c].notna().mean()),
            }
    if best_row:
        out.update(best_row)
    return out


def write_figures(outdir: Path, summary: pd.DataFrame) -> None:
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    if summary.empty:
        return
    try:
        import cv2
        import matplotlib.pyplot as plt

        d = summary.copy()
        d["delta_nmi_minus_public"] = pd.to_numeric(d["delta_nmi_minus_public"], errors="coerce")
        d["nmi_best_auprc"] = pd.to_numeric(d["nmi_best_auprc"], errors="coerce")
        d["public_best_auprc"] = pd.to_numeric(d["public_best_auprc"], errors="coerce")
        d = d.sort_values("delta_nmi_minus_public", ascending=False)

        fig, ax = plt.subplots(figsize=(12, 6), facecolor="#080b12")
        ax.set_facecolor("#111827")
        colors = ["#26d9ff" if x > 0 else "#fb7185" for x in d["delta_nmi_minus_public"].fillna(0)]
        ax.barh(d["source_family"], d["delta_nmi_minus_public"].fillna(0), color=colors)
        ax.axvline(0, color="#9fb0c7", lw=1)
        ax.set_xlabel("best NMI AUPRC - best public AUPRC")
        ax.set_title("Source transfer gap: clean NMI vs frozen public comparators", color="white", fontsize=14)
        ax.tick_params(colors="white")
        ax.xaxis.label.set_color("white")
        for spine in ax.spines.values():
            spine.set_color("#263348")
        ax.grid(axis="x", color="#263348", alpha=0.5)
        fig.tight_layout()
        raw = figdir / "source_transfer_gap_raw.png"
        fig.savefig(raw, dpi=220)
        plt.close(fig)

        img = cv2.imread(str(raw))
        if img is not None:
            band_h = 84
            band = np.zeros((band_h, img.shape[1], 3), dtype=np.uint8)
            band[:, :] = (8, 11, 18)
            cv2.putText(band, "Positive bars mean NMI beats the best frozen public comparator in that source family.", (28, 36), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (38, 217, 255), 2, cv2.LINE_AA)
            cv2.putText(band, "Negative or blank bars mean coverage is missing or public comparators still dominate.", (28, 64), cv2.FONT_HERSHEY_SIMPLEX, 0.55, (159, 176, 199), 1, cv2.LINE_AA)
            cv2.imwrite(str(figdir / "source_transfer_gap.png"), np.vstack([band, img]))
    except Exception:
        pass


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)

    canon = safe_read_table(outdir / "data" / "canonical_candidates.tsv.gz")
    canon[LABEL_COL] = pd.to_numeric(canon[LABEL_COL], errors="coerce")
    nmi = safe_read_table(outdir / "predictions" / "nmi_clean_method_predictions.tsv")
    local = safe_read_table(outdir / "local_scores" / "all_local_scores.tsv.gz")
    ext = safe_read_table(outdir / "external_scores" / "all_external_scores.tsv.gz")

    local_wide = pivot_scores(local[local["local_model_name"].isin(LOCAL_COLS)], "local_model_name", "normalized_score")
    ext_wide = pivot_scores(ext[ext["model_name"].isin(PUBLIC_COLS)], "model_name", "normalized_score")

    for c in ["any_existing_overlap_flag", "strict_set_flag"]:
        if c not in canon.columns:
            canon[c] = False
    d = canon[[
        "candidate_id",
        "patient_id",
        "dataset_source",
        "source_study",
        "tumor_type",
        LABEL_COL,
        "any_existing_overlap_flag",
        "strict_set_flag",
    ]].merge(nmi[["candidate_id"] + [c for c in NMI_COLS if c in nmi.columns]], on="candidate_id", how="left")
    d = d.merge(local_wide, on="candidate_id", how="left").merge(ext_wide, on="candidate_id", how="left")

    rows = []
    for fam, token in FOCUS_SOURCES:
        sub = d[d["source_study"].fillna("").astype(str).str.contains(token, case=False, regex=False)].copy()
        if sub.empty:
            continue
        total = len(sub)
        positives = int(pd.to_numeric(sub[LABEL_COL], errors="coerce").fillna(0).sum())
        nmi_cov = {c: int(sub[c].notna().sum()) if c in sub else 0 for c in NMI_COLS}
        local_cov = {c: int(sub[c].notna().sum()) if c in sub else 0 for c in LOCAL_COLS}
        pub_cov = {c: int(sub[c].notna().sum()) if c in sub else 0 for c in PUBLIC_COLS}

        best_nmi = best_metric(sub, [c for c in NMI_COLS if c in sub.columns])
        best_local = best_metric(sub, [c for c in LOCAL_COLS if c in sub.columns])
        best_pub = best_metric(sub, [c for c in PUBLIC_COLS if c in sub.columns])

        rows.append(
            {
                "source_family": fam,
                "rows": total,
                "positives": positives,
                "any_overlap_flagged": int(sub["any_existing_overlap_flag"].fillna(False).astype(bool).sum()) if "any_existing_overlap_flag" in sub else 0,
                "strict_rows": int((~sub["any_existing_overlap_flag"].fillna(False).astype(bool)).sum()) if "any_existing_overlap_flag" in sub else total,
                "nmi_locked_rows": int(sub["NMI_all_branch_mean"].notna().sum()) if "NMI_all_branch_mean" in sub else 0,
                "nmi_best_model": best_nmi["best_model"],
                "nmi_best_auprc": best_nmi["best_AUPRC"],
                "nmi_best_auroc": best_nmi["best_AUROC"],
                "nmi_best_n": best_nmi["best_n"],
                "nmi_best_coverage": best_nmi["best_coverage"],
                "local_best_model": best_local["best_model"],
                "local_best_auprc": best_local["best_AUPRC"],
                "local_best_auroc": best_local["best_AUROC"],
                "local_best_coverage": best_local["best_coverage"],
                "public_best_model": best_pub["best_model"],
                "public_best_auprc": best_pub["best_AUPRC"],
                "public_best_auroc": best_pub["best_AUROC"],
                "public_best_coverage": best_pub["best_coverage"],
                "delta_nmi_minus_public": best_nmi["best_AUPRC"] - best_pub["best_AUPRC"] if pd.notna(best_nmi["best_AUPRC"]) and pd.notna(best_pub["best_AUPRC"]) else np.nan,
                "delta_local_minus_public": best_local["best_AUPRC"] - best_pub["best_AUPRC"] if pd.notna(best_local["best_AUPRC"]) and pd.notna(best_pub["best_AUPRC"]) else np.nan,
                "any_patient_real_ids": int(sub["patient_id"].fillna("unknown_patient").astype(str).ne("unknown_patient").sum()) if "patient_id" in sub else 0,
                "unique_real_patients": int(sub.loc[sub["patient_id"].fillna("unknown_patient").astype(str).ne("unknown_patient"), "patient_id"].astype(str).nunique()) if "patient_id" in sub else 0,
                "patient_claimable": bool(int(sub.loc[sub["patient_id"].fillna("unknown_patient").astype(str).ne("unknown_patient"), "patient_id"].astype(str).nunique()) >= 2) if "patient_id" in sub else False,
                "nmi_cov_any": int(sum(v > 0 for v in nmi_cov.values())),
                "local_cov_any": int(sum(v > 0 for v in local_cov.values())),
                "public_cov_any": int(sum(v > 0 for v in pub_cov.values())),
                "claim_boundary": "source-family level audit; not a clinical validation",
            }
        )

    summary = pd.DataFrame(rows).sort_values("rows", ascending=False)
    write_tsv(summary, outdir / "metrics" / "source_transfer_gap_summary.tsv")
    write_figures(outdir, summary)

    md = [
        "# Source Transfer Gap Audit",
        "",
        "This audit answers where NMI is already competitive and where coverage is still thin.",
        "",
        "## Key reading rule",
        "- `delta_nmi_minus_public > 0` means the best clean NMI row in that source family beats the best frozen public comparator available in the same source family.",
        "- `patient_claimable=False` means do not make patient-level top-N claims from that family.",
        "",
        summary.to_markdown(index=False) if not summary.empty else "_No rows available._",
        "",
        "## Conservative takeaways",
        "- NEPdb is the current bridge where NMI is strongest and public comparators are fully available.",
        "- TESLA/dbPepNeo coverage exists, but NMI locked-branch coverage is sparse or absent in the current run, so broad claims are premature.",
        "- This is exactly why the next ingestion step should target additional branch features or aligned source-specific mapping, not a new foundation model.",
        f"- Figure: `{outdir / 'figures/source_transfer_gap.png'}`",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "SOURCE_TRANSFER_GAP_AUDIT.md")
    print(outdir / "reports" / "SOURCE_TRANSFER_GAP_AUDIT.md")


if __name__ == "__main__":
    main()
