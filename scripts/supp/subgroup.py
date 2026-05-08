"""
Supp #4 — Subgroup AUC analysis.

Joins predictions TSV with sample_master_v3.tsv (TCGA-THCA clinical) on
12-char patient barcode (submitter_id[:12]).
Subgroups:
  - AJCC stage group
  - driver_anchor (BRAF / RAS / unknown / other)
  - molecular_subtype (BRAF_like / RAS_like)
  - histology_subtype (cPTC / FVPTC / other)
  - tds_group (high / mid / low)
  - sex
  - tissue source site (TCGA-XX-YYYY → XX)

AUC + 1000-bootstrap 95% CI per subgroup with N>=5 and >=2 classes.

Saves:
  analysis_supp/subgroup_aucs.tsv
  analysis_supp/subgroup_aucs.json
  analysis_supp/figures/figS3_subgroup_forest.{png,pdf}
"""
from __future__ import annotations
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

PRED = Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/clam_per_slide_predictions.tsv")
MAN = Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/phase2_tcga_clam/slide_manifest.tsv")
SAMPLE_MASTER = Path("project/metadata/sample_master_v3.tsv")
OUT = Path("project/results/p2_image_dm1_v2_foundation_clam_2026_05_07/analysis_supp")
FIGS = OUT / "figures"

N_BOOT = 1000
SEED = 42


def boot_auc(y, p, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    n = len(y)
    if n < 2 or len(np.unique(y)) < 2:
        return None, None, None
    point = float(roc_auc_score(y, p))
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb = y[idx]
        pb = p[idx]
        if len(np.unique(yb)) < 2:
            continue
        aucs.append(roc_auc_score(yb, pb))
    if not aucs:
        return point, None, None
    return point, float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def main() -> None:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    pred = pd.read_csv(PRED, sep="\t")
    man = pd.read_csv(MAN, sep="\t")
    sm = pd.read_csv(SAMPLE_MASTER, sep="\t")

    df = pred.merge(man[["file_id", "submitter_id", "case_id"]],
                    left_on="slide", right_on="file_id", how="left")
    df["patient"] = df["submitter_id"].str[:12]
    df["tss"] = df["submitter_id"].str.slice(5, 7)

    # Map to sample_master clinical features. Dedup by patient (mean float, first cat).
    tcga = sm[sm["sample_id"].astype(str).str.startswith("TCGA-")].copy()
    tcga["patient"] = tcga["sample_id"].str[:12]
    cols_needed = [
        "patient", "ajcc_stage_group", "driver_anchor",
        "molecular_subtype", "histology_subtype", "tds_group", "sex",
    ]
    tcga = tcga[cols_needed].drop_duplicates("patient", keep="first")
    df = df.merge(tcga, on="patient", how="left")

    # Stage simplified into T-like buckets:
    def stage_bucket(s):
        if pd.isna(s):
            return "unknown"
        s = str(s)
        if "Stage IV" in s:
            return "Stage IV"
        if "Stage III" in s:
            return "Stage III"
        if "Stage II" in s:
            return "Stage II"
        if "Stage I" in s:
            return "Stage I"
        return "other"

    df["stage_bucket"] = df["ajcc_stage_group"].apply(stage_bucket)

    # Driver simplified
    def driver_bucket(s):
        if pd.isna(s):
            return "unknown"
        s = str(s)
        if s == "BRAF":
            return "BRAF"
        if s == "RAS":
            return "RAS"
        if s == "unknown":
            return "unknown"
        return "other"

    df["driver_bucket"] = df["driver_anchor"].apply(driver_bucket)

    cohort_specs = [
        ("overall", "all"),
        ("stage_bucket", None),
        ("driver_bucket", None),
        ("molecular_subtype", None),
        ("histology_subtype", None),
        ("tds_group", None),
        ("sex", None),
        ("tss", None),
    ]

    rows = []
    skipped = []
    for col, only in cohort_specs:
        if col == "overall":
            sub = df
            y = sub["label"].to_numpy(int)
            p = sub["prob_DM1"].to_numpy(float)
            auc, lo, hi = boot_auc(y, p)
            rows.append(dict(group="overall", level="all", n=int(len(sub)),
                             n_pos=int((y == 1).sum()), n_neg=int((y == 0).sum()),
                             auc=auc, ci_lo=lo, ci_hi=hi))
            continue
        if col not in df.columns:
            skipped.append(f"{col}: column missing")
            continue
        for lvl, sub in df.groupby(col, dropna=False):
            level = "NA" if pd.isna(lvl) else str(lvl)
            y = sub["label"].to_numpy(int)
            p = sub["prob_DM1"].to_numpy(float)
            n = len(sub)
            n_pos = int((y == 1).sum())
            n_neg = int((y == 0).sum())
            if n < 5 or n_pos == 0 or n_neg == 0:
                rows.append(dict(group=col, level=level, n=int(n),
                                 n_pos=n_pos, n_neg=n_neg,
                                 auc=None, ci_lo=None, ci_hi=None))
                continue
            auc, lo, hi = boot_auc(y, p)
            rows.append(dict(group=col, level=level, n=int(n),
                             n_pos=n_pos, n_neg=n_neg,
                             auc=auc, ci_lo=lo, ci_hi=hi))

    res_df = pd.DataFrame(rows)
    res_df.to_csv(OUT / "subgroup_aucs.tsv", sep="\t", index=False)
    res_df.to_json(OUT / "subgroup_aucs.json", orient="records", indent=2)

    # Forest plot — only AUC-computable rows
    plot_df = res_df[res_df["auc"].notna()].copy()
    plot_df["label_str"] = plot_df["group"] + " | " + plot_df["level"] + \
        " (n=" + plot_df["n"].astype(str) + ")"
    # sort by group, then by level
    plot_df = plot_df.sort_values(["group", "level"]).reset_index(drop=True)
    fig_h = max(3.5, 0.30 * len(plot_df) + 1.5)
    fig, ax = plt.subplots(figsize=(7.5, fig_h))
    y_pos = np.arange(len(plot_df))[::-1]
    ax.errorbar(plot_df["auc"], y_pos,
                xerr=[plot_df["auc"] - plot_df["ci_lo"],
                      plot_df["ci_hi"] - plot_df["auc"]],
                fmt="o", color="#1f77b4", ecolor="#888", capsize=3, markersize=5)
    ax.axvline(0.5, color="gray", linestyle="--", linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df["label_str"], fontsize=8)
    ax.set_xlabel("Subgroup AUC (95% CI, 1000-bootstrap)")
    ax.set_xlim(0.0, 1.05)
    ax.set_title("ViT-L CLAM — Subgroup AUC stratification (TCGA-THCA)")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGS / "figS3_subgroup_forest.png", dpi=300)
    fig.savefig(FIGS / "figS3_subgroup_forest.pdf")
    plt.close(fig)

    summary = {
        "n_total_predictions": int(len(df)),
        "n_subgroup_rows": int(len(res_df)),
        "n_with_auc": int(len(plot_df)),
        "n_bootstrap": int(N_BOOT),
        "skipped": skipped,
        "columns_used": ["stage_bucket", "driver_bucket", "molecular_subtype",
                         "histology_subtype", "tds_group", "sex", "tss"],
        "runtime_sec": float(time.time() - t0),
    }
    (OUT / "subgroup_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[subgroup] rows={len(res_df)}  with_auc={len(plot_df)}  skipped={skipped}")


if __name__ == "__main__":
    main()
