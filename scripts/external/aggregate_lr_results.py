"""
Aggregate per-sample LR results across all 16 GSE250521 samples → per-condition stats.

Output:
  project/results/pantheonos_demo/_AGGREGATE/
    per_sample_summary.csv    — one row per sample (n_groups, n_sig, top LR)
    lr_condition_means.csv    — for each LR pair, mean across samples per condition
    lr_condition_pvals.csv    — same, but p-values
    lr_condition_recurrence.csv — fraction of samples per condition where each pair is significant
    lr_top_per_condition.csv  — top recurrent LR pairs per stage
    aggregate_dotplot.png/pdf — heatmap of recurrence
    cross_stage_delta.csv     — delta between PT/PTC/LPTC/ATC for top pairs
    SUMMARY.md                — narrative
"""
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("project/results/pantheonos_demo")
OUT = ROOT / "_AGGREGATE"

STAGE_OF = {
    "N-1": "PT", "N-2": "PT", "N-3": "PT", "N-4": "PT",
    "PTC-1": "PTC", "PTC-2": "PTC", "PTC-3": "PTC", "PTC-4": "PTC",
    "LPTC-1": "LPTC", "LPTC-2": "LPTC", "LPTC-3": "LPTC", "LPTC-4": "LPTC",
    "ATC-1": "ATC", "ATC-2": "ATC", "ATC-3": "ATC", "ATC-4": "ATC",
}
STAGE_ORDER = ["PT", "PTC", "LPTC", "ATC"]


def load_sample(sample_dir: Path) -> dict | None:
    meta = sample_dir / "run_meta.json"
    means = sample_dir / "lr_means.csv"
    pvals = sample_dir / "lr_pvalues.csv"
    if not (meta.exists() and means.exists() and pvals.exists()):
        return None
    meta_d = json.loads(meta.read_text())
    if meta_d.get("lr_skipped", False) or not means.exists() or not pvals.exists():
        return {
            "sample": sample_dir.name,
            "stage": STAGE_OF.get(sample_dir.name, "?"),
            "meta": meta_d,
            "means": None,
            "pvals": None,
            "skipped": True,
        }
    try:
        return {
            "sample": sample_dir.name,
            "stage": STAGE_OF.get(sample_dir.name, "?"),
            "meta": meta_d,
            "means": pd.read_csv(means, header=[0,1], index_col=[0,1]),
            "pvals": pd.read_csv(pvals, header=[0,1], index_col=[0,1]),
            "skipped": False,
        }
    except Exception as e:
        print(f"  WARN parse {sample_dir.name}: {e}")
        return None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    samples = []
    for d in sorted(ROOT.iterdir()):
        if d.is_dir() and not d.name.startswith("_"):
            s = load_sample(d)
            if s:
                samples.append(s)
    print(f"Loaded {len(samples)} samples")
    by_stage = defaultdict(list)
    for s in samples:
        by_stage[s["stage"]].append(s)

    # 1) Per-sample summary
    per_sample = []
    for s in samples:
        m = s["meta"]
        per_sample.append({
            "sample": s["sample"],
            "stage": s["stage"],
            "n_spots": m["n_sp_spots"],
            "n_sc_cells": m["n_sc_cells"],
            "n_groups": m["n_groups"],
            "n_significant_lr": m["n_significant_lr"],
            "elapsed_sec": m["elapsed_sec"],
        })
    pd.DataFrame(per_sample).to_csv(OUT / "per_sample_summary.csv", index=False)

    # 2) LR pair recurrence per stage — fraction of samples where p<0.05 (any cluster pair)
    samples_with_lr = [s for s in samples if not s.get("skipped", False)]
    all_lr_pairs = set()
    for s in samples_with_lr:
        all_lr_pairs.update(s["pvals"].index.tolist())
    all_lr_pairs = sorted(all_lr_pairs)
    print(f"  {len(all_lr_pairs)} unique LR pairs across {len(samples_with_lr)} samples (LR-skipped: {len(samples)-len(samples_with_lr)})")

    rec_records = []
    for lr in all_lr_pairs:
        row = {"lr_pair": f"{lr[0]}__{lr[1]}" if isinstance(lr, tuple) else str(lr)}
        for stage in STAGE_ORDER:
            stage_samples_all = by_stage.get(stage, [])
            stage_samples_lr = [s for s in stage_samples_all if not s.get("skipped", False)]
            if not stage_samples_lr:
                row[f"{stage}_n"] = len(stage_samples_all)
                row[f"{stage}_recurrence"] = np.nan
                continue
            sig_count = 0
            for s in stage_samples_lr:
                if lr in s["pvals"].index:
                    vals = s["pvals"].loc[lr].values
                    vals = vals[~pd.isna(vals)].astype(float)
                    if (vals < 0.05).any():
                        sig_count += 1
            row[f"{stage}_n"] = len(stage_samples_all)
            row[f"{stage}_recurrence"] = sig_count / len(stage_samples_lr)
        rec_records.append(row)
    rec_df = pd.DataFrame(rec_records)
    rec_df.to_csv(OUT / "lr_condition_recurrence.csv", index=False)

    # 3) Top recurrent LR pairs per stage (recurrence ≥ 0.75)
    top_per_stage = []
    for stage in STAGE_ORDER:
        col = f"{stage}_recurrence"
        if col not in rec_df.columns:
            continue
        sub = rec_df[rec_df[col] >= 0.75].sort_values(col, ascending=False).head(50)
        for _, r in sub.iterrows():
            top_per_stage.append({
                "stage": stage,
                "lr_pair": r["lr_pair"],
                "recurrence": r[col],
                "n_samples": r[f"{stage}_n"],
                **{f"{s}_recurrence": r.get(f"{s}_recurrence", np.nan) for s in STAGE_ORDER if s != stage},
            })
    pd.DataFrame(top_per_stage).to_csv(OUT / "lr_top_per_condition.csv", index=False)

    # 4) Stage-specific (gain/loss) delta — pairs in stage X but not in others
    stage_specific = []
    for stage in STAGE_ORDER:
        for _, r in rec_df.iterrows():
            others = [r.get(f"{s}_recurrence", np.nan) for s in STAGE_ORDER if s != stage]
            others = [x for x in others if not np.isnan(x)]
            if len(others) == 0:
                continue
            this = r.get(f"{stage}_recurrence", np.nan)
            if np.isnan(this):
                continue
            delta = this - max(others)
            if delta >= 0.5:
                stage_specific.append({
                    "stage": stage,
                    "lr_pair": r["lr_pair"],
                    "this_recurrence": this,
                    "max_other": max(others),
                    "delta": delta,
                })
    pd.DataFrame(stage_specific).sort_values(["stage", "delta"], ascending=[True, False]).to_csv(
        OUT / "cross_stage_delta.csv", index=False
    )

    # 5) Heatmap: top 50 most variable LR pairs (max - min across stages)
    rec_df["variability"] = rec_df[[f"{s}_recurrence" for s in STAGE_ORDER]].max(axis=1) - \
                            rec_df[[f"{s}_recurrence" for s in STAGE_ORDER]].min(axis=1)
    top_var = rec_df.sort_values("variability", ascending=False).head(50)
    heat = top_var[[f"{s}_recurrence" for s in STAGE_ORDER]].copy()
    heat.index = top_var["lr_pair"].apply(lambda x: x[:50])
    heat.columns = STAGE_ORDER

    import matplotlib.pyplot as plt
    import seaborn as sns
    fig, ax = plt.subplots(figsize=(7, 14))
    sns.heatmap(heat, cmap="rocket_r", vmin=0, vmax=1, ax=ax, cbar_kws={"label": "Recurrence (fraction sig)"})
    ax.set_title("Top 50 most stage-variable LR pairs\n(GSE250521 Visium × Lu 2023 sc → MOSCOT → Squidpy)")
    fig.tight_layout()
    fig.savefig(OUT / "aggregate_heatmap.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "aggregate_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    # 6) Stage-trend bar plot (n_significant_lr)
    fig, ax = plt.subplots(figsize=(8, 5))
    ps = pd.DataFrame(per_sample)
    sns.boxplot(data=ps, x="stage", y="n_significant_lr", order=STAGE_ORDER, ax=ax,
                palette="Blues", showfliers=False)
    sns.stripplot(data=ps, x="stage", y="n_significant_lr", order=STAGE_ORDER, ax=ax,
                  color="black", size=6, jitter=0.15)
    ax.set_yscale("log")
    ax.set_ylabel("# significant LR pairs (log)")
    ax.set_title("LR diversity across thyroid dediff axis")
    fig.tight_layout()
    fig.savefig(OUT / "stage_trend_n_significant.png", dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "stage_trend_n_significant.pdf", bbox_inches="tight")
    plt.close(fig)

    print(f"[done] {OUT}")
    print(f"  per_sample_summary.csv      — {len(per_sample)} rows")
    print(f"  lr_condition_recurrence.csv — {len(rec_df)} LR pairs × {len(STAGE_ORDER)} stages")
    print(f"  lr_top_per_condition.csv    — {len(top_per_stage)} top pairs")
    print(f"  cross_stage_delta.csv       — {len(stage_specific)} stage-specific pairs")
    print(f"  aggregate_heatmap.{{png,pdf}}")
    print(f"  stage_trend_n_significant.{{png,pdf}}")


if __name__ == "__main__":
    main()
