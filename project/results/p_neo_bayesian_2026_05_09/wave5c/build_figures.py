"""Wave 5C — build comparison figure + WAVE5C_REPORT.md.

Inputs:
  peft_results.tsv           -- this run's 4 PEFT methods
  ../wave4a/wave4a_results.tsv (optional)  -- Wave 4A LoRA baseline
  ../wave1/...  (optional)   -- Wave 1 baseline AUROC

Outputs:
  fig_peft_comparison.png/pdf
  WAVE5C_REPORT.md
"""
from __future__ import annotations
import argparse, json, os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def fmt(v):
    if v is None or (isinstance(v, float) and (np.isnan(v))):
        return "NA"
    return f"{v:.3f}"


def load_no_overlap_aurocs(results_path):
    if not Path(results_path).exists():
        return {}
    df = pd.read_csv(results_path, sep="\t")
    sub = df[df["testset"] == "ITSNdb_no_overlap"]
    return dict(zip(sub["method"], sub["AUROC"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="peft_results.tsv")
    ap.add_argument("--wave4a", default="../wave4a/wave4a_results.tsv")
    ap.add_argument("--out_dir", default=".")
    args = ap.parse_args()
    out = Path(args.out_dir)

    df = pd.read_csv(args.results, sep="\t")

    # ----- bar chart: ITSNdb_no_overlap AUROC across PEFT methods + refs -----
    no_overlap = df[df["testset"] == "ITSNdb_no_overlap"].copy()
    methods = no_overlap["method"].tolist()
    aurocs = no_overlap["AUROC"].tolist()
    lo = no_overlap["AUROC_lo95"].tolist()
    hi = no_overlap["AUROC_hi95"].tolist()
    n_params = no_overlap["n_trainable_params"].tolist()

    # add reference lines
    wave4a_aurocs = load_no_overlap_aurocs(args.wave4a)
    wave4a_lora = wave4a_aurocs.get("lora", None)
    WAVE1_BASELINE = 0.41
    MHCFLURRY = 0.668

    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(len(methods))
    bars = ax.bar(x, aurocs, color=["#888", "#1f77b4", "#2ca02c", "#d62728"][:len(methods)],
                  alpha=0.85, edgecolor="black")
    yerr_lo = np.array(aurocs) - np.array(lo)
    yerr_hi = np.array(hi) - np.array(aurocs)
    ax.errorbar(x, aurocs, yerr=[yerr_lo, yerr_hi], fmt="none", color="black", capsize=4, lw=1)
    for i, (a, n) in enumerate(zip(aurocs, n_params)):
        ax.text(i, max(0, a) + 0.02, f"{a:.3f}\n({n:,} p)", ha="center", va="bottom", fontsize=9)
    ax.axhline(WAVE1_BASELINE, color="grey", linestyle="--", lw=1.2,
               label=f"Wave 1 baseline ({WAVE1_BASELINE})")
    ax.axhline(MHCFLURRY, color="orange", linestyle="--", lw=1.2,
               label=f"MHCflurry ({MHCFLURRY})")
    if wave4a_lora is not None:
        ax.axhline(wave4a_lora, color="purple", linestyle=":", lw=1.4,
                   label=f"Wave 4A LoRA ({wave4a_lora:.3f})")
    ax.axhline(0.5, color="black", linestyle=":", lw=0.8, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=15)
    ax.set_ylabel("AUROC (ITSNdb_no_overlap, n=106)")
    ax.set_title("Wave 5C — Modern PEFT methods on ESM2-150M (q/k/v adapters)")
    ax.set_ylim(0, max(0.9, max(aurocs + [WAVE1_BASELINE, MHCFLURRY]) + 0.1))
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "fig_peft_comparison.png", dpi=150)
    fig.savefig(out / "fig_peft_comparison.pdf")
    plt.close(fig)
    print(f"saved fig_peft_comparison.png/pdf")

    # ----- WAVE5C_REPORT.md -----
    lines = []
    lines.append("# Wave 5C — Modern PEFT methods on ESM2-150M (DoRA / VeRA / AdaLoRA vs LoRA)\n")
    lines.append("## Setup\n")
    lines.append("- Base LM: `facebook/esm2_t30_150M_UR50D` (640 hidden, 30 layers).")
    lines.append("- PEFT target modules: `q/k/v` projections in attention.")
    lines.append("- Head: mean-pool peptide + HLA encodings, concat → 256 → 256 → 1.")
    lines.append("- Loss: focal-BCE (γ=2, label_smooth=0.05), AdamW, lr=3e-4, bs=32, 20 epochs.")
    lines.append("- Sampler: source-balanced WeightedRandomSampler over the 4 train sources.")
    lines.append("- Train pool: n=2396 (CEDAR + TESLA_mmc4 + NEPdb + TESLA_mmc7_validation).")
    lines.append("- Bootstrap CI: 1000 resamples.")
    lines.append("- Headline eval: ITSNdb_no_overlap (n=106, in_master=False).\n")

    lines.append("## Trainable parameters\n")
    lines.append("| Method | rank | base PEFT params | total trainable (base+head) |")
    lines.append("|---|---:|---:|---:|")
    seen = set()
    for _, r in df.drop_duplicates("method").iterrows():
        if r["method"] in seen:
            continue
        seen.add(r["method"])
        lines.append(f"| {r['method']} | {r['rank']} | {int(r['n_base_peft_params']):,} | {int(r['n_trainable_params']):,} |")
    lines.append("")

    lines.append("## Results — AUROC across testsets (with 95% bootstrap CI)\n")
    pivot = df.pivot_table(index="testset", columns="method", values="AUROC", aggfunc="first")
    pivot_lo = df.pivot_table(index="testset", columns="method", values="AUROC_lo95", aggfunc="first")
    pivot_hi = df.pivot_table(index="testset", columns="method", values="AUROC_hi95", aggfunc="first")
    lines.append("| testset | " + " | ".join(pivot.columns) + " |")
    lines.append("|---" * (len(pivot.columns) + 1) + "|")
    for ts in pivot.index:
        cells = [ts]
        for m in pivot.columns:
            a = pivot.loc[ts, m]
            l = pivot_lo.loc[ts, m] if ts in pivot_lo.index else None
            h = pivot_hi.loc[ts, m] if ts in pivot_hi.index else None
            if pd.isna(a):
                cells.append("NA")
            else:
                cells.append(f"{a:.3f} [{l:.3f},{h:.3f}]" if not pd.isna(l) else f"{a:.3f}")
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")

    # Headline
    lines.append("## Headline — ITSNdb_no_overlap (n=106)\n")
    no_overlap_sorted = no_overlap.sort_values("AUROC", ascending=False)
    lines.append("| method | AUROC | 95% CI | n_trainable | Δ Wave1(0.41) | Δ MHCflurry(0.668) |")
    lines.append("|---|---:|:---:|---:|---:|---:|")
    for _, r in no_overlap_sorted.iterrows():
        a = r["AUROC"]
        ci = f"[{r['AUROC_lo95']:.3f},{r['AUROC_hi95']:.3f}]"
        lines.append(f"| {r['method']} | {a:.3f} | {ci} | {int(r['n_trainable_params']):,} | "
                     f"{a - WAVE1_BASELINE:+.3f} | {a - MHCFLURRY:+.3f} |")
    lines.append("")

    if wave4a_lora is not None:
        lines.append(f"**Wave 4A plain LoRA reference (no_overlap)**: AUROC = {wave4a_lora:.3f}\n")

    # Param-efficiency
    no_overlap_sorted["param_eff"] = no_overlap_sorted["AUROC"] / no_overlap_sorted["n_trainable_params"]
    eff_sorted = no_overlap_sorted.sort_values("param_eff", ascending=False)
    lines.append("## Parameter efficiency (AUROC / n_trainable)\n")
    lines.append("| method | AUROC | n_trainable | AUROC / param (×1e6) |")
    lines.append("|---|---:|---:|---:|")
    for _, r in eff_sorted.iterrows():
        lines.append(f"| {r['method']} | {r['AUROC']:.3f} | {int(r['n_trainable_params']):,} | "
                     f"{1e6 * r['param_eff']:.3f} |")
    lines.append("")

    # Honesty section
    lines.append("## Honest reading\n")
    best_row = no_overlap_sorted.iloc[0]
    best_method = best_row["method"]
    best_auroc = best_row["AUROC"]
    lora_row = no_overlap[no_overlap["method"] == "lora_baseline"]
    lora_auroc = lora_row["AUROC"].iloc[0] if len(lora_row) > 0 else None
    if lora_auroc is not None:
        delta_lora = best_auroc - lora_auroc
        better = "outperforms" if delta_lora > 0.005 else ("matches" if abs(delta_lora) <= 0.005 else "lags")
        lines.append(f"- Best PEFT variant on ITSNdb_no_overlap: **{best_method}** (AUROC={best_auroc:.3f}); "
                     f"plain LoRA-baseline run alongside = {lora_auroc:.3f}, Δ = {delta_lora:+.3f} → {best_method} {better} plain LoRA.")
    if best_auroc < MHCFLURRY:
        lines.append(f"- All four PEFT methods are **below MHCflurry ({MHCFLURRY})** on no_overlap. "
                     "PEFT-on-ESM2 alone does not close the gap on this small (n=2396) train pool — "
                     "consistent with the Wave 1/4 finding that the bottleneck is data, not the adapter family.")
    elif best_auroc > MHCFLURRY:
        lines.append(f"- **Best PEFT method clears MHCflurry** ({best_auroc:.3f} > {MHCFLURRY}). "
                     "First time this generation of methods passes the supervised binding-affinity prior.")

    lines.append("")
    lines.append("## Files\n")
    lines.append("- `train_peft.py` — main script")
    lines.append("- `peft_results.tsv` — long-format AUROC table")
    lines.append("- `predictions_lora_baseline.tsv`, `predictions_dora.tsv`, `predictions_vera.tsv`, `predictions_adalora.tsv`")
    lines.append("- `fig_peft_comparison.png/pdf`")
    lines.append("- `timings.json`")
    lines.append("")

    (out / "WAVE5C_REPORT.md").write_text("\n".join(lines))
    print("saved WAVE5C_REPORT.md")


if __name__ == "__main__":
    main()
