"""
PoC figures: Moffitt KM, neoQ histogram, immune-module heatmap, vaccine ranking.
Writes PNG and SVG to /data/pdac_poc/figures/.
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

RAW = Path("/data/pdac_poc/raw")
RES = Path("/data/pdac_poc/results")
PROC = Path("/data/pdac_poc/processed")
FIG = Path("/data/pdac_poc/figures")
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
})

EMERALD = "#35d39d"; ROSE = "#ef5f79"; VIOLET = "#9b7cff"
TEAL = "#32b8c6"; AMBER = "#f2b84b"; INK = "#0a101b"


def fig_km():
    moff = pd.read_csv(RES / "moffitt_calls.tsv", sep="\t", index_col=0)
    clin = pd.read_csv(RAW / "clinical.tsv", sep="\t").set_index("sampleId")
    os_m = pd.to_numeric(clin.get("OS_MONTHS"), errors="coerce")
    os_s = clin.get("OS_STATUS", pd.Series(dtype=str))
    ev = os_s.fillna("").str.startswith("1").astype(int)
    df = moff.join(pd.DataFrame({"os_months": os_m, "event": ev}))
    df = df.dropna(subset=["os_months"])
    df = df[df["os_months"] > 0]

    fig, ax = plt.subplots(figsize=(5.5, 4.0), dpi=140)
    for label, col in [("basal-like", ROSE), ("classical", EMERALD)]:
        sub = df[df["moffitt_call"] == label]
        kmf = KaplanMeierFitter()
        kmf.fit(sub["os_months"], sub["event"], label=f"{label} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=col, linewidth=2)
    ax.set_xlabel("Months from diagnosis"); ax.set_ylabel("Overall survival")
    ax.set_title("TCGA-PAAD · Moffitt Basal vs Classical (PoC)")
    ax.set_xlim(0, 80)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "fig_km_moffitt.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "fig_km_moffitt.svg", bbox_inches="tight")
    plt.close(fig)


def fig_neoq_hist():
    df = pd.read_csv(RES / "neoantigen_quality_table.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.5), dpi=140)
    axes[0].hist(df["NeoQ"], bins=30, color=VIOLET, alpha=0.85)
    axes[0].set_xlabel("NeoQ = R × D"); axes[0].set_ylabel("Peptides")
    axes[0].set_title("NeoQ distribution (Balachandran-PoC)")
    by_gene = df.groupby("gene")["NeoQ"].agg(["max", "count", "median"])
    by_gene = by_gene.sort_values("max", ascending=True)
    axes[1].barh(by_gene.index, by_gene["max"], color=TEAL, alpha=0.9)
    axes[1].set_xlabel("Top NeoQ"); axes[1].set_title("Top NeoQ by gene")
    for i, (g, row) in enumerate(by_gene.iterrows()):
        axes[1].text(row["max"] + 0.002, i, f"n={int(row['count'])}",
                     va="center", fontsize=8, color="#555")
    fig.tight_layout()
    fig.savefig(FIG / "fig_neoq.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "fig_neoq.svg", bbox_inches="tight")
    plt.close(fig)


def fig_immune_heatmap():
    s = json.load(open(PROC / "pdac_immune_readiness_summary.json"))
    cohen = s["cohens_d_basal_vs_classical"]
    cohen_kras = s["cohens_d_kras_mut_vs_wt"]
    keys = sorted(cohen.keys())
    arr = np.array([[cohen[k], cohen_kras.get(k, np.nan)] for k in keys])
    fig, ax = plt.subplots(figsize=(5.5, 4.0), dpi=140)
    im = ax.imshow(arr, cmap="RdBu_r", vmin=-1.0, vmax=1.0, aspect="auto")
    ax.set_yticks(range(len(keys))); ax.set_yticklabels(keys, fontsize=9)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["basal − classical", "KRASmut − wt"],
                                              rotation=20, ha="right")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            v = arr[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                        fontsize=8, color="#0a0a0a")
    ax.set_title("PDAC TME modules (Cohen's d)")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.04, label="Cohen's d")
    fig.tight_layout()
    fig.savefig(FIG / "fig_immune_heatmap.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "fig_immune_heatmap.svg", bbox_inches="tight")
    plt.close(fig)


def fig_vaccine_priority():
    df = pd.read_csv(RES / "vaccine_targets.tsv", sep="\t")
    df = df.sort_values("priority_score", ascending=True)
    colors = {"public_neoantigen": EMERALD, "tumor_associated_self": AMBER,
              "private_neoantigen": VIOLET}
    cols = [colors.get(c, "#888") for c in df["class"]]
    fig, ax = plt.subplots(figsize=(7.0, 4.6), dpi=140)
    ax.barh(df["target"], df["priority_score"], color=cols, alpha=0.92)
    for i, (t, p, prev) in enumerate(zip(df["target"], df["priority_score"],
                                         df["tcga_prevalence_pct"])):
        ax.text(p + 0.005, i, f"prev={prev}%" if pd.notna(prev) else "—",
                va="center", fontsize=8, color="#555")
    ax.set_xlabel("Priority score (0.5·prev + 0.3·KORHLA + 0.2·NeoQ)")
    ax.set_title("PDAC vaccine target ranking (PoC)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.92)
               for c in [EMERALD, AMBER, VIOLET]]
    ax.legend(handles, ["public neoantigen", "tumor-associated self",
                         "private neoantigen"], loc="lower right",
              frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_vaccine_priority.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "fig_vaccine_priority.svg", bbox_inches="tight")
    plt.close(fig)


def fig_hla_coverage():
    s = json.load(open(PROC / "vaccine_priority_summary.json"))
    pops = ["Korean", "Pan-Asian", "European"]
    covs = [s["korean_off_the_shelf_population_coverage_pct"],
            None,  # placeholder; future
            s["european_off_the_shelf_population_coverage_pct"]]
    fig, ax = plt.subplots(figsize=(4.8, 3.2), dpi=140)
    bars = ax.bar([p for p, c in zip(pops, covs) if c is not None],
                  [c for c in covs if c is not None],
                  color=[ROSE, AMBER, EMERALD][:sum(c is not None for c in covs)],
                  alpha=0.92)
    for b in bars:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1,
                f"{b.get_height():.1f}%", ha="center", fontsize=9)
    ax.set_ylim(0, max(c for c in covs if c is not None) + 8)
    ax.set_ylabel("Off-the-shelf cassette pop. coverage (%)")
    ax.set_title("Off-the-shelf KRAS+TP53+TAA cassette · pop. coverage")
    fig.tight_layout()
    fig.savefig(FIG / "fig_hla_coverage.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "fig_hla_coverage.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_km()
    fig_neoq_hist()
    fig_immune_heatmap()
    fig_vaccine_priority()
    fig_hla_coverage()
    print("[08] figures →", FIG)
    for f in sorted(FIG.glob("*.png")):
        print(" ", f.name)
