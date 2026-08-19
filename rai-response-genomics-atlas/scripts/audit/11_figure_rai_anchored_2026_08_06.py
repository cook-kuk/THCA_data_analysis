#!/usr/bin/env python3
"""AUDIT 11 — figure for the RAI-anchored reanalysis.

(a) Kaplan-Meier by median split on the RAI-anchored clock, annotated with BOTH the
    naive log-rank p and the cut-point-corrected p, so the retraction is visible.
(b) Permutation null of the maximally selected log-rank statistic, with the observed
    value marked. This is the panel that shows why P = 0.043 does not survive.
(c) Forest: cohort rules, adjustments, sensitivity, and the treatment x score interaction.

Output
  figures/audit/figure_rai_anchored_reanalysis_2026_08_06.{png,pdf}
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "audit"
SEED = 20260806
NPERM = 5000

a9spec = importlib.util.spec_from_file_location(
    "a9", ROOT / "scripts" / "audit" / "09_rai_anchored_robustness_2026_08_06.py")
_a9 = importlib.util.module_from_spec(a9spec)
a9spec.loader.exec_module(_a9)

a10spec = importlib.util.spec_from_file_location(
    "a10", ROOT / "scripts" / "audit" / "10_explicit_i131_label_2026_08_06.py")
_a10 = importlib.util.module_from_spec(a10spec)
a10spec.loader.exec_module(_a10)

RED, GREEN, INK, GREY = "#9c2b25", "#2f6f4f", "#15151a", "#8a8a92"
AMBER, BLUE = "#b5761f", "#37618e"


def km_curve(t, e):
    order = np.argsort(t)
    t, e = np.asarray(t)[order], np.asarray(e)[order]
    times, surv, s, n = [0.0], [1.0], 1.0, len(t)
    for i, ti in enumerate(t):
        if e[i] == 1:
            s *= (1 - 1 / (n - i))
        times.append(ti)
        surv.append(s)
    return np.array(times), np.array(surv)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    d = _a10.build()
    a = _a10.anchored(d[d["radfile_mci"]])
    t, e, x = a["time"].values, a["event"].values, a["RAI_8"].values
    n, nev = len(a), int(e.sum())

    # -------- recompute the cut-point sweep and permutation null
    lo_q, hi_q = np.quantile(x, [0.10, 0.90])
    cand = np.unique(x[(x >= lo_q) & (x <= hi_q)])
    _o = np.argsort(t, kind="mergesort")
    ev_pos = np.where(e[_o] == 1)[0]
    n_risk = (n - np.arange(n))[ev_pos].astype(float)

    def max_chi(xv):
        G = (xv[_o][:, None] > cand[None, :])
        n1 = np.cumsum(G[::-1], axis=0)[::-1]
        keep = (G.sum(axis=0) >= 15) & ((~G).sum(axis=0) >= 15)
        n1e = n1[ev_pos][:, keep].astype(float)
        p = n1e / n_risk[:, None]
        O1 = G[ev_pos][:, keep].sum(axis=0).astype(float)
        V = (p * (1 - p)).sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            chi = np.where(V > 0, (O1 - p.sum(axis=0)) ** 2 / V, np.nan)
        j = int(np.nanargmax(chi))
        return float(chi[j]), float(cand[keep][j])

    obs_chi, obs_cut = max_chi(x)
    perm = np.array([max_chi(rng.permutation(x))[0] for _ in range(NPERM)])
    p_corr = float((perm >= obs_chi).mean())

    med = float(np.median(x))
    lo, hi = x <= med, x > med
    chi_med, p_med = _a9.logrank(t[lo], e[lo], t[hi], e[hi])

    # Explicit axes rectangles rather than a gridspec: panel (c) needs a reserved
    # margin on BOTH sides (row labels left, numeric column right) and a shared
    # wspace cannot give that.
    fig = plt.figure(figsize=(18.4, 6.0))
    AX_A = fig.add_axes([0.045, 0.135, 0.245, 0.60])
    AX_B = fig.add_axes([0.365, 0.135, 0.215, 0.60])
    AX_C = fig.add_axes([0.728, 0.135, 0.122, 0.60])

    # ---------------------------------------------------------------- (a) KM
    ax = AX_A
    for m, col, lab in [(lo, RED, "score ≤ median"), (hi, GREEN, "score > median")]:
        tt, ss = km_curve(t[m] / 30.44, e[m])
        ax.step(tt, 100 * ss, where="post", color=col, lw=2.1,
                label=f"{lab}  (n={int(m.sum())}, events={int(e[m].sum())})")
    ax.set_xlabel("months since the first radioiodine course")
    ax.set_ylabel("recurrence-free (%)")
    ax.set_ylim(60, 101)
    ax.set_xlim(0, min(t.max() / 30.44, 120))
    ax.legend(loc="lower left", fontsize=8.2, frameon=False,
              bbox_to_anchor=(-0.01, -0.015))
    ax.set_title("a · Kaplan–Meier on the RAI-anchored clock\n"
                 f"time origin = first radioiodine course (n={n}, {nev} events)",
                 fontsize=10.2, loc="left", pad=10)
    ax.text(0.97, 0.97,
            f"naive log-rank P = {p_med:.3f}\n"
            f"cut-point-corrected P = {p_corr:.3f}\n"
            "→ does not survive correction",
            transform=ax.transAxes, ha="right", va="top", fontsize=9.0,
            linespacing=1.5, color=INK,
            bbox=dict(boxstyle="round,pad=0.45", fc="#fdf3f2", ec=RED, lw=1.0))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # ------------------------------------------------- (b) permutation null
    ax = AX_B
    ax.hist(perm, bins=48, color="#d8d3c8", edgecolor="white", lw=0.4)
    ax.axvline(obs_chi, color=RED, lw=2.2)
    ax.axvline(chi_med, color=AMBER, lw=1.6, ls="--")
    ax.set_xlabel("maximally selected log-rank χ²  (5,000 permutations of the score)")
    ax.set_ylabel("permutations")
    ax.set_xlim(0, min(perm.max() * 1.03, 22))
    ax.set_title("b · why P = 0.043 dies\n"
                 "the median is 1 of ~180 candidate cut-points",
                 fontsize=10.2, loc="left", pad=10)
    ymax = ax.get_ylim()[1]
    ax.annotate(f"observed max χ² = {obs_chi:.2f}  (cut = {obs_cut:.3f})",
                xy=(obs_chi, ymax * 0.58),
                xytext=(0.98, 0.93), textcoords="axes fraction",
                fontsize=8.8, color=RED, ha="right", va="top",
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.1,
                                connectionstyle="arc3,rad=0.18"))
    ax.annotate(f"median-split χ² = {chi_med:.2f}",
                xy=(chi_med, ymax * 0.30),
                xytext=(0.98, 0.76), textcoords="axes fraction",
                fontsize=8.6, color=AMBER, ha="right", va="top",
                arrowprops=dict(arrowstyle="->", color=AMBER, lw=1.0,
                                connectionstyle="arc3,rad=0.22"))
    ax.text(0.98, 0.56,
            f"{100*p_corr:.1f}% of permutations reach\n"
            f"the observed maximum  →  P = {p_corr:.3f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.8,
            linespacing=1.45,
            bbox=dict(boxstyle="round,pad=0.42", fc="#f7f4ee", ec=GREY, lw=0.8))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # ------------------------------------------------------------ (c) forest
    ax = AX_C
    rob = pd.read_csv(ROOT / "results" / "tables"
                      / "audit09_rai_anchored_robustness_2026_08_06.tsv", sep="\t")
    ex = pd.read_csv(ROOT / "results" / "tables"
                     / "audit10_explicit_i131_results_2026_08_06.tsv", sep="\t")

    def pick(df, key, col="analysis"):
        m = df[df[col].str.contains(key, regex=False, na=False)]
        return m.iloc[0] if len(m) else None

    items = []
    for lab, row, kind in [
        ("permissive · any mCi course", pick(ex, "permissive"), "cohort"),
        ("strict · not explicitly denied", pick(ex, "strict"), "cohort"),
        ("explicit-only · I-131 confirmed", pick(ex, "explicit-only"), "cohort"),
        ("+ stage", pick(rob, "Cox score + stage"), "adj"),
        ("+ age", pick(rob, "Cox score + age"), "adj"),
        ("+ cumulative dose", pick(rob, "Cox score + cumulative dose"), "adj"),
        ("exclude external-beam", pick(rob, "exclude external-beam"), "adj"),
        ("within explicit I-131", pick(ex, "within explicit I-131"), "strat"),
        ("within explicit NO I-131", pick(ex, "within explicit NO I-131"), "strat"),
        ("score × treatment INTERACTION", pick(ex, "INTERACTION"), "inter"),
    ]:
        if row is None:
            continue
        if "hr" in row.index and pd.notna(row.get("hr")):
            hr, cl, ch, p = row["hr"], row["ci_low"], row["ci_high"], row["p"]
            nn, ee = row["n"], row["events"]
        else:
            st = str(row.get("statistic", ""))
            try:
                hr = float(st.split("HR=")[1].split(" ")[0])
                cl = float(st.split("[")[1].split(",")[0])
                ch = float(st.split(",")[-1].rstrip("]"))
            except Exception:
                continue
            p, nn, ee = row["p_value"], row["n"], row["events"]
        items.append((lab, hr, cl, ch, p, int(nn), int(ee), kind))

    cols = {"cohort": INK, "adj": GREY, "strat": BLUE, "inter": RED}
    ys = np.arange(len(items))[::-1]
    LOX, HIX = 0.06, 8.0
    for y, (lab, hr, cl, ch, p, nn, ee, kind) in zip(ys, items):
        c = cols[kind]
        clip_l, clip_h = max(cl, LOX * 1.02), min(ch, HIX * 0.98)
        ax.plot([clip_l, clip_h], [y, y], color=c, lw=1.7, solid_capstyle="round")
        if ch > HIX:                       # mark a truncated upper bound
            ax.plot([clip_h], [y], marker=">", color=c, ms=5, mec="none")
        ax.plot([hr], [y], "o", color=c, ms=7.5 if kind == "inter" else 6,
                mec="white", mew=1.0, zorder=3)
        # row label OUTSIDE the axes on the left, numbers OUTSIDE on the right
        ax.text(-0.06, y, lab, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=8.6,
                fontweight="bold" if kind == "inter" else "normal", color=c)
        ax.text(1.09, y,
                f"{hr:4.2f} [{cl:4.2f},{ch:5.2f}]  P={p:.3f}  {ee:>2d}/{nn:<3d}",
                transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=7.7, color=GREY if kind == "adj" else (
                    RED if kind == "inter" else INK),
                family="DejaVu Sans Mono")
    ax.axvline(1, color=INK, lw=0.9, ls=":")
    ax.set_xscale("log")
    ax.set_xlim(LOX, HIX)
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4])
    ax.set_xticklabels(["0.1", "0.25", "0.5", "1", "2", "4"])
    ax.set_ylim(-1.35, len(items) - 0.35)
    ax.set_yticks([])
    ax.set_xlabel("hazard ratio (log scale)  ·  95% CI")
    ax.set_title("c · every cohort rule, adjustment and stratum gives the same answer\n"
                 "nothing crosses significance; the interaction is estimable but null",
                 fontsize=10.2, loc="left", pad=10, x=-0.62)
    ax.legend(handles=[Line2D([], [], color=cols[k], marker="o", ls="-", ms=5.5, label=v)
                       for k, v in [("cohort", "cohort definition"),
                                    ("adj", "adjustment / sensitivity"),
                                    ("strat", "treatment stratum"),
                                    ("inter", "interaction term")]],
              loc="upper center", bbox_to_anchor=(0.30, -0.145), fontsize=7.9,
              frameon=False, ncol=4, columnspacing=1.5, handletextpad=0.5)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)

    fig.suptitle("Radioiodine-anchored reanalysis — the eight-gene score and recurrence "
                 "after radioiodine (TCGA-THCA)",
                 fontsize=13.4, fontweight="bold", x=0.045, ha="left", y=0.965)
    fig.text(0.045, 0.905,
             "Time origin moved from diagnosis to the first radioiodine course; cohort "
             "restricted to patients who actually received radioiodine. "
             "Prompted by an external reviewer objection.",
             fontsize=9.3, color="#43434b", ha="left")

    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"figure_rai_anchored_reanalysis_2026_08_06.{ext}",
                    dpi=200, facecolor="white")
    print("wrote", OUT / "figure_rai_anchored_reanalysis_2026_08_06.png")
    print(f"panel b: observed max chi2={obs_chi:.3f}, corrected P={p_corr:.4f}")
    print(f"forest rows: {len(items)}")


if __name__ == "__main__":
    main()
