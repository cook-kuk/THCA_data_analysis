"""K2 (PRJEB11591) HLA distribution boundary check vs AFND South Korea baseline.

Two questions:
1) Is the K2 cohort representative of the Korean reference panel (selection-bias defense)?
2) Per-allele binomial test for the AFND-anchored alleles (DPB1*05:01, DRB1*04:05, DRB1*15:01, B*46:01).

Inputs:
- /data/thca/_repo_offload/arcasHLA/*.genotype.json (n≈260 K2 samples)
- project/manuscript_p2_brief/lit_enrich_2026_05_02/data/afnd_alleles.json

Outputs (project/results/dm1_robustness_v2026_05_08/hla_boundary/):
- k2_allele_freq.tsv          per locus, per 4-digit allele
- k2_vs_afnd_korea.tsv        binomial test for the 4 AFND-anchored alleles
- hla_boundary_summary.json
- hla_boundary_panel.png      bar+CI overlay for the 4 alleles
"""
from __future__ import annotations

import glob
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ARCAS_DIR = Path("/data/thca/_repo_offload/arcasHLA")
AFND_PATH = Path(
    "/home/seungho/personal/THCA_data_analysis/project/manuscript_p2_brief/lit_enrich_2026_05_02/data/afnd_alleles.json"
)
OUT_DIR = Path(
    "/home/seungho/personal/THCA_data_analysis/project/results/dm1_robustness_v2026_05_08/hla_boundary"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

LOCI = ["A", "B", "C", "DPB1", "DQB1", "DRB1"]


def to_4digit(allele: str) -> str:
    """A*11:303 -> A*11:303; B*39:24:01 -> B*39:24."""
    parts = allele.split(":")
    return ":".join(parts[:2])


def load_k2() -> dict[str, list[str]]:
    """Return {locus: [allele_4digit, ...]} pooled across samples (2 alleles per sample per locus)."""
    out: dict[str, list[str]] = {locus: [] for locus in LOCI}
    files = sorted(glob.glob(str(ARCAS_DIR / "*.genotype.json")))
    print(f"[K2] {len(files)} arcasHLA genotype files")
    for f in files:
        try:
            with open(f) as fh:
                g = json.load(fh)
        except json.JSONDecodeError:
            continue
        for locus in LOCI:
            for a in g.get(locus, []) or []:
                if a:
                    out[locus].append(to_4digit(a))
    return out


def kor_baseline_freq(afnd: dict, allele: str) -> tuple[float, int] | None:
    """Sample-size-weighted Korean allele freq."""
    pops = afnd.get(allele, {}).get("Korea", []) or []
    if not pops:
        return None
    total_n = 0
    total_freq_n = 0.0
    for p in pops:
        try:
            n = int(str(p["sample_size"]).replace(",", ""))
            f = float(p["allele_freq"])
        except (KeyError, ValueError):
            continue
        total_n += n
        total_freq_n += f * n
    if total_n == 0:
        return None
    return total_freq_n / total_n, total_n


def binom_two_sided(x: int, n: int, p0: float) -> float:
    """Two-sided binomial test (scipy)."""
    return stats.binomtest(x, n=n, p=p0, alternative="two-sided").pvalue


def wilson_ci(x: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    z = stats.norm.ppf(1 - alpha / 2)
    p = x / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = (z / denom) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, centre - half), min(1.0, centre + half))


def main() -> None:
    print("[1/4] loading K2 genotypes")
    k2 = load_k2()
    n_samples = max(len(v) // 2 for v in k2.values())
    print(f"      n_samples ≈ {n_samples}")

    print("[2/4] computing per-allele K2 frequencies")
    rows = []
    for locus, alleles in k2.items():
        n_chr = len(alleles)
        c = Counter(alleles)
        for a, k in c.most_common():
            f = k / n_chr if n_chr else 0.0
            lo, hi = wilson_ci(k, n_chr)
            rows.append(
                {
                    "locus": locus,
                    "allele": a,
                    "k": k,
                    "n_chr": n_chr,
                    "freq": f,
                    "ci_lo": lo,
                    "ci_hi": hi,
                }
            )
    out_path = OUT_DIR / "k2_allele_freq.tsv"
    with open(out_path, "w") as fh:
        fh.write("locus\tallele\tk\tn_chr\tfreq\tci_lo\tci_hi\n")
        for r in rows:
            fh.write(
                f"{r['locus']}\t{r['allele']}\t{r['k']}\t{r['n_chr']}\t{r['freq']:.6f}\t{r['ci_lo']:.6f}\t{r['ci_hi']:.6f}\n"
            )
    print(f"      wrote {out_path}")

    print("[3/4] AFND-anchored allele tests")
    afnd = json.load(open(AFND_PATH))
    anchored = ["DPB1*05:01", "DRB1*04:05", "DRB1*15:01", "B*46:01"]
    summary_rows = []
    panel_data = []
    for allele in anchored:
        locus = allele.split("*")[0]
        n_chr = len(k2[locus])
        k = sum(1 for a in k2[locus] if a == allele)
        k2_freq = k / n_chr if n_chr else 0.0
        k2_lo, k2_hi = wilson_ci(k, n_chr)
        kor = kor_baseline_freq(afnd, allele)
        if kor is None:
            print(f"      {allele}: no Korean baseline, skip")
            continue
        kor_p, kor_n = kor
        p_two = binom_two_sided(k, n_chr, kor_p) if n_chr else float("nan")
        ratio = k2_freq / kor_p if kor_p > 0 else float("nan")
        summary_rows.append(
            {
                "allele": allele,
                "k2_n_chr": n_chr,
                "k2_carriers_chr": k,
                "k2_freq": k2_freq,
                "k2_ci_lo": k2_lo,
                "k2_ci_hi": k2_hi,
                "afnd_korea_freq": kor_p,
                "afnd_korea_n": kor_n,
                "freq_ratio_k2_over_afnd": ratio,
                "binom_two_sided_p": p_two,
            }
        )
        panel_data.append(
            (allele, k2_freq, k2_lo, k2_hi, kor_p, p_two)
        )
    test_path = OUT_DIR / "k2_vs_afnd_korea.tsv"
    with open(test_path, "w") as fh:
        fh.write(
            "allele\tk2_n_chr\tk2_carriers_chr\tk2_freq\tk2_ci_lo\tk2_ci_hi\tafnd_korea_freq\tafnd_korea_n\tfreq_ratio_k2_over_afnd\tbinom_two_sided_p\n"
        )
        for r in summary_rows:
            fh.write(
                f"{r['allele']}\t{r['k2_n_chr']}\t{r['k2_carriers_chr']}\t{r['k2_freq']:.6f}\t{r['k2_ci_lo']:.6f}\t{r['k2_ci_hi']:.6f}\t{r['afnd_korea_freq']:.6f}\t{r['afnd_korea_n']}\t{r['freq_ratio_k2_over_afnd']:.4f}\t{r['binom_two_sided_p']:.4e}\n"
            )
    print(f"      wrote {test_path}")

    print("[4/4] panel figure")
    fig, ax = plt.subplots(figsize=(7, 4))
    xs = np.arange(len(panel_data))
    k2_freqs = [d[1] for d in panel_data]
    err_lo = [d[1] - d[2] for d in panel_data]
    err_hi = [d[3] - d[1] for d in panel_data]
    afnd_freqs = [d[4] for d in panel_data]
    pvals = [d[5] for d in panel_data]
    ax.bar(xs - 0.18, k2_freqs, width=0.32, color="#2c6fbb", label="K2 (PRJEB11591) RNA-imputed")
    ax.errorbar(xs - 0.18, k2_freqs, yerr=[err_lo, err_hi], fmt="none", ecolor="black", capsize=3)
    ax.bar(xs + 0.18, afnd_freqs, width=0.32, color="#bbbbbb", label="AFND South Korea baseline")
    for i, (k2f, p) in enumerate(zip(k2_freqs, pvals)):
        sig = "n.s." if p > 0.05 else ("*" if p > 0.01 else ("**" if p > 0.001 else "***"))
        ax.text(i, max(k2f, afnd_freqs[i]) + 0.02, sig, ha="center", fontsize=10)
    ax.set_xticks(xs)
    ax.set_xticklabels([d[0] for d in panel_data], rotation=20, ha="right")
    ax.set_ylabel("Allele frequency")
    ax.set_title(f"K2 cohort vs AFND Korean baseline (n_samples≈{n_samples})")
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "hla_boundary_panel.png", dpi=180)
    fig.savefig(OUT_DIR / "hla_boundary_panel.pdf")
    plt.close(fig)

    summary = {
        "n_samples_arcasHLA_K2": int(n_samples),
        "loci_tested": LOCI,
        "n_unique_alleles_per_locus": {
            locus: len(set(k2[locus])) for locus in LOCI
        },
        "anchored_alleles": [
            {**r, "binom_two_sided_p": float(r["binom_two_sided_p"])}
            for r in summary_rows
        ],
        "interpretation": {
            "all_anchored_within_2x_of_baseline": all(
                0.5 <= r["freq_ratio_k2_over_afnd"] <= 2.0 for r in summary_rows
            ),
            "any_significant_deviation_p_lt_0p05": any(
                r["binom_two_sided_p"] < 0.05 for r in summary_rows
            ),
        },
        "generated_at": "2026-05-08 hla_boundary",
    }
    with open(OUT_DIR / "hla_boundary_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
