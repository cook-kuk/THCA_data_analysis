#!/usr/bin/env python3
"""E-MTAB-12837 / E-MTAB-12900 — 8-gene panel vs redifferentiation + I-131 response.

These two ArrayExpress submissions are the only public per-patient datasets that pair a
tumour molecular profile with the response to a MAPK-inhibitor redifferentiation protocol
followed by radioiodine — i.e. an intervention whose entire purpose is to restore iodine
handling. Response is recorded per sample as RECIST and PERCIST responder / non-responder,
and the submissions also carry `tumor cell content`, the covariate that invalidated the
GSE151179 result.

  E-MTAB-12837  n=13  BRAF V600E arm     (ENA PRJEB61187)
  E-MTAB-12900  n=8   RAS-mutant arm     (ENA PRJEB61819)
  Assay: HTG EdgeSeq targeted RNA-seq, single-end.  Trial paper PMID 37074727.

Because the assay is a targeted panel, the first thing this script reports is how many of
the eight panel genes are actually measurable. If a gene is absent from the HTG probe set
it will show zero counts across every sample, and the panel cannot be scored here — that
is a real possible outcome and is reported rather than worked around.

Outputs:
  results/tables/meraiode_panel_per_sample_2026_08_06.tsv
  results/tables/meraiode_response_test_2026_08_06.tsv
  results/figures/figure_meraiode_redifferentiation_2026_08_06.{png,pdf}
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RAW = Path("/data/rai_atlas/raw/E-MTAB")
FASTQ = RAW / "fastq"
QUANT = RAW / "quant"
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
for d in (QUANT, FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

KALLISTO = Path("/tmp/claude-1000/-home-seungho-personal-THCA-data-analysis/"
                "b9611370-1293-405b-bb35-984bcf91262c/scratchpad/kallisto/kallisto")
# Index restricted to the eight panel genes (93 transcripts). The full gencode index on
# disk was built by an incompatible kallisto version, and a targeted assay only needs the
# panel loci anyway. Abundance is therefore normalised as counts per million TOTAL reads
# in the run (from the ENA read_count), not as kallisto TPM, because TPM computed inside
# an eight-gene index is a within-panel composition and not comparable across samples.
INDEX = RAW / "panel8.idx"
PANEL_FA = Path("/data/thca/reference_kallisto/gencode.v44.8gene.fa")

PANEL_8 = ["SLC5A5", "TPO", "TG", "TSHR", "PAX8", "NKX2-1", "FOXE1", "DIO1"]
STAMP = "2026_08_06"
SEED = 20260806
# HTG EdgeSeq is a probe-based assay: fragments are a fixed length, so single-end
# kallisto needs an explicit fragment-length mean/sd rather than an inferred one.
FRAG_MEAN, FRAG_SD = 120, 20


def load_sdrf():
    frames = []
    for acc in ("E-MTAB-12837", "E-MTAB-12900"):
        d = pd.read_csv(RAW / f"{acc}.sdrf.txt", sep="\t")
        cols = {c.lower().strip(): c for c in d.columns}

        def col(*keys):
            for k in keys:
                for lc, orig in cols.items():
                    if k in lc:
                        return orig
            return None

        out = pd.DataFrame({
            "study": acc,
            "run": d[col("comment[ena_run]", "ena_run")] if col("comment[ena_run]", "ena_run")
                   else np.nan,
            "individual": d[col("characteristics[individual]")],
            "recist": d[col("characteristics[recist]")].astype(str).str.strip().str.lower(),
            "percist": d[col("characteristics[percist]")].astype(str).str.strip().str.lower(),
            "clinical": d[col("characteristics[clinical history]")].astype(str).str.strip().str.lower(),
            "purity": pd.to_numeric(
                d[col("characteristics[tumor cell content]")].astype(str)
                .str.extract(r"([0-9]*\.?[0-9]+)")[0], errors="coerce"),
            "file": d[col("submitted_file_name")],
        })
        frames.append(out)
    s = pd.concat(frames, ignore_index=True)
    if s["run"].isna().all():
        ena = pd.read_csv(RAW / "ena_runs.tsv", sep="\t")
        ena = ena[ena.run_accession != "run_accession"]
        s = s.reset_index(drop=True)
        s["run"] = ena["run_accession"].reset_index(drop=True)[:len(s)]
    return s


def quantify(runs):
    if not KALLISTO.exists():
        raise SystemExit(f"kallisto binary not found at {KALLISTO}")
    if not INDEX.exists():
        print("building panel index ...")
        subprocess.run([str(KALLISTO), "index", "-i", str(INDEX), str(PANEL_FA)],
                       check=True, capture_output=True)
    for run in runs:
        out = QUANT / run
        if (out / "abundance.tsv").exists():
            continue
        fq = FASTQ / f"{run}.fastq.gz"
        if not fq.exists():
            print(f"  MISSING fastq for {run} — skipped")
            continue
        out.mkdir(parents=True, exist_ok=True)
        cmd = [str(KALLISTO), "quant", "-i", str(INDEX), "-o", str(out),
               "--single", "-l", str(FRAG_MEAN), "-s", str(FRAG_SD), "-t", "4", str(fq)]
        print("  quantifying", run)
        subprocess.run(cmd, check=True, capture_output=True)


def total_reads():
    """Total reads per run, from the ENA report — the denominator for CPM."""
    ena = pd.read_csv(RAW / "ena_runs.tsv", sep="\t")
    ena = ena[ena["run_accession"] != "run_accession"]
    if "read_count" not in ena.columns:
        return {}
    return dict(zip(ena["run_accession"], pd.to_numeric(ena["read_count"], errors="coerce")))


def gene_cpm(runs):
    """Panel-gene counts per million total sequenced reads."""
    tot = total_reads()
    cols = {}
    for run in runs:
        f = QUANT / run / "abundance.tsv"
        if not f.exists():
            continue
        a = pd.read_csv(f, sep="\t")
        a["gene"] = a["target_id"].str.split("|").str[5]
        counts = a.groupby("gene")["est_counts"].sum()
        denom = tot.get(run)
        if not denom or not np.isfinite(denom):
            denom = counts.sum()
        cols[run] = counts / denom * 1e6
    return pd.DataFrame(cols)


def main():
    s = load_sdrf()
    print(f"samples in SDRF: {len(s)}  ({s.study.value_counts().to_dict()})")
    print("response labels:", s.clinical.value_counts().to_dict())

    quantify(s["run"].dropna().tolist())
    g = gene_cpm(s["run"].dropna().tolist())
    if g.empty:
        raise SystemExit("no quantifications produced — check fastq download")
    print(f"quantified samples: {g.shape[1]}")

    present = [x for x in PANEL_8 if x in g.index]
    detected = [x for x in present if (g.loc[x] > 0).sum() >= max(3, 0.5 * g.shape[1])]
    print(f"panel genes in reference: {len(present)}/8")
    print(f"panel genes DETECTED on the HTG probe set (nonzero in >=50% samples): "
          f"{len(detected)}/8 -> {detected}")
    for x in present:
        print(f"    {x:8s} median CPM {g.loc[x].median():10.2f}  "
              f"nonzero {int((g.loc[x] > 0).sum())}/{g.shape[1]}")

    if len(detected) < 4:
        print("\nSTOP: the HTG EdgeSeq panel does not cover enough of the 8 genes to score "
              "the signature. Reporting coverage only; no association test is run.")
        pd.DataFrame({"gene": present,
                      "median_cpm": [float(g.loc[x].median()) for x in present],
                      "n_nonzero": [int((g.loc[x] > 0).sum()) for x in present],
                      "n_samples": g.shape[1],
                      "detected": [x in detected for x in present]}).to_csv(
            TAB / f"meraiode_panel_coverage_{STAMP}.tsv", sep="\t", index=False)
        return

    lg = np.log2(g.loc[detected] + 1)
    z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1, ddof=1) + 1e-9, axis=0)
    panel_z = z.mean(axis=0)

    s = s.set_index("run")
    s["panel_z"] = panel_z.reindex(s.index)
    s = s.dropna(subset=["panel_z"])
    s.to_csv(TAB / f"meraiode_panel_per_sample_{STAMP}.tsv", sep="\t")

    rows = []
    for label, col in [("RECIST", "recist"), ("PERCIST", "percist"),
                       ("clinical history", "clinical")]:
        r = s[s[col].str.contains("respon", na=False)]
        pos = r.loc[~r[col].str.contains("non", na=False), "panel_z"].values
        neg = r.loc[r[col].str.contains("non", na=False), "panel_z"].values
        if len(pos) < 2 or len(neg) < 2:
            continue
        na, nb = len(pos), len(neg)
        sp = np.sqrt(((na - 1) * pos.var(ddof=1) + (nb - 1) * neg.var(ddof=1)) / (na + nb - 2))
        dd = (pos.mean() - neg.mean()) / sp if sp > 0 else np.nan
        u = stats.mannwhitneyu(pos, neg, alternative="two-sided")
        rows.append(dict(endpoint=label, n_responder=na, n_nonresponder=nb,
                         median_responder=float(np.median(pos)),
                         median_nonresponder=float(np.median(neg)),
                         cohens_d=float(dd), auc=float(u.statistic / (na * nb)),
                         p=float(u.pvalue)))
    res = pd.DataFrame(rows)
    res.to_csv(TAB / f"meraiode_response_test_{STAMP}.tsv", sep="\t", index=False)
    print("\n", res.to_string(index=False) if len(res) else "no testable endpoint")

    if not len(res):
        return
    fig, ax = plt.subplots(figsize=(6.4, 5.2), facecolor="white", layout="constrained")
    r = s[s["clinical"].str.contains("respon", na=False)]
    pos = r.loc[~r["clinical"].str.contains("non", na=False), "panel_z"].values
    neg = r.loc[r["clinical"].str.contains("non", na=False), "panel_z"].values
    bp = ax.boxplot([neg, pos], tick_labels=["Non-responder", "Responder"],
                    showfliers=False, patch_artist=True, widths=0.55)
    for patch, c in zip(bp["boxes"], ["#9c4742", "#2f6f4f"]):
        patch.set_facecolor(c); patch.set_alpha(0.28); patch.set_edgecolor(c)
    rng = np.random.default_rng(SEED)
    for i, (vals, c) in enumerate(zip([neg, pos], ["#9c4742", "#2f6f4f"]), start=1):
        ax.scatter(rng.normal(i, 0.06, len(vals)), vals, s=42, color=c, alpha=0.85,
                   edgecolor="white", linewidth=0.6, zorder=3)
        ax.text(i, 0.02, f"n = {len(vals)}", ha="center", va="bottom", fontsize=9,
                color="#444", transform=ax.get_xaxis_transform(),
                bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85))
    row = res[res.endpoint == "clinical history"]
    ttl = (f"d = {row.cohens_d.iloc[0]:+.2f} · P = {row.p.iloc[0]:.3g}"
           if len(row) else "")
    ax.axhline(0, color="#666", lw=0.6, ls="--")
    ax.set_ylabel(f"{len(detected)}-gene panel z (of 8 measurable on HTG panel)")
    ax.set_title(f"Redifferentiation + I-131 response\nE-MTAB-12837/12900 · {ttl}",
                 fontsize=11, loc="left", fontweight="bold")
    for sp_ in ("top", "right"):
        ax.spines[sp_].set_visible(False)
    fig.savefig(FIG / f"figure_meraiode_redifferentiation_{STAMP}.png", dpi=180,
                bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"figure_meraiode_redifferentiation_{STAMP}.pdf",
                bbox_inches="tight", facecolor="white")
    print(f"wrote {FIG / f'figure_meraiode_redifferentiation_{STAMP}.png'}")


if __name__ == "__main__":
    main()
