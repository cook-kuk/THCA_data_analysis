#!/usr/bin/env python3
"""
v14 strengthening pass: external cohort scan.

Goal: Test if BRAF V600E -> high TACSTD2 (TROP2) expression replicates in an
independent thyroid cohort beyond TCGA-THCA + the GEO sets already used by
v8 / v9 / v13 / v14 (GSE27155, GSE29265, GSE33630, GSE60542).

Cohort selection rationale (see report):
- ICGC-THCA-CN, THCA-SA, THCA-US: open S3 bucket release_28 contains only
  donor / sample / specimen / ssm_open / cnsm / meth_array / mirna_seq / pexp /
  exp_seq for THCA-US. THCA-CN and THCA-SA carry only mutation calls (no
  expression). THCA-US is the same TCGA-THCA cohort already used in v14, so it
  is not an independent replication. ICGC raw RNA-seq is DACO-controlled.
- CPTAC: PDC GraphQL listing of every CPTAC{2,3,KF,Other,PTRC} study confirms
  there is NO CPTAC thyroid / THCA proteomic study (no thyroid in any
  primary_site or submitter_id_name).
- Therefore use a fresh public GEO cohort: GSE58545 (Bauer et al., 2014),
  HG-U133A microarray, 27 PTC tumors with explicit BRAF/RET/RAS status in GEO
  metadata + 18 normal thyroids. Not used by any prior v8/v9/v13/v14 analysis
  (grep over results/ confirms novelty).

This script:
  1. parses the series_matrix metadata to extract BRAF/RET/RAS status per GSM
  2. extracts the TACSTD2 probe (202286_s_at on GPL96) from the
     RMA-normalised expression matrix embedded in the series_matrix file
  3. computes a sample z-score over PTC tumours only (matching the v14
     definition: tumour-level claim, exclude normals from z-score)
  4. compares mean TACSTD2 z(PTC) BRAF V600E (+) vs BRAF (-): RET-fusion or
     RAS-mutant or, lacking other data, genotype-undescribed PTC. Welch t test
     and Mann-Whitney both reported.
  5. writes results/v14_strengthening/external_cohort_tacstd2_braf.tsv
"""
from __future__ import annotations

import gzip
import math
import os
import re
import statistics
import sys
from pathlib import Path

PROJ = Path("/home/seungho/personal/THCA_data_analysis/project")
SM = Path("/tmp/gse58545/GSE58545_series_matrix.txt.gz")
OUT_DIR = PROJ / "results" / "v14_strengthening"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_TSV = OUT_DIR / "external_cohort_tacstd2_braf.tsv"

# GPL96 probe for TACSTD2 / TROP2
TACSTD2_PROBE = "202286_s_at"


def parse_series_matrix(path: Path):
    """Return (sample_meta_dict, header_gsms, probe_row) for TACSTD2 probe."""
    sample_titles = []
    sample_gsms = []
    sample_source = []
    sample_chars: list[list[str]] = []  # list of rows, each is list per sample
    probe_row = None
    in_table = False
    table_header = None
    with gzip.open(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                if table_header is None:
                    fields = [x.strip().strip('"') for x in line.split("\t")]
                    table_header = fields
                    continue
                # only fetch the row we care about
                if line.startswith(f'"{TACSTD2_PROBE}"\t') or line.startswith(f"{TACSTD2_PROBE}\t"):
                    fields = line.split("\t")
                    probe_row = [fields[0].strip('"')] + [
                        float(x) if x not in ("", '""') else math.nan for x in fields[1:]
                    ]
                continue
            # metadata
            if line.startswith("!Sample_title"):
                sample_titles = [x.strip('"') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                sample_gsms = [x.strip('"') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                sample_source = [x.strip('"') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                sample_chars.append([x.strip('"') for x in line.split("\t")[1:]])
    return sample_titles, sample_gsms, sample_source, sample_chars, table_header, probe_row


def extract_braf_status(title: str, source: str, chars_for_sample: list[str]) -> str:
    """Return BRAF / RET / RAS / OTHER / NORMAL.

    GSE58545 records BRAF/RET/RAS status as a Sample_characteristics_ch1 field
    of the form
        "braf/ret/ras status; ... : BRAF(+)"   (or RET(+), RAS(+))
    For normal-thyroid samples the corresponding characteristics row is
    something else entirely (Sex, age, ...).  We therefore look at the
    characteristics field whose key text starts with 'braf/ret/ras status'
    and parse only its value (the part after the final ':').
    The PTC vs normal label is reliably encoded in Sample_source_name_ch1.
    """
    src = (source or "").lower()
    if "normal" in src:
        return "NORMAL"
    # Find the characteristics field that records the genotype call.
    geno_call = ""
    for f in chars_for_sample:
        fl = f.lower()
        if fl.startswith("braf/ret/ras status"):
            # value is after the LAST ':'
            geno_call = f.rsplit(":", 1)[-1].strip().upper()
            break
    # also look in the title as a fallback (titles encode 'PTC - BRAF(+) - NIS131')
    title_u = (title or "").upper()
    if not geno_call:
        if "BRAF(+)" in title_u:
            geno_call = "BRAF(+)"
        elif "RET(+)" in title_u:
            geno_call = "RET(+)"
        elif "RAS(+)" in title_u:
            geno_call = "RAS(+)"
    if geno_call.startswith("BRAF"):
        return "BRAF_V600E"
    if geno_call.startswith("RET"):
        return "RET_FUSION"
    if geno_call.startswith("RAS"):
        return "RAS_MUT"
    return "OTHER"


def welch_t(x: list[float], y: list[float]):
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return float("nan"), float("nan"), float("nan")
    mx, my = statistics.mean(x), statistics.mean(y)
    vx, vy = statistics.variance(x), statistics.variance(y)
    se = math.sqrt(vx / nx + vy / ny)
    if se == 0:
        return float("nan"), float("nan"), float("nan")
    t = (mx - my) / se
    df_num = (vx / nx + vy / ny) ** 2
    df_den = (vx / nx) ** 2 / max(nx - 1, 1) + (vy / ny) ** 2 / max(ny - 1, 1)
    df = df_num / df_den if df_den > 0 else float("nan")
    # two-sided p via stdlib regularised incomplete beta lookup absent;
    # use a survival-function approximation via math.erfc on z when df>=30,
    # else a Welch-Satterthwaite + t to z bridge.
    p = student_t_two_sided(t, df)
    return t, df, p


def student_t_two_sided(t: float, df: float) -> float:
    """Two-sided Student-t p-value using a series for the regularised
    incomplete beta function. Pure stdlib.
    P(|T| >= |t|) = I_x(df/2, 1/2) with x = df/(df + t^2)  (NR 6.4.9)."""
    if math.isnan(t) or math.isnan(df) or df <= 0:
        return float("nan")
    x = df / (df + t * t)
    p = reg_incomplete_beta(x, df / 2.0, 0.5)
    return max(0.0, min(1.0, p))


def reg_incomplete_beta(x: float, a: float, b: float) -> float:
    """Regularised incomplete beta I_x(a,b) via continued fraction.
    Numerical Recipes 6.4."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    bt = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log(1.0 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * betacf(x, a, b) / a
    return 1.0 - bt * betacf(1.0 - x, b, a) / b


def betacf(x: float, a: float, b: float, max_iter: int = 200, eps: float = 3e-7) -> float:
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        del_ = d * c
        h *= del_
        if abs(del_ - 1.0) < eps:
            return h
    return h


def mann_whitney_u(x: list[float], y: list[float]):
    """Two-sided Mann-Whitney U with normal approximation, no ties correction
    needed for this small sample (returns p, U)."""
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return float("nan"), float("nan")
    pooled = [(v, 0) for v in x] + [(v, 1) for v in y]
    pooled.sort()
    # rank with ties = average
    ranks = [0.0] * len(pooled)
    i = 0
    while i < len(pooled):
        j = i
        while j + 1 < len(pooled) and pooled[j + 1][0] == pooled[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0  # 1-based mean rank
        for k in range(i, j + 1):
            ranks[k] = avg
        i = j + 1
    Rx = sum(r for r, (_v, g) in zip(ranks, pooled) if g == 0)
    U1 = Rx - nx * (nx + 1) / 2.0
    U2 = nx * ny - U1
    U = min(U1, U2)
    mu = nx * ny / 2.0
    sigma = math.sqrt(nx * ny * (nx + ny + 1) / 12.0)
    if sigma == 0:
        return U, float("nan")
    z = (U - mu) / sigma
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return U, p


def main():
    if not SM.exists():
        sys.exit(f"missing series_matrix at {SM}")
    titles, gsms, sources, chars, header, probe_row = parse_series_matrix(SM)
    if probe_row is None:
        sys.exit(f"probe {TACSTD2_PROBE} not found")

    # align expression to GSMs
    if header[0] != "ID_REF":
        sys.exit(f"unexpected table header: {header[:3]}")
    expr_gsms = header[1:]
    expr_vals = probe_row[1:]
    if expr_gsms != gsms:
        # try to map by gsm
        idx = {g: i for i, g in enumerate(expr_gsms)}
        expr_vals = [expr_vals[idx[g]] for g in gsms]

    # build per-sample meta
    n = len(gsms)
    chars_per_sample = [[row[i] for row in chars] for i in range(n)]
    statuses = [
        extract_braf_status(titles[i], sources[i], chars_per_sample[i])
        for i in range(n)
    ]

    rows = []
    for i in range(n):
        rows.append(
            {
                "gsm": gsms[i],
                "title": titles[i],
                "source": sources[i],
                "status": statuses[i],
                "tacstd2_log2_rma": expr_vals[i],
            }
        )

    # z-score within tumours only (PTC); v14 reports tumour-level z
    tumour_status = {"BRAF_V600E", "RET_FUSION", "RAS_MUT", "OTHER"}
    tumour_vals = [r["tacstd2_log2_rma"] for r in rows if r["status"] in tumour_status]
    mu_t = statistics.mean(tumour_vals)
    sd_t = statistics.pstdev(tumour_vals) or 1.0
    for r in rows:
        if r["status"] in tumour_status:
            r["tacstd2_z_in_PTC"] = (r["tacstd2_log2_rma"] - mu_t) / sd_t
        else:
            r["tacstd2_z_in_PTC"] = float("nan")

    # also a z-score vs all samples (PTC + normal) for sanity / contrast
    all_vals = [r["tacstd2_log2_rma"] for r in rows]
    mu_a = statistics.mean(all_vals)
    sd_a = statistics.pstdev(all_vals) or 1.0
    for r in rows:
        r["tacstd2_z_in_all"] = (r["tacstd2_log2_rma"] - mu_a) / sd_a

    # comparison: BRAF V600E (+) vs all other PTC (BRAF-WT umbrella = RET+RAS+OTHER)
    braf_vals = [r["tacstd2_log2_rma"] for r in rows if r["status"] == "BRAF_V600E"]
    wt_vals = [
        r["tacstd2_log2_rma"]
        for r in rows
        if r["status"] in {"RET_FUSION", "RAS_MUT", "OTHER"}
    ]
    normal_vals = [r["tacstd2_log2_rma"] for r in rows if r["status"] == "NORMAL"]

    t, df, p_t = welch_t(braf_vals, wt_vals)
    U, p_mw = mann_whitney_u(braf_vals, wt_vals)
    t_n, df_n, p_t_n = welch_t(braf_vals, normal_vals)

    # also compute z-mean per group
    def z_mean(rows_filter, key="tacstd2_z_in_PTC"):
        vs = [r[key] for r in rows if rows_filter(r) and not math.isnan(r[key])]
        return statistics.mean(vs) if vs else float("nan")

    summary = {
        "cohort": "GSE58545 (Bauer 2014)",
        "platform": "GPL96 / HG-U133A (RMA log2)",
        "probe": TACSTD2_PROBE,
        "n_total": n,
        "n_BRAF_V600E_PTC": len(braf_vals),
        "n_BRAF_WT_PTC": len(wt_vals),
        "n_NORMAL": len(normal_vals),
        "mean_log2_TACSTD2_BRAF_V600E": statistics.mean(braf_vals) if braf_vals else float("nan"),
        "mean_log2_TACSTD2_BRAF_WT_PTC": statistics.mean(wt_vals) if wt_vals else float("nan"),
        "mean_log2_TACSTD2_NORMAL": statistics.mean(normal_vals) if normal_vals else float("nan"),
        "mean_z_PTC_BRAF_V600E": z_mean(lambda r: r["status"] == "BRAF_V600E"),
        "mean_z_PTC_BRAF_WT": z_mean(lambda r: r["status"] in {"RET_FUSION", "RAS_MUT", "OTHER"}),
        "welch_t": t,
        "welch_df": df,
        "welch_two_sided_p": p_t,
        "mannwhitney_U": U,
        "mannwhitney_two_sided_p": p_mw,
        "welch_t_BRAFvsNORMAL": t_n,
        "welch_p_BRAFvsNORMAL": p_t_n,
        "v14_predicted_direction": "BRAF V600E > BRAF-WT (mean TACSTD2 z >0)",
        "replicates_v14_direction": (
            (statistics.mean(braf_vals) > statistics.mean(wt_vals))
            if braf_vals and wt_vals
            else None
        ),
    }

    # write per-sample tsv + summary block at the bottom
    with OUT_TSV.open("w") as f:
        f.write(
            "gsm\ttitle\tstatus\ttacstd2_log2_rma\ttacstd2_z_in_PTC\ttacstd2_z_in_all\n"
        )
        for r in rows:
            f.write(
                "\t".join(
                    [
                        r["gsm"],
                        r["title"],
                        r["status"],
                        f"{r['tacstd2_log2_rma']:.4f}",
                        ("nan" if math.isnan(r["tacstd2_z_in_PTC"]) else f"{r['tacstd2_z_in_PTC']:.4f}"),
                        f"{r['tacstd2_z_in_all']:.4f}",
                    ]
                )
                + "\n"
            )
        f.write("# SUMMARY\n")
        for k, v in summary.items():
            f.write(f"#   {k}\t{v}\n")

    # also dump a small JSON for downstream tooling
    import json
    (OUT_DIR / "external_cohort_tacstd2_braf_summary.json").write_text(
        json.dumps(summary, indent=2, default=str)
    )

    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
