#!/usr/bin/env python3
"""Benchmark runnable external TCR experts across fast challenge panels."""

from __future__ import annotations

import subprocess
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

from common import OUT, SEED


TCR_OUT = OUT / "tcr_extension"
DISC_OUT = TCR_OUT / "model_discovery"
OUTDIR = DISC_OUT / "external_tcr_expert_benchmark"
PMTNET = Path("/data/thca/tcr_model_repos/pMTnet")
TEPCAM = Path("/data/thca/tcr_model_repos/TEPCAM")
PANEL_POS_N = int(os.environ.get("CROSS_NEO_TCR_PANEL_POS_N", "100"))


AA_RE = r"[ACDEFGHIKLMNPQRSTVWY]+"
PEP_RE = r"[ACDEFGHIKLMNPQRSTVWY]{8,15}"


def run(cmd: list[str], cwd: Path, stdout_path: Path, stderr_path: Path) -> int:
    proc = subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True, check=False)
    stdout_path.write_text(proc.stdout)
    stderr_path.write_text(proc.stderr)
    return int(proc.returncode)


def clean_hla_for_pmtnet(hla: object) -> str:
    return str(hla or "").replace("HLA-", "").strip()


def metric_row(y: pd.Series, score: pd.Series) -> dict[str, float]:
    yy = y.astype(int).to_numpy()
    ss = pd.to_numeric(score, errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(ss)
    yy, ss = yy[mask], ss[mask]
    out = {
        "n": int(len(yy)),
        "n_pos": int(yy.sum()) if len(yy) else 0,
        "n_neg": int(len(yy) - yy.sum()) if len(yy) else 0,
        "positive_mean": float(np.mean(ss[yy == 1])) if np.any(yy == 1) else np.nan,
        "negative_mean": float(np.mean(ss[yy == 0])) if np.any(yy == 0) else np.nan,
    }
    if len(yy) and len(np.unique(yy)) == 2:
        out["auprc"] = float(average_precision_score(yy, ss))
        out["auroc"] = float(roc_auc_score(yy, ss))
    else:
        out["auprc"] = np.nan
        out["auroc"] = np.nan
    return out


def choose_panel(pos: pd.DataFrame, name: str, query: str | None, n: int, seed: int) -> tuple[pd.DataFrame, dict[str, int]]:
    sub = pos.query(query).copy() if query else pos.copy()
    before = len(sub)
    # Avoid letting a single dominant peptide-HLA clone family define the panel.
    sub["_pmhc"] = sub["peptide"].astype(str) + "|" + sub["hla_4digit"].astype(str)
    sub = sub.sort_values(["source_dataset", "_pmhc", "cdr3_beta"]).drop_duplicates(["cdr3_beta", "peptide", "hla_4digit"])
    if len(sub) > n:
        per_group = max(3, int(np.ceil(n / max(1, sub["_pmhc"].nunique()))))
        capped = pd.concat(
            [g.sample(n=min(per_group, len(g)), replace=False, random_state=seed) for _, g in sub.groupby("_pmhc")],
            ignore_index=True,
        )
        sub = capped.sample(n=min(n, len(capped)), random_state=seed)
    sub = sub.drop(columns=["_pmhc"], errors="ignore").reset_index(drop=True)
    return sub, {"panel": name, "eligible_positive_rows": int(before), "sampled_positive_rows": int(len(sub))}


def make_decoys(pos: pd.DataFrame, global_pos_triplets: set[tuple[str, str, str]], global_pos_pairs: set[tuple[str, str]], seed: int) -> tuple[pd.DataFrame, dict[str, int]]:
    rng = np.random.default_rng(seed)
    dec = pos.copy()
    original = dec["cdr3_beta"].astype(str).to_numpy()
    shuffled = original.copy()
    for _ in range(100):
        rng.shuffle(shuffled)
        same = shuffled == original
        triplet_conflict = [
            (cdr3, pep, hla) in global_pos_triplets
            for cdr3, pep, hla in zip(shuffled, dec["peptide"].astype(str), dec["hla_4digit"].astype(str))
        ]
        pair_conflict = [
            (cdr3, pep) in global_pos_pairs
            for cdr3, pep in zip(shuffled, dec["peptide"].astype(str))
        ]
        bad = same | np.array(triplet_conflict, dtype=bool) | np.array(pair_conflict, dtype=bool)
        if not bad.any():
            break
    dec["cdr3_beta"] = shuffled
    dec["decoy_triplet_conflict"] = [
        (cdr3, pep, hla) in global_pos_triplets
        for cdr3, pep, hla in zip(dec["cdr3_beta"].astype(str), dec["peptide"].astype(str), dec["hla_4digit"].astype(str))
    ]
    dec["decoy_pair_conflict"] = [
        (cdr3, pep) in global_pos_pairs for cdr3, pep in zip(dec["cdr3_beta"].astype(str), dec["peptide"].astype(str))
    ]
    dec["decoy_self_same_tcr"] = dec["cdr3_beta"].astype(str).to_numpy() == original
    keep = ~(dec["decoy_triplet_conflict"] | dec["decoy_pair_conflict"] | dec["decoy_self_same_tcr"])
    report = {
        "candidate_decoys": int(len(dec)),
        "kept_decoys": int(keep.sum()),
        "dropped_conflicting_decoys": int((~keep).sum()),
    }
    return dec.loc[keep].copy(), report


def build_panels() -> tuple[pd.DataFrame, pd.DataFrame]:
    tcr = pd.read_parquet(
        TCR_OUT / "tcr_registry.parquet",
        columns=[
            "row_id",
            "source_dataset",
            "peptide",
            "hla_4digit",
            "mhc_class",
            "cdr3_alpha",
            "cdr3_beta",
            "binding_label_binary",
            "paired_tcr_available",
            "is_cancer_context",
            "is_pathogen_context",
        ],
    )
    label = pd.to_numeric(tcr["binding_label_binary"], errors="coerce")
    pos = tcr[
        label.eq(1)
        & tcr["cdr3_beta"].fillna("").astype(str).str.fullmatch(AA_RE)
        & tcr["peptide"].fillna("").astype(str).str.fullmatch(PEP_RE)
        & tcr["hla_4digit"].fillna("").astype(str).str.contains(r"HLA-[ABC]\*", regex=True)
    ].copy()
    pos = pos.drop_duplicates(["source_dataset", "cdr3_beta", "peptide", "hla_4digit"]).reset_index(drop=True)
    pos["is_cancer_context"] = pos["is_cancer_context"].fillna(False).astype(bool)
    pos["is_pathogen_context"] = pos["is_pathogen_context"].fillna(False).astype(bool)
    pos["paired_tcr_available"] = pos["paired_tcr_available"].fillna(False).astype(bool)

    global_pos_triplets = set(zip(pos["cdr3_beta"].astype(str), pos["peptide"].astype(str), pos["hla_4digit"].astype(str)))
    global_pos_pairs = set(zip(pos["cdr3_beta"].astype(str), pos["peptide"].astype(str)))
    specs = [
        ("all_public", None, PANEL_POS_N),
        ("cancer_context", "is_cancer_context == True", PANEL_POS_N),
        ("paired_tcr", "paired_tcr_available == True", PANEL_POS_N),
        ("vdjdb", "source_dataset == 'VDJdb'", PANEL_POS_N),
        ("mcpas", "source_dataset == 'McPAS-TCR'", PANEL_POS_N),
        ("tenx", "source_dataset.str.contains('10x', na=False)", PANEL_POS_N),
    ]
    panel_frames: list[pd.DataFrame] = []
    report_rows: list[dict[str, object]] = []
    for i, (name, query, n) in enumerate(specs):
        panel_pos, rep = choose_panel(pos, name, query, n, SEED + i)
        dec, dec_rep = make_decoys(panel_pos, global_pos_triplets, global_pos_pairs, SEED + 1000 + i)
        panel = pd.concat(
            [
                panel_pos.assign(panel=name, label=1, pair_type="positive"),
                dec.assign(panel=name, label=0, pair_type="same_panel_shuffled_tcr_decoy"),
            ],
            ignore_index=True,
        )
        panel_frames.append(panel)
        report_rows.append({**rep, **dec_rep, "final_panel_rows": int(len(panel))})

    pairs = pd.concat(panel_frames, ignore_index=True)
    pairs = pairs.reset_index(drop=True)
    pairs["eval_id"] = [f"bench_{i:06d}" for i in range(len(pairs))]
    pairs["hla_for_pmtnet"] = pairs["hla_4digit"].map(clean_hla_for_pmtnet)
    pairs["pmhc_group"] = pairs["peptide"].astype(str) + "|" + pairs["hla_4digit"].astype(str)
    pairs["tcr_group"] = pairs["cdr3_beta"].astype(str)
    return pairs, pd.DataFrame(report_rows)


def run_pmtnet(scored: pd.DataFrame) -> pd.DataFrame:
    uniq = (
        scored[["cdr3_beta", "peptide", "hla_for_pmtnet"]]
        .drop_duplicates()
        .rename(columns={"cdr3_beta": "CDR3", "peptide": "Antigen", "hla_for_pmtnet": "HLA"})
    )
    input_path = OUTDIR / "pmtnet_input_unique.csv"
    uniq.to_csv(input_path, index=False)
    pred_path = OUTDIR / "pmtnet_run/prediction.csv"
    if pred_path.exists() and pred_path.stat().st_size > 0:
        pred = pd.read_csv(pred_path)
    else:
        rc = run(
            [
                "python",
                str(PMTNET / "pMTnet.py"),
                "-input",
                str(input_path),
                "-library",
                str(PMTNET / "library"),
                "-output",
                str(OUTDIR / "pmtnet_run"),
                "-output_log",
                str(OUTDIR / "pmtnet_output.log"),
            ],
            PMTNET,
            OUTDIR / "pmtnet_stdout.txt",
            OUTDIR / "pmtnet_stderr.txt",
        )
        if rc != 0:
            scored["pmtnet_score"] = np.nan
            return scored
        pred = pd.read_csv(pred_path)
    pred = pred.rename(columns={"CDR3": "cdr3_beta", "Antigen": "peptide", "HLA": "hla_for_pmtnet"})
    pred["pmtnet_rank"] = pd.to_numeric(pred["Rank"], errors="coerce")
    pred["pmtnet_score"] = 1.0 - pred["pmtnet_rank"]
    return scored.merge(pred[["cdr3_beta", "peptide", "hla_for_pmtnet", "pmtnet_rank", "pmtnet_score"]], on=["cdr3_beta", "peptide", "hla_for_pmtnet"], how="left")


def run_tepcam(scored: pd.DataFrame) -> pd.DataFrame:
    uniq = scored[["cdr3_beta", "peptide"]].drop_duplicates().reset_index(drop=True)
    tepcam_input = OUTDIR / "tepcam_input_unique.csv"
    pd.DataFrame({"TCR": uniq["cdr3_beta"], "epitope": uniq["peptide"], "Label": 0}).to_csv(tepcam_input, index=False)
    out_csv = OUTDIR / "tepcam_output_unique.csv"
    metrics_txt = OUTDIR / "tepcam_metrics_raw.txt"
    for p in [out_csv, metrics_txt]:
        if p.exists():
            p.unlink()
    rc = run(
        [
            "python",
            "scripts/test.py",
            "--file_path",
            str(tepcam_input),
            "--model_path",
            "./ckpts/tepcam_test.pt",
            "--output_file",
            str(out_csv),
            "--metric_file",
            str(metrics_txt),
            "--batch_size",
            "256",
        ],
        TEPCAM,
        OUTDIR / "tepcam_stdout.txt",
        OUTDIR / "tepcam_stderr.txt",
    )
    if rc != 0:
        scored["tepcam_score"] = np.nan
        return scored
    pred = pd.read_csv(out_csv)
    if len(pred) != len(uniq):
        scored["tepcam_score"] = np.nan
        return scored
    pred = pd.DataFrame(
        {
            "cdr3_beta": uniq["cdr3_beta"].to_numpy(),
            "peptide": uniq["peptide"].to_numpy(),
            "tepcam_score": pd.to_numeric(pred["prediction"], errors="coerce").to_numpy(),
        }
    )
    return scored.merge(pred[["cdr3_beta", "peptide", "tepcam_score"]], on=["cdr3_beta", "peptide"], how="left")


def add_scores(scored: pd.DataFrame) -> pd.DataFrame:
    scored["ensemble_mean_pmt_tepcam"] = scored[["pmtnet_score", "tepcam_score"]].mean(axis=1, skipna=True)
    both = scored[["pmtnet_score", "tepcam_score"]].notna().all(axis=1)
    scored["ensemble_logreg_groupcv"] = np.nan
    groups = scored.loc[both, "pmhc_group"].astype(str)
    if both.sum() >= 50 and groups.nunique() >= 3 and scored.loc[both, "label"].nunique() == 2:
        n_splits = min(5, groups.nunique())
        X = scored.loc[both, ["pmtnet_score", "tepcam_score"]].to_numpy(dtype=float)
        y = scored.loc[both, "label"].astype(int).to_numpy()
        oof = np.zeros(len(y), dtype=float)
        for train_idx, test_idx in GroupKFold(n_splits=n_splits).split(X, y, groups=groups):
            clf = LogisticRegression(max_iter=1000, class_weight="balanced")
            clf.fit(X[train_idx], y[train_idx])
            oof[test_idx] = clf.predict_proba(X[test_idx])[:, 1]
        scored.loc[both, "ensemble_logreg_groupcv"] = oof
    return scored


def summarize(scored: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    models = [
        ("pMTnet", "pmtnet_score", "TCRbeta + peptide + HLA"),
        ("TEPCAM", "tepcam_score", "TCRbeta + peptide"),
        ("mean_pMTnet_TEPCAM", "ensemble_mean_pmt_tepcam", "mean of available pMTnet/TEPCAM scores"),
        ("logreg_groupcv_pMTnet_TEPCAM", "ensemble_logreg_groupcv", "out-of-fold peptide-HLA GroupKFold ensemble"),
    ]
    rows = []
    for panel, sub in scored.groupby("panel"):
        for model, col, inputs in models:
            row = {"panel": panel, "model": model, "score_column": col, "inputs": inputs}
            row.update(metric_row(sub["label"], sub[col]))
            rows.append(row)
    for model, col, inputs in models:
        row = {"panel": "pooled", "model": model, "score_column": col, "inputs": inputs}
        row.update(metric_row(scored["label"], scored[col]))
        rows.append(row)
    metrics = pd.DataFrame(rows)

    corr_rows = []
    score_cols = ["pmtnet_score", "tepcam_score", "ensemble_mean_pmt_tepcam"]
    for panel, sub in scored.groupby("panel"):
        for i, a in enumerate(score_cols):
            for b in score_cols[i + 1 :]:
                valid = sub[[a, b]].apply(pd.to_numeric, errors="coerce").notna().all(axis=1)
                rho = spearmanr(sub.loc[valid, a], sub.loc[valid, b]).statistic if valid.sum() > 2 else np.nan
                corr_rows.append({"panel": panel, "score_a": a, "score_b": b, "n": int(valid.sum()), "spearman_rho": float(rho)})
    corrs = pd.DataFrame(corr_rows)

    source_rows = []
    both = scored[["pmtnet_score", "tepcam_score"]].notna().all(axis=1)
    usable = scored[both].copy()
    sources = [s for s, n in usable["source_dataset"].value_counts().items() if n >= 50]
    for source in sources:
        train = usable["source_dataset"].ne(source)
        test = usable["source_dataset"].eq(source)
        if train.sum() < 50 or test.sum() < 20 or usable.loc[train, "label"].nunique() < 2 or usable.loc[test, "label"].nunique() < 2:
            continue
        clf = LogisticRegression(max_iter=1000, class_weight="balanced")
        clf.fit(usable.loc[train, ["pmtnet_score", "tepcam_score"]].to_numpy(float), usable.loc[train, "label"].astype(int).to_numpy())
        pred = clf.predict_proba(usable.loc[test, ["pmtnet_score", "tepcam_score"]].to_numpy(float))[:, 1]
        row = {"heldout_source": source, "model": "source_heldout_logreg_pMTnet_TEPCAM"}
        row.update(metric_row(usable.loc[test, "label"], pd.Series(pred, index=usable.loc[test].index)))
        source_rows.append(row)
    source_metrics = pd.DataFrame(source_rows)
    return metrics, corrs, source_metrics


def write_report(scored: pd.DataFrame, panel_report: pd.DataFrame, metrics: pd.DataFrame, source_metrics: pd.DataFrame) -> None:
    pooled = metrics[metrics["panel"].eq("pooled")].copy()
    top = pooled.sort_values("auprc", ascending=False)
    lines = [
        "# External TCR Expert Challenge Benchmark",
        "",
        "Fast multi-panel benchmark using positive public TCR-pMHC rows and same-panel shuffled-TCR decoys.",
        "",
        "Claim boundary: this is a model-selection diagnostic. Decoys are synthetic and public-source overlap is possible; do not use as external SOTA evidence.",
        "",
        "## Pooled Result",
        "",
        "| Model | n | AUPRC | AUROC | Positive mean | Negative mean |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in top.itertuples(index=False):
        lines.append(f"| {r.model} | {r.n} | {r.auprc:.4f} | {r.auroc:.4f} | {r.positive_mean:.4f} | {r.negative_mean:.4f} |")
    lines.extend(["", "## Panel Construction", "", "| Panel | Eligible positives | Sampled positives | Kept decoys | Dropped decoys | Final rows |", "|---|---:|---:|---:|---:|---:|"])
    for r in panel_report.itertuples(index=False):
        lines.append(f"| {r.panel} | {r.eligible_positive_rows} | {r.sampled_positive_rows} | {r.kept_decoys} | {r.dropped_conflicting_decoys} | {r.final_panel_rows} |")
    lines.extend(["", "## Panel-Level Best Raw Expert", "", "| Panel | Best raw expert | AUPRC | AUROC |", "|---|---|---:|---:|"])
    raw = metrics[metrics["model"].isin(["pMTnet", "TEPCAM"])]
    for panel, sub in raw.groupby("panel"):
        best = sub.sort_values("auprc", ascending=False).iloc[0]
        lines.append(f"| {panel} | {best['model']} | {best['auprc']:.4f} | {best['auroc']:.4f} |")
    if not source_metrics.empty:
        lines.extend(["", "## Source-Heldout Ensemble", "", "| Heldout source | n | AUPRC | AUROC |", "|---|---:|---:|---:|"])
        for r in source_metrics.itertuples(index=False):
            lines.append(f"| {r.heldout_source} | {r.n} | {r.auprc:.4f} | {r.auroc:.4f} |")
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- Promote pMTnet as the first HLA-aware optional TCR expert.",
            "- Promote TEPCAM as a second sequence-only optional TCR expert.",
            "- Use the mean pMTnet/TEPCAM score as the first no-training diagnostic ensemble; only use logistic ensembles inside strict folds.",
            "- Keep all outputs out of the main pMHC model unless paired TCR data are present.",
        ]
    )
    (OUTDIR / "04_external_tcr_expert_benchmark.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "pmtnet_run").mkdir(parents=True, exist_ok=True)
    pairs, panel_report = build_panels()
    pairs.to_csv(OUTDIR / "external_tcr_expert_benchmark_pairs.tsv", sep="\t", index=False, na_rep="NA")
    panel_report.to_csv(OUTDIR / "external_tcr_expert_panel_report.tsv", sep="\t", index=False)
    scored = run_pmtnet(pairs)
    scored = run_tepcam(scored)
    scored = add_scores(scored)
    scored.to_csv(OUTDIR / "external_tcr_expert_benchmark_scores.tsv", sep="\t", index=False, na_rep="NA")
    metrics, corrs, source_metrics = summarize(scored)
    metrics.to_csv(OUTDIR / "external_tcr_expert_benchmark_metrics.tsv", sep="\t", index=False, na_rep="NA")
    corrs.to_csv(OUTDIR / "external_tcr_expert_score_correlations.tsv", sep="\t", index=False, na_rep="NA")
    source_metrics.to_csv(OUTDIR / "external_tcr_expert_source_heldout_metrics.tsv", sep="\t", index=False, na_rep="NA")
    write_report(scored, panel_report, metrics, source_metrics)
    pooled = metrics[metrics["panel"].eq("pooled")].sort_values("auprc", ascending=False)
    print(pooled[["model", "n", "n_pos", "n_neg", "auprc", "auroc"]].to_string(index=False))
    print(f"[external-tcr-expert-benchmark] rows={len(scored)} out={OUTDIR}")


if __name__ == "__main__":
    main()
