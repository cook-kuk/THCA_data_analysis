#!/usr/bin/env python3
"""Compare runnable external TCR models on the same fast pilot pairs."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

from common import OUT


TCR_OUT = OUT / "tcr_extension"
DISC_OUT = TCR_OUT / "model_discovery"
PILOT_OUT = DISC_OUT / "external_tcr_model_pilot_compare"
PMTNET_OUT = DISC_OUT / "pmtnet_cross_neo_pilot"
PANPEP = Path("/data/thca/tcr_model_repos/PanPep")
TEPCAM = Path("/data/thca/tcr_model_repos/TEPCAM")
ERGO = Path("/data/thca/tcr_model_repos/ERGO-II")


def run(cmd: list[str], cwd: Path, stdout_path: Path, stderr_path: Path) -> int:
    proc = subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True, check=False)
    stdout_path.write_text(proc.stdout)
    stderr_path.write_text(proc.stderr)
    return int(proc.returncode)


def score_metrics(y: pd.Series, score: pd.Series) -> dict[str, float]:
    yy = y.astype(int).to_numpy()
    ss = pd.to_numeric(score, errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(ss)
    yy, ss = yy[mask], ss[mask]
    out: dict[str, float] = {
        "n": int(len(yy)),
        "n_pos": int(yy.sum()) if len(yy) else 0,
        "n_neg": int(len(yy) - yy.sum()) if len(yy) else 0,
        "positive_score_mean": float(np.mean(ss[yy == 1])) if np.any(yy == 1) else np.nan,
        "decoy_score_mean": float(np.mean(ss[yy == 0])) if np.any(yy == 0) else np.nan,
    }
    if len(yy) and len(np.unique(yy)) == 2:
        out["auprc"] = float(average_precision_score(yy, ss))
        out["auroc"] = float(roc_auc_score(yy, ss))
    else:
        out["auprc"] = np.nan
        out["auroc"] = np.nan
    return out


def ergo_mhc(hla: object) -> str:
    text = str(hla or "")
    if "*" not in text:
        return "UNK"
    gene, allele = text.split("*", 1)
    allele = allele.split(":", 1)[0]
    return f"{gene}*{allele}"


def main() -> None:
    PILOT_OUT.mkdir(parents=True, exist_ok=True)

    pairs = pd.read_csv(PMTNET_OUT / "pmtnet_pilot_pairs.tsv", sep="\t")
    pairs = pairs.reset_index(drop=True).assign(eval_id=lambda d: [f"pilot_{i:05d}" for i in range(len(d))])
    pairs["label"] = pairs["label"].astype(int)

    pmtnet = pd.read_csv(PMTNET_OUT / "pmtnet_pilot_scored.tsv", sep="\t")
    pmtnet = pmtnet[["cdr3_beta", "peptide", "hla_4digit", "pmtnet_score"]].copy()
    scored = pairs.merge(pmtnet, on=["cdr3_beta", "peptide", "hla_4digit"], how="left")

    # PanPep is peptide+TCRbeta only; keep input sorted to match PanPep's grouped output order.
    panpep_in = (
        scored[["eval_id", "peptide", "cdr3_beta"]]
        .rename(columns={"peptide": "Peptide", "cdr3_beta": "CDR3"})
        .sort_values(["Peptide", "CDR3", "eval_id"])
        .reset_index(drop=True)
    )
    panpep_input = PILOT_OUT / "panpep_input_sorted.csv"
    panpep_output = PILOT_OUT / "panpep_output.csv"
    panpep_in[["Peptide", "CDR3"]].to_csv(panpep_input, index=False)
    panpep_rc = run(
        [
            "python",
            "PanPep.py",
            "--learning_setting",
            "zero-shot",
            "--input",
            str(panpep_input),
            "--output",
            str(panpep_output),
        ],
        PANPEP,
        PILOT_OUT / "panpep_stdout.txt",
        PILOT_OUT / "panpep_stderr.txt",
    )
    if panpep_rc == 0:
        panpep_pred = pd.read_csv(panpep_output)
        if len(panpep_pred) == len(panpep_in):
            panpep_pred["eval_id"] = panpep_in["eval_id"].to_numpy()
            scored = scored.merge(
                panpep_pred[["eval_id", "Score"]].rename(columns={"Score": "panpep_score"}),
                on="eval_id",
                how="left",
            )
        else:
            scored["panpep_score"] = np.nan
    else:
        scored["panpep_score"] = np.nan

    tepcam_in = scored[["cdr3_beta", "peptide", "label"]].rename(
        columns={"cdr3_beta": "TCR", "peptide": "epitope", "label": "Label"}
    )
    tepcam_input = PILOT_OUT / "tepcam_input.csv"
    tepcam_output = PILOT_OUT / "tepcam_output.csv"
    tepcam_metrics_txt = PILOT_OUT / "tepcam_metrics_raw.txt"
    tepcam_in.to_csv(tepcam_input, index=False)
    tepcam_rc = run(
        [
            "python",
            "scripts/test.py",
            "--file_path",
            str(tepcam_input),
            "--model_path",
            "./ckpts/tepcam_test.pt",
            "--output_file",
            str(tepcam_output),
            "--metric_file",
            str(tepcam_metrics_txt),
            "--batch_size",
            "128",
        ],
        TEPCAM,
        PILOT_OUT / "tepcam_stdout.txt",
        PILOT_OUT / "tepcam_stderr.txt",
    )
    if tepcam_rc == 0:
        tepcam_pred = pd.read_csv(tepcam_output)
        if len(tepcam_pred) == len(scored):
            scored["tepcam_score"] = pd.to_numeric(tepcam_pred["prediction"], errors="coerce").to_numpy()
        else:
            scored["tepcam_score"] = np.nan
    else:
        scored["tepcam_score"] = np.nan

    ergo_input = PILOT_OUT / "ergo_input.csv"
    ergo_output = PILOT_OUT / "ergo_vdjdb_output.csv"
    ergo_df = pd.DataFrame(
        {
            "TRA": np.nan,
            "TRB": scored["cdr3_beta"],
            "TRAV": np.nan,
            "TRAJ": np.nan,
            "TRBV": np.nan,
            "TRBJ": np.nan,
            "T-Cell-Type": np.where(scored["mhc_class"].astype(str).str.upper().eq("II"), "MHCII", "MHCI"),
            "Peptide": scored["peptide"],
            "MHC": scored["hla_4digit"].map(ergo_mhc),
        }
    )
    ergo_df.to_csv(ergo_input, index=False)
    ergo_rc = run(
        ["python", "Predict.py", "vdjdb", str(ergo_input), str(ergo_output)],
        ERGO,
        PILOT_OUT / "ergo_stdout.txt",
        PILOT_OUT / "ergo_stderr.txt",
    )
    if ergo_rc == 0:
        ergo_pred = pd.read_csv(ergo_output)
        if len(ergo_pred) == len(scored):
            scored["ergo_vdjdb_score"] = pd.to_numeric(ergo_pred["Score"], errors="coerce").to_numpy()
        else:
            scored["ergo_vdjdb_score"] = np.nan
    else:
        scored["ergo_vdjdb_score"] = np.nan

    models = [
        ("pMTnet", "pmtnet_score", "TCRbeta + peptide + HLA", "P0 runnable; pilot claim-limited"),
        ("PanPep", "panpep_score", "TCRbeta + peptide", "P1 runnable; no HLA in score"),
        ("TEPCAM", "tepcam_score", "TCRbeta + peptide", "P1 runnable; no HLA in score"),
        ("ERGO-II-vdjdb", "ergo_vdjdb_score", "TCRbeta + peptide + coarse MHC", "P1 runnable; adapter is beta-only/unknown-VJ"),
    ]
    rows = []
    for model, col, inputs, claim in models:
        row = {"model": model, "score_column": col, "inputs": inputs, "claim_boundary": claim}
        row.update(score_metrics(scored["label"], scored[col]))
        rows.append(row)
    comp = pd.DataFrame(rows)
    comp.to_csv(PILOT_OUT / "external_tcr_model_pilot_comparison.tsv", sep="\t", index=False, na_rep="NA")

    corr_rows = []
    score_cols = [m[1] for m in models]
    for i, a in enumerate(score_cols):
        for b in score_cols[i + 1 :]:
            valid = scored[[a, b]].apply(pd.to_numeric, errors="coerce").notna().all(axis=1)
            rho = spearmanr(scored.loc[valid, a], scored.loc[valid, b]).statistic if valid.sum() > 2 else np.nan
            corr_rows.append({"score_a": a, "score_b": b, "n": int(valid.sum()), "spearman_rho": float(rho)})
    pd.DataFrame(corr_rows).to_csv(PILOT_OUT / "external_tcr_model_score_correlations.tsv", sep="\t", index=False, na_rep="NA")
    scored.to_csv(PILOT_OUT / "external_tcr_model_pilot_scores.tsv", sep="\t", index=False, na_rep="NA")

    lines = [
        "# External TCR Model Pilot Comparison",
        "",
        "Same positive public TCR-pMHC rows versus same peptide-HLA shuffled-TCR decoys used for the pMTnet pilot.",
        "",
        "| Model | Inputs | n | AUPRC | AUROC | Positive mean | Decoy mean | Claim boundary |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['model']} | {row['inputs']} | {row['n']} | {row['auprc']:.4f} | {row['auroc']:.4f} | "
            f"{row['positive_score_mean']:.4f} | {row['decoy_score_mean']:.4f} | {row['claim_boundary']} |"
        )
    lines.extend(
        [
            "",
            "Interpretation boundary: these are fast pilot numbers on synthetic shuffled decoys, not a clean external benchmark.",
            "pMTnet has an input advantage because it uses HLA; TEPCAM and PanPep are TCRbeta-peptide models, and the ERGO-II pilot used only coarse MHC with unknown V/J and no alpha chain.",
        ]
    )
    (PILOT_OUT / "03_external_tcr_model_pilot_comparison.md").write_text("\n".join(lines) + "\n")
    print(comp.to_string(index=False))
    print(f"[external-model-pilot] out={PILOT_OUT}")


if __name__ == "__main__":
    main()
