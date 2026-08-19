#!/usr/bin/env python3
"""Run a fast pMTnet pilot on local TCR-pMHC positives versus shuffled decoys."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from common import OUT, SEED


TCR_OUT = OUT / "tcr_extension"
DISC_OUT = TCR_OUT / "model_discovery"
PMTNET = Path("/data/thca/tcr_model_repos/pMTnet")


def clean_hla(hla: object) -> str:
    text = str(hla or "").strip()
    text = text.replace("HLA-", "")
    return text


def main() -> None:
    out = DISC_OUT / "pmtnet_cross_neo_pilot"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    tcr = pd.read_parquet(
        TCR_OUT / "tcr_registry.parquet",
        columns=[
            "row_id",
            "source_dataset",
            "peptide",
            "hla_4digit",
            "mhc_class",
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
        & tcr["cdr3_beta"].fillna("").astype(str).str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")
        & tcr["peptide"].fillna("").astype(str).str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{8,15}")
        & tcr["hla_4digit"].fillna("").astype(str).str.contains(r"HLA-[ABC]\*", regex=True)
    ].copy()
    pos = pos.drop_duplicates(["cdr3_beta", "peptide", "hla_4digit"])
    if len(pos) > 200:
        pos = pos.sample(n=200, random_state=SEED)
    dec = pos.copy()
    shuffled = dec["cdr3_beta"].to_numpy().copy()
    for _ in range(20):
        rng.shuffle(shuffled)
        if not np.any(shuffled == dec["cdr3_beta"].to_numpy()):
            break
    dec["cdr3_beta"] = shuffled
    dec["row_id"] = dec["row_id"].astype(str) + "_shuffled_decoy"

    eval_df = pd.concat(
        [
            pos.assign(label=1, pair_type="positive"),
            dec.assign(label=0, pair_type="same_pmhc_shuffled_tcr_decoy"),
        ],
        ignore_index=True,
    )
    eval_df["hla_for_pmtnet"] = eval_df["hla_4digit"].map(clean_hla)
    eval_df = (
        eval_df.sort_values("label", ascending=False)
        .drop_duplicates(["cdr3_beta", "peptide", "hla_for_pmtnet"], keep="first")
        .sample(frac=1.0, random_state=SEED)
        .reset_index(drop=True)
    )
    input_df = pd.DataFrame(
        {
            "CDR3": eval_df["cdr3_beta"].astype(str),
            "Antigen": eval_df["peptide"].astype(str),
            "HLA": eval_df["hla_for_pmtnet"],
        }
    )
    eval_df.to_csv(out / "pmtnet_pilot_pairs.tsv", sep="\t", index=False, na_rep="NA")
    input_df.to_csv(out / "pmtnet_input.csv", index=False)

    cmd = [
        "python",
        str(PMTNET / "pMTnet.py"),
        "-input",
        str(out / "pmtnet_input.csv"),
        "-library",
        str(PMTNET / "library"),
        "-output",
        str(out),
        "-output_log",
        str(out / "pmtnet_output.log"),
    ]
    proc = subprocess.run(cmd, cwd=str(PMTNET), text=True, capture_output=True, check=False)
    (out / "pmtnet_stdout.txt").write_text(proc.stdout)
    (out / "pmtnet_stderr.txt").write_text(proc.stderr)
    if proc.returncode != 0:
        raise SystemExit(f"pMTnet failed with code {proc.returncode}; see {out}")

    pred = pd.read_csv(out / "prediction.csv")
    scored = eval_df.merge(
        pred.rename(columns={"CDR3": "cdr3_beta", "Antigen": "peptide", "HLA": "hla_for_pmtnet"}),
        on=["cdr3_beta", "peptide", "hla_for_pmtnet"],
        how="left",
    )
    scored["pmtnet_rank"] = pd.to_numeric(scored["Rank"], errors="coerce")
    scored["pmtnet_score"] = 1.0 - scored["pmtnet_rank"]
    scored.to_csv(out / "pmtnet_pilot_scored.tsv", sep="\t", index=False, na_rep="NA")
    valid = scored["pmtnet_score"].notna()
    y = scored.loc[valid, "label"].astype(int).to_numpy()
    s = scored.loc[valid, "pmtnet_score"].astype(float).to_numpy()
    metrics = {
        "model": "pMTnet",
        "pilot": "positive_vs_same_pmhc_shuffled_tcr_decoy",
        "n": int(valid.sum()),
        "n_pos": int(y.sum()),
        "n_neg": int(len(y) - y.sum()),
        "auprc": float(average_precision_score(y, s)) if len(set(y)) > 1 else np.nan,
        "auroc": float(roc_auc_score(y, s)) if len(set(y)) > 1 else np.nan,
        "positive_score_mean": float(scored.loc[valid & scored["label"].eq(1), "pmtnet_score"].mean()),
        "decoy_score_mean": float(scored.loc[valid & scored["label"].eq(0), "pmtnet_score"].mean()),
        "positive_rank_mean": float(scored.loc[valid & scored["label"].eq(1), "pmtnet_rank"].mean()),
        "decoy_rank_mean": float(scored.loc[valid & scored["label"].eq(0), "pmtnet_rank"].mean()),
    }
    pd.DataFrame([metrics]).to_csv(out / "pmtnet_pilot_metrics.tsv", sep="\t", index=False)
    lines = [
        "# pMTnet CROSS-Neo-TCR Pilot",
        "",
        "Fast pilot using local pMTnet on positive public TCR-pMHC rows versus same peptide-HLA shuffled-TCR decoys.",
        "",
        "| Metric | Value |",
        "|---|---:|",
    ]
    for k, v in metrics.items():
        if k in {"model", "pilot"}:
            lines.append(f"| {k} | {v} |")
        else:
            lines.append(f"| {k} | {v:.4f} |" if isinstance(v, float) else f"| {k} | {v} |")
    lines.extend(
        [
            "",
            "Interpretation: lower pMTnet rank means stronger predicted binding, so `pmtnet_score = 1 - rank` is used for AUPRC/AUROC.",
            "This is still a pilot because shuffled decoys are synthetic and public-source leakage has not been fully eliminated.",
        ]
    )
    (out / "pmtnet_pilot_report.md").write_text("\n".join(lines) + "\n")
    print(f"[pmtnet-pilot] n={metrics['n']} auprc={metrics['auprc']:.4f} auroc={metrics['auroc']:.4f} out={out}")


if __name__ == "__main__":
    main()
