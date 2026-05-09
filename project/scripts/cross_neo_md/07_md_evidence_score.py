#!/usr/bin/env python3
"""Rule-based MD evidence scoring for CROSS-Neo structural audit."""

from __future__ import annotations

import numpy as np
import pandas as pd

from common_md import OUT, ensure_dirs, write_markdown_table


def clip01(x):
    if pd.isna(x):
        return 0.0
    return float(max(0.0, min(1.0, x)))


def label(score, complete_fraction):
    if complete_fraction < 0.1:
        return "MD_INSUFFICIENT_RUNTIME"
    if score >= 0.82:
        return "MD_VERY_STRONG"
    if score >= 0.68:
        return "MD_STRONG"
    if score >= 0.50:
        return "MD_MODERATE"
    if score >= 0.30:
        return "MD_WEAK"
    return "MD_FAIL"


def main() -> None:
    ensure_dirs()
    qc = pd.read_csv(OUT / "qc/md_qc_summary.tsv", sep="\t") if (OUT / "qc/md_qc_summary.tsv").exists() else pd.DataFrame()
    pmhc_ts = pd.read_csv(OUT / "contacts/pmhc_contact_timeseries.tsv", sep="\t") if (OUT / "contacts/pmhc_contact_timeseries.tsv").exists() else pd.DataFrame()
    tcr_ts = pd.read_csv(OUT / "contacts/tcr_peptide_contact_timeseries.tsv", sep="\t") if (OUT / "contacts/tcr_peptide_contact_timeseries.tsv").exists() else pd.DataFrame()
    repl = pd.read_csv(OUT / "replicates/replicate_consistency_summary.tsv", sep="\t") if (OUT / "replicates/replicate_consistency_summary.tsv").exists() else pd.DataFrame()
    rows = []
    for _, r in qc.iterrows():
        if not bool(r.get("trajectory_available", False)):
            continue
        run_id = r["run_id"]
        pmhc_q = np.nan
        if not pmhc_ts.empty:
            sub = pmhc_ts[pmhc_ts["run_id"].eq(run_id)]
            if not sub.empty:
                pmhc_q = float(sub["native_contact_fraction_q"].tail(max(1, len(sub) // 4)).mean())
        tcr_contacts = np.nan
        if not tcr_ts.empty:
            sub = tcr_ts[tcr_ts["run_id"].eq(run_id)]
            if not sub.empty:
                tcr_contacts = float(sub["tcr_peptide_heavy_atom_contacts"].tail(max(1, len(sub) // 4)).mean())
        peptide_rmsd = pd.to_numeric(pd.Series([r.get("peptide_rmsd_final_nm")]), errors="coerce").iloc[0]
        drift = pd.to_numeric(pd.Series([r.get("peptide_com_drift_final_nm")]), errors="coerce").iloc[0]
        complete = float(r.get("trajectory_complete_fraction", 0) or 0)
        p_stability = 0.45 * clip01(1 - peptide_rmsd / 0.5) + 0.30 * clip01(1 - drift / 0.6) + 0.25 * clip01(pmhc_q)
        tcr_score = clip01(tcr_contacts / 50.0) if not pd.isna(tcr_contacts) else 0.0
        sim_qc = 0.0 if bool(r.get("nan_coordinates", False)) else clip01(complete)
        rep_score = 0.2
        if not repl.empty:
            rr = repl[repl["candidate"].eq(r["candidate"])]
            if not rr.empty:
                cat = str(rr.iloc[0].get("interpretation_category", ""))
                rep_score = {"robust_stable": 1.0, "replicate_inconsistent": 0.4, "unstable_peptide": 0.1, "insufficient_replicates": 0.2}.get(cat, 0.2)
        counter = 0.0
        score = 0.38 * p_stability + 0.27 * tcr_score + 0.20 * sim_qc + 0.10 * rep_score + 0.05 * counter
        rows.append(
            {
                "run_id": run_id,
                "candidate": r["candidate"],
                "condition": r["condition"],
                "pMHC_stability_score": p_stability,
                "TCR_recognition_score": tcr_score,
                "counterfactual_specificity_score": counter,
                "replicate_confidence_score": rep_score,
                "simulation_qc_score": sim_qc,
                "MD_evidence_score": score,
                "MD_evidence_label": label(score, complete),
                "runtime_fraction": complete,
                "peptide_rmsd_final_nm": peptide_rmsd,
                "pmhc_native_q_tail": pmhc_q,
                "tcr_peptide_contacts_tail": tcr_contacts,
                "claim_boundary": "Structural audit only; not immunogenicity proof.",
            }
        )
    scores = pd.DataFrame(rows)
    scores.to_csv(OUT / "md_evidence_scores.tsv", sep="\t", index=False)
    for cand, sub in scores.groupby("candidate") if not scores.empty else []:
        path = OUT / f"md_evidence_score_card_{cand.replace('/', '_').replace('*','').replace(':','')}.md"
        lines = [
            f"# MD Evidence Score Card: {cand}",
            "",
            write_markdown_table(sub, ["run_id", "condition", "MD_evidence_score", "MD_evidence_label", "pMHC_stability_score", "TCR_recognition_score", "simulation_qc_score", "runtime_fraction"], max_rows=20),
            "",
            "Boundary: this score is rule-based and interpretable. It is not trained and cannot prove immunogenicity.",
        ]
        path.write_text("\n".join(lines) + "\n")
    lines = [
        "# MD Evidence Score Report",
        "",
        "Weights are transparent heuristics: pMHC stability, TCR contacts, runtime/QC, replicate consistency, and counterfactual specificity. Counterfactual specificity is zero until WT/decoy trajectories exist.",
        "",
        write_markdown_table(scores, ["run_id", "candidate", "MD_evidence_score", "MD_evidence_label", "pMHC_stability_score", "TCR_recognition_score", "simulation_qc_score"], max_rows=100),
    ]
    (OUT / "07_md_evidence_score_report.md").write_text("\n".join(lines) + "\n")
    print(f"[md-score] rows={len(scores)} out={OUT}")
    if not scores.empty:
        print(scores[["run_id", "candidate", "MD_evidence_score", "MD_evidence_label"]].to_string(index=False))


if __name__ == "__main__":
    main()
