#!/usr/bin/env python3
"""Prepare the next mutant/WT/decoy/control MD batch design.

This does not claim WT specificity. Rows with missing or inferred WT sequence
are explicitly marked as blocked or manual-confirmation required.
"""

from __future__ import annotations

import random
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "project/results/cross_neo_md_audit_2026_05_10"
CROSS = REPO / "project/results/cross_neo_v2_sota_sprint_2026_05_09"
DESIGN = OUT / "counterfactual_design"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def safe_id(text: object) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in str(text)).strip("_")


def preserve_anchor_scramble(peptide: str, seed: int = 17) -> str:
    """Scramble non-anchor positions while preserving P2 and C-terminal anchors."""
    pep = str(peptide).upper()
    if len(pep) < 5:
        return pep[::-1]
    rng = random.Random(seed + sum(ord(c) for c in pep))
    fixed = {1, len(pep) - 1}
    mutable_idx = [i for i in range(len(pep)) if i not in fixed]
    chars = [pep[i] for i in mutable_idx]
    for _ in range(100):
        rng.shuffle(chars)
        decoy = list(pep)
        for i, ch in zip(mutable_idx, chars):
            decoy[i] = ch
        out = "".join(decoy)
        if out != pep:
            return out
    return pep[::-1]


def tentative_wt(peptide: str, hla: str, antigen_sources: str) -> tuple[str, str]:
    sources = str(antigen_sources).upper()
    if peptide == "GADGVGKSAL" and "KRAS" in sources:
        return "GAGGVGKSAL", "inferred_kras_g12d_like_requires_manual_confirmation"
    return "", "missing_wt_sequence_blocked"


def top_template(ready: pd.DataFrame, peptide: str, hla: str) -> tuple[str, str]:
    if ready.empty:
        return "", ""
    sub = ready[(ready["target_peptide"].astype(str) == peptide) & (ready["hla_4digit"].astype(str) == hla)].copy()
    if sub.empty:
        return "", ""
    if "contacts_tcr_peptide" in sub.columns:
        sub = sub.sort_values("contacts_tcr_peptide", ascending=False)
    row = sub.iloc[0]
    return str(row.get("pdb_id", "")), str(row.get("extracted_pdb", ""))


def add_row(rows: list[dict], **kwargs) -> None:
    defaults = {
        "batch_id": "",
        "priority_rank": "",
        "row_id": "",
        "peptide": "",
        "hla_4digit": "",
        "control_type": "",
        "complex_kind": "",
        "sequence": "",
        "sequence_status": "",
        "template_pdb_id": "",
        "template_pdb_path": "",
        "replicates": 3,
        "ns_per_replicate": 10,
        "total_ns": 30,
        "readiness": "",
        "blocked_by": "",
        "claim_status": "structural_audit_only_not_immunogenicity_proof",
        "recommended_after": "",
        "rationale": "",
    }
    defaults.update(kwargs)
    rows.append(defaults)


def main() -> None:
    DESIGN.mkdir(parents=True, exist_ok=True)
    queue = read_tsv(CROSS / "tcr_extension/md_escalation/md_escalation_queue_top20.tsv")
    wetlab = read_tsv(CROSS / "tcr_extension/tcr_wetlab_candidate_prioritization_unique_pmhc.tsv")
    ready = read_tsv(CROSS / "tcr_extension/md_escalation/p0_structures/p0_md_pilot_ready_complexes.tsv")
    if queue.empty:
        raise SystemExit("Missing md escalation queue")

    wetlab_key = {}
    if not wetlab.empty:
        for _, r in wetlab.iterrows():
            wetlab_key[(str(r.get("peptide", "")), str(r.get("hla_4digit", "")))] = r.to_dict()

    rows: list[dict] = []
    control_rows: list[dict] = []
    for rank, (_, q) in enumerate(queue.head(12).iterrows(), start=1):
        peptide = str(q.get("peptide", ""))
        hla = str(q.get("hla_4digit", ""))
        row_id = str(q.get("row_id", ""))
        tier = str(q.get("md_tier", ""))
        wt_info = wetlab_key.get((peptide, hla), {})
        wt_seq = "" if pd.isna(wt_info.get("wildtype_peptide", "")) else str(wt_info.get("wildtype_peptide", ""))
        if not wt_seq or wt_seq.lower() == "nan":
            wt_seq, wt_status = tentative_wt(peptide, hla, wt_info.get("example_tcr_antigen_sources", ""))
        else:
            wt_status = "registry_wt_sequence_available"
        decoy = preserve_anchor_scramble(peptide)
        pdb_id, template = top_template(ready, peptide, hla)
        has_tcr_template = bool(template)
        base = {
            "priority_rank": rank,
            "row_id": row_id,
            "peptide": peptide,
            "hla_4digit": hla,
            "template_pdb_id": pdb_id,
            "template_pdb_path": template,
        }
        p0 = tier.startswith("P0")

        add_row(
            rows,
            **base,
            batch_id=f"P{rank:02d}_{safe_id(peptide)}_{safe_id(hla)}_mutant_pmhc_3x10ns",
            control_type="mutant",
            complex_kind="pMHC",
            sequence=peptide,
            sequence_status="observed_mutant_or_candidate_peptide",
            readiness="ready_template_or_structure_generation_needed",
            blocked_by="" if p0 else "requires_pmhc_structure_generation",
            rationale="baseline peptide-MHC stability audit",
        )
        if p0:
            add_row(
                rows,
                **base,
                batch_id=f"P{rank:02d}_{safe_id(peptide)}_{safe_id(hla)}_mutant_tcr_pmhc_3x10ns",
                control_type="mutant",
                complex_kind="TCR-pMHC",
                sequence=peptide,
                sequence_status="observed_mutant_or_candidate_peptide",
                readiness="ready_existing_tcr_pmhc_template" if has_tcr_template else "blocked",
                blocked_by="" if has_tcr_template else "no_tcr_pmhc_template",
                rationale="recognition-interface stability audit",
            )

        for complex_kind in (["pMHC", "TCR-pMHC"] if p0 else ["pMHC"]):
            add_row(
                rows,
                **base,
                batch_id=f"P{rank:02d}_{safe_id(peptide)}_{safe_id(hla)}_wt_{safe_id(complex_kind)}_3x10ns",
                control_type="wildtype",
                complex_kind=complex_kind,
                sequence=wt_seq,
                sequence_status=wt_status,
                readiness="manual_confirmation_required" if wt_seq else "blocked",
                blocked_by="confirm_wt_sequence_before_simulation" if wt_seq else "missing_wt_sequence",
                rationale="required for mutant-WT specificity and WT cross-reactivity risk",
            )

        for complex_kind in (["pMHC", "TCR-pMHC"] if p0 else ["pMHC"]):
            add_row(
                rows,
                **base,
                batch_id=f"P{rank:02d}_{safe_id(peptide)}_{safe_id(hla)}_anchor_preserved_decoy_{safe_id(complex_kind)}_3x10ns",
                control_type="anchor_preserved_scrambled_decoy",
                complex_kind=complex_kind,
                sequence=decoy,
                sequence_status="algorithmic_decoy_preserving_p2_and_cterminal_anchor",
                readiness="structure_generation_required",
                blocked_by="prepare_decoy_structure_from_template",
                rationale="specificity control; instability or lost contacts should separate decoy from candidate",
            )

        for complex_kind in (["pMHC", "TCR-pMHC"] if p0 else ["pMHC"]):
            add_row(
                rows,
                **base,
                batch_id=f"P{rank:02d}_{safe_id(peptide)}_{safe_id(hla)}_same_hla_positive_control_{safe_id(complex_kind)}_3x10ns",
                control_type="same_or_similar_hla_positive_control",
                complex_kind=complex_kind,
                sequence="",
                sequence_status="not_selected",
                readiness="blocked",
                blocked_by="select_nonidentical_known_positive_same_or_similar_hla_control",
                rationale="runtime and analysis sanity control; should not be used as cancer specificity proof",
            )

        control_rows.append(
            {
                "row_id": row_id,
                "peptide": peptide,
                "hla_4digit": hla,
                "candidate_sequence": peptide,
                "tentative_wt_sequence": wt_seq,
                "wt_status": wt_status,
                "anchor_preserved_decoy": decoy,
                "decoy_rule": "preserve P2 and C-terminal anchor, shuffle remaining residues",
                "template_pdb_id": pdb_id,
                "template_pdb_path": template,
            }
        )

    manifest = pd.DataFrame(rows)
    controls = pd.DataFrame(control_rows)
    manifest.to_csv(DESIGN / "counterfactual_md_batch_manifest.tsv", sep="\t", index=False)
    controls.to_csv(DESIGN / "candidate_control_sequences.tsv", sep="\t", index=False)
    with (DESIGN / "candidate_decoy_controls.fasta").open("w") as handle:
        for _, r in controls.iterrows():
            handle.write(f">{r['row_id']}|{r['peptide']}|{r['hla_4digit']}|candidate\n{r['candidate_sequence']}\n")
            if r["tentative_wt_sequence"]:
                handle.write(f">{r['row_id']}|{r['peptide']}|{r['hla_4digit']}|wt_manual_confirm\n{r['tentative_wt_sequence']}\n")
            handle.write(f">{r['row_id']}|{r['peptide']}|{r['hla_4digit']}|anchor_preserved_decoy\n{r['anchor_preserved_decoy']}\n")

    ready_rows = int((manifest["readiness"].astype(str).str.contains("ready", na=False)).sum())
    blocked_rows = int((manifest["readiness"].astype(str).str.contains("blocked", na=False)).sum())
    report = [
        "# Counterfactual MD Batch Design",
        "",
        "## Boundary",
        "",
        "This is a simulation-design manifest. It does not establish mutant-specific TCR recognition. WT rows with inferred sequences are marked for manual confirmation before any run.",
        "",
        "## Summary",
        "",
        f"- candidate/control jobs: {len(manifest)}",
        f"- ready or near-ready rows: {ready_rows}",
        f"- blocked rows: {blocked_rows}",
        f"- candidate sequence rows: {len(controls)}",
        "",
        "## Immediate P0 Plan",
        "",
        "- `HMTEVVRHC / HLA-A*02:01`: finish current 10 ns, sync DCD, rerun full analysis, then add WT/decoy controls once WT is confirmed.",
        "- `GADGVGKSAL / HLA-C*08:02`: add WT confirmation, anchor-preserved decoy, and same-HLA positive-control selection before long MD.",
        "",
        "## RunPod Submission Template",
        "",
        "Use this only after structure preparation has produced PDB inputs for each manifest row.",
        "",
        "```bash",
        "python run_openmm_pilot.py --input INPUT.pdb --outdir OUTDIR --mode explicit --platform CUDA --ns 10 --report-steps 50000 --minimize-iterations 1000 --timestep-fs 2.0 --temperature-k 300",
        "```",
        "",
        "## Outputs",
        "",
        "- `counterfactual_md_batch_manifest.tsv`",
        "- `candidate_control_sequences.tsv`",
        "- `candidate_decoy_controls.fasta`",
    ]
    (DESIGN / "counterfactual_md_batch_plan.md").write_text("\n".join(report) + "\n")
    print(f"[md-counterfactual-design] jobs={len(manifest)} ready_like={ready_rows} blocked={blocked_rows}")
    print(DESIGN / "counterfactual_md_batch_manifest.tsv")


if __name__ == "__main__":
    main()
