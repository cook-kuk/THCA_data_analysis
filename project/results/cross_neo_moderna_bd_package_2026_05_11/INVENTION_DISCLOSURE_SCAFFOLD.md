# Invention Disclosure Scaffold

Not legal advice. Use this as a structured intake for counsel or technology transfer before any external outreach.

## Working Title

Knowledge-graph-guided evolutionary selection of individualized cancer vaccine neoantigens with audit-controlled wetlab translation.

## Problem

Individualized mRNA vaccine payloads have limited antigen slots. Candidate pools can contain many plausible peptides, but existing presentation or immunogenicity scores do not fully solve false-positive burden, explainability, patient-level slot allocation, or wetlab-control design.

## Technical Solution

- Use a knowledge graph to encode biological, assay, HLA, presentation, immunogenicity, self-similarity, TCR, and translation constraints.
- Use a genetic/evolutionary controller to search candidate scoring/selection policies.
- Maintain separated production-priority and claim-safe scoring lanes.
- Emit per-candidate reason codes, uncertainty/audit flags, and top-N wetlab-control maps.

## Candidate Claim Areas

- KG-guided evolutionary optimization for neoantigen ranking under a fixed payload slot budget.
- Dual-lane prioritization: production/experiment-priority score plus claim-safe fallback score.
- Blinded validation workflow with frozen score hashes before label unblinding.
- Automated mutant/WT/decoy assay-control selection tied to candidate ranking.
- Patient-level top-N payload selection with HLA/source diversity constraints and uncertainty flags.

## Evidence Snapshot

- Same-board 2,715-candidate KG-GA: AUPRC 0.933, AUROC 0.943.
- Frozen validation-like subset: AUPRC 0.978.
- Current strict no-overlap reliability caveat: BAR-Neo_confidence AUPRC 0.857 vs KG-GA-controller 0.624.

## Do Not Disclose Before Filing

- Source code, exact feature list, scoring weights, graph schema details, candidate-level score tables, training corpus joins, or endpoint labels.

