# Blinded Data Request

Please provide a blinded candidate table with the following where available:

- patient_id
- mutant_peptide
- optional WT peptide
- HLA allele
- expression / RNA support
- variant context
- payload slot budget / intended rank window
- sponsor baseline score if it can remain blinded until lock
- eventual labels withheld until score freeze

Optional controls:
- decoys
- self-like peptides
- known non-binders
- wetlab-negative historical candidates

Deliverable after freeze:
- per-candidate score
- per-patient top-N ranking
- uncertainty flags
- reason codes
- assay-control annotations

