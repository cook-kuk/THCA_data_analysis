# GA/RL -> BAR-Neo-NG Bridge Method Card

## What It Does

This bridge takes GA/RL-discovered or GA/RL-ranked neoantigen candidates from BioDarwin and KG-GA outputs, maps them to CLEAN-NeoBench candidates by candidate ID or peptide-HLA pair, and applies BAR-Neo-NG claim gates.

## Claim Boundary

GA/RL fitness is treated as a discovery and prioritization signal. Clean claims require BAR-Neo-NG gate survival, leakage/source/HLA audit, public-training overlap resolution, and downstream metadata or assay evidence.

## Reviewer-Safe Use

- GA/RL high score plus BAR-Neo-NG T1: pHLA assay-design candidate.
- GA/RL high score plus BAR-Neo-NG blocked: useful failure/audit case, not a clean claim.
- Industrial public rows: post-freeze/manual-QA only, never clean training rows.
- MHC-II scout rows: separate benchmark required.
