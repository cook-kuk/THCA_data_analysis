# PRIME 2.1 scoring on ITSNdb (leakage-stratified)

## Install
- PRIME 2.1 cloned from https://github.com/GfellerLab/PRIME
- Bundled `lib/PRIME.x` was Mac arm64; recompiled for Linux x86_64 with `g++ -O3 PRIME.cc -o PRIME.x`.
- MixMHCpred 3.0 cloned from https://github.com/GfellerLab/MixMHCpred (PRIME dependency).
- Python deps: logomaker installed via pip; numpy/pandas/scipy already present.
- MAFFT installed via `apt-get install -y mafft` (used by MixMHCpred for new-allele alignment).

## HLA format conversion
- Input bundle: `HLA-A*02:01`
- PRIME format: `A0201` (drop `HLA-`, `*`, `:`).

## Allele coverage
- Unique alleles in ITSNdb subset: 35
- Unsupported by PRIME 2.1: none
- n_peptides skipped due to unsupported HLA: 0

## Scoring
- For each peptide, PRIME ran with `-a <singleAllele>` so that `%Rank_bestAllele` = that peptide's allele rank.
- score_prime = -%Rank (lower %Rank = more immunogenic, so we flip sign so HIGHER = more immunogenic).

## Runtime
- 49.0s wall-clock for 319 (peptide, HLA) pairs across 35 alleles.

## Results
- **ITSNdb_combined** (n=319, n_pos=136): AUROC = 0.567 [0.501, 0.628]
- **ITSNdb_no_overlap** (n=106, n_pos=33): AUROC = 0.572 [0.452, 0.686]
- **ITSNdb_in_master** (n=213, n_pos=103): AUROC = 0.564 [0.487, 0.636]
