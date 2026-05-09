# MHCflurry Training Caveat - 2026-05-09

## Verdict

Treat MHCflurry as a **public-pretrained upper-bound comparator**, not a clean external baseline.

It is not proven to be "cheating" on ITSNdb, but we also cannot exclude training-set overlap from the files available locally. The no-overlap flag used in the wave analyses only means no overlap with **our in-house master table**, not no overlap with MHCflurry's own training data.

## Local provenance

Observed local install:

- Package imported in repo venv: `mhcflurry 2.2.1`.
- Loaded presentation predictor supports 14,883 alleles.
- `provenance_string`: `generated on Thu Jun 11 13:37:18 2020`.
- Downloaded model: `models_class1_presentation.20200611.tar.bz2`.
- Model `info.txt`: trained on `Thu Jun 11 13:37:18 2020` with package `mhcflurry 1.7.0`.

Local files:

- `/home/seungho/.local/share/mhcflurry/4/2.2.0/models_class1_presentation/DOWNLOAD_INFO.csv`
- `/home/seungho/.local/share/mhcflurry/4/2.2.0/models_class1_presentation/models/info.txt`
- `/home/seungho/.local/share/mhcflurry/4/2.2.0/models_class1_presentation/models/affinity_predictor/info.txt`
- `/home/seungho/.local/share/mhcflurry/4/2.2.0/models_class1_presentation/models/processing_predictor_without_flanks/info.txt`

## Source-level read

Official MHCflurry documentation states pretrained models are fit to mass-spec-identified MHC-I ligands and peptide/MHC affinity measurements deposited in IEDB plus a few other sources, and that the presentation predictor integrates binding-affinity and antigen-processing predictors.

Implication: if ITSNdb peptides or closely related measurements appear in IEDB/MS-derived training data before June 2020, MHCflurry's `ITSNdb_no_overlap` score may be partially contaminated.

## Curation policy

Use this language:

- OK: "MHCflurry public-pretrained comparator"
- OK: "MHCflurry upper-bound comparator with unresolved training overlap"
- OK: "MHCflurry has no inflation against our master-table overlap flag"

Avoid this language:

- "clean external MHCflurry baseline"
- "true no-overlap MHCflurry performance"
- "MHCflurry proves field ceiling"

## Practical consequence

For paper claims, the safest local anchor is `Structure_LR`:

- It has known local training/evaluation provenance.
- It remains competitive on no-overlap.
- In the strict Wave 8 artifact audit after removing TCR/self exact hits, `Structure_LR` is the cleaner leader.

MHCflurry should remain in the figure, but visually annotated as public-pretrained / training-corpus-unresolved.
