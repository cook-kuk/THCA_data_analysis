# Public Training Corpora Drop Directory

Place row-level public-tool training corpora here as TSV files.

Required minimum columns:
- `peptide`
- `hla` or `hla_allele_4digit` when available
- `public_tool`
- `source_dataset` or `source_name` when available
- `train_split`, `label`, `assay_type`, and `publication_or_url` when available

After adding files, rerun:

```bash
python scripts/clean_neobench/run_clean_neobench_pipeline.py \
  --repo-root . \
  --output-root project/results/clean_neobench_barneo_2026_05_09
```

Public pretrained tools remain caveated until row-level exact peptide, exact peptide-HLA, and near-peptide overlap audit passes.
