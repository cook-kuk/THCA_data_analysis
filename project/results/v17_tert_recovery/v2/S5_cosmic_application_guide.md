# COSMIC v100 manual application guide

COSMIC bulk access requires a free academic license. To recover sample-level
TERT promoter mutations for thyroid:

1. Register at https://cancer.sanger.ac.uk/cosmic/register (academic email)
2. Download `CosmicMutantExport.tsv.gz` and `CosmicSample.tsv.gz` from
   `https://cancer.sanger.ac.uk/cosmic/download` (v100, GRCh38)
3. After download, run:

```bash
zcat CosmicMutantExport.tsv.gz | awk -F'\t' '$1=="TERT"' > TERT_all.tsv
# Filter to thyroid using sample mapping
zcat CosmicSample.tsv.gz | awk -F'\t' '$5 ~ /thyroid/i' > thyroid_samples.tsv
```

4. Join `TERT_all.tsv` to `thyroid_samples.tsv` on sample_id.
5. Filter rows where `Mutation CDS == c.-124C>T` or `c.-146C>T` (= C228T/C250T promoter hotspots).

Once the file is on disk, point this script at it via env var
`COSMIC_TSV_PATH=/path/to/TERT_all.tsv` and rerun.
