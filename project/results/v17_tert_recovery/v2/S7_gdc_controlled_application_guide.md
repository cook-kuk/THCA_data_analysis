# dbGaP / GDC controlled-access application guide for TCGA-THCA TERT recovery

The TCGA-THCA WES MAFs in open access **do not capture TERT promoter** because
the promoter falls outside the WES exonic baits. To recover TERT promoter
status for TCGA-THCA, you need WGS BAMs (controlled access).

## Steps

1. **dbGaP application** — study accession `phs000178` (TCGA), data use limitation
   "General Research Use" or "Health/Medical/Biomedical".
2. Get an eRA Commons account if you don't have one (PI or sponsoring PI).
3. Submit a Data Access Request (DAR) on dbGaP web.
4. After approval (typically 2–6 weeks):
   - Download GDC token from the GDC portal.
   - Use `gdc-client` to download THCA WGS BAMs:
     ```bash
     gdc-client download -t gdc-user-token.txt -m thca_wgs_manifest.tsv
     ```
   - Generate a manifest from `parsed/S7_gdc_controlled_thca_seq_files.tsv` filtering to
     `experimental_strategy == "WGS"` and `data_format == "BAM"`.
5. **Variant calling at TERT promoter region**:
   - Region: chr5:1,294,500-1,295,800 (GRCh38) or chr5:1,295,000-1,296,000 (GRCh37).
   - Use `mutect2` or `varscan2` with low-VAF settings (TERT promoter is often
     subclonal, ~5-15% VAF).
   - Recommended: `gatk Mutect2 -R Homo_sapiens_assembly38.fasta -I tumor.bam -L TERT_promoter.bed`

## Realistic timeline
- DAR submission to data download: 4–8 weeks.
- Variant calling: 1–2 days for ~50 BAMs on a modest cluster.

This source is **flagged for future work**, not actionable in the current sprint.
