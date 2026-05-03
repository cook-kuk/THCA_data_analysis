#!/bin/bash
# Pod D — K2 STAR re-quantification (D7-P3 closure)
# Run inside RunPod CPU pod after SSH connect
#
# Strategy: 5M-read partial download from ENA + STAR + featureCounts
#           on K2 PRJEB11591 cohort to enable raw-score 3-cohort meta.
# Expected runtime: ~30 hr on 16-core CPU (or 8-12 hr on 32-core)
# Expected cost: $5-10
set -e

WORK=/workspace/k2_star
mkdir -p $WORK $WORK/fastq $WORK/star_index $WORK/aligned $WORK/results
cd $WORK

# ============================================================
# Phase 1 — Install tools
# ============================================================
echo "=== [Phase 1] Tool install ==="
DEBIAN_FRONTEND=noninteractive apt-get update -qq
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    curl wget pigz parallel git unzip \
    samtools bedtools bowtie2 \
    python3-pip 2>&1 | tail -3

# STAR aligner
wget -q https://github.com/alexdobin/STAR/raw/master/bin/Linux_x86_64_static/STAR -O /usr/local/bin/STAR
chmod +x /usr/local/bin/STAR
STAR --version

# Subread (featureCounts)
wget -q https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-Linux-x86_64.tar.gz/download -O /tmp/subread.tgz
tar xzf /tmp/subread.tgz -C /tmp/
cp /tmp/subread-2.0.6-Linux-x86_64/bin/featureCounts /usr/local/bin/
featureCounts -v 2>&1 | head -1

pip install -q pandas numpy

# ============================================================
# Phase 2 — Reference download + STAR index build
# ============================================================
echo "=== [Phase 2] Reference setup ==="
cd $WORK/star_index

# GENCODE v44 primary assembly + basic annotation
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz
pigz -d *.gz

# STAR index — sparse (sjdb 100, sparse 1)
echo "  Building STAR genome index (~30-45 min)..."
STAR --runMode genomeGenerate \
     --genomeDir $WORK/star_index \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v44.basic.annotation.gtf \
     --sjdbOverhang 99 \
     --runThreadN 16 \
     --genomeSAsparseD 2 \
     --limitGenomeGenerateRAM 60000000000 \
     2>&1 | tail -5
echo "  STAR index built"

# ============================================================
# Phase 3 — K2 ENA partial download (PoC: 30 samples)
# ============================================================
echo "=== [Phase 3] K2 ENA partial download (30 samples PoC) ==="
cd $WORK/fastq

# Get K2 SRR list (first 30 samples for PoC)
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB11591&result=read_run&fields=run_accession,sample_alias,fastq_ftp,read_count&format=tsv" -o /tmp/k2_runs.tsv
head -31 /tmp/k2_runs.tsv > /tmp/k2_runs_30.tsv
wc -l /tmp/k2_runs_30.tsv

# Parallel download with 5M-read crop
fetch_one() {
    SRR=$1
    URL1=$(echo "$2" | tr ';' '\n' | head -1)
    URL2=$(echo "$2" | tr ';' '\n' | tail -1)
    cd $WORK/fastq
    if [ -f "${SRR}_1.fq.gz" ] && [ -f "${SRR}_2.fq.gz" ]; then
        echo "$(date +%H:%M:%S) skip $SRR (already done)"
        return
    fi
    echo "$(date +%H:%M:%S) start $SRR"
    curl -s "ftp://$URL1" | zcat 2>/dev/null | head -20000000 | pigz -p 2 > "${SRR}_1.fq.gz"
    curl -s "ftp://$URL2" | zcat 2>/dev/null | head -20000000 | pigz -p 2 > "${SRR}_2.fq.gz"
    echo "$(date +%H:%M:%S) done $SRR"
}
export -f fetch_one
export WORK

tail -n +2 /tmp/k2_runs_30.tsv | awk -F'\t' '{print $1, $3}' | parallel -j 8 --colsep ' ' fetch_one {1} {2}

ls *.fq.gz | wc -l

# ============================================================
# Phase 4 — STAR alignment
# ============================================================
echo "=== [Phase 4] STAR alignment ==="
cd $WORK/aligned

align_one() {
    SRR=$1
    cd $WORK/aligned
    if [ -f "${SRR}/Aligned.sortedByCoord.out.bam" ]; then return; fi
    mkdir -p $SRR
    STAR --runMode alignReads \
         --genomeDir $WORK/star_index \
         --readFilesIn $WORK/fastq/${SRR}_1.fq.gz $WORK/fastq/${SRR}_2.fq.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix ${SRR}/ \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN 4 \
         --outSAMunmapped Within 2>&1 | tail -3
}
export -f align_one
export WORK

tail -n +2 /tmp/k2_runs_30.tsv | awk '{print $1}' | parallel -j 4 align_one {}

ls $WORK/aligned/*/Aligned.sortedByCoord.out.bam | wc -l

# ============================================================
# Phase 5 — featureCounts
# ============================================================
echo "=== [Phase 5] featureCounts ==="
cd $WORK/results

featureCounts -a $WORK/star_index/gencode.v44.basic.annotation.gtf \
              -o $WORK/results/k2_counts.tsv \
              -p --countReadPairs -T 16 \
              $WORK/aligned/*/Aligned.sortedByCoord.out.bam \
              2>&1 | tail -10

# Extract 8-gene counts
GENE_8="SLC5A5\|TPO\|TG\|TSHR\|PAX8\|NKX2-1\|FOXE1\|DIO1"
grep -E "^Geneid|gene_name=($GENE_8)" $WORK/results/k2_counts.tsv > $WORK/results/k2_8gene_counts.tsv

ls -la $WORK/results/

echo "=== ALL DONE ==="
df -h $WORK
