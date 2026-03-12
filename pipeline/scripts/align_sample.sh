#!/bin/bash
set -eo pipefail
export HOME=/home/royweinstain
export LD_LIBRARY_PATH=""
source /home/royweinstain/miniconda3/etc/profile.d/conda.sh
conda activate bioinfo

PIPELINE="/home/royweinstain/pipeline"
INDEX="$PIPELINE/reference_genomes/argan/hisat2_index/Sspinosum_hap1"
FASTQ_DIR="$PIPELINE/datasets/argan_196_vs_YM3/raw/fastq"
BAM_DIR="$PIPELINE/datasets/argan_196_vs_YM3/raw/bam"
mkdir -p "$BAM_DIR"

NAME="$1"
R1="$FASTQ_DIR/${NAME}_1.fastq.gz"
R2="$FASTQ_DIR/${NAME}_2.fastq.gz"
BAM="$BAM_DIR/${NAME}.bam"

if [ -f "$BAM" ] && [ -f "${BAM}.bai" ]; then
    echo "SKIP $NAME (BAM already exists)"
    exit 0
fi

echo "=== Aligning $NAME ==="
echo "R1: $R1"
echo "R2: $R2"
echo "Start: $(date)"

hisat2 -x "$INDEX" \
       -1 "$R1" -2 "$R2" \
       --dta -p 8 \
       --summary-file "$BAM_DIR/${NAME}.hisat2.log" \
    | samtools sort -@ 4 -o "$BAM" -

samtools index "$BAM"

echo "=== Done $NAME ==="
echo "End: $(date)"
cat "$BAM_DIR/${NAME}.hisat2.log"
ls -lh "$BAM"
