#!/bin/bash
# Align 6 apple T0 samples with HISAT2 and count with featureCounts
set -eo pipefail
export HOME=/home/royweinstain
export LD_LIBRARY_PATH=""
source /home/royweinstain/miniconda3/etc/profile.d/conda.sh
conda activate bioinfo

PIPELINE="/mnt/c/Users/royweinstain/Desktop/Claude/RNAseq hard and easy/pipeline"
INDEX="$PIPELINE/reference_genomes/apple/hisat2_index/apple_GDDH13"
GTF="$PIPELINE/reference_genomes/apple/Malus_domestica_golden.ASM211411v1.62.gtf"
FASTQ_DIR="$PIPELINE/datasets/apple_J_vs_A/raw/fastq"
BAM_DIR="$PIPELINE/datasets/apple_J_vs_A/raw/bam"
PROCESSED="$PIPELINE/datasets/apple_J_vs_A/processed"

mkdir -p "$BAM_DIR" "$PROCESSED"

SAMPLES=("J_0_rep1" "J_0_rep2" "J_0_rep3" "A_0_rep1" "A_0_rep2" "A_0_rep3")

echo "=== Aligning apple samples with HISAT2 ==="
echo "Index: $INDEX"
echo "GTF: $GTF"
echo ""

for SAMPLE in "${SAMPLES[@]}"; do
    R1="$FASTQ_DIR/${SAMPLE}_1.fastq.gz"
    R2="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"
    BAM="$BAM_DIR/${SAMPLE}.bam"
    LOG="$BAM_DIR/${SAMPLE}.hisat2.log"

    if [ -f "$BAM" ]; then
        echo "--- $SAMPLE: BAM exists, skipping alignment ---"
        continue
    fi

    echo "--- Aligning $SAMPLE ---"
    echo "Start: $(date)"

    hisat2 -x "$INDEX" \
        -1 "$R1" -2 "$R2" \
        --dta -p 8 \
        2>"$LOG" \
        | samtools sort -@ 4 -o "$BAM"

    samtools index "$BAM"

    echo "  Alignment rate: $(grep 'overall alignment rate' "$LOG")"
    echo "  Done: $(date)"
done

echo ""
echo "=== Running featureCounts ==="
echo "Start: $(date)"

BAM_FILES=""
for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILES="$BAM_FILES $BAM_DIR/${SAMPLE}.bam"
done

featureCounts \
    -a "$GTF" \
    -o "$PROCESSED/apple_counts_raw.txt" \
    -T 8 \
    -p --countReadPairs \
    -t exon \
    -g gene_id \
    $BAM_FILES

echo ""
echo "=== featureCounts complete ==="
echo "End: $(date)"

# Format counts matrix
echo ""
echo "=== Formatting counts matrix ==="
python3 "$PIPELINE/scripts/apple_format_counts.py" \
    "$PROCESSED/apple_counts_raw.txt" \
    "$PROCESSED/apple_counts.csv"

echo "Done. Output: $PROCESSED/apple_counts.csv"
