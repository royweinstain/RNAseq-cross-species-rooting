#!/bin/bash
set -eo pipefail
export HOME=/home/royweinstain
export LD_LIBRARY_PATH=""
source /home/royweinstain/miniconda3/etc/profile.d/conda.sh
conda activate bioinfo

PIPELINE="/home/royweinstain/pipeline"
GTF="$PIPELINE/reference_genomes/argan/Sspinosum_hap1.gtf"
BAM_DIR="$PIPELINE/datasets/argan_196_vs_YM3/raw/bam"
OUT_DIR="$PIPELINE/datasets/argan_196_vs_YM3/processed"
mkdir -p "$OUT_DIR"

echo "=== Running featureCounts ==="
echo "GTF: $GTF"
echo "BAM dir: $BAM_DIR"
echo "Start: $(date)"

featureCounts \
    -a "$GTF" \
    -o "$OUT_DIR/argan_counts_raw.txt" \
    -T 8 \
    -p --countReadPairs \
    -t exon -g gene_id \
    "$BAM_DIR/YM3-T0-1.bam" \
    "$BAM_DIR/YM3-T0-2.bam" \
    "$BAM_DIR/YM3-T0-3.bam" \
    "$BAM_DIR/196-T0-1.bam" \
    "$BAM_DIR/196-T0-2.bam" \
    "$BAM_DIR/196-T0-3.bam"

echo ""
echo "=== featureCounts complete ==="
echo "End: $(date)"
echo ""
echo "Output:"
ls -lh "$OUT_DIR/argan_counts_raw.txt"
echo ""
echo "First few lines:"
head -5 "$OUT_DIR/argan_counts_raw.txt"
echo ""
echo "Gene count:"
tail -n +3 "$OUT_DIR/argan_counts_raw.txt" | wc -l
