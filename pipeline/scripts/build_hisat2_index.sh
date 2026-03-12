#!/bin/bash
# Build HISAT2 index for argan chromosome-level genome
set -eo pipefail
export HOME=/home/royweinstain
export LD_LIBRARY_PATH=""
source /home/royweinstain/miniconda3/etc/profile.d/conda.sh
conda activate bioinfo

PIPELINE="/mnt/c/Users/royweinstain/Desktop/Claude/RNAseq hard and easy/pipeline"
GENOME="$PIPELINE/reference_genomes/argan/Sspinosum_hap1_genome.fa"
INDEX_DIR="$PIPELINE/reference_genomes/argan/hisat2_index"
INDEX_PREFIX="$INDEX_DIR/Sspinosum_hap1"

mkdir -p "$INDEX_DIR"

echo "=== Building HISAT2 index ==="
echo "Genome: $GENOME"
echo "Index prefix: $INDEX_PREFIX"
echo "Start time: $(date)"

hisat2-build -p 8 "$GENOME" "$INDEX_PREFIX"

echo "=== Index build complete ==="
echo "End time: $(date)"
ls -lh "$INDEX_DIR/"
