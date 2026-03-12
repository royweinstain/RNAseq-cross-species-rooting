#!/bin/bash
# Build HISAT2 index for Malus domestica GDDH13 v1.1
set -eo pipefail
export HOME=/home/royweinstain
export LD_LIBRARY_PATH=""
source /home/royweinstain/miniconda3/etc/profile.d/conda.sh
conda activate bioinfo

PIPELINE="/mnt/c/Users/royweinstain/Desktop/Claude/RNAseq hard and easy/pipeline"
GENOME="$PIPELINE/reference_genomes/apple/Malus_domestica_golden.ASM211411v1.dna.toplevel.fa"
INDEX_DIR="$PIPELINE/reference_genomes/apple/hisat2_index"
INDEX_PREFIX="$INDEX_DIR/apple_GDDH13"

mkdir -p "$INDEX_DIR"

echo "=== Building HISAT2 index for apple GDDH13 v1.1 ==="
echo "Genome: $GENOME"
echo "Index prefix: $INDEX_PREFIX"
echo "Start: $(date)"

hisat2-build -p 8 "$GENOME" "$INDEX_PREFIX"

echo "End: $(date)"
echo "=== Index build complete ==="
ls -lh "$INDEX_DIR"/
