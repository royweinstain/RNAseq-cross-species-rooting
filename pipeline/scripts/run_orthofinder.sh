#!/bin/bash
# Run OrthoFinder on 3 proteomes (argan, poplar, arabidopsis)
set -eo pipefail
export HOME=/home/royweinstain
export LD_LIBRARY_PATH=""
source /home/royweinstain/miniconda3/etc/profile.d/conda.sh
conda activate bioinfo

PIPELINE="/mnt/c/Users/royweinstain/Desktop/Claude/RNAseq hard and easy/pipeline"
PROTEOMES="$PIPELINE/orthogroups/proteomes"

echo "=== Running OrthoFinder ==="
echo "Input: $PROTEOMES"
echo "Proteomes:"
for f in "$PROTEOMES"/*.fa; do
    echo "  $(basename $f): $(grep -c '^>' $f) proteins"
done
echo "Start time: $(date)"

orthofinder -f "$PROTEOMES" -t 8 -a 4

echo ""
echo "=== OrthoFinder complete ==="
echo "End time: $(date)"

# Find and copy results
RESULTS_DIR=$(ls -dt "$PROTEOMES/OrthoFinder/Results_"* 2>/dev/null | head -1)
if [ -n "$RESULTS_DIR" ]; then
    echo "Results: $RESULTS_DIR"
    echo "Orthogroups file:"
    ls -lh "$RESULTS_DIR/Orthogroups/Orthogroups.tsv"
    echo "Stats:"
    cat "$RESULTS_DIR/Comparative_Genomics_Statistics/Statistics_Overall.tsv" 2>/dev/null || true
    # Copy to pipeline location
    cp "$RESULTS_DIR/Orthogroups/Orthogroups.tsv" "$PIPELINE/orthogroups/Orthogroups_v2.tsv"
    echo "Copied to: $PIPELINE/orthogroups/Orthogroups_v2.tsv"
fi
