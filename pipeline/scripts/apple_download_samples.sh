#!/bin/bash
# Download 6 T0 RNA-seq samples for apple juvenile vs adult comparison
# BioProject: PRJNA392909
# Li et al. (2021) The Plant Journal 107(6):1663-1680
# DOI: 10.1111/tpj.15406
set -eo pipefail
export HOME=/home/royweinstain
export LD_LIBRARY_PATH=""
source /home/royweinstain/miniconda3/etc/profile.d/conda.sh
conda activate bioinfo

PIPELINE="/mnt/c/Users/royweinstain/Desktop/Claude/RNAseq hard and easy/pipeline"
OUTDIR="$PIPELINE/datasets/apple_J_vs_A/raw/fastq"
mkdir -p "$OUTDIR"

echo "=== Downloading Apple T0 samples (6 samples) ==="
echo "Juvenile (Mx-J, easy-to-root): basal sucker cuttings"
echo "Adult (Mx-A, hard-to-root): canopy shoot cuttings"
echo ""

# Juvenile T0 (easy-to-root)
declare -A SAMPLES=(
    ["J_0_rep1"]="SRR15322003"
    ["J_0_rep2"]="SRR15322002"
    ["J_0_rep3"]="SRR15322030"
    ["A_0_rep1"]="SRR15494548"
    ["A_0_rep2"]="SRR15494559"
    ["A_0_rep3"]="SRR15494546"
)

for SAMPLE in "${!SAMPLES[@]}"; do
    SRR="${SAMPLES[$SAMPLE]}"
    echo "--- Downloading $SAMPLE ($SRR) ---"
    if [ -f "$OUTDIR/${SAMPLE}_1.fastq.gz" ] && [ -f "$OUTDIR/${SAMPLE}_2.fastq.gz" ]; then
        echo "  Already exists, skipping."
        continue
    fi
    fasterq-dump "$SRR" -O "$OUTDIR" -e 4 --split-files
    # Rename to sample names
    mv "$OUTDIR/${SRR}_1.fastq" "$OUTDIR/${SAMPLE}_1.fastq"
    mv "$OUTDIR/${SRR}_2.fastq" "$OUTDIR/${SAMPLE}_2.fastq"
    # Compress
    gzip "$OUTDIR/${SAMPLE}_1.fastq" &
    gzip "$OUTDIR/${SAMPLE}_2.fastq" &
    wait
    echo "  Done: ${SAMPLE}_1.fastq.gz, ${SAMPLE}_2.fastq.gz"
done

echo ""
echo "=== All downloads complete ==="
ls -lh "$OUTDIR"/*.fastq.gz
