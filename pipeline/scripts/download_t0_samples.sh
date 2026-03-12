#!/bin/bash
# Download T0 (ctrl) samples from SRA for argan cross-species comparison
set -eo pipefail
export HOME=/home/royweinstain
export LD_LIBRARY_PATH=""
source /home/royweinstain/miniconda3/etc/profile.d/conda.sh
conda activate bioinfo

OUTDIR="/mnt/c/Users/royweinstain/Desktop/Claude/RNAseq hard and easy/pipeline/datasets/argan_196_vs_YM3/raw/fastq"
mkdir -p "$OUTDIR"

# T0 samples: ctrl = before any treatment
declare -A SAMPLES=(
    ["SRR22045439"]="YM3-T0-1"
    ["SRR22045441"]="YM3-T0-2"
    ["SRR22045440"]="YM3-T0-3"
    ["SRR22045436"]="196-T0-1"
    ["SRR22045437"]="196-T0-2"
    ["SRR22045438"]="196-T0-3"
)

echo "=== Downloading T0 samples from PRJNA863910 ==="
echo "Output: $OUTDIR"
echo "Start time: $(date)"

for SRR in "${!SAMPLES[@]}"; do
    NAME="${SAMPLES[$SRR]}"

    if [ -f "$OUTDIR/${NAME}_1.fastq.gz" ] && [ -f "$OUTDIR/${NAME}_2.fastq.gz" ]; then
        echo "SKIP $NAME (already exists)"
        continue
    fi

    echo ""
    echo "--- Downloading $SRR -> $NAME ---"
    fasterq-dump "$SRR" --outdir "$OUTDIR" --split-files --threads 4 --progress

    echo "Compressing $NAME ..."
    gzip -c "$OUTDIR/${SRR}_1.fastq" > "$OUTDIR/${NAME}_1.fastq.gz"
    gzip -c "$OUTDIR/${SRR}_2.fastq" > "$OUTDIR/${NAME}_2.fastq.gz"
    rm "$OUTDIR/${SRR}_1.fastq" "$OUTDIR/${SRR}_2.fastq"

    echo "DONE $NAME"
    ls -lh "$OUTDIR/${NAME}"_*.fastq.gz
done

echo ""
echo "=== All T0 downloads complete ==="
echo "End time: $(date)"
ls -lh "$OUTDIR/"
