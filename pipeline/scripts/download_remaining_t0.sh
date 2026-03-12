#!/bin/bash
set -eo pipefail
export HOME=/home/royweinstain
export LD_LIBRARY_PATH=""
source /home/royweinstain/miniconda3/etc/profile.d/conda.sh
conda activate bioinfo

OUTDIR="/mnt/c/Users/royweinstain/Desktop/Claude/RNAseq hard and easy/pipeline/datasets/argan_196_vs_YM3/raw/fastq"

download() {
    local SRR=$1 NAME=$2
    if [ -f "$OUTDIR/${NAME}_1.fastq.gz" ] && [ -f "$OUTDIR/${NAME}_2.fastq.gz" ]; then
        echo "SKIP $NAME"; return 0
    fi
    echo "=== Downloading $SRR -> $NAME === $(date)"
    fasterq-dump "$SRR" --outdir "$OUTDIR" --split-files --threads 4 --progress
    echo "Compressing..."
    gzip -c "$OUTDIR/${SRR}_1.fastq" > "$OUTDIR/${NAME}_1.fastq.gz"
    gzip -c "$OUTDIR/${SRR}_2.fastq" > "$OUTDIR/${NAME}_2.fastq.gz"
    rm "$OUTDIR/${SRR}_1.fastq" "$OUTDIR/${SRR}_2.fastq"
    echo "DONE $NAME $(date)"
}

download SRR22045441 YM3-T0-2
download SRR22045440 YM3-T0-3
download SRR22045436 196-T0-1
download SRR22045437 196-T0-2
download SRR22045438 196-T0-3

echo "=== ALL T0 DOWNLOADS COMPLETE ==="
ls -lh "$OUTDIR/"
