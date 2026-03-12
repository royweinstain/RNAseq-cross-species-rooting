#!/bin/bash
# Download Malus domestica GDDH13 v1.1 reference genome from Ensembl Plants (release 62)
# Assembly: ASM211411v1 (GCA_002114115.1)
# Doubled haploid of Golden Delicious, cultivar GDDH13 (X9273)
# Daccord et al. (2017) Nature Genetics 49:1099-1106
set -eo pipefail

PIPELINE="/mnt/c/Users/royweinstain/Desktop/Claude/RNAseq hard and easy/pipeline"
OUTDIR="$PIPELINE/reference_genomes/apple"
mkdir -p "$OUTDIR"

BASE_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62"

echo "=== Downloading Malus domestica GDDH13 v1.1 genome ==="

# Genome FASTA
echo "Downloading genome FASTA..."
wget -O "$OUTDIR/Malus_domestica_golden.ASM211411v1.dna.toplevel.fa.gz" \
    "$BASE_URL/fasta/malus_domestica_golden/dna/Malus_domestica_golden.ASM211411v1.dna.toplevel.fa.gz"

# GFF3 annotation
echo "Downloading GFF3 annotation..."
wget -O "$OUTDIR/Malus_domestica_golden.ASM211411v1.62.gff3.gz" \
    "$BASE_URL/gff3/malus_domestica_golden/Malus_domestica_golden.ASM211411v1.62.gff3.gz"

# Protein FASTA (for OrthoFinder)
echo "Downloading protein FASTA..."
wget -O "$OUTDIR/Malus_domestica_golden.ASM211411v1.pep.all.fa.gz" \
    "$BASE_URL/fasta/malus_domestica_golden/pep/Malus_domestica_golden.ASM211411v1.pep.all.fa.gz"

# Decompress genome and annotation for HISAT2/featureCounts
echo "Decompressing..."
gunzip -k "$OUTDIR/Malus_domestica_golden.ASM211411v1.dna.toplevel.fa.gz"
gunzip -k "$OUTDIR/Malus_domestica_golden.ASM211411v1.62.gff3.gz"
gunzip -k "$OUTDIR/Malus_domestica_golden.ASM211411v1.pep.all.fa.gz"

# Convert GFF3 to GTF for featureCounts (gffread from cufflinks or conda)
echo "Converting GFF3 to GTF..."
gffread "$OUTDIR/Malus_domestica_golden.ASM211411v1.62.gff3" \
    -T -o "$OUTDIR/Malus_domestica_golden.ASM211411v1.62.gtf"

echo ""
echo "=== Download complete ==="
echo "Genome: $OUTDIR/Malus_domestica_golden.ASM211411v1.dna.toplevel.fa"
echo "GFF3:   $OUTDIR/Malus_domestica_golden.ASM211411v1.62.gff3"
echo "GTF:    $OUTDIR/Malus_domestica_golden.ASM211411v1.62.gtf"
echo "Protein: $OUTDIR/Malus_domestica_golden.ASM211411v1.pep.all.fa"
echo "Proteins: $(grep -c '^>' "$OUTDIR/Malus_domestica_golden.ASM211411v1.pep.all.fa")"
