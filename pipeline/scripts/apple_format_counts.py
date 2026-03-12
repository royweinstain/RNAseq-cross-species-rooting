#!/usr/bin/env python3
"""Format featureCounts output into a clean counts matrix for DESeq2.

Usage: python apple_format_counts.py apple_counts_raw.txt apple_counts.csv
"""

import sys
import pandas as pd

if len(sys.argv) < 3:
    print("Usage: python apple_format_counts.py <input_raw.txt> <output.csv>")
    sys.exit(1)

raw_file = sys.argv[1]
out_file = sys.argv[2]

# Read featureCounts output (skip first comment line)
df = pd.read_csv(raw_file, sep="\t", comment="#")

# First column is Geneid, columns 6+ are counts (columns 1-5 are Chr, Start, End, Strand, Length)
count_cols = df.columns[6:]

# Build clean dataframe
counts = df[["Geneid"]].copy()
counts.rename(columns={"Geneid": "gene_id"}, inplace=True)

# Rename BAM path columns to sample names
for col in count_cols:
    # Extract sample name from BAM path (e.g., .../J_0_rep1.bam -> J_0_rep1)
    sample_name = col.split("/")[-1].replace(".bam", "")
    counts[sample_name] = df[col]

counts.to_csv(out_file, index=False)
print(f"Formatted counts: {len(counts)} genes x {len(count_cols)} samples -> {out_file}")
