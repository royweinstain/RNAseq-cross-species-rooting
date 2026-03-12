#!/usr/bin/env python3
"""Convert featureCounts output to pipeline-ready CSV."""
import csv
import sys
import os

infile = sys.argv[1]
outfile = sys.argv[2] if len(sys.argv) > 2 else infile.replace("_raw.txt", ".csv")

# featureCounts output: skip first line (command), then header + data
# Columns: Geneid, Chr, Start, End, Strand, Length, then BAM columns
with open(infile) as f:
    lines = f.readlines()

# Skip the command line (starts with #)
data_lines = [l for l in lines if not l.startswith('#')]
reader = csv.reader(data_lines, delimiter='\t')
header = next(reader)

# Extract sample names from BAM paths (last column names)
# featureCounts uses full path as column name, extract just the sample name
sample_cols = []
for i, col in enumerate(header[6:], start=6):
    # Extract sample name from path like /path/to/YM3-T0-1.bam
    name = os.path.basename(col).replace('.bam', '')
    sample_cols.append((i, name))

print(f"Found {len(sample_cols)} samples: {[s[1] for s in sample_cols]}")

# Write CSV with gene_id + sample columns
with open(outfile, 'w', newline='') as out:
    writer = csv.writer(out)
    # Header
    writer.writerow(['gene_id'] + [s[1] for s in sample_cols])

    count = 0
    for row in reader:
        gene_id = row[0]
        counts = [row[i] for i, _ in sample_cols]
        # Skip genes with all zeros
        if all(c == '0' for c in counts):
            continue
        writer.writerow([gene_id] + counts)
        count += 1

print(f"Wrote {count} genes (non-zero) to {outfile}")
