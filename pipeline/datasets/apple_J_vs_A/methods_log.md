# Methods Log — Apple Juvenile vs Adult Dataset

## 1. Data Source

### 1.1 RNA-seq Data
- **Publication:** Li et al. (2021) "An HD-ZIP transcription factor, MxHB13, integrates auxin-regulated and juvenility-determined control of adventitious rooting in *Malus xiaojinensis*." *The Plant Journal* 107(6): 1663-1680.
- **DOI:** 10.1111/tpj.15406
- **BioProject:** PRJNA392909
- **Species:** *Malus xiaojinensis* (syn. *Malus baccata* var. *xiaojinensis*)
- **Institution:** College of Horticulture, China Agricultural University, Beijing, China
- **Total samples in BioProject:** 54 (2 phases × 2 treatments × 5 timepoints × 3 replicates)
- **Samples used here:** 6 (T0 only, no treatment)
- **Sequencing platform:** Illumina HiSeq X Ten, paired-end

### 1.2 Experimental Design (from original paper)
- **Juvenile phase (Mx-J):** Semi-lignified leafy cuttings (8-10 cm) from basal suckers. High miR156 expression. Rooting rate ~77% with IBA treatment.
- **Adult phase (Mx-A):** Semi-lignified leafy cuttings from canopy of reproductively mature trees. High SPL26 expression. Rooting rate ~11% with IBA treatment.
- **T0 samples:** Cuttings before any IBA treatment (baseline expression)
- **Biological replicates:** 3 per condition, each replicate pooled from ≥50 individual cuttings

### 1.3 SRA Accessions
| SRA Run | Sample Name | Phase | Description |
|---------|-------------|-------|-------------|
| SRR15322003 | J_0_rep1 | Juvenile | Mx-J, T0, replicate 1 |
| SRR15322002 | J_0_rep2 | Juvenile | Mx-J, T0, replicate 2 |
| SRR15322030 | J_0_rep3 | Juvenile | Mx-J, T0, replicate 3 |
| SRR15494548 | A_0_rep1 | Adult | Mx-A, T0, replicate 1 |
| SRR15494559 | A_0_rep2 | Adult | Mx-A, T0, replicate 2 |
| SRR15494546 | A_0_rep3 | Adult | Mx-A, T0, replicate 3 |

### 1.4 Reference Genome
- **Species:** *Malus domestica* (domesticated apple, closest high-quality genome to *M. xiaojinensis*)
- **Assembly:** GDDH13 v1.1 (ASM211411v1, GCA_002114115.1)
- **Source:** Ensembl Plants release 62
- **Publication:** Daccord N, Celton JM, Linsmith G, et al. (2017) "High-quality de novo assembly of the apple genome and methylome dynamics of early fruit development." *Nature Genetics* 49: 1099-1106.
- **Type:** Doubled haploid of Golden Delicious, cultivar GDDH13 (X9273)
- **Chromosomes:** 17
- **Protein-coding genes:** ~40,624
- **Gene ID format:** MD01G0000000 (chromosome-based)
- **Note:** The original paper (Li et al., 2021) used the older Velasco et al. (2010) genome with MDP0000... gene IDs. We use the newer GDDH13 assembly for improved gene models.

---

## 2. Software Versions

| Tool | Version | Purpose |
|------|---------|---------|
| sra-tools (fasterq-dump) | 3.2.1 | SRA download |
| HISAT2 | 2.2.2 | Splice-aware read alignment |
| samtools | 1.23 | BAM processing |
| subread (featureCounts) | 2.1.1 | Gene-level read counting |
| gffread | (conda) | GFF3 to GTF conversion |
| DESeq2 | (R/Bioconductor) | Differential expression |

---

## 3. Processing Steps

### 3.1 Genome Download
- **Script:** `pipeline/scripts/apple_download_genome.sh`
- **Files downloaded:**
  - Genome FASTA: `Malus_domestica_golden.ASM211411v1.dna.toplevel.fa`
  - GFF3 annotation: `Malus_domestica_golden.ASM211411v1.62.gff3`
  - Protein FASTA: `Malus_domestica_golden.ASM211411v1.pep.all.fa` (40,624 proteins)
  - GTF converted from GFF3 via gffread
- **Location:** `pipeline/reference_genomes/apple/`

### 3.2 SRA Download
- **Script:** `pipeline/scripts/apple_download_samples.sh`
- **Tool:** sra-tools v3.2.1, fasterq-dump
- **Output:** `pipeline/datasets/apple_J_vs_A/raw/fastq/` (12 files: 6 samples × 2 paired-end)

### 3.3 HISAT2 Index Build
- **Script:** `pipeline/scripts/apple_build_index.sh`
- **Command:** `hisat2-build -p 8 <genome.fa> <index_prefix>`
- **Output:** `pipeline/reference_genomes/apple/hisat2_index/apple_GDDH13.*`

### 3.4 Read Alignment
- **Script:** `pipeline/scripts/apple_align_and_count.sh`
- **Command:** `hisat2 -x INDEX -1 R1 -2 R2 --dta -p 8 | samtools sort -@ 4 -o BAM`
- **Output:** `pipeline/datasets/apple_J_vs_A/raw/bam/*.bam`
- **Alignment rates:** [to be filled after running]

### 3.5 featureCounts
- **Script:** (same as 3.4)
- **Command:** `featureCounts -a GTF -o OUT -T 8 -p --countReadPairs -t exon -g gene_id`
- **Output:** `pipeline/datasets/apple_J_vs_A/processed/apple_counts_raw.txt`
- **Formatted output:** `pipeline/datasets/apple_J_vs_A/processed/apple_counts.csv`
- **Annotation stats:** [to be filled after running]

### 3.6 Differential Expression
- **Method:** DESeq2 (raw counts, 3 vs 3)
- **Contrast:** Adult (hard-to-root) vs Juvenile (easy-to-root)
- **Convention:** Positive log2FC = higher expression in adult (hard-to-root)
- **Thresholds:** padj < 0.05, |log2FC| >= 1.0
- **Output:** `pipeline/datasets/apple_J_vs_A/processed/de_results.csv`
- **Results:** [to be filled after running]

---

## 4. Key Differences from Original Paper Analysis

| Aspect | Li et al. (2021) | Our analysis |
|--------|-----------------|--------------|
| Reference genome | Velasco et al. (2010), MDP IDs | GDDH13 v1.1, MD IDs |
| Transcript assembly | StringTie | Not used (gene-level counting) |
| Quantification | StringTie → prepDE.py (counts) | featureCounts (counts) |
| DE analysis | DESeq2 | DESeq2 (same) |
| Comparisons | All timepoints, ±IBA, J vs A | T0 only: J vs A |
| Goal | Identify MxHB13 regulatory module | Find conserved rooting genes across species |

---

## 5. Run Execution Order

```bash
# In WSL, with conda env bioinfo activated:
bash pipeline/scripts/apple_download_genome.sh    # Step 1: Get reference genome
bash pipeline/scripts/apple_download_samples.sh   # Step 2: Get raw reads
bash pipeline/scripts/apple_build_index.sh         # Step 3: Build HISAT2 index
bash pipeline/scripts/apple_align_and_count.sh     # Step 4: Align + count

# Back in Windows/Python:
# DE analysis will be run through the main pipeline when config is updated
```

---

## 6. Citations

- Li X et al. (2021) *The Plant Journal* 107(6): 1663-1680. DOI: 10.1111/tpj.15406
- Daccord N et al. (2017) *Nature Genetics* 49: 1099-1106. (GDDH13 genome)
- Kim D et al. (2019) *Nature Biotechnology* 37: 907-915. (HISAT2)
- Danecek P et al. (2021) *GigaScience* 10: giab008. (samtools)
- Liao Y et al. (2014) *Bioinformatics* 30: 923-930. (featureCounts)
- Love MI et al. (2014) *Genome Biology* 15: 550. (DESeq2)
