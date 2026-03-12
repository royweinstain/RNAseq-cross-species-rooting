# Methods Log — RNA-seq Cross-Species Rooting Analysis

## Purpose
This log documents all data sources, software versions, processing parameters, and commands
used in the re-alignment of argan RNA-seq reads and re-analysis of cross-species orthogroups.
It serves as the basis for a publication Methods section.

---

## 1. Data Sources

### 1.1 Raw RNA-seq Reads
- **Project:** PRJNA863910
- **Publication:** Frontiers in Plant Science, 2022 — Argan adventitious rooting transcriptomics
- **Samples:** 30 total (2 clones x 2 tissues x 3 timepoints x 3 biological replicates)
  - Clones: ARS1 (= YM3, easy-to-root) and ARS7 (= 196, hard-to-root)
  - Tissues: Rooting zone (R) and Grafting zone (G)
  - Timepoints: T0, 24h, 120h
- **Download date:** [pending]
- **Tool:** sra-tools fasterq-dump

### 1.2 Reference Genomes

#### Argan (Sideroxylon spinosum) — Chromosome-level Assembly (2025)
- **Source:** ENA PRJEB88017 (haplotig 1)
- **Assembly:** Chromosome-level, 11 T2T chromosomes
- **BUSCO completeness:** 97.8%
- **Protein-coding genes:** ~28,720 per haplotype
- **Download date:** [pending]
- **Files:** genome FASTA, GFF3 annotation, protein FASTA

#### Poplar (Populus trichocarpa)
- **Source:** Phytozome / Ensembl Plants
- **Assembly version:** v3.x (to match Potri. gene IDs)
- **Proteins:** ~42,000
- **Purpose:** OrthoFinder proteome only (expression data uses existing datasets)
- **Download date:** [pending]

#### Arabidopsis (Arabidopsis thaliana)
- **Source:** Ensembl Plants / TAIR10
- **Proteins:** ~27,000
- **Purpose:** OrthoFinder anchor species only (no expression data)
- **Download date:** [pending]

---

## 2. Software Versions

| Tool | Version | Purpose |
|------|---------|---------|
| Miniconda | 3 (latest) | Package manager |
| HISAT2 | 2.2.2 | Splice-aware read alignment |
| samtools | 1.23 | BAM processing |
| subread (featureCounts) | 2.1.1 | Read counting per gene |
| OrthoFinder | (installed via bioconda) | Orthogroup inference |
| sra-tools | 3.2.1 | SRA download |
| FastQC | 0.12.1 | Read quality assessment |
| DESeq2 | [pending] | Differential expression (counts) |
| R | [pending] | Statistical computing |

---

## 3. Processing Steps

### 3.1 WSL + Conda Setup
- **Date:** [pending]
- **Commands:**
```bash
# [will be filled as commands are run]
```

### 3.2 Genome Downloads
- **Date:** [pending]
- **Commands:**
```bash
# [will be filled as commands are run]
```

### 3.3 SRA Downloads
- **Date:** 2026-03-11 to 2026-03-12
- **SRA accessions (T0 samples only):**
| SRA Run | Sample Name | Clone | Timepoint | Replicate |
|---------|-------------|-------|-----------|-----------|
| SRR22045442 | YM3-T0-1 | ARS1 (easy) | T0 | 1 |
| SRR22045441 | YM3-T0-2 | ARS1 (easy) | T0 | 2 |
| SRR22045440 | YM3-T0-3 | ARS1 (easy) | T0 | 3 |
| SRR22045436 | 196-T0-1 | ARS7 (hard) | T0 | 1 |
| SRR22045437 | 196-T0-2 | ARS7 (hard) | T0 | 2 |
| SRR22045438 | 196-T0-3 | ARS7 (hard) | T0 | 3 |

### 3.4 HISAT2 Indexing
- **Date:** [pending]
- **Commands:**
```bash
# [will be filled as commands are run]
```

### 3.5 Read Alignment
- **Date:** 2026-03-12
- **Parameters:** `hisat2 -x INDEX -1 R1 -2 R2 --dta -p 8 | samtools sort -@ 4 -o BAM`
- **Mapping rates:**
| Sample | Total read pairs | Overall alignment rate |
|--------|-----------------|----------------------|
| YM3-T0-1 | 8,811,230 | 86.84% |
| YM3-T0-2 | 14,789,072 | 86.41% |
| YM3-T0-3 | 14,813,474 | 85.83% |
| 196-T0-1 | 18,205,400 | 88.65% |
| 196-T0-2 | 21,085,600 | 88.39% |
| 196-T0-3 | 19,155,914 | 87.91% |

### 3.6 featureCounts
- **Date:** 2026-03-12
- **Parameters:** `featureCounts -a GTF -o OUT -T 8 -p --countReadPairs -t exon -g gene_id`
- **Annotation:** 164,386 features (exons), 24,721 meta-features (genes), 13 chromosomes/contigs
- **Assignment rates:**
| Sample | Total alignments | Assigned | Rate |
|--------|-----------------|----------|------|
| YM3-T0-1 | 15,001,531 | 9,421,870 | 62.8% |
| YM3-T0-2 | 25,866,876 | 15,353,754 | 59.4% |
| YM3-T0-3 | 25,366,841 | 12,638,163 | 49.8% |
| 196-T0-1 | 32,231,124 | 9,421,564 | 29.2% |
| 196-T0-2 | 30,423,545 | 12,960,514 | 42.6% |
| 196-T0-3 | 21,772,154 | 13,925,751 | 64.0% |
- **Genes with non-zero counts:** 23,480 / 24,721

### 3.7 OrthoFinder
- **Date:** 2026-03-11
- **Input proteomes:**
  - argan: Sspinosum_hap1 proteins (from chromosome-level genome)
  - poplar: Ptrichocarpa v3.x proteins (Ensembl Plants)
  - arabidopsis: TAIR10 proteins (Ensembl Plants)
- **Output:** 20,157 orthogroups total
  - With argan+poplar: 13,941
  - With all 3 species: 12,470

### 3.8 Pipeline Re-run
- **Date:** 2026-03-12
- **Config changes:**
  - `orthofinder_custom_file`: `Orthogroups.csv` → `Orthogroups_v2.tsv`
  - Argan `expression_file`: old XLOC_ FPKM Excel → new counts CSV from featureCounts
  - Argan `data_format`: `fpkm` (limma-voom) → `counts` (DESeq2)
  - Argan `gene_id_prefixes`: `["XLOC_"]` → `["g"]`
  - `orthofinder_column` values updated to match new OrthoFinder output
- **DE analysis:** DESeq2 for argan (20,707 genes after filtering), limma-voom for both poplar datasets
- **Results:**
  - Argan OG mapping rate: **92.9%** (was 19.6%)
  - Poplar OP42_vs_T89 OG mapping: 65.7%
  - Poplar G_vs_B OG mapping: 68.1%
  - OGs with data in ≥2 datasets: **218** (was 87)
  - Conserved OGs (score ≥ 0.8): **92** (was 38)
    - up_in_hard: 43 | up_in_easy: 49
  - Annotated conserved OGs (with Arabidopsis orthologs): **81**

---

## 4. Quality Metrics Summary

### Previous Run (2026-03-09, XLOC_ IDs)
- Argan OG mapping rate: 19.6%
- Poplar OG mapping rate: 57.3%
- Orthogroups with DE in both species: 87
- Conserved OGs (score >= 0.8): 38

### Current Run (2026-03-12, chromosome-level genome)
- Argan OG mapping rate: 92.9%
- Poplar OP42_vs_T89 OG mapping rate: 65.7%
- Poplar G_vs_B OG mapping rate: 68.1%
- Orthogroups with data in ≥2 datasets: 218
- Conserved OGs (score ≥ 0.8): 92 (43 up_in_hard, 49 up_in_easy)

---

## 5. File Checksums
```
# [will be filled with md5sum outputs]
```
