# Apple Juvenile vs Adult Dataset (apple_J_vs_A)

## Source Publication
Li X, Shen F, Xu X, Zheng Q, Wang Y, Wu T, Li W, Qiu C, Xu X, Han Z, Zhang X. 2021.
An HD-ZIP transcription factor, MxHB13, integrates auxin-regulated and juvenility-determined
control of adventitious rooting in *Malus xiaojinensis*.
*The Plant Journal* 107(6): 1663-1680. DOI: 10.1111/tpj.15406

## BioProject
PRJNA392909 (NCBI SRA)

## Species
*Malus xiaojinensis* (syn. *Malus baccata* var. *xiaojinensis*) — wild apple rootstock

## Experimental Design
- **Juvenile (Mx-J):** Leafy stem cuttings from basal suckers (easy-to-root, ~77% rooting rate)
- **Adult (Mx-A):** Leafy stem cuttings from canopy of reproductively mature trees (hard-to-root, ~11%)
- **Timepoint:** T0 (before any auxin treatment)
- **Replicates:** 3 biological replicates per condition (50 cuttings pooled per replicate)
- **Sequencing:** Illumina HiSeq X Ten, paired-end, ~11-15M read pairs per sample

## Samples Used (T0 only — 6 of 54 total in BioProject)
| SRA Run | Sample Name | Phase | Rooting |
|---------|-------------|-------|---------|
| SRR15322003 | J_0_rep1 | Juvenile | Easy |
| SRR15322002 | J_0_rep2 | Juvenile | Easy |
| SRR15322030 | J_0_rep3 | Juvenile | Easy |
| SRR15494548 | A_0_rep1 | Adult | Hard |
| SRR15494559 | A_0_rep2 | Adult | Hard |
| SRR15494546 | A_0_rep3 | Adult | Hard |

## Reference Genome
*Malus domestica* GDDH13 v1.1 (doubled haploid of Golden Delicious)
- Assembly: ASM211411v1 (GCA_002114115.1)
- Source: Ensembl Plants release 62
- Daccord et al. (2017) *Nature Genetics* 49: 1099-1106
- ~40,624 protein-coding genes; 17 chromosomes
- Gene ID format: MD01G0000000

Note: The original paper used the older Velasco et al. (2010) apple genome (MDP0000... IDs).
We use the newer GDDH13 assembly for better gene models and consistency with current databases.

## Processing Pipeline
1. Download raw reads: `scripts/apple_download_samples.sh`
2. Download genome: `scripts/apple_download_genome.sh`
3. Build HISAT2 index: `scripts/apple_build_index.sh`
4. Align + count: `scripts/apple_align_and_count.sh`
   - HISAT2 v2.2.2 (--dta, 8 threads) → samtools v1.23 sort
   - featureCounts (subread v2.1.1): paired-end, exon-level, gene_id aggregation
5. DE analysis: DESeq2 via pipeline step 1 (juvenile vs adult, positive log2FC = higher in adult/hard)

## Biological Context
This dataset compares developmental stages (juvenile vs adult) rather than genotypes. The juvenile
phase has high miR156 expression and adventitious rooting competence. The adult phase has high
SPL26 expression, which represses MxHB13, reducing auxin transport and rooting capacity. This is
biologically related but distinct from the genotype-based easy-vs-hard comparisons in the argan
and poplar datasets.

## Full BioProject
The full PRJNA392909 contains 54 samples: 2 phases (J/A) × 2 treatments (IBA/water) ×
5 timepoints (0, 6, 12, 24, 48h) × 3 replicates. The IBA-treated and later timepoint samples
are available for future analysis of auxin-responsive rooting genes.
