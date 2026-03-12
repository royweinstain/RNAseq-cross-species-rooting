# Methods

## Plant material and RNA-seq datasets

Three independent transcriptomic datasets were used to identify conserved gene expression patterns associated with adventitious rooting capacity across tree species.

**Argan dataset (argan_196_vs_YM3).** RNA-seq data for *Argania spinosa* (syn. *Sideroxylon spinosum*) were obtained from the BioProject PRJNA863910 (Tzeela et al., 2022). Two clonally propagated genotypes differing in rooting ability — ARS1 (clone YM3; easy-to-root) and ARS7 (clone 196; hard-to-root) — were compared at the T0 timepoint (prior to auxin treatment), with three biological replicates per genotype. Raw paired-end reads were downloaded from NCBI SRA using sra-tools v3.2.1 (fasterq-dump).

**Poplar OP42 vs T89 dataset (poplar_OP42_vs_T89).** Cambium tissue-specific transcriptome data for hybrid poplar OP42 (easy-to-root) and hybrid aspen T89 (hard-to-root) were obtained from supplementary tables of Ranjan et al. (2022). Log2-transformed FPKM expression values at the T0 timepoint were used, with three biological replicates per genotype.

**Poplar G vs B dataset (poplar_G_vs_B).** Transcriptome data from an F1 progeny population derived from a cross between *Populus deltoides* 'Danhong' and *P. simonii* 'Tongliao1' were obtained from supplementary tables of Sun et al. (2019). From 434 phenotyped genotypes, three fine-rooting (G1–G3; easy-to-root) and three poor-rooting (B1–B3; hard-to-root) genotypes were selected for RNA-seq, with three biological replicates each (9 vs 9 samples). The dataset comprised FPKM values for 10,172 pre-filtered differentially expressed genes.

## Read alignment and quantification

For the argan dataset, raw paired-end reads from six T0 samples were aligned to the chromosome-level phased assembly of *Sideroxylon spinosum* haplotig 1 (Mateus et al., 2025; ENA accession PRJEB88017). This assembly comprises 11 telomere-to-telomere chromosomes with BUSCO completeness exceeding 97.8% and approximately 28,720 predicted protein-coding genes. Reads were aligned using HISAT2 v2.2.2 (Kim et al., 2019) with the `--dta` flag for downstream transcript assembly compatibility, using 8 threads. Alignments were piped directly to samtools v1.23 (Danecek et al., 2021) for coordinate sorting. Overall alignment rates ranged from 85.8% to 88.7% across the six samples.

Gene-level read counts were generated using featureCounts from the Subread package v2.1.1 (Liao et al., 2014) in paired-end mode (`-p --countReadPairs`), counting at the exon level with gene-level aggregation (`-t exon -g gene_id`). The reference annotation contained 164,386 exon features across 24,721 genes on 13 chromosomes and contigs. Assignment rates ranged from 29.2% to 64.0%, and 23,480 genes had non-zero counts across the six samples.

The two poplar datasets were not re-processed from raw reads; expression values were used as provided in the original publications.

## Differential expression analysis

Differential expression (DE) was assessed independently for each dataset using R (R Core Team, 2024). In all analyses, the contrast was defined as hard-to-root versus easy-to-root, such that positive log2 fold-change (log2FC) values indicate higher expression in hard-to-root genotypes.

For the argan dataset (raw counts), DE analysis was performed using DESeq2 (Love et al., 2014). Genes with fewer than 10 counts in fewer than 3 samples were excluded prior to analysis. The remaining 20,707 genes were tested using the default median-of-ratios normalization and Wald test, with Benjamini–Hochberg (BH) adjustment for multiple testing.

For the poplar OP42 vs T89 dataset (log2 FPKM values), the limma package (Ritchie et al., 2015) was applied without voom transformation, as the data were already log2-transformed. A total of 41,335 genes were tested.

For the poplar G vs B dataset (FPKM values), the limma-voom pipeline was applied with TMM normalization via edgeR (Robinson et al., 2010). A total of 10,172 genes were tested.

Genes were classified as significantly differentially expressed if the BH-adjusted p-value was less than 0.05 and the absolute log2FC was at least 1.0. Under these criteria, 1,983, 2,356, and 648 significant DE genes were identified in the argan, poplar OP42 vs T89, and poplar G vs B datasets, respectively.

## Ortholog inference

Orthologous gene groups were inferred using OrthoFinder v3.1.3 (Emms and Kelly, 2019) with DIAMOND for sequence similarity search (BLOSUM62 scoring matrix), FAMSA for multiple sequence alignment, and FastTree for gene tree inference. Input proteomes comprised *Arabidopsis thaliana* TAIR10 (~27,000 proteins; Lamesch et al., 2012), *Sideroxylon spinosum* haplotig 1 (~28,720 proteins; Mateus et al., 2025), and *Populus trichocarpa* v3.x (~42,000 proteins; Goodstein et al., 2012). *Arabidopsis thaliana* was included as an anchor species for functional annotation only and did not contribute expression data.

OrthoFinder identified 20,157 orthogroups in total, of which 13,941 contained representatives from both argan and poplar, and 12,470 contained representatives from all three species.

## Cross-species conservation scoring

Differentially expressed genes from each dataset were mapped to orthogroups using normalized gene identifiers with transcript and isoform suffixes stripped. Orthogroup mapping rates were 92.9% for argan, 65.7% for poplar OP42 vs T89, and 68.1% for poplar G vs B.

For each orthogroup in each dataset, the direction of differential expression was classified by majority vote among the significant genes within that orthogroup (up_in_hard if the majority had positive log2FC, up_in_easy if negative). A conservation score was then computed across datasets as the fraction of datasets agreeing on the modal direction. An orthogroup was classified as conserved if the conservation score was at least 0.80 and expression data were present in at least two of the three datasets.

Of 218 orthogroups with significant DE data in at least two datasets, 92 were classified as conserved (43 up_in_hard, 49 up_in_easy). All 92 conserved orthogroups had perfect cross-dataset agreement (conservation score = 1.0).

## Functional annotation

Arabidopsis orthologs from the 92 conserved orthogroups were annotated via the MyGene.info v3 REST API (Xin et al., 2016). Gene symbol, gene name, Gene Ontology biological process, and GO molecular function terms were retrieved. Of the 92 conserved orthogroups, 81 had Arabidopsis orthologs with functional annotations available.

## Visualization

Volcano plots, conservation heatmaps, and direction agreement charts were generated using matplotlib (Hunter, 2007) and seaborn (Waskom, 2021) in Python.

---

# References

Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. 2021. Twelve years of SAMtools and BCFtools. *GigaScience* 10(2): giab008.

Emms DM, Kelly S. 2019. OrthoFinder: phylogenetic orthology inference for comparative genomics. *Genome Biology* 20: 238.

Goodstein DM, Shu S, Howson R, Neupane R, Hayes RD, Fazo J, Mitros T, Dirks W, Hellsten U, Putnam N, Rokhsar DS. 2012. Phytozome: a comparative platform for green plant genomics. *Nucleic Acids Research* 40(D1): D1178–D1186.

Hunter JD. 2007. Matplotlib: a 2D graphics environment. *Computing in Science & Engineering* 9(3): 90–95.

Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. 2019. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology* 37(8): 907–915.

Lamesch P, Berardini TZ, Li D, Swarbreck D, Wilks C, Sasidharan R, Muller R, Dreher K, Alexander DL, Garcia-Hernandez M, Karthikeyan AS, Lee CH, Nelson WD, Ploetz L, Singh S, Wensel A, Huala E. 2012. The Arabidopsis Information Resource (TAIR): improved gene annotation and new tools. *Nucleic Acids Research* 40(D1): D1202–D1210.

Liao Y, Smyth GK, Shi W. 2014. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics* 30(7): 923–930.

Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15(12): 550.

Mateus CS, Oliveira T, Nunes J, Wallberg A, Sousa P, Dias D. 2025. Chromosome-level phased genome assembly of *Sideroxylon spinosum* (Sapotaceae). *Scientific Data* 12: 1430.

R Core Team. 2024. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria.

Ranjan A, Persson S, Gorzás A, Bhalerao RP. 2022. Cell-type-specific gene expression profiling reveals roles for vascular cambium cell programs in secondary growth. *Journal of Experimental Botany* 73(12): 4046–4064.

Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research* 43(7): e47.

Robinson MD, McCarthy DJ, Smyth GK. 2010. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics* 26(1): 139–140.

Sun P, Jia H, Zhang Y, Li J, Lu M, Hu J. 2019. Deciphering genetic architecture of adventitious root and related shoot traits in poplar using QTL mapping and RNA-seq data. *International Journal of Molecular Sciences* 20(24): 6114.

Tzeela S, Guillen-Palma A, Pery S, Cohen H, Tel-Zur N, Fait A, Galili S, Levy A, Hochberg U, David-Schwartz R, Sadot E. 2022. Transcriptomics of morphological and biochemical changes in shoot cuttings of argan tree during adventitious root induction. *Frontiers in Plant Science* 13: 1002703.

Waskom ML. 2021. seaborn: statistical data visualization. *Journal of Open Source Software* 6(60): 3021.

Xin J, Mark A, Afrasiabi C, Tsueng G, Jber M, Bharber K, Head S, Lu D, Shafir A, Shepherd O, Grossmann R, Wu C, Su AI. 2016. High-performance web services for querying gene and variant annotation. *Genome Biology* 17: 91.
