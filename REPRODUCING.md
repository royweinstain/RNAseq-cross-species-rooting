# Reproducing the Analysis

This guide explains how to reproduce the full cross-species rooting analysis from this repository.

## Quick Start (steps 2–5 only)

If you only want to reproduce the cross-species analysis starting from pre-computed DE results (included in this repo), you can skip read alignment and differential expression:

```bash
# 1. Install Python dependencies
pip install -r requirements.txt

# 2. Install R packages
Rscript install_r_packages.R

# 3. Run the pipeline from step 2 onward
cd pipeline/
python run_pipeline.py config.yaml --steps 2,3,4,5
```

This uses the DE results already in `pipeline/datasets/*/processed/de_results.csv` and the orthogroup table in `pipeline/orthogroups/Orthogroups_v2.tsv`.

## Full Reproduction (all steps)

### Prerequisites

- Python 3.10+
- R 4.x with Bioconductor
- WSL (Windows) or Linux for read alignment steps
- Conda (for bioinformatics tools)

### Step 0: Install dependencies

```bash
# Python
pip install -r requirements.txt

# R
Rscript install_r_packages.R

# WSL/Linux bioinformatics tools (only for re-aligning argan reads)
conda env create -f environment_wsl.yml
conda activate bioinfo
```

### Step 1: Obtain expression data

The poplar expression data files are included in this repo:
- `pipeline/datasets/poplar_OP42_vs_T89/raw/20230807 poplar.xlsx` — from Ranjan et al. (2022) supplementary
- `pipeline/datasets/poplar_G_vs_B/raw/20230823 poplus different lines.xlsx` — from Sun et al. (2019) supplementary

The argan counts matrix is also included:
- `pipeline/datasets/argan_196_vs_YM3/processed/argan_counts.csv`

**If you want to regenerate argan counts from raw reads** (optional):

1. Download raw reads from NCBI SRA (BioProject PRJNA863910):
   - See `pipeline/datasets/argan_196_vs_YM3/sra_samples.tsv` for accessions
   - Use `pipeline/scripts/download_t0_samples.sh`

2. Download the *S. spinosum* reference genome (Mateus et al., 2025):
   - ENA accession PRJEB88017, haplotig 1
   - Place genome FASTA and GFF3 in `pipeline/reference_genomes/argan/`

3. Build index and align (in WSL):
   ```bash
   conda activate bioinfo
   bash pipeline/scripts/build_hisat2_index.sh
   bash pipeline/scripts/align_sample.sh   # for each sample
   bash pipeline/scripts/run_featurecounts.sh
   python pipeline/scripts/format_counts.py
   ```

### Step 2: Run differential expression (step 1)

```bash
cd pipeline/
python run_pipeline.py config.yaml --steps 1
```

This runs DESeq2 on argan counts and limma/limma-voom on the poplar datasets, producing `de_results.csv` in each dataset's `processed/` folder.

### Step 3: Run the rest of the pipeline (steps 2–5)

```bash
python run_pipeline.py config.yaml --steps 2,3,4,5
```

This will:
- Parse orthogroups and map genes (step 2)
- Compute cross-species conservation scores (step 3)
- Generate figures and run report (step 4)
- Annotate conserved orthogroups via MyGene.info (step 5)

Output goes to `pipeline/output/<timestamp>/`.

### Step 4: Re-running OrthoFinder (optional)

The orthogroup table (`Orthogroups_v2.tsv`) is included in this repo. To regenerate it:

1. Download proteomes:
   - *A. thaliana* TAIR10 proteins from Ensembl Plants
   - *S. spinosum* hap1 proteins from ENA PRJEB88017
   - *P. trichocarpa* v3.x proteins from Ensembl Plants or Phytozome

2. Place `.fa` files in `pipeline/orthogroups/proteomes/`

3. Run in WSL:
   ```bash
   conda activate bioinfo
   bash pipeline/scripts/run_orthofinder.sh
   ```

## Expected Results

With the included data, steps 2–5 should produce:
- 218 orthogroups with DE data in >=2 datasets
- **92 conserved orthogroups** (43 up_in_hard, 49 up_in_easy)
- All 92 with conservation score = 1.0
- 81 annotated with Arabidopsis GO terms

## Configuration

All thresholds are in `pipeline/config.yaml`:
- `padj_threshold: 0.05`
- `log2fc_threshold: 1.0`
- `conservation_threshold: 0.80`
- `min_datasets_for_conservation: 2`
