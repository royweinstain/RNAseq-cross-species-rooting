# RNA-seq Cross-Species Rooting Analysis Pipeline

A modular bioinformatics pipeline to identify conserved gene expression patterns distinguishing easy-to-root from hard-to-root tree genotypes across species.

## Overview

This pipeline integrates RNA-seq differential expression data from three independent datasets (one argan, two poplar) and uses OrthoFinder-based ortholog inference to find conserved rooting-associated genes. Arabidopsis is included as an anchor species for functional annotation only.

## Datasets

| Dataset | Species | Comparison | Format | Samples | Source |
|---------|---------|-----------|--------|---------|--------|
| argan_196_vs_YM3 | *Argania spinosa* | ARS7 (hard) vs ARS1 (easy) | Raw counts | 3 vs 3 | Tzeela et al. (2022), BioProject PRJNA863910 |
| poplar_OP42_vs_T89 | *Populus* spp. | T89 (hard) vs OP42 (easy) | Log2 FPKM | 3 vs 3 | Ranjan et al. (2022) |
| poplar_G_vs_B | *Populus* F1 progeny | B (hard) vs G (easy) | FPKM | 9 vs 9 | Sun et al. (2019) |

## Pipeline Steps

1. **Differential expression** (`steps/01_de_analysis.R`) — DESeq2 for counts, limma/limma-voom for FPKM
2. **Orthogroup parsing** (`steps/02_parse_orthogroups.py`) — Map genes to OrthoFinder orthogroups
3. **Cross-species conservation** (`steps/03_cross_species.py`) — Majority-vote scoring across datasets
4. **Reporting** (`steps/04_report.py`) — Volcano plots, heatmaps, summary statistics
5. **Functional annotation** (`steps/05_annotate.py`) — Arabidopsis GO terms via MyGene.info

## Usage

```bash
cd pipeline/
python run_pipeline.py config.yaml              # full run
python run_pipeline.py config.yaml --dry-run     # preview only
python run_pipeline.py config.yaml --steps 2,3,4 # skip R step
```

## Reproducing the Analysis

All expression data, DE results, and orthogroup tables needed to reproduce steps 2–5 are included in this repository. See **[REPRODUCING.md](REPRODUCING.md)** for detailed instructions.

**Quick start:**
```bash
pip install -r requirements.txt
Rscript install_r_packages.R
cd pipeline/
python run_pipeline.py config.yaml --steps 2,3,4,5
```

### Data not included (too large for git)

- **Raw reads**: NCBI SRA BioProject PRJNA863910 (use `scripts/download_t0_samples.sh`)
- **Reference genome**: *S. spinosum* chromosome-level assembly, ENA PRJEB88017
- **Proteomes**: Ensembl Plants (Arabidopsis TAIR10, *P. trichocarpa* v3.x) and ENA (*S. spinosum*)

These are only needed if re-running read alignment or OrthoFinder from scratch.

## WSL Scripts (`pipeline/scripts/`)

These bash scripts were used in WSL (Windows Subsystem for Linux) with a conda `bioinfo` environment:

| Script | Purpose |
|--------|---------|
| `download_t0_samples.sh` | Download argan T0 SRA reads |
| `build_hisat2_index.sh` | Build HISAT2 genome index |
| `align_sample.sh` | Align reads (HISAT2 + samtools sort) |
| `run_featurecounts.sh` | Gene-level read counting |
| `run_orthofinder.sh` | OrthoFinder 3-species run |
| `format_counts.py` | Format count matrix for pipeline |

## Key Results (2026-03-12 run)

- 92 conserved orthogroups (43 up in hard-to-root, 49 up in easy-to-root)
- All 92 with perfect cross-dataset agreement (conservation score = 1.0)
- 81 annotated with Arabidopsis GO terms

## References

- Tzeela et al. (2022) *Front. Plant Sci.* 13: 1002703
- Ranjan et al. (2022) *J. Exp. Bot.* 73(12): 4046–4064
- Sun et al. (2019) *Int. J. Mol. Sci.* 20(24): 6114
- Mateus et al. (2025) *Sci. Data* 12: 1430

## License

This project is provided for academic and research purposes.
