# Install R packages required for step 01 (differential expression)
# Run once: Rscript install_r_packages.R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c("DESeq2", "edgeR", "limma"), ask = FALSE)

install.packages(c("yaml", "readxl"), repos = "https://cloud.r-project.org")
