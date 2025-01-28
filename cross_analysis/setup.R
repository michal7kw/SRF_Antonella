#!/usr/bin/env Rscript

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required packages
BiocManager::install(c(
    "GenomicRanges",
    "ChIPseeker",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db",
    "ComplexHeatmap",
    "circlize",
    "GenomeInfoDb",
    "rtracklayer",
    "RColorBrewer"
), force = TRUE)

# Install CRAN packages
install.packages(c(
    "ggplot2"
))
