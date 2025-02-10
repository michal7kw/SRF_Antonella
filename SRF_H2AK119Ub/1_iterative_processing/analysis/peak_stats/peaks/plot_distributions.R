library(ggplot2)
library(gridExtra)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
data_dir <- args[2]

# Read and process length distribution data
length_dist <- read.table(file.path(data_dir, paste0(sample, "_length_dist_sorted.txt")))
colnames(length_dist) <- c("Length", "Count")

# Read and process length vs fold enrichment data
length_fold <- read.table(file.path(data_dir, paste0(sample, "_length_vs_fold.txt")))
colnames(length_fold) <- c("Length", "FoldEnrichment")

# Create length distribution plot
p1 <- ggplot(length_dist, aes(x=Length, y=Count)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    labs(title=paste("Peak Length Distribution -", sample),
         x="Peak Length (bp)", y="Number of Peaks") +
    scale_x_continuous(breaks=seq(0, max(length_dist$Length), by=2000))

# Create length vs fold enrichment scatter plot
p2 <- ggplot(length_fold, aes(x=Length, y=FoldEnrichment)) +
    geom_point(alpha=0.3) +
    theme_minimal() +
    labs(title=paste("Peak Length vs Fold Enrichment -", sample),
         x="Peak Length (bp)", y="Fold Enrichment") +
    scale_x_continuous(breaks=seq(0, max(length_fold$Length), by=2000))

# Combine plots into single PDF
pdf(file.path(data_dir, paste0(sample, "_peak_analysis.pdf")), width=12, height=8)
grid.arrange(p1, p2, ncol=2)
dev.off()

# Calculate and save length quantiles
quantiles <- quantile(length_fold$Length, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
write.table(quantiles, file=file.path(data_dir, paste0(sample, "_length_quantiles.txt")),
            quote=FALSE, col.names=FALSE)
