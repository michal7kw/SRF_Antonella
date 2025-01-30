# Load required libraries
suppressPackageStartupMessages({
    library(DiffBind)
    library(tidyverse)
    library(rtracklayer)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(GenomicRanges)
    library(clusterProfiler)
})

# Utility functions
log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
}

validate_files <- function(files) {
    missing <- files[!file.exists(files)]
    if (length(missing) > 0) {
        stop("Missing files:\n", paste(missing, collapse="\n"))
    }
    TRUE
}

create_dirs <- function(dirs) {
    for (dir in dirs) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
}

# Function to standardize chromosome names in peak files
standardize_peak_files <- function(peak_files) {
    for (file in peak_files) {
        log_message(sprintf("Standardizing chromosome names in %s", file))
        
        # Read peak file with exact MACS2 peak format specifications
        peaks <- try({
            read.table(file, header=FALSE, 
                      colClasses=c("character", "numeric", "numeric", "character",
                                 "numeric", "character", "numeric", "numeric", 
                                 "numeric", "numeric"),
                      col.names=c("chr", "start", "end", "name", "score", 
                                "strand", "signalValue", "pValue", "qValue", "peak"))
        })
        
        # If 10 columns fail (narrowPeak), try 9 columns (broadPeak)
        if (inherits(peaks, "try-error")) {
            peaks <- read.table(file, header=FALSE,
                              colClasses=c("character", "numeric", "numeric", "character",
                                         "numeric", "character", "numeric", "numeric", "numeric"),
                              col.names=c("chr", "start", "end", "name", "score",
                                        "strand", "signalValue", "pValue", "qValue"))
        }
        
        # Add 'chr' prefix if missing
        if (!all(grepl("^chr", peaks$chr))) {
            peaks$chr <- paste0("chr", peaks$chr)
            # Handle special case where "chr" might have been added to "chrM"
            peaks$chr <- gsub("^chrchr", "chr", peaks$chr)
        }
        
        # Keep only standard chromosomes
        standard_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
        peaks <- peaks[peaks$chr %in% standard_chroms,]
        
        # Write back
        write.table(peaks, file, sep="\t", quote=FALSE, 
                   row.names=FALSE, col.names=FALSE)
        
        log_message(sprintf("Processed %d peaks in %s", nrow(peaks), file))
    }
}

#' Perform differential binding analysis
#' @param peak_type Character, either "broad" or "narrow"
perform_diffbind <- function(peak_type = c("broad", "narrow")) {
    peak_type <- match.arg(peak_type)
    log_message(sprintf("Starting differential binding analysis for %s peaks", peak_type))
    
    # Create output directories
    dirs <- c(
        file.path("analysis", paste0("diffbind_", peak_type)),
        file.path("analysis", paste0("annotation_", peak_type)),
        "logs"
    )
    create_dirs(dirs)
    
    # Create sample sheet
    samples <- data.frame(
        SampleID = c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)),
        Factor = rep("H2AK119Ub", 6),
        Condition = rep(c("YAF", "GFP"), each=3),
        Replicate = rep(1:3, 2),
        bamReads = file.path("analysis/aligned", 
                            c(paste0("YAF_", 1:3, ".dedup.bam"),
                              paste0("GFP_", 1:3, ".dedup.bam"))),
        Peaks = file.path("analysis/peaks",
                         paste0(c(paste0("YAF_", 1:3), paste0("GFP_", 1:3)),
                               "_", peak_type, "_peaks.", 
                               ifelse(peak_type == "broad", "broadPeak", "narrowPeak"))),
        PeakCaller = rep(peak_type, 6)
    )
    
    # Validate input files
    validate_files(c(samples$bamReads, samples$Peaks))
    
    # Standardize chromosome names in peak files
    standardize_peak_files(samples$Peaks)
    
    # Create DiffBind object with specific configuration
    log_message("Creating DiffBind object...")
    dba_data <- dba(sampleSheet=samples,
                    minOverlap=1,  # Reduced from 2 to be less stringent
                    peakFormat="bed",
                    peakCaller=peak_type,
                    config=data.frame(AnalysisMethod="max",
                                    fragmentSize=150,
                                    doBlacklist=TRUE,
                                    RunParallel=FALSE))
    
    # Print diagnostic information
    log_message(sprintf("Number of binding sites before counting: %d", 
                       length(dba.show(dba_data)$Intervals)))
    
    # Count reads with less stringent parameters
    log_message("Counting reads...")
    dba_data <- dba.count(dba_data, 
                         bUseSummarizeOverlaps=TRUE,
                         minCount=1,       # Minimum read count
                         score=DBA_SCORE_READS)  # Use raw read counts
    
    # Print diagnostic information after counting
    count_info <- dba.show(dba_data)
    log_message(sprintf("Number of binding sites after counting: %d", 
                       length(count_info$Intervals)))
    log_message("Read counts per sample:")
    tryCatch({
        print(count_info[,c("Reads.counted", "FRiP")])
    }, error = function(e) {
        log_message("Unable to print detailed count information", "WARNING")
        print(count_info)  # Print full count info instead
    })
    
    # Normalize and perform differential analysis with adjusted parameters
    log_message("Performing differential analysis...")
    dba_data <- dba.normalize(dba_data, 
                             normalize=DBA_NORM_TMM,
                             background=FALSE)  # Don't subtract background
    
    dba_data <- dba.contrast(dba_data, 
                            categories=DBA_CONDITION,
                            minMembers=2)
    
    dba_data <- dba.analyze(dba_data, 
                           method=DBA_DESEQ2,
                           bSubControl=FALSE,  # Don't subtract control reads
                           bFullLibrarySize=TRUE)  # Use full library size
    
    # Extract and save results with less stringent thresholds
    log_message("Extracting and saving results...")
    dba_results <- dba.report(dba_data,
                             th=1,            # No p-value threshold
                             bCalled=FALSE,   # Include all sites
                             bNormalized=TRUE)
    
    # Save complete analysis
    saveRDS(dba_data, file.path("analysis", paste0("diffbind_", peak_type), "all_peaks.rds"))
    saveRDS(dba_results, file.path("analysis", paste0("diffbind_", peak_type), "significant_peaks.rds"))
    
    # Convert to data frame and save
    df_results <- as.data.frame(dba_results)
    write.csv(df_results, 
             file.path("analysis", paste0("diffbind_", peak_type), "differential_peaks.csv"), 
             row.names=FALSE)
    
    log_message(sprintf("Found %d differential binding sites", nrow(df_results)))
    
    # Annotate peaks
    log_message("Annotating peaks...")
    peaks_gr <- dba_results
    seqlevelsStyle(peaks_gr) <- "UCSC"  # Ensure UCSC-style chromosome names
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peakAnno <- annotatePeak(peaks_gr, TxDb=txdb,
                            tssRegion=c(-3000, 3000),
                            verbose=FALSE)
    
    # Save annotation results
    log_message("Saving annotation results...")
    saveRDS(peakAnno, file.path("analysis", paste0("annotation_", peak_type), "peak_annotation.rds"))
    write.csv(as.data.frame(peakAnno),
             file.path("analysis", paste0("annotation_", peak_type), "peak_annotation.csv"),
             row.names=FALSE)
    
    # Generate annotation plots
    log_message("Generating annotation plots...")
    pdf(file.path("analysis", paste0("annotation_", peak_type), "annotation_plots.pdf"))
    plotAnnoPie(peakAnno)
    plotDistToTSS(peakAnno)
    dev.off()
    
    # # Perform GO enrichment analysis
    # log_message("Performing GO enrichment analysis...")
    # genes <- unique(peakAnno@anno$geneId)
    # ego <- enrichGO(gene = genes,
    #                OrgDb = org.Hs.eg.db,
    #                keyType = "ENTREZID",
    #                ont = "BP",
    #                pAdjustMethod = "BH",
    #                pvalueCutoff = 0.05,
    #                qvalueCutoff = 0.2)
    
    # # Save GO analysis results
    # saveRDS(ego, file.path("analysis", paste0("annotation_", peak_type), "go_enrichment.rds"))
    # write.csv(as.data.frame(ego),
    #          file.path("analysis", paste0("annotation_", peak_type), "go_enrichment.csv"),
    #          row.names=FALSE)
    
    gc()
    
    return(dba_results)
}

# Run analysis for both peak types
# broad_results <- perform_diffbind("broad")
narrow_results <- perform_diffbind("narrow")