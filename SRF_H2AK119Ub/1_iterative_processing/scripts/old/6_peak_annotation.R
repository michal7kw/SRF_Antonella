# Load required libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(rtracklayer)
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(DiffBind)

# Function to convert DBA object to GRanges
dba_to_granges <- function(dba_obj) {
    tryCatch({
        # Print DiffBind object info
        print("DiffBind object info:")
        print(dba_obj)
        
        # First try to get all peaks without any filtering
        peaks_df <- dba.report(dba_obj, th=1, bUsePval=FALSE)
        
        if (is.null(peaks_df) || nrow(peaks_df) == 0) {
            # Try getting peaks directly from the DBA object
            peaks_df <- as.data.frame(dba_obj$binding)
            
            if (is.null(peaks_df) || nrow(peaks_df) == 0) {
                # Try getting peaks from the peaks slot
                peaks_df <- as.data.frame(dba_obj$peaks[[1]])
                
                if (is.null(peaks_df) || nrow(peaks_df) == 0) {
                    stop("No peaks found in DBA object after trying multiple methods")
                }
            }
        }
        
        # Convert peaks_df to data.frame if it's not already
        peaks_df <- as.data.frame(peaks_df)
        
        # Print debugging information
        print("Peaks data summary:")
        print(paste("Number of peaks:", nrow(peaks_df)))
        print("Column names:")
        print(colnames(peaks_df))
        print("First few rows:")
        print(head(peaks_df))
        
        # Map column names to standard format if needed
        col_mapping <- list(
            Chr = grep("^chr$|^seqnames$|^Chr$", colnames(peaks_df), value = TRUE)[1],
            Start = grep("^start$|^Start$", colnames(peaks_df), value = TRUE)[1],
            End = grep("^end$|^End$", colnames(peaks_df), value = TRUE)[1]
        )
        
        # Check if required columns exist
        if (is.na(col_mapping$Chr) || is.na(col_mapping$Start) || is.na(col_mapping$End)) {
            stop("Required columns (Chr/seqnames, Start/start, End/end) not found in peaks data")
        }
        
        # Create GRanges object with minimal required fields
        gr <- GRanges(
            seqnames = as.character(peaks_df[[col_mapping$Chr]]),
            ranges = IRanges(
                start = as.numeric(as.character(peaks_df[[col_mapping$Start]])), 
                end = as.numeric(as.character(peaks_df[[col_mapping$End]]))
            ),
            strand = "*"
        )
        
        # Add metadata columns if they exist
        metadata_cols <- c(
            "Fold" = grep("^Fold$|^fold$", colnames(peaks_df), value = TRUE)[1],
            "FDR" = grep("^FDR$|^fdr$", colnames(peaks_df), value = TRUE)[1],
            "p.value" = grep("^p.value$|^pvalue$|^PValue$", colnames(peaks_df), value = TRUE)[1],
            "Conc" = grep("^Conc$|^conc$", colnames(peaks_df), value = TRUE)[1]
        )
        
        for (col_name in names(metadata_cols)) {
            col <- metadata_cols[[col_name]]
            if (!is.na(col) && col %in% colnames(peaks_df)) {
                mcols(gr)[[col_name]] <- as.numeric(as.character(peaks_df[[col]]))
            }
        }
        
        # Print GRanges summary
        print("Created GRanges object:")
        print(gr)
        
        return(gr)
    }, error = function(e) {
        stop(paste("Error in dba_to_granges:", e$message))
    })
}

# Function to standardize chromosome names
standardize_chromosomes <- function(gr) {
    tryCatch({
        # Print initial state
        print("Initial seqlevels:")
        print(seqlevels(gr))
        
        # First, check if chromosomes have 'chr' prefix
        current_chroms <- seqlevels(gr)
        
        # Add 'chr' prefix if missing
        if (!all(grepl("^chr", current_chroms))) {
            new_levels <- current_chroms
            for (i in seq_along(current_chroms)) {
                if (!grepl("^chr", current_chroms[i])) {
                    if (current_chroms[i] == "MT") {
                        new_levels[i] <- "chrM"
                    } else if (current_chroms[i] == "23") {
                        new_levels[i] <- "chrX"
                    } else if (current_chroms[i] == "24") {
                        new_levels[i] <- "chrY"
                    } else {
                        new_levels[i] <- paste0("chr", current_chroms[i])
                    }
                }
            }
            
            # Create mapping
            names(new_levels) <- current_chroms
            
            # Print mapping for debugging
            print("Chromosome name mapping:")
            print(new_levels)
            
            # Rename the seqlevels
            seqlevels(gr) <- new_levels
        }
        
        # Define standard chromosomes
        standard_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
        
        # Print chromosome stats
        print("Chromosome statistics before filtering:")
        print(table(seqnames(gr)))
        
        # Keep only standard chromosomes
        gr <- gr[seqnames(gr) %in% standard_chroms]
        
        if (length(gr) == 0) {
            stop("No peaks remaining after chromosome filtering")
        }
        
        # Set standard chromosome levels
        seqlevels(gr) <- standard_chroms[standard_chroms %in% seqlevels(gr)]
        
        # Print final state
        print("Final seqlevels:")
        print(seqlevels(gr))
        print("Final chromosome statistics:")
        print(table(seqnames(gr)))
        
        # Ensure genome info is set
        genome(gr) <- "hg38"
        
        return(gr)
    }, error = function(e) {
        stop(paste("Error in standardize_chromosomes:", e$message))
    })
}

# Function to perform annotation analysis
perform_annotation <- function(peak_type) {
    tryCatch({
        # Create output directories
        base_dir <- file.path("analysis", paste0("annotation_", peak_type))
        dir.create(file.path(base_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(base_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
        
        # Load differential binding results
        print(paste("Loading", peak_type, "peak differential binding results..."))
        diffbind_dir <- file.path("analysis", paste0("diffbind_", peak_type))
        
        # Check if files exist
        all_peaks_file <- file.path(diffbind_dir, "all_peaks.rds")
        sig_peaks_file <- file.path(diffbind_dir, "significant_peaks.rds")
        
        if (!file.exists(all_peaks_file)) {
            stop(paste("File not found:", all_peaks_file))
        }
        if (!file.exists(sig_peaks_file)) {
            stop(paste("File not found:", sig_peaks_file))
        }
        
        all_peaks <- readRDS(all_peaks_file)
        sig_peaks <- readRDS(sig_peaks_file)
        
        # Print DiffBind object information
        print("DiffBind object summary:")
        print(all_peaks)
        print(sig_peaks)
        
        # Convert DBA objects to GRanges
        print("Converting DBA objects to GRanges...")
        all_peaks_gr <- dba_to_granges(all_peaks)
        sig_peaks_gr <- dba_to_granges(sig_peaks)
        
        # Standardize chromosome names
        print("Standardizing chromosome names...")
        all_peaks_gr <- standardize_chromosomes(all_peaks_gr)
        sig_peaks_gr <- standardize_chromosomes(sig_peaks_gr)
        
        # Get TxDb object
        print("Loading TxDb...")
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        
        # Annotate peaks
        print("Annotating peaks...")
        peak_anno <- annotatePeak(sig_peaks_gr, 
                                TxDb = txdb,
                                tssRegion = c(-3000, 3000),
                                verbose = FALSE)
        
        # Save annotation results
        print("Saving annotation results...")
        saveRDS(peak_anno, file.path(base_dir, "peak_annotation.rds"))
        write.csv(as.data.frame(peak_anno), 
                 file.path(base_dir, "tables", "peak_annotation.csv"),
                 row.names = FALSE)
        
        # Create annotation plots
        print("Creating annotation plots...")
        pdf(file.path(base_dir, "figures", "annotation_plots.pdf"))
        plotAnnoPie(peak_anno)
        plotDistToTSS(peak_anno)
        dev.off()
        
        # Perform GO enrichment analysis
        print("Performing GO enrichment analysis...")
        genes <- peak_anno@anno$geneId
        ego <- enrichGO(gene = genes,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)
        
        if (!is.null(ego) && nrow(ego) > 0) {
            # Save GO results
            write.csv(as.data.frame(ego), 
                     file.path(base_dir, "tables", "go_enrichment.csv"),
                     row.names = FALSE)
            
            # Create GO plots
            pdf(file.path(base_dir, "figures", "go_enrichment_plots.pdf"))
            print(dotplot(ego, showCategory = 20))
            print(emapplot(ego, showCategory = 50))
            print(cnetplot(ego, showCategory = 10))
            dev.off()
        } else {
            print("No significant GO terms found")
        }
        
        print(paste("Annotation analysis completed for", peak_type, "peaks"))
        return(peak_anno)
        
    }, error = function(e) {
        stop(paste("Error in perform_annotation:", e$message))
    })
}

# Main execution
print("Starting peak annotation analysis...")

# Process narrow peaks
print("\nProcessing narrow peaks...")
narrow_anno <- perform_annotation("narrow")

# Process broad peaks
print("\nProcessing broad peaks...")
broad_anno <- perform_annotation("broad")

print("Peak annotation analysis completed successfully")
