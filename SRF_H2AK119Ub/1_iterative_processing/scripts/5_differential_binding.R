source("scripts/utils.R")

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
    
    # Create DiffBind object and analyze
    dba_data <- dba(sampleSheet=samples)
    dba_data <- dba.count(dba_data)
    dba_data <- dba.normalize(dba_data)
    dba_data <- dba.contrast(dba_data, categories=DBA_CONDITION)
    dba_data <- dba.analyze(dba_data)
    
    # Extract and save results
    results <- dba.report(dba_data)
    saveRDS(dba_data, file.path("analysis", paste0("diffbind_", peak_type), "all_peaks.rds"))
    saveRDS(results, file.path("analysis", paste0("diffbind_", peak_type), "significant_peaks.rds"))
    
    # Create and save plots
    plots_dir <- file.path("analysis", paste0("plots_", peak_type))
    create_dirs(plots_dir)
    
    df_results <- as.data.frame(results)
    save_plot(create_ma_plot(df_results), 
             file.path(plots_dir, "ma_plot.pdf"))
    save_plot(create_volcano_plot(df_results), 
             file.path(plots_dir, "volcano_plot.pdf"))
    
    # Save summary statistics
    stats <- calculate_summary_stats(df_results)
    write.csv(stats, 
              file.path("analysis", paste0("diffbind_", peak_type), "summary_stats.csv"),
              row.names = FALSE)
    
    log_message(sprintf("Completed differential binding analysis for %s peaks", peak_type))
    return(list(dba = dba_data, results = results))
}

# Run analysis for both peak types
broad_results <- perform_diffbind("broad")
narrow_results <- perform_diffbind("narrow") 