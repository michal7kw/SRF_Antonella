# Function to plot overlap analysis between V5 and H2AK119Ub peaks
plot_overlap_analysis <- function(v5_peaks, h2a_peaks, condition) {
  # Ensure both inputs are GRanges objects
  if (!is(v5_peaks, "GRanges") || !is(h2a_peaks, "GRanges")) {
    log_message("Both inputs must be GRanges objects", "ERROR")
    return(NULL)
  }
  
  # Find overlaps
  overlaps <- findOverlaps(v5_peaks, h2a_peaks)
  v5_with_h2a <- unique(queryHits(overlaps))
  h2a_with_v5 <- unique(subjectHits(overlaps))
  
  # Calculate statistics
  n_v5_peaks <- length(v5_peaks)
  n_h2a_peaks <- length(h2a_peaks)
  n_overlaps <- length(v5_with_h2a)
  
  # Create data for summary plot
  summary_data <- data.frame(
    Category = c("V5 Peaks", paste0(condition, " H2AK119Ub Peaks")),
    Total = c(n_v5_peaks, n_h2a_peaks),
    Overlapping = c(n_overlaps, length(h2a_with_v5)),
    Percent = c(100 * n_overlaps / n_v5_peaks, 
                100 * length(h2a_with_v5) / n_h2a_peaks)
  )
  
  # Create summary plot
  summary_plot <- ggplot(summary_data, aes(x = Category, y = Total)) +
    geom_bar(stat = "identity", fill = "lightblue") +
    geom_bar(aes(y = Overlapping), stat = "identity", fill = "darkblue") +
    geom_text(aes(label = sprintf("%.1f%%", Percent), y = Total + max(Total) * 0.05), 
              size = 4) +
    theme_bw() +
    labs(title = paste0("Peak Overlap Summary - ", condition),
         x = "", y = "Number of Peaks") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Create width distribution plot for overlapping peaks
  width_dist_plot <- NULL
  if (n_overlaps > 0) {
    overlapping_v5_peaks <- v5_peaks[v5_with_h2a]
    overlapping_h2a_peaks <- h2a_peaks[h2a_with_v5]
    
    width_data <- data.frame(
      Width = c(width(overlapping_v5_peaks), width(overlapping_h2a_peaks)),
      Type = c(rep("V5 Peaks", length(overlapping_v5_peaks)),
               rep(paste0(condition, " H2AK119Ub"), length(overlapping_h2a_peaks)))
    )
    
    width_dist_plot <- ggplot(width_data, aes(x = Width, fill = Type)) +
      geom_density(alpha = 0.5) +
      scale_x_log10(labels = comma) +
      theme_bw() +
      labs(title = paste0("Width Distribution of Overlapping Peaks - ", condition),
           x = "Peak Width (bp)", y = "Density") +
      theme(legend.position = "bottom")
  }
  
  # Return plots
  return(list(
    summary = summary_plot,
    width_dist = width_dist_plot
  ))
}
