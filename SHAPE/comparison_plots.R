# =========================================================================
# SETUP: Load Libraries and Define File Paths
# =========================================================================

# Install packages if you don't have them
# install.packages(c("data.table", "ggplot2"))

library(data.table)
library(ggplot2)

# Define paths to your two replicate files
rep1_file <- "/scratch/groups/nicolemm/rodell/SHAPE/InVitro_Rep1_20241209/shapemapper_fullpool_HosseinRecs/Rep1_NoPUS_profile_annotated.txt"
rep2_file <- "/scratch/groups/nicolemm/rodell/SHAPE/InVitro_Rep2_20250731/shapemapper_fullpool_HosseinRecs/Rep2_NoPUS_profile_annotated.txt"

# Define a directory to save the plots
output_dir <- "/scratch/groups/nicolemm/rodell/SHAPE/InVitro_2rep_Analysis_fullpool/comparison_plots"
dir.create(output_dir, showWarnings = FALSE) # Create directory if it doesn't exist

# =========================================================================
# STEP 1: Read in Profile Data
# =========================================================================
cat("Reading data...\n")
rep1_dt <- fread(rep1_file)
rep2_dt <- fread(rep2_file)
cat("Data loaded successfully.\n")

# =========================================================================
# STEP 2: Reusable Scatter Plot Function
# =========================================================================

#' Create a standardized scatter plot comparing Norm_profile between two replicates.
#'
#' @param merged_data A data.table that has been merged, containing columns
#'   'Norm_profile.x' and 'Norm_profile.y'.
#' @param plot_title The main title for the plot.
#' @param file_prefix A succinct name for the output .png and .pdf files.
#' @param output_path The directory where plots will be saved.
create_scatter_plot <- function(merged_data, plot_title, file_prefix, output_path) {
  
  # Ensure there's data to plot
  if (nrow(merged_data) < 2) {
    cat("Skipping plot:", plot_title, "- Not enough data points.\n")
    return(NULL)
  }
  
  # Calculate correlation and R-squared
  correlation_test <- cor.test(merged_data$Norm_profile.x, merged_data$Norm_profile.y)
  r_value <- correlation_test$estimate
  r_squared <- r_value^2
  
  # Create a label string for the plot
  stats_label <- sprintf("R = %.3f\nR² = %.3f", r_value, r_squared)
  
  # Determine axis limits to maintain a 1:1 aspect ratio
  max_val <- max(c(merged_data$Norm_profile.x, merged_data$Norm_profile.y), na.rm = TRUE)
  min_val <- min(c(merged_data$Norm_profile.x, merged_data$Norm_profile.y), na.rm = TRUE)
  
  p <- ggplot(merged_data, aes(x = Norm_profile.x, y = Norm_profile.y)) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", formula = 'y ~ x') +
    coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
    labs(
      title = plot_title,
      subtitle = paste(nrow(merged_data), "points"),
      x = "Replicate 1 Norm_profile",
      y = "Replicate 2 Norm_profile"
    ) +
    annotate(
      "text", x = min_val, y = max_val, label = stats_label, 
      hjust = 0, vjust = 1, size = 6, fontface = "bold"
    ) +
    theme_minimal(base_family = "sans", base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      aspect.ratio = 1
    )
  
  # Save plots
  ggsave(
    filename = file.path(output_path, paste0(file_prefix, ".png")),
    plot = p, width = 6, height = 6, dpi = 300, create.dir = TRUE
  )
  ggsave(
    filename = file.path(output_path, paste0(file_prefix, ".pdf")),
    plot = p, width = 6, height = 6, create.dir = TRUE
  )
  
  cat("Saved plot:", file_prefix, "\n")
  return(p)
}


# =========================================================================
# PLOT 1: All Data Comparison
# =========================================================================
cat("\n--- Processing all data ---\n")
# Merge by unique nucleotide-RNA combination
merged_all_dt <- merge(rep1_dt, rep2_dt, by = c("RNA_name", "Nucleotide"))
create_scatter_plot(merged_all_dt, "Replicate Comparison: All Data", "rep_comparison_all_data", output_dir)


# =========================================================================
# STEP 3: Filter by QC_flag and Re-plot
# =========================================================================
cat("\n--- Processing QC-filtered data ---\n")
rep1_qc_filtered <- rep1_dt[QC_flag == "good"]
rep2_qc_filtered <- rep2_dt[QC_flag == "good"]

merged_qc_dt <- merge(rep1_qc_filtered, rep2_qc_filtered, by = c("RNA_name", "Nucleotide"))
create_scatter_plot(merged_qc_dt, "Replicate Comparison: QC Filtered (Good Only)", "rep_comparison_qc_filtered", output_dir)


# =========================================================================
# STEP 4: Calculate Read Depth Bins
# =========================================================================
cat("\n--- Calculating read depth bins ---\n")

calculate_depth_bins <- function(dt) {
  # Calculate average depths for each RNA_name
  depth_summary <- dt[, .(
    avg_mod = mean(Modified_effective_depth, na.rm = TRUE),
    avg_unt = mean(Untreated_effective_depth, na.rm = TRUE)
  ), by = RNA_name]
  
  # Find the lower of the two averages
  depth_summary[, min_avg_effective_depth := pmin(avg_mod, avg_unt)]
  
  # Assign bins based on the lower average depth
  depth_summary[, read_bin := fcase(
    min_avg_effective_depth >= 5000, "great",
    min_avg_effective_depth >= 2500, "good",
    min_avg_effective_depth >= 1000, "decent",
    min_avg_effective_depth >= 500, "questionable",
    default = "bad"
  )]
  
  # Merge the bin information back into the original data table
  dt_binned <- merge(dt, depth_summary[, .(RNA_name, read_bin)], by = "RNA_name")
  
  # Print counts
  cat("Read depth bin counts (number of unique RNA_names per bin):\n")
  print(table(depth_summary$read_bin))
  
  return(dt_binned)
}

rep1_binned <- calculate_depth_bins(rep1_dt)
rep2_binned <- calculate_depth_bins(rep2_dt)


# =========================================================================
# STEP 5: Create Plots Based on Read Depth Bins (Corrected)
# =========================================================================
cat("\n--- Generating plots based on read depth bins ---\n")

# Merge the binned data tables for easier filtering
merged_binned_dt <- merge(rep1_binned, rep2_binned, by = c("RNA_name", "Nucleotide"), suffixes = c(".rep1", ".rep2"))

# --- Rename columns before the loops ---
# This is more efficient and avoids modification-by-reference issues inside the loop.
setnames(merged_binned_dt, 
         old = c("Norm_profile.rep1", "Norm_profile.rep2"), 
         new = c("Norm_profile.x", "Norm_profile.y"))


# --- Plots for individual bins (where both reps are in the same bin) ---
# bin_levels <- c("great", "good", "decent", "questionable", "bad")
# for (bin in bin_levels) {
#   # Now we just filter. No renaming is needed inside the loop.
#   filtered_for_bin <- merged_binned_dt[read_bin.rep1 == bin & read_bin.rep2 == bin]
  
#   create_scatter_plot(
#     filtered_for_bin,
#     plot_title = paste("Replicate Comparison: Read Bin '", bin, "'"),
#     file_prefix = paste0("rep_comparison_bin_", bin),
#     output_path = output_dir
#   )
# }


# --- Plots for cumulative bins ---
cat("\n--- Generating plots for cumulative read depth bins ---\n")
cumulative_bins <- list(
  great = c("great"),
  great_good = c("great", "good"),
  great_good_decent = c("great", "good", "decent"),
  great_good_decent_quest = c("great", "good", "decent", "questionable")
)

for (i in seq_along(cumulative_bins)) {
  bin_set <- cumulative_bins[[i]]
  set_name <- names(cumulative_bins)[i]
  
  # Filter the main table which already has the correct column names.
  filtered_for_set <- merged_binned_dt[read_bin.rep1 %in% bin_set & read_bin.rep2 %in% bin_set]
  
  create_scatter_plot(
    filtered_for_set,
    plot_title = paste("Replicate Comparison: Bins '", set_name, "'"),
    file_prefix = paste0("rep_comparison_cumulative_", set_name),
    output_path = output_dir
  )
}

# =========================================================================
# STEP 6: Re-create Plots with Outlier Removal
# =========================================================================

cat("\n--- Regenerating plots with a tailored outlier removal method ---\n")

#' Creates plots after removing outliers using a method tailored for reactivity data.
#'
#' This function first identifies a sensible outlier threshold by calculating the IQR
#' only on the more reactive data points (defaulting to those > 1.0), then applies
#' this threshold to the entire dataset before plotting.
#'
#' @param merged_data A data.table containing 'Norm_profile.x' and 'Norm_profile.y'.
#' @param plot_title The main title for the plot.
#' @param file_prefix A succinct name for the output files.
#' @param output_path The directory where plots will be saved.
#' @param pre_filter_threshold A value to isolate the 'region of interest' for IQR calculation.
create_scatter_plot_with_outlier_removal <- function(merged_data, plot_title, file_prefix, output_path, pre_filter_threshold = 1.0) {
  
  # --- Tailored Outlier Detection ---
  # Isolate the "region of interest" where at least one point is reactive
  pre_filter_threshold <- 2.0
  high_value_data <- merged_data[Norm_profile.x > pre_filter_threshold | Norm_profile.y > pre_filter_threshold]
  
  upper_fence <- Inf # Default to no upper limit
  
  if (nrow(high_value_data) > 10) { 
    # Calculate IQR and upper fence ONLY on these high-value points
    high_values <- c(high_value_data$Norm_profile.x, high_value_data$Norm_profile.y)
    q1 <- quantile(high_values, 0.25, na.rm = TRUE)
    q3 <- quantile(high_values, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    upper_fence <- q3 + (3.0 * iqr)
  } else {
    # Fallback for cases with very few reactive points: use a reasonable hard cap.
    upper_fence <- 10.0
    cat(paste("For plot '", plot_title, "': Not enough reactive points for IQR, using fallback fence of", upper_fence, "\n"))
  }
  
  # Filter the ORIGINAL data where BOTH points are below the new fence
  data_filtered <- merged_data[Norm_profile.x <= upper_fence & Norm_profile.y <= upper_fence]
  
  num_removed <- nrow(merged_data) - nrow(data_filtered)
  cat(paste("For plot '", plot_title, "':\n", sep=""))
  cat(paste("Removing", num_removed, "points.\n"))
  
  # --- Call the original plotting function with the filtered data ---
  create_scatter_plot(
    merged_data = data_filtered,
    plot_title = paste(plot_title, "(Outliers Removed)"),
    file_prefix = paste0(file_prefix, "_outliers_removed"),
    output_path = output_path
  )
}

# --- Re-plot 1: All Data (Outliers Removed) ---
create_scatter_plot_with_outlier_removal(
  merged_all_dt, 
  "Replicate Comparison: All Data", 
  "rep_comparison_all_data", 
  output_dir
)

# --- Re-plot 2: QC Filtered (Outliers Removed) ---
create_scatter_plot_with_outlier_removal(
  merged_qc_dt, 
  "Replicate Comparison: QC Filtered (Good Only)", 
  "rep_comparison_qc_filtered", 
  output_dir
)


# --- Re-plot 3: Cumulative Bins (Outliers Removed) ---
cat("\n--- Regenerating cumulative bin plots with outlier removal ---\n")
for (i in seq_along(cumulative_bins)) {
  bin_set <- cumulative_bins[[i]]
  set_name <- names(cumulative_bins)[i]
  
  filtered_for_set <- merged_binned_dt[read_bin.rep1 %in% bin_set & read_bin.rep2 %in% bin_set]
  
  # Call the new outlier-removing plot function
  create_scatter_plot_with_outlier_removal(
    filtered_for_set,
    plot_title = paste("Replicate Comparison: Bins '", set_name, "'"),
    file_prefix = paste0("rep_comparison_cumulative_", set_name),
    output_path = output_dir
  )
}

cat("\nOutlier removal analysis complete. All plots saved to:", output_dir, "\n")


cat("\nAnalysis complete. All plots saved to:", output_dir, "\n")