#!/usr/bin/env Rscript

# ===================================================================================
# Description:
#   Compares SHAPE reactivity profiles from two files. It generates
#   scatter plots of Norm_profile values for matching nucleotides, calculates
#   Pearson correlation, and creates a second plot excluding outliers.
#
# Arguments:
#   --nanopore <path>      Path to the annotated profile file for Nanopore files.
#   --illumina <path>      Path to the annotated profile file for Illumina files.
#   --out-prefix <str>     Prefix for the output plot files (e.g., sample name).
# ===================================================================================

# --- Setup ---
required_packages <- c("argparse", "data.table", "ggplot2", "ggrepel", "pROC")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(pROC))

# --- Argument Parsing ---
parser <- ArgumentParser(description = "Compare SHAPE reactivity profiles between two profile files.")
parser$add_argument("--nanopore", required = TRUE, help = "Annotated profile file for nanopore.")
parser$add_argument("--illumina", required = TRUE, help = "Annotated profile file for illumina.")
parser$add_argument("--out-prefix", required = TRUE, help = "Prefix for output plot files.")

args <- parser$parse_args()

# --- Read and Prepare Data ---
cat("Reading and merging replicate data...\n")
dt1 <- fread(args$nanopore)
dt2 <- fread(args$illumina)

# --- QC Flag Concordance Analysis ---
cat("Analyzing QC flag concordance for common RNA molecules...\n")

# Extract unique RNA_name and its QC_flag from each dataset.
qc1 <- unique(dt1[, .(RNA_name, QC_flag)])
qc2 <- unique(dt2[, .(RNA_name, QC_flag)])

# Merge the two datasets on RNA_name to find common molecules.
merged_qc <- merge(qc1, qc2, by = "RNA_name", suffixes = c(".nanopore", ".illumina"))

# Categorize the QC flag relationship for each common RNA.
merged_qc[, category := fcase(
  QC_flag.nanopore == 'good' & QC_flag.illumina == 'good', "Good in Both",
  QC_flag.nanopore == 'poor' & QC_flag.illumina == 'poor', "Poor in Both",
  default = "Discordant (One Good, One Poor)"
)]

# Count the number of RNA molecules in each category.
qc_summary <- merged_qc[, .N, by = category]

# Rename columns for clarity in the final output.
setnames(qc_summary, c("category", "N"), c("QC_Concordance", "RNA_Count"))

# Print the summary table to the console.
cat("Summary of QC Flag Concordance:\n")
print(qc_summary)
cat("\n")

# ---- Plot Comparison ----

# Filter for good quality data and select necessary columns
dt1_good <- dt1[QC_flag == 'good', .(RNA_name, Nucleotide, Norm_profile)]
dt2_good <- dt2[QC_flag == 'good', .(RNA_name, Nucleotide, Norm_profile)]

# Merge the two files on matching RNA and nucleotide positions
merged_dt <- merge(dt1_good, dt2_good, by = c("RNA_name", "Nucleotide"), suffixes = c(".nanopore", ".illumina"))
merged_dt <- na.omit(merged_dt) # Remove rows where a match wasn't found in both

cat("\n[DIAGNOSTIC] First 6 rows of the merged data:\n")
print(head(merged_dt))
cat("\n")

cat(" -> Found", nrow(merged_dt), "matching high-quality nucleotide positions between files.\n")

if (nrow(merged_dt) < 2) {
  stop("Not enough common data points to generate a correlation plot.")
}

# --- Plotting Function ---
theme_base_custom <- function(base_size = 16) {
  theme_minimal(base_size = base_size, base_family = "sans") %+replace%
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1.2), margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = rel(0.8), margin = margin(b = 5)),
      axis.title = element_text(size = rel(1.0)),
      axis.text = element_text(size = rel(0.8)),
      plot.margin = margin(5, 5, 5, 5)
    )
}
theme_scatter <- function(base_size = 16, base_family = "sans") {
  theme_base_custom(base_size) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
}
generate_plot <- function(data, plot_title, filename) {
  # Calculate Pearson correlation
  corr <- cor(data$Norm_profile.nanopore, data$Norm_profile.illumina, method = "pearson")
  corr_label <- paste0("Pearson's r = ", round(corr, 3))

  p <- ggplot(data, aes(x = Norm_profile.nanopore, y = Norm_profile.illumina)) +
    geom_point(alpha = 0.2, shape = 16) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    annotate("text", x = -Inf, y = Inf, label = corr_label, hjust = -0.2, vjust = 1.5, size = 5) +
    labs(
      title = plot_title,
      x = "SHAPE Reactivity (Nanopore)",
      y = "SHAPE Reactivity (Illumina)"
    ) +
    theme_scatter(base_size = 12, base_family = "sans") +
    coord_fixed() # Ensure 1:1 aspect ratio

  ggsave(filename, plot = p, width = 3, height = 3, dpi = 300)
  cat(" -> Plot saved to:", filename, "\n")
}

# --- Generate Plots ---
sample_name <- basename(args$out_prefix)

# 1. Plot with all data
generate_plot(merged_dt, paste(sample_name, "Replicate Comparison (All Data)"), paste0(args$out_prefix, "_scatter.png"))
generate_plot(merged_dt, paste(sample_name, "Replicate Comparison (All Data)"), paste0(args$out_prefix, "_scatter.pdf"))

# 2. Plot with outliers removed
# Define outliers as values > 3 SD from the mean of all Norm_profile values
all_norm_values <- c(merged_dt$Norm_profile.nanopore, merged_dt$Norm_profile.illumina)
mean_norm <- mean(all_norm_values, na.rm = TRUE)
sd_norm <- sd(all_norm_values, na.rm = TRUE)

# Define outlier bounds
lower_bound <- mean_norm - (3 * sd_norm)
upper_bound <- mean_norm + (3 * sd_norm)

# Filter to keep rows where both values are within the non-outlier range
non_outlier_dt <- merged_dt[
    (Norm_profile.nanopore >= lower_bound & Norm_profile.nanopore <= upper_bound) &
    (Norm_profile.illumina >= lower_bound & Norm_profile.illumina <= upper_bound)
]

if (nrow(non_outlier_dt) > 2) {
    generate_plot(non_outlier_dt, paste(sample_name, "Replicate Comparison (Outliers Removed)"), paste0(args$out_prefix, "_scatter_no_outliers.png"))
    generate_plot(non_outlier_dt, paste(sample_name, "Replicate Comparison (Outliers Removed)"), paste0(args$out_prefix, "_scatter_no_outliers.pdf"))
} else {
    cat(" -> Not enough data to generate plot with outliers removed.\n")
}

# --- ROC Curve Analysis at Multiple Thresholds ---
cat("Generating faceted ROC curve for thresholds 0.3, 0.5, 0.7...\n")

# Define the thresholds to test
thresholds_to_test <- c(0.3, 0.5, 0.7)

# Use lapply to iterate over thresholds, calculating ROC data for each
roc_data_list <- lapply(thresholds_to_test, function(th) {
  # Binarize the data using the current threshold 'th'
  roc_dt <- copy(merged_dt)
  roc_dt[, illumina_class := ifelse(Norm_profile.illumina > th, "reactive", "unreactive")]
  roc_dt[, illumina_class := factor(illumina_class, levels = c("unreactive", "reactive"))]

  # Calculate the ROC curve object
  roc_obj <- roc(response = roc_dt$illumina_class, predictor = roc_dt$Norm_profile.nanopore, quiet = TRUE)

  # Extract coordinates and AUC, then return as a list of data.tables
  coords_dt <- data.table(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  )
  auc_val <- as.numeric(auc(roc_obj))

  return(list(
    coords = coords_dt,
    auc = auc_val,
    threshold_label = paste("Illumina Threshold >", th)
  ))
})

# Combine coordinates for plotting
plot_coords <- rbindlist(lapply(roc_data_list, `[[`, "coords"), idcol = "group_id")
# Create a data.frame for AUC labels
auc_labels <- data.table(
  group_id = seq_along(roc_data_list),
  label = sapply(roc_data_list, function(x) paste0("AUC = ", round(x$auc, 3)))
)

# Assign threshold labels for faceting
threshold_labels_vec <- sapply(roc_data_list, `[[`, "threshold_label")
plot_coords[, threshold_level := factor(threshold_labels_vec[group_id], levels = threshold_labels_vec)]
auc_labels[, threshold_level := factor(threshold_labels_vec[group_id], levels = threshold_labels_vec)]
cat(" -> AUC values | Threshold > 0.3:", round(roc_data_list[[1]]$auc, 3), "| > 0.5:", round(roc_data_list[[2]]$auc, 3), "| > 0.7:", round(roc_data_list[[3]]$auc, 3), "\n")

# Create the faceted plot
faceted_roc_plot <- ggplot(plot_coords, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "steelblue", size = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  # Use the separate data.frame for text annotation
  geom_text(data = auc_labels, aes(x = 0.75, y = 0.25, label = label), size = 5) +
  # Create a separate panel for each threshold
  facet_wrap(~threshold_level) +
  labs(
    title = paste(sample_name, "ROC Curve Comparison"),
    subtitle = "Nanopore predicting binarized Illumina SHAPE reactivity",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_scatter(base_size = 12) +
  coord_fixed()

# Save the faceted plot
roc_filename <- paste0(args$out_prefix, "_ROC")
ggsave(paste0(roc_filename, ".png"), plot = faceted_roc_plot, width = 6, height = 3, dpi = 300)
ggsave(paste0(roc_filename, ".pdf"), plot = faceted_roc_plot, width = 6, height = 3)
cat(" -> Faceted ROC plot saved to:", paste0(roc_filename, ".png/.pdf"), "\n")


cat("Replicate comparison complete.\n")

# --- DIAGNOSTIC: Visualize Predictor Distribution ---
cat("[DIAGNOSTIC] Generating predictor distribution plot...\n")
# Create the binary class just for the 0.5 threshold for this example
diagnostic_dt <- copy(merged_dt)
diagnostic_dt[, illumina_class := ifelse(Norm_profile.illumina > 0.5, "reactive", "unreactive")]

dist_plot <- ggplot(diagnostic_dt, aes(x = Norm_profile.nanopore, fill = illumina_class)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("reactive" = "#e41a1c", "unreactive" = "#377eb8")) +
  labs(
    title = "[DIAGNOSTIC] Nanopore Score Distribution by Illumina Class",
    subtitle = "Ground truth threshold = 0.5",
    x = "Nanopore SHAPE Reactivity",
    fill = "Illumina Class"
  ) +
  theme_minimal(base_size = 14)

dist_filename <- paste0(args$out_prefix, "_DIAGNOSTIC_distribution.png")
ggsave(dist_filename, plot = dist_plot, width = 8, height = 5, dpi = 300)
cat(" -> Distribution plot saved to:", dist_filename, "\n")
# --- END DIAGNOSTIC ---