#!/usr/bin/env Rscript

# ===================================================================================
# Description:
#   Compares SHAPE reactivity profiles from two replicate files. It generates
#   scatter plots of Norm_profile values for matching nucleotides, calculates
#   Pearson correlation, and creates a second plot excluding outliers.
#
# Arguments:
#   --rep1 <path>      Path to the annotated profile file for replicate 1.
#   --rep2 <path>      Path to the annotated profile file for replicate 2.
#   --out-prefix <str> Prefix for the output plot files (e.g., sample name).
# ===================================================================================

# --- Setup ---
required_packages <- c("argparse", "data.table", "ggplot2", "ggrepel")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

# --- Argument Parsing ---
parser <- ArgumentParser(description = "Compare SHAPE reactivity profiles between two replicates.")
parser$add_argument("--rep1", required = TRUE, help = "Annotated profile file for replicate 1.")
parser$add_argument("--rep2", required = TRUE, help = "Annotated profile file for replicate 2.")
parser$add_argument("--out-prefix", required = TRUE, help = "Prefix for output plot files.")

args <- parser$parse_args()

# --- Read and Prepare Data ---
cat("Reading and merging replicate data...\n")
dt1 <- fread(args$rep1)
dt2 <- fread(args$rep2)

# Filter for good quality data and select necessary columns
dt1_good <- dt1[QC_flag == 'good', .(RNA_name, Nucleotide, Norm_profile)]
dt2_good <- dt2[QC_flag == 'good', .(RNA_name, Nucleotide, Norm_profile)]

# Merge the two replicates on matching RNA and nucleotide positions
merged_dt <- merge(dt1_good, dt2_good, by = c("RNA_name", "Nucleotide"), suffixes = c(".rep1", ".rep2"))
merged_dt <- na.omit(merged_dt) # Remove rows where a match wasn't found in both
cat(" -> Found", nrow(merged_dt), "matching high-quality nucleotide positions between replicates.\n")

if (nrow(merged_dt) < 2) {
  stop("Not enough common data points to generate a correlation plot.")
}

# --- Plotting Function ---
generate_plot <- function(data, plot_title, filename) {
  # Calculate Pearson correlation
  corr <- cor(data$Norm_profile.rep1, data$Norm_profile.rep2, method = "pearson")
  corr_label <- paste0("Pearson's r = ", round(corr, 3))

  p <- ggplot(data, aes(x = Norm_profile.rep1, y = Norm_profile.rep2)) +
    geom_point(alpha = 0.2, shape = 16) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    annotate("text", x = -Inf, y = Inf, label = corr_label, hjust = -0.2, vjust = 1.5, size = 5) +
    labs(
      title = plot_title,
      x = "SHAPE Reactivity (Replicate 1)",
      y = "SHAPE Reactivity (Replicate 2)"
    ) +
    theme_bw(base_size = 14, base_family = "sans") +
    theme(
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_fixed() # Ensure 1:1 aspect ratio

  ggsave(filename, plot = p, width = 5, height = 5, dpi = 300)
  cat(" -> Plot saved to:", filename, "\n")
}

# --- Generate Plots ---
sample_name <- basename(args$out_prefix)

# 1. Plot with all data
generate_plot(merged_dt, paste(sample_name, "Replicate Comparison (All Data)"), paste0(args$out_prefix, "_replicate_comparison.png"))

# 2. Plot with outliers removed
# Define outliers as values > 3 SD from the mean of all Norm_profile values
all_norm_values <- c(merged_dt$Norm_profile.rep1, merged_dt$Norm_profile.rep2)
mean_norm <- mean(all_norm_values, na.rm = TRUE)
sd_norm <- sd(all_norm_values, na.rm = TRUE)

# Define outlier bounds
lower_bound <- mean_norm - (3 * sd_norm)
upper_bound <- mean_norm + (3 * sd_norm)

# Filter to keep rows where both values are within the non-outlier range
non_outlier_dt <- merged_dt[
    (Norm_profile.rep1 >= lower_bound & Norm_profile.rep1 <= upper_bound) &
    (Norm_profile.rep2 >= lower_bound & Norm_profile.rep2 <= upper_bound)
]

if (nrow(non_outlier_dt) > 2) {
    generate_plot(non_outlier_dt, paste(sample_name, "Replicate Comparison (Outliers Removed)"), paste0(args$out_prefix, "_replicate_comparison_no_outliers.png"))
} else {
    cat(" -> Not enough data to generate plot with outliers removed.\n")
}

cat("Replicate comparison complete.\n")