#!/usr/bin/env Rscript

# =========================================================================
# Description:
#   This script annotates a combined SHAPE profile file with quality control (QC)
#   flags. It takes a list of RNA names identified as having poor quality and
#   adds a 'QC_flag' column to the main profile data, marking the corresponding
#   entries as 'poor'.
#
#   This script will automatically install required packages ('argparse', 'data.table')
#   from CRAN if they are not found in the current R environment.
#
# Arguments:
#   -p, --profile <path>  Path to the combined profile.txt file from SHAPEprofile_combine.R.
#                         (Required)
#   -l, --list    <path>  Path to the text file containing the list of poor-quality
#                         RNA names (one name per line). (Required)
#   -o, --outfile <path>  Path for the final, annotated output file. (Required)
#
# Example Usage:
#   Rscript annotate_profiles.R \
#       --profile /path/to/combined_profile.txt \
#       --list /path/to/poor_quality_rnas.txt \
#       --outfile /path/to/annotated_profile.txt
# =========================================================================

# --- Automatic Package Installation ---
# Define the list of packages required for this script
required_packages <- c("argparse", "data.table")

# Check if each package is installed. If not, install it.
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Package '", pkg, "' not found. Installing...\n", sep = "")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# --- Load necessary libraries ---
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))


# --- Argument Parsing ---
parser <- ArgumentParser(description = "Annotate SHAPE profiles with QC flags.")
parser$add_argument("-p", "--profile", type = "character", required = TRUE,
                    help = "Path to the combined profile.txt file.")
parser$add_argument("-l", "--list", type = "character", required = TRUE,
                    help = "Path to the list of poor-quality RNA names.")
parser$add_argument("-o", "--outfile", type = "character", required = TRUE,
                    help = "Path for the final annotated output file.")

args <- parser$parse_args()

# --- Print parameters for user verification ---
cat("Starting SHAPE profile annotation...\n")
cat("----------------------------------------\n")
cat("Profile File:", args$profile, "\n")
cat("RNA List File:", args$list, "\n")
cat("Output File:", args$outfile, "\n")
cat("----------------------------------------\n")

# --- Validate Inputs ---
if (!file.exists(args$profile)) {
  stop(paste("Profile file not found:", args$profile))
}
if (!file.exists(args$list)) {
  stop(paste("RNA list file not found:", args$list))
}

# --- Read Data ---
cat("Reading profile data...\n")
profile_data <- fread(args$profile, sep = "\t")
cat("Reading list of poor-quality RNAs...\n")
poor_quality_rnas <- fread(args$list, header = FALSE, sep = "\n", col.names = c("RNA_name"))


# --- Feature Engineering: Extract RNA Name from 'file_name' ---
cat("Extracting RNA names from 'file_name' column...\n")
# This regex captures the text between "..._chunk_[number]_" and "_profile.txt"
profile_data[, RNA_name := sub("^.*_chunk_\\d+_(.*)_profile\\.txt$", "\\1", file_name)]
cat(" -> A new 'RNA_name' column has been created.\n")


# --- Annotation ---
cat("Annotating data based on the new 'RNA_name' column...\n")
profile_data[, QC_flag := "good"] # Default to good
profile_data[RNA_name %in% poor_quality_rnas$RNA_name, QC_flag := "poor"] # Update to poor


# --- Write Output ---
output_dir <- dirname(args$outfile)
if (!dir.exists(output_dir)) {
  cat("Creating output directory:", output_dir, "\n")
  dir.create(output_dir, recursive = TRUE)
}

cat("Writing annotated data to:", args$outfile, "\n")
setcolorder(profile_data, neworder = c("RNA_name", "QC_flag", "file_name"))
fwrite(profile_data, args$outfile, sep = "\t", row.names = FALSE, quote = FALSE)


# --- Summary ---
qc_summary <- profile_data[, .N, by = QC_flag]
cat("\nAnnotation complete.\n")
cat("Summary of QC flags:\n")
print(qc_summary)