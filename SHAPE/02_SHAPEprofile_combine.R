#!/usr/bin/env Rscript

# =========================================================================
# Description:
#   This script automates the process of combining SHAPE-Mapper profile.txt files.
#   It searches a parent directory for subdirectories matching a given pattern,
#   finds all "*_profile.txt" files within them, and aggregates them into a
#   single, comprehensive output file. A 'source_file' column is added to
#   track the origin of each row.
#
# Arguments:
#   -d, --dir     <path>    Path to the parent directory containing the SHAPE-Mapper
#                           output subdirectories. (Required)
#   -p, --pattern <string>  A text pattern to identify the target subdirectories
#                           (e.g., "MySample_Rep1_chunk_"). (Required)
#   -o, --outfile <path>    Path for the final combined output file. (Required)
#
# Example Usage:
#   Rscript SHAPEprofile_combine_v2.R \
#       --dir /path/to/shapemapper_outputs \
#       --pattern "Rep1_NoPus_chunk_" \
#       --outfile /path/to/shapemapper_outputs/Rep1_NoPus_combined_profile.txt
# =========================================================================

# --- Automatic Package Installation ---
# Define the list of packages required for this script
required_packages <- c("argparse", "data.table")

# Check if each package is installed. If not, install it.
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Package '", pkg, "' not found. Installing...\n", sep = "")
    # Install the package from a default CRAN mirror
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# --- Load necessary libraries ---
# Suppress startup messages for cleaner output
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))

# --- Argument Parsing ---
parser <- ArgumentParser(description = "Combine SHAPE-Mapper profile.txt files from multiple directories.")
parser$add_argument("-d", "--dir", type = "character", required = TRUE,
                    help = "Parent directory to search for analysis subdirectories.")
parser$add_argument("-p", "--pattern", type = "character", required = TRUE,
                    help = "Pattern to identify target subdirectories (e.g., 'Sample_Rep1_').")
parser$add_argument("-o", "--outfile", type = "character", required = TRUE,
                    help = "Path for the final combined output file.")

# In case of error, Rscript prints help and exits
args <- parser$parse_args()

# --- Print parameters for user verification ---
cat("Starting SHAPE profile aggregation...\n")
cat("----------------------------------------\n")
cat("Parent Directory:", args$dir, "\n")
cat("Subdirectory Pattern:", args$pattern, "\n")
cat("Output File:", args$outfile, "\n")
cat("----------------------------------------\n")

# --- Find target directories ---
if (!dir.exists(args$dir)) {
  stop(paste("Parent directory not found:", args$dir))
}

# List all subdirectories in the parent directory
subdirs <- list.dirs(path = args$dir, full.names = TRUE, recursive = FALSE)

# Filter for directories that match the specified pattern
target_dirs <- grep(args$pattern, subdirs, value = TRUE)

if (length(target_dirs) == 0) {
  stop(paste("No subdirectories found in", args$dir, "matching the pattern:", args$pattern))
}
cat("Found", length(target_dirs), "matching directories to process.\n\n")

# --- Find and process profile.txt files ---
all_data_frames <- list()

# Loop through each directory
for (dir in target_dirs) {
  
  cat("Processing directory:", dir, "\n")
  
  # Get a list of all files that end with profile.txt in the directory
  file_list <- list.files(path = dir, pattern = "profile\\.txt$", full.names = TRUE)

  # Loop through each file found
  for (file_path in file_list) {
    # Read the tab-separated file into a data frame
    df <- fread(file_path, sep = "\t", header = TRUE)

    # Add a new column with the filename (without path) as the first column
    df_file_name <- basename(file_path)
    df$file_name <- df_file_name  # Add the filename as a new column

    # Append the data frame to the list
    all_data_frames[[length(all_data_frames) + 1]] <- df
    cat("Appended data from file:", df_file_name, "\n")  # Debug-Trace
  }
}

# Combine all data frames into one
if (length(all_data_frames) > 0) {

  combined_data <- rbindlist(all_data_frames, use.names = TRUE, fill = TRUE)

  # Ensure the output directory exists
  output_dir <- dirname(args$outfile)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Write the combined data to the output file
  fwrite(combined_data, args$outfile, sep = "\t", row.names = FALSE, quote = FALSE)

  cat("Data successfully concatenated and written to", args$outfile, "\n")
} else {
  cat("No files found matching the pattern across all directories.\n")
}