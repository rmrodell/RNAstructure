#!/usr/bin/env Rscript

# =========================================================================
# Description:
#   This script combines individual SHAPE profile files and annotates them
#   with quality control (QC) flags. It searches a directory for profile files,
#   aggregates them, extracts RNA names, and marks any RNAs found in a "failed list"
#   with a 'poor' QC flag.
#
# Arguments:
#   -d, --dir         <path>    Directory containing the individual profile.txt files.
#   -p, --pattern     <string>  Regex pattern to find profile files (e.g., "Pool2_Rep1_.*_profile\\.txt").
#   --prefix          <string>  The sample prefix used in filenames (e.g., "Pool2_Rep1").
#   --failed-list     <path>    Path to the text file of failed RNA names.
#   -o, --outfile     <path>    Path for the final combined and annotated output file.
# =========================================================================

# --- Automatic Package Installation ---
required_packages <- c("argparse", "data.table")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))

# --- Argument Parsing ---
parser <- ArgumentParser(description = "Combine and annotate SHAPE profile files.")
parser$add_argument("-d", "--dir", required = TRUE, help = "Directory with profile files.")
parser$add_argument("-p", "--pattern", required = TRUE, help = "Regex pattern to find profile files.")
parser$add_argument("--prefix", required = TRUE, help = "Sample prefix for parsing filenames.")
parser$add_argument("--failed-list", required = TRUE, help = "Path to the list of failed RNA names.")
parser$add_argument("-o", "--outfile", required = TRUE, help = "Path for the final annotated output file.")

args <- parser$parse_args()

# --- Print parameters ---
cat("Starting SHAPE profile combination and annotation...\n")
cat("----------------------------------------\n")
cat("Input Directory:", args$dir, "\n")
cat("File Pattern:", args$pattern, "\n")
cat("Sample Prefix:", args$prefix, "\n")
cat("Failed List:", args$failed_list, "\n")
cat("Output File:", args$outfile, "\n")
cat("----------------------------------------\n")

# --- Find and Process Profile Files ---
if (!dir.exists(args$dir)) stop(paste("Input directory not found:", args$dir))

profile_files <- list.files(path = args$dir, pattern = args$pattern, full.names = TRUE)

if (length(profile_files) == 0) {
  stop(paste("No profile files found in", args$dir, "matching pattern:", args$pattern))
}
cat("Found", length(profile_files), "profile files to combine.\n")

# Read all files into a single data.table
combined_data <- rbindlist(lapply(profile_files, function(f) {
  dt <- fread(f, sep = "\t")
  dt[, file_name := basename(f)]
  return(dt)
}), use.names = TRUE, fill = TRUE)

# --- Extract RNA Name ---
cat("Extracting RNA names from 'file_name' column...\n")
# Regex uses the prefix to extract the target name
# Example: "Pool2_Rep1_TARGETNAME_profile.txt" -> "TARGETNAME"
regex_pattern <- paste0("^", args$prefix, "_(.*)_profile\\.txt$")
combined_data[, RNA_name := sub(regex_pattern, "\\1", file_name)]

# --- Annotate with QC Flags ---
cat("Annotating data with QC flags...\n")
if (file.exists(args$failed_list) && file.info(args$failed_list)$size > 0) {
  failed_rnas <- fread(args$failed_list, header = FALSE, col.names = "RNA_name")
  cat(" -> Found", nrow(failed_rnas), "failed RNAs to flag.\n")
  
  combined_data[, QC_flag := "good"]
  combined_data[RNA_name %in% failed_rnas$RNA_name, QC_flag := "poor"]
} else {
  cat(" -> Failed RNA list not found or is empty. Flagging all as 'good'.\n")
  combined_data[, QC_flag := "good"]
}

# --- Write Output ---
output_dir <- dirname(args$outfile)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

setcolorder(combined_data, neworder = c("RNA_name", "QC_flag", "file_name"))
fwrite(combined_data, args$outfile, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Combined and annotated data written to:", args$outfile, "\n")

# --- Summary ---
qc_summary <- combined_data[, .N, by = QC_flag]
cat("\nAnnotation complete.\n")
cat("Summary of QC flags:\n")
print(qc_summary)