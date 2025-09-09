#!/usr/bin/env Rscript

library(data.table)

# --- Run to Extract Poor Quality RNAs ---
# Run from the command line for each sample
# Optimized to work on naming structure of pool1
# grep "possible poor quality reactivity profiles" ${LOG_DIR_REP1}/Rep1_NoPus_chunk_*_shapemapper.log | \
# awk -F':' '{sub(/^[ ]+/, "", $NF); print $NF}' | \
# awk '{gsub(/, /,"\n"); print}' | \
# sort -u > "${OUTPUT_LIST_REP1}"

# --- Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: ./annotate_profiles.R <path_to_profile.txt> <path_to_rna_list.txt> <output_file_name.txt>")
}

# Assign arguments to variables
profile_file <- args[1]
rna_list_file <- args[2]
output_file <- args[3]

cat("Starting annotation process...\n")
cat(" -> Profile file:", profile_file, "\n")
cat(" -> RNA list file:", rna_list_file, "\n")

# --- Read Data ---
cat("Reading profile data...\n")
profile_data <- fread(profile_file, sep = "\t")

cat("Reading list of poor-quality RNAs...\n")
poor_quality_rnas <- fread(rna_list_file, header = FALSE, sep = "\n", col.names = c("rna_name_to_flag"))

# --- Feature Engineering: Extract RNA Name from 'file_name' ---
cat("Extracting RNA names from 'file_name' column...\n")

# This regex captures the text between "..._chunk_[number]_" and "_profile.txt"
# For example, in "Rep2_NoPus_chunk_1_ABCC2_chr10_99808770_profile.txt",
# it will extract "ABCC2_chr10_99808770".
profile_data[, RNA_name := sub("^.*_chunk_\\d+_(.*)_profile\\.txt$", "\\1", file_name)]

cat(" -> A new 'RNA_name' column has been created.\n")

# --- Annotation ---
cat("Annotating data based on the new 'RNA_name' column...\n")

# 1. Initialize a new column 'QC_flag' with a default 'good' value.
profile_data[, QC_flag := "good"]

# 2. Update the 'QC_flag' to 'poor' for rows where the newly extracted RNA_name
#    is present in our list of poor quality RNAs.
profile_data[RNA_name %in% poor_quality_rnas$rna_name_to_flag, QC_flag := "poor"]

# --- Write Output ---
cat("Writing annotated data to:", output_file, "\n")
# We'll move the new columns to the front for easier viewing
setcolorder(profile_data, neworder = c("RNA_name", "QC_flag", "file_name"))
fwrite(profile_data, output_file, sep = "\t", quote = TRUE)

# --- Summary ---
qc_summary <- profile_data[, .N, by = QC_flag]
cat("\nAnnotation complete.\n")
cat("Summary of QC flags:\n")
print(qc_summary)