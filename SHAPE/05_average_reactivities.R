#!/usr/bin/env Rscript

# =========================================================================
# Description:
#   This script reads one or more annotated SHAPE profile files, combines them,
#   and calculates an inverse-variance weighted average reactivity for each
#   nucleotide position across replicates. It filters for high-quality data
#   (QC_flag == 'good', finite Norm_profile) before averaging.
#
#   The final output is a set of individual .shape files (one for each RNA),
#   which is a standard format for RNA structure prediction software.
#
# Arguments:
#   -i, --infiles <paths>  One or more paths to annotated profile files, separated
#                          by spaces. (Required)
#   -o, --outdir  <path>   Path to the directory where the output .shape files
#                          will be saved. (Required)
#
# Example Usage:
#   Rscript average_reactivities.R \
#       --infiles rep1_annotated.txt rep2_annotated.txt \
#       --outdir /path/to/final_shape_files
# =========================================================================

# --- Automatic Package Installation ---
required_packages <- c("argparse", "data.table")
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
parser <- ArgumentParser(description = "Calculate weighted average reactivities and generate .shape files.")
parser$add_argument("-i", "--infiles", nargs = '+', type = "character", required = TRUE,
                    help = "One or more space-separated paths to annotated profile files.")
parser$add_argument("-o", "--outdir", type = "character", required = TRUE,
                    help = "Path to the output directory for .shape files.")

args <- parser$parse_args()

# --- Print parameters for user verification ---
cat("Starting SHAPE profile averaging and .shape file generation...\n")
cat("------------------------------------------------------------\n")
cat("Input Files:\n")
for(f in args$infiles) cat("  -", f, "\n")
cat("Output Directory:", args$outdir, "\n")
cat("------------------------------------------------------------\n")

# =========================================================================
# STEP 1: Load and Combine Data
# =========================================================================

for(f in args$infiles) {
  if (!file.exists(f)) stop(paste("Input file not found:", f))
}

# Read all files into a list of data.tables
file_list <- lapply(args$infiles, fread)
names(file_list) <- sub("\\.txt$", "", basename(args$infiles))

# Combine the list into a single data.table with a 'replicate' column
combined_dt <- rbindlist(file_list, idcol = "replicate")
cat(" -> Combined", length(args$infiles), "files into a single table with", nrow(combined_dt), "rows.\n\n")

# =========================================================================
# STEP 2: Filter Data for High-Quality Values
# =========================================================================
initial_rows <- nrow(combined_dt)

filtered_dt <- combined_dt[
  QC_flag == "good" & 
    is.finite(Norm_profile) & 
    Norm_profile != -999
]

num_removed <- initial_rows - nrow(filtered_dt)
cat(paste("Filtering: Removed", num_removed, "rows due to 'poor' QC, NA, or -999 values.\n"))

# =========================================================================
# STEP 3: Calculate Inverse-Variance Weighted Average Reactivity
# =========================================================================
averaged_reactivities <- filtered_dt[, .(
  avg_reactivity = {
    if (.N > 1) {
      weights <- 1 / (Norm_stderr^2)
      if (any(!is.finite(weights))) {
        weights[!is.finite(weights)] <- 1e9 
      }
      sum(Norm_profile * weights) / sum(weights)
    } else {
      Norm_profile
    }
  }
), by = .(RNA_name, Nucleotide)]

cat(paste("Calculated weighted averages for", nrow(averaged_reactivities), "nucleotide positions.\n"))

# =========================================================================
# STEP 4: Prepare Full-Length .shape Files (Handle Missing Nucleotides)
# =========================================================================
rna_lengths <- combined_dt[, .(max_nt = max(Nucleotide)), by = RNA_name]
full_template <- rna_lengths[, .(Nucleotide = 1:max_nt), by = RNA_name]

final_shapes <- merge(full_template, averaged_reactivities, by = c("RNA_name", "Nucleotide"), all.x = TRUE)
final_shapes[is.na(avg_reactivity), avg_reactivity := -999]

cat("Generated full-length profiles and filled missing nucleotides with -999.\n")

# =========================================================================
# STEP 5: Write .shape Files for Each RNA
# =========================================================================
if (!dir.exists(args$outdir)) {
  cat(" -> Creating output directory:", args$outdir, "\n")
  dir.create(args$outdir, recursive = TRUE)
}

shape_list <- split(final_shapes, by = "RNA_name")

num_files_written <- 0
for (rna in names(shape_list)) {
  output_file_path <- file.path(args$outdir, paste0(rna, ".shape"))
  data_to_write <- shape_list[[rna]][, .(Nucleotide, avg_reactivity)]
  
  fwrite(data_to_write, 
         file = output_file_path, 
         sep = "\t", 
         col.names = FALSE)
  
  num_files_written <- num_files_written + 1
}

cat(paste("\nProcess complete. Wrote", num_files_written, "'.shape' files to:", args$outdir, "\n"))