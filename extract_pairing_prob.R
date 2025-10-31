#!/usr/bin/env Rscript

# =========================================================================
# Description:
#   This script processes RNAfold dot plot output files (*dp.ps) to calculate a
#   per-nucleotide score. It reads all *dp.ps files from an input directory,
#   parses the pairing probability data, and for each nucleotide, sums the
#   squared probabilities of it being in a pair.
#
#   The final output is a matrix where rows are RNA names and columns are
#   nucleotide positions, containing the calculated scores.
#
# Arguments:
#   -i, --indir   <path>   Path to the directory containing the RNAfold *dp.ps
#                          output files. (Required)
#   -o, --outfile <path>   Path for the final summary CSV file. (Required)
#   -l, --max-len <int>    Maximum sequence length to consider. This determines the
#                          number of columns in the output matrix. (Default: 200)
#
# Example Usage:
#   Rscript extract_pairing_probabilities.R \
#       --indir /path/to/rnafold_output \
#       --outfile /path/to/pairing_probabilities.csv \
#       --max-len 300
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
parser <- ArgumentParser(description = "Extract and summarize pairing probabilities from RNAfold *dp.ps files.")
parser$add_argument("-i", "--indir", type = "character", required = TRUE,
                    help = "Directory containing *dp.ps files.")
parser$add_argument("-o", "--outfile", type = "character", required = TRUE,
                    help = "Path for the output summary CSV file.")
parser$add_argument("-l", "--max-len", type = "integer", default = 130,
                    help = "Maximum sequence length for the output matrix. [Default: %(default)s]")

args <- parser$parse_args()

# --- Print parameters for user verification ---
cat("Starting pairing probability extraction...\n")
cat("----------------------------------------\n")
cat("Input Directory:", args$indir, "\n")
cat("Output File:", args$outfile, "\n")
cat("Max Sequence Length:", args$max_len, "\n")
cat("----------------------------------------\n")

# --- Validate Inputs and Find Files ---
if (!dir.exists(args$indir)) {
  stop(paste("Input directory not found:", args$indir))
}

file_list <- list.files(path = args$indir, pattern = "_dp\\.ps$", full.names = TRUE)

if (length(file_list) == 0) {
  stop(paste("No '*_dp.ps' files found in directory:", args$indir))
}
cat("Found", length(file_list), "files to process.\n\n")

# --- Core Logic: Process files and aggregate data ---

all_results <- list()

for (file_path in file_list) {
  # cat("Processing:", basename(file_path), "\n")
  
  # Read file and extract lines with pairing probability data ("ubox" lines)
  lines <- readLines(file_path)
  ubox_lines <- grep("ubox$", lines, value = TRUE)
  
  if (length(ubox_lines) == 0) {
    cat(" -> Warning: No 'ubox' lines found. Skipping file.\n")
    next
  }
  
  # Efficiently parse all relevant lines into a data.table
  parsed_data <- lapply(ubox_lines, function(line) {
    parts <- as.numeric(strsplit(line, "\\s+")[[1]][1:3])
    # Returns a list for each line: c(i, j, p)
    if(length(parts) == 3) list(i = parts[1], j = parts[2], p = parts[3]) else NULL
  })
  
  # Remove any NULL entries from failed parsing
  parsed_data <- rbindlist(Filter(Negate(is.null), parsed_data))
  
  if (nrow(parsed_data) == 0) {
    cat(" -> Warning: Failed to parse 'ubox' lines. Skipping file.\n")
    next
  }

  # For each pair (i, j) with probability p, we need to add p^2 to both i and j.
  # This is done by creating two tables and binding them.
  dt1 <- parsed_data[, .(nucleotide = i, value = p^2)]
  dt2 <- parsed_data[, .(nucleotide = j, value = p^2)]
  long_format_dt <- rbind(dt1, dt2)

  # Sum the squared probabilities for each nucleotide within this file
  summed_probs <- long_format_dt[nucleotide <= args$max_len, 
                                 .(sum_sq_prob = sum(value)), 
                                 by = nucleotide]

  # Add the RNA name (derived from the filename)
  rna_name <- sub("_dp\\.ps$", "", basename(file_path))
  summed_probs[, rna_name := rna_name]

  all_results <- append(all_results, list(summed_probs))
}

cat("\nAll files processed. Aggregating results...\n")

# --- Final Aggregation and Reshaping ---
# Combine all results into one long-format table
master_table <- rbindlist(all_results)

# Reshape the data from long to wide format using dcast
# Rows = rna_name, Columns = nucleotide, Values = sum_sq_prob
wide_format <- dcast(master_table, rna_name ~ nucleotide, value.var = "sum_sq_prob", fill = 0)

# Ensure all columns from 1 to max_len exist
all_cols <- as.character(1:args$max_len)
missing_cols <- setdiff(all_cols, names(wide_format))
if (length(missing_cols) > 0) {
  wide_format[, (missing_cols) := 0]
}

# Order the columns numerically
setcolorder(wide_format, c("rna_name", all_cols))

# --- Write Output ---
cat("Writing final matrix to:", args$outfile, "\n")
output_dir <- dirname(args$outfile)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

fwrite(wide_format, args$outfile, sep = ",", row.names = FALSE)

cat("\nSuccess! Processing complete.\n")

