#!/usr/bin/env Rscript

# =========================================================================
# Description:
#   This script parses all .fold files from a specified directory, extracting
#   structural and energetic information for each RNA. It handles the complex format
#   of RNAfold output, including MFE, Centroid, and MEA structures and their
#   associated energies and scores.
#
#   It produces a comprehensive summary CSV file of all successfully parsed data and
#   an optional separate log for any files that failed to parse.
#
# Arguments:
#   -i, --indir      <path>   Path to the directory containing the .fold files.
#                             (Required)
#   -o, --outfile    <path>   Path for the final summary CSV file. (Required)
#   --failed-log <path>   Optional path to save a log of files that could
#                             not be parsed.
#
# Example Usage:
#   Rscript extract_fold_summary.R \
#       --indir /path/to/rnafold_output \
#       --outfile /path/to/rnafold_summary.csv \
#       --failed-log /path/to/rnafold_failures.csv
# =========================================================================

# --- Automatic Package Installation ---
required_packages <- c("argparse", "data.table", "stringr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Package '", pkg, "' not found. Installing...\n", sep = "")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# --- Load necessary libraries ---
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))


# --- Command-Line Argument Parsing ---
parser <- ArgumentParser(description = "Parse RNAfold .fold files into a summary CSV.")
parser$add_argument("-i", "--indir", type = "character", required = TRUE,
                    help = "Input directory containing .fold files.")
parser$add_argument("-o", "--outfile", type = "character", required = TRUE,
                    help = "Output path for the summary CSV file.")
parser$add_argument("--failed-log", type = "character",
                    help = "Optional path for a log of failed files.")

args <- parser$parse_args()

# --- Print parameters for user verification ---
cat("Starting RNAfold .fold file parsing...\n")
cat("----------------------------------------\n")
cat("Input Directory:", args$indir, "\n")
cat("Output File:", args$outfile, "\n")
if (!is.null(args$failed_log)) {
    cat("Failed Files Log:", args$failed_log, "\n")
}
cat("----------------------------------------\n")

# --- Validate Inputs and Find Files ---
if (!dir.exists(args$indir)) {
  stop(paste("Input directory not found:", args$indir))
}

file_list <- list.files(path = args$indir, pattern = "\\.fold$", full.names = TRUE)

if (length(file_list) == 0) {
  stop(paste("No '.fold' files found in directory:", args$indir))
}
cat("Found", length(file_list), "'.fold' files to process.\n\n")

# --- Core Extraction Function ---

extract_fold_info <- function(file_path) {
  tryCatch({
    lines <- readLines(file_path)
    
    if (length(lines) < 7) {
      stop("File does not contain enough lines for full parsing.")
    }
    
    # Define a robust regex for parsing numbers, including scientific notation
    num_regex <- "[-]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?"

    # --- Line-by-Line Extraction (Corrected based on new format) ---
    
    # Lines 1 & 2: ID and Sequence
    sequence_id <- sub(">", "", lines[1])
    sequence <- lines[2]
    
    # Line 3: MFE (Minimum Free Energy)
    mfe_structure <- str_extract(lines[3], "^[().]+")
    mfe_energy <- as.numeric(str_extract(lines[3], sprintf("(?<=\\()\\s*(%s)(?=\\s*\\)$)", num_regex)))
    
    # Line 4: Probabilistic Dot-Plot Representation
    dot_plot_representation <- str_extract(lines[4], "^[().|,{}]+")
    dot_plot_energy <- as.numeric(str_extract(lines[4], sprintf("(?<=\\[)\\s*(%s)(?=\\s*\\]$)", num_regex)))
    
    # Line 5: Centroid Structure
    centroid_structure <- str_extract(lines[5], "^[().]+")
    centroid_energy <- as.numeric(str_extract(lines[5], sprintf("(?<=\\{)\\s*(%s)", num_regex)))
    centroid_distance <- as.numeric(str_extract(lines[5], sprintf("(?<=d=)\\s*(%s)(?=\\s*\\})", num_regex)))
    
    # Line 6: MEA (Maximum Expected Accuracy) Structure
    mea_structure <- str_extract(lines[6], "^[().]+")
    mea_energy <- as.numeric(str_extract(lines[6], sprintf("(?<=\\{)\\s*(%s)(?=\\s+MEA=)", num_regex)))
    mea_score <- as.numeric(str_extract(lines[6], sprintf("(?<=MEA=)\\s*(%s)(?=\\s*\\})", num_regex)))
    
    # Line 7: Ensemble Statistics
    numbers_in_line7 <- str_extract_all(lines[7], num_regex)[[1]]
    mfe_frequency <- as.numeric(numbers_in_line7[1])
    ensemble_diversity <- as.numeric(numbers_in_line7[2])
    
    # --- Validation and Result Assembly ---
    essential_values <- c(mfe_structure, mfe_energy, centroid_structure, centroid_energy, 
                          mea_structure, mea_energy, mea_score, 
                          mfe_frequency, ensemble_diversity)
    if (any(is.na(essential_values))) {
      stop("Failed to extract one or more essential values. Check file format.")
    }
    
    result <- data.table(
      sequence_id = sequence_id,
      sequence = sequence,
      mfe_structure = mfe_structure,
      mfe_energy = mfe_energy,
      dot_plot_representation = dot_plot_representation,
      dot_plot_energy = dot_plot_energy,
      centroid_structure = centroid_structure,
      centroid_energy = centroid_energy,
      centroid_distance = centroid_distance,
      mea_structure = mea_structure,
      mea_energy = mea_energy,
      mea_score = mea_score,
      mfe_frequency = mfe_frequency,
      ensemble_diversity = ensemble_diversity
    )
    
    return(result)
    
  }, error = function(e) {
    return(list(file_name = basename(file_path), error_message = e$message))
  })
}

# --- Main Processing Loop ---
print(paste("Processing", length(file_list), "files..."))
results_list <- lapply(file_list, extract_fold_info)

is_error <- sapply(results_list, function(x) "error_message" %in% names(x))
failed_files_list <- results_list[is_error]
successful_results <- results_list[!is_error]

master_data_table <- rbindlist(successful_results, use.names = TRUE, fill = TRUE)
cat(" -> Successfully processed", nrow(master_data_table), "files.\n")
cat(" -> Encountered errors in", length(failed_files_list), "files.\n\n")



# --- Output and Summary ---
# Ensure output directory exists
output_dir <- dirname(args$outfile)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Write summary of successfully parsed files
cat("Writing summary file to:", args$outfile, "\n")
fwrite(master_data_table, args$outfile, sep = ",", row.names = FALSE)

# Write log of failed files, if requested and if any failed
if (!is.null(args$failed_log) && length(failed_files_list) > 0) {
    failed_dt <- rbindlist(failed_files_list, use.names = TRUE, fill = TRUE)
    cat("Writing failed files log to:", args$failed_log, "\n")
    fwrite(failed_dt, args$failed_log, sep = ",", row.names = FALSE)
}

print("Processing complete")