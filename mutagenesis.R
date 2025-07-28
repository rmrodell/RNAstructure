#!/usr/bin/env Rscript

cat("mutagenesis.R started\n")
cat("Command line arguments:", paste(commandArgs(trailingOnly = TRUE), collapse = " "), "\n")
flush.console()

# --- Load packages ---
suppressPackageStartupMessages({
  library(optparse)
})

# --- Define command-line options ---
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input CSV from motif_matcher", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "mutants_output.csv",
              help = "Output CSV with mutant sequences [default = %default]", metavar = "character"),
  make_option("--protect_position", type = "integer", default = 60,
              help = "1-based position to exclude from mutation [default = %default]")
)

# --- Parse options ---
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

cat("Parsed options:\n")
print(opt)
flush.console()

input_csv <- opt$input
output_csv <- opt$output
protect_pos <- opt$protect_position

# input_csv <- "/scratch/users/rodell/motifmatcher/20250727/output_test17/individual_results/pos59_off0-1_no_u1_p1_2-4_u2_3-4_p2_2-3/analysis_results.csv"
# output_csv <- "/scratch/users/rodell/motifmatcher/20250727/output_test17/individual_results/pos59_off0-1_no_u1_p1_2-4_u2_3-4_p2_2-3/mutagenesistest.csv"
# protect_pos <- 60

cat("After parsing:\n")
cat("Input file:", input_csv, "\n")
cat("Output file:", output_csv, "\n")
cat("Protect position:", protect_pos, "\n")
flush.console()

if (is.null(input_csv)) {
  stop("Input CSV is required. Use --input to specify the file.")
}

cat("Checking input file:\n")
cat("File exists:", file.exists(input_csv), "\n")
if (file.exists(input_csv)) {
  data <- read.csv(input_csv)
  cat("Rows in input:", nrow(data), "\n")
  cat("Columns in input:", paste(colnames(data), collapse = ", "), "\n")
} else {
  cat("Input file does not exist\n")
}
flush.console()

# --- Load input ---
results <- read.csv(input_csv, stringsAsFactors = FALSE)

# --- Helper functions ---
# Function to convert dot-bracket to base-pair indices
db_to_pairs <- function(dot_bracket) {
  db <- strsplit(dot_bracket, "")[[1]]
  stack <- c()
  pairs <- rep(0, length(db))

  for (i in seq_along(db)) {
    if (db[i] == "(") {
      stack <- c(i, stack)
    } else if (db[i] == ")") {
      if (length(stack) == 0) stop("Unmatched closing parenthesis!")
      j <- stack[1]
      stack <- stack[-1]
      pairs[i] <- j
      pairs[j] <- i
    }
  }
  return(pairs)
}

# Dunctions to generate mutations
generate_mutations <- function(seq_vec, pairs, start, end, protect_pos) {
  if (is.na(start) || is.na(end)) {
    return(list(disruption = NA, compensatory = NA))
  }
  
  indices <- start:end
  indices <- indices[indices != protect_pos]
  pair_indices <- pairs[indices]
  valid_indices <- indices[pairs[indices] != 0 & pairs[indices] != protect_pos]

  if (length(valid_indices) > 0) {
    seq_disrupt <- seq_vec
    seq_disrupt[valid_indices] <- seq_vec[pairs[valid_indices]]
    disruption <- paste0(seq_disrupt, collapse = "")

    seq_comp <- seq_disrupt
    seq_comp[pairs[valid_indices]] <- seq_vec[valid_indices]
    compensatory <- paste0(seq_comp, collapse = "")

    return(list(disruption = disruption, compensatory = compensatory))
  } else {
    return(list(disruption = NA, compensatory = NA))
  }
}

# --- Initialize columns for mutant sequences ---
results$paired1_disruption <- NA
results$paired1_compensatory <- NA
results$paired2_disruption <- NA
results$paired2_compensatory <- NA

# Check if paired2 columns exist and contain non-NA values
has_paired2 <- all(c("paired2_start", "paired2_end") %in% names(results)) &&
               any(!is.na(results$paired2_start)) && any(!is.na(results$paired2_end))
if (has_paired2) {
  cat("Paired2 columns found with non-NA values. Processing both paired1 and paired2.\n")
} else {
  cat("Paired2 columns not found or contain only NA values. Processing only paired1.\n")
}

# --- Process each sequence ---
for (row in 1:nrow(results)) {
  seq_vec <- strsplit(results$fasta_sequence[row], "")[[1]]
  dbn <- results$dot_bracket[row]
  pairs <- db_to_pairs(dbn)

  # Paired1
  paired1_mutations <- generate_mutations(
    seq_vec, pairs, 
    results$paired1_start[row], results$paired1_end[row], 
    protect_pos
  )
  results$paired1_disruption[row] <- paired1_mutations$disruption
  results$paired1_compensatory[row] <- paired1_mutations$compensatory

  # Paired2 (if exists and contains non-NA values)
  if (has_paired2 && !is.na(results$paired2_start[row]) && !is.na(results$paired2_end[row])) {
    paired2_mutations <- generate_mutations(
      seq_vec, pairs, 
      results$paired2_start[row], results$paired2_end[row], 
      protect_pos
    )
    results$paired2_disruption[row] <- paired2_mutations$disruption
    results$paired2_compensatory[row] <- paired2_mutations$compensatory
  }
}

# --- Save mutated sequences ---
write.csv(results, output_csv, row.names = FALSE)