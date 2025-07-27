#!/usr/bin/env Rscript

.libPaths("/home/users/rodell/R/x86_64-pc-linux-gnu-library/4.2")
library(future)
library(future.apply)
library(future.callr)
library(fs)
library(optparse)
library(dplyr)
library(stringr)

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- dirname(sub("--file=", "", args[grep("--file=", args)]))
  return(normalizePath(script_path))
}

SCRIPT_DIR <- get_script_path()

# Define parameter ranges
PARAMETER_RANGES <- list(
  input_position = 59, 
  offset_min = 0:0,
  offset_max = 1:1,
  include_unpaired1 = c(FALSE, TRUE),
  paired1_min = 2:2,
  paired1_max = 3:4,
  unpaired2_min = 3:3,
  unpaired2_max = 4:4,
  include_paired2 = c(FALSE, TRUE)
)

UNPAIRED1_RANGES <- list(
  unpaired1_min = 1:2,
  unpaired1_max = 3:3
)

PAIRED2_RANGES <- list(
  paired2_min = 2:2,
  paired2_max = 3:3
)

# Capture warnings
warning_messages <- character()
withCallingHandlers(
  {


# --- Define options ---
option_list <- list(
  make_option(c("-i", "--input"), type="character", default="input_pool1.csv",
              help="Input CSV file [default= %default]"),
  make_option(c("-f", "--fold_dir"), type="character", default="/scratch/users/rodell/RNAfold_psipos",
              help="Directory containing .fold files [default= %default]"),
  make_option(c("-o", "--output_dir"), type="character", default="parametersweep_output",
              help="Output directory for results [default= %default]"),
  make_option(c("-d", "--dataset"), type="character", default="pool1_psipos_v2.csv",
              help="Dataset file for evaluation [default= %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# --- Assign options ---
input_file <- opt$input
fold_dir <- opt$fold_dir
output_dir <- opt$output_dir
dataset_file <- opt$dataset


# Create output directory:
dir_create(output_dir)

# Load the dataset
pool1_df <- read.csv(dataset_file)

# --- Helper Functions ---

# Set up parallel strategy
parallel_strategy <- function(workers) plan(callr, workers = workers)

# Determine valid parameter combinations
generate_parameter_combinations <- function() {
  base_combinations <- expand.grid(PARAMETER_RANGES)
  
  # Generate all possible combinations for unpaired1 and paired2
  unpaired1_combinations <- expand.grid(UNPAIRED1_RANGES)
  paired2_combinations <- expand.grid(PAIRED2_RANGES)
  
  # Merge all combinations
  all_combinations <- merge(base_combinations, unpaired1_combinations)
  all_combinations <- merge(all_combinations, paired2_combinations)
  
  # Set values to NA where the region is not included
  all_combinations$unpaired1_min[!all_combinations$include_unpaired1] <- NA
  all_combinations$unpaired1_max[!all_combinations$include_unpaired1] <- NA
  all_combinations$paired2_min[!all_combinations$include_paired2] <- NA
  all_combinations$paired2_max[!all_combinations$include_paired2] <- NA
  
  valid_combinations <- all_combinations[
    (is.na(all_combinations$unpaired1_min) | all_combinations$unpaired1_min <= all_combinations$unpaired1_max) &
    all_combinations$paired1_min <= all_combinations$paired1_max &
    all_combinations$unpaired2_min <= all_combinations$unpaired2_max &
    all_combinations$offset_min <= all_combinations$offset_max &
    (is.na(all_combinations$paired2_min) | all_combinations$paired2_min <= all_combinations$paired2_max),
  ]
  
  # Define the desired column order
  column_order <- c(
    "input_position", "offset_min", "offset_max",
    "include_unpaired1", "unpaired1_min", "unpaired1_max",
    "paired1_min", "paired1_max",
    "unpaired2_min", "unpaired2_max",
    "include_paired2", "paired2_min", "paired2_max"
  )
  
  # Reorder the columns
  valid_combinations <- valid_combinations[, column_order]
  
  return(unique(valid_combinations))
}


# Function to calculate F1 score
calculate_f1_score <- function(actual, predicted) {
  true_positives <- sum(actual == 1 & predicted == 1)
  false_positives <- sum(actual == 0 & predicted == 1)
  false_negatives <- sum(actual == 1 & predicted == 0)
  
  precision <- ifelse(true_positives + false_positives > 0,
                      true_positives / (true_positives + false_positives),
                      0)
  recall <- ifelse(true_positives + false_negatives > 0,
                   true_positives / (true_positives + false_negatives),
                   0)
  
  f1_score <- ifelse(precision + recall > 0,
                     2 * (precision * recall) / (precision + recall),
                     0)
  
  return(f1_score)
}


calculate_sequence_f1_scores <- function(pool1_df) {
  unuar_both_f1 <- calculate_f1_score(pool1_df$both, pool1_df$UNUAR)
  uguag_both_f1 <- calculate_f1_score(pool1_df$both, pool1_df$UGUAG)
  
  return(sprintf("Sequence-only F1 Scores:\nUNUAR vs Both: %.4f\nUGUAG vs Both: %.4f", 
                 unuar_both_f1, uguag_both_f1))
}

create_run_id <- function(input_position, offset_min, offset_max, 
                          include_unpaired1, unpaired1_min, unpaired1_max,
                          paired1_min, paired1_max, 
                          unpaired2_min, unpaired2_max,
                          include_paired2, paired2_min, paired2_max) {
  paste0(
    "pos", input_position, "_",
    "off", offset_min, "-", offset_max, "_",
    if(include_unpaired1) paste0("u1_", unpaired1_min, "-", unpaired1_max, "_") else "no_u1_",
    "p1_", paired1_min, "-", paired1_max, "_",
    "u2_", unpaired2_min, "-", unpaired2_max,
    if(include_paired2) paste0("_p2_", paired2_min, "-", paired2_max) else "_no_p2"
  )
}

run_motif_matcher <- function(motif_matcher_path, output_file, params) {
  cmd_args <- c(
    motif_matcher_path,
    "--min_paired1", params$paired1_min,
    "--max_paired1", params$paired1_max,
    "--min_unpaired2", params$unpaired2_min,
    "--max_unpaired2", params$unpaired2_max,
    "--input_position", params$input_position,
    "--offset_min", params$offset_min,
    "--offset_max", params$offset_max,
    "--output", output_file,
    "--input", params$input_file,
    "--fold_dir", params$fold_dir
  )

  if(params$include_unpaired1) {
    cmd_args <- c(cmd_args,
      "--min_unpaired1", params$unpaired1_min,
      "--max_unpaired1", params$unpaired1_max,
      "--include_unpaired1", "TRUE"
    )
  } else {
    cmd_args <- c(cmd_args, "--include_unpaired1", "FALSE")
  }

  if(params$include_paired2) {
    cmd_args <- c(cmd_args,
      "--min_paired2", params$paired2_min,
      "--max_paired2", params$paired2_max,
      "--include_paired2", "TRUE"
    )
  } else {
    cmd_args <- c(cmd_args, "--include_paired2", "FALSE")
  }
  
  system2("Rscript", args = cmd_args, stdout = TRUE, stderr = TRUE)
}

calculate_statistics <- function(df, filtered_pool1_df) {
  structure_motif <- rep(1, nrow(df))
  unuar_and_structure <- ifelse(filtered_pool1_df$UNUAR == 1 & structure_motif == 1, 1, 0)
  uguag_and_structure <- ifelse(filtered_pool1_df$UGUAG == 1 & structure_motif == 1, 1, 0)
  
  list(
    total_sequences = nrow(pool1_df),
    total_matches = nrow(df),
    both_count = sum(filtered_pool1_df$both == 1),
    incell_count = sum(filtered_pool1_df$incell == 1),
    invitro_count = sum(filtered_pool1_df$invitro == 1),
    f1_both_structure = calculate_f1_score(filtered_pool1_df$both, structure_motif),
    f1_both_unuar_and_structure = calculate_f1_score(filtered_pool1_df$both, unuar_and_structure),
    f1_both_uguag_and_structure = calculate_f1_score(filtered_pool1_df$both, uguag_and_structure)
  )
}

run_analysis <- function(params) {
  run_id <- create_run_id(params$input_position, params$offset_min, params$offset_max, 
                          params$include_unpaired1, params$unpaired1_min, params$unpaired1_max,
                          params$paired1_min, params$paired1_max, 
                          params$unpaired2_min, params$unpaired2_max,
                          params$include_paired2, params$paired2_min, params$paired2_max)
  
  run_dir <- file.path(output_dir, "individual_results", run_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  
  output_file <- file.path(run_dir, "motif_matches.csv")
  motif_matcher_path <- file.path(SCRIPT_DIR, "motif_matcher_v3.R")
  
  if (!file.exists(motif_matcher_path)) {
    stop(paste("motif_matcher_v3.R not found at", motif_matcher_path))
  }

  cmd_output <- run_motif_matcher(motif_matcher_path, output_file, params)
  
  if (!file.exists(output_file)) {
    cat("Error: Output file not created for run:", run_id, "\n")
    cat("Command output:\n", paste(cmd_output, collapse = "\n"), "\n")
    return(NULL)
  }

  df <- read.csv(output_file)
  
  filtered_pool1_df <- pool1_df %>%
    filter(sapply(chr, function(x) any(str_detect(df$filename, fixed(x)))))
  
  stats <- calculate_statistics(df, filtered_pool1_df)
  
  analysis_results <- c(stats, params)
  
  analysis_results_file <- file.path(run_dir, "analysis_results.csv")
  tryCatch({
    write.csv(data.frame(analysis_results), analysis_results_file, row.names = FALSE)
    cat("Analysis results saved to:", analysis_results_file, "\n")
  }, error = function(e) {
    cat("Error saving analysis results:", e$message, "\n")
  })

  return(analysis_results)
}

run_parameter_sweep <- function() {
  # Set up parallel processing
  parallel_strategy(future::availableCores() - 1)

  # Generate valid parameter combinations
  valid_combinations <- generate_parameter_combinations()

  cat("Number of valid combinations to test:", nrow(valid_combinations), "\n")
  
  # Run analysis in parallel
  results <- future_lapply(1:nrow(valid_combinations), function(i) {
    params <- valid_combinations[i,]
    params$input_file <- input_file
    params$fold_dir <- fold_dir
    run_analysis(params)
  }, future.seed = TRUE)

  # Process results
  processed_results <- process_results(results)

  return(processed_results)
}

process_results <- function(results) {
  # Remove NULL results
  results <- results[!sapply(results, is.null)]

  # Check if we have any valid results
  if (length(results) == 0) {
    cat("Error: No valid results were obtained.\n")
    return(NULL)
  }

  # Combine results
  all_results <- do.call(rbind, lapply(list.files(file.path(output_dir, "individual_results"), full.names = TRUE), function(dir) {
    read.csv(file.path(dir, "analysis_results.csv"))
  }))

  # Save all results
  write.csv(all_results, file.path(output_dir, "all_analysis_results.csv"), row.names = FALSE)

  cat("Number of valid results:", nrow(all_results), "\n")

  return(all_results)
}

main <- function() {
  # Print the output directory
  cat("Results will be saved in:", output_dir, "\n")

  # Calculate and print sequence F1 scores
  cat("Calculating sequence-only F1 scores...\n")
  sequence_f1_scores <- calculate_sequence_f1_scores(pool1_df)
  cat(sequence_f1_scores, "\n\n")

  # Add a timestamp for the start of the process
  cat("Process started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  # Run parameter sweep
  cat("Starting parameter sweep. Results will be saved in:", output_dir, "\n")
  results <- run_parameter_sweep()
  cat("Parameter sweep completed. Results are saved in:", output_dir, "\n")

  # Print summary of output files
  cat("\nSummary of output files:\n")
  cat("All analysis results:", file.path(output_dir, "all_analysis_results.csv"), "\n")
  cat("Individual analysis results for each parameter combination are also saved in this directory.\n")

  # Add a timestamp for the end of the process
  cat("Process finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  return(results)
}

# --- Run the main function ---

results <- main()

  },
  warning = function(w) {
    warning_messages <<- c(warning_messages, w$message)
  }
)

# Print captured warnings at the end
if (length(warning_messages) > 0) {
  cat("\nWarnings encountered during execution:\n")
  for (i in seq_along(warning_messages)) {
    cat(i, ": ", warning_messages[i], "\n", sep = "")
  }
} else {
  cat("\nNo warnings encountered during execution.\n")
}