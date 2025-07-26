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
  input_position = 59, # Use c(59, 60) to test multiple positions
  offset_min = 0:0,
  offset_max = 1:2,
  unpaired1_min = 1:2,
  unpaired1_max = 3:3,
  paired1_min = 2:2,
  paired1_max = 3:4,
  unpaired2_min = 3:3,
  unpaired2_max = 4:4,
  include_paired2 = c(FALSE, TRUE)
)

PAIRED2_RANGES <- list(
  paired2_min = 2:2,
  paired2_max = 3:3
)

# Capture warnings
warning_messages <- character()
withCallingHandlers(
  {


# Set up parallel strategy
parallel_strategy <- function(workers) plan(callr, workers = workers)

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

# Use the parsed options
input_file <- opt$input
fold_dir <- opt$fold_dir
output_dir <- opt$output_dir
dataset_file <- opt$dataset


# Create output directory:
dir_create(output_dir)

# Load the dataset
pool1_df <- read.csv(dataset_file)

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
  
  return(list(
    UNUAR_both = unuar_both_f1,
    UGUAG_both = uguag_both_f1
  ))
}

# Print the output directory
cat("Results will be saved in:", output_dir, "\n")

# Add a timestamp for the start of the process
cat("Process started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

run_analysis <- function(unpaired1_min, unpaired1_max, 
                         paired1_min, paired1_max, 
                         unpaired2_min, unpaired2_max, 
                         paired2_min = NULL, paired2_max = NULL, 
                         include_paired2 = FALSE,
                         input_position, offset_min, offset_max) {
  
  # Create a unique identifier for this run
  run_id <- paste0(
    "pos", input_position, "_",
    "off", offset_min, "-", offset_max, "_",
    "u1_", unpaired1_min, "-", unpaired1_max, "_",
    "p1_", paired1_min, "-", paired1_max, "_",
    "u2_", unpaired2_min, "-", unpaired2_max,
    if(include_paired2) paste0("_p2_", paired2_min, "-", paired2_max) else "",
    "_", if(include_paired2) "with_p2" else "no_p2"
  )
  
  # Create a subfolder for the output of this run
  run_dir <- path(output_dir, "individual_results", run_id)
  dir_create(run_dir)
  
  output_file <- path(run_dir, "motif_matches.csv")
  
  motif_matcher_path <- file.path(SCRIPT_DIR, "motif_matcher_v2.R")
  
  if (!file.exists(motif_matcher_path)) {
    stop(paste("motif_matcher_v2.R not found at", motif_matcher_path))
  }

  # Construct the command to run motif_matcher_v2.R
  cmd_args <- c(
    motif_matcher_path,
    "--min_unpaired1", unpaired1_min,
    "--max_unpaired1", unpaired1_max,
    "--min_paired1", paired1_min,
    "--max_paired1", paired1_max,
    "--min_unpaired2", unpaired2_min,
    "--max_unpaired2", unpaired2_max,
    if(include_paired2) c(
      "--min_paired2", paired2_min,
      "--max_paired2", paired2_max,
      "--include_paired2", "TRUE"
    ) else c("--include_paired2", "FALSE"),
    "--input_position", input_position,
    "--offset_min", offset_min,
    "--offset_max", offset_max,
    "--output", output_file,
    "--input", input_file,
    "--fold_dir", fold_dir
  )
  
  # Run the command and capture output
  cmd_output <- system2("Rscript", args = cmd_args, stdout = TRUE, stderr = TRUE)
  
  # Check if the output file exists
  if (!file.exists(output_file)) {
    cat("Error: Output file not created for run:", run_id, "\n")
    cat("Command output:\n", paste(cmd_output, collapse = "\n"), "\n")
    return(NULL)
  }

  # Read the results
  df <- read.csv(output_file)
  
  # Filter pool1_df to match the sequences in df
  filtered_pool1_df <- pool1_df %>%
    filter(sapply(chr, function(x) any(str_detect(df$filename, fixed(x)))))
  
  # Calculate various prediction scores
  structure_motif <- rep(1, nrow(df))  # All rows in df represent found motifs
  unuar_and_structure <- ifelse(filtered_pool1_df$UNUAR == 1 & structure_motif == 1, 1, 0)
  uguag_and_structure <- ifelse(filtered_pool1_df$UGUAG == 1 & structure_motif == 1, 1, 0)
  
  # Calculate basic statistics
  total_sequences <- nrow(pool1_df)
  total_matches <- nrow(df)
  both_count <- sum(filtered_pool1_df$both == 1)
  incell_count <- sum(filtered_pool1_df$incell == 1)
  invitro_count <- sum(filtered_pool1_df$invitro == 1)
  
  # Calculate F1 scores for each prediction method
  analysis_results <- data.frame(
    total_sequences = total_sequences,
    total_matches = total_matches,
    both_count = both_count,
    incell_count = incell_count,
    invitro_count = invitro_count,
    f1_both_structure = calculate_f1_score(filtered_pool1_df$both, structure_motif),
    f1_both_unuar_and_structure = calculate_f1_score(filtered_pool1_df$both, unuar_and_structure),
    f1_both_uguag_and_structure = calculate_f1_score(filtered_pool1_df$both, uguag_and_structure),
    input_position = input_position,
    offset_min = offset_min,
    offset_max = offset_max,
    unpaired1_min = unpaired1_min,
    unpaired1_max = unpaired1_max,
    paired1_min = paired1_min,
    paired1_max = paired1_max,
    unpaired2_min = unpaired2_min,
    unpaired2_max = unpaired2_max,
    paired2_min = ifelse(include_paired2, paired2_min, NA),
    paired2_max = ifelse(include_paired2, paired2_max, NA),
    include_paired2 = include_paired2
  )
  
  # Save analysis results in the run-specific folder
  analysis_results_file <- path(run_dir, "analysis_results.csv")
  tryCatch({
    write.csv(analysis_results, analysis_results_file, row.names = FALSE)
    cat("Analysis results saved to:", analysis_results_file, "\n")
  }, error = function(e) {
    cat("Error saving analysis results:", e$message, "\n")
  })

  #cat("Completed run:", run_id, "\n")
  return(analysis_results)
}

# Function to run limited combinations in parallel on multiple cores
run_limited_combinations <- function() {
  # Set up parallel processing
  parallel_strategy(future::availableCores() - 1)

  # Calculate F1 scores for UNUAR and UGUAG
  sequence_f1_scores <- calculate_sequence_f1_scores(pool1_df)

  ## NEW

  # Create base parameter combinations without paired2
  parameter_combinations <- expand.grid(PARAMETER_RANGES)

  # Create paired2 combinations
  paired2_combinations <- expand.grid(PAIRED2_RANGES)

  # For combinations where include_paired2 is TRUE, add paired2 parameters
  paired2_expanded <- merge(
    parameter_combinations[parameter_combinations$include_paired2, ],
    paired2_combinations
  )

  # For combinations where include_paired2 is FALSE, add NA for paired2 parameters
  non_paired2 <- parameter_combinations[!parameter_combinations$include_paired2, ]
  non_paired2$paired2_min <- NA
  non_paired2$paired2_max <- NA

  # Combine all combinations
  all_combinations <- rbind(paired2_expanded, non_paired2)

  # Filter out invalid combinations and ensure uniqueness
  valid_combinations <- unique(all_combinations[
    all_combinations$unpaired1_min <= all_combinations$unpaired1_max &
    all_combinations$paired1_min <= all_combinations$paired1_max &
    all_combinations$unpaired2_min <= all_combinations$unpaired2_max &
    all_combinations$offset_min <= all_combinations$offset_max &
    (!all_combinations$include_paired2 | 
     (all_combinations$paired2_min <= all_combinations$paired2_max)),
  ])

  cat("Number of valid combinations to test:", nrow(valid_combinations), "\n")
  
  # Run analysis in parallel
  results <- future_lapply(1:nrow(valid_combinations), function(i) {
    params <- valid_combinations[i,]
    run_analysis(
      params$unpaired1_min, params$unpaired1_max,
      params$paired1_min, params$paired1_max,
      params$unpaired2_min, params$unpaired2_max,
      params$paired2_min, params$paired2_max,
      !is.na(params$paired2_min),
      params$input_position,
      params$offset_min,
      params$offset_max
    )
  }, future.seed = TRUE)

  # Remove NULL results
  results <- results[!sapply(results, is.null)]

  # Check if we have any valid results
  if (length(results) == 0) {
    cat("Error: No valid results were obtained.\n")
    return(NULL)
  }

  # Combine results
  all_results <- do.call(rbind, lapply(list.files(path(output_dir, "individual_results"), full.names = TRUE), function(dir) {
  read.csv(path(dir, "analysis_results.csv"))
  }))

  # Save all results
  write.csv(all_results, path(output_dir, "all_analysis_results.csv"), row.names = FALSE)

  cat("Number of valid results:", nrow(all_results), "\n")

  return(all_results)
}


# Run limited combinations
cat("Starting limited parameter sweep. Results will be saved in:", output_dir, "\n")
run_limited_combinations()
cat("Limited parameter sweep completed. Results are saved in:", output_dir, "\n")

# Add this new print statement
cat("\nSummary of output files:\n")
cat("All analysis results:", path(output_dir, "all_analysis_results.csv"), "\n")
cat("Individual analysis results for each parameter combination are also saved in this directory.\n")

# Add a timestamp for the start of the process
cat("Process finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


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