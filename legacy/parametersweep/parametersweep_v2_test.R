#!/usr/bin/env Rscript

.libPaths("/home/users/rodell/R/x86_64-pc-linux-gnu-library/4.2")
library(future)
library(future.apply)
library(future.callr)
library(fs)
library(optparse)
library(dplyr)
library(stringr)

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
  make_option(c("-d", "--dataset"), type="character", default="pool1_psipos_info.csv",
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
  
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  
  f1_score <- 2 * (precision * recall) / (precision + recall)
  return(f1_score)
}

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
  unuar_incell_f1 <- calculate_f1_score(pool1_df$incell, pool1_df$UNUAR)
  unuar_invitro_f1 <- calculate_f1_score(pool1_df$invitro, pool1_df$UNUAR)
  uguag_incell_f1 <- calculate_f1_score(pool1_df$incell, pool1_df$UGUAG)
  uguag_invitro_f1 <- calculate_f1_score(pool1_df$invitro, pool1_df$UGUAG)
  
  return(list(
    UNUAR_incell = unuar_incell_f1,
    UNUAR_invitro = unuar_invitro_f1,
    UGUAG_incell = uguag_incell_f1,
    UGUAG_invitro = uguag_invitro_f1
  ))
}

# Print the output directory
cat("Results will be saved in:", output_dir, "\n")

# Add a timestamp for the start of the process
cat("Process started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


# Function to run the analysis with specific parameters
run_analysis <- function(unpaired1_min, unpaired1_max, 
                         paired1_min, paired1_max, 
                         unpaired2_min, unpaired2_max, 
                         paired2_min = NULL, paired2_max = NULL, 
                         include_paired2 = FALSE) {
  
  # Create a unique identifier for this run
  run_id <- paste0(
    "u1_", unpaired1_min, "-", unpaired1_max, "_",
    "p1_", paired1_min, "-", paired1_max, "_",
    "u2_", unpaired2_min, "-", unpaired2_max,
    if(include_paired2) paste0("_p2_", paired2_min, "-", paired2_max) else "",
    "_", if(include_paired2) "with_p2" else "no_p2"
  )
  
  output_file <- path(output_dir, paste0(run_id, ".csv"))
  
  # Construct the command to run motif_matcher_v2.R
  cmd_args <- c(
    "motif_matcher_v2.R",
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
    "--output", output_file,
    "--input", input_file,
    "--fold_dir", fold_dir
  )
  
  # Run the command
  system2("Rscript", args = cmd_args)
  
  # Check if the output file exists
  if (!file.exists(output_file)) {
    cat("Error: Output file not created for run:", run_id, "\n")
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
  incell_count <- sum(filtered_pool1_df$incell == 1)
  invitro_count <- sum(filtered_pool1_df$invitro == 1)
  
  # Calculate F1 scores for each prediction method
  analysis_results <- data.frame(
    total_sequences = total_sequences,
    total_matches = total_matches,
    incell_count = incell_count,
    invitro_count = invitro_count,
    f1_incell_structure = calculate_f1_score(filtered_pool1_df$incell, structure_motif),
    f1_invitro_structure = calculate_f1_score(filtered_pool1_df$invitro, structure_motif),
    f1_incell_unuar_and_structure = calculate_f1_score(filtered_pool1_df$incell, unuar_and_structure),
    f1_invitro_unuar_and_structure = calculate_f1_score(filtered_pool1_df$invitro, unuar_and_structure),
    f1_incell_uguag_and_structure = calculate_f1_score(filtered_pool1_df$incell, uguag_and_structure),
    f1_invitro_uguag_and_structure = calculate_f1_score(filtered_pool1_df$invitro, uguag_and_structure)
  )
  
  
  # Add the parameters to the results
  analysis_results$unpaired1_min <- unpaired1_min
  analysis_results$unpaired1_max <- unpaired1_max
  analysis_results$paired1_min <- paired1_min
  analysis_results$paired1_max <- paired1_max
  analysis_results$unpaired2_min <- unpaired2_min
  analysis_results$unpaired2_max <- unpaired2_max
  analysis_results$paired2_min <- ifelse(include_paired2, paired2_min, NA)
  analysis_results$paired2_max <- ifelse(include_paired2, paired2_max, NA)
  analysis_results$include_paired2 <- include_paired2
  
  cat("Completed run:", run_id, "\n")
  return(analysis_results)
}

# Function to run limited combinations in parallel on multiple cores
run_limited_combinations <- function() {
  # Set up parallel processing
  parallel_strategy(future::availableCores() - 1)

  # Calculate F1 scores for UNUAR and UGUAG
  sequence_f1_scores <- calculate_sequence_f1_scores(pool1_df)

  # Create base parameter combinations without paired2
  parameter_combinations <- expand.grid(
    unpaired1_min = 1:2, unpaired1_max = 2:5,
    paired1_min = 1:1, paired1_max = 5:7,
    unpaired2_min = 3:3, unpaired2_max = 7:7,
    include_paired2 = c(FALSE, TRUE)
  )

  # Create paired2 combinations
  paired2_combinations <- expand.grid(
    paired2_min = 3:3,
    paired2_max = 7:7
  )

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
    all_combinations$unpaired1_min < all_combinations$unpaired1_max &
    all_combinations$paired1_min < all_combinations$paired1_max &
    all_combinations$unpaired2_min < all_combinations$unpaired2_max &
    (!all_combinations$include_paired2 | 
     (all_combinations$paired2_min < all_combinations$paired2_max)),
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
      !is.na(params$paired2_min)
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
  all_results <- do.call(rbind, results)

  # Save all results
  write.csv(all_results, path(output_dir, "all_analysis_results.csv"), row.names = FALSE)

  cat("Number of valid results:", nrow(all_results), "\n")
  cat("Columns in all_results:", paste(colnames(all_results), collapse=", "), "\n")

  
  # Find best parameters for each method
  best_structure_incell <- all_results[which.max(all_results$f1_incell_structure), ]
  best_structure_invitro <- all_results[which.max(all_results$f1_invitro_structure), ]
  best_unuar_structure_incell <- all_results[which.max(all_results$f1_incell_unuar_and_structure), ]
  best_unuar_structure_invitro <- all_results[which.max(all_results$f1_invitro_unuar_and_structure), ]
  best_uguag_structure_incell <- all_results[which.max(all_results$f1_incell_uguag_and_structure), ]
  best_uguag_structure_invitro <- all_results[which.max(all_results$f1_invitro_uguag_and_structure), ]

  # Print sequence F1 scores for debugging
  cat("UNUAR F1 scores - incell:", sequence_f1_scores$UNUAR_incell, 
      "invitro:", sequence_f1_scores$UNUAR_invitro, "\n")
  cat("UGUAG F1 scores - incell:", sequence_f1_scores$UGUAG_incell, 
      "invitro:", sequence_f1_scores$UGUAG_invitro, "\n")

  # Prepare summary of best results
  best_results <- data.frame(
    Method = c("UNUAR", "UGUAG", "Structure", "UNUAR+Structure", "UGUAG+Structure"),
    Incell_F1 = c(sequence_f1_scores$UNUAR_incell, 
                  sequence_f1_scores$UGUAG_incell,
                  max(all_results$f1_incell_structure, na.rm = TRUE),
                  max(all_results$f1_incell_unuar_and_structure, na.rm = TRUE),
                  max(all_results$f1_incell_uguag_and_structure, na.rm = TRUE)),
    Invitro_F1 = c(sequence_f1_scores$UNUAR_invitro, 
                   sequence_f1_scores$UGUAG_invitro,
                   max(all_results$f1_invitro_structure, na.rm = TRUE),
                   max(all_results$f1_invitro_unuar_and_structure, na.rm = TRUE),
                   max(all_results$f1_invitro_uguag_and_structure, na.rm = TRUE))
  )


  # Print and save best results
  cat("\nBest F1 scores for each method:\n")
  print(best_results)
  write.csv(best_results, path(output_dir, "best_f1_scores.csv"), row.names = FALSE)

  # Save best parameters for each method
  write.csv(best_structure_incell, path(output_dir, "best_structure_incell_params.csv"), row.names = FALSE)
  write.csv(best_structure_invitro, path(output_dir, "best_structure_invitro_params.csv"), row.names = FALSE)
  write.csv(best_unuar_structure_incell, path(output_dir, "best_unuar_structure_incell_params.csv"), row.names = FALSE)
  write.csv(best_unuar_structure_invitro, path(output_dir, "best_unuar_structure_invitro_params.csv"), row.names = FALSE)
  write.csv(best_uguag_structure_incell, path(output_dir, "best_uguag_structure_incell_params.csv"), row.names = FALSE)
  write.csv(best_uguag_structure_invitro, path(output_dir, "best_uguag_structure_invitro_params.csv"), row.names = FALSE)

  # Find top 20 parameter sets for incell and invitro for structure
  top20_incell <- all_results[order(all_results$f1_incell_structure, decreasing = TRUE), ][1:20, ]
  top20_invitro <- all_results[order(all_results$f1_invitro_structure, decreasing = TRUE), ][1:20, ]
  
  # Save top 20 parameter sets
  write.csv(top20_incell, path(output_dir, "top20_parameters_incell_structure.csv"), row.names = FALSE)
  write.csv(top20_invitro, path(output_dir, "top20_parameters_invitro_structure.csv"), row.names = FALSE)
  
  return(list(all_results = all_results,
    best_results = best_results,
    best_structure_incell = best_structure_incell,
    best_structure_invitro = best_structure_invitro,
    best_unuar_structure_incell = best_unuar_structure_incell,
    best_unuar_structure_invitro = best_unuar_structure_invitro,
    best_uguag_structure_incell = best_uguag_structure_incell,
    best_uguag_structure_invitro = best_uguag_structure_invitro,
    incell = top20_incell, invitro = top20_invitro))
}


# Run limited combinations
cat("Starting limited parameter sweep. Results will be saved in:", output_dir, "\n")
run_limited_combinations()
cat("Limited parameter sweep completed. Results are saved in:", output_dir, "\n")

# Add this new print statement
cat("\nSummary of output files:\n")
cat("1. All analysis results:", path(output_dir, "all_analysis_results.csv"), "\n")
cat("2. Best F1 scores for each method:", path(output_dir, "best_f1_scores.csv"), "\n")
cat("3. Best parameters for structure (incell):", path(output_dir, "best_structure_incell_params.csv"), "\n")
cat("4. Best parameters for structure (invitro):", path(output_dir, "best_structure_invitro_params.csv"), "\n")
cat("5. Best parameters for UNUAR+structure (incell):", path(output_dir, "best_unuar_structure_incell_params.csv"), "\n")
cat("6. Best parameters for UNUAR+structure (invitro):", path(output_dir, "best_unuar_structure_invitro_params.csv"), "\n")
cat("7. Best parameters for UGUAG+structure (incell):", path(output_dir, "best_uguag_structure_incell_params.csv"), "\n")
cat("8. Best parameters for UGUAG+structure (invitro):", path(output_dir, "best_uguag_structure_invitro_params.csv"), "\n")
cat("9. Top 20 parameters for incell, structure only:", path(output_dir, "top20_parameters_incell.csv"), "\n")
cat("10. Top 20 parameters for invitro, structure only:", path(output_dir, "top20_parameters_invitro.csv"), "\n")
cat("11. Individual analysis results for each parameter combination are also saved in this directory.\n")

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