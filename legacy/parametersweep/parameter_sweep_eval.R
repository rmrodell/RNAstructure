#!/usr/bin/env Rscript

library(fs)
library(optparse)
library(dplyr)
library(stringr)

# Capture warnings
warning_messages <- character()
withCallingHandlers(
  {

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

# Function to analyze dataframe
analyze_dataframe <- function(df, pool1_df) {
  filtered_df <- pool1_df %>%
    filter(sapply(chr, function(x) any(str_detect(df$filename, fixed(x)))))
  
  total_incell <- sum(pool1_df$incell == 1)
  total_invitro <- sum(pool1_df$invitro == 1)
  
  results <- filtered_df %>%
    summarise(
      total_matches = n(),
      incell_count = sum(incell == 1),
      invitro_count = sum(invitro == 1),
      incell_perc_overmotif = (sum(incell == 1) / total_matches) * 100,
      invitro_perc_overmotif = (sum(invitro == 1) / total_matches) * 100,
      incell_perc_overmodified = (sum(incell == 1) / total_incell) * 100,
      invitro_perc_overmodified = (sum(invitro == 1) / total_invitro) * 100
    )
  
  results$f1_incell <- calculate_f1_score(pool1_df$incell, 
                                          ifelse(filtered_df$incell == 1, 1, 0))
  results$f1_invitro <- calculate_f1_score(pool1_df$invitro, 
                                           ifelse(filtered_df$invitro == 1, 1, 0))
  
  return(results)
}

# Function to calculate total combinations
calculate_total_combinations <- function(include_paired2) {
  unpaired1_combinations <- length(0:3) * length(2:7)
  paired1_combinations <- length(0:3) * length(2:7)
  unpaired2_combinations <- length(0:3) * length(2:7)
  
  total_without_paired2 <- unpaired1_combinations * paired1_combinations * unpaired2_combinations
  
  if (include_paired2) {
    paired2_combinations <- length(0:3) * length(2:7)
    total_with_paired2 <- total_without_paired2 * paired2_combinations
    return(total_without_paired2 + total_with_paired2)
  } else {
    return(total_without_paired2)
  }
}

# Calculate and print the total number of combinations
total_combinations <- calculate_total_combinations(TRUE)
cat("Total number of combinations to be run:", total_combinations, "\n")

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
  
  # Analyze the results
  df <- read.csv(output_file)
  analysis_results <- analyze_dataframe(df, pool1_df)
  
  # Ensure all expected columns are present
  expected_columns <- c("total_matches", "incell_count", "invitro_count", 
                        "incell_perc_overmotif", "invitro_perc_overmotif", 
                        "incell_perc_overmodified", "invitro_perc_overmodified",
                        "f1_incell", "f1_invitro")
  
  for (col in expected_columns) {
    if (!(col %in% names(analysis_results))) {
      analysis_results[[col]] <- NA
    }
  }
  
  # Add the parameters to the results, ensuring they're always present
  analysis_results$unpaired1_min <- unpaired1_min
  analysis_results$unpaired1_max <- unpaired1_max
  analysis_results$paired1_min <- paired1_min
  analysis_results$paired1_max <- paired1_max
  analysis_results$unpaired2_min <- unpaired2_min
  analysis_results$unpaired2_max <- unpaired2_max
  analysis_results$paired2_min <- ifelse(include_paired2, paired2_min, NA)
  analysis_results$paired2_max <- ifelse(include_paired2, paired2_max, NA)
  analysis_results$include_paired2 <- include_paired2
  
  # Ensure consistent column order
  all_columns <- c(expected_columns, 
                   "unpaired1_min", "unpaired1_max", 
                   "paired1_min", "paired1_max", 
                   "unpaired2_min", "unpaired2_max", 
                   "paired2_min", "paired2_max", 
                   "include_paired2")
  
  analysis_results <- analysis_results[, all_columns]
  
  cat("Completed run:", run_id, "\n")
  return(analysis_results)
}

# Function to run limited combinations
run_limited_combinations <- function() {
  all_results <- data.frame()
  
  for (unpaired1_min in 0:3) {
    for (unpaired1_max in 2:7) {
      for (paired1_min in 0:3) {
        for (paired1_max in 2:7) {
          for (unpaired2_min in 0:3) {
            for (unpaired2_max in 2:7) {
              # Test without paired2
              result <- run_analysis(unpaired1_min, unpaired1_max, 
                                     paired1_min, paired1_max, 
                                     unpaired2_min, unpaired2_max, 
                                     include_paired2 = FALSE)
              all_results <- rbind(all_results, result)
              
              # Test with paired2
              for (paired2_min in 0:3) {
                for (paired2_max in 2:7) {
                  result <- run_analysis(unpaired1_min, unpaired1_max, 
                                         paired1_min, paired1_max, 
                                         unpaired2_min, unpaired2_max, 
                                         paired2_min, paired2_max, 
                                         include_paired2 = TRUE)
                  all_results <- rbind(all_results, result)
                }
              }
            }
          }
        }
      }
    }
  }

  # Save all results
  write.csv(all_results, path(output_dir, "all_analysis_results.csv"), row.names = FALSE)
  
  # Find top 20 parameter sets for incell and invitro
  top20_incell <- all_results[order(all_results$f1_incell, decreasing = TRUE), ][1:20, ]
  top20_invitro <- all_results[order(all_results$f1_invitro, decreasing = TRUE), ][1:20, ]
  
  cat("Top 20 parameters for incell:\n")
  print(top20_incell)
  cat("\nTop 20 parameters for invitro:\n")
  print(top20_invitro)
  
  # Save top 20 parameter sets
  write.csv(top20_incell, path(output_dir, "top20_parameters_incell.csv"), row.names = FALSE)
  write.csv(top20_invitro, path(output_dir, "top20_parameters_invitro.csv"), row.names = FALSE)
  
  return(list(incell = top20_incell, invitro = top20_invitro))
}


# Run limited combinations
cat("Starting limited parameter sweep. Results will be saved in:", output_dir, "\n")
run_limited_combinations()
cat("Limited parameter sweep completed. Results are saved in:", output_dir, "\n")

# Add this new print statement
cat("\nSummary of output files:\n")
cat("1. All analysis results:", path(output_dir, "all_analysis_results.csv"), "\n")
cat("2. Top 20 parameters for incell:", path(output_dir, "top20_parameters_incell.csv"), "\n")
cat("3. Top 20 parameters for invitro:", path(output_dir, "top20_parameters_invitro.csv"), "\n")
cat("4. Individual analysis results for each parameter combination are also saved in this directory.\n")

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