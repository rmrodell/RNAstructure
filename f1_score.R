#!/usr/bin/env Rscript

# Load required libraries
library(readr)

# Define the F1 score calculation function
calculate_f1_score <- function(actual, predicted) {
  true_positives <- sum(actual == 1 & predicted == 1)
  false_positives <- sum(actual == 0 & predicted == 1)
  false_negatives <- sum(actual == 1 & predicted == 0)
  
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  
  f1_score <- 2 * (precision * recall) / (precision + recall)
  return(f1_score)
}

# Main function to calculate reference F1 scores
calculate_reference_f1_scores <- function(file_path) {
  # Read the CSV file
  data <- read_csv(file_path)
  
  # Calculate F1 scores
  unuar_incell_f1 <- calculate_f1_score(data$incell, data$UNUAR)
  unuar_invitro_f1 <- calculate_f1_score(data$invitro, data$UNUAR)
  uguag_incell_f1 <- calculate_f1_score(data$incell, data$UGUAG)
  uguag_invitro_f1 <- calculate_f1_score(data$invitro, data$UGUAG)
  
  # Print results
  cat("Reference F1 Scores:\n")
  cat("UNUAR for incell:", unuar_incell_f1, "\n")
  cat("UNUAR for invitro:", unuar_invitro_f1, "\n")
  cat("UGUAG for incell:", uguag_incell_f1, "\n")
  cat("UGUAG for invitro:", uguag_invitro_f1, "\n")
  
  # Return results as a list
  return(list(
    UNUAR_incell = unuar_incell_f1,
    UNUAR_invitro = unuar_invitro_f1,
    UGUAG_incell = uguag_incell_f1,
    UGUAG_invitro = uguag_invitro_f1
  ))
}

# Check if the script is being run directly
if (!interactive()) {
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) != 1) {
    stop("Usage: Rscript calculate_reference_f1.R <path_to_csv_file>")
  }
  
  file_path <- args[1]
  
  # Calculate and print F1 scores
  results <- calculate_reference_f1_scores(file_path)
  
  # Save results to CSV
  write_csv(results$results_df, "reference_f1_scores.csv")
  cat("Results saved to reference_f1_scores.csv\n")

  # Optionally, save results to a file
  saveRDS(results, file = "reference_f1_scores.rds")
}