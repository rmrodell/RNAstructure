#!/usr/bin/env Rscript

library(fs)

library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type="character", default="input_pool1.csv",
              help="Input CSV file [default= %default]"),
  make_option(c("-f", "--fold_dir"), type="character", default="/scratch/users/rodell/RNAfold_psipos",
              help="Directory containing .fold files [default= %default]"),
  make_option(c("-o", "--output_dir"), type="character", default="testsweep_output",
              help="Output directory for results [default= %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Use the parsed options
input_file <- opt$input
fold_dir <- opt$fold_dir
output_dir <- opt$output_dir


# Create output directory:
dir_create(output_dir)


# Function to calculate total combinations
calculate_total_combinations <- function(include_paired2) {
  unpaired1_combinations <- length(0:1) * length(2:3)
  paired1_combinations <- length(0:1) * length(2:3)
  unpaired2_combinations <- length(0:1) * length(2:3)
  
  total_without_paired2 <- unpaired1_combinations * paired1_combinations * unpaired2_combinations
  
  if (include_paired2) {
    paired2_combinations <- length(0:1) * length(2:3)
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
    # Add any other necessary arguments for motif_matcher_v2.R here
    "--input", input_file,
    "--fold_dir", fold_dir
  )
  
  # Run the command and capture the output
  result <- system2("Rscript", args = cmd_args, stdout = TRUE, stderr = TRUE)
  
  # Log the output
  writeLines(result, path(output_dir, paste0(run_id, "_log.txt")))
  
  cat("Completed run:", run_id, "\n")
}

# Function to run limited combinations
run_limited_combinations <- function() {
  for (unpaired1_min in 0:1) {
    for (unpaired1_max in 2:3) {
      for (paired1_min in 0:1) {
        for (paired1_max in 2:3) {
          for (unpaired2_min in 0:1) {
            for (unpaired2_max in 2:3) {
              # Test without paired2
              run_analysis(unpaired1_min, unpaired1_max, 
                           paired1_min, paired1_max, 
                           unpaired2_min, unpaired2_max, 
                           include_paired2 = FALSE)
              
              # Test with paired2
              for (paired2_min in 0:1) {
                for (paired2_max in 2:3) {
                  run_analysis(unpaired1_min, unpaired1_max, 
                               paired1_min, paired1_max, 
                               unpaired2_min, unpaired2_max, 
                               paired2_min, paired2_max, 
                               include_paired2 = TRUE)
                }
              }
            }
          }
        }
      }
    }
  }
}

# Run limited combinations
cat("Starting limited parameter sweep. Results will be saved in:", output_dir, "\n")
run_limited_combinations()
cat("Limited parameter sweep completed. Results are saved in:", output_dir, "\n")


# identify best outputs

calculate_f1_score <- function(actual, predicted) {
  true_positives <- sum(actual == 1 & predicted == 1)
  false_positives <- sum(actual == 0 & predicted == 1)
  false_negatives <- sum(actual == 1 & predicted == 0)
  
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  
  f1_score <- 2 * (precision * recall) / (precision + recall)
  return(f1_score)
}