#!/usr/bin/env Rscript

library(fs)

# Create a temporary directory for outputs
temp_dir <- tempdir()
output_dir <- path(temp_dir, "motif_search_results")
dir_create(output_dir)

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
  
  # Construct the command to run your main script
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
    "--output", output_file
  )
  
  # Run the command and capture the output
  result <- system2("Rscript", args = cmd_args, stdout = TRUE, stderr = TRUE)
  
  # Optionally, you can log the output
  writeLines(result, path(output_dir, paste0(run_id, "_log.txt")))
  
  cat("Completed run:", run_id, "\n")
}

# Function to run all combinations
run_all_combinations <- function() {
  for (unpaired1_min in 0:3) {
    for (unpaired1_max in max(2, unpaired1_min):7) {
      for (paired1_min in 0:3) {
        for (paired1_max in max(2, paired1_min):7) {
          for (unpaired2_min in 0:3) {
            for (unpaired2_max in max(2, unpaired2_min):7) {
              # Test without paired2
              run_analysis(unpaired1_min, unpaired1_max, 
                           paired1_min, paired1_max, 
                           unpaired2_min, unpaired2_max, 
                           include_paired2 = FALSE)
              
              # Test with paired2
              for (paired2_min in 0:3) {
                for (paired2_max in max(2, paired2_min):7) {
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

# Run all combinations
cat("Starting parameter sweep. Results will be saved in:", output_dir, "\n")
run_all_combinations()
cat("Parameter sweep completed. Results are saved in:", output_dir, "\n")