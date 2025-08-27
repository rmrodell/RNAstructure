# --- Configuration ---
# Set this to the parent directory containing all your result folders.
# "." means it will execute in the working directory
base_dir <- "."

# --- Main Processing Logic ---

# Get a list of all immediate subdirectories
sub_dirs <- list.dirs(path = base_dir, recursive = FALSE, full.names = TRUE)

cat("Scanning", length(sub_dirs), "directories...\n\n")

# Loop through each subdirectory
for (dir_path in sub_dirs) {
  checkpoint_file <- file.path(dir_path, "checkpoint.rds")
  final_file <- file.path(dir_path, "all_analysis_results.csv")
  
  # Process only if checkpoint exists AND the final file does NOT exist
  if (file.exists(checkpoint_file) && !file.exists(final_file)) {
    
    cat("Processing:", basename(dir_path), "...\n")
    
    # Load the list of results
    results_list <- readRDS(checkpoint_file)
    
    # Filter out NULLs (failed or incomplete jobs)
    valid_results <- results_list[!sapply(results_list, is.null)]
    
    # Proceed only if there are valid results to combine
    if (length(valid_results) > 0) {
      # Combine the list of data frames into one
      all_results <- do.call(rbind, valid_results)
      
      # Write the final CSV file
      write.csv(all_results, final_file, row.names = FALSE)
      cat(" -> Created all_analysis_results.csv\n")
    } else {
      cat(" -> Checkpoint found, but contained no valid results. Skipping.\n")
    }
  }
}

cat("\nScan complete. All applicable directories have been processed.\n")