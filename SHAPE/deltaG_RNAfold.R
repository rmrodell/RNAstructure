#!/usr/bin/env Rscript

library(data.table)
library(stringr)
library(optparse)

# this extracts all numerical values from a folder of .fold files from RNAfold and returns a csv file

# Parse command-line arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input directory containing .fold files", metavar="DIRECTORY"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output directory for the CSV files", metavar="DIRECTORY")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Input and output directories must be specified.", call.=FALSE)
}

print("Started")

input_directory <- opt$input
output_directory <- opt$output

# Ensure output directory exists
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

file_list <- list.files(path = input_directory, pattern = "\\.fold$", full.names = TRUE)

# Define the columns we want to extract
columns <- c("sequence_id", "sequence", "mfe_structure", "mfe_energy", 
             "centroid_energy", "ensemble_energy", "ensemble_distance",
             "mea_energy", "mea_score", "mfe_frequency", "ensemble_diversity")

# Initialize the master data frame and failed files data frame
master_data_frame <- data.frame(matrix(ncol = length(columns), nrow = 0))
colnames(master_data_frame) <- columns
failed_files_df <- data.frame(file_name = character(), error_message = character(), stringsAsFactors = FALSE)

extract_fold_info <- function(file_path) {
  tryCatch({
    lines <- readLines(file_path)
    
    if (length(lines) < 7) {
      stop("File does not contain enough lines")
    }
    
    sequence_id <- sub(">", "", lines[1])
    sequence <- lines[2]
    mfe_structure <- str_extract(lines[3], "^[().]+")
    
    # Updated all number extractions to handle potential spaces
    mfe_energy <- as.numeric(str_extract(lines[3], "(?<=\\()\\s*[-]?[0-9.]+(?=\\s*\\)$)"))
    centroid_energy <- as.numeric(str_extract(lines[4], "(?<=\\[)\\s*[-]?[0-9.]+(?=\\s*\\]$)"))
    ensemble_energy <- as.numeric(str_extract(lines[5], "(?<=\\{)\\s*[-]?[0-9.]+(?=\\s+d=)"))
    ensemble_distance <- as.numeric(str_extract(lines[5], "(?<=d=)\\s*[0-9.]+(?=\\s*\\})"))
    mea_energy <- as.numeric(str_extract(lines[6], "(?<=\\{)\\s*[-]?[0-9.]+(?=\\s+MEA=)"))
    mea_score <- as.numeric(str_extract(lines[6], "(?<=MEA=)\\s*[0-9.]+(?=\\s*\\})"))
    
    stats <- str_match(lines[7], "frequency of mfe structure in ensemble\\s*([0-9.]+)\\s*;\\s*ensemble diversity\\s*([0-9.]+)")
    mfe_frequency <- as.numeric(stats[2])
    ensemble_diversity <- as.numeric(stats[3])
    
    if (any(is.na(c(mfe_energy, centroid_energy, ensemble_energy, ensemble_distance, mea_energy, mea_score, mfe_frequency, ensemble_diversity)))) {
      stop("Failed to extract one or more numeric values")
    }
    
    result <- c(sequence_id, sequence, mfe_structure, mfe_energy, centroid_energy,
                ensemble_energy, ensemble_distance, mea_energy, mea_score,
                mfe_frequency, ensemble_diversity)
    names(result) <- columns  # Use the 'columns' vector defined earlier
    return(result)
  }, error = function(e) {
    return(list(error = TRUE, message = e$message))
  })
}

number_iterations <- 0
for (file_path in file_list) {
  number_iterations <- number_iterations + 1
  print(paste("Processing file", number_iterations, "of", length(file_list), ":", basename(file_path)))
  
  fold_info <- extract_fold_info(file_path)
  if (!is.list(fold_info) || !fold_info$error) {
    master_data_frame <- do.call(rbind, list(master_data_frame, fold_info))
  } else {
    failed_files_df <- rbind(failed_files_df, data.frame(file_name = basename(file_path), error_message = fold_info$message))
  }
}

colnames(master_data_frame) <- columns

# Write successful extractions to CSV
output_file <- file.path(output_directory, "RNAfold_deltaG.csv")
write.csv(master_data_frame, output_file, row.names = FALSE)
print(paste("Success! Output written to", output_file))

# Write failed files to a separate CSV
failed_files_output <- file.path(output_directory, "deltaG_failed_files.csv")
write.csv(failed_files_df, failed_files_output, row.names = FALSE)
print(paste("Failed files information written to", failed_files_output))

# Print summary
print(paste("Total files processed:", length(file_list)))
print(paste("Successfully processed files:", nrow(master_data_frame)))
print(paste("Failed files:", nrow(failed_files_df)))