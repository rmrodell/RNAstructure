#!/usr/bin/env Rscript

library(data.table)  # For fast reading and binding of data

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check that enough arguments are provided
if (length(args) != 2) {
  stop("Please provide a comma-separated list of directories and an output file name.")
}

# Parse command line arguments
input_directories <- unlist(strsplit(args[1], ","))  # Split the directories by comma
output_file <- args[2]  # Desired output file name

# Initialize an empty list to store data frames
all_data_frames <- list()

# Loop through each directory
for (input_directory in input_directories) {
  # Get a list of all files that end with profile.txt in the directory
  file_list <- list.files(path = input_directory, pattern = "profile\\.txt$", full.names = TRUE)

  # Loop through each file found
  for (file_path in file_list) {
    # Read the tab-separated file into a data frame
    df <- fread(file_path, sep = "\t", header = TRUE)

    # Add a new column with the filename (without path) as the first column
    df_file_name <- basename(file_path)
    df$file_name <- df_file_name  # Add the filename as a new column

    # Append the data frame to the list
    all_data_frames[[length(all_data_frames) + 1]] <- df
  }
}

# Combine all data frames into one
# If there are data frames in the list
if (length(all_data_frames) > 0) {
  combined_data <- rbindlist(all_data_frames, use.names = TRUE, fill = TRUE)

  # Write the combined data to the output file
  fwrite(combined_data, output_file, sep = "\t", quote = TRUE, row.names = FALSE)
  cat("Data successfully concatenated and written to", output_file, "\n")
} else {
  cat("No files found matching the pattern.\n")
}