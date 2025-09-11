#!/usr/bin/env Rscript

library(data.table)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values
default_num_columns <- 130
default_output_file <- "SHAPE_reactivities.csv"

# Check for sufficient arguments
if (length(args) < 3) {
  stop("Please provide directories (comma-separated), number of columns, and output CSV filename.")
}

# Parse command line arguments
input_directories <- unlist(strsplit(args[1], ","))  # Split by comma and convert to a vector
num_columns <- as.integer(args[2])  # Convert to integer
output_file <- args[3]  # File name for output

print("Started")
print(paste("Input Directories:", input_directories))
print(paste("Number of Columns:", num_columns))
print(paste("Output CSV Name:", output_file))


# Initialize the master data frame
master_data_frame <- data.frame(matrix(0, nrow = 0, ncol = num_columns))
colnames(master_data_frame) <- as.character(1:num_columns)

number_iterations <- 0 

for (input_directory in input_directories) {
  # Ensure we get a list of .shape files from each specified directory
  file_list <- list.files(path = input_directory, pattern = "\\.shape$", full.names = TRUE)
  
  for (file_path in file_list) {
    number_iterations <- number_iterations + 1
    print(paste("Processing file number:", number_iterations))
    row_name <- basename(file_path)  # Set row name to be the file name
    print(paste("Processing file:", row_name))
    
    # Read the content of the .shape file
    lines <- readLines(file_path)

    # Initialize a vector to store values for the current row based on columns
    current_values <- numeric(num_columns)

    # Process each line in the file
    for (line in lines) {
      # Split the line by tab
      parts <- unlist(strsplit(line, "\t"))  # Explicitly using tab as the delimiter 
      
      if (length(parts) == 2) {
        col1 <- as.integer(parts[1])  # First column
        value <- as.numeric(parts[2])  # Second column
        
        if (!is.na(value) && col1 <= num_columns && col1 > 0) {
          # Store the value in the corresponding position
          current_values[col1] <- current_values[col1] + value  # Direct accumulation
        } else {
          print(paste("Warning: Invalid data in file:", basename(file_path), "Column 1:", col1, "Value:", value))
        }
      } else {
        print(paste("Warning: Unexpected format in line:", line))
      }
    }

    # Add current values to the master data frame
    master_data_frame[row_name, ] <- current_values  # Directly set the values for the row
    
    # Ensure the row is added to the master data frame if it was not present (avoid NA issues)
    if (!(row_name %in% rownames(master_data_frame))) {
      master_data_frame[row_name, ] <- NA  # Adding a new row and filling it with NA initially
    }
  }
}

# Writing the master data frame to a CSV file
write.csv(master_data_frame, output_file, row.names = TRUE)
print("Success!")