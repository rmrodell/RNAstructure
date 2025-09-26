#!/usr/bin/env Rscript

library(data.table)

print("Started")

wdir = getwd()

input_directory <- paste(wdir)

file_list <- list.files(path = input_directory, pattern = "dp\\.ps$", full.names = TRUE)


num_columns <- 150

# Initialize the master data frame
master_data_frame <- data.frame(matrix(0, nrow = length(file_list), ncol = num_columns))
colnames(master_data_frame) <- as.character(1:num_columns)
rownames(master_data_frame) <- basename(file_list)  # Set row names to the file names


# This function extracts the first three numbers from a line (which is what we need for extracting relevant info from the dp.ps files)
extract_numbers <- function(line) {
  # Split the line into parts by whitespace
  parts <- strsplit(line, "\\s+")[[1]]
  
  if (length(parts) == 4) {
    # Extract and return the first three numbers
    return(as.numeric(parts[1:3]))
  } else {
    return(NULL)
  }
}

number_iterations <- 0 
for (file_path in file_list) {
  number_iterations <- number_iterations + 1
  print(paste("Number of files examined:", number_iterations))
  row_name <- basename(file_path)

  # Read the content of the .ps file into a vector of lines
  lines <- readLines(file_path)
  
  # Filter lines that end with 'lbox'
  lbox_lines <- grep("lbox$", lines, value = TRUE)
  
  for (line in lbox_lines) {
    numbers <- extract_numbers(line)
    
    if (!is.null(numbers)) {
      col1 <- as.integer(numbers[1])
      col2 <- as.integer(numbers[2])
      value <- numbers[3]
      
      
      if (col1 <= num_columns && col1 > 0) {
        master_data_frame[row_name, as.character(col1)] <- master_data_frame[row_name, as.character(col1)] + (value^2)
      }
      if (col2 <= num_columns && col2 > 0) {
        master_data_frame[row_name, as.character(col2)] <- master_data_frame[row_name, as.character(col2)] + (value^2)
      }
    }
  }
}

write.csv(master_data_frame, paste(wdir, "/", "MFE_pairingprob.csv", sep=""),row.names = T)
print("Success!")

