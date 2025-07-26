#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "https://cloud.r-project.org")
}

# --- Load libraries ---
suppressPackageStartupMessages({
  library(optparse)
})

# --- Define options ---
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "CSV file with sequence IDs"),
  make_option(c("-f", "--fold_dir"), type = "character", help = "Directory containing .fold files"),
  make_option(c("-o", "--output"), type = "character", help = "Output CSV file"),
  make_option("--input_position", type = "integer", default = 59, help = "1-indexed position to start searching"),
  make_option("--offset_min", type = "integer", default = 0),
  make_option("--offset_max", type = "integer", default = 2),
  make_option("--min_unpaired1", type = "integer", default = 0),  # New option
  make_option("--max_unpaired1", type = "integer", default = 5),  # New option
  make_option("--min_paired1", type = "integer", default = 3),
  make_option("--max_paired1", type = "integer", default = 6),
  make_option("--min_unpaired2", type = "integer", default = 3),  # Changed from min_unpaired
  make_option("--max_unpaired2", type = "integer", default = 7),  # Changed from max_unpaired
  make_option("--min_paired2", type = "integer", default = 3),
  make_option("--max_paired2", type = "integer", default = 10),
  make_option("--include_unpaired1", type = "logical", default = TRUE, action = "store_false", help = "Include unpaired1 region"),
  make_option("--include_paired2", type = "logical", default = TRUE, action = "store_false", help = "Include paired2 region")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Assign options ---
input_position <- opt$input_position
offset_range <- opt$offset_min:opt$offset_max
min_unpaired1 <- opt$min_unpaired1  # New
max_unpaired1 <- opt$max_unpaired1  # New
min_paired1 <- opt$min_paired1
max_paired1 <- opt$max_paired1
min_paired2 <- opt$min_paired2
max_paired2 <- opt$max_paired2
min_unpaired2 <- opt$min_unpaired2  # Changed from min_unpaired
max_unpaired2 <- opt$max_unpaired2  # Changed from max_unpaired
include_unpaired1 <- opt$include_unpaired1
include_paired2 <- opt$include_paired2

sequence_names_file <- opt$input
fold_files_dir <- opt$fold_dir
output_csv <- opt$output

# --- Read sequence IDs ---
sequence_df <- read.csv(sequence_names_file, stringsAsFactors = FALSE)
sequence_ids <- sequence_df$id

# --- Get matching .fold files ---
all_fold_files <- list.files(fold_files_dir, pattern = "\\.fold$", full.names = TRUE)
filtered_files <- all_fold_files[sapply(all_fold_files, function(file) {
  any(sapply(sequence_ids, function(id) grepl(id, basename(file), fixed = TRUE)))
})]

# --- Helper: motif matcher ---
parse_dot_bracket <- function(dot_bracket) {
  chars <- strsplit(dot_bracket, "")[[1]]
  if (!all(chars %in% c(".", "(", ")"))) {
    stop("Invalid characters in dot-bracket notation")
  }
  return(chars)
}

is_valid_paired_region <- function(chars, start, end) {
  all(chars[start:end] %in% c("(", ")"))
}

is_valid_unpaired_region <- function(chars, start, end) {
  all(chars[start:end] == ".")
}

generate_result_list <- function(i, unpaired1_len, paired1_len, unpaired2_len, paired2_len, include_unpaired1, include_paired2) {
  result <- list(
    paired1_start = i + if(include_unpaired1) unpaired1_len else 0,
    paired1_end = i + if(include_unpaired1) unpaired1_len else 0 + paired1_len - 1,
    paired1_length = paired1_len,
    unpaired2_start = i + if(include_unpaired1) unpaired1_len else 0 + paired1_len,
    unpaired2_end = i + if(include_unpaired1) unpaired1_len else 0 + paired1_len + unpaired2_len - 1,
    unpaired2_length = unpaired2_len
  )
  
  if (include_unpaired1) {
    result <- c(list(
      unpaired1_start = i,
      unpaired1_end = i + unpaired1_len - 1,
      unpaired1_length = unpaired1_len
    ), result)
  } else {
    result <- c(list(
      unpaired1_start = NA,
      unpaired1_end = NA,
      unpaired1_length = NA
    ), result)
  }
  
  if (include_paired2) {
    result <- c(result, list(
      paired2_start = i + if(include_unpaired1) unpaired1_len else 0 + paired1_len + unpaired2_len,
      paired2_end = i + if(include_unpaired1) unpaired1_len else 0 + paired1_len + unpaired2_len + paired2_len - 1,
      paired2_length = paired2_len
    ))
  } else {
    result <- c(result, list(
      paired2_start = NA,
      paired2_end = NA,
      paired2_length = NA
    ))
  }
  
  return(result)
}

has_motif_at <- function(dot_bracket, i, min_unpaired1, max_unpaired1, min_paired1, max_paired1, 
                         min_unpaired2, max_unpaired2, min_paired2, max_paired2, 
                         include_unpaired1, include_paired2) {
  chars <- parse_dot_bracket(dot_bracket)
  n <- length(chars)
  max_total_len <- (if(include_unpaired1) max_unpaired1 else 0) + max_paired1 + max_unpaired2 + (if(include_paired2) max_paired2 else 0)
  if (i + max_total_len - 1 > n) return(NULL)

  unpaired1_range <- if(include_unpaired1) max_unpaired1:min_unpaired1 else 0:0
  paired2_range <- if(include_paired2) max_paired2:min_paired2 else 0:0

  for (unpaired1_len in unpaired1_range) {
    if (include_unpaired1) {
      if (i + unpaired1_len - 1 > n) next
      if (!is_valid_unpaired_region(chars, i, i + unpaired1_len - 1)) next
    }

    for (paired1_len in max_paired1:min_paired1) {
      paired1_start <- i + (if(include_unpaired1) unpaired1_len else 0)
      paired1_end <- paired1_start + paired1_len - 1
      if (paired1_end > n) next
      if (!is_valid_paired_region(chars, paired1_start, paired1_end)) next

      for (unpaired2_len in max_unpaired2:min_unpaired2) {
        unpaired2_start <- paired1_end + 1
        unpaired2_end <- unpaired2_start + unpaired2_len - 1
        if (unpaired2_end > n) next
        if (!is_valid_unpaired_region(chars, unpaired2_start, unpaired2_end)) next

        for (paired2_len in paired2_range) {
          if (include_paired2) {
            paired2_start <- unpaired2_end + 1
            paired2_end <- paired2_start + paired2_len - 1
            if (paired2_end > n) next
            if (!is_valid_paired_region(chars, paired2_start, paired2_end)) next
          }
          
          return(generate_result_list(i, unpaired1_len, paired1_len, unpaired2_len, paired2_len, include_unpaired1, include_paired2))
        }
      }
    }
  }
  return(NULL)
}

has_motif <- function(dot_bracket, base_position, offsets,
                      min_unpaired1, max_unpaired1, min_paired1, max_paired1, 
                      min_unpaired2, max_unpaired2, min_paired2, max_paired2,
                      include_unpaired1, include_paired2) {
  for (offset in offsets) {
    i <- base_position + offset
    result <- has_motif_at(dot_bracket, i, min_unpaired1, max_unpaired1, min_paired1, max_paired1, 
                           min_unpaired2, max_unpaired2, min_paired2, max_paired2,
                           include_unpaired1, include_paired2)
    if (!is.null(result)) {
      return(result)
    }
  }
  return(NULL)
}

# --- Analyze .fold files ---
results <- data.frame()

for (file in filtered_files) {
  lines <- readLines(file, warn = FALSE)
  if (length(lines) >= 3) {
    fasta_seq <- lines[2]
    dot_bracket <- strsplit(lines[3], " ")[[1]][1]

    match_info <- has_motif(
      dot_bracket,
      input_position,
      offset_range,
      min_unpaired1, max_unpaired1,
      min_paired1, max_paired1,
      min_unpaired2, max_unpaired2,
      min_paired2, max_paired2,
      include_unpaired1, include_paired2
    )

    if (!is.null(match_info)) {
      results <- rbind(results, data.frame(
        filename = basename(file),
        unpaired1_start = match_info$unpaired1_start,
        unpaired1_end = match_info$unpaired1_end,
        unpaired1_length = match_info$unpaired1_length,
        paired1_start = match_info$paired1_start,
        paired1_end = match_info$paired1_end,
        paired1_length = match_info$paired1_length,
        unpaired2_start = match_info$unpaired2_start,
        unpaired2_end = match_info$unpaired2_end,
        unpaired2_length = match_info$unpaired2_length,
        paired2_start = match_info$paired2_start,
        paired2_end = match_info$paired2_end,
        paired2_length = match_info$paired2_length,
        fasta_sequence = fasta_seq,
        dot_bracket = dot_bracket,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# --- Save results ---
write.csv(results, output_csv, row.names = FALSE)