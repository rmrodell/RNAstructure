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
  make_option("--offset_max", type = "integer", default = 1),
  make_option("--min_unpaired1", type = "integer", default = 0),
  make_option("--max_unpaired1", type = "integer", default = 1),
  make_option("--min_paired1", type = "integer", default = 2),
  make_option("--max_paired1", type = "integer", default = 6),
  make_option("--min_unpaired2", type = "integer", default = 2),
  make_option("--max_unpaired2", type = "integer", default = 6),
  make_option("--include_paired2", type = "logical", default = TRUE, 
            action = "store_false", 
            help = "Include paired2 region in the search (default: TRUE)"),
  make_option("--min_paired2", type = "integer", default = 3),
  make_option("--max_paired2", type = "integer", default = 10)
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Assign options ---
input_position <- opt$input_position
offset_range <- opt$offset_min:opt$offset_max
min_paired1 <- opt$min_paired1
max_paired1 <- opt$max_paired1
include_paired2 <- opt$include_paired2
min_paired2 <- opt$min_paired2
max_paired2 <- opt$max_paired2
min_unpaired1 <- opt$min_unpaired1
max_unpaired1 <- opt$max_unpaired1
min_unpaired2 <- opt$min_unpaired2
max_unpaired2 <- opt$max_unpaired2


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
has_motif_at <- function(dot_bracket, i, min_u1, max_u1, min_p1, max_p1, min_u2, max_u2, min_p2, max_p2, include_p2) {
  chars <- strsplit(dot_bracket, "")[[1]]
  n <- length(chars)
  max_total_len <- max_u1 + max_p1 + max_u2 + (if(include_p2) max_p2 else 0)
  if (i + max_total_len - 1 > n) return(NULL)

  for (u1_len in max_u1:min_u1) {
    u1_start <- i
    u1_end <- u1_start + u1_len - 1
    if (u1_end > n) next
    if (!all(chars[u1_start:u1_end] == ".")) next

    for (p1_len in max_p1:min_p1) {
      p1_start <- u1_end + 1
      p1_end <- p1_start + p1_len - 1
      if (p1_end > n) next
      if (!all(chars[p1_start:p1_end] %in% c("(", ")"))) next

      for (u2_len in max_u2:min_u2) {
        u2_start <- p1_end + 1
        u2_end <- u2_start + u2_len - 1
        if (u2_end > n) next
        if (!all(chars[u2_start:u2_end] == ".")) next

        # if paired2 is NOT included
        if (!include_p2) {
          return(list(
            unpaired1_start = u1_start,
            unpaired1_end = u1_end,
            unpaired1_length = u1_len,
            paired1_start = p1_start,
            paired1_end = p1_end,
            paired1_length = p1_len,
            unpaired2_start = u2_start,
            unpaired2_end = u2_end,
            unpaired2_length = u2_len
          ))
        }

        # if paired2 is included
        for (p2_len in max_p2:min_p2) {
          p2_start <- u2_end + 1
          p2_end <- p2_start + p2_len - 1
          if (p2_end > n) next
          if (all(chars[p2_start:p2_end] %in% c("(", ")"))) {
            return(list(
              unpaired1_start = u1_start,
              unpaired1_end = u1_end,
              unpaired1_length = u1_len,
              paired1_start = p1_start,
              paired1_end = p1_end,
              paired1_length = p1_len,
              unpaired2_start = u2_start,
              unpaired2_end = u2_end,
              unpaired2_length = u2_len,
              paired2_start = p2_start,
              paired2_end = p2_end,
              paired2_length = p2_len
            ))
          }
        }
      }
    }
  }
  return(NULL)
}

has_motif <- function(dot_bracket, base_position, offsets,
                      min_u1, max_u1, min_p1, max_p1, min_u2, max_u2, min_p2, max_p2, include_p2) {
  for (offset in offsets) {
    i <- base_position + offset
    result <- has_motif_at(dot_bracket, i, min_u1, max_u1, min_p1, max_p1, min_u2, max_u2, min_p2, max_p2, include_p2)
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
      include_paired2
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
        paired2_start = if(include_paired2) match_info$paired2_start else NA,
        paired2_end = if(include_paired2) match_info$paired2_end else NA,
        paired2_length = if(include_paired2) match_info$paired2_length else NA,
        fasta_sequence = fasta_seq,
        dot_bracket = dot_bracket,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# --- Save results ---
write.csv(results, output_csv, row.names = FALSE)
