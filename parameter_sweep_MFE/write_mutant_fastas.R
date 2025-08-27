#!/usr/bin/env Rscript

# --- Argument parsing ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript write_mutant_fastas.R <mutant_csv> <output_dir>")
}
mutant_csv <- args[1]
output_dir <- args[2]

# --- Load data ---
data <- read.csv(mutant_csv, stringsAsFactors = FALSE)

# --- Ensure output directory exists ---
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Extract base IDs from filename column ---
if (!"filename" %in% names(data)) {
  stop("Input CSV must contain a 'filename' column")
}
base_ids <- sub("\\.fold$", "", data$filename)

# --- Determine available mutant types ---
mutant_types <- c("paired1_disruption", "paired1_compensatory")
if ("paired2_disruption" %in% names(data) && any(!is.na(data$paired2_disruption))) {
  mutant_types <- c(mutant_types, "paired2_disruption", "paired2_compensatory")
}

if ("combined_disruption" %in% names(data) && any(!is.na(data$combined_disruption))) {
  mutant_types <- c(mutant_types, "combined_disruption", "combined_compensatory")
}

# --- Write one FASTA file per mutant type ---
for (type in mutant_types) {
  fasta_path <- file.path(output_dir, paste0(type, ".fa"))
  con <- file(fasta_path, "w")

  for (i in seq_len(nrow(data))) {
    seq <- data[[type]][i]
    if (!is.na(seq) && nzchar(seq)) {
      writeLines(c(paste0(">", base_ids[i]), seq), con)
    }
  }

  close(con)
}

cat("Created FASTA files for mutant types:", paste(mutant_types, collapse = ", "), "\n")