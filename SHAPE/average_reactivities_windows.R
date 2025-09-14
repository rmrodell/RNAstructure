#!/usr/bin/env Rscript

# =========================================================================
# SETUP: Load Libraries and Define Paths
# =========================================================================
library(data.table)
library(dplyr)
library(tidyr)
library(Biostrings)

# --- Input & Output Paths ---
rep1_file <- "/scratch/groups/nicolemm/rodell/SHAPE/InVitro_Rep1_20241209/shapemapper_fullpool_highqual/Rep1_NoPUS_profile_annotated.txt"
rep2_file <- "/scratch/groups/nicolemm/rodell/SHAPE/InVitro_Rep2_20250731/shapemapper_fullpool_highqual/Rep2_NoPUS_profile_annotated.txt"

position_info_file <- "/scratch/groups/nicolemm/rodell/SHAPE/InVitro_2repHighQ_Analysis/deltaG_windows/Pool1_DiffLengths.csv"

base_output_dir  <- "/scratch/groups/nicolemm/rodell/SHAPE/InVitro_2repHighQ_Analysis/deltaG_windows"
shape_output_dir <- file.path(base_output_dir, "windowed_averaged_shape_files")
fasta_output_path <- file.path(base_output_dir, "pool1psipos_windows.fasta")

dir.create(shape_output_dir, showWarnings = FALSE, recursive = TRUE)
cat("Starting SHAPE profile averaging and windowing...\n")

# =========================================================================
# STEP 1: Load and Combine Data
# =========================================================================
rep1_dt <- fread(rep1_file)
rep2_dt <- fread(rep2_file)

combined_dt <- rbindlist(list(rep1 = rep1_dt, rep2 = rep2_dt), idcol = "replicate")
cat("Data loaded and combined successfully.\n")

# =========================================================================
# STEP 2: Filter Data for High-Quality Values
# =========================================================================
initial_rows <- nrow(combined_dt)

filtered_dt <- combined_dt[
  QC_flag == "good" & 
    is.finite(Norm_profile) & 
    Norm_profile != -999
]

num_removed <- initial_rows - nrow(filtered_dt)
cat(paste("Filtering: Removed", num_removed, "rows due to 'poor' QC, NA, or -999 values.\n"))

# =========================================================================
# STEP 3: Calculate Inverse-Variance Weighted Average Reactivity
# =========================================================================
averaged_reactivities <- filtered_dt[, .(
  avg_reactivity = {
    if (.N == 2) {
      weights <- 1 / (Norm_stderr^2)
      if (any(!is.finite(weights))) {
        weights[!is.finite(weights)] <- 1e9 
      }
      sum(Norm_profile * weights) / sum(weights)
    } else {
      Norm_profile
    }
  }
), by = .(RNA_name, Nucleotide)]

cat(paste("Calculated weighted averages for", nrow(averaged_reactivities), "nucleotide positions.\n"))

# =========================================================================
# STEP 4: Prepare Full-Length .shape Files (Handle Missing Nucleotides)
# =========================================================================
rna_lengths <- combined_dt[, .(max_nt = max(Nucleotide)), by = RNA_name]
full_template <- rna_lengths[, .(Nucleotide = 1:max_nt), by = RNA_name]

final_shapes_dt <- merge(full_template, averaged_reactivities, by = c("RNA_name", "Nucleotide"), all.x = TRUE)
final_shapes_dt[is.na(avg_reactivity), avg_reactivity := -999]

cat("Generated full-length profiles and filled missing nucleotides with -999.\n")


# =========================================================================
# STEP 5: Prepare Window Definitions and Join with SHAPE data
# =========================================================================

position_info_all <- read.csv(position_info_file)

# --- Sanitize names to ensure they match between the two datasets ---
position_info_all$RNA_name <- gsub(",", "_", position_info_all$name)

nested_shapes <- final_shapes_dt %>%
  group_by(RNA_name) %>%
  nest(shape_data = c(Nucleotide, avg_reactivity))

window_definitions <- position_info_all %>%
  as_tibble() %>%
  left_join(nested_shapes, by = "RNA_name") %>%
  filter(!sapply(shape_data, is.null)) # Remove rows with no matching SHAPE data

cat("Joined sequence definitions with SHAPE data.\n")

# =========================================================================
# STEP 6: Generate All Windowed Sequences and Sliced SHAPE Profiles
# =========================================================================

window_sizes <- c(60, 50, 40, 30, 20)
window_types <- c("upstream", "centered", "downstream", "random")

all_windows_data <- window_definitions %>%
  crossing(size = window_sizes, type = window_types) %>%
  rowwise() %>%
  mutate(
    # --- Define window coordinates ---
    start_pos = case_when(
      type == "upstream"   ~ original_psipos - size + 1,
      type == "centered"   ~ original_psipos - (size / 2),
      type == "downstream" ~ original_psipos,
      type == "random"     ~ sample.int(original_seq_length - size + 1, 1)
    ),
    end_pos = start_pos + size - 1,
    start_pos = pmax(1, start_pos),
    end_pos = pmin(original_seq_length, end_pos),

    # --- Create the windowed sequence and new FASTA name ---
    window_seq = substr(original_seq, start_pos, end_pos),
    fasta_name = paste0(RNA_name, "_", type, size),

    # --- Slice the SHAPE data to match the window ---
    window_shape_data = list(
      shape_data %>%
        filter(Nucleotide >= start_pos, Nucleotide <= end_pos) %>%
        # Re-index the Nucleotide column to start from 1 for the new file
        mutate(Nucleotide = 1:n())
    )
  ) %>%
  ungroup() %>%
  # Keep only the necessary final columns
  select(fasta_name, window_seq, window_shape_data)

cat("Generated all windowed sequences and sliced SHAPE profiles.\n")

# =========================================================================
# STEP 7: Write All FASTA and .shape Files
# =========================================================================

# --- Write the FASTA file ---
fasta_sequences <- DNAStringSet(all_windows_data$window_seq)
names(fasta_sequences) <- all_windows_data$fasta_name
writeXStringSet(fasta_sequences, filepath = fasta_output_path)
cat(paste("Wrote", length(fasta_sequences), "sequences to", fasta_output_path, "\n"))

# --- Write the per-window .shape files ---
num_shape_files_written <- 0

for (i in 1:nrow(all_windows_data)) {
  fasta_name <- all_windows_data$fasta_name[i]
  shape_data <- all_windows_data$window_shape_data[[i]]
  
  # Ensure the data has the two required columns for .shape format
  data_to_write <- shape_data[, c("Nucleotide", "avg_reactivity")]
  
  output_file_path <- file.path(shape_output_dir, paste0(fasta_name, ".shape"))
  
  fwrite(data_to_write, 
         file = output_file_path, 
         sep = "\t", 
         col.names = FALSE)
  
  num_shape_files_written <- num_shape_files_written + 1
}

cat(paste("Wrote", num_shape_files_written, "windowed '.shape' files to:", shape_output_dir, "\n"))