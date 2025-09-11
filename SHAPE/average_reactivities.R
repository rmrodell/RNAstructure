#!/usr/bin/env Rscript

# =========================================================================
# SETUP: Load Libraries and Define Paths
# =========================================================================
library(data.table)

# --- Input & Output Paths ---
rep1_file <- "/scratch/groups/nicolemm/rodell/SHAPE/InVitro_Rep1_20241209/shapemapper_fullpool_HosseinRecs/Rep1_NoPUS_profile_annotated.txt"
rep2_file <- "/scratch/groups/nicolemm/rodell/SHAPE/InVitro_Rep2_20250731/shapemapper_fullpool_HosseinRecs/Rep2_NoPUS_profile_annotated.txt"
output_dir <- "/scratch/groups/nicolemm/rodell/SHAPE/InVitro_2rep_Analysis_fullpool/averaged_shape_files"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("Starting SHAPE profile averaging...\n")

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

final_shapes <- merge(full_template, averaged_reactivities, by = c("RNA_name", "Nucleotide"), all.x = TRUE)
final_shapes[is.na(avg_reactivity), avg_reactivity := -999]

cat("Generated full-length profiles and filled missing nucleotides with -999.\n")

# =========================================================================
# STEP 5: Write .shape Files for Each RNA
# =========================================================================
shape_list <- split(final_shapes, by = "RNA_name")

num_files_written <- 0
for (rna in names(shape_list)) {
  output_file_path <- file.path(output_dir, paste0(rna, ".shape"))
  data_to_write <- shape_list[[rna]][, .(Nucleotide, avg_reactivity)]
  
  fwrite(data_to_write, 
         file = output_file_path, 
         sep = "\t", 
         col.names = FALSE)
  
  num_files_written <- num_files_written + 1
}

cat(paste("\nProcess complete. Wrote", num_files_written, "'.shape' files to:", output_dir, "\n"))