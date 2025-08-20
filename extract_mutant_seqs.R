library(dplyr)
library(readr)
library(purrr)

# --- 1. CONFIGURATION ---
input_dir <- "/scratch/users/rodell/motifmatcher/20250816/analysis/bestmoremut" 
output_file <- file.path(input_dir, "bestmoremut.csv")

mutation_types <- c(
  "paired1_disruption", "paired1_compensatory",
  "paired2_disruption", "paired2_compensatory",
  "combined_disruption", "combined_compensatory"
)


# --- 2. LOAD DATA ---
cat("Loading data...\n")

# input csv with all the mutated sequences
main_data_path <- file.path(input_dir, "mutagenesis.csv")
if (!file.exists(main_data_path)) {
  stop("Error: 'mutagenesis.csv' not found in: ", input_dir)
}
main_data <- read_csv(main_data_path, show_col_types = FALSE)

# csvs corresponding to motifmatcher for the mutants
found_filenames <- mutation_types %>%
  set_names() %>%
  map(~{
    filepath <- file.path(input_dir, paste0(.x, "_motifmatcher2.csv"))
    if (file.exists(filepath)) {
      read_csv(filepath, show_col_types = FALSE, col_select = "filename") %>% pull(filename)
    } else {
      cat("Warning: Check file not found:", filepath, "\n")
      character(0)
    }
  })


# --- 3. PERFORM ANALYSIS ---
cat("Performing analysis...\n")

analysis_results <- main_data %>%
  
  mutate(chr = sub("\\.fold$", "", filename), .after = filename) %>%
  
  mutate(
    # --- Disruption Columns ---
    # Rule: If filename is in the corresponding check file, it fails (NA). Otherwise, it succeeds (U->T).
    paired1_disruption    = if_else(filename %in% found_filenames$paired1_disruption,    NA_character_, gsub("U", "T", paired1_disruption)),
    paired2_disruption    = if_else(filename %in% found_filenames$paired2_disruption,    NA_character_, gsub("U", "T", paired2_disruption)),
    combined_disruption = if_else(filename %in% found_filenames$combined_disruption, NA_character_, gsub("U", "T", combined_disruption)),
    
    # --- Compensatory Columns ---
    # Rule: If its disruption counterpart (above) is NOT NA AND its own filename is in its check file, it succeeds. Otherwise, it fails (NA).
    paired1_compensatory    = if_else(!is.na(paired1_disruption)    & filename %in% found_filenames$paired1_compensatory,    gsub("U", "T", paired1_compensatory),    NA_character_),
    paired2_compensatory    = if_else(!is.na(paired2_disruption)    & filename %in% found_filenames$paired2_compensatory,    gsub("U", "T", paired2_compensatory),    NA_character_),
    combined_compensatory = if_else(!is.na(combined_disruption) & filename %in% found_filenames$combined_compensatory, gsub("U", "T", combined_compensatory), NA_character_)
  )


# --- 4. CALCULATE AND DISPLAY COUNTS ---
cat("Calculating final counts...\n\n")

final_counts <- analysis_results %>%
  summarise(across(all_of(mutation_types), ~sum(!is.na(.))))

print(as.data.frame(final_counts))
cat("\n")

  
# --- 5. WRITE RESULTS ---
cat("Writing results to:", output_file, "\n")
write_csv(analysis_results, output_file)

cat("Analysis complete.\n")