#!/usr/bin/env Rscript

.libPaths("/home/users/rodell/R/x86_64-pc-linux-gnu-library/4.2")
library(future)
library(future.apply)
library(future.callr)
library(fs)
library(optparse)
library(dplyr)
library(stringr)

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- dirname(sub("--file=", "", args[grep("--file=", args)]))
  return(normalizePath(script_path))
}

SCRIPT_DIR <- get_script_path()
# SCRIPT_DIR <- "/scratch/users/rodell/motifmatcher/RNAstructure"

# Helper function to parse parameter strings from the command line
parse_param_string <- function(s) {
  if (is.null(s)) return(NULL)
  if (grepl(":", s)) return(eval(parse(text = s)))
  vals <- strsplit(s, ",")[[1]]
  if (all(toupper(vals) %in% c("TRUE", "FALSE"))) return(as.logical(vals))
  num_vals <- suppressWarnings(as.numeric(vals))
  if (!any(is.na(num_vals))) return(num_vals)
  return(vals)
}

# Define default parameter ranges (used if not provided on command line)
DEFAULT_PARAMETER_RANGES <- list(
  input_position = 59, 
  offset_min = 0:0,
  offset_max = 1:2,
  include_unpaired1 = c(FALSE, TRUE),
  paired1_min = 1:3,
  paired1_max = 4:9,
  unpaired2_min = 1:3,
  unpaired2_max = 4:9,
  include_paired2 = c(FALSE, TRUE)
)

DEFAULT_UNPAIRED1_RANGES <- list(
  unpaired1_min = 1:3,
  unpaired1_max = 4:9
)

DEFAULT_PAIRED2_RANGES <- list(
  paired2_min = 1:3,
  paired2_max = 4:9
)

# PARAMETER_RANGES <- list(
#   input_position = 59, 
#   offset_min = 0:0,
#   offset_max = 1:1,
#   include_unpaired1 = c(FALSE, TRUE),
#   paired1_min = 1:1,
#   paired1_max = 4:4,
#   unpaired2_min = 3:3,
#   unpaired2_max = 4:4,
#   include_paired2 = c(FALSE, TRUE)
# )

# UNPAIRED1_RANGES <- list(
#   unpaired1_min = 1:1,
#   unpaired1_max = 4:4
# )

# PAIRED2_RANGES <- list(
#   paired2_min = 1:1,
#   paired2_max = 4:4
# )

# Capture warnings
warning_messages <- character()
withCallingHandlers(
  {


# --- Define options ---
option_list <- list(
  make_option(c("-i", "--input"), type="character", default="input_pool1.csv",
              help="Input CSV file [default= %default]"),
  make_option(c("-f", "--fold_dir"), type="character", default="/scratch/users/rodell/RNAfold_psipos",
              help="Directory containing .fold files [default= %default]"),
  make_option(c("-o", "--output_dir"), type="character", default="parametersweep_output",
              help="Output directory for results [default= %default]"),
  make_option(c("-d", "--dataset"), type="character", default="pool1_psipos_v2.csv",
              help="Dataset file for evaluation [default= %default]"),

  # Options for parameter ranges (e.g., "1:5", "TRUE,FALSE")
  make_option("--param-input_position", type="character", default=NULL),
  make_option("--param-offset_min", type="character", default=NULL),
  make_option("--param-offset_max", type="character", default=NULL),
  make_option("--param-include_unpaired1", type="character", default=NULL),
  make_option("--param-paired1_min", type="character", default=NULL),
  make_option("--param-paired1_max", type="character", default=NULL),
  make_option("--param-unpaired2_min", type="character", default=NULL),
  make_option("--param-unpaired2_max", type="character", default=NULL),
  make_option("--param-include_paired2", type="character", default=NULL),
  make_option("--param-unpaired1_min", type="character", default=NULL),
  make_option("--param-unpaired1_max", type="character", default=NULL),
  make_option("--param-paired2_min", type="character", default=NULL),
  make_option("--param-paired2_max", type="character", default=NULL)


)

opt <- parse_args(OptionParser(option_list=option_list))

# --- Assign options ---
input_file <- opt$input
fold_dir <- opt$fold_dir
output_dir <- opt$output_dir
dataset_file <- opt$dataset

# input_file <- "input_pool1.csv"
# fold_dir <- "/scratch/users/rodell/RNAfold_psipos"
# output_dir <- "20250816/output_test4"
# dataset_file <- "pool1_psipos_v2.csv"

# Build parameter lists by overriding defaults with any provided command-line arguments.
PARAMETER_RANGES <- DEFAULT_PARAMETER_RANGES
UNPAIRED1_RANGES <- DEFAULT_UNPAIRED1_RANGES
PAIRED2_RANGES <- DEFAULT_PAIRED2_RANGES

update_param <- function(param_list, param_name, opt_value) {
  if (!is.null(opt_value)) {
    parsed_value <- parse_param_string(opt_value)
    param_list[[gsub("param-", "", param_name)]] <- parsed_value
  }
  return(param_list)
}

# Update parameter lists from command-line options
PARAMETER_RANGES <- update_param(PARAMETER_RANGES, "param-input_position", opt$`param-input_position`)
PARAMETER_RANGES <- update_param(PARAMETER_RANGES, "param-offset_min", opt$`param-offset_min`)
PARAMETER_RANGES <- update_param(PARAMETER_RANGES, "param-offset_max", opt$`param-offset_max`)
PARAMETER_RANGES <- update_param(PARAMETER_RANGES, "param-include_unpaired1", opt$`param-include_unpaired1`)
PARAMETER_RANGES <- update_param(PARAMETER_RANGES, "param-paired1_min", opt$`param-paired1_min`)
PARAMETER_RANGES <- update_param(PARAMETER_RANGES, "param-paired1_max", opt$`param-paired1_max`)
PARAMETER_RANGES <- update_param(PARAMETER_RANGES, "param-unpaired2_min", opt$`param-unpaired2_min`)
PARAMETER_RANGES <- update_param(PARAMETER_RANGES, "param-unpaired2_max", opt$`param-unpaired2_max`)
PARAMETER_RANGES <- update_param(PARAMETER_RANGES, "param-include_paired2", opt$`param-include_paired2`)

UNPAIRED1_RANGES <- update_param(UNPAIRED1_RANGES, "param-unpaired1_min", opt$`param-unpaired1_min`)
UNPAIRED1_RANGES <- update_param(UNPAIRED1_RANGES, "param-unpaired1_max", opt$`param-unpaired1_max`)

PAIRED2_RANGES <- update_param(PAIRED2_RANGES, "param-paired2_min", opt$`param-paired2_min`)
PAIRED2_RANGES <- update_param(PAIRED2_RANGES, "param-paired2_max", opt$`param-paired2_max`)

# Create output directory:
dir_create(output_dir)

# Load the dataset
pool1_df <- read.csv(dataset_file)

# --- Helper Functions ---

# Set up parallel strategy
parallel_strategy <- function(workers) plan(callr, workers = workers)

# Determine valid parameter combinations
generate_parameter_combinations <- function() {
  base_combinations <- expand.grid(PARAMETER_RANGES)
  
  # Generate all possible combinations for unpaired1 and paired2
  unpaired1_combinations <- expand.grid(UNPAIRED1_RANGES)
  paired2_combinations <- expand.grid(PAIRED2_RANGES)
  
  # Merge all combinations
  all_combinations <- merge(base_combinations, unpaired1_combinations)
  all_combinations <- merge(all_combinations, paired2_combinations)
  
  # Set values to NA where the region is not included
  all_combinations$unpaired1_min[!all_combinations$include_unpaired1] <- NA
  all_combinations$unpaired1_max[!all_combinations$include_unpaired1] <- NA
  all_combinations$paired2_min[!all_combinations$include_paired2] <- NA
  all_combinations$paired2_max[!all_combinations$include_paired2] <- NA
  
  valid_combinations <- all_combinations[
    (is.na(all_combinations$unpaired1_min) | all_combinations$unpaired1_min <= all_combinations$unpaired1_max) &
    all_combinations$paired1_min <= all_combinations$paired1_max &
    all_combinations$unpaired2_min <= all_combinations$unpaired2_max &
    all_combinations$offset_min <= all_combinations$offset_max &
    (is.na(all_combinations$paired2_min) | all_combinations$paired2_min <= all_combinations$paired2_max),
  ]
  
  # Define the desired column order
  column_order <- c(
    "input_position", "offset_min", "offset_max",
    "include_unpaired1", "unpaired1_min", "unpaired1_max",
    "paired1_min", "paired1_max",
    "unpaired2_min", "unpaired2_max",
    "include_paired2", "paired2_min", "paired2_max"
  )
  
  # Reorder the columns
  valid_combinations <- valid_combinations[, column_order]
  
  return(unique(valid_combinations))
}

# Function to calculate F1 score
calculate_f1_score <- function(actual, predicted) {
  true_positives <- sum(actual == 1 & predicted == 1)
  false_positives <- sum(actual == 0 & predicted == 1)
  false_negatives <- sum(actual == 1 & predicted == 0)
  
  precision <- ifelse(true_positives + false_positives > 0,
                      true_positives / (true_positives + false_positives),
                      0)
  recall <- ifelse(true_positives + false_negatives > 0,
                   true_positives / (true_positives + false_negatives),
                   0)
  
  f1_score <- ifelse(precision + recall > 0,
                     2 * (precision * recall) / (precision + recall),
                     0)
  
  return(f1_score)
}

calculate_sequence_f1_scores <- function(pool1_df) {
  unuar_both_f1 <- calculate_f1_score(pool1_df$both, pool1_df$UNUAR)
  uguag_both_f1 <- calculate_f1_score(pool1_df$both, pool1_df$UGUAG)
  
  return(sprintf("Sequence-only F1 Scores:\nUNUAR vs Both: %.4f\nUGUAG vs Both: %.4f", 
                 unuar_both_f1, uguag_both_f1))
}

create_run_id <- function(input_position, offset_min, offset_max, 
                          include_unpaired1, unpaired1_min, unpaired1_max,
                          paired1_min, paired1_max, 
                          unpaired2_min, unpaired2_max,
                          include_paired2, paired2_min, paired2_max) {
  paste0(
    "pos", input_position, "_",
    "off", offset_min, "-", offset_max, "_",
    if(include_unpaired1) paste0("u1_", unpaired1_min, "-", unpaired1_max, "_") else "no_u1_",
    "p1_", paired1_min, "-", paired1_max, "_",
    "u2_", unpaired2_min, "-", unpaired2_max,
    if(include_paired2) paste0("_p2_", paired2_min, "-", paired2_max) else "_no_p2"
  )
}

run_motif_matcher <- function(motif_matcher_path, output_file, params, fold_dir) {
  
  cmd_args <- c(
    motif_matcher_path,
    "--min_paired1", params$paired1_min,
    "--max_paired1", params$paired1_max,
    "--min_unpaired2", params$unpaired2_min,
    "--max_unpaired2", params$unpaired2_max,
    "--input_position", params$input_position,
    "--offset_min", params$offset_min,
    "--offset_max", params$offset_max,
    "--output", output_file,
    "--input", params$input_file,
    "--fold_dir", params$fold_dir
  )

  if(params$include_unpaired1) {
    cmd_args <- c(cmd_args,
      "--min_unpaired1", params$unpaired1_min,
      "--max_unpaired1", params$unpaired1_max,
      "--include_unpaired1", "TRUE"
    )
  } else {
    cmd_args <- c(cmd_args, "--include_unpaired1", "FALSE")
  }

  if(params$include_paired2) {
    cmd_args <- c(cmd_args,
      "--min_paired2", params$paired2_min,
      "--max_paired2", params$paired2_max,
      "--include_paired2", "TRUE"
    )
  } else {
    cmd_args <- c(cmd_args, "--include_paired2", "FALSE")
  }
  
  system2("Rscript", args = cmd_args, stdout = TRUE, stderr = TRUE)
}

calculate_statistics <- function(df, filtered_pool1_df) {
  
  if (nrow(df) == 0) {
    return(list(
      total_sequences = nrow(pool1_df),
      total_matches = 0,
      both_count = 0,
      incell_count = 0,
      invitro_count = 0,
      f1_both_structure = 0,
      f1_both_unuar_and_structure = 0,
      f1_both_uguag_and_structure = 0
    ))
  }
  
  structure_motif <- rep(1, nrow(df))
  unuar_and_structure <- ifelse(filtered_pool1_df$UNUAR == 1 & structure_motif == 1, 1, 0)
  uguag_and_structure <- ifelse(filtered_pool1_df$UGUAG == 1 & structure_motif == 1, 1, 0)
  
  list(
    total_sequences = nrow(pool1_df),
    total_matches = nrow(df),
    both_count = sum(filtered_pool1_df$both == 1),
    incell_count = sum(filtered_pool1_df$incell == 1),
    invitro_count = sum(filtered_pool1_df$invitro == 1),
    f1_both_structure = calculate_f1_score(filtered_pool1_df$both, structure_motif),
    f1_both_unuar_and_structure = calculate_f1_score(filtered_pool1_df$both, unuar_and_structure),
    f1_both_uguag_and_structure = calculate_f1_score(filtered_pool1_df$both, uguag_and_structure)
  )
}

filter_and_count <- function(df, pool_df, match_column = "chr", count_column = "both") {
  # Ensure df has a filename column
  if (!"filename" %in% names(df)) {
    stop("df must have a filename column")
  }
  
  # Filter pool_df based on the match_column
  filtered_pool_df <- pool_df %>%
    filter(sapply(!!sym(match_column), function(x) any(str_detect(df$filename, fixed(x)))))
  
  # Calculate statistics
  list(
    total_input = nrow(df),
    total_matches = nrow(filtered_pool_df),
    count = sum(filtered_pool_df[[count_column]] == 1, na.rm = TRUE)
  )
}

compare_initial_to_disruption <- function(initial_df, rnafold_df) {
  # Check if 'filename' column exists in both dataframes
  if (!"filename" %in% names(initial_df) || !"filename" %in% names(rnafold_df)) {
    warning("'filename' column not found in one or both dataframes")
    return(NA)
  }
  # Count sequences in initial_df not in disruption files
  disrupted_count <- sum(!initial_df$filename %in% rnafold_df$filename)
  
  return(disrupted_count)
}

# Implementation for mutant generation, run on motif_matches.csv for each individual parameter file
generate_mutants <- function(motif_matches_file, protect_pos = 60, run_id) {
    # Check if the file exists
    if (!file.exists(motif_matches_file)) {
      stop("motif_matches.csv not found: ", motif_matches_file)
    }

    test_output <- file.path(output_dir, "individual_results", run_id, "mutants.csv")
    #cat("Test output file:", test_output, "\n")
    
    # Construct the path to mutagenesis.R
    mutagenesis_path <- file.path(SCRIPT_DIR, "mutagenesis.R")
    
    # Check if the script exists
    if (!file.exists(mutagenesis_path)) {
      stop("mutagenesis.R not found in the script directory: ", SCRIPT_DIR)
    }
    
    # Call the mutagenesis.R script
    result <- system2("Rscript", 
                      args = c(mutagenesis_path, 
                               "-i", motif_matches_file, 
                               #"-o", temp_output,
                               "-o", test_output,
                               "--protect_position", protect_pos),
                      stdout = TRUE,
                      stderr = TRUE)
    
    
    # Check if the output file is missing or empty
    if (!file.exists(test_output) || file.info(test_output)$size == 0) {
      cat("Warning: Mutant output file", test_output, "is missing or empty for run_id:", run_id, "\n")
      return(data.frame())
    }

    # Read the results
    #mutant_results <- read.csv(temp_output, stringsAsFactors = FALSE)
    mutant_results <- read.csv(test_output, stringsAsFactors = FALSE)

    # # Clean up temporary file
    # file.remove(temp_output)

    return(mutant_results)
    
}

create_mutant_fastas <- function(input, output_dir) {
  #cat("Creating mutant FASTAs\n")
  #cat("Debug: Received output_dir:", output_dir, "\n")
  
  if (is.null(output_dir) || !dir.exists(output_dir)) {
    stop("Invalid or missing output directory: ", output_dir)
  }
  
  # Create a directory for FASTA files if it doesn't exist
  fasta_dir <- file.path(output_dir, "mutant_fastas")
  dir.create(fasta_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (is.character(input) && file.exists(input)) {
    # If input is a file path, read it
    mutant_data <- read.csv(input, stringsAsFactors = FALSE)
  } else if (is.data.frame(input)) {
    # If input is already a data frame, use it directly
    mutant_data <- input
  } else {
    stop("Invalid input: must be either a file path or a data frame")
  }
  
  # Construct the path to write_mutant_fastas.R
  write_mutant_fastas_path <- file.path(SCRIPT_DIR, "write_mutant_fastas.R")
  
  # Check if the script exists
  if (!file.exists(write_mutant_fastas_path)) {
    stop("write_mutant_fastas.R not found in the script directory: ", SCRIPT_DIR)
  }
  
  # Create a temporary CSV file from the mutant data
  temp_csv <- tempfile(fileext = ".csv")
  write.csv(mutant_data, temp_csv, row.names = FALSE)
  
  # Call the write_mutant_fastas.R script
  result <- tryCatch({
    system2("Rscript", 
            args = c(write_mutant_fastas_path, temp_csv, fasta_dir),
            stdout = TRUE,
            stderr = TRUE)
  }, error = function(e) {
    cat("Error occurred while running write_mutant_fastas.R:", e$message, "\n")
    return(NULL)
  })
  
  # Clean up temporary file
  file.remove(temp_csv)
  
  # Check for errors
  if (is.null(result)) {
    stop("Failed to execute write_mutant_fastas.R")
  } else if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    stop("Error in write_mutant_fastas.R: ", paste(result, collapse = "\n"))
  }
  
  # Get list of created FASTA files
  fasta_files <- list.files(fasta_dir, pattern = "\\.fa$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    cat("Warning: No FASTA files were created in", fasta_dir, "\n")
  } else {
    #cat("Created", length(fasta_files), "FASTA files in", fasta_dir, "\n")
  }
  
  return(fasta_files)
}

run_rnafold <- function(run_id, output_dir, script_dir) {
  # Construct paths
  mutant_fastas_dir <- file.path(output_dir, "individual_results", run_id, "mutant_fastas")
  rnafold_results_dir <- file.path(output_dir, "individual_results", run_id, "rnafold_results")
  rnafold_script_path <- file.path(script_dir, "run_RNAfold.sh")
  
  # Create rnafold_results directory if it doesn't exist
  if (!dir.exists(rnafold_results_dir)) {
    dir.create(rnafold_results_dir, recursive = TRUE)
  }
  
  # Get list of all .fa files in mutant_fastas directory
  fa_files <- list.files(mutant_fastas_dir, pattern = "\\.fa$", full.names = TRUE)
  
  # Run RNAfold for each .fa file
  for (fa_file in fa_files) {
    base_name <- tools::file_path_sans_ext(basename(fa_file))
    output_subdir <- file.path(rnafold_results_dir, base_name)
    
    # Create subdirectory for this mutant's results
    if (!dir.exists(output_subdir)) {
      dir.create(output_subdir)
    }
    
    # Run RNAfold bash script
    system2(command = "bash",
            args = c(rnafold_script_path, fa_file, output_subdir),
            stdout = "TRUE",
            stderr = "TRUE")
    
    
    # Check for .fold files
    fold_files <- list.files(output_subdir, pattern = "\\.fold$")
    if (length(fold_files) > 0) {
      #cat(sprintf("RNAfold successfully generated %d .fold file(s) for %s\n", 
                  #length(fold_files), base_name))
    } else {
      cat(sprintf("Warning: No .fold files found for %s\n", base_name))
    }
  }
  
  # Final check
  all_fold_files <- list.files(rnafold_results_dir, pattern = "\\.fold$", recursive = TRUE)
  #cat(sprintf("Total .fold files generated: %d\n", length(all_fold_files)))
}

find_motifs <- function(params, run_dir, output_filename = "motif_matches.csv", fold_dir = NULL) {
  output_file <- file.path(run_dir, output_filename)
  motif_matcher_path <- file.path(SCRIPT_DIR, "motif_matcher_v3.R")
  
  if (!file.exists(motif_matcher_path)) {
    stop(paste("motif_matcher_v3.R not found at", motif_matcher_path))
  }

  if (!is.null(fold_dir)) {
    params$fold_dir <- fold_dir
  }

  cmd_output <- run_motif_matcher(motif_matcher_path, output_file, params)
  
  # Check if the file is missing OR empty before trying to read it
  if (!file.exists(output_file) || file.info(output_file)$size <= 3) {
    cat("Warning: Output file", output_file, "is missing or empty for run:", basename(run_dir), "\n")
    cat("Command output was:\n", paste(cmd_output, collapse = "\n"), "\n")
    # Return an empty dataframe to avoid crashing downstream processes
    return(data.frame()) 
  }


  df <- read.csv(output_file)
  return(df)
}

run_initial_motif_analysis <- function(params, run_dir) {
  initial_result <- find_motifs(params, run_dir)
  if (is.null(initial_result)) {
    cat("Error: Initial analysis returned NULL\n")
    return(NULL)
  }
  
  filtered_pool1_df <- pool1_df %>%
    filter(sapply(chr, function(x) any(str_detect(initial_result$filename, fixed(x)))))
  
  initial_stats <- calculate_statistics(initial_result, filtered_pool1_df)
  
  return(list(result = initial_result, stats = initial_stats))
}

run_rnafold_and_motif_matching <- function(params, run_id, run_dir, output_dir, script_dir, initial_results) {
  
  #cat("Run RNAfold on mutants.. \n")
  
  rnafold_result <- run_rnafold(run_id = run_id, output_dir = output_dir, script_dir = script_dir)
  rnafold_results_dir <- file.path(output_dir, "individual_results", run_id, "rnafold_results")
  
  motif_match_results <- list()
  rnafold_stats <- list()
  
  for (subfolder in list.dirs(rnafold_results_dir, recursive = FALSE)) {
    base_name <- basename(subfolder)
    
    #cat("Check motifs for", base_name, "\n")
    
    dir.create(subfolder, recursive = TRUE, showWarnings = FALSE)
    
    rnafold_df <- find_motifs(params, run_dir, 
                              output_filename = paste0(base_name, "_motif_matches.csv"),
                              fold_dir = subfolder)
    
    
    motif_match_results[[base_name]] <- rnafold_df
    
    if (grepl("disruption", base_name)) {
      disrupted_count <- compare_initial_to_disruption(initial_results, rnafold_df)
      rnafold_stats[[base_name]] <- list(disrupted_count = disrupted_count)
    } else {
      rnafold_stats[[base_name]] <- list(disrupted_count = NA)
    }
  }
  
  # Compare mutant sequences
  sequence_comparison <- compare_mutant_sequences(motif_match_results)
  
  return(list(results = motif_match_results, stats = rnafold_stats, sequence_comparison = sequence_comparison))
}

# Helper function to safely extract filenames from a results data frame
extract_filenames <- function(results_list, key) {
  df <- results_list[[key]]
  
  # Check if df is a valid, non-empty data frame with the 'filename' column
  if (!is.null(df) && is.data.frame(df) && "filename" %in% names(df) && nrow(df) > 0) {
    return(df$filename)
  } else {
    # Otherwise, return an empty character vector (represents an empty set)
    return(character(0))
  }
}

# Identify sequences in disrupted and compensated and modified
compare_mutant_sequences <- function(motif_match_results) {
  # Safely extract sequence names using the new helper function
  paired1_comp_seq <- extract_filenames(motif_match_results, "paired1_compensatory")
  paired1_disr_seq <- extract_filenames(motif_match_results, "paired1_disruption")
  paired2_comp_seq <- extract_filenames(motif_match_results, "paired2_compensatory")
  paired2_disr_seq <- extract_filenames(motif_match_results, "paired2_disruption")
  combined_comp_seq <- extract_filenames(motif_match_results, "combined_compensatory")
  combined_disr_seq <- extract_filenames(motif_match_results, "combined_disruption")
  
  count_unique <- function(comp_seq, disr_seq) {
    if (is.null(comp_seq) || is.null(disr_seq)) return(list(diff = NA, count = NA))
    
    unique_set <- setdiff(comp_seq, disr_seq)
    
    # If unique_set is empty, we can return 0 counts without calling filter_and_count
    if (length(unique_set) == 0) {
      return(list(diff = 0, count = 0))
    }
    
    count_result <- filter_and_count(data.frame(filename = unique_set), pool1_df)
    list(diff = count_result$total_input, count = count_result$count)
  }
  
  paired1_result <- count_unique(paired1_comp_seq, paired1_disr_seq)
  paired2_result <- count_unique(paired2_comp_seq, paired2_disr_seq)
  combined_result <- count_unique(combined_comp_seq, combined_disr_seq)
  
  return(list(
    paired1_disr_comp = paired1_result$diff,
    paired1_disr_comp_mod = paired1_result$count,
    paired2_disr_comp = paired2_result$diff,
    paired2_disr_comp_mod = paired2_result$count,
    combined_disr_comp = combined_result$diff,
    combined_disr_comp_mod = combined_result$count
  ))
}

# functions to combine results
define_all_possible_columns <- function() {
  base_columns <- c(
    "input_position", "offset_min", "offset_max", 
    "include_unpaired1", "unpaired1_min", "unpaired1_max",
    "paired1_min", "paired1_max", 
    "unpaired2_min", "unpaired2_max",
    "include_paired2", "paired2_min", "paired2_max"
  )
  
  stat_columns <- c("total_sequences", "total_matches", "both_count", "incell_count", "invitro_count",
                    "f1_both_structure", "f1_both_unuar_and_structure", "f1_both_uguag_and_structure")
  
  initial_stat_columns <- paste0("initial_", stat_columns)
  
  filter_count_columns <- c("total_input", "total_matches", "count")

  rnafold_types <- c("paired1_disruption", "paired2_disruption")
  rnafold_columns <- paste0("rnafold_", rnafold_types, "_disrupted_count")
  
  comparison_columns <- c("paired1_disr_comp", "paired1_disr_comp_mod", "paired2_disr_comp", "paired2_disr_comp_mod", "combined_disr_comp", "combined_disr_comp_mod")

  return(c(base_columns, initial_stat_columns, rnafold_columns, comparison_columns))
}

generate_summary_df <- function(params, initial_analysis, rnafold_analysis, all_possible_columns) {
  summary_df <- data.frame(matrix(NA, nrow = 1, ncol = length(all_possible_columns)))
  colnames(summary_df) <- all_possible_columns
  
  # Fill in parameter values
  for (param_name in names(params)) {
    if (param_name %in% colnames(summary_df)) {
      summary_df[[param_name]] <- params[[param_name]]
    }
  }
  
  # Fill in initial analysis results
  for (stat_name in names(initial_analysis$stats)) {
    col_name <- paste0("initial_", stat_name)
    if (col_name %in% colnames(summary_df)) {
      summary_df[[col_name]] <- initial_analysis$stats[[stat_name]]
    }
  }
  
  # Fill in RNAfold analysis results
  for (name in names(rnafold_analysis$stats)) {
    col_name <- paste0("rnafold_", name, "_disrupted_count")
    if (col_name %in% colnames(summary_df)) {
      summary_df[[col_name]] <- rnafold_analysis$stats[[name]]$disrupted_count
    }
  }


  # Fill in sequence comparison results
  comparison_columns <- c("paired1_disr_comp", "paired1_disr_comp_mod", "paired2_disr_comp", "paired2_disr_comp_mod", "combined_disr_comp", "combined_disr_comp_mod")
  for (comp_name in comparison_columns) {
    if (comp_name %in% colnames(summary_df)) {
      summary_df[[comp_name]] <- rnafold_analysis$sequence_comparison[[comp_name]]
    }
  }
  
  return(summary_df)
}

# main running loop
run_parameter_sweep <- function(chunk_size = 4) {
  valid_combinations <- generate_parameter_combinations()
  total_combinations <- nrow(valid_combinations)
  all_possible_columns <- define_all_possible_columns()
  
  cat("Total number of combinations to test:", total_combinations, "\n")

  num_cores <- future::availableCores() - 1
  plan(multisession, workers = num_cores)

  results_list <- vector("list", total_combinations)
  
  for (chunk_start in seq(1, total_combinations, by = chunk_size)) {
  #for (chunk_start in seq(371, total_combinations, by = chunk_size)) {
    chunk_end <- min(chunk_start + chunk_size - 1, total_combinations)
    cat("Processing combinations", chunk_start, "to", chunk_end, "\n")
    
    ## KF: tells future_lapply to run sequentially, making debugging easier
    #plan(sequential)
    
    chunk_results <- future_lapply(chunk_start:chunk_end, function(i) {
      params <- valid_combinations[i,]
      params$input_file <- input_file
      params$fold_dir <- fold_dir
      
      run_id <- create_run_id(params$input_position, params$offset_min, params$offset_max,
                              params$include_unpaired1, params$unpaired1_min, params$unpaired1_max,
                              params$paired1_min, params$paired1_max,
                              params$unpaired2_min, params$unpaired2_max,
                              params$include_paired2, params$paired2_min, params$paired2_max)
      run_dir <- file.path(output_dir, "individual_results", run_id)
      
      dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
      #cat("Processing ", run_id, "\n")

      #cat("Running initial motif matcher analysis...\n")
      initial_analysis <- run_initial_motif_analysis(params, run_dir)
      if (is.null(initial_analysis)) return(NULL)
      
      #cat("Generating mutants...\n")
      motif_matches_file <- file.path(run_dir, "motif_matches.csv")
      if (!file.exists(motif_matches_file)) {
        cat("Error: motif_matches.csv does not exist\n")
        return(NULL)
      }
      
      mutant_result <- generate_mutants(motif_matches_file, run_id = run_id, protect_pos = params$input_position + 1)
      if (is.null(mutant_result) || !is.data.frame(mutant_result)) {
        cat("Error: generate_mutants failed\n")
        return(NULL)
      }

      #cat("Creating mutant FASTA files...\n")
      fasta_files <- create_mutant_fastas(input = mutant_result, output_dir = run_dir)
      
      #cat("Starting RNAfold and subsequent motif matching...\n")
      rnafold_analysis <- run_rnafold_and_motif_matching(params, run_id, run_dir, output_dir, SCRIPT_DIR, initial_analysis$result)
      
      #cat("Generating summary dataframe...\n")
      summary_df <- generate_summary_df(params, initial_analysis, rnafold_analysis, all_possible_columns)

      result <- list(
        params = params,
        initial_result = initial_analysis$result,
        mutant_result = mutant_result,
        fasta_files = fasta_files,
        rnafold_analysis = rnafold_analysis,
        summary = summary_df
      )
      
      run_dir_files <- list.files(run_dir, include.dirs = TRUE, full.names = TRUE, recursive = TRUE)
      unlink(run_dir_files, recursive = TRUE)
      
      saveRDS(result, file.path(run_dir, "intermediate_result.rds"))
      #write.csv(summary_df, file.path(run_dir, "summary.csv"), row.names = FALSE)

      #cat("Completed processing combination", i, "\n\n")
      return(summary_df)
    
    }, future.seed = TRUE)
    
    results_list[chunk_start:chunk_end] <- chunk_results
    
    # Checkpoint: save progress
    saveRDS(results_list, file.path(output_dir, "checkpoint.rds"))
  }
  
  # all_results <- do.call(rbind, results_list)
  
  valid_results <- results_list[!sapply(results_list, is.null)]
  
  # Check if there are any valid results to bind
  if (length(valid_results) == 0) {
    cat("Warning: No valid results were generated to create a final summary file.\n")
    return(data.frame()) # Return an empty data frame
  }
  
  all_results <- do.call(rbind, valid_results)
  
  write.csv(all_results, file.path(output_dir, "all_analysis_results.csv"), row.names = FALSE)
  
  return(all_results)
}


main <- function() {
  # Print the output directory
  cat("Results will be saved in:", output_dir, "\n")

  # Calculate and print sequence F1 scores
  cat("Calculating sequence-only F1 scores...\n")
  sequence_f1_scores <- calculate_sequence_f1_scores(pool1_df)
  cat(sequence_f1_scores, "\n\n")

  # Add a timestamp for the start of the process
  cat("Process started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  # Run parameter sweep
  cat("Starting parameter sweep. Results will be saved in:", output_dir, "\n")
  results <- run_parameter_sweep(chunk_size = 128)
  cat("Parameter sweep completed. Results are saved in:", output_dir, "\n")

  # Print summary of output files
  cat("\nSummary of output files:\n")
  cat("All analysis results:", file.path(output_dir, "all_analysis_results.csv"), "\n")
  cat("Individual analysis results for each parameter combination are also saved in this directory.\n")

  # Add a timestamp for the end of the process
  cat("Process finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  return(results)
}



# --- Run the main function ---

results <- main()



  },
  warning = function(w) {
    warning_messages <<- c(warning_messages, w$message)
  }
)


# Print captured warnings at the end
if (length(warning_messages) > 0) {
  cat("\nWarnings encountered during execution:\n")
  for (i in seq_along(warning_messages)) {
    cat(i, ": ", warning_messages[i], "\n", sep = "")
  }
} else {
  cat("\nNo warnings encountered during execution.\n")
}