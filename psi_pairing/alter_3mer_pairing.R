#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
#
# RNA Pairing Status Flip Pipeline
#
# Description:
# This script designs mutations in RNA sequences to "flip" the pairing status
# of a user-defined 3-nucleotide target region. It can change a paired region
# to an unpaired one, and vice-versa.
#
# Workflow:
# 1. Reads RNA sequences and their structures from `.fold` files.
# 2. For each RNA, it determines the pairing status of the target region.
# 3. Applies a specific mutation strategy (e.g., break pairs, create a hairpin)
#    to achieve the desired structural flip.
# 4. Generates a FASTA file of the mutated sequences.
# 5. Calls an external script (`run_RNAfold.sh`) to predict the new structures.
# 6. Analyzes the results to determine if the mutations were successful.
# 7. Outputs a detailed CSV report and a summary text file.
#
# Dependencies:
# - R packages: optparse, Biostrings, dplyr, purrr, this.path
# - External software: ViennaRNA package (specifically `RNAfold`).
# - A companion script `run_RNAfold.sh` located in the same directory.
#
# -----------------------------------------------------------------------------


# ---- 1. Load Libraries ----
# Suppress startup messages for a cleaner console output.
suppressPackageStartupMessages(library(optparse))       # For parsing command-line arguments.
suppressPackageStartupMessages(library(Biostrings))     # For handling DNA/RNA sequences.
suppressPackageStartupMessages(library(dplyr))          # For data manipulation (e.g., mutate, left_join).
suppressPackageStartupMessages(library(purrr))          # For functional programming (e.g., map, compact).
suppressPackageStartupMessages(library(this.path))      # To reliably find the script's own directory.


# ---- 2. Parse Command-Line Arguments ----
# Define the command-line options the script will accept.
option_list <- list(
  make_option(c("-f", "--fold_dir"), type = "character", default = NULL, 
              help = "Required. Path to the input directory containing original .fold files.", metavar = "DIR"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL, 
              help = "Required. Path to the directory where all output files will be saved.", metavar = "DIR"),
  make_option(c("-t", "--target_pos"), type = "integer", default = NULL, 
              help = "Required. The 1-based CENTER of a 3-nucleotide block to flip (+/- 1).", metavar = "INT"),
  make_option(c("-p", "--protect_pos"), type = "integer", default = NULL,
              help = "Optional. The 1-based CENTER of a 5-nucleotide block to protect from mutation (+/- 2).", metavar = "INT")
)

# Create the parser and parse the arguments.
opt_parser <- OptionParser(option_list = option_list, description = "Mutates RNA to flip the pairing status of a 3-mer region.")
opt <- parse_args(opt_parser)

# Validate that all required arguments have been provided. If not, print help and stop.
if (is.null(opt$fold_dir) || is.null(opt$output_dir) || is.null(opt$target_pos)) {
  print_help(opt_parser)
  stop("All three arguments (--fold_dir, --output_dir, --target_pos) must be supplied.", call. = FALSE)
}


# ---- 3. Configure Paths and Parameters ----
# This section centralizes all dynamic paths and key parameters for easy management.
script_dir <- this.path::this.dir() # Get the directory where this script is located.

CONFIG <- list(
  # Define the 3-nucleotide target region based on the user-provided center position.
  target_positions = (opt$target_pos - 1):(opt$target_pos + 1),
  
  # Define the 5-nucleotide region to protect from modification
  protected_positions = if (!is.null(opt$protect_pos)) {
      (opt$protect_pos - 2):(opt$protect_pos + 2)
    } else {
      NULL
    },

  # Input/Output paths
  org_fold_dir    = opt$fold_dir,
  output_dir      = opt$output_dir,
  
  # Path to the helper script for running RNAfold.
  rna_fold_script = file.path(script_dir, "run_RNAfold.sh"),
  
  # Paths for intermediate files generated during the workflow.
  intermediate_dir   = file.path(opt$output_dir, "intermediate"),
  mutated_fasta_path = file.path(opt$output_dir, "intermediate", "mutated_sequences.fa"),
  mut_fold_dir       = file.path(opt$output_dir, "intermediate", "mutated_structures"),
  
  # Paths for the final output reports.
  output_csv         = file.path(opt$output_dir, "pairing_status_flip_analysis.csv"),
  output_summary     = file.path(opt$output_dir, "summary_report.txt")
)


# ---- 4. Helper and Core Logic Functions ----

#' Read a single .fold file.
#'
#' Assumes a ViennaRNA-style format where line 1 is the name, line 2 is the
#' sequence, and line 6 contains the dot-bracket structure.
#'
#' @param file_path Path to the .fold file.
#' @return A data frame with name, seq, and dot_bracket, or NULL if the file is malformed.
read_fold_file <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  # Check if the file has enough lines to be valid.
  if (length(lines) >= 3) {
    name <- sub(">", "", lines[1])
    seq <- lines[2]
    # The dot-bracket structure is the first element on line 6, space-separated.
    dot_bracket <- strsplit(lines[6], " ")[[1]][1]
    return(data.frame(name = name, seq = seq, dot_bracket = dot_bracket, stringsAsFactors = FALSE))
  } else {
    # Return NULL for invalid files, to be filtered out later.
    return(NULL)
  }
}

#' Read all .fold files in a directory.
#'
#' @param fold_dir Path to the directory containing .fold files.
#' @return A single data frame combining data from all valid .fold files.
read_all_fold_files <- function(fold_dir) {
  fold_files <- list.files(fold_dir, pattern = "\\.fold$", full.names = TRUE)
  if (length(fold_files) == 0) stop(paste("No .fold files found in:", fold_dir))
  
  # Apply read_fold_file to each file and bind the results into one data frame.
  do.call(rbind, lapply(fold_files, read_fold_file))
}

#' Create and save a FASTA file from a data frame.
#'
#' @param df The data frame containing sequence data.
#' @param seq_column The name of the column with the sequences.
#' @param name_column The name of the column with the sequence names.
#' @param fasta_path The file path to save the FASTA file.
create_and_save_fasta <- function(df, seq_column, name_column, fasta_path) {
  if (nrow(df) == 0) return()
  
  # Create a named character vector (required by Biostrings).
  named_sequences <- setNames(df[[seq_column]], df[[name_column]])
  
  # Biostrings' DNAStringSet requires 'T' instead of 'U'. RNAfold handles this correctly.
  dna_sequences <- gsub("U", "T", named_sequences)
  
  # Write the sequences to a FASTA file. width=10000 ensures each sequence is on a single line.
  writeXStringSet(DNAStringSet(dna_sequences), fasta_path, width = 10000)
}

#' Get the pairing status of a single nucleotide.
#'
#' @param dot_bracket The full dot-bracket structure string.
#' @param position The 1-based position to check.
#' @return "unpaired" if the character is '.', "paired" otherwise.
get_pairing_status <- function(dot_bracket, position) {
  if (substr(dot_bracket, position, position) == ".") "unpaired" else "paired"
}

#' Get the overall pairing status of a region.
#'
#' @param dot_bracket The full dot-bracket structure string.
#' @param positions A vector of 1-based positions to check.
#' @return "all_paired", "all_unpaired", or "mixed".
get_region_status <- function(dot_bracket, positions) {
  statuses <- sapply(positions, function(p) get_pairing_status(dot_bracket, p))
  if (all(statuses == "paired")) {
    "all_paired"
  } else if (all(statuses == "unpaired")) {
    "all_unpaired"
  } else {
    "mixed"
  }
}

#' Find the pairing partner of a nucleotide.
#'
#' This function correctly handles nested structures by using a counter.
#'
#' @param dot_bracket The full dot-bracket structure string.
#' @param target_pos The 1-based position of the nucleotide whose partner is sought.
#' @return The 1-based position of the partner, or NA if unpaired or not found.
find_pairing_partner <- function(dot_bracket, target_pos) {
  char_at_target <- substr(dot_bracket, target_pos, target_pos)
  if (char_at_target == '.') return(NA)
  
  chars <- strsplit(dot_bracket, "")[[1]]
  balance <- 1
  
  if (char_at_target == '(') { # Search forward for the matching ')'
    for (i in (target_pos + 1):length(chars)) {
      if (chars[i] == '(') balance <- balance + 1
      if (chars[i] == ')') balance <- balance - 1
      if (balance == 0) return(i)
    }
  } else if (char_at_target == ')') { # Search backward for the matching '('
    for (i in (target_pos - 1):1) {
      if (chars[i] == ')') balance <- balance + 1
      if (chars[i] == '(') balance <- balance - 1
      if (balance == 0) return(i)
    }
  }
  return(NA) # Should not be reached in a valid structure
}

#' Get the Watson-Crick complement of an RNA base.
#'
#' @param base A single character: "A", "U", "C", or "G".
#' @return The complementary base, or NA for invalid input.
get_complement <- function(base) {
  case_when(
    base == "A" ~ "U",
    base == "U" ~ "A",
    base == "G" ~ "C",
    base == "C" ~ "G",
    TRUE ~ NA_character_
  )
}

# ---- 5. Mutation Strategy Functions ----

#' Strategy 1: Mutate to break existing base pairs.
#'
#' This function aims to make the target region unpaired. It finds the pairing
#' partner of each nucleotide in the target region and mutates the *partner*
#' to disrupt the bond (e.g., for an A-U pair, it mutates U to A).
#'
#' @param org_data A data frame row with original sequence and structure.
#' @param config The configuration list.
#' @return A list containing the mutated sequence and mutation info, or NULL if no mutations were made.
mutate_to_break_pairs <- function(org_data, config) {
  seq_vec <- strsplit(org_data$seq, "")[[1]]
  mutation_infos <- list()
  
  for (pos in config$target_positions) {
    if (get_pairing_status(org_data$dot_bracket, pos) == "paired") {
      partner_pos <- find_pairing_partner(org_data$dot_bracket, pos)
      if (!is.na(partner_pos)) {
        # Check if the partner position is in the protected region.
        if (!is.null(config$protected_positions) && partner_pos %in% config$protected_positions) {
          # If it is protected, skip this mutation and continue to the next.
          next
        }
        original_base_at_partner <- seq_vec[partner_pos]
        # Mutate the partner base to be the same as the target base, ensuring a mismatch.
        new_base_at_partner <- seq_vec[pos]
        seq_vec[partner_pos] <- new_base_at_partner
        # Record the mutation in standard notation (e.g., "U105A").
        mutation_infos[[length(mutation_infos) + 1]] <- sprintf("%s%d%s", original_base_at_partner, partner_pos, new_base_at_partner)
      }
    }
  }
  
  if (length(mutation_infos) == 0) return(NULL) # No mutations were applicable.
  
  return(list(
    sequence = paste(seq_vec, collapse = ""),
    info_string = paste(unlist(mutation_infos), collapse = ", "),
    strategy_used = "break_pairs"
  ))
}

#' Strategy 2a: Attempt to extend a nearby stem to pair the target region.
#'
#' This is the "stitch up" or "zip up" method. It looks for a pre-existing
#' stem near the target region and introduces complementary mutations in the
#' intervening unpaired nucleotides to extend the stem over the target region.
#'
#' @param org_data A data frame row with original sequence and structure.
#' @param config The configuration list.
#' @return A list containing the mutated sequence and info, or NULL if it fails.
attempt_extend_stem <- function(org_data, config) {
  seq_len <- nchar(org_data$seq)
  SEARCH_WINDOW <- 5 # How far (in nucleotides) to look for a nearby stem.
  
  # Search outwards from the target region for the nearest paired base.
  for (i in 1:SEARCH_WINDOW) {
    # Check upstream and downstream simultaneously.
    search_positions <- c(min(config$target_positions) - i, max(config$target_positions) + i)
    
    for (nearest_paired_base in search_positions) {
      # --- Condition Checks ---
      # Skip if out of bounds or if the found position is not actually paired.
      if (nearest_paired_base < 1 || nearest_paired_base > seq_len || 
          get_pairing_status(org_data$dot_bracket, nearest_paired_base) == "unpaired") next
      
      partner_of_nearest <- find_pairing_partner(org_data$dot_bracket, nearest_paired_base)
      if (is.na(partner_of_nearest)) next
      
      # --- Define Mutation Regions ---
      if (nearest_paired_base < min(config$target_positions)) { # Found an upstream stem.
        # Define the two strands to be "zipped up".
        positions_to_pair_1 <- (nearest_paired_base + 1):max(config$target_positions)
        positions_to_pair_2 <- rev((partner_of_nearest - length(positions_to_pair_1)):(partner_of_nearest - 1))
      } else { # Found a downstream stem.
        positions_to_pair_1 <- min(config$target_positions):(nearest_paired_base - 1)
        positions_to_pair_2 <- (partner_of_nearest + 1):(partner_of_nearest + length(positions_to_pair_1))
      }
      
      # --- Conflict Checks ---
      full_zip_region <- c(positions_to_pair_1, positions_to_pair_2)
      # Fail if any position is out of bounds, or if any position is already paired (conflict).
      if (any(full_zip_region < 1) || any(full_zip_region > seq_len) || 
          any(sapply(full_zip_region, function(p) get_pairing_status(org_data$dot_bracket, p) != "unpaired"))) next
      
      # Check for conflict with the protected region before attempting mutation.
      if (!is.null(config$protected_positions) && any(positions_to_pair_2 %in% config$protected_positions)) {
        # If any mutation site is protected, this attempt is invalid. Try the next search position.
        next
      }

      # --- Perform Mutation ---
      seq_vec <- strsplit(org_data$seq, "")[[1]]
      mutation_infos <- character(length(positions_to_pair_1))
      possible <- TRUE
      
      # Loop through the positions to be paired and apply complementary mutations.
      for (j in seq_along(positions_to_pair_1)) {
        pos1 <- positions_to_pair_1[j]
        pos2 <- positions_to_pair_2[j]
        complement <- get_complement(seq_vec[pos1])
        if (is.na(complement)) {
          possible <- FALSE
          break # Stop if a base is not A, U, C, G (e.g., 'N').
        }
        original_base <- seq_vec[pos2]
        seq_vec[pos2] <- complement
        mutation_infos[j] <- sprintf("%s%d%s", original_base, pos2, complement)
      }
      
      if (!possible) next # If mutation failed, try the next search position.
      
      # If successful, return the result immediately.
      return(list(
        sequence = paste(seq_vec, collapse = ""),
        info_string = paste(mutation_infos, collapse = ", "),
        strategy_used = "extend_stem"
      ))
    }
  }
  return(NULL) # Return NULL if no suitable stem was found to extend.
}


#' Strategy 2b: Mutate to create a new de novo hairpin.
#'
#' This is a fallback strategy if `attempt_extend_stem` fails. It creates a
#' new stem-loop by mutating a downstream region to be complementary to the
#' 5-nucleotide target region, separated by a short loop.
#'
#' @param org_data A data frame row with original sequence and structure.
#' @param config The configuration list.
#' @return A list containing the mutated sequence and info, or NULL if it fails.
mutate_to_create_hairpin <- function(org_data, config) {
  seq_len <- nchar(org_data$seq)
  loop_size <- 4 # A standard small hairpin loop size.
  
  # Define the region that will become the other side of the new stem.
  target_stem_positions <- config$target_positions
  partner_start_pos <- max(target_stem_positions) + loop_size + 1
  partner_end_pos <- partner_start_pos + length(target_stem_positions) - 1
  partner_stem_positions <- partner_start_pos:partner_end_pos
  
  # --- Conflict Checks ---
  # Fail if the new partner region would go past the end of the sequence.
  if (partner_end_pos > seq_len) return(NULL)
  
  # Fail if any nucleotide in the entire proposed hairpin region is already paired.
  full_hairpin_region <- c(target_stem_positions, partner_stem_positions)
  if (any(sapply(full_hairpin_region, function(p) get_pairing_status(org_data$dot_bracket, p) != "unpaired"))) {
    return(NULL)
  }
  
  # Check if the region to be mutated overlaps with the protected region.
  if (!is.null(config$protected_positions) && any(partner_stem_positions %in% config$protected_positions)) {
    # If there is an overlap, this strategy is not possible.
    return(NULL)
  }

  # --- Perform Mutation ---
  seq_vec <- strsplit(org_data$seq, "")[[1]]
  mutation_infos <- character(length(target_stem_positions))
  possible <- TRUE
  
  # The partner stem positions must be iterated in reverse to match correctly.
  rev_partner_stem_positions <- rev(partner_stem_positions)
  
  for (i in seq_along(target_stem_positions)) {
    pos1 <- target_stem_positions[i]
    pos2 <- rev_partner_stem_positions[i]
    
    complement <- get_complement(seq_vec[pos1])
    if (is.na(complement)) {
      possible <- FALSE
      break
    }
    
    original_base <- seq_vec[pos2]
    seq_vec[pos2] <- complement
    mutation_infos[i] <- sprintf("%s%d%s", original_base, pos2, complement)
  }
  
  if (!possible) return(NULL)

  return(list(
    sequence = paste(seq_vec, collapse = ""),
    info_string = paste(mutation_infos, collapse = ", "),
    strategy_used = "seed_hairpin_3bp"
  ))
}

#' Controller function to manage strategies for CREATING pairs.
#'
#' This function orchestrates the pairing strategies in order of priority.
#' It first tries the more minimal "extend_stem" method. If that fails,
#' it falls back to creating a new hairpin from scratch.
#'
#' @param org_data A data frame row with original sequence and structure.
#' @param config The configuration list.
#' @return A list containing the mutated sequence and info, or NULL if all strategies fail.
mutate_to_create_pairs <- function(org_data, config) {
  # Priority 1: Try the minimal "stitch up" method first.
  result <- attempt_extend_stem(org_data, config)
  
  # If that was successful, return the result immediately.
  if (!is.null(result)) {
    return(result)
  }
  
  # Priority 2: If extending failed, fall back to creating a new hairpin.
  result <- mutate_to_create_hairpin(org_data, config)
  
  return(result)
}


# ---- 6. Main Workflow ----

main_workflow <- function(config) {
  
  # Create output directories if they don't exist.
  dir.create(config$intermediate_dir, showWarnings = FALSE, recursive = TRUE)
  
  # --- Step 1: Read data and plan mutations ---
  cat("Step 1: Reading input files and planning mutations...\n")
  all_fold_files <- list.files(config$org_fold_dir, pattern = "\\.fold$", full.names = TRUE)
  
  # Use purrr::map to iterate through each file and generate a list of mutation plans.
  mutation_plan_list <- purrr::map(all_fold_files, function(file) {
    org_data <- read_fold_file(file)
    # Skip if file is invalid or sequence is too short for the target region.
    if (is.null(org_data) || nchar(org_data$seq) < max(config$target_positions)) return(NULL)
    
    original_status <- get_region_status(org_data$dot_bracket, config$target_positions)
    mutation_results_list <- list()
    
    # Apply strategies based on the original pairing status.
    if (original_status == "all_unpaired") {
      # Goal: Make paired.
      mutation_results_list[[1]] <- mutate_to_create_pairs(org_data, config)
    } else if (original_status == "all_paired") {
      # Goal: Make unpaired.
      mutation_results_list[[1]] <- mutate_to_break_pairs(org_data, config)
    } else if (original_status == "mixed") {
      # Goal: Try both making it fully paired and fully unpaired.
      mutation_results_list[[1]] <- mutate_to_break_pairs(org_data, config)
      mutation_results_list[[2]] <- mutate_to_create_pairs(org_data, config)
    }
    
    # Filter out any strategies that failed (returned NULL).
    valid_results <- purrr::compact(mutation_results_list)
    if (length(valid_results) == 0) return(NULL)
    
    # Format the valid results into a data frame for each successful strategy.
    purrr::map(valid_results, function(res) {
      data.frame(
        name = org_data$name,
        original_seq = org_data$seq,
        original_dot_bracket = org_data$dot_bracket,
        target_pos_range = paste(min(config$target_positions), max(config$target_positions), sep = "-"),
        original_status = original_status,
        strategy = res$strategy_used,
        mutated_seq = res$sequence,
        mutation_info = res$info_string,
        stringsAsFactors = FALSE
      )
    })
  }) %>% purrr::flatten() # Flatten the nested list of data frames into a single list.
  
  # Combine all individual data frames into one large one.
  mutation_plan_df <- bind_rows(mutation_plan_list)
  if (nrow(mutation_plan_df) == 0) {
    cat("No sequences were eligible for mutation or no strategies were successful. Exiting.\n")
    return()
  }
  
  # --- Step 2: Run RNAfold on mutated sequences ---
  cat("Step 2: Generating FASTA and running RNAfold on", nrow(mutation_plan_df), "mutated sequences...\n")
  create_and_save_fasta(mutation_plan_df, "mutated_seq", "name", config$mutated_fasta_path)
  
  # Execute the external shell script.
  # Using system2 is safer and provides more control than system().
  status <- system2("bash", args = c(config$rna_fold_script, config$mutated_fasta_path, config$mut_fold_dir), stdout = TRUE, stderr = TRUE)
  
  # Check if the command failed.
  if (inherits(status, "try-error") || (is.numeric(status) && status != 0)) {
     cat("Error: The 'run_RNAfold.sh' script failed to execute. Please check that RNAfold is installed and the script has execute permissions.\n")
     cat("Error details:\n", status, "\n")
     return()
  }

  # --- Step 3: Read and analyze RNAfold results ---
  cat("Step 3: Analyzing RNAfold output...\n")
  mutated_fold_data <- read_all_fold_files(config$mut_fold_dir)
  if (nrow(mutated_fold_data) == 0) {
    cat("Warning: RNAfold did not produce any output .fold files. Cannot complete analysis.\n")
    return()
  }
  
  # Join the original mutation plan with the new structures predicted by RNAfold.
  final_results <- mutation_plan_df %>%
    left_join(mutated_fold_data %>% select(name, dot_bracket_mut = dot_bracket), by = "name")
  
  # Determine if each mutation attempt was successful.
  final_analysis <- final_results %>%
    rowwise() %>% # Process row-by-row for the next calculation.
    mutate(new_status = get_region_status(dot_bracket_mut, config$target_positions)) %>%
    ungroup() %>% # Switch back to standard data frame processing.
    mutate(
      is_successful = case_when(
        # A "break" strategy is successful if it results in an all_unpaired region.
        (original_status == "all_paired" | original_status == "mixed") & 
          new_status == "all_unpaired" & 
          strategy == "break_pairs" ~ TRUE,
        # A "create" strategy is successful if it results in an all_paired region.
        (original_status == "all_unpaired" | original_status == "mixed") & 
          new_status == "all_paired" & 
          strategy %in% c("extend_stem", "seed_hairpin_3bp") ~ TRUE,
        # All other outcomes are considered failures.
        TRUE ~ FALSE
      )
    )
  
  # --- Step 4: Save final reports ---
  cat("Step 4: Saving final reports...\n")
  write.csv(final_analysis, config$output_csv, row.names = FALSE)
  
  # Create a summary text for the report.
  summary_text <- c(
    "===== 3-mer Pairing Status Flip: Summary Report =====",
    "",
    paste("Target Region (Center", opt$target_pos, "):", min(config$target_positions), "to", max(config$target_positions)),
    paste("Total unique input sequences processed:", length(all_fold_files)),
    paste("Total mutation attempts generated:", nrow(final_analysis)),
    "---",
    paste("Total successful flips:", sum(final_analysis$is_successful, na.rm = TRUE)),
    paste("Success Rate:", scales::percent(sum(final_analysis$is_successful, na.rm = TRUE) / nrow(final_analysis))),
    "==================================================="
  )
  writeLines(summary_text, config$output_summary)
  
  cat("Workflow complete. Outputs are in:", config$output_dir, "\n")
}

# --- Execute the main function ---
main_workflow(CONFIG)