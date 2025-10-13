#!/usr/bin/env Rscript

# --- 1. Initialization and Setup ---

# Load required packages
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(Biostrings)
})

# --- Helper function to find the script's directory ---
get_script_dir <- function() {
  # Get the command arguments
  args <- commandArgs(trailingOnly = FALSE)
  # Find the --file argument
  file_arg <- args[grep("^--file=", args)]
  if (length(file_arg) == 0) {
    # Fallback for interactive sessions (like in RStudio)
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      return(dirname(rstudioapi::getActiveDocumentContext()$path))
    } else {
      # Fallback to current working directory if not in a script or RStudio
      return(getwd())
    }
  }
  # Extract path from --file argument
  return(dirname(sub("--file=", "", file_arg)))
}


# --- Helper function for parsing .fold files ---
parse_fold_file <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  if (length(lines) < 6) return(NULL)

  sequence_name <- sub(">", "", lines[1]) 
  sequence <- lines[2]
  
  mea_line <- lines[6]
  mea_structure <- str_trim(str_extract(mea_line, "^[().]+"))
  
  tibble(
    sequence_name = sequence_name, 
    original_filename = basename(filepath),
    sequence = sequence,
    original_mea = mea_structure
  )
}

# --- Helper function to convert dot-bracket to a pairing vector ---
db_to_pairs <- function(dot_bracket) {
  db <- strsplit(dot_bracket, "")[[1]]
  stack <- c()
  pairs <- rep(0, length(db))
  
  for (i in seq_along(db)) {
    if (db[i] == "(") {
      stack <- c(i, stack)
    } else if (db[i] == ")") {
      if (length(stack) > 0) {
        j <- stack[1]
        stack <- stack[-1]
        pairs[i] <- j
        pairs[j] <- i
      }
    }
  }
  return(pairs)
}

# --- Helper function to calculate structural identity ---
calculate_identity <- function(str1, str2) {
  if (is.na(str1) || is.na(str2) || nchar(str1) != nchar(str2)) {
    return(NA_real_)
  }
  s1 <- strsplit(str1, "")[[1]]
  s2 <- strsplit(str2, "")[[1]]
  
  hamming_dist <- sum(s1 != s2)
  identity <- (nchar(str1) - hamming_dist) / nchar(str1)
  return(identity)
}


# --- Core Mutagenesis Function ---
generate_mutations <- function(sequence, pairs_vec, window_indices, protect_indices) {
  
  seq_vec <- strsplit(sequence, "")[[1]]
  
  # Find indices within the window that are paired and not protected
  indices_to_mutate <- window_indices[
    window_indices %in% which(pairs_vec > 0) & 
    !window_indices %in% protect_indices
  ]
  
  # Find their partners, ensuring partners are also not protected
  partner_indices <- pairs_vec[indices_to_mutate]
  valid_pairs <- partner_indices > 0 & !partner_indices %in% protect_indices
  
  indices_to_mutate <- indices_to_mutate[valid_pairs]
  partner_indices <- partner_indices[valid_pairs]
  
  if (length(indices_to_mutate) == 0) {
    return(tibble(
      disrupt_seq = NA_character_, disrupt_mutations = NA_character_,
      comp_seq = NA_character_, comp_mutations = NA_character_
    ))
  }
  
  # Generate Disruption Mutation
  disrupt_vec <- seq_vec
  disrupt_vec[partner_indices] <- seq_vec[indices_to_mutate] # Partner becomes same as original
  disrupt_seq <- paste0(disrupt_vec, collapse = "")
  
  disrupt_mutations <- paste0(
    seq_vec[partner_indices], partner_indices, seq_vec[indices_to_mutate],
    collapse = ";"
  )

  # Generate Compensatory Mutation (from the DISRUPTED sequence)
  comp_vec <- disrupt_vec
  # Restore pairing by mutating the original base
  comp_vec[indices_to_mutate] <- case_when(
    disrupt_vec[partner_indices] == "A" ~ "U",
    disrupt_vec[partner_indices] == "U" ~ "A",
    disrupt_vec[partner_indices] == "G" ~ "C",
    disrupt_vec[partner_indices] == "C" ~ "G",
    TRUE ~ "N" # Fallback
  )
  comp_seq <- paste0(comp_vec, collapse = "")
  
  comp_mutations_part1 <- disrupt_mutations
  comp_mutations_part2 <- paste0(
    seq_vec[indices_to_mutate], indices_to_mutate, comp_vec[indices_to_mutate],
    collapse = ";"
  )
  comp_mutations <- paste(comp_mutations_part1, comp_mutations_part2, sep=";")
  
  return(tibble(
    disrupt_seq, disrupt_mutations,
    comp_seq, comp_mutations
  ))
}

# --- Main Execution Block ---
main <- function() {
  cat("--- Mutagenesis and Evaluation Script Started ---\n\n")

  # --- 2. Parse Command-Line Arguments ---
  option_list <- list(
    make_option(c("--input_dir"), type = "character", help = "Directory of .fold files"),
    make_option(c("--output_dir"), type = "character", help = "Directory for final CSV results"),
    make_option(c("--position_zero"), type = "integer", help = "Central reference position (1-based)")
  )
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)

  if (is.null(opt$input_dir) || is.null(opt$output_dir) || is.null(opt$position_zero)) {
    stop("All arguments (--input_dir, --output_dir, --position_zero) are required.", call. = FALSE)
  }

  dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # --- Setup paths and parameters ---
  script_dir <- get_script_dir()
  run_rnafold_script <- file.path(script_dir, "run_RNAfold.sh")
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  temp_dir <- file.path(opt$output_dir, paste0("mutagenesis_run_", timestamp))
  dir.create(temp_dir, showWarnings = FALSE)

  cat(paste("Input directory:", opt$input_dir, "\n"))
  cat(paste("Output directory:", opt$output_dir, "\n"))
  cat(paste("Temporary directory:", temp_dir, "\n"))
  cat(paste("Position Zero:", opt$position_zero, "\n"))
  cat(paste("RNAfold script path:", run_rnafold_script, "\n\n"))

  # --- Define mutation windows and protection zone ---
  protect_indices <- (opt$position_zero - 2):(opt$position_zero + 2)
  upstream_window <- (opt$position_zero - 10):(opt$position_zero - 1)
  downstream_window <- (opt$position_zero + 1):(opt$position_zero + 15)
  total_window <- c(upstream_window, downstream_window)
  
  cat("Defined Regions (1-based indices):\n")
  cat(paste("  - Protection zone:", paste(protect_indices, collapse=", "), "\n"))
  cat(paste("  - Upstream window:", min(upstream_window), "to", max(upstream_window), "\n"))
  cat(paste("  - Downstream window:", min(downstream_window), "to", max(downstream_window), "\n\n"))
  
  # --- 3. Process Input Files and Generate All Mutants ---
  cat("--- Phase 1: Parsing input files and generating mutants ---\n")
  
  fold_files <- list.files(opt$input_dir, pattern = "\\.fold$", full.names = TRUE)
  if (length(fold_files) == 0) stop("No .fold files found in the input directory.")
  
  original_data <- map_dfr(fold_files, parse_fold_file) %>%
    mutate(
      original_seq_len = nchar(sequence),
      original_pairs = map(original_mea, db_to_pairs)
    )
  
  all_mutants <- list()
  
  for (i in 1:nrow(original_data)) {
    row <- original_data[i, ]
    cat(paste("Processing file:", row$original_filename, "\n"))
    
    # Generate mutations for each window
    mut_up <- generate_mutations(row$sequence, row$original_pairs[[1]], upstream_window, protect_indices)
    mut_down <- generate_mutations(row$sequence, row$original_pairs[[1]], downstream_window, protect_indices)
    mut_total <- generate_mutations(row$sequence, row$original_pairs[[1]], total_window, protect_indices)
    
    # Structure data for binding
    results_tbl <- bind_rows(
      tibble(type = "upstream_disruption", seq = mut_up$disrupt_seq, mutations = mut_up$disrupt_mutations, window = list(upstream_window)),
      tibble(type = "upstream_compensation", seq = mut_up$comp_seq, mutations = mut_up$comp_mutations, window = list(upstream_window)),
      tibble(type = "downstream_disruption", seq = mut_down$disrupt_seq, mutations = mut_down$disrupt_mutations, window = list(downstream_window)),
      tibble(type = "downstream_compensation", seq = mut_down$comp_seq, mutations = mut_down$comp_mutations, window = list(downstream_window)),
      tibble(type = "total_disruption", seq = mut_total$disrupt_seq, mutations = mut_total$disrupt_mutations, window = list(total_window)),
      tibble(type = "total_compensation", seq = mut_total$comp_seq, mutations = mut_total$comp_mutations, window = list(total_window))
    ) %>%
       mutate(
        original_filename = row$original_filename, 
        original_sequence_name = row$sequence_name, 
        .before = 1
      ) %>%
      filter(!is.na(seq)) # Remove cases where no mutations could be made
      
    all_mutants[[i]] <- results_tbl
  }
  
  mutants_df <- bind_rows(all_mutants)
  if (nrow(mutants_df) == 0) {
      cat("\nNo valid mutations could be generated for any input files. Exiting.\n")
      return()
  }

  mutants_df <- mutants_df %>%
    mutate(fasta_header = paste(str_remove(original_filename, "\\.fold$"), type, sep = "_"))
  
  cat(paste("\nGenerated a total of", nrow(mutants_df), "mutant sequences.\n\n"))

  # --- 4. Run RNAfold in Batch ---
  cat("--- Phase 2: Running RNAfold on all mutants ---\n")
  
  # Write all sequences to a single multi-FASTA file
  multi_fasta_path <- file.path(temp_dir, "all_mutants.fa")
  sequences_set <- RNAStringSet(mutants_df$seq)
  names(sequences_set) <- mutants_df$fasta_header
  writeXStringSet(sequences_set, multi_fasta_path, width = 200)
  cat(paste("Wrote all mutants to:", multi_fasta_path, "\n"))
  
  # Run the shell script
  cat("Executing run_RNAfold.sh. This may take some time...\n")
  rnafold_out_dir <- file.path(temp_dir, "rnafold_outputs")
  dir.create(rnafold_out_dir)

  system2_command <- run_rnafold_script
  system2_args <- c(multi_fasta_path, rnafold_out_dir)
  
  cat("Command:", system2_command, paste(system2_args, collapse = " "), "\n")
  
  exit_code <- system2(system2_command, args = system2_args, stdout = TRUE, stderr = TRUE)
  
  if (any(grepl("Error", exit_code, ignore.case = TRUE))) {
    cat("\n--- WARNING: RNAfold script reported an error. --- \n")
    print(exit_code)
    cat("--- Proceeding with evaluation, but results may be incomplete. ---\n\n")
  } else {
    cat("RNAfold script completed.\n\n")
  }

  # --- 5. Evaluate Mutation Success ---
  cat("--- Phase 3: Evaluating mutation success ---\n")

  # Merge back original data for comparison
  results <- mutants_df %>%
    left_join(original_data, by = "original_filename")

  # Parse new .fold files
  get_mutant_mea <- function(header) {
    mutant_fold_path <- file.path(rnafold_out_dir, paste0(header, ".fold"))
    if (!file.exists(mutant_fold_path)) {
      cat(paste("  - WARNING: Expected output file not found:", basename(mutant_fold_path), "\n"))
      return(NA_character_)
    }
    lines <- readLines(mutant_fold_path, n = 6, warn = FALSE)
    if (length(lines) < 6) return(NA_character_)
    return(str_trim(str_extract(lines[6], "^[().]+")))
  }
  
  results <- results %>%
    mutate(mutant_mea = map_chr(fasta_header, get_mutant_mea))

  # --- Calculate metrics ---
  cat("Calculating success metrics...\n")
  results <- results %>%
    mutate(
      # Metric 1: Global structural identity
      global_identity = map2_dbl(original_mea, mutant_mea, calculate_identity),
      
      # Metric 2: Local structural identity (for compensation)
      local_identity = pmap_dbl(list(original_mea, mutant_mea, window), function(o, m, w) {
        if (is.na(o) || is.na(m)) return(NA_real_)
        calculate_identity(substr(o, min(w), max(w)), substr(m, min(w), max(w)))
      }),
      
      # Metric 3: Disruption efficacy
    #   disruption_efficacy = pmap_dbl(list(original_mea, mutant_mea, window, original_pairs), function(o, m, w, p) {
    #     if (is.na(o) || is.na(m)) return(NA_real_)
        
    #     # Paired bases in the original structure within the window
    #     orig_paired_in_window <- w[substr(o, w, w) %in% c("(", ")")]
    #     if (length(orig_paired_in_window) == 0) return(1.0) # Or NA? Let's say 100% success if nothing to disrupt.
        
    #     # Of those, which are now unpaired in the mutant structure?
    #     now_unpaired <- sum(substr(m, orig_paired_in_window, orig_paired_in_window) == ".")
        
    #     return(now_unpaired / length(orig_paired_in_window))
    #   }),

      # Calculate the percentage of the window that is unpaired
      window_unpaired_percent = pmap_dbl(list(mutant_mea, window), function(m, w) {
        if (is.na(m)) return(NA_real_)
        
        # Extract the structure of the target window
        window_structure <- substr(m, min(w), max(w))
        if (nchar(window_structure) == 0) return(NA_real_) # Avoid division by zero
        
        # Count the number of dots (unpaired bases)
        unpaired_count <- str_count(window_structure, "\\.")
        
        return(unpaired_count / nchar(window_structure))
      }),

      # --- Extract the dot-bracket notation for the target window ---
      original_window_mea = pmap_chr(list(original_mea, window), function(mea, w) {
        if (is.na(mea)) return(NA_character_)
        substr(mea, min(w), max(w))
      }),
      mutant_window_mea = pmap_chr(list(mutant_mea, window), function(mea, w) {
        if (is.na(mea)) return(NA_character_)
        substr(mea, min(w), max(w))
      })
    )

  # --- Apply success criteria ---
  results <- results %>%
    mutate(
      is_disruption = str_detect(type, "disruption"),
      base_type = str_extract(type, "upstream|downstream|total")
    ) %>%
    mutate(
      success_criteria = case_when(
        is_disruption & base_type == "upstream"   ~ window_unpaired_percent >= 0.80 & global_identity >= 0.80,
        is_disruption & base_type == "downstream" ~ window_unpaired_percent >= 0.80 & global_identity >= 0.75,
        is_disruption & base_type == "total"      ~ window_unpaired_percent >= 0.80 & global_identity >= 0.60,
        !is_disruption ~ local_identity >= 0.90 & global_identity >= 0.80,
        TRUE ~ FALSE
      )
    )

  # Handle conditional success for compensation
  disruption_successes <- results %>%
    filter(is_disruption & success_criteria) %>%
    select(original_filename, base_type, disruption_success = success_criteria)

  results <- results %>%
    left_join(disruption_successes, by = c("original_filename", "base_type")) %>%
    mutate(
      success = case_when(
        is_disruption ~ success_criteria,
        !is_disruption & is.na(disruption_success) ~ FALSE, # Corresponding disruption failed
        !is_disruption & disruption_success == TRUE ~ success_criteria, # Disruption succeeded, so use compensation criteria
        TRUE ~ FALSE
      )
    ) %>%
    select(-is_disruption, -base_type, -success_criteria, -disruption_success)

  # --- 6. Save Final Outputs ---
  cat("\n--- Phase 4: Saving final results ---\n")
  
  # Select and reorder columns for clarity
  final_results <- results %>%
    select(
      original_filename,
      type,
      success,
      mutations,
      window_unpaired_percent,
      local_identity,
      global_identity,
      fasta_header,
      seq,
      original_window_mea,
      mutant_window_mea,
      original_mea,
      mutant_mea,
      sequence
    )
  
  successful_results <- final_results %>% filter(success == TRUE)

  # Write to CSV
  all_results_path <- file.path(opt$output_dir, "all_mutagenesis_results.csv")
  successful_results_path <- file.path(opt$output_dir, "successful_mutagenesis_results.csv")
  
  write.csv(final_results, all_results_path, row.names = FALSE)
  write.csv(successful_results, successful_results_path, row.names = FALSE)
  
  cat(paste("Saved all results to:", all_results_path, "\n"))
  cat(paste("Saved successful results to:", successful_results_path, "\n"))
  cat(paste("\nIntermediate files are located in:", temp_dir, "\n"))
  cat("--- Script Finished ---\n")
}

# Run the main function
main()