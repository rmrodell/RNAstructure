#!/usr/bin/env Rscript

# --- 1. Initialization and Setup ---

# --- Auto-install required packages ---
# This block checks for missing packages and installs them if necessary.

# List all required packages
cran_packages <- c("optparse", "dplyr", "stringr", "purrr", "tibble", "stringdist")
bioc_packages <- c("Biostrings")

# Function to install a list of packages from CRAN
install_cran <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(paste("Installing missing CRAN package:", pkg, "\n"))
      install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    }
  }
}

# Function to install a list of packages from Bioconductor
install_bioc <- function(pkg_list) {
  # First, ensure BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  # Now install the packages
  for (pkg in pkg_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat(paste("Installing missing Bioconductor package:", pkg, "\n"))
      BiocManager::install(pkg, ask = FALSE)
    }
  }
}

# Run the installation functions
install_cran(cran_packages)
install_bioc(bioc_packages)

# Load required packages
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(Biostrings)
  library(stringdist)
})

# --- Helper function to find the script's directory ---
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)]
  if (length(file_arg) == 0) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      return(dirname(rstudioapi::getActiveDocumentContext()$path))
    } else {
      return(getwd())
    }
  }
  return(dirname(sub("--file=", "", file_arg)))
}

# --- Helper function for parsing .fold files (MFE structure and energy) ---
parse_fold_file <- function(filepath, verbose = FALSE) {
  lines <- readLines(filepath, warn = FALSE)
  if (length(lines) < 3) return(NULL)

  sequence_name <- sub(">", "", lines[1]) 
  sequence <- lines[2]
  
  mfe_line <- lines[3]
  mfe_match <- str_match(mfe_line, "\\s+\\(([^)]+)\\)$")
  original_mfe <- as.numeric(mfe_match[1, 2])
  original_db <- str_trim(str_extract(mfe_line, ".*(?=\\s+\\()"))

  if (verbose) {
    cat(sprintf("\n[DEBUG] Parsing file: %s\n", basename(filepath)))
    cat(sprintf("  - Raw MFE line: '%s'\n", mfe_line))
    cat(sprintf("  - Extracted DB: '%s'\n", ifelse(is.na(original_db), "NA", original_db)))
    cat(sprintf("  - Extracted MFE: %s\n", ifelse(is.na(original_mfe), "NA", original_mfe)))
  }

  tibble(
    sequence_name = sequence_name, 
    original_filename = basename(filepath),
    sequence = sequence,
    original_db = original_db,
    original_mfe = original_mfe
  )
}

# --- Helper function to convert dot-bracket to a pairing vector ---
db_to_pairs <- function(dot_bracket, verbose = FALSE) {
  if (verbose) {
    cat(sprintf("    [DB_TO_PAIRS] Received input of class: %s\n", class(dot_bracket)))
    cat(sprintf("    [DB_TO_PAIRS] Input string: '%s'\n", dot_bracket))
  }

  if (is.na(dot_bracket) || nchar(dot_bracket) == 0) {
    if(verbose) cat("    [DB_TO_PAIRS] Input is NA or empty. Returning NA.\n")
    return(NA)
  }

  if (is.na(dot_bracket)) return(NA)
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

  if (verbose) {
    cat(sprintf("    [DEBUG C] 'db_to_pairs' is returning a vector. Number of paired bases found: %d\n", sum(pairs > 0)))
  }

  return(pairs)
}


# --- Core MFE Mutagenesis Function ---
generate_mfe_mutations <- function(sequence, pairs_vec, window_indices, protect_indices, verbose = FALSE) {

  if (verbose) {
    cat(sprintf("    [DEBUG E] 'generate_mfe_mutations' received 'pairs_vec'. Class: %s, Length: %d, Head: %s\n", 
                class(pairs_vec), length(pairs_vec), paste(head(pairs_vec), collapse=",")))
  }
  
  if(is.null(pairs_vec) || is.na(pairs_vec[1])) {
    if(verbose) cat("    [DEBUG] Received invalid pairs vector (NULL or NA). Cannot generate mutations.\n")
    return(tibble())
  }
  
  seq_vec <- strsplit(sequence, "")[[1]]
  
  # Identify all unique pairs (i, j) where i < j
  all_pairs_i <- which(pairs_vec > 0 & seq_along(pairs_vec) < pairs_vec)
  all_pairs_j <- pairs_vec[all_pairs_i]
  
  if (verbose) {
    cat(sprintf("    [DEBUG] Found %d total pairs in the original structure.\n", length(all_pairs_i)))
    if (length(all_pairs_i) > 0) {
       pair_strings <- paste0("(", all_pairs_i, ", ", all_pairs_j, ")")
       cat(paste0("      - Pairs: ", paste(pair_strings, collapse = ", "), "\n"))
    }
  }
  
  # --- Filter pairs step-by-step for debugging ---
  is_in_window <- (all_pairs_i %in% window_indices) | (all_pairs_j %in% window_indices)
  is_protected <- (all_pairs_i %in% protect_indices) | (all_pairs_j %in% protect_indices)
  valid_pair_indices <- is_in_window & !is_protected

  target_i <- all_pairs_i[valid_pair_indices]
  target_j <- all_pairs_j[valid_pair_indices]
  
  if (verbose) {
      cat(sprintf("    [DEBUG] Filtering pairs:\n"))
      cat(sprintf("      - Window: %d to %d\n", min(window_indices), max(window_indices)))
      cat(sprintf("      - Protection Zone: %d to %d\n", min(protect_indices), max(protect_indices)))
      cat(sprintf("      - Found %d pairs with at least one base in the window.\n", sum(is_in_window)))
      cat(sprintf("      - Found %d pairs with at least one base in the protection zone (will be excluded).\n", sum(is_protected)))
      cat(sprintf("      - Result: %d valid pairs to mutate in this window.\n", length(target_i)))
      if (length(target_i) > 0) {
          valid_pair_strings <- paste0("(", target_i, ", ", target_j, ")")
          cat(paste0("        - Valid pairs: ", paste(valid_pair_strings, collapse=", "), "\n"))
      }
  }
  
  if (length(target_i) == 0) {
    return(tibble()) # Return empty tibble if no pairs to mutate
  }
  
  # --- Separate pairs into GC/AU types ---
  base_i <- seq_vec[target_i]
  base_j <- seq_vec[target_j]
  
  is_gc_pair <- (base_i == "G" & base_j == "C") | (base_i == "C" & base_j == "G")
  is_au_pair <- (base_i == "A" & base_j == "U") | (base_i == "U" & base_j == "A")
  
  if (verbose) {
      cat(sprintf("    [DEBUG] Categorizing %d valid pairs:\n", length(target_i)))
      cat(sprintf("      - Found %d G-C pairs.\n", sum(is_gc_pair)))
      cat(sprintf("      - Found %d A-U pairs.\n", sum(is_au_pair)))
  }

  # --- Generate Destabilizing Mutations (GC -> AU) ---
  destab_seq <- NA_character_
  destab_muts <- NA_character_
  if (any(is_gc_pair)) {
    gc_i <- target_i[is_gc_pair]
    gc_j <- target_j[is_gc_pair]
    
    mut_vec <- seq_vec
    mut_vec[gc_i] <- ifelse(seq_vec[gc_i] == "G", "A", "U")
    mut_vec[gc_j] <- ifelse(seq_vec[gc_j] == "C", "U", "A")
    destab_seq <- paste0(mut_vec, collapse = "")
    destab_muts <- paste0(seq_vec[c(gc_i, gc_j)], c(gc_i, gc_j), mut_vec[c(gc_i, gc_j)], collapse = ";")
    if (verbose) cat(sprintf("    [DEBUG] Generated DESTABILIZING mutant with changes: %s\n", destab_muts))
  }
  
  # --- Generate Stabilizing Mutations (AU -> GC) ---
  stab_seq <- NA_character_
  stab_muts <- NA_character_
  if (any(is_au_pair)) {
    au_i <- target_i[is_au_pair]
    au_j <- target_j[is_au_pair]

    mut_vec <- seq_vec
    mut_vec[au_i] <- ifelse(seq_vec[au_i] == "A", "G", "C")
    mut_vec[au_j] <- ifelse(seq_vec[au_j] == "U", "C", "G")
    stab_seq <- paste0(mut_vec, collapse = "")
    stab_muts <- paste0(seq_vec[c(au_i, au_j)], c(au_i, au_j), mut_vec[c(au_i, au_j)], collapse = ";")
    if (verbose) cat(sprintf("    [DEBUG] Generated STABILIZING mutant with changes: %s\n", stab_muts))
  }

  final_mutants <- bind_rows(
    tibble(type_suffix = "destabilizing", mutant_seq = destab_seq, mutations = destab_muts),
    tibble(type_suffix = "stabilizing", mutant_seq = stab_seq, mutations = stab_muts)
  ) %>% filter(!is.na(mutant_seq))
  
  if (verbose) cat(sprintf("    [DEBUG] Returning %d mutant(s) from this function call.\n", nrow(final_mutants)))

  return(final_mutants)
}


# --- Main Execution Block ---
main <- function() {
  cat("--- MFE Mutagenesis and Evaluation Script Started ---\n\n")

  # --- 2. Parse Command-Line Arguments ---
  option_list <- list(
    make_option(c("--input_dir"), type = "character", help = "Directory of .fold files"),
    make_option(c("--output_dir"), type = "character", help = "Directory for final CSV results"),
    make_option(c("--position_zero"), type = "integer", help = "Central reference position (1-based)"),
    make_option(c("--verbose"), action="store_true", default=FALSE, help="Print detailed debugging information")
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
  temp_dir <- file.path(opt$output_dir, paste0("mfe_mutagenesis_run_", timestamp))
  dir.create(temp_dir, showWarnings = FALSE)

  # --- Define mutation windows and protection zone ---
  protect_indices <- (opt$position_zero - 2):(opt$position_zero + 2)
  upstream_window <- (opt$position_zero - 10):(opt$position_zero - 1)
  downstream_window <- (opt$position_zero + 1):(opt$position_zero + 15)
  total_window <- c(upstream_window, downstream_window)
  local_identity_window_size <- 15
  
  # --- 3. Process Input Files and Generate All Mutants ---
  cat("--- Phase 1: Parsing input files and generating mutants ---\n")
  
  fold_files <- list.files(opt$input_dir, pattern = "\\.fold$", full.names = TRUE)
  if (length(fold_files) == 0) stop("No .fold files found in the input directory.")
  
  original_data <- map_dfr(fold_files, ~parse_fold_file(.x, verbose = opt$verbose)) %>%
    mutate(original_pairs = map(original_db, ~db_to_pairs(.x, verbose = opt$verbose)))
  
  mutants_list <- pmap(original_data, function(...) {
    row <- list(...)
    cat(sprintf("\n[INFO] Processing file: %s\n", row$original_filename))
    
    if (opt$verbose) cat("  [INFO] Analyzing UPSTREAM window...\n")
    mut_up <- generate_mfe_mutations(row$sequence, row$original_pairs, upstream_window, protect_indices, opt$verbose) %>% mutate(window_name="upstream")
    
    if (opt$verbose) cat("  [INFO] Analyzing DOWNSTREAM window...\n")
    mut_down <- generate_mfe_mutations(row$sequence, row$original_pairs, downstream_window, protect_indices, opt$verbose) %>% mutate(window_name="downstream")
    
    if (opt$verbose) cat("  [INFO] Analyzing TOTAL window...\n")
    mut_total <- generate_mfe_mutations(row$sequence, row$original_pairs, total_window, protect_indices, opt$verbose) %>% mutate(window_name="total")
    
    all_window_muts <- bind_rows(mut_up, mut_down, mut_total)
    
    if (nrow(all_window_muts) == 0) {
      if(opt$verbose) cat("[INFO] No valid mutations were generated for this file.\n")
      return(all_window_muts)
    }
    
    all_window_muts %>%
      mutate(
        type = paste(window_name, type_suffix, sep = "_"),
        original_filename = row$original_filename,
        .before = 1
      ) %>%
      select(-window_name, -type_suffix)
  })
  
  mutants_df <- bind_rows(mutants_list)
  if (nrow(mutants_df) == 0) {
      cat("\nNo valid mutations could be generated for any input files. Exiting.\n")
      return()
  }

  mutants_df <- mutants_df %>%
    mutate(fasta_header = paste(str_remove(original_filename, "\\.fold$"), type, sep = "_"))
  
  cat(paste("\nGenerated a total of", nrow(mutants_df), "mutant sequences.\n\n"))

  # --- 4. Run RNAfold in Batch ---
  # ... (rest of the script is unchanged) ...
  cat("--- Phase 2: Running RNAfold on all mutants ---\n")
  
  multi_fasta_path <- file.path(temp_dir, "all_mutants.fa")
  sequences_set <- RNAStringSet(mutants_df$mutant_seq)
  names(sequences_set) <- mutants_df$fasta_header
  writeXStringSet(sequences_set, multi_fasta_path, width = 200)
  
  rnafold_out_dir <- file.path(temp_dir, "rnafold_outputs")
  dir.create(rnafold_out_dir)
  system2(run_rnafold_script, args = c(multi_fasta_path, rnafold_out_dir), stdout = !opt$verbose, stderr = !opt$verbose)
  cat("RNAfold script completed.\n\n")

  # --- 5. Evaluate Mutation Success ---
  cat("--- Phase 3: Evaluating mutations ---\n")

  results <- mutants_df %>% left_join(original_data, by = "original_filename")

  get_mutant_fold_data <- function(header) {
    path <- file.path(rnafold_out_dir, paste0(header, ".fold"))
    if (!file.exists(path)) return(tibble(mutant_db = NA_character_, mutant_mfe = NA_real_))
    parse_fold_file(path) %>% select(mutant_db = original_db, mutant_mfe = original_mfe)
  }
  
  mutant_fold_data <- map_dfr(results$fasta_header, get_mutant_fold_data)
  results <- bind_cols(results, mutant_fold_data) %>% filter(!is.na(mutant_db))
  
  # --- Calculate metrics ---
  cat("Calculating MFE change and structural identity metrics...\n")
  results <- results %>%
    mutate(
      mfe_change = mutant_mfe - original_mfe,
      overall_identity = 1 - (stringdist(original_db, mutant_db, method = "hamming") / nchar(original_db)),
      
      local_win_start = pmax(1, opt$position_zero - local_identity_window_size),
      local_win_end = pmin(nchar(original_db), opt$position_zero + local_identity_window_size),
      local_orig_db = str_sub(original_db, local_win_start, local_win_end),
      local_mut_db = str_sub(mutant_db, local_win_start, local_win_end),
      local_identity = 1 - (stringdist(local_orig_db, local_mut_db, method = "hamming") / nchar(local_orig_db)),
      
      motif_win_start = pmax(1, opt$position_zero - 2),
      motif_win_end = pmin(nchar(original_db), opt$position_zero + 2),
      motif_orig_db = str_sub(original_db, motif_win_start, motif_win_end),
      motif_mut_db = str_sub(mutant_db, motif_win_start, motif_win_end),
      motif_identity = 1 - (stringdist(motif_orig_db, motif_mut_db, method = "hamming") / nchar(motif_orig_db))
    )

  # --- Filter for significant MFE change ---
  significant_results <- results %>% filter(abs(mfe_change) > 5)

  # --- 6. Save Final Outputs ---
  cat("\n--- Phase 4: Saving final results ---\n")
  
  final_cols <- c(
    "original_filename", "type", "mutations", "mfe_change", 
    "overall_identity", "local_identity", "motif_identity",
    "original_mfe", "mutant_mfe", "original_db", "mutant_db",
    "sequence", "mutant_seq", "fasta_header"
  )
  
  all_results_path <- file.path(opt$output_dir, "all_mfe_mutants.csv")
  significant_results_path <- file.path(opt$output_dir, "significant_mfe_mutants.csv")
  
  write.csv(results %>% select(any_of(final_cols)), all_results_path, row.names = FALSE)
  write.csv(significant_results %>% select(any_of(final_cols)), significant_results_path, row.names = FALSE)
  
  cat(paste("Saved all results to:", all_results_path, "\n"))
  cat(paste("Saved significant MFE change (|>5|) results to:", significant_results_path, "\n"))
  cat(paste("\nIntermediate files are located in:", temp_dir, "\n"))
  cat("--- Script Finished ---\n")
}

# Run the main function
main()