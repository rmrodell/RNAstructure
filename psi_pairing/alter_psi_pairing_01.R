#!/usr/bin/env Rscript

# ---- Libraries ----
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(this.path))

# ---- Argument Parsing ----
option_list <- list(
  make_option(c("-f", "--fold_dir"), type = "character", default = NULL,
              help = "Required. Path to the input directory containing original .fold files.", metavar = "DIR"),
  
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Required. Path to the directory where all output files will be saved.", metavar = "DIR"),
  
  make_option(c("-t", "--target_pos"), type = "integer", default = NULL,
              help = "Required. The 1-based start of a 2-nucleotide block to flip.", metavar = "INT")
)

opt_parser <- OptionParser(option_list = option_list, description = "Mutates RNA to flip the pairing status of a target nucleotide.")
opt <- parse_args(opt_parser)

# Check for required arguments
if (is.null(opt$fold_dir) || is.null(opt$output_dir) || is.null(opt$target_pos)) {
  print_help(opt_parser)
  stop("All three arguments (--fold_dir, --output_dir, --target_pos) must be supplied.", call. = FALSE)
}

# ---- Interactive Debugging Setup ----
# If running in RStudio, manually define the 'opt' list with your test arguments.
# cat("Running in interactive mode. Using manual 'opt' list for debugging.\n")
# 
# opt <- list(
#   fold_dir = "/scratch/users/rodell/RNAfold_psipos",
#   output_dir = "/scratch/users/rodell/psi_pairing/makeandbreakpair01/test1",
#   target_pos = 59
# )


# ---- Dynamic Configuration ----
# Similar structure to the original script.
script_dir <- this.path::this.dir()
CONFIG <- list(
  target_pos1 = opt$target_pos,
  target_pos2 = opt$target_pos + 1,
  
  org_fold_dir    = opt$fold_dir,
  
  # Paths to helper scripts
  rna_fold_script = file.path(script_dir, "run_RNAfold.sh"), 
  
  # Paths for intermediate and final files
  intermediate_dir   = file.path(opt$output_dir, "intermediate"),
  mutated_fasta_path = file.path(opt$output_dir, "intermediate", "mutated_sequences.fa"),
  mut_fold_dir       = file.path(opt$output_dir, "intermediate", "mutated_structures"),
  
  output_csv         = file.path(opt$output_dir, "pairing_status_flip_analysis.csv"),
  output_summary     = file.path(opt$output_dir, "summary_report.txt")
)



# ---- File Operations ----
read_fold_file <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  if (length(lines) >= 3) {
    name <- sub(">", "", lines[1])
    seq <- lines[2]
    dot_bracket <- strsplit(lines[6], " ")[[1]][1]
    return(data.frame(name = name, seq = seq, dot_bracket = dot_bracket, stringsAsFactors = FALSE))
  } else {
    warning(paste("File", file_path, "does not have the expected format"))
    return(NULL)
  }
}

read_all_fold_files <- function(fold_dir) {
  fold_files <- list.files(fold_dir, pattern = "\\.fold$", full.names = TRUE)
  if (length(fold_files) == 0) {
    stop(paste("No .fold files found in the specified directory:", fold_dir))
  }
  fold_data <- do.call(rbind, lapply(fold_files, read_fold_file))
  return(fold_data)
}

create_and_save_fasta <- function(df, seq_column, name_column, fasta_path) {
  if (nrow(df) == 0) {
    cat("Warning: No sequences to write to FASTA file:", fasta_path, "\n")
    return()
  }
  named_sequences <- setNames(df[[seq_column]], df[[name_column]])
  dna_sequences <- gsub("U", "T", named_sequences)
  mutated_sequences <- DNAStringSet(dna_sequences)
  writeXStringSet(mutated_sequences, fasta_path, width = 10000)
}


# ---- Core Logic Functions ----

get_pairing_status <- function(dot_bracket, position) {
  char <- substr(dot_bracket, position, position)
  if (char == ".") {
    return("unpaired")
  } else if (char %in% c("(", ")", "[", "]", "{", "}")) {
    return("paired")
  } else {
    return("unknown")
  }
}

get_combined_pairing_status <- function(dot_bracket, pos1, pos2) {
  s1 <- get_pairing_status(dot_bracket, pos1)
  s2 <- get_pairing_status(dot_bracket, pos2)
  paste(s1, s2, sep = "-")
}

# Finds the partner of a paired nucleotide.
find_pairing_partner <- function(dot_bracket, target_pos) {
  char_at_target <- substr(dot_bracket, target_pos, target_pos)
  
  if (char_at_target == '.') return(NA)
  
  chars <- strsplit(dot_bracket, "")[[1]]
  
  if (char_at_target == '(') {
    balance <- 1
    for (i in (target_pos + 1):length(chars)) {
      if (chars[i] == '(') balance <- balance + 1
      if (chars[i] == ')') balance <- balance - 1
      if (balance == 0) return(i)
    }
  } else if (char_at_target == ')') {
    balance <- 1
    for (i in (target_pos - 1):1) {
      if (chars[i] == ')') balance <- balance + 1
      if (chars[i] == '(') balance <- balance - 1
      if (balance == 0) return(i)
    }
  }
  # Note: This logic can be extended for other bracket types like [] or {} if needed.
  
  return(NA) # Return NA if no partner is found (malformed dot-bracket)
}

# Scenario 1: Make paired/mixed sites unpaired.
mutate_to_break_pairs <- function(org_data, config) {
  seq_vec <- strsplit(org_data$seq, "")[[1]]
  positions_to_check <- c(config$target_pos1, config$target_pos2)
  mutation_infos <- list()
  
  for (pos in positions_to_check) {
    if (get_pairing_status(org_data$dot_bracket, pos) == "paired") {
      partner_pos <- find_pairing_partner(org_data$dot_bracket, pos)
      if (!is.na(partner_pos)) {
        original_base <- seq_vec[partner_pos]
        target_base_for_disruption <- seq_vec[pos]
        seq_vec[partner_pos] <- target_base_for_disruption
        mutation_infos[[length(mutation_infos) + 1]] <- sprintf("%s%d%s", original_base, partner_pos, target_base_for_disruption)
      }
    }
  }
  
  if (length(mutation_infos) == 0) return(NULL)
  
  return(list(
    sequence = paste(seq_vec, collapse = ""),
    info_string = paste(unlist(mutation_infos), collapse = ", "),
    strategy_used = "break_pairs"
  ))
}

# Scenario 2: Make an unpaired site paired. This is the hardest part.
get_complement <- function(base) {
  case_when(
    base == "A" ~ "U",
    base == "U" ~ "A",
    base == "G" ~ "C",
    base == "C" ~ "G",
    TRUE ~ NA_character_
  )
}

# Strategy 1: Zip
# Tries to pair a target by extending the nearest stem by 3 bp.
attempt_extend_stem <- function(org_data, config) {
  seq_len <- nchar(org_data$seq)
  SEARCH_WINDOW <- 5
  
  for (i in 1:SEARCH_WINDOW) {
    for (nearest_paired_base in c(config$target_pos1 - i, config$target_pos2 + i)) {
      if (nearest_paired_base < 1 || nearest_paired_base > seq_len || get_pairing_status(org_data$dot_bracket, nearest_paired_base) == "unpaired") next
      
      partner_of_nearest <- find_pairing_partner(org_data$dot_bracket, nearest_paired_base)
      if(is.na(partner_of_nearest)) next
      
      if (config$target_pos1 > nearest_paired_base) { # Upstream stem
        positions_to_pair_1 <- (nearest_paired_base + 1):config$target_pos2
        positions_to_pair_2 <- rev((partner_of_nearest - length(positions_to_pair_1)):(partner_of_nearest - 1))
      } else { # Downstream stem
        positions_to_pair_1 <- rev(config$target_pos1:(nearest_paired_base - 1))
        positions_to_pair_2 <- (partner_of_nearest + 1):(partner_of_nearest + length(positions_to_pair_1))
      }
      
      all_new_positions <- c(positions_to_pair_1, positions_to_pair_2)
      if (any(all_new_positions < 1) || any(all_new_positions > seq_len) || any(sapply(all_new_positions, function(p) get_pairing_status(org_data$dot_bracket, p) != "unpaired"))) next
      
      seq_vec <- strsplit(org_data$seq, "")[[1]]
      mutation_infos <- map2_chr(positions_to_pair_1, positions_to_pair_2, function(pos1, pos2) {
        complement <- get_complement(seq_vec[pos1])
        original_base_to_mutate <- seq_vec[pos2]
        if(is.na(complement)) return(NA_character_)
        seq_vec[pos2] <<- complement
        sprintf("%s%d%s", original_base_to_mutate, pos2, complement)
      })
      
      if(any(is.na(mutation_infos))) next
      
      return(list(
        sequence = paste(seq_vec, collapse = ""),
        info_string = paste(mutation_infos, collapse = ", "),
        strategy_used = "extend_stem"
      ))
    }
  }
  return(NULL)
}

# Strategy 2: Seed
# Tries to pair a target by creating a new hairpin.
# Forms a 2-bp stem creating a stable tetraloop ((..)).
attempt_seed_hairpin <- function(org_data, config) {
  seq_len <- nchar(org_data$seq)
  
  # Use the new, consistent variable names
  pos_to_mutate_1 <- config$target_pos1 + 4 
  pos_to_mutate_2 <- config$target_pos1 + 3
  
  # Sanity Checks
  if (pos_to_mutate_1 > seq_len) return(NULL)
  
  required_unpaired <- config$target_pos1:(config$target_pos1 + 4)
  for (pos in required_unpaired) {
    if (substr(org_data$dot_bracket, pos, pos) != ".") return(NULL)
  }
  
  # Perform Mutation
  seq_vec <- strsplit(org_data$seq, "")[[1]]
  
  # Use target_pos1 and target_pos2 explicitly
  base_at_target1 <- seq_vec[config$target_pos1]
  base_at_target2 <- seq_vec[config$target_pos2] # This is target_pos1 + 1
  
  complement_1 <- get_complement(base_at_target1)
  complement_2 <- get_complement(base_at_target2)
  
  # Check for non-standard bases that can't be complemented
  if (is.na(complement_1) || is.na(complement_2)) return(NULL)
  
  seq_vec[pos_to_mutate_1] <- complement_1
  seq_vec[pos_to_mutate_2] <- complement_2
  
  # Format Output
  mutated_seq <- paste(seq_vec, collapse = "")
  mutation_info <- sprintf("%s%d%s, %s%d%s", 
                           substr(org_data$seq, pos_to_mutate_1, pos_to_mutate_1), pos_to_mutate_1, complement_1,
                           substr(org_data$seq, pos_to_mutate_2, pos_to_mutate_2), pos_to_mutate_2, complement_2)
  
  return(list(
    sequence = mutated_seq,
    info_string = mutation_info,
    strategy_used = "seed_hairpin"
  ))
}

find_and_mutate_to_create_pair <- function(org_data, config) {
  
  # Strategy 1: Try to extend the nearest stem (high priority)
  result <- attempt_extend_stem(org_data, config)
  
  if (!is.null(result)) {
    return(result)
  }
  
  # Strategy 2: If extending failed, try to seed a new hairpin (fallback)
  result <- attempt_seed_hairpin(org_data, config)
  
  return(result)
}

# ---- Analysis Functions ----

# Compares original and new status and determines success.
analyze_flip_results <- function(final_df, config) {
  if(nrow(final_df) == 0) return(final_df)
  
  final_df %>%
    rowwise() %>%
    mutate(
      new_status = get_combined_pairing_status(dot_bracket_mut, config$target_pos1, config$target_pos2),
      is_successful = case_when(
        grepl("paired", original_status) & new_status == "unpaired-unpaired" ~ TRUE,
        original_status == "unpaired-unpaired" & new_status == "paired-paired" ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    ungroup()
}

# ---- Main Workflow ----
main_workflow <- function(config) {
  dir.create(config$intermediate_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat("Step 1: Planning mutations...\n")
  all_fold_files <- list.files(config$org_fold_dir, pattern = "\\.fold$", full.names = TRUE)
  
  mutation_plan_list <- purrr::map(all_fold_files, function(file) {
    org_data <- read_fold_file(file)
    if(is.null(org_data) || nchar(org_data$seq) < config$target_pos2) return(NULL)
    
    original_status <- get_combined_pairing_status(org_data$dot_bracket, config$target_pos1, config$target_pos2)
    
    mutation_result <- NULL
    if (original_status == "unpaired-unpaired") {
      mutation_result <- find_and_mutate_to_create_pair(org_data, config)
    } else if (grepl("paired", original_status)) {
      mutation_result <- mutate_to_break_pairs(org_data, config)
    }
    
    if (!is.null(mutation_result)) {
      return(data.frame(
        name = org_data$name,
        original_seq = org_data$seq,
        original_dot_bracket = org_data$dot_bracket,
        target_pos = paste(config$target_pos1, config$target_pos2, sep = "-"),
        original_status = original_status,
        strategy = mutation_result$strategy_used,
        mutated_seq = mutation_result$sequence,
        mutation_info = mutation_result$info_string,
        stringsAsFactors = FALSE
      ))
    }
    return(NULL)
  })
  
  mutation_plan_df <- bind_rows(mutation_plan_list)
  if (nrow(mutation_plan_df) == 0) {
    cat("No sequences were eligible for mutation. Exiting.\n")
    return()
  }
  
  cat("Step 2: Running RNAfold on", nrow(mutation_plan_df), "mutated sequences...\n")
  create_and_save_fasta(mutation_plan_df, "mutated_seq", "name", config$mutated_fasta_path)
  system2("bash", args = c(config$rna_fold_script, config$mutated_fasta_path, config$mut_fold_dir), stdout = TRUE, stderr = TRUE)
  
  cat("Step 3: Analyzing results...\n")
  mutated_fold_data <- read_all_fold_files(config$mut_fold_dir)
  if (nrow(mutated_fold_data) == 0) {
    cat("RNAfold did not produce any output files. Cannot complete analysis.\n")
    return()
  }
  
  final_results <- mutation_plan_df %>%
    left_join(mutated_fold_data %>% select(name, dot_bracket_mut = dot_bracket), by = "name")
  
  final_analysis <- analyze_flip_results(final_results, config)
  
  cat("Step 4: Saving final reports...\n")
  
  # --- NEW: Enhanced Reporting Calculations ---
  # Total attempted for each initial status
  total_uu <- sum(final_analysis$original_status == "unpaired-unpaired")
  total_pp <- sum(final_analysis$original_status == "paired-paired")
  total_pu <- sum(final_analysis$original_status == "paired-unpaired")
  total_up <- sum(final_analysis$original_status == "unpaired-paired")
  
  # Successful flips for each initial status
  successful_uu <- sum(final_analysis$original_status == "unpaired-unpaired" & final_analysis$is_successful, na.rm = TRUE)
  successful_pp <- sum(final_analysis$original_status == "paired-paired" & final_analysis$is_successful, na.rm = TRUE)
  successful_pu <- sum(final_analysis$original_status == "paired-unpaired" & final_analysis$is_successful, na.rm = TRUE)
  successful_up <- sum(final_analysis$original_status == "unpaired-paired" & final_analysis$is_successful, na.rm = TRUE)
  
  write.csv(final_analysis, config$output_csv, row.names = FALSE)
  
  summary_text <- c(
    "===== Pairing Status Flip: Summary Report =====",
    "",
    paste("Target Positions:", config$target_pos1, "and", config$target_pos2),
    paste("Total sequences processed:", length(all_fold_files)),
    paste("Total sequences with mutation attempted:", nrow(final_analysis)),
    "---",
    "Detailed Flip Analysis by Initial State:",
    "",
    "  Initial State: unpaired-unpaired",
    paste("    - Goal: Create a paired-paired block."),
    paste("    - Sequences Attempted:", total_uu),
    paste("    - Successfully Flipped to paired-paired:", successful_uu),
    "",
    "  Initial State: paired-paired",
    paste("    - Goal: Create an unpaired-unpaired block."),
    paste("    - Sequences Attempted:", total_pp),
    paste("    - Successfully Flipped to unpaired-unpaired:", successful_pp),
    "",
    "  Initial State: paired-unpaired",
    paste("    - Goal: Create an unpaired-unpaired block."),
    paste("    - Sequences Attempted:", total_pu),
    paste("    - Successfully Flipped to unpaired-unpaired:", successful_pu),
    "",
    "  Initial State: unpaired-paired",
    paste("    - Goal: Create an unpaired-unpaired block."),
    paste("    - Sequences Attempted:", total_up),
    paste("    - Successfully Flipped to unpaired-unpaired:", successful_up),
    "---",
    paste("Total successful flips (any type):", sum(final_analysis$is_successful, na.rm = TRUE)),
    "============================================="
  )
  writeLines(summary_text, config$output_summary)
  
  cat("Workflow complete. Outputs are in:", config$output_dir, "\n")
}

# ---- Execution ----
main_workflow(CONFIG)
