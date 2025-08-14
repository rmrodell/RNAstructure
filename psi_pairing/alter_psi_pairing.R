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
              help = "Required. The 1-based nucleotide position whose pairing status you want to flip.", metavar = "INT")
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
#   output_dir = "/scratch/users/rodell/psi_pairing/makeandbreakpair/test1",
#   target_pos = 59
# )


# ---- Dynamic Configuration ----
# Similar structure to the original script.
script_dir <- this.path::this.dir()
CONFIG <- list(
  target_pos = opt$target_pos,
  
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
    dot_bracket <- strsplit(lines[3], " ")[[1]][1]
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

# Scenario 1: Make a paired site unpaired.
mutate_to_break_pair <- function(sequence, partner_pos, target_base) {
  seq_vec <- strsplit(sequence, "")[[1]]
  original_base <- seq_vec[partner_pos]
  
  mutated_base <- target_base
  
  seq_vec[partner_pos] <- mutated_base
  
  mutated_seq <- paste(seq_vec, collapse = "")
  mutation_info <- sprintf("%s%d%s", original_base, partner_pos, mutated_base)
  
  return(list(
    sequence = mutated_seq,
    info_string = mutation_info,
    strategy_used = "break_pair" # Add this line
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
#' Tries to pair a target by extending the nearest stem by 3 bp.
attempt_extend_stem <- function(org_data, config) {
  seq_len <- nchar(org_data$seq)
  target_pos <- config$target_pos
  SEARCH_WINDOW <- 5 # How far to look for a paired region.
  MAX_GAP_SIZE <- 3  # Max unpaired bases between stem and target.
  
  # Search for the nearest paired base
  for (i in 1:SEARCH_WINDOW) {
    positions_to_check <- c(target_pos - i, target_pos + i)
    
    for (nearest_paired_base in positions_to_check) {
      if (nearest_paired_base < 1 || nearest_paired_base > seq_len) next
      
      if (substr(org_data$dot_bracket, nearest_paired_base, nearest_paired_base) != ".") {
        
        partner_of_nearest <- find_pairing_partner(org_data$dot_bracket, nearest_paired_base)
        if(is.na(partner_of_nearest)) next
        
        # --- NEW: Generalize for gaps ---
        gap_size <- abs(target_pos - nearest_paired_base) - 1
        if (gap_size < 0 || gap_size > MAX_GAP_SIZE) next
        
        # Determine the full range of positions to pair
        if (target_pos > nearest_paired_base) { # Upstream stem: ...))...t...
          # Positions on the target side of the new helix
          positions_to_pair_1 <- (nearest_paired_base + 1):target_pos
          # Positions on the opposite strand that will be mutated
          positions_to_pair_2 <- rev((partner_of_nearest - gap_size - 1):(partner_of_nearest - 1))
        } else { # Downstream stem: ...t...((...
          positions_to_pair_1 <- rev((target_pos):(nearest_paired_base - 1))
          positions_to_pair_2 <- (partner_of_nearest + 1):(partner_of_nearest + 1 + gap_size)
        }
        
        # --- Sanity Checks ---
        all_new_positions <- c(positions_to_pair_1, positions_to_pair_2)
        if (any(all_new_positions < 1) || any(all_new_positions > seq_len)) next
        
        is_unpaired <- sapply(all_new_positions, function(p) substr(org_data$dot_bracket, p, p) == ".")
        if (!all(is_unpaired)) next
        
        # --- If we get here, the extension is valid. Perform all mutations. ---
        seq_vec <- strsplit(org_data$seq, "")[[1]]
        mutation_infos <- list()
        
        for (j in 1:length(positions_to_pair_1)) {
          pos1 <- positions_to_pair_1[j]
          pos2 <- positions_to_pair_2[j]
          
          base_to_get_comp_from <- seq_vec[pos1]
          original_base_to_mutate <- seq_vec[pos2]
          
          complement <- get_complement(base_to_get_comp_from)
          if(is.na(complement)) next # Skip if we have an 'N' or other weird base
          
          # Apply mutation
          seq_vec[pos2] <- complement
          
          # Record info
          mutation_infos[[j]] <- sprintf("%s%d%s", original_base_to_mutate, pos2, complement)
        }
        
        # Return the result
        return(list(
          sequence = paste(seq_vec, collapse = ""),
          info_string = paste(unlist(mutation_infos), collapse = ", "),
          strategy_used = "extend_stem"
        ))
      }
    }
  }
  
  # If the loop finishes without finding a valid extension
  return(NULL)
}

# Strategy 2: Seed
# Tries to pair a target by creating a new hairpin.
# Forms a 2-bp stem creating a stable tetraloop ((..)).
attempt_seed_hairpin <- function(org_data, config) {
  seq_len <- nchar(org_data$seq)
  target_pos <- config$target_pos
  
  # Define positions for a ((..)) tetraloop relative to target_pos
  # target_pos will be the first '(', target_pos+1 the second '('
  # target_pos+3 the first ')', target_pos+4 the second ')'
  pos_to_mutate_1 <- target_pos + 4 # Partner for target_pos
  pos_to_mutate_2 <- target_pos + 3 # Partner for target_pos + 1
  
  # --- Sanity Checks ---
  # 1. Ensure all positions are within the sequence bounds.
  if (pos_to_mutate_1 > seq_len) return(NULL)
  
  # 2. Ensure all involved positions are currently unpaired.
  required_unpaired <- c(target_pos, target_pos + 1, target_pos + 2, target_pos + 3, target_pos + 4)
  for (pos in required_unpaired) {
    if (substr(org_data$dot_bracket, pos, pos) != ".") return(NULL)
  }
  
  # --- Perform Mutation ---
  seq_vec <- strsplit(org_data$seq, "")[[1]]
  
  # Mutation 1: Pair target_pos with pos_to_mutate_1
  base_at_target <- seq_vec[target_pos]
  base_to_mutate_1 <- seq_vec[pos_to_mutate_1]
  complement_1 <- get_complement(base_at_target)
  seq_vec[pos_to_mutate_1] <- complement_1
  
  # Mutation 2: Pair target_pos + 1 with pos_to_mutate_2
  base_at_adjacent <- seq_vec[target_pos + 1]
  base_to_mutate_2 <- seq_vec[pos_to_mutate_2]
  complement_2 <- get_complement(base_at_adjacent)
  seq_vec[pos_to_mutate_2] <- complement_2
  
  # --- Format Output ---
  mutated_seq <- paste(seq_vec, collapse = "")
  mutation_info <- sprintf("%s%d%s, %s%d%s", 
                           base_to_mutate_1, pos_to_mutate_1, complement_1,
                           base_to_mutate_2, pos_to_mutate_2, complement_2)
  
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
analyze_flip_results <- function(final_df, target_pos) {
  if(nrow(final_df) == 0) return(final_df)
  
  final_df %>%
    rowwise() %>%
    mutate(
      new_status = get_pairing_status(dot_bracket_mut, target_pos),
      is_successful = case_when(
        # The only success condition we are currently testing
        original_status == "paired" & new_status == "unpaired" ~ TRUE,
        # Future success condition
        original_status == "unpaired" & new_status == "paired" ~ TRUE, 
        TRUE ~ FALSE
      )
    ) %>%
    ungroup()
}

# ---- Main Workflow ----

main_workflow <- function(config) {
  # 1. Create output directories
  dir.create(config$intermediate_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 2. Process each input file to determine strategy and generate mutations
  cat("Step 1: Planning mutations...\n")
  all_fold_files <- list.files(config$org_fold_dir, pattern = "\\.fold$", full.names = TRUE)
  
  mutation_plan_list <- purrr::map(all_fold_files, function(file) {
    # Read original data
    org_data <- read_fold_file(file)
    if(is.null(org_data)) return(NULL)
    
    # Determine initial status
    original_status <- get_pairing_status(org_data$dot_bracket, config$target_pos)
    
    mutation_result <- NULL
    
    # --- Main Logic Gate ---
    # Apply strategy based on status
    if (original_status == "paired") {
      partner_pos <- find_pairing_partner(org_data$dot_bracket, config$target_pos)
      
      if (is.na(partner_pos)) {
        cat(sprintf("  - WARNING: Could not find partner for paired site %d in %s. Skipping.\n", config$target_pos, org_data$name))
        return(NULL)
      }
      
      # 1. Get the base at the target position
      target_base <- substr(org_data$seq, config$target_pos, config$target_pos)
      
      # 2. Pass it to the updated function
      mutation_result <- mutate_to_break_pair(org_data$seq, partner_pos, target_base)

      #strategy <- "break_pair"
      
      
    } else { # status is "unpaired"
      mutation_result <- find_and_mutate_to_create_pair(org_data, config)
    }
    
    # If a mutation was possible, return a row for our results table
    if (!is.null(mutation_result)) {
      return(
        data.frame(
          name = org_data$name,
          original_seq = org_data$seq,
          original_dot_bracket = org_data$dot_bracket,
          target_pos = config$target_pos,
          original_status = original_status,
          strategy = mutation_result$strategy_used,
          mutated_seq = mutation_result$sequence,
          mutation_info = mutation_result$info_string,
          stringsAsFactors = FALSE
        )
      )
    } else {
      return(NULL) # This sequence could not be mutated
    }
  })
  
  # Combine all successful mutation plans into one data frame
  mutation_plan_df <- bind_rows(mutation_plan_list)
  
  if (nrow(mutation_plan_df) == 0) {
    cat("No sequences were eligible for mutation. Exiting.\n")
    return()
  }
  
  # 3. Save mutated sequences and run RNAfold
  cat("Step 2: Running RNAfold on", nrow(mutation_plan_df), "mutated sequences...\n")
  create_and_save_fasta(mutation_plan_df, "mutated_seq", "name", config$mutated_fasta_path)
  
  system2("bash", 
          args = c(config$rna_fold_script, 
                   config$mutated_fasta_path, 
                   config$mut_fold_dir),
          stdout = TRUE, stderr = TRUE)
  
  # 4. Load new structures and perform final analysis
  cat("Step 3: Analyzing results...\n")
  mutated_fold_data <- read_all_fold_files(config$mut_fold_dir)
  
  if (nrow(mutated_fold_data) == 0) {
    cat("RNAfold did not produce any output files. Cannot complete analysis.\n")
    return()
  }
  
  # Merge original plan with new results
  final_results <- mutation_plan_df %>%
    left_join(mutated_fold_data %>% select(name, dot_bracket_mut = dot_bracket), by = "name")
  
  # Check if the flip was successful
  final_analysis <- analyze_flip_results(final_results)
  
  # 5. Save outputs
  cat("Step 4: Saving final reports...\n")
  write.csv(final_analysis, config$output_csv, row.names = FALSE)
  
  # Generate a simple summary report
  break_pair_results <- filter(final_analysis, strategy == "break_pair")
  extend_stem_results <- filter(final_analysis, strategy == "extend_stem")
  seed_hairpin_results <- filter(final_analysis, strategy == "seed_hairpin")
  
  summary_text <- c(
    "===== Pairing Status Flip: Summary Report =====",
    "",
    paste("Target Position:", config$target_pos),
    paste("Total files processed:", length(all_fold_files)),
    "---",
    "Strategy 1: Unpair a Paired Site",
    paste("  - Sequences attempted:", nrow(break_pair_results)),
    paste("  - Successfully flipped to unpaired:", sum(break_pair_results$is_successful, na.rm = TRUE)),
    "---",
    "Strategy 2: Pair an Unpaired Site",
    "  Sub-strategy: Extend Nearest Stem",
    paste("    - Sequences attempted:", nrow(extend_stem_results)),
    paste("    - Successfully flipped to paired:", sum(extend_stem_results$is_successful, na.rm = TRUE)),
    "  Sub-strategy: Seed New Hairpin (Fallback)",
    paste("    - Sequences attempted:", nrow(seed_hairpin_results)),
    paste("    - Successfully flipped to paired:", sum(seed_hairpin_results$is_successful, na.rm = TRUE)),
    "---",
    paste("Total mutations attempted:", nrow(final_analysis)),
    paste("Total successful flips (any type):", sum(final_analysis$is_successful, na.rm = TRUE)),
    "============================================="
  )
  writeLines(summary_text, config$output_summary)
  
  cat("Workflow complete. Outputs are in:", config$output_dir, "\n")
}

# ---- Execution ----
main_workflow(CONFIG)
