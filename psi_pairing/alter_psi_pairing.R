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
cat("Running in interactive mode. Using manual 'opt' list for debugging.\n")

opt <- list(
  fold_dir = "/scratch/users/rodell/RNAfold_psipos",
  output_dir = "/scratch/users/rodell/psi_pairing/breakpair/repulsive59",
  target_pos = 59
)


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
# mutate_to_break_pair <- function(sequence, partner_pos, target_base) {
#   seq_vec <- strsplit(sequence, "")[[1]]
#   original_base <- seq_vec[partner_pos]
#   
#   # # Mutation strategy: purine -> purine, pyrimidine -> pyrimidine
#   # mutated_base <- case_when(
#   #   original_base == "A" ~ "G",
#   #   original_base == "G" ~ "A",
#   #   original_base == "C" ~ "U",
#   #   original_base == "U" ~ "C",
#   #   TRUE ~ NA_character_ # Handle unexpected characters like 'N'
#   # )
#   
#   # New strategy: The mutated base is now a copy of the target base.
#   mutated_base <- target_base 
#   
#   # If the original base was not A, U, G, or C, we can't mutate it.
#   if (is.na(mutated_base)) {
#     warning(sprintf("Cannot mutate unrecognized base '%s' at position %d.", original_base, partner_pos))
#     return(NULL) # Return NULL to indicate failure
#   }
#   
#   seq_vec[partner_pos] <- mutated_base
#   
#   mutated_seq <- paste(seq_vec, collapse = "")
#   mutation_info <- sprintf("%s%d%s", original_base, partner_pos, mutated_base)
#   
#   return(list(
#     sequence = mutated_seq,
#     info_string = mutation_info
#   ))
# }

mutate_to_break_pair <- function(sequence, partner_pos, target_base) {
  seq_vec <- strsplit(sequence, "")[[1]]
  original_base <- seq_vec[partner_pos]
  
  # New strategy: The mutated base is now a copy of the target base.
  mutated_base <- target_base
  
  seq_vec[partner_pos] <- mutated_base
  
  mutated_seq <- paste(seq_vec, collapse = "")
  mutation_info <- sprintf("%s%d%s", original_base, partner_pos, mutated_base)
  
  return(list(
    sequence = mutated_seq,
    info_string = mutation_info
  ))
}

# Scenario 2: Make an unpaired site paired. This is the hardest part.
# This function would contain the "zip it up" logic.
find_and_mutate_to_create_pair <- function(sequence, dot_bracket, target_pos) {
  # Strategy:
  # 1. Identify the base at target_pos (e.g., 'G').
  # 2. Find the desired complementary base (e.g., 'C').
  # 3. Search for a suitable partner position. A good heuristic is to look for the
  #    nearest existing stem and try to extend it.
  # 4. For example, if target_pos is in a loop, find the closest base of the
  #    enclosing stem. Let's say the stem is from i to j. If target_pos is k,
  #    and k is near i, a good candidate partner is near j.
  # 5. Let's say we identify position `p` as the best candidate partner.
  # 6. Mutate the base at `p` to the complementary base (e.g., 'C').
  # 7. **Crucial:** To "zip it up" and make this new pair stable, you might also
  #    need to introduce a supporting mutation. For example, mutate `p-1` to
  #    be complementary to `target_pos+1`. This creates a stable 2-bp stem.
  # Returns: list(mutated_seq, mutation_info) or NULL if no good candidate is found.
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
      # 
      # mutation_result <- mutate_to_break_pair(org_data$seq, partner_pos)
      strategy <- "break_pair"
      
      
    } else { # status is "unpaired"
      #mutation_result <- find_and_mutate_to_create_pair(org_data$seq, org_data$dot_bracket, config$target_pos)
      
      # --- PLACEHOLDER for the harder case ---
      #cat(sprintf("  - INFO: Target position %d in %s is already unpaired. Logic to create a pair is not yet implemented. Skipping.\n", config$target_pos, org_data$name))
      return(NULL) # Skip this sequence for now
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
          strategy = strategy,
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
  total_processed <- length(all_fold_files)
  total_attempted <- nrow(mutation_plan_df)
  total_successful <- sum(final_analysis$is_successful, na.rm = TRUE)
  
  summary_text <- c(
    "===== Pairing Status Flip: Summary Report =====",
    "",
    paste("Target Position:", config$target_pos),
    paste("Total files processed:", total_processed),
    "",
    "--- Strategy: Unpair a Paired Site ---",
    paste("Sequences where target was initially paired:", total_attempted),
    paste("Successfully flipped to unpaired:", total_successful),
    paste("Failed to flip:", total_attempted - total_successful),
    "",
    "--- Strategy: Pair an Unpaired Site ---",
    "Logic for this strategy is not yet implemented.",
    "============================================="
  )
  writeLines(summary_text, config$output_summary)
  
  cat("Workflow complete. Outputs are in:", config$output_dir, "\n")
}

# ---- Execution ----
main_workflow(CONFIG)
