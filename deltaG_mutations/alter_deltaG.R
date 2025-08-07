#!/usr/bin/env Rscript

# ---- Libraries ----
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(this.path))

# ---- Argument Parsing with optparse ----
# 1. Define the list of options the script will accept.
option_list <- list(
  make_option(c("-f", "--fold_dir"), type = "character", default = NULL,
              help = "Required. Path to the input directory containing original .fold files.", metavar = "DIR"),
  
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Required. Path to the directory where all output files will be saved.", metavar = "DIR"),
  
  make_option(c("-p", "--psipos"), type = "integer", default = 59,
              help = "Position of the pseudouridine site. [default: %default]", metavar = "INT"),
  
  make_option(c("-w", "--window"), type = "integer", default = 5,
              help = "Window size around psipos to look for G-C pairs. [default: %default]", metavar = "INT"),
  
  make_option(c("-d", "--deltaG_threshold"), type = "double", default = 5.0,
              help = "Threshold for significant deltaG change. [default: %default]", metavar = "NUM"),
  
  make_option(c("-m", "--mode"), type = "character", default = "destabilize",
              help = paste("Operational mode:",
                           "'destabilize' (GC->AU/GU) aims for a less negative (higher) delta G.",
                           "'stabilize' (AU->GC/GU) aims for a more negative (lower) delta G.",
                           "[default: %default]"), metavar = "MODE"),
  
  make_option(c("-s", "--mutation_strategy"), type = "character", default = "AU",
              help = paste("Mutation strategy. For 'destabilize' mode: 'AU' or 'GU'.",
                           "For 'stabilize' mode: 'GC' or 'GU'.",
                           "[default: %default]"), metavar = "STRATEGY")
)

# 2. Create the parser and parse the arguments.
opt_parser <- OptionParser(option_list = option_list, description = "Perform mutation analysis to increase RNA deltaG and compare structures.")
opt <- parse_args(opt_parser)

# 3. Manually check for required arguments, as optparse doesn't have a 'required=TRUE' flag.
if (is.null(opt$fold_dir) || is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("Both --fold_dir and --output_dir must be supplied.", call. = FALSE)
}

# Validate the mode
if (!opt$mode %in% c("destabilize", "stabilize")) {
  stop("Invalid --mode. Must be 'destabilize' or 'stabilize'.", call. = FALSE)
}

# Validate the strategy based on the mode
if (opt$mode == "destabilize" && !opt$mutation_strategy %in% c("AU", "GU")) {
  stop("For 'destabilize' mode, --mutation_strategy must be 'AU' or 'GU'.", call. = FALSE)
}
if (opt$mode == "stabilize" && !opt$mutation_strategy %in% c("GC", "GU")) {
  stop("For 'stabilize' mode, --mutation_strategy must be 'GC' or 'GU'.", call. = FALSE)
}

# ---- Interactive Debugging Setup ----
# If running in RStudio, manually define the 'opt' list with your test arguments.
# cat("Running in interactive mode. Using manual 'opt' list for debugging.\n")
# 
# opt <- list(
#   fold_dir = "/scratch/users/rodell/RNAfold_psipos",
#   output_dir = "/scratch/users/rodell/deltaG/alter_deltaG/both/stabilizeGC_window7",
#   psipos = 59,
#   window = 7, # This was specified in your command
#   deltaG_threshold = 5.0,
#   mode = "stabilize",
#   mutation_strategy = "GC"
# )


# ---- Dynamic Configuration ----
# Build the CONFIG list dynamically using parsed arguments and the script's location.
script_dir <- this.path::this.dir()

CONFIG <- list(
  # Analysis Parameters from command line
  psipos = opt$psipos,
  window = opt$window,
  deltaG_threshold = opt$deltaG_threshold,
  mode = opt$mode,
  mutation_strategy = opt$mutation_strategy,
  
  # Input paths based on user input
  org_fold_dir    = opt$fold_dir,
  org_deltaG_path = file.path(opt$fold_dir, "RNAfold_deltaG.csv"), # Assumes deltaG file is in fold_dir
  
  # Paths to helper scripts, located relative to this script
  rna_fold_script = file.path(script_dir, "run_RNAfold.sh"),
  delta_g_script  = file.path(script_dir, "deltaG_RNAfold.R"),
  
  # Intermediate file paths, organized within the main output directory
  intermediate_dir   = file.path(opt$output_dir, "intermediate"),
  mutated_fasta_path = file.path(opt$output_dir, "intermediate", "mutated_sequences.fa"),
  mut_fold_dir       = file.path(opt$output_dir, "intermediate", "mutated_structures"),
  inc_deltaG_path    = file.path(opt$output_dir, "intermediate", "RNAfold_deltaG.csv"),
  
  # Final output file paths
  output_csv         = file.path(opt$output_dir, "comprehensive_mutation_analysis.csv"),
  output_fasta       = file.path(opt$output_dir, "increased_deltaG_all.fa"),
  output_summary     = file.path(opt$output_dir, "summary_statistics.txt")
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

# ---- Sequence Mutations ----

identify_gc_pairs <- function(sequence, dot_bracket, psipos, window) {
  seq_length <- nchar(sequence)
  start <- max(1, psipos - window)
  end <- min(seq_length, psipos + window)
  
  gc_pairs <- list()
  stack <- list()
  
  protected_range <- (psipos-2):(psipos+2)
  
  for (i in 1:seq_length) {
    if (substr(dot_bracket, i, i) == "(") {
      stack <- c(stack, i)
    } else if (substr(dot_bracket, i, i) == ")") {
      if (length(stack) > 0) {
        open <- stack[[length(stack)]]
        stack <- stack[-length(stack)]
        
        if ((open >= start && open <= end) || (i >= start && i <= end)) {
          base1 <- substr(sequence, open, open)
          base2 <- substr(sequence, i, i)
          if ((base1 == "G" && base2 == "C") || (base1 == "C" && base2 == "G")) {
            if (!(open %in% protected_range) && !(i %in% protected_range)) {
              gc_pairs[[length(gc_pairs) + 1]] <- c(open, i)
            }
          }
        }
      }
    }
  }
  
  return(gc_pairs)
}

identify_au_pairs <- function(sequence, dot_bracket, psipos, window) {
  seq_length <- nchar(sequence)
  start <- max(1, psipos - window)
  end <- min(seq_length, psipos + window)
  
  au_pairs <- list()
  stack <- list()
  protected_range <- (psipos - 2):(psipos + 2)
  
  for (i in 1:seq_length) {
    if (substr(dot_bracket, i, i) == "(") {
      stack <- c(stack, i)
    } else if (substr(dot_bracket, i, i) == ")") {
      if (length(stack) > 0) {
        open <- stack[[length(stack)]]
        stack <- stack[-length(stack)]
        
        if ((open >= start && open <= end) || (i >= start && i <= end)) {
          base1 <- substr(sequence, open, open)
          base2 <- substr(sequence, i, i)
          if ((base1 == "A" && base2 == "U") || (base1 == "U" && base2 == "A")) {
            if (!(open %in% protected_range) && !(i %in% protected_range)) {
              au_pairs[[length(au_pairs) + 1]] <- c(open, i)
            }
          }
        }
      }
    }
  }
  return(au_pairs)
}

mutate_sequence <- function(sequence, pairs_to_mutate, mode, strategy) {
  seq_vec <- strsplit(sequence, "")[[1]]
  mutations <- c()
  
  # Logic Gate for Mode
  if (mode == "destabilize") {
    # This is the existing logic for GC -> AU/GU
    if (strategy == "AU") {
      for (pair in pairs_to_mutate) {
        if (seq_vec[pair[1]] == "G") { seq_vec[pair[1]] <- "A" } else { seq_vec[pair[1]] <- "U" }
        mutations <- c(mutations, pair[1])
        if (seq_vec[pair[2]] == "G") { seq_vec[pair[2]] <- "A" } else { seq_vec[pair[2]] <- "U" }
        mutations <- c(mutations, pair[2])
      }
    } else if (strategy == "GU") {
      for (pair in pairs_to_mutate) {
        if (seq_vec[pair[1]] == "C") {
          seq_vec[pair[1]] <- "U"; mutations <- c(mutations, pair[1])
        } else if (seq_vec[pair[2]] == "C") {
          seq_vec[pair[2]] <- "U"; mutations <- c(mutations, pair[2])
        }
      }
    }
  } else if (mode == "stabilize") {
    # This is the new logic for AU -> GC/GU
    if (strategy == "GC") {
      for (pair in pairs_to_mutate) {
        if (seq_vec[pair[1]] == "A") { seq_vec[pair[1]] <- "G" } else { seq_vec[pair[1]] <- "C" }
        mutations <- c(mutations, pair[1])
        if (seq_vec[pair[2]] == "A") { seq_vec[pair[2]] <- "G" } else { seq_vec[pair[2]] <- "C" }
        mutations <- c(mutations, pair[2])
      }
    } else if (strategy == "GU") {
      for (pair in pairs_to_mutate) {
        if (seq_vec[pair[1]] == "A") {
          seq_vec[pair[1]] <- "G"; mutations <- c(mutations, pair[1])
        } else if (seq_vec[pair[2]] == "A") {
          seq_vec[pair[2]] <- "G"; mutations <- c(mutations, pair[2])
        }
      }
    }
  }
  
  return(list(sequence = paste(seq_vec, collapse = ""), mutations = mutations))
}

format_locations <- function(gc_pairs, mutations, sequence) {
  if (length(gc_pairs) == 0 || length(mutations) == 0) {
    return(list(basepair = "NA", mutation = "NA"))
  }
  
  basepair_str <- paste(sapply(gc_pairs, function(pair) paste(pair, collapse = ":")), collapse = ", ")
  
  mutation_str <- paste(sapply(mutations, function(pos) {
    original <- substr(sequence, pos, pos)
    mutated <- if (original == "G") "A" else if (original == "C") "U" else original
    sprintf("%s%d%s", original, pos, mutated)
  }), collapse = ", ")
  
  return(list(basepair = basepair_str, mutation = mutation_str))
}

determine_location <- function(mutations, psipos) {
  if (length(mutations) == 0) return("NA")
  
  upstream <- sum(mutations < psipos)
  downstream <- sum(mutations > psipos)
  
  if (upstream > 0 && downstream > 0) {
    return("both")
  } else if (upstream > 0) {
    return("upstream")
  } else if (downstream > 0) {
    return("downstream")
  } else {
    return("NA")
  }
}

# ---- Analysis Functions ----

read_and_merge_deltaG <- function(results_df, org_deltaG_path, inc_deltaG_path, deltaG_threshold) {
  # Check for existence of input files
  if (!file.exists(org_deltaG_path)) stop("Original deltaG file not found at: ", org_deltaG_path)
  if (!file.exists(inc_deltaG_path)) stop("Increased deltaG file not found at: ", inc_deltaG_path)
  
  org_deltaG_all <- read.csv(org_deltaG_path)
  inc_deltaG_all <- read.csv(inc_deltaG_path)
  
  results_df %>%
    left_join(org_deltaG_all %>% select(sequence_id, org_deltaG = mfe_energy), by = c("name" = "sequence_id")) %>%
    left_join(inc_deltaG_all %>% select(sequence_id, mut_deltaG = mfe_energy), by = c("name" = "sequence_id")) %>%
    mutate(
      deltaG_change = mut_deltaG - org_deltaG,
      deltaG_category = case_when(
        deltaG_change > deltaG_threshold ~ "increase",
        deltaG_change < -deltaG_threshold ~ "decrease",
        TRUE ~ "no change"
      )
    )
}

compare_dot_brackets <- function(org_dot_bracket, mut_dot_bracket) {
  if (is.na(org_dot_bracket) || is.na(mut_dot_bracket)) {
    return(list(num_differences = NA, fraction_different = NA, different_pairings = list()))
  }
  if (nchar(org_dot_bracket) != nchar(mut_dot_bracket)) {
    stop("Dot-bracket notations have different lengths")
  }
  
  org_chars <- strsplit(org_dot_bracket, "")[[1]]
  mut_chars <- strsplit(mut_dot_bracket, "")[[1]]
  
  different_positions <- which(org_chars != mut_chars)
  differences <- length(different_positions)
  
  different_pairings <- data.frame(
    position = different_positions,
    original = org_chars[different_positions],
    mutated = mut_chars[different_positions],
    stringsAsFactors = FALSE
  )
  
  return(list(
    num_differences = differences,
    fraction_different = differences / nchar(org_dot_bracket),
    different_pairings = different_pairings
  ))
}

compare_fold_structures <- function(org_fold_dir, mut_fold_dir) {
  org_fold_data <- read_all_fold_files(org_fold_dir)
  mut_fold_data <- read_all_fold_files(mut_fold_dir)
  
  merged_data <- org_fold_data %>%
    left_join(mut_fold_data, by = "name", suffix = c("_org", "_mut"))
  
  comparison_results <- merged_data %>%
    rowwise() %>%
    mutate(
      comparison = list(compare_dot_brackets(dot_bracket_org, dot_bracket_mut)),
      num_differences = comparison$num_differences,
      fraction_different = comparison$fraction_different,
      different_pairings = list(comparison$different_pairings)
    ) %>%
    ungroup() %>%
    select(-comparison)
  
  return(comparison_results)
}

analyze_pairings <- function(different_pairings, psipos, mutations) {
  
  if (is.null(different_pairings) || nrow(different_pairings) == 0) {
    return(list(
      in_5mer = FALSE,
      same_as_mutations = FALSE,
      num_in_5mer = 0,
      num_same_as_mutations = 0
    ))
  }
  
  diff_positions <- different_pairings$position
  five_mer <- (psipos - 2):(psipos + 2)
  
  in_5mer <- any(diff_positions %in% five_mer)
  num_in_5mer <- sum(diff_positions %in% five_mer)
  
  same_as_mutations <- any(diff_positions %in% mutations)
  num_same_as_mutations <- sum(diff_positions %in% mutations)
  
  return(list(
    in_5mer = in_5mer,
    same_as_mutations = same_as_mutations,
    num_in_5mer = num_in_5mer,
    num_same_as_mutations = num_same_as_mutations
  ))
}

generate_summary_stats <- function(comprehensive_results, config) {
  # --- Section 0: Run Parameters ---
  run_parameters_text <- c(
    "===== Comprehensive Mutation Analysis Summary =====",
    "",
    "Run Parameters:",
    paste("  - Goal (Mode):", toupper(config$mode)),
    paste("  - Mutation Strategy:", toupper(config$mutation_strategy)),
    paste("  - Window Size:", config$window),
    "---",
    ""
  )
  
  # --- Section 1: DeltaG Categories (with clear labels) ---
  deltaG_categories <- table(comprehensive_results$deltaG_category)
  increase_count <- ifelse("increase" %in% names(deltaG_categories), deltaG_categories["increase"], 0)
  decrease_count <- ifelse("decrease" %in% names(deltaG_categories), deltaG_categories["decrease"], 0)
  no_change_count <- ifelse("no change" %in% names(deltaG_categories), deltaG_categories["no change"], 0)
  
  deltaG_text <- c(
    "1. Overall Delta G Change Distribution:",
    paste("   - Destabilized (delta G became less negative/higher):", increase_count),
    paste("   - Stabilized (delta G became more negative/lower):", decrease_count),
    paste("   - No significant change:", no_change_count),
    ""
  )
  
  # --- Section 2: Mutation Location Summary ---
  mutated_df <- comprehensive_results %>% filter(num_mutations > 0)
  upstream_count <- sum(mutated_df$location == "upstream", na.rm = TRUE)
  downstream_count <- sum(mutated_df$location == "downstream", na.rm = TRUE)
  both_count <- sum(mutated_df$location == "both", na.rm = TRUE)
  
  location_text <- c(
    "2. Mutation Location Summary (for all mutated sequences):",
    paste("   - Upstream only:", upstream_count),
    paste("   - Downstream only:", downstream_count),
    paste("   - Both upstream and downstream:", both_count),
    ""
  )
  
  # --- Section 3: Detailed Analysis of Successful Sequences ---
  success_df <- data.frame()
  header_text <- ""
  
  if (config$mode == "destabilize") {
    success_df <- comprehensive_results %>% filter(deltaG_category == "increase")
    header_text <- "3. Detailed Analysis of Successfully DESTABILIZED Sequences (delta G increased):"
  } else { # mode == "stabilize"
    success_df <- comprehensive_results %>% filter(deltaG_category == "decrease")
    header_text <- "3. Detailed Analysis of Successfully STABILIZED Sequences (delta G decreased):"
  }
  
  if (nrow(success_df) > 0) {
    fraction_diff_less_0.05 <- sum(success_df$fraction_different < 0.05, na.rm = TRUE)
    fraction_diff_less_0.1 <- sum(success_df$fraction_different < 0.1, na.rm = TRUE)
    fraction_diff_less_0.15 <- sum(success_df$fraction_different < 0.15, na.rm = TRUE)
    no_diff_pair_in_5mer <- sum(!success_df$diff_pair_in_5mer, na.rm = TRUE)
    no_diff_pair_same_as_mutations <- sum(!success_df$diff_pair_same_as_mutations, na.rm = TRUE)
    all_conditions_met <- sum(
      success_df$fraction_different < 0.15 &
        !success_df$diff_pair_in_5mer &
        !success_df$diff_pair_same_as_mutations,
      na.rm = TRUE
    )
    
    success_analysis_text <- c(
      header_text,
      paste("   Total successful sequences:", nrow(success_df)),
      paste("   a. Structure conservation (fraction different):"),
      paste("      - Less than 5% different:", fraction_diff_less_0.05),
      paste("      - Less than 10% different:", fraction_diff_less_0.1),
      paste("      - Less than 15% different:", fraction_diff_less_0.15),
      paste("   b. No structural changes in 5-mer protected region:", no_diff_pair_in_5mer),
      paste("   c. No structural changes at the exact mutation sites:", no_diff_pair_same_as_mutations),
      "",
      paste("   => Sequences meeting key conditions (fraction < 0.15 AND b AND c):", all_conditions_met)
    )
  } else {
    success_analysis_text <- c(
      header_text,
      paste("   No sequences were successfully", toupper(config$mode), "d according to the threshold.")
    )
  }
  
  # --- Combine and Write ---
  summary_text <- c(
    run_parameters_text,
    deltaG_text,
    location_text,
    success_analysis_text,
    "",
    "================================================="
  )
  
  output_file <- config$output_summary
  writeLines(summary_text, output_file)
  cat("Summary statistics have been written to", output_file, "\n")
}

# ---- Main Workflow ----

process_fold_file <- function(file_path, config) {
  fold_data <- read_fold_file(file_path)
  if (is.null(fold_data)) return(NULL)
  
  # Logic Gate: Choose which pairs to identify based on the mode
  if (config$mode == "destabilize") {
    pairs_to_mutate <- identify_gc_pairs(fold_data$seq, fold_data$dot_bracket, config$psipos, config$window)
  } else {# mode == "stabilize"
    pairs_to_mutate <- identify_au_pairs(fold_data$seq, fold_data$dot_bracket, config$psipos, config$window)
  }
  
  # Pass all necessary parameters to the updated mutate_sequence function
  mutation_result <- mutate_sequence(fold_data$seq, pairs_to_mutate, config$mode, config$mutation_strategy)
  
  locations <- format_locations(pairs_to_mutate, mutation_result$mutations, fold_data$seq)
  location <- determine_location(mutation_result$mutations, config$psipos)
  
  data.frame(
    name = fold_data$name,
    seq = fold_data$seq,
    psipos = config$psipos,
    mutated_seq = mutation_result$sequence,
    mutations = I(list(mutation_result$mutations)),
    location_basepair = locations$basepair,
    location_mutation = locations$mutation,
    location = location,
    num_mutations = length(mutation_result$mutations),
    stringsAsFactors = FALSE
  )
}

main_workflow <- function(config) {
  # Create output directories if they don't already exist.
  dir.create(config$intermediate_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(config$mut_fold_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Process original fold files to generate mutations
  cat("Processing original fold files...\n")
  all_fold_files <- list.files(config$org_fold_dir, pattern = "\\.fold$", full.names = TRUE)
  results <- do.call(rbind, lapply(all_fold_files, function(file) process_fold_file(file, config)))
  
  # Save mutated sequences to an intermediate FASTA file
  create_and_save_fasta(results, "mutated_seq", "name", config$mutated_fasta_path)
  
  # Run RNAfold on the mutated sequences
  cat("Running RNAfold on mutated sequences...\n")
  system2("bash", 
          args = c(config$rna_fold_script, 
                   config$mutated_fasta_path, 
                   config$mut_fold_dir),
          stdout = TRUE, stderr = TRUE)
  
  # Run deltaG calculation on the new .fold files
  cat("Extracting deltaG for mutated sequences...\n")
  # Note the arguments for the deltaG script: -i is the input dir, -o is the output dir
  system2("Rscript",
          args = c(config$delta_g_script, 
                   "-i", config$mut_fold_dir,
                   "-o", config$intermediate_dir),
          stdout = TRUE, stderr = TRUE)
  
  cat("Comparing original and mutated structures...\n")
  # Merge deltaG data
  results_with_deltaG <- read_and_merge_deltaG(results, config$org_deltaG_path, config$inc_deltaG_path, config$deltaG_threshold)
  
  # Compare fold structures
  fold_comparison <- compare_fold_structures(config$org_fold_dir, config$mut_fold_dir)
  
  # Merge all results
  comprehensive_results <- results_with_deltaG %>%
    left_join(fold_comparison, by = "name") %>%
    rowwise() %>%
    mutate(
      pairing_analysis = list(analyze_pairings(different_pairings, psipos, mutations)),
      diff_pair_in_5mer = pairing_analysis$in_5mer,
      diff_pair_same_as_mutations = pairing_analysis$same_as_mutations,
      num_diff_pair_in_5mer = pairing_analysis$num_in_5mer,
      num_diff_pair_same_as_mutations = pairing_analysis$num_same_as_mutations
    ) %>%
    ungroup() %>%
    select(-pairing_analysis)
  
  # Convert list columns to string representation for CSV output
  comprehensive_results <- comprehensive_results %>%
    mutate(across(where(is.list), ~sapply(.x, function(x) paste(unlist(x), collapse = ","))))
  
  # Save final results
  cat("Saving final results...\n")
  write.csv(comprehensive_results, config$output_csv, row.names = FALSE)
  
  # Save a final FASTA file of sequences that showed increased deltaG
  # increased_deltaG_results <- filter(comprehensive_results, deltaG_category == "increase")
  # create_and_save_fasta(increased_deltaG_results, "mutated_seq", "name", config$output_fasta)
  if (config$mode == "destabilize") {
    successful_results <- filter(comprehensive_results, deltaG_category == "increase")
  } else { # mode == "stabilize"
    successful_results <- filter(comprehensive_results, deltaG_category == "decrease")
  }
  create_and_save_fasta(successful_results, "mutated_seq", "name", config$output_fasta)
  
  # Generate and save the final summary statistics
  generate_summary_stats(comprehensive_results, config)
  
  cat("Workflow complete. Outputs are in:", dirname(config$output_csv), "\n")
  #return(comprehensive_results)
}


# ---- Execution ----
# The main workflow is called with the CONFIG list built from optparse arguments.
main_workflow(CONFIG)
