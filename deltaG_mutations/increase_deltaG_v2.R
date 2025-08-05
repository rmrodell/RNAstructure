# ---- Libraries ----
library(Biostrings)
library(dplyr)
library(purrr)

# ---- Configuration ----
CONFIG <- list(
  psipos = 59,
  window = 10,
  deltaG_threshold = 5,
  fold_dir = "/scratch/users/rodell/RNAfold_psipos",
  org_deltaG_path = "/scratch/users/rodell/RNAfold_psipos/RNAfold_deltaG.csv",
  inc_deltaG_path = "/scratch/users/rodell/deltaG/increase_deltaG/v2/RNAfold_deltaG.csv",
  mut_fold_dir = "/scratch/users/rodell/deltaG/increase_deltaG/v2/RNAfold_incG",
  output_csv = "comprehensive_mutation_analysis.csv",
  output_fasta = "increased_deltaG_all.fa",  # Changed from "increased_deltaG_actually.fa"
  rna_fold_script = "/scratch/users/rodell/motifmatcher/RNAstructure/run_RNAfold.sh",
  rna_fold_output_folder = "RNAfold_incG",
  delta_g_script = "/scratch/users/rodell/motifmatcher/RNAstructure/deltaG_mutations/deltaG_RNAfold.R",
  delta_g_input = "/scratch/users/rodell/deltaG/increase_deltaG/v2/RNAfold_incG",
  delta_g_output = "/scratch/users/rodell/deltaG/increase_deltaG/v2"
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
  fold_data <- do.call(rbind, lapply(fold_files, read_fold_file))
  return(fold_data)
}

create_and_save_fasta <- function(df, seq_column, name_column, fasta_path) {
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

mutate_sequence <- function(sequence, gc_pairs) {
  seq_vec <- strsplit(sequence, "")[[1]]
  mutations <- c()
  
  for (pair in gc_pairs) {
    if (seq_vec[pair[1]] == "G") {
      seq_vec[pair[1]] <- "A"
      mutations <- c(mutations, pair[1])
    } else if (seq_vec[pair[1]] == "C") {
      seq_vec[pair[1]] <- "U"
      mutations <- c(mutations, pair[1])
    }
    
    if (seq_vec[pair[2]] == "G") {
      seq_vec[pair[2]] <- "A"
      mutations <- c(mutations, pair[2])
    } else if (seq_vec[pair[2]] == "C") {
      seq_vec[pair[2]] <- "U"
      mutations <- c(mutations, pair[2])
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
  
  if (length(different_pairings) == 0) {
    cat("different_pairings is empty\n")
    return(list(
      in_5mer = FALSE,
      same_as_mutations = FALSE,
      num_in_5mer = 0,
      num_same_as_mutations = 0
    ))
  }
  
  # Extract the positions from different_pairings
  diff_positions <- different_pairings[[1]]
  
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

generate_summary_stats <- function(comprehensive_results, output_file) {
  # Calculate summary statistics
  deltaG_categories <- table(comprehensive_results$deltaG_category)
  
  increased_deltaG_df <- comprehensive_results[comprehensive_results$deltaG_category == "increase", ]
  
  fraction_diff_less_0.05 <- sum(increased_deltaG_df$fraction_different < 0.05)
  fraction_diff_less_0.1 <- sum(increased_deltaG_df$fraction_different < 0.1)
  fraction_diff_less_0.15 <- sum(increased_deltaG_df$fraction_different < 0.15)
  
  no_diff_pair_in_5mer <- sum(!increased_deltaG_df$diff_pair_in_5mer)
  no_diff_pair_same_as_mutations <- sum(!increased_deltaG_df$diff_pair_same_as_mutations)
  
  all_conditions_met <- sum(
    increased_deltaG_df$fraction_different < 0.15 &
      !increased_deltaG_df$diff_pair_in_5mer &
      !increased_deltaG_df$diff_pair_same_as_mutations
  )
  
  # Create summary text
  summary_text <- c(
    "Summary Statistics:",
    "",
    "1. DeltaG categories:",
    paste("   Increased deltaG:", deltaG_categories["increase"]),
    paste("   Decreased deltaG:", deltaG_categories["decrease"]),
    paste("   No change in deltaG:", deltaG_categories["no change"]),
    "",
    "2. For sequences with increased deltaG:",
    paste("   Total sequences:", nrow(increased_deltaG_df)),
    paste("   a. Fraction different < 0.05:", fraction_diff_less_0.05),
    paste("      Fraction different < 0.10:", fraction_diff_less_0.1),
    paste("      Fraction different < 0.15:", fraction_diff_less_0.15),
    paste("   b. No differences in 5-mer region:", no_diff_pair_in_5mer),
    paste("   c. No differences at mutation sites:", no_diff_pair_same_as_mutations),
    paste("   Sequences meeting all conditions (fraction < 0.15 & b & c):", all_conditions_met)
  )
  
  # Write summary to file
  writeLines(summary_text, output_file)
  
  cat("Summary statistics have been written to", output_file, "\n")
}
# ---- Main Workflow ----
process_fold_file <- function(file_path, config) {
  fold_data <- read_fold_file(file_path)
  if (is.null(fold_data)) return(NULL)
  
  gc_pairs <- identify_gc_pairs(fold_data$seq, fold_data$dot_bracket, config$psipos, config$window)
  mutation_result <- mutate_sequence(fold_data$seq, gc_pairs)
  locations <- format_locations(gc_pairs, mutation_result$mutations, fold_data$seq)
  location <- determine_location(mutation_result$mutations, config$psipos)
  
  data.frame(
    name = fold_data$name,
    seq = fold_data$seq,
    psipos = config$psipos,
    mutated_seq = mutation_result$sequence,
    mutations = I(list(mutation_result$mutations)),  # Store mutations as a list column
    location_basepair = locations$basepair,
    location_mutation = locations$mutation,
    location = location,
    num_mutations = length(mutation_result$mutations),
    stringsAsFactors = FALSE
  )
}

main_workflow <- function(config) {
  # Process fold files
  cat("Process fold files..\n")
  all_fold_files <- list.files(config$fold_dir, pattern = "\\.fold$", full.names = TRUE)
  results <- do.call(rbind, lapply(all_fold_files, function(file) process_fold_file(file, config)))
  
  # Save initial FASTA file
  create_and_save_fasta(results, "mutated_seq", "name", config$output_fasta)
  
  # Run RNAfold
  cat("Run RNAfold..\n")
  rna_fold_result <- system2("bash", 
                             args = c(config$rna_fold_script, config$output_fasta, 
                                      config$rna_fold_output_folder),
                             stdout = TRUE,
                             stderr = TRUE)

  # Run deltaG_RNAfold.R
  cat("Extract deltaG..\n")
  delta_g_result <- system2("Rscript",
                            args = c(config$delta_g_script, 
                                     "-i", config$delta_g_input,
                                     "-o", config$delta_g_output),
                            stdout = TRUE,
                            stderr = TRUE)

  cat("Compare original and mutated structures..\n")
  # Merge deltaG data
  results_with_deltaG <- read_and_merge_deltaG(results, config$org_deltaG_path, config$inc_deltaG_path, config$deltaG_threshold)
  
  # Compare fold structures
  fold_comparison <- compare_fold_structures(config$fold_dir, config$mut_fold_dir)
  
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
    ungroup()%>%
    select(-pairing_analysis)
  
  # Convert list columns to string representation
  comprehensive_results <- comprehensive_results %>%
    mutate(across(where(is.list), ~sapply(.x, paste, collapse = ",")))
  
  # Save results
  write.csv(comprehensive_results, config$output_csv, row.names = FALSE)
  create_and_save_fasta(comprehensive_results %>% filter(deltaG_category == "increase"), "mutated_seq", "name", config$output_fasta)
  
  generate_summary_stats(comprehensive_results, "summary_statistics.txt")
  
  return(comprehensive_results)
}

# ---- Execution ----
results <- main_workflow(CONFIG)

