library(Biostrings)
library(dplyr)
library(purrr)

# ---- Helper functions ----

# Function to read and parse .fold files
read_fold_file <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  if (length(lines) >= 3) {
    name <- sub(">", "", lines[1])  # Remove ">" if present
    seq <- lines[2]
    dot_bracket <- strsplit(lines[3], " ")[[1]][1]
    return(data.frame(name = name, seq = seq, dot_bracket = dot_bracket, stringsAsFactors = FALSE))
  } else {
    warning(paste("File", file_path, "does not have the expected format"))
    return(NULL)
  }
}

# Function to read all .fold files in a directory
read_all_fold_files <- function(fold_dir) {
  fold_files <- list.files(fold_dir, pattern = "\\.fold$", full.names = TRUE)
  fold_data <- map_dfr(fold_files, read_fold_file)
  return(fold_data)
}

# Function to identify G:C base pairs in the region of interest
identify_gc_pairs <- function(sequence, dot_bracket, psipos, window = 20) {
  seq_length <- nchar(sequence)
  start <- max(1, psipos - window)
  end <- min(seq_length, psipos + window)
  
  gc_pairs <- list()
  stack <- list()
  
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
            gc_pairs[[length(gc_pairs) + 1]] <- c(open, i)
          }
        }
      }
    }
  }
  
  # Validation
  for (pair in gc_pairs) {
    base1 <- substr(sequence, pair[1], pair[1])
    base2 <- substr(sequence, pair[2], pair[2])
    if (!((base1 == "G" && base2 == "C") || (base1 == "C" && base2 == "G"))) {
      stop(sprintf("Invalid G:C pair found at positions %d and %d", pair[1], pair[2]))
    }
  }
  
  return(gc_pairs)
}

# Function to mutate C to T in identified base pairs
mutate_sequence <- function(sequence, gc_pairs, psipos) {
  seq_vec <- strsplit(sequence, "")[[1]]
  mutations <- c()
  
  for (pair in gc_pairs) {
    if (abs(pair[1] - psipos) >= 3 && abs(pair[2] - psipos) >= 3) {
      if (seq_vec[pair[1]] == "C") {
        seq_vec[pair[1]] <- "U"
        mutations <- c(mutations, pair[1])
      }
      if (seq_vec[pair[2]] == "C") {
        seq_vec[pair[2]] <- "U"
        mutations <- c(mutations, pair[2])
      }
    }
  }
  
  return(list(sequence = paste(seq_vec, collapse = ""), mutations = mutations))
}

# Function to format location strings
format_locations <- function(gc_pairs, mutations) {
  if (length(gc_pairs) == 0 || length(mutations) == 0) {
    return(list(basepair = "NA", mutation = "NA"))
  }
  
  basepair_str <- paste(sapply(gc_pairs, function(pair) paste(pair, collapse = ":")), collapse = ", ")
  mutation_str <- paste(mutations, collapse = ", ")
  
  return(list(basepair = basepair_str, mutation = mutation_str))
}

# Function to determine mutation location relative to psipos
determine_location <- function(mutations, psipos) {
  if (length(mutations) == 0) return("NA")
  
  if (all(mutations < psipos)) {
    return("upstream")
  } else if (all(mutations > psipos)) {
    return("downstream")
  } else {
    return("both")
  }
}

# Main running function for making mutations
#process_fold_file <- function(file_path, psipos = 59) {
  fold_data <- read_fold_file(file_path)
  if (is.null(fold_data)) return(NULL)
  
  gc_pairs <- identify_gc_pairs(fold_data$seq, fold_data$dot_bracket, psipos)
  mutation_result <- mutate_sequence(fold_data$seq, gc_pairs, psipos)
  locations <- format_locations(gc_pairs, mutation_result$mutations)
  location <- determine_location(mutation_result$mutations, psipos)
  
  incG_seq_dna <- gsub("U", "T", mutation_result$sequence)
  
  return(data.frame(
    name = fold_data$name,
    seq = fold_data$seq,
    psipos = psipos,
    incG_seq = mutation_result$sequence,
    incG_seq_dna = incG_seq_dna,
    location_basepair = locations$basepair,
    location_mutation = locations$mutation,
    location = location,
    num_mutations = length(mutation_result$mutations),
    stringsAsFactors = FALSE
  ))
}
# Main running function for making mutations
process_fold_file <- function(file_path, psipos = 59) {
  fold_data <- read_fold_file(file_path)
  if (is.null(fold_data)) return(NULL)
  
  gc_pairs <- identify_gc_pairs(fold_data$seq, fold_data$dot_bracket, psipos)
  mutation_result <- mutate_sequence(fold_data$seq, gc_pairs, psipos)
  locations <- format_locations(gc_pairs, mutation_result$mutations)
  location <- determine_location(mutation_result$mutations, psipos)
  
  result_df <- data.frame(
    name = fold_data$name,
    seq = fold_data$seq,
    psipos = psipos,
    incG_seq = mutation_result$sequence,
    location_basepair = locations$basepair,
    location_mutation = locations$mutation,
    location = location,
    num_mutations = length(mutation_result$mutations),
    stringsAsFactors = FALSE
  )
  
  # Use create_and_save_fasta to handle DNA conversion
  fasta_file <- paste0(tools::file_path_sans_ext(file_path), "_mutated.fa")
  create_and_save_fasta(result_df, "incG_seq", "name", "increased_deltaG_all.fa")
  
  return(result_df)
}

# Function to read and merge deltaG data
read_and_merge_deltaG <- function(results_df, org_deltaG_path, inc_deltaG_path, deltaG_threshold = 5) {
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

# Function to filter and save increased deltaG sequences
filter_and_save_increased_deltaG <- function(df, csv_path) {
  increased_deltaG <- df %>% filter(deltaG_category == "increase")
  write.csv(increased_deltaG, csv_path, row.names = FALSE)
  increased_deltaG
}

# Function to create and save FASTA file
create_and_save_fasta <- function(df, seq_column, name_column, fasta_path) {
  named_sequences <- setNames(df[[seq_column]], df[[name_column]])
  dna_sequences <- gsub("U", "T", named_sequences)
  mutated_sequences <- DNAStringSet(dna_sequences)
  
  stopifnot(all(names(mutated_sequences) == df[[name_column]]))
  
  writeXStringSet(mutated_sequences, fasta_path, width = 10000)
}

# Function to print summary
print_summary <- function(df) {
  df %>%
    summarise(
      Total_sequences = n(),
      Increased_deltaG = sum(deltaG_category == "increase"),
      Decreased_deltaG = sum(deltaG_category == "decrease"),
      No_change_deltaG = sum(deltaG_category == "no change")
    ) %>%
    print()
}

# Main function to run the analysis
analyze_mutations <- function(results_df, org_deltaG_path, inc_deltaG_path, org_fold_dir, mut_fold_dir, csv_output_path, fasta_output_path, deltaG_threshold = 5) {
  # Read and merge deltaG data
  results_df_with_deltaG <- read_and_merge_deltaG(results_df, org_deltaG_path, inc_deltaG_path, deltaG_threshold)
  
  # Compare fold structures
  fold_comparison <- compare_fold_structures(org_fold_dir, mut_fold_dir)
  
  # Merge deltaG and fold structure comparison results
  comprehensive_results <- results_df_with_deltaG %>%
    left_join(fold_comparison, by = "name") %>%
    mutate(
      deltaG_change = mut_deltaG - org_deltaG,
      deltaG_category = case_when(
        deltaG_change > deltaG_threshold ~ "increase",
        deltaG_change < -deltaG_threshold ~ "decrease",
        TRUE ~ "no change"
      )
    ) %>%
    select(
      name, 
      num_mutations, 
      org_deltaG, 
      mut_deltaG, 
      deltaG_change, 
      deltaG_category,
      num_differences,
      fraction_different,
      seq_org,
      seq_mut,
      dot_bracket_org,
      dot_bracket_mut,
      different_pairings
    )
  
  # Filter for sequences where deltaG actually increased
  increased_deltaG_actually <- comprehensive_results %>%
    filter(deltaG_category == "increase")
  
  # Save results
  write.csv(comprehensive_results, csv_output_path, row.names = FALSE)
  create_and_save_fasta(increased_deltaG_actually, "seq_mut", "name", fasta_output_path)
  
  # Print summary
  print_summary(comprehensive_results)
  
  # Return the comprehensive results
  return(comprehensive_results)
}

# Function to compare dot-bracket notations
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

# Workflow to compare dot-brackets
compare_fold_structures <- function(org_fold_dir, mut_fold_dir) {
  # Read original and mutated fold data
  org_fold_data <- read_all_fold_files(org_fold_dir)
  mut_fold_data <- read_all_fold_files(mut_fold_dir)
  
  # Merge the datasets
  merged_data <- org_fold_data %>%
    left_join(mut_fold_data, by = "name", suffix = c("_org", "_mut"))
  
  # Compare dot-brackets
  comparison_results <- merged_data %>%
    mutate(
      comparison = map2(dot_bracket_org, dot_bracket_mut, compare_dot_brackets),
      num_differences = map_dbl(comparison, "num_differences"),
      fraction_different = map_dbl(comparison, "fraction_different"),
      different_pairings = map(comparison, "different_pairings")
    )
  
  return(comparison_results)
}

# ---- Execution ----

# Main script
fold_dir <- "/scratch/users/rodell/RNAfold_psipos"
all_fold_files <- list.files(fold_dir, pattern = "\\.fold$", full.names = TRUE)


# Execute full script
results <- lapply(all_fold_files, process_fold_file)
results_df <- do.call(rbind, results)


# Write output files
write.csv(results_df, "increased_deltaG_all.csv", row.names = FALSE)


# TO DO: implement as system2() calls
# bash /scratch/users/rodell/motifmatcher/RNAstructure/run_RNAfold.sh increased_deltaG_all.fa RNAfold_incG
# Rscript /scratch/users/rodell/motifmatcher/RNAstructure/deltaG_mutations/deltaG_RNAfold.R -i /scratch/users/rodell/deltaG/increase_deltaG/RNAfold_incG -o /scratch/users/rodell/deltaG/increase_deltaG


# Analyze mutations
# analyze_mutations(
#   results_df,
#   "/scratch/users/rodell/RNAfold_psipos/RNAfold_deltaG.csv",
#   "/scratch/users/rodell/deltaG/increase_deltaG/RNAfold_deltaG.csv",
#   "increased_deltaG_actually.csv",
#   "increased_deltaG_actually.fa",
#   deltaG_threshold = 5
# )
# 
# # Compare dot brakcets
# # Usage
# org_fold_dir <- "/scratch/users/rodell/RNAfold_psipos"
# mut_fold_dir <- "/scratch/users/rodell/deltaG/increase_deltaG/RNAfold_incG"
# 
# comparison_results <- compare_fold_structures(org_fold_dir, mut_fold_dir)

# Usage
results <- analyze_mutations(
  results_df = results_df,  # Your original results dataframe
  org_deltaG_path = "/scratch/users/rodell/RNAfold_psipos/RNAfold_deltaG.csv",
  inc_deltaG_path = "/scratch/users/rodell/deltaG/increase_deltaG/RNAfold_deltaG.csv",
  org_fold_dir = "/scratch/users/rodell/RNAfold_psipos",
  mut_fold_dir = "/scratch/users/rodell/deltaG/increase_deltaG/RNAfold_incG",
  csv_output_path = "comprehensive_mutation_analysis.csv",
  fasta_output_path = "increased_deltaG_actually.fa",
  deltaG_threshold = 5
)
