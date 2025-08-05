library(Biostrings)

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

process_fold_file <- function(file_path, psipos = 59) {
  fold_data <- read_fold_file(file_path)
  if (is.null(fold_data)) return(NULL)
  
  gc_pairs <- identify_gc_pairs(fold_data$seq, fold_data$dot_bracket, psipos)
  mutation_result <- mutate_sequence(fold_data$seq, gc_pairs, psipos)
  locations <- format_locations(gc_pairs, mutation_result$mutations)
  location <- determine_location(mutation_result$mutations, psipos)
  
  return(data.frame(
    name = fold_data$name,
    seq = fold_data$seq,
    psipos = psipos,
    incG_seq = mutation_result$sequence,
    location_basepair = locations$basepair,
    location_mutation = locations$mutation,
    location = location,
    num_mutations = length(mutation_result$mutations),
    stringsAsFactors = FALSE
  ))
}

process_fold_file <- function(file_path, psipos = 59) {
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
# ---- Execution ----

# Main script
fold_dir <- "/scratch/users/rodell/RNAfold_psipos"
all_fold_files <- list.files(fold_dir, pattern = "\\.fold$", full.names = TRUE)


# Execute full script
results <- lapply(all_fold_files, process_fold_file)
results_df <- do.call(rbind, results)

# Create a named DNAStringSet
named_sequences <- DNAStringSet(results_df$incG_seq_dna)
names(named_sequences) <- results_df$name

# Write output files
write.csv(results_df, "increased_deltaG_all.csv", row.names = FALSE)
writeXStringSet(named_sequences, "increased_deltaG_all.fa", width = 10000)


# # testing script
# # Test function by function
# # Test read_fold_file
# cat("Testing read_fold_file:\n")
# sample_file <- all_fold_files[2]  # Take the first file as a sample
# print(read_fold_file(sample_file))
# fold_data <- read_fold_file(sample_file)
# 
# # Test identify_gc_pairs
# cat("Testing identify_gc_pairs:\n")
# gc_pairs <- identify_gc_pairs(fold_data$seq, fold_data$dot_bracket, psipos = 59, window = 20)
# 
# cat("Identified G:C pairs:\n")
# for (i in seq_along(gc_pairs)) {
#   pair <- gc_pairs[[i]]
#   base1 <- substr(fold_data$seq, pair[1], pair[1])
#   base2 <- substr(fold_data$seq, pair[2], pair[2])
#   cat(sprintf("Pair %d: %s-%s at positions %d-%d\n",
#               i, base1, base2, pair[1], pair[2]))
# }
# 
# print(gc_pairs)
# 
# # Test mutate_sequence
# cat("\nTesting mutate_sequence:\n")
# mutation_result <- mutate_sequence(fold_data$seq, gc_pairs, psipos = 59)
# print(mutation_result)
# 
# # Test format_locations
# cat("\nTesting format_locations:\n")
# locations <- format_locations(gc_pairs, mutation_result$mutations)
# print(locations)
# 
# # Test determine_location
# cat("\nTesting determine_location:\n")
# location <- determine_location(mutation_result$mutations, psipos = 59)
# print(location)
# 
# # Test process_fold_file
# cat("\nTesting process_fold_file:\n")
# result <- process_fold_file(sample_file, psipos = 59)
# print(result)
# 
# # Test the full process on a few files
# cat("\nTesting full process on multiple files:\n")
# test_files <- all_fold_files[1:5]  # Test with the first 5 files
# results <- lapply(test_files, process_fold_file)
# results_df <- do.call(rbind, results)
# print(results_df)
