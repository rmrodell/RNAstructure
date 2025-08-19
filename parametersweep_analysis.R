# ---- Libraries ----
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)

# ---- Values ----
# total number of modified sites
incell_total <- 178
invitro_total <- 217

# fl scores
f1_UNUAR_both <- 0.5544
f1_UGUAG_both <- 0.6689 
f1_UNUAR_incell <- 0.555555555555556
f1_UGUAG_invitro <- 0.754990925589837
f1_UGUAG_incell <- 0.652818991097923
f1_UNUAR_invitro <- 0.671875
f1_ML_incell <- 0.666
f1_ML_invitro <- 0.54

# ---- Functions ----
# determine length of a range
calc_range_robust <- function(min_val, max_val) {
  range_size <- max_val - min_val + 1
  
  # If the calculated range is NA, replace it with 1 (the neutral multiplier)
  range_size[is.na(range_size)] <- 1
  
  return(range_size)
}

create_scatter_plot <- function(data, x_var, y_var, title, x_label, y_label,
                                unuar_line = NULL,
                                uguag_line = NULL,
                                y_lims = NULL,
                                log2_x_scale = FALSE) {
  
  # 1. Start with the base plot
  plot <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(alpha = 0.6, shape = 16) + # Use a solid shape for better visibility
    labs(title = title, x = x_label, y = y_label) +
    theme_minimal(base_size = 18, base_family = "sans") +
    theme(
      plot.title = element_text(size = 21.6, hjust = 0.5, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  # 2. Conditionally add UNUAR line and annotation if a value is provided
  if (!is.null(unuar_line)) {
    plot <- plot +
      geom_hline(yintercept = unuar_line, color = "red", linetype = "dashed") +
      annotate("text", x = Inf, y = unuar_line, label = "UNUAR", 
               color = "red", hjust = 1.1, vjust = -0.5, size = 4)
  }
  
  # 3. Conditionally add UGUAG line and annotation
  if (!is.null(uguag_line)) {
    plot <- plot +
      geom_hline(yintercept = uguag_line, color = "blue", linetype = "dashed") +
      annotate("text", x = Inf, y = uguag_line, label = "UGUAG", 
               color = "blue", hjust = 1.1, vjust = -0.5, size = 4)
  }
  
  # 4. Conditionally set y-axis limits (using the safer coord_cartesian)
  if (!is.null(y_lims)) {
    plot <- plot + coord_cartesian(ylim = y_lims)
  }
  
  # 5. Conditionally transform the x-axis
  if (log2_x_scale) {
    plot <- plot + scale_x_continuous(trans = 'log2')
  }
  
  return(plot)
}


# ---- Datasets ----
base_dir <- "/scratch/users/rodell/motifmatcher/20250816"

# 1. Get only the immediate subdirectories (this is very fast)
target_dirs <- list.dirs(path = base_dir, recursive = FALSE)

# 2. Construct the full, expected path for the CSV in each directory
potential_file_paths <- file.path(target_dirs, "all_analysis_results.csv")

# 3. Keep only the paths for files that actually exist
file_paths_to_load <- potential_file_paths[file.exists(potential_file_paths)]

# 4. Read and combine
data <- map_dfr(file_paths_to_load, read.csv, .id = "source_file")

# finish processing the data

# create complexity score
data <- data %>%
  mutate(
    complexity_score =
      # For optional 'unpaired1': calculates range or returns 1 if NA
      calc_range_robust(unpaired1_min, unpaired1_max) *
      
      # For required 'paired1': always calculates range
      calc_range_robust(paired1_min, paired1_max) *
      
      # For required 'unpaired2': always calculates range
      calc_range_robust(unpaired2_min, unpaired2_max) *
      
      # For optional 'paired2': calculates range or returns 1 if NA
      calc_range_robust(paired2_min, paired2_max)
  )

write.csv(data, file = "all_analysis_results.csv", row.names = FALSE)

# ---- Plotting ----

# sample plot
# create_scatter_plot(
#   data = data,
#   x_var = "complexity_score",
#   y_var = "initial_f1_both_structure", # Choosing this F1 score for the example
#   
#   # Titles and labels for clarity
#   title = "F1 Score (Structure) vs. Parameter Complexity",
#   x_label = "Parameter Complexity Score (log2 Scale)",
#   y_label = "F1 Score",
#   
#   # Provide the F1 baseline values to draw the reference lines
#   unuar_line = f1_UNUAR_both,
#   uguag_line = f1_UGUAG_both,
#   
#   # Set the y-axis limits to the standard F1 score range [0, 1]
#   y_lims = c(0, 1),
#   
#   # Enable log scaling for the x-axis, as complexity can vary widely
#   log2_x_scale = TRUE
# )

# plot f1 scores versus complexity and both count

# Define the F1 score columns to plot on the y-axis
f1_columns <- c(
  "initial_f1_both_structure",
  "initial_f1_both_unuar_and_structure",
  "initial_f1_both_uguag_and_structure"
)

# Define the variables for the x-axis
x_variables <- c(
  "complexity_score",
  "initial_both_count"
)

# Create a list to store the plots
plot_list <- list()

# Use a nested loop to create and save a plot for each combination
for (y_col in f1_columns) {
  for (x_col in x_variables) {
    
    # Create nice, dynamic labels
    pretty_y_label <- gsub("_", " ", gsub("initial_f1_both_", "", y_col))
    pretty_x_label <- gsub("_", " ", x_col)
    
    plot_title <- paste(tools::toTitleCase(pretty_y_label), "vs.", tools::toTitleCase(pretty_x_label))
    
    # Check if the x-axis should be log-scaled
    use_log_scale <- (x_col == "complexity_score")
    
    # Generate the plot using our unified function
    p <- create_scatter_plot(
      data = data,
      x_var = x_col,
      y_var = y_col,
      title = plot_title,
      x_label = tools::toTitleCase(pretty_x_label),
      y_label = "F1 Score",
      unuar_line = f1_UNUAR_both,
      uguag_line = f1_UGUAG_both,
      y_lims = c(0, 1),
      log2_x_scale = use_log_scale
    )
    
    # Store the plot in the list, using a descriptive name
    plot_name <- paste(y_col, "vs", x_col, sep = "_")
    plot_list[[plot_name]] <- p
    
    # 1. Create a clean filename
    file_name <- paste0(y_col, "_vs_", x_col, ".png")
    
    # 2. Use ggsave to save the plot object 'p'
    ggsave(filename = file_name, plot = p,
      width = 8, height = 8, units = "in")
  }
}

# plot complexity versus both count
p <- create_scatter_plot(
  data = data,
  x_var = "complexity_score",
  y_var = "initial_both_count",
  
  # Provide descriptive titles and labels
  title = "Initial Match Count vs. Parameter Complexity",
  x_label = "Parameter Complexity Score (log2 Scale)",
  y_label = "Initial Sequence Matches (Both)",
  
  # Activate the log scale for the x-axis
  log2_x_scale = TRUE
)

ggsave(filename = "initial_both_count_vs_complexity_score.png",
  plot = p, width = 8, height = 8, units = "in")

# ---- Specific Queries ----

# # What is the lowest complexity score that returns the highest count score?
# 
# # Step 1: Find the maximum value of 'initial_both_count'
# max_count <- max(data$initial_both_count, na.rm = TRUE)
# cat("The highest count score found in the data is:", max_count, "\n")
# 
# # Step 2: Filter the data to find all rows achieving the max count,
# #         then find the minimum complexity within that group and filter again.
# 
# optimal_rows <- data %>%
#   # First, get all rows that have the highest possible count
#   filter(initial_both_count == max_count)
# 
# # Display the resulting row(s)
# print(optimal_rows)


# What is the initial both count for the highest f1 score?

# Define the F1 score columns you want to analyze
f1_columns_to_check <- c(
  "initial_f1_both_structure",
  "initial_f1_both_unuar_and_structure",
  "initial_f1_both_uguag_and_structure"
)

# Loop through each F1 column to find the single best row
for (f1_col in f1_columns_to_check) {

  cat("=================================================================\n")
  cat("ANALYSIS FOR:", f1_col, "\n")
  cat("(Tie-breaker: lowest complexity score)\n")
  cat("=================================================================\n")

  # a. Find the maximum F1 score in the column
  max_f1_value <- data %>%
    pull(.data[[f1_col]]) %>%
    max(na.rm = TRUE)

  # b. Filter to get ALL rows that have this maximum F1 score
  top_f1_rows <- data %>%
    filter(.data[[f1_col]] == max_f1_value)

  # c. Count how many rows are tied for the top F1 score
  num_tied_rows <- nrow(top_f1_rows)

  # d. From those tied rows, select the single best one using the tie-breaker
  single_best_row <- top_f1_rows %>%
    arrange(complexity_score) %>%
    slice(1)

  # e. Extract key values from that single best row
  best_f1_score <- single_best_row[[f1_col]]
  corresponding_count <- single_best_row$initial_both_count
  tiebreaker_complexity <- single_best_row$complexity_score

  # f. Print the results, including the new count
  cat("The highest F1 score found is:", round(best_f1_score, 4), "\n")
  cat("Number of parameter combinations yielding this F1 score:", num_tied_rows, "\n")
  cat("The corresponding 'initial_both_count' for the best row (after tie-break) is:", corresponding_count, "\n")
  cat("The complexity score of this best row is:", tiebreaker_complexity, "\n\n")

  # g. Print the full row for complete details
  cat("Full details of the best row (after tie-break):\n")
  print(single_best_row)
  cat("\n")
}

# check if the best row has the maximal number of mutations

# Define the columns you want to inspect when there's a tie for the F1 score.
columns_to_check <- c(
  "rnafold_paired1_disruption_disrupted_count",
  "rnafold_paired2_disruption_disrupted_count",
  "paired1_disr_comp",
  "paired1_disr_comp_mod",
  "paired2_disr_comp",
  "paired2_disr_comp_mod",
  "combined_disr_comp",
  "combined_disr_comp_mod"
)

# Loop through each F1 score column you want to analyze
for (f1_col in f1_columns_to_check) {
  
  cat("=================================================================\n")
  cat("ANALYSIS FOR:", f1_col, "\n")
  cat("(Tie-breaker: lowest complexity score)\n")
  cat("=================================================================\n")
  
  # a. Find the maximum F1 score in the column
  max_f1_value <- data %>%
    pull(.data[[f1_col]]) %>%
    max(na.rm = TRUE)
  
  # b. Filter to get ALL rows that have this maximum F1 score
  top_f1_rows <- data %>%
    filter(.data[[f1_col]] == max_f1_value)
  
  # c. Count how many rows are tied for the top F1 score
  num_tied_rows <- nrow(top_f1_rows)
  
  cat("The highest F1 score found is:", round(max_f1_value, 4), "\n")
  cat("Number of parameter combinations yielding this F1 score:", num_tied_rows, "\n\n")
  
  # Only perform this detailed analysis if there is more than one top row.
  if (num_tied_rows > 1) {
    
    # 1. Calculate the range (min and max) for the specified columns across all tied rows.
    cat("--- Range of Values Across All", num_tied_rows, "Tied Rows ---\n")
    
    value_ranges <- top_f1_rows %>%
      summarise(across(all_of(columns_to_check), list(min = min, max = max), .names = "{.col}_{.fn}"))
    
    # Print the ranges in a readable format
    for (col_name in columns_to_check) {
      min_val <- value_ranges[[paste0(col_name, "_min")]]
      max_val <- value_ranges[[paste0(col_name, "_max")]]
      cat(sprintf("  %-40s Min: %-10s Max: %s\n", col_name, round(min_val, 4), round(max_val, 4)))
    }
    cat("\n")
    
    # 2. Check if the max values are present in the least complex subset.
    cat("--- Complexity Check ---\n")
    
    # Identify the least complex subset (rows with the minimum complexity score among the ties)
    min_complexity_in_ties <- min(top_f1_rows$complexity_score, na.rm = TRUE)
    least_complex_subset <- top_f1_rows %>%
      filter(complexity_score == min_complexity_in_ties)
    
    cat("The lowest complexity score among these tied rows is:", min_complexity_in_ties, "\n")
    
    all_max_in_least_complex <- TRUE
    for (col_name in columns_to_check) {
      overall_max_value <- value_ranges[[paste0(col_name, "_max")]]
      
      # Check if this overall max value exists in the least complex subset for this column
      is_max_in_subset <- overall_max_value %in% least_complex_subset[[col_name]]
      
      if (!is_max_in_subset) {
        all_max_in_least_complex <- FALSE
        cat(sprintf("  WARNING for '%s': Max value (%s) is NOT found in the least complex subset.\n", 
                    col_name, round(overall_max_value, 4)))
      }
    }
    
    if (all_max_in_least_complex) {
      cat("  SUCCESS: For all checked columns, the max value is available in the least complex subset.\n")
    }
    cat("\n")
  }

  # From all tied rows, select the single best one using the complexity tie-breaker
  single_best_row <- top_f1_rows %>%
    arrange(complexity_score) %>%
    slice(1)
  
  # Extract key values from that single best row
  corresponding_count <- single_best_row$initial_both_count
  tiebreaker_complexity <- single_best_row$complexity_score
  
  # Print the final result for the chosen row
  cat("--- Best Row Selected (After Tie-break) ---\n")
  cat("The corresponding 'initial_both_count' is:", corresponding_count, "\n")
  cat("The complexity score of this row is:", tiebreaker_complexity, "\n\n")
  
  cat("Full details of the best row:\n")
  print(as.data.frame(single_best_row)) # as.data.frame for cleaner printing
  cat("\n")
}


# find best f1 rows with maximal number of mutations

# Define the columns to use
key_metric_cols <- c(
  "combined_disr_comp_mod",
  "paired2_disr_comp_mod",
  "paired1_disr_comp_mod"
)

# Loop through each primary F1 score column to analyze
for (f1_col in f1_columns_to_check) {
  
  cat("=================================================================\n")
  cat("ANALYSIS FOR:", f1_col, "\n")
  cat("=================================================================\n")
  
  # Find all rows that share the highest F1 score for the current column
  max_f1_value <- max(data[[f1_col]], na.rm = TRUE)
  top_f1_rows <- data %>% filter(.data[[f1_col]] == max_f1_value)
  
  cat("Highest F1 score is:", round(max_f1_value, 4), 
      "found in", nrow(top_f1_rows), "rows.\n\n")
  
  # --- Find the specific rows of interest from the top-performing set ---
  
  # 1. The row with the LEAST complexity
  row_min_complexity <- top_f1_rows %>%
    arrange(complexity_score) %>%
    slice(1) %>%
    mutate(reason_for_selection = "Least Complexity")
  
  # 2. Find the rows that maximize each of the other key metrics
  rows_max_metrics <- purrr::map_df(key_metric_cols, function(metric_col) {
    top_f1_rows %>%
      arrange(desc(.data[[metric_col]])) %>%
      slice(1) %>%
      mutate(reason_for_selection = paste("Max", metric_col))
  })
  
  # --- Combine, de-duplicate, and display the results ---
  
  # Bind all selected rows together
  all_selected_rows <- bind_rows(
    row_min_complexity,
    rows_max_metrics
  )
  
  # Group by all columns except our added "reason" to find unique rows.
  # Then, paste all the reasons a unique row was selected into a single string.
  final_rows_to_print <- all_selected_rows %>%
    group_by(across(-reason_for_selection)) %>%
    summarise(
      reason_for_selection = paste(reason_for_selection, collapse = "; "), 
      .groups = "drop"
    ) %>%
    select(reason_for_selection, everything()) # Move the new reason column to the front
  
  cat("--- Key Rows from Top-Performing Set ---\n")
  print(as.data.frame(final_rows_to_print))
  cat("\n")
}


