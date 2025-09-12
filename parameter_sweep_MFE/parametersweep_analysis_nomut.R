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
base_dir <- "." # current working directory

# 1. Get only the immediate subdirectories (this is very fast)
target_dirs <- list.dirs(path = base_dir, recursive = FALSE)

# 2. Construct the full, expected path for the CSV in each directory
potential_file_paths <- file.path(target_dirs, "all_analysis_results.csv")

# 3. Keep only the paths for files that actually exist
file_paths_to_load <- potential_file_paths[file.exists(potential_file_paths)]

# 4. Read and combine
data <- map_dfr(file_paths_to_load, read.csv, .id = "source_file")

#data <- read.csv("all_analysis_results.csv")

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

write.csv(data, file = "analysis/all_results.csv", row.names = FALSE)

# ---- Plotting ----

# plot f1 scores versus complexity and both count

# Define the F1 score columns to plot on the y-axis
f1_columns <- c(
  "f1_both_structure",
  "f1_both_unuar_and_structure",
  "f1_both_uguag_and_structure"
)

# Define the variables for the x-axis
x_variables <- c(
  "complexity_score",
  "both_count"
)

# Create a list to store the plots
plot_list <- list()

# Use a nested loop to create and save a plot for each combination
for (y_col in f1_columns) {
  for (x_col in x_variables) {
    
    # Create nice, dynamic labels
    pretty_y_label <- gsub("_", " ", gsub("f1_both_", "", y_col))
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
    file_name <- paste0("analysis/", y_col, "_vs_", x_col, ".png")
    
    # 2. Use ggsave to save the plot object 'p'
    ggsave(filename = file_name, plot = p,
      width = 8, height = 8, units = "in")
  }
}

# plot complexity versus both count
p <- create_scatter_plot(
  data = data,
  x_var = "complexity_score",
  y_var = "both_count",
  
  # Provide descriptive titles and labels
  title = "Initial Match Count vs. Parameter Complexity",
  x_label = "Parameter Complexity Score (log2 Scale)",
  y_label = "Initial Sequence Matches (Both)",
  
  # Activate the log scale for the x-axis
  log2_x_scale = TRUE
)

ggsave(filename = "analysis/both_count_vs_complexity_score.png",
  plot = p, width = 8, height = 8, units = "in")

# ---- Specific Queries ----


# What is the initial both count for the highest f1 score?

# Define the F1 score columns you want to analyze
f1_columns_to_check <- c(
  "f1_both_structure",
  "f1_both_unuar_and_structure",
  "f1_both_uguag_and_structure"
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
  corresponding_count <- single_best_row$both_count
  tiebreaker_complexity <- single_best_row$complexity_score

  # f. Print the results, including the new count
  cat("The highest F1 score found is:", round(best_f1_score, 4), "\n")
  cat("Number of parameter combinations yielding this F1 score:", num_tied_rows, "\n")
  cat("The corresponding 'both_count' for the best row (after tie-break) is:", corresponding_count, "\n")
  cat("The complexity score of this best row is:", tiebreaker_complexity, "\n\n")

  # g. Print the full row for complete details
  cat("Full details of the best row (after tie-break):\n")
  print(single_best_row)
  cat("\n")
}

# ---- Find a Row ----

# bestf1 from 130 nt Extend & Trim parameter sweep:
bestf1_trimmed <- data %>%
  filter(
    offset_min == 0,
    offset_max == 1,
    unpaired1_min == 1,
    unpaired1_max == 4,
    paired1_min == 3,
    paired1_max == 5,
    unpaired2_min == 2,
    unpaired2_max == 5,
    paired2_min == 1,
    paired2_max == 4
  )

cat("Extend and trim bestf1: \n")
print(bestf1_trimmed)