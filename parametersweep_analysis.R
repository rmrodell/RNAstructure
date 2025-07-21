library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

create_f1_plot <- function(data, x_var, y_var, title, x_label, y_label, 
                           unuar_line, uguag_line, 
                           unuar_f1, uguag_f1,
                           log2_x_scale = FALSE) {
  plot <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point() +
    geom_hline(yintercept = unuar_line, color = "red", linetype = "dashed") +
    geom_hline(yintercept = uguag_line, color = "blue", linetype = "dashed") +
    annotate("text", x = max(data[[x_var]]), y = unuar_line, label = "UNUAR", 
             color = "red", hjust = 1, vjust = -0.5, size = 4) +
    annotate("text", x = max(data[[x_var]]), y = uguag_line, label = "UGUAG", 
             color = "blue", hjust = 1, vjust = -0.5, size = 4) +
    labs(title = title,
         x = x_label,
         y = y_label) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal(base_size = 18, base_family = "sans") +
    theme(plot.title = element_text(size = 21.6, hjust = 0.5),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  
  if (log2_x_scale) {
    plot <- plot + scale_x_continuous(trans = 'log2')
  }
  
  return(plot)
}

data <- read.csv("/scratch/users/rodell/motifmatcher/parametersweep_paralleled_output/all_analysis_results.csv")

# total number of modified sites
incell_total <- 178
invitro_total <- 217

# fl scores
f1_UNUAR_incell <- 0.555555555555556
f1_UGUAG_invitro <- 0.754990925589837
f1_UGUAG_incell <- 0.652818991097923
f1_UNUAR_invitro <- 0.671875
f1_ML_incell <- 0.666
f1_ML_invitro <- 0.54

# plot f1 versus matches
p1 <- create_f1_plot(data, "incell_count", "f1_incell_uguag_and_structure", 
                     "F1 In Cell Structure + UGUAG vs In Cell Count",
                     "In Cell Count", "F1 In Cell Structure + UGUAG",
                     f1_UNUAR_incell, f1_UGUAG_incell)

p2 <- create_f1_plot(data, "incell_count", "f1_incell_structure", 
                     "F1 In Cell Structure vs In Cell Count",
                     "In Cell Count", "F1 In Cell Structure",
                     f1_UNUAR_incell, f1_UGUAG_incell)

p3 <- create_f1_plot(data, "invitro_count", "f1_invitro_uguag_and_structure", 
                     "F1 In Vitro Structure + UGUAG vs In Vitro Count",
                     "In Vitro Count", "F1 In Vitro Structure + UGUAG",
                     f1_UNUAR_invitro, f1_UGUAG_invitro)

p4 <- create_f1_plot(data, "invitro_count", "f1_invitro_structure", 
                     "F1 In Vitro Structure vs In Vitro Count",
                     "In Vitro Count", "F1 In Vitro Structure",
                     f1_UNUAR_invitro, f1_UGUAG_invitro)

# save the plots to pngs
plot_list <- list(p1, p2, p3, p4)
plot_names <- c("output_plots/plot1.png", "output_plots/plot2.png", "output_plots/plot3.png", "output_plots/plot4.png")

for (i in 1:4) {
  ggsave(plot_names[i], plot = plot_list[[i]], width = 7, height = 7, units = "in", dpi = 300)
}


# create complexity score

# Function to calculate range
calc_range <- function(min_val, max_val) {
  return(max_val - min_val + 1)  # Adding 1 to include both min and max in the range
}

# Calculate complexity score
data$complexity_score <- 
  calc_range(data$unpaired1_min, data$unpaired1_max) *
  calc_range(data$paired1_min, data$paired1_max) *
  calc_range(data$unpaired2_min, data$unpaired2_max)

# If include_paired2 is TRUE, multiply by paired2 range as well
data$complexity_score <- ifelse(data$include_paired2, 
                                data$complexity_score * calc_range(data$paired2_min, data$paired2_max),
                                data$complexity_score)


# plot complexity versus f1
create_f1_plot(data, "complexity_score", "f1_incell_structure", 
               "F1 In Cell Structure vs Complexity",
               "log2(Complexity Score)", "F1 In Cell Structure",
               f1_UNUAR_incell, f1_UGUAG_incell,
               log2_x_scale = TRUE)

create_f1_plot(data, "complexity_score", "f1_incell_uguag_and_structure", 
               "F1 In Cell Structure + UGUAG vs Complexity",
               "log2(Complexity Score)", "F1 In Cell Structure + UGUAG",
               f1_UNUAR_incell, f1_UGUAG_incell,
               log2_x_scale = TRUE)

create_f1_plot(data, "complexity_score", "f1_invitro_uguag_and_structure", 
               "F1 In Vitro Structure + UGUAG vs Complexity",
               "log2(Complexity Score)", "F1 In Vitro Structure + UGUAG",
               f1_UNUAR_invitro, f1_UGUAG_invitro,
               log2_x_scale = TRUE)

create_f1_plot(data, "complexity_score", "f1_invitro_structure", 
               "F1 In Vitro Structure vs Complexity",
               "log2(Complexity Score)", "F1 In Vitro Structure",
               f1_UNUAR_invitro, f1_UGUAG_invitro,
               log2_x_scale = TRUE)

# filter for count > 30
# sort to find lowest complexity with highest f1
# adjust the complexity-f1 function to handle multiple rows with the same complexity and f1

find_best_scores <- function(data, f1_column, count_column, count_threshold = 30) {
  # Filter the data based on the count
  filtered_data <- data[data[[count_column]] > count_threshold, ]
  
  # If no data remains after filtering, return an empty dataframe
  if (nrow(filtered_data) == 0) {
    return(data.frame())
  }
  
  # Sort the filtered data
  sorted_data <- filtered_data[order(-filtered_data[[f1_column]], filtered_data$complexity_score), ]
  
  # Get the best F1 score and complexity
  best_f1 <- sorted_data[[f1_column]][1]
  best_complexity <- sorted_data$complexity_score[1]
  
  # Find all rows that match the best F1 score and complexity
  best_rows <- sorted_data[sorted_data[[f1_column]] == best_f1 & 
                             sorted_data$complexity_score == best_complexity, ]
  
  # Return the dataframe of best rows
  return(best_rows)
}

results_incell <- find_best_scores(data, "f1_incell_structure", "incell_count", 30)

# List of F1 columns
f1_columns <- c("f1_incell_structure", "f1_invitro_structure", 
                "f1_incell_unuar_and_structure", "f1_invitro_unuar_and_structure", 
                "f1_incell_uguag_and_structure", "f1_invitro_uguag_and_structure")

# Initialize an empty list to store results
results_list <- list()

# Loop through each F1 column
for (f1_col in f1_columns) {
  # Determine whether it's incell or invitro
  count_col <- ifelse(grepl("incell", f1_col), "incell_count", "invitro_count")
  
  # Find best scores
  best_scores <- find_best_scores(data, f1_col, count_col, count_threshold = 30)
  
  # Add a column to indicate which F1 metric this is
  best_scores$f1_metric <- f1_col
  
  # Add to results list
  results_list[[f1_col]] <- best_scores
}

# Combine all results into a single dataframe
all_results <- bind_rows(results_list)

write.csv(all_results, "bestf1_30seq_all.csv", row.names = FALSE)

# find row with best average f1 across all f1 in the row
# calculate the average F1 score
all_results <- all_results %>%
  rowwise() %>%
  mutate(avg_f1_score = mean(c_across(all_of(f1_columns)))) %>%
  ungroup()


# plot the f1 scores
# Prepare the data
# filter the rows based on the f1_metric, selecting the best for invitro and incell
filtered_results <- all_results %>%
  filter(f1_metric %in% c("f1_incell_structure", "f1_invitro_unuar_and_structure"))

# create the plot data
plot_data <- filtered_results %>%
  select(f1_metric, f1_incell_structure, f1_invitro_structure, 
         f1_incell_unuar_and_structure, f1_invitro_unuar_and_structure, 
         f1_incell_uguag_and_structure, f1_invitro_uguag_and_structure) %>%
  pivot_longer(cols = -f1_metric, names_to = "score_type", values_to = "f1_score") %>%
  mutate(
    type = case_when(
      grepl("incell", score_type) ~ "incell",
      grepl("invitro", score_type) ~ "invitro",
      TRUE ~ NA_character_
    ),
    metric = case_when(
      grepl("unuar_and_structure", score_type) ~ "UNUAR + Structure",
      grepl("uguag_and_structure", score_type) ~ "UGUAG + Structure",
      grepl("structure$", score_type) ~ "Structure",
      TRUE ~ NA_character_
    )
  )

new_rows <- tibble(
  f1_metric = c("all","all","all","all", "ML", "ML"),
  metric = c("UNUAR", "UNUAR", "UGUAG", "UGUAG", "ML", "ML"),
  f1_score = c(f1_UNUAR_incell, f1_UNUAR_invitro, f1_UGUAG_incell, f1_UGUAG_invitro, f1_ML_incell, f1_ML_invitro),
  type = c("incell", "invitro", "incell", "invitro", "incell", "invitro")
)

plot_data <- bind_rows(plot_data, new_rows)

# Function to create the bar plot
create_bar_plot <- function(data, plot_type, title) {
  
  # Define the metric order within the function
  metric_order <- c("UNUAR", "UGUAG", "Structure", "UNUAR + Structure", "UGUAG + Structure", "ML")
  
  ggplot(data, aes(x = metric, y = f1_score, fill = f1_metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = sprintf("%.3f", f1_score)), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 5) +
    labs(title = title,
         x = "Metric",
         y = "F1 Score",
         fill = "Parameters") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_discrete(limits = metric_order) +
    theme_minimal(base_size = 18, base_family = "sans") +
    theme(
      plot.title = element_text(size = 21.6, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom",
      legend.box = "horiztontal",
      #legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key.size = unit(1, "lines"),  
      legend.text = element_text(size = 12), 
      legend.title = element_text(size = 12)  
    ) +
    scale_fill_brewer(palette = "Set2")
}

    
# Create incell plot
create_bar_plot(
  filter(plot_data, type == "incell"),
  "incell",  "F1 Scores for In Cellulo Prediction"
)

# Create invitro plot
create_bar_plot(
  filter(plot_data, type == "invitro"),
  "invitro", "F1 Scores for In Vitro Prediction"
)  
