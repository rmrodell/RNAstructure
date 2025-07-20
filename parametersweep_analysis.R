library(dplyr)
library(tidyr)


data <- read.csv("parametersweep_paralleled_output/all_analysis_results.csv")

# total number of modified sites
incell_total <- 178
invitro_total <- 217

#fl scores
f1_UNUAR_incell <- 0.555555555555556
f1_UGUAG_invitro <- 0.754990925589837
f1_UGUAG_incell <- 0.652818991097923
f1_UNUAR_invitro <- 0.671875


# function to analyze all rows with top value in a column
analyze_top_values <- function(data, column_name) {
  top_values <- data %>%
    arrange(desc(!!sym(column_name))) %>%
    mutate(rank = row_number(),
           max_value = max(!!sym(column_name))) %>%
    filter(!!sym(column_name) == max_value) %>%
    select(-rank, -max_value)
  
  # Get the top value
  top_value <- top_values %>% pull(!!sym(column_name)) %>% max(na.rm = TRUE)
  
  range_summary <- top_values %>%
    summarize(
      across(
        c("total_matches", 
          "incell_count", "invitro_count",
          "unpaired1_min", "unpaired1_max", 
          "paired1_min", "paired1_max", 
          "unpaired2_min", "unpaired2_max", 
          "paired2_min", "paired2_max"),
        list(
          min = ~min(.x, na.rm = TRUE),
          max = ~max(.x, na.rm = TRUE),
          na_count = ~sum(is.na(.x)),
          non_na_count = ~sum(!is.na(.x))
        ),
        .names = "{.col}_{.fn}"
      ),
      include_paired2_true = sum(include_paired2),
      include_paired2_false = sum(!include_paired2),
      include_paired2_total = n()
    )
  
  numeric_summary <- range_summary %>%
    select(-starts_with("include_paired2"), -ends_with("na_count"), -ends_with("non_na_count")) %>%
    pivot_longer(
      everything(),
      names_to = c("column", ".value"),
      names_pattern = "(.+)_(min|max)"
    )
  
  logical_summary <- range_summary %>%
    select(starts_with("include_paired2")) %>%
    mutate(
      include_paired2_true_percent = include_paired2_true / include_paired2_total * 100,
      include_paired2_false_percent = include_paired2_false / include_paired2_total * 100
    )
  
  # Print summaries
  cat(sprintf("Top value for %s: %s\n\n", column_name, top_value))
  cat("Numeric Summary:\n")
  print(numeric_summary, n = Inf)
  cat("\nLogical Summary:\n")
  print(logical_summary)
  
  # Return results
  return(list(
    top_values = top_values,
    numeric_summary = numeric_summary,
    logical_summary = logical_summary
  ))
}

# top scores just based on count
top_incellcount <- analyze_top_values(data, "incell_count")
top_invitrocount <- analyze_top_values(data, "invitro_count")

# top scores based on f1 for just the structure
top_incell_f1structure <- analyze_top_values(data, "f1_incell_structure")
top_invitro_f1structure <- analyze_top_values(data, "f1_invitro_structure")

# top scores based on f1 for structure + UNUAR
top_incell_f1UNUARstructure <- analyze_top_values(data, "f1_incell_unuar_and_structure")
top_invitro_f1UNUARstructure <- analyze_top_values(data, "f1_invitro_unuar_and_structure")

# top scores based on f1 for structure + UGUAG
top_incell_f1UGUAGstructure <- analyze_top_values(data, "f1_incell_uguag_and_structure")
top_invitro_f1UGUAGstructure <- analyze_top_values(data, "f1_invitro_uguag_and_structure")



# exclude values where f1 = 1
library(dplyr)
library(tidyr)

analyze_top_values_exclude_one <- function(data, column_name) {
  top_values <- data %>%
    filter(!!sym(column_name) < 1) %>%  # Exclude values of 1
    arrange(desc(!!sym(column_name))) %>%
    mutate(rank = row_number(),
           max_value = max(!!sym(column_name), na.rm = TRUE)) %>%
    filter(!!sym(column_name) == max_value) %>%
    select(-rank, -max_value)
  
  # Get the top value
  top_value <- top_values %>% pull(!!sym(column_name)) %>% max(na.rm = TRUE)
  
  range_summary <- top_values %>%
    summarize(
      across(
        c("total_matches", "unpaired1_min", "unpaired1_max", 
          "paired1_min", "paired1_max", 
          "unpaired2_min", "unpaired2_max", 
          "paired2_min", "paired2_max"),
        list(
          min = ~min(.x, na.rm = TRUE),
          max = ~max(.x, na.rm = TRUE),
          na_count = ~sum(is.na(.x)),
          non_na_count = ~sum(!is.na(.x))
        ),
        .names = "{.col}_{.fn}"
      ),
      include_paired2_true = sum(include_paired2),
      include_paired2_false = sum(!include_paired2),
      include_paired2_total = n()
    )
  
  numeric_summary <- range_summary %>%
    select(-starts_with("include_paired2"), -ends_with("na_count"), -ends_with("non_na_count")) %>%
    pivot_longer(
      everything(),
      names_to = c("column", ".value"),
      names_pattern = "(.+)_(min|max)"
    )
  
  logical_summary <- range_summary %>%
    select(starts_with("include_paired2")) %>%
    mutate(
      include_paired2_true_percent = include_paired2_true / include_paired2_total * 100,
      include_paired2_false_percent = include_paired2_false / include_paired2_total * 100
    )
  
  # Print summaries
  cat(sprintf("Top value for %s (excluding 1): %s\n\n", column_name, top_value))
  cat("Numeric Summary:\n")
  print(numeric_summary, n = Inf)
  cat("\nLogical Summary:\n")
  print(logical_summary)
  
  # Return results
  return(list(
    top_value = top_value,
    top_values = top_values,
    numeric_summary = numeric_summary,
    logical_summary = logical_summary
  ))
}

# top scores based on f1 for structure + UNUAR
top_not1_incell_f1UNUARstructure <- analyze_top_values_exclude_one(data, "f1_incell_unuar_and_structure")
top_not1_invitro_f1UNUARstructure <- analyze_top_values_exclude_one(data, "f1_invitro_unuar_and_structure")

# top scores based on f1 for structure + UGUAG
top_not1_incell_f1UGUAGstructure <- analyze_top_values_exclude_one(data, "f1_incell_uguag_and_structure")
top_not1_invitro_f1UGUAGstructure <- analyze_top_values_exclude_one(data, "f1_invitro_uguag_and_structure")
