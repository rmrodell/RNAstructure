

# Define parameter ranges
PARAMETER_RANGES <- list(
  input_position = 59, 
  offset_min = 0:0,
  offset_max = 1:2,
  include_unpaired1 = c(FALSE, TRUE),
  paired1_min = 1:3,
  paired1_max = 4:9,
  unpaired2_min = 1:3,
  unpaired2_max = 4:9,
  include_paired2 = c(FALSE, TRUE)
)

UNPAIRED1_RANGES <- list(
  unpaired1_min = 1:3,
  unpaired1_max = 4:9
)

PAIRED2_RANGES <- list(
  paired2_min = 1:2,
  paired2_max = 4:9
)

generate_parameter_combinations <- function() {
  base_combinations <- expand.grid(PARAMETER_RANGES)
  
  # Generate all possible combinations for unpaired1 and paired2
  unpaired1_combinations <- expand.grid(UNPAIRED1_RANGES)
  paired2_combinations <- expand.grid(PAIRED2_RANGES)
  
  # Merge all combinations
  all_combinations <- merge(base_combinations, unpaired1_combinations)
  all_combinations <- merge(all_combinations, paired2_combinations)
  
  # Set values to NA where the region is not included
  all_combinations$unpaired1_min[!all_combinations$include_unpaired1] <- NA
  all_combinations$unpaired1_max[!all_combinations$include_unpaired1] <- NA
  all_combinations$paired2_min[!all_combinations$include_paired2] <- NA
  all_combinations$paired2_max[!all_combinations$include_paired2] <- NA
  
  valid_combinations <- all_combinations[
    (is.na(all_combinations$unpaired1_min) | all_combinations$unpaired1_min <= all_combinations$unpaired1_max) &
      all_combinations$paired1_min <= all_combinations$paired1_max &
      all_combinations$unpaired2_min <= all_combinations$unpaired2_max &
      all_combinations$offset_min <= all_combinations$offset_max &
      (is.na(all_combinations$paired2_min) | all_combinations$paired2_min <= all_combinations$paired2_max),
  ]
  
  # Define the desired column order
  column_order <- c(
    "input_position", "offset_min", "offset_max",
    "include_unpaired1", "unpaired1_min", "unpaired1_max",
    "paired1_min", "paired1_max",
    "unpaired2_min", "unpaired2_max",
    "include_paired2", "paired2_min", "paired2_max"
  )
  
  # Reorder the columns
  valid_combinations <- valid_combinations[, column_order]
  
  return(unique(valid_combinations))
}


combinations <- generate_parameter_combinations()

nrow(combinations)

#set 1
PARAMETER_RANGES <- list(
  input_position = 59, 
  offset_min = 0:0,
  offset_max = 1:2,
  include_unpaired1 = c(FALSE, TRUE),
  paired1_min = 1:3,
  paired1_max = 4:9,
  unpaired2_min = 1:3,
  unpaired2_max = 4:4,
  include_paired2 = c(FALSE, TRUE)
)

UNPAIRED1_RANGES <- list(
  unpaired1_min = 1:3,
  unpaired1_max = 4:9
)

PAIRED2_RANGES <- list(
  paired2_min = 1:2,
  paired2_max = 4:9
)

#set2
PARAMETER_RANGES <- list(
  input_position = 59, 
  offset_min = 0:0,
  offset_max = 1:2,
  include_unpaired1 = c(FALSE, TRUE),
  paired1_min = 1:3,
  paired1_max = 4:9,
  unpaired2_min = 1:3,
  unpaired2_max = 5:5,
  include_paired2 = c(FALSE, TRUE)
)

UNPAIRED1_RANGES <- list(
  unpaired1_min = 1:3,
  unpaired1_max = 4:9
)

PAIRED2_RANGES <- list(
  paired2_min = 1:2,
  paired2_max = 4:9
)

#set3
PARAMETER_RANGES <- list(
  input_position = 59, 
  offset_min = 0:0,
  offset_max = 1:2,
  include_unpaired1 = c(FALSE, TRUE),
  paired1_min = 1:3,
  paired1_max = 4:9,
  unpaired2_min = 1:3,
  unpaired2_max = 6:6,
  include_paired2 = c(FALSE, TRUE)
)

UNPAIRED1_RANGES <- list(
  unpaired1_min = 1:3,
  unpaired1_max = 4:9
)

PAIRED2_RANGES <- list(
  paired2_min = 1:2,
  paired2_max = 4:9
)

#set4
PARAMETER_RANGES <- list(
  input_position = 59, 
  offset_min = 0:0,
  offset_max = 1:2,
  include_unpaired1 = c(FALSE, TRUE),
  paired1_min = 1:3,
  paired1_max = 4:9,
  unpaired2_min = 1:3,
  unpaired2_max = 7:7,
  include_paired2 = c(FALSE, TRUE)
)

UNPAIRED1_RANGES <- list(
  unpaired1_min = 1:3,
  unpaired1_max = 4:9
)

PAIRED2_RANGES <- list(
  paired2_min = 1:2,
  paired2_max = 4:9
)

#set5
PARAMETER_RANGES <- list(
  input_position = 59, 
  offset_min = 0:0,
  offset_max = 1:2,
  include_unpaired1 = c(FALSE, TRUE),
  paired1_min = 1:3,
  paired1_max = 4:9,
  unpaired2_min = 1:3,
  unpaired2_max = 8:8,
  include_paired2 = c(FALSE, TRUE)
)

UNPAIRED1_RANGES <- list(
  unpaired1_min = 1:3,
  unpaired1_max = 4:9
)

PAIRED2_RANGES <- list(
  paired2_min = 1:2,
  paired2_max = 4:9
)

#set6
PARAMETER_RANGES <- list(
  input_position = 59, 
  offset_min = 0:0,
  offset_max = 1:2,
  include_unpaired1 = c(FALSE, TRUE),
  paired1_min = 1:3,
  paired1_max = 4:9,
  unpaired2_min = 1:3,
  unpaired2_max = 9:9,
  include_paired2 = c(FALSE, TRUE)
)

UNPAIRED1_RANGES <- list(
  unpaired1_min = 1:3,
  unpaired1_max = 4:9
)

PAIRED2_RANGES <- list(
  paired2_min = 1:2,
  paired2_max = 4:9
)
