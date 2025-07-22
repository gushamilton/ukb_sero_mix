# Test script for UK Biobank serology analysis
library(readr)
library(dplyr)

# Test data loading function
load_antibody_data <- function(field_id, instance = "i0", data_file = "serology_export_name.tsv") {
  field_name <- paste0(field_id, "_", instance)
  
  if (!file.exists(data_file)) {
    stop("Data file not found: ", data_file)
  }
  
  data <- read_tsv(data_file, 
                   col_select = c("eid", !!field_name),
                   show_col_types = FALSE,
                   na = c("", "NA", "False", "True")) %>%
    rename(value = !!field_name) %>%
    filter(!is.na(value)) %>%
    mutate(value = as.numeric(value)) %>%
    filter(value > 0)
  
  cat("Loaded", nrow(data), "observations for", field_name, "\n")
  return(data)
}

# Test with a few key antibodies
test_fields <- c("p23000", "p23001", "p23049", "p23048", "p23026", "p23039")

cat("Testing UK Biobank data loading...\n")

for (field in test_fields) {
  tryCatch({
    data <- load_antibody_data(field)
    cat("Field", field, ":", nrow(data), "observations, range:", 
        round(range(data$value), 3), "\n")
  }, error = function(e) {
    cat("Error with field", field, ":", e$message, "\n")
  })
} 