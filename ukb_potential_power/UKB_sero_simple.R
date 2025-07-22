# UK Biobank Serology - Simple, Linear Script
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

# --- 1. Load the data ---
# Load the file that has the full titles as headers.
# `read_tsv` will likely create warnings due to the messy format, which is okay.
cat("Loading data from 'serology_export_title.tsv'...\n")
data <- read_tsv("serology_export_title.tsv", show_col_types = FALSE)
cat("Data loaded.\n")

# --- 2. Define a simple map for the columns we want to plot ---
# We map a short, easy name to the exact, long column name from the file.
col_map <- c(
  "hsv1"      = "1gG antigen for Herpes Simplex virus-1 | Instance 0",
  "hsv2"      = "2mgG unique antigen for Herpes Simplex virus-2 | Instance 0",
  "ebv_vca"   = "VCA p18 antigen for Epstein-Barr Virus | Instance 0",
  "vzv"       = "gE / gI antigen for Varicella Zoster Virus  | Instance 0",
  "cmv_pp150" = "pp150 Nter antigen for Human Cytomegalovirus | Instance 0",
  "hp_caga"   = "CagA antigen for Helicobacter pylori | Instance 0",
  "hp_vaca"   = "VacA antigen for Helicobacter pylori | Instance 0",
  "hp_groel"  = "GroEL antigen for Helicobacter pylori | Instance 0",
  "hp_urea"   = "UreA antigen for Helicobacter pylori | Instance 0",
  "toxo_p22"  = "p22 antigen for Toxoplasma gondii | Instance 0"
)

# --- 3. Rename columns and clean the data ---
# We will only keep and process the columns defined above.
# `any_of()` is important because it won't fail if a column name from our map isn't found.
cat("Selecting and renaming columns...\n")
data_clean <- data %>%
  select(any_of(col_map)) %>%
  rename(any_of(col_map)) %>%
  # Now, convert every selected column to numeric.
  # This is the critical step that fixes the 'discrete value' error.
  mutate(across(everything(), as.numeric))

cat("Data cleaned and converted to numeric.\n")

# --- 4. Define thresholds for our simple names ---
thresholds <- c(
  hsv1 = 150, hsv2 = 150, ebv_vca = 250, vzv = 100,
  cmv_pp150 = 100, hp_caga = 400, hp_vaca = 100,
  hp_groel = 80, hp_urea = 130, toxo_p22 = 100
)

# --- 5. Plotting Loop ---
# This can be copy-pasted into an interactive terminal.
all_plots <- list()
cat("Starting to generate plots...\n")

for (col_simple_name in names(data_clean)) {
  
  plot_data <- data_clean %>%
    select(value = !!sym(col_simple_name)) %>%
    filter(!is.na(value) & value > 0)
  
  if (nrow(plot_data) == 0) {
    cat("  -> No valid data for '", col_simple_name, "', skipping.\n", sep="")
    next
  }
  
  cat("  -> Plotting '", col_simple_name, "' with ", nrow(plot_data), " records.\n", sep="")
  
  current_threshold <- thresholds[col_simple_name]
  
  p_linear <- ggplot(plot_data, aes(x = value)) +
    geom_histogram(bins = 50, fill = "navyblue", alpha = 0.7) +
    geom_vline(xintercept = current_threshold, color = "red", size = 1) +
    labs(title = paste(col_simple_name, "(Linear)"), y = NULL) +
    theme_bw()
    
  p_log <- ggplot(plot_data, aes(x = value)) +
    geom_histogram(bins = 50, fill = "darkorange", alpha = 0.7) +
    geom_vline(xintercept = current_threshold, color = "red", size = 1) +
    scale_x_log10() +
    labs(title = paste(col_simple_name, "(Log)"), y = NULL) +
    theme_bw()
    
  all_plots[[col_simple_name]] <- p_linear / p_log
}

# --- 6. Save the final combined plot ---
if (length(all_plots) > 0) {
  cat("Combining", length(all_plots), "plots and saving to PDF...\n")
  combined_plots <- wrap_plots(all_plots, ncol = 2)
  ggsave("serology_histograms.pdf", combined_plots, width = 12, height = length(all_plots) * 2.5, limitsize = FALSE)
  cat("Finished! Plots saved to 'serology_histograms.pdf'\n")
} else {
  cat("Warning: No plots were generated.\n")
} 