################################################################################
#   UK BIOBANK SEROLOGY ANALYSIS — 2025‑01‑21
################################################################################
#  ∙ Iterates over key antibodies in UK Biobank with genetic signals
#  ∙ Plots distributions on both linear and log scales
#  ∙ Includes seropositivity thresholds from Table 1 of the referenced study
#  ∙ Creates comprehensive visualization for remote analysis
################################################################################

## ---- packages ----------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(ggplot2, dplyr, tibble, patchwork, scales, readr, purrr)

## ---- antibody field mapping and thresholds -----------------------------------
# Based on Table 1 from https://pubmed.ncbi.nlm.nih.gov/33204752/
# Key antibodies with genetic signals and their seropositivity thresholds
antibody_info <- tribble(
  ~field_id, ~antibody_name, ~threshold, ~description, ~antigen_target,
  # Herpesviridae
  "p23049", "HSV1 IgG", 150, "Herpes Simplex Virus 1 IgG", "mgG-1",
  "p23048", "HSV2 IgG", 150, "Herpes Simplex Virus 2 IgG", "mgG-2",
  "p23001", "EBV VCA IgG", 250, "Epstein-Barr Virus VCA IgG", "VCA p18",
  "p23001", "EBV EBNA1 IgG", 250, "Epstein-Barr Virus EBNA1 IgG", "EBNA-1",
  "p23001", "EBV ZEBRA IgG", 100, "Epstein-Barr Virus ZEBRA IgG", "ZEBRA",
  "p23001", "EBV EA-D IgG", 100, "Epstein-Barr Virus EA-D IgG", "EA-D",
  "p23000", "CMV pp150 IgG", 100, "Cytomegalovirus pp150 IgG", "pp150",
  "p23000", "CMV pp52 IgG", 150, "Cytomegalovirus pp52 IgG", "pp52",
  "p23000", "CMV pp28 IgG", 200, "Cytomegalovirus pp28 IgG", "pp28",
  "p23026", "VZV IgG", 100, "Varicella Zoster Virus IgG", "Glycoproteins E and I",
  # Bacteria - H.pylori
  "p23031", "H.pylori CagA IgG", 400, "Helicobacter pylori CagA IgG", "CagA",
  "p23030", "H.pylori VacA IgG", 100, "Helicobacter pylori VacA IgG", "VacA",
  "p23004", "H.pylori GroEL IgG", 80, "Helicobacter pylori GroEL IgG", "GroEL",
  "p23006", "H.pylori UreA IgG", 130, "Helicobacter pylori UreA IgG", "UreA",
  "p23043", "H.pylori OMP IgG", 170, "Helicobacter pylori OMP IgG", "OMP",
  "p23043", "H.pylori Catalase IgG", 180, "Helicobacter pylori Catalase IgG", "Catalase",
  # Parasite
  "p23039", "Toxo p22 IgG", 100, "Toxoplasma gondii p22 IgG", "p22",
  "p23039", "Toxo sag1 IgG", 160, "Toxoplasma gondii sag1 IgG", "sag1",
  # Additional H.pylori antigens (may be in different fields)
  "p23042", "H.pylori NapA IgG", 80, "Helicobacter pylori NapA IgG", "NapA",
  "p23016", "H.pylori HP0305 IgG", 80, "H.pylori HP0305 IgG", "HP0305",
  "p23017", "H.pylori HP1564 IgG", 80, "H.pylori HP1564 IgG", "HP1564",
  "p23025", "H.pylori HP0596 IgG", 80, "H.pylori HP0596 IgG", "HP0596",
  "p23024", "H.pylori HP0875 IgG", 80, "H.pylori HP0875 IgG", "HP0875",
  "p23023", "H.pylori HP0231 IgG", 80, "H.pylori HP0231 IgG", "HP0231",
  # IgA versions (thresholds may be different)
  "p23018", "H.pylori IgA", 80, "Helicobacter pylori IgA", "General IgA",
  "p23022", "H.pylori HP1564 IgA", 80, "H.pylori HP1564 IgA", "HP1564 IgA",
  "p23010", "H.pylori HP0596 IgA", 80, "H.pylori HP0596 IgA", "HP0596 IgA",
  "p23011", "H.pylori HP0875 IgA", 80, "H.pylori HP0875 IgA", "HP0875 IgA",
  "p23027", "H.pylori HP0231 IgA", 80, "H.pylori HP0231 IgA", "HP0231 IgA",
  "p23015", "H.pylori HP0305 IgA", 80, "H.pylori HP0305 IgA", "HP0305 IgA",
  "p23029", "H.pylori VacA IgA", 80, "H.pylori VacA IgA", "VacA IgA",
  "p23032", "H.pylori CagA IgA", 80, "H.pylori CagA IgA", "CagA IgA",
  "p23014", "H.pylori UreA IgA", 80, "H.pylori UreA IgA", "UreA IgA",
  "p23028", "H.pylori GroEL IgA", 80, "H.pylori GroEL IgA", "GroEL IgA",
  "p23019", "H.pylori NapA IgA", 80, "H.pylori NapA IgA", "NapA IgA"
)

# Focus on the key antibodies with genetic signals
# Note: Some antibodies may be measured in the same field but different antigens
key_antibodies <- antibody_info %>%
  distinct(antibody_name, antigen_target, .keep_all = TRUE) %>%
  arrange(antibody_name, antigen_target)

print("Key antibodies with genetic signals:")
print(key_antibodies %>% select(antibody_name, antigen_target, field_id, threshold))

## ---- data loading function ---------------------------------------------------
load_antibody_data <- function(field_id, instance = "i0", data_file = "serology_export_name.tsv") {
  # Load data from UK Biobank serology export file
  field_name <- paste0(field_id, "_", instance)
  
  # Read the data file
  if (!file.exists(data_file)) {
    stop("Data file not found: ", data_file)
  }
  
  # Read only the columns we need for efficiency
  data <- read_tsv(data_file, 
                   col_select = c("eid", !!field_name),
                   show_col_types = FALSE,
                   na = c("", "NA", "False", "True")) %>%
    rename(value = !!field_name) %>%
    filter(!is.na(value)) %>%
    mutate(value = as.numeric(value))
  
  # Remove any non-positive values (antibody levels should be positive)
  data <- data %>% filter(value > 0)
  
  cat("  Loaded", nrow(data), "observations for", field_name, "\n")
  
  return(data)
}

## ---- plotting functions -----------------------------------------------------
create_antibody_plot <- function(data, antibody_name, antigen_target, threshold) {
  # Create linear scale plot
  p1 <- ggplot(data, aes(x = value)) +
    geom_histogram(aes(y = ..density..), bins = 50, 
                   fill = "#4E79A7", alpha = 0.7, color = "white") +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1) +
    labs(title = paste("Linear Scale:", antibody_name),
         subtitle = paste("Antigen:", antigen_target),
         x = "Antibody Level (MFI)", y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8))
  
  # Create log scale plot
  p2 <- ggplot(data, aes(x = value)) +
    geom_histogram(aes(y = ..density..), bins = 50,
                   fill = "#F28E2B", alpha = 0.7, color = "white") +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1) +
    scale_x_log10(labels = comma) +
    labs(title = paste("Log Scale:", antibody_name),
         subtitle = paste("Antigen:", antigen_target),
         x = "Antibody Level (MFI, log scale)", y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8))
  
  # Combine plots
  p1 / p2
}

## ---- main analysis function -------------------------------------------------
analyze_ukb_serology <- function(data_file = "serology_export_name.tsv") {
  cat("Starting UK Biobank serology analysis...\n")
  cat("Using data file:", data_file, "\n")
  
  # Check if data file exists
  if (!file.exists(data_file)) {
    stop("Data file not found: ", data_file, "\nPlease ensure the file is in the current directory.")
  }
  
  # Create list to store plots
  all_plots <- list()
  successful_analyses <- 0
  
  # Analyze each key antibody
  for (i in 1:nrow(key_antibodies)) {
    antibody <- key_antibodies[i, ]
    cat("Processing:", antibody$antibody_name, "(", antibody$field_id, ")\n")
    
    # Try to load data with error handling
    tryCatch({
      # Load data
      data <- load_antibody_data(antibody$field_id, data_file = data_file)
      
      # Skip if no data
      if (nrow(data) == 0) {
        cat("  No data available for", antibody$field_id, "- skipping\n")
        next
      }
      
      # Create plot
      plot <- create_antibody_plot(data, antibody$antibody_name, antibody$antigen_target, antibody$threshold)
      
      # Store plot
      all_plots[[successful_analyses + 1]] <- plot
      
      # Add antibody info to plot
      all_plots[[successful_analyses + 1]] <- all_plots[[successful_analyses + 1]] + 
        plot_annotation(
          title = paste("Field:", antibody$field_id, "-", antibody$antibody_name),
          subtitle = paste("Antigen:", antibody$antigen_target, "| Threshold:", antibody$threshold, "| N:", nrow(data)),
          theme = theme(plot.title = element_text(size = 12, face = "bold"))
        )
      
      successful_analyses <- successful_analyses + 1
      
    }, error = function(e) {
      cat("  Error processing", antibody$field_id, ":", e$message, "\n")
    })
  }
  
  cat("Successfully analyzed", successful_analyses, "out of", nrow(key_antibodies), "antibodies\n")
  
  if (successful_analyses == 0) {
    stop("No antibodies could be analyzed. Please check the data file and field IDs.")
  }
  
  # Create comprehensive plot
  cat("Creating comprehensive visualization...\n")
  
  # Calculate layout
  n_plots <- length(all_plots)
  if (n_plots == 0) {
    stop("No plots were created. Please check the data and field IDs.")
  }
  n_cols <- 3
  n_rows <- ceiling(n_plots / n_cols)
  
  # Combine all plots
  combined_plot <- wrap_plots(all_plots, ncol = n_cols, nrow = n_rows) +
    plot_annotation(
      title = "UK Biobank Antibody Distributions",
      subtitle = "Linear and Log Scale Views with Seropositivity Thresholds",
      caption = "Red dashed lines indicate seropositivity thresholds from Table 1",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        plot.caption = element_text(size = 10, hjust = 0)
      )
    )
  
  # Save plot
  cat("Saving plot to 'ukb_antibody_distributions.pdf'...\n")
  ggsave("ukb_antibody_distributions.pdf", 
         combined_plot, 
         width = 18, 
         height = 6 * n_rows, 
         limitsize = FALSE)
  
  # Create summary statistics
  cat("Creating summary statistics...\n")
  summary_stats <- tibble()
  
  for (i in 1:nrow(key_antibodies)) {
    antibody <- key_antibodies[i, ]
    
    tryCatch({
      data <- load_antibody_data(antibody$field_id, data_file = data_file)
      
      if (nrow(data) > 0) {
        stats <- tibble(
          antibody_name = antibody$antibody_name,
          antigen_target = antibody$antigen_target,
          field_id = antibody$field_id,
          threshold = antibody$threshold,
          n_total = nrow(data),
          n_seropositive = sum(data$value > antibody$threshold),
          seroprevalence = n_seropositive / n_total,
          mean_level = mean(data$value),
          median_level = median(data$value),
          sd_level = sd(data$value),
          min_level = min(data$value),
          max_level = max(data$value)
        )
        summary_stats <- bind_rows(summary_stats, stats)
      }
    }, error = function(e) {
      cat("  Error creating stats for", antibody$field_id, ":", e$message, "\n")
    })
  }
  
  # Save summary statistics
  write_csv(summary_stats, "ukb_antibody_summary.csv")
  
  cat("Analysis complete!\n")
  cat("Files created:\n")
  cat("- ukb_antibody_distributions.pdf: Comprehensive visualization\n")
  cat("- ukb_antibody_summary.csv: Summary statistics\n")
  
  return(list(
    plots = all_plots,
    summary = summary_stats,
    combined_plot = combined_plot
  ))
}

## ---- run analysis ----------------------------------------------------------
# Run the analysis
results <- analyze_ukb_serology()

# Print summary
print("Summary of key antibodies analyzed:")
print(results$summary %>% 
      select(antibody_name, antigen_target, field_id, seroprevalence, mean_level, threshold) %>%
      arrange(desc(seroprevalence))) 