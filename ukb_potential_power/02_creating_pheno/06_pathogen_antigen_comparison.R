################################################################################
#   UKB SEROLOGY: PATHOGEN-ANTIGEN COMPARISON & GWAS OUTPUTS - SCRIPT 08
################################################################################
#  ▸ Strategy: Compare individual antigen probabilities with pathogen-level
#    Bayesian regression probabilities. Generate scatter plots and output
#    multiple serostatus measures for GWAS analysis.
#  ▸ Outputs: Butler-Laporte/UKB hard seroprevalence, Bayesian hard seroprevalence,
#    soft probabilities, and appropriate weights for GWAS.
################################################################################

# --- 1. PACKAGES --------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, glue, patchwork, scales)

# --- 2. INPUT ----------------------------------------------------------------
pheno_path <- "quickdraws_input/phenotypes_master_postqc.tsv"
if (!file.exists(pheno_path)) stop("Cannot find ", pheno_path, ". Run script 02 first.")
master <- read_tsv(pheno_path, show_col_types = FALSE, guess_max = 5000)

# Load Bayesian regression results
bayes_path <- "quickdraws_input/phenotypes_pathogen_bayesReg.tsv"
if (!file.exists(bayes_path)) stop("Cannot find ", bayes_path, ". Run script 06 first.")
bayes_results <- read_tsv(bayes_path, show_col_types = FALSE)

cat("Loaded", nrow(master), "individuals\n")

# --- 3. Pathogen → antigen map ----------------------------------------------
pathogen_map <- list(
  EBV = c("ebv_vca", "ebv_ebna1", "ebv_zebra", "ebv_ead"),
  CMV = c("cmv_pp150", "cmv_pp52", "cmv_pp28"),
  CT  = c("ct_pgp3",  "ct_mompa", "ct_mompd", "ct_tarpf1", "ct_tarpf2"),
  HP  = c("hp_caga",  "hp_vaca", "hp_omp", "hp_groel", "hp_catalase", "hp_urea")
)

# UKB rule definitions --------------------------------------------------------
ukb_rules <- list(
  EBV = list(type = "count", n = 2),
  CMV = list(type = "count", n = 2),
  HP  = list(type = "count", n = 2),
  CT  = list(type = "pgp3")
)

make_ukb_label <- function(df, antigens, rule){
  if (rule$type == "count") {
    cols <- paste0(antigens, "_sero_hard")
    hard_matrix <- df[cols]
    hard_matrix[is.na(hard_matrix)] <- 0
    as.integer(rowSums(hard_matrix, na.rm = TRUE) >= rule$n)
  } else { # pgp3 rule
    ifelse(is.na(df$ct_pgp3_sero_hard), 0, df$ct_pgp3_sero_hard)
  }
}

# --- 4. Create comparison plots and outputs ----------------------------------
all_plots <- list()
gwas_outputs <- list()

for (pth in names(pathogen_map)) {
  antigens <- pathogen_map[[pth]]
  soft_cols <- paste0(antigens, "_sero_soft")
  hard_cols <- paste0(antigens, "_sero_hard")
  
  # Get data for this pathogen
  pathogen_data <- master %>%
    select(FID, IID, all_of(c(soft_cols, hard_cols))) %>%
    left_join(
      bayes_results %>% select(FID, IID, 
                              !!glue("{tolower(pth)}_sero_soft"),
                              !!glue("{tolower(pth)}_sero_hard_50")),
      by = c("FID", "IID")
    ) %>%
    drop_na(!!glue("{tolower(pth)}_sero_soft"))  # Keep only those with Bayesian results
  
  if (nrow(pathogen_data) == 0) {
    warning("No data for ", pth, " - skipping")
    next
  }
  
  # Create UKB hard labels
  pathogen_data$ukb_hard <- make_ukb_label(pathogen_data, antigens, ukb_rules[[pth]])
  
  # Rename Bayesian columns for clarity
  pathogen_data <- pathogen_data %>%
    rename(
      bayes_soft = !!glue("{tolower(pth)}_sero_soft"),
      bayes_hard = !!glue("{tolower(pth)}_sero_hard_50")
    )
  
  # Create scatter plots for each antigen
  antigen_plots <- list()
  
  for (i in seq_along(antigens)) {
    antigen <- antigens[i]
    soft_col <- paste0(antigen, "_sero_soft")
    
    # Calculate correlation
    cor_val <- cor(pathogen_data[[soft_col]], pathogen_data$bayes_soft, use = "complete.obs")
    
    p <- pathogen_data %>%
      ggplot(aes(x = !!sym(soft_col), y = bayes_soft)) +
      geom_point(alpha = 0.6, size = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(
        title = glue("{pth}: {antigen} vs Pathogen-Level"),
        subtitle = glue("r = {round(cor_val, 3)}"),
        x = glue("{antigen} Soft Probability"),
        y = "Pathogen-Level Soft Probability"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 9))
    
    antigen_plots[[antigen]] <- p
  }
  
  # Combine antigen plots
  if (length(antigen_plots) <= 4) {
    combined_plot <- wrap_plots(antigen_plots, ncol = 2)
  } else {
    combined_plot <- wrap_plots(antigen_plots, ncol = 3)
  }
  
  all_plots[[pth]] <- combined_plot + 
    plot_annotation(title = glue("{pth}: Individual Antigens vs Pathogen-Level Probabilities"))
  
  # Calculate summary statistics
  ukb_seroprev <- mean(pathogen_data$ukb_hard, na.rm = TRUE)
  bayes_seroprev <- mean(pathogen_data$bayes_hard, na.rm = TRUE)
  mean_soft <- mean(pathogen_data$bayes_soft, na.rm = TRUE)
  
  # Create GWAS outputs
  gwas_outputs[[pth]] <- pathogen_data %>%
    select(FID, IID) %>%
    mutate(
      # Butler-Laporte/UKB hard seroprevalence
      !!glue("{tolower(pth)}_ukb_hard") := pathogen_data$ukb_hard,
      
      # Best Bayesian seroprevalence (hard)
      !!glue("{tolower(pth)}_bayes_hard") := pathogen_data$bayes_hard,
      
      # Soft probability
      !!glue("{tolower(pth)}_bayes_soft") := pathogen_data$bayes_soft,
      
      # Weights for GWAS - confidence-based weights
      # For soft phenotypes: use soft probability as weight
      !!glue("{tolower(pth)}_w_soft") := pathogen_data$bayes_soft,
      
      # For hard phenotypes: use confidence (distance from 0.5) as weight
      !!glue("{tolower(pth)}_w_hard") := 2 * abs(pathogen_data$bayes_soft - 0.5)
    )
  
  # Print summary
  cat(glue("\n--- {pth} Summary ---\n"))
  cat(glue("UKB Hard Seroprevalence: {scales::percent(ukb_seroprev, accuracy = 0.1)}\n"))
  cat(glue("Bayesian Hard Seroprevalence: {scales::percent(bayes_seroprev, accuracy = 0.1)}\n"))
  cat(glue("Mean Soft Probability: {round(mean_soft, 3)}\n"))
  cat(glue("Sample size: {nrow(pathogen_data)}\n"))
}

# --- 5. Save plots -----------------------------------------------------------
if (!dir.exists("results")) dir.create("results")

# Save individual pathogen plots
for (pth in names(all_plots)) {
  ggsave(
    glue("results/{tolower(pth)}_antigen_vs_pathogen_scatters.pdf"),
    all_plots[[pth]],
    width = 12, height = 8, dpi = 300
  )
}

# Combine all plots into one file
if (length(all_plots) > 0) {
  all_combined <- wrap_plots(all_plots, ncol = 1) +
    plot_annotation(title = "All Pathogens: Antigen vs Pathogen-Level Comparisons")
  
  ggsave("results/all_pathogens_antigen_comparisons.pdf", all_combined,
         width = 14, height = 10, dpi = 300)
}

# --- 6. Create GWAS output files ---------------------------------------------
if (!dir.exists("quickdraws_input")) dir.create("quickdraws_input")

# Combine all GWAS outputs
gwas_combined <- reduce(gwas_outputs, full_join, by = c("FID", "IID"))

# Write main phenotype file with all measures
write_tsv(gwas_combined, "quickdraws_input/phenotypes_gwas_comparison.tsv")
cat("\nSaved GWAS comparison phenotypes: phenotypes_gwas_comparison.tsv\n")

# Write separate files for each pathogen and measure type
for (pth in names(gwas_outputs)) {
  low <- tolower(pth)
  res <- gwas_outputs[[pth]]
  
  # UKB hard phenotypes
  write_tsv(res %>% select(FID, IID, Phenotype = !!glue("{low}_ukb_hard")),
            glue("quickdraws_input/{low}_ukb_hard.pheno"))
  
  # Bayesian hard phenotypes  
  write_tsv(res %>% select(FID, IID, Phenotype = !!glue("{low}_bayes_hard")),
            glue("quickdraws_input/{low}_bayes_hard.pheno"))
  
  # Bayesian soft phenotypes
  write_tsv(res %>% select(FID, IID, Phenotype = !!glue("{low}_bayes_soft")),
            glue("quickdraws_input/{low}_bayes_soft.pheno"))
  
  # Weight files - confidence-based weights
  write_tsv(res %>% select(FID, IID, Weight = !!glue("{low}_w_soft")),
            glue("quickdraws_input/{low}_p_soft.weights"))
  write_tsv(res %>% select(FID, IID, Weight = !!glue("{low}_w_hard")),
            glue("quickdraws_input/{low}_sero_hard.weights"))
}

# Create summary table
summary_table <- map_dfr(names(gwas_outputs), function(pth) {
  res <- gwas_outputs[[pth]]
  low <- tolower(pth)
  
  tibble(
    Pathogen = pth,
    N_Total = nrow(res),
    UKB_Hard_Seroprev = mean(res[[glue("{low}_ukb_hard")]], na.rm = TRUE),
    Bayes_Hard_Seroprev = mean(res[[glue("{low}_bayes_hard")]], na.rm = TRUE),
    Mean_Soft_Prob = mean(res[[glue("{low}_bayes_soft")]], na.rm = TRUE),
    Mean_Soft_Weight = mean(res[[glue("{low}_w_soft")]], na.rm = TRUE),
    Mean_Hard_Weight = mean(res[[glue("{low}_w_hard")]], na.rm = TRUE)
  )
})

write_tsv(summary_table, "results/gwas_comparison_summary.tsv")
cat("Summary table saved to results/gwas_comparison_summary.tsv\n")

cat("\nScript 08 complete!\n")
cat("Generated:\n")
cat("- Scatter plots comparing individual antigens vs pathogen-level probabilities\n")
cat("- GWAS phenotype files with UKB hard, Bayesian hard, and soft measures\n")
cat("- Confidence-based weight files for GWAS analysis\n")
cat("- Summary table with seroprevalence comparisons\n") 