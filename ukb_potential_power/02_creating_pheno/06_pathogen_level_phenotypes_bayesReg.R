################################################################################
#   UKB SEROLOGY: PATHOGEN-LEVEL PHENOTYPES (BAYESIAN REGRESSION + COMPARISON) - SCRIPT 06
################################################################################
#  ▸ Strategy: Use the single-antigen posterior probabilities (p_soft) generated
#    by script 02 as predictors in a Bayesian logistic regression that
#    approximates the UK Biobank rule (≥2 antigens positive, or pgp3 for CT).
#  ▸ This learns data-driven weights for each antigen while remaining robust and
#    fast – logistic regression is easy for MCMC to sample.
#  ▸ Also generates comparison plots and GWAS outputs with confidence-based weights.
################################################################################

# --- 1. PACKAGES --------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, glue, brms, future, here, bayestestR, patchwork, scales)

plan(multisession, workers = max(1, availableCores() - 1))
options(mc.cores = max(1, availableCores() - 1))

# --- 2. INPUT ----------------------------------------------------------------
pheno_path <- "quickdraws_input/phenotypes_master_postqc.tsv"
if (!file.exists(pheno_path)) stop("Cannot find ", pheno_path, ". Run script 02 first.")
master <- read_tsv(pheno_path, show_col_types = FALSE, guess_max = 5000)
cat("Loaded", nrow(master), "individuals\n")

# --- 3. Pathogen → antigen map ----------------------------------------------
pathogen_map <- list(
  EBV = c("ebv_vca", "ebv_ebna1", "ebv_zebra", "ebv_ead"),
  CMV = c("cmv_pp150", "cmv_pp52", "cmv_pp28"),
  CT  = c("ct_pgp3",  "ct_mompa", "ct_mompd", "ct_tarpf1", "ct_tarpf2"),
  HP  = c("hp_caga",  "hp_vaca", "hp_omp", "hp_groel", "hp_catalase", "hp_urea")
)

# Butler-Laporte et al. 2023 MFI thresholds for exact comparison --------------
butler_laporte_thresholds <- list(
  # EBV antigens
  ebv_vca_p18 = 250,
  ebv_ebna1 = 250, 
  ebv_zebra = 100,
  ebv_ea_d = 100,
  
  # CMV antigens
  cmv_pp150 = 100,
  cmv_pp52 = 150,
  cmv_pp28 = 200,
  
  # H. pylori antigens
  hp_caga = 400,
  hp_vaca = 100,
  hp_omp = 170,
  hp_groel = 80,
  hp_catalase = 180,
  hp_urea = 130,
  
  # C. trachomatis antigens
  ct_mompa = 100,
  ct_mompd = 100,
  ct_tarpf1 = 100,
  ct_tarpf2 = 100,
  ct_porb = 80,
  ct_pgp3 = 200
)

# UKB rule as noisy target (Butler-Laporte et al. 2023) ------------------------
ukb_rules <- list(
  EBV = list(type = "count", n = 2),  # Positive for 2 or more antigens
  CMV = list(type = "count", n = 2),  # Positive for 2 or more antigens
  HP  = list(type = "count_exclude_caga", n = 2),  # Positive for 2 or more antigens, except CagA
  CT  = list(type = "pgp3_or_complex")  # Positive for pGP3 OR negative for pGP3 but positive for 2 out of 5 remaining antigens
)

# Function to create hard cutoffs using Butler-Laporte thresholds
make_butler_laporte_hard <- function(df, antigens, pathogen) {
  hard_cols <- character(length(antigens))
  
  for (i in seq_along(antigens)) {
    antigen <- antigens[i]
    mfi_col <- paste0(antigen, "_mfi")
    threshold <- butler_laporte_thresholds[[antigen]]
    
    if (is.null(threshold)) {
      stop("Missing threshold for antigen: ", antigen)
    }
    
    hard_cols[i] <- paste0(antigen, "_bl_hard")
    df[[hard_cols[i]]] <- as.integer(df[[mfi_col]] >= threshold)
  }
  
  # Apply Butler-Laporte rules
  if (pathogen == "EBV") {
    # Positive for 2 or more antigens
    df[[paste0(tolower(pathogen), "_bl_hard")]] <- 
      as.integer(rowSums(df[hard_cols], na.rm = TRUE) >= 2)
  } else if (pathogen == "CMV") {
    # Positive for 2 or more antigens
    df[[paste0(tolower(pathogen), "_bl_hard")]] <- 
      as.integer(rowSums(df[hard_cols], na.rm = TRUE) >= 2)
  } else if (pathogen == "HP") {
    # Positive for 2 or more antigens, except CagA
    cols_no_caga <- hard_cols[!grepl("caga", hard_cols)]
    df[[paste0(tolower(pathogen), "_bl_hard")]] <- 
      as.integer(rowSums(df[cols_no_caga], na.rm = TRUE) >= 2)
  } else if (pathogen == "CT") {
    # Positive for pGP3 OR negative for pGP3 but positive for 2 out of 5 remaining antigens
    pgp3_col <- "ct_pgp3_bl_hard"
    other_cols <- hard_cols[hard_cols != "ct_pgp3_bl_hard"]
    
    pgp3_positive <- df[[pgp3_col]] == 1
    pgp3_negative_but_others <- (df[[pgp3_col]] == 0) & 
                               (rowSums(df[other_cols], na.rm = TRUE) >= 2)
    
    df[[paste0(tolower(pathogen), "_bl_hard")]] <- 
      as.integer(pgp3_positive | pgp3_negative_but_others)
  }
  
  return(df)
}

make_label <- function(df, antigens, rule){
  if (rule$type == "count") {
    # Standard count rule: positive if ≥n antigens are positive
    cols <- paste0(antigens, "_sero_hard")
    as.integer(rowSums(df[cols], na.rm = TRUE) >= rule$n)
  } else if (rule$type == "count_exclude_caga") {
    # HP rule: positive if ≥2 antigens are positive, excluding CagA
    cols <- paste0(antigens, "_sero_hard")
    # Remove CagA from the count
    cols_no_caga <- cols[!grepl("caga", cols)]
    as.integer(rowSums(df[cols_no_caga], na.rm = TRUE) >= rule$n)
  } else if (rule$type == "pgp3_or_complex") {
    # CT rule: positive for pGP3 OR negative for pGP3 but positive for 2 out of 5 remaining antigens
    pgp3_col <- "ct_pgp3_sero_hard"
    other_cols <- paste0(antigens[antigens != "ct_pgp3"], "_sero_hard")
    
    # Positive if pGP3 is positive
    pgp3_positive <- df[[pgp3_col]] == 1
    
    # OR if pGP3 is negative but ≥2 other antigens are positive
    pgp3_negative_but_others <- (df[[pgp3_col]] == 0) & 
                               (rowSums(df[other_cols], na.rm = TRUE) >= 2)
    
    as.integer(pgp3_positive | pgp3_negative_but_others)
  } else {
    stop("Unknown rule type: ", rule$type)
  }
}

# --- 4. Bayesian Regression + Comparison --------------------------------------
results <- list()
model_summaries <- list()
all_plots <- list()
gwas_outputs <- list()

for (pth in names(pathogen_map)) {
  antigens <- pathogen_map[[pth]]
  soft_cols <- paste0(antigens, "_sero_soft")
  hard_cols <- paste0(antigens, "_sero_hard")

  df <- master %>%
    select(FID, IID, all_of(c(soft_cols, hard_cols))) %>%
    drop_na()
  if (nrow(df) < 100) {
    warning("Skipping ", pth, ": less than 100 complete cases")
    next
  }

  df$z <- make_label(df, antigens, ukb_rules[[pth]])

  # Create Butler-Laporte hard cutoffs for exact comparison
  df <- make_butler_laporte_hard(df, antigens, pth)

  # Use full dataset (no subsampling) -----------------------------------------
  df_model <- df

  # Build formula -----------------------------------------------------------
  rhs <- paste(soft_cols, collapse = " + ")
  brm_formula <- bf(as.formula(paste("z ~", rhs)), family = bernoulli(link = "logit"))

  cat(glue("\nFitting Bayesian logistic model for {pth} on {nrow(df_model)} individuals (full dataset)...\n"))

  fit <- brm(brm_formula, data = df_model,
             chains = 2, iter = 2000, warmup = 500,
             refresh = 0,
             prior = set_prior("normal(0,1)", class = "b"))

  # Store model summary
  model_summaries[[pth]] <- summary(fit)

  # Predict on *all* complete-case rows with credible intervals
  pred_full <- fitted(fit, newdata = df, scale = "response", probs = c(0.025, 0.25, 0.75, 0.975))
  pi_hat <- pred_full[, "Estimate"]
  pi_lower <- pred_full[, "Q2.5"]
  pi_upper <- pred_full[, "Q97.5"]
  pi_ci_width <- pi_upper - pi_lower

  # Multiple serostatus thresholds
  thresholds <- c(0.3, 0.5, 0.7, 0.9)
  serostatus_cols <- list()
  
  for (thresh in thresholds) {
    col_name <- glue("{tolower(pth)}_sero_hard_{thresh*100}")
    serostatus_cols[[col_name]] <- as.integer(pi_hat >= thresh)
  }

  # Calculate seroprevalence using Bayesian posterior samples
  # Get posterior samples for all predictions
  post_samples <- posterior_epred(fit, newdata = df, draws = 1000)
  
  # Calculate seroprevalence for each posterior sample
  seroprev_samples <- apply(post_samples, 1, function(x) mean(x >= 0.5))
  
  # Get credible intervals from posterior distribution
  seroprev_ci <- quantile(seroprev_samples, c(0.025, 0.975))
  seroprevalence <- mean(seroprev_samples)

  results[[pth]] <- tibble(
    FID = df$FID,
    IID = df$IID,
    !!glue("{tolower(pth)}_sero_soft") := pi_hat,
    !!glue("{tolower(pth)}_sero_soft_lower") := pi_lower,
    !!glue("{tolower(pth)}_sero_soft_upper") := pi_upper,
    !!glue("{tolower(pth)}_sero_soft_ci_width") := pi_ci_width,
    !!!serostatus_cols,
    !!glue("{tolower(pth)}_w_soft")    := pi_hat,
    !!glue("{tolower(pth)}_w_hard")    := 2 * abs(pi_hat - 0.5)
  )
  
  # Print summary statistics
  cat(glue("\n--- {pth} Summary ---\n"))
  cat(glue("Seroprevalence (≥0.5): {scales::percent(seroprevalence, accuracy = 0.1)} ({scales::percent(seroprev_ci[1], accuracy = 0.1)} - {scales::percent(seroprev_ci[2], accuracy = 0.1)})\n"))
  cat(glue("Mean CI width: {round(mean(pi_ci_width), 3)}\n"))
  cat(glue("Coefficients:\n"))
  coef_summary <- summary(fit)$fixed
  for (i in 2:nrow(coef_summary)) {  # Skip intercept
    antigen_name <- rownames(coef_summary)[i]
    antigen_clean <- str_remove(antigen_name, "_sero_soft")
    cat(glue("  {antigen_clean}: {round(coef_summary[i, 'Estimate'], 3)} ({round(coef_summary[i, 'l-95% CI'], 3)}, {round(coef_summary[i, 'u-95% CI'], 3)})\n"))
  }

  # --- Create comparison plots and GWAS outputs ------------------------------
  # Create UKB hard labels for comparison
  df$ukb_hard <- make_label(df, antigens, ukb_rules[[pth]])
  
  # Rename Bayesian columns for clarity
  df <- df %>%
    mutate(
      bayes_soft = pi_hat,
      bayes_hard = as.integer(pi_hat >= 0.5)
    )
  
  # Create scatter plots for each antigen
  antigen_plots <- list()
  
  for (i in seq_along(antigens)) {
    antigen <- antigens[i]
    soft_col <- paste0(antigen, "_sero_soft")
    
    # Calculate correlation
    cor_val <- cor(df[[soft_col]], df$bayes_soft, use = "complete.obs")
    
    p <- df %>%
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
  
  # Create GWAS outputs
  gwas_outputs[[pth]] <- df %>%
    select(FID, IID) %>%
    mutate(
      # Butler-Laporte exact hard seroprevalence (for comparison with paper)
      !!glue("{tolower(pth)}_bl_hard") := df[[paste0(tolower(pth), "_bl_hard")]],
      
      # UKB rule hard seroprevalence (noisy target)
      !!glue("{tolower(pth)}_ukb_hard") := df$ukb_hard,
      
      # Best Bayesian seroprevalence (hard)
      !!glue("{tolower(pth)}_bayes_hard") := df$bayes_hard,
      
      # Soft probability
      !!glue("{tolower(pth)}_bayes_soft") := df$bayes_soft,
      
      # Weights for GWAS - confidence-based weights
      # For soft phenotypes: use soft probability as weight
      !!glue("{tolower(pth)}_w_soft") := df$bayes_soft,
      
      # For hard phenotypes: use confidence (distance from 0.5) as weight
      !!glue("{tolower(pth)}_w_hard") := 2 * abs(df$bayes_soft - 0.5)
    )
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

# --- 6. Combine results and write outputs ------------------------------------
combined <- reduce(results, full_join, by = c("FID", "IID"))

if (!dir.exists("quickdraws_input")) dir.create("quickdraws_input")

write_tsv(combined, "quickdraws_input/phenotypes_pathogen_bayesReg.tsv")
cat("\nSaved pathogen-level Bayesian-regression phenotypes: phenotypes_pathogen_bayesReg.tsv\n")

# Write weight files per pathogen --------------------------------------------
for (pth in names(results)) {
  low <- tolower(pth)
  res <- results[[pth]]
  write_tsv(res %>% select(FID, IID, Weight = !!glue("{low}_w_soft")),
            glue("quickdraws_input/{low}_p_soft.weights"))
  write_tsv(res %>% select(FID, IID, Weight = !!glue("{low}_w_hard")),
            glue("quickdraws_input/{low}_sero_hard.weights"))
}

# Write seroprevalence summary
seroprev_summary <- map_dfr(names(results), function(pth) {
  res <- results[[pth]]
  soft_col <- glue("{tolower(pth)}_sero_soft")
  
  tibble(
    Pathogen = pth,
    N_Total = nrow(res),
    N_Positive_50 = sum(res[[glue("{tolower(pth)}_sero_hard_50")]]),
    N_Positive_70 = sum(res[[glue("{tolower(pth)}_sero_hard_70")]]),
    Seroprevalence_50 = N_Positive_50 / N_Total,
    Seroprevalence_70 = N_Positive_70 / N_Total,
    Mean_CI_Width = mean(res[[glue("{tolower(pth)}_sero_soft_ci_width")]])
  )
})

write_tsv(seroprev_summary, "results/bayesian_regression_seroprevalence_summary.tsv")

# --- 7. Create GWAS comparison outputs ---------------------------------------
# Combine all GWAS outputs
gwas_combined <- reduce(gwas_outputs, full_join, by = c("FID", "IID"))

# Write main phenotype file with all measures
write_tsv(gwas_combined, "quickdraws_input/phenotypes_gwas_comparison.tsv")
cat("\nSaved GWAS comparison phenotypes: phenotypes_gwas_comparison.tsv\n")

# Write separate files for each pathogen and measure type
for (pth in names(gwas_outputs)) {
  low <- tolower(pth)
  res <- gwas_outputs[[pth]]
  
  # Butler-Laporte hard phenotypes (exact paper comparison)
  write_tsv(res %>% select(FID, IID, !!glue("{low}_bl_hard")),
            glue("quickdraws_input/{low}_bl_hard.pheno"))
  
  # UKB hard phenotypes
  write_tsv(res %>% select(FID, IID, !!glue("{low}_ukb_hard")),
            glue("quickdraws_input/{low}_ukb_hard.pheno"))
  
  # Bayesian hard phenotypes  
  write_tsv(res %>% select(FID, IID, !!glue("{low}_bayes_hard")),
            glue("quickdraws_input/{low}_bayes_hard.pheno"))
  
  # Bayesian soft phenotypes
  write_tsv(res %>% select(FID, IID, !!glue("{low}_bayes_soft")),
            glue("quickdraws_input/{low}_bayes_soft.pheno"))
}

# Create summary table
summary_table <- map_dfr(names(gwas_outputs), function(pth) {
  res <- gwas_outputs[[pth]]
  low <- tolower(pth)
  
  tibble(
    Pathogen = pth,
    N_Total = nrow(res),
    Butler_Laporte_Seroprev = mean(res[[glue("{low}_bl_hard")]], na.rm = TRUE),
    UKB_Hard_Seroprev = mean(res[[glue("{low}_ukb_hard")]], na.rm = TRUE),
    Bayes_Hard_Seroprev = mean(res[[glue("{low}_bayes_hard")]], na.rm = TRUE),
    Mean_Soft_Prob = mean(res[[glue("{low}_bayes_soft")]], na.rm = TRUE),
    Mean_Soft_Weight = mean(res[[glue("{low}_w_soft")]], na.rm = TRUE),
    Mean_Hard_Weight = mean(res[[glue("{low}_w_hard")]], na.rm = TRUE)
  )
})

write_tsv(summary_table, "results/gwas_comparison_summary.tsv")

cat("Seroprevalence summary saved to results/bayesian_regression_seroprevalence_summary.tsv\n")
cat("GWAS comparison summary saved to results/gwas_comparison_summary.tsv\n")
cat("Weight files written.\n")
cat("Scatter plots saved to results/\n")
cat("Script 06 complete!\n") 