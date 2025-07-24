################################################################################
#   UKB SEROLOGY: PATHOGEN-LEVEL PHENOTYPES (BAYESIAN REGRESSION) - SCRIPT 06
################################################################################
#  ▸ Strategy: use the single-antigen posterior probabilities (p_soft) generated
#    by script 02 as predictors in a Bayesian logistic regression that
#    approximates the UK Biobank rule (≥2 antigens positive, or pgp3 for CT).
#  ▸ This learns data-driven weights for each antigen while remaining robust and
#    fast – logistic regression is easy for MCMC to sample.
#  ▸ For development the script subsamples ≤2,000 complete cases per pathogen.
#    Comment out the slice_sample() line to run on the full dataset.
################################################################################

# --- 1. PACKAGES --------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, glue, brms, future, here, bayestestR, patchwork)

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

# UKB rule as noisy target -----------------------------------------------------
ukb_rules <- list(
  EBV = list(type = "count", n = 2),
  CMV = list(type = "count", n = 2),
  HP  = list(type = "count", n = 2),
  CT  = list(type = "pgp3")
)

make_label <- function(df, antigens, rule){
  if (rule$type == "count") {
    cols <- paste0(antigens, "_sero_hard")
    as.integer(rowSums(df[cols], na.rm = TRUE) >= rule$n)
  } else { # pgp3 rule
    df$ct_pgp3_sero_hard
  }
}

# --- 4. Loop over pathogens ---------------------------------------------------
results <- list()
model_summaries <- list()
diagnostic_plots <- list()

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

  # Subsample for development speed -----------------------------------------
  set.seed(42)
  if (nrow(df) > 2000) df_model <- slice_sample(df, n = 200) else df_model <- df

  # Build formula -----------------------------------------------------------
  rhs <- paste(soft_cols, collapse = " + ")
  brm_formula <- bf(as.formula(paste("z ~", rhs)), family = bernoulli(link = "logit"))

  cat("\nFitting Bayesian logistic model for", pth, "on", nrow(df_model), "individuals...\n")

  fit <- brm(brm_formula, data = df_model,
             chains = 2, iter = 2000, warmup = 500,
             refresh = 0,
             prior = set_prior("normal(0,1)", class = "b"))

  # Store model summary
  model_summaries[[pth]] <- summary(fit)
  
  # Create diagnostic plots
  p1 <- pp_check(fit, type = "dens_overlay", ndraws = 100) + 
    ggtitle(glue("{pth}: Posterior Predictive Check"))
  
  # Get coefficient plots - brms plot() returns a list, so we need to handle it properly
  coef_plots <- plot(fit, variable = "b_", regex = TRUE)
  if (is.list(coef_plots)) {
    p2 <- coef_plots[[1]] + ggtitle(glue("{pth}: Coefficient Posteriors"))
  } else {
    p2 <- coef_plots + ggtitle(glue("{pth}: Coefficient Posteriors"))
  }
  
  diagnostic_plots[[pth]] <- p1 + p2 + plot_layout(ncol = 2)

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
  cat(glue("Seroprevalence (≥0.5): {seroprevalence:.1%} ({seroprev_ci[1]:.1%} - {seroprev_ci[2]:.1%})\n"))
  cat(glue("Mean CI width: {mean(pi_ci_width):.3f}\n"))
  cat(glue("Coefficients:\n"))
  coef_summary <- summary(fit)$fixed
  for (i in 2:nrow(coef_summary)) {  # Skip intercept
    antigen_name <- rownames(coef_summary)[i]
    antigen_clean <- str_remove(antigen_name, "_sero_soft")
    cat(glue("  {antigen_clean}: {coef_summary[i, 'Estimate']:.3f} ({coef_summary[i, 'l-95% CI']:.3f}, {coef_summary[i, 'u-95% CI']:.3f})\n"))
  }
}

# --- 5. Create comprehensive diagnostic plots --------------------------------
cat("\nCreating diagnostic plots...\n")

# Combined diagnostic plot
all_diagnostics <- wrap_plots(diagnostic_plots, ncol = 2)
ggsave("results/bayesian_regression_diagnostics.pdf", all_diagnostics, 
       width = 12, height = 8, dpi = 300)

# Distribution plots for each pathogen
dist_plots <- list()
for (pth in names(results)) {
  p <- results[[pth]] %>%
    ggplot(aes(x = !!sym(glue("{tolower(pth)}_sero_soft")))) +
    geom_histogram(bins = 30, alpha = 0.7, fill = "steelblue") +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
    labs(title = glue("{pth}: Soft Probability Distribution"),
         x = "Soft Probability", y = "Count") +
    theme_minimal()
  
  dist_plots[[pth]] <- p
}

combined_dist <- wrap_plots(dist_plots, ncol = 2)
ggsave("results/bayesian_regression_distributions.pdf", combined_dist, 
       width = 12, height = 8, dpi = 300)

# Seroprevalence comparison plot
seroprev_data <- map_dfr(names(results), function(pth) {
  tibble(
    Pathogen = pth,
    Seroprevalence = sum(results[[pth]][[glue("{tolower(pth)}_sero_hard_50")]]) / nrow(results[[pth]]),
    Method = "Bayesian Regression"
  )
})

seroprev_plot <- seroprev_data %>%
  ggplot(aes(x = Pathogen, y = Seroprevalence, fill = Method)) +
  geom_col(alpha = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Seroprevalence Estimates by Pathogen",
       y = "Seroprevalence", x = "Pathogen") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/bayesian_regression_seroprevalence.pdf", seroprev_plot, 
       width = 8, height = 6, dpi = 300)

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

cat("\nDiagnostic plots saved to results/\n")
cat("Seroprevalence summary saved to results/bayesian_regression_seroprevalence_summary.tsv\n")
cat("Weight files written.\nScript 06 complete.\n") 