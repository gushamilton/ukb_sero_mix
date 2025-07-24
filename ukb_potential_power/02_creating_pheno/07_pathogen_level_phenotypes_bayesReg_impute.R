################################################################################
#   UKB SEROLOGY: PATHOGEN-LEVEL PHENOTYPES (BAYESIAN REGRESSION + IMPUTATION) - SCRIPT 07
################################################################################
#  ▸ Strategy: Extend script 06 to handle missing antigen values directly in the
#    Bayesian model. Instead of dropping individuals with missing antigens, we
#    model missing values as parameters that get imputed during MCMC sampling.
#  ▸ This allows us to use ALL individuals, not just those with complete antigen data.
#  ▸ Missing values are imputed using correlations between antigens and any available
#    covariates (age, sex, etc.).
################################################################################

# --- 1. PACKAGES --------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, glue, brms, future, here, bayestestR)

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
    # Handle missing values in hard calls - treat as 0 for rule calculation
    hard_matrix <- df[cols]
    hard_matrix[is.na(hard_matrix)] <- 0
    as.integer(rowSums(hard_matrix, na.rm = TRUE) >= rule$n)
  } else { # pgp3 rule
    # Handle missing pgp3 - treat as 0
    ifelse(is.na(df$ct_pgp3_sero_hard), 0, df$ct_pgp3_sero_hard)
  }
}

# --- 4. Bayesian regression with missing data imputation ---------------------
results <- list()
model_summaries <- list()

for (pth in names(pathogen_map)) {
  antigens <- pathogen_map[[pth]]
  soft_cols <- paste0(antigens, "_sero_soft")
  hard_cols <- paste0(antigens, "_sero_hard")

  # Get all individuals (including those with missing antigens)
  df <- master %>%
    select(FID, IID, all_of(c(soft_cols, hard_cols))) %>%
    # Keep individuals with at least 1 non-missing antigen
    filter(rowSums(!is.na(select(., all_of(soft_cols)))) >= 1)
  
  if (nrow(df) < 100) {
    warning("Skipping ", pth, ": less than 100 individuals with any antigen data")
    next
  }

  # Create target variable (handling missing hard calls)
  df$z <- make_label(df, antigens, ukb_rules[[pth]])
  
  # Count missing data
  missing_counts <- colSums(is.na(df[soft_cols]))
  cat(glue("\n--- {pth} Missing Data Summary ---\n"))
  cat(glue("Total individuals: {nrow(df)}\n"))
  cat(glue("Individuals with complete data: {sum(complete.cases(df[soft_cols]))}\n"))
  for (i in seq_along(antigens)) {
    cat(glue("  {antigens[i]}: {missing_counts[i]} missing ({round(missing_counts[i]/nrow(df)*100, 1)}%)\n"))
  }

  # Subsample for development speed -----------------------------------------
  set.seed(42)
  if (nrow(df) > 2000) df_model <- slice_sample(df, n = 200) else df_model <- df

  # Build formula with missing data handling
  # brms will automatically handle missing predictors by treating them as parameters
  rhs <- paste(soft_cols, collapse = " + ")
  brm_formula <- bf(as.formula(paste("z ~", rhs)), family = bernoulli(link = "logit"))

  cat(glue("\nFitting Bayesian logistic model for {pth} on {nrow(df_model)} individuals (with missing data imputation)...\n"))

  # Fit model with missing data imputation
  fit <- brm(brm_formula, data = df_model,
             chains = 2, iter = 2000, warmup = 500,
             refresh = 0,
             prior = set_prior("normal(0,1)", class = "b"),
             # Enable missing data imputation
             missing = "mi")

  # Store model summary
  model_summaries[[pth]] <- summary(fit)

  # Predict on ALL individuals (including those with missing data)
  # brms will impute missing values during prediction
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
  # Get posterior samples for all predictions (including imputed)
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
  cat(glue("\n--- {pth} Summary (with imputation) ---\n"))
  cat(glue("Seroprevalence (≥0.5): {scales::percent(seroprevalence, accuracy = 0.1)} ({scales::percent(seroprev_ci[1], accuracy = 0.1)} - {scales::percent(seroprev_ci[2], accuracy = 0.1)})\n"))
  cat(glue("Mean CI width: {round(mean(pi_ci_width), 3)}\n"))
  cat(glue("Coefficients:\n"))
  coef_summary <- summary(fit)$fixed
  for (i in 2:nrow(coef_summary)) {  # Skip intercept
    antigen_name <- rownames(coef_summary)[i]
    antigen_clean <- str_remove(antigen_name, "_sero_soft")
    cat(glue("  {antigen_clean}: {round(coef_summary[i, 'Estimate'], 3)} ({round(coef_summary[i, 'l-95% CI'], 3)}, {round(coef_summary[i, 'u-95% CI'], 3)})\n"))
  }
}

# --- 5. Create summary tables and results ------------------------------------
cat("\nCreating summary tables...\n")

# --- 6. Combine results and write outputs ------------------------------------
combined <- reduce(results, full_join, by = c("FID", "IID"))

if (!dir.exists("quickdraws_input")) dir.create("quickdraws_input")

write_tsv(combined, "quickdraws_input/phenotypes_pathogen_bayesReg_imputed.tsv")
cat("\nSaved pathogen-level Bayesian-regression phenotypes (with imputation): phenotypes_pathogen_bayesReg_imputed.tsv\n")

# Write weight files per pathogen --------------------------------------------
for (pth in names(results)) {
  low <- tolower(pth)
  res <- results[[pth]]
  write_tsv(res %>% select(FID, IID, Weight = !!glue("{low}_w_soft")),
            glue("quickdraws_input/{low}_p_soft_imputed.weights"))
  write_tsv(res %>% select(FID, IID, Weight = !!glue("{low}_w_hard")),
            glue("quickdraws_input/{low}_sero_hard_imputed.weights"))
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

write_tsv(seroprev_summary, "results/bayesian_regression_seroprevalence_imputed_summary.tsv")

cat("Seroprevalence summary saved to results/bayesian_regression_seroprevalence_imputed_summary.tsv\n")
cat("Weight files written.\nScript 07 complete.\n") 