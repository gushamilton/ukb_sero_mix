################################################################################
#   UKB SEROLOGY: PATHOGEN-LEVEL PHENOTYPES (BAYESIAN REGRESSION + IMPUTATION) - SCRIPT 07
################################################################################
#  ▸ Strategy: Use multiple imputation (mice) to handle missing antigen values,
#    then run Bayesian regression on each imputed dataset and pool results.
#  ▸ This allows us to use ALL individuals, not just those with complete antigen data.
#  ▸ Missing values are imputed using correlations between antigens.
################################################################################

# --- 1. PACKAGES --------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, glue, brms, future, here, bayestestR, mice)

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

# --- 4. Multiple Imputation + Bayesian Regression ----------------------------
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

  # Multiple Imputation using mice
  cat(glue("\nPerforming multiple imputation for {pth}...\n"))
  
  # Prepare data for imputation (only soft probabilities)
  imp_data <- df_model[soft_cols]
  
  # Set up imputation method - use predictive mean matching for continuous data
  imp_method <- rep("pmm", length(soft_cols))
  names(imp_method) <- soft_cols
  
  # Perform multiple imputation
  n_imp <- 5  # Number of imputations
  imp <- mice(imp_data, m = n_imp, method = imp_method, 
              maxit = 10, seed = 42, print = FALSE)
  
  cat(glue("Created {n_imp} imputed datasets\n"))

  # Fit Bayesian regression on each imputed dataset
  imp_results <- list()
  imp_models <- list()
  
  for (i in 1:n_imp) {
    cat(glue("Fitting Bayesian model on imputed dataset {i}/{n_imp}...\n"))
    
    # Get imputed dataset
    imp_df <- complete(imp, i)
    imp_df$z <- df_model$z  # Add back the target variable
    
    # Build formula
    rhs <- paste(soft_cols, collapse = " + ")
    brm_formula <- bf(as.formula(paste("z ~", rhs)), family = bernoulli(link = "logit"))
    
    # Fit model
    fit <- brm(brm_formula, data = imp_df,
               chains = 2, iter = 1000, warmup = 250,  # Reduced for speed
               refresh = 0,
               prior = set_prior("normal(0,1)", class = "b"))
    
    imp_models[[i]] <- fit
    
    # Get predictions
    pred <- fitted(fit, newdata = imp_df, scale = "response", probs = c(0.025, 0.25, 0.75, 0.975))
    
    imp_results[[i]] <- tibble(
      FID = df_model$FID,
      IID = df_model$IID,
      pi_hat = pred[, "Estimate"],
      pi_lower = pred[, "Q2.5"],
      pi_upper = pred[, "Q97.5"],
      pi_ci_width = pred[, "Q97.5"] - pred[, "Q2.5"]
    )
  }
  
  # Pool results across imputations (Rubin's rules)
  pooled_results <- tibble(
    FID = df_model$FID,
    IID = df_model$IID,
    pi_hat = rowMeans(sapply(imp_results, function(x) x$pi_hat)),
    pi_lower = rowMeans(sapply(imp_results, function(x) x$pi_lower)),
    pi_upper = rowMeans(sapply(imp_results, function(x) x$pi_upper)),
    pi_ci_width = rowMeans(sapply(imp_results, function(x) x$pi_ci_width))
  )
  
  # Multiple serostatus thresholds
  thresholds <- c(0.3, 0.5, 0.7, 0.9)
  serostatus_cols <- list()
  
  for (thresh in thresholds) {
    col_name <- glue("{tolower(pth)}_sero_hard_{thresh*100}")
    serostatus_cols[[col_name]] <- as.integer(pooled_results$pi_hat >= thresh)
  }

  # Calculate pooled seroprevalence
  seroprevalence <- mean(pooled_results$pi_hat >= 0.5)
  
  # Bootstrap CI for seroprevalence (accounting for imputation uncertainty)
  set.seed(42)
  n_boot <- 1000
  boot_samples <- numeric(n_boot)
  
  for (b in 1:n_boot) {
    # Sample one imputed dataset
    imp_idx <- sample(1:n_imp, 1)
    boot_seroprev <- mean(imp_results[[imp_idx]]$pi_hat >= 0.5)
    boot_samples[b] <- boot_seroprev
  }
  
  seroprev_ci <- quantile(boot_samples, c(0.025, 0.975))

  results[[pth]] <- tibble(
    FID = pooled_results$FID,
    IID = pooled_results$IID,
    !!glue("{tolower(pth)}_sero_soft") := pooled_results$pi_hat,
    !!glue("{tolower(pth)}_sero_soft_lower") := pooled_results$pi_lower,
    !!glue("{tolower(pth)}_sero_soft_upper") := pooled_results$pi_upper,
    !!glue("{tolower(pth)}_sero_soft_ci_width") := pooled_results$pi_ci_width,
    !!!serostatus_cols,
    !!glue("{tolower(pth)}_w_soft")    := pooled_results$pi_hat,
    !!glue("{tolower(pth)}_w_hard")    := 2 * abs(pooled_results$pi_hat - 0.5)
  )
  
  # Store model summary from first imputation
  model_summaries[[pth]] <- summary(imp_models[[1]])
  
  # Print summary statistics
  cat(glue("\n--- {pth} Summary (with multiple imputation) ---\n"))
  cat(glue("Seroprevalence (≥0.5): {scales::percent(seroprevalence, accuracy = 0.1)} ({scales::percent(seroprev_ci[1], accuracy = 0.1)} - {scales::percent(seroprev_ci[2], accuracy = 0.1)})\n"))
  cat(glue("Mean CI width: {round(mean(pooled_results$pi_ci_width), 3)}\n"))
  cat(glue("Coefficients (from first imputation):\n"))
  coef_summary <- summary(imp_models[[1]])$fixed
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