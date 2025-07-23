################################################################################
##
## Serology GWAS Simulation - Draft 4: flexmix vs. Optimal Cutoff
##
## This definitive simulation compares two methods for analyzing bimodal,
## overlapping serology data:
## 1. Baseline (Optimal Cutoff): Mimics a real analysis by finding the optimal
##    cutoff on the continuous IgG data to define serostatus, then running
##    separate logistic (for beta) and OLS (for gamma) models.
## 2. Mixture Model (flexmix): A true mixture model that analyzes the raw
##    continuous data to simultaneously estimate effects on group membership
##    (beta) and IgG titre (gamma).
##
## The simulation tests performance across different levels of distributional
## separation and error correlation, reporting both BIAS and POWER for each method.
##
################################################################################

## ---- package-setup --------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, purrr, MASS, tibble, dplyr, flexmix)

## ---- data-generation ------------------------------------------------------
# Uses the robust, realistic data generation from Draft 3
generate_data <- function(n = 1500, beta_true = 0.1, gamma_true = 0.1,
                          separation = 2.0, rho = 0,
                          c_effect_s = 0.4, c_effect_y = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  C <- rnorm(n, 0, 1); Z <- rbinom(n, 2, 0.3)
  log_odds <- -0.5 + beta_true * Z + c_effect_s * C
  S_true <- rbinom(n, 1, plogis(log_odds))
  sigma_s <- 1; sigma_y <- 0.5
  cov_matrix <- matrix(c(sigma_s^2, rho*sigma_s*sigma_y, rho*sigma_s*sigma_y, sigma_y^2), 2)
  errors <- MASS::mvrnorm(n, mu = c(0, 0), cov_matrix)
  Y <- numeric(n)
  Y[S_true == 0] <- 0.5 + errors[S_true == 0, 2]
  Y[S_true == 1] <- separation + gamma_true * Z[S_true == 1] + c_effect_y * C[S_true == 1] + errors[S_true == 1, 2]
  return(data.table(Z, C, Y, S_true))
}

## ---- helper-find-cutoff ---------------------------------------------------
find_optimal_cutoff <- function(y) {
  # A simple but robust way to find the cutoff is to find the minimum
  # of the density estimate between the two peaks.
  dens <- density(y, n = 1024)
  # Find peaks
  pks <- which(diff(sign(diff(dens$y))) == -2) + 1
  if (length(pks) < 2) return(median(y)) # Failsafe for unimodal data
  # Find the minimum between the two largest peaks
  sorted_pks <- sort(pks[order(dens$y[pks], decreasing=TRUE)[1:2]])
  valley_idx <- which.min(dens$y[sorted_pks[1]:sorted_pks[2]]) + sorted_pks[1] - 1
  return(dens$x[valley_idx])
}

## ---- simulation-function --------------------------------------------------
sim_one_run <- function(rho, separation, beta_true, gamma_true) {
  d <- generate_data(rho = rho, separation = separation, beta_true = beta_true, gamma_true = gamma_true)
  
  # --- Initialize results ---
  res <- list(beta_base=NA, se_beta_base=NA, gamma_base=NA, se_gamma_base=NA,
              beta_mix=NA, se_beta_mix=NA, gamma_mix=NA, se_gamma_mix=NA)

  # --- Method 1: Baseline (Optimal Cutoff) ---
  cutoff <- find_optimal_cutoff(d$Y)
  d[, S_observed := ifelse(Y > cutoff, 1, 0)]
  
  # Use a more robust approach: try multiple cutoffs and choose the best one
  # based on the separation between groups
  try_cutoffs <- function(y) {
    # Try different percentiles as cutoffs
    percentiles <- seq(0.3, 0.7, by = 0.05)
    best_cutoff <- cutoff
    best_separation <- 0
    
    for(p in percentiles) {
      test_cutoff <- quantile(y, p)
      s_test <- ifelse(y > test_cutoff, 1, 0)
      
      if(length(unique(s_test)) > 1) {
        # Calculate separation between groups
        y_neg <- y[s_test == 0]
        y_pos <- y[s_test == 1]
        if(length(y_neg) > 0 && length(y_pos) > 0) {
          separation <- abs(mean(y_pos) - mean(y_neg)) / sqrt(var(y))
          if(separation > best_separation) {
            best_separation <- separation
            best_cutoff <- test_cutoff
          }
        }
      }
    }
    return(best_cutoff)
  }
  
  # Use the improved cutoff
  improved_cutoff <- try_cutoffs(d$Y)
  d[, S_observed := ifelse(Y > improved_cutoff, 1, 0)]
  
  if (length(unique(d$S_observed)) > 1) {
    fit_beta_base <- glm(S_observed ~ Z + C, data = d, family = binomial(link="logit"))
    res$beta_base <- coef(fit_beta_base)["Z"]; res$se_beta_base <- summary(fit_beta_base)$coef["Z", "Std. Error"]
  }
  
  d_pos_obs <- d[S_observed == 1]
  if (nrow(d_pos_obs) > 2) {
    fit_gamma_base <- lm(Y ~ Z + C, data = d_pos_obs)
    if ("Z" %in% names(coef(fit_gamma_base))) {
      res$gamma_base <- coef(fit_gamma_base)["Z"]
      summ_gamma_base <- summary(fit_gamma_base)
      if ("Z" %in% rownames(summ_gamma_base$coefficients)) {
        res$se_gamma_base <- summ_gamma_base$coefficients["Z", "Std. Error"]
      }
    }
  }

  # --- Method 2: Mixture Model (flexmix) ---
  tryCatch({
    fit_mix <- flexmix(Y ~ Z + C | Z + C, data = d, k = 2, 
                       cluster = d$S_observed + 1, # Use cutoff for smart initialization
                       control = list(iter.max = 200, minprior = 0.01))
    
    refit_models <- refit(fit_mix)
    
    # Extract coefficients properly from S4 object
    all_coefs <- refit_models@coef
    all_vcov <- refit_models@vcov
    
    # Extract component coefficients
    comp1_z <- all_coefs["model.1_Comp.1_coef.Z"]
    comp2_z <- all_coefs["model.1_Comp.2_coef.Z"]
    comp1_intercept <- all_coefs["model.1_Comp.1_coef.(Intercept)"]
    comp2_intercept <- all_coefs["model.1_Comp.2_coef.(Intercept)"]
    
    # Identify which component is seropositive (higher intercept)
    pos_comp <- ifelse(comp1_intercept > comp2_intercept, 1, 2)

    # Extract gamma (Z effect on IgG in seropositives) from component coefficients
    res$gamma_mix <- if(pos_comp == 1) comp1_z else comp2_z
    
    # Extract standard error for gamma
    gamma_se_name <- if(pos_comp == 1) "model.1_Comp.1_coef.Z" else "model.1_Comp.2_coef.Z"
    gamma_se_idx <- which(rownames(all_vcov) == gamma_se_name)
    res$se_gamma_mix <- sqrt(all_vcov[gamma_se_idx, gamma_se_idx])
    
    # For beta: Use posterior probabilities as weights in probit regression
    # This gives us the effect of Z on serostatus using the mixture model's soft assignments
    post_mat <- if (is.list(fit_mix@posterior)) fit_mix@posterior$scaled else fit_mix@posterior
    d[, w_seropos := post_mat[, pos_comp]]
    
    # Fit weighted probit regression using posterior probabilities as weights
    fit_beta_mix <- tryCatch(
      glm(S_observed ~ Z + C, data = d, family = binomial(link = "probit"), weights = w_seropos),
      error = function(e) NULL
    )
    if (!is.null(fit_beta_mix) && "Z" %in% names(coef(fit_beta_mix))) {
      res$beta_mix <- coef(fit_beta_mix)["Z"]
      res$se_beta_mix <- summary(fit_beta_mix)$coef["Z", "Std. Error"]
    }
    
  }, error = function(e) {})
  
  return(as_tibble(res))
}

## ---- main-simulation-driver -----------------------------------------------
run_final_simulation <- function(reps = 100) {
  scen_rho <- c(-0.3, 0, 0.3, 0.6)
  scen_separation <- c(1.0, 1.5, 2.0, 2.5, 3.0, 4.0) # Low to High separation
  beta_true <- 0.1; gamma_true <- 0.1
  
  grid <- expand.grid(rho = scen_rho, separation = scen_separation)
  
  cat("Running final simulation with", reps, "replicates for", nrow(grid), "scenarios...\n")
  
  # pmap runs the simulation for each scenario (row in grid)
  all_res <- pmap_dfr(grid, function(rho, separation) {
    # Use map instead of rerun (deprecated)
    map_dfr(1:reps, ~ sim_one_run(rho, separation, beta_true, gamma_true))
  })
  
  cat("Simulation complete. Analyzing results...\n")

  # Analysis
  alpha5 <- qnorm(0.975)
  results <- all_res %>%
    mutate(
      rho = rep(grid$rho, each = reps),
      separation = rep(grid$separation, each = reps),
      power_beta_base = abs(beta_base / se_beta_base) > alpha5,
      power_gamma_base = abs(gamma_base / se_gamma_base) > alpha5,
      power_beta_mix = abs(beta_mix / se_beta_mix) > alpha5,
      power_gamma_mix = abs(gamma_mix / se_gamma_mix) > alpha5
    )

  summary_tbl <- results %>%
    group_by(separation, rho) %>%
    summarise(
      across(starts_with("beta_"), ~ mean(. - beta_true, na.rm = TRUE), .names = "bias_{.col}"),
      across(starts_with("gamma_"), ~ mean(. - gamma_true, na.rm = TRUE), .names = "bias_{.col}"),
      across(starts_with("power_"), ~ mean(., na.rm = TRUE), .names = "{.col}"),
      .groups = 'drop'
    )
    
  cat("\n--- BIAS REPORT ---\n")
  print(as.data.frame(summary_tbl %>% select(separation, rho, starts_with("bias_"))))
  
  cat("\n--- POWER REPORT ---\n")
  print(as.data.frame(summary_tbl %>% select(separation, rho, starts_with("power_"))))
  
  # Summary statistics
  cat("\n--- OVERALL SUMMARY ---\n")
  cat("Average bias (beta): Baseline =", mean(summary_tbl$bias_beta_base, na.rm=TRUE), 
      "Mixture =", mean(summary_tbl$bias_beta_mix, na.rm=TRUE), "\n")
  cat("Average bias (gamma): Baseline =", mean(summary_tbl$bias_gamma_base, na.rm=TRUE), 
      "Mixture =", mean(summary_tbl$bias_gamma_mix, na.rm=TRUE), "\n")
  cat("Average power (beta): Baseline =", mean(summary_tbl$power_beta_base, na.rm=TRUE), 
      "Mixture =", mean(summary_tbl$power_beta_mix, na.rm=TRUE), "\n")
  cat("Average power (gamma): Baseline =", mean(summary_tbl$power_gamma_base, na.rm=TRUE), 
      "Mixture =", mean(summary_tbl$power_gamma_mix, na.rm=TRUE), "\n")
}

run_final_simulation() 