################################################################################
#   GWAS POWER SIMULATION: REALISTIC MIXTURE MODEL COMPARISON
################################################################################
#  ∙ This script simulates GWAS data to provide a REALISTIC comparison of the
#    statistical power AND BIAS of different analysis methods.
#
#  ∙ It simulates data with a known "ground truth" serostatus, then fits a
#    mixture model to the continuous data to derive "soft" probabilities,
#    mimicking a real-world analysis workflow.
#
#  ∙ Methods Compared:
#    - Gold Standard: Uses the true, latent serostatus.
#    - Hard Cutoff: The standard approach using a fixed MFI threshold.
#    - Mixture Model: A realistic approach using posterior probabilities
#      estimated by fitting a skew-t mixture model in each simulation.
################################################################################

# --- 1. Packages & Setup ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tibble, purrr, sn, mixsmsn, furrr, future)

# --- Set up parallel processing ---
# Use multicore on Mac/Linux for efficiency. For Windows, use plan(multisession).
plan(multicore)

# --- 2. Simulation Parameters ---
sim_params <- list(
  cmv_pp150 = list(pii = c(0.45, 0.55), mu = c(3.5, 6.5), sigma2 = c(1.0, 1.2), shape = c(1.5, -2), hard_cutoff = 100),
  hsv1 = list(pii = c(0.3, 0.7), mu = c(3.0, 7.0), sigma2 = c(1.5, 1.5), shape = c(1.5, -2.5), hard_cutoff = 150)
)

# --- 3. Core Simulation & Analysis Functions ---

#' Mixture model fitting function (from serology_power_analysis.R)
fit_skew_mix <- function(y) {
  tryCatch(
    mixsmsn::smsn.mix(y, g = 2, family = "Skew.t", nu = 4, group = FALSE, calc.im = FALSE, obs.prob = TRUE),
    error = function(e) { NULL }
  )
}

#' Generate a simulated dataset
generate_sim_data <- function(n_samples, params, beta, gamma, maf) {
  # Simulate two independent SNPs
  snp_beta <- rbinom(n_samples, 2, maf)
  snp_gamma <- rbinom(n_samples, 2, maf)
  
  # Simulate "true" serostatus based on snp_beta's effect
  prob_positive <- params$pii[2] + (snp_beta * beta)
  S_true <- rbinom(n_samples, 1, prob_positive)
  
  # Simulate MFI values, adding the effect of snp_gamma to the positive mean
  n_pos <- sum(S_true == 1)
  n_neg <- sum(S_true == 0)
  
  # Create MFI vectors, one for each true status
  mfi_vec <- numeric(n_samples)
  
  # For true positives, the mean is shifted by the gamma effect
  mfi_vec[S_true == 1] <- exp(rsn(n_pos, xi = params$mu[2] + (snp_gamma[S_true == 1] * gamma), omega = sqrt(params$sigma2[2]), alpha = params$shape[2]))
  mfi_vec[S_true == 0] <- exp(rsn(n_neg, xi = params$mu[1], omega = sqrt(params$sigma2[1]), alpha = params$shape[1]))
  
  tibble(snp_beta, snp_gamma, S_true, mfi = mfi_vec) %>%
    mutate(S_hard = ifelse(mfi >= params$hard_cutoff, 1, 0))
}

#' Run power and bias analysis for a single simulated dataset
run_one_analysis <- function(sim_data) {
  
  # --- NEW: Fit mixture model to get realistic soft probabilities ---
  fit <- fit_skew_mix(log(sim_data$mfi[sim_data$mfi > 0]))
  p_soft_estimated <- NA
  if (!is.null(fit)) {
    pos_comp <- which.max(fit$mu)
    p_soft_estimated <- fit$obs.prob[, pos_comp]
  }
  sim_data$p_soft_estimated <- p_soft_estimated

  # --- Analysis for BETA effect (on seropositivity) ---
  fit_beta_true <- glm(S_true ~ snp_beta, data = sim_data, family = binomial())
  fit_beta_hard <- glm(S_hard ~ snp_beta, data = sim_data, family = binomial())
  fit_beta_soft <- lm(p_soft_estimated ~ snp_beta, data = sim_data) # Use estimated probs

  # --- Analysis for GAMMA effect (on IgG Level) ---
  fit_gamma_true <- lm(log(mfi) ~ snp_gamma, data = filter(sim_data, S_true == 1))
  fit_gamma_hard <- lm(log(mfi) ~ snp_gamma, data = filter(sim_data, S_hard == 1))
  fit_gamma_soft <- lm(log(mfi) ~ snp_gamma, data = sim_data, weights = p_soft_estimated) # Use estimated probs

  # Safely extract coefficient estimates and p-values
  extract_coef <- function(model, coef_name) {
    if (is.null(model) || any(is.na(coef(model)))) return(c(NA, NA))
    s <- summary(model)$coefficients
    if (coef_name %in% rownames(s)) {
      p_val_col <- if ("Pr(>|z|)" %in% colnames(s)) "Pr(>|z|)" else "Pr(>|t|)"
      return(s[coef_name, c("Estimate", p_val_col)])
    }
    return(c(NA, NA))
  }
  
  beta_res <- rbind(extract_coef(fit_beta_true, "snp_beta"), extract_coef(fit_beta_hard, "snp_beta"), extract_coef(fit_beta_soft, "snp_beta"))
  gamma_res <- rbind(extract_coef(fit_gamma_true, "snp_gamma"), extract_coef(fit_gamma_hard, "snp_gamma"), extract_coef(fit_gamma_soft, "snp_gamma"))
  
  tibble(
    p_beta_true = beta_res[1,2], p_beta_hard = beta_res[2,2], p_beta_soft = beta_res[3,2],
    p_gamma_true = gamma_res[1,2], p_gamma_hard = gamma_res[2,2], p_gamma_soft = gamma_res[3,2],
    est_beta_true = beta_res[1,1], est_beta_hard = beta_res[2,1], est_beta_soft = beta_res[3,1],
    est_gamma_true = gamma_res[1,1], est_gamma_hard = gamma_res[2,1], est_gamma_soft = gamma_res[3,1]
  )
}

# --- 4. Main Simulation Driver ---
run_power_simulation <- function(n_sims, n_samples, antibody_name, beta, gamma, maf) {
  cat("\n--- Running REALISTIC Power & Bias Simulation for:", antibody_name, "---\n")
  cat("Sims:", n_sims, "| Samples:", n_samples, "| Beta:", beta, "| Gamma:", gamma, "| MAF:", maf, "\n")
  
  params <- sim_params[[antibody_name]]
  
  results <- future_map_dfr(1:n_sims, ~{
    # The cat statement from the original loop was removed as it's not safe in parallel.
    # A progress bar is enabled in future_map_dfr instead.
    sim_data <- generate_sim_data(n_samples, params, beta, gamma, maf)
    run_one_analysis(sim_data)
  }, .options = furrr_options(seed = TRUE), .progress = TRUE)
  
  alpha <- 0.05
  
  power_beta <- results %>% summarise(Analysis = "Seropositivity (Beta)", Power_Gold_Standard = mean(p_beta_true < alpha, na.rm=T), Power_Hard_Cutoff = mean(p_beta_hard < alpha, na.rm=T), Power_Mixture_Model = mean(p_beta_soft < alpha, na.rm=T))
  power_gamma <- results %>% summarise(Analysis = "IgG Level (Gamma)", Power_Gold_Standard = mean(p_gamma_true < alpha, na.rm=T), Power_Hard_Cutoff = mean(p_gamma_hard < alpha, na.rm=T), Power_Mixture_Model = mean(p_gamma_soft < alpha, na.rm=T))
  power_summary <- bind_rows(power_beta, power_gamma) %>% mutate(Antibody = antibody_name, .before = 1)
  
  bias_beta <- results %>% summarise(Analysis = "Seropositivity (Beta)", True_Effect = beta, Mean_Est_Gold_Standard = mean(est_beta_true, na.rm=T), Mean_Est_Hard_Cutoff = mean(est_beta_hard, na.rm=T), `Mean_Est_Mixture_Model (Prob Scale)` = mean(est_beta_soft, na.rm=T))
  bias_gamma <- results %>% summarise(Analysis = "IgG Level (Gamma)", True_Effect = gamma, Mean_Est_Gold_Standard = mean(est_gamma_true, na.rm=T), Mean_Est_Hard_Cutoff = mean(est_gamma_hard, na.rm=T), Mean_Est_Mixture_Model = mean(est_gamma_soft, na.rm=T))
  bias_summary <- bind_rows(bias_beta, bias_gamma) %>% mutate(Antibody = antibody_name, .before = 1)
  
  return(list(power = power_summary, bias = bias_summary))
}

# --- 5. Execute Simulations ---
N_SIMULATIONS <- 20
N_SAMPLES <- 4000
BETA_EFFECT <- 0.05 # Effect on serostatus risk (prob scale)
GAMMA_EFFECT <- 0.06 # Effect on log-MFI value
SNP_MAF <- 0.2

# Run sims and extract results
results_beta_cmv <- run_power_simulation(N_SIMULATIONS, N_SAMPLES, "cmv_pp150", BETA_EFFECT, 0, SNP_MAF)
results_beta_hsv1 <- run_power_simulation(N_SIMULATIONS, N_SAMPLES, "hsv1", BETA_EFFECT, 0, SNP_MAF)
results_gamma_cmv <- run_power_simulation(N_SIMULATIONS, N_SAMPLES, "cmv_pp150", 0, GAMMA_EFFECT, SNP_MAF)
results_gamma_hsv1 <- run_power_simulation(N_SIMULATIONS, N_SAMPLES, "hsv1", 0, GAMMA_EFFECT, SNP_MAF)

# --- 6. Report Final Results ---
cat("\n\n--- POWER Results: SNP effect on SEROPOSITIVITY (Beta) ---\n")
print(bind_rows(results_beta_cmv$power, results_beta_hsv1$power))

cat("\n\n--- BIAS Results: SNP effect on SEROPOSITIVITY (Beta) ---\n")
print(bind_rows(results_beta_cmv$bias, results_beta_hsv1$bias))

cat("\n\n--- POWER Results: SNP effect on IGG LEVEL (Gamma) ---\n")
print(bind_rows(results_gamma_cmv$power, results_gamma_hsv1$power))

cat("\n\n--- BIAS Results: SNP effect on IGG LEVEL (Gamma) ---\n")
print(bind_rows(results_gamma_cmv$bias, results_gamma_hsv1$bias))

cat("\n* Power: The percentage of simulations where the SNP's effect was detected (p < 0.05).\n")
cat("* Bias: Compares the mean estimated effect size to the true effect size.\n") 

