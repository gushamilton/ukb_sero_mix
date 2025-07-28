################################################################################
#   GWAS POWER SIMULATION (V5): FOCUSED SCENARIOS
################################################################################
#  ∙ This script's SOLE responsibility is to run the comprehensive simulation
#    for three specific scenarios and save the raw, unsummarized results.
#
#  ∙ SCENARIOS:
#    1. High Separation: Clear distinction between seropositive/negative groups.
#    2. Low Separation: Significant overlap between groups.
#    3. Truncated: The "Low Separation" data with left-truncation to mimic
#       a limit of detection.
################################################################################

# --- 1. Packages & Setup ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, tibble, purrr, sn, mixsmsn, furrr, future, readr)

# --- Set up parallel processing ---
plan(multicore)

# --- 2. Simulation Parameters ---
base_params <- list(
    high_separation = list(pii = c(0.4, 0.6), mu = c(4.0, 8.0), sigma2 = c(1.5, 1.2), shape = c(1.5, -3)),
    low_separation  = list(pii = c(0.4, 0.6), mu = c(4.0, 6.2), sigma2 = c(1.5, 1.2), shape = c(1.5, -3))
)

find_optimal_cutoff <- function(params) {
    set.seed(12345)
    n_test <- 100000
    component <- rbinom(n_test, 1, params$pii[2])
    mfi_log <- ifelse(component == 1,
        rsn(n_test, xi = params$mu[2], omega = sqrt(params$sigma2[2]), alpha = params$shape[2]),
        rsn(n_test, xi = params$mu[1], omega = sqrt(params$sigma2[1]), alpha = params$shape[1])
    )
    mfi <- exp(mfi_log)
    cutoffs <- quantile(mfi, seq(0.05, 0.95, 0.01))
    accuracies <- sapply(cutoffs, function(c) mean(ifelse(mfi >= c, 1, 0) == component))
    optimal_cutoff <- cutoffs[which.max(accuracies)]
    cat("Optimal cutoff for antibody:", round(optimal_cutoff, 2), "with accuracy:", round(max(accuracies), 3), "\n")
    return(optimal_cutoff)
}

cat("--- Calculating Optimal Hard Cutoffs ---\n")
sim_params <- map(base_params, ~ c(.x, hard_cutoff = find_optimal_cutoff(.x)))

# --- 3. Core Simulation & Analysis Functions ---
generate_sim_data <- function(n_samples, params, beta, gamma, maf) {
    snp <- rbinom(n_samples, 2, maf)
    base_log_odds <- qlogis(params$pii[2])
    prob_positive <- plogis(base_log_odds + (snp * beta))
    S_true <- rbinom(n_samples, 1, prob_positive)
    mfi_vec <- numeric(n_samples)
    n_pos <- sum(S_true == 1); n_neg <- sum(S_true == 0)
    if (n_pos > 0) mfi_vec[S_true == 1] <- exp(rsn(n_pos, xi = params$mu[2] + (snp[S_true == 1] * gamma), omega = sqrt(params$sigma2[2]), alpha = params$shape[2]))
    if (n_neg > 0) mfi_vec[S_true == 0] <- exp(rsn(n_neg, xi = params$mu[1], omega = sqrt(params$sigma2[1]), alpha = params$shape[1]))
    tibble(snp, S_true, mfi = mfi_vec) %>% mutate(S_hard = ifelse(mfi >= params$hard_cutoff, 1, 0))
}

run_one_analysis <- function(sim_data) {
    fit_skew_mix <- function(y) tryCatch(mixsmsn::smsn.mix(y, g = 2, family = "Skew.t", nu = 4, group = FALSE, calc.im = FALSE, obs.prob = TRUE), error = function(e) NULL)
    fit <- fit_skew_mix(log(sim_data$mfi[sim_data$mfi > 0 & is.finite(log(sim_data$mfi))]))
    p_soft_estimated <- if (!is.null(fit)) fit$obs.prob[, which.max(fit$mu)] else NA
    sim_data$p_soft_estimated <- p_soft_estimated
    fit_beta_true <- glm(S_true ~ snp, data = sim_data, family = binomial())
    fit_beta_hard <- glm(S_hard ~ snp, data = sim_data, family = binomial())
    fit_beta_soft <- if (any(!is.na(p_soft_estimated))) glm(p_soft_estimated ~ snp, data = sim_data, family = quasibinomial()) else NULL
    fit_gamma_true <- if (sum(sim_data$S_true == 1) > 2) lm(log(mfi) ~ snp, data = filter(sim_data, S_true == 1)) else NULL
    fit_gamma_hard <- if (sum(sim_data$S_hard == 1) > 2) lm(log(mfi) ~ snp, data = filter(sim_data, S_hard == 1)) else NULL
    fit_gamma_soft <- if (any(!is.na(p_soft_estimated))) lm(log(mfi) ~ snp, data = sim_data, weights = p_soft_estimated) else NULL
    extract_coef <- function(model, coef_name = "snp") {
        if (is.null(model) || any(is.na(coef(model)))) return(c(NA, NA, NA))
        s <- summary(model)$coefficients
        if (coef_name %in% rownames(s)) {
            p_val_col <- if ("Pr(>|z|)" %in% colnames(s)) "Pr(>|z|)" else "Pr(>|t|)"
            return(s[coef_name, c("Estimate", p_val_col, "Std. Error")])
        }
        return(c(NA, NA, NA))
    }
    beta_res <- rbind(extract_coef(fit_beta_true), extract_coef(fit_beta_hard), extract_coef(fit_beta_soft))
    gamma_res <- rbind(extract_coef(fit_gamma_true), extract_coef(fit_gamma_hard), extract_coef(fit_gamma_soft))
    tibble(
        p_beta_true = beta_res[1,2], p_beta_hard = beta_res[2,2], p_beta_soft = beta_res[3,2],
        p_gamma_true = gamma_res[1,2], p_gamma_hard = gamma_res[2,2], p_gamma_soft = gamma_res[3,2],
        est_beta_true = beta_res[1,1], est_beta_hard = beta_res[2,1], est_beta_soft = beta_res[3,1],
        est_gamma_true = gamma_res[1,1], est_gamma_hard = gamma_res[2,1], est_gamma_soft = gamma_res[3,1],
        se_beta_true = beta_res[1,3], se_beta_hard = beta_res[2,3], se_beta_soft = beta_res[3,3],
        se_gamma_true = gamma_res[1,3], se_gamma_hard = gamma_res[2,3], se_gamma_soft = gamma_res[3,3]
    )
}

# --- 4. Main Simulation Execution ---
N_SIMULATIONS <- 100
N_SAMPLES <- 2000
SNP_MAF <- 0.2
beta_effects <- seq(-0.25, 0.25, length.out = 5)
gamma_effects <- seq(-0.1, 0.1, length.out = 5)

# Define the three scenarios to run
scenarios_to_run <- tibble(
    scenario_name = c("high_separation", "low_separation", "low_separation_truncated"),
    base_params_name = c("high_separation", "low_separation", "low_separation"),
    is_truncated = c(FALSE, FALSE, TRUE)
)

# Create the full grid of all simulation runs
full_sim_grid <- scenarios_to_run %>%
    crossing(
        beta = beta_effects,
        gamma = gamma_effects,
        sim_id = 1:N_SIMULATIONS
    )

cat("--- Starting Comprehensive Simulation ---\n")
raw_results <- future_map_dfr(1:nrow(full_sim_grid), ~{
    set.seed(.x)
    params_row <- full_sim_grid[.x, ]
    antibody_params <- sim_params[[params_row$base_params_name]]
    
    sim_data <- generate_sim_data(N_SAMPLES, antibody_params, params_row$beta, params_row$gamma, SNP_MAF)
    
    # Apply truncation if specified for the scenario
    if (params_row$is_truncated) {
        lod_threshold <- exp(antibody_params$mu[1] - 0.5 * sqrt(antibody_params$sigma2[1]))
        sim_data <- sim_data %>% mutate(mfi = if_else(mfi < lod_threshold, lod_threshold, mfi))
    }
    
    analysis_res <- run_one_analysis(sim_data)
    
    # Return a single tibble with parameters and results
    bind_cols(
        select(params_row, scenario_name, beta, gamma, sim_id),
        analysis_res
    )
}, .options = furrr_options(seed = NULL), .progress = FALSE)

raw_results
# --- 5. Save RAW Results ---
cat("\n--- Simulation Complete. Saving RAW Results... ---\n")

write_tsv(raw_results, "raw_simulation_results.tsv")
saveRDS(sim_params, "simulation_params.rds")

cat("RAW results saved to raw_simulation_results.tsv\n")
