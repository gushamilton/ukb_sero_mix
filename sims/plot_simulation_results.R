################################################################################
#   GWAS POWER SIMULATION (V8): FOUR-METHOD COMPARISON
################################################################################
#  ∙ This script runs a comprehensive simulation to compare four distinct
#    analysis methods and saves the raw, unsummarized results.
#
#  ∙ UPDATE:
#    - Increased the severity of the left-truncation for the 'truncated'
#      scenario to make it a more distinct and challenging test case.
################################################################################

# --- 1. Packages & Setup ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tibble, purrr, sn, mixsmsn, furrr, future, readr, tidyr, pracma)

# --- Set up parallel processing ---
plan(multicore)

# --- 2. Simulation Parameters ---
base_params <- list(
    high_separation = list(pii = c(0.4, 0.6), mu = c(4.0, 8.0), sigma2 = c(1.5, 1.2), shape = c(1.5, -3)),
    low_separation  = list(pii = c(0.4, 0.6), mu = c(4.0, 6.2), sigma2 = c(1.5, 1.2), shape = c(1.5, -3))
)

# --- 3. Core Simulation & Analysis Functions ---

# Helper function to generate data
generate_sim_data <- function(n_samples, params, beta, gamma, maf) {
    snp <- rbinom(n_samples, 2, maf)
    base_log_odds <- qlogis(params$pii[2])
    prob_positive <- plogis(base_log_odds + (snp * beta))
    S_true <- rbinom(n_samples, 1, prob_positive)
    mfi_vec <- numeric(n_samples)
    n_pos <- sum(S_true == 1); n_neg <- sum(S_true == 0)
    if (n_pos > 0) mfi_vec[S_true == 1] <- exp(rsn(n_pos, xi = params$mu[2] + (snp[S_true == 1] * gamma), omega = sqrt(params$sigma2[2]), alpha = params$shape[2]))
    if (n_neg > 0) mfi_vec[S_true == 0] <- exp(rsn(n_neg, xi = params$mu[1], omega = sqrt(params$sigma2[1]), alpha = params$shape[1]))
    tibble(snp, S_true, mfi = mfi_vec)
}

# Helper function for Method 3: Find the valley between two peaks
find_density_valley <- function(x) {
    d <- density(x, na.rm = TRUE)
    pks <- findpeaks(d$y)
    if (nrow(pks) < 2) return(NA) # Not bimodal, can't find a valley
    pks <- pks[order(pks[,1], decreasing = TRUE), ] # Sort by peak height
    
    valley_idx <- which.min(d$y[pks[2,2]:pks[1,2]]) + pks[2,2] - 1
    return(d$x[valley_idx])
}

# Helper function for Method 4: Find a cutoff from a small, noisy panel
find_external_cutoff <- function(params, n_panel = 300) {
    panel_data <- generate_sim_data(n_panel, params, 0, 0, 0.2)
    
    # Find cutoff that maximizes accuracy IN THE SMALL PANEL
    possible_cutoffs <- quantile(panel_data$mfi, seq(0.05, 0.95, 0.01))
    accuracies <- sapply(possible_cutoffs, function(c) {
        mean((panel_data$mfi >= c) == panel_data$S_true)
    })
    
    return(possible_cutoffs[which.max(accuracies)])
}

# Main analysis function, now runs all 4 methods
run_one_analysis <- function(sim_data) {
    # Method 2: Mixture Model
    fit_skew_mix <- function(y) tryCatch(mixsmsn::smsn.mix(y, g = 2, family = "Skew.t", nu = 4, group = FALSE, calc.im = FALSE, obs.prob = TRUE), error = function(e) NULL)
    fit <- fit_skew_mix(log(sim_data$mfi[sim_data$mfi > 0 & is.finite(log(sim_data$mfi))]))
    p_soft_estimated <- if (!is.null(fit)) fit$obs.prob[, which.max(fit$mu)] else NA
    
    # Run GLMs for all methods
    fit_beta_true <- glm(S_true ~ snp, data = sim_data, family = binomial())
    fit_beta_soft <- if (any(!is.na(p_soft_estimated))) glm(p_soft_estimated ~ snp, data = sim_data, family = quasibinomial()) else NULL
    fit_beta_valley <- glm(S_valley ~ snp, data = sim_data, family = binomial())
    fit_beta_external <- glm(S_external ~ snp, data = sim_data, family = binomial())

    fit_gamma_true <- if (sum(sim_data$S_true) > 2) lm(log(mfi) ~ snp, data = filter(sim_data, S_true == 1)) else NULL
    fit_gamma_soft <- if (any(!is.na(p_soft_estimated))) lm(log(mfi) ~ snp, data = sim_data, weights = p_soft_estimated) else NULL
    fit_gamma_valley <- if (sum(sim_data$S_valley) > 2) lm(log(mfi) ~ snp, data = filter(sim_data, S_valley == 1)) else NULL
    fit_gamma_external <- if (sum(sim_data$S_external) > 2) lm(log(mfi) ~ snp, data = filter(sim_data, S_external == 1)) else NULL

    # Helper to safely extract coefficients
    extract_coef <- function(model, coef_name = "snp") {
        if (is.null(model) || any(is.na(coef(model)))) return(c(NA, NA, NA))
        s <- summary(model)$coefficients
        p_val_col <- if ("Pr(>|z|)" %in% colnames(s)) "Pr(>|z|)" else "Pr(>|t|)"
        if (coef_name %in% rownames(s)) return(s[coef_name, c("Estimate", p_val_col, "Std. Error")])
        return(c(NA, NA, NA))
    }
    
    res <- list(
        beta_true = extract_coef(fit_beta_true), beta_soft = extract_coef(fit_beta_soft),
        beta_valley = extract_coef(fit_beta_valley), beta_external = extract_coef(fit_beta_external),
        gamma_true = extract_coef(fit_gamma_true), gamma_soft = extract_coef(fit_gamma_soft),
        gamma_valley = extract_coef(fit_gamma_valley), gamma_external = extract_coef(fit_gamma_external)
    )
    
    # Combine results into a single-row tibble
    tibble(
        est_beta_true = res$beta_true[1], p_beta_true = res$beta_true[2], se_beta_true = res$beta_true[3],
        est_beta_soft = res$beta_soft[1], p_beta_soft = res$beta_soft[2], se_beta_soft = res$beta_soft[3],
        est_beta_valley = res$beta_valley[1], p_beta_valley = res$beta_valley[2], se_beta_valley = res$beta_valley[3],
        est_beta_external = res$beta_external[1], p_beta_external = res$beta_external[2], se_beta_external = res$beta_external[3],
        est_gamma_true = res$gamma_true[1], p_gamma_true = res$gamma_true[2], se_gamma_true = res$gamma_true[3],
        est_gamma_soft = res$gamma_soft[1], p_gamma_soft = res$gamma_soft[2], se_gamma_soft = res$gamma_soft[3],
        est_gamma_valley = res$gamma_valley[1], p_gamma_valley = res$gamma_valley[2], se_gamma_valley = res$gamma_valley[3],
        est_gamma_external = res$gamma_external[1], p_gamma_external = res$gamma_external[2], se_gamma_external = res$gamma_external[3]
    )
}

# --- 4. Main Simulation Execution ---
N_SIMULATIONS <- 500
N_SAMPLES <- 2000
SNP_MAF <- 0.2
beta_effects <- seq(-0.25, 0.25, length.out = 5)
gamma_effects <- seq(-0.1, 0.1, length.out = 5)

scenarios_to_run <- tibble(
    scenario_name = c("high_separation", "low_separation", "low_separation_truncated"),
    base_params_name = c("high_separation", "low_separation", "low_separation"),
    is_truncated = c(FALSE, FALSE, TRUE)
)

full_sim_grid <- scenarios_to_run %>%
    crossing(beta = beta_effects, gamma = gamma_effects, sim_id = 1:N_SIMULATIONS)

cat("--- Starting Comprehensive Simulation ---\n")
raw_results <- future_map_dfr(1:nrow(full_sim_grid), ~{
    set.seed(.x)
    params_row <- full_sim_grid[.x, ]
    antibody_params <- base_params[[params_row$base_params_name]]
    
    # Generate main data
    sim_data <- generate_sim_data(N_SAMPLES, antibody_params, params_row$beta, params_row$gamma, SNP_MAF)
    
    # Apply truncation if specified
    if (params_row$is_truncated) {
        # Increased truncation from 0.5 to 1.0 standard deviations below the mean
        lod_threshold <- exp(antibody_params$mu[1] - 1.0 * sqrt(antibody_params$sigma2[1]))
        sim_data <- sim_data %>% mutate(mfi = if_else(mfi < lod_threshold, lod_threshold, mfi))
    }
    
    # Derive cutoffs for hard-cutoff methods
    valley_cutoff <- find_density_valley(log(sim_data$mfi))
    external_cutoff <- find_external_cutoff(antibody_params)
    
    # Classify individuals using the hard cutoffs
    sim_data <- sim_data %>%
        mutate(
            S_valley = if (is.na(valley_cutoff)) {
                           rep(0, n()) # If no valley, classify all as negative
                       } else {
                           as.integer(log(mfi) >= valley_cutoff)
                       },
            S_external = as.integer(mfi >= external_cutoff)
        )
    
    analysis_res <- run_one_analysis(sim_data)
    
    bind_cols(select(params_row, scenario_name, beta, gamma, sim_id), analysis_res)
}, .options = furrr_options(seed = NULL), .progress = TRUE)

# --- 5. Save RAW Results ---
cat("\n--- Simulation Complete. Saving RAW Results... ---\n")
write_tsv(raw_results, "raw_simulation_results.tsv")
saveRDS(base_params, "simulation_params.rds") # Save base params for plotting
cat("RAW results saved to raw_simulation_results.tsv\n")
