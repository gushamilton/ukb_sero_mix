################################################################################
#   GWAS POWER SIMULATION (FINAL VERSION): FOUR-METHOD COMPARISON
################################################################################
#  ∙ This script runs the final, comprehensive simulation to compare four
#    distinct analysis methods across four challenging, realistic scenarios.
#  ∙ It saves the raw, unsummarized results for analysis by the plotting script.
#
#  ∙ SCENARIOS TESTED:
#    1. Standard Bimodal: A clean, baseline scenario.
#    2. Extreme Negative Skew: A highly skewed negative distribution.
#    3. Bimodal Gamma Overlap: Two complex, overlapping skewed distributions.
#    4. Bimodal Gamma Truncated: The above, with a hard floor on low values.
#
#  ∙ METHODS TESTED:
#    1. Gold Standard: Theoretical best-case using true serostatus.
#    2. Skew-T Mixture Model (Soft): The primary adaptive method being tested.
#    3. GMM Hard Cutoff: A robust, unsupervised hard-cutoff method.
#    4. Noisy External Cutoff: A realistic hard cutoff from a small panel.
################################################################################

# --- 1. Packages & Setup ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tibble, purrr, sn, mclust, furrr, future, readr, tidyr)

# --- Set up parallel processing ---
plan(multicore)

# --- 2. Simulation Parameters ---
# Final parameters from the visualization script
base_params <- list(
    standard_bimodal = list(type = "skew-normal", pii = c(0.4, 0.6), 
                               mu_neg = 4.0, sigma2_neg = 1.5, shape_neg = 1.5,
                               mu_pos = 7.5, sigma2_pos = 1.2, shape_pos = -3),
                               
    extreme_negative_skew = list(type = "gamma-skew-normal", pii = c(0.4, 0.6),
                                gamma_shape_neg = 0.8, gamma_rate_neg = 0.1, 
                                mu_pos = 5.0, sigma2_pos = 1.5, shape_pos = -3),

    bimodal_gamma_overlap = list(type = "gamma-gamma", pii = c(0.4, 0.6),
                         gamma_shape_neg = 1.5, gamma_rate_neg = 0.3,
                         gamma_shape_pos = 7.5, gamma_rate_pos = 0.6),
                         
    bimodal_gamma_truncated = list(type = "gamma-gamma-truncated", pii = c(0.4, 0.6),
                         gamma_shape_neg = 1.5, gamma_rate_neg = 0.3,
                         gamma_shape_pos = 7.5, gamma_rate_pos = 0.6,
                         truncation_floor = 1.5)
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
    if (params$type == "skew-normal") {
        if (n_pos > 0) mfi_vec[S_true == 1] <- exp(rsn(n_pos, xi = params$mu_pos + (snp[S_true == 1] * gamma), omega = sqrt(params$sigma2_pos), alpha = params$shape_pos))
        if (n_neg > 0) mfi_vec[S_true == 0] <- exp(rsn(n_neg, xi = params$mu_neg, omega = sqrt(params$sigma2_neg), alpha = params$shape_neg))
    } else if (params$type == "gamma-skew-normal") {
        if (n_pos > 0) mfi_vec[S_true == 1] <- exp(rsn(n_pos, xi = params$mu_pos + (snp[S_true == 1] * gamma), omega = sqrt(params$sigma2_pos), alpha = params$shape_pos))
        if (n_neg > 0) mfi_vec[S_true == 0] <- rgamma(n_neg, shape = params$gamma_shape_neg, rate = params$gamma_rate_neg)
    } else {
        if (n_pos > 0) mfi_vec[S_true == 1] <- rgamma(n_pos, shape = params$gamma_shape_pos, rate = params$gamma_rate_pos - (snp[S_true == 1] * gamma))
        if (n_neg > 0) mfi_vec[S_true == 0] <- rgamma(n_neg, shape = params$gamma_shape_neg, rate = params$gamma_rate_neg)
    }
    if (params$type == "gamma-gamma-truncated") {
        mfi_vec[mfi_vec < params$truncation_floor] <- params$truncation_floor
    }
    tibble(snp, S_true, mfi = mfi_vec)
}

# Helper function for Method 3: Find cutoff using a Gaussian Mixture Model
find_gmm_cutoff <- function(x) {
    gmm <- tryCatch(Mclust(x, G=2, verbose=FALSE, modelNames="V"), error = function(e) NULL)
    if (is.null(gmm)) return(NA)
    mu1 <- gmm$parameters$mean[1]; sig1 <- sqrt(gmm$parameters$variance$sigmasq[1])
    mu2 <- gmm$parameters$mean[2]; sig2 <- sqrt(gmm$parameters$variance$sigmasq[2])
    f <- function(x) dnorm(x, mu1, sig1) - dnorm(x, mu2, sig2)
    intersect <- tryCatch(uniroot(f, lower=min(mu1, mu2), upper=max(mu1, mu2))$root, error=function(e) NA)
    return(intersect)
}

# Helper function for Method 4: Find a cutoff from a small, noisy panel
find_external_cutoff <- function(params, n_panel = 300) {
    panel_data <- generate_sim_data(n_panel, params, 0, 0, 0.2)
    possible_cutoffs <- quantile(panel_data$mfi, seq(0.05, 0.95, 0.01), na.rm=TRUE)
    accuracies <- sapply(possible_cutoffs, function(c) mean((panel_data$mfi >= c) == panel_data$S_true, na.rm=T))
    return(possible_cutoffs[which.max(accuracies)])
}

# Main analysis function
run_one_analysis <- function(sim_data) {
    fit_skew_mix <- function(y) tryCatch(mixsmsn::smsn.mix(y, g = 2, family = "Skew.t", nu = 4, group = FALSE, calc.im = FALSE, obs.prob = TRUE), error = function(e) NULL)
    fit <- fit_skew_mix(log(sim_data$mfi[sim_data$mfi > 0 & is.finite(log(sim_data$mfi))]))
    p_soft_estimated <- if (!is.null(fit)) fit$obs.prob[, which.max(fit$mu)] else NA
    
    fit_beta_true <- glm(S_true ~ snp, data = sim_data, family = binomial())
    fit_beta_soft <- if (any(!is.na(p_soft_estimated))) glm(p_soft_estimated ~ snp, data = sim_data, family = quasibinomial()) else NULL
    fit_beta_gmm <- glm(S_gmm ~ snp, data = sim_data, family = binomial())
    fit_beta_external <- glm(S_external ~ snp, data = sim_data, family = binomial())

    fit_gamma_true <- if (sum(sim_data$S_true) > 2) lm(log(mfi) ~ snp, data = filter(sim_data, S_true == 1)) else NULL
    fit_gamma_soft <- if (any(!is.na(p_soft_estimated))) lm(log(mfi) ~ snp, data = sim_data, weights = p_soft_estimated) else NULL
    fit_gamma_gmm <- if (sum(sim_data$S_gmm) > 2) lm(log(mfi) ~ snp, data = filter(sim_data, S_gmm == 1)) else NULL
    fit_gamma_external <- if (sum(sim_data$S_external) > 2) lm(log(mfi) ~ snp, data = filter(sim_data, S_external == 1)) else NULL

    extract_coef <- function(model, coef_name = "snp") {
        if (is.null(model) || any(is.na(coef(model)))) return(c(NA, NA, NA))
        s <- summary(model)$coefficients
        p_val_col <- if ("Pr(>|z|)" %in% colnames(s)) "Pr(>|z|)" else "Pr(>|t|)"
        if (coef_name %in% rownames(s)) return(s[coef_name, c("Estimate", p_val_col, "Std. Error")])
        return(c(NA, NA, NA))
    }
    
    res <- list(
        beta_true = extract_coef(fit_beta_true), beta_soft = extract_coef(fit_beta_soft),
        beta_gmm = extract_coef(fit_beta_gmm), beta_external = extract_coef(fit_beta_external),
        gamma_true = extract_coef(fit_gamma_true), gamma_soft = extract_coef(fit_gamma_soft),
        gamma_gmm = extract_coef(fit_gamma_gmm), gamma_external = extract_coef(fit_gamma_external)
    )
    
    tibble(
        est_beta_true = res$beta_true[1], p_beta_true = res$beta_true[2], se_beta_true = res$beta_true[3],
        est_beta_soft = res$beta_soft[1], p_beta_soft = res$beta_soft[2], se_beta_soft = res$beta_soft[3],
        est_beta_gmm = res$beta_gmm[1], p_beta_gmm = res$beta_gmm[2], se_beta_gmm = res$beta_gmm[3],
        est_beta_external = res$beta_external[1], p_beta_external = res$beta_external[2], se_beta_external = res$beta_external[3],
        est_gamma_true = res$gamma_true[1], p_gamma_true = res$gamma_true[2], se_gamma_true = res$gamma_true[3],
        est_gamma_soft = res$gamma_soft[1], p_gamma_soft = res$gamma_soft[2], se_gamma_soft = res$gamma_soft[3],
        est_gamma_gmm = res$gamma_gmm[1], p_gamma_gmm = res$gamma_gmm[2], se_gamma_gmm = res$gamma_gmm[3],
        est_gamma_external = res$gamma_external[1], p_gamma_external = res$gamma_external[2], se_gamma_external = res$gamma_external[3]
    )
}

# --- 4. Main Simulation Execution ---
N_SIMULATIONS <- 500
N_SAMPLES <- 2000
SNP_MAF <- 0.2
beta_effects <- seq(-0.25, 0.25, length.out = 5)
gamma_effects <- seq(-0.1, 0.1, length.out = 5)

full_sim_grid <- crossing(
    scenario_name = names(base_params),
    beta = beta_effects,
    gamma = gamma_effects,
    sim_id = 1:N_SIMULATIONS
)

cat("--- Starting Comprehensive Simulation ---\n")
raw_results <- future_map_dfr(1:nrow(full_sim_grid), ~{
    set.seed(.x)
    params_row <- full_sim_grid[.x, ]
    antibody_params <- base_params[[params_row$scenario_name]]
    
    sim_data <- generate_sim_data(N_SAMPLES, antibody_params, params_row$beta, params_row$gamma, SNP_MAF)
    
    gmm_cutoff <- find_gmm_cutoff(log(sim_data$mfi))
    external_cutoff <- find_external_cutoff(antibody_params)
    
    sim_data <- sim_data %>%
        mutate(
            S_gmm = if (is.na(gmm_cutoff)) rep(0, n()) else as.integer(log(mfi) >= gmm_cutoff),
            S_external = as.integer(mfi >= external_cutoff)
        )
    
    analysis_res <- run_one_analysis(sim_data)
    
    bind_cols(select(params_row, -sim_id), analysis_res)
}, .options = furrr_options(seed = NULL), .progress = TRUE)

# --- 5. Save RAW Results ---
cat("\n--- Simulation Complete. Saving RAW Results... ---\n")
write_tsv(raw_results, "raw_simulation_results.tsv")
saveRDS(base_params, "simulation_params.rds")
cat("RAW results saved to raw_simulation_results.tsv\n")
