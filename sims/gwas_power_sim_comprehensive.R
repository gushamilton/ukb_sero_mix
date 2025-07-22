################################################################################
#   GWAS POWER SIMULATION: COMPREHENSIVE ANALYSIS OF BIAS AND POWER
################################################################################
#  ∙ This script extends the previous simulation to systematically explore how
#    varying genetic effect sizes (beta and gamma) influence statistical power
#    and estimation bias across different analysis methods.
#
#  ∙ It iterates over a grid of beta (effect on seropositivity) and gamma
#    (effect on antibody level) values to provide a comprehensive picture.
#
#  ∙ The primary output is a set of plots visualizing:
#    1. The power gain of using a mixture model over a hard cutoff.
#    2. The absolute bias in effect size estimation for each method.
################################################################################

# --- 1. Packages & Setup ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tibble, purrr, sn, mixsmsn, furrr, future, ggplot2, tidyr)

# --- Set up parallel processing ---
# Use multicore on Mac/Linux for efficiency. For Windows, use plan(multisession).
plan(multicore)

# --- 2. Simulation Parameters (from gwas_power_simulation.R) ---
# MODIFIED: Adjusted mu parameters for seropositive groups to better match
# the bimodal distribution seen in real-world data (e.g., the provided plot).
# The seropositive peaks are now at higher MFI values.

# Function to find optimal hard cutoff for fair comparison
find_optimal_cutoff <- function(params) {
  # Generate a large sample to estimate optimal cutoff
  set.seed(12345)  # Fixed seed for reproducibility
  n_test <- 100000
  component <- rbinom(n_test, 1, params$pii[2])
  mfi_log <- ifelse(component == 1,
                    rsn(n_test, xi = params$mu[2], omega = sqrt(params$sigma2[2]), alpha = params$shape[2]),
                    rsn(n_test, xi = params$mu[1], omega = sqrt(params$sigma2[1]), alpha = params$shape[1]))
  mfi <- exp(mfi_log)
  
  # Find cutoff that maximizes accuracy
  cutoffs <- quantile(mfi, seq(0.05, 0.95, 0.01))
  accuracies <- sapply(cutoffs, function(c) {
    predicted <- ifelse(mfi >= c, 1, 0)
    mean(predicted == component)
  })
  
  optimal_cutoff <- cutoffs[which.max(accuracies)]
  cat("Optimal cutoff:", round(optimal_cutoff, 2), "with accuracy:", round(max(accuracies), 3), "\n")
  return(optimal_cutoff)
}

# Base parameters (without hard cutoffs initially)
base_params <- list(
  cmv_pp150 = list(pii = c(0.45, 0.55), mu = c(3.7, 8.2), sigma2 = c(1.0, 0.8), shape = c(2, -4)),
  hsv1 = list(pii = c(0.3, 0.7), mu = c(3.0, 8.7), sigma2 = c(1.5, 1.0), shape = c(1.5, -3)),
  hsv1_overlap = list(pii = c(0.3, 0.7), mu = c(3.0, 4.5), sigma2 = c(1.5, 1.0), shape = c(1.5, -3))
)

# Calculate optimal cutoffs for each scenario
cat("--- Calculating Optimal Hard Cutoffs ---\n")
sim_params <- list()
for (name in names(base_params)) {
  cat("Calculating optimal cutoff for", name, "...\n")
  optimal_cutoff <- find_optimal_cutoff(base_params[[name]])
  sim_params[[name]] <- c(base_params[[name]], hard_cutoff = optimal_cutoff)
}

# --- 3. Core Simulation & Analysis Functions (from gwas_power_simulation.R) ---

#' Mixture model fitting function
fit_skew_mix <- function(y) {
  tryCatch(
    mixsmsn::smsn.mix(y, g = 2, family = "Skew.t", nu = 4, group = FALSE, calc.im = FALSE, obs.prob = TRUE),
    error = function(e) { NULL }
  )
}

#' Generate a simulated dataset
generate_sim_data <- function(n_samples, params, beta, gamma, maf) {
  snp <- rbinom(n_samples, 2, maf)

  # MODIFIED: Simulate a logistic effect for beta on the log-odds scale.
  # 'beta' now represents the per-allele change in the log-odds of being seropositive.
  base_log_odds <- qlogis(params$pii[2])
  prob_positive <- plogis(base_log_odds + (snp * beta))
  S_true <- rbinom(n_samples, 1, prob_positive)

  n_pos <- sum(S_true == 1)
  n_neg <- sum(S_true == 0)

  mfi_vec <- numeric(n_samples)

  # For true positives, the mean is shifted by the gamma effect
  if (n_pos > 0) {
    mfi_vec[S_true == 1] <- exp(rsn(n_pos, xi = params$mu[2] + (snp[S_true == 1] * gamma), omega = sqrt(params$sigma2[2]), alpha = params$shape[2]))
  }
  # For true negatives, we use the first component of the mixture
  if (n_neg > 0) {
    mfi_vec[S_true == 0] <- exp(rsn(n_neg, xi = params$mu[1], omega = sqrt(params$sigma2[1]), alpha = params$shape[1]))
  }

  tibble(snp, S_true, mfi = mfi_vec) %>%
    mutate(S_hard = ifelse(mfi >= params$hard_cutoff, 1, 0))
}

#' Run power and bias analysis for a single simulated dataset
run_one_analysis <- function(sim_data) {
  fit <- fit_skew_mix(log(sim_data$mfi[sim_data$mfi > 0]))
  p_soft_estimated <- NA
  if (!is.null(fit)) {
    pos_comp <- which.max(fit$mu)
    p_soft_estimated <- fit$obs.prob[, pos_comp]
  }
  sim_data$p_soft_estimated <- p_soft_estimated

  fit_beta_true <- glm(S_true ~ snp, data = sim_data, family = binomial())
  fit_beta_hard <- glm(S_hard ~ snp, data = sim_data, family = binomial())
  
  # MODIFIED: Use a quasibinomial GLM for the mixture model probabilities.
  # This estimates the effect on the log-odds scale, making it comparable to other methods.
  fit_beta_soft <- glm(p_soft_estimated ~ snp, data = sim_data, family = quasibinomial())

  fit_gamma_true <- lm(log(mfi) ~ snp, data = filter(sim_data, S_true == 1))
  fit_gamma_hard <- lm(log(mfi) ~ snp, data = filter(sim_data, S_hard == 1))
  fit_gamma_soft <- lm(log(mfi) ~ snp, data = sim_data, weights = p_soft_estimated)

  extract_coef <- function(model, coef_name) {
    if (is.null(model) || any(is.na(coef(model)))) return(c(NA, NA, NA))
    s <- summary(model)$coefficients
    if (coef_name %in% rownames(s)) {
      p_val_col <- if ("Pr(>|z|)" %in% colnames(s)) "Pr(>|z|)" else "Pr(>|t|)"
      return(s[coef_name, c("Estimate", p_val_col, "Std. Error")])
    }
    return(c(NA, NA, NA))
  }
  
  beta_res <- rbind(extract_coef(fit_beta_true, "snp"), extract_coef(fit_beta_hard, "snp"), extract_coef(fit_beta_soft, "snp"))
  gamma_res <- rbind(extract_coef(fit_gamma_true, "snp"), extract_coef(fit_gamma_hard, "snp"), extract_coef(fit_gamma_soft, "snp"))
  
  tibble(
    p_beta_true = beta_res[1,2], p_beta_hard = beta_res[2,2], p_beta_soft = beta_res[3,2],
    p_gamma_true = gamma_res[1,2], p_gamma_hard = gamma_res[2,2], p_gamma_soft = gamma_res[3,2],
    est_beta_true = beta_res[1,1], est_beta_hard = beta_res[2,1], est_beta_soft = beta_res[3,1],
    est_gamma_true = gamma_res[1,1], est_gamma_hard = gamma_res[2,1], est_gamma_soft = gamma_res[3,1],
    se_beta_true = beta_res[1,3], se_beta_hard = beta_res[2,3], se_beta_soft = beta_res[3,3],
    se_gamma_true = gamma_res[1,3], se_gamma_hard = gamma_res[2,3], se_gamma_soft = gamma_res[3,3]
  )
}

# --- 4. Main Simulation Driver (Modified for comprehensive grid search) ---
run_power_simulation <- function(n_sims, n_samples, antibody_name, beta, gamma, maf) {
  params <- sim_params[[antibody_name]]
  
  results <- future_map_dfr(1:n_sims, ~{
    sim_data <- generate_sim_data(n_samples, params, beta, gamma, maf)
    run_one_analysis(sim_data)
  }, .options = furrr_options(seed = TRUE))
  
  alpha <- 0.05
  
  # --- NEW: Calculate Precision (Mean Standard Error) ---
  precision_beta <- results %>% summarise(Analysis = "Seropositivity (Beta)", Mean_SE_Gold_Standard = mean(se_beta_true, na.rm=T), Mean_SE_Hard_Cutoff = mean(se_beta_hard, na.rm=T), Mean_SE_Mixture_Model = mean(se_beta_soft, na.rm=T))
  precision_gamma <- results %>% summarise(Analysis = "IgG Level (Gamma)", Mean_SE_Gold_Standard = mean(se_gamma_true, na.rm=T), Mean_SE_Hard_Cutoff = mean(se_gamma_hard, na.rm=T), Mean_SE_Mixture_Model = mean(se_gamma_soft, na.rm=T))
  precision_summary <- bind_rows(precision_beta, precision_gamma) %>% 
    mutate(Antibody = antibody_name, Beta_Effect = beta, Gamma_Effect = gamma, .before = 1)
  
  bias_beta <- results %>% summarise(Analysis = "Seropositivity (Beta)", True_Effect = beta, Mean_Est_Gold_Standard = mean(est_beta_true, na.rm=T), Mean_Est_Hard_Cutoff = mean(est_beta_hard, na.rm=T), Mean_Est_Mixture_Model = mean(est_beta_soft, na.rm=T))
  bias_gamma <- results %>% summarise(Analysis = "IgG Level (Gamma)", True_Effect = gamma, Mean_Est_Gold_Standard = mean(est_gamma_true, na.rm=T), Mean_Est_Hard_Cutoff = mean(est_gamma_hard, na.rm=T), Mean_Est_Mixture_Model = mean(est_gamma_soft, na.rm=T))
  bias_summary <- bind_rows(bias_beta, bias_gamma) %>% 
    mutate(Antibody = antibody_name, Beta_Effect = beta, Gamma_Effect = gamma, .before = 1)
  
  return(list(precision = precision_summary, bias = bias_summary))
}

# --- 5. Define and Execute Comprehensive Simulation ---
N_SIMULATIONS <- 20
N_SAMPLES <- 2000
SNP_MAF <- 0.2

# MODIFIED: Adjusted beta effect sizes for the new logistic (log-odds) scale.
# This range corresponds to per-allele odds ratios from ~0.78 to ~1.28.
# Now using a 5x5 grid to ensure 0 is included for the null effect.
beta_effects <- seq(-0.25, 0.25, length.out = 5)
gamma_effects <- seq(-0.1, 0.1, length.out = 5)

sim_grid <- expand_grid(
  antibody_name = c("cmv_pp150", "hsv1", "hsv1_overlap"),
  beta = beta_effects,
  gamma = gamma_effects
)

cat("--- Starting Comprehensive Simulation Grid ---\n")
cat("Total scenarios to run:", nrow(sim_grid), "\n")
cat("Simulations per scenario:", N_SIMULATIONS, "\n")




# Use future_map to iterate over the grid and show progress
all_results <- future_map(1:nrow(sim_grid), ~{
  params <- sim_grid[.x, ]
  run_power_simulation(
    n_sims = N_SIMULATIONS,
    n_samples = N_SAMPLES,
    antibody_name = params$antibody_name,
    beta = params$beta,
    gamma = params$gamma,
    maf = SNP_MAF
  )
}, .options = furrr_options(seed = TRUE), .progress = TRUE)

# --- 6. Process and Consolidate Results ---
precision_results <- map_dfr(all_results, "precision")
bias_results <- map_dfr(all_results, "bias")

cat("\n--- Simulation Complete. Generating Plots... ---\n")

# --- 7. Generate Plots (Revised) ---
simulation_results <- bind_rows(precision_results, bias_results)


# --- Plotting Function for PRECISION (Mean SE) ---
plot_precision <- function(antibody, effect_type, data) {
  
  plot_data <- data %>%
    filter(Antibody == antibody, Analysis == effect_type) %>%
    pivot_longer(
      cols = starts_with("Mean_SE_"),
      names_to = "Method",
      names_prefix = "Mean_SE_",
      values_to = "Mean_SE"
    ) %>%
    mutate(Method = gsub("_", " ", Method))
  
  if (effect_type == "Seropositivity (Beta)") {
    p <- ggplot(plot_data, aes(x = Beta_Effect, y = Mean_SE, color = Method, group = Method)) +
      facet_wrap(~Gamma_Effect, labeller = label_bquote(gamma == .(round(Gamma_Effect, 2)))) +
      labs(x = "True Beta Effect Size")
  } else {
    p <- ggplot(plot_data, aes(x = Gamma_Effect, y = Mean_SE, color = Method, group = Method)) +
      facet_wrap(~Beta_Effect, labeller = label_bquote(beta == .(round(Beta_Effect, 2)))) +
      labs(x = "True Gamma Effect Size")
  }
  
  p + geom_line() + geom_point() +
    labs(
      title = paste(toupper(antibody), ": Precision (Mean Standard Error) in", effect_type, "Estimation"),
      subtitle = "Lower SE indicates higher precision. Faceted by the other SNP effect.",
      y = "Mean Standard Error of Estimate",
      color = "Analysis Method"
    ) +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "lightblue"))
}

# --- Plotting Function for BIAS ---
plot_bias <- function(antibody, effect_type, data) {
  
  subtitle_text <- ""

  if (effect_type == "Seropositivity (Beta)") {
    plot_data <- data %>%
      filter(Antibody == antibody, Analysis == "Seropositivity (Beta)") %>%
      mutate(
        Bias_Gold_Standard = Mean_Est_Gold_Standard - True_Effect,
        Bias_Hard_Cutoff = Mean_Est_Hard_Cutoff - True_Effect,
        Bias_Mixture_Model = Mean_Est_Mixture_Model - True_Effect
      ) %>%
      select(Antibody, Beta_Effect, Gamma_Effect, Bias_Gold_Standard, Bias_Hard_Cutoff, Bias_Mixture_Model)
    
    x_var <- "Beta_Effect"
    facet_var <- "Gamma_Effect"
    x_lab <- "True Beta Effect Size"
    facet_lab <- "gamma"
    subtitle_text <- "Faceted by SNP effect on IgG Level (gamma)"
    
  } else { # IgG Level (Gamma)
    plot_data <- data %>%
      filter(Antibody == antibody, Analysis == "IgG Level (Gamma)") %>%
      mutate(
        Bias_Gold_Standard = Mean_Est_Gold_Standard - True_Effect,
        Bias_Hard_Cutoff = Mean_Est_Hard_Cutoff - True_Effect,
        Bias_Mixture_Model = Mean_Est_Mixture_Model - True_Effect
      ) %>%
      select(Antibody, Beta_Effect, Gamma_Effect, Bias_Gold_Standard, Bias_Hard_Cutoff, Bias_Mixture_Model)
      
    x_var <- "Gamma_Effect"
    facet_var <- "Beta_Effect"
    x_lab <- "True Gamma Effect Size"
    facet_lab <- "beta"
    subtitle_text <- "Faceted by SNP effect on Seropositivity (beta)"
  }
  
  plot_data <- plot_data %>%
    pivot_longer(
      cols = starts_with("Bias_"),
      names_to = "Method",
      names_prefix = "Bias_",
      values_to = "Bias"
    ) %>%
    mutate(Method = gsub("_", " ", Method))
  
  # Create a labeller expression that can capture the function's environment
  labeller_expr <- bquote(.(as.name(facet_lab)) == round(.(as.name(facet_var)), 2))
  
  ggplot(plot_data, aes(x = .data[[x_var]], y = Bias, color = Method, group = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line() + geom_point() +
    facet_wrap(vars(.data[[facet_var]]), labeller = as_labeller(labeller_expr)) +
    labs(
      title = paste(toupper(antibody), ": Bias in", effect_type, "Estimation"),
      subtitle = subtitle_text,
      x = x_lab,
      y = "Bias (Estimated - True)",
      color = "Analysis Method"
    ) +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill="lightcoral"))
}


# --- 8. Save Plots to PDF ---
pdf("comprehensive_power_bias_plots.pdf", width = 12, height = 8)


generate_sim_data(1000, sim_params$hsv1, 0.05, 0.06, 0.2) %>%
  ggplot(aes(x = mfi)) +
  geom_histogram() +
  scale_x_log10() +
  ggtitle("Simulated MFI Distribution for HSV-1")



generate_sim_data(1000, sim_params$cmv_pp150, 0.05, 0.06, 0.2) %>%
  ggplot(aes(x = mfi)) +
  geom_histogram() +
  scale_x_log10() +
  ggtitle("Simulated MFI Distribution for CMV")

# CMV Plots
print(plot_precision("cmv_pp150", "Seropositivity (Beta)", precision_results))
print(plot_bias("cmv_pp150", "Seropositivity (Beta)", bias_results))
print(plot_precision("cmv_pp150", "IgG Level (Gamma)", precision_results))
print(plot_bias("cmv_pp150", "IgG Level (Gamma)", bias_results))

# HSV1 Plots
print(plot_precision("hsv1", "Seropositivity (Beta)", precision_results))
print(plot_bias("hsv1", "Seropositivity (Beta)", bias_results))
print(plot_precision("hsv1", "IgG Level (Gamma)", precision_results))
print(plot_bias("hsv1", "IgG Level (Gamma)", bias_results))

# HSV1 Overlap Plots
print(plot_precision("hsv1_overlap", "Seropositivity (Beta)", precision_results))
print(plot_bias("hsv1_overlap", "Seropositivity (Beta)", bias_results))
print(plot_precision("hsv1_overlap", "IgG Level (Gamma)", precision_results))
print(plot_bias("hsv1_overlap", "IgG Level (Gamma)", bias_results))

dev.off()

cat("\n--- Plots saved to comprehensive_power_bias_plots.pdf ---\n") 


