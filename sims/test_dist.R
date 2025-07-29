################################################################################
#   VISUALIZATION SCRIPT FOR SIMULATION SCENARIOS
################################################################################
#  ∙ This script's only purpose is to generate a large sample from each of the
#    four final, refined simulation scenarios and plot their distributions.
#  ∙ This allows for a visual check before running the full simulation.
################################################################################

# --- 1. Packages & Setup ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tibble, purrr, sn, ggplot2)

# --- 2. Simulation Parameters ---
# UPDATE: Parameters adjusted to slightly increase overlap in all scenarios.
base_params <- list(
    standard_bimodal = list(type = "skew-normal", pii = c(0.4, 0.6), 
                               mu_neg = 4.0, sigma2_neg = 1.5, shape_neg = 1.5,
                               mu_pos = 7.5, sigma2_pos = 1.2, shape_pos = -3), # mu_pos was 7.8
                               
    extreme_negative_skew = list(type = "gamma-skew-normal", pii = c(0.4, 0.6),
                                gamma_shape_neg = 0.8, gamma_rate_neg = 0.1, 
                                mu_pos = 5.0, sigma2_pos = 1.5, shape_pos = -3), # mu_pos was 5.3

    bimodal_gamma_overlap = list(type = "gamma-gamma", pii = c(0.4, 0.6),
                         gamma_shape_neg = 1.5, gamma_rate_neg = 0.3,
                         gamma_shape_pos = 7.5, gamma_rate_pos = 0.6),  # rate_pos was 0.5
                         
    bimodal_gamma_truncated = list(type = "gamma-gamma-truncated", pii = c(0.4, 0.6),
                         gamma_shape_neg = 1.5, gamma_rate_neg = 0.3,
                         gamma_shape_pos = 7.5, gamma_rate_pos = 0.6, # rate_pos was 0.5
                         truncation_floor = 1.5)
)

# --- 3. Data Generation Function ---
generate_sim_data <- function(n_samples, params, beta = 0, gamma = 0, maf = 0) {
    snp <- rbinom(n_samples, 2, maf)
    base_log_odds <- qlogis(params$pii[2])
    prob_positive <- plogis(base_log_odds + (snp * beta))
    S_true <- rbinom(n_samples, 1, prob_positive)
    mfi_vec <- numeric(n_samples)
    n_pos <- sum(S_true == 1); n_neg <- sum(S_true == 0)

    # Generate MFI based on scenario type
    if (params$type == "skew-normal") {
        if (n_pos > 0) mfi_vec[S_true == 1] <- exp(rsn(n_pos, xi = params$mu_pos + (snp[S_true == 1] * gamma), omega = sqrt(params$sigma2_pos), alpha = params$shape_pos))
        if (n_neg > 0) mfi_vec[S_true == 0] <- exp(rsn(n_neg, xi = params$mu_neg, omega = sqrt(params$sigma2_neg), alpha = params$shape_neg))
    } else if (params$type == "gamma-skew-normal") {
        if (n_pos > 0) mfi_vec[S_true == 1] <- exp(rsn(n_pos, xi = params$mu_pos + (snp[S_true == 1] * gamma), omega = sqrt(params$sigma2_pos), alpha = params$shape_pos))
        if (n_neg > 0) mfi_vec[S_true == 0] <- rgamma(n_neg, shape = params$gamma_shape_neg, rate = params$gamma_rate_neg)
    } else { # Handles both gamma-gamma and gamma-gamma-truncated
        # Note: The 'gamma' effect is modeled as a shift on the rate parameter here
        if (n_pos > 0) mfi_vec[S_true == 1] <- rgamma(n_pos, shape = params$gamma_shape_pos, rate = params$gamma_rate_pos - (snp[S_true == 1] * gamma))
        if (n_neg > 0) mfi_vec[S_true == 0] <- rgamma(n_neg, shape = params$gamma_shape_neg, rate = params$gamma_rate_neg)
    }
    
    # Apply truncation for the specific scenario
    if (params$type == "gamma-gamma-truncated") {
        mfi_vec[mfi_vec < params$truncation_floor] <- params$truncation_floor
    }
    
    tibble(snp, S_true, mfi = mfi_vec)
}

# --- 4. Generate and Plot Data ---

# Generate a large sample (N=20,000) for each scenario
plot_data <- map_dfr(names(base_params), ~{
    generate_sim_data(20000, base_params[[.x]]) %>%
        mutate(
            Scenario = tools::toTitleCase(gsub("_", " ", .x)),
            S_true_label = if_else(S_true == 1, "Seropositive", "Seronegative")
        )
})

# Create the plot
distribution_plot <- ggplot(plot_data, aes(x = mfi, fill = S_true_label)) +
    geom_histogram(bins = 150, alpha = 0.7, position = "identity") +
    facet_wrap(~Scenario, scales = "free", ncol = 2) +
    scale_x_log10(labels = scales::comma, n.breaks = 7) +
    labs(
        title = "Sample MFI Distributions for Each Proposed Scenario",
        subtitle = "Data generated with no genetic effect (beta=0, gamma=0)",
        x = "MFI (log scale)",
        y = "Count",
        fill = "True Serostatus"
    ) +
    scale_fill_manual(values = c("Seronegative" = "#F8766D", "Seropositive" = "#00BFC4")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", strip.background = element_rect(fill="lightblue"))

# Print the plot to the screen
print(distribution_plot)

# Optionally, save the plot to a file
 ggsave("scenario_distributions.png", plot = distribution_plot, width = 12, height = 10)

