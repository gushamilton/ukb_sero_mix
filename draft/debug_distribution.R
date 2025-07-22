# --- Packages ---
library(ggplot2)
library(dplyr)
library(sn)

# --- Simulation Parameters ---
# These parameters are tuned to visually match the provided reference plot for cmv_pp150.
sim_params <- list(
  cmv_pp150 = list(
    pii = c(0.45, 0.55),      # Mixing proportions (neg, pos)
    mu = c(3.7, 8.2),         # Location parameters (log scale)
    sigma2 = c(1.0, 0.8),     # Scale parameters (variance)
    shape = c(2, -4),         # Skewness parameters
    hard_cutoff = 100
  )
)

# --- Data Generation Function ---
generate_sim_data <- function(n_samples, params) {
  # Simplified for debugging: ignore SNP effects
  S_true <- rbinom(n_samples, 1, params$pii[2])
  
  n_pos <- sum(S_true == 1)
  n_neg <- sum(S_true == 0)
  
  mfi_vec <- numeric(n_samples)
  
  log_mfi_neg <- rsn(n_neg, xi = params$mu[1], omega = sqrt(params$sigma2[1]), alpha = params$shape[1])
  log_mfi_pos <- rsn(n_pos, xi = params$mu[2], omega = sqrt(params$sigma2[2]), alpha = params$shape[2])
  
  mfi_vec[S_true == 0] <- exp(log_mfi_neg)
  mfi_vec[S_true == 1] <- exp(log_mfi_pos)
  
  tibble(S_true, mfi = mfi_vec, log_mfi = log(mfi_vec))
}

# --- Generate and Plot Data ---
set.seed(42)
N_SAMPLES <- 5000
params <- sim_params$cmv_pp150
sim_data <- generate_sim_data(N_SAMPLES, params)

# Function to get the scaled density of a component
scaled_dsn <- function(x, p, xi, omega, alpha) {
  p * dsn(x, xi = xi, omega = omega, alpha = alpha)
}

# --- Create the Plot ---
p <- ggplot(sim_data, aes(x = log_mfi)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "grey", color = "lightgrey") +
  stat_function(
    fun = scaled_dsn,
    args = list(p = params$pii[1], xi = params$mu[1], omega = sqrt(params$sigma2[1]), alpha = params$shape[1]),
    geom = "area", fill = "blue", alpha = 0.5
  ) +
  stat_function(
    fun = scaled_dsn,
    args = list(p = params$pii[2], xi = params$mu[2], omega = sqrt(params$sigma2[2]), alpha = params$shape[2]),
    geom = "area", fill = "red", alpha = 0.5
  ) +
  geom_vline(xintercept = log(params$hard_cutoff), linetype = "dashed", color = "black", linewidth = 1) +
  scale_x_continuous(
    name = "MFI (log scale)",
    breaks = log(c(1, 10, 100, 1000, 10000)),
    labels = c("1", "10", "100", "1 000", "10 000")
  ) +
  labs(
    title = "Simulated CMV pp150 MFI Distribution (Debug)",
    subtitle = paste("Threshold =", params$hard_cutoff, "(dashed line)"),
    y = "Density"
  ) +
  theme_bw()

# Save the plot
ggsave("debug_mfi_distribution.pdf", plot = p, width = 8, height = 6)

cat("Debug plot saved to debug_mfi_distribution.pdf\n")
cat("Please check this plot. If it looks correct, I will update the main simulation script with these parameters.\n") 