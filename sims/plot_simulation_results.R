################################################################################
#   GWAS POWER SIMULATION: ANALYSIS & PLOTTING SCRIPT (FINAL VERSION)
################################################################################
#  ∙ This script loads RAW simulation results from the final simulation script
#    and performs all summarization and plotting for the four-method,
#    four-scenario comparison.
################################################################################

# --- 1. Packages & Setup ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tibble, purrr, ggplot2, tidyr, readr, patchwork, scales, sn, grid, gridExtra)

# --- 2. Load RAW Simulation Results ---
cat("Loading RAW simulation results...\n")
raw_results <- read_tsv("raw_simulation_results.tsv.gz", show_col_types = FALSE)
sim_params <- readRDS("simulation_params.rds")

# --- 3. Plotting & Table Functions ---

# Function to generate sample data (needed for histogram plot)
generate_sim_data <- function(n_samples, params, beta = 0, gamma = 0, maf = 0) {
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

#' Plot sample MFI distributions for each scenario
plot_scenario_distributions <- function(sim_params) {
    cat("Generating sample MFI distribution plots...\n")
    plot_data <- map_dfr(names(sim_params), ~{
        generate_sim_data(20000, sim_params[[.x]]) %>%
            mutate(Scenario = tools::toTitleCase(gsub("_", " ", .x)), 
                   S_true_label = if_else(S_true == 1, "Seropositive", "Seronegative"))
    })
    ggplot(plot_data, aes(x = mfi, fill = S_true_label)) +
        geom_histogram(bins = 150, alpha = 0.7, position = "identity") + facet_wrap(~Scenario, scales = "free", ncol=2) +
        scale_x_log10(labels = scales::comma, n.breaks=6) +
        labs(title = "Sample MFI Distributions for Each Scenario", subtitle = "Data generated with no genetic effect (beta=0, gamma=0)", x = "MFI (log scale)", y = "Count", fill = "True Serostatus") +
        scale_fill_manual(values = c("Seronegative" = "#F8766D", "Seropositive" = "#00BFC4")) + theme_bw() + theme(legend.position = "bottom")
}

#' Plot a metric (Bias or Mean SE) with shaded error ribbons
plot_simulation_metric <- function(scen_name, raw_data, metric) {
    metric_label <- if (metric == "Bias") "Bias (Estimated - True)" else "Mean Standard Error"
    target_cols <- if (metric == "Bias") names(select(raw_data, starts_with("est_"))) else names(select(raw_data, starts_with("se_")))

    summary_data <- raw_data %>%
        filter(scenario_name == {{scen_name}}) %>%
        pivot_longer(cols = all_of(target_cols), names_to = "Metric", values_to = "Value") %>%
        mutate(Value = case_when(
            metric == "Bias" & grepl("beta", Metric) ~ Value - beta,
            metric == "Bias" & grepl("gamma", Metric) ~ Value - gamma,
            TRUE ~ Value
        )) %>%
        group_by(beta, gamma, Metric) %>%
        summarise(mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE), .groups = "drop") %>%
        mutate(
            Effect = if_else(grepl("beta", Metric), "Beta", "Gamma"),
            Method = factor(case_when(
                grepl("gmm", Metric) ~ "GMM Cutoff",
                grepl("external", Metric) ~ "Noisy External Cutoff",
                grepl("soft", Metric) ~ "Mixture Model",
                grepl("true", Metric) ~ "Gold Standard"
            ), levels = c("Gold Standard", "Mixture Model", "Noisy External Cutoff", "GMM Cutoff"))
        )

    p_beta <- ggplot(filter(summary_data, Effect == "Beta"), aes(x = beta, y = mean, color = Method, fill = Method)) +
        geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, linetype = 0) + geom_line() + geom_point() +
        facet_wrap(~gamma, labeller = label_bquote(gamma == .(round(gamma, 2)))) + labs(x = "True Beta Effect Size", y = metric_label)
    p_gamma <- ggplot(filter(summary_data, Effect == "Gamma"), aes(x = gamma, y = mean, color = Method, fill = Method)) +
        geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, linetype = 0) + geom_line() + geom_point() +
        facet_wrap(~beta, labeller = label_bquote(beta == .(round(beta, 2)))) + labs(x = "True Gamma Effect Size", y = metric_label)
    if (metric == "Bias") { p_beta <- p_beta + geom_hline(yintercept = 0, linetype = "dashed", color = "grey50"); p_gamma <- p_gamma + geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") }
    (p_beta + theme_bw() + theme(legend.position = "none")) / (p_gamma + theme_bw() + theme(legend.position = "bottom")) + plot_layout(guides = "collect") +
    plot_annotation(title = paste(toupper(gsub("_", " ", scen_name)), ":", metric_label), caption = "Error ribbons show Mean +/- 1 SD.")
}

#' Plot statistical power
plot_power <- function(scen_name, raw_data) {
    power_data <- raw_data %>%
        filter(scenario_name == {{scen_name}}) %>%
        select(beta, gamma, starts_with("p_")) %>%
        pivot_longer(cols = starts_with("p_"), names_to = "Metric", values_to = "PValue") %>%
        group_by(beta, gamma, Metric) %>%
        summarise(Power = mean(PValue < 0.05, na.rm = TRUE), .groups = "drop") %>%
        mutate(
            Effect = if_else(grepl("beta", Metric), "Beta", "Gamma"),
            Method = factor(case_when(
                grepl("gmm", Metric) ~ "GMM Cutoff",
                grepl("external", Metric) ~ "Noisy External Cutoff",
                grepl("soft", Metric) ~ "Mixture Model",
                grepl("true", Metric) ~ "Gold Standard"
            ), levels = c("Gold Standard", "Mixture Model", "Noisy External Cutoff", "GMM Cutoff"))
        )
    p_beta <- ggplot(filter(power_data, Effect == "Beta"), aes(x = beta, y = Power, color = Method, group = Method)) +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + geom_hline(yintercept = 0.80, linetype = "dashed", color = "blue") +
        geom_line() + geom_point() + facet_wrap(~gamma, labeller = label_bquote(gamma == .(round(gamma, 2)))) +
        scale_y_continuous(labels = scales::percent) + labs(x = "True Beta Effect Size", y = "Power (or Type I Error)")
    p_gamma <- ggplot(filter(power_data, Effect == "Gamma"), aes(x = gamma, y = Power, color = Method, group = Method)) +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + geom_hline(yintercept = 0.80, linetype = "dashed", color = "blue") +
        geom_line() + geom_point() + facet_wrap(~beta, labeller = label_bquote(beta == .(round(beta, 2)))) +
        scale_y_continuous(labels = scales::percent) + labs(x = "True Gamma Effect Size", y = "Power (or Type I Error)")
    (p_beta + theme_bw() + theme(legend.position = "none")) / (p_gamma + theme_bw() + theme(legend.position = "bottom")) + plot_layout(guides = "collect") +
    plot_annotation(title = paste(toupper(gsub("_", " ", scen_name)), ": Statistical Power"), caption = "Red line = 5% Type I Error rate. Blue line = 80% Power.")
}

#' Create and draw a summary table
draw_summary_table <- function(table_data, title) {
    grob <- tableGrob(table_data, rows = NULL, theme = ttheme_minimal(base_size = 8))
    grid.newpage(); grid.draw(grob); grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
}

# --- 4. Save All Plots and Tables to PDF ---
pdf("simulation_analysis_final.pdf", width = 12, height = 10)

# Page 1: Sample MFI Distributions
print(plot_scenario_distributions(sim_params))

# Page 2: Type I Error Table
t1e_summary <- raw_results %>%
    filter(beta == 0, gamma == 0) %>%
    select(scenario_name, starts_with("p_")) %>%
    group_by(Scenario = scenario_name) %>%
    summarise(across(everything(), ~mean(.x < 0.05, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(cols = -Scenario, names_to = "Metric", values_to = "T1E_Rate") %>%
    mutate(
        Effect = if_else(grepl("beta", Metric), "Beta", "Gamma"),
        Method = factor(case_when(
            grepl("gmm", Metric) ~ "GMM Cutoff",
            grepl("external", Metric) ~ "Noisy External",
            grepl("soft", Metric) ~ "Mixture Model",
            grepl("true", Metric) ~ "Gold Standard"
        )),
        T1E_Rate = scales::percent(T1E_Rate, accuracy = 0.1)
    ) %>%
    select(Scenario, Effect, Method, T1E_Rate) %>%
    pivot_wider(names_from = Method, values_from = T1E_Rate) %>%
    mutate(Scenario = tools::toTitleCase(gsub("_", " ", Scenario)))
draw_summary_table(t1e_summary, "Type I Error Rate (at Alpha = 0.05)")

# Subsequent Pages: Detailed plots for each scenario
scenarios <- unique(raw_results$scenario_name)
for (scen in scenarios) {
    print(plot_simulation_metric(scen, raw_results, "Bias"))
    print(plot_simulation_metric(scen, raw_results, "Mean_SE"))
    print(plot_power(scen, raw_results))
}

dev.off()
cat("\n--- All plots and tables saved to simulation_analysis_final.pdf ---\n")
