################################################################################
#   GWAS POWER SIMULATION: ANALYSIS & PLOTTING SCRIPT (V5)
################################################################################
#  ∙ This script loads RAW simulation results and performs all summarization
#    and plotting for the focused scenarios.
################################################################################

# --- 1. Packages & Setup ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tibble, purrr, ggplot2, tidyr, readr, patchwork, scales, sn, grid, gridExtra, Hmisc)

# --- 2. Load RAW Simulation Results ---
cat("Loading RAW simulation results...\n")
raw_results <- read_tsv("raw_simulation_results.tsv", show_col_types = FALSE)
sim_params <- readRDS("simulation_params.rds")

# --- 3. Plotting & Table Functions ---

#' Plot a metric (Bias or Mean SE) with shaded error ribbons
plot_simulation_metric <- function(scen_name, raw_data, metric) {
    
    metric_label <- if (metric == "Bias") "Bias (Estimated - True)" else "Mean Standard Error"
    
    # Select the correct columns for the chosen metric
    if (metric == "Bias") {
        est_cols <- c("est_beta_true", "est_beta_hard", "est_beta_soft", "est_gamma_true", "est_gamma_hard", "est_gamma_soft")
    } else { # Mean_SE
        est_cols <- c("se_beta_true", "se_beta_hard", "se_beta_soft", "se_gamma_true", "se_gamma_hard", "se_gamma_soft")
    }

    # Summarize the raw data *inside* the function
    summary_data <- raw_data %>%
        filter(scenario_name == {{scen_name}}) %>%
        pivot_longer(cols = all_of(est_cols), names_to = "Metric", values_to = "Value") %>%
        mutate(
            Value = if(metric == "Bias" & grepl("beta", Metric)) Value - beta else Value,
            Value = if(metric == "Bias" & grepl("gamma", Metric)) Value - gamma else Value
        ) %>%
        group_by(beta, gamma, Metric) %>%
        summarise(mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE), .groups = "drop") %>%
        mutate(
            Effect = if_else(grepl("beta", Metric), "Beta", "Gamma"),
            Method = factor(case_when(
                grepl("hard", Metric) ~ "Hard Cutoff",
                grepl("soft", Metric) ~ "Mixture Model",
                grepl("true", Metric) ~ "Gold Standard"
            ), levels = c("Gold Standard", "Mixture Model", "Hard Cutoff"))
        )

    # Plot for Beta effect
    p_beta <- ggplot(filter(summary_data, Effect == "Beta"), aes(x = beta, y = mean, color = Method, fill = Method)) +
        geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, linetype = 0) +
        geom_line() + geom_point() +
        facet_wrap(~gamma, labeller = label_bquote(gamma == .(round(gamma, 2)))) +
        labs(x = "True Beta Effect Size", y = metric_label)
        
    # Plot for Gamma effect
    p_gamma <- ggplot(filter(summary_data, Effect == "Gamma"), aes(x = gamma, y = mean, color = Method, fill = Method)) +
        geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, linetype = 0) +
        geom_line() + geom_point() +
        facet_wrap(~beta, labeller = label_bquote(beta == .(round(beta, 2)))) +
        labs(x = "True Gamma Effect Size", y = metric_label)
        
    if (metric == "Bias") {
        p_beta <- p_beta + geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")
        p_gamma <- p_gamma + geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")
    }

    (p_beta / p_gamma) + plot_layout(guides = "collect") & theme_bw() & theme(legend.position = "bottom") &
    plot_annotation(title = paste(toupper(gsub("_", " ", scen_name)), ":", metric_label))
}

#' Create and draw a summary table
draw_summary_table <- function(table_data, title) {
    grob <- tableGrob(table_data, rows = NULL, theme = ttheme_minimal(base_size = 9))
    grid.newpage()
    grid.draw(grob)
    grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))
}

# --- 4. Save All Plots and Tables to PDF ---
pdf("simulation_analysis_final.pdf", width = 12, height = 10)

# Page 1: Type I Error Table
t1e_summary <- raw_results %>%
    filter(beta == 0, gamma == 0) %>%
    select(scenario_name, starts_with("p_")) %>%
    group_by(Scenario = scenario_name) %>%
    summarise(across(everything(), ~mean(.x < 0.05, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(cols = -Scenario, names_to = "Metric", values_to = "T1E_Rate") %>%
    mutate(
        Effect = if_else(grepl("beta", Metric), "Beta", "Gamma"),
        Method = Hmisc::capitalize(sub(".*_(.*)", "\\1", Metric)),
        T1E_Rate = scales::percent(T1E_Rate, accuracy = 0.1)
    ) %>%
    select(Scenario, Effect, Method, T1E_Rate) %>%
    pivot_wider(names_from = Method, values_from = T1E_Rate) %>%
    mutate(Scenario = Hmisc::capitalize(gsub("_", " ", Scenario)))

draw_summary_table(t1e_summary, "Type I Error Rate (at Alpha = 0.05)")

# Page 2 onwards: Detailed plots
scenarios <- unique(raw_results$scenario_name)
for (scen in scenarios) {
    print(plot_simulation_metric(scen, raw_results, "Bias"))
    print(plot_simulation_metric(scen, raw_results, "Mean_SE"))
}

dev.off()
cat("\n--- All plots and tables saved to simulation_analysis_final.pdf ---\n")
