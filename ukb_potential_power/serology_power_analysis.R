################################################################################
#   UKB SEROLOGY: COMPREHENSIVE POWER ANALYSIS WITH HISTOGRAMS
################################################################################
#  ∙ This script provides a complete analysis of UKB serology data including:
#    1. Plain histograms with seropositivity cutoffs
#    2. Mixture model fitted histograms with component densities
#    3. Power analysis with inflation factors (1/λ²)
#  ∙ Uses a 2-component skew-t mixture model to estimate classification accuracy
#  ∙ Calculates sample size inflation factor for each antibody
################################################################################


# --- 1. Packages & Setup ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(mixsmsn, sn, dplyr, tibble, readr, purrr, ggplot2, patchwork)

stopifnot("smsn.mix" %in% ls("package:mixsmsn"), "pst" %in% ls("package:sn"))
`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- 2. Define All Antigens from Butler-Laporte et al. ---
# This tribble maps a short, usable name to the exact (messy) column header
# from the title file and includes the official MFI threshold.
antigen_map <- tribble(
  ~short_name, ~full_name, ~threshold,
  # Herpesviridae
  "hsv1", "1gG antigen for Herpes Simplex virus-1 | Instance 0", 150,
  "hsv2", "2mgG unique antigen for Herpes Simplex virus-2 | Instance 0", 150,
  "ebv_vca", "VCA p18 antigen for Epstein-Barr Virus | Instance 0", 250,
  "ebv_ebna1", "EBNA-1 antigen for Epstein-Barr Virus | Instance 0", 250,
  "ebv_zebra", "ZEBRA antigen for Epstein-Barr Virus | Instance 0", 100,
  "ebv_ead", "EA-D antigen for Epstein-Barr Virus | Instance 0", 100,
  "cmv_pp150", "pp150 Nter antigen for Human Cytomegalovirus | Instance 0", 100,
  "cmv_pp52", "pp 52 antigen for Human Cytomegalovirus | Instance 0", 150,
  "cmv_pp28", "pp 28 antigen for Human Cytomegalovirus | Instance 0", 200,
  "vzv", "gE / gI antigen for Varicella Zoster Virus  | Instance 0", 100,
  "hhv6_ie1a", "IE1A antigen for Human Herpesvirus-6 | Instance 0", 100,
  "hhv6_ie1b", "IE1B antigen for Human Herpesvirus-6 | Instance 0", 100,
  "hhv6_p101k", "p101 k antigen for Human Herpesvirus-6 | Instance 0", 100,
  "hhv7_u14", "U14 antigen for Human Herpesvirus-7 | Instance 0", 100,
  # Polyomaviridae
  "bkv_vp1", "BK VP1 antigen for Human Polyomavirus BKV | Instance 0", 250,
  "jcv_vp1", "JC VP1 antigen for Human Polyomavirus JCV | Instance 0", 250,
  "mcv_vp1", "MC VP1 antigen for Merkel Cell Polyomavirus | Instance 0", 250,
  # Bacteria
  "ct_mompa", "momp A antigen for Chlamydia trachomatis | Instance 0", 100,
  "ct_mompd", "momp D antigen for Chlamydia trachomatis | Instance 0", 100,
  "ct_tarpf1", "tarp-D F1 antigen for Chlamydia trachomatis | Instance 0", 100,
  "ct_tarpf2", "tarp-D F2 antigen for Chlamydia trachomatis | Instance 0", 100,
  "ct_porb", "PorB antigen for Chlamydia trachomatis | Instance 0", 80,
  "ct_pgp3", "pGP3 antigen for Chlamydia trachomatis | Instance 0", 200,
  # H. pylori
  "hp_caga", "CagA antigen for Helicobacter pylori | Instance 0", 400,
  "hp_vaca", "VacA antigen for Helicobacter pylori | Instance 0", 100,
  "hp_omp", "OMP antigen for Helicobacter pylori | Instance 0", 170,
  "hp_groel", "GroEL antigen for Helicobacter pylori | Instance 0", 80,
  "hp_catalase", "Catalase antigen for Helicobacter pylori | Instance 0", 180,
  "hp_urea", "UreA antigen for Helicobacter pylori | Instance 0", 130,
  # T. gondii
  "toxo_p22", "p22 antigen for Toxoplasma gondii | Instance 0", 100,
  "toxo_sag1", "sag1 antigen for Toxoplasma gondii | Instance 0", 160,
  # Additional antigens (added 2025-07-22)
  "kshv_k8_1", "K8.1 antigen for Kaposi's Sarcoma-Associated Herpesvirus | Instance 0", 175,
  "kshv_lana", "LANA antigen for Kaposi's Sarcoma-Associated Herpesvirus | Instance 0", 100,
  # Hepatitis B
  "hbv_hbc", "HBc antigen for Hepatitis B Virus | Instance 0", 100,
  "hbv_hbe", "HBe antigen for Hepatitis B Virus | Instance 0", 150,
  # Hepatitis C
  "hcv_core", "Core antigen for Hepatitis C Virus | Instance 0", 150,
  "hcv_ns3", "NS3 antigen for Hepatitis C Virus | Instance 0", 150,
  # HIV / HTLV
  "hiv_env", "HIV-1 env antigen for Human Immunodeficiency Virus | Instance 0", 150,
  "hiv_gag", "HIV-1 gag antigen for Human Immunodeficiency Virus | Instance 0", 600,
  "htlv_env", "HTLV-1 env antigen for Human T-Lymphotropic Virus 1 | Instance 0", 150,
  "htlv_gag", "HTLV-1 gag antigen for Human T-Lymphotropic Virus 1 | Instance 0", 1500,
  # Papillomavirus HPV16/18
  "hpv16_l1", "L1 antigen for Human Papillomavirus type-16 | Instance 0", 175,
  "hpv16_e6", "E6 antigen for Human Papillomavirus type-16 | Instance 0", 120,
  "hpv16_e7", "E7 antigen for Human Papillomavirus type-16 | Instance 0", 150,
  "hpv18_l1", "L1 antigen for Human Papillomavirus type-18 | Instance 0", 175
)

# --- 3. Load and Clean Data ---
cat("Loading data from 'serology_export_title.tsv'...\n")
data <- read_tsv("serology_export_title.tsv", show_col_types = FALSE, guess_max = 5000)

# Select and rename only the columns we need
cols_to_keep <- antigen_map$full_name
names(cols_to_keep) <- antigen_map$short_name

data_clean <- data %>%
  select(any_of(cols_to_keep)) %>%
  rename(any_of(cols_to_keep)) %>%
  mutate(across(everything(), as.numeric))
cat("Data loaded and cleaned for", ncol(data_clean), "antigens.\n")


# --- 4. Mixture Model & Plotting Utility Functions ---

#' Mixture model fitting function
fit_skew_mix <- function(y, fam = "Skew.t") {
  # Simplified for just skew-t
  tryCatch(
    mixsmsn::smsn.mix(
      y, g = 2,
      family   = fam,
      nu       = 4,
      group    = FALSE,
      calc.im  = FALSE,
      obs.prob = TRUE
    ) %>%
      { attr(., "family_used") <- fam; . },
    error = function(e) {
      warning("Failed to fit mixture model: ", e$message)
      return(NULL)
    }
  )
}

#' Calculate lambda statistics for power analysis
cdf_skew <- function(q, mu, sigma2, shape, fam, nu = 4) {
  omega <- sqrt(sigma2)
  sn::pst(q, xi = mu, omega = omega, alpha = shape, nu = nu)
}

lambda_stats <- function(fit, pos, cut_log, nu_default = 4) {
  fam <- attr(fit, "family_used") %||% "Skew.t"
  neg <- setdiff(1:2, pos)

  # Robustly get degrees of freedom, defaulting if needed
  get_nu <- function(k) {
    if (!is.null(fit$nu) && !is.na(fit$nu[k]) && fit$nu[k] > 0) fit$nu[k] else nu_default
  }
  nu_pos <- get_nu(pos)
  nu_neg <- get_nu(neg)

  # Calculate Se and Sp based on the fitted distributions and the hard cutoff
  Se <- 1 - cdf_skew(cut_log, fit$mu[pos], fit$sigma2[pos], fit$shape[pos], fam, nu_pos)
  Sp <-     cdf_skew(cut_log, fit$mu[neg], fit$sigma2[neg], fit$shape[neg], fam, nu_neg)

  # Calculate lambda and the inflation factor
  lambda  <- Se + Sp - 1
  tibble(
    Sensitivity = Se,
    Specificity = Sp,
    Lambda = lambda,
    Inflation_Factor = 1 / lambda^2
  )
}

#' Create plain histogram with cutoff
plot_plain_histogram <- function(data_vec, hard_cutoff, title) {
  plot_df <- tibble(mfi = data_vec)
  
  p <- ggplot(plot_df, aes(x = mfi)) +
    geom_histogram(aes(y = after_stat(density)), bins = 75, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = hard_cutoff, color = "red", linetype = "dashed", size = 1) +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = title,
      subtitle = paste0("Threshold = ", hard_cutoff, " (red dashed line)"),
      x = "MFI (log scale)",
      y = "Density"
    ) +
    theme_light()
    
  return(p)
}

#' Create mixture model diagnostic plot
plot_mixture_diagnosis <- function(data_vec, fit_model, hard_cutoff, title) {
  if (is.null(fit_model)) return(NULL)
  
  pos_comp <- which.max(fit_model$mu)
  neg_comp <- setdiff(1:2, pos_comp)
  
  # Create a dataframe for plotting
  plot_df <- tibble(mfi = data_vec)
  
  # Generate points for the density curves
  x_range <- range(log(plot_df$mfi[plot_df$mfi > 0]), na.rm = TRUE)
  x_points <- seq(x_range[1], x_range[2], length.out = 500)
  
  # Calculate weighted densities for each component
  density_pos <- dsn(x_points, fit_model$mu[pos_comp], fit_model$sigma2[pos_comp], fit_model$shape[pos_comp]) * fit_model$pii[pos_comp]
  density_neg <- dsn(x_points, fit_model$mu[neg_comp], fit_model$sigma2[neg_comp], fit_model$shape[neg_comp]) * fit_model$pii[neg_comp]
  
  density_df <- tibble(
    x = exp(x_points),
    positive = density_pos,
    negative = density_neg
  )
  
  # Create the plot
  p <- ggplot(plot_df, aes(x = mfi)) +
    geom_histogram(aes(y = after_stat(density)), bins = 75, fill = "grey", alpha = 0.6) +
    geom_area(data = density_df, aes(x = x, y = positive), fill = "red", alpha = 0.4) +
    geom_area(data = density_df, aes(x = x, y = negative), fill = "blue", alpha = 0.4) +
    geom_vline(xintercept = hard_cutoff, color = "black", linetype = "dashed", size = 1) +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = title,
      subtitle = paste0("Threshold = ", hard_cutoff, " (dashed line)"),
      x = "MFI (log scale)",
      y = "Density"
    ) +
    theme_light()
    
  return(p)
}


# --- 5. Main Analysis Loop ---
cat("\n--- Starting Comprehensive Analysis ---\n")
all_results <- list()
all_plain_plots <- list()
all_mixture_plots <- list()

# Iterate through the antigen map
for (i in 1:nrow(antigen_map)) {
  ab_name <- antigen_map$short_name[i]
  hard_cutoff <- antigen_map$threshold[i]
  
  # Check if the column exists in the cleaned data
  if (!ab_name %in% names(data_clean)) {
    cat("  -> Antigen '", ab_name, "' not found in data. Skipping.\n", sep="")
    next
  }
  
  cat("Analyzing:", ab_name, "\n")
  
  titres <- data_clean %>% pull(!!sym(ab_name)) %>% na.omit()
  
  if (length(titres[titres > 0]) < 200) { # Increased threshold for robustness
    cat("  -> Not enough valid data points. Skipping.\n")
    next
  }
  
  # Create plain histogram
  plain_plot <- plot_plain_histogram(titres, hard_cutoff, ab_name)
  all_plain_plots[[ab_name]] <- plain_plot
  
  # Fit mixture model and create diagnostic plot
  titres_log <- log(titres[titres > 0])
  fit <- fit_skew_mix(titres_log)
  
  mixture_plot <- plot_mixture_diagnosis(titres, fit, hard_cutoff, ab_name)
  if (!is.null(mixture_plot)) all_mixture_plots[[ab_name]] <- mixture_plot
  
  if (is.null(fit)) {
    cat("  -> Mixture model failed to converge. Skipping stats.\n")
    next
  }
  
  # Calculate power statistics
  pos_component <- which.max(fit$mu)
  stats <- lambda_stats(fit, pos_component, log(hard_cutoff))
  all_results[[ab_name]] <- stats %>% mutate(Antibody = ab_name, .before = 1)
}


# --- 6. Generate and Save Plots ---

# Save plain histograms
if (length(all_plain_plots) > 0) {
  cat("\nSaving", length(all_plain_plots), "plain histograms to PDF...\n")
  plots_per_page <- 8
  sorted_plots <- all_plain_plots[order(names(all_plain_plots))]
  
  pdf("ukb_antibody_distributions_final.pdf", width = 14, height = 14)
  for (i in seq(1, length(sorted_plots), by = plots_per_page)) {
    page_plot_list <- sorted_plots[i:min(i + plots_per_page - 1, length(sorted_plots))]
    combined_page <- wrap_plots(page_plot_list, ncol = 2)
    print(combined_page)
  }
  dev.off()
  cat("Plain histograms saved to 'ukb_antibody_distributions_final.pdf'\n")
}

# Save mixture model diagnostic plots
if (length(all_mixture_plots) > 0) {
  cat("\nSaving", length(all_mixture_plots), "mixture model diagnostic plots to PDF...\n")
  sorted_mixture_plots <- all_mixture_plots[order(names(all_mixture_plots))]
  
  pdf("serology_diagnostic_plots.pdf", width = 14, height = 14)
  for (i in seq(1, length(sorted_mixture_plots), by = plots_per_page)) {
    page_plot_list <- sorted_mixture_plots[i:min(i + plots_per_page - 1, length(sorted_mixture_plots))]
    combined_page <- wrap_plots(page_plot_list, ncol = 2)
    print(combined_page)
  }
  dev.off()
  cat("Mixture model diagnostic plots saved to 'serology_diagnostic_plots.pdf'\n")
}


# --- 7. Report Final Results ---
if (length(all_results) > 0) {
  summary_table <- bind_rows(all_results) %>% arrange(desc(Inflation_Factor))
  cat("\n\n--- Power Analysis Summary ---\n\n")
  print(summary_table, n = 50)
  cat("\n* Inflation_Factor: Estimated power loss from using a hard cutoff.\n")
  cat("* Lambda: Classification accuracy (Se + Sp - 1).\n")
  cat("* Higher Inflation_Factor = Greater power gain from mixture model.\n")
} else {
  cat("\nAnalysis complete, but no results were generated.\n")
}

cat("\n--- Analysis Complete ---\n")
cat("Files created:\n")
cat("- ukb_antibody_distributions_final.pdf: Plain histograms with cutoffs\n")
cat("- serology_diagnostic_plots.pdf: Mixture model fitted histograms\n")
cat("- Power analysis results printed above\n") 

write_tsv(summary_table, "power_analysis_results.tsv")