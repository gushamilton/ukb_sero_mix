################################################################################
#   UKB SEROLOGY: COMPREHENSIVE ANTIGEN DIAGNOSTIC PLOTS
################################################################################
#  ∙ This script generates detailed diagnostic plots for each antigen to:
#    1. Validate the mixture model fits
#    2. Compare different phenotype definitions
#    3. Understand the distribution of weights and probabilities
#    4. Assess the quality of the phenotype generation process
################################################################################

# --- 1. SETUP & CONFIGURATION ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, mixsmsn, sn, here, glue, pheatmap, furrr, data.table, patchwork, cowplot)

# Setup parallel processing
plan(multisession, workers = availableCores() - 1)

# Ensure the output directory exists
if (!dir.exists("antigen_diagnostics")) dir.create("antigen_diagnostics")

# Define the core antigens (same as in main script)
core_antigens <- c(
    "cmv_pp150", "cmv_pp28", "cmv_pp52",
    "hsv1", "hsv2",
    "ebv_vca", "ebv_ebna1", "ebv_zebra", "ebv_ead",
    "hp_omp", "hp_urea", "hp_caga", "hp_vaca",
    "toxo_p22",
    "bkv_vp1", "jcv_vp1", "mcv_vp1",
    "hhv6_ie1a", "hhv6_ie1b", "hhv7_u14",
    "ct_pgp3", "ct_mompd", "ct_tarpf2"
)

# Antigen mapping (same as in main script)
antigen_map <- tribble(
  ~short_name, ~full_name, ~threshold,
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
  "bkv_vp1", "BK VP1 antigen for Human Polyomavirus BKV | Instance 0", 250,
  "jcv_vp1", "JC VP1 antigen for Human Polyomavirus JCV | Instance 0", 250,
  "mcv_vp1", "MC VP1 antigen for Merkel Cell Polyomavirus | Instance 0", 250,
  "ct_mompa", "momp A antigen for Chlamydia trachomatis | Instance 0", 100,
  "ct_mompd", "momp D antigen for Chlamydia trachomatis | Instance 0", 100,
  "ct_tarpf1", "tarp-D F1 antigen for Chlamydia trachomatis | Instance 0", 100,
  "ct_tarpf2", "tarp-D F2 antigen for Chlamydia trachomatis | Instance 0", 100,
  "ct_porb", "PorB antigen for Chlamydia trachomatis | Instance 0", 80,
  "ct_pgp3", "pGP3 antigen for Chlamydia trachomatis | Instance 0", 200,
  "hp_caga", "CagA antigen for Helicobacter pylori | Instance 0", 400,
  "hp_vaca", "VacA antigen for Helicobacter pylori | Instance 0", 100,
  "hp_omp", "OMP antigen for Helicobacter pylori | Instance 0", 170,
  "hp_groel", "GroEL antigen for Helicobacter pylori | Instance 0", 80,
  "hp_catalase", "Catalase antigen for Helicobacter pylori | Instance 0", 180,
  "hp_urea", "UreA antigen for Helicobacter pylori | Instance 0", 130,
  "toxo_p22", "p22 antigen for Toxoplasma gondii | Instance 0", 100,
  "toxo_sag1", "sag1 antigen for Toxoplasma gondii | Instance 0", 160,
  "kshv_k8_1", "K8.1 antigen for Kaposi's Sarcoma-Associated Herpesvirus | Instance 0", 175,
  "kshv_lana", "LANA antigen for Kaposi's Sarcoma-Associated Herpesvirus | Instance 0", 100,
  "hbv_hbc", "HBc antigen for Hepatitis B Virus | Instance 0", 100,
  "hbv_hbe", "HBe antigen for Hepatitis B Virus | Instance 0", 150,
  "hcv_core", "Core antigen for Hepatitis C Virus | Instance 0", 150,
  "hcv_ns3", "NS3 antigen for Hepatitis C Virus | Instance 0", 150,
  "hiv_env", "HIV-1 env antigen for Human Immunodeficiency Virus | Instance 0", 150,
  "hiv_gag", "HIV-1 gag antigen for Human Immunodeficiency Virus | Instance 0", 600,
  "htlv_env", "HTLV-1 env antigen for Human T-Lymphotropic Virus 1 | Instance 0", 150,
  "htlv_gag", "HTLV-1 gag antigen for Human T-Lymphotropic Virus 1 | Instance 0", 1500,
  "hpv16_l1", "L1 antigen for Human Papillomavirus type-16 | Instance 0", 175,
  "hpv16_e6", "E6 antigen for Human Papillomavirus type-16 | Instance 0", 120,
  "hpv16_e7", "E7 antigen for Human Papillomavirus type-16 | Instance 0", 150,
  "hpv18_l1", "L1 antigen for Human Papillomavirus type-18 | Instance 0", 175
) %>%
filter(short_name %in% core_antigens)

# --- 2. LOAD DATA ---
cat("Loading phenotype data...\n")

# Load the master phenotype table
if (file.exists("phenotypes_master.tsv")) {
  master_pheno_table <- read_tsv("phenotypes_master.tsv", show_col_types = FALSE)
} else {
  stop("Phenotype file not found. Please run 02_create_serology_phenotypes.R first.")
}

# Load mixture model fits
if (file.exists("mixture_model_fits.rds")) {
  all_model_fits <- readRDS("mixture_model_fits.rds")
  cat("Loaded mixture model fits for", length(all_model_fits), "antigens.\n")
} else {
  warning("Mixture model fits not found. Plot 1 will not show fitted densities.")
  all_model_fits <- NULL
}

# Load raw data for comparison
if (file.exists("serology_export_title.tsv")) {
  raw_data <- read_tsv("serology_export_title.tsv", show_col_types = FALSE, guess_max = 5000)
  cols_to_keep <- setNames(antigen_map$full_name, antigen_map$short_name)
  data_clean <- raw_data %>%
    select(FID = `Participant ID`, IID = `Participant ID`, any_of(cols_to_keep)) %>%
    rename(any_of(cols_to_keep)) %>%
    mutate(across(-c(FID, IID), as.numeric)) %>%
    filter(!if_all(-c(FID, IID), is.na))
} else {
  stop("Raw serology data not found.")
}

cat("Data loaded for", ncol(data_clean) - 2, "antigens and", nrow(data_clean), "samples.\n")

# Load the single, combined covariate file
covar_table <- NULL
if (file.exists("quickdraws_input/covariates_all.tsv")) {
  covar_table <- read_tsv("quickdraws_input/covariates_all.tsv", show_col_types = FALSE)
  cat("Loaded combined covariates from quickdraws_input/covariates_all.tsv\n")
} else {
  warning("Combined covariate file not found. PC and other covariate-based plots will be skipped.")
}

# --- 3. DIAGNOSTIC PLOTTING FUNCTIONS ---

#' Fit mixture model (same as in main script)
fit_skew_mix <- function(y) {
  y_log <- log(y[y > 0 & is.finite(y)])
  if (length(y_log) < 200) return(NULL)
  
  tryCatch({
    mixsmsn::smsn.mix(
      y_log, g = 2, family = "Skew.t",
      nu = 4, group = FALSE, calc.im = FALSE, obs.prob = TRUE
    )
  }, error = function(e) {
    warning("Failed to fit mixture model: ", e$message)
    return(NULL)
  })
}

#' Plot 1: Original distribution with hard cutoff and mixture model fits
plot_original_distribution <- function(antigen_name, raw_data, threshold, model_fits = NULL) {
  mfi_values <- raw_data[[antigen_name]]
  mfi_values <- mfi_values[mfi_values > 0 & is.finite(mfi_values)]
  
  p <- ggplot(tibble(mfi = mfi_values), aes(x = mfi)) +
    geom_histogram(aes(y = after_stat(density)), bins = 75, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 1) +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = glue("Original MFI Distribution: {antigen_name}"),
      subtitle = glue("Hard cutoff threshold = {threshold}"),
      x = "MFI (log scale)",
      y = "Density"
    ) +
    theme_light()
  
  # Add mixture model fits if available
  if (!is.null(model_fits) && antigen_name %in% names(model_fits)) {
    fit <- model_fits[[antigen_name]]
    pos_comp <- which.max(fit$mu)
    neg_comp <- setdiff(1:2, pos_comp)

    # Generate density curves using skew-normal (dsn) for plotting stability
    x_range <- range(log(mfi_values), na.rm = TRUE)
    x_points <- seq(x_range[1], x_range[2], length.out = 500)
    
    # Calculate component densities
    density_pos <- sn::dsn(x_points, xi = fit$mu[pos_comp], omega = sqrt(fit$sigma2[pos_comp]), 
                          alpha = fit$shape[pos_comp]) * fit$pii[pos_comp]
    density_neg <- sn::dsn(x_points, xi = fit$mu[neg_comp], omega = sqrt(fit$sigma2[neg_comp]), 
                          alpha = fit$shape[neg_comp]) * fit$pii[neg_comp]
    
    density_df <- tibble(
      x = exp(x_points),
      positive = density_pos,
      negative = density_neg,
      total = positive + negative
    )
    
    p <- p +
      geom_line(data = density_df, aes(x = x, y = positive), color = "red", linewidth = 1) +
      geom_line(data = density_df, aes(x = x, y = negative), color = "blue", linewidth = 1)
  }
  
  return(p)
}

#' Plot 2: Distribution colored by mixture model probability
plot_probability_colored_distribution <- function(antigen_name, pheno_data) {
  igg_col  <- glue("{antigen_name}_IgG_raw")
  soft_col <- glue("{antigen_name}_sero_soft")

  # Check required columns
  if (!all(c(igg_col, soft_col) %in% names(pheno_data))) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                              label = "Required columns missing for this antigen") +
             theme_void())
  }

  df <- pheno_data %>%
    select(igg = all_of(igg_col), prob = all_of(soft_col)) %>%
    filter(is.finite(igg) & !is.na(prob))

  p <- ggplot(df, aes(x = exp(igg), y = prob, colour = prob)) +
    geom_point(alpha = 0.4, size = 0.6, stroke = 0) +
    scale_colour_viridis_c(name = "Seropositive\nProbability", option = "plasma") +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(title = glue("IgG vs Seropositive Probability: {antigen_name}"),
         x = "MFI (log scale)", y = "Seropositive Probability") +
    theme_light()

  return(p)
}

#' Plot 3: Distribution of high-confidence seropositives
plot_high_confidence_seropositives <- function(antigen_name, pheno_data, prob_threshold = 0.96, mfi_cutoff = NA) {
  igg_col <- glue("{antigen_name}_IgG_raw")
  soft_col <- glue("{antigen_name}_sero_soft")
  
  if (!all(c(igg_col, soft_col) %in% names(pheno_data))) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Required columns missing for this antigen") + theme_void())
  }
  
  plot_data <- pheno_data %>%
    select(igg = all_of(igg_col), prob = all_of(soft_col)) %>%
    filter(!is.na(igg) & !is.na(prob))
  
  high_conf_pos <- plot_data %>% filter(prob > prob_threshold)
  high_conf_neg <- plot_data %>% filter(prob < (1 - prob_threshold))
  ambiguous     <- plot_data %>% filter(prob >= (1 - prob_threshold) & prob <= prob_threshold)
  
  p <- ggplot() +
    geom_histogram(data = high_conf_pos, aes(x = exp(igg), fill = "High-conf Seropositive"), 
                   bins = 50, alpha = 0.7) +
    geom_histogram(data = high_conf_neg, aes(x = exp(igg), fill = "High-conf Seronegative"), 
                   bins = 50, alpha = 0.7) +
    geom_histogram(data = ambiguous, aes(x = exp(igg), fill = "Ambiguous"), 
                   bins = 50, alpha = 0.5) +
    scale_fill_manual(values = c("High-conf Seropositive" = "red", 
                                 "High-conf Seronegative" = "blue", 
                                 "Ambiguous" = "grey")) +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = glue("High-Confidence Seropositive Distribution: {antigen_name}"),
      subtitle = glue("Prob threshold = {prob_threshold}"),
      x = "MFI (log scale)",
      y = "Count",
      fill = "Classification"
    ) +
    theme_light()
  
  # Optionally add MFI cutoff vertical line if provided
  if (is.finite(mfi_cutoff)) {
    p <- p + geom_vline(xintercept = mfi_cutoff, colour = "black", linetype = "dashed", linewidth = 1)
  }
  
  return(p)
}

#' Plot 4: Distribution of phenotypes by classification
plot_phenotype_distributions <- function(antigen_name, pheno_data) {
  igg_col <- glue("{antigen_name}_IgG_raw")
  soft_col <- glue("{antigen_name}_sero_soft")
  hard_col <- glue("{antigen_name}_sero_hard")
  
  if (!all(c(igg_col, soft_col, hard_col) %in% names(pheno_data))) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available") + theme_void())
  }
  
  plot_data <- pheno_data %>%
    select(igg = all_of(igg_col), soft = all_of(soft_col), hard = all_of(hard_col)) %>%
    filter(!is.na(igg) & !is.na(soft) & !is.na(hard)) %>%
    mutate(
      classification = case_when(
        hard == 1 & soft > 0.8 ~ "Hard+Soft Positive",
        hard == 0 & soft < 0.2 ~ "Hard+Soft Negative", 
        hard == 1 & soft <= 0.8 ~ "Hard Positive, Soft Low",
        hard == 0 & soft >= 0.2 ~ "Hard Negative, Soft High",
        TRUE ~ "Other"
      )
    )
  
  p <- ggplot(plot_data, aes(x = exp(igg), fill = classification)) +
    geom_histogram(bins = 75, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c(
      "Hard+Soft Positive" = "red",
      "Hard+Soft Negative" = "blue",
      "Hard Positive, Soft Low" = "orange",
      "Hard Negative, Soft High" = "purple",
      "Other" = "grey"
    )) +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = glue("Phenotype Distribution by Classification: {antigen_name}"),
      subtitle = "Comparing hard vs soft classifications",
      x = "MFI (log scale)",
      y = "Count",
      fill = "Classification"
    ) +
    theme_light()
  
  return(p)
}

#' Plot 5: Scatter plot of IgG vs probability
plot_igg_vs_probability <- function(antigen_name, pheno_data, threshold) {
  igg_col <- glue("{antigen_name}_IgG_raw")
  soft_col <- glue("{antigen_name}_sero_soft")
  
  if (!all(c(igg_col, soft_col) %in% names(pheno_data))) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available") + theme_void())
  }
  
  plot_data <- pheno_data %>%
    select(igg = all_of(igg_col), prob = all_of(soft_col)) %>%
    filter(!is.na(igg) & !is.na(prob))
  
  p <- ggplot(plot_data, aes(x = exp(igg), y = prob)) +
    geom_point(alpha = 0.6, size = 0.5) +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = 0.5, color = "blue", linetype = "dashed", linewidth = 1) +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = glue("IgG Level vs Seropositive Probability: {antigen_name}"),
      subtitle = glue("Red line = hard threshold, Blue line = 50% probability"),
      x = "MFI (log scale)",
      y = "Seropositive Probability"
    ) +
    theme_light()
  
  return(p)
}

#' Plot 6: Weight distributions
plot_weight_distributions <- function(antigen_name, pheno_data) {
  w_hard_col <- glue("{antigen_name}_w_hard")
  w_igg_col  <- glue("{antigen_name}_w_igg")
  igg_col    <- glue("{antigen_name}_IgG_raw")

  if (!all(c(w_hard_col, w_igg_col, igg_col) %in% names(pheno_data))) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available") + theme_void())
  }

  plot_data <- pheno_data %>%
    select(igg = all_of(igg_col), w_hard = all_of(w_hard_col), w_igg = all_of(w_igg_col)) %>%
    filter(!is.na(igg) & !is.na(w_hard) & !is.na(w_igg))

  # Density curves: unweighted vs weighted (IgG weights)
  density_unweighted <- ggplot(plot_data, aes(x = exp(igg))) +
    geom_density(color = "grey40", linewidth = 1) +
    labs(y = "Density", x = "MFI (log scale)") +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    theme_light()

  density_weighted <- ggplot(plot_data, aes(x = exp(igg), weight = w_igg)) +
    geom_density(color = "blue", linewidth = 1) +
    labs(y = "Weighted Density", x = "MFI (log scale)") +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    theme_light()

  # Histogram of weight values
  weight_hist <- ggplot(plot_data %>% pivot_longer(c(w_hard, w_igg), names_to = "type", values_to = "weight"),
                        aes(x = weight, fill = type)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c("w_hard" = "red", "w_igg" = "blue"),
                      labels = c("Hard Classification Weight", "IgG Analysis Weight")) +
    labs(x = "Weight Value", y = "Count", fill = "Weight Type") +
    theme_light()

  combined <- (density_unweighted + density_weighted) / weight_hist +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(title = glue("IgG Distributions & Weights: {antigen_name}"))

  return(combined)
}

#' Plot 7: Hard vs Soft classification comparison
plot_classification_comparison <- function(antigen_name, pheno_data, threshold) {
  soft_col <- glue("{antigen_name}_sero_soft")
  hard_col <- glue("{antigen_name}_sero_hard")
  
  if (!all(c(soft_col, hard_col) %in% names(pheno_data))) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available") + theme_void())
  }
  
  plot_data <- pheno_data %>%
    select(soft = all_of(soft_col), hard = all_of(hard_col)) %>%
    filter(!is.na(soft) & !is.na(hard))
  
  # Create confusion matrix style plot
  p <- ggplot(plot_data, aes(x = hard, y = soft)) +
    geom_jitter(alpha = 0.6, size = 0.5) +
    geom_hline(yintercept = 0.5, color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = 0.5, color = "blue", linetype = "dashed", linewidth = 1) +
    scale_x_continuous(breaks = c(0, 1), labels = c("Seronegative", "Seropositive")) +
    labs(
      title = glue("Hard vs Soft Classification: {antigen_name}"),
      subtitle = "Red line = 50% soft threshold, Blue line = hard threshold",
      x = "Hard Classification (0/1)",
      y = "Soft Probability (0-1)"
    ) +
    theme_light()
  
  return(p)
}

#' Plot 8: Summary statistics
plot_summary_stats <- function(antigen_name, pheno_data, raw_data, threshold) {
  igg_col <- glue("{antigen_name}_IgG_raw")
  soft_col <- glue("{antigen_name}_sero_soft")
  hard_col <- glue("{antigen_name}_sero_hard")
  
  if (!all(c(igg_col, soft_col, hard_col) %in% names(pheno_data))) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available") + theme_void())
  }
  
  # Calculate summary statistics
  stats_data <- pheno_data %>%
    select(igg = all_of(igg_col), soft = all_of(soft_col), hard = all_of(hard_col)) %>%
    filter(!is.na(igg) & !is.na(soft) & !is.na(hard))
  
  n_total <- nrow(stats_data)
  n_hard_pos <- sum(stats_data$hard == 1)
  n_soft_high <- sum(stats_data$soft > 0.9)
  n_soft_low <- sum(stats_data$soft < 0.1)
  n_ambiguous <- sum(stats_data$soft >= 0.1 & stats_data$soft <= 0.9)
  
  # Create summary text
  summary_text <- glue(
    "Total samples: {n_total}\n",
    "Hard seropositive: {n_hard_pos} ({round(n_hard_pos/n_total*100, 1)}%)\n",
    "Soft high-conf (>0.9): {n_soft_high} ({round(n_soft_high/n_total*100, 1)}%)\n",
    "Soft low-conf (<0.1): {n_soft_low} ({round(n_soft_low/n_total*100, 1)}%)\n",
    "Ambiguous (0.1-0.9): {n_ambiguous} ({round(n_ambiguous/n_total*100, 1)}%)"
  )
  
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = summary_text, 
             hjust = 0.5, vjust = 0.5, size = 4) +
    labs(title = glue("Summary Statistics: {antigen_name}")) +
    theme_void()
  
  return(p)
}

#' Plot 5: Scatter IgG vs PC1_assay_noise coloured by serostatus group
plot_pc1_vs_igg <- function(antigen_name, pheno_data, covar_data) {
  if (is.null(covar_data) || !"PC1_assay_noise" %in% names(covar_data)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "PC1_assay_noise data not available") + theme_void())
  }
  igg_col  <- glue("{antigen_name}_IgG_raw")
  soft_col <- glue("{antigen_name}_sero_soft")
  
  if (!all(c("IID", igg_col, soft_col) %in% names(pheno_data))) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Data not available") + theme_void())
  }
  
  merged <- pheno_data %>%
    select(IID, igg = all_of(igg_col), prob = all_of(soft_col)) %>%
    left_join(covar_data, by = "IID") %>%
    filter(!is.na(igg) & !is.na(prob) & !is.na(PC1_assay_noise)) %>%
    mutate(group = case_when(
      prob > 0.95 ~ "Likely Seropositive (>0.95)",
      prob < 0.05 ~ "Likely Seronegative (<0.05)",
      TRUE        ~ "Ambiguous (0.05–0.95)"
    ))
  
  p <- ggplot(merged, aes(x = exp(igg), y = PC1_assay_noise, colour = group)) +
    geom_point(alpha = 0.5, size = 0.6, stroke = 0) +
    scale_colour_manual(values = c(
      "Likely Seropositive (>0.95)" = "red",
      "Likely Seronegative (<0.05)" = "blue",
      "Ambiguous (0.05–0.95)" = "grey"
    )) +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = glue("IgG vs PC1_assay_noise: {antigen_name}"),
      x = "MFI (log scale)",
      y = "PC1_assay_noise",
      colour = "Group"
    ) +
    theme_light()
  return(p)
}

#' Plot 8: PC1_assay_noise vs hard/soft classifications with trend lines
plot_pc1_trends <- function(antigen_name, pheno_data, covar_data) {
  if (is.null(covar_data) || !all(c("PC1_assay_noise", "PC1_seropositive") %in% names(covar_data))) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "PC1 components not available") + theme_void())
  }
  
  soft_col <- glue("{antigen_name}_sero_soft")
  igg_col  <- glue("{antigen_name}_IgG_raw")
  
  if (!all(c("IID", soft_col, igg_col) %in% names(pheno_data))) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Required columns missing for this antigen") + theme_void())
  }
  
  # Check if we have both PC1 components
  if (!all(c("PC1_assay_noise", "PC1_seropositive") %in% names(covar_data))) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "PC1 components not available") + theme_void())
  }
  
  merged <- pheno_data %>%
    select(IID, soft = all_of(soft_col), igg = all_of(igg_col)) %>%
    left_join(covar_data, by = "IID") %>%
    filter(!is.na(soft) & !is.na(igg) & !is.na(PC1_assay_noise) & !is.na(PC1_seropositive))
  
  # Split into high and low probability groups
  high_prob <- merged %>% filter(soft > 0.95)
  low_prob  <- merged %>% filter(soft < 0.05)
  
  # Create four facets: PC1_assay_noise and PC1_seropositive for each group
  plot_data <- bind_rows(
    low_prob %>% mutate(group = "Low Prob (<0.05)", pc_type = "PC1_assay_noise", pc_value = PC1_assay_noise),
    low_prob %>% mutate(group = "Low Prob (<0.05)", pc_type = "PC1_seropositive", pc_value = PC1_seropositive),
    high_prob %>% mutate(group = "High Prob (>0.95)", pc_type = "PC1_assay_noise", pc_value = PC1_assay_noise),
    high_prob %>% mutate(group = "High Prob (>0.95)", pc_type = "PC1_seropositive", pc_value = PC1_seropositive)
  )
  
  p <- ggplot(plot_data, aes(x = exp(igg), y = pc_value)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
    facet_grid(pc_type ~ group, scales = "free") +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = glue("PC1 Components vs IgG Level: {antigen_name}"),
      subtitle = "Top: PC1_assay_noise (seroneg), Bottom: PC1_seropositive (seropos)",
      x = "MFI (log scale)",
      y = "PC1 Value"
    ) +
    theme_light()
  
  return(p)
}

#' Plot 6: Histogram of seropositive probabilities
plot_probability_histogram <- function(antigen_name, pheno_data) {
  soft_col <- glue("{antigen_name}_sero_soft")
  
  if (!all(c(soft_col) %in% names(pheno_data))) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Required columns missing for this antigen") + theme_void())
  }
  
  plot_data <- pheno_data %>%
    select(prob = all_of(soft_col)) %>%
    filter(!is.na(prob))
  
  p <- ggplot(plot_data, aes(x = prob)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = glue("Distribution of Seropositive Probabilities: {antigen_name}"),
      subtitle = "Red line = 50% threshold",
      x = "Seropositive Probability",
      y = "Count"
    ) +
    theme_light()
  
  return(p)
}

# --- Adjusted MAIN PLOTTING LOOP (single comprehensive PDF) ---
cat("\n--- Generating Comprehensive Diagnostic Plots (single PDF with all antigens) ---\n")

# Create a list to store all plots
all_plots <- list()

for (antigen in core_antigens) {
  cat("Processing:", antigen, "\n")
  
  threshold <- antigen_map$threshold[antigen_map$short_name == antigen]
  
  # Generate base stats for subtitle
  soft_col <- glue("{antigen}_sero_soft")
  if (soft_col %in% names(master_pheno_table)) {
    probs <- master_pheno_table[[soft_col]]
    n_total <- sum(!is.na(probs))
    n_pos   <- sum(probs > 0.95, na.rm = TRUE)
    n_neg   <- sum(probs < 0.05, na.rm = TRUE)
    n_ambig <- n_total - n_pos - n_neg
    subtitle_text <- glue("N={n_total} | >0.95={n_pos} | <0.05={n_neg} | Ambig={n_ambig}")
  } else {
    subtitle_text <- "Sample counts unavailable"
  }
  
  # Generate plots
  p1 <- plot_original_distribution(antigen, data_clean, threshold, all_model_fits)
  p2 <- plot_probability_colored_distribution(antigen, master_pheno_table)
  p3 <- plot_high_confidence_seropositives(antigen, master_pheno_table, 0.96, threshold)
  p4 <- plot_phenotype_distributions(antigen, master_pheno_table)
  p5 <- plot_igg_vs_probability(antigen, master_pheno_table, threshold)
  p6 <- plot_probability_histogram(antigen, master_pheno_table)
  p7 <- plot_classification_comparison(antigen, master_pheno_table, threshold)
  p8 <- plot_pc1_trends(antigen, master_pheno_table, covar_table)

  combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6) / (p7 + p8) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title   = glue("Diagnostics: {antigen}"),
      subtitle = subtitle_text,
      theme   = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  # Store the combined plot
  all_plots[[antigen]] <- combined_plot
  cat("  -> Processed", antigen, "\n")
}

# Save as single comprehensive PDF with one antigen per page
cat("\n--- Creating comprehensive PDF with one antigen per page ---\n")

# Use pdf() device to create multi-page PDF
pdf("antigen_diagnostics/comprehensive_antigen_diagnostics.pdf", 
    width = 16, height = 18, onefile = TRUE)

for (antigen in names(all_plots)) {
  print(all_plots[[antigen]])
}

dev.off()

cat("  -> Saved comprehensive PDF: antigen_diagnostics/comprehensive_antigen_diagnostics.pdf\n")
cat("  -> Total antigens:", length(all_plots), "\n")
cat("  -> One antigen per page\n")

# --- 5. CREATE SUMMARY REPORT ---
cat("\n--- Creating Summary Report ---\n")

# Generate summary statistics for all antigens
summary_stats <- map_dfr(core_antigens, function(antigen) {
  igg_col <- glue("{antigen}_IgG_raw")
  soft_col <- glue("{antigen}_sero_soft")
  hard_col <- glue("{antigen}_sero_hard")
  
  if (!all(c(igg_col, soft_col, hard_col) %in% names(master_pheno_table))) {
    return(tibble(antigen = antigen, status = "Data not available"))
  }
  
  stats_data <- master_pheno_table %>%
    select(igg = all_of(igg_col), soft = all_of(soft_col), hard = all_of(hard_col)) %>%
    filter(!is.na(igg) & !is.na(soft) & !is.na(hard))
  
  n_total <- nrow(stats_data)
  n_hard_pos <- sum(stats_data$hard == 1)
  n_soft_high <- sum(stats_data$soft > 0.9)
  n_soft_low <- sum(stats_data$soft < 0.1)
  n_ambiguous <- sum(stats_data$soft >= 0.1 & stats_data$soft <= 0.9)
  
  tibble(
    antigen = antigen,
    n_total = n_total,
    n_hard_seropositive = n_hard_pos,
    pct_hard_seropositive = round(n_hard_pos/n_total*100, 1),
    n_soft_high_conf = n_soft_high,
    pct_soft_high_conf = round(n_soft_high/n_total*100, 1),
    n_soft_low_conf = n_soft_low,
    pct_soft_low_conf = round(n_soft_low/n_total*100, 1),
    n_ambiguous = n_ambiguous,
    pct_ambiguous = round(n_ambiguous/n_total*100, 1)
  )
})

write_tsv(summary_stats, "antigen_diagnostics/summary_statistics.tsv")
cat("  -> Saved summary statistics: antigen_diagnostics/summary_statistics.tsv\n")

# Create a summary plot
summary_plot <- summary_stats %>%
  filter(!is.na(n_total)) %>%
  select(antigen, pct_hard_seropositive, pct_soft_high_conf, pct_ambiguous) %>%
  pivot_longer(-antigen, names_to = "metric", values_to = "percentage") %>%
  mutate(metric = case_when(
    metric == "pct_hard_seropositive" ~ "Hard Seropositive",
    metric == "pct_soft_high_conf" ~ "Soft High Conf",
    metric == "pct_ambiguous" ~ "Ambiguous"
  )) %>%
  ggplot(aes(x = reorder(antigen, percentage), y = percentage, fill = metric)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Hard Seropositive" = "red", 
                              "Soft High Conf" = "blue", 
                              "Ambiguous" = "grey")) +
  coord_flip() +
  labs(
    title = "Summary of Classification Results Across All Antigens",
    x = "Antigen",
    y = "Percentage (%)",
    fill = "Classification Type"
  ) +
  theme_light() +
  theme(axis.text.y = element_text(size = 8))

ggsave("antigen_diagnostics/summary_classification_plot.pdf", summary_plot, 
       width = 12, height = 10)
cat("  -> Saved summary classification plot: antigen_diagnostics/summary_classification_plot.pdf\n")

cat("\n--- Diagnostic Analysis Complete ---\n")
cat("All plots saved in: antigen_diagnostics/\n") 