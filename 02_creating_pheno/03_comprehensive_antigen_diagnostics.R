################################################################################
#
#   UKB SEROLOGY: COMPREHENSIVE ANTIGEN & PATHOGEN DIAGNOSTIC PLOTS
#   VERSION: 9.0 (Final - With PDF size reduction via point shape)
#
#   ▸ This script is compatible with the output of '01_final_master_sero_phenotype.R'.
#   ▸ It preserves all original plots, adapting them to the new data columns.
#   ▸ It uses `shape = "."` in dense scatter plots to reduce PDF file size
#     without requiring external dependencies like Cairo.
#   ▸ It includes the requested pathogen-level diagnostics at the end.
#
################################################################################

# --- 1. SETUP & CONFIGURATION ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, mixsmsn, sn, here, glue, pheatmap, furrr, data.table, patchwork, cowplot, scales)

# Setup parallel processing
plan(multisession, workers = availableCores() - 1)

# Ensure output directories exist
if (!dir.exists("antigen_diagnostics")) dir.create("antigen_diagnostics")
if (!dir.exists("pathogen_diagnostics")) dir.create("pathogen_diagnostics")

# --- GLOBAL DEFINITIONS (from user's script) ---
core_antigens <- c(
    "cmv_pp150", "cmv_pp28", "cmv_pp52",
    "hsv1", "hsv2",
    "ebv_vca", "ebv_ebna1", "ebv_zebra", "ebv_ead",
    "hp_omp", "hp_urea", "hp_caga", "hp_vaca", "hp_groel", "hp_catalase",
    "toxo_p22",
    "bkv_vp1", "jcv_vp1", "mcv_vp1",
    "hhv6_ie1a", "hhv6_ie1b", "hhv7_u14",
    "ct_pgp3", "ct_mompd", "ct_tarpf2", "ct_mompa", "ct_tarpf1",
    "kshv_lana"
)
pathogens <- c("ebv", "cmv", "ct", "hp") # For new section

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
  "kshv_lana", "LANA antigen for Kaposi's Sarcoma-Associated Herpesvirus | Instance 0", 100
) %>% filter(short_name %in% core_antigens)

# --- 2. DATA LOADING & RECONSTRUCTION ---
cat("--- Loading and Reconstructing Data from Script 01 Outputs ---\n")

if (!dir.exists("quickdraws_input") || !file.exists("quickdraws_input/analysis_manifest.tsv")) {
  stop("Outputs from '01_final_master_sero_phenotype.R' not found. Please run it first.")
}
phenos_quant <- read_tsv("quickdraws_input/phenotypes_quantitative.tsv", show_col_types = FALSE)
phenos_binary <- read_tsv("quickdraws_input/phenotypes_binary.tsv", show_col_types = FALSE)
manifest <- read_tsv("quickdraws_input/analysis_manifest.tsv", show_col_types = FALSE)

master_pheno_table <- full_join(phenos_quant, phenos_binary, by = c("FID", "IID"))
cat("Loaded and merged quantitative and binary phenotype files.\n")

weights_files <- manifest %>% filter(!is.na(weights_file)) %>% distinct(weights_file) %>% pull(weights_file)
p_soft_data_list <- map(weights_files, function(w_file) {
  antigen_name <- str_remove(w_file, "_p_soft\\.weights")
  new_col_name <- glue("{antigen_name}_p_soft")
  read_tsv(file.path("quickdraws_input", w_file), show_col_types = FALSE) %>%
    rename(!!new_col_name := Weight)
})
master_pheno_table <- reduce(p_soft_data_list, full_join, by = c("FID", "IID")) %>%
  right_join(master_pheno_table, by = c("FID", "IID"))
cat("Reconstructed full phenotype table with", ncol(master_pheno_table), "columns.\n")

if (file.exists("mixture_model_fits.rds")) {
  all_model_fits <- readRDS("mixture_model_fits.rds")
  cat("Loaded mixture model fits for", length(all_model_fits), "antigens.\n")
} else {
  warning("Mixture model fits not found. Plot 1 will not show fitted densities.")
  all_model_fits <- NULL
}

latent_factors_data <- NULL
if (all(c("latent_factor_1", "latent_factor_2", "latent_factor_igg") %in% names(master_pheno_table))) {
  latent_factors_data <- master_pheno_table %>% 
    select(FID, IID, latent_factor_1, latent_factor_2, latent_factor_igg)
  cat("Extracted latent factors from phenotype data for plotting\n")
} else {
  warning("Latent factors not found in phenotype data. PC-based plots will be skipped.")
}

# --- 3. DIAGNOSTIC PLOTTING FUNCTIONS (ADAPTED) ---

#' Plot 1: Original distribution with hard cutoff and mixture model fits
plot_original_distribution <- function(antigen_name, threshold, model_fits) {
  igg_col <- glue("{antigen_name}_IgG_raw")
  if(!igg_col %in% names(master_pheno_table)) return(ggplot() + annotate("text", x=0.5, y=0.5, label="Data missing") + theme_void())
  
  df <- master_pheno_table %>% select(igg_raw = all_of(igg_col)) %>% drop_na()
  
  p <- ggplot(df, aes(x = exp(igg_raw))) +
    geom_histogram(aes(y = after_stat(density)), bins = 75, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 1) +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = glue("MFI Distribution: {antigen_name}"),
      subtitle = glue("BL Hard Threshold = {threshold}"),
      x = "MFI (log scale)", y = "Density"
    ) + theme_light()

  if (!is.null(model_fits) && antigen_name %in% names(model_fits) && requireNamespace("sn", quietly=TRUE)) {
    fit <- model_fits[[antigen_name]]
    pos_comp <- which.max(fit$mu); neg_comp <- which.min(fit$mu)
    x_range <- range(df$igg_raw, na.rm = TRUE); x_points <- seq(x_range[1], x_range[2], length.out = 500)
    density_pos <- sn::dsn(x_points, xi=fit$mu[pos_comp], omega=sqrt(fit$sigma2[pos_comp]), alpha=fit$shape[pos_comp]) * fit$pii[pos_comp]
    density_neg <- sn::dsn(x_points, xi=fit$mu[neg_comp], omega=sqrt(fit$sigma2[neg_comp]), alpha=fit$shape[neg_comp]) * fit$pii[neg_comp]
    density_df <- tibble(x = exp(x_points), positive = density_pos, negative = density_neg)
    p <- p +
      geom_line(data = density_df, aes(x=x, y=positive), color="red", linewidth=1) +
      geom_line(data = density_df, aes(x=x, y=negative), color="blue", linewidth=1)
  }
  return(p)
}

#' Plot 2: Distribution colored by mixture model probability
plot_probability_colored_distribution <- function(antigen_name, pheno_data) {
  igg_col <- glue("{antigen_name}_IgG_raw"); soft_col <- glue("{antigen_name}_p_soft")
  if (!all(c(igg_col, soft_col) %in% names(pheno_data))) return(ggplot() + annotate("text",x=0.5,y=0.5,label="Data missing") + theme_void())
  df <- pheno_data %>% select(igg = all_of(igg_col), prob = all_of(soft_col)) %>% filter(is.finite(igg) & !is.na(prob))
  ggplot(df, aes(x = exp(igg), y = prob, colour = prob)) +
    geom_point(alpha = 0.5, shape = ".") +
    scale_colour_viridis_c(name = "Seropositive\nProbability", option = "plasma") +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(title = glue("IgG vs Seropositive Probability: {antigen_name}"), x = "MFI (log scale)", y = "Seropositive Probability (p_soft)") +
    theme_light()
}

#' Plot 3: Distribution of high-confidence seropositives
plot_high_confidence_seropositives <- function(antigen_name, pheno_data, prob_threshold = 0.96) {
  igg_col <- glue("{antigen_name}_IgG_raw"); soft_col <- glue("{antigen_name}_p_soft")
  if (!all(c(igg_col, soft_col) %in% names(pheno_data))) return(ggplot() + annotate("text",x=0.5,y=0.5,label="Data missing") + theme_void())
  
  plot_data <- pheno_data %>% select(igg = all_of(igg_col), prob = all_of(soft_col)) %>% filter(!is.na(igg) & !is.na(prob))
  high_conf_pos <- plot_data %>% filter(prob > prob_threshold)
  high_conf_neg <- plot_data %>% filter(prob < (1 - prob_threshold))
  ambiguous <- plot_data %>% filter(prob >= (1 - prob_threshold) & prob <= prob_threshold)
  
  ggplot() +
    geom_histogram(data = high_conf_pos, aes(x = exp(igg), fill = "High-conf Positive"), bins = 50, alpha = 0.7) +
    geom_histogram(data = high_conf_neg, aes(x = exp(igg), fill = "High-conf Negative"), bins = 50, alpha = 0.7) +
    geom_histogram(data = ambiguous, aes(x = exp(igg), fill = "Ambiguous"), bins = 50, alpha = 0.5) +
    scale_fill_manual(values = c("High-conf Positive" = "red", "High-conf Negative" = "blue", "Ambiguous" = "grey")) +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = glue("Distribution by Confidence: {antigen_name}"), subtitle = glue("Prob threshold = {prob_threshold}"),
      x = "MFI (log scale)", y = "Count", fill = "Classification (from p_soft)"
    ) + theme_light()
}

#' Plot 4: Distribution of phenotypes by classification
plot_phenotype_distributions <- function(antigen_name, pheno_data) {
  igg_col <- glue("{antigen_name}_IgG_raw"); bl_col <- glue("{antigen_name}_sero_hard_BL"); mix_col <- glue("{antigen_name}_sero_hard_mix")
  if (!all(c(igg_col, bl_col, mix_col) %in% names(pheno_data))) return(ggplot() + annotate("text",x=0.5,y=0.5,label="Data missing") + theme_void())
  
  plot_data <- pheno_data %>%
    select(igg = all_of(igg_col), bl = all_of(bl_col), mix = all_of(mix_col)) %>%
    filter(!is.na(igg) & !is.na(bl) & !is.na(mix)) %>%
    mutate(
      classification = case_when(
        bl == 1 & mix == 1 ~ "Concordant Positive",
        bl == 0 & mix == 0 ~ "Concordant Negative",
        bl == 1 & mix == 0 ~ "Discordant (BL+, Mix-)",
        bl == 0 & mix == 1 ~ "Discordant (BL-, Mix+)",
        TRUE ~ "Other"
      )
    )
  
  ggplot(plot_data, aes(x = exp(igg), fill = classification)) +
    geom_histogram(bins = 75, alpha = 0.8, position = "identity") +
    scale_fill_manual(values = c("Concordant Positive" = "red", "Concordant Negative" = "blue", "Discordant (BL+, Mix-)" = "orange", "Discordant (BL-, Mix+)" = "purple")) +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = glue("Phenotype Distribution by Classification: {antigen_name}"), subtitle = "Comparing BL vs. Mixture-Model Hard Calls",
      x = "MFI (log scale)", y = "Count", fill = "Classification"
    ) + theme_light()
}

#' Plot 5: Scatter plot of IgG vs probability
plot_igg_vs_probability <- function(antigen_name, pheno_data, threshold) {
  igg_col <- glue("{antigen_name}_IgG_raw"); soft_col <- glue("{antigen_name}_p_soft")
  if (!all(c(igg_col, soft_col) %in% names(pheno_data))) return(ggplot() + annotate("text",x=0.5,y=0.5,label="Data missing") + theme_void())
  
  plot_data <- pheno_data %>% select(igg = all_of(igg_col), prob = all_of(soft_col)) %>% filter(!is.na(igg) & !is.na(prob))
  
  ggplot(plot_data, aes(x = exp(igg), y = prob)) +
    geom_point(alpha = 0.6, shape = ".") +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = 0.5, color = "blue", linetype = "dashed", linewidth = 1) +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = glue("IgG Level vs. Seropositive Probability: {antigen_name}"), subtitle = "Red line = BL threshold, Blue line = 50% probability",
      x = "MFI (log scale)", y = "Seropositive Probability (p_soft)"
    ) + theme_light()
}

#' Plot 6: Weight distributions (Placeholder)
plot_weight_distributions <- function(antigen_name, pheno_data) {
    ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Plot not generated: Weight columns (e.g., _w_hard, _w_igg)\n are not created by the current version of script 01.", size = 3) +
        labs(title = glue("Weight Distributions (Data Unavailable): {antigen_name}")) +
        theme_void()
}

#' Plot 7: Hard vs Soft classification comparison
plot_classification_comparison <- function(antigen_name, pheno_data, threshold) {
  soft_col <- glue("{antigen_name}_p_soft"); bl_col <- glue("{antigen_name}_sero_hard_BL"); mix_col <- glue("{antigen_name}_sero_hard_mix")
  if (!all(c(soft_col, bl_col, mix_col) %in% names(pheno_data))) return(ggplot() + annotate("text",x=0.5,y=0.5,label="Data missing") + theme_void())

  plot_data <- pheno_data %>%
    select(p_soft = all_of(soft_col), bl = all_of(bl_col), mix = all_of(mix_col)) %>%
    filter(!is.na(p_soft) & !is.na(bl) & !is.na(mix))

  ggplot(plot_data, aes(y = p_soft)) +
    geom_jitter(aes(x = bl, color = "BL Hard Call"), width = 0.2, alpha = 0.3, shape = ".") +
    geom_jitter(aes(x = mix + 0.1, color = "Mix Hard Call"), width = 0.2, alpha = 0.3, shape = ".") + # offset for visibility
    scale_color_manual(values = c("BL Hard Call" = "red", "Mix Hard Call" = "blue")) +
    scale_x_continuous(breaks = c(0, 1), labels = c("Seronegative", "Seropositive")) +
    labs(
      title = glue("Hard Calls vs. Soft Probability: {antigen_name}"), subtitle = "Comparing BL and Mixture-Model hard calls against p_soft",
      x = "Hard Classification (0/1)", y = "Soft Probability (p_soft)", color = "Hard Call Type"
    ) + theme_light()
}

#' Plot 8: Summary statistics
plot_summary_stats <- function(antigen_name, pheno_data) {
  p_soft_col <- glue("{antigen_name}_p_soft"); bl_col <- glue("{antigen_name}_sero_hard_BL"); mix_col <- glue("{antigen_name}_sero_hard_mix")
  if (!all(c(p_soft_col, bl_col, mix_col) %in% names(pheno_data))) return(ggplot() + annotate("text",x=0.5,y=0.5,label="Data missing") + theme_void())
  
  stats_data <- pheno_data %>% select(p_soft = all_of(p_soft_col), bl = all_of(bl_col), mix = all_of(mix_col)) %>% drop_na()
  
  n_total <- nrow(stats_data)
  n_bl_pos <- sum(stats_data$bl == 1); pct_bl_pos <- round(n_bl_pos/n_total*100, 1)
  n_mix_pos <- sum(stats_data$mix == 1); pct_mix_pos <- round(n_mix_pos/n_total*100, 1)
  n_soft_high <- sum(stats_data$p_soft > 0.95); pct_soft_high <- round(n_soft_high/n_total*100, 1)
  n_ambiguous <- sum(stats_data$p_soft >= 0.05 & stats_data$p_soft <= 0.95); pct_ambiguous <- round(n_ambiguous/n_total*100, 1)

  summary_text <- glue(
    "Total samples: {n_total}\n\n",
    "BL Hard Seropositive: {n_bl_pos} ({pct_bl_pos}%)\n",
    "Mix Hard Seropositive: {n_mix_pos} ({pct_mix_pos}%)\n\n",
    "p_soft > 0.95: {n_soft_high} ({pct_soft_high}%)\n",
    "p_soft 0.05-0.95 (Ambiguous): {n_ambiguous} ({pct_ambiguous}%)"
  )
  ggplot() + annotate("text", x=0.5,y=0.5,label=summary_text,hjust=0.5,vjust=0.5,size=4) + labs(title=glue("Summary Statistics: {antigen_name}")) + theme_void()
}

#' Plot 9: Scatter IgG vs latent_factor_1 coloured by serostatus group
plot_pc1_vs_igg <- function(antigen_name, pheno_data, latent_data) {
  if (is.null(latent_data)) return(ggplot() + annotate("text",x=0.5,y=0.5,label="Latent factor data missing") + theme_void())
  igg_col <- glue("{antigen_name}_IgG_raw"); soft_col <- glue("{antigen_name}_p_soft")
  if (!all(c("IID", igg_col, soft_col) %in% names(pheno_data))) return(ggplot() + annotate("text",x=0.5,y=0.5,label="Data missing") + theme_void())

  merged <- pheno_data %>% select(IID, igg = all_of(igg_col), prob = all_of(soft_col)) %>%
    left_join(latent_data, by = "IID") %>% filter(!is.na(igg) & !is.na(prob) & !is.na(latent_factor_1)) %>%
    mutate(group = case_when(
      prob > 0.95 ~ "Likely Seropositive (>0.95)", prob < 0.05 ~ "Likely Seronegative (<0.05)", TRUE ~ "Ambiguous (0.05–0.95)"
    ))
  ggplot(merged, aes(x = exp(igg), y = latent_factor_1, colour = group)) +
    geom_point(alpha = 0.5, shape = ".") +
    scale_colour_manual(values = c("Likely Seropositive (>0.95)" = "red", "Likely Seronegative (<0.05)" = "blue", "Ambiguous (0.05–0.95)" = "grey")) +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(title = glue("IgG vs latent_factor_1: {antigen_name}"), x = "MFI (log scale)", y = "latent_factor_1", colour = "Group") +
    theme_light()
}

#' Plot 10: latent_factor_1 vs hard/soft classifications with trend lines
plot_pc1_trends <- function(antigen_name, pheno_data, latent_data) {
  if (is.null(latent_data)) return(ggplot() + annotate("text",x=0.5,y=0.5,label="Latent factor data missing") + theme_void())
  soft_col <- glue("{antigen_name}_p_soft"); igg_col <- glue("{antigen_name}_IgG_raw")
  if (!all(c("IID", soft_col, igg_col) %in% names(pheno_data))) return(ggplot() + annotate("text",x=0.5,y=0.5,label="Data missing") + theme_void())
  
  merged <- pheno_data %>% select(IID, soft = all_of(soft_col), igg = all_of(igg_col)) %>%
    left_join(latent_data, by = "IID") %>% filter(!is.na(soft) & !is.na(igg) & !is.na(latent_factor_1) & !is.na(latent_factor_2))
  
  high_prob <- merged %>% filter(soft > 0.95); low_prob <- merged %>% filter(soft < 0.05)
  
  plot_data <- bind_rows(
    low_prob %>% mutate(group = "Low Prob (<0.05)", pc_type = "latent_factor_1", pc_value = latent_factor_1),
    low_prob %>% mutate(group = "Low Prob (<0.05)", pc_type = "latent_factor_2", pc_value = latent_factor_2),
    high_prob %>% mutate(group = "High Prob (>0.95)", pc_type = "latent_factor_1", pc_value = latent_factor_1),
    high_prob %>% mutate(group = "High Prob (>0.95)", pc_type = "latent_factor_2", pc_value = latent_factor_2)
  )
  
  ggplot(plot_data, aes(x = exp(igg), y = pc_value)) +
    geom_point(alpha = 0.3, shape = ".") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
    facet_grid(pc_type ~ group, scales = "free") +
    scale_x_log10(labels = scales::label_number(accuracy = 1)) +
    labs(
      title = glue("Latent Factor Components vs IgG Level: {antigen_name}"), subtitle = "Top: latent_factor_1 (seroneg), Bottom: latent_factor_2 (seropos)",
      x = "MFI (log scale)", y = "Latent Factor Value"
    ) + theme_light()
}

#' Plot 11: Histogram of seropositive probabilities
plot_probability_histogram <- function(antigen_name, pheno_data) {
  soft_col <- glue("{antigen_name}_p_soft")
  if (!soft_col %in% names(pheno_data)) return(ggplot() + annotate("text",x=0.5,y=0.5,label="Data missing") + theme_void())
  
  plot_data <- pheno_data %>% select(prob = all_of(soft_col)) %>% filter(!is.na(prob))
  
  ggplot(plot_data, aes(x = prob)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", linewidth = 1) +
    labs(
      title = glue("Distribution of Seropositive Probabilities: {antigen_name}"), subtitle = "Red line = 50% threshold",
      x = "Seropositive Probability (p_soft)", y = "Count"
    ) + theme_light()
}

# --- 4. MAIN ANTIGEN PLOTTING LOOP (from user's script) ---
cat("\n--- Generating Comprehensive Diagnostic Plots (single PDF with all antigens) ---\n")

# Use a temporary directory for individual plots to avoid issues with parallel saving
temp_plot_dir <- "antigen_diagnostics/temp_plots"
if (!dir.exists(temp_plot_dir)) dir.create(temp_plot_dir, recursive = TRUE)

future_walk(core_antigens, function(antigen) {
  cat("Processing:", antigen, "\n")
  if (!glue("{antigen}_p_soft") %in% names(master_pheno_table)) {
      cat("  -> Skipping", antigen, "- essential data columns not found.\n")
      return(NULL)
  }
  
  threshold <- antigen_map$threshold[antigen_map$short_name == antigen]
  soft_col <- glue("{antigen}_p_soft")
  probs <- master_pheno_table[[soft_col]]; n_total <- sum(!is.na(probs)); n_pos <- sum(probs > 0.95, na.rm=T); n_neg <- sum(probs < 0.05, na.rm=T); n_ambig <- n_total - n_pos - n_neg
  subtitle_text <- glue("N={n_total} | >0.95={n_pos} | <0.05={n_neg} | Ambig={n_ambig}")

  # The plot list from the original user script is preserved.
  p1 <- plot_original_distribution(antigen, threshold, all_model_fits)
  p2 <- plot_probability_colored_distribution(antigen, master_pheno_table)
  p3 <- plot_high_confidence_seropositives(antigen, master_pheno_table, 0.96)
  p4 <- plot_phenotype_distributions(antigen, master_pheno_table)
  p5 <- plot_igg_vs_probability(antigen, master_pheno_table, threshold)
  p6 <- plot_probability_histogram(antigen, master_pheno_table)
  p7 <- plot_classification_comparison(antigen, master_pheno_table, threshold)
  p8 <- plot_pc1_trends(antigen, master_pheno_table, latent_factors_data)
  
  combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6) / (p7 + p8) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = glue("Comprehensive Diagnostics: {antigen}"),
      subtitle = subtitle_text,
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  # Save individual plot object to temporary file
  saveRDS(combined_plot, file.path(temp_plot_dir, glue("plot_{antigen}.rds")))
  
  cat("  -> Processed", antigen, "\n")
}, .options = furrr_options(seed = TRUE))


# Assemble the final PDF from the saved plot objects
cat("\nAssembling final PDF from generated plots...\n")
plot_files <- list.files(temp_plot_dir, pattern = "\\.rds$", full.names = TRUE)
# Ensure correct order
ordered_antigens <- map_chr(plot_files, ~str_extract(basename(.x), "(?<=plot_).*(?=\\.rds$)"))
plot_files <- plot_files[order(match(ordered_antigens, core_antigens))]

all_plots <- map(plot_files, readRDS)

pdf("antigen_diagnostics/comprehensive_antigen_diagnostics.pdf", width = 18, height = 22, onefile = TRUE)
for (p in all_plots) { print(p) }
dev.off()

# Clean up temporary directory
unlink(temp_plot_dir, recursive = TRUE)
cat("  -> Saved comprehensive antigen PDF\n")


# --- 5. CREATE SUMMARY REPORT (ADAPTED) ---
cat("\n--- Creating Summary Report ---\n")

summary_stats <- map_dfr(core_antigens, function(antigen) {
  p_soft_col <- glue("{antigen}_p_soft"); bl_col <- glue("{antigen}_sero_hard_BL"); mix_col <- glue("{antigen}_sero_hard_mix")
  if (!all(c(p_soft_col, bl_col, mix_col) %in% names(master_pheno_table))) return(NULL)
  
  stats_data <- master_pheno_table %>% select(p_soft = all_of(p_soft_col), bl = all_of(bl_col), mix = all_of(mix_col)) %>% drop_na()
  n_total <- nrow(stats_data)
  
  tibble(
    antigen = antigen,
    n_total = n_total,
    n_bl_pos = sum(stats_data$bl == 1), pct_bl_pos = round(n_bl_pos/n_total*100, 1),
    n_mix_pos = sum(stats_data$mix == 1), pct_mix_pos = round(n_mix_pos/n_total*100, 1),
    n_soft_high = sum(stats_data$p_soft > 0.95), pct_soft_high = round(n_soft_high/n_total*100, 1),
    n_ambiguous = sum(stats_data$p_soft >= 0.05 & stats_data$p_soft <= 0.95), pct_ambiguous = round(n_ambiguous/n_total*100, 1)
  )
})
write_tsv(summary_stats, "antigen_diagnostics/summary_statistics.tsv")
cat("  -> Saved summary statistics.\n")


# --- 6. PATHOGEN-LEVEL DIAGNOSTIC ANALYSIS (NEW SECTION) ---
cat("\n--- Generating Pathogen-Level Diagnostic Report (0-1 Scale) ---\n")

pdf_report_path <- "pathogen_diagnostics/pathogen_diagnostics_summary_report.pdf"
pdf(pdf_report_path, width = 10, height = 8, onefile = TRUE)

for (pth in pathogens) {
  cat("  -> Adding pathogen to report:", toupper(pth), "\n")
  p_soft_std_col <- glue("{pth}_psoft_std"); hard_bl_col <- glue("{pth}_sero_hard_BL"); hard_mix_col <- glue("{pth}_sero_hard_mix")
  pheno_cols <- c(p_soft_std_col, hard_bl_col, hard_mix_col)

  if (!all(pheno_cols %in% names(master_pheno_table))) {
    cat("     ! Skipping pathogen", toupper(pth), "- columns not found.\n")
    plot.new(); title(main = paste("Data for", toupper(pth), "not available."))
    next
  }
  pathogen_df <- master_pheno_table %>%
    select(all_of(pheno_cols)) %>%
    set_names(c("p_soft_std", "BL_hard", "Mix_hard")) %>%
    drop_na()

  if (nrow(pathogen_df) < 100) {
    cat("     ! Skipping pathogen", toupper(pth), "- not enough data.\n")
    plot.new(); title(main = paste("Not enough data for", toupper(pth), "diagnostics."))
    next
  }

  pathogen_df <- pathogen_df %>% mutate(p_soft_0_1 = plogis(p_soft_std))

  corr_df <- pathogen_df %>% select(p_soft_0_1, BL_hard, Mix_hard)
  corr_matrix <- cor(corr_df)
  pheatmap(
    corr_matrix, display_numbers = TRUE, number_format = "%.3f",
    color = colorRampPalette(c("#4575b4", "white", "#d73027"))(50),
    main = glue("Correlation of Phenotype Definitions for {toupper(pth)}"),
    fontsize = 12, border_color = "grey60"
  )
  
  prob_dist_plot <- ggplot(pathogen_df, aes(x = p_soft_0_1)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "#7570b3", alpha = 0.8, color = "white") +
    geom_density(color = "#1b9e77", linewidth = 1.2) +
    labs(
      title = glue("Distribution of Scaled Seropositivity Probability for {toupper(pth)}"),
      subtitle = "Standardized scores transformed to a 0-1 scale using the logistic function",
      x = "Scaled Seropositivity Probability (0-1 Scale)", y = "Density"
    ) + theme_minimal()
  print(prob_dist_plot)
}
dev.off()
cat("--- Pathogen-Level Diagnostics Complete. Report saved to:", pdf_report_path, "---\n")

cat("\n\n--- UNIFIED DIAGNOSTICS SCRIPT COMPLETE ---\n\n")