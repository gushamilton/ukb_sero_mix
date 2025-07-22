################################################################################
#   UKB SEROLOGY: CREATE PHENOTYPE & COVARIATE FILES FOR QUICKDRAWS
################################################################################
#  ∙ This script processes UKB serology data to generate phenotype, covariate,
#    and weights files suitable for a comprehensive GWAS analysis using Quickdraws.
#  ∙ It implements both hard-cutoff and mixture-model-based "soft" phenotypes.
#  ∙ It also generates a novel covariate, PC1_neg, representing the principal
#    component of seronegativity across all selected antigens.
################################################################################


# --- 1. SETUP & CONFIGURATION ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, mixsmsn, sn, here, glue, pheatmap)

# Ensure the output directory exists
if (!dir.exists("quickdraws_input")) dir.create("quickdraws_input")

# Define the final core list of antigens to be processed
# This list is based on visual inspection of plots and lambda statistics
core_antigens <- c(
    "cmv_pp150", "cmv_pp28", "cmv_pp52",
    "hsv1", "hsv2",
    "ebv_vca", "ebv_ebna1", "ebv_zebra", "ebv_ead",
    "hp_omp", "hp_urea", "hp_caga", "hp_vaca",
    "toxo_p22",
    "bkv_vp1", "jcv_vp1", "mcv_vp1",
    "hhv6_ie1a", "hhv6_ie1b", "hhv7_u14",
    "ct_pgp3", "ct_mompd", "ct_tarpf2",
    "kshv_lana"
)
cat("Defined", length(core_antigens), "core antigens for analysis.\n")

# This maps a short, usable name to the exact (messy) column header
# from the title file and includes the official MFI threshold.
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


# --- 2. LOAD & PREPARE DATA ---
cat("Loading and cleaning serology data...\n")
# The raw data file is expected to be in the parent directory
raw_data <- read_tsv("serology_export_title.tsv", show_col_types = FALSE, guess_max = 5000)

# Select and rename only the columns we need, including FID/IID. Drop rows with missing values across all
cols_to_keep <- setNames(antigen_map$full_name, antigen_map$short_name)
data_clean <- raw_data %>%
  select(FID, IID, any_of(cols_to_keep)) %>%
  rename(any_of(cols_to_keep)) %>%
  mutate(across(-c(FID, IID), as.numeric)) %>%
  drop_na()
  
cat("Data loaded for", ncol(data_clean) - 2, "antigens and", nrow(data_clean), "samples.\n")


# --- 3. CORE FUNCTIONS ---

#' Fit a 2-component skew-t mixture model
#'
#' @param y Numeric vector of MFI values.
#' @return A fitted smsn.mix object or NULL on error.
fit_skew_mix <- function(y) {
  y_log <- log(y[y > 0 & is.finite(y)])
  if (length(y_log) < 200) return(NULL) # Not enough data to fit
  
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

#' Process a single antigen to generate all phenotype and weight columns
#'
#' @param antigen_short_name The short name of the antigen (e.g., "cmv_pp150").
#' @param mfi_data A tibble with FID, IID, and a column for the MFI values.
#' @param threshold The official hard seropositivity threshold.
#' @return A tibble with FID, IID, and all generated phenotype/weight columns.
process_antigen <- function(antigen_short_name, mfi_data, threshold) {
  cat("  -> Processing:", antigen_short_name, "\n")
  
  # Fit the mixture model
  fit <- fit_skew_mix(mfi_data[[antigen_short_name]])
  if (is.null(fit)) {
    warning("Mixture model failed for ", antigen_short_name, ". Skipping.")
    return(NULL)
  }
  
  # Identify positive component and get posterior probabilities
  pos_comp <- which.max(fit$mu)
  p_soft_df <- tibble(
      y_log = fit$y,
      p_soft = fit$obs.prob[, pos_comp]
  )
  
  # Join posteriors back to the original data
  mfi_data %>%
    mutate(y_log = log(.data[[antigen_short_name]])) %>%
    left_join(p_soft_df, by = "y_log") %>%
    select(-y_log) %>%
    # Generate the 4 phenotype/weight columns
    mutate(
      !!glue("{antigen_short_name}_IgG_raw") := log(.data[[antigen_short_name]]),
      !!glue("{antigen_short_name}_sero_hard") := if_else(.data[[antigen_short_name]] >= threshold, 1, 0),
      !!glue("{antigen_short_name}_sero_soft") := p_soft,
      !!glue("{antigen_short_name}_p_neg") := 1 - p_soft, # For PC1_neg
      # Weights for downstream analysis
      !!glue("{antigen_short_name}_w_hard") := 2 * abs(p_soft - 0.5),
      !!glue("{antigen_short_name}_w_igg") := p_soft
    )
}


# --- 4. MAIN ANALYSIS LOOP ---
cat("\n--- Starting Phenotype Generation ---\n")
all_pheno_list <- list(select(data_clean, FID, IID))

for (i in 1:nrow(antigen_map)) {
  ab_name <- antigen_map$short_name[i]
  threshold <- antigen_map$threshold[i]
  
  if (!ab_name %in% names(data_clean)) next
  
  antigen_phenos <- process_antigen(ab_name, data_clean %>% select(FID, IID, !!sym(ab_name)), threshold)
  
  if (!is.null(antigen_phenos)) {
    all_pheno_list[[ab_name]] <- select(antigen_phenos, -starts_with("p_soft"), -!!sym(ab_name))
  }
}

# Combine all generated phenotypes into a single master tibble
master_pheno_table <- reduce(all_pheno_list, left_join, by = c("FID", "IID"))
cat("Master phenotype table generated with", ncol(master_pheno_table) - 2, "columns.\n")


# --- 5. GENERATE PC1_neg COVARIATE ---
cat("\n--- Generating PC1_neg Covariate ---\n")

# Create the matrix of seronegative probabilities
p_neg_matrix <- master_pheno_table %>%
  select(FID, IID, ends_with("_p_neg")) %>%
  column_to_rownames("IID") %>%
  select(-FID) %>%
  as.matrix()

# Impute missing values (e.g., from failed fits) with the column mean
for(j in 1:ncol(p_neg_matrix)){
    p_neg_matrix[is.na(p_neg_matrix[,j]), j] <- mean(p_neg_matrix[,j], na.rm = TRUE)
}

# Run PCA
pca_results <- prcomp(p_neg_matrix, center = TRUE, scale. = TRUE)
pc1_neg_df <- pca_results$x %>%
  as_tibble(rownames = "IID") %>%
  select(IID, PC1_neg = PC1) %>%
  mutate(IID = as.integer(IID))

cat("PC1_neg explains", round(summary(pca_results)$importance[2,1] * 100, 2), "% of variance.\n")


# --- 6. CREATE QUICKDRAWS INPUT FILES ---
cat("\n--- Writing Quickdraws Input Files ---\n")

# 6a. Master Covariate File
# NOTE: This assumes `age_sex.txt` and `ukb_pcs.txt` are in the parent directory.
# You may need to adjust paths for your environment.
covar_base <- read_tsv("../03_quickdraw_gwas/age_sex.txt")
pcs_base <- read_tsv("../pcs.txt") # Placeholder for PCs 1-20

master_covar_table <- covar_base %>%
  left_join(pcs_base, by = c("FID", "IID")) %>%
  left_join(pc1_neg_df, by = "IID") %>%
  select(FID, IID, age, sex, starts_with("PC"), PC1_neg)

write_tsv(master_covar_table, "quickdraws_input/covariates_master.tsv")
cat("  -> Wrote master covariate file: quickdraws_input/covariates_master.tsv\n")


# 6b. Phenotype and Weights Files
pheno_cols <- master_pheno_table %>% select(-ends_with("_w_hard"), -ends_with("_w_igg"), -ends_with("_p_neg"))

# Split quantitative (IgG + sero_soft) and binary (sero_hard) phenotypes
phenos_quant <- pheno_cols %>% select(FID, IID, matches("_IgG_raw$|_sero_soft$"))
phenos_binary <- pheno_cols %>% select(FID, IID, matches("_sero_hard$"))

write_tsv(phenos_quant, "quickdraws_input/phenotypes_quantitative.tsv")
write_tsv(phenos_binary, "quickdraws_input/phenotypes_binary.tsv")
cat("  -> Wrote separate phenotype files: quantitative and binary.\n")

# (Optionally keep the combined file for reference)
write_tsv(pheno_cols, "quickdraws_input/phenotypes_master.tsv")
cat("  -> Wrote master phenotype file (all traits): quickdraws_input/phenotypes_master.tsv\n")

# Create a manifest of all analyses to run
analysis_manifest <- list()

for (ab in core_antigens) {
    # IgG Raw (no weights)
    analysis_manifest[[length(analysis_manifest) + 1]] <- tibble(phenotype_name = glue("{ab}_IgG_raw"), analysis_type = "linear", weights_file = NA)
    
    # IgG Weighted
    weights_igg <- master_pheno_table %>% select(FID, IID, Weight = !!glue("{ab}_w_igg"))
    write_tsv(weights_igg, glue("quickdraws_input/{ab}_IgG_wgt.weights"), col_names = TRUE)
    analysis_manifest[[length(analysis_manifest) + 1]] <- tibble(phenotype_name = glue("{ab}_IgG_raw"), analysis_type = "linear", weights_file = glue("{ab}_IgG_wgt.weights"))

    # Sero Hard (with weights)
    weights_hard <- master_pheno_table %>% select(FID, IID, Weight = !!glue("{ab}_w_hard"))
    write_tsv(weights_hard, glue("quickdraws_input/{ab}_sero_hard.weights"), col_names = TRUE)
    analysis_manifest[[length(analysis_manifest) + 1]] <- tibble(phenotype_name = glue("{ab}_sero_hard"), analysis_type = "logistic", weights_file = glue("{ab}_sero_hard.weights"))
    
    # Sero Soft (no weights)
    analysis_manifest[[length(analysis_manifest) + 1]] <- tibble(phenotype_name = glue("{ab}_sero_soft"), analysis_type = "linear", weights_file = NA)
}

manifest_df <- bind_rows(analysis_manifest)
write_tsv(manifest_df, "quickdraws_input/analysis_manifest.tsv")
cat("  -> Wrote", nrow(manifest_df), "weight files and analysis manifest.\n")


# --- 7. QC & REPORTING on PC1_neg ---
cat("\n--- QC for PC1_neg ---\n")
cor_matrix <- cor(
    master_pheno_table %>% select(ends_with("_IgG_raw")),
    pc1_neg_df$PC1_neg,
    use = "pairwise.complete.obs"
)

cor_df <- as_tibble(cor_matrix, rownames = "Antigen") %>%
    rename(Correlation = V1) %>%
    mutate(Antigen = str_remove(Antigen, "_IgG_raw")) %>%
    arrange(desc(abs(Correlation)))

cat("Top 5 correlations with PC1_neg:\n")
print(head(cor_df, 5))

p <- ggplot(cor_df, aes(x = Correlation, y = reorder(Antigen, Correlation))) +
    geom_col(fill = "steelblue") +
    labs(
        title = "Correlation of PC1_neg with Raw log-MFI of Each Antigen",
        subtitle = "PC1_neg captures shared variance in seronegativity",
        x = "Pearson Correlation", y = "Antigen"
    ) +
    theme_light()

ggsave("quickdraws_input/pc1_neg_correlation_plot.pdf", p, width = 8, height = 10)
cat("  -> Saved PC1_neg correlation plot.\n")

cat("\n--- Phenotype Generation Complete ---\n") 