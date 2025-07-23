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
pacman::p_load(tidyverse, mixsmsn, sn, here, glue, pheatmap, furrr, data.table)

# Setup parallel processing
plan(multisession, workers = availableCores() - 1)

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
  select(FID = `Participant ID`, IID = `Participant ID`, any_of(cols_to_keep)) %>%
  rename(any_of(cols_to_keep)) %>%
  mutate(across(-c(FID, IID), as.numeric)) %>%
  filter(!if_all(-c(FID, IID), is.na))

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
#' @return A list containing the phenotype tibble and the mixture model fit, or NULL on failure.
process_antigen <- function(antigen_short_name, mfi_data, threshold) {
  cat("  -> Processing:", antigen_short_name, "\n")
  titres <- mfi_data[[antigen_short_name]]
  
  # Indices of valid (positive) MFI values for model fitting
  idx_valid <- which(titres > 0 & is.finite(titres))
  if (length(idx_valid) < 200) {
    warning("Not enough valid data points for ", antigen_short_name, ". Skipping.")
    return(NULL)
  }
  
  # Fit model on log-transformed valid titres
  fit <- fit_skew_mix(titres[idx_valid])
  if (is.null(fit)) {
    warning("Mixture model failed for ", antigen_short_name, ". Skipping.")
    return(NULL)
  }
  
  pos_comp <- which.max(fit$mu)
  
  # Initialise a raw probability vector
p_soft_vec_raw <- rep(NA_real_, length(titres))
p_soft_vec_raw[idx_valid] <- fit$obs.prob[, pos_comp]

# --- START: Monotonicity Fix ---
# This fix addresses the biologically implausible "dip" in probability at very high MFI values.
# It ensures that the probability of being seropositive cannot decrease as MFI increases.
prob_df <- tibble(
    idx = idx_valid,
    mfi = titres[idx_valid],
    prob_raw = p_soft_vec_raw[idx_valid]
) %>%
    arrange(mfi) %>% # 1. Sort by MFI value
    mutate(prob_mono = cummax(prob_raw)) # 2. Calculate the cumulative maximum probability

# Create the final, corrected probability vector
p_soft_vec <- rep(NA_real_, length(titres))
p_soft_vec[prob_df$idx] <- prob_df$prob_mono
# --- END: Monotonicity Fix ---
  
  # Build the phenotype tibble
  pheno_tibble <- tibble(
    FID = mfi_data$FID,
    IID = mfi_data$IID,
    !!glue("{antigen_short_name}_IgG_raw")   := log(titres),
    !!glue("{antigen_short_name}_sero_hard") := if_else(titres >= threshold, 1, 0, missing = 0),
    !!glue("{antigen_short_name}_sero_soft") := p_soft_vec,
    !!glue("{antigen_short_name}_p_neg")     := 1 - p_soft_vec,
    !!glue("{antigen_short_name}_w_hard")    := 2 * abs(p_soft_vec - 0.5),
    !!glue("{antigen_short_name}_w_igg")     := p_soft_vec
  )
  
  # Return both phenotype data and model fit
  list(
    phenotypes = pheno_tibble,
    model_fit = fit,
    antigen_name = antigen_short_name,
    threshold = threshold,
    idx_valid = idx_valid
  )
}


# --- 4. MAIN ANALYSIS LOOP (PARALLELIZED) ---
cat("\n--- Starting Phenotype Generation ---\n")

# Prepare a list of arguments for pmap
args_for_pmap <- list(
  antigen_short_name = antigen_map$short_name,
  mfi_data = map(antigen_map$short_name, ~select(data_clean, FID, IID, all_of(.x))),
  threshold = antigen_map$threshold
)

# Filter for antigens actually present in the data
is_present <- antigen_map$short_name %in% names(data_clean)
args_for_pmap <- map(args_for_pmap, ~ .x[is_present])

# Process all antigens in parallel using pmap
# It returns a list of lists, one for each successfully processed antigen
list_of_results <- future_pmap(args_for_pmap, process_antigen, .options = furrr_options(seed = TRUE))

# Filter out NULLs from failed fits
valid_results <- compact(list_of_results)

# Extract phenotype tibbles and model fits
list_of_pheno_tibbles <- map(valid_results, ~.x$phenotypes)
all_model_fits <- map(valid_results, ~.x$model_fit)
names(all_model_fits) <- map(valid_results, ~.x$antigen_name)

# Save all model fits as a single RDS file
saveRDS(all_model_fits, "mixture_model_fits.rds")
cat("  -> Saved mixture model fits: mixture_model_fits.rds\n")

# Filter out NULLs from failed fits and prepend the base FID/IID tibble
# The base tibble ensures all samples are kept, even if some have no phenotype data
all_valid_tibbles <- c(
    list(select(data_clean, FID, IID)), 
    list_of_pheno_tibbles
)

# Combine all generated phenotypes into a single master tibble by joining on FID/IID
master_pheno_table <- reduce(all_valid_tibbles, left_join, by = c("FID", "IID"))

cat("Master phenotype table generated with", ncol(master_pheno_table) - 2, "columns.\n")


# --- 5. GENERATE PC1_assay_noise COVARIATE (IMPROVED LOGIC) ---
cat("\n--- Generating PC1_assay_noise Covariate ---\n")

# Start with a matrix of the raw IgG levels for all core antigens
igg_levels_matrix <- master_pheno_table %>% 
    select(FID, IID, ends_with("_IgG_raw"))

# For each antigen, set the IgG level to NA for individuals who are NOT confidently seronegative
# A "confident negative" is defined as p_soft < 0.1
for (ab_name in core_antigens) {
    igg_col <- glue("{ab_name}_IgG_raw")
    soft_col <- glue("{ab_name}_sero_soft")

    if (all(c(igg_col, soft_col) %in% names(master_pheno_table))) {
        # Get indices of individuals who are NOT confident negatives
        idx_to_mask <- which(master_pheno_table[[soft_col]] >= 0.1)
        if (length(idx_to_mask) > 0) {
            igg_levels_matrix[idx_to_mask, igg_col] <- NA
        }
    }
}

# Now, prepare the matrix for PCA: numeric, with IID as rownames
pca_input_matrix <- igg_levels_matrix %>%
  select(-any_of(c("FID", "IID"))) %>%
  # drop antigens that failed and dont exist in the master table
  select(any_of(glue("{core_antigens}_IgG_raw"))) %>% 
  as.matrix()

# Impute the remaining NAs (column-wise) using the mean of the confident negatives for that antigen
for(j in 1:ncol(pca_input_matrix)){
    col_mean <- mean(pca_input_matrix[,j], na.rm = TRUE)
    if (is.finite(col_mean)) {
      pca_input_matrix[is.na(pca_input_matrix[,j]), j] <- col_mean
    } else {
      # If all values are NA (shouldn't happen), fill with 0
      pca_input_matrix[,j] <- 0
    }
}

# Check for columns with zero variance and remove them before PCA
col_vars <- apply(pca_input_matrix, 2, var, na.rm = TRUE)
constant_cols <- which(col_vars == 0 | is.na(col_vars))
if (length(constant_cols) > 0) {
    cat("  -> Removing", length(constant_cols), "constant columns before PCA:\n")
    cat("     ", paste(colnames(pca_input_matrix)[constant_cols], collapse = ", "), "\n")
    pca_input_matrix <- pca_input_matrix[, -constant_cols, drop = FALSE]
}

# Check if we have enough columns for PCA
if (ncol(pca_input_matrix) < 2) {
    cat("  -> Warning: Not enough variable columns for PCA. Creating dummy PC1_assay_noise.\n")
    pc1_assay_noise_df <- tibble(IID = data_clean$IID, PC1_assay_noise = 0)
    variance_explained_neg <- 0
} else {
    # Run PCA for seronegative IgG levels
    pca_results_neg <- prcomp(pca_input_matrix, center = TRUE, scale. = TRUE)
    pc1_assay_noise_df <- pca_results_neg$x %>%
      as_tibble() %>%
      select(PC1_assay_noise = PC1) %>%
      bind_cols(select(data_clean, IID), .)
    
    variance_explained_neg <- summary(pca_results_neg)$importance["Proportion of Variance", "PC1"]
}
cat("PC1_assay_noise explains", round(variance_explained_neg * 100, 2), "% of variance in seronegative IgG levels.\n")

# --- 5b. GENERATE PC1_seropositive COVARIATE (FOR COMPARISON) ---
cat("\n--- Generating PC1_seropositive Covariate ---\n")

# Start with a matrix of the raw IgG levels for all core antigens
igg_levels_matrix_pos <- master_pheno_table %>% 
    select(FID, IID, ends_with("_IgG_raw"))

# For each antigen, set the IgG level to NA for individuals who are NOT confidently seropositive
# A "confident positive" is defined as p_soft > 0.9
for (ab_name in core_antigens) {
    igg_col <- glue("{ab_name}_IgG_raw")
    soft_col <- glue("{ab_name}_sero_soft")

    if (all(c(igg_col, soft_col) %in% names(master_pheno_table))) {
        # Get indices of individuals who are NOT confident positives
        idx_to_mask <- which(master_pheno_table[[soft_col]] <= 0.9)
        if (length(idx_to_mask) > 0) {
            igg_levels_matrix_pos[idx_to_mask, igg_col] <- NA
        }
    }
}

# Now, prepare the matrix for PCA: numeric, with IID as rownames
pca_input_matrix_pos <- igg_levels_matrix_pos %>%
  select(-any_of(c("FID", "IID"))) %>%
  # drop antigens that failed and dont exist in the master table
  select(any_of(glue("{core_antigens}_IgG_raw"))) %>% 
  as.matrix()

# Impute the remaining NAs (column-wise) using the mean of the confident positives for that antigen
for(j in 1:ncol(pca_input_matrix_pos)){
    col_mean <- mean(pca_input_matrix_pos[,j], na.rm = TRUE)
    if (is.finite(col_mean)) {
      pca_input_matrix_pos[is.na(pca_input_matrix_pos[,j]), j] <- col_mean
    } else {
      # If all values are NA (shouldn't happen), fill with 0
      pca_input_matrix_pos[,j] <- 0
    }
}

# Check for columns with zero variance and remove them before PCA
col_vars_pos <- apply(pca_input_matrix_pos, 2, var, na.rm = TRUE)
constant_cols_pos <- which(col_vars_pos == 0 | is.na(col_vars_pos))
if (length(constant_cols_pos) > 0) {
    cat("  -> Removing", length(constant_cols_pos), "constant columns before seropositive PCA:\n")
    cat("     ", paste(colnames(pca_input_matrix_pos)[constant_cols_pos], collapse = ", "), "\n")
    pca_input_matrix_pos <- pca_input_matrix_pos[, -constant_cols_pos, drop = FALSE]
}

# Check if we have enough columns for PCA
if (ncol(pca_input_matrix_pos) < 2) {
    cat("  -> Warning: Not enough variable columns for seropositive PCA. Creating dummy PC1_seropositive.\n")
    pc1_seropositive_df <- tibble(IID = data_clean$IID, PC1_seropositive = 0)
    variance_explained_pos <- 0
} else {
    # Run PCA for seropositive IgG levels
    pca_results_pos <- prcomp(pca_input_matrix_pos, center = TRUE, scale. = TRUE)
    pc1_seropositive_df <- pca_results_pos$x %>%
      as_tibble() %>%
      select(PC1_seropositive = PC1) %>%
      bind_cols(select(data_clean, IID), .)
    
    variance_explained_pos <- summary(pca_results_pos)$importance["Proportion of Variance", "PC1"]
}
cat("PC1_seropositive explains", round(variance_explained_pos * 100, 2), "% of variance in seropositive IgG levels.\n")


# --- 6. CREATE QUICKDRAWS INPUT FILES ---
cat("\n--- Writing Quickdraws Input Files ---\n")

# 6a. Master Covariate File ---
cat("Reading covariate files with a post-hoc fix for mixed delimiters...\n")

# Reading age/sex should be fine
covar_base <- read_delim("age_sex.txt")

# Read the messy PC file without headers
pcs_fixed <- fread("covariates_clean.txt", sep = " ")

# --- Post-hoc cleanup based on the file structure ---
# 
# --- Perform the join ---
# Merge but **exclude** PC1 covariates for the first GWAS runs
master_covar_table <- merge(covar_base, pcs_fixed, by = c("FID", "IID")) %>%
  as_tibble() %>%
  select(FID, IID, age, sex, all_of(paste0("PC", 1:20)))

write_tsv(master_covar_table, "quickdraws_input/covariates_master.tsv")
cat("  -> Wrote master covariate file (no extra PC1 covars): quickdraws_input/covariates_master.tsv\n")


# 6b. Phenotype and Weights Files
pheno_cols <- master_pheno_table %>% select(-ends_with("_w_hard"), -ends_with("_w_igg"), -ends_with("_p_neg"))

# Split quantitative (IgG + sero_soft + IgG_seropos_only) and binary (sero_hard) phenotypes
phenos_quant <- pheno_cols %>% select(FID, IID, matches("_IgG_raw$|_IgG_seropos_only$|_sero_soft$"))
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
    # Create Butler-Laporte style IgG phenotype: IgG levels only in seropositives (NA for seronegatives)
    igg_seropos_only <- master_pheno_table %>%
        mutate(!!glue("{ab}_IgG_seropos_only") := if_else(!!sym(glue("{ab}_sero_hard")) == 1, !!sym(glue("{ab}_IgG_raw")), NA_real_))
    
    # Update master_pheno_table with the new filtered phenotype
    master_pheno_table <- master_pheno_table %>%
        mutate(!!glue("{ab}_IgG_seropos_only") := igg_seropos_only[[glue("{ab}_IgG_seropos_only")]])
    
    # Analysis A: Classic binary serostatus (unweighted)
    analysis_manifest[[length(analysis_manifest) + 1]] <- tibble(phenotype_name = glue("{ab}_sero_hard"), analysis_type = "logistic", weights_file = NA)
    
    # Analysis B: Butler-Laporte IgG analysis (IgG only in seropositives, unweighted)
    analysis_manifest[[length(analysis_manifest) + 1]] <- tibble(phenotype_name = glue("{ab}_IgG_seropos_only"), analysis_type = "linear", weights_file = NA)
    
    # Analysis C: IgG Raw with p_soft weights (probabilistic weighting)
    weights_igg <- master_pheno_table %>% select(FID, IID, Weight = !!glue("{ab}_w_igg"))
    write_tsv(weights_igg, glue("quickdraws_input/{ab}_IgG_wgt.weights"), col_names = TRUE)
    analysis_manifest[[length(analysis_manifest) + 1]] <- tibble(phenotype_name = glue("{ab}_IgG_raw"), analysis_type = "linear", weights_file = glue("{ab}_IgG_wgt.weights"))
    
    # Analysis D1: Sero Hard with p_soft weights
    weights_psoft <- master_pheno_table %>% select(FID, IID, Weight = !!glue("{ab}_sero_soft"))
    write_tsv(weights_psoft, glue("quickdraws_input/{ab}_p_soft.weights"), col_names = TRUE)
    analysis_manifest[[length(analysis_manifest) + 1]] <- tibble(phenotype_name = glue("{ab}_sero_hard"), analysis_type = "logistic", weights_file = glue("{ab}_p_soft.weights"))
    
    # Analysis D2: Sero Hard with hard weights (emphasize confident cases/controls)
    weights_hard <- master_pheno_table %>% select(FID, IID, Weight = !!glue("{ab}_w_hard"))
    write_tsv(weights_hard, glue("quickdraws_input/{ab}_sero_hard.weights"), col_names = TRUE)
    analysis_manifest[[length(analysis_manifest) + 1]] <- tibble(phenotype_name = glue("{ab}_sero_hard"), analysis_type = "logistic", weights_file = glue("{ab}_sero_hard.weights"))
    
    # Analysis E: Sero Soft (quantitative, unweighted)
    analysis_manifest[[length(analysis_manifest) + 1]] <- tibble(phenotype_name = glue("{ab}_sero_soft"), analysis_type = "linear", weights_file = NA)
}

manifest_df <- bind_rows(analysis_manifest)
write_tsv(manifest_df, "quickdraws_input/analysis_manifest.tsv")
cat("  -> Wrote", nrow(manifest_df), "weight files and analysis manifest.\n")

# --- 6c. Create per-analysis minimal phenotype / covariate files for Step-2 ---
cat("  -> Generating per-analysis phenotype files for Step-2 runs...\n")

step2_dir <- "quickdraws_input/step2"
if (!dir.exists(step2_dir)) dir.create(step2_dir)

# Write one covariate file (identical for all analyses)
write_tsv(master_covar_table, file.path(step2_dir, "covariates_master.tsv"))

# Create minimal phenotype files - NOW after all phenotypes are created
for (row_idx in seq_len(nrow(manifest_df))) {
  pheno_name <- manifest_df$phenotype_name[row_idx]
  pheno_file <- file.path(step2_dir, glue("{pheno_name}.phen.tsv"))
  # Remove the file.exists check since we want to create/overwrite all files
  tmp <- master_pheno_table %>% select(FID, IID, all_of(pheno_name))
  write_tsv(tmp, pheno_file)
}
cat("  -> Wrote minimal phenotype files to", step2_dir, "\n")


# --- 7. QC & REPORTING on PC1_assay_noise ---
cat("\n--- QC for PC1_assay_noise ---\n")
cor_matrix <- cor(
    master_pheno_table %>% select(ends_with("_IgG_raw")),
    pc1_assay_noise_df$PC1_assay_noise,
    use = "pairwise.complete.obs"
)

cor_df <- as_tibble(cor_matrix, rownames = "Antigen") %>%
    rename(Correlation = V1) %>%
    mutate(Antigen = str_remove(Antigen, "_IgG_raw")) %>%
    arrange(desc(abs(Correlation)))

cat("Top 5 correlations with PC1_assay_noise:\n")
print(head(cor_df, 5))

p <- ggplot(cor_df, aes(x = Correlation, y = reorder(Antigen, Correlation))) +
    geom_col(fill = "steelblue") +
    labs(
        title = "Correlation of PC1_assay_noise with Raw log-MFI of Each Antigen",
        subtitle = "PC1 from IgG levels of seronegatives, representing shared assay noise",
        x = "Pearson Correlation", y = "Antigen"
    ) +
    theme_light()

ggsave("quickdraws_input/pc1_assay_noise_correlation_plot.pdf", p, width = 8, height = 10)
cat("  -> Saved PC1_assay_noise correlation plot.\n")

cat("\n--- Phenotype Generation Complete ---\n") 


# --- 8. WITHIN-ANTIGEN PHENOTYPE CORRELATIONS ---
cat("\n--- Correlating Phenotype Definitions ---\n")

# A helper function to safely calculate raw correlations for one antigen
calculate_correlations <- function(ab_name, data) {
  # Define column names
  igg_col <- glue("{ab_name}_IgG_raw")
  soft_col <- glue("{ab_name}_sero_soft")
  hard_col <- glue("{ab_name}_sero_hard")
  
  # Ensure all necessary columns exist
  required_cols <- c(igg_col, soft_col, hard_col)
  if (!all(required_cols %in% names(data))) {
    return(NULL)
  }
  
  # Subset to non-missing data for this antigen
  subset_data <- data %>%
    select(all_of(required_cols)) %>%
    na.omit()
    
  if (nrow(subset_data) < 2) return(NULL)

  # Calculate standard Pearson correlations
  cor_igg_soft <- cor(subset_data[[igg_col]], subset_data[[soft_col]])
  cor_hard_soft <- cor(subset_data[[hard_col]], subset_data[[soft_col]])
  
  tibble(
    antigen = ab_name,
    cor_igg_vs_soft = cor_igg_soft,
    cor_hard_vs_soft = cor_hard_soft
  )
}

# Apply the function to all core antigens and bind the results
phenotype_correlations <- map_dfr(core_antigens, calculate_correlations, data = master_pheno_table)

# Save the correlation results
write_tsv(phenotype_correlations, "quickdraws_input/phenotype_correlations.tsv")
cat("  -> Wrote phenotype correlation summary: quickdraws_input/phenotype_correlations.tsv\n")

# Optional: Create a visualization of the correlations
cor_long <- phenotype_correlations %>%
  pivot_longer(
    cols = -antigen,
    names_to = "correlation_type",
    values_to = "correlation_value"
  )

p_cor <- ggplot(cor_long, aes(x = correlation_type, y = correlation_value, fill = correlation_type)) +
  geom_col() +
  facet_wrap(~antigen, ncol = 4) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(
    title = "Within-Antigen Phenotype Correlations (Raw Data)",
    subtitle = "Comparing different phenotype definitions for each antigen",
    x = "Correlation Type",
    y = "Pearson Correlation"
  )
  
ggsave("quickdraws_input/phenotype_correlation_plot.pdf", p_cor, width = 12, height = 15)
cat("  -> Saved phenotype correlation plot.\n")


cat("\n--- Analysis Script Complete ---\n") 

