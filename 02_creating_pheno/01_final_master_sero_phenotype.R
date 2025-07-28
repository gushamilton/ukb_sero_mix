################################################################################
#
#   UKB SEROLOGY: UNIFIED PHENOTYPE GENERATION SCRIPT
#
#   ▸ Strategy: This script consolidates the full phenotype generation pipeline,
#     from raw MFI data to GWAS-ready input files. It implements the
#     five-analysis strategy for both individual antigens and aggregated pathogens.
#
################################################################################

# --- 1. SETUP & CONFIGURATION -------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
    tidyverse, glue, brms, future, here, bayestestR, patchwork, scales,
    mixsmsn, sn, furrr, data.table
)

# Setup parallel processing
plan(multisession, workers = max(1, availableCores() - 1))
options(mc.cores = max(1, availableCores() - 1))

# Ensure output directories exist
if (!dir.exists("quickdraws_input")) dir.create("quickdraws_input")
if (!dir.exists("results")) dir.create("results")


# --- 2. GLOBAL DEFINITIONS ----------------------------------------------------

# Define the final core list of antigens to be processed
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
cat("Defined", length(core_antigens), "core antigens for analysis.\n")

# This maps a short, usable name to the exact column header and official MFI threshold.
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
) %>%
filter(short_name %in% unique(core_antigens))

# Pathogen to antigen map for aggregation step
pathogen_map <- list(
  EBV = c("ebv_vca", "ebv_ebna1", "ebv_zebra", "ebv_ead"),
  CMV = c("cmv_pp150", "cmv_pp52", "cmv_pp28"),
  CT  = c("ct_pgp3",  "ct_mompa", "ct_mompd", "ct_tarpf1", "ct_tarpf2"),
  HP  = c("hp_caga",  "hp_vaca", "hp_omp", "hp_groel", "hp_catalase", "hp_urea")
)

# UKB rules (Butler-Laporte et al. 2023) for pathogen definition
ukb_rules <- list(
  EBV = list(type = "count", n = 2),
  CMV = list(type = "count", n = 2),
  HP  = list(type = "count_exclude_caga", n = 2),
  CT  = list(type = "pgp3_or_complex")
)


# --- 3. DATA LOADING ----------------------------------------------------------
cat("Loading and cleaning serology data...\n")
# The raw data file is expected to be in the parent directory
raw_data_path <- "serology_export_title.tsv"
if (!file.exists(raw_data_path)) stop("Raw data file not found: ", raw_data_path)
raw_data <- read_tsv(raw_data_path, show_col_types = FALSE, guess_max = 5000)

# Select and rename only the columns we need
cols_to_keep <- setNames(antigen_map$full_name, antigen_map$short_name)
data_clean <- raw_data %>%
  select(FID = `Participant ID`, IID = `Participant ID`, any_of(cols_to_keep)) %>%
  rename(any_of(cols_to_keep)) %>%
  mutate(across(-c(FID, IID), as.numeric)) %>%
  filter(!if_all(-c(FID, IID), is.na))

cat("Data loaded for", ncol(data_clean) - 2, "antigens and", nrow(data_clean), "samples.\n")


# --- 4. CORE FUNCTIONS --------------------------------------------------------

#' Fit a 2-component skew-t mixture model
fit_skew_mix <- function(y) {
  y_log <- log(y[y > 0 & is.finite(y)])
  if (length(y_log) < 200) return(NULL)
  tryCatch({
    smsn.mix(y_log, g = 2, family = "Skew.t", nu = 4, group = FALSE, calc.im = FALSE, obs.prob = TRUE)
  }, error = function(e) { NULL })
}

#' Calculate the optimal threshold from p_soft using Youden's J statistic
calculate_youden_threshold <- function(p_soft, truth_hard) {
  df <- tibble(p_soft, truth_hard) %>% drop_na()
  if (nrow(df) < 100 || n_distinct(df$truth_hard) < 2) return(0.5) # Default if data is sparse

  thresholds <- seq(0.05, 0.95, by = 0.01)
  youden_scores <- map_dbl(thresholds, ~{
    pred_hard <- as.integer(df$p_soft >= .x)
    tp <- sum(pred_hard == 1 & df$truth_hard == 1)
    tn <- sum(pred_hard == 0 & df$truth_hard == 0)
    fp <- sum(pred_hard == 1 & df$truth_hard == 0)
    fn <- sum(pred_hard == 0 & df$truth_hard == 1)
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    youden <- sensitivity + specificity - 1
    if (is.nan(youden)) 0 else youden
  })
  
  best_threshold <- thresholds[which.max(youden_scores)]
  return(best_threshold)
}


#' Create a pathogen-level hard label based on UK Biobank rules
make_ukb_rule_label <- function(df, antigens, rule, col_suffix = "_sero_hard_BL"){
  # Ensure suffix is applied to antigen names
  antigens_with_suffix <- paste0(antigens, col_suffix)

  if (rule$type == "count") {
    existing_cols <- intersect(antigens_with_suffix, names(df))
    as.integer(rowSums(df[existing_cols], na.rm = TRUE) >= rule$n)
  } else if (rule$type == "count_exclude_caga") {
    cols_no_caga <- antigens_with_suffix[!grepl("caga", antigens_with_suffix)]
    existing_cols <- intersect(cols_no_caga, names(df))
    as.integer(rowSums(df[existing_cols], na.rm = TRUE) >= rule$n)
  } else if (rule$type == "pgp3_or_complex") {
    pgp3_col <- paste0("ct_pgp3", col_suffix)
    other_cols <- paste0(antigens[antigens != "ct_pgp3"], col_suffix)
    existing_other_cols <- intersect(other_cols, names(df))
    if (!pgp3_col %in% names(df)) return(rep(NA_integer_, nrow(df)))

    pgp3_positive <- df[[pgp3_col]] == 1
    pgp3_negative_but_others <- (df[[pgp3_col]] == 0) & (rowSums(df[existing_other_cols], na.rm = TRUE) >= 2)
    as.integer(pgp3_positive | pgp3_negative_but_others)
  } else {
    stop("Unknown rule type: ", rule$type)
  }
}


# --- 5. ANTIGEN-LEVEL PHENOTYPE GENERATION --------------------------------------
cat("\n--- Starting Antigen-Level Phenotype Generation ---\n")

process_antigen <- function(antigen_short_name, mfi_data, threshold) {
  cat("  -> Processing:", antigen_short_name, "\n")
  titres <- mfi_data[[antigen_short_name]]
  
  idx_valid <- which(titres > 0 & is.finite(titres))
  if (length(idx_valid) < 200) {
      warning("Not enough valid data for ", antigen_short_name, ". Skipping.")
      return(NULL)
  }
  
  fit <- fit_skew_mix(titres[idx_valid])
  if (is.null(fit)) {
      warning("Mixture model failed for ", antigen_short_name, ". Skipping.")
      return(NULL)
  }
  
  pos_comp <- which.max(fit$mu)
  
  # Generate raw p_soft
  p_soft_vec_raw <- rep(NA_real_, length(titres))
  p_soft_vec_raw[idx_valid] <- fit$obs.prob[, pos_comp]
  
  # Apply monotonicity fix
  prob_df <- tibble(idx = idx_valid, mfi = titres[idx_valid], prob_raw = p_soft_vec_raw[idx_valid]) %>%
      arrange(mfi) %>%
      mutate(prob_mono = cummax(prob_raw))
  p_soft_vec <- rep(NA_real_, length(titres))
  p_soft_vec[prob_df$idx] <- prob_df$prob_mono
  
  # --- Generate the 5 required phenotype sets ---
  sero_hard_BL <- if_else(titres >= threshold, 1, 0, missing = 0)
  igg_seropos_BL <- if_else(sero_hard_BL == 1, log(titres), NA_real_)
  youden_thresh <- calculate_youden_threshold(p_soft_vec, sero_hard_BL)
  sero_hard_mix <- if_else(p_soft_vec >= youden_thresh, 1, 0)
  p_soft_std <- scale(p_soft_vec)[,1]
  
  pheno_tibble <- tibble(
    FID = mfi_data$FID, IID = mfi_data$IID,
    !!glue("{antigen_short_name}_IgG_raw") := log(titres),
    !!glue("{antigen_short_name}_p_soft") := p_soft_vec,
    !!glue("{antigen_short_name}_sero_hard_BL") := sero_hard_BL,
    !!glue("{antigen_short_name}_IgG_seropos_BL") := igg_seropos_BL,
    !!glue("{antigen_short_name}_sero_hard_mix") := sero_hard_mix,
    !!glue("{antigen_short_name}_psoft_std") := p_soft_std
  )
  
  # Return both phenotypes and the model fit
  list(phenotypes = pheno_tibble, model_fit = fit, antigen_name = antigen_short_name)
}

# Process all antigens in parallel
args_for_pmap <- list(
  antigen_short_name = antigen_map$short_name,
  mfi_data = map(antigen_map$short_name, ~select(data_clean, FID, IID, all_of(.x))),
  threshold = antigen_map$threshold
)

is_present <- antigen_map$short_name %in% names(data_clean)
args_for_pmap <- map(args_for_pmap, ~ .x[is_present])

# pmap now returns a list of lists (phenotypes + model_fit)
list_of_results <- future_pmap(args_for_pmap, process_antigen, .options = furrr_options(seed = TRUE))

# Filter out NULLs from failed fits and extract components
valid_results <- compact(list_of_results)
list_of_pheno_tibbles <- map(valid_results, "phenotypes")
all_model_fits <- map(valid_results, "model_fit")
names(all_model_fits) <- map_chr(valid_results, "antigen_name")

saveRDS(all_model_fits, "mixture_model_fits.rds")
cat("  -> Saved mixture model fits for diagnostics: mixture_model_fits.rds\n")

# Combine all generated phenotypes into a single master table
all_valid_tibbles <- c(list(select(data_clean, FID, IID)), list_of_pheno_tibbles)
antigen_pheno_table <- reduce(all_valid_tibbles, left_join, by = c("FID", "IID"))

cat("Antigen-level phenotype table generated with", ncol(antigen_pheno_table) - 2, "columns.\n")


# --- 6. PATHOGEN-LEVEL PHENOTYPE GENERATION -------------------------------------
cat("\n--- Starting Pathogen-Level Phenotype Generation ---\n")

pathogen_results_list <- list()

for (pth in names(pathogen_map)) {
  antigens <- pathogen_map[[pth]]
  soft_cols <- paste0(antigens, "_p_soft")
  
  required_cols <- c("FID", "IID", soft_cols)
  
  if (pth == "HP") {
      required_cols <- required_cols[!grepl("caga", required_cols)]
  }
  
  existing_cols <- intersect(required_cols, names(antigen_pheno_table))
  if(length(existing_cols) < 3) {
      warning("Skipping ", pth, ": not enough antigen data.")
      next
  }

  df <- antigen_pheno_table %>%
    select(all_of(existing_cols)) %>%
    drop_na()

  if (nrow(df) < 100) {
    warning("Skipping ", pth, ": less than 100 complete cases")
    next
  }

  df_bl <- antigen_pheno_table %>%
      select(FID, IID, ends_with("_sero_hard_BL")) %>%
      filter(IID %in% df$IID)
  
  df$z <- make_ukb_rule_label(df_bl, antigens, ukb_rules[[pth]], col_suffix = "_sero_hard_BL")
  
  if(n_distinct(df$z) < 2) {
    warning("Skipping ", pth, ": outcome variable 'z' is unary.")
    next
  }

  existing_soft_cols <- intersect(soft_cols, names(df))
  rhs <- paste(existing_soft_cols, collapse = " + ")
  brm_formula <- bf(as.formula(paste("z ~", rhs)), family = bernoulli(link = "logit"))

  cat(glue("\nFitting Bayesian logistic model for {pth} on {nrow(df)} individuals...\n"))
  fit <- brm(brm_formula, data = df, chains = 2, iter = 2000, warmup = 500, refresh = 0)

  pred_full <- fitted(fit, newdata = df, scale = "response")
  pathogen_p_soft <- pred_full[, "Estimate"]
  
  pathogen_sero_hard_BL <- make_ukb_rule_label(df_bl, antigens, ukb_rules[[pth]], col_suffix = "_sero_hard_BL")
  youden_thresh_pathogen <- calculate_youden_threshold(pathogen_p_soft, pathogen_sero_hard_BL)
  pathogen_sero_hard_mix <- if_else(pathogen_p_soft >= youden_thresh_pathogen, 1, 0)
  pathogen_psoft_std <- scale(pathogen_p_soft)[,1]
  
  pathogen_results_list[[pth]] <- tibble(
      FID = df$FID, IID = df$IID,
      !!glue("{tolower(pth)}_p_soft") := pathogen_p_soft,
      !!glue("{tolower(pth)}_sero_hard_BL") := pathogen_sero_hard_BL,
      !!glue("{tolower(pth)}_sero_hard_mix") := pathogen_sero_hard_mix,
      !!glue("{tolower(pth)}_psoft_std") := pathogen_psoft_std
  )
}

if (length(pathogen_results_list) > 0) {
    pathogen_pheno_table <- reduce(pathogen_results_list, full_join, by = c("FID", "IID"))
    cat("Pathogen-level phenotype table generated.\n")
    master_pheno_table <- full_join(antigen_pheno_table, pathogen_pheno_table, by = c("FID", "IID"))
} else {
    master_pheno_table <- antigen_pheno_table
    cat("No pathogen-level phenotypes were generated.\n")
}


# --- 7. COVARIATE GENERATION ----------------------------------------------------
cat("\n--- Generating Covariates ---\n")

# 7a. Phenotype-derived latent factors
igg_levels_matrix <- master_pheno_table %>% select(FID, IID, ends_with("_IgG_raw"))

# latent_factor_1
igg_neg_matrix <- igg_levels_matrix
for (ab_name in core_antigens) {
    igg_col <- glue("{ab_name}_IgG_raw"); soft_col <- glue("{ab_name}_p_soft")
    if (all(c(igg_col, soft_col) %in% names(master_pheno_table))) {
        idx_to_mask <- which(master_pheno_table[[soft_col]] >= 0.1)
        if (length(idx_to_mask) > 0) igg_neg_matrix[idx_to_mask, igg_col] <- NA
    }
}
pca_input_neg <- igg_neg_matrix %>% select(-FID, -IID) %>% as.matrix()
for(j in 1:ncol(pca_input_neg)) pca_input_neg[is.na(pca_input_neg[,j]), j] <- mean(pca_input_neg[,j], na.rm = TRUE)
pca_input_neg[is.na(pca_input_neg)] <- 0
col_vars_neg <- apply(pca_input_neg, 2, var, na.rm = TRUE)
constant_cols_neg <- which(col_vars_neg == 0 | is.na(col_vars_neg))
if (length(constant_cols_neg) > 0) {
    cat("  -> Removing", length(constant_cols_neg), "constant columns before seronegative PCA.\n")
    pca_input_neg <- pca_input_neg[, -constant_cols_neg, drop = FALSE]
}
pca_results_neg <- prcomp(pca_input_neg, center = TRUE, scale. = TRUE)
master_pheno_table$latent_factor_1 <- pca_results_neg$x[,1]

# latent_factor_2
igg_pos_matrix <- igg_levels_matrix
for (ab_name in core_antigens) {
    igg_col <- glue("{ab_name}_IgG_raw"); soft_col <- glue("{ab_name}_p_soft")
    if (all(c(igg_col, soft_col) %in% names(master_pheno_table))) {
        idx_to_mask <- which(master_pheno_table[[soft_col]] <= 0.9)
        if (length(idx_to_mask) > 0) igg_pos_matrix[idx_to_mask, igg_col] <- NA
    }
}
pca_input_pos <- igg_pos_matrix %>% select(-FID, -IID) %>% as.matrix()
for(j in 1:ncol(pca_input_pos)) pca_input_pos[is.na(pca_input_pos[,j]), j] <- mean(pca_input_pos[,j], na.rm = TRUE)
pca_input_pos[is.na(pca_input_pos)] <- 0
col_vars_pos <- apply(pca_input_pos, 2, var, na.rm = TRUE)
constant_cols_pos <- which(col_vars_pos == 0 | is.na(col_vars_pos))
if (length(constant_cols_pos) > 0) {
    cat("  -> Removing", length(constant_cols_pos), "constant columns before seropositive PCA.\n")
    pca_input_pos <- pca_input_pos[, -constant_cols_pos, drop = FALSE]
}
pca_results_pos <- prcomp(pca_input_pos, center = TRUE, scale. = TRUE)
master_pheno_table$latent_factor_2 <- pca_results_pos$x[,1]

# latent_factor_igg
pca_input_global <- igg_levels_matrix %>% select(-FID, -IID) %>% as.matrix()
for(j in 1:ncol(pca_input_global)) pca_input_global[is.na(pca_input_global[,j]), j] <- mean(pca_input_global[,j], na.rm = TRUE)
pca_input_global[is.na(pca_input_global)] <- 0
col_vars_global <- apply(pca_input_global, 2, var, na.rm = TRUE)
constant_cols_global <- which(col_vars_global == 0 | is.na(col_vars_global))
if (length(constant_cols_global) > 0) {
    cat("  -> Removing", length(constant_cols_global), "constant columns before global IgG PCA.\n")
    pca_input_global <- pca_input_global[, -constant_cols_global, drop = FALSE]
}
pca_results_global <- prcomp(pca_input_global, center = TRUE, scale. = TRUE)
master_pheno_table$latent_factor_igg <- pca_results_global$x[,1]

# 7b. External covariates
age_sex_data <- read_delim("age_sex.txt", show_col_types = FALSE)
genetic_pcs <- fread("covariates_clean.txt", sep = " ")

base_covar <- age_sex_data %>%
  left_join(genetic_pcs, by = c("FID", "IID")) %>%
  select(FID, IID, age, sex, all_of(paste0("PC", 1:20)))
write_tsv(base_covar, "quickdraws_input/covariates_base.tsv")
cat("  -> Wrote base covariates (age, sex, genetic PCs 1-20)\n")

all_covar <- base_covar %>%
  left_join(master_pheno_table %>% select(FID, IID, latent_factor_1, latent_factor_2, latent_factor_igg), by = c("FID", "IID"))
write_tsv(all_covar, "quickdraws_input/covariates_all.tsv")
cat("  -> Wrote combined covariates (all)\n")


# --- 8. BUILD FULL MANIFEST & WRITE COMPLETE FILES FOR VISUALIZATION ------------
cat("\n--- Building Full Manifest for All Phenotypes ---\n")

manifest_list <- list()

# 8a. Add Antigen Phenotypes to Manifest
successful_antigens <- names(all_model_fits)
for (ab in successful_antigens) {
    manifest_list[[length(manifest_list) + 1]] <- tibble(id = "BL_SER0", phenotype_name = glue("{ab}_sero_hard_BL"), analysis_type = "logistic", weights_file = NA)
    manifest_list[[length(manifest_list) + 1]] <- tibble(id = "MIX_SER0", phenotype_name = glue("{ab}_sero_hard_mix"), analysis_type = "logistic", weights_file = NA)
    manifest_list[[length(manifest_list) + 1]] <- tibble(id = "BL_IGG_SERPOS", phenotype_name = glue("{ab}_IgG_seropos_BL"), analysis_type = "linear", weights_file = NA)
    manifest_list[[length(manifest_list) + 1]] <- tibble(id = "P_SOFT", phenotype_name = glue("{ab}_psoft_std"), analysis_type = "linear", weights_file = NA)
    manifest_list[[length(manifest_list) + 1]] <- tibble(id = "IGG_WGT", phenotype_name = glue("{ab}_IgG_raw"), analysis_type = "linear", weights_file = glue("{ab}_p_soft.weights"))
}

# 8b. Add Pathogen Phenotypes to Manifest
successful_pathogens <- tolower(names(pathogen_results_list))
for (pth in successful_pathogens) {
    manifest_list[[length(manifest_list) + 1]] <- tibble(id = "BL_SER0", phenotype_name = glue("{pth}_sero_hard_BL"), analysis_type = "logistic", weights_file = NA)
    manifest_list[[length(manifest_list) + 1]] <- tibble(id = "MIX_SER0", phenotype_name = glue("{pth}_sero_hard_mix"), analysis_type = "logistic", weights_file = NA)
    manifest_list[[length(manifest_list) + 1]] <- tibble(id = "P_SOFT", phenotype_name = glue("{pth}_psoft_std"), analysis_type = "linear", weights_file = NA)
}

# 8c. Add Latent Factors to Manifest
all_latent_factors <- c("latent_factor_1", "latent_factor_2", "latent_factor_igg")
for (lt in all_latent_factors) {
    manifest_list[[length(manifest_list) + 1]] <- tibble(id = "LATENT_FACTOR", phenotype_name = lt, analysis_type = "linear", weights_file = NA)
}

manifest_full <- bind_rows(manifest_list)

# 8d. Perform Phenotype QC on the full set of phenotypes
phenotypes_to_check <- unique(manifest_full$phenotype_name)
excluded_traits <- future_map_dfr(phenotypes_to_check, function(pheno_name) {
  if (!pheno_name %in% names(master_pheno_table)) { return(tibble(phenotype_name = pheno_name, reason = "Column not found")) }
  pheno_vec_clean <- na.omit(master_pheno_table[[pheno_name]])
  if (length(pheno_vec_clean) == 0) return(tibble(phenotype_name = pheno_name, reason = "All values NA"))
  if (n_distinct(pheno_vec_clean) < 2) return(tibble(phenotype_name = pheno_name, reason = "Unary trait"))
  is_binary <- all(pheno_vec_clean %in% c(0, 1))
  if (is_binary && min(table(pheno_vec_clean)) < 100) return(tibble(phenotype_name = pheno_name, reason = "Minor cat < 100"))
  if (!is_binary && sd(pheno_vec_clean) < 0.01) return(tibble(phenotype_name = pheno_name, reason = "SD < 0.01"))
  return(NULL)
}, .options = furrr_options(seed = TRUE))

if(nrow(excluded_traits) > 0) {
  write_tsv(excluded_traits, "quickdraws_input/excluded_traits_qc.tsv")
  cat("  -> QC: Excluded", nrow(excluded_traits), "traits. See 'excluded_traits_qc.tsv'.\n")
} else {
  cat("  -> QC: All traits passed variance and unary checks.\n")
}

good_trait_names <- setdiff(phenotypes_to_check, excluded_traits$phenotype_name)
manifest_full <- manifest_full %>% filter(phenotype_name %in% good_trait_names)

# 8e. Write COMPLETE, UNFILTERED files for plotting and visualization
cat("\n--- Writing Complete Output Files for Visualization ---\n")
phenos_quant_full <- master_pheno_table %>% select(FID, IID, any_of(manifest_full %>% filter(analysis_type == "linear") %>% pull(phenotype_name) %>% unique()))
phenos_binary_full <- master_pheno_table %>% select(FID, IID, any_of(manifest_full %>% filter(analysis_type == "logistic") %>% pull(phenotype_name) %>% unique()))
write_tsv(phenos_quant_full, "quickdraws_input/phenotypes_quantitative.tsv")
write_tsv(phenos_binary_full, "quickdraws_input/phenotypes_binary.tsv")
write_tsv(manifest_full, "quickdraws_input/analysis_manifest.tsv")
cat("  -> Wrote complete phenotype files and manifest to quickdraws_input/\n")

for (w_file in na.omit(unique(manifest_full$weights_file))) {
    ab_name <- str_remove(w_file, "_p_soft\\.weights")
    weight_col <- glue("{ab_name}_p_soft")
    if (weight_col %in% names(master_pheno_table)) {
        weights_data <- master_pheno_table %>% select(FID, IID, Weight = all_of(weight_col))
        write_tsv(weights_data, file.path("quickdraws_input", w_file))
    }
}
cat("  -> Wrote all corresponding weight files.\n")


# --- 9. FILTER & WRITE GWAS-READY FILES -----------------------------------------
cat("\n--- Creating Filtered, GWAS-Ready Files ---\n")

# 9a. Apply Pathogen-Level Filtering to Manifest for GWAS
cat("  -> Applying pathogen-level filtering to create GWAS analysis set...\n")
suppressed_antigens <- unlist(pathogen_map, use.names = FALSE)
phenotypes_to_exclude <- c(
    glue("{suppressed_antigens}_sero_hard_BL"),
    glue("{suppressed_antigens}_sero_hard_mix"),
    glue("{suppressed_antigens}_psoft_std")
)
manifest_gwas <- manifest_full %>% 
    filter(!phenotype_name %in% phenotypes_to_exclude)
cat("  -> Created a filtered manifest for GWAS with", nrow(manifest_gwas), "analyses.\n")

# 9b. Filter for Genotyped Participants
gwas_ready_dir <- "quickdraws_input/gwas_ready_files"
if (!dir.exists(gwas_ready_dir)) dir.create(gwas_ready_dir)
fam_file1 <- "genotype_antibody_only.fam"
fam_file2 <- "chr1_imputed_antibody.fam"

if (!file.exists(fam_file1) || !file.exists(fam_file2)) {
    warning("One or both .fam files not found. Skipping filtering step for GWAS-ready files.")
} else {
    ids_geno1 <- fread(fam_file1, header = FALSE, select = 2, data.table = FALSE)[[1]]
    ids_geno2 <- fread(fam_file2, header = FALSE, select = 2, data.table = FALSE)[[1]]
    
    genotyped_ids <- intersect(ids_geno1, ids_geno2)
    final_ids_for_gwas <- intersect(genotyped_ids, master_pheno_table$IID)
    cat(glue("  -> Found {length(final_ids_for_gwas)} participants with both genotype and phenotype data.\n"))
    cat(glue("  -> Writing filtered files to: {gwas_ready_dir}\n"))

    # Write filtered manifest
    write_tsv(manifest_gwas, file.path(gwas_ready_dir, "analysis_manifest.tsv"))

    # Create filtered phenotype tables based on the GWAS manifest
    phenos_quant_gwas <- master_pheno_table %>% select(FID, IID, any_of(manifest_gwas %>% filter(analysis_type == "linear") %>% pull(phenotype_name) %>% unique()))
    phenos_binary_gwas <- master_pheno_table %>% select(FID, IID, any_of(manifest_gwas %>% filter(analysis_type == "logistic") %>% pull(phenotype_name) %>% unique()))

    # Filter by genotype ID and write
    write_tsv(phenos_quant_gwas %>% filter(IID %in% final_ids_for_gwas), file.path(gwas_ready_dir, "phenotypes_quantitative.tsv"))
    write_tsv(phenos_binary_gwas %>% filter(IID %in% final_ids_for_gwas), file.path(gwas_ready_dir, "phenotypes_binary.tsv"))

    # Filter and write covariate and weight files
    write_tsv(all_covar %>% filter(IID %in% final_ids_for_gwas), file.path(gwas_ready_dir, "covariates_all.tsv"))
    write_tsv(base_covar %>% filter(IID %in% final_ids_for_gwas), file.path(gwas_ready_dir, "covariates_base.tsv"))

    for (w_file in na.omit(unique(manifest_gwas$weights_file))) {
        original_path <- file.path("quickdraws_input", w_file)
        if (file.exists(original_path)) {
            weights_data <- read_tsv(original_path, show_col_types = FALSE) %>%
                filter(IID %in% final_ids_for_gwas)
            write_tsv(weights_data, file.path(gwas_ready_dir, w_file))
        }
    }
    cat(glue("  -> Filtering complete. GWAS-ready files are in {gwas_ready_dir}\n"))
}


# --- 10. QC REPORTING & VISUALIZATION -------------------------------------------
cat("\n--- Generating Final QC Reports ---\n")

latent_cor <- cor(master_pheno_table %>% select(all_of(all_latent_factors)), use = "pairwise.complete.obs")
write_tsv(as_tibble(latent_cor, rownames = "Latent_Factor"), "results/latent_factor_correlations.tsv")
cat("  -> Saved latent factor correlations.\n")

# Summary of seroprevalence
qc_passed_prefixes <- unique(str_remove(manifest_full$phenotype_name, "_IgG_raw|_p_soft|_sero_hard_BL|_IgG_seropos_BL|_sero_hard_mix|_psoft_std"))
summary_table <- map_dfr(qc_passed_prefixes, function(trait_prefix) {
  bl_col <- glue("{trait_prefix}_sero_hard_BL")
  mix_col <- glue("{trait_prefix}_sero_hard_mix")
  if (!all(c(bl_col, mix_col) %in% names(master_pheno_table))) return(NULL)
  tibble(
    trait = trait_prefix,
    seroprev_BL = mean(master_pheno_table[[bl_col]], na.rm = TRUE),
    seroprev_mix = mean(master_pheno_table[[mix_col]], na.rm = TRUE)
  )
})
write_tsv(summary_table, "results/seroprevalence_comparison_summary.tsv")
cat("  -> Saved seroprevalence comparison summary.\n")


cat("\n\n--- UNIFIED PHENOTYPE SCRIPT COMPLETE ---\n\n")