################################################################################
#   UKB SEROLOGY: PATHOGEN-LEVEL SOFT PHENOTYPES (SCRIPT 04)
################################################################################
#  ∙ Combines multiple antigen measurements for each pathogen into a single
#    posterior infection probability using a 2-component multivariate Gaussian
#    mixture (full covariance) fitted with mclust.
#  ∙ Per-antigen mixture probabilities from Script 02 are used to initialise
#    the EM algorithm, preserving previous information and accelerating
#    convergence.
#  ∙ Outputs (per pathogen):
#      ‑ <pathogen>_sero_soft   – posterior P(infected)
#      ‑ <pathogen>_sero_hard   – MAP call (prob ≥ 0.5)
#      ‑ <pathogen>_w_soft      – weight = p_sero_soft
#      ‑ <pathogen>_w_hard      – weight = 2 × |p – 0.5|
#  ∙ Saves fitted models to: pathogen_mixture_model_fits.rds
#  ∙ Writes phenotype & weight files into quickdraws_input/, mirroring Script 02
################################################################################

# --- 1. SETUP & PACKAGES -------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, mclust, furrr, glue, data.table, here)

plan(multisession, workers = availableCores() - 1)

# Ensure output directory exists ------------------------------------------------
if (!dir.exists("quickdraws_input")) dir.create("quickdraws_input")

# --- 2. INPUT DATA -------------------------------------------------------------
cat("Loading master phenotype table (Script 02 output)…\n")
pheno_file <- "quickdraws_input/phenotypes_master_postqc.tsv"
if (!file.exists(pheno_file)) {
  stop("Required file not found: ", pheno_file,
       "\n-> Run 02_create_serology_phenotypes.R first.")
}
master_pheno <- read_tsv(pheno_file, show_col_types = FALSE, guess_max = 5000)
cat("Loaded", nrow(master_pheno), "participants.\n")

# --- 3. PATHOGEN → ANTIGEN MAP -------------------------------------------------
pathogen_map <- list(
  EBV   = c("ebv_vca",  "ebv_ebna1", "ebv_zebra", "ebv_ead"),
  CMV   = c("cmv_pp150", "cmv_pp52",  "cmv_pp28"),
  HHV6  = c("hhv6_ie1a", "hhv6_ie1b", "hhv6_p101k"),
  # HBV and HCV removed (too rare for reliable mixture fitting)
  # HPV16 removed (not available in phenotype data)
  CT    = c("ct_pgp3",  "ct_mompa",  "ct_mompd", "ct_tarpf1", "ct_tarpf2", "ct_porb"),
  HP    = c("hp_caga",  "hp_vaca",  "hp_omp", "hp_groel", "hp_catalase", "hp_urea")
)

# --- 4. HELPER FUNCTIONS -------------------------------------------------------

# Get the *_IgG_raw column names for a set of antigens -------------------------
get_igg_cols <- function(antigens) glue("{antigens}_IgG_raw")

# Get per-antigen soft probability columns -------------------------------------
get_soft_cols <- function(antigens) glue("{antigens}_sero_soft")

# Fit pathogen-level mixture ----------------------------------------------------
fit_pathogen_mixture <- function(pathogen, antigens, data_tbl) {
  cat("\n-- Fitting", pathogen, "(", length(antigens), "antigens)…\n")

  # Keep only antigens present in the data
  igg_cols  <- intersect(get_igg_cols(antigens), names(data_tbl))
  soft_cols <- intersect(get_soft_cols(antigens), names(data_tbl))

  if (length(igg_cols) < 2) {
    warning("  -> Skipping ", pathogen, ": need ≥2 antigen columns, found ", length(igg_cols))
    return(NULL)
  }

  # Assemble data matrix (log-MFI); mean-impute NAs for mixture fitting
  mat <- data_tbl %>% select(all_of(igg_cols)) %>% as.matrix()
  for (j in seq_len(ncol(mat))) {
    col_mean <- mean(mat[, j], na.rm = TRUE)
    if (is.finite(col_mean)) {
      mat[is.na(mat[, j]), j] <- col_mean
    }
  }

  # Initialise with per-antigen soft probabilities (if available)
  init_class <- NULL
  if (length(soft_cols) == length(igg_cols)) {
    p_soft_mat <- data_tbl %>% select(all_of(soft_cols)) %>% as.matrix()
    row_means  <- rowMeans(p_soft_mat, na.rm = TRUE)
    # NA rows fall back to 0 (assume seronegative)
    row_means[is.na(row_means)] <- 0
    init_class <- ifelse(row_means > 0.5, 2, 1)  # mclust classes are 1/2
  }

  # Fit 2-component full-covariance Gaussian mixture
  m <- Mclust(mat, G = 2, modelNames = "VVV",
              initialization = list(classification = init_class))

  # Identify positive component = larger average log-MFI across antigens
  comp_means <- m$parameters$mean  # matrix(antigen × component)
  pos_comp <- which.max(colSums(comp_means))

  p_soft <- m$z[, pos_comp]
  p_soft[is.na(p_soft)] <- 0  # should not happen

  tibble(
    FID = data_tbl$FID,
    IID = data_tbl$IID,
    !!glue("{tolower(pathogen)}_sero_soft") := p_soft,
    !!glue("{tolower(pathogen)}_sero_hard") := if_else(p_soft >= 0.5, 1, 0),
    !!glue("{tolower(pathogen)}_w_soft")    := p_soft,
    !!glue("{tolower(pathogen)}_w_hard")    := 2 * abs(p_soft - 0.5)
  ) %>%
    list(model_fit = m)
}

# --- 5. MAIN LOOP --------------------------------------------------------------
cat("\n--- Starting pathogen-level mixture fitting ---\n")

pathogen_results <- future_map(names(pathogen_map), function(p) {
  fit_pathogen_mixture(p, pathogen_map[[p]], master_pheno)
}, .options = furrr_options(seed = TRUE))

# Remove NULLs (skipped pathogens)
valid_res <- pathogen_results[!map_lgl(pathogen_results, is.null)]

if (length(valid_res) == 0) {
  stop("No pathogen models were successfully fitted. Check that antigen columns exist in the phenotype data.")
}

# Extract phenotype tibbles and model objects -----------------------------------
pheno_tibbles <- map(valid_res, 1)
model_fits    <- map(valid_res, 2)
names(model_fits) <- tolower(names(pathogen_map))[!map_lgl(pathogen_results, is.null)]

# Save model fits ----------------------------------------------------------------
saveRDS(model_fits, "pathogen_mixture_model_fits.rds")
cat("  -> Saved mixture models to pathogen_mixture_model_fits.rds\n")

# --- 6. MERGE & WRITE OUTPUTS --------------------------------------------------
cat("\n--- Writing phenotype & weight files ---\n")

# Merge all pathogen phenotypes
pathogen_pheno_tbl <- reduce(pheno_tibbles, left_join, by = c("FID", "IID"))

# Write combined phenotype table
write_tsv(pathogen_pheno_tbl, "quickdraws_input/phenotypes_pathogen.tsv")
cat("  -> Wrote quickdraws_input/phenotypes_pathogen.tsv\n")

# Weight files (soft + hard) ----------------------------------------------------
for (pathogen in names(pathogen_map)) {
  soft_col <- glue("{tolower(pathogen)}_w_soft")
  hard_col <- glue("{tolower(pathogen)}_w_hard")
  if (all(c(soft_col, hard_col) %in% names(pathogen_pheno_tbl))) {
    soft_w <- pathogen_pheno_tbl %>% select(FID, IID, Weight = !!soft_col)
    write_tsv(soft_w, glue("quickdraws_input/{tolower(pathogen)}_p_soft.weights"))

    hard_w <- pathogen_pheno_tbl %>% select(FID, IID, Weight = !!hard_col)
    write_tsv(hard_w, glue("quickdraws_input/{tolower(pathogen)}_sero_hard.weights"))
  }
}
cat("  -> Wrote per-pathogen weight files.\n")

# --- 7. SEROPREVALENCE SUMMARY ----------------------------------------------
cat("\n--- Seroprevalence estimates (hard calls) ---\n")
seroprev <- map_dbl(names(pathogen_map), function(p) {
  col <- glue("{tolower(p)}_sero_hard")
  if (col %in% names(pathogen_pheno_tbl)) {
    mean(pathogen_pheno_tbl[[col]] == 1, na.rm = TRUE)
  } else NA_real_
})

seroprev_tbl <- tibble(
  pathogen = names(pathogen_map),
  seroprevalence = round(seroprev * 100, 2)
)

print(seroprev_tbl)

write_tsv(seroprev_tbl, "quickdraws_input/pathogen_seroprevalence.tsv")
cat("  -> Wrote seroprevalence summary to quickdraws_input/pathogen_seroprevalence.tsv\n")

cat("\n--- Script 04 complete ---\n") 