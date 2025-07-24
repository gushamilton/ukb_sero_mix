################################################################################
#   UKB SEROLOGY: PATHOGEN-LEVEL PHENOTYPES (BAYESIAN REGRESSION) - SCRIPT 06
################################################################################
#  ▸ Strategy: use the single-antigen posterior probabilities (p_soft) generated
#    by script 02 as predictors in a Bayesian logistic regression that
#    approximates the UK Biobank rule (≥2 antigens positive, or pgp3 for CT).
#  ▸ This learns data-driven weights for each antigen while remaining robust and
#    fast – logistic regression is easy for MCMC to sample.
#  ▸ For development the script subsamples ≤2,000 complete cases per pathogen.
#    Comment out the slice_sample() line to run on the full dataset.
################################################################################

# --- 1. PACKAGES --------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, glue, brms, future, here)

plan(multisession, workers = max(1, availableCores() - 1))
options(mc.cores = max(1, availableCores() - 1))

# --- 2. INPUT ----------------------------------------------------------------
pheno_path <- "quickdraws_input/phenotypes_master_postqc.tsv"
if (!file.exists(pheno_path)) stop("Cannot find ", pheno_path, ". Run script 02 first.")
master <- read_tsv(pheno_path, show_col_types = FALSE, guess_max = 5000)
cat("Loaded", nrow(master), "individuals\n")

# --- 3. Pathogen → antigen map ----------------------------------------------
pathogen_map <- list(
  EBV = c("ebv_vca", "ebv_ebna1", "ebv_zebra", "ebv_ead"),
  CMV = c("cmv_pp150", "cmv_pp52", "cmv_pp28"),
  CT  = c("ct_pgp3",  "ct_mompa", "ct_mompd", "ct_tarpf1", "ct_tarpf2"),
  HP  = c("hp_caga",  "hp_vaca", "hp_omp", "hp_groel", "hp_catalase", "hp_urea")
)

# UKB rule as noisy target -----------------------------------------------------
ukb_rules <- list(
  EBV = list(type = "count", n = 2),
  CMV = list(type = "count", n = 2),
  HP  = list(type = "count", n = 2),
  CT  = list(type = "pgp3")
)

make_label <- function(df, antigens, rule){
  if (rule$type == "count") {
    cols <- paste0(antigens, "_sero_hard")
    as.integer(rowSums(df[cols], na.rm = TRUE) >= rule$n)
  } else { # pgp3 rule
    df$ct_pgp3_sero_hard
  }
}

# --- 4. Loop over pathogens ---------------------------------------------------
results <- list()

for (pth in names(pathogen_map)) {
  antigens <- pathogen_map[[pth]]
  soft_cols <- paste0(antigens, "_sero_soft")
  hard_cols <- paste0(antigens, "_sero_hard")

  df <- master %>%
    select(FID, IID, all_of(c(soft_cols, hard_cols))) %>%
    drop_na()
  if (nrow(df) < 100) {
    warning("Skipping ", pth, ": less than 100 complete cases")
    next
  }

  df$z <- make_label(df, antigens, ukb_rules[[pth]])

  # Subsample for development speed -----------------------------------------
  set.seed(42)
  if (nrow(df) > 2000) df_model <- slice_sample(df, n = 2000) else df_model <- df

  # Build formula -----------------------------------------------------------
  rhs <- paste(soft_cols, collapse = " + ")
  brm_formula <- bf(as.formula(paste("z ~", rhs)), family = bernoulli(link = "logit"))

  cat("\nFitting Bayesian logistic model for", pth, "on", nrow(df_model), "individuals...\n")

  fit <- brm(brm_formula, data = df_model,
             chains = 2, iter = 2000, warmup = 500,
             refresh = 0,
             prior = set_prior("normal(0,1)", class = "b"))

  # Predict on *all* complete-case rows (not just the subsample)
  pi_hat <- fitted(fit, newdata = df, scale = "response")[, "Estimate"]

  results[[pth]] <- tibble(
    FID = df$FID,
    IID = df$IID,
    !!glue("{tolower(pth)}_sero_soft") := pi_hat,
    !!glue("{tolower(pth)}_sero_hard") := as.integer(pi_hat >= 0.5),
    !!glue("{tolower(pth)}_w_soft")    := pi_hat,
    !!glue("{tolower(pth)}_w_hard")    := 2 * abs(pi_hat - 0.5)
  )
}

combined <- reduce(results, full_join, by = c("FID", "IID"))

# --- 5. Write outputs ---------------------------------------------------------
if (!dir.exists("quickdraws_input")) dir.create("quickdraws_input")

write_tsv(combined, "quickdraws_input/phenotypes_pathogen_bayesReg.tsv")
cat("\nSaved pathogen-level Bayesian-regression phenotypes: phenotypes_pathogen_bayesReg.tsv\n")

# Write weight files per pathogen --------------------------------------------
for (pth in names(results)) {
  low <- tolower(pth)
  res <- results[[pth]]
  write_tsv(res %>% select(FID, IID, Weight = !!glue("{low}_w_soft")),
            glue("quickdraws_input/{low}_p_soft.weights"))
  write_tsv(res %>% select(FID, IID, Weight = !!glue("{low}_w_hard")),
            glue("quickdraws_input/{low}_sero_hard.weights"))
}

cat("Weight files written.\nScript 06 complete.\n") 