################################################################################
#   UKB SEROLOGY: PATHOGEN-LEVEL PHENOTYPES (BAYESIAN) - SCRIPT 05
################################################################################
#  ∙ This script provides a Bayesian alternative to the frequentist mixture
#    models, using the `brms` package (a frontend for Stan).
#  ∙ For each pathogen, it fits a univariate 2-component Gaussian mixture model
#    to each of its constituent antigens. This is the Bayesian equivalent of
#    the model in script 02 and is highly robust.
#  ∙ The resulting posterior probabilities of seropositivity for each antigen
#    are then combined using a simple mean to get a final pathogen-level call.
#  ∙ NOTE: This script is configured to run on a small subsample (n=200) for
#    speed. To run on the full dataset, comment out the subsampling line.
################################################################################

# --- 1. SETUP & PACKAGES -------------------------------------------------------
# brms may require a one-time setup of its backend, CmdStanR.
# If it's the first time, run: install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, glue, data.table, here, brms, future)

# Enable parallel processing for brms
plan(multisession, workers = availableCores() - 1)
options(mc.cores = availableCores() - 1)


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
  CT    = c("ct_pgp3",  "ct_mompa",  "ct_mompd", "ct_tarpf1", "ct_tarpf2"),
  HP    = c("hp_caga",  "hp_vaca",  "hp_omp", "hp_groel", "hp_catalase", "hp_urea")
)

# --- 4. MAIN BAYESIAN ANALYSIS LOOP --------------------------------------------
cat("\n--- Starting Bayesian pathogen-level modeling ---\n")

# Process each pathogen
pathogen_pheno_list <- map(names(pathogen_map), function(pathogen_name) {
  
  cat("\n-- Processing", pathogen_name, "...\n")
  antigens <- pathogen_map[[pathogen_name]]
  igg_cols <- glue("{antigens}_IgG_raw")
  
  # Find individuals with complete data for this pathogen's antigens
  complete_cases <- master_pheno %>%
    select(all_of(c("FID", "IID", igg_cols))) %>%
    na.omit()
  
  if(nrow(complete_cases) < 200) {
      warning("Skipping ", pathogen_name, ", not enough complete cases.")
      return(NULL)
  }

  # ============================================================================
  # === FOR DEVELOPMENT: Subsample to 200 individuals for speed.             ===
  # === COMMENT OUT THE NEXT LINE TO RUN ON THE FULL DATASET.                ===
  # ============================================================================
  set.seed(123) # for reproducibility
  analysis_data <- sample_n(complete_cases, 200)
  cat("  -> Using a subsample of", nrow(analysis_data), "for rapid testing.\n")


  # This will store the posterior probabilities for each antigen
  posterior_probs <- tibble(FID = analysis_data$FID, IID = analysis_data$IID)
  
  # Loop through each antigen for the current pathogen
  for (antigen in antigens) {
    col_name <- glue("{antigen}_IgG_raw")
    cat("    -> Fitting Bayesian mixture model for:", antigen, "\n")
    
    # Define the 2-component mixture model formula
    model_formula <- bf(!!sym(col_name) ~ 1, family = mixture(gaussian, gaussian))
    
    # Fit the model using brms
    # Using 'silent = 2' to suppress Stan compilation messages
    fit <- brm(model_formula, 
               data = analysis_data, 
               chains = 2, iter = 1000, warmup = 300,
               silent = 2, refresh = 0)
    
    # Extract posterior predictions (gives probabilities for each component)
    preds <- posterior_epred(fit)
    
    # Identify the "high titre" component (the one with the larger mean)
    mu <- as.data.frame(fit)$b_mu2 # mu1 is fixed to 0
    pos_comp <- if (mean(mu) > 0) 2 else 1
    
    # Get the mean posterior probability of belonging to the positive component
    p_soft <- colMeans(preds[,,pos_comp])
    
    posterior_probs[[glue("{antigen}_sero_soft")]] <- p_soft
  }
  
  # Combine the per-antigen probabilities using the mean
  prob_cols_to_combine <- glue("{antigens}_sero_soft")
  p_soft_combined <- rowMeans(posterior_probs[prob_cols_to_combine], na.rm = TRUE)
  
  # Create the final phenotype tibble for this pathogen
  tibble(
    FID = analysis_data$FID,
    IID = analysis_data$IID,
    !!glue("{tolower(pathogen_name)}_sero_soft") := p_soft_combined,
    !!glue("{tolower(pathogen_name)}_sero_hard") := if_else(p_soft_combined >= 0.5, 1, 0)
  )
})

# --- 5. FINALIZE AND SAVE OUTPUTS --------------------------------------------
cat("\n--- Finalizing and saving outputs ---\n")

# Filter out any NULL results from skipped pathogens and merge
valid_phenos <- compact(pathogen_pheno_list)
if(length(valid_phenos) == 0) stop("No pathogens were successfully modeled.")
final_pheno_tbl <- reduce(valid_phenos, full_join, by = c("FID", "IID"))

# Write the final phenotype file
write_tsv(final_pheno_tbl, "quickdraws_input/phenotypes_pathogen_bayesian.tsv")
cat("  -> Wrote Bayesian phenotypes to quickdraws_input/phenotypes_pathogen_bayesian.tsv\n")

cat("\n--- Bayesian Script 05 complete ---\n") 