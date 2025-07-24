################################################################################
#   UKB SEROLOGY: PATHOGEN-LEVEL SOFT PHENOTYPES (SCRIPT 04) - NOISY-OR
################################################################################
#  ∙ Combines multiple per-antigen soft probabilities (from Script 02) into a
#    single pathogen-level probability using a "Noisy-OR" model.
#  ∙ This approach avoids fitting a new multivariate mixture model, making it
#    more robust and sidestepping potential convergence issues with mclust.
#  ∙ It assumes conditional independence of antigens given the true latent status.
#  ∙ Outputs (per pathogen):
#      ‑ <pathogen>_sero_soft   – Combined P(infected)
#      ‑ <pathogen>_sero_hard   – MAP call (prob ≥ 0.5)
#      ‑ <pathogen>_w_soft      – weight = p_sero_soft
#      ‑ <pathogen>_w_hard      – weight = 2 × |p – 0.5|
#  ∙ Writes phenotype & weight files into quickdraws_input/
################################################################################

# --- 1. SETUP & PACKAGES -------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, glue, data.table, here, binom)

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
  CT    = c("ct_pgp3",  "ct_mompa",  "ct_mompd", "ct_tarpf1", "ct_tarpf2"),
  HP    = c("hp_caga",  "hp_vaca",  "hp_omp", "hp_groel", "hp_catalase", "hp_urea")
)

# --- 4. COMBINATION LOGIC -------------------------------------------------------

#' Combine per-antigen probabilities using a Noisy-OR model.
#' P(Pathogen+) = 1 - P(all antigens negative) = 1 - product(1 - P(antigen_i positive))
#' @param p_vector A numeric vector of per-antigen probabilities for one person.
#' @return A single combined probability.
noisy_or <- function(p_vector) {
  1 - prod(1 - p_vector, na.rm = TRUE)
}

# --- 5. MAIN LOOP --------------------------------------------------------------
cat("\n--- Starting pathogen-level probability combination ---\n")

pheno_tibbles <- map(names(pathogen_map), function(pathogen_name) {
  cat("-- Processing", pathogen_name, "...\n")
  
  antigens <- pathogen_map[[pathogen_name]]
  soft_cols <- glue("{antigens}_sero_soft")
  
  # Keep only antigens present in the master phenotype table
  cols_to_use <- intersect(soft_cols, names(master_pheno))
  
  if (length(cols_to_use) < 1) {
    warning("  -> Skipping ", pathogen_name, ": no 'sero_soft' columns found.")
    return(NULL)
  }
  
  # Extract the matrix of per-antigen probabilities
  prob_matrix <- as.matrix(master_pheno[cols_to_use])
  
  # Apply the noisy-OR function row-wise
  p_soft_combined <- apply(prob_matrix, 1, noisy_or)
  
  # Create the new phenotype tibble
  tibble(
    FID = master_pheno$FID,
    IID = master_pheno$IID,
    !!glue("{tolower(pathogen_name)}_sero_soft") := p_soft_combined,
    !!glue("{tolower(pathogen_name)}_sero_hard") := if_else(p_soft_combined >= 0.5, 1, 0),
    !!glue("{tolower(pathogen_name)}_w_soft")    := p_soft_combined,
    !!glue("{tolower(pathogen_name)}_w_hard")    := 2 * abs(p_soft_combined - 0.5)
  )
})

# Filter out any NULL results from skipped pathogens
pheno_tibbles <- compact(pheno_tibbles)

if (length(pheno_tibbles) == 0) {
  stop("No pathogen phenotypes were generated.")
}

# --- 6. MERGE & WRITE OUTPUTS --------------------------------------------------
cat("\n--- Writing phenotype & weight files ---\n")

# Merge all pathogen phenotypes into a single table
pathogen_pheno_tbl <- reduce(pheno_tibbles, left_join, by = c("FID", "IID"))

# Write combined phenotype table
write_tsv(pathogen_pheno_tbl, "quickdraws_input/phenotypes_pathogen.tsv")
cat("  -> Wrote quickdraws_input/phenotypes_pathogen.tsv\n")

# Write individual weight files (soft + hard) for each pathogen
for (pathogen_tibble in pheno_tibbles) {
  pathogen_name_lower <- str_remove(names(pathogen_tibble)[3], "_sero_soft")
  
  soft_col <- glue("{pathogen_name_lower}_w_soft")
  hard_col <- glue("{pathogen_name_lower}_w_hard")

  if (all(c(soft_col, hard_col) %in% names(pathogen_tibble))) {
    soft_w <- pathogen_tibble %>% select(FID, IID, Weight = !!sym(soft_col))
    write_tsv(soft_w, glue("quickdraws_input/{pathogen_name_lower}_p_soft.weights"))

    hard_w <- pathogen_tibble %>% select(FID, IID, Weight = !!sym(hard_col))
    write_tsv(hard_w, glue("quickdraws_input/{pathogen_name_lower}_sero_hard.weights"))
  }
}
cat("  -> Wrote per-pathogen weight files.\n")


# --- 7. SEROPREVALENCE SUMMARY ----------------------------------------------
cat("\n--- Seroprevalence estimates (p_hard with 95% CI) ---\n")

# Use a loop to calculate seroprevalence and CI for each pathogen
successful_pathogens <- map_chr(pheno_tibbles, ~str_remove(names(.x)[3], "_sero_soft"))

seroprev_list <- map(successful_pathogens, function(p_name_lower) {
  hard_call_col <- glue("{p_name_lower}_sero_hard")
  
  counts <- table(pathogen_pheno_tbl[[hard_call_col]])
  k <- if ("1" %in% names(counts)) counts["1"] else 0
  n <- sum(counts)
  ci <- binom.confint(k, n, methods = "wilson")
  
  tibble(
    pathogen = toupper(p_name_lower),
    n_total = n,
    n_positive = k,
    seroprevalence_pct = round(ci$mean * 100, 2),
    ci_lower_pct = round(ci$lower * 100, 2),
    ci_upper_pct = round(ci$upper * 100, 2)
  )
})

seroprev_tbl <- bind_rows(seroprev_list)

print(seroprev_tbl)

write_tsv(seroprev_tbl, "quickdraws_input/pathogen_seroprevalence.tsv")
cat("  -> Wrote seroprevalence summary to quickdraws_input/pathogen_seroprevalence.tsv\n")

cat("\n--- Script 04 complete ---\n") 