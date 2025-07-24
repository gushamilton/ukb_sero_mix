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
pacman::p_load(tidyverse, glue, data.table, here, binom, cowplot, gridExtra)

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
  CT    = c("ct_pgp3",  "ct_mompa",  "ct_mompd", "ct_tarpf1", "ct_tarpf2"),
  HP    = c("hp_caga",  "hp_vaca",  "hp_omp", "hp_groel", "hp_catalase", "hp_urea")
)

# --- 4. COMBINATION LOGIC -------------------------------------------------------

#' Combine per-antigen probabilities using the mean.
#' This is a more conservative approach than Noisy-OR.
#' @param p_vector A numeric vector of per-antigen probabilities for one person.
#' @return A single combined probability (the mean).
combine_probs <- function(p_vector) {
  mean(p_vector, na.rm = TRUE)
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
  
  # Apply the combination function row-wise
  p_soft_combined <- apply(prob_matrix, 1, combine_probs)
  
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

# --- 8. DIAGNOSTICS: CORRELATIONS AND PLOTS -----------------------------------
cat("\n--- Generating diagnostic correlations and plots ---\n")

# To apply UKB rules, we need the original single-antigen hard calls
master_antigen_hard_calls <- master_pheno %>%
  select(FID, IID, any_of(glue("{pathogen_map %>% unlist %>% unique}_sero_hard")))

# Join back to get all necessary columns in one table
diag_tbl <- left_join(pathogen_pheno_tbl, master_antigen_hard_calls, by = c("FID", "IID"))


# --- 8a. UKB Rule-Based Serostatus Definition ---
ukb_rules <- list(
  EBV = list(type = "count", n = 2),
  CMV = list(type = "count", n = 2),
  HP  = list(type = "count", n = 2),
  CT  = list(type = "pgp3_positive", antigen = "ct_pgp3_sero_hard")
)

for (pathogen_name in names(pathogen_map)) {
    if (pathogen_name %in% names(ukb_rules)) {
        rule <- ukb_rules[[pathogen_name]]
        p_lower <- tolower(pathogen_name)
        ukb_col_name <- glue("{p_lower}_sero_ukb_rule")
        
        if (rule$type == "count") {
            antigens <- pathogen_map[[pathogen_name]]
            hard_cols <- glue("{antigens}_sero_hard")
            cols_to_use <- intersect(hard_cols, names(diag_tbl))
            if(length(cols_to_use) > 0) {
                antigen_sum <- rowSums(diag_tbl[cols_to_use], na.rm = TRUE)
                diag_tbl[[ukb_col_name]] <- ifelse(antigen_sum >= rule$n, 1, 0)
            }
        } else if (rule$type == "pgp3_positive") {
            if (rule$antigen %in% names(diag_tbl)) {
                diag_tbl[[ukb_col_name]] <- diag_tbl[[rule$antigen]]
            }
        }
    }
}


# --- 8b. Correlation Matrices ---
cat("\nCorrelations between combined probability and constituent antigen probabilities:\n")
for (pathogen_name in successful_pathogens) {
    p_upper <- toupper(pathogen_name)
    
    combined_col <- glue("{pathogen_name}_sero_soft")
    antigen_cols <- glue("{pathogen_map[[p_upper]]}_sero_soft")
    
    # Need to join to get single-antigen soft probs
    cols_to_cor_data <- master_pheno %>%
        select(any_of(c("FID", "IID", antigen_cols))) %>%
        left_join(pathogen_pheno_tbl %>% select(FID, IID, !!sym(combined_col)), by = c("FID", "IID"))

    cols_to_cor <- intersect(c(combined_col, antigen_cols), names(cols_to_cor_data))
    
    if (length(cols_to_cor) > 1) {
        cat("\n---", p_upper, "---\n")
        cor_matrix <- cor(cols_to_cor_data[cols_to_cor], use = "pairwise.complete.obs")
        print(round(cor_matrix, 3))
    }
}


# --- 8c. Scatter Plots (Combined vs. Individual) ---
cat("\nGenerating scatter plots (Combined vs. Individual Antigens)...\n")
scatter_plot_list <- map(successful_pathogens, function(p_name_lower) {
    p_upper <- toupper(p_name_lower)
    combined_col <- glue("{p_name_lower}_sero_soft")
    antigen_cols <- glue("{pathogen_map[[p_upper]]}_sero_soft")
    
    cols_to_plot <- intersect(antigen_cols, names(master_pheno))
    
    plot_data <- master_pheno %>%
        select(FID, IID, all_of(cols_to_plot)) %>%
        left_join(diag_tbl %>% select(FID, IID, !!sym(combined_col)), by = c("FID", "IID"))

    if (length(cols_to_plot) > 0) {
        plot_data %>%
            pivot_longer(cols = all_of(cols_to_plot), names_to = "antigen", values_to = "antigen_prob") %>%
            ggplot(aes(x = antigen_prob, y = .data[[combined_col]])) +
            geom_point(alpha = 0.05, color = "darkblue", shape = 16) +
            facet_wrap(~antigen, ncol = 3) +
            labs(title = glue("Combined vs. Individual Antigen Probabilities for {p_upper}"),
                 x = "Individual Antigen 'sero_soft' Probability",
                 y = "Combined (Mean) 'sero_soft' Probability") +
            theme_bw(base_size = 14)
    } else { NULL }
})

ggsave("quickdraws_input/pathogen_vs_antigen_scatters.pdf", 
       marrangeGrob(grobs = compact(scatter_plot_list), nrow=1, ncol=1), 
       width = 12, height = 8)
cat("  -> Saved scatter plots to quickdraws_input/pathogen_vs_antigen_scatters.pdf\n")


# --- 8d. Violin Plots (Soft vs. Hard Calls) ---
cat("\nGenerating violin plots (Soft vs. Hard Calls)...\n")
comparison_plot_list <- map(successful_pathogens, function(p_name_lower) {
    p_upper <- toupper(p_name_lower)
    soft_col <- glue("{p_name_lower}_sero_soft")
    hard_col <- glue("{p_name_lower}_sero_hard")
    ukb_col  <- glue("{p_name_lower}_sero_ukb_rule")

    if(!ukb_col %in% names(diag_tbl)) return(NULL)

    plot_data <- diag_tbl %>%
      select(all_of(c(soft_col, hard_col, ukb_col))) %>%
      rename(Soft = 1, `Mean Prob >= 0.5` = 2, `UKB Rule` = 3) %>%
      pivot_longer(-Soft, names_to = "method", values_to = "status")

    ggplot(plot_data, aes(x = as.factor(status), y = Soft, fill = method)) +
        geom_violin(trim = FALSE, alpha = 0.7) +
        geom_boxplot(width=0.1, fill="white", alpha=0.7) +
        facet_wrap(~method, scales = "free_x") +
        labs(title = glue("Soft Probability Distribution by Hard Call Method for {p_upper}"),
             x = "Hard Serostatus (0 or 1)",
             y = "Combined 'sero_soft' Probability") +
        theme_light(base_size = 14) +
        theme(legend.position = "none", strip.text = element_text(face="bold"))
})

ggsave("quickdraws_input/pathogen_soft_vs_hard_comparison.pdf", 
       marrangeGrob(grobs = compact(comparison_plot_list), nrow=1, ncol=1), 
       width = 10, height = 6)
cat("  -> Saved comparison plots to quickdraws_input/pathogen_soft_vs_hard_comparison.pdf\n")


# --- 9. PLOT DISTRIBUTIONS ----------------------------------------------------
cat("\n--- Generating distribution plots ---\n")

plot_list <- map(successful_pathogens, function(p_name_lower) {
  soft_col <- glue("{p_name_lower}_sero_soft")
  
  ggplot(pathogen_pheno_tbl, aes(x = .data[[soft_col]])) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.8) +
    labs(
      title = glue("Distribution of Combined Probability for {toupper(p_name_lower)}"),
      x = "Combined 'sero_soft' Probability",
      y = "Count"
    ) +
    theme_light()
})

# Arrange plots into a grid and save
n_pathogens <- length(plot_list)
n_cols <- if (n_pathogens <= 4) 2 else 3
plot_grid <- cowplot::plot_grid(plotlist = plot_list, ncol = n_cols)

ggsave("quickdraws_input/pathogen_probability_distributions.pdf", plot_grid, 
       width = n_cols * 5, height = ceiling(n_pathogens / n_cols) * 4)
cat("  -> Wrote distribution plots to quickdraws_input/pathogen_probability_distributions.pdf\n")


cat("\n--- Script 04 complete ---\n") 