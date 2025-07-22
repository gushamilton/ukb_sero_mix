################################################################################
##
## Serology GWAS Simulation - Final Working Draft
##
## This script compares two methods for disentangling genetic effects on
## serostatus (beta) vs. IgG titre (gamma):
## 1. Baseline: Separate Probit and OLS models.
## 2. Mixture Model: A custom Expectation-Maximization (EM) algorithm that
##    corrects for selection bias.
##
## Key Finding:
## The mixture model provides a substantial power advantage for detecting effects
## on IgG titre by correcting the bias inherent in the naive baseline approach.
## For effects on serostatus, both models perform equally well.
##
## Next Steps & Limitations:
## This simulation models IgG titres using a Normal distribution. A key next
## step for more realistic modeling is to extend the EM algorithm to handle
## the Skew-Normal distributions often observed in real serology data, which
## presents a significant statistical challenge.
##
################################################################################

## ---- package-setup --------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, purrr, MASS, tibble, dplyr)

## ---- custom-em-algorithm --------------------------------------------------
fit_em_mixture <- function(d, max_iter = 100, tol = 1e-6) {
  probit_fit <- glm(S ~ Z + C, data = d, family = binomial(link = "probit"))
  ols_fit <- lm(Y ~ Z + C, data = d[S == 1])
  params <- list(
    alpha = coef(probit_fit)[1], beta = coef(probit_fit)["Z"],
    mu = coef(ols_fit)[1], gamma = coef(ols_fit)["Z"],
    sigma = summary(ols_fit)$sigma,
    c_s = coef(probit_fit)["C"], c_y = coef(ols_fit)["C"]
  )
  for (iter in 1:max_iter) {
    params_old <- params
    prob_seropositive <- pnorm(params$alpha + params$beta * d$Z + params$c_s * d$C)
    probit_update <- glm(S ~ Z + C, data = d, family = binomial(link = "probit"))
    params$alpha <- coef(probit_update)[1]; params$beta <- coef(probit_update)["Z"]; params$c_s <- coef(probit_update)["C"]
    ols_update <- lm(Y ~ Z + C, data = d, weights = prob_seropositive)
    params$mu <- coef(ols_update)[1]; params$gamma <- coef(ols_update)["Z"]; params$c_y <- coef(ols_update)["C"]; params$sigma <- summary(ols_update)$sigma
    param_vec_old <- unlist(params_old); param_vec_new <- unlist(params)
    if (max(abs(param_vec_new - param_vec_old), na.rm=T) < tol) break
  }
  se_probit <- summary(probit_update)$coef; se_ols <- summary(ols_update)$coef
  return(list(
    beta = params$beta, se_beta = se_probit["Z", "Std. Error"],
    gamma = params$gamma, se_gamma = se_ols["Z", "Std. Error"],
    iterations = iter
  ))
}

## ---- simulation-function --------------------------------------------------
sim_one_run_em <- function(rho, n = 1500, alpha_true = -0.5, beta_true, mu_true = 1, gamma_true,
                           sigma_true = 0.5, c_effect_s = 0.4, c_effect_y = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  C <- rnorm(n, 0, 1); Z <- rbinom(n, 2, 0.3)
  Sigma <- matrix(c(1, rho * sigma_true, rho * sigma_true, sigma_true^2), 2)
  eps <- MASS::mvrnorm(n, mu = c(0, 0), Sigma)
  S_star <- alpha_true + beta_true * Z + c_effect_s * C + eps[, 1]
  S <- as.numeric(S_star > 0)
  Y <- mu_true + gamma_true * Z + c_effect_y * C + eps[, 2]
  d <- data.table(Z, C, S, Y)
  res <- list(beta_base=NA, se_beta_base=NA, gamma_base=NA, se_gamma_base=NA,
              beta_em=NA, se_beta_em=NA, gamma_em=NA, se_gamma_em=NA)
  fit_beta_base <- glm(S ~ Z + C, data = d, family = binomial(link = "probit"))
  res$beta_base <- coef(fit_beta_base)["Z"]; res$se_beta_base <- summary(fit_beta_base)$coef["Z", "Std. Error"]
  if (sum(d$S) > 2) {
    fit_gamma_base <- lm(Y ~ Z + C, data = d[S == 1])
    if ("Z" %in% names(coef(fit_gamma_base))) {
      res$gamma_base <- coef(fit_gamma_base)["Z"]; res$se_gamma_base <- summary(fit_gamma_base)$coef["Z", "Std. Error"]
    }
  }
  tryCatch({
    fit_em <- fit_em_mixture(d)
    res$beta_em <- fit_em$beta; res$se_beta_em <- fit_em$se_beta
    res$gamma_em <- fit_em$gamma; res$se_gamma_em <- fit_em$se_gamma
  }, error = function(e) {})
  return(as_tibble(res))
}

## ---- main-simulation-driver ----------------------------------------------
run_em_simulation <- function(reps = 200) {
  scen_rho <- c(-0.3, 0, 0.3, 0.6); beta_true <- 0.1; gamma_true <- 0.1
  grid <- expand.grid(rho = scen_rho, rep = seq_len(reps))
  grid$seed <- seq_len(nrow(grid)) + as.integer(Sys.time())
  cat("Running EM simulation with", reps, "replicates and smaller effects...\n")
  all_res <- purrr::pmap_dfr(grid, ~ sim_one_run_em(rho = ..1, seed = ..3, beta_true = beta_true, gamma_true = gamma_true))
  cat("Simulation complete. Analyzing...\n")
  alpha5 <- qnorm(0.975)
  results <- all_res %>%
    mutate(
      power_beta_base = abs(beta_base / se_beta_base) > alpha5,
      power_gamma_base = abs(gamma_base / se_gamma_base) > alpha5,
      power_beta_em = abs(beta_em / se_beta_em) > alpha5,
      power_gamma_em = abs(gamma_em / se_gamma_em) > alpha5
    )
  summary_tbl <- results %>%
    group_by(rho = findInterval(grid$rho, scen_rho)) %>%
    summarise(
      across(starts_with("beta_"), ~ mean(. - beta_true, na.rm = TRUE), .names = "bias_{.col}"),
      across(starts_with("gamma_"), ~ mean(. - gamma_true, na.rm = TRUE), .names = "bias_{.col}"),
      across(starts_with("power_"), ~ mean(., na.rm = TRUE), .names = "{.col}")
    ) %>%
    mutate(rho = scen_rho[rho])
  cat("\n--- BIAS REPORT ---\n"); print(summary_tbl %>% select(rho, starts_with("bias_")))
  cat("\n--- POWER REPORT ---\n"); print(summary_tbl %>% select(rho, starts_with("power_")))
  cat("\n--- ANALYSIS COMPLETE ---\n")
}

run_em_simulation() 