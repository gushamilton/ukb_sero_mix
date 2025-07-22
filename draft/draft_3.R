################################################################################
##
## Serology GWAS Simulation - Draft 3: Investigating Distributional Separation
##
################################################################################

## ---- package-setup --------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, purrr, MASS, tibble, dplyr)

## ---- data-generation-v3 ---------------------------------------------------
generate_data_v3 <- function(n = 1500, beta_true = 0.1, gamma_true = 0.1,
                             separation = 2.0, rho = 0,
                             c_effect_s = 0.4, c_effect_y = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  C <- rnorm(n, 0, 1); Z <- rbinom(n, 2, 0.3)
  
  # Latent serostatus is a function of Z and C
  log_odds <- -0.5 + beta_true * Z + c_effect_s * C
  S_true <- rbinom(n, 1, plogis(log_odds))
  
  # Simulate correlated errors for the underlying continuous variables
  sigma_s <- 1; sigma_y <- 0.5
  cov_matrix <- matrix(c(sigma_s^2, rho * sigma_s * sigma_y, rho * sigma_s * sigma_y, sigma_y^2), 2)
  errors <- MASS::mvrnorm(n, mu = c(0, 0), cov_matrix)
  
  # Generate Y from two different Normal distributions
  Y <- numeric(n)
  Y[S_true == 0] <- 0.5 + errors[S_true == 0, 2]
  Y[S_true == 1] <- separation + gamma_true * Z[S_true == 1] + c_effect_y * C[S_true == 1] + errors[S_true == 1, 2]
  
  return(data.table(Z, C, Y, S_true))
}

## ---- custom-em-algorithm --------------------------------------------------
fit_em_mixture_v3 <- function(d, max_iter = 100, tol = 1e-6) {
  # Initialize with kmeans guess
  initial_groups <- kmeans(d$Y, 2)$cluster
  d[, prob_pos := initial_groups - 1]
  
  params <- list() # Store params
  
  for (iter in 1:max_iter) {
    params_old <- params
    
    # E-Step: Update soft probabilities
    fit_prob <- glm(prob_pos ~ Z + C, data=d, family=binomial(link="probit"))
    prob_seropositive <- predict(fit_prob, type = "response")
    
    dens_neg <- dnorm(d$Y, mean = mean(d$Y[d$prob_pos < 0.5]), sd = sd(d$Y[d$prob_pos < 0.5]))
    dens_pos <- dnorm(d$Y, mean = mean(d$Y[d$prob_pos > 0.5]), sd = sd(d$Y[d$prob_pos > 0.5]))
    
    d$prob_pos <- (prob_seropositive * dens_pos) / (prob_seropositive * dens_pos + (1-prob_seropositive) * dens_neg)
    d$prob_pos[is.na(d$prob_pos)] <- 0.5

    # M-Step: Update model parameters
    probit_update <- glm(prob_pos ~ Z + C, data = d, family = binomial(link = "probit"))
    ols_update <- lm(Y ~ Z + C, data = d, weights = d$prob_pos)
    
    params$beta <- coef(probit_update)["Z"]; params$se_beta <- summary(probit_update)$coef["Z", "Std. Error"]
    params$gamma <- coef(ols_update)["Z"]; params$se_gamma <- summary(ols_update)$coef["Z", "Std. Error"]
    
    if (iter > 1 && max(abs(unlist(params) - unlist(params_old)), na.rm=T) < tol) break
  }
  return(params)
}

## ---- diagnostic-plots -----------------------------------------------------
generate_diagnostic_plots <- function(separation_levels) {
  cat("Generating diagnostic plots to 'distribution_plots.pdf'...\n")
  pdf("distribution_plots.pdf", width = 10, height = 4); par(mfrow = c(1, 3))
  for (sep in separation_levels) {
    d <- generate_data_v3(separation = sep, seed = 123)
    plot(density(d$Y[d$S_true==0]), col='blue', xlim=range(d$Y), main=paste("Separation =", sep), xlab="IgG Titre (Y)", ylim=c(0,1.2))
    lines(density(d$Y[d$S_true==1]), col='red')
    legend("topright", legend=c("Seronegative", "Seropositive"), col=c("blue", "red"), lty=1, bty="n")
  }
  dev.off()
}

## ---- main-simulation ------------------------------------------------------
run_simulation_v3 <- function(reps = 50) {
  scen_rho <- c(-0.3, 0, 0.3, 0.6)
  scen_separation <- c(1.5, 2.0, 3.0)
  beta_true <- 0.1; gamma_true <- 0.1
  
  generate_diagnostic_plots(scen_separation)
  
  grid <- expand.grid(rho = scen_rho, separation = scen_separation, rep = seq_len(reps))
  
  cat("Running simulation for", nrow(grid), "scenarios...\n")
  
  all_res <- pmap_dfr(list(i = 1:nrow(grid)), function(i) {
    d <- generate_data_v3(rho = grid$rho[i], separation = grid$separation[i], seed=i)
    
    # Baseline Model
    beta_base <- coef(glm(S_true ~ Z + C, data = d, family = binomial(link = "probit")))["Z"]
    gamma_base <- NA
    if(sum(d$S_true) > 2) gamma_base <- coef(lm(Y ~ Z + C, data = d[S_true == 1]))["Z"]
    
    # EM Mixture Model
    em_fit <- tryCatch(fit_em_mixture_v3(d), error=function(e) NULL)
    
    tibble(rho=grid$rho[i], separation=grid$separation[i], 
           beta_base=beta_base, gamma_base=gamma_base, 
           beta_em=em_fit$beta, gamma_em=em_fit$gamma)
  })
  
  # Analysis
  summary_tbl <- all_res %>%
    group_by(rho, separation) %>%
    summarise(
      bias_beta_base = mean(beta_base - beta_true, na.rm=T),
      bias_gamma_base = mean(gamma_base - gamma_true, na.rm=T),
      bias_beta_em = mean(beta_em - beta_true, na.rm=T),
      bias_gamma_em = mean(gamma_em - gamma_true, na.rm=T),
      .groups = 'drop'
    )
    
  cat("\n--- BIAS REPORT ---\n")
  print(as.data.frame(summary_tbl))
}

run_simulation_v3() 