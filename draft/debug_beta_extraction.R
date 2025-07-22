# Debug script to test different beta extraction methods from flexmix
library(pacman)
p_load(data.table, purrr, MASS, tibble, dplyr, flexmix)

# Data generation function
generate_data <- function(n = 1000, beta = 0.1, gamma = 0.1, rho = 0.3, separation = 2.0) {
  # Generate covariates
  Z <- rbinom(n, 1, 0.5)  # Binary genetic variant
  C <- rnorm(n, 0, 1)     # Continuous covariate
  
  # Generate correlated errors
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  errors <- mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  epsilon_S <- errors[, 1]
  epsilon_Y <- errors[, 2]
  
  # Generate true serostatus (probit model)
  S_star <- 0.5 + beta * Z + 0.3 * C + epsilon_S
  S_true <- ifelse(S_star > 0, 1, 0)
  
  # Generate IgG titre (normal model for seropositives)
  mu_pos <- 3.0 + gamma * Z + 0.2 * C + separation  # Higher mean for seropositives
  mu_neg <- 1.0 + 0.1 * C                           # Lower mean for seronegatives
  
  Y <- ifelse(S_true == 1, 
              mu_pos + epsilon_Y,
              mu_neg + epsilon_Y)
  
  # Create observed serostatus based on cutoff
  cutoff <- quantile(Y, 0.6)  # 60th percentile as cutoff
  S_observed <- ifelse(Y > cutoff, 1, 0)
  
  data.table(Z = Z, C = C, Y = Y, S_true = S_true, S_observed = S_observed)
}

# Test different beta extraction methods
test_beta_extraction <- function() {
  cat("=== TESTING BETA EXTRACTION METHODS ===\n")
  
  # Generate data with known beta
  true_beta <- 0.1
  d <- generate_data(n = 1000, beta = true_beta, gamma = 0.1, rho = 0.3, separation = 2.0)
  
  cat("True beta:", true_beta, "\n")
  cat("True seropositive proportion:", mean(d$S_true), "\n")
  cat("Observed seropositive proportion:", mean(d$S_observed), "\n\n")
  
  # Method 1: True logistic regression (gold standard)
  true_logit <- glm(S_true ~ Z + C, data = d, family = binomial(link = "logit"))
  true_beta_logit <- coef(true_logit)["Z"]
  cat("Method 1 - True logistic (logit scale):", true_beta_logit, "\n")
  
  # Method 2: True probit regression
  true_probit <- glm(S_true ~ Z + C, data = d, family = binomial(link = "probit"))
  true_beta_probit <- coef(true_probit)["Z"]
  cat("Method 2 - True probit (probit scale):", true_beta_probit, "\n")
  
  # Method 3: Observed logistic regression
  obs_logit <- glm(S_observed ~ Z + C, data = d, family = binomial(link = "logit"))
  obs_beta_logit <- coef(obs_logit)["Z"]
  cat("Method 3 - Observed logistic (logit scale):", obs_beta_logit, "\n")
  
  # Method 4: Observed probit regression
  obs_probit <- glm(S_observed ~ Z + C, data = d, family = binomial(link = "probit"))
  obs_beta_probit <- coef(obs_probit)["Z"]
  cat("Method 4 - Observed probit (probit scale):", obs_beta_probit, "\n")
  
  # Method 5: Flexmix concomitant coefficient (raw)
  tryCatch({
    fit_mix <- flexmix(Y ~ Z + C | Z + C, data = d, k = 2,
                       cluster = d$S_observed + 1,
                       control = list(iter.max = 200, minprior = 0.01))
    
    refit_models <- refit(fit_mix)
    all_coefs <- refit_models@coef
    
    # Identify seropositive component
    comp1_intercept <- all_coefs["model.1_Comp.1_coef.(Intercept)"]
    comp2_intercept <- all_coefs["model.1_Comp.2_coef.(Intercept)"]
    pos_comp <- ifelse(comp1_intercept > comp2_intercept, 1, 2)
    
    concomitant_coef <- all_coefs["concomitant_Comp.2.alpha"]
    cat("Method 5 - Flexmix concomitant (raw):", concomitant_coef, "\n")
    
    # Method 6: Flexmix concomitant with probit scaling (1.6)
    cat("Method 6 - Flexmix concomitant (probit scaled):", concomitant_coef / 1.6, "\n")
    
    # Method 7: Flexmix concomitant with inverse scaling (0.625)
    cat("Method 7 - Flexmix concomitant (inverse scaled):", concomitant_coef * 0.625, "\n")
    
    # Method 8: Use posterior probabilities to fit weighted regression
    post_mat <- if (is.list(fit_mix@posterior)) fit_mix@posterior$scaled else fit_mix@posterior
    d[, p_seropos := post_mat[, pos_comp]]
    
    # Fit logistic regression on posterior probabilities
    post_logit <- glm(p_seropos ~ Z + C, data = d, family = binomial(link = "logit"))
    post_beta_logit <- coef(post_logit)["Z"]
    cat("Method 8 - Posterior-weighted logistic:", post_beta_logit, "\n")
    
    # Method 9: Use posterior probabilities to fit probit regression
    post_probit <- glm(p_seropos ~ Z + C, data = d, family = binomial(link = "probit"))
    post_beta_probit <- coef(post_probit)["Z"]
    cat("Method 9 - Posterior-weighted probit:", post_beta_probit, "\n")
    
    # Method 10: Use posterior probabilities as weights in logistic regression
    d[, w_seropos := p_seropos]
    weighted_logit <- glm(S_observed ~ Z + C, data = d, family = binomial(link = "logit"), weights = w_seropos)
    weighted_beta_logit <- coef(weighted_logit)["Z"]
    cat("Method 10 - Weighted logistic (posterior as weights):", weighted_beta_logit, "\n")
    
    # Method 11: Use posterior probabilities as weights in probit regression
    weighted_probit <- glm(S_observed ~ Z + C, data = d, family = binomial(link = "probit"), weights = w_seropos)
    weighted_beta_probit <- coef(weighted_probit)["Z"]
    cat("Method 11 - Weighted probit (posterior as weights):", weighted_beta_probit, "\n")
    
    # Calculate biases
    cat("\n=== BIAS ANALYSIS ===\n")
    cat("True beta (probit):", true_beta, "\n")
    cat("Bias for Method 5 (raw concomitant):", concomitant_coef - true_beta, "\n")
    cat("Bias for Method 6 (probit scaled):", (concomitant_coef / 1.6) - true_beta, "\n")
    cat("Bias for Method 7 (inverse scaled):", (concomitant_coef * 0.625) - true_beta, "\n")
    cat("Bias for Method 8 (posterior logistic):", post_beta_logit - true_beta, "\n")
    cat("Bias for Method 9 (posterior probit):", post_beta_probit - true_beta, "\n")
    cat("Bias for Method 10 (weighted logistic):", weighted_beta_logit - true_beta, "\n")
    cat("Bias for Method 11 (weighted probit):", weighted_beta_probit - true_beta, "\n")
    
    # Find best method
    biases <- c(
      "raw" = concomitant_coef - true_beta,
      "probit_scaled" = (concomitant_coef / 1.6) - true_beta,
      "inverse_scaled" = (concomitant_coef * 0.625) - true_beta,
      "posterior_logit" = post_beta_logit - true_beta,
      "posterior_probit" = post_beta_probit - true_beta,
      "weighted_logit" = weighted_beta_logit - true_beta,
      "weighted_probit" = weighted_beta_probit - true_beta
    )
    
    best_method <- names(biases)[which.min(abs(biases))]
    cat("\nBest method (lowest absolute bias):", best_method, "-> Bias:", biases[best_method], "\n")
    
  }, error = function(e) {
    cat("Flexmix failed:", e$message, "\n")
  })
}

# Run the test
test_beta_extraction() 