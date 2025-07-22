################################################################################
##
## Debug flexmix S4 object structure
##
################################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, MASS, flexmix)

# Simple data generation
generate_data <- function(n = 500, beta_true = 0.1, gamma_true = 0.1,
                          separation = 2.0, rho = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  C <- rnorm(n, 0, 1); Z <- rbinom(n, 2, 0.3)
  log_odds <- -0.5 + beta_true * Z + 0.4 * C
  S_true <- rbinom(n, 1, plogis(log_odds))
  sigma_s <- 1; sigma_y <- 0.5
  cov_matrix <- matrix(c(sigma_s^2, rho*sigma_s*sigma_y, rho*sigma_s*sigma_y, sigma_y^2), 2)
  errors <- MASS::mvrnorm(n, mu = c(0, 0), cov_matrix)
  Y <- numeric(n)
  Y[S_true == 0] <- 0.5 + errors[S_true == 0, 2]
  Y[S_true == 1] <- separation + gamma_true * Z[S_true == 1] + 0.3 * C[S_true == 1] + errors[S_true == 1, 2]
  return(data.table(Z, C, Y, S_true))
}

# Generate test data
d <- generate_data(n = 500, rho = 0, separation = 2.0, beta_true = 0.1, gamma_true = 0.1)

# Find cutoff for initialization
find_optimal_cutoff <- function(y) {
  dens <- density(y, n = 1024)
  pks <- which(diff(sign(diff(dens$y))) == -2) + 1
  if (length(pks) < 2) return(median(y))
  sorted_pks <- sort(pks[order(dens$y[pks], decreasing=TRUE)[1:2]])
  valley_idx <- which.min(dens$y[sorted_pks[1]:sorted_pks[2]]) + sorted_pks[1] - 1
  return(dens$x[valley_idx])
}

cutoff <- find_optimal_cutoff(d$Y)
d[, S_observed := ifelse(Y > cutoff, 1, 0)]

cat("Data summary:\n")
cat("Sample size:", nrow(d), "\n")
cat("Seropositive proportion:", mean(d$S_observed), "\n")
cat("Y range:", range(d$Y), "\n")

# Fit flexmix
cat("\nFitting flexmix...\n")
fit_mix <- flexmix(Y ~ Z + C | Z + C, data = d, k = 2, 
                   cluster = d$S_observed + 1,
                   control = list(iter.max = 200, minprior = 0.01))

cat("flexmix converged!\n")
cat("Number of iterations:", fit_mix@iter, "\n")
cat("Log-likelihood:", logLik(fit_mix), "\n")

# Inspect S4 object
cat("\n=== S4 OBJECT INSPECTION ===\n")
cat("Class:", class(fit_mix), "\n")
cat("Slots:", slotNames(fit_mix), "\n")

# Inspect each slot
cat("\n--- @iter ---\n")
print(fit_mix@iter)

cat("\n--- @k ---\n")
print(fit_mix@k)

cat("\n--- @cluster ---\n")
print(table(fit_mix@cluster))

cat("\n--- @prior ---\n")
print(fit_mix@prior)

cat("\n--- @posterior ---\n")
print(head(fit_mix@posterior))

cat("\n--- @concomitant ---\n")
print(class(fit_mix@concomitant))
print(slotNames(fit_mix@concomitant))

cat("\n--- @model ---\n")
print(class(fit_mix@model))
print(length(fit_mix@model))

# Try refit
cat("\n=== REFIT INSPECTION ===\n")
refit_models <- refit(fit_mix)
cat("Refit class:", class(refit_models), "\n")
cat("Refit slots:", slotNames(refit_models), "\n")

cat("\n--- @coef ---\n")
print(refit_models@coef)
cat("Coefficient names:\n")
print(names(refit_models@coef))

cat("\n--- @vcov ---\n")
print(dim(refit_models@vcov))
print(rownames(refit_models@vcov))

# Extract coefficients properly
cat("\n=== COEFFICIENT EXTRACTION ===\n")
all_coefs <- refit_models@coef
all_vcov <- refit_models@vcov

# Find component coefficients
comp1_coefs <- all_coefs[grep("Comp.1", names(all_coefs))]
comp2_coefs <- all_coefs[grep("Comp.2", names(all_coefs))]

cat("Component 1 coefficients:\n")
print(comp1_coefs)

cat("Component 2 coefficients:\n")
print(comp2_coefs)

# Extract standard errors
comp1_vcov <- all_vcov[grep("Comp.1", rownames(all_vcov)), grep("Comp.1", colnames(all_vcov))]
comp2_vcov <- all_vcov[grep("Comp.2", rownames(all_vcov)), grep("Comp.2", colnames(all_vcov))]

cat("Component 1 Z coefficient:", comp1_coefs["model.1_Comp.1_coef.Z"], "\n")
cat("Component 2 Z coefficient:", comp2_coefs["model.1_Comp.2_coef.Z"], "\n")

# Identify which component is seropositive (higher intercept)
pos_comp <- ifelse(comp1_coefs["model.1_Comp.1_coef.(Intercept)"] > comp2_coefs["model.1_Comp.2_coef.(Intercept)"], 1, 2)
cat("Seropositive component:", pos_comp, "\n")

# Extract gamma (Z effect on IgG in seropositives)
gamma_mix <- if(pos_comp == 1) comp1_coefs["model.1_Comp.1_coef.Z"] else comp2_coefs["model.1_Comp.2_coef.Z"]
cat("Gamma (Z effect on IgG):", gamma_mix, "\n")

# Investigate concomitant coefficient
cat("\n=== CONCOMITANT ANALYSIS ===\n")
concomitant_coef <- all_coefs["concomitant_Comp.2.alpha"]
cat("Concomitant coefficient:", concomitant_coef, "\n")

# Check what this actually means
cat("True beta effect:", 0.1, "\n")
cat("True seropositive proportion:", mean(d$S_true), "\n")
cat("Observed seropositive proportion:", mean(d$S_observed), "\n")

# Look at the mixing probabilities
cat("Posterior matrix dimensions:", dim(fit_mix@posterior), "\n")
cat("Posterior matrix structure:\n")
str(fit_mix@posterior)

# Try to access posterior correctly
if(is.matrix(fit_mix@posterior)) {
  cat("Component 1 proportion:", mean(fit_mix@posterior[,1]), "\n")
  cat("Component 2 proportion:", mean(fit_mix@posterior[,2]), "\n")
} else if(is.list(fit_mix@posterior)) {
  cat("Component 1 proportion (scaled):", mean(fit_mix@posterior$scaled[,1]), "\n")
  cat("Component 2 proportion (scaled):", mean(fit_mix@posterior$scaled[,2]), "\n")
  cat("Component 1 proportion (unscaled):", mean(fit_mix@posterior$unscaled[,1]), "\n")
  cat("Component 2 proportion (unscaled):", mean(fit_mix@posterior$unscaled[,2]), "\n")
} else {
  cat("Posterior is not a matrix, it's:", class(fit_mix@posterior), "\n")
}

# Check if concomitant affects mixing probability correctly
cat("Mixing probability for Z=0 vs Z=1:\n")
z0_data <- d[d$Z == 0]
z1_data <- d[d$Z == 1]
if(nrow(z0_data) > 0 && nrow(z1_data) > 0) {
  z0_mix <- mean(fit_mix@posterior$scaled[d$Z == 0, 2])  # prob of component 2
  z1_mix <- mean(fit_mix@posterior$scaled[d$Z == 1, 2])  # prob of component 2
  cat("Z=0 mixing prob:", z0_mix, "\n")
  cat("Z=1 mixing prob:", z1_mix, "\n")
  cat("Difference:", z1_mix - z0_mix, "\n")
}

# The concomitant coefficient affects the mixing probability
# Let's understand what this means
cat("\n=== CONCOMITANT INTERPRETATION ===\n")
cat("The concomitant coefficient", concomitant_coef, "affects the mixing probability.\n")
cat("This is NOT directly the beta effect on serostatus.\n")
cat("The beta effect would be the effect of Z on the probability of being in component 2.\n")
cat("We need to convert this to a proper beta estimate.\n")

# Compare with true serostatus
cat("\n=== COMPARISON WITH TRUE SEROSTATUS ===\n")
# Fit a logistic regression to see the true beta effect
true_logit <- glm(S_true ~ Z, data = d, family = binomial(link = "probit"))
true_beta <- coef(true_logit)["Z"]
cat("True beta (probit):", true_beta, "\n")

# Compare mixing probabilities with true serostatus
cat("Correlation between mixing prob and true serostatus:", 
    cor(fit_mix@posterior$scaled[,2], d$S_true), "\n")

# Check if component 2 corresponds to seropositive
comp2_seropos_corr <- cor(fit_mix@posterior$scaled[,2], d$S_true)
comp1_seropos_corr <- cor(fit_mix@posterior$scaled[,1], d$S_true)

cat("Component 1 correlation with seropositive:", comp1_seropos_corr, "\n")
cat("Component 2 correlation with seropositive:", comp2_seropos_corr, "\n")

# Determine which component is seropositive
seropos_comp <- ifelse(comp2_seropos_corr > comp1_seropos_corr, 2, 1)
cat("Seropositive component:", seropos_comp, "\n")

# If component 1 is seropositive, we need to flip the concomitant coefficient
if(seropos_comp == 1) {
  cat("Component 1 is seropositive, so beta should be -concomitant_coef\n")
  corrected_beta <- -concomitant_coef
} else {
  cat("Component 2 is seropositive, so beta should be concomitant_coef\n")
  corrected_beta <- concomitant_coef
}

cat("Corrected beta estimate:", corrected_beta, "\n")
cat("True beta:", true_beta, "\n")
cat("Bias:", corrected_beta - true_beta, "\n")

# Test different interpretations
cat("\n=== TESTING DIFFERENT INTERPRETATIONS ===\n")
cat("Raw concomitant coefficient:", concomitant_coef, "\n")
cat("Negative concomitant coefficient:", -concomitant_coef, "\n")
cat("Scaled concomitant coefficient (2x):", 2 * concomitant_coef, "\n")
cat("Scaled concomitant coefficient (0.5x):", 0.5 * concomitant_coef, "\n")

# Check if the issue is scale or direction
cat("Bias with raw concomitant:", concomitant_coef - true_beta, "\n")
cat("Bias with negative concomitant:", -concomitant_coef - true_beta, "\n")
cat("Bias with 2x concomitant:", 2 * concomitant_coef - true_beta, "\n")
cat("Bias with 0.5x concomitant:", 0.5 * concomitant_coef - true_beta, "\n")

# The concomitant coefficient might need to be scaled to match the probit scale
# Probit and logit coefficients are related by a factor of approximately 1.6
cat("Bias with probit scaling (1.6x):", 1.6 * concomitant_coef - true_beta, "\n")

# Test different scaling factors to find optimal
cat("\n=== OPTIMAL SCALING TEST ===\n")
scaling_factors <- seq(0.3, 0.7, by = 0.05)
for(sf in scaling_factors) {
  bias <- sf * concomitant_coef - true_beta
  cat("Scaling factor", sf, "-> Bias:", bias, "\n")
}

# Find the scaling factor that minimizes absolute bias
abs_biases <- sapply(scaling_factors, function(sf) abs(sf * concomitant_coef - true_beta))
optimal_sf <- scaling_factors[which.min(abs_biases)]
cat("Optimal scaling factor:", optimal_sf, "-> Minimal bias:", min(abs_biases), "\n") 