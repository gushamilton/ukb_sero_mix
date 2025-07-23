################################################################################
#   ANALYSE A SINGLE VECTOR OF ANTIBODY TITRES — 2025‑07‑21
################################################################################
#  ∙ Accepts a numeric vector *Y* of antibody measurements (arbitrary scale)
#  ∙ Optionally a numeric covariate *Z* (same length) for RE estimation
#  ∙ Fits a 2‑component skew‑t mixture (mixsmsn + sn)
#  ∙ Provides soft posteriors p_soft, hard calls at mixture mid‑point
#  ∙ Computes Se, Sp, λ, λ², sample‑size inflation, optional RE
#  ∙ Draws histogram (two peaks) + posterior S‑curve
#  ∙ Returns a list with the fitted object, metrics, cutoff, and the plot
################################################################################

## ---- packages ----------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(mixsmsn, sn, ggplot2, dplyr, tibble, patchwork, scales)

stopifnot("smsn.mix" %in% ls("package:mixsmsn"),
          all(c("psn", "pst") %in% ls("package:sn")))

## ---- utilities ----------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b   # null‑coalesce

fit_skew_mix <- function(y, fam = c("Skew.t", "Skew.norm")) {
  fam <- match.arg(fam)
  mixsmsn::smsn.mix(
    y, g = 2,
    family   = fam,
    nu       = if (fam == "Skew.t") 4 else NULL,
    group    = FALSE,
    calc.im  = FALSE,
    obs.prob = TRUE                # ensure posterior probabilities are returned
  ) %>%
    { attr(., "family_used") <- fam; . }
}

cdf_skew <- function(q, mu, sigma2, shape, fam, nu = 4) {
  omega <- sqrt(sigma2)
  if (fam == "Skew.t") sn::pst(q, xi = mu, omega = omega, alpha = shape, nu = nu)
  else                  sn::psn(q, xi = mu, omega = omega, alpha = shape)
}

lambda_stats <- function(fit, pos, cut_log, nu_default = 4) {
  fam <- attr(fit, "family_used") %||% "Skew.t"
  neg <- setdiff(1:2, pos)
  get_nu <- function(k) if (!is.null(fit$nu) && !is.na(fit$nu[k]) && fit$nu[k] > 0) fit$nu[k] else nu_default
  nu_pos <- get_nu(pos); nu_neg <- get_nu(neg)
  Se <- 1 - cdf_skew(cut_log, fit$mu[pos], fit$sigma2[pos], fit$shape[pos], fam, nu_pos)
  Sp <-     cdf_skew(cut_log, fit$mu[neg], fit$sigma2[neg], fit$shape[neg], fam, nu_neg)
  λ  <- Se + Sp - 1
  tibble(Se, Sp, lambda = λ, lambda2 = λ^2, infl = 1/λ^2)
}

compute_stats <- function(df, λ) {
  if (!("Z" %in% names(df)))
    return(tibble(beta_soft = NA_real_, beta_hard = NA_real_, RE = NA_real_))
  df <- mutate(df, Zc = scale(Z, center = TRUE, scale = FALSE)[, 1])
  hard <- glm(S_hard ~ Zc, family = binomial(), data = df)
  soft <- lm (p_soft ~ Zc, data = df)
  var_hard    <- var(df$S_hard)
  sigma2_soft <- summary(soft)$sigma^2
  RE <- (1 / λ^2) * (var_hard / sigma2_soft)
  tibble(beta_soft = coef(soft)["Zc"],
         beta_hard = coef(hard)["Zc"],
         RE = RE)
}

diag_plot <- function(d, cut) {
  p1 <- ggplot(d, aes(Y)) +
    geom_histogram(aes(y = ..density..), bins = 60,
                   fill = "#99c", alpha = .4) +
    geom_vline(xintercept = cut, col = "red") +
    scale_x_continuous(labels = comma) + labs(y = "Density")

  p2 <- ggplot(d, aes(log(Y), p_soft)) +
    geom_point(alpha = .25, size = .8, col = "#006799") +
    geom_vline(xintercept = log(cut), col = "red") +
    scale_y_continuous(lims = c(0, 1)) + labs(y = "p_soft")

  p1 / p2
}

## ---- main user‑facing function ----------------------------------------------
#' Analyse antibody titres with a 2‑component skew‑t mixture
#' @param Y   numeric vector of strictly positive titres
#' @param Z   optional numeric covariate (same length). Needed for RE.
#' @param fam "Skew.t" (default) or "Skew.norm"
#' @param plot logical. Print diagnostic plot?
#' @return    list(fit, cutoff, metrics, data, plot)
#' @export
analyze_antibody <- function(Y, Z = NULL, fam = "Skew.t", plot = TRUE) {
  stopifnot(is.numeric(Y), all(is.finite(Y)), all(Y > 0))
  df <- tibble(Y = as.numeric(Y))
  if (!is.null(Z)) {
    stopifnot(is.numeric(Z), length(Z) == length(Y))
    df$Z <- as.numeric(Z)
  }

  fit <- fit_skew_mix(log(df$Y), fam = fam)
  pos <- which.max(fit$mu)
  post <- fit$obs.prob %||% fit$pro
  df$p_soft <- post[, pos]

  cut <- exp(weighted.mean(fit$mu, fit$pii))    # mixture mid‑point
  df$S_hard <- as.integer(df$Y > cut)

  lam <- lambda_stats(fit, pos, log(cut))
  stats <- bind_cols(lam, compute_stats(df, lam$lambda))

  plt <- diag_plot(df, cut) + plot_annotation(tag_levels = "A")
  if (isTRUE(plot)) print(plt)

  invisible(list(fit = fit, cutoff = cut, metrics = stats, data = df, plot = plt))
}

################################################################################
# Example (commented):
# library(readr); titres <- read_csv("my_titres.csv")$titre
# res <- analyze_antibody(titres)
# res$metrics
################################################################################

d <- rnorm(100000, mean = 100, sd = 10)
res <- analyze_antibody(d)
res$metrics