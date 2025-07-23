################################################################################
#   SOFT‑vs‑HARD SEROSTATUS POWER CHECKER  (CRAN‑safe 2025‑07‑21)
################################################################################

# ---- packages ----------------------------------------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(mixsmsn, sn, ggplot2, dplyr, tibble, patchwork, scales, purrr)

stopifnot("smsn.mix" %in% ls("package:mixsmsn"),
          all(c("psn","pst") %in% ls("package:sn")))

# ---- utilities ----------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b   # null‑coalesce (defined early)

set.seed(1234)

# ---- helper functions ---------------------------------------------------------
fit_skew_mix <- function(y, fam = c("Skew.t", "Skew.norm")) {
  fam <- match.arg(fam)
  fit <- mixsmsn::smsn.mix(
    y, g = 2,
    family   = fam,
    nu       = if (fam == "Skew.t") 4 else NULL,
    group    = FALSE,
    calc.im  = FALSE,
    obs.prob = TRUE       # ensure posterior probabilities are returned
  )
  attr(fit, "family_used") <- fam   # stash for later
  fit
}

cdf_skew <- function(q, mu, sigma2, shape, fam, nu = 4) {
  omega <- sqrt(sigma2)
  if (fam == "Skew.t")
    sn::pst(q, xi = mu, omega = omega, alpha = shape, nu = nu)
  else
    sn::psn(q, xi = mu, omega = omega, alpha = shape)
}

lambda_stats <- function(fit, pos, cut_log, nu_default = 4) {
  fam <- attr(fit, "family_used") %||% "Skew.t"
  neg <- setdiff(1:2, pos)
  # robust d.f. extraction
  get_nu <- function(k) {
    if (!is.null(fit$nu) && !is.na(fit$nu[k]) && fit$nu[k] > 0)
      fit$nu[k]
    else
      nu_default
  }
  nu_pos <- get_nu(pos)
  nu_neg <- get_nu(neg)

  Se <- 1 - cdf_skew(cut_log, fit$mu[pos], fit$sigma2[pos],
                     fit$shape[pos], fam, nu_pos)
  Sp <-     cdf_skew(cut_log, fit$mu[neg], fit$sigma2[neg],
                     fit$shape[neg], fam, nu_neg)
  λ  <- Se + Sp - 1
  tibble(Se, Sp, lambda = λ,
         lambda2 = λ^2, infl = 1/λ^2)
}

soft_hard_RE <- function(df, λ) {
  df <- mutate(df, Zc = scale(Z, center = TRUE, scale = FALSE)[, 1])
  hard <- glm(S_hard ~ Zc, family = binomial, data = df)
  soft <- lm (p_soft ~ Zc, data = df)
  var_hard    <- var(df$S_hard)
  sigma2_soft <- summary(soft)$sigma^2
  tibble(beta_soft = coef(soft)["Zc"],
         beta_hard = coef(hard)["Zc"],
         RE = (1/λ^2) * (var_hard / sigma2_soft))
}

sim_titre <- function(n, prev, d, sigma, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  true_inf <- rbinom(n, 1, prev)
  y <- rnorm(n, 0, sigma)
  y[true_inf == 1] <- rnorm(sum(true_inf == 1), d, sigma)
  tibble(Y = exp(y), Z = rnorm(n))
}

diag_plot <- function(d, cut) {
  p1 <- ggplot(d, aes(Y)) +
    geom_histogram(aes(y = ..density..), bins = 60,
                   fill = "#99c", alpha = .4) +
    geom_vline(xintercept = cut, col = "red") +
    scale_x_log10(labels = comma) +
    labs(y = "Density")

  p2 <- ggplot(d, aes(log(Y), p_soft)) +
    geom_point(alpha = .25, size = .8, col = "#006799") +
    geom_vline(xintercept = log(cut), col = "red") +
    labs(y = "p_soft")

  p1 / p2
}

# ---- simulation archetypes ----------------------------------------------------
arch <- tribble(
  ~name, ~prev, ~d, ~sigma,
  "HighSep",   0.95, 2.5, 1.0,
  "HSV1_like", 0.70, 1.2, 1.0,
  "HSV2_like", 0.16, 1.0, 1.0
)

res   <- list()
plots <- list()

for (i in seq_len(nrow(arch))) {
  cfg <- arch[i, ]
  d   <- sim_titre(8000, cfg$prev, cfg$d, cfg$sigma, seed = 100 + i)

  fit <- fit_skew_mix(log(d$Y), fam = "Skew.t")
  pos <- which.max(fit$mu)
  post <- fit$obs.prob %||% fit$pro
  d$p_soft <- post[, pos]

  # mixture cut‑off uses component weights
  cut <- exp(weighted.mean(fit$mu, fit$pii))
  d$S_hard <- as.integer(d$Y > cut)

  met <- lambda_stats(fit, pos, log(cut))
  cmp <- soft_hard_RE(d, met$lambda)

  res[[i]]   <- bind_cols(cfg, met, cmp)
  plots[[i]] <- diag_plot(d, cut) + plot_annotation(
    title = sprintf("%s  (prev=%.2f, d/σ=%.2f)", cfg$name, cfg$prev, cfg$d/cfg$sigma))
}

summary_tbl <- bind_rows(res) %>%
  select(name, prev, d, lambda, infl, RE, beta_soft, beta_hard)
print(summary_tbl)

print(wrap_plots(plots) + plot_annotation(tag_levels = 'A'))
