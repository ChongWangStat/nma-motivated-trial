# Selection-adjusted confidence intervals under the second-stage normal model.
# Fixed-effect A-B contrast, independent old and new estimates, and fixed known
# (or stage-2 conditioned-on) variances. No external packages required.
# Matches equation (selectiveci) in the manuscript.

method1_cdf <- function(eta, t, w, criterion = c("lt", "band"),
                        alpha_sel = 0.10, band = c(0.05, 0.10)) {
  criterion <- match.arg(criterion)
  stopifnot(w > 0, w < 1, is.finite(t))
  if (criterion == "lt") {
    lower <- qnorm(1 - alpha_sel / 2); upper <- Inf
  } else {
    lower <- qnorm(1 - band[2] / 2)
    upper <- qnorm(1 - band[1] / 2)
  }
  mu_old <- sqrt(w) * eta
  s_new <- sqrt(1 - w)
  f <- function(u) {
    pnorm((t - sqrt(w) * u - (1 - w) * eta) / s_new) *
      dnorm(u, mean = mu_old, sd = 1)
  }
  num <- integrate(f, lower, upper, subdivisions = 500L)$value +
         integrate(f, -upper, -lower, subdivisions = 500L)$value
  den <- (pnorm(upper - mu_old) - pnorm(lower - mu_old)) +
         (pnorm(-lower - mu_old) - pnorm(-upper - mu_old))
  if (!is.finite(den) || den <= 0) stop("Selection probability too small for stable inversion")
  num / den
}

method1_ci <- function(z_old, v_old, y_new, v_new,
                       level = 0.95, criterion = "lt", ...) {
  stopifnot(v_old > 0, v_new > 0, level > 0, level < 1)
  w <- v_new / (v_old + v_new)
  vstar <- 1 / (1/v_old + 1/v_new)
  znew <- y_new / sqrt(v_new)
  t <- sqrt(w) * z_old + sqrt(1 - w) * znew
  alpha <- 1 - level
  cdf <- function(eta) method1_cdf(eta, t, w, criterion, ...)
  # The CDF is decreasing in eta; roots here give the equal-tail interval.
  # These brackets cover conventional z-statistics; enlarge if necessary.
  lo <- uniroot(function(eta) cdf(eta) - (1 - alpha / 2), c(-12, 12))$root
  hi <- uniroot(function(eta) cdf(eta) - alpha / 2, c(-12, 12))$root
  lor <- sqrt(vstar) * c(lower = lo, upper = hi)
  list(log_odds_ratio = lor, odds_ratio = exp(lor),
       t = t, w = w, level = level, criterion = criterion)
}

# Example using rounded manuscript statistics (Table 7):
# ynew <- log(.4/.6) - log(.72/.28)
# vnew <- 1/(25*.4*.6) + 1/(25*.72*.28)
# method1_ci(-1.77, 0.140, ynew, vnew, criterion = "lt")
