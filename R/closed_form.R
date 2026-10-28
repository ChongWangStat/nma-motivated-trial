## Closed-form inference for a new A-B trial whose conduct was motivated by an existing NMA.
## Notation (manuscript Section 2.2.3): z_old, v_old = z statistic and variance of the existing-network
## A-B estimate; y_new, v_new = log odds ratio and variance of the new trial; w = v_new/(v_old + v_new).

z_sel <- function(alpha_sel) qnorm(1 - alpha_sel/2)

## Selection region R for z_old: "lt" -> existing p < alpha_sel; "band" -> band[1] < existing p < band[2]
selection_region <- function(criterion = c("lt", "band"), alpha_sel = 0.10, band = c(0.05, 0.10)){
  criterion <- match.arg(criterion)
  if (criterion == "lt") c(lo = z_sel(alpha_sel), hi = Inf) else c(lo = z_sel(band[2]), hi = z_sel(band[1]))
}

## g_c(u) = Pr(|sqrt(w) u + sqrt(1-w) V| >= c), V ~ N(delta_new, 1)
g_tail <- function(u, c, w, delta_new = 0){
  if (w >= 1) return(as.numeric(abs(u) >= c))
  s <- sqrt(1 - w); m <- sqrt(w)*u + s*delta_new
  (1 - pnorm((c - m)/s)) + pnorm((-c - m)/s)
}

## q(c; w, R) = Pr(|T| >= c | z_old in R), with z_old ~ N(delta_old, 1)
q_sel <- function(c, w, R, delta_old = 0, delta_new = 0){
  f <- function(u) g_tail(u, c, w, delta_new)*dnorm(u - delta_old)
  num <- integrate(f, R["lo"], R["hi"], subdivisions = 500L)$value +
         integrate(f, -R["hi"], -R["lo"], subdivisions = 500L)$value
  den <- (pnorm(R["hi"] - delta_old) - pnorm(R["lo"] - delta_old)) +
         (pnorm(-R["lo"] - delta_old) - pnorm(-R["hi"] - delta_old))
  unname(num/den)
}

## Main function: unadjusted ("naive"), alone, Method 1 (eq. 3) and Method 2 (eq. 5) p-values
nma_update_test <- function(z_old, v_old, y_new, v_new, criterion = "lt", alpha_sel = 0.10, band = c(0.05, 0.10)){
  w <- unname(v_new/(v_old + v_new))
  est <- (z_old*sqrt(v_old)/v_old + y_new/v_new)/(1/v_old + 1/v_new)
  t <- est*sqrt(1/v_old + 1/v_new)                 # = sqrt(w) z_old + sqrt(1-w) z_new  (Theorem 1)
  R <- selection_region(criterion, alpha_sel, band)
  c(w = w, t = t, p_naive = 2*pnorm(-abs(t)), p_alone = 2*pnorm(-abs(y_new/sqrt(v_new))),
    p_method1 = min(1, q_sel(abs(t), w, R)),
    p_method2 = g_tail(z_old, abs(t), w))
}

## Critical value c_alpha(w, R) of Method 1: reject iff |t| >= c_alpha (Corollary 1).
## Under the rule "p < alpha_sel" c_alpha is nondecreasing in w (Sidak 1968); see the Appendix of the paper.
crit_value <- function(w, criterion = "lt", alpha = 0.05, ...){
  R <- selection_region(criterion, ...)
  if (w <= 0) return(qnorm(1 - alpha/2))
  uniroot(function(c) q_sel(c, w, R) - alpha, c(0.5, 8), tol = 1e-10)$root
}

## Type I error of the unadjusted analysis (new trial added to the network, selection ignored;
## called "naive" in the code) under selection (Proposition 1); depends only on w and R
naive_type1 <- function(w, criterion = "lt", alpha = 0.05, ...)
  q_sel(qnorm(1 - alpha/2), w, selection_region(criterion, ...))

## Conditional critical value of Method 2 given z_old = u: |t| >= c_u  <=>  p_method2 <= alpha
crit_value_m2 <- function(u, w, alpha = 0.05)
  uniroot(function(c) g_tail(u, c, w) - alpha, c(0, abs(sqrt(w)*u) + 10*sqrt(1 - w) + 10), tol = 1e-10)$root

## Conditional power of a planned new trial given the observed z_old (Proposition 3).
## d = planning log odds ratio (A vs B), pB = risk under B, n = animals per arm.
## Method 1 uses the selection rule given by criterion/alpha_sel/band; Method 2 is computed exactly
## (method2_approx gives the one-sided approximation Phi(sign(z_old)*Delta - z_{1-alpha})).
cond_power <- function(n, d, z_old, v_old, pB, criterion = "lt", alpha = 0.05, alpha_sel = 0.10, band = c(0.05, 0.10)){
  pA <- plogis(qlogis(pB) + d)
  v_new <- 1/(n*pA*(1 - pA)) + 1/(n*pB*(1 - pB)); w <- v_new/(v_old + v_new); dn <- d/sqrt(v_new)
  ca <- crit_value(w, criterion, alpha, alpha_sel = alpha_sel, band = band)
  cu <- crit_value_m2(z_old, w, alpha)
  c(n = n, w = w,
    alone   = (1 - pnorm(qnorm(1 - alpha/2) - dn)) + pnorm(-qnorm(1 - alpha/2) - dn),
    method1 = g_tail(z_old, ca, w, dn),
    method2 = g_tail(z_old, cu, w, dn),
    method2_approx = pnorm(sign(z_old)*dn - qnorm(1 - alpha)),
    c_alpha = ca)
}

## Sample size: smallest n (per arm) beyond which conditional power stays >= target, searched on 5..n_max.
## Conditional power of Method 1 need not increase with n when |z_old| is large, so the equation
## power(n) = target can have several roots. Returns NA if the target is not reached by n_max.
n_for_power <- function(which, power = 0.8, n_max = 3000, ...){
  grid <- seq(5, n_max, by = 5)
  pw <- sapply(grid, function(n) cond_power(n, ...)[which])
  if (pw[length(pw)] < power) return(NA_integer_)
  below <- which(pw < power)
  if (length(below) == 0) return(5L)
  lo <- grid[max(below)]                       # last grid point below target; all later grid points meet it
  fine <- (lo + 1):(lo + 5)
  pf <- sapply(fine, function(n) cond_power(n, ...)[which])
  as.integer(fine[which(pf >= power)[1]])
}
