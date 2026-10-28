## Resampling versions of Methods 1 and 2 (manuscript Sections 2.2-2.3), as described in the text.
## Correction relative to github.com/lamckeen/Error-Control-NMA: in Method 1 the resampled data are kept
## when the A-B p-value of the EXISTING network (new trial removed) meets the selection rule; the original
## scripts applied the rule to the network that includes the new trial.

resample_events <- function(d, risk){
  for (k in 1:3){
    tr <- d[[paste0("Arm.", k)]]; nn <- d[[paste0("Total.number.in.arm.", k)]]
    ok <- which(tr != "" & !is.na(nn)); r <- d[[paste0("Number.of.Event.in.arm.", k)]]
    r[ok] <- rbinom(length(ok), nn[ok], risk[tr[ok]])
    r[ok][r[ok] == 0] <- 0.5; full <- ok[r[ok] == nn[ok]]; r[full] <- nn[full] - 0.5
    d[[paste0("Number.of.Event.in.arm.", k)]] <- r
  }
  d
}

## Null risks: NMA-implied risks from the network including the new trial, with A and B set to their
## subject-weighted average (p_A = p_B under H0).
null_risks <- function(d_all, fit_all, A, B){
  risk <- nma_risks(d_all, fit_all); nt <- attr(risk, "n_by_trt")
  pool <- (risk[A]*nt[A] + risk[B]*nt[B])/(nt[A] + nt[B])
  risk[c(A, B)] <- pool
  risk
}

resample_method1 <- function(d_all, new_id = 99, A = "Ceftiofur hydrochloride", B = "Trimethoprim",
                             criterion = "lt", alpha_sel = 0.10, band = c(0.05, 0.10), n_iter = 10000){
  fit <- nma_fe(d_all, A, B); risk <- null_risks(d_all, fit, A, B)
  keep_rule <- if (criterion == "lt") function(p) p < alpha_sel else function(p) p > band[1] & p < band[2]
  tstar <- numeric(n_iter); keep <- logical(n_iter)
  for (b in seq_len(n_iter)){
    db <- resample_events(d_all, risk)
    tstar[b] <- nma_fe(db, A, B)$z
    keep[b]  <- keep_rule(nma_fe(db[db$Study.number != new_id, ], A, B)$p)   # existing network only
  }
  c(p = mean(abs(tstar[keep]) >= abs(fit$z)), n_kept = sum(keep))
}

resample_method2 <- function(d_all, new_id = 99, A = "Ceftiofur hydrochloride", B = "Trimethoprim", n_iter = 10000){
  fit <- nma_fe(d_all, A, B); risk <- null_risks(d_all, fit, A, B)
  old <- nma_fe(d_all[d_all$Study.number != new_id, ], A, B)
  nw <- d_all[d_all$Study.number == new_id, ]
  nA <- nw$Total.number.in.arm.1; nB <- nw$Total.number.in.arm.2
  tstar <- replicate(n_iter, {
    r <- rbinom(2, c(nA, nB), risk[c(A, B)]); n <- c(nA, nB)
    r[r == 0] <- 0.5; r[r == n] <- n[r == n] - 0.5; ph <- r/n
    y <- log(ph[1]/(1 - ph[1])) - log(ph[2]/(1 - ph[2])); v <- sum(1/(n*ph*(1 - ph)))
    est <- (old$est/old$se^2 + y/v)/(1/old$se^2 + 1/v); est*sqrt(1/old$se^2 + 1/v)   # Theorem 1
  })
  c(p = mean(abs(tstar) >= abs(fit$z)))
}
