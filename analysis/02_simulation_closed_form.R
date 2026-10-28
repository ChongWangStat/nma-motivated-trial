## Simulation with the closed-form Methods 1 and 2 (manuscript Tables 3-4). Run from the repository root.
## About 3 minutes with the default numbers of simulated networks.
source("R/nma_fe.R"); source("R/closed_form.R")
set.seed(20261028)
N_H0 <- 100000; N_H1 <- 40000
d <- read_brd("data/BRD_Arm_Network.csv"); f <- nma_fe(d); risk <- nma_risks(d, f)
A <- "Ceftiofur hydrochloride"; B <- "Trimethoprim"
## design objects for fast refitting
arms <- do.call(rbind, lapply(seq_len(nrow(d)), function(i){
  k <- d$Number.of.arms[i]; tr <- c(d$Arm.1[i], d$Arm.2[i], d$Arm.3[i])[1:k]
  nn <- c(d$Total.number.in.arm.1[i], d$Total.number.in.arm.2[i], d$Total.number.in.arm.3[i])[1:k]
  data.frame(study = i, pos = 1:k, trt = tr, n = nn)}))
fit_X <- local({ # design matrix and contrast from the observed data
  main_trts <- setdiff(unique(c(d$Arm.1, d$Arm.2, d$Arm.3)), ""); base <- f$baseline; nb <- f$nb
  rows <- do.call(rbind, lapply(seq_len(nrow(d)), function(i){ k <- d$Number.of.arms[i]
    data.frame(study = i, num = if (k == 2) 1 else c(1, 3), ref = 2)}))
  X <- t(apply(rows, 1, function(r){ i <- r["study"]; x <- rep(0, length(nb))
    trt <- c(d$Arm.1[i], d$Arm.2[i], d$Arm.3[i])[r["num"]]
    if (d$Arm.2[i] == base) x[nb == trt] <- 1 else { x[nb == d$Arm.2[i]] <- -1; x[nb == trt] <- 1 }; x }))
  cvec <- rep(0, length(nb)); cvec[nb == A] <- 1; cvec[nb == B] <- -1
  idx <- function(st, pos) which(arms$study == st & arms$pos == pos)
  list(X = X, cvec = cvec, num = mapply(idx, rows$study, rows$num), ref = mapply(idx, rows$study, rows$ref),
       mrow = lapply(which(d$Number.of.arms == 3), function(i) which(rows$study == i)))
})
sim_old <- function(pA, pB){
  pr <- risk[arms$trt]; pr[arms$trt == A] <- pA; pr[arms$trt == B] <- pB
  r <- rbinom(nrow(arms), arms$n, pr); n <- arms$n
  r[r == 0] <- 0.5; r[r == n] <- n[r == n] - 0.5
  ph <- r/n; lo <- log(ph/(1 - ph)); iv <- 1/(n*ph*(1 - ph))
  Y <- lo[fit_X$num] - lo[fit_X$ref]; S <- diag(iv[fit_X$num] + iv[fit_X$ref])
  for (m in fit_X$mrow) S[m[1], m[2]] <- S[m[2], m[1]] <- iv[fit_X$ref[m[1]]]
  Si <- solve(S); V <- solve(t(fit_X$X) %*% Si %*% fit_X$X); mu <- V %*% (t(fit_X$X) %*% Si %*% Y)
  v <- drop(t(fit_X$cvec) %*% V %*% fit_X$cvec); c(z = sum(fit_X$cvec*mu)/sqrt(v), v = v)
}
sim_new <- function(n, pA, pB){
  r <- rbinom(2, n/2, c(pA, pB)); m <- rep(n/2, 2); r[r == 0] <- 0.5; r[r == m] <- m[r == m] - 0.5
  ph <- r/m; c(y = log(ph[1]/(1 - ph[1])) - log(ph[2]/(1 - ph[2])), v = sum(1/(m*ph*(1 - ph))))
}
run <- function(nsim, pA, pB, ns = c(30, 50, 100, 150, 200)){
  old <- t(replicate(nsim, sim_old(pA, pB))); pold <- 2*pnorm(-abs(old[, "z"]))
  out <- list()
  for (crit in c("lt", "band")){
    keep <- if (crit == "lt") which(pold < .1) else which(pold > .05 & pold < .1)
    for (n in ns){
      res <- t(sapply(keep, function(j){ nw <- sim_new(n, pA, pB)
        nma_update_test(unname(old[j, "z"]), unname(old[j, "v"]), unname(nw["y"]), unname(nw["v"]), crit)}))
      out[[paste(crit, n)]] <- data.frame(criterion = crit, n = n, selected = length(keep), mean_w = mean(res[, "w"]),
        method1 = mean(res[, "p_method1"] < .05), method2 = mean(res[, "p_method2"] < .05),
        alone = mean(res[, "p_alone"] < .05), naive = mean(res[, "p_naive"] < .05))
    }
  }
  do.call(rbind, out)
}
nt <- attr(risk, "n_by_trt"); p0 <- unname((risk[A]*nt[A] + risk[B]*nt[B])/(nt[A] + nt[B]))
H0 <- run(N_H0, p0, p0);   print(H0, digits = 3); write.csv(H0, "results/table3_type1_closed_form.csv", row.names = FALSE)
H1 <- run(N_H1, .35, .55); print(H1, digits = 3); write.csv(H1, "results/table4_power_closed_form.csv", row.names = FALSE)
