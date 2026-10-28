## Pairwise meta-analysis as a special case (manuscript Sections 2.4.4 and 3.3). Run from the repository root.
## The existing evidence is the only direct gamithromycin (GAMI) vs florfenicol (FLOR) trial in the BRD network
## (Van Donkersgoed et al. 2013; 50/297 vs 60/305 affected), i.e. a pairwise meta-analysis of a single trial,
## as in Simulation Study 2 of Ye et al. (2023). A new GAMI-FLOR trial is conducted only if the existing
## p-value meets the selection rule, and is then analysed with the closed-form Methods 1 and 2.
## About 1 minute.
source("R/closed_form.R")
set.seed(20261029)
nG <- 297; nF <- 305
lor <- function(r, n){                       # log odds ratio (arm 1 vs arm 2) and variance, 0.5 correction
  r[r == 0] <- 0.5; r[r == n] <- n[r == n] - 0.5; ph <- r/n
  c(y = log(ph[1]/(1 - ph[1])) - log(ph[2]/(1 - ph[2])), v = sum(1/(n*ph*(1 - ph))))
}
obs <- lor(c(50, 60), c(nG, nF))
cat(sprintf("Observed direct GAMI vs FLOR: log OR %.3f, SE %.3f, z %.2f, p %.3f\n",
            obs["y"], sqrt(obs["v"]), obs["y"]/sqrt(obs["v"]), 2*pnorm(-abs(obs["y"]/sqrt(obs["v"])))))
run <- function(nsim, pG, pF, ns = c(30, 50, 100, 150, 200)){
  old <- t(replicate(nsim, lor(rbinom(2, c(nG, nF), c(pG, pF)), c(nG, nF))))
  zold <- old[, "y"]/sqrt(old[, "v"]); pold <- 2*pnorm(-abs(zold))
  out <- list()
  for (crit in c("lt", "band")){
    keep <- if (crit == "lt") which(pold < .1) else which(pold > .05 & pold < .1)
    for (n in ns){
      res <- t(sapply(keep, function(j){
        nw <- lor(rbinom(2, n/2, c(pG, pF)), rep(n/2, 2))
        nma_update_test(zold[j], unname(old[j, "v"]), unname(nw["y"]), unname(nw["v"]), crit)}))
      out[[paste(crit, n)]] <- data.frame(criterion = crit, n = n, selected = length(keep), mean_w = mean(res[, "w"]),
        method1 = mean(res[, "p_method1"] < .05), method2 = mean(res[, "p_method2"] < .05),
        alone = mean(res[, "p_alone"] < .05), naive = mean(res[, "p_naive"] < .05))
    }
  }
  do.call(rbind, out)
}
p0 <- (50 + 60)/(nG + nF)                       # pooled risk under H0
H0 <- run(100000, p0, p0)
H1 <- run(40000, 0.12, 0.20)                    # alternative: GAMI 12% vs FLOR 20% (OR 0.55)
print(H0, digits = 3); print(H1, digits = 3)
write.csv(H0, "results/sim_pairwise_type1.csv", row.names = FALSE)
write.csv(H1, "results/sim_pairwise_power.csv", row.names = FALSE)
