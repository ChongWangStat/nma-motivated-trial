## Cross-check: corrected resampling Methods 1 and 2 versus the closed forms on simulated selected data sets
## (selection rule p < .1, null hypothesis). About 9 minutes with the defaults. Run from the repository root.
source("R/nma_fe.R"); source("R/closed_form.R"); source("R/resampling.R")
set.seed(42)
N_SETS <- 30; N_ITER <- 2000; n_new <- 50
d <- read_brd("data/BRD_Arm_Network.csv"); f <- nma_fe(d); risk <- nma_risks(d, f)
A <- "Ceftiofur hydrochloride"; B <- "Trimethoprim"; nt <- attr(risk, "n_by_trt")
p0 <- unname((risk[A]*nt[A] + risk[B]*nt[B])/(nt[A] + nt[B])); r0 <- risk; r0[c(A, B)] <- p0
res <- NULL
while (is.null(res) || nrow(res) < N_SETS){
  dd <- resample_events(d, r0); old <- nma_fe(dd)
  if (old$p >= .1) next
  r <- rbinom(2, n_new/2, p0); r[r == 0] <- .5; r[r == n_new/2] <- n_new/2 - .5
  da <- add_trial(dd, A, B, r[1], n_new/2, r[2], n_new/2)
  ph <- r/(n_new/2); y <- log(ph[1]/(1 - ph[1])) - log(ph[2]/(1 - ph[2])); v <- sum(1/(n_new/2*ph*(1 - ph)))
  cf <- nma_update_test(old$z, old$se^2, y, v, "lt")
  res <- rbind(res, c(cf, rs1 = unname(resample_method1(da, n_iter = N_ITER)["p"]),
                      rs2 = unname(resample_method2(da, n_iter = N_ITER)["p"])))
}
res <- as.data.frame(res)
cat(sprintf("Method 1: correlation %.3f, mean |difference| %.3f; Method 2: correlation %.3f, mean |difference| %.3f\n",
  cor(res$p_method1, res$rs1), mean(abs(res$p_method1 - res$rs1)), cor(res$p_method2, res$rs2), mean(abs(res$p_method2 - res$rs2))))
write.csv(res, "results/check_resampling_vs_closed_form.csv", row.names = FALSE)
