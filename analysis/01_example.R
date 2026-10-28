## Illustrative example (manuscript Section 3.6.1, Table 3.5). Run from the repository root.
source("R/nma_fe.R"); source("R/closed_form.R"); source("R/resampling.R")
set.seed(1)
d   <- read_brd("data/BRD_Arm_Network.csv")
old <- nma_fe(d)
cat(sprintf("Existing network, CEFTH vs TRIM: logOR = %.3f, SE = %.3f, z = %.3f, p = %.4f\n", old$est, old$se, old$z, old$p))
## new trial: CEFTH 10/25 vs TRIM 18/25 affected
d_all <- add_trial(d, "Ceftiofur hydrochloride", "Trimethoprim", 10, 25, 18, 25)
y_new <- log(.4/.6) - log(.72/.28); v_new <- 1/(25*.4*.6) + 1/(25*.72*.28)
cf_lt   <- nma_update_test(old$z, old$se^2, y_new, v_new, "lt")
cf_band <- nma_update_test(old$z, old$se^2, y_new, v_new, "band")
refit <- nma_fe(d_all)
stopifnot(abs(refit$z - cf_lt["t"]) < 1e-8)       # Theorem 1: closed form equals refitting the network
rs1 <- resample_method1(d_all, criterion = "lt", n_iter = 10000)
rs2 <- resample_method2(d_all, n_iter = 10000)
out <- data.frame(
  analysis = c("New trial alone", "Existing network plus new trial, ignoring selection",
               "Method 1, closed form, rule p < .1", "Method 1, resampling, rule p < .1",
               "Method 1, closed form, rule .05 < p < .1", "Method 2, closed form", "Method 2, resampling"),
  p_value = round(c(cf_lt["p_alone"], cf_lt["p_naive"], cf_lt["p_method1"], rs1["p"],
                    cf_band["p_method1"], cf_lt["p_method2"], rs2["p"]), 4))
print(out); cat("w =", round(cf_lt["w"], 4), " t =", round(cf_lt["t"], 4), " resamples kept (Method 1):", rs1["n_kept"], "\n")
write.csv(out, "results/table3_5_example.csv", row.names = FALSE)
