## Critical values, type I error of the unadjusted ("naive") analysis, conditional type I error of Method 1, and the design example
## (manuscript Section 2.2.4, Tables 6-7). Run from the repository root.
source("R/nma_fe.R"); source("R/closed_form.R")
ws <- c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95)
tab <- data.frame(w = ws,
  c_lt = sapply(ws, crit_value, criterion = "lt"), c_band = sapply(ws, crit_value, criterion = "band"),
  naive_lt = sapply(ws, naive_type1, criterion = "lt"), naive_band = sapply(ws, naive_type1, criterion = "band"))
print(round(tab, 3)); write.csv(round(tab, 4), "results/table6_critical_values.csv", row.names = FALSE)
cat("Limits as w -> 1: c_lt =", round(qnorm(1 - .0025), 3), " c_band =", round(qnorm(1 - .02625), 3), "\n")
o <- optimize(function(w) -naive_type1(w, "band"), c(.01, .99))
cat(sprintf("Maximum naive type I error, rule .05<p<.1: %.4f at w = %.3f\n", -o$objective, o$minimum))
## conditional type I error of Method 1 given |z_old| = u (rule p < .1)
for (w in c(.5, .72)){ ca <- crit_value(w, "lt")
  cat(sprintf("w = %.2f, c = %.3f:", w, ca)); for (u in c(1.70, 2, 2.5, 3)) cat(sprintf("  u=%.2f %.4f", u, g_tail(u, ca, w))); cat("\n") }
## design example: n per arm for 80% conditional power given the observed BRD network
d <- read_brd("data/BRD_Arm_Network.csv"); f <- nma_fe(d); risk <- nma_risks(d, f)
pB <- unname(risk["Trimethoprim"])
des <- do.call(rbind, lapply(c(0.52, 0.60), function(or){
  args <- list(d = log(or), z_old = f$z, v_old = f$se^2, pB = pB)
  data.frame(planning_OR = or,
             alone = do.call(n_for_power, c(list("alone"), args)),
             method1 = do.call(n_for_power, c(list("method1"), args)),
             method2 = do.call(n_for_power, c(list("method2"), args)))}))
print(des); write.csv(des, "results/table7_sample_size.csv", row.names = FALSE)
## approximation error of the one-sided formula for Method 2 power
for (w_n in c(25, 120, 400)) print(round(cond_power(w_n, log(.52), f$z, f$se^2, pB), 5))
